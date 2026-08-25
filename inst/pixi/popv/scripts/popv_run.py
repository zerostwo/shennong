#!/usr/bin/env python
"""Run PopV consensus cell-type annotation for Shennong."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
import sys


import numpy as np
import pandas as pd
from scipy import io, sparse


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _read_adata(path: Path):
    import anndata as ad

    matrix = sparse.csr_matrix(io.mmread(path / "matrix.mtx")).transpose().tocsr()
    obs = pd.read_csv(path / "obs.csv", index_col=0)
    var = pd.read_csv(path / "var.csv", index_col=0)
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    if "feature_id" in adata.var.columns:
        adata.var_names = adata.var["feature_id"].astype(str).to_numpy()
    return adata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    config = _read_json(Path(args.config))

    import anndata as ad
    import scanpy as sc

    seed = int(config.get("seed", 0))
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
    except Exception:  # pragma: no cover - torch is a hard popv dependency
        pass
    try:
        import scvi

        scvi.settings.seed = seed
    except Exception:
        pass

    query_adata = _read_adata(input_dir / "query")
    ref_adata = _read_adata(input_dir / "reference")
    for name, obj in (("query", query_adata), ("reference", ref_adata)):
        if not sparse.issparse(obj.X):
            raise MemoryError(f"The {name} expression matrix must remain sparse.")
        if obj.X.min() < 0:
            raise ValueError(
                f"The {name} matrix must contain raw, non-negative counts; "
                "PopV normalizes internally."
            )

    label_key = config["label_key"]
    if label_key not in ref_adata.obs.columns:
        raise KeyError(
            f"Reference metadata column '{label_key}' is required by PopV."
        )
    ref_adata.obs[label_key] = ref_adata.obs[label_key].astype(str)

    common_genes = np.intersect1d(
        ref_adata.var_names.astype(str), query_adata.var_names.astype(str)
    )
    min_genes = int(config.get("min_shared_genes", 2))
    if len(common_genes) < min_genes:
        raise ValueError(
            f"Query and reference share only {len(common_genes)} features; "
            f"PopV requires at least {min_genes}."
        )
    query_adata = query_adata[:, common_genes].copy()
    ref_adata = ref_adata[:, common_genes].copy()

    import popv

    save_path_trained_models = str(Path(output_dir) / "trained_models")
    hvg_config = config.get("hvg")
    process_query_kwargs = dict(
        ref_labels_key=label_key,
        # ref_batch_key and cl_obo_folder are positional-required upstream.
        ref_batch_key=config.get("ref_batch_key"),
        cl_obo_folder=config.get("cl_obo_folder", False),
        query_batch_key=config.get("query_batch_key"),
        prediction_mode=config.get("prediction_mode", "retrain"),
        unknown_celltype_label=config.get("unknown_celltype_label", "unknown"),
        n_samples_per_label=int(config.get("n_samples_per_label", 300)),
        relabel_reference_cells=bool(config.get("relabel_reference_cells", False)),
        hvg=int(hvg_config) if hvg_config not in (None, False) else None,
        save_path_trained_models=save_path_trained_models,
    )
    if process_query_kwargs["query_batch_key"] is None:
        del process_query_kwargs["query_batch_key"]
    processed = popv.preprocessing.Process_Query(query_adata, ref_adata, **process_query_kwargs)

    methods = config.get("methods")
    popv.annotation.annotate_data(
        processed.adata,
        methods=list(methods) if methods else None,
        save_path=str(output_dir),
    )

    manifest = {
        "method": "popv",
        "popv_version": getattr(popv, "__version__", "unknown"),
        "prediction_mode": process_query_kwargs["prediction_mode"],
        "methods_requested": methods,
        "ontology_enabled": bool(process_query_kwargs["cl_obo_folder"]),
        "seed": seed,
        "n_query_cells_input": int(len(pd.read_csv(input_dir / "query" / "obs.csv", index_col=0))),
        "predictions_path": str(output_dir / "predictions.csv"),
    }
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    print("PopV annotation completed.")


if __name__ == "__main__":
    main()
