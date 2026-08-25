#!/usr/bin/env python
"""PopV conformance oracle: direct upstream API call, no Shennong wrappers.

Reads the same MatrixMarket inputs the candidate exports and invokes
popv.preprocessing.Process_Query + popv.annotation.annotate_data directly.
Writes predictions.csv in the same layout as the Shennong runner so a plain
frame comparison proves scientific fidelity of the wrapper path.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

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

    seed = int(config["seed"])
    random.seed(seed)
    np.random.seed(seed)
    import torch

    torch.manual_seed(seed)
    import scvi

    scvi.settings.seed = seed

    import popv

    query_adata = _read_adata(input_dir / "query")
    ref_adata = _read_adata(input_dir / "reference")
    label_key = config["label_key"]
    ref_adata.obs[label_key] = ref_adata.obs[label_key].astype(str)

    common_genes = np.intersect1d(
        ref_adata.var_names.astype(str), query_adata.var_names.astype(str)
    )
    query_adata = query_adata[:, common_genes].copy()
    ref_adata = ref_adata[:, common_genes].copy()

    hvg_config = config.get("hvg")
    query_batch_kwargs = (
        {} if config.get("query_batch_key") is None
        else {"query_batch_key": config["query_batch_key"]}
    )
    processed = popv.preprocessing.Process_Query(
        query_adata,
        ref_adata,
        ref_labels_key=label_key,
        ref_batch_key=config.get("ref_batch_key"),
        cl_obo_folder=False,
        prediction_mode="retrain",
        unknown_celltype_label="unknown",
        n_samples_per_label=int(config.get("n_samples_per_label", 300)),
        hvg=int(hvg_config) if hvg_config not in (None, False) else None,
        save_path_trained_models=str(output_dir / "trained_models"),
        **query_batch_kwargs,
    )
    methods = config.get("methods")
    popv.annotation.annotate_data(
        processed.adata,
        methods=list(methods) if methods else None,
        save_path=str(output_dir),
    )
    print("Oracle PopV run completed.")


if __name__ == "__main__":
    main()
