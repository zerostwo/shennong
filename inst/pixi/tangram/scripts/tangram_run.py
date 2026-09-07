#!/usr/bin/env python
"""Run Tangram mapping from Shennong-exported Seurat objects."""

from __future__ import annotations

import argparse
import importlib.metadata as metadata
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _read_adata(path: Path) -> ad.AnnData:
    matrix = sparse.csr_matrix(io.mmread(path / "matrix.mtx")).transpose().tocsr()
    obs = pd.read_csv(path / "obs.csv", index_col=0)
    var = pd.read_csv(path / "var.csv", index_col=0)
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    if "feature_id" in adata.var.columns:
        adata.var_names = adata.var["feature_id"].astype(str).to_numpy()
    spatial_path = path / "spatial.csv"
    if spatial_path.exists():
        adata.obsm["spatial"] = pd.read_csv(spatial_path, index_col=0).loc[adata.obs_names].to_numpy()
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

    import tangram as tg

    spatial = _read_adata(input_dir / "query")
    reference = _read_adata(input_dir / "reference")
    genes = list(reference.var_names.intersection(spatial.var_names))
    if len(genes) == 0:
        raise ValueError("No shared genes between reference and spatial objects for Tangram.")
    max_mapping_gb = float(config.get("max_mapping_gb", 2.0))
    estimated_mapping_gb = reference.n_obs * spatial.n_obs * 4 / (1024**3)
    if not np.isfinite(max_mapping_gb) or max_mapping_gb <= 0:
        raise ValueError("max_mapping_gb must be one positive finite number.")
    if estimated_mapping_gb > max_mapping_gb:
        raise MemoryError(
            f"Tangram mapping is estimated at {estimated_mapping_gb:.2f} GiB "
            f"before working copies, exceeding max_mapping_gb={max_mapping_gb}."
        )
    tg.pp_adatas(reference, spatial, genes=genes)
    ad_map = tg.map_cells_to_space(reference, spatial, mode=config.get("mode", "cells"))
    mapping_path = output_dir / "mapping.csv"
    mapping = ad_map.X.toarray() if sparse.issparse(ad_map.X) else np.asarray(ad_map.X)
    mapping = np.asarray(mapping, dtype=np.float32)
    if mapping.shape != (reference.n_obs, spatial.n_obs):
        raise ValueError("Tangram mapping dimensions do not match reference and spatial cell identities.")
    if not np.isfinite(mapping).all() or np.any(mapping < 0):
        raise ValueError("Tangram mapping must contain finite, non-negative probabilities.")
    if not np.allclose(mapping.sum(axis=1), 1.0, rtol=1e-4, atol=1e-6):
        raise ValueError("Each Tangram reference-cell mapping row must sum to 1 within tolerance.")
    pd.DataFrame(mapping, index=reference.obs_names, columns=spatial.obs_names).to_csv(mapping_path)

    cell_type_key = config.get("cell_type_key")
    obs_path = None
    if cell_type_key and cell_type_key in reference.obs:
        annotations = pd.get_dummies(reference.obs[cell_type_key])
        projected = pd.DataFrame(mapping.T @ annotations.to_numpy(), index=spatial.obs_names, columns=annotations.columns)
        obs_path = output_dir / "obs.csv"
        projected.to_csv(obs_path)
    h5ad_path = output_dir / "tangram_map.h5ad"
    if bool(config.get("write_h5ad", False)):
        ad_map.write_h5ad(h5ad_path)
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "method": "tangram",
                "mapping_path": str(mapping_path),
                "projected_obs_path": str(obs_path) if obs_path else None,
                "output_h5ad": str(h5ad_path) if h5ad_path.exists() else None,
                "n_cells": int(spatial.n_obs),
                "n_reference_cells": int(reference.n_obs),
                "n_shared_genes": len(genes),
                "mapping_dtype": str(mapping.dtype),
                "estimated_mapping_gb": estimated_mapping_gb,
                "max_mapping_gb": max_mapping_gb,
                "tangram_version": metadata.version("tangram-sc"),
            },
            handle,
            indent=2,
        )


if __name__ == "__main__":
    main()
