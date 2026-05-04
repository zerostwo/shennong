#!/usr/bin/env python
"""Run scArches/scPoli-style object workflows for Shennong."""

from __future__ import annotations

import argparse
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
    return adata


def _write_embedding(adata: ad.AnnData, key: str, path: Path) -> None:
    arr = np.asarray(adata.obsm[key])
    pd.DataFrame(arr, index=adata.obs_names).to_csv(path)


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

    import scanpy as sc
    import scarches  # noqa: F401

    adata = _read_adata(input_dir / "query")
    batch_key = config.get("batch_key")
    labels_key = config.get("labels_key")
    n_pcs = int(config.get("n_pcs", 30))

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars), batch_key=batch_key)
    if "highly_variable" in adata.var:
        adata = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(n_pcs, max(1, adata.n_obs - 1), max(1, adata.n_vars - 1)))
    latent_key = config.get("latent_key", "X_scarches")
    adata.obsm[latent_key] = adata.obsm["X_pca"]
    if labels_key and labels_key in adata.obs:
        adata.obs[f"{config.get('method', 'scarches')}_label"] = adata.obs[labels_key].astype(str)

    adata.obs.to_csv(output_dir / "obs.csv")
    _write_embedding(adata, latent_key, output_dir / "latent.csv")
    if config.get("write_h5ad", True):
        adata.write_h5ad(output_dir / "scarches.h5ad")
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "method": config.get("method", "scarches"),
                "latent_path": str(output_dir / "latent.csv"),
                "obs_path": str(output_dir / "obs.csv"),
                "output_h5ad": str(output_dir / "scarches.h5ad"),
                "note": "Default Shennong runner exports a PCA latent after validating the scarches environment; pass a custom script for a trained reference mapping workflow.",
            },
            handle,
            indent=2,
        )


if __name__ == "__main__":
    main()
