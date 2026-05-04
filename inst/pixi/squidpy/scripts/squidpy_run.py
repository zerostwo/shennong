#!/usr/bin/env python
"""Run Squidpy spatial graph workflows from a Shennong-exported Seurat object."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
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
    spatial = pd.read_csv(path / "spatial.csv", index_col=0)
    adata.obsm["spatial"] = spatial.loc[adata.obs_names].to_numpy()
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

    import squidpy as sq

    adata = _read_adata(input_dir / "query")
    sq.gr.spatial_neighbors(adata, coord_type=config.get("coord_type", "generic"))
    cluster_key = config.get("cluster_key")
    if cluster_key and cluster_key in adata.obs:
        sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
    adata.obs.to_csv(output_dir / "obs.csv")
    adata.write_h5ad(output_dir / "squidpy.h5ad")
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump({"method": "squidpy", "output_h5ad": str(output_dir / "squidpy.h5ad")}, handle, indent=2)


if __name__ == "__main__":
    main()
