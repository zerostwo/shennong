#!/usr/bin/env python
"""Run Squidpy spatial graph workflows from a Shennong-exported Seurat object."""

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
    connectivities = adata.obsp["spatial_connectivities"].tocoo()
    graph_path = output_dir / "spatial_graph.csv"
    pd.DataFrame(
        {
            "source": adata.obs_names.to_numpy()[connectivities.row],
            "target": adata.obs_names.to_numpy()[connectivities.col],
            "weight": connectivities.data,
        }
    ).to_csv(graph_path, index=False)

    cluster_key = config.get("cluster_key")
    enrichment_path = None
    if cluster_key and cluster_key in adata.obs:
        adata.obs[cluster_key] = adata.obs[cluster_key].astype(str).astype("category")
        sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
        result = adata.uns.get(f"{cluster_key}_nhood_enrichment", {})
        zscore = result.get("zscore")
        if zscore is not None:
            categories = adata.obs[cluster_key].cat.categories.astype(str).to_numpy()
            row, col = np.indices(np.asarray(zscore).shape)
            enrichment = pd.DataFrame(
                {
                    "group_1": categories[row.ravel()],
                    "group_2": categories[col.ravel()],
                    "zscore": np.asarray(zscore).ravel(),
                }
            )
            counts = result.get("count")
            if counts is not None and np.asarray(counts).shape == np.asarray(zscore).shape:
                enrichment["count"] = np.asarray(counts).ravel()
            enrichment_path = output_dir / "neighborhood_enrichment.csv"
            enrichment.to_csv(enrichment_path, index=False)

    h5ad_path = output_dir / "squidpy.h5ad"
    if bool(config.get("write_h5ad", False)):
        adata.write_h5ad(h5ad_path)
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "method": "squidpy",
                "spatial_graph_path": str(graph_path),
                "neighborhood_enrichment_path": str(enrichment_path) if enrichment_path else None,
                "output_h5ad": str(h5ad_path) if h5ad_path.exists() else None,
                "n_cells": int(adata.n_obs),
                "n_features": int(adata.n_vars),
                "squidpy_version": sq.__version__,
            },
            handle,
            indent=2,
        )


if __name__ == "__main__":
    main()
