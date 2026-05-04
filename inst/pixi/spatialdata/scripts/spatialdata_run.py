#!/usr/bin/env python
"""Create a SpatialData table from a Shennong-exported Seurat object."""

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

    from spatialdata import SpatialData

    query_dir = input_dir / "query"
    matrix = sparse.csr_matrix(io.mmread(query_dir / "matrix.mtx")).transpose().tocsr()
    obs = pd.read_csv(query_dir / "obs.csv", index_col=0)
    var = pd.read_csv(query_dir / "var.csv", index_col=0)
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    if "feature_id" in adata.var.columns:
        adata.var_names = adata.var["feature_id"].astype(str).to_numpy()
    spatial_path = query_dir / "spatial.csv"
    if spatial_path.exists():
        adata.obsm["spatial"] = pd.read_csv(spatial_path, index_col=0).loc[adata.obs_names].to_numpy()

    sdata = SpatialData(tables={config.get("table_name", "table"): adata})
    zarr_path = output_dir / config.get("zarr_name", "spatialdata.zarr")
    sdata.write(zarr_path)
    adata.obs.to_csv(output_dir / "obs.csv")
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump({"method": "spatialdata", "zarr_path": str(zarr_path)}, handle, indent=2)


if __name__ == "__main__":
    main()
