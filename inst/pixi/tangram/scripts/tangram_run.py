#!/usr/bin/env python
"""Run Tangram mapping from Shennong-exported Seurat objects."""

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
    tg.pp_adatas(reference, spatial, genes=genes)
    ad_map = tg.map_cells_to_space(reference, spatial, mode=config.get("mode", "cells"))
    pd.DataFrame(np.asarray(ad_map.X), index=reference.obs_names, columns=spatial.obs_names).to_csv(output_dir / "mapping.csv")

    cell_type_key = config.get("cell_type_key")
    if cell_type_key and cell_type_key in reference.obs:
        annotations = pd.get_dummies(reference.obs[cell_type_key])
        projected = pd.DataFrame(np.asarray(ad_map.X).T @ annotations.to_numpy(), index=spatial.obs_names, columns=annotations.columns)
        projected.to_csv(output_dir / "obs.csv")
    else:
        spatial.obs.to_csv(output_dir / "obs.csv")
    ad_map.write_h5ad(output_dir / "tangram_map.h5ad")
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump({"method": "tangram", "mapping_path": str(output_dir / "mapping.csv"), "output_h5ad": str(output_dir / "tangram_map.h5ad")}, handle, indent=2)


if __name__ == "__main__":
    main()
