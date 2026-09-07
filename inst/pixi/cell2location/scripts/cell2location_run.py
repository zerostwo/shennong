#!/usr/bin/env python
"""Run cell2location from Shennong-exported spatial Seurat data."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys


import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _read_unique_id_table(path: Path, id_column: str) -> pd.DataFrame:
    table = pd.read_csv(path, dtype=str)
    if id_column not in table.columns:
        raise ValueError(f"{path.name} must contain a {id_column} column.")
    raw_ids = table.pop(id_column)
    if raw_ids.isna().any():
        raise ValueError(f"{path.name} {id_column} values must be unique and non-empty.")
    identifiers = raw_ids.astype(str).str.strip()
    if identifiers.eq("").any() or identifiers.duplicated().any():
        raise ValueError(f"{path.name} {id_column} values must be unique and non-empty.")
    table.index = identifiers.to_numpy()
    return table


def _validate_raw_counts(matrix: sparse.csr_matrix) -> None:
    values = matrix.data
    if values.size and (not np.isfinite(values).all() or np.min(values) < 0):
        raise ValueError("cell2location input must contain finite, non-negative raw counts.")
    if values.size and np.any(np.abs(values - np.rint(values)) > 1e-8):
        raise ValueError(
            "cell2location input must contain integer-like raw counts; values are never rounded silently."
        )


def _read_adata(path: Path) -> ad.AnnData:
    matrix = sparse.csr_matrix(io.mmread(path / "matrix.mtx")).transpose().tocsr()
    obs = _read_unique_id_table(path / "obs.csv", "cell_id")
    var = _read_unique_id_table(path / "var.csv", "feature_id")
    if matrix.shape != (obs.shape[0], var.shape[0]):
        raise ValueError("Exported expression matrix dimensions do not match obs.csv and var.csv.")
    _validate_raw_counts(matrix)
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.obs_names = obs.index.astype(str)
    adata.var_names = var.index.astype(str)
    return adata


def _read_reference_signatures(path: Path) -> pd.DataFrame:
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle), [])
    state_header = [value.strip() for value in header[1:]]
    if (
        len(header) < 2
        or any(not value for value in state_header)
        or len(state_header) != len(set(state_header))
    ):
        raise ValueError(
            "reference_signatures cell-state identifiers must be unique and non-empty."
        )
    raw = pd.read_csv(path)
    if raw.shape[0] < 1 or raw.shape[1] < 2:
        raise ValueError(
            "reference_signatures must contain feature IDs and at least one cell-state column."
        )
    raw_ids = raw.iloc[:, 0]
    if raw_ids.isna().any():
        raise ValueError("reference_signatures feature identifiers must be unique and non-empty.")
    feature_ids = raw_ids.astype(str).str.strip()
    state_ids = pd.Index(state_header)
    if (
        feature_ids.eq("").any()
        or feature_ids.duplicated().any()
        or (state_ids.str.len() == 0).any()
    ):
        raise ValueError(
            "reference_signatures feature and cell-state identifiers must be unique and non-empty."
        )
    try:
        values = raw.iloc[:, 1:].to_numpy(dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError("reference_signatures values must be numeric.") from error
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("reference_signatures values must be finite and non-negative.")
    return pd.DataFrame(values, index=feature_ids.to_numpy(), columns=state_ids)


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

    adata = _read_adata(input_dir / "query")
    input_n_features = int(adata.n_vars)
    signatures_path = config.get("reference_signatures")
    if not signatures_path:
        raise ValueError("cell2location object workflow requires `reference_signatures`.")
    signatures = _read_reference_signatures(Path(signatures_path))
    common = adata.var_names.intersection(signatures.index.astype(str))
    if len(common) == 0:
        raise ValueError("No shared genes between spatial object and reference signatures.")
    adata = adata[:, common].copy()
    signatures = signatures.loc[common, :]

    import torch

    import cell2location
    from cell2location.models import Cell2location

    cell2location.models.Cell2location.setup_anndata(adata=adata)
    model = Cell2location(
        adata,
        cell_state_df=signatures,
        N_cells_per_location=config.get("n_cells_per_location", 30),
        detection_alpha=config.get("detection_alpha", 20),
    )
    model.train(max_epochs=int(config.get("max_epochs", 30000)), batch_size=config.get("batch_size"))
    adata = model.export_posterior(adata, sample_kwargs={"num_samples": int(config.get("num_samples", 1000))})

    posterior = adata.obsm.get("q05_cell_abundance_w_sf")
    if posterior is None:
        posterior = adata.obsm.get("means_cell_abundance_w_sf")
    if posterior is None:
        raise RuntimeError("cell2location completed without posterior cell-abundance estimates.")
    abundance = pd.DataFrame(posterior, index=adata.obs_names)
    abundance.to_csv(output_dir / "cell_abundance.csv")
    abundance.to_csv(output_dir / "obs.csv")
    if config.get("write_h5ad", False):
        adata.write_h5ad(output_dir / "cell2location.h5ad")
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "method": "cell2location",
                "output_h5ad": str(output_dir / "cell2location.h5ad") if (output_dir / "cell2location.h5ad").exists() else None,
                "n_cells": int(adata.n_obs),
                "n_features": input_n_features,
                "n_shared_features": int(adata.n_vars),
                "cell2location_version": cell2location.__version__,
                "torch_version": torch.__version__,
            },
            handle,
            indent=2,
        )


if __name__ == "__main__":
    main()
