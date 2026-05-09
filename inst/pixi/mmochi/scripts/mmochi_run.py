#!/usr/bin/env python
"""Run MMoCHi ADT landmark registration for Shennong."""

from __future__ import annotations

import argparse
import importlib.metadata as metadata
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import mmochi as mmc


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _drop_none(mapping: dict | None) -> dict:
    if mapping is None or not isinstance(mapping, dict):
        return {}
    return {key: value for key, value in mapping.items() if value is not None}


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return value.item()
    if not isinstance(value, (list, tuple, dict, np.ndarray)) and pd.isna(value):
        return None
    return value


def _version(package: str) -> str | None:
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return getattr(mmc, "__version__", None)


def _read_input(input_dir: Path, data_key: str) -> ad.AnnData:
    protein = pd.read_csv(input_dir / "protein.csv")
    obs = pd.read_csv(input_dir / "obs.csv")
    if "cell_id" not in protein.columns:
        raise ValueError("protein.csv must contain a cell_id column.")
    if "cell_id" not in obs.columns:
        raise ValueError("obs.csv must contain a cell_id column.")

    protein = protein.set_index("cell_id")
    obs = obs.set_index("cell_id").reindex(protein.index)
    if obs.shape[1] > 0 and obs.isna().all(axis=1).any():
        raise ValueError("obs.csv does not contain metadata for the exported cells.")

    var = pd.DataFrame(index=["shennong_dummy"])
    adata = ad.AnnData(
        X=np.zeros((protein.shape[0], 1), dtype=np.float32),
        obs=obs,
        var=var,
    )
    adata.obs_names = protein.index.astype(str)
    adata.var_names = var.index.astype(str)
    adata.obsm[data_key] = protein.astype(float)
    return adata


def _coerce_inclusion_mask(adata: ad.AnnData, inclusion_mask: Any) -> Any:
    if inclusion_mask is None:
        return None
    if isinstance(inclusion_mask, str):
        if inclusion_mask not in adata.obs.columns:
            raise ValueError(f"inclusion_mask column {inclusion_mask!r} was not found in obs.csv.")
        return inclusion_mask
    mask = np.asarray(inclusion_mask)
    if mask.shape[0] != adata.n_obs:
        raise ValueError("inclusion_mask must have one value per exported cell.")
    return mask.astype(bool)


def run(input_dir: Path, output_dir: Path, config: dict) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    data_key = config.get("data_key") or "protein"
    key_added = config.get("key_added") or "landmark_protein"
    batch_key = config.get("batch_key")
    if batch_key is None:
        raise ValueError("MMoCHi landmark registration requires a batch_key.")

    adata = _read_input(input_dir=input_dir, data_key=data_key)
    if batch_key not in adata.obs.columns:
        raise ValueError(f"Batch key {batch_key!r} was not found in obs.csv.")

    inclusion_mask = _coerce_inclusion_mask(adata, config.get("inclusion_mask"))
    landmark_args = _drop_none(config.get("landmark_args"))
    mmc.landmark_register_adts(
        adata,
        batch_key=batch_key,
        data_key=data_key,
        key_added=key_added,
        show=config.get("show", False),
        single_peaks=config.get("single_peaks") or [],
        marker_bandwidths=config.get("marker_bandwidths") or {},
        peak_overrides=config.get("peak_overrides") or {},
        inclusion_mask=inclusion_mask,
        **landmark_args,
    )

    corrected = adata.obsm[key_added]
    if not isinstance(corrected, pd.DataFrame):
        corrected = pd.DataFrame(corrected, index=adata.obs_names, columns=adata.obsm[data_key].columns)
    corrected.to_csv(output_dir / "landmark_protein.csv")
    adata.obs.to_csv(output_dir / "obs.csv")

    peaks = {
        key_added: adata.uns.get(f"{key_added}_peaks", {}),
        data_key: adata.uns.get(f"{data_key}_peaks", {}),
    }
    with (output_dir / "peaks.json").open("w", encoding="utf-8") as handle:
        json.dump(_jsonable(peaks), handle, indent=2)

    manifest = {
        "method": "mmochi",
        "batch_key": batch_key,
        "data_key": data_key,
        "key_added": key_added,
        "corrected_csv": str(output_dir / "landmark_protein.csv"),
        "obs_csv": str(output_dir / "obs.csv"),
        "peaks_json": str(output_dir / "peaks.json"),
        "n_cells": int(adata.n_obs),
        "n_proteins": int(corrected.shape[1]),
        "mmochi_version": _version("mmochi"),
        "scanpy_version": sc.__version__,
    }
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    run(
        input_dir=Path(args.input_dir),
        output_dir=Path(args.output_dir),
        config=_read_json(Path(args.config)),
    )


if __name__ == "__main__":
    main()
