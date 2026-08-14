#!/usr/bin/env python
"""Build a batch-balanced KNN graph from a Shennong PCA export."""

from __future__ import annotations

import argparse
import importlib.metadata as metadata
import json
from pathlib import Path
import sys

_SHENNONG_PIXI_SHARED = Path(__file__).resolve().parents[2] / "_shared"
sys.path.insert(0, str(_SHENNONG_PIXI_SHARED))
from shennong_autozyme import activate_autozyme

AUTOZYME_STATUS = activate_autozyme(["scanpy"])

import bbknn.matrix
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _drop_none(mapping: dict | None) -> dict:
    if not isinstance(mapping, dict):
        return {}
    return {key: value for key, value in mapping.items() if value is not None}


def _version(package: str) -> str | None:
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def run(input_dir: Path, output_dir: Path, config: dict) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    pca = pd.read_csv(input_dir / "pca.csv", index_col=0)
    obs = pd.read_csv(input_dir / "obs.csv", dtype=str)
    if not {"cell_id", "batch"}.issubset(obs.columns):
        raise ValueError("obs.csv must contain cell_id and batch columns.")
    obs = obs.set_index("cell_id").reindex(pca.index)
    if obs["batch"].isna().any():
        raise ValueError("BBKNN metadata is missing batch values for exported cells.")

    seed = int(config.get("seed") or 717)
    np.random.seed(seed)
    pca_values = pca.to_numpy(dtype=float, copy=False)
    bbknn_args = _drop_none(config.get("bbknn_args"))
    bbknn_args.setdefault("pynndescent_random_state", seed)
    bbknn_args.setdefault("use_faiss", False)
    distances, connectivities, parameters = bbknn.matrix.bbknn(
        pca_values,
        obs["batch"].astype(str).to_numpy(),
        **bbknn_args,
    )
    distances = sp.csr_matrix(distances)
    connectivities = sp.csr_matrix(connectivities)
    scipy.io.mmwrite(output_dir / "distances.mtx", distances)
    scipy.io.mmwrite(output_dir / "connectivities.mtx", connectivities)

    adata = ad.AnnData(X=np.zeros((pca.shape[0], 1), dtype=np.float32), obs=obs.copy())
    adata.obsm["X_pca"] = pca_values
    adata.obsp["distances"] = distances
    adata.obsp["connectivities"] = connectivities
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": parameters,
    }
    umap_args = _drop_none(config.get("umap_args"))
    aliases = {
        "seed.use": "random_state",
        "min.dist": "min_dist",
        "n.components": "n_components",
        "n.epochs": "maxiter",
    }
    umap_args = {aliases.get(key, key): value for key, value in umap_args.items()}
    for unsupported in ("graph", "reduction", "reduction.name", "reduction.key", "umap.method", "verbose"):
        umap_args.pop(unsupported, None)
    umap_args.setdefault("random_state", seed)
    sc.tl.umap(adata, **umap_args)
    pd.DataFrame(
        adata.obsm["X_umap"],
        index=pca.index,
        columns=[f"UMAP_{idx + 1}" for idx in range(adata.obsm["X_umap"].shape[1])],
    ).to_csv(output_dir / "umap.csv")

    manifest = {
        "method": "bbknn",
        "graph_name": config.get("graph_name") or "bbknn_snn",
        "n_cells": int(pca.shape[0]),
        "n_pcs": int(pca.shape[1]),
        "n_batches": int(obs["batch"].nunique()),
        "umap_csv": str(output_dir / "umap.csv"),
        "parameters": parameters,
        "bbknn_version": _version("bbknn"),
        "autozyme": AUTOZYME_STATUS,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(
            manifest,
            indent=2,
            default=lambda value: value.item() if isinstance(value, np.generic) else str(value),
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(Path(args.input_dir), Path(args.output_dir), _read_json(Path(args.config)))


if __name__ == "__main__":
    main()
