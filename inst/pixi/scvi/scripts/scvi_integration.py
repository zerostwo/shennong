#!/usr/bin/env python
"""Run scVI/scANVI for Shennong's R integration wrapper."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import scvi


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _drop_none(mapping: dict | None) -> dict:
    if mapping is None:
        return {}
    if not isinstance(mapping, dict):
        return {}
    return {key: value for key, value in mapping.items() if value is not None}


def _read_input(input_dir: Path) -> ad.AnnData:
    counts = scipy.io.mmread(input_dir / "counts.mtx")
    counts = sp.csr_matrix(counts).transpose().tocsr()

    features = pd.read_csv(input_dir / "features.csv")["feature_id"].astype(str).to_numpy()
    cells = pd.read_csv(input_dir / "cells.csv")["cell_id"].astype(str).to_numpy()
    obs = pd.read_csv(input_dir / "obs.csv", dtype=str)
    if "cell_id" not in obs.columns:
        raise ValueError("obs.csv must contain a cell_id column.")
    obs = obs.set_index("cell_id").reindex(cells)

    var = pd.DataFrame(index=features)
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.obs_names = cells
    adata.var_names = features
    adata.var_names_make_unique()

    protein_counts_path = input_dir / "protein_counts.mtx"
    proteins_path = input_dir / "proteins.csv"
    if protein_counts_path.exists() and proteins_path.exists():
        protein_counts = scipy.io.mmread(protein_counts_path)
        protein_counts = sp.csr_matrix(protein_counts).transpose().tocsr()
        proteins = pd.read_csv(proteins_path)["protein_id"].astype(str).to_numpy()
        if protein_counts.shape[0] != adata.n_obs:
            raise ValueError("protein_counts.mtx must have one row per cell after transposition.")
        adata.obsm["protein_expression"] = pd.DataFrame(
            protein_counts.toarray(),
            index=adata.obs_names,
            columns=proteins,
        )
    return adata


def _write_latent(latent: np.ndarray, cells: pd.Index, output_dir: Path, prefix: str) -> Path:
    latent_df = pd.DataFrame(
        latent,
        index=cells,
        columns=[f"{prefix}_{idx + 1}" for idx in range(latent.shape[1])],
    )
    path = output_dir / "latent.csv"
    latent_df.to_csv(path)
    return path


def _prepare_labels(adata: ad.AnnData, labels_key: str, unlabeled_category: str) -> None:
    labels = adata.obs[labels_key].astype("object")
    labels = labels.where(labels.notna() & (labels.astype(str) != ""), unlabeled_category)
    adata.obs[labels_key] = labels.astype(str)


def run(input_dir: Path, output_dir: Path, config: dict) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    method = config.get("method", "scvi")
    batch_key = config.get("batch_key")
    labels_key = config.get("labels_key")
    unlabeled_category = config.get("unlabeled_category", "Unknown")
    n_latent = int(config.get("n_latent") or 30)
    seed = int(config.get("seed") or 717)

    scvi.settings.seed = seed
    adata = _read_input(input_dir)
    if batch_key is not None and batch_key not in adata.obs.columns:
        raise ValueError(f"Batch key {batch_key!r} was not found in obs.csv.")

    setup_kwargs = {}
    if batch_key is not None:
        setup_kwargs["batch_key"] = batch_key

    obs_out = pd.DataFrame(index=adata.obs_names)
    if method == "totalvi":
        protein_obsm_key = config.get("protein_obsm_key") or "protein_expression"
        if "protein_expression" not in adata.obsm:
            raise ValueError("totalVI requires protein_counts.mtx and proteins.csv in the input directory.")
        if protein_obsm_key != "protein_expression":
            adata.obsm[protein_obsm_key] = adata.obsm["protein_expression"]
        totalvi_setup_kwargs = dict(setup_kwargs)
        totalvi_setup_kwargs["protein_expression_obsm_key"] = protein_obsm_key
        scvi.model.TOTALVI.setup_anndata(adata, **totalvi_setup_kwargs)
        totalvi_model_args = _drop_none(config.get("totalvi_model_args"))
        totalvi_model_args.setdefault("n_latent", n_latent)
        backend = scvi.model.TOTALVI(adata, **totalvi_model_args)
        totalvi_train_args = _drop_none(config.get("totalvi_train_args"))
        if config.get("max_epochs") is not None:
            totalvi_train_args.setdefault("max_epochs", int(config["max_epochs"]))
        backend.train(**totalvi_train_args)
        latent = backend.get_latent_representation()
        latent_path = _write_latent(latent, adata.obs_names, output_dir, method.upper())
        obs_path = output_dir / "obs.csv"
        obs_out.to_csv(obs_path)
        h5ad_path = output_dir / "integrated.h5ad"
        if bool(config.get("write_h5ad", True)):
            adata.obsm[f"X_{method}"] = latent
            adata.write_h5ad(h5ad_path)
        manifest = {
            "method": method,
            "batch_key": batch_key,
            "labels_key": labels_key,
            "latent_csv": str(latent_path),
            "obs_csv": str(obs_path),
            "output_h5ad": str(h5ad_path) if h5ad_path.exists() else None,
            "n_cells": int(adata.n_obs),
            "n_features": int(adata.n_vars),
            "n_proteins": int(adata.obsm[protein_obsm_key].shape[1]),
            "scvi_version": scvi.__version__,
            "scanpy_version": sc.__version__,
        }
        with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2)
        return

    scvi.model.SCVI.setup_anndata(adata, **setup_kwargs)

    model_args = _drop_none(config.get("model_args"))
    model_args.setdefault("n_latent", n_latent)
    model = scvi.model.SCVI(adata, **model_args)

    train_args = _drop_none(config.get("train_args"))
    if config.get("max_epochs") is not None:
        train_args.setdefault("max_epochs", int(config["max_epochs"]))
    model.train(**train_args)

    backend = model
    if method == "scanvi":
        if labels_key is None or labels_key not in adata.obs.columns:
            raise ValueError("scANVI requires a labels_key present in obs.csv.")
        _prepare_labels(adata, labels_key, unlabeled_category)
        scanvi_model_args = _drop_none(config.get("scanvi_model_args"))
        backend = scvi.model.SCANVI.from_scvi_model(
            model,
            labels_key=labels_key,
            unlabeled_category=unlabeled_category,
            **scanvi_model_args,
        )
        scanvi_train_args = _drop_none(config.get("scanvi_train_args"))
        if config.get("scanvi_max_epochs") is not None:
            scanvi_train_args.setdefault("max_epochs", int(config["scanvi_max_epochs"]))
        backend.train(**scanvi_train_args)
        obs_out["scanvi_prediction"] = backend.predict()

    latent = backend.get_latent_representation()
    latent_path = _write_latent(latent, adata.obs_names, output_dir, method.upper())
    obs_path = output_dir / "obs.csv"
    obs_out.to_csv(obs_path)

    h5ad_path = output_dir / "integrated.h5ad"
    if bool(config.get("write_h5ad", True)):
        adata.obsm[f"X_{method}"] = latent
        if method == "scanvi" and "scanvi_prediction" in obs_out.columns:
            adata.obs["scanvi_prediction"] = obs_out["scanvi_prediction"]
        adata.write_h5ad(h5ad_path)

    manifest = {
        "method": method,
        "batch_key": batch_key,
        "labels_key": labels_key,
        "latent_csv": str(latent_path),
        "obs_csv": str(obs_path),
        "output_h5ad": str(h5ad_path) if h5ad_path.exists() else None,
        "n_cells": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "scvi_version": scvi.__version__,
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
