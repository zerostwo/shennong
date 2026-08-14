#!/usr/bin/env python
"""Train a real scPoli latent integration model for Shennong."""

from __future__ import annotations

import argparse
import importlib.metadata as metadata
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import torch
from scarches.models.scpoli import scPoli


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _drop_none(mapping: dict | None) -> dict:
    if not isinstance(mapping, dict):
        return {}
    return {key: value for key, value in mapping.items() if value is not None}


def _read_input(input_dir: Path) -> ad.AnnData:
    if (input_dir / "counts.mtx").exists():
        matrix_path = input_dir / "counts.mtx"
        obs_path = input_dir / "obs.csv"
        var_path = input_dir / "features.csv"
    else:
        query_dir = input_dir / "query"
        matrix_path = query_dir / "matrix.mtx"
        obs_path = query_dir / "obs.csv"
        var_path = query_dir / "var.csv"

    matrix = sp.csr_matrix(scipy.io.mmread(matrix_path)).transpose().tocsr()
    obs = pd.read_csv(obs_path, dtype=str)
    if "cell_id" not in obs.columns:
        raise ValueError("obs.csv must contain a cell_id column.")
    obs = obs.set_index("cell_id")
    var = pd.read_csv(var_path, dtype=str)
    feature_column = "feature_id" if "feature_id" in var.columns else var.columns[0]
    features = var[feature_column].astype(str).to_numpy()
    if matrix.shape != (obs.shape[0], len(features)):
        raise ValueError("Exported expression matrix dimensions do not match obs.csv and feature metadata.")
    if matrix.data.size and (not np.isfinite(matrix.data).all() or np.min(matrix.data) < 0):
        raise ValueError("scPoli input must contain finite, non-negative expression values.")
    adata = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=features))
    adata.obs_names = obs.index.astype(str)
    adata.var_names_make_unique()
    if not sp.issparse(adata.X):
        raise MemoryError("The complete scPoli expression matrix must remain sparse.")
    return adata


def _version(package: str) -> str | None:
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def _get_latent_in_sparse_batches(model: scPoli, adata: ad.AnnData, batch_size: int) -> np.ndarray:
    """Encode sparse input without ever densifying the complete expression matrix."""
    device = next(model.model.parameters()).device
    condition_codes = []
    for condition_key in model.condition_keys_:
        values = adata.obs[condition_key].astype(str).to_numpy()
        encoder = model.model.condition_encoders[condition_key]
        unknown = set(values).difference(encoder)
        if unknown:
            raise ValueError(f"Incorrect conditions for {condition_key!r}: {sorted(unknown)}")
        condition_codes.append(np.asarray([encoder[value] for value in values], dtype=np.int64))
    conditions = torch.as_tensor(np.column_stack(condition_codes), device=device)

    latent_batches = []
    for start in range(0, adata.n_obs, batch_size):
        stop = min(start + batch_size, adata.n_obs)
        expression_batch = adata.X[start:stop, :]
        if sp.issparse(expression_batch):
            # scPoli's neural network consumes dense tensors, but this conversion
            # is bounded to one minibatch; adata.X remains sparse throughout.
            expression_batch = expression_batch.toarray()
        expression_batch = torch.as_tensor(expression_batch, device=device, dtype=torch.float32)
        latent = model.model.get_latent(expression_batch, conditions[start:stop, :], mean=True)
        latent_batches.append(latent.detach().cpu().numpy())
    return np.concatenate(latent_batches, axis=0)


def run(input_dir: Path, output_dir: Path, config: dict) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = _read_input(input_dir)
    batch_key = config.get("batch_key")
    labels_key = config.get("labels_key")
    if not batch_key or batch_key not in adata.obs.columns:
        raise ValueError("scPoli integration requires a batch_key present in obs.csv.")
    adata.obs[batch_key] = adata.obs[batch_key].astype(str).astype("category")
    if labels_key:
        if labels_key not in adata.obs.columns:
            raise ValueError(f"Labels key {labels_key!r} was not found in obs.csv.")
        adata.obs[labels_key] = adata.obs[labels_key].astype(str).astype("category")

    seed = int(config.get("seed") or 717)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    model_args = _drop_none(config.get("model_args"))
    model_args.setdefault("condition_keys", batch_key)
    model_args.setdefault("cell_type_keys", labels_key)
    model_args.setdefault("latent_dim", int(config.get("n_latent") or 10))
    model_args.setdefault("embedding_dims", int(config.get("embedding_dims") or 5))
    model_args.setdefault("recon_loss", "nb")
    model = scPoli(adata=adata, **model_args)

    n_epochs = max(1, int(config.get("n_epochs") or config.get("max_epochs") or 100))
    pretraining_epochs = int(config.get("pretraining_epochs") or int(n_epochs * 0.9))
    pretraining_epochs = max(0, min(pretraining_epochs, n_epochs))
    train_args = _drop_none(config.get("train_args"))
    train_args.setdefault("n_epochs", n_epochs)
    train_args.setdefault("pretraining_epochs", pretraining_epochs)
    if not labels_key:
        train_args.setdefault("prototype_training", False)
        train_args.setdefault("unlabeled_prototype_training", False)
    model.train(**train_args)
    model.model.eval()

    latent_batch_size = max(1, int(config.get("latent_batch_size") or 2048))
    latent = _get_latent_in_sparse_batches(model, adata, latent_batch_size)
    latent_path = output_dir / "latent.csv"
    pd.DataFrame(
        latent,
        index=adata.obs_names,
        columns=[f"SCPOLI_{idx + 1}" for idx in range(latent.shape[1])],
    ).to_csv(latent_path)
    obs_out = pd.DataFrame(index=adata.obs_names)
    obs_path = output_dir / "obs.csv"
    obs_out.to_csv(obs_path)

    model_dir = output_dir / "model"
    if bool(config.get("save_model", True)):
        model.save(model_dir, overwrite=True, save_anndata=False)
    h5ad_path = output_dir / "integrated.h5ad"
    if bool(config.get("write_h5ad", True)):
        adata.obsm["X_scpoli"] = latent
        adata.write_h5ad(h5ad_path)

    manifest = {
        "method": "scpoli",
        "batch_key": batch_key,
        "labels_key": labels_key,
        "source_layer": config.get("source_layer") or config.get("layer"),
        "latent_csv": str(latent_path),
        "obs_csv": str(obs_path),
        "output_h5ad": str(h5ad_path) if h5ad_path.exists() else None,
        "model_dir": str(model_dir) if model_dir.exists() else None,
        "n_cells": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "scarches_version": _version("scarches"),
        "torch_version": _version("torch"),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(Path(args.input_dir), Path(args.output_dir), _read_json(Path(args.config)))


if __name__ == "__main__":
    main()
