#!/usr/bin/env python
"""Run cell2location from Shennong-exported spatial Seurat data."""

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

    import cell2location
    from cell2location.models import Cell2location

    adata = _read_adata(input_dir / "query")
    signatures_path = config.get("reference_signatures")
    if not signatures_path:
        raise ValueError("cell2location object workflow requires `reference_signatures`.")
    signatures = pd.read_csv(signatures_path, index_col=0)
    common = adata.var_names.intersection(signatures.index.astype(str))
    if len(common) == 0:
        raise ValueError("No shared genes between spatial object and reference signatures.")
    adata = adata[:, common].copy()
    signatures = signatures.loc[common, :]

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
    if posterior is not None:
        abundance = pd.DataFrame(posterior, index=adata.obs_names)
        abundance.to_csv(output_dir / "cell_abundance.csv")
        abundance.to_csv(output_dir / "obs.csv")
    else:
        adata.obs.to_csv(output_dir / "obs.csv")
    if config.get("write_h5ad", True):
        adata.write_h5ad(output_dir / "cell2location.h5ad")
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump({"method": "cell2location", "output_h5ad": str(output_dir / "cell2location.h5ad")}, handle, indent=2)


if __name__ == "__main__":
    main()
