# seurat_merge changelog

## 2026-08-01

- Added the isolated `SeuratObject::merge.Assay5` accelerator.
- Restricted registration to one method; v3, SCT, subclasses, malformed calls,
  and custom backend attributes remain on captured upstream behavior.
- Added exact contract coverage for metadata, LogMaps, raw sparse attributes,
  input immutability, lazy dots, opt-out, and fallback behavior.
