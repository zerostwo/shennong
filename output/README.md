# Shennong visual guide

These artifacts summarize the current object-centered workflow and plotting
surface. They are generated from the live package registries rather than a
manually copied method inventory.

- `diagrams/shennong-object-workflow.drawio`: editable diagrams.net source.
- `diagrams/shennong-object-workflow.svg`: browser-ready workflow preview.
- `pdf/shennong-object-workflow.pdf`: print-ready workflow diagram.
- `pdf/shennong-cheatsheet.pdf`: two-page bilingual API cheat sheet.
- `shennong-visual-guide-data.json`: generation snapshot with package version,
  public exports, analysis methods, and result-plot views.

Regenerate from the repository root with:

```sh
uv run --with reportlab --with pdfplumber --with pypdf --with pymupdf \
  --with fonttools python scripts/generate_shennong_visual_guide.py
```

The generator validates its curated function groups against the current
namespace before writing output. The committed visual artifacts are excluded
from the package source tarball through `.Rbuildignore`.
