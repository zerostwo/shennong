#!/usr/bin/env python3
"""Assemble the local scientific-figure review after atlas-figure-layouts.R.

Requires PyMuPDF and Pillow, only for local PDF inspection/contact sheets.
No reference images are modified or included in the installed R package.
"""
from pathlib import Path
import csv
import html
import json

import pymupdf
from PIL import Image, ImageDraw

root = Path("data-local/atlas-layout-review")
rows = list(csv.DictReader(Path("inst/extdata/zhang-atlas-layouts.csv").open()))
unique = {}
for row in rows:
    unique.setdefault(row["style"], row)
document = pymupdf.open()
font_sizes, outside = set(), []
for style in unique:
    pdf = pymupdf.open(root / style / f"{style}.pdf")
    document.insert_pdf(pdf)
    for page in pdf:
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                for span in line["spans"]:
                    font_sizes.add(round(span["size"], 2))
                    if not page.rect.contains(pymupdf.Rect(span["bbox"])):
                        outside.append([style, span["text"]])
        page.get_pixmap(matrix=pymupdf.Matrix(1.7, 1.7)).save(
            root / style / f"{style}-pdf-review.png")
document.save(root / "atlas-layout-gallery.pdf")
selected = ["split_triangle", "circular", "annotated", "raincloud", "split_violin",
            "paired_box", "forest", "upset", "bipartite", "clone_tree", "chord", "rewiring"]
canvas = Image.new("RGB", (1800, 1800), "white")
draw = ImageDraw.Draw(canvas)
for i, style in enumerate(selected):
    picture = Image.open(root / style / f"{style}-pdf-review.png").convert("RGB")
    picture.thumbnail((580, 420))
    x, y = (i % 3) * 600 + (600 - picture.width) // 2, (i // 3) * 450 + 25
    canvas.paste(picture, (x, y))
    draw.text(((i % 3) * 600 + 18, (i // 3) * 450 + 5), style, fill="#173247")
canvas.save(root / "gallery-overview.png")
report = dict(pdf_pages=len(document), font_sizes_pt=sorted(font_sizes),
              out_of_bounds_text=outside)
(root / "pdf-review.json").write_text(json.dumps(report, indent=2))

parts = ['''<!doctype html><meta charset="utf-8"><title>Shennong atlas layouts</title>
<style>body{max-width:1450px;margin:30px auto;padding:20px;font:16px/1.65 system-ui;
color:#172c3b;background:#f4f6f8}header{margin-bottom:30px}main{display:grid;
grid-template-columns:repeat(auto-fit,minmax(450px,1fr));gap:24px}article{background:white;
border:1px solid #d9e0e5;border-radius:8px;padding:18px}img{width:100%;height:auto}
code{background:#edf2f7;padding:3px 7px}a{color:#176794}input{width:85%;padding:12px;
font-size:16px}h1{font-size:30px}h2{font-size:20px}small{color:#506170}</style>
<header><h1>Shennong：23 个论文图形结构样式</h1>
<p>按结构复刻；每张图可打开 PDF、PNG 和源数据清单。以下演示为合成数据，
不代表论文生物学结果。Nature Fig. 2c 原始 408 个 Ro/e 值另列。</p>
<p><a href="atlas-layout-gallery.pdf">全部布局 PDF</a> ·
<a href="nature-original-values/nature-original-values.pdf">Nature 原始数值复画</a> ·
<a href="layout-catalog.csv">面板与函数对照表</a></p>
<input placeholder="筛选：style / function / paper" oninput="document.querySelectorAll('article').forEach(x=>x.hidden=!x.textContent.toLowerCase().includes(this.value.toLowerCase()))"></header><main>''']
for style, row in unique.items():
    esc = html.escape
    parts.append(f'''<article><h2>{esc(style)}</h2>
<code>{esc(row['function_name'])}(..., style="{esc(style)}")</code>
<p>{esc(row['reference'])}</p><a href="{style}/{style}.pdf">
<img src="{style}/{style}.png" loading="lazy"></a>
<p>{esc(row['description'])}</p><small>{esc(row['data_requirements'])}</small>
<p><a href="{style}/{style}.pdf">PDF</a> · <a href="{style}/{style}.png">PNG</a> ·
<a href="{style}/{style}_manifest.json">数据清单</a></p></article>''')
parts.append('''</main><h2>原图核对入口（本地参考，不随包分发）</h2><p>
<a href="../zhang-atlas-figures/b-f2.jpg">B Fig. 2：UpSet</a> ·
<a href="../zhang-atlas-figures/b-f4.jpg">B Fig. 4：双三角</a> ·
<a href="../zhang-atlas-figures/cm-p7.png">Nature Fig. 4：环形／分裂小提琴</a> ·
<a href="../zhang-atlas-figures/fib-f6.jpg">Fibroblast Fig. 6：弦图／配对</a></p>
<p>尚未覆盖：ATAC、FACS、显微图像、Venn、圆形系统发育树及所有补充图。
T 细胞论文目前仅有图注核对。Raincloud 是基于原图视觉元素的改编。</p>''')
(root / "index.html").write_text("\n".join(parts))
print(json.dumps(report))
assert not outside, "PDF text extends beyond the page"
assert font_sizes == {8.0}, "Unexpected PDF typography"
