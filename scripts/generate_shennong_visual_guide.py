#!/usr/bin/env python3
"""Generate the Shennong object workflow diagram and package cheat sheet.

Run from the repository root with:

    uv run --with reportlab --with pypdf --with pdfplumber \
      --with pymupdf --with fonttools python scripts/generate_shennong_visual_guide.py

The script queries the live package registries before drawing. Curated module
cards define the explanatory layer; validation fails if a named public function
is no longer exported.
"""

from __future__ import annotations

import html
import json
import math
import re
import subprocess
import textwrap
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from reportlab.lib.colors import HexColor, white
from reportlab.lib.pagesizes import A3, A4, landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
OUT_PDF = ROOT / "output" / "pdf"
OUT_DIAGRAM = ROOT / "output" / "diagrams"

FONT_SOURCE = ROOT / "tmp" / "pdfs" / "NotoSansSC-wght.ttf"
FONT_REGULAR = ROOT / "tmp" / "pdfs" / "NotoSansSC-Regular.ttf"
FONT_BOLD = ROOT / "tmp" / "pdfs" / "NotoSansSC-Bold.ttf"
FONT_MONO = "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"
FONT_URL = "https://github.com/google/fonts/raw/main/ofl/notosanssc/NotoSansSC%5Bwght%5D.ttf"

COLORS = {
    "ink": "#16324F",
    "muted": "#52667A",
    "paper": "#F7FAFC",
    "line": "#CBD5E1",
    "core": "#3751A6",
    "io": "#0F8B8D",
    "cluster": "#2F80ED",
    "annotation": "#7B61A8",
    "program": "#D97706",
    "composition": "#D6456D",
    "dynamic": "#2E8B57",
    "spatial": "#23967F",
    "bulk": "#5B677A",
    "result": "#374151",
    "plot": "#006D77",
    "export": "#7C3AED",
    "method": "#64748B",
    "white": "#FFFFFF",
}


@dataclass(frozen=True)
class Node:
    key: str
    x: float
    y: float
    width: float
    height: float
    title: str
    lines: tuple[str, ...]
    fill: str
    stroke: str | None = None
    title_color: str = "#FFFFFF"
    text_color: str = "#FFFFFF"
    dashed: bool = False
    container: bool = False


@dataclass(frozen=True)
class Edge:
    source: str
    target: str
    label: str = ""
    dashed: bool = False
    color: str = "#64748B"


MODULES = [
    {
        "title": "Data & I/O 数据与读写",
        "color": COLORS["io"],
        "functions": [
            ("sn_list_10x_paths()", "发现 10x/STARsolo/spatial 输入"),
            ("sn_read() / sn_write()", "qs2, h5ad, BPCells 与表格 I/O"),
            ("sn_set_layer_backend()", "dgCMatrix <-> BPCells layer"),
        ],
    },
    {
        "title": "Preprocess & QC 预处理",
        "color": COLORS["io"],
        "functions": [
            ("sn_initialize_seurat_object()", "创建中心 Seurat/Shennong object"),
            ("sn_add_qc_metrics()", "按 assay/layer 重算 QC"),
            ("sn_filter_cells() / sn_filter_genes()", "过滤细胞和基因"),
            ("sn_normalize_data()", "log, scran 或 SCTransform"),
            ("sn_find_doublets()", "scDblFinder / Scrublet"),
            ("sn_remove_ambient_contamination()", "decontX/decontPro/SoupX"),
        ],
    },
    {
        "title": "Cluster & Integrate 聚类整合",
        "color": COLORS["cluster"],
        "functions": [
            ("sn_run_cluster()", "聚类、参数网格、批次整合"),
            ("sn_run_multimodal()", "WNN/totalVI/Coralysis/MMoCHi"),
            ("sn_transfer_labels()", "Seurat/Coralysis/scANVI/scArches"),
            ("sn_assess_integration()", "整合质量评估"),
        ],
    },
    {
        "title": "Annotate, DE & Pathways 注释与差异",
        "color": COLORS["annotation"],
        "functions": [
            ("sn_run_annotation()", "SingleR/CellTypist/PopV/..."),
            ("sn_review_annotation()", "低置信度证据审查"),
            ("sn_find_de()", "markers/contrast/pseudobulk/bulk"),
            ("sn_run_enrichment()", "ORA/GSEA"),
        ],
    },
    {
        "title": "Programs, GRN & Metabolism 程序调控",
        "color": COLORS["program"],
        "functions": [
            ("sn_score_programs() / sn_test_programs()", "UCell/AUCell/GSVA/ssGSEA/mean"),
            ("sn_discover_programs()", "NMF/cNMF/Hotspot"),
            ("sn_run_grn()", "GENIE3/SCENIC/GRNBoost2"),
            ("sn_run_metabolism()", "gene sets/scMetabolism/scFEA/Compass"),
        ],
    },
    {
        "title": "Composition, DA, State & CNV 组成与状态",
        "color": COLORS["composition"],
        "functions": [
            ("sn_calculate_composition()", "计数与比例"),
            ("sn_compare_composition()", "样本作为统计重复"),
            ("sn_test_abundance()", "Propeller/Milo/scCODA/permutation"),
            ("sn_prioritize_states()", "Augur/RareQ/Scissor"),
            ("sn_run_cnv()", "inferCNVpy/CopyKAT"),
        ],
    },
    {
        "title": "Communication & Dynamics 通讯与动态",
        "color": COLORS["dynamic"],
        "functions": [
            ("sn_run_cell_communication()", "LIANA/CellChat/NicheNet/..."),
            ("sn_run_trajectory()", "Slingshot/Monocle3/Palantir"),
            ("sn_run_velocity()", "scVelo/RegVelo"),
            ("sn_run_fate()", "CellRank"),
        ],
    },
    {
        "title": "Spatial 空间组学",
        "color": COLORS["spatial"],
        "functions": [
            ("sn_run_spatial()", "统一 task dispatcher"),
            ("sn_find_spatial_features()", "Moran's I/nnSVG/SPARK-X"),
            ("sn_find_spatial_domains()", "BANKSY/stLearn/BayesSpace/..."),
            ("sn_run_spatial_neighborhood()", "邻域富集与共现"),
            ("sn_run_spatial_deconvolution()", "cell2location"),
        ],
    },
    {
        "title": "Bulk transcriptomics Bulk 转录组",
        "color": COLORS["bulk"],
        "functions": [
            ("sn_run_bulk()", "QC/DE/pathway/network/survival dispatcher"),
            ("sn_assess_bulk_qc()", "PCA、相关性与异常值"),
            ("sn_score_bulk_pathways()", "sample-level pathway scores"),
            ("sn_run_wgcna() / sn_run_survival()", "网络与生存"),
            ("sn_run_clinical_association()", "临床表型关联"),
        ],
    },
]


def register_fonts() -> None:
    from fontTools.ttLib import TTFont as VariableFont
    from fontTools.varLib.instancer import instantiateVariableFont

    FONT_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    if not FONT_SOURCE.exists():
        urllib.request.urlretrieve(FONT_URL, FONT_SOURCE)
    for path, weight in ((FONT_REGULAR, 400), (FONT_BOLD, 700)):
        if path.exists():
            continue
        font = VariableFont(FONT_SOURCE)
        instantiateVariableFont(font, {"wght": weight}, inplace=True)
        font.save(path)
        font.close()
    pdfmetrics.registerFont(TTFont("NotoCJK", str(FONT_REGULAR)))
    pdfmetrics.registerFont(TTFont("NotoCJKBold", str(FONT_BOLD)))
    pdfmetrics.registerFont(TTFont("GuideMono", FONT_MONO))


def collect_registry() -> dict:
    expression = r'''
      suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
      methods <- sn_list_methods()
      methods <- methods[, c("task", "name", "default", "implemented", "runtime")]
      plots <- sn_list_plot_methods()
      exports <- sort(getNamespaceExports("Shennong"))
      payload <- list(
        version = as.character(utils::packageVersion("Shennong")),
        exports = exports,
        methods = methods,
        plots = plots
      )
      cat(jsonlite::toJSON(payload, dataframe = "rows", auto_unbox = TRUE, null = "null"))
    '''
    result = subprocess.run(
        ["Rscript", "-e", expression],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(result.stdout)
    exports = set(payload["exports"])
    named = {
        name
        for module in MODULES
        for function, _ in module["functions"]
        for name in re.findall(r"sn_[A-Za-z0-9_]+(?=\(\))", function)
    }
    missing = sorted(name for name in named if name not in exports)
    if missing:
        raise RuntimeError("Curated functions are not exported: " + ", ".join(missing))
    plot_pairs = [(row["analysis_type"], row["view"]) for row in payload["plots"]]
    if not plot_pairs:
        raise RuntimeError("The result-plot registry is empty.")
    if len(plot_pairs) != len(set(plot_pairs)):
        raise RuntimeError("The result-plot registry contains duplicate type/view pairs.")
    defaults = [row["analysis_type"] for row in payload["plots"] if row["default"]]
    plot_types = {row["analysis_type"] for row in payload["plots"]}
    if set(defaults) != plot_types or len(defaults) != len(plot_types):
        raise RuntimeError("Each result type must expose exactly one default plot view.")
    return payload


def write_registry_snapshot(payload: dict) -> None:
    snapshot = {
        "generated_from": "live Shennong namespace and registries",
        "package_version": payload["version"],
        "exports": payload["exports"],
        "analysis_method_tasks": payload["methods"],
        "plot_views": payload["plots"],
    }
    path = ROOT / "output" / "shennong-visual-guide-data.json"
    path.write_text(json.dumps(snapshot, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def workflow_nodes() -> list[Node]:
    return [
        Node(
            "input", 40, 170, 285, 180, "Input 输入",
            ("10x / STARsolo / H5AD / qs2", "matrix / list / SummarizedExperiment", "sn_read()", "sn_list_10x_paths()"),
            COLORS["io"],
        ),
        Node(
            "qc", 40, 400, 285, 245, "Prepare & QC 预处理",
            ("sn_initialize_seurat_object()", "sn_add_qc_metrics()", "sn_filter_cells() / sn_filter_genes()", "sn_normalize_data()", "sn_find_doublets()", "sn_remove_ambient_contamination()"),
            COLORS["io"],
        ),
        Node(
            "object", 385, 355, 340, 300, "Shennong object",
            (
                "Seurat-centered mutable analysis object",
                "Assays + layers: RNA / ADT / spatial",
                "meta.data: cell + sample provenance",
                "reductions / graphs / images",
                "misc$shennong$results[type][result_id]",
                "commands + reproducibility provenance",
            ),
            COLORS["core"], stroke="#263A82",
        ),
        Node(
            "method_registry", 790, 55, 610, 90, "Method registry 方法发现",
            ("sn_list_methods(task)  |  sn_get_method_status(method, task)", "23 tasks; runtime, availability, defaults, citations"),
            COLORS["method"], dashed=True,
        ),
        Node(
            "analysis_container", 770, 170, 650, 855, "Analysis modules 分析模块",
            ("Each module consumes the object or a validated standalone input.",),
            "#EEF4FB", stroke="#94A3B8", title_color=COLORS["ink"], text_color=COLORS["muted"], container=True,
        ),
        Node(
            "cluster", 800, 255, 285, 205, "Cluster & Integrate",
            ("sn_run_cluster()", "sn_run_multimodal()", "Harmony / CCA / RPCA", "scVI / scANVI / scPoli", "WNN / totalVI / Coralysis"),
            COLORS["cluster"],
        ),
        Node(
            "annotation", 1105, 255, 285, 205, "Annotation, DE & Pathway",
            ("sn_run_annotation()", "sn_find_de()", "sn_run_enrichment()", "review confidence + ontology"),
            COLORS["annotation"],
        ),
        Node(
            "program", 800, 485, 285, 205, "Programs & Regulation",
            ("sn_score_programs()", "sn_test_programs()", "sn_discover_programs()", "sn_run_grn()", "sn_run_metabolism()"),
            COLORS["program"],
        ),
        Node(
            "composition", 1105, 485, 285, 205, "Composition, DA & CNV",
            ("sn_calculate_composition()", "sn_compare_composition()", "sn_test_abundance()", "sn_prioritize_states()", "sn_run_cnv()"),
            COLORS["composition"],
        ),
        Node(
            "dynamic", 800, 715, 285, 205, "Communication & Dynamics",
            ("sn_run_cell_communication()", "sn_run_trajectory()", "sn_run_velocity()", "sn_run_fate()"),
            COLORS["dynamic"],
        ),
        Node(
            "spatial", 1105, 715, 285, 205, "Spatial",
            ("sn_run_spatial(task=)", "SVG / domain", "neighborhood / mapping", "deconvolution / communication"),
            COLORS["spatial"],
        ),
        Node(
            "bulk", 385, 760, 340, 220, "Standalone bulk lane",
            ("matrix / list / SummarizedExperiment", "sn_run_bulk(workflow=)", "sn_assess_bulk_qc() / sn_find_de()", "sn_score_bulk_pathways() / sn_run_wgcna()", "sn_run_survival() / sn_run_clinical_association()"),
            COLORS["bulk"],
        ),
        Node(
            "result", 1480, 250, 365, 220, "Versioned result 结果契约",
            ("schema_version = 2.0.0", "analysis_type + result_id", "tables$primary + diagnostics", "provenance + backend", "sn_list_results() / sn_get_result()", "sn_store_result() / sn_validate_result()"),
            COLORS["result"],
        ),
        Node(
            "plot", 1480, 540, 365, 245, "Visualization 可视化",
            ("sn_list_plot_methods()", "26 analysis types / 78 views", "sn_plot_result(object,", "  analysis_type, result_id, view)", "direct plots: dim / feature / composition", "distribution / association"),
            COLORS["plot"],
        ),
        Node(
            "export", 1480, 850, 365, 175, "Validate & Export 输出",
            ("sn_apply_figure_profile()", "sn_validate_figure()", "sn_export_figure()", "sn_export_figure_bundle()"),
            COLORS["export"],
        ),
    ]


def workflow_edges() -> list[Edge]:
    return [
        Edge("input", "object", "read / initialize"),
        Edge("qc", "object", "update object"),
        Edge("object", "analysis_container", "analyze"),
        Edge("method_registry", "analysis_container", "select backend", dashed=True),
        Edge("analysis_container", "result", "store or return"),
        Edge("bulk", "result", "validated result"),
        Edge("result", "plot", "resolve type + id"),
        Edge("object", "plot", "metadata / stored result", dashed=True, color="#3751A6"),
        Edge("plot", "export", "profile + validate"),
    ]


def add_drawio_cell(root: ET.Element, node: Node, index: int) -> None:
    title_color = node.title_color
    body_color = node.text_color
    line_html = "<br>".join(html.escape(line) for line in node.lines)
    value = (
        f'<div style="font-size:16px;color:{title_color};"><b>{html.escape(node.title)}</b></div>'
        f'<div style="font-size:12px;line-height:1.45;color:{body_color};margin-top:8px;">{line_html}</div>'
    )
    stroke = node.stroke or node.fill
    opacity = "opacity=45;" if node.container else ""
    dashed = "dashed=1;dashPattern=8 6;" if node.dashed else ""
    style = (
        "rounded=1;whiteSpace=wrap;html=1;align=left;verticalAlign=top;spacing=16;"
        f"fillColor={node.fill};strokeColor={stroke};strokeWidth=2;{dashed}{opacity}"
        "fontFamily=Noto Sans CJK SC;shadow=0;"
    )
    cell = ET.SubElement(
        root,
        "mxCell",
        {"id": node.key, "value": value, "style": style, "vertex": "1", "parent": "1"},
    )
    ET.SubElement(
        cell,
        "mxGeometry",
        {"x": str(node.x), "y": str(node.y), "width": str(node.width), "height": str(node.height), "as": "geometry"},
    )


def write_drawio(nodes: list[Node], edges: list[Edge]) -> Path:
    mxfile = ET.Element(
        "mxfile",
        {"host": "app.diagrams.net", "modified": "2026-09-06T00:00:00.000Z", "agent": "Shennong generator", "version": "24.7.17", "type": "device"},
    )
    diagram = ET.SubElement(mxfile, "diagram", {"id": "shennong-object-workflow", "name": "Shennong workflow"})
    graph = ET.SubElement(
        diagram,
        "mxGraphModel",
        {
            "dx": "1900", "dy": "1100", "grid": "1", "gridSize": "10", "guides": "1", "tooltips": "1",
            "connect": "1", "arrows": "1", "fold": "1", "page": "1", "pageScale": "1",
            "pageWidth": "1900", "pageHeight": "1100", "math": "0", "shadow": "0",
        },
    )
    root = ET.SubElement(graph, "root")
    ET.SubElement(root, "mxCell", {"id": "0"})
    ET.SubElement(root, "mxCell", {"id": "1", "parent": "0"})
    header_value = (
        '<div style="font-size:24px;color:#FFFFFF;"><b>SHENNONG OBJECT-CENTERED WORKFLOW</b></div>'
        '<div style="font-size:14px;color:#DCE8F5;margin-top:6px;">输入 -&gt; 对象 -&gt; 分析模块 -&gt; versioned result -&gt; 可视化与交付</div>'
    )
    header = ET.SubElement(
        root,
        "mxCell",
        {
            "id": "header", "value": header_value,
            "style": "rounded=1;whiteSpace=wrap;html=1;align=left;verticalAlign=middle;spacingLeft=24;fillColor=#16324F;strokeColor=#16324F;fontFamily=Noto Sans CJK SC;",
            "vertex": "1", "parent": "1",
        },
    )
    ET.SubElement(header, "mxGeometry", {"x": "28", "y": "22", "width": "1844", "height": "86", "as": "geometry"})
    for index, node in enumerate(nodes, start=2):
        add_drawio_cell(root, node, index)
    for index, edge in enumerate(edges, start=100):
        dashed = "dashed=1;dashPattern=8 6;" if edge.dashed else ""
        style = (
            "edgeStyle=orthogonalEdgeStyle;rounded=1;orthogonalLoop=1;jettySize=auto;html=1;"
            f"strokeColor={edge.color};strokeWidth=3;endArrow=block;endFill=1;{dashed}"
            "fontFamily=Noto Sans CJK SC;fontSize=12;labelBackgroundColor=#FFFFFF;"
        )
        cell = ET.SubElement(
            root,
            "mxCell",
            {
                "id": f"edge-{index}", "value": html.escape(edge.label), "style": style,
                "edge": "1", "parent": "1", "source": edge.source, "target": edge.target,
            },
        )
        ET.SubElement(cell, "mxGeometry", {"relative": "1", "as": "geometry"})
    ET.indent(mxfile, space="  ")
    path = OUT_DIAGRAM / "shennong-object-workflow.drawio"
    ET.ElementTree(mxfile).write(path, encoding="utf-8", xml_declaration=True)
    return path


def svg_text(parent: ET.Element, x: float, y: float, text: str, size: int, fill: str, weight: str = "normal") -> None:
    el = ET.SubElement(
        parent,
        "text",
        {
            "x": str(x), "y": str(y), "font-size": str(size), "font-weight": weight,
            "fill": fill, "font-family": "Noto Sans CJK SC, Noto Sans, sans-serif",
        },
    )
    el.text = text


def svg_wrapped(parent: ET.Element, x: float, y: float, text: str, width: float, size: int, fill: str, line_height: float) -> float:
    max_chars = max(12, int(width / (size * 0.58)))
    lines = textwrap.wrap(text, width=max_chars, break_long_words=False, break_on_hyphens=False) or [""]
    for line in lines:
        svg_text(parent, x, y, line, size, fill)
        y += line_height
    return y


def write_svg(nodes: list[Node], edges: list[Edge]) -> Path:
    svg = ET.Element(
        "svg",
        {
            "xmlns": "http://www.w3.org/2000/svg", "width": "1900", "height": "1100",
            "viewBox": "0 0 1900 1100", "role": "img", "aria-labelledby": "title desc",
        },
    )
    title = ET.SubElement(svg, "title", {"id": "title"})
    title.text = "Shennong object-centered workflow"
    desc = ET.SubElement(svg, "desc", {"id": "desc"})
    desc.text = "Inputs and quality control feed a Shennong object, analysis modules produce versioned results, and the plot registry drives validated figure export."
    defs = ET.SubElement(svg, "defs")
    marker = ET.SubElement(defs, "marker", {"id": "arrow", "markerWidth": "10", "markerHeight": "10", "refX": "9", "refY": "3", "orient": "auto", "markerUnits": "strokeWidth"})
    ET.SubElement(marker, "path", {"d": "M0,0 L0,6 L9,3 z", "fill": COLORS["method"]})
    ET.SubElement(svg, "rect", {"x": "0", "y": "0", "width": "1900", "height": "1100", "fill": COLORS["paper"]})
    ET.SubElement(svg, "rect", {"x": "28", "y": "22", "width": "1844", "height": "86", "rx": "22", "fill": COLORS["ink"]})
    svg_text(svg, 58, 63, "SHENNONG OBJECT-CENTERED WORKFLOW", 28, COLORS["white"], "700")
    svg_text(svg, 58, 91, "输入 -> 对象 -> 分析模块 -> versioned result -> 可视化与交付", 18, "#DCE8F5")

    by_key = {node.key: node for node in nodes}
    for edge in edges:
        source = by_key[edge.source]
        target = by_key[edge.target]
        sx = source.x + source.width
        sy = source.y + source.height / 2
        tx = target.x
        ty = target.y + target.height / 2
        if target.y > source.y + source.height and abs(tx - sx) < 120:
            sx = source.x + source.width / 2
            sy = source.y + source.height
            tx = target.x + target.width / 2
            ty = target.y
        path = f"M {sx} {sy} C {(sx + tx) / 2} {sy}, {(sx + tx) / 2} {ty}, {tx} {ty}"
        attrs = {"d": path, "fill": "none", "stroke": edge.color, "stroke-width": "3", "marker-end": "url(#arrow)"}
        if edge.dashed:
            attrs["stroke-dasharray"] = "10 8"
        ET.SubElement(svg, "path", attrs)
        if edge.label:
            lx = (sx + tx) / 2
            ly = (sy + ty) / 2 - 7
            ET.SubElement(svg, "rect", {"x": str(lx - 55), "y": str(ly - 15), "width": "110", "height": "24", "rx": "10", "fill": "#FFFFFF", "opacity": "0.9"})
            label = ET.SubElement(svg, "text", {"x": str(lx), "y": str(ly + 2), "text-anchor": "middle", "font-size": "13", "fill": edge.color, "font-family": "Noto Sans CJK SC, sans-serif"})
            label.text = edge.label

    for node in nodes:
        stroke = node.stroke or node.fill
        attrs = {
            "x": str(node.x), "y": str(node.y), "width": str(node.width), "height": str(node.height),
            "rx": "18", "fill": node.fill, "stroke": stroke, "stroke-width": "2.5",
        }
        if node.container:
            attrs["fill-opacity"] = "0.62"
        if node.dashed:
            attrs["stroke-dasharray"] = "10 8"
        ET.SubElement(svg, "rect", attrs)
        title_y = node.y + 34
        svg_text(svg, node.x + 18, title_y, node.title, 18 if not node.container else 20, node.title_color, "700")
        current_y = title_y + 31
        for line in node.lines:
            current_y = svg_wrapped(svg, node.x + 20, current_y, line, node.width - 40, 14, node.text_color, 21)
            current_y += 3

    svg_text(svg, 40, 1080, "Generated from Shennong registries | editable source: .drawio | canonical plotting: sn_plot_result()", 14, COLORS["muted"])
    path = OUT_DIAGRAM / "shennong-object-workflow.svg"
    ET.ElementTree(svg).write(path, encoding="utf-8", xml_declaration=True)
    return path


def hex_color(value: str):
    return HexColor(value)


def wrap_for_pdf(text: str, font: str, size: float, max_width: float) -> list[str]:
    words = text.split()
    if not words:
        return [""]
    lines: list[str] = []
    current = words[0]
    for word in words[1:]:
        candidate = current + " " + word
        if pdfmetrics.stringWidth(candidate, font, size) <= max_width:
            current = candidate
        else:
            lines.append(current)
            current = word
    lines.append(current)
    return lines


def draw_pdf_node(c: canvas.Canvas, node: Node, scale: float, offset_x: float, offset_y: float, design_height: float) -> None:
    x = offset_x + node.x * scale
    y = offset_y + (design_height - node.y - node.height) * scale
    width = node.width * scale
    height = node.height * scale
    c.setFillColor(hex_color(node.fill))
    c.setStrokeColor(hex_color(node.stroke or node.fill))
    c.setLineWidth(1.2)
    if node.dashed:
        c.setDash(5, 4)
    else:
        c.setDash()
    c.roundRect(x, y, width, height, 10 * scale, fill=1, stroke=1)
    if node.container:
        c.setFillAlpha(1)
    c.setFillColor(hex_color(node.title_color))
    c.setFont("NotoCJKBold", max(8, 18 * scale))
    c.drawString(x + 16 * scale, y + height - 29 * scale, node.title)
    c.setFillColor(hex_color(node.text_color))
    c.setFont("NotoCJK", max(6.5, 13 * scale))
    cursor_y = y + height - 58 * scale
    body_size = max(6.5, 13 * scale)
    line_height = body_size * 1.32
    for line in node.lines:
        for wrapped in wrap_for_pdf(line, "NotoCJK", body_size, width - 32 * scale):
            c.drawString(x + 17 * scale, cursor_y, wrapped)
            cursor_y -= line_height
        cursor_y -= body_size * 0.28


def draw_pdf_edge(c: canvas.Canvas, edge: Edge, by_key: dict[str, Node], scale: float, offset_x: float, offset_y: float, design_height: float) -> None:
    source = by_key[edge.source]
    target = by_key[edge.target]
    sx = offset_x + (source.x + source.width) * scale
    sy = offset_y + (design_height - source.y - source.height / 2) * scale
    tx = offset_x + target.x * scale
    ty = offset_y + (design_height - target.y - target.height / 2) * scale
    if target.y > source.y + source.height and abs(target.x - source.x - source.width) < 120:
        sx = offset_x + (source.x + source.width / 2) * scale
        sy = offset_y + (design_height - source.y - source.height) * scale
        tx = offset_x + (target.x + target.width / 2) * scale
        ty = offset_y + (design_height - target.y) * scale
    c.setStrokeColor(hex_color(edge.color))
    c.setFillColor(hex_color(edge.color))
    c.setLineWidth(1.6)
    c.setDash(5, 4) if edge.dashed else c.setDash()
    middle_x = (sx + tx) / 2
    path = c.beginPath()
    path.moveTo(sx, sy)
    path.curveTo(middle_x, sy, middle_x, ty, tx, ty)
    c.drawPath(path, stroke=1, fill=0)
    angle = math.atan2(ty - ty, tx - middle_x) if tx != middle_x else 0
    arrow = 7
    c.line(tx, ty, tx - arrow, ty + arrow / 2)
    c.line(tx, ty, tx - arrow, ty - arrow / 2)
    if edge.label:
        label_size = 6.5
        label_width = pdfmetrics.stringWidth(edge.label, "NotoCJK", label_size) + 10
        lx = (sx + tx) / 2
        ly = (sy + ty) / 2
        c.setFillColor(white)
        c.roundRect(lx - label_width / 2, ly - 5, label_width, 12, 4, fill=1, stroke=0)
        c.setFillColor(hex_color(edge.color))
        c.setFont("NotoCJK", label_size)
        c.drawCentredString(lx, ly - 1, edge.label)


def write_workflow_pdf(nodes: list[Node], edges: list[Edge]) -> Path:
    path = OUT_PDF / "shennong-object-workflow.pdf"
    page_width, page_height = landscape(A3)
    design_width, design_height = 1900, 1100
    margin = 20
    scale = min((page_width - 2 * margin) / design_width, (page_height - 2 * margin) / design_height)
    offset_x = (page_width - design_width * scale) / 2
    offset_y = (page_height - design_height * scale) / 2
    c = canvas.Canvas(str(path), pagesize=(page_width, page_height), pageCompression=1)
    c.setTitle("Shennong object-centered workflow")
    c.setAuthor("Shennong")
    c.setFillColor(hex_color(COLORS["paper"]))
    c.rect(0, 0, page_width, page_height, fill=1, stroke=0)
    c.setFillColor(hex_color(COLORS["ink"]))
    c.roundRect(offset_x + 28 * scale, offset_y + (design_height - 108) * scale, 1844 * scale, 86 * scale, 14, fill=1, stroke=0)
    c.setFillColor(white)
    # Keep pure-Latin header text on PDF base fonts. Some Poppler builds can
    # omit repeated glyphs from a multi-page subsetted CJK font.
    c.setFont("Helvetica-Bold", 18)
    c.drawString(offset_x + 58 * scale, offset_y + (design_height - 62) * scale, "SHENNONG OBJECT-CENTERED WORKFLOW")
    c.setFillColor(hex_color("#DCE8F5"))
    c.setFont("NotoCJK", 10)
    c.drawString(offset_x + 58 * scale, offset_y + (design_height - 91) * scale, "输入 -> 对象 -> 分析模块 -> versioned result -> 可视化与交付")
    by_key = {node.key: node for node in nodes}
    for edge in edges:
        draw_pdf_edge(c, edge, by_key, scale, offset_x, offset_y, design_height)
    for node in nodes:
        draw_pdf_node(c, node, scale, offset_x, offset_y, design_height)
    c.setFillColor(hex_color(COLORS["muted"]))
    c.setFont("NotoCJK", 7)
    c.drawString(offset_x + 40 * scale, offset_y + 12 * scale, "Generated from live Shennong registries | editable source: shennong-object-workflow.drawio")
    c.showPage()
    c.save()
    return path


def draw_header(c: canvas.Canvas, page_number: int, subtitle: str, version: str) -> None:
    page_width, page_height = landscape(A4)
    c.setFillColor(hex_color(COLORS["ink"]))
    c.rect(0, page_height - 50, page_width, 50, fill=1, stroke=0)
    c.setFillColor(white)
    c.setFont("NotoCJKBold", 18)
    c.drawString(22, page_height - 29, "SHENNONG CHEAT SHEET")
    c.setFont("NotoCJK", 8.5)
    c.drawString(285, page_height - 28, subtitle)
    c.setFont("Courier", 7)
    c.drawRightString(page_width - 22, page_height - 19, f"v{version}")
    c.drawRightString(page_width - 22, page_height - 34, f"PAGE {page_number} / 2")


def draw_footer(c: canvas.Canvas, page_number: int) -> None:
    page_width, _ = landscape(A4)
    c.setStrokeColor(hex_color(COLORS["line"]))
    c.line(20, 17, page_width - 20, 17)
    c.setFillColor(hex_color(COLORS["muted"]))
    c.setFont("NotoCJK", 6.5)
    c.drawString(22, 7, "Discover first: sn_list_methods() + sn_list_plot_methods() | Results: analysis_type + result_id")
    c.drawRightString(page_width - 22, 7, f"Shennong object workflow | {page_number}")


def draw_code_strip(c: canvas.Canvas, y_top: float) -> float:
    page_width, _ = landscape(A4)
    x = 20
    width = page_width - 40
    height = 66
    y = y_top - height
    c.setFillColor(hex_color("#EAF0F8"))
    c.setStrokeColor(hex_color("#C8D4E3"))
    c.roundRect(x, y, width, height, 8, fill=1, stroke=1)
    c.setFillColor(hex_color(COLORS["ink"]))
    c.setFont("NotoCJKBold", 8.5)
    c.drawString(x + 11, y + height - 16, "QUICK START  最短主线")
    snippets = [
        "x <- sn_read(path)",
        "x <- sn_initialize_seurat_object(x, project = 'atlas')",
        "x <- sn_run_cluster(x, integration_method = 'harmony')",
        "x <- sn_run_annotation(x, group_by = 'cluster', result_id = 'labels')",
        "sn_list_results(x)",
        "sn_plot_result(x, analysis_type = 'de', result_id = 'markers', view = 'volcano')",
    ]
    c.setFont("GuideMono", 6.2)
    c.setFillColor(hex_color("#263B52"))
    col_width = (width - 30) / 2
    for index, snippet in enumerate(snippets):
        col = index // 3
        row = index % 3
        c.drawString(x + 11 + col * (col_width + 8), y + height - 31 - row * 11, snippet)
    return y - 8


def measure_text_lines(text: str, font: str, size: float, width: float) -> list[str]:
    return wrap_for_pdf(text, font, size, width)


def draw_panel(
    c: canvas.Canvas,
    x: float,
    y_top: float,
    width: float,
    title: str,
    color: str,
    items: Iterable[tuple[str, str]],
    footer: str | None = None,
    body_size: float = 6.35,
) -> float:
    items = list(items)
    inner_width = width - 16
    item_layout: list[tuple[list[str], list[str]]] = []
    content_height = 0.0
    for function, description in items:
        fn_lines = measure_text_lines(function, "GuideMono", body_size, inner_width)
        desc_lines = measure_text_lines(description, "NotoCJK", body_size, inner_width)
        item_layout.append((fn_lines, desc_lines))
        content_height += len(fn_lines) * 8.2 + len(desc_lines) * 8.0 + 4.0
    footer_lines = measure_text_lines(footer, "NotoCJK", 5.8, inner_width) if footer else []
    content_height += len(footer_lines) * 7.2 + (5 if footer_lines else 0)
    height = 27 + content_height + 8
    y = y_top - height
    c.setFillColor(white)
    c.setStrokeColor(hex_color("#D8E1EA"))
    c.roundRect(x, y, width, height, 7, fill=1, stroke=1)
    c.setFillColor(hex_color(color))
    c.roundRect(x, y + height - 24, width, 24, 7, fill=1, stroke=0)
    c.rect(x, y + height - 24, width, 12, fill=1, stroke=0)
    c.setFillColor(white)
    c.setFont("NotoCJKBold", 8)
    c.drawString(x + 8, y + height - 16, title)
    cursor = y + height - 35
    for fn_lines, desc_lines in item_layout:
        c.setFillColor(hex_color(color))
        c.setFont("GuideMono", body_size)
        for line in fn_lines:
            c.drawString(x + 8, cursor, line)
            cursor -= 8.2
        c.setFillColor(hex_color(COLORS["muted"]))
        c.setFont("NotoCJK", body_size)
        for line in desc_lines:
            c.drawString(x + 8, cursor, line)
            cursor -= 8.0
        cursor -= 4.0
    if footer_lines:
        c.setStrokeColor(hex_color("#E2E8F0"))
        c.line(x + 8, cursor + 2, x + width - 8, cursor + 2)
        cursor -= 6
        c.setFillColor(hex_color(COLORS["muted"]))
        c.setFont("NotoCJK", 5.8)
        for line in footer_lines:
            c.drawString(x + 8, cursor, line)
            cursor -= 7.2
    return y - 8


def draw_page_one(c: canvas.Canvas, version: str) -> None:
    page_width, page_height = landscape(A4)
    c.setFillColor(hex_color(COLORS["paper"]))
    c.rect(0, 0, page_width, page_height, fill=1, stroke=0)
    draw_header(c, 1, "对象、预处理、整合、注释与结果契约", version)
    top = draw_code_strip(c, page_height - 58)
    margin, gap = 20, 8
    col_width = (page_width - 2 * margin - 2 * gap) / 3
    xs = [margin, margin + col_width + gap, margin + 2 * (col_width + gap)]

    y = top
    y = draw_panel(c, xs[0], y, col_width, MODULES[0]["title"], MODULES[0]["color"], MODULES[0]["functions"])
    y = draw_panel(c, xs[0], y, col_width, MODULES[1]["title"], MODULES[1]["color"], MODULES[1]["functions"], footer="关键：assay/layer 明确；样本和受试者 provenance 不要静默合并。")
    y = draw_panel(
        c, xs[0], y, col_width, "Verb grammar 命名语法", COLORS["method"],
        [
            ("sn_list_* / sn_get_*", "列举摘要 / 读取实体"),
            ("sn_run_* / sn_find_* / sn_score_*", "执行工作流 / 推断 / 打分"),
            ("sn_check_* / sn_validate_*", "诊断环境 / 断言契约"),
            ("sn_assess_* / sn_calculate_* / sn_plot_*", "QC 分析 / 指标 / 图形"),
        ],
        body_size=5.7,
    )

    y = top
    y = draw_panel(c, xs[1], y, col_width, MODULES[2]["title"], MODULES[2]["color"], MODULES[2]["functions"], footer="先用 sn_list_methods(task) 查看默认后端和当前可用性。")
    y = draw_panel(c, xs[1], y, col_width, MODULES[3]["title"], MODULES[3]["color"], MODULES[3]["functions"])
    y = draw_panel(
        c, xs[1], y, col_width, "Method controls 后端控制", COLORS["method"],
        [
            ("sn_list_methods(task)", "方法、默认值、runtime 与 availability"),
            ("sn_get_method_status(method, task)", "依赖、安装动作、输入和引用"),
            ("sn_get_integration_control_template()", "整合后端的完整 control surface"),
        ],
        body_size=5.7,
    )

    y = top
    y = draw_panel(
        c, xs[2], y, col_width, "Versioned result 统一结果", COLORS["result"],
        [
            ("sn_list_results(object)", "查看 analysis_type / result_id"),
            ("sn_get_result(object, type, result_id)", "读取已存结果，不直接摸 misc"),
            ("sn_validate_result(result)", "校验 schema 和 primary table"),
            ("sn_store_result(object, type, result_id, result)", "写入规范存储"),
            ("sn_audit_results(object)", "审计旧结果和非注册 artifact"),
        ],
        footer="主表位于 tables$primary；结果身份由 analysis_type + result_id 决定。",
    )
    y = draw_panel(
        c, xs[2], y, col_width, "Object anatomy 对象结构", COLORS["core"],
        [
            ("Assays / layers", "RNA, ADT, spatial; counts/data/custom layers"),
            ("meta.data", "cell, sample, donor, condition 和 labels"),
            ("reductions / graphs / images", "UMAP/PCA/latent/邻接图/空间坐标"),
            ("misc$shennong$results", "版本化结果而不是临时 list"),
        ],
    )
    y = draw_panel(
        c, xs[2], y, col_width, "Reproducibility 复现要点", COLORS["cluster"],
        [
            ("result_id", "为每个分析结果命名，不覆盖旧结果"),
            ("assay + layer", "声明表达来源，尤其 counts 与 corrected layers"),
            ("sample_by + seed", "保留统计重复并固定随机性"),
        ],
        body_size=5.7,
    )
    draw_footer(c, 1)


def draw_page_two(c: canvas.Canvas, version: str, plots: list[dict]) -> None:
    page_width, page_height = landscape(A4)
    c.setFillColor(hex_color(COLORS["paper"]))
    c.rect(0, 0, page_width, page_height, fill=1, stroke=0)
    draw_header(c, 2, "扩展分析、统一可视化与出版输出", version)
    top = page_height - 60
    margin, gap = 20, 8
    col_width = (page_width - 2 * margin - 2 * gap) / 3
    xs = [margin, margin + col_width + gap, margin + 2 * (col_width + gap)]

    y = top
    y = draw_panel(c, xs[0], y, col_width, MODULES[5]["title"], MODULES[5]["color"], MODULES[5]["functions"], footer="sample-level 比较先按 sample 聚合；cell 不是生物学重复。")
    y = draw_panel(c, xs[0], y, col_width, MODULES[4]["title"], MODULES[4]["color"], MODULES[4]["functions"])
    y = draw_panel(
        c, xs[0], y, col_width, "Accepted inputs 输入契约", COLORS["core"],
        [
            ("Seurat", "单细胞、多模态与空间主入口"),
            ("matrix / list / SummarizedExperiment", "standalone bulk 入口"),
            ("sn_analysis_result", "直接进入结果可视化"),
            ("data.frame", "兼容 composition、distribution 和 association"),
        ],
        body_size=5.7,
    )

    y = top
    y = draw_panel(c, xs[1], y, col_width, MODULES[6]["title"], MODULES[6]["color"], MODULES[6]["functions"])
    y = draw_panel(c, xs[1], y, col_width, MODULES[7]["title"], MODULES[7]["color"], MODULES[7]["functions"])
    y = draw_panel(c, xs[1], y, col_width, MODULES[8]["title"], MODULES[8]["color"], MODULES[8]["functions"], body_size=5.9)

    plot_by_type: dict[str, list[str]] = {}
    defaults: dict[str, str] = {}
    for row in plots:
        plot_by_type.setdefault(row["analysis_type"], []).append(row["view"])
        if row["default"]:
            defaults[row["analysis_type"]] = row["view"]
    selected = ["de", "enrichment", "annotation", "cell_communication", "trajectory", "cnv", "metabolism", "bulk_qc"]
    plot_items = []
    for analysis_type in selected:
        views = ", ".join(plot_by_type[analysis_type])
        plot_items.append((analysis_type, f"{views} | default: {defaults[analysis_type]}"))

    y = top
    y = draw_panel(
        c, xs[2], y, col_width, "Canonical plotting 统一绘图", COLORS["plot"],
        [
            ("sn_list_plot_methods(type)", "发现 accepted_input, default, required_parameters"),
            ("sn_plot_result(object, ...)", "统一 object / analysis_type / result_id / view"),
            ("sn_plot_distribution()", "violin/box/histogram/density/ridge"),
            ("sn_plot_association()", "cell/sample numeric association scatter"),
        ],
        footer="Seurat 中有多个同类结果时显式给 result_id；只有 unique/default 才自动选择。",
    )
    y = draw_panel(c, xs[2], y, col_width, "Common result views 常用视图", COLORS["annotation"], plot_items, body_size=5.35)
    y = draw_panel(
        c, xs[2], y, col_width, "Publication 出版输出", COLORS["export"],
        [
            ("sn_list_figure_profiles()", "screen/column/page/slide"),
            ("sn_apply_figure_profile()", "应用字体、线宽和版面约束"),
            ("sn_validate_figure()", "出版前结构化检查"),
            ("sn_export_figure()", "PDF/SVG/TIFF/PNG"),
            ("sn_export_figure_bundle()", "figure + source data + manifest"),
        ],
        body_size=5.8,
    )
    draw_footer(c, 2)


def write_cheatsheet(payload: dict) -> Path:
    path = OUT_PDF / "shennong-cheatsheet.pdf"
    c = canvas.Canvas(str(path), pagesize=landscape(A4), pageCompression=1)
    c.setTitle("Shennong cheat sheet")
    c.setAuthor("Shennong")
    c.setSubject("Object-centered workflow, analysis methods, visualization, and publication output")
    draw_page_one(c, payload["version"])
    c.showPage()
    draw_page_two(c, payload["version"], payload["plots"])
    c.showPage()
    c.save()
    return path


def main() -> None:
    OUT_PDF.mkdir(parents=True, exist_ok=True)
    OUT_DIAGRAM.mkdir(parents=True, exist_ok=True)
    register_fonts()
    payload = collect_registry()
    write_registry_snapshot(payload)
    nodes = workflow_nodes()
    edges = workflow_edges()
    artifacts = [
        write_drawio(nodes, edges),
        write_svg(nodes, edges),
        write_workflow_pdf(nodes, edges),
        write_cheatsheet(payload),
    ]
    for artifact in artifacts:
        print(artifact.relative_to(ROOT))


if __name__ == "__main__":
    main()
