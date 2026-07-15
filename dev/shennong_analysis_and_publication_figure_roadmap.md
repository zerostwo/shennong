# Shennong 分析能力补全与发表级可视化实施规范

> 面向 Codex 的开发路线图  
> 版本：2026-07-15  
> 审查对象：`zerostwo/shennong`  
> 对照项目：`YuLab-SMU/sclet`、`mengxu98/scop`、`omicverse/omicverse`、`tidyplots`

------------------------------------------------------------------------

## 0. 文档范围

本阶段只完成两条主线：

1.  **补齐 Shennong 的组学分析能力**  
    参考 sclet、scop 和 OmicVerse 的分析方法，审查 Shennong
    当前已经具备的能力、仍然缺失的主分析环节，以及应该优先增加的后端。

2.  **建立 Shennong 的发表级图形系统**  
    在现有 `sn_plot_*()`、自动点大小、调色板、栅格化和 panel
    尺寸基础上，建设一套根据数据规模、分组数量、panel
    数量和输出目标自动确定图形尺寸与排版的 Figure Engine。

本阶段**不优先实现**：

- MCP server；
- REST API；
- 完整 Agent 产品；
- 新的统一生物对象；
- Shennong DB 的服务端功能；
- Web UI。

但是新增分析函数必须保持 **Agent-ready**：

- 输入参数可序列化；
- 输出结构明确；
- 无交互式菜单；
- 可记录 provenance；
- 可以映射为未来的 JSON Schema。

------------------------------------------------------------------------

# 1. 核心结论

Shennong 不应通过“包装最多的方法”与 scop 或 OmicVerse 竞争。

Shennong 的差异化应当是：

> **为真实生物信息学项目提供有明确默认方法、完整分析闭环、统一结果语义和发表级输出的
> R 工作流。**

建议形成以下产品特征：

``` text
可靠的数据接入
        +
经过选择的分析主线
        +
统一结果与诊断
        +
发表级智能可视化
        +
可选的 R/Python 后端
```

而不是：

``` text
数百个互不关联的 RunXXX() wrapper
```

sclet 值得借鉴的是：

- 以真实研究问题组织分析 mainline；
- analysis-state contract；
- 每个结果可查询、可审计；
- Python 环境隔离；
- workflow shell 将独立方法串联成完整路径。

scop 值得借鉴的是：

- 单细胞和空间方法覆盖广；
- Seurat v5 实战兼容；
- 多种 annotation、trajectory、communication 和 spatial backend；
- 大量针对结果的专用图形。

OmicVerse 值得借鉴的是：

- bulk、single-cell、spatial 和 multi-omics 的模块化分区；
- function registry；
- lazy loading；
- annotation、trajectory、GRN、metabolism、CNV、drug response
  和报告能力；
- Python 重型后端的可选加载。

Shennong
不应复制它们的全部方法，而应选择少量高价值默认后端，并为每一条主线提供：

``` text
validate → run → diagnose → store → plot → interpret
```

------------------------------------------------------------------------

# 2. 审查依据

本次比较基于 2026-07-15 可见的仓库接口，而不是仅依据项目宣传。

## 2.1 Shennong 当前接口

当前 `NAMESPACE` 已经导出以下主要能力：

- 数据读取与写出：10x、10x
  Spatial、H5、H5AD、BPCells、qs/qs2、STARsolo；
- 数据下载、Zenodo 索引及 Shennong Data Server；
- 细胞和基因过滤；
- doublet detection；
- ambient RNA correction；
- cell-cycle scoring；
- Seurat 聚类和多种 integration backend；
- scVI、scANVI、scArches、scPoli；
- integration assessment；
- CellTypist 和 label transfer；
- marker/DE、pseudobulk、enrichment；
- composition、RO/E、Milo；
- regulatory activity；
- CellChat、LIANA/NicheNet/CellPhoneDB 类通信分析；
- BayesPrism/CIBERSORTx 类 bulk deconvolution；
- InferCNVpy；
- cell2location、Tangram、Squidpy、stLearn、SpatialData；
- scDesign3 simulation；
- ellmer-based interpretation；
- 基础发表风格图形。

这说明 Shennong
已经不是一个早期空壳，下一阶段应当补齐主分析链条，而不是重写已有功能。

## 2.2 sclet 的方法组织

sclet 基于 `SingleCellExperiment`，将分析组织为 11 条 mainline：

1.  数据接入与互操作；
2.  核心单细胞流程；
3.  整合；
4.  参考映射和注释；
5.  trajectory、velocity 和 fate；
6.  gene program、regulon 和 mechanistic interpretation；
7.  DE、pseudobulk 和 enrichment；
8.  cell-cell communication；
9.  perturbation priority 和 rare-cell detection；
10. spatial deconvolution、colocalization 和 niche；
11. multimodal expansion。

它还提供：

- `AuditAnalysisChain()`；
- `PipelineSummary()`；
- `CommandLog()`；
- `get_*()` / `has_*()` 结果访问器；
- `SummarizeContextForLLM()`；
- workflow-level 函数。

值得借鉴的是“状态和主线”，而不是它的 Seurat 风格函数命名。

## 2.3 scop 的方法覆盖

scop 的优势是极广的单细胞和空间方法覆盖，包括：

- 多种 QC、doublet、decontamination；
- 多种降维和整合；
- CellTypist、SingleR、SciBet、scmap、Symphony 等注释；
- CytoTRACE/FitDevo；
- Slingshot、Monocle、PAGA、Palantir、CellRank、WOT、scVelo；
- GSVA、GSEA、metabolism、scFEA；
- DoRothEA、SCENIC/SCENIC+ 和多种 GRN；
- CellChat、CellPhoneDB、LIANA、NicheNet、MultiNicheNet；
- Milo、propeller、scCODA；
- Augur、Scissor、RareQ、in-silico perturbation；
- 广泛的 spatial domain、SVG、deconvolution、neighborhood、integration
  和 communication 方法；
- bulk DE、deconvolution 和组成分析。

它的风险也很明显：公共 API 极大，方法选择成本和依赖维护成本高。

## 2.4 OmicVerse 的方法覆盖

OmicVerse 当前将功能分为：

- bulk；
- single；
- space；
- bulk2single；
- metabol；
- protein；
- genetics；
- AIRR；
- epigenomics；
- agent/MCP/report。

其中单细胞模块包含：

- SCSA、MetaTiME、GPT annotation、CellVote 和 ontology mapping；
- VIA、Monocle-like trajectory、pseudotime fate、velocity 和 CellRank；
- SCENIC、GRN、NMF/cNMF、Hotspot；
- Milo、Augur、perturbation 和 dynamic features；
- CNV；
- metabolism 和 metabolite communication；
- CellPhoneDB、LIANA；
- MOFA、GLUE、SIMBA、TOSICA；
- drug-response prediction。

bulk 模块包含：

- DESeq2-like DE；
- WGCNA；
- GSEA；
- PPI/STRING；
- TCGA；
- batch correction；
- deconvolution。

spatial 模块包含：

- STAGATE、CellCharter 和 CAST；
- STAligner；
- SVG/Moran’s I；
- spatial trajectory；
- Tangram；
- spatial deconvolution；
- H&E 到表达预测。

------------------------------------------------------------------------

# 3. Shennong 当前能力矩阵

图例：

- **完整**：已有清晰用户接口和结果工作流；
- **部分**：已有后端或单点函数，但尚未形成完整分析主线；
- **缺失**：当前公共 API 中没有明确实现；
- **延后**：不应在当前里程碑优先实现。

| 分析能力 | Shennong | sclet | scop | OmicVerse | Shennong 决策 |
|----|---:|---:|---:|---:|----|
| 多格式数据接入 | 完整 | 完整 | 部分 | 完整 | 保持并强化验证 |
| 单细胞 QC | 完整 | 完整 | 完整 | 完整 | 增加统一 QC report |
| doublet / ambient RNA | 完整 | 完整 | 完整 | 部分 | 保持 |
| cell-cycle | 完整 | 部分 | 完整 | 完整 | 增加连续 cell-cycle backend 可选项 |
| clustering | 完整 | 完整 | 完整 | 完整 | 保持 |
| integration | 完整但偏后端化 | Harmony/scVI | 极广 | 广 | 增加统一 method registry 和诊断 |
| integration assessment | 完整 | 状态审计 | LISI 等 | 部分 | Shennong 优势，继续强化 |
| annotation：CellTypist/label transfer | 部分完整 | 完整 | 完整 | 完整 | 补齐 consensus、confidence、ontology |
| SingleR/scmap/Symphony | 缺失 | 有 | 有 | 部分 | P1 增加 |
| marker / cell-level DE | 完整 | 完整 | 完整 | 完整 | 保持 |
| pseudobulk DE | 完整 | 完整 | 完整 | 有 | 补齐复杂 design 和 repeated measures |
| standalone bulk DE | 部分 | 不主打 | 完整 | 完整 | P2 增加完整 bulk workflow |
| enrichment / GSEA | 完整 | 完整 | 完整 | 完整 | 统一结果与图形 |
| per-cell gene-set scoring | 缺失/不完整 | AUCell/GSVA/UCell | 完整 | AUCell 等 | P1 增加 |
| TF/pathway activity | 部分完整 | program workflow | 完整 | 完整 | 保留 decoupleR，补齐统一 scoring |
| GRN / SCENIC | 缺失 | 有 | 完整 | 完整 | P1 增加可选后端 |
| NMF/cNMF/Hotspot programs | 缺失 | 部分 | NMF | 完整 | P2 选择性增加 |
| trajectory / pseudotime | 缺失 | 完整 | 完整 | 完整 | P1 必须增加 |
| RNA velocity | 缺失 | 完整 | 完整 | 完整 | P2 Python 后端 |
| fate probability / CellRank | 缺失 | 完整 | 完整 | 完整 | P2 Python 后端 |
| dynamic genes / tradeSeq | 缺失 | 部分 | 完整 | 完整 | P1/P2 增加 |
| composition table | 完整 | 完整 | 完整 | 完整 | 保持 |
| differential abundance：Milo | 完整 | 有 | 有 | 有 | 保持 |
| propeller / scCODA | 缺失 | 部分 | 有 | 部分 | P1 增加 propeller，scCODA 可选 |
| perturbation priority / Augur | 缺失 | 有 | 有 | 有 | P2 增加 |
| phenotype-associated states / Scissor | 缺失 | 部分 | 有 | 部分 | P2 可选 |
| rare-cell detection | 部分 | 有 | RareQ | 部分 | 统一输出，增加 RareQ 可选 |
| CellChat / LIANA / NicheNet / CPDB | 部分完整 | 统一 CCI | 完整 | 完整 | 统一 result schema 和 consensus |
| MultiNicheNet | 缺失 | 无明确主线 | 有 | 部分 | P2 增加 |
| InferCNV | 完整后端 | 无主线 | 有 | 有 | 增加诊断和比较 |
| metabolism | 缺失 | gene programs | scMetabolism/scFEA | 完整 | P2 增加 |
| metacell | 缺失 | SuperCell | 有 | 有 | P3 可选 |
| spatial deconvolution/mapping | 部分完整 | 有 | 完整 | 完整 | 保持并统一 |
| spatial domain | 部分后端 | niche | 完整 | 完整 | P1 spatial 增加 |
| spatially variable genes | 部分后端 | 有 | 完整 | 完整 | P1 spatial 增加 |
| spatial neighborhood | 部分后端 | 有 | 完整 | 完整 | P1 spatial 增加 |
| spatial integration/alignment | 缺少统一主线 | 部分 | 完整 | 完整 | P2 spatial 增加 |
| spatial communication | 缺少统一主线 | CCI | 完整 | 部分 | P2 spatial 增加 |
| CITE-seq | 部分：WNN/totalVI/Coralysis | 预留 | 完整 | 多模态 | 整理成独立主线 |
| scATAC / RNA+ATAC | 缺失 | 预留 | 完整 | GLUE/ATAC | P3 |
| MOFA / general multi-omics | 缺失 | 预留 | 部分 | 完整 | P3 |
| WGCNA | 缺失 | 无 | 部分 | 完整 | P2 bulk |
| PPI/STRING | 缺失 | 无 | 部分 | 完整 | P3 bulk |
| survival / clinical association | 缺失 | 无 | 部分 | TCGA | P2 bulk |
| proteomics/metabolomics workflow | 缺失 | 无 | 无主线 | 有模块 | P4，暂不宣称已支持 |
| publication-aware figure sizing | 部分 | 基础图形 | 图形种类多 | 图形种类多 | Shennong 核心差异化 |
| journal export / figure QA | 缺失 | 缺失 | 缺失 | 缺失 | Shennong 核心差异化 |

------------------------------------------------------------------------

# 4. 第一条主线：补齐 Shennong 分析能力

## 4.1 开发原则

### 4.1.1 一个研究问题对应一个稳定入口

不要为每个后端创建一个同等地位的顶层函数。

推荐：

``` r

sn_run_annotation(
  object,
  method = "singleR"
)

sn_run_trajectory(
  object,
  method = "slingshot"
)

sn_score_programs(
  object,
  method = "ucell"
)
```

不推荐继续扩张：

``` r

sn_run_singleR()
sn_run_scmap()
sn_run_symphony()
sn_run_slingshot()
sn_run_monocle3()
sn_run_palantir()
```

底层专用函数可以作为 internal adapter：

``` text
.sn_annotation_singleR()
.sn_annotation_celltypist()
.sn_trajectory_slingshot()
.sn_trajectory_cellrank()
```

### 4.1.2 默认后端必须有明确选择

Shennong 的价值是帮助用户做合理选择，而不是把后端列表交还给用户。

建议默认：

| 问题 | 默认后端 | 备选后端 |
|----|----|----|
| 单细胞 annotation | SingleR + marker evidence | CellTypist、Seurat transfer、scANVI |
| reference mapping | Seurat transfer | Symphony、scmap、scArches |
| per-cell signature scoring | UCell | AUCell、singscore |
| bulk/pseudobulk pathway scoring | GSVA | ssGSEA、PLAGE、z-score |
| trajectory | Slingshot | Monocle3、Palantir |
| dynamic genes | tradeSeq | GAM/loess |
| RNA velocity | scVelo | velocyto-derived backends |
| fate | CellRank | pseudotime terminal-state approximation |
| differential abundance | propeller | Milo、scCODA |
| CCI consensus | LIANA | CellChat、CellPhoneDB |
| ligand-target inference | NicheNet | MultiNicheNet |
| TF/pathway activity | decoupleR + CollecTRI/PROGENy | DoRothEA |
| GRN | SCENIC/pySCENIC optional | GENIE3/GRNBoost2 |
| spatial SVG | nnSVG 或 Moran’s I | SPARK-X |
| spatial domain | BANKSY | BayesSpace、stLearn |
| spatial deconvolution | cell2location | RCTD、CARD |
| bulk DE | edgeR/DESeq2 按 design 选择 | limma-voom、dream |

### 4.1.3 不切换 Python 环境管理方案

sclet 使用 `basilisk` 的思想值得借鉴，但 Shennong 已经决定使用 pixi 管理
Python 后端。

本阶段应保持：

``` text
R public API
    ↓
structured input files
    ↓
pixi environment
    ↓
Python script
    ↓
structured artifacts
```

不要同时维护：

- pixi；
- basilisk；
- reticulate-managed conda；
- 用户全局 Python。

### 4.1.4 重型依赖继续放在 Suggests 或 pixi

核心安装必须保持轻量。

每个可选方法必须提供：

``` r

sn_method_status("cellrank")
```

返回：

``` text
Method: cellrank
Available: FALSE
Runtime: pixi
Environment: trajectory
Install action: sn_prepare_pixi_environment("trajectory")
```

------------------------------------------------------------------------

## 4.2 先建立最小分析结果契约

不要求本阶段创建新的 ShennongObject，但必须统一新增方法的结果结构。

建议内部结构：

``` r

result <- list(
  schema_version = "1.0",
  analysis_type = "trajectory",
  name = "cd8_slingshot",
  method = "slingshot",
  backend = "slingshot",
  input = list(
    assay = "RNA",
    layer = "data",
    reduction = "umap",
    cells = 12543L,
    features = 2000L
  ),
  parameters = list(),
  tables = list(),
  embeddings = list(),
  graphs = list(),
  models = list(),
  diagnostics = list(),
  warnings = character(),
  provenance = list(
    package_versions = list(),
    random_seed = 1L,
    timestamp = ""
  )
)
```

继续存入当前集中注册的 `object@misc` result collections。

必须新增或统一：

``` r

sn_store_result(object, type, name, result)
sn_get_result(object, type, name)
sn_list_results(object, type = NULL)
sn_delete_result(object, type, name)
sn_validate_result(result)
```

现有的：

``` r

sn_get_de_result()
sn_get_enrichment_result()
sn_get_milo_result()
```

可以暂时保留，但底层必须读取同一 schema。

------------------------------------------------------------------------

# 5. P1：补齐单细胞分析主线

P1 完成后，Shennong 应当能够从 count matrix
一直完成到机制解释和论文图，而不需要用户切换到另一个综合包。

## 5.1 Annotation 与 reference mapping

### 当前基础

已有：

- [`sn_run_celltypist()`](https://songqi.org/shennong/dev/reference/sn_run_celltypist.md)；
- [`sn_transfer_labels()`](https://songqi.org/shennong/dev/reference/sn_transfer_labels.md)；
- scANVI；
- scArches；
- scPoli；
- annotation evidence；
- LLM interpretation。

### 需要新增

``` r

sn_run_annotation(
  object,
  group_by = "seurat_clusters",
  method = c(
    "consensus",
    "singleR",
    "celltypist",
    "seurat",
    "symphony",
    "scmap",
    "scanvi"
  ),
  reference = NULL,
  tissue = NULL,
  disease = NULL,
  species = NULL,
  ontology = TRUE,
  store_name = "annotation"
)
```

### `method = "consensus"` 的建议逻辑

``` text
cluster markers
    +
reference prediction
    +
canonical marker database
    +
negative marker conflicts
    +
tissue/disease context
    ↓
candidate labels
    ↓
confidence calibration
    ↓
hierarchical label
```

标准输出至少包括：

- cell-level prediction；
- cluster-level prediction；
- level 1/2/3 label；
- prediction score；
- second-best label；
- margin；
- supporting markers；
- conflicting markers；
- reference coverage；
- ontology identifier；
- low-confidence flag。

### 必须增加的函数

``` r

sn_annotation_consensus()
sn_annotation_confidence()
sn_map_cell_ontology()
sn_review_annotation()
sn_plot_annotation_confidence()
sn_plot_annotation_markers()
sn_plot_annotation_confusion()
```

### 验收标准

- 支持 cluster-level 和 cell-level；
- 无 reference 时可基于 markers 工作；
- 有 reference 时保留每种 backend 的原始预测；
- 不允许 LLM 直接覆盖计算结果；
- LLM 只能对 evidence packet 进行解释和候选排序；
- 所有最终标签均可追溯。

------------------------------------------------------------------------

## 5.2 Gene set、program 和 regulatory analysis

### 当前基础

已有：

- signature registry；
- enrichment；
- regulatory activity；
- decoupleR/DoRothEA/PROGENy 依赖。

### 当前缺口

缺少统一的：

- per-cell scoring；
- per-sample scoring；
- program discovery；
- program comparison；
- regulon network；
- program-level plotting。

### 新增统一入口

``` r

sn_score_programs(
  object,
  signatures,
  method = c("ucell", "aucell", "gsva", "ssgsea", "mean"),
  assay = NULL,
  layer = "data",
  name = NULL
)
```

``` r

sn_discover_programs(
  object,
  method = c("nmf", "cnmf", "hotspot"),
  n_programs = "auto",
  group_by = NULL,
  name = NULL
)
```

``` r

sn_run_grn(
  object,
  method = c("scenic", "pyscenic", "genie3", "grnboost2"),
  name = NULL
)
```

### 方法优先级

P1：

- UCell；
- GSVA；
- AUCell；
- decoupleR activity；
- program differential testing。

P2：

- NMF/cNMF；
- Hotspot；
- SCENIC/pySCENIC；
- regulon specificity；
- dynamic program analysis。

### 输出

- score matrix；
- per-cell/per-sample metadata；
- program gene weights；
- program activity contrasts；
- regulon activity；
- network edges；
- stability diagnostics。

------------------------------------------------------------------------

## 5.3 Trajectory、pseudotime 和 dynamic genes

这是当前最大的单细胞功能缺口。

### 新增统一入口

``` r

sn_run_trajectory(
  object,
  method = c("slingshot", "monocle3", "palantir"),
  reduction = NULL,
  cluster_by = NULL,
  start = NULL,
  end = NULL,
  lineages = NULL,
  store_name = "trajectory"
)
```

### P1 默认实现

- Slingshot；
- pseudotime；
- lineage assignment；
- tradeSeq dynamic genes；
- lineage-specific smooth curves；
- branch comparison。

### P2 Python 扩展

``` r

sn_run_velocity(
  object,
  method = "scvelo",
  ...
)

sn_run_fate(
  object,
  method = "cellrank",
  ...
)
```

### 必须提供的结果

- pseudotime；
- lineage probability；
- principal curves/graph；
- terminal states；
- dynamic genes；
- branch-specific genes；
- fitted expression trends；
- diagnostics。

### 必须提供的图

- embedding + trajectory；
- pseudotime embedding；
- lineage probability；
- dynamic heatmap；
- gene trend；
- branch comparison；
- velocity stream；
- terminal-state probability。

------------------------------------------------------------------------

## 5.4 Differential abundance 与 perturbation priority

### 当前基础

已有：

- composition；
- compare composition；
- RO/E；
- Milo；
- rare-cell detection。

### 需要增加

``` r

sn_test_abundance(
  object,
  method = c("propeller", "milo", "scCODA", "permutation"),
  sample_by,
  condition_by,
  cell_type_by,
  design = NULL
)
```

推荐：

- `propeller` 作为 cell-type proportion 的默认样本级方法；
- `Milo` 用于 neighborhood-level differential abundance；
- `scCODA` 作为可选 Python compositional backend；
- permutation 作为简单验证方法。

### Perturbation

``` r

sn_prioritize_states(
  object,
  method = c("augur", "scissor", "rareQ"),
  phenotype,
  sample_by = NULL,
  ...
)
```

输出：

- state ranking；
- AUC/score；
- permutation null；
- phenotype association；
- uncertainty；
- sample contribution。

------------------------------------------------------------------------

## 5.5 Cell-cell communication 统一化

Shennong 已经具备多个后端，但当前需要从“能运行”升级为“结果可比较”。

### 推荐入口

``` r

sn_run_cell_communication(
  object,
  method = c("liana", "cellchat", "cellphonedb", "nichenet", "multinichenet"),
  sender = NULL,
  receiver = NULL,
  condition_by = NULL,
  sample_by = NULL,
  consensus = TRUE,
  store_name = "communication"
)
```

### 统一结果 schema

至少统一：

- source；
- target；
- ligand；
- receptor；
- score；
- p/q value；
- rank；
- method；
- condition；
- sample；
- pathway；
- target genes；
- evidence source。

### 必须增加

- method concordance；
- consensus ranking；
- condition comparison；
- sample-aware pseudobulk communication；
- NicheNet ligand-target matrix；
- MultiNicheNet；
- spatial distance constraint 的预留字段。

------------------------------------------------------------------------

## 5.6 CNV、malignancy 和 metabolism

### CNV

已有 InferCNVpy backend，但需要补齐：

``` r

sn_run_cnv(
  object,
  method = c("infercnvpy", "copykat"),
  reference_cells,
  genome = NULL,
  store_name = "cnv"
)
```

新增：

- CNV quality diagnostics；
- malignant score；
- subclone assignment；
- sample-aware CNV；
- chromosome heatmap；
- CNV UMAP；
- CNV-expression association。

### Metabolism

``` r

sn_run_metabolism(
  object,
  method = c("geneset", "scmetabolism", "scfea", "compass"),
  ...
)
```

P2 先实现：

- curated metabolic signatures；
- UCell/GSVA-based activity；
- differential metabolic activity。

scFEA/Compass 作为可选重型后端，不应作为默认。

------------------------------------------------------------------------

# 6. P2：补齐空间转录组主线

Shennong 已有多个 Python
backend，但需要一个不依赖用户理解后端细节的空间主线。

## 6.1 统一入口

``` r

sn_run_spatial(
  object,
  task = c(
    "qc",
    "svg",
    "domain",
    "neighborhood",
    "deconvolution",
    "mapping",
    "integration",
    "communication"
  ),
  method = "auto",
  ...
)
```

也可以保留更明确的用户函数：

``` r

sn_find_spatial_features()
sn_find_spatial_domains()
sn_run_spatial_neighborhood()
sn_run_spatial_deconvolution()
sn_run_spatial_mapping()
sn_integrate_spatial()
sn_run_spatial_communication()
```

## 6.2 优先方法

### SVG

- Moran’s I / Squidpy；
- nnSVG；
- SPARK-X 可选。

### Spatial domain

- BANKSY 默认；
- BayesSpace；
- stLearn；
- CellCharter 可选。

### Neighborhood

- Squidpy graph；
- neighborhood enrichment；
- co-occurrence；
- ligand-receptor distance filtering。

### Deconvolution

- cell2location；
- RCTD；
- CARD；
- Tangram 归类为 mapping，而非严格 deconvolution。

### Spatial integration

- sample alignment；
- shared latent space；
- STAligner 可选；
- spatial graph-aware batch comparison。

### Spatial communication

- SpatialCellChat；
- spatially constrained LIANA；
- distance-aware NicheNet。

## 6.3 空间结果必须统一

不能只保存 backend 文件。

统一保存：

- coordinates；
- spatial graph；
- domains；
- SVG table；
- deconvolution proportions；
- cell mapping；
- neighborhood statistics；
- spatial interactions；
- image transform 和 scale factor；
- backend artifacts。

------------------------------------------------------------------------

# 7. P2：建立完整 bulk transcriptomics 主线

Shennong 的 DESCRIPTION 宣称支持多组学，但当前 bulk
主线尚不完整。应先完成 bulk RNA-seq，再考虑 proteomics/metabolomics。

## 7.1 数据结构

支持：

- count matrix + sample metadata；
- `SummarizedExperiment`；
- Shennong DB lazy table materialization；
- TCGA/GTEx/Toil 类表达矩阵。

## 7.2 统一入口

``` r

sn_run_bulk(
  object,
  workflow = c("qc", "de", "pathway", "network", "survival"),
  ...
)
```

明确函数：

``` r

sn_assess_bulk_qc()
sn_find_bulk_de()
sn_score_bulk_pathways()
sn_run_wgcna()
sn_run_survival()
sn_run_clinical_association()
```

## 7.3 DE 方法决策

``` text
raw integer counts
    ├── simple design → edgeR or DESeq2
    ├── complex design → edgeR GLM / DESeq2
    ├── repeated measures → dream
    └── normalized continuous expression → limma
```

函数应该接受：

``` r

sn_find_bulk_de(
  object,
  design = ~ batch + condition,
  contrast = c("condition", "tumor", "normal"),
  method = "auto"
)
```

## 7.4 必须补齐

- library size / expression distribution；
- sample PCA；
- sample correlation；
- outlier detection；
- design matrix validation；
- contrast validation；
- independent filtering；
- shrinkage；
- batch covariates；
- repeated measures；
- GSEA/GSVA；
- WGCNA；
- phenotype/module association；
- survival/Cox；
- immune deconvolution；
- pathway/TF activity；
- publication result bundle。

------------------------------------------------------------------------

# 8. P3：多模态与其他组学

只有在 scRNA、spatial 和 bulk 三条主线稳定后再实施。

## 8.1 CITE-seq

当前已有 WNN、totalVI 和 Coralysis 相关决策，下一步应独立成文档化主线：

``` r

sn_run_multimodal(
  object,
  modality = "cite_seq",
  method = c("wnn", "totalvi", "coralysis")
)
```

需要：

- RNA/ADT QC；
- CLR/DSB normalization；
- ADT background；
- joint embedding；
- modality weight；
- multimodal marker；
- protein/RNA concordance。

## 8.2 scATAC 和 RNA+ATAC

建议 adapter-first：

- Signac；
- ArchR；
- chromVAR；
- WNN；
- GLUE；
- MultiMAP。

不要自己重写 peak calling 或 motif engine。

## 8.3 General multi-omics

采用：

- `MultiAssayExperiment`；
- MOFA；
- SNF/DIABLO 可选；
- sample mapping；
- factor activity；
- clinical association。

## 8.4 暂时不要实现

- 大而全的 proteomics engine；
- 大而全的 metabolomics engine；
- 所有 OmicVerse genetics/AIRR/epigenomics 子模块；
- scop 中全部 niche 方法；
- 所有降维算法；
- 所有空间包 wrapper。

------------------------------------------------------------------------

# 9. 推荐的分析模块目录

``` text
R/
├── analysis-registry.R
├── analysis-result.R
├── analysis-provenance.R
│
├── annotation.R
├── annotation-consensus.R
├── annotation-ontology.R
│
├── program-scoring.R
├── program-discovery.R
├── regulatory-activity.R
├── grn.R
│
├── trajectory.R
├── trajectory-dynamics.R
├── velocity.R
├── fate.R
│
├── differential-abundance.R
├── perturbation-priority.R
├── rare-cells.R
│
├── communication.R
├── communication-consensus.R
│
├── cnv.R
├── metabolism.R
│
├── spatial-svg.R
├── spatial-domain.R
├── spatial-neighborhood.R
├── spatial-deconvolution.R
├── spatial-integration.R
├── spatial-communication.R
│
├── bulk-qc.R
├── bulk-de.R
├── bulk-pathway.R
├── bulk-network.R
├── bulk-survival.R
│
└── multimodal.R
```

Python adapters：

``` text
inst/pixi/
├── annotation/
├── trajectory/
├── grn/
├── communication/
├── cnv/
├── metabolism/
├── spatial/
└── multimodal/
```

方法注册：

``` text
inst/methods/
├── annotation.yml
├── trajectory.yml
├── program_scoring.yml
├── differential_abundance.yml
├── communication.yml
├── spatial.yml
└── bulk.yml
```

示例 YAML：

``` yaml
name: slingshot
task: trajectory
runtime: r
package: slingshot
input:
  object:
    - Seurat
    - SingleCellExperiment
requires:
  reduction: true
  cluster_by: true
outputs:
  - pseudotime
  - lineage_weights
  - curves
supports:
  branching: true
  large_data: true
```

------------------------------------------------------------------------

# 10. 第二条主线：Shennong Publication Figure Engine

## 10.1 当前基础

Shennong 当前已经具备：

- [`sn_plot_dim()`](https://songqi.org/shennong/dev/reference/sn_plot_dim.md)；
- [`sn_plot_feature()`](https://songqi.org/shennong/dev/reference/sn_plot_feature.md)；
- [`sn_plot_heatmap()`](https://songqi.org/shennong/dev/reference/sn_plot_heatmap.md)；
- [`sn_plot_dot()`](https://songqi.org/shennong/dev/reference/sn_plot_dot.md)；
- [`sn_plot_violin()`](https://songqi.org/shennong/dev/reference/sn_plot_violin.md)；
- [`sn_plot_boxplot()`](https://songqi.org/shennong/dev/reference/sn_plot_boxplot.md)；
- [`sn_plot_barplot()`](https://songqi.org/shennong/dev/reference/sn_plot_barplot.md)；
- [`sn_plot_composition()`](https://songqi.org/shennong/dev/reference/sn_plot_composition.md)；
- [`sn_plot_milo()`](https://songqi.org/shennong/dev/reference/sn_plot_milo.md)；
- named palette registry；
- 连续和离散调色板；
- 自动 UMAP 点大小；
- 大量点图层栅格化；
- patchwork 处理；
- compact legend；
- panel width / height；
- feature density；
- stored-result reuse。

这是很好的基础，但仍缺少：

1.  统一的图形尺寸计算器；
2.  统一的 export；
3.  统一的期刊/版面 profile；
4.  图形 QA；
5.  结果专用图形；
6.  多 panel 自动布局；
7.  label/legend overflow 预测；
8.  可视回归测试；
9.  figure manifest；
10. 一键导出主图、补图和源数据。

------------------------------------------------------------------------

## 10.2 设计目标

用户应当能够：

``` r

p <- sn_plot_dim(
  object,
  group_by = "cell_type"
)

sn_save_figure(
  p,
  "fig1a.pdf",
  profile = "single_column"
)
```

无需手工反复尝试：

``` r

width = 7
height = 5
pt.size = 0.1
legend.position = ...
```

图形系统必须根据：

- plot type；
- 数据点数量；
- category 数量；
- feature 数量；
- panel 数量；
- facet 数量；
- label 最大长度；
- legend item 数量；
- annotation track 数量；
- 输出 profile；

自动给出：

- width；
- height；
- point size；
- alpha；
- font size；
- line width；
- legend position；
- rasterization；
- DPI；
- panel layout；
- pagination。

------------------------------------------------------------------------

# 11. 不要立即引入复杂 Figure 类

第一阶段应保持所有绘图函数返回原生对象：

- `ggplot`；
- `patchwork`；
- `ComplexHeatmap`。

在对象上附加：

``` r

attr(p, "shennong_figure_spec")
```

例如：

``` r

list(
  plot_type = "embedding",
  data_summary = list(
    n_points = 5162897,
    n_groups = 18,
    n_panels = 1
  ),
  recommended = list(
    width_mm = 90,
    height_mm = 78,
    point_size = 0.04,
    alpha = 0.7,
    rasterize = TRUE,
    raster_dpi = 600
  ),
  profile = "single_column"
)
```

新增：

``` r

sn_figure_spec(p)
sn_recommend_figure_size(p, profile = NULL)
sn_apply_figure_profile(p, profile)
sn_validate_figure(p)
sn_save_figure(p, filename, profile = NULL)
sn_export_figure_bundle(p, path)
```

这样仍然支持：

``` r

p + ggplot2::theme(...)
```

------------------------------------------------------------------------

# 12. Figure profile

先提供通用 profile，不声称自动满足所有期刊实时要求。

``` yaml
screen:
  width_mm: 160
  height_mm: 100
  base_font_pt: 10
  raster_dpi: 150

single_column:
  width_mm: 85
  max_height_mm: 120
  base_font_pt: 7
  raster_dpi: 600

one_half_column:
  width_mm: 120
  max_height_mm: 160
  base_font_pt: 7
  raster_dpi: 600

double_column:
  width_mm: 180
  max_height_mm: 220
  base_font_pt: 7
  raster_dpi: 600

supplement_page:
  width_mm: 180
  max_height_mm: 240
  base_font_pt: 7
  raster_dpi: 450

slide_16_9:
  width_mm: 254
  height_mm: 143
  base_font_pt: 16
  raster_dpi: 180
```

文件位置：

``` text
inst/figure-profiles/default.yml
inst/figure-profiles/journals/
```

可以增加 journal preset，但必须在文档中提醒用户核对当前 author
guidelines。

------------------------------------------------------------------------

# 13. 自动尺寸规则

以下公式作为第一版启发式规则，随后通过 golden figures 校准。

## 13.1 UMAP / embedding

输入：

- `n_cells`；
- `n_panels`；
- `n_groups`；
- `max_label_chars`；
- 是否显示 cluster label；
- 是否 split。

建议：

``` r

point_size <- clamp(
  1.2 * (n_cells / 1000)^(-0.35),
  min = 0.03,
  max = 1.6
)

rasterize <- n_cells * n_panels > 50000
```

panel：

``` text
单 panel：55–75 mm 正方形
多 panel：每 panel 至少 45–55 mm
```

legend：

- `n_groups <= 8`：右侧；

- `9–18`：右侧多列或底部；

- `>18`：底部多列；

- 长 label 优先底部；

- > 30 类别时发出“颜色区分不足”警告。

## 13.2 Dot plot

``` text
width =
  left_label_margin
  + n_groups × 4–6 mm
  + legend_width

height =
  top_margin
  + n_features × 3.5–5 mm
  + facet_strip_height
```

规则：

- gene label 长度决定左边距；

- grouped feature list 使用 free-width facet；

- > 40 genes 自动建议分页；

- > 25 groups 自动减小 dot scale 并旋转标签；

- 保证最终字号不低于 profile 最小字号。

## 13.3 Heatmap

输入：

- rows；
- columns/groups；
- annotation tracks；
- label lengths；
- cell-level 或 average。

规则：

- cell-level 大矩阵默认 rasterize body；

- row/column labels、annotation 和 legend 保持 vector；

- > 100 rows 自动分页或只标记 selected labels；

- average heatmap 按 group 数计算宽度；

- > 50 groups 时优先横向分页；

- output PDF 中禁止数百万个 vector rect。

## 13.4 Violin / box / bar

``` text
width =
  max(
    profile minimum,
    n_categories × category_width
    + label_margin
    + comparison_margin
  )
```

自动行为：

- category \<= 6：垂直；
- category 7–15：扩大宽度或旋转 30–45°；
- category \>15：建议 horizontal；
- 显著性 bracket 增加顶部空间；
- paired data 自动预留连接线；
- sample-level points 优先显示，cell-level points 仅作背景。

## 13.5 Composition

根据 sample 数和 cell-type 数决定：

- bar width；
- legend 行列；
- palette warning；
- 是否显示 label；
- 是否转为 heatmap；
- sample \>30 时建议横向或 heatmap。

## 13.6 Network / CCI

根据 node/edge 数：

- \<=15 nodes：完整 label；

- 16–40：缩小 label、过滤低权重 edge；

- > 40：默认聚合到 pathway 或 cell type；

- chord/network 不允许在默认状态展示数千条 edge；

- 自动返回过滤摘要。

## 13.7 Spatial

根据：

- spot/cell 数；
- tissue aspect ratio；
- image pixel dimensions；
- facet 数；
- scale bar；
- crop。

空间图必须保持真实坐标比例，不使用固定正方形强制拉伸。

------------------------------------------------------------------------

# 14. 必须新增的图形类型

## 14.1 Core/QC

``` r

sn_plot_qc()
sn_plot_qc_thresholds()
sn_plot_doublets()
sn_plot_ambient_correction()
sn_plot_hvg()
sn_plot_elbow()
sn_plot_cluster_tree()
sn_plot_resolution_sweep()
sn_plot_integration()
```

## 14.2 Annotation

``` r

sn_plot_annotation_confidence()
sn_plot_annotation_markers()
sn_plot_annotation_confusion()
sn_plot_reference_projection()
```

## 14.3 DE 和 abundance

``` r

sn_plot_de(
  result,
  type = c("volcano", "ma", "effect", "heatmap")
)

sn_plot_abundance(
  result,
  type = c("composition", "effect", "neighborhood", "beeswarm")
)
```

## 14.4 Enrichment 和 programs

``` r

sn_plot_enrichment(
  result,
  type = c("dot", "bar", "ridge", "network", "emap")
)

sn_plot_gsea()
sn_plot_program_activity()
sn_plot_program_heatmap()
sn_plot_regulon()
```

## 14.5 Trajectory

``` r

sn_plot_trajectory()
sn_plot_pseudotime()
sn_plot_dynamic_heatmap()
sn_plot_gene_trend()
sn_plot_velocity()
sn_plot_fate()
```

## 14.6 Communication

``` r

sn_plot_communication(
  result,
  type = c("bubble", "heatmap", "network", "chord", "river")
)

sn_plot_ligand_target()
sn_plot_communication_comparison()
```

## 14.7 Spatial

``` r

sn_plot_spatial()
sn_plot_spatial_feature()
sn_plot_spatial_domain()
sn_plot_spatial_svg()
sn_plot_spatial_neighborhood()
sn_plot_spatial_deconvolution()
sn_plot_spatial_communication()
```

## 14.8 Bulk

``` r

sn_plot_bulk_qc()
sn_plot_bulk_pca()
sn_plot_sample_correlation()
sn_plot_bulk_de()
sn_plot_wgcna()
sn_plot_survival()
```

------------------------------------------------------------------------

# 15. 统一图形语法

Shennong 应保持显式、可发现的函数，同时减少重复参数。

所有 `sn_plot_*()` 尽量共享：

``` r

profile = NULL
width = "auto"
height = "auto"
font_size = "auto"
rasterize = "auto"
raster_dpi = "auto"
palette = "auto"
legend = "auto"
labels = "auto"
```

`"auto"` 必须是真正由 figure engine 解析，而不是每个函数各写一套逻辑。

内部统一：

``` r

spec <- .sn_build_figure_spec(
  plot_type,
  data_summary,
  profile,
  overrides
)

p <- .sn_render_...
p <- .sn_apply_figure_spec(p, spec)

attr(p, "shennong_figure_spec") <- spec
```

------------------------------------------------------------------------

# 16. 导出系统

## 16.1 `sn_save_figure()`

``` r

sn_save_figure(
  plot,
  filename,
  profile = NULL,
  width = "auto",
  height = "auto",
  units = "mm",
  dpi = "auto",
  background = "white",
  embed_fonts = TRUE,
  validate = TRUE
)
```

支持：

- PDF；
- SVG；
- TIFF；
- PNG。

默认：

- PDF/SVG：文字和几何保持 vector；
- 高密度点/heatmap body：局部 raster；
- TIFF：600 dpi 用于线稿/混合图；
- PNG：屏幕和文档预览；
- 背景显式；
- 不依赖当前 RStudio device 尺寸。

## 16.2 Figure bundle

``` r

sn_export_figure_bundle(
  plot,
  path = "figures/Figure_2A",
  formats = c("pdf", "png"),
  include_data = TRUE,
  include_spec = TRUE,
  include_session = TRUE
)
```

输出：

``` text
Figure_2A/
├── Figure_2A.pdf
├── Figure_2A.png
├── Figure_2A_data.csv
├── Figure_2A_spec.yml
├── Figure_2A_session.txt
└── Figure_2A_manifest.json
```

------------------------------------------------------------------------

# 17. Figure QA

新增：

``` r

sn_validate_figure(plot, profile = NULL)
```

检查：

- 最终尺寸；
- 最小字号；
- legend overflow；
- label clipping；
- panel 太小；
- 类别过多；
- 调色板重复；
- 色觉缺陷风险；
- raster DPI；
- PDF 体积；
- 透明背景；
- 非嵌入字体；
- 极端 aspect ratio；
- 空 panel；
- 缺失值被静默删除；
- heatmap label 不可读；
- 网络 edge 过多。

返回：

``` text
Figure validation: WARNING

✓ Width: 85 mm
✓ Raster layer: 600 dpi
✓ Minimum text size: 7 pt
! 24 legend categories may be difficult to distinguish
! Longest label is 31 characters
Suggested action:
  use legend = "bottom", ncol = 4
```

------------------------------------------------------------------------

# 18. 可视化测试策略

## 18.1 单元测试

测试自动尺寸规则：

``` r

expect_gt(
  sn_recommend_figure_size(dot_40_genes)$height,
  sn_recommend_figure_size(dot_10_genes)$height
)
```

测试：

- panel 数增加时尺寸不减小；
- label 更长时 margin 不减小；
- 点数超过阈值时启用 raster；
- 用户覆盖值优先于 auto；
- profile 尺寸限制有效；
- unsupported format 清晰报错。

## 18.2 可视回归

使用：

- `vdiffr`；
- rendered gtable 检查；
- golden PBMC figures；
- spatial aspect-ratio figures；
- large synthetic object。

必须覆盖：

- 500 cells；
- 5,000 cells；
- 50,000 cells；
- 500,000 cells；
- 5,000,000 cells 的 spec calculation。

最大数据测试不必真实渲染所有点，可以对 spec calculator 使用模拟
metadata。

## 18.3 导出测试

验证：

- PDF 页面尺寸；
- TIFF/PNG 像素尺寸；
- DPI；
- 文件可打开；
- alpha/background；
- multipage；
- figure bundle manifest；
- source data 完整。

------------------------------------------------------------------------

# 19. 推荐可视化文件结构

``` text
R/
├── figure-spec.R
├── figure-profile.R
├── figure-layout.R
├── figure-style.R
├── figure-export.R
├── figure-validate.R
├── figure-bundle.R
│
├── plot-core.R
├── plot-qc.R
├── plot-annotation.R
├── plot-de.R
├── plot-enrichment.R
├── plot-program.R
├── plot-trajectory.R
├── plot-communication.R
├── plot-spatial.R
├── plot-bulk.R
└── palettes.R

inst/figure-profiles/
├── default.yml
└── journals/

tests/testthat/
├── test-figure-spec.R
├── test-figure-profile.R
├── test-figure-export.R
├── test-figure-validation.R
└── test-plot-*.R
```

当前过大的 `R/visualization.R`
应按功能逐步拆分，但每次拆分必须保持行为不变并有回归测试。

------------------------------------------------------------------------

# 20. 里程碑与 PR 拆分

不要在一个 PR 中实现全部内容。

## Milestone A：分析契约与 registry

### PR A1

``` text
feat(analysis): add unified method registry
```

- [`sn_list_methods()`](https://songqi.org/shennong/dev/reference/sn_list_methods.md)；
- [`sn_method_status()`](https://songqi.org/shennong/dev/reference/sn_method_status.md)；
- YAML registry；
- optional dependency diagnostics。

### PR A2

``` text
feat(results): standardize analysis result schema
```

- result validator；
- generic store/get/list；
- provenance；
- 迁移一个现有 DE result 作为示例。

------------------------------------------------------------------------

## Milestone B：单细胞主线补齐

### PR B1

``` text
feat(annotation): add consensus annotation workflow
```

- SingleR；
- confidence；
- marker evidence；
- ontology；
- plot。

### PR B2

``` text
feat(programs): add unified gene-set scoring
```

- UCell；
- AUCell；
- GSVA；
- result and plots。

### PR B3

``` text
feat(trajectory): add Slingshot and dynamic gene workflow
```

- Slingshot；
- tradeSeq；
- trajectory plots。

### PR B4

``` text
feat(abundance): add propeller and unified DA interface
```

- propeller；
- Milo adapter；
- result plots。

### PR B5

``` text
refactor(communication): standardize communication results
```

- backend harmonization；
- consensus；
- comparison plots。

### PR B6

``` text
feat(grn): add optional SCENIC workflow
```

### PR B7

``` text
feat(dynamics): add scVelo and CellRank pixi backend
```

------------------------------------------------------------------------

## Milestone C：空间与 bulk

### PR C1

``` text
feat(spatial): add spatial feature and domain workflows
```

### PR C2

``` text
feat(spatial): add neighborhood and communication workflows
```

### PR C3

``` text
feat(bulk): add bulk QC and differential expression
```

### PR C4

``` text
feat(bulk): add WGCNA and clinical association
```

------------------------------------------------------------------------

## Milestone D：Figure Engine

### PR D1

``` text
feat(figures): add figure profiles and auto-size engine
```

### PR D2

``` text
feat(figures): add publication export and validation
```

### PR D3

``` text
refactor(plots): migrate existing plots to figure specs
```

迁移顺序：

1.  dim；
2.  feature；
3.  dot；
4.  heatmap；
5.  violin/box/bar；
6.  composition；
7.  Milo。

### PR D4

``` text
feat(plots): add result-aware DE and enrichment figures
```

### PR D5

``` text
feat(plots): add trajectory and communication figures
```

### PR D6

``` text
test(figures): add visual regression suite
```

------------------------------------------------------------------------

# 21. 每个新分析方法的完成定义

一个方法只有同时满足以下条件才算完成：

有清晰的顶层 Shennong API；

有 input validation；

有合理默认 method；

optional dependency 不会破坏 package load；

返回或存储统一 result schema；

保存参数、版本、seed 和 warnings；

有 getter；

有至少一个诊断图；

有至少一个 publication plot；

有 testthat 测试；

有最小示例；

有完整 vignette；

`R CMD check` 通过；

Linux/macOS/Windows 中可运行的 R 方法有 CI；

Python backend 有 CPU smoke test；

大数据路径不会无意 densify sparse matrix；

API 参数可以转换为 JSON。

------------------------------------------------------------------------

# 22. Figure Engine 的完成定义

所有现有 `sn_plot_*()` 返回原生可编辑对象；

每张图具有 `shennong_figure_spec`；

[`sn_save_figure()`](https://songqi.org/shennong/dev/reference/sn_save_figure.md)
不依赖当前 device；

支持 mm/in/cm；

PDF/SVG vector 元素保持 vector；

大点层和 heatmap 可局部 raster；

自动尺寸随数据复杂度单调变化；

用户覆盖参数优先；

有 single/double column profile；

有 figure validation；

有 multipage；

有 source-data export；

有 visual regression；

所有 vignette 图由同一导出系统生成；

图例、字体和 panel 在最终物理尺寸下可读。

------------------------------------------------------------------------

# 23. Shennong 与三个对照项目的最终差异

## 相对于 sclet

Shennong 不必复制它的全部 `Run*()` 风格。

应吸收：

- state-aware results；
- workflow mainline；
- provenance；
- AI evidence summary。

Shennong 的差异：

- 更强的数据服务连接；
- 更强的 Seurat 大数据实战；
- pixi Python runtime；
- 更完整的 bulk/spatial 计划；
- publication figure engine。

## 相对于 scop

Shennong 不应比拼函数数量。

应吸收：

- annotation/trajectory/spatial 的核心方法；
- 结果专用图形；
- 多种 optional backend。

Shennong 的差异：

- 更少但更有主次的方法；
- 自动 backend recommendation；
- 统一结果 schema；
- 明确 diagnostic；
- 发表尺寸感知；
- figure QA；
- 数据和结果可回写 Shennong 生态。

## 相对于 OmicVerse

Shennong 不应在当前阶段复制所有 Python 深度学习模块。

应吸收：

- 模块化；
- function registry；
- lazy optional dependencies；
- bulk/single/spatial 分区；
- analysis report。

Shennong 的差异：

- R/Bioconductor/Seurat 用户优先；
- 统计 design 优先；
- 结果可审计；
- 默认方法透明；
- publication-ready output 为核心产品，而不是附属 plotting module。

------------------------------------------------------------------------

# 24. 最终开发顺序

严格按照以下顺序：

``` text
1. 分析 method registry 和 result contract
2. annotation consensus
3. per-cell program scoring
4. trajectory + dynamic genes
5. differential abundance
6. communication result standardization
7. spatial SVG/domain/neighborhood
8. bulk QC/DE/pathway
9. Figure Engine core
10. 迁移现有图形
11. 新分析结果图形
12. velocity/fate/GRN 等重型后端
13. multimodal
```

Figure Engine 可以与分析方法并行开发，但必须先稳定：

``` text
figure profile
auto size
export
validation
```

再大量添加新的 plot types。

------------------------------------------------------------------------

# 25. 对 Codex 的直接执行要求

1.  不一次性增加大量依赖。
2.  不把 Python 包加入 R 的强制依赖。
3.  不使用 `eval(parse())`。
4.  不在顶层 API 暴露任意 Python 函数名。
5.  不无意将 sparse matrix 转成 dense matrix。
6.  不在每个函数重复实现 result storage。
7.  不在每个 plot 函数重复实现尺寸逻辑。
8.  不删除现有测试以让新实现通过。
9.  不降低现有 Seurat 兼容性。
10. 每个 PR 使用 conventional commit。
11. 每个新增功能必须更新：
    - `NAMESPACE`；
    - roxygen；
    - `_pkgdown.yml`；
    - vignette；
    - tests；
    - `NEWS.md`；
    - `docs/codex/Decisions.md`。
12. 每个 backend 必须说明：
    - 为什么选它；
    - 什么情况下不应使用；
    - 输入要求；
    - 输出结构；
    - CPU/GPU；
    - 依赖安装；
    - 引用。

------------------------------------------------------------------------

# 26. 参考来源

审查日期：2026-07-15。

## Shennong

- Repository: <https://github.com/zerostwo/shennong>
- Current exported API:
  <https://github.com/zerostwo/shennong/blob/main/NAMESPACE>
- Package dependencies:
  <https://github.com/zerostwo/shennong/blob/main/DESCRIPTION>
- Visualization vignette:
  <https://github.com/zerostwo/shennong/blob/main/vignettes/visualization.Rmd>
- Modernization decisions:
  <https://github.com/zerostwo/shennong/blob/main/docs/codex/Decisions.md>

## sclet

- Repository: <https://github.com/YuLab-SMU/sclet>
- README/mainlines:
  <https://github.com/YuLab-SMU/sclet/blob/devel/README.md>
- Exported API:
  <https://github.com/YuLab-SMU/sclet/blob/devel/NAMESPACE>
- Visualization implementation:
  <https://github.com/YuLab-SMU/sclet/blob/devel/R/visualization.R>

## scop

- Repository: <https://github.com/mengxu98/scop>
- README/capability overview:
  <https://github.com/mengxu98/scop/blob/main/README.md>
- Exported API: <https://github.com/mengxu98/scop/blob/main/NAMESPACE>

## OmicVerse

- Repository: <https://github.com/omicverse/omicverse>
- README/platform overview:
  <https://github.com/omicverse/omicverse/blob/master/README.md>
- Top-level modules:
  <https://github.com/omicverse/omicverse/blob/master/omicverse/__init__.py>
- Single-cell module:
  <https://github.com/omicverse/omicverse/blob/master/omicverse/single/__init__.py>
- Bulk module:
  <https://github.com/omicverse/omicverse/blob/master/omicverse/bulk/__init__.py>
- Spatial module:
  <https://github.com/omicverse/omicverse/blob/master/omicverse/space/__init__.py>

## tidyplots

- Website: <https://tidyplots.org/>
- Documentation: <https://jbengler.github.io/tidyplots/>
- Advanced plotting, rasterization and multipage layouts:
  <https://jbengler.github.io/tidyplots/articles/Advanced-plotting.html>

------------------------------------------------------------------------

# 27. 一句话目标

> Shennong
> 应成为一个能够完成真实组学分析闭环，并默认生成可直接进入论文初稿图形的
> R 生物信息学框架。
