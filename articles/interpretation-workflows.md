# Interpretation workflow

This article shows how the new interpretation layer fits on top of the
current Shennong analysis workflow. The key idea is that analysis stays
deterministic, while interpretation reuses stored DE and enrichment
results to build evidence bundles, prompts, and optional LLM-style
summaries.

The workflow below demonstrates:

- marker discovery with
  [`sn_find_de()`](https://songqi.org/shennong/reference/sn_find_de.md)
- enrichment storage with
  [`sn_store_enrichment()`](https://songqi.org/shennong/reference/sn_store_enrichment.md)
- evidence preparation for annotation, DE, and writing tasks
- prompt generation with
  [`sn_build_prompt()`](https://songqi.org/shennong/reference/sn_build_prompt.md)
- provider-based interpretation with
  [`sn_interpret_annotation()`](https://songqi.org/shennong/reference/sn_interpret_annotation.md)

``` r
library(Shennong)
library(dplyr)
library(knitr)
library(Seurat)
```

## Marker discovery remains the analysis layer

Interpretation starts from stored DE results rather than recomputing
statistics. Here,
[`sn_find_de()`](https://songqi.org/shennong/reference/sn_find_de.md)
has already saved cluster markers into `object@misc$de_results`.

``` r
knitr::kable(
  pbmc@misc$de_results$cluster_markers$table |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(order_by = avg_log2FC, n = 3, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(cluster, gene, avg_log2FC, p_val_adj)
)
```

| cluster | gene            | avg_log2FC | p_val_adj |
|:--------|:----------------|-----------:|----------:|
| 0       | BPI             |   5.812083 |         0 |
| 0       | ENSG00000289381 |   5.667320 |         0 |
| 0       | MTARC1          |   5.624992 |         0 |
| 1       | IATPR           |   4.782931 |         0 |
| 1       | TTC39C-AS1      |   4.436556 |         0 |
| 1       | PI16            |   4.155404 |         0 |
| 2       | EDAR            |   4.702752 |         0 |
| 2       | ANKRD55         |   4.373791 |         0 |
| 2       | ADTRP           |   4.130501 |         0 |
| 3       | TCL1A           |  10.568722 |         0 |
| 3       | ENSG00000257275 |  10.050273 |         0 |
| 3       | SCN3A           |   7.897144 |         0 |
| 4       | LINC02446       |   7.370057 |         0 |
| 4       | CD8B            |   5.570102 |         0 |
| 4       | ENSG00000310107 |   5.000580 |         0 |
| 5       | CLEC10A         |   6.518629 |         0 |
| 5       | FCER1A          |   4.873550 |         0 |
| 5       | CYP2S1          |   4.847983 |         0 |
| 6       | SLC4A10         |   9.112196 |         0 |
| 6       | ENSG00000228033 |   7.882580 |         0 |
| 6       | IL23R           |   7.568894 |         0 |
| 7       | IGHA1           |  10.217845 |         0 |
| 7       | IGHG2           |   9.501537 |         0 |
| 7       | IGHG3           |   8.979619 |         0 |
| 8       | ENSG00000288970 |   9.465656 |         0 |
| 8       | ENSG00000291157 |   8.923955 |         0 |
| 8       | MYOM2           |   8.903953 |         0 |
| 9       | MT-ND6          |   3.780937 |         1 |
| 9       | MT-CYB          |   3.628350 |         0 |
| 9       | MT-CO3          |   3.523565 |         0 |
| 10      | LYPD2           |  11.502831 |         0 |
| 10      | ENSG00000301038 |   9.496434 |         0 |
| 10      | UICLM           |   9.381918 |         0 |
| 11      | CLDN5           |  13.204275 |         0 |
| 11      | PF4V1           |  13.066932 |         0 |
| 11      | CMTM5           |  13.059807 |         0 |

## Store enrichment for later interpretation

The interpretation layer expects enrichment results to be stored in
`object@misc$enrichment_results`, just like DE results already live in
`object@misc$de_results`.

``` r
knitr::kable(pbmc@misc$enrichment_results$cluster0_gsea$table)
```

| ID           | Description                                                                                                                                      |        NES |  p.adjust |
|:-------------|:-------------------------------------------------------------------------------------------------------------------------------------------------|-----------:|----------:|
| <GO:0009617> | response to bacterium                                                                                                                            |  1.7311779 | 0.0000002 |
| <GO:0006954> | inflammatory response                                                                                                                            |  1.6105454 | 0.0000002 |
| <GO:0006955> | immune response                                                                                                                                  |  1.4839565 | 0.0000002 |
| <GO:0002682> | regulation of immune system process                                                                                                              |  1.4732294 | 0.0000007 |
| <GO:0006952> | defense response                                                                                                                                 |  1.4607568 | 0.0000008 |
| <GO:0051240> | positive regulation of multicellular organismal process                                                                                          |  1.4871140 | 0.0000019 |
| <GO:0042742> | defense response to bacterium                                                                                                                    |  1.8313202 | 0.0000135 |
| <GO:0002684> | positive regulation of immune system process                                                                                                     |  1.4805177 | 0.0000135 |
| <GO:0032103> | positive regulation of response to external stimulus                                                                                             |  1.5547745 | 0.0000151 |
| <GO:0032101> | regulation of response to external stimulus                                                                                                      |  1.4716948 | 0.0000176 |
| <GO:0050729> | positive regulation of inflammatory response                                                                                                     |  1.8997862 | 0.0000198 |
| <GO:0031349> | positive regulation of defense response                                                                                                          |  1.5900966 | 0.0000209 |
| <GO:0050727> | regulation of inflammatory response                                                                                                              |  1.7205855 | 0.0000221 |
| <GO:0051707> | response to other organism                                                                                                                       |  1.4197927 | 0.0000244 |
| <GO:0009607> | response to biotic stimulus                                                                                                                      |  1.4181958 | 0.0000244 |
| <GO:0002274> | myeloid leukocyte activation                                                                                                                     |  1.7638875 | 0.0000320 |
| <GO:0050776> | regulation of immune response                                                                                                                    |  1.4752673 | 0.0000320 |
| <GO:0043207> | response to external biotic stimulus                                                                                                             |  1.4276942 | 0.0000508 |
| <GO:0071216> | cellular response to biotic stimulus                                                                                                             |  1.7425289 | 0.0000603 |
| <GO:0001817> | regulation of cytokine production                                                                                                                |  1.5048487 | 0.0000603 |
| <GO:2000026> | regulation of multicellular organismal development                                                                                               |  1.4717772 | 0.0000638 |
| <GO:0001816> | cytokine production                                                                                                                              |  1.5000333 | 0.0001101 |
| <GO:0051241> | negative regulation of multicellular organismal process                                                                                          |  1.4694270 | 0.0001368 |
| <GO:0006909> | phagocytosis                                                                                                                                     |  1.7263507 | 0.0001462 |
| <GO:0002237> | response to molecule of bacterial origin                                                                                                         |  1.6563569 | 0.0001576 |
| <GO:0051094> | positive regulation of developmental process                                                                                                     |  1.4484662 | 0.0001576 |
| <GO:0071222> | cellular response to lipopolysaccharide                                                                                                          |  1.7687924 | 0.0001611 |
| <GO:0098542> | defense response to other organism                                                                                                               |  1.4217183 | 0.0001611 |
| <GO:0071219> | cellular response to molecule of bacterial origin                                                                                                |  1.7474425 | 0.0001982 |
| <GO:0031347> | regulation of defense response                                                                                                                   |  1.4653261 | 0.0001982 |
| <GO:0040011> | locomotion                                                                                                                                       |  1.4446355 | 0.0001982 |
| <GO:0140546> | defense response to symbiont                                                                                                                     |  1.4396349 | 0.0002119 |
| <GO:0042776> | proton motive force-driven mitochondrial ATP synthesis                                                                                           | -2.5150280 | 0.0002216 |
| <GO:0002764> | immune response-regulating signaling pathway                                                                                                     |  1.5056810 | 0.0002690 |
| <GO:0098657> | import into cell                                                                                                                                 |  1.4543130 | 0.0002901 |
| <GO:0032611> | interleukin-1 beta production                                                                                                                    |  1.8269296 | 0.0004308 |
| <GO:0032651> | regulation of interleukin-1 beta production                                                                                                      |  1.8269296 | 0.0004308 |
| <GO:0033993> | response to lipid                                                                                                                                |  1.4640098 | 0.0004308 |
| <GO:0050793> | regulation of developmental process                                                                                                              |  1.3273177 | 0.0004526 |
| <GO:0050778> | positive regulation of immune response                                                                                                           |  1.4399659 | 0.0005009 |
| <GO:0045087> | innate immune response                                                                                                                           |  1.4051820 | 0.0005246 |
| <GO:0006935> | chemotaxis                                                                                                                                       |  1.6087118 | 0.0005429 |
| <GO:0042330> | taxis                                                                                                                                            |  1.6087118 | 0.0005429 |
| <GO:0071396> | cellular response to lipid                                                                                                                       |  1.5333662 | 0.0006392 |
| <GO:0044419> | biological process involved in interspecies interaction between organisms                                                                        |  1.3442300 | 0.0006593 |
| <GO:0048513> | animal organ development                                                                                                                         |  1.3339201 | 0.0006593 |
| <GO:0032496> | response to lipopolysaccharide                                                                                                                   |  1.6350102 | 0.0007314 |
| <GO:0006959> | humoral immune response                                                                                                                          |  1.8048292 | 0.0007960 |
| <GO:0006897> | endocytosis                                                                                                                                      |  1.4376989 | 0.0008369 |
| <GO:0001818> | negative regulation of cytokine production                                                                                                       |  1.6234124 | 0.0015579 |
| <GO:0009888> | tissue development                                                                                                                               |  1.3568760 | 0.0015579 |
| <GO:0098581> | detection of external biotic stimulus                                                                                                            |  1.9215674 | 0.0015675 |
| <GO:0036230> | granulocyte activation                                                                                                                           |  1.8970159 | 0.0019143 |
| <GO:0032731> | positive regulation of interleukin-1 beta production                                                                                             |  1.8240141 | 0.0019506 |
| <GO:0032612> | interleukin-1 production                                                                                                                         |  1.7696306 | 0.0023228 |
| <GO:0032652> | regulation of interleukin-1 production                                                                                                           |  1.7696306 | 0.0023228 |
| <GO:0007186> | G protein-coupled receptor signaling pathway                                                                                                     |  1.5102554 | 0.0025663 |
| <GO:0002768> | immune response-regulating cell surface receptor signaling pathway                                                                               |  1.5430602 | 0.0025989 |
| <GO:0009595> | detection of biotic stimulus                                                                                                                     |  1.9127146 | 0.0031257 |
| <GO:0045597> | positive regulation of cell differentiation                                                                                                      |  1.4305462 | 0.0032585 |
| <GO:0097529> | myeloid leukocyte migration                                                                                                                      |  1.6705455 | 0.0035965 |
| <GO:0043065> | positive regulation of apoptotic process                                                                                                         |  1.4932671 | 0.0038755 |
| <GO:0050764> | regulation of phagocytosis                                                                                                                       |  1.7280451 | 0.0041623 |
| <GO:0060326> | cell chemotaxis                                                                                                                                  |  1.5846432 | 0.0044073 |
| <GO:0001819> | positive regulation of cytokine production                                                                                                       |  1.4569330 | 0.0046150 |
| <GO:0002833> | positive regulation of response to biotic stimulus                                                                                               |  1.4473689 | 0.0046812 |
| <GO:0043068> | positive regulation of programmed cell death                                                                                                     |  1.4553635 | 0.0047009 |
| <GO:0007155> | cell adhesion                                                                                                                                    |  1.3455047 | 0.0050729 |
| <GO:0045595> | regulation of cell differentiation                                                                                                               |  1.3391475 | 0.0050863 |
| <GO:0019730> | antimicrobial humoral response                                                                                                                   |  1.8799339 | 0.0052951 |
| <GO:0050766> | positive regulation of phagocytosis                                                                                                              |  1.7726523 | 0.0054774 |
| <GO:0042119> | neutrophil activation                                                                                                                            |  1.8463663 | 0.0058850 |
| <GO:0045089> | positive regulation of innate immune response                                                                                                    |  1.4398017 | 0.0059581 |
| <GO:0071621> | granulocyte chemotaxis                                                                                                                           |  1.7518978 | 0.0059826 |
| <GO:0045807> | positive regulation of endocytosis                                                                                                               |  1.6735320 | 0.0066123 |
| <GO:0007200> | phospholipase C-activating G protein-coupled receptor signaling pathway                                                                          |  1.8018553 | 0.0071723 |
| <GO:0032732> | positive regulation of interleukin-1 production                                                                                                  |  1.7449046 | 0.0075291 |
| <GO:0002252> | immune effector process                                                                                                                          |  1.3903214 | 0.0075291 |
| <GO:0030595> | leukocyte chemotaxis                                                                                                                             |  1.5978051 | 0.0089051 |
| <GO:0030100> | regulation of endocytosis                                                                                                                        |  1.5094922 | 0.0089051 |
| <GO:0050865> | regulation of cell activation                                                                                                                    |  1.3980308 | 0.0089051 |
| <GO:0002253> | activation of immune response                                                                                                                    |  1.3884949 | 0.0089051 |
| <GO:0032637> | interleukin-8 production                                                                                                                         |  1.7161460 | 0.0091516 |
| <GO:0032677> | regulation of interleukin-8 production                                                                                                           |  1.7161460 | 0.0091516 |
| <GO:0050830> | defense response to Gram-positive bacterium                                                                                                      |  1.7250597 | 0.0093551 |
| <GO:0002757> | immune response-activating signaling pathway                                                                                                     |  1.3994554 | 0.0093551 |
| <GO:0045321> | leukocyte activation                                                                                                                             |  1.3394480 | 0.0093551 |
| <GO:0002694> | regulation of leukocyte activation                                                                                                               |  1.4127773 | 0.0114410 |
| <GO:1903036> | positive regulation of response to wounding                                                                                                      |  1.7948439 | 0.0115761 |
| <GO:0042127> | regulation of cell population proliferation                                                                                                      |  1.3154341 | 0.0115761 |
| <GO:0030888> | regulation of B cell proliferation                                                                                                               |  1.7921079 | 0.0117822 |
| <GO:0002831> | regulation of response to biotic stimulus                                                                                                        |  1.3772310 | 0.0127784 |
| <GO:0030593> | neutrophil chemotaxis                                                                                                                            |  1.7144721 | 0.0130954 |
| <GO:0040012> | regulation of locomotion                                                                                                                         |  1.3506459 | 0.0133928 |
| <GO:0001775> | cell activation                                                                                                                                  |  1.3167525 | 0.0133928 |
| <GO:0002824> | positive regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains |  1.6926417 | 0.0134336 |
| <GO:1901701> | cellular response to oxygen-containing compound                                                                                                  |  1.3103909 | 0.0134336 |
| <GO:2000145> | regulation of cell motility                                                                                                                      |  1.3413096 | 0.0144665 |
| <GO:0010574> | regulation of vascular endothelial growth factor production                                                                                      |  1.7751085 | 0.0146559 |
| <GO:0002548> | monocyte chemotaxis                                                                                                                              |  1.7575009 | 0.0146559 |
| <GO:0042116> | macrophage activation                                                                                                                            |  1.6499793 | 0.0146559 |
| <GO:0002224> | toll-like receptor signaling pathway                                                                                                             |  1.6468745 | 0.0146559 |
| <GO:0097530> | granulocyte migration                                                                                                                            |  1.6300622 | 0.0155426 |
| <GO:0060284> | regulation of cell development                                                                                                                   |  1.3730263 | 0.0155426 |
| <GO:0043408> | regulation of MAPK cascade                                                                                                                       |  1.3946913 | 0.0169641 |
| <GO:0045088> | regulation of innate immune response                                                                                                             |  1.3621389 | 0.0174899 |
| <GO:1901700> | response to oxygen-containing compound                                                                                                           |  1.2825818 | 0.0193648 |
| <GO:0035872> | nucleotide-binding domain, leucine rich repeat containing receptor signaling pathway                                                             |  1.7907397 | 0.0204251 |
| <GO:0031214> | biomineral tissue development                                                                                                                    |  1.6225668 | 0.0204251 |
| <GO:0030890> | positive regulation of B cell proliferation                                                                                                      |  1.7682800 | 0.0208995 |
| <GO:0045730> | respiratory burst                                                                                                                                |  1.7239171 | 0.0208995 |
| <GO:0032640> | tumor necrosis factor production                                                                                                                 |  1.5498738 | 0.0208995 |
| <GO:0032680> | regulation of tumor necrosis factor production                                                                                                   |  1.5498738 | 0.0208995 |
| <GO:0050900> | leukocyte migration                                                                                                                              |  1.4599614 | 0.0208995 |
| <GO:0010628> | positive regulation of gene expression                                                                                                           |  1.2805391 | 0.0208995 |
| <GO:0014002> | astrocyte development                                                                                                                            |  1.7401563 | 0.0210238 |
| <GO:1902105> | regulation of leukocyte differentiation                                                                                                          |  1.4794830 | 0.0210238 |
| <GO:0030282> | bone mineralization                                                                                                                              |  1.6596120 | 0.0213866 |
| <GO:0071706> | tumor necrosis factor superfamily cytokine production                                                                                            |  1.5291733 | 0.0217159 |
| <GO:1903555> | regulation of tumor necrosis factor superfamily cytokine production                                                                              |  1.5291733 | 0.0217159 |
| <GO:0043030> | regulation of macrophage activation                                                                                                              |  1.7372249 | 0.0218318 |
| <GO:0019884> | antigen processing and presentation of exogenous antigen                                                                                         |  1.7793077 | 0.0227945 |
| <GO:0003018> | vascular process in circulatory system                                                                                                           |  1.5921963 | 0.0228386 |
| <GO:0006809> | nitric oxide biosynthetic process                                                                                                                |  1.6596647 | 0.0229099 |
| <GO:0008283> | cell population proliferation                                                                                                                    |  1.2657385 | 0.0229099 |
| <GO:1990266> | neutrophil migration                                                                                                                             |  1.6575644 | 0.0246459 |
| <GO:1903706> | regulation of hemopoiesis                                                                                                                        |  1.4306509 | 0.0247377 |
| <GO:0015849> | organic acid transport                                                                                                                           |  1.5492122 | 0.0267895 |
| <GO:0046942> | carboxylic acid transport                                                                                                                        |  1.5492122 | 0.0267895 |
| <GO:0070887> | cellular response to chemical stimulus                                                                                                           |  1.2387703 | 0.0267895 |
| <GO:0048646> | anatomical structure formation involved in morphogenesis                                                                                         |  1.3247154 | 0.0269568 |
| <GO:0002250> | adaptive immune response                                                                                                                         |  1.4242513 | 0.0271216 |
| <GO:0030316> | osteoclast differentiation                                                                                                                       |  1.5761878 | 0.0276777 |
| <GO:0002683> | negative regulation of immune system process                                                                                                     |  1.4132612 | 0.0286836 |
| <GO:0035295> | tube development                                                                                                                                 |  1.3447217 | 0.0286836 |
| <GO:0061844> | antimicrobial humoral immune response mediated by antimicrobial peptide                                                                          |  1.7540889 | 0.0303851 |
| <GO:0050867> | positive regulation of cell activation                                                                                                           |  1.4133418 | 0.0314647 |
| <GO:0031663> | lipopolysaccharide-mediated signaling pathway                                                                                                    |  1.6457332 | 0.0324136 |
| <GO:0030855> | epithelial cell differentiation                                                                                                                  |  1.3993928 | 0.0324136 |
| <GO:0002460> | adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains                        |  1.4858017 | 0.0325300 |
| <GO:0034308> | primary alcohol metabolic process                                                                                                                |  1.7036020 | 0.0327578 |
| <GO:0035239> | tube morphogenesis                                                                                                                               |  1.3525886 | 0.0327578 |
| <GO:1903707> | negative regulation of hemopoiesis                                                                                                               |  1.6238173 | 0.0330676 |
| <GO:0030334> | regulation of cell migration                                                                                                                     |  1.3124991 | 0.0336558 |
| <GO:0032635> | interleukin-6 production                                                                                                                         |  1.5170580 | 0.0356295 |
| <GO:0032675> | regulation of interleukin-6 production                                                                                                           |  1.5170580 | 0.0356295 |
| <GO:0002218> | activation of innate immune response                                                                                                             |  1.3658585 | 0.0356295 |
| <GO:0050832> | defense response to fungus                                                                                                                       |  1.6722821 | 0.0356900 |
| <GO:0007204> | positive regulation of cytosolic calcium ion concentration                                                                                       |  1.6661843 | 0.0375275 |
| <GO:0080134> | regulation of response to stress                                                                                                                 |  1.2490878 | 0.0375275 |
| <GO:0031640> | killing of cells of another organism                                                                                                             |  1.7257654 | 0.0376510 |
| <GO:0141061> | disruption of cell in another organism                                                                                                           |  1.7257654 | 0.0376510 |
| <GO:0045766> | positive regulation of angiogenesis                                                                                                              |  1.5763666 | 0.0376510 |
| <GO:1904018> | positive regulation of vasculature development                                                                                                   |  1.5763666 | 0.0376510 |
| <GO:0043410> | positive regulation of MAPK cascade                                                                                                              |  1.4266714 | 0.0381057 |
| <GO:0002429> | immune response-activating cell surface receptor signaling pathway                                                                               |  1.4170521 | 0.0381057 |
| <GO:0002697> | regulation of immune effector process                                                                                                            |  1.4111633 | 0.0381057 |
| <GO:0034097> | response to cytokine                                                                                                                             |  1.2818026 | 0.0381057 |
| <GO:0002821> | positive regulation of adaptive immune response                                                                                                  |  1.5716253 | 0.0389699 |
| <GO:0016477> | cell migration                                                                                                                                   |  1.2533005 | 0.0389699 |
| <GO:0015718> | monocarboxylic acid transport                                                                                                                    |  1.6044321 | 0.0393772 |
| <GO:0043588> | skin development                                                                                                                                 |  1.5272329 | 0.0393785 |
| <GO:0071345> | cellular response to cytokine stimulus                                                                                                           |  1.2938392 | 0.0393785 |
| <GO:0010232> | vascular transport                                                                                                                               |  1.6786193 | 0.0401850 |
| <GO:0150104> | transport across blood-brain barrier                                                                                                             |  1.6786193 | 0.0401850 |
| <GO:0042088> | T-helper 1 type immune response                                                                                                                  |  1.6587302 | 0.0401850 |
| <GO:0048514> | blood vessel morphogenesis                                                                                                                       |  1.3712769 | 0.0401850 |
| <GO:0002521> | leukocyte differentiation                                                                                                                        |  1.3366073 | 0.0401850 |
| <GO:1901652> | response to peptide                                                                                                                              |  1.2797657 | 0.0401850 |
| <GO:0002444> | myeloid leukocyte mediated immunity                                                                                                              |  1.5384006 | 0.0409746 |
| <GO:0002699> | positive regulation of immune effector process                                                                                                   |  1.4826276 | 0.0409746 |
| <GO:0019882> | antigen processing and presentation                                                                                                              |  1.5927145 | 0.0414720 |
| <GO:0051966> | regulation of synaptic transmission, glutamatergic                                                                                               |  1.6687274 | 0.0421299 |
| <GO:0050829> | defense response to Gram-negative bacterium                                                                                                      |  1.6932121 | 0.0423562 |
| <GO:0006120> | mitochondrial electron transport, NADH to ubiquinone                                                                                             | -2.0998283 | 0.0427655 |
| <GO:1902106> | negative regulation of leukocyte differentiation                                                                                                 |  1.5984243 | 0.0427655 |
| <GO:0002526> | acute inflammatory response                                                                                                                      |  1.6394696 | 0.0427948 |
| <GO:2001235> | positive regulation of apoptotic signaling pathway                                                                                               |  1.5488255 | 0.0427948 |
| <GO:0032715> | negative regulation of interleukin-6 production                                                                                                  |  1.7045426 | 0.0434628 |
| <GO:2001057> | reactive nitrogen species metabolic process                                                                                                      |  1.6009042 | 0.0445992 |
| <GO:0043409> | negative regulation of MAPK cascade                                                                                                              |  1.5846964 | 0.0445992 |
| <GO:0002696> | positive regulation of leukocyte activation                                                                                                      |  1.4331068 | 0.0445992 |
| <GO:0002758> | innate immune response-activating signaling pathway                                                                                              |  1.3699578 | 0.0446689 |
| <GO:0051606> | detection of stimulus                                                                                                                            |  1.5664101 | 0.0461606 |
| <GO:0098754> | detoxification                                                                                                                                   |  1.5868723 | 0.0464884 |
| <GO:0034121> | regulation of toll-like receptor signaling pathway                                                                                               |  1.6942734 | 0.0468131 |
| <GO:0061900> | glial cell activation                                                                                                                            |  1.5776321 | 0.0484337 |
| <GO:0038094> | Fc-gamma receptor signaling pathway                                                                                                              |  1.6601777 | 0.0485703 |
| <GO:0002825> | regulation of T-helper 1 type immune response                                                                                                    |  1.6454623 | 0.0485703 |
| <GO:0009913> | epidermal cell differentiation                                                                                                                   |  1.5504135 | 0.0485703 |
| <GO:0002822> | regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains          |  1.5128276 | 0.0485703 |
| <GO:0051050> | positive regulation of transport                                                                                                                 |  1.3004694 | 0.0485703 |
| <GO:0042092> | type 2 immune response                                                                                                                           |  1.6910163 | 0.0489883 |
| <GO:0002828> | regulation of type 2 immune response                                                                                                             |  1.6434698 | 0.0504721 |
| <GO:0046209> | nitric oxide metabolic process                                                                                                                   |  1.6024404 | 0.0504721 |
| <GO:0141060> | disruption of anatomical structure in another organism                                                                                           |  1.6811001 | 0.0504973 |
| <GO:0009620> | response to fungus                                                                                                                               |  1.6204966 | 0.0504973 |
| <GO:0038093> | Fc receptor signaling pathway                                                                                                                    |  1.6472371 | 0.0541426 |
| <GO:0010720> | positive regulation of cell development                                                                                                          |  1.3784924 | 0.0541426 |
| <GO:0002221> | pattern recognition receptor signaling pathway                                                                                                   |  1.3594016 | 0.0546016 |
| <GO:0062207> | regulation of pattern recognition receptor signaling pathway                                                                                     |  1.4531097 | 0.0560651 |
| <GO:0002573> | myeloid leukocyte differentiation                                                                                                                |  1.4292119 | 0.0560651 |
| <GO:0030856> | regulation of epithelial cell differentiation                                                                                                    |  1.5769538 | 0.0569299 |
| <GO:0098739> | import across plasma membrane                                                                                                                    |  1.5341619 | 0.0575615 |
| <GO:0015711> | organic anion transport                                                                                                                          |  1.4243154 | 0.0580688 |
| <GO:0001568> | blood vessel development                                                                                                                         |  1.3213467 | 0.0584583 |
| <GO:0008284> | positive regulation of cell population proliferation                                                                                             |  1.3037362 | 0.0584583 |
| <GO:0002446> | neutrophil mediated immunity                                                                                                                     |  1.6331789 | 0.0592446 |
| <GO:0010573> | vascular endothelial growth factor production                                                                                                    |  1.6644986 | 0.0614934 |
| <GO:0009653> | anatomical structure morphogenesis                                                                                                               |  1.1991559 | 0.0614934 |
| <GO:0051093> | negative regulation of developmental process                                                                                                     |  1.3160880 | 0.0620126 |
| <GO:0007188> | adenylate cyclase-modulating G protein-coupled receptor signaling pathway                                                                        |  1.5927877 | 0.0629821 |
| <GO:0010721> | negative regulation of cell development                                                                                                          |  1.4519174 | 0.0629821 |
| <GO:0044283> | small molecule biosynthetic process                                                                                                              |  1.3716102 | 0.0629821 |
| <GO:0002761> | regulation of myeloid leukocyte differentiation                                                                                                  |  1.5292340 | 0.0632076 |
| <GO:0030097> | hemopoiesis                                                                                                                                      |  1.2537145 | 0.0632076 |
| <GO:0003008> | system process                                                                                                                                   |  1.2471546 | 0.0632076 |
| <GO:0060627> | regulation of vesicle-mediated transport                                                                                                         |  1.3205428 | 0.0633411 |
| <GO:0010647> | positive regulation of cell communication                                                                                                        |  1.1949940 | 0.0633411 |
| <GO:0071347> | cellular response to interleukin-1                                                                                                               |  1.5962764 | 0.0638182 |
| <GO:0045429> | positive regulation of nitric oxide biosynthetic process                                                                                         |  1.6379795 | 0.0640335 |
| <GO:0008544> | epidermis development                                                                                                                            |  1.4843759 | 0.0640335 |
| <GO:0009967> | positive regulation of signal transduction                                                                                                       |  1.2203635 | 0.0640335 |
| <GO:0072359> | circulatory system development                                                                                                                   |  1.2782293 | 0.0645624 |
| <GO:0045670> | regulation of osteoclast differentiation                                                                                                         |  1.5942878 | 0.0649896 |
| <GO:0032543> | mitochondrial translation                                                                                                                        | -1.9604158 | 0.0651602 |
| <GO:0015908> | fatty acid transport                                                                                                                             |  1.5541463 | 0.0651602 |
| <GO:0014911> | positive regulation of smooth muscle cell migration                                                                                              |  1.6096600 | 0.0651876 |
| <GO:0006898> | receptor-mediated endocytosis                                                                                                                    |  1.3768014 | 0.0700264 |
| <GO:0048708> | astrocyte differentiation                                                                                                                        |  1.5757223 | 0.0714289 |
| <GO:0032755> | positive regulation of interleukin-6 production                                                                                                  |  1.5091168 | 0.0714289 |
| <GO:0001944> | vasculature development                                                                                                                          |  1.3152407 | 0.0714289 |
| <GO:0046903> | secretion                                                                                                                                        |  1.2629945 | 0.0714289 |
| <GO:0045765> | regulation of angiogenesis                                                                                                                       |  1.4567895 | 0.0718935 |
| <GO:0002292> | T cell differentiation involved in immune response                                                                                               |  1.5387126 | 0.0719469 |
| <GO:0098609> | cell-cell adhesion                                                                                                                               |  1.3046571 | 0.0729517 |
| <GO:0002283> | neutrophil activation involved in immune response                                                                                                |  1.6206026 | 0.0734333 |
| <GO:0002366> | leukocyte activation involved in immune response                                                                                                 |  1.3843168 | 0.0734333 |
| <GO:0048870> | cell motility                                                                                                                                    |  1.2305335 | 0.0736825 |
| <GO:0014909> | smooth muscle cell migration                                                                                                                     |  1.5948386 | 0.0743758 |
| <GO:0014910> | regulation of smooth muscle cell migration                                                                                                       |  1.5948386 | 0.0743758 |
| <GO:0014812> | muscle cell migration                                                                                                                            |  1.5846747 | 0.0743758 |
| <GO:0050871> | positive regulation of B cell activation                                                                                                         |  1.5803547 | 0.0743758 |
| <GO:0002695> | negative regulation of leukocyte activation                                                                                                      |  1.5043589 | 0.0743758 |
| <GO:0040017> | positive regulation of locomotion                                                                                                                |  1.3183062 | 0.0743758 |
| <GO:0043122> | regulation of canonical NF-kappaB signal transduction                                                                                            |  1.3721532 | 0.0747851 |
| <GO:0042100> | B cell proliferation                                                                                                                             |  1.5392315 | 0.0748563 |
| <GO:0032760> | positive regulation of tumor necrosis factor production                                                                                          |  1.4900726 | 0.0759205 |
| <GO:0001659> | temperature homeostasis                                                                                                                          |  1.4604346 | 0.0764074 |
| <GO:0002263> | cell activation involved in immune response                                                                                                      |  1.3831113 | 0.0775238 |
| <GO:0007193> | adenylate cyclase-inhibiting G protein-coupled receptor signaling pathway                                                                        |  1.5858432 | 0.0779741 |
| <GO:0006801> | superoxide metabolic process                                                                                                                     |  1.5306577 | 0.0822336 |
| <GO:0043062> | extracellular structure organization                                                                                                             |  1.4909941 | 0.0822336 |
| <GO:1901342> | regulation of vasculature development                                                                                                            |  1.4379271 | 0.0822336 |
| <GO:1903557> | positive regulation of tumor necrosis factor superfamily cytokine production                                                                     |  1.4649860 | 0.0846910 |
| <GO:0042554> | superoxide anion generation                                                                                                                      |  1.5810174 | 0.0871439 |
| <GO:0000165> | MAPK cascade                                                                                                                                     |  1.2927656 | 0.0871439 |
| <GO:0070373> | negative regulation of ERK1 and ERK2 cascade                                                                                                     |  1.6207599 | 0.0913549 |
| <GO:0048002> | antigen processing and presentation of peptide antigen                                                                                           |  1.5930644 | 0.0913549 |
| <GO:0032757> | positive regulation of interleukin-8 production                                                                                                  |  1.5490381 | 0.0913549 |
| <GO:0001525> | angiogenesis                                                                                                                                     |  1.3430031 | 0.0913549 |
| <GO:0009611> | response to wounding                                                                                                                             |  1.3008949 | 0.0913549 |
| <GO:0048878> | chemical homeostasis                                                                                                                             |  1.2555479 | 0.0924763 |
| <GO:0030216> | keratinocyte differentiation                                                                                                                     |  1.5088297 | 0.0929200 |
| <GO:0046651> | lymphocyte proliferation                                                                                                                         |  1.3799851 | 0.0929200 |
| <GO:0007417> | central nervous system development                                                                                                               |  1.2684381 | 0.0929200 |
| <GO:0032720> | negative regulation of tumor necrosis factor production                                                                                          |  1.5512549 | 0.0936296 |
| <GO:1903556> | negative regulation of tumor necrosis factor superfamily cytokine production                                                                     |  1.5512549 | 0.0936296 |
| <GO:0046328> | regulation of JNK cascade                                                                                                                        |  1.4497716 | 0.0949187 |
| <GO:1904407> | positive regulation of nitric oxide metabolic process                                                                                            |  1.5815176 | 0.0957250 |
| <GO:0042592> | homeostatic process                                                                                                                              |  1.1976949 | 0.0957250 |
| <GO:0007600> | sensory perception                                                                                                                               |  1.4437191 | 0.0958225 |
| <GO:0007254> | JNK cascade                                                                                                                                      |  1.4380816 | 0.0958225 |
| <GO:0007249> | canonical NF-kappaB signal transduction                                                                                                          |  1.3447652 | 0.0958225 |
| <GO:0030335> | positive regulation of cell migration                                                                                                            |  1.3111094 | 0.0958225 |
| <GO:0032946> | positive regulation of mononuclear cell proliferation                                                                                            |  1.4690072 | 0.0961176 |
| <GO:0050671> | positive regulation of lymphocyte proliferation                                                                                                  |  1.4690072 | 0.0961176 |
| <GO:0019221> | cytokine-mediated signaling pathway                                                                                                              |  1.3141256 | 0.0961176 |
| <GO:0045637> | regulation of myeloid cell differentiation                                                                                                       |  1.4214954 | 0.0965002 |
| <GO:2000147> | positive regulation of cell motility                                                                                                             |  1.3215508 | 0.0965002 |
| <GO:0023056> | positive regulation of signaling                                                                                                                 |  1.1880363 | 0.0965002 |
| <GO:0070371> | ERK1 and ERK2 cascade                                                                                                                            |  1.4004769 | 0.0985173 |
| <GO:0043123> | positive regulation of canonical NF-kappaB signal transduction                                                                                   |  1.3750862 | 0.0985173 |
| <GO:0141085> | regulation of inflammasome-mediated signaling pathway                                                                                            |  1.5553827 | 0.1005377 |
| <GO:1900225> | regulation of NLRP3 inflammasome complex assembly                                                                                                |  1.5553827 | 0.1005377 |
| <GO:0032944> | regulation of mononuclear cell proliferation                                                                                                     |  1.4051265 | 0.1027585 |
| <GO:0050670> | regulation of lymphocyte proliferation                                                                                                           |  1.4051265 | 0.1027585 |
| <GO:0030500> | regulation of bone mineralization                                                                                                                |  1.5393526 | 0.1049920 |
| <GO:0002763> | positive regulation of myeloid leukocyte differentiation                                                                                         |  1.5379946 | 0.1049920 |
| <GO:0050920> | regulation of chemotaxis                                                                                                                         |  1.4431963 | 0.1058964 |
| <GO:0043277> | apoptotic cell clearance                                                                                                                         |  1.6046580 | 0.1065678 |
| <GO:0048145> | regulation of fibroblast proliferation                                                                                                           |  1.5733071 | 0.1075170 |
| <GO:0070374> | positive regulation of ERK1 and ERK2 cascade                                                                                                     |  1.4368481 | 0.1075170 |
| <GO:0070663> | regulation of leukocyte proliferation                                                                                                            |  1.3909797 | 0.1082823 |
| <GO:1905039> | carboxylic acid transmembrane transport                                                                                                          |  1.5100463 | 0.1083482 |
| <GO:0070269> | pyroptotic inflammatory response                                                                                                                 |  1.5436497 | 0.1088463 |
| <GO:0050864> | regulation of B cell activation                                                                                                                  |  1.4718458 | 0.1088463 |
| <GO:0006805> | xenobiotic metabolic process                                                                                                                     |  1.5647943 | 0.1091168 |
| <GO:0070372> | regulation of ERK1 and ERK2 cascade                                                                                                              |  1.4400488 | 0.1091168 |
| <GO:0032943> | mononuclear cell proliferation                                                                                                                   |  1.3744021 | 0.1100713 |
| <GO:0150076> | neuroinflammatory response                                                                                                                       |  1.4904061 | 0.1102563 |
| <GO:0032602> | chemokine production                                                                                                                             |  1.4583786 | 0.1145453 |
| <GO:0032642> | regulation of chemokine production                                                                                                               |  1.4583786 | 0.1145453 |
| <GO:0032743> | positive regulation of interleukin-2 production                                                                                                  |  1.5898742 | 0.1159185 |
| <GO:0090022> | regulation of neutrophil chemotaxis                                                                                                              |  1.5340569 | 0.1165493 |
| <GO:0007584> | response to nutrient                                                                                                                             |  1.4407674 | 0.1165493 |
| <GO:0016053> | organic acid biosynthetic process                                                                                                                |  1.4316359 | 0.1165493 |
| <GO:0046394> | carboxylic acid biosynthetic process                                                                                                             |  1.4316359 | 0.1165493 |
| <GO:0050878> | regulation of body fluid levels                                                                                                                  |  1.3759658 | 0.1180017 |
| <GO:0008285> | negative regulation of cell population proliferation                                                                                             |  1.3212841 | 0.1180017 |
| <GO:0045948> | positive regulation of translational initiation                                                                                                  | -1.6729340 | 0.1189852 |
| <GO:0072001> | renal system development                                                                                                                         |  1.4434655 | 0.1189852 |
| <GO:0002269> | leukocyte activation involved in inflammatory response                                                                                           |  1.5215868 | 0.1196845 |
| <GO:0032615> | interleukin-12 production                                                                                                                        |  1.5210345 | 0.1196845 |
| <GO:0032655> | regulation of interleukin-12 production                                                                                                          |  1.5210345 | 0.1196845 |
| <GO:0030198> | extracellular matrix organization                                                                                                                |  1.4858482 | 0.1196845 |
| <GO:0045229> | external encapsulating structure organization                                                                                                    |  1.4858482 | 0.1196845 |
| <GO:0045576> | mast cell activation                                                                                                                             |  1.4720557 | 0.1213409 |
| <GO:0031664> | regulation of lipopolysaccharide-mediated signaling pathway                                                                                      |  1.5796040 | 0.1223744 |
| <GO:0003013> | circulatory system process                                                                                                                       |  1.3623217 | 0.1229649 |
| <GO:0140053> | mitochondrial gene expression                                                                                                                    | -1.7238042 | 0.1230649 |
| <GO:0048066> | developmental pigmentation                                                                                                                       |  1.5837841 | 0.1230649 |
| <GO:0045444> | fat cell differentiation                                                                                                                         |  1.4256719 | 0.1230649 |
| <GO:0070661> | leukocyte proliferation                                                                                                                          |  1.3462462 | 0.1230649 |
| <GO:0007162> | negative regulation of cell adhesion                                                                                                             |  1.3827585 | 0.1236302 |
| <GO:0032370> | positive regulation of lipid transport                                                                                                           |  1.5503595 | 0.1240748 |
| <GO:0045428> | regulation of nitric oxide biosynthetic process                                                                                                  |  1.5010696 | 0.1241152 |
| <GO:0050866> | negative regulation of cell activation                                                                                                           |  1.4225232 | 0.1274782 |
| <GO:0032692> | negative regulation of interleukin-1 production                                                                                                  |  1.5700024 | 0.1280848 |
| <GO:0006518> | peptide metabolic process                                                                                                                        |  1.5470160 | 0.1280848 |
| <GO:0045778> | positive regulation of ossification                                                                                                              |  1.5180252 | 0.1280848 |
| <GO:0051251> | positive regulation of lymphocyte activation                                                                                                     |  1.3489030 | 0.1280848 |
| <GO:0042310> | vasoconstriction                                                                                                                                 |  1.5347942 | 0.1307644 |
| <GO:0045624> | positive regulation of T-helper cell differentiation                                                                                             |  1.5654432 | 0.1307899 |
| <GO:0002819> | regulation of adaptive immune response                                                                                                           |  1.4001766 | 0.1316829 |
| <GO:0071826> | protein-RNA complex organization                                                                                                                 | -1.4708629 | 0.1318103 |
| <GO:0051090> | regulation of DNA-binding transcription factor activity                                                                                          |  1.3899415 | 0.1335267 |
| <GO:0045622> | regulation of T-helper cell differentiation                                                                                                      |  1.5241399 | 0.1348242 |
| <GO:0002534> | cytokine production involved in inflammatory response                                                                                            |  1.4758784 | 0.1373082 |
| <GO:1900015> | regulation of cytokine production involved in inflammatory response                                                                              |  1.4758784 | 0.1373082 |
| <GO:0097028> | dendritic cell differentiation                                                                                                                   |  1.5178993 | 0.1408616 |
| <GO:0009593> | detection of chemical stimulus                                                                                                                   |  1.5350995 | 0.1411782 |
| <GO:0051249> | regulation of lymphocyte activation                                                                                                              |  1.3021010 | 0.1421009 |
| <GO:0048871> | multicellular organismal-level homeostasis                                                                                                       |  1.2528070 | 0.1421009 |
| <GO:0035150> | regulation of tube size                                                                                                                          |  1.4847069 | 0.1424376 |
| <GO:0035296> | regulation of tube diameter                                                                                                                      |  1.4847069 | 0.1424376 |
| <GO:0097746> | blood vessel diameter maintenance                                                                                                                |  1.4847069 | 0.1424376 |
| <GO:0033273> | response to vitamin                                                                                                                              |  1.4970773 | 0.1426045 |
| <GO:0030970> | retrograde protein transport, ER to cytosol                                                                                                      | -1.7159847 | 0.1437815 |
| <GO:1903513> | endoplasmic reticulum to cytosol transport                                                                                                       | -1.7159847 | 0.1437815 |
| <GO:0070646> | protein modification by small protein removal                                                                                                    | -1.5918988 | 0.1437815 |
| <GO:0070498> | interleukin-1-mediated signaling pathway                                                                                                         |  1.5588319 | 0.1437815 |
| <GO:0001523> | retinoid metabolic process                                                                                                                       |  1.5479096 | 0.1437815 |
| <GO:0002437> | inflammatory response to antigenic stimulus                                                                                                      |  1.4913032 | 0.1437815 |
| <GO:0030225> | macrophage differentiation                                                                                                                       |  1.4884445 | 0.1437815 |
| <GO:0001774> | microglial cell activation                                                                                                                       |  1.4872765 | 0.1437815 |
| <GO:0070665> | positive regulation of leukocyte proliferation                                                                                                   |  1.4382788 | 0.1437815 |
| <GO:0060429> | epithelium development                                                                                                                           |  1.2436749 | 0.1437815 |
| <GO:0061138> | morphogenesis of a branching epithelium                                                                                                          |  1.4791501 | 0.1443923 |
| <GO:0031623> | receptor internalization                                                                                                                         |  1.4045101 | 0.1447613 |
| <GO:0002443> | leukocyte mediated immunity                                                                                                                      |  1.3016519 | 0.1464178 |
| <GO:1990748> | cellular detoxification                                                                                                                          |  1.4644670 | 0.1464892 |
| <GO:0006690> | icosanoid metabolic process                                                                                                                      |  1.4496928 | 0.1464892 |
| <GO:0051130> | positive regulation of cellular component organization                                                                                           |  1.2126566 | 0.1473947 |
| <GO:0042476> | odontogenesis                                                                                                                                    |  1.4806166 | 0.1493785 |
| <GO:0002532> | production of molecular mediator involved in inflammatory response                                                                               |  1.4219384 | 0.1512347 |
| <GO:0072593> | reactive oxygen species metabolic process                                                                                                        |  1.3314160 | 0.1526870 |
| <GO:1905954> | positive regulation of lipid localization                                                                                                        |  1.5152818 | 0.1528120 |
| <GO:0032609> | type II interferon production                                                                                                                    |  1.4939523 | 0.1528120 |
| <GO:0032649> | regulation of type II interferon production                                                                                                      |  1.4939523 | 0.1528120 |
| <GO:0032623> | interleukin-2 production                                                                                                                         |  1.4749848 | 0.1530816 |
| <GO:0032663> | regulation of interleukin-2 production                                                                                                           |  1.4749848 | 0.1530816 |
| <GO:0071674> | mononuclear cell migration                                                                                                                       |  1.3729717 | 0.1546472 |
| <GO:0030278> | regulation of ossification                                                                                                                       |  1.4926171 | 0.1551647 |
| <GO:0030099> | myeloid cell differentiation                                                                                                                     |  1.2654099 | 0.1579078 |
| <GO:1903825> | organic acid transmembrane transport                                                                                                             |  1.4506214 | 0.1582785 |
| <GO:0060760> | positive regulation of response to cytokine stimulus                                                                                             |  1.4612963 | 0.1589759 |
| <GO:0002700> | regulation of production of molecular mediator of immune response                                                                                |  1.3893121 | 0.1590751 |
| <GO:0008033> | tRNA processing                                                                                                                                  | -1.5796426 | 0.1610282 |
| <GO:0006968> | cellular defense response                                                                                                                        |  1.4921033 | 0.1610282 |
| <GO:0007599> | hemostasis                                                                                                                                       |  1.3728563 | 0.1610282 |
| <GO:0006082> | organic acid metabolic process                                                                                                                   |  1.2664146 | 0.1610282 |
| <GO:0043436> | oxoacid metabolic process                                                                                                                        |  1.2664146 | 0.1610282 |
| <GO:0072006> | nephron development                                                                                                                              |  1.4765113 | 0.1615749 |
| <GO:0001822> | kidney development                                                                                                                               |  1.4151271 | 0.1634438 |
| <GO:1903131> | mononuclear cell differentiation                                                                                                                 |  1.2706556 | 0.1634438 |
| <GO:0030155> | regulation of cell adhesion                                                                                                                      |  1.2330078 | 0.1634438 |
| <GO:0000245> | spliceosomal complex assembly                                                                                                                    | -1.6863894 | 0.1637724 |
| <GO:0050731> | positive regulation of peptidyl-tyrosine phosphorylation                                                                                         |  1.5289300 | 0.1646293 |
| <GO:0046173> | polyol biosynthetic process                                                                                                                      |  1.5007361 | 0.1678849 |
| <GO:0032368> | regulation of lipid transport                                                                                                                    |  1.4500782 | 0.1678849 |
| <GO:0070167> | regulation of biomineral tissue development                                                                                                      |  1.4477150 | 0.1678849 |
| <GO:2000379> | positive regulation of reactive oxygen species metabolic process                                                                                 |  1.4321146 | 0.1678849 |
| <GO:0062208> | positive regulation of pattern recognition receptor signaling pathway                                                                            |  1.4067719 | 0.1678849 |
| <GO:0001906> | cell killing                                                                                                                                     |  1.3736525 | 0.1678849 |
| <GO:0019752> | carboxylic acid metabolic process                                                                                                                |  1.2498692 | 0.1678849 |
| <GO:0016125> | sterol metabolic process                                                                                                                         |  1.4555348 | 0.1683078 |
| <GO:0034113> | heterotypic cell-cell adhesion                                                                                                                   |  1.4621477 | 0.1694251 |
| <GO:0061028> | establishment of endothelial barrier                                                                                                             |  1.4663967 | 0.1714773 |
| <GO:1905952> | regulation of lipid localization                                                                                                                 |  1.4323851 | 0.1714773 |
| <GO:0002275> | myeloid cell activation involved in immune response                                                                                              |  1.4208947 | 0.1714773 |
| <GO:0002286> | T cell activation involved in immune response                                                                                                    |  1.3870017 | 0.1714773 |
| <GO:0030501> | positive regulation of bone mineralization                                                                                                       |  1.5290547 | 0.1733457 |
| <GO:0022618> | protein-RNA complex assembly                                                                                                                     | -1.4873591 | 0.1735078 |
| <GO:0008347> | glial cell migration                                                                                                                             |  1.4483300 | 0.1736878 |
| <GO:0007585> | respiratory gaseous exchange by respiratory system                                                                                               |  1.4836832 | 0.1740120 |
| <GO:0038034> | signal transduction in absence of ligand                                                                                                         |  1.4798928 | 0.1740120 |
| <GO:0097192> | extrinsic apoptotic signaling pathway in absence of ligand                                                                                       |  1.4798928 | 0.1740120 |
| <GO:0097237> | cellular response to toxic substance                                                                                                             |  1.4372019 | 0.1753763 |
| <GO:0032691> | negative regulation of interleukin-1 beta production                                                                                             |  1.5065054 | 0.1805818 |
| <GO:0035249> | synaptic transmission, glutamatergic                                                                                                             |  1.4857484 | 0.1805818 |
| <GO:0034219> | carbohydrate transmembrane transport                                                                                                             |  1.4297984 | 0.1805818 |
| <GO:0051452> | intracellular pH reduction                                                                                                                       |  1.4287834 | 0.1805818 |
| <GO:0002688> | regulation of leukocyte chemotaxis                                                                                                               |  1.4157340 | 0.1805818 |
| <GO:0120162> | positive regulation of cold-induced thermogenesis                                                                                                |  1.4576423 | 0.1818018 |
| <GO:0090025> | regulation of monocyte chemotaxis                                                                                                                |  1.5008103 | 0.1820656 |
| <GO:2001244> | positive regulation of intrinsic apoptotic signaling pathway                                                                                     |  1.5007119 | 0.1820656 |
| <GO:0034314> | Arp2/3 complex-mediated actin nucleation                                                                                                         | -1.5486012 | 0.1829531 |
| <GO:0022603> | regulation of anatomical structure morphogenesis                                                                                                 |  1.2361575 | 0.1829531 |
| <GO:0006911> | phagocytosis, engulfment                                                                                                                         |  1.4638114 | 0.1829876 |
| <GO:0019751> | polyol metabolic process                                                                                                                         |  1.4696733 | 0.1834219 |
| <GO:0080164> | regulation of nitric oxide metabolic process                                                                                                     |  1.4492466 | 0.1834219 |
| <GO:0045806> | negative regulation of endocytosis                                                                                                               |  1.4334785 | 0.1834219 |
| <GO:0044546> | NLRP3 inflammasome complex assembly                                                                                                              |  1.4367252 | 0.1844816 |
| <GO:0001763> | morphogenesis of a branching structure                                                                                                           |  1.4588341 | 0.1856369 |
| <GO:0045834> | positive regulation of lipid metabolic process                                                                                                   |  1.4575951 | 0.1856369 |
| <GO:0055082> | intracellular chemical homeostasis                                                                                                               |  1.2359889 | 0.1857159 |
| <GO:1900016> | negative regulation of cytokine production involved in inflammatory response                                                                     |  1.5138550 | 0.1866769 |
| <GO:0050931> | pigment cell differentiation                                                                                                                     |  1.4951333 | 0.1867495 |
| <GO:1902107> | positive regulation of leukocyte differentiation                                                                                                 |  1.3698945 | 0.1867495 |
| <GO:1903708> | positive regulation of hemopoiesis                                                                                                               |  1.3698945 | 0.1867495 |
| <GO:0009887> | animal organ morphogenesis                                                                                                                       |  1.2666572 | 0.1867495 |
| <GO:0051145> | smooth muscle cell differentiation                                                                                                               |  1.4811835 | 0.1872740 |
| <GO:0042093> | T-helper cell differentiation                                                                                                                    |  1.4536943 | 0.1872740 |
| <GO:0046456> | icosanoid biosynthetic process                                                                                                                   |  1.4525739 | 0.1872740 |
| <GO:0034341> | response to type II interferon                                                                                                                   |  1.3876276 | 0.1872740 |
| <GO:0006066> | alcohol metabolic process                                                                                                                        |  1.3367007 | 0.1872740 |
| <GO:0141124> | intracellular signaling cassette                                                                                                                 |  1.1404736 | 0.1872740 |
| <GO:0045672> | positive regulation of osteoclast differentiation                                                                                                |  1.4910359 | 0.1897858 |
| <GO:0032102> | negative regulation of response to external stimulus                                                                                             |  1.3139144 | 0.1906154 |
| <GO:0002367> | cytokine production involved in immune response                                                                                                  |  1.3853596 | 0.1915162 |
| <GO:0002718> | regulation of cytokine production involved in immune response                                                                                    |  1.3853596 | 0.1915162 |
| <GO:0002702> | positive regulation of production of molecular mediator of immune response                                                                       |  1.3849622 | 0.1915162 |
| <GO:0006890> | retrograde vesicle-mediated transport, Golgi to endoplasmic reticulum                                                                            | -1.6698186 | 0.1941876 |
| <GO:0032613> | interleukin-10 production                                                                                                                        |  1.4346038 | 0.1941876 |
| <GO:0032653> | regulation of interleukin-10 production                                                                                                          |  1.4346038 | 0.1941876 |
| <GO:0070169> | positive regulation of biomineral tissue development                                                                                             |  1.4636816 | 0.1943475 |
| <GO:0048659> | smooth muscle cell proliferation                                                                                                                 |  1.4184494 | 0.1943475 |
| <GO:0048660> | regulation of smooth muscle cell proliferation                                                                                                   |  1.4184494 | 0.1943475 |
| <GO:0031960> | response to corticosteroid                                                                                                                       |  1.4151711 | 0.1943475 |
| <GO:0035265> | organ growth                                                                                                                                     |  1.4139824 | 0.1943475 |
| <GO:0051046> | regulation of secretion                                                                                                                          |  1.2671736 | 0.1943475 |
| <GO:0006081> | aldehyde metabolic process                                                                                                                       |  1.5080258 | 0.1948243 |
| <GO:0046928> | regulation of neurotransmitter secretion                                                                                                         |  1.4653296 | 0.1948243 |
| <GO:0002705> | positive regulation of leukocyte mediated immunity                                                                                               |  1.4299393 | 0.1948243 |
| <GO:0046466> | membrane lipid catabolic process                                                                                                                 |  1.4280067 | 0.1948243 |
| <GO:0032963> | collagen metabolic process                                                                                                                       |  1.4273977 | 0.1948243 |
| <GO:0072073> | kidney epithelium development                                                                                                                    |  1.5040005 | 0.1953867 |
| <GO:0150077> | regulation of neuroinflammatory response                                                                                                         |  1.4852537 | 0.1953867 |
| <GO:0106106> | cold-induced thermogenesis                                                                                                                       |  1.3818925 | 0.1953867 |
| <GO:0097191> | extrinsic apoptotic signaling pathway                                                                                                            |  1.3123009 | 0.1953867 |
| <GO:0046649> | lymphocyte activation                                                                                                                            |  1.1945464 | 0.1953867 |
| <GO:0071346> | cellular response to type II interferon                                                                                                          |  1.4005938 | 0.1997606 |
| <GO:0050728> | negative regulation of inflammatory response                                                                                                     |  1.3662271 | 0.1997606 |
| <GO:0051962> | positive regulation of nervous system development                                                                                                |  1.3486464 | 0.1997606 |
| <GO:0006903> | vesicle targeting                                                                                                                                | -1.5199891 | 0.2013310 |
| <GO:0045596> | negative regulation of cell differentiation                                                                                                      |  1.2655497 | 0.2033351 |
| <GO:0065009> | regulation of molecular function                                                                                                                 |  1.1980867 | 0.2033351 |
| <GO:0045669> | positive regulation of osteoblast differentiation                                                                                                |  1.4530452 | 0.2034863 |
| <GO:0008643> | carbohydrate transport                                                                                                                           |  1.3993664 | 0.2034863 |
| <GO:0002708> | positive regulation of lymphocyte mediated immunity                                                                                              |  1.4123178 | 0.2042849 |
| <GO:0016579> | protein deubiquitination                                                                                                                         | -1.5856466 | 0.2044537 |
| <GO:0045638> | negative regulation of myeloid cell differentiation                                                                                              |  1.4171672 | 0.2044537 |
| <GO:0051250> | negative regulation of lymphocyte activation                                                                                                     |  1.3855387 | 0.2054322 |
| <GO:0051453> | regulation of intracellular pH                                                                                                                   |  1.3735360 | 0.2054322 |
| <GO:2000377> | regulation of reactive oxygen species metabolic process                                                                                          |  1.3546862 | 0.2054322 |
| <GO:0048771> | tissue remodeling                                                                                                                                |  1.3512652 | 0.2054322 |
| <GO:0051049> | regulation of transport                                                                                                                          |  1.1609444 | 0.2054322 |
| <GO:0006654> | phosphatidic acid biosynthetic process                                                                                                           |  1.4678569 | 0.2086177 |
| <GO:0046473> | phosphatidic acid metabolic process                                                                                                              |  1.4678569 | 0.2086177 |
| <GO:0033002> | muscle cell proliferation                                                                                                                        |  1.3676980 | 0.2086177 |
| <GO:0006885> | regulation of pH                                                                                                                                 |  1.3676873 | 0.2086177 |
| <GO:0070555> | response to interleukin-1                                                                                                                        |  1.3935639 | 0.2088003 |
| <GO:0007595> | lactation                                                                                                                                        |  1.4616096 | 0.2091498 |
| <GO:0035751> | regulation of lysosomal lumen pH                                                                                                                 |  1.4067735 | 0.2104726 |
| <GO:0045744> | negative regulation of G protein-coupled receptor signaling pathway                                                                              |  1.4163098 | 0.2112519 |
| <GO:0002285> | lymphocyte activation involved in immune response                                                                                                |  1.3096132 | 0.2112519 |
| <GO:0001501> | skeletal system development                                                                                                                      |  1.3062935 | 0.2112519 |
| <GO:0046329> | negative regulation of JNK cascade                                                                                                               |  1.4683796 | 0.2134420 |
| <GO:0032620> | interleukin-17 production                                                                                                                        |  1.4335473 | 0.2138130 |
| <GO:0032660> | regulation of interleukin-17 production                                                                                                          |  1.4335473 | 0.2138130 |
| <GO:0048144> | fibroblast proliferation                                                                                                                         |  1.4108443 | 0.2138130 |
| <GO:0002287> | alpha-beta T cell activation involved in immune response                                                                                         |  1.3934077 | 0.2138130 |
| <GO:0002293> | alpha-beta T cell differentiation involved in immune response                                                                                    |  1.3934077 | 0.2138130 |
| <GO:0002294> | CD4-positive, alpha-beta T cell differentiation involved in immune response                                                                      |  1.3934077 | 0.2138130 |
| <GO:0008360> | regulation of cell shape                                                                                                                         |  1.3457107 | 0.2138130 |
| <GO:0022407> | regulation of cell-cell adhesion                                                                                                                 |  1.2595718 | 0.2138130 |
| <GO:0007035> | vacuolar acidification                                                                                                                           |  1.4082599 | 0.2151158 |
| <GO:0032729> | positive regulation of type II interferon production                                                                                             |  1.4074432 | 0.2151158 |
| <GO:0002762> | negative regulation of myeloid leukocyte differentiation                                                                                         |  1.4383957 | 0.2152543 |
| <GO:0006721> | terpenoid metabolic process                                                                                                                      |  1.4347964 | 0.2152543 |
| <GO:0016101> | diterpenoid metabolic process                                                                                                                    |  1.4347964 | 0.2152543 |
| <GO:1904646> | cellular response to amyloid-beta                                                                                                                |  1.4277765 | 0.2152543 |
| <GO:0002220> | innate immune response activating cell surface receptor signaling pathway                                                                        |  1.3548272 | 0.2152543 |
| <GO:0050801> | monoatomic ion homeostasis                                                                                                                       |  1.2339923 | 0.2152543 |
| <GO:0043124> | negative regulation of canonical NF-kappaB signal transduction                                                                                   |  1.4008037 | 0.2156934 |
| <GO:1903008> | organelle disassembly                                                                                                                            | -1.6384310 | 0.2189068 |
| <GO:1903034> | regulation of response to wounding                                                                                                               |  1.3890002 | 0.2207505 |
| <GO:0006826> | iron ion transport                                                                                                                               |  1.4544293 | 0.2209746 |
| <GO:0097006> | regulation of plasma lipoprotein particle levels                                                                                                 |  1.4087565 | 0.2224351 |
| <GO:1902533> | positive regulation of intracellular signal transduction                                                                                         |  1.1686474 | 0.2224351 |
| <GO:0009636> | response to toxic substance                                                                                                                      |  1.3389996 | 0.2224451 |
| <GO:0022408> | negative regulation of cell-cell adhesion                                                                                                        |  1.3363908 | 0.2260660 |
| <GO:0002440> | production of molecular mediator of immune response                                                                                              |  1.3195398 | 0.2260660 |
| <GO:0050877> | nervous system process                                                                                                                           |  1.2223328 | 0.2260660 |
| <GO:0032495> | response to muramyl dipeptide                                                                                                                    |  1.4309302 | 0.2264691 |
| <GO:1905146> | lysosomal protein catabolic process                                                                                                              |  1.4516774 | 0.2269586 |
| <GO:0071526> | semaphorin-plexin signaling pathway                                                                                                              |  1.4514607 | 0.2269586 |
| <GO:0048143> | astrocyte activation                                                                                                                             |  1.4505812 | 0.2269586 |
| <GO:0072348> | sulfur compound transport                                                                                                                        |  1.4530087 | 0.2274042 |
| <GO:0072009> | nephron epithelium development                                                                                                                   |  1.4481559 | 0.2274042 |
| <GO:0002886> | regulation of myeloid leukocyte mediated immunity                                                                                                |  1.3870576 | 0.2290195 |
| <GO:0016064> | immunoglobulin mediated immune response                                                                                                          |  1.3772900 | 0.2292965 |
| <GO:0019724> | B cell mediated immunity                                                                                                                         |  1.3772900 | 0.2292965 |
| <GO:0019229> | regulation of vasoconstriction                                                                                                                   |  1.4470127 | 0.2310448 |
| <GO:1903530> | regulation of secretion by cell                                                                                                                  |  1.2365109 | 0.2311681 |
| <GO:0007632> | visual behavior                                                                                                                                  | -1.5726529 | 0.2327568 |
| <GO:0008542> | visual learning                                                                                                                                  | -1.5726529 | 0.2327568 |
| <GO:0140352> | export from cell                                                                                                                                 |  1.1893684 | 0.2348540 |
| <GO:0042440> | pigment metabolic process                                                                                                                        |  1.4440189 | 0.2390589 |
| <GO:0002703> | regulation of leukocyte mediated immunity                                                                                                        |  1.2829646 | 0.2390589 |
| <GO:0043954> | cellular component maintenance                                                                                                                   |  1.4175356 | 0.2405778 |
| <GO:0008203> | cholesterol metabolic process                                                                                                                    |  1.3905573 | 0.2419177 |
| <GO:0009063> | amino acid catabolic process                                                                                                                     |  1.4385238 | 0.2421037 |
| <GO:0021782> | glial cell development                                                                                                                           |  1.3668284 | 0.2421678 |
| <GO:0034614> | cellular response to reactive oxygen species                                                                                                     |  1.3656721 | 0.2421678 |
| <GO:0002752> | cell surface pattern recognition receptor signaling pathway                                                                                      |  1.3571663 | 0.2421678 |
| <GO:0141087> | positive regulation of inflammasome-mediated signaling pathway                                                                                   |  1.4213073 | 0.2444361 |
| <GO:1900227> | positive regulation of NLRP3 inflammasome complex assembly                                                                                       |  1.4213073 | 0.2444361 |
| <GO:0045604> | regulation of epidermal cell differentiation                                                                                                     |  1.3953024 | 0.2470173 |
| <GO:0120254> | olefinic compound metabolic process                                                                                                              |  1.3799912 | 0.2470173 |
| <GO:0001676> | long-chain fatty acid metabolic process                                                                                                          |  1.3846793 | 0.2478131 |
| <GO:0007596> | blood coagulation                                                                                                                                |  1.3219095 | 0.2487681 |
| <GO:0019722> | calcium-mediated signaling                                                                                                                       |  1.3172692 | 0.2527151 |
| <GO:0039531> | regulation of cytoplasmic pattern recognition receptor signaling pathway                                                                         |  1.3171008 | 0.2527151 |
| <GO:0055123> | digestive system development                                                                                                                     |  1.4170145 | 0.2537514 |
| <GO:0061041> | regulation of wound healing                                                                                                                      |  1.3937193 | 0.2561565 |
| <GO:0008645> | hexose transmembrane transport                                                                                                                   |  1.3687582 | 0.2561565 |
| <GO:0015749> | monosaccharide transmembrane transport                                                                                                           |  1.3687582 | 0.2561565 |
| <GO:1904659> | D-glucose transmembrane transport                                                                                                                |  1.3687582 | 0.2561565 |
| <GO:0032722> | positive regulation of chemokine production                                                                                                      |  1.3650084 | 0.2561565 |
| <GO:0002720> | positive regulation of cytokine production involved in immune response                                                                           |  1.3631639 | 0.2561565 |
| <GO:0120161> | regulation of cold-induced thermogenesis                                                                                                         |  1.3573744 | 0.2561565 |
| <GO:1904645> | response to amyloid-beta                                                                                                                         |  1.4116924 | 0.2571915 |
| <GO:0044057> | regulation of system process                                                                                                                     |  1.2629747 | 0.2571915 |
| <GO:0033194> | response to hydroperoxide                                                                                                                        |  1.4289497 | 0.2574341 |
| <GO:0015804> | neutral amino acid transport                                                                                                                     |  1.4140918 | 0.2574341 |
| <GO:0006633> | fatty acid biosynthetic process                                                                                                                  |  1.3735520 | 0.2574341 |
| <GO:0048678> | response to axon injury                                                                                                                          |  1.4123352 | 0.2615723 |
| <GO:0140467> | integrated stress response signaling                                                                                                             |  1.3821067 | 0.2649083 |
| <GO:0002183> | cytoplasmic translational initiation                                                                                                             | -1.5206992 | 0.2668669 |
| <GO:0015909> | long-chain fatty acid transport                                                                                                                  |  1.4063459 | 0.2668669 |
| <GO:0045815> | transcription initiation-coupled chromatin remodeling                                                                                            |  1.4023664 | 0.2668669 |
| <GO:0042886> | amide transport                                                                                                                                  |  1.2947316 | 0.2674149 |
| <GO:0048662> | negative regulation of smooth muscle cell proliferation                                                                                          |  1.4236593 | 0.2684304 |
| <GO:0045620> | negative regulation of lymphocyte differentiation                                                                                                |  1.4215711 | 0.2684304 |
| <GO:0000768> | syncytium formation by plasma membrane fusion                                                                                                    |  1.4189308 | 0.2684304 |
| <GO:0006949> | syncytium formation                                                                                                                              |  1.4189308 | 0.2684304 |
| <GO:0140253> | cell-cell fusion                                                                                                                                 |  1.4189308 | 0.2684304 |
| <GO:0031638> | zymogen activation                                                                                                                               |  1.4055199 | 0.2684304 |
| <GO:1900017> | positive regulation of cytokine production involved in inflammatory response                                                                     |  1.3999656 | 0.2684304 |
| <GO:0001508> | action potential                                                                                                                                 |  1.3977529 | 0.2684304 |
| <GO:0046545> | development of primary female sexual characteristics                                                                                             |  1.3970317 | 0.2684304 |
| <GO:0034143> | regulation of toll-like receptor 4 signaling pathway                                                                                             |  1.3951395 | 0.2684304 |
| <GO:0071622> | regulation of granulocyte chemotaxis                                                                                                             |  1.3949541 | 0.2684304 |
| <GO:1902652> | secondary alcohol metabolic process                                                                                                              |  1.3582754 | 0.2684304 |
| <GO:0006639> | acylglycerol metabolic process                                                                                                                   |  1.3450049 | 0.2684304 |
| <GO:0051092> | positive regulation of NF-kappaB transcription factor activity                                                                                   |  1.3420982 | 0.2684304 |
| <GO:0051047> | positive regulation of secretion                                                                                                                 |  1.2651167 | 0.2684304 |
| <GO:0045785> | positive regulation of cell adhesion                                                                                                             |  1.2167550 | 0.2684304 |
| <GO:0055080> | monoatomic cation homeostasis                                                                                                                    |  1.2093006 | 0.2684304 |
| <GO:0046889> | positive regulation of lipid biosynthetic process                                                                                                |  1.3911820 | 0.2696224 |
| <GO:0006879> | intracellular iron ion homeostasis                                                                                                               |  1.4013841 | 0.2702653 |
| <GO:0044281> | small molecule metabolic process                                                                                                                 |  1.1408713 | 0.2715777 |
| <GO:0050921> | positive regulation of chemotaxis                                                                                                                |  1.3371763 | 0.2722012 |
| <GO:0045682> | regulation of epidermis development                                                                                                              |  1.3880653 | 0.2744813 |
| <GO:0032527> | protein exit from endoplasmic reticulum                                                                                                          | -1.3793236 | 0.2744813 |
| <GO:1990845> | adaptive thermogenesis                                                                                                                           |  1.3084763 | 0.2744813 |
| <GO:0042113> | B cell activation                                                                                                                                |  1.2543985 | 0.2744813 |
| <GO:0010817> | regulation of hormone levels                                                                                                                     |  1.2411421 | 0.2744813 |
| <GO:0000041> | transition metal ion transport                                                                                                                   |  1.3698374 | 0.2772656 |
| <GO:0050678> | regulation of epithelial cell proliferation                                                                                                      |  1.3003083 | 0.2773923 |
| <GO:0009791> | post-embryonic development                                                                                                                       |  1.3793794 | 0.2783083 |
| <GO:0048705> | skeletal system morphogenesis                                                                                                                    |  1.3521430 | 0.2789715 |
| <GO:0030641> | regulation of cellular pH                                                                                                                        |  1.3277416 | 0.2802551 |
| <GO:0050817> | coagulation                                                                                                                                      |  1.2964559 | 0.2802551 |
| <GO:0031123> | RNA 3’-end processing                                                                                                                            | -1.4520502 | 0.2843647 |
| <GO:0006865> | amino acid transport                                                                                                                             |  1.3412837 | 0.2843647 |
| <GO:0032945> | negative regulation of mononuclear cell proliferation                                                                                            |  1.3527660 | 0.2858556 |
| <GO:0070509> | calcium ion import                                                                                                                               |  1.4083651 | 0.2867910 |
| <GO:0032940> | secretion by cell                                                                                                                                |  1.1750265 | 0.2902854 |
| <GO:0070252> | actin-mediated cell contraction                                                                                                                  |  1.3794406 | 0.2906917 |
| <GO:0043473> | pigmentation                                                                                                                                     |  1.3312048 | 0.2906917 |
| <GO:0001503> | ossification                                                                                                                                     |  1.2412417 | 0.2932200 |
| <GO:0003007> | heart morphogenesis                                                                                                                              |  1.3236214 | 0.2936153 |
| <GO:0009395> | phospholipid catabolic process                                                                                                                   |  1.3783884 | 0.2938439 |
| <GO:0140632> | canonical inflammasome complex assembly                                                                                                          |  1.3330628 | 0.2942303 |
| <GO:0141084> | inflammasome-mediated signaling pathway                                                                                                          |  1.3330628 | 0.2942303 |
| <GO:0099024> | plasma membrane invagination                                                                                                                     |  1.3719802 | 0.2961927 |
| <GO:0046460> | neutral lipid biosynthetic process                                                                                                               |  1.3716124 | 0.2961927 |
| <GO:0046463> | acylglycerol biosynthetic process                                                                                                                |  1.3716124 | 0.2961927 |
| <GO:0072330> | monocarboxylic acid biosynthetic process                                                                                                         |  1.3233587 | 0.2961927 |
| <GO:0006836> | neurotransmitter transport                                                                                                                       |  1.3265310 | 0.2984337 |
| <GO:0097242> | amyloid-beta clearance                                                                                                                           |  1.3758575 | 0.3014001 |
| <GO:2000191> | regulation of fatty acid transport                                                                                                               |  1.3742680 | 0.3014001 |
| <GO:0007566> | embryo implantation                                                                                                                              |  1.3741671 | 0.3014001 |
| <GO:0032481> | positive regulation of type I interferon production                                                                                              |  1.3189351 | 0.3014001 |
| <GO:0036336> | dendritic cell migration                                                                                                                         |  1.3923643 | 0.3031368 |
| <GO:0070741> | response to interleukin-6                                                                                                                        |  1.3909513 | 0.3031368 |
| <GO:0071354> | cellular response to interleukin-6                                                                                                               |  1.3909513 | 0.3031368 |
| <GO:1902692> | regulation of neuroblast proliferation                                                                                                           |  1.3878799 | 0.3031368 |
| <GO:0046849> | bone remodeling                                                                                                                                  |  1.3439241 | 0.3031368 |
| <GO:0070664> | negative regulation of leukocyte proliferation                                                                                                   |  1.3389221 | 0.3031368 |
| <GO:0045598> | regulation of fat cell differentiation                                                                                                           |  1.3257868 | 0.3031368 |
| <GO:0031348> | negative regulation of defense response                                                                                                          |  1.2541568 | 0.3031368 |
| <GO:0022409> | positive regulation of cell-cell adhesion                                                                                                        |  1.2418063 | 0.3031368 |
| <GO:2000628> | regulation of miRNA metabolic process                                                                                                            |  1.3478497 | 0.3033972 |
| <GO:1903532> | positive regulation of secretion by cell                                                                                                         |  1.2495624 | 0.3033972 |
| <GO:0019432> | triglyceride biosynthetic process                                                                                                                |  1.3650441 | 0.3039894 |
| <GO:0007267> | cell-cell signaling                                                                                                                              |  1.1690622 | 0.3062937 |
| <GO:0007269> | neurotransmitter secretion                                                                                                                       |  1.3280195 | 0.3068572 |
| <GO:0099643> | signal release from synapse                                                                                                                      |  1.3280195 | 0.3068572 |
| <GO:0006638> | neutral lipid metabolic process                                                                                                                  |  1.3218582 | 0.3099469 |
| <GO:0032740> | positive regulation of interleukin-17 production                                                                                                 |  1.3849989 | 0.3126728 |
| <GO:0002027> | regulation of heart rate                                                                                                                         |  1.3811983 | 0.3148826 |
| <GO:0150146> | cell junction disassembly                                                                                                                        |  1.3794963 | 0.3148826 |
| <GO:0048638> | regulation of developmental growth                                                                                                               |  1.2723102 | 0.3148826 |
| <GO:0006873> | intracellular monoatomic ion homeostasis                                                                                                         |  1.1982738 | 0.3148826 |
| <GO:0032733> | positive regulation of interleukin-10 production                                                                                                 |  1.3673492 | 0.3159056 |
| <GO:0051149> | positive regulation of muscle cell differentiation                                                                                               |  1.3665798 | 0.3159056 |
| <GO:0032835> | glomerulus development                                                                                                                           |  1.3591291 | 0.3159056 |
| <GO:0048736> | appendage development                                                                                                                            |  1.3339796 | 0.3159056 |
| <GO:0060173> | limb development                                                                                                                                 |  1.3339796 | 0.3159056 |
| <GO:0042063> | gliogenesis                                                                                                                                      |  1.2453570 | 0.3160535 |
| <GO:2000630> | positive regulation of miRNA metabolic process                                                                                                   |  1.3788174 | 0.3192641 |
| <GO:0007498> | mesoderm development                                                                                                                             |  1.3459401 | 0.3214986 |
| <GO:0030217> | T cell differentiation                                                                                                                           |  1.2239313 | 0.3227592 |
| <GO:0002686> | negative regulation of leukocyte migration                                                                                                       |  1.3776364 | 0.3256615 |
| <GO:0002685> | regulation of leukocyte migration                                                                                                                |  1.2730630 | 0.3256615 |
| <GO:0032787> | monocarboxylic acid metabolic process                                                                                                            |  1.1956901 | 0.3256615 |
| <GO:0140962> | multicellular organismal-level chemical homeostasis                                                                                              |  1.3627644 | 0.3271747 |
| <GO:0016485> | protein processing                                                                                                                               |  1.2677185 | 0.3283925 |
| <GO:0035107> | appendage morphogenesis                                                                                                                          |  1.3346742 | 0.3289326 |
| <GO:0035108> | limb morphogenesis                                                                                                                               |  1.3346742 | 0.3289326 |
| <GO:0046330> | positive regulation of JNK cascade                                                                                                               |  1.3228261 | 0.3289326 |
| <GO:0002690> | positive regulation of leukocyte chemotaxis                                                                                                      |  1.3115753 | 0.3289326 |
| <GO:0008202> | steroid metabolic process                                                                                                                        |  1.2779059 | 0.3302466 |
| <GO:0034728> | nucleosome organization                                                                                                                          |  1.3454045 | 0.3306837 |
| <GO:0045646> | regulation of erythrocyte differentiation                                                                                                        |  1.3504012 | 0.3325451 |
| <GO:0042981> | regulation of apoptotic process                                                                                                                  |  1.1275485 | 0.3364272 |
| <GO:0090322> | regulation of superoxide metabolic process                                                                                                       |  1.3570022 | 0.3387493 |
| <GO:0060840> | artery development                                                                                                                               |  1.3513947 | 0.3395791 |
| <GO:0006939> | smooth muscle contraction                                                                                                                        |  1.3464939 | 0.3395791 |
| <GO:0008585> | female gonad development                                                                                                                         |  1.3416860 | 0.3395791 |
| <GO:0010876> | lipid localization                                                                                                                               |  1.2120826 | 0.3418476 |
| <GO:0018108> | peptidyl-tyrosine phosphorylation                                                                                                                |  1.3075091 | 0.3425336 |
| <GO:0018212> | peptidyl-tyrosine modification                                                                                                                   |  1.3075091 | 0.3425336 |
| <GO:0048754> | branching morphogenesis of an epithelial tube                                                                                                    |  1.3535836 | 0.3433752 |
| <GO:0060419> | heart growth                                                                                                                                     |  1.3459441 | 0.3446474 |
| <GO:0051260> | protein homooligomerization                                                                                                                      |  1.2913047 | 0.3446474 |
| <GO:0032728> | positive regulation of interferon-beta production                                                                                                |  1.3471416 | 0.3456736 |
| <GO:0030879> | mammary gland development                                                                                                                        |  1.3026230 | 0.3456736 |
| <GO:0007189> | adenylate cyclase-activating G protein-coupled receptor signaling pathway                                                                        |  1.3406915 | 0.3457279 |
| <GO:0007039> | protein catabolic process in the vacuole                                                                                                         |  1.3668954 | 0.3463526 |
| <GO:0030279> | negative regulation of ossification                                                                                                              |  1.3656229 | 0.3463526 |
| <GO:0051385> | response to mineralocorticoid                                                                                                                    |  1.3643362 | 0.3463526 |
| <GO:0042475> | odontogenesis of dentin-containing tooth                                                                                                         |  1.3453377 | 0.3468773 |
| <GO:1902622> | regulation of neutrophil migration                                                                                                               |  1.3394035 | 0.3483040 |
| <GO:0032735> | positive regulation of interleukin-12 production                                                                                                 |  1.3355113 | 0.3483040 |
| <GO:0042220> | response to cocaine                                                                                                                              | -1.2929941 | 0.3486450 |
| <GO:0032753> | positive regulation of interleukin-4 production                                                                                                  |  1.3571536 | 0.3496423 |
| <GO:0019216> | regulation of lipid metabolic process                                                                                                            |  1.2566682 | 0.3496423 |
| <GO:0051147> | regulation of muscle cell differentiation                                                                                                        |  1.3118580 | 0.3504197 |
| <GO:0032411> | positive regulation of transporter activity                                                                                                      |  1.3511075 | 0.3519994 |
| <GO:0014855> | striated muscle cell proliferation                                                                                                               |  1.3501970 | 0.3519994 |
| <GO:0008037> | cell recognition                                                                                                                                 |  1.3035136 | 0.3519994 |
| <GO:1903522> | regulation of blood circulation                                                                                                                  |  1.2848177 | 0.3519994 |
| <GO:0046631> | alpha-beta T cell activation                                                                                                                     |  1.2618038 | 0.3566348 |
| <GO:0007159> | leukocyte cell-cell adhesion                                                                                                                     |  1.2091485 | 0.3566348 |
| <GO:0042060> | wound healing                                                                                                                                    |  1.2007919 | 0.3566348 |
| <GO:0002687> | positive regulation of leukocyte migration                                                                                                       |  1.2800865 | 0.3568323 |
| <GO:0007586> | digestion                                                                                                                                        |  1.3288043 | 0.3599789 |
| <GO:0022600> | digestive system process                                                                                                                         |  1.3288043 | 0.3599789 |
| <GO:0002709> | regulation of T cell mediated immunity                                                                                                           |  1.3287145 | 0.3599789 |
| <GO:0002861> | regulation of inflammatory response to antigenic stimulus                                                                                        |  1.3571698 | 0.3604894 |
| <GO:0046323> | D-glucose import                                                                                                                                 |  1.3293243 | 0.3604894 |
| <GO:0051384> | response to glucocorticoid                                                                                                                       |  1.2963687 | 0.3604894 |
| <GO:0043542> | endothelial cell migration                                                                                                                       |  1.2668823 | 0.3604894 |
| <GO:1901570> | fatty acid derivative biosynthetic process                                                                                                       |  1.3550532 | 0.3626866 |
| <GO:0006869> | lipid transport                                                                                                                                  |  1.2028258 | 0.3626866 |
| <GO:0031099> | regeneration                                                                                                                                     |  1.2747928 | 0.3630885 |
| <GO:0010543> | regulation of platelet activation                                                                                                                |  1.3310840 | 0.3676354 |
| <GO:0021799> | cerebral cortex radially oriented cell migration                                                                                                 |  1.3499937 | 0.3678048 |
| <GO:0050672> | negative regulation of lymphocyte proliferation                                                                                                  |  1.3157144 | 0.3697533 |
| <GO:0043300> | regulation of leukocyte degranulation                                                                                                            |  1.3078377 | 0.3794936 |
| <GO:1903039> | positive regulation of leukocyte cell-cell adhesion                                                                                              |  1.2175384 | 0.3819014 |
| <GO:0014037> | Schwann cell differentiation                                                                                                                     | -1.2560865 | 0.3836495 |
| <GO:0048661> | positive regulation of smooth muscle cell proliferation                                                                                          |  1.3168708 | 0.3859572 |
| <GO:0070613> | regulation of protein processing                                                                                                                 |  1.3261750 | 0.3860497 |
| <GO:1903317> | regulation of protein maturation                                                                                                                 |  1.3261750 | 0.3860497 |
| <GO:0071466> | cellular response to xenobiotic stimulus                                                                                                         |  1.2937683 | 0.3860497 |
| <GO:1901654> | response to ketone                                                                                                                               |  1.2491601 | 0.3868813 |
| <GO:0050868> | negative regulation of T cell activation                                                                                                         |  1.2778933 | 0.3869596 |
| <GO:0042445> | hormone metabolic process                                                                                                                        |  1.2756800 | 0.3869596 |
| <GO:0043603> | amide metabolic process                                                                                                                          |  1.1751718 | 0.3881852 |
| <GO:1903409> | reactive oxygen species biosynthetic process                                                                                                     |  1.3153646 | 0.3892021 |
| <GO:0001912> | positive regulation of leukocyte mediated cytotoxicity                                                                                           |  1.3338025 | 0.3906347 |
| <GO:0031343> | positive regulation of cell killing                                                                                                              |  1.3338025 | 0.3906347 |
| <GO:0007042> | lysosomal lumen acidification                                                                                                                    |  1.3261800 | 0.3906347 |
| <GO:0016192> | vesicle-mediated transport                                                                                                                       |  1.1017506 | 0.3906347 |
| <GO:0051091> | positive regulation of DNA-binding transcription factor activity                                                                                 |  1.2891499 | 0.3909645 |
| <GO:0050679> | positive regulation of epithelial cell proliferation                                                                                             |  1.2644440 | 0.3965937 |
| <GO:0042632> | cholesterol homeostasis                                                                                                                          |  1.3176240 | 0.3969771 |
| <GO:0055092> | sterol homeostasis                                                                                                                               |  1.3176240 | 0.3969771 |
| <GO:0048732> | gland development                                                                                                                                |  1.2167782 | 0.3969771 |
| <GO:0002449> | lymphocyte mediated immunity                                                                                                                     |  1.2160905 | 0.3969771 |
| <GO:0070542> | response to fatty acid                                                                                                                           |  1.3154820 | 0.3996569 |
| <GO:0048565> | digestive tract development                                                                                                                      |  1.3137133 | 0.3996569 |
| <GO:0045580> | regulation of T cell differentiation                                                                                                             |  1.2513064 | 0.3996589 |
| <GO:0019218> | regulation of steroid metabolic process                                                                                                          |  1.3125508 | 0.4023219 |
| <GO:0046660> | female sex differentiation                                                                                                                       |  1.3208824 | 0.4034339 |
| <GO:1903305> | regulation of regulated secretory pathway                                                                                                        |  1.2725765 | 0.4056638 |
| <GO:0032814> | regulation of natural killer cell activation                                                                                                     |  1.3255133 | 0.4059623 |
| <GO:0009581> | detection of external stimulus                                                                                                                   |  1.3253601 | 0.4059623 |
| <GO:0009582> | detection of abiotic stimulus                                                                                                                    |  1.3253601 | 0.4059623 |
| <GO:0033559> | unsaturated fatty acid metabolic process                                                                                                         |  1.3063093 | 0.4059623 |
| <GO:0070979> | protein K11-linked ubiquitination                                                                                                                | -1.4024875 | 0.4071270 |
| <GO:0046636> | negative regulation of alpha-beta T cell activation                                                                                              |  1.3325371 | 0.4120268 |
| <GO:0006694> | steroid biosynthetic process                                                                                                                     |  1.2834816 | 0.4120268 |
| <GO:0045639> | positive regulation of myeloid cell differentiation                                                                                              |  1.2692659 | 0.4120268 |
| <GO:0055094> | response to lipoprotein particle                                                                                                                 |  1.3263126 | 0.4126909 |
| <GO:0071404> | cellular response to low-density lipoprotein particle stimulus                                                                                   |  1.3263126 | 0.4126909 |
| <GO:0045619> | regulation of lymphocyte differentiation                                                                                                         |  1.2409901 | 0.4126909 |
| <GO:0043010> | camera-type eye development                                                                                                                      |  1.2405730 | 0.4129613 |
| <GO:0015695> | organic cation transport                                                                                                                         |  1.2931789 | 0.4129858 |
| <GO:0046683> | response to organophosphorus                                                                                                                     |  1.2747931 | 0.4134662 |
| <GO:0042110> | T cell activation                                                                                                                                |  1.1523713 | 0.4134662 |
| <GO:0043067> | regulation of programmed cell death                                                                                                              |  1.1095007 | 0.4134662 |
| <GO:0051588> | regulation of neurotransmitter transport                                                                                                         |  1.3130410 | 0.4145061 |
| <GO:0010827> | regulation of D-glucose transmembrane transport                                                                                                  |  1.3060010 | 0.4145061 |
| <GO:1903038> | negative regulation of leukocyte cell-cell adhesion                                                                                              |  1.2720658 | 0.4145061 |
| <GO:0060420> | regulation of heart growth                                                                                                                       |  1.3221652 | 0.4146454 |
| <GO:0045621> | positive regulation of lymphocyte differentiation                                                                                                |  1.2586469 | 0.4146454 |
| <GO:0046620> | regulation of organ growth                                                                                                                       |  1.2937679 | 0.4164160 |
| <GO:0007423> | sensory organ development                                                                                                                        |  1.2047608 | 0.4166350 |
| <GO:0060562> | epithelial tube morphogenesis                                                                                                                    |  1.2356975 | 0.4169169 |
| <GO:0050777> | negative regulation of immune response                                                                                                           |  1.2382779 | 0.4190262 |
| <GO:0030307> | positive regulation of cell growth                                                                                                               |  1.2591497 | 0.4201523 |
| <GO:0050000> | chromosome localization                                                                                                                          | -1.2858302 | 0.4215117 |
| <GO:0032633> | interleukin-4 production                                                                                                                         |  1.3274496 | 0.4222219 |
| <GO:0032673> | regulation of interleukin-4 production                                                                                                           |  1.3274496 | 0.4222219 |
| <GO:0032371> | regulation of sterol transport                                                                                                                   |  1.3182238 | 0.4222219 |
| <GO:0032374> | regulation of cholesterol transport                                                                                                              |  1.3182238 | 0.4222219 |
| <GO:0033173> | calcineurin-NFAT signaling cascade                                                                                                               |  1.3178277 | 0.4222219 |
| <GO:0019674> | NAD+ metabolic process                                                                                                                           |  1.3167643 | 0.4222219 |
| <GO:0141137> | positive regulation of gene expression, epigenetic                                                                                               |  1.3064120 | 0.4222219 |
| <GO:0010934> | macrophage cytokine production                                                                                                                   |  1.3063314 | 0.4222219 |
| <GO:0010935> | regulation of macrophage cytokine production                                                                                                     |  1.3063314 | 0.4222219 |
| <GO:0009743> | response to carbohydrate                                                                                                                         |  1.2502956 | 0.4222219 |
| <GO:0032924> | activin receptor signaling pathway                                                                                                               | -1.2388630 | 0.4222219 |
| <GO:0030003> | intracellular monoatomic cation homeostasis                                                                                                      |  1.1705680 | 0.4222219 |
| <GO:0043372> | positive regulation of CD4-positive, alpha-beta T cell differentiation                                                                           |  1.3076442 | 0.4244534 |
| <GO:2000278> | regulation of DNA biosynthetic process                                                                                                           |  1.3126624 | 0.4256049 |
| <GO:1903037> | regulation of leukocyte cell-cell adhesion                                                                                                       |  1.1839950 | 0.4268987 |
| <GO:0001885> | endothelial cell development                                                                                                                     |  1.2879165 | 0.4280688 |
| <GO:0150063> | visual system development                                                                                                                        |  1.2299589 | 0.4290811 |
| <GO:0030336> | negative regulation of cell migration                                                                                                            |  1.2131302 | 0.4290811 |
| <GO:0003333> | amino acid transmembrane transport                                                                                                               |  1.2895659 | 0.4310485 |
| <GO:0050790> | regulation of catalytic activity                                                                                                                 |  1.1545971 | 0.4310485 |
| <GO:0002064> | epithelial cell development                                                                                                                      |  1.2305781 | 0.4312543 |
| <GO:0061614> | miRNA transcription                                                                                                                              |  1.2962753 | 0.4342443 |
| <GO:1902893> | regulation of miRNA transcription                                                                                                                |  1.2962753 | 0.4342443 |
| <GO:0006641> | triglyceride metabolic process                                                                                                                   |  1.2613135 | 0.4347254 |
| <GO:0007405> | neuroblast proliferation                                                                                                                         |  1.3014243 | 0.4360734 |
| <GO:0051153> | regulation of striated muscle cell differentiation                                                                                               |  1.2997938 | 0.4360734 |
| <GO:0070302> | regulation of stress-activated protein kinase signaling cascade                                                                                  |  1.2946852 | 0.4360734 |
| <GO:0040013> | negative regulation of locomotion                                                                                                                |  1.2044608 | 0.4360734 |
| <GO:2000146> | negative regulation of cell motility                                                                                                             |  1.2044608 | 0.4360734 |
| <GO:0006979> | response to oxidative stress                                                                                                                     |  1.1756106 | 0.4360734 |
| <GO:0010634> | positive regulation of epithelial cell migration                                                                                                 |  1.2938783 | 0.4361758 |
| <GO:0060348> | bone development                                                                                                                                 |  1.2563971 | 0.4361758 |
| <GO:0002712> | regulation of B cell mediated immunity                                                                                                           |  1.3149812 | 0.4362221 |
| <GO:0002889> | regulation of immunoglobulin mediated immune response                                                                                            |  1.3149812 | 0.4362221 |
| <GO:0031016> | pancreas development                                                                                                                             |  1.3078269 | 0.4362221 |
| <GO:0002888> | positive regulation of myeloid leukocyte mediated immunity                                                                                       |  1.3072854 | 0.4362221 |
| <GO:0043299> | leukocyte degranulation                                                                                                                          |  1.2540690 | 0.4362221 |
| <GO:2001237> | negative regulation of extrinsic apoptotic signaling pathway                                                                                     |  1.2538573 | 0.4362221 |
| <GO:0002711> | positive regulation of T cell mediated immunity                                                                                                  |  1.3021998 | 0.4366997 |
| <GO:0050870> | positive regulation of T cell activation                                                                                                         |  1.1991875 | 0.4366997 |
| <GO:1902895> | positive regulation of miRNA transcription                                                                                                       |  1.3049233 | 0.4389746 |
| <GO:0006691> | leukotriene metabolic process                                                                                                                    |  1.3040666 | 0.4389746 |
| <GO:2001257> | regulation of cation channel activity                                                                                                            |  1.3000485 | 0.4389746 |
| <GO:0002090> | regulation of receptor internalization                                                                                                           |  1.2768111 | 0.4389746 |
| <GO:0001776> | leukocyte homeostasis                                                                                                                            |  1.2446730 | 0.4389746 |
| <GO:0042129> | regulation of T cell proliferation                                                                                                               |  1.2195644 | 0.4389746 |
| <GO:0034381> | plasma lipoprotein particle clearance                                                                                                            |  1.2966629 | 0.4419698 |
| <GO:0034383> | low-density lipoprotein particle clearance                                                                                                       |  1.2966629 | 0.4419698 |
| <GO:0034110> | regulation of homotypic cell-cell adhesion                                                                                                       |  1.2963182 | 0.4419698 |
| <GO:0032088> | negative regulation of NF-kappaB transcription factor activity                                                                                   |  1.3121413 | 0.4428353 |
| <GO:0048844> | artery morphogenesis                                                                                                                             |  1.3023327 | 0.4467102 |
| <GO:0008630> | intrinsic apoptotic signaling pathway in response to DNA damage                                                                                  |  1.2712175 | 0.4467102 |
| <GO:0031124> | mRNA 3’-end processing                                                                                                                           | -1.2179487 | 0.4483928 |
| <GO:0046427> | positive regulation of receptor signaling pathway via JAK-STAT                                                                                   |  1.3037989 | 0.4489603 |
| <GO:1904894> | positive regulation of receptor signaling pathway via STAT                                                                                       |  1.3037989 | 0.4489603 |
| <GO:0006720> | isoprenoid metabolic process                                                                                                                     |  1.2875635 | 0.4489603 |
| <GO:0003206> | cardiac chamber morphogenesis                                                                                                                    |  1.2844103 | 0.4489603 |
| <GO:0046324> | regulation of D-glucose import                                                                                                                   |  1.2815192 | 0.4489603 |
| <GO:0071402> | cellular response to lipoprotein particle stimulus                                                                                               |  1.2796265 | 0.4489603 |
| <GO:0043370> | regulation of CD4-positive, alpha-beta T cell differentiation                                                                                    |  1.2698025 | 0.4489603 |
| <GO:0098742> | cell-cell adhesion via plasma-membrane adhesion molecules                                                                                        |  1.2661280 | 0.4489603 |
| <GO:0042102> | positive regulation of T cell proliferation                                                                                                      |  1.2486264 | 0.4489603 |
| <GO:0008277> | regulation of G protein-coupled receptor signaling pathway                                                                                       |  1.2370046 | 0.4489603 |
| <GO:0035710> | CD4-positive, alpha-beta T cell activation                                                                                                       |  1.2292472 | 0.4489603 |
| <GO:0050769> | positive regulation of neurogenesis                                                                                                              |  1.2275600 | 0.4489603 |
| <GO:0032479> | regulation of type I interferon production                                                                                                       |  1.2268951 | 0.4489603 |
| <GO:0032606> | type I interferon production                                                                                                                     |  1.2268951 | 0.4489603 |
| <GO:2001236> | regulation of extrinsic apoptotic signaling pathway                                                                                              |  1.2202024 | 0.4489603 |
| <GO:0010001> | glial cell differentiation                                                                                                                       |  1.2166964 | 0.4489603 |
| <GO:0051051> | negative regulation of transport                                                                                                                 |  1.1631438 | 0.4489603 |
| <GO:0007507> | heart development                                                                                                                                |  1.1600003 | 0.4489603 |
| <GO:0019725> | cellular homeostasis                                                                                                                             |  1.1243345 | 0.4489603 |
| <GO:0007399> | nervous system development                                                                                                                       |  1.0885937 | 0.4511498 |
| <GO:0070098> | chemokine-mediated signaling pathway                                                                                                             |  1.2889285 | 0.4520878 |
| <GO:0043405> | regulation of MAP kinase activity                                                                                                                |  1.2886673 | 0.4520878 |
| <GO:0060907> | positive regulation of macrophage cytokine production                                                                                            |  1.2838340 | 0.4520878 |
| <GO:0001570> | vasculogenesis                                                                                                                                   |  1.2819076 | 0.4520878 |
| <GO:0046637> | regulation of alpha-beta T cell differentiation                                                                                                  |  1.2420504 | 0.4520878 |
| <GO:0050654> | chondroitin sulfate proteoglycan metabolic process                                                                                               |  1.2947618 | 0.4525433 |
| <GO:0033363> | secretory granule organization                                                                                                                   |  1.2741120 | 0.4535717 |
| <GO:0001909> | leukocyte mediated cytotoxicity                                                                                                                  |  1.2268093 | 0.4535717 |
| <GO:0048880> | sensory system development                                                                                                                       |  1.2174459 | 0.4535717 |
| <GO:0051156> | glucose 6-phosphate metabolic process                                                                                                            |  1.2790994 | 0.4541790 |
| <GO:0001961> | positive regulation of cytokine-mediated signaling pathway                                                                                       |  1.2771409 | 0.4543407 |
| <GO:0048873> | homeostasis of number of cells within a tissue                                                                                                   |  1.2730762 | 0.4549315 |
| <GO:0014074> | response to purine-containing compound                                                                                                           |  1.2465369 | 0.4553320 |
| <GO:0140895> | cell surface toll-like receptor signaling pathway                                                                                                |  1.2457119 | 0.4553320 |
| <GO:0008217> | regulation of blood pressure                                                                                                                     |  1.2313676 | 0.4570132 |
| <GO:0007420> | brain development                                                                                                                                |  1.1521127 | 0.4582234 |
| <GO:0048259> | regulation of receptor-mediated endocytosis                                                                                                      |  1.2240823 | 0.4585194 |
| <GO:0007565> | female pregnancy                                                                                                                                 |  1.2185743 | 0.4585194 |
| <GO:0010569> | regulation of double-strand break repair via homologous recombination                                                                            | -1.4011868 | 0.4601170 |
| <GO:0002313> | mature B cell differentiation involved in immune response                                                                                        |  1.2902071 | 0.4601170 |
| <GO:0048762> | mesenchymal cell differentiation                                                                                                                 |  1.2166089 | 0.4601170 |
| <GO:0001932> | regulation of protein phosphorylation                                                                                                            |  1.1706993 | 0.4601170 |
| <GO:0044093> | positive regulation of molecular function                                                                                                        |  1.1507367 | 0.4601170 |
| <GO:0007169> | cell surface receptor protein tyrosine kinase signaling pathway                                                                                  |  1.1409415 | 0.4601170 |
| <GO:0048639> | positive regulation of developmental growth                                                                                                      |  1.2495988 | 0.4611021 |
| <GO:0008015> | blood circulation                                                                                                                                |  1.1855108 | 0.4621253 |
| <GO:0050714> | positive regulation of protein secretion                                                                                                         |  1.2247552 | 0.4629748 |
| <GO:0032890> | regulation of organic acid transport                                                                                                             |  1.2672716 | 0.4654161 |
| <GO:0042098> | T cell proliferation                                                                                                                             |  1.2070422 | 0.4671220 |
| <GO:0046777> | protein autophosphorylation                                                                                                                      |  1.2205322 | 0.4675151 |
| <GO:0046634> | regulation of alpha-beta T cell activation                                                                                                       |  1.2127633 | 0.4705170 |
| <GO:0045453> | bone resorption                                                                                                                                  |  1.2661903 | 0.4713940 |
| <GO:1902235> | regulation of endoplasmic reticulum stress-induced intrinsic apoptotic signaling pathway                                                         |  1.2842048 | 0.4716311 |
| <GO:0046890> | regulation of lipid biosynthetic process                                                                                                         |  1.2266200 | 0.4719517 |
| <GO:0002065> | columnar/cuboidal epithelial cell differentiation                                                                                                |  1.2713750 | 0.4729697 |
| <GO:0050673> | epithelial cell proliferation                                                                                                                    |  1.1833177 | 0.4729697 |
| <GO:0030193> | regulation of blood coagulation                                                                                                                  |  1.2736260 | 0.4733488 |
| <GO:1900046> | regulation of hemostasis                                                                                                                         |  1.2736260 | 0.4733488 |
| <GO:0034502> | protein localization to chromosome                                                                                                               |  1.2692003 | 0.4733488 |
| <GO:1903307> | positive regulation of regulated secretory pathway                                                                                               |  1.2648139 | 0.4733488 |
| <GO:0072175> | epithelial tube formation                                                                                                                        |  1.2563023 | 0.4733488 |
| <GO:0009719> | response to endogenous stimulus                                                                                                                  |  1.0978647 | 0.4733488 |
| <GO:0035148> | tube formation                                                                                                                                   |  1.2461761 | 0.4760283 |
| <GO:0060322> | head development                                                                                                                                 |  1.1367498 | 0.4767698 |
| <GO:0043414> | macromolecule methylation                                                                                                                        | -1.3016511 | 0.4771137 |
| <GO:0045582> | positive regulation of T cell differentiation                                                                                                    |  1.2098422 | 0.4771137 |
| <GO:0008016> | regulation of heart contraction                                                                                                                  |  1.2458870 | 0.4773813 |
| <GO:0014033> | neural crest cell differentiation                                                                                                                |  1.2867344 | 0.4776289 |
| <GO:0070935> | 3’-UTR-mediated mRNA stabilization                                                                                                               |  1.2664068 | 0.4787575 |
| <GO:1990823> | response to leukemia inhibitory factor                                                                                                           |  1.2419707 | 0.4787879 |
| <GO:1990868> | response to chemokine                                                                                                                            |  1.2608490 | 0.4805055 |
| <GO:1990869> | cellular response to chemokine                                                                                                                   |  1.2608490 | 0.4805055 |
| <GO:0097720> | calcineurin-mediated signaling                                                                                                                   |  1.2847951 | 0.4807621 |
| <GO:0006887> | exocytosis                                                                                                                                       |  1.1517699 | 0.4807621 |
| <GO:0007215> | glutamate receptor signaling pathway                                                                                                             |  1.2704893 | 0.4809794 |
| <GO:2000177> | regulation of neural precursor cell proliferation                                                                                                |  1.2587721 | 0.4809794 |
| <GO:0010906> | regulation of glucose metabolic process                                                                                                          |  1.2476474 | 0.4809794 |
| <GO:0045471> | response to ethanol                                                                                                                              |  1.2377087 | 0.4809794 |
| <GO:0032928> | regulation of superoxide anion generation                                                                                                        |  1.2801390 | 0.4827707 |
| <GO:0032930> | positive regulation of superoxide anion generation                                                                                               |  1.2801390 | 0.4827707 |
| <GO:0055017> | cardiac muscle tissue growth                                                                                                                     |  1.2682905 | 0.4827707 |
| <GO:0045601> | regulation of endothelial cell differentiation                                                                                                   |  1.2613964 | 0.4827707 |
| <GO:0016525> | negative regulation of angiogenesis                                                                                                              |  1.2542109 | 0.4846099 |
| <GO:1901343> | negative regulation of vasculature development                                                                                                   |  1.2542109 | 0.4846099 |
| <GO:2000181> | negative regulation of blood vessel morphogenesis                                                                                                |  1.2542109 | 0.4846099 |
| <GO:0017157> | regulation of exocytosis                                                                                                                         |  1.2002488 | 0.4846099 |
| <GO:0033627> | cell adhesion mediated by integrin                                                                                                               |  1.2656500 | 0.4848588 |
| <GO:0009409> | response to cold                                                                                                                                 |  1.2640805 | 0.4848588 |
| <GO:2000514> | regulation of CD4-positive, alpha-beta T cell activation                                                                                         |  1.2459588 | 0.4848588 |
| <GO:0043367> | CD4-positive, alpha-beta T cell differentiation                                                                                                  |  1.2238524 | 0.4848588 |
| <GO:0002456> | T cell mediated immunity                                                                                                                         |  1.2232390 | 0.4848588 |
| <GO:0062197> | cellular response to chemical stress                                                                                                             |  1.1628598 | 0.4848588 |
| <GO:0030098> | lymphocyte differentiation                                                                                                                       |  1.1405507 | 0.4848588 |
| <GO:0051216> | cartilage development                                                                                                                            |  1.2524699 | 0.4852321 |
| <GO:0046638> | positive regulation of alpha-beta T cell differentiation                                                                                         |  1.2482558 | 0.4852321 |
| <GO:0001654> | eye development                                                                                                                                  |  1.1959919 | 0.4895683 |
| <GO:0050730> | regulation of peptidyl-tyrosine phosphorylation                                                                                                  |  1.2539379 | 0.4901295 |
| <GO:0048872> | homeostasis of number of cells                                                                                                                   |  1.1446557 | 0.4909091 |
| <GO:0071320> | cellular response to cAMP                                                                                                                        |  1.2764053 | 0.4925044 |
| <GO:0008589> | regulation of smoothened signaling pathway                                                                                                       |  1.2751297 | 0.4925044 |
| <GO:0048536> | spleen development                                                                                                                               |  1.2745514 | 0.4925044 |
| <GO:0070301> | cellular response to hydrogen peroxide                                                                                                           |  1.2516557 | 0.4925044 |
| <GO:0007259> | cell surface receptor signaling pathway via JAK-STAT                                                                                             |  1.2229054 | 0.4925044 |
| <GO:0097696> | cell surface receptor signaling pathway via STAT                                                                                                 |  1.2229054 | 0.4925044 |
| <GO:0098771> | inorganic ion homeostasis                                                                                                                        |  1.1425770 | 0.4925044 |
| <GO:0060485> | mesenchyme development                                                                                                                           |  1.1998612 | 0.4928389 |
| <GO:0001934> | positive regulation of protein phosphorylation                                                                                                   |  1.1960770 | 0.4928389 |
| <GO:0002431> | Fc receptor mediated stimulatory signaling pathway                                                                                               |  1.2596127 | 0.4929571 |
| <GO:0044282> | small molecule catabolic process                                                                                                                 |  1.1851743 | 0.4941498 |
| <GO:0014032> | neural crest cell development                                                                                                                    |  1.2588997 | 0.4954461 |
| <GO:0098773> | skin epidermis development                                                                                                                       |  1.2609872 | 0.4973689 |
| <GO:1902475> | L-alpha-amino acid transmembrane transport                                                                                                       |  1.2537221 | 0.4973689 |
| <GO:0002312> | B cell activation involved in immune response                                                                                                    |  1.2467590 | 0.4973689 |
| <GO:0007589> | body fluid secretion                                                                                                                             |  1.2466131 | 0.4973689 |
| <GO:0032608> | interferon-beta production                                                                                                                       |  1.2334217 | 0.4973689 |
| <GO:0032648> | regulation of interferon-beta production                                                                                                         |  1.2334217 | 0.4973689 |
| <GO:0071763> | nuclear membrane organization                                                                                                                    | -1.2131595 | 0.4973689 |
| <GO:0023061> | signal release                                                                                                                                   |  1.1494106 | 0.4973689 |
| <GO:0035914> | skeletal muscle cell differentiation                                                                                                             |  1.2682920 | 0.4999439 |
| <GO:0002062> | chondrocyte differentiation                                                                                                                      |  1.2558901 | 0.5001143 |
| <GO:0048008> | platelet-derived growth factor receptor signaling pathway                                                                                        |  1.2539900 | 0.5001143 |
| <GO:0030149> | sphingolipid catabolic process                                                                                                                   |  1.2534068 | 0.5001143 |
| <GO:0050810> | regulation of steroid biosynthetic process                                                                                                       |  1.2486351 | 0.5001143 |
| <GO:0046165> | alcohol biosynthetic process                                                                                                                     |  1.2306667 | 0.5001143 |
| <GO:0016054> | organic acid catabolic process                                                                                                                   |  1.2152163 | 0.5001143 |
| <GO:0046395> | carboxylic acid catabolic process                                                                                                                |  1.2152163 | 0.5001143 |
| <GO:2000106> | regulation of leukocyte apoptotic process                                                                                                        |  1.2113906 | 0.5001143 |
| <GO:0015833> | peptide transport                                                                                                                                |  1.1991535 | 0.5001143 |
| <GO:0000302> | response to reactive oxygen species                                                                                                              |  1.1950357 | 0.5001143 |
| <GO:0061448> | connective tissue development                                                                                                                    |  1.1927960 | 0.5001143 |
| <GO:0051960> | regulation of nervous system development                                                                                                         |  1.1676082 | 0.5001143 |
| <GO:0010038> | response to metal ion                                                                                                                            |  1.1654568 | 0.5001143 |
| <GO:0050873> | brown fat cell differentiation                                                                                                                   |  1.2671034 | 0.5014925 |
| <GO:0034204> | lipid translocation                                                                                                                              |  1.2550532 | 0.5036030 |
| <GO:0045332> | phospholipid translocation                                                                                                                       |  1.2550532 | 0.5036030 |
| <GO:0032409> | regulation of transporter activity                                                                                                               |  1.2319022 | 0.5041597 |
| <GO:0040007> | growth                                                                                                                                           |  1.1112722 | 0.5041597 |
| <GO:0061082> | myeloid leukocyte cytokine production                                                                                                            |  1.2349257 | 0.5056411 |
| <GO:0055088> | lipid homeostasis                                                                                                                                |  1.2103971 | 0.5079086 |
| <GO:0032414> | positive regulation of ion transmembrane transporter activity                                                                                    |  1.2492656 | 0.5100115 |
| <GO:0008625> | extrinsic apoptotic signaling pathway via death domain receptors                                                                                 |  1.2064990 | 0.5111656 |
| <GO:0030168> | platelet activation                                                                                                                              |  1.2087043 | 0.5113393 |
| <GO:0046632> | alpha-beta T cell differentiation                                                                                                                |  1.1891953 | 0.5113393 |
| <GO:0048864> | stem cell development                                                                                                                            |  1.2578123 | 0.5121503 |
| <GO:0002335> | mature B cell differentiation                                                                                                                    |  1.2534434 | 0.5121503 |
| <GO:0072538> | T-helper 17 type immune response                                                                                                                 |  1.2480401 | 0.5121503 |
| <GO:0045446> | endothelial cell differentiation                                                                                                                 |  1.2063465 | 0.5121503 |
| <GO:0009410> | response to xenobiotic stimulus                                                                                                                  |  1.1592253 | 0.5121503 |
| <GO:0060759> | regulation of response to cytokine stimulus                                                                                                      |  1.1827349 | 0.5148621 |
| <GO:0042327> | positive regulation of phosphorylation                                                                                                           |  1.1750968 | 0.5148621 |
| <GO:1901568> | fatty acid derivative metabolic process                                                                                                          |  1.2427718 | 0.5154902 |
| <GO:0016042> | lipid catabolic process                                                                                                                          |  1.1679631 | 0.5154902 |
| <GO:0006631> | fatty acid metabolic process                                                                                                                     |  1.1582311 | 0.5154902 |
| <GO:0042542> | response to hydrogen peroxide                                                                                                                    |  1.2010193 | 0.5208206 |
| <GO:0045667> | regulation of osteoblast differentiation                                                                                                         |  1.2043048 | 0.5225091 |
| <GO:0045859> | regulation of protein kinase activity                                                                                                            |  1.1768667 | 0.5242848 |
| <GO:0060537> | muscle tissue development                                                                                                                        |  1.1604876 | 0.5242848 |
| <GO:0043433> | negative regulation of DNA-binding transcription factor activity                                                                                 |  1.2362790 | 0.5265226 |
| <GO:0061351> | neural precursor cell proliferation                                                                                                              |  1.2026915 | 0.5275125 |
| <GO:0038066> | p38MAPK cascade                                                                                                                                  |  1.2097449 | 0.5279461 |
| <GO:0003205> | cardiac chamber development                                                                                                                      |  1.2270821 | 0.5284569 |
| <GO:0021885> | forebrain cell migration                                                                                                                         |  1.2241455 | 0.5284569 |
| <GO:0022029> | telencephalon cell migration                                                                                                                     |  1.2241455 | 0.5284569 |
| <GO:2000516> | positive regulation of CD4-positive, alpha-beta T cell activation                                                                                |  1.2339715 | 0.5305489 |
| <GO:0046633> | alpha-beta T cell proliferation                                                                                                                  |  1.2285094 | 0.5305489 |
| <GO:0045055> | regulated exocytosis                                                                                                                             |  1.1717201 | 0.5305489 |
| <GO:0071692> | protein localization to extracellular region                                                                                                     |  1.1506016 | 0.5319114 |
| <GO:0043304> | regulation of mast cell degranulation                                                                                                            |  1.2185278 | 0.5337539 |
| <GO:0002347> | response to tumor cell                                                                                                                           |  1.2318972 | 0.5377242 |
| <GO:0043094> | metabolic compound salvage                                                                                                                       |  1.2249332 | 0.5377859 |
| <GO:0032872> | regulation of stress-activated MAPK cascade                                                                                                      |  1.2155295 | 0.5377859 |
| <GO:0070228> | regulation of lymphocyte apoptotic process                                                                                                       |  1.2240778 | 0.5383576 |
| <GO:0050886> | endocrine process                                                                                                                                |  1.2230721 | 0.5383576 |
| <GO:0050818> | regulation of coagulation                                                                                                                        |  1.2208870 | 0.5384494 |
| <GO:0097305> | response to alcohol                                                                                                                              |  1.1633665 | 0.5384494 |
| <GO:0031100> | animal organ regeneration                                                                                                                        |  1.2280999 | 0.5396248 |
| <GO:0007548> | sex differentiation                                                                                                                              |  1.1825280 | 0.5431122 |
| <GO:0045137> | development of primary sexual characteristics                                                                                                    |  1.1793204 | 0.5435858 |
| <GO:0044706> | multi-multicellular organism process                                                                                                             |  1.1742152 | 0.5435858 |
| <GO:0009914> | hormone transport                                                                                                                                |  1.1738373 | 0.5435858 |
| <GO:0046879> | hormone secretion                                                                                                                                |  1.1738373 | 0.5435858 |
| <GO:0018205> | peptidyl-lysine modification                                                                                                                     | -1.1274298 | 0.5435858 |
| <GO:0055021> | regulation of cardiac muscle tissue growth                                                                                                       |  1.2318028 | 0.5457188 |
| <GO:2000243> | positive regulation of reproductive process                                                                                                      |  1.2240143 | 0.5479963 |
| <GO:0010324> | membrane invagination                                                                                                                            |  1.2188956 | 0.5479963 |
| <GO:0050772> | positive regulation of axonogenesis                                                                                                              |  1.2181617 | 0.5479963 |
| <GO:0012501> | programmed cell death                                                                                                                            |  1.0673154 | 0.5479963 |
| <GO:0045600> | positive regulation of fat cell differentiation                                                                                                  |  1.2216720 | 0.5482115 |
| <GO:0061081> | positive regulation of myeloid leukocyte cytokine production involved in immune response                                                         |  1.2192366 | 0.5482115 |
| <GO:0035019> | somatic stem cell population maintenance                                                                                                         |  1.2279918 | 0.5505081 |
| <GO:1901605> | alpha-amino acid metabolic process                                                                                                               |  1.2128763 | 0.5511619 |
| <GO:0016126> | sterol biosynthetic process                                                                                                                      |  1.2270039 | 0.5517235 |
| <GO:0006739> | NADP+ metabolic process                                                                                                                          |  1.2148277 | 0.5517235 |
| <GO:0051336> | regulation of hydrolase activity                                                                                                                 |  1.1697118 | 0.5528988 |
| <GO:0007167> | enzyme-linked receptor protein signaling pathway                                                                                                 |  1.0929587 | 0.5564961 |
| <GO:0071677> | positive regulation of mononuclear cell migration                                                                                                |  1.1839484 | 0.5581745 |
| <GO:0016486> | peptide hormone processing                                                                                                                       |  1.2206131 | 0.5621860 |
| <GO:0170039> | proteinogenic amino acid metabolic process                                                                                                       |  1.2128737 | 0.5621860 |
| <GO:0071897> | DNA biosynthetic process                                                                                                                         |  1.2090666 | 0.5621860 |
| <GO:0046635> | positive regulation of alpha-beta T cell activation                                                                                              |  1.1895713 | 0.5621860 |
| <GO:0008088> | axo-dendritic transport                                                                                                                          | -1.1072371 | 0.5621860 |
| <GO:0043489> | RNA stabilization                                                                                                                                |  1.1894441 | 0.5634975 |
| <GO:0070232> | regulation of T cell apoptotic process                                                                                                           |  1.2075121 | 0.5640394 |
| <GO:0007229> | integrin-mediated signaling pathway                                                                                                              |  1.1749604 | 0.5640394 |
| <GO:0006915> | apoptotic process                                                                                                                                |  1.0635489 | 0.5640394 |
| <GO:0045577> | regulation of B cell differentiation                                                                                                             |  1.2123621 | 0.5644330 |
| <GO:1905517> | macrophage migration                                                                                                                             |  1.2163373 | 0.5651086 |
| <GO:0071887> | leukocyte apoptotic process                                                                                                                      |  1.1738117 | 0.5651086 |
| <GO:0019217> | regulation of fatty acid metabolic process                                                                                                       |  1.2092568 | 0.5664336 |
| <GO:0008219> | cell death                                                                                                                                       |  1.0641418 | 0.5664336 |
| <GO:0034599> | cellular response to oxidative stress                                                                                                            |  1.1459179 | 0.5664457 |
| <GO:0006513> | protein monoubiquitination                                                                                                                       | -1.2412725 | 0.5671993 |
| <GO:0002637> | regulation of immunoglobulin production                                                                                                          |  1.2072787 | 0.5677340 |
| <GO:0051591> | response to cAMP                                                                                                                                 |  1.2046012 | 0.5677340 |
| <GO:0071496> | cellular response to external stimulus                                                                                                           |  1.1996210 | 0.5677340 |
| <GO:0001890> | placenta development                                                                                                                             |  1.1780062 | 0.5677340 |
| <GO:0009306> | protein secretion                                                                                                                                |  1.1391839 | 0.5677340 |
| <GO:0035592> | establishment of protein localization to extracellular region                                                                                    |  1.1391839 | 0.5677340 |
| <GO:0045927> | positive regulation of growth                                                                                                                    |  1.1689384 | 0.5677511 |
| <GO:0052547> | regulation of peptidase activity                                                                                                                 |  1.2030882 | 0.5694647 |
| <GO:0007266> | Rho protein signal transduction                                                                                                                  |  1.1710478 | 0.5703577 |
| <GO:0048589> | developmental growth                                                                                                                             |  1.1162518 | 0.5726421 |
| <GO:1903426> | regulation of reactive oxygen species biosynthetic process                                                                                       |  1.2032673 | 0.5747293 |
| <GO:0043086> | negative regulation of catalytic activity                                                                                                        |  1.1646866 | 0.5747293 |
| <GO:0071675> | regulation of mononuclear cell migration                                                                                                         |  1.1642986 | 0.5755496 |
| <GO:0014706> | striated muscle tissue development                                                                                                               |  1.1458441 | 0.5755496 |
| <GO:0050808> | synapse organization                                                                                                                             |  1.1323688 | 0.5755496 |
| <GO:0055013> | cardiac muscle cell development                                                                                                                  |  1.1969421 | 0.5798959 |
| <GO:0022604> | regulation of cell morphogenesis                                                                                                                 |  1.1417603 | 0.5798959 |
| <GO:0045648> | positive regulation of erythrocyte differentiation                                                                                               |  1.1990344 | 0.5828895 |
| <GO:0003151> | outflow tract morphogenesis                                                                                                                      |  1.2072572 | 0.5833096 |
| <GO:0090594> | inflammatory response to wounding                                                                                                                |  1.2050849 | 0.5833096 |
| <GO:0001704> | formation of primary germ layer                                                                                                                  |  1.1930607 | 0.5833096 |
| <GO:1902108> | regulation of mitochondrial membrane permeability involved in apoptotic process                                                                  |  1.1923037 | 0.5833096 |
| <GO:1902110> | positive regulation of mitochondrial membrane permeability involved in apoptotic process                                                         |  1.1923037 | 0.5833096 |
| <GO:0050848> | regulation of calcium-mediated signaling                                                                                                         |  1.1792844 | 0.5833096 |
| <GO:0051963> | regulation of synapse assembly                                                                                                                   |  1.1753604 | 0.5833096 |
| <GO:0043604> | amide biosynthetic process                                                                                                                       |  1.1613129 | 0.5833096 |
| <GO:0006029> | proteoglycan metabolic process                                                                                                                   |  1.1670737 | 0.5833161 |
| <GO:0009612> | response to mechanical stimulus                                                                                                                  |  1.1651483 | 0.5835285 |
| <GO:1902369> | negative regulation of RNA catabolic process                                                                                                     |  1.1643917 | 0.5847285 |
| <GO:0050863> | regulation of T cell activation                                                                                                                  |  1.1194208 | 0.5868553 |
| <GO:0010717> | regulation of epithelial to mesenchymal transition                                                                                               |  1.1985170 | 0.5871671 |
| <GO:0035337> | fatty-acyl-CoA metabolic process                                                                                                                 |  1.1936410 | 0.5871671 |
| <GO:0140448> | signaling receptor ligand precursor processing                                                                                                   |  1.1926576 | 0.5871671 |
| <GO:0051053> | negative regulation of DNA metabolic process                                                                                                     |  1.1919273 | 0.5871671 |
| <GO:0044703> | multi-organism reproductive process                                                                                                              |  1.1636134 | 0.5871671 |
| <GO:1903573> | negative regulation of response to endoplasmic reticulum stress                                                                                  |  1.1877210 | 0.5874468 |
| <GO:0007040> | lysosome organization                                                                                                                            |  1.1415794 | 0.5874468 |
| <GO:0080171> | lytic vacuole organization                                                                                                                       |  1.1415794 | 0.5874468 |
| <GO:0061564> | axon development                                                                                                                                 |  1.1228046 | 0.5874468 |
| <GO:0010828> | positive regulation of D-glucose transmembrane transport                                                                                         |  1.1865418 | 0.5912215 |
| <GO:0015698> | inorganic anion transport                                                                                                                        |  1.1889114 | 0.5948188 |
| <GO:0008406> | gonad development                                                                                                                                |  1.1589720 | 0.5948188 |
| <GO:0006814> | sodium ion transport                                                                                                                             |  1.1573196 | 0.5963567 |
| <GO:0032879> | regulation of localization                                                                                                                       |  1.0551933 | 0.5964878 |
| <GO:0048585> | negative regulation of response to stimulus                                                                                                      |  1.0589271 | 0.5984025 |
| <GO:2000107> | negative regulation of leukocyte apoptotic process                                                                                               |  1.1889589 | 0.5991717 |
| <GO:0002706> | regulation of lymphocyte mediated immunity                                                                                                       |  1.1567844 | 0.6003529 |
| <GO:2001238> | positive regulation of extrinsic apoptotic signaling pathway                                                                                     |  1.1917351 | 0.6005988 |
| <GO:0030203> | glycosaminoglycan metabolic process                                                                                                              |  1.1820039 | 0.6020364 |
| <GO:0007416> | synapse assembly                                                                                                                                 |  1.1557888 | 0.6029485 |
| <GO:2000008> | regulation of protein localization to cell surface                                                                                               |  1.1955972 | 0.6031099 |
| <GO:0010718> | positive regulation of epithelial to mesenchymal transition                                                                                      |  1.1860169 | 0.6031099 |
| <GO:0097035> | regulation of membrane lipid distribution                                                                                                        |  1.1811894 | 0.6031099 |
| <GO:0030900> | forebrain development                                                                                                                            |  1.1284802 | 0.6038182 |
| <GO:0042593> | glucose homeostasis                                                                                                                              |  1.1508062 | 0.6038624 |
| <GO:0042325> | regulation of phosphorylation                                                                                                                    |  1.1088753 | 0.6039862 |
| <GO:1902531> | regulation of intracellular signal transduction                                                                                                  |  1.0530896 | 0.6039862 |
| <GO:0051865> | protein autoubiquitination                                                                                                                       | -1.1604708 | 0.6049908 |
| <GO:0097306> | cellular response to alcohol                                                                                                                     |  1.1618419 | 0.6053709 |
| <GO:1902373> | negative regulation of mRNA catabolic process                                                                                                    |  1.1616461 | 0.6053709 |
| <GO:0046883> | regulation of hormone secretion                                                                                                                  |  1.1497928 | 0.6053709 |
| <GO:0044092> | negative regulation of molecular function                                                                                                        |  1.1182990 | 0.6055159 |
| <GO:0010638> | positive regulation of organelle organization                                                                                                    |  1.1091154 | 0.6055159 |
| <GO:0120163> | negative regulation of cold-induced thermogenesis                                                                                                |  1.1850629 | 0.6060784 |
| <GO:0006643> | membrane lipid metabolic process                                                                                                                 |  1.1392959 | 0.6074380 |
| <GO:0045445> | myoblast differentiation                                                                                                                         |  1.1597949 | 0.6074931 |
| <GO:0001843> | neural tube closure                                                                                                                              |  1.1924211 | 0.6090298 |
| <GO:0060606> | tube closure                                                                                                                                     |  1.1924211 | 0.6090298 |
| <GO:0030851> | granulocyte differentiation                                                                                                                      |  1.1900120 | 0.6090298 |
| <GO:0002260> | lymphocyte homeostasis                                                                                                                           |  1.1820986 | 0.6090298 |
| <GO:0015837> | amine transport                                                                                                                                  |  1.1857631 | 0.6166430 |
| <GO:0006998> | nuclear envelope organization                                                                                                                    | -1.2327636 | 0.6187744 |
| <GO:0015918> | sterol transport                                                                                                                                 |  1.1518726 | 0.6187744 |
| <GO:0051174> | regulation of phosphorus metabolic process                                                                                                       |  1.0895048 | 0.6187744 |
| <GO:0007613> | memory                                                                                                                                           |  1.1547300 | 0.6188320 |
| <GO:1902914> | regulation of protein polyubiquitination                                                                                                         |  1.1818642 | 0.6189932 |
| <GO:1901890> | positive regulation of cell junction assembly                                                                                                    |  1.1834810 | 0.6197286 |
| <GO:0055006> | cardiac cell development                                                                                                                         |  1.1768267 | 0.6197286 |
| <GO:0043029> | T cell homeostasis                                                                                                                               |  1.1759315 | 0.6197286 |
| <GO:2000316> | regulation of T-helper 17 type immune response                                                                                                   |  1.1738810 | 0.6197286 |
| <GO:0019827> | stem cell population maintenance                                                                                                                 |  1.1707143 | 0.6197286 |
| <GO:0098727> | maintenance of cell number                                                                                                                       |  1.1707143 | 0.6197286 |
| <GO:0048255> | mRNA stabilization                                                                                                                               |  1.1705542 | 0.6197286 |
| <GO:0001838> | embryonic epithelial tube formation                                                                                                              |  1.1657408 | 0.6197286 |
| <GO:0016331> | morphogenesis of embryonic epithelium                                                                                                            |  1.1657408 | 0.6197286 |
| <GO:0045165> | cell fate commitment                                                                                                                             |  1.1459030 | 0.6197286 |
| <GO:0002698> | negative regulation of immune effector process                                                                                                   |  1.1424438 | 0.6197286 |
| <GO:0033500> | carbohydrate homeostasis                                                                                                                         |  1.1363118 | 0.6197286 |
| <GO:0031667> | response to nutrient levels                                                                                                                      |  1.0862604 | 0.6197286 |
| <GO:0045010> | actin nucleation                                                                                                                                 | -1.0635309 | 0.6197286 |
| <GO:0045663> | positive regulation of myoblast differentiation                                                                                                  |  1.1673752 | 0.6198380 |
| <GO:0032506> | cytokinetic process                                                                                                                              | -1.1683712 | 0.6213762 |
| <GO:0006637> | acyl-CoA metabolic process                                                                                                                       |  1.1679183 | 0.6225211 |
| <GO:0035383> | thioester metabolic process                                                                                                                      |  1.1679183 | 0.6225211 |
| <GO:0006693> | prostaglandin metabolic process                                                                                                                  |  1.1653570 | 0.6225211 |
| <GO:0035051> | cardiocyte differentiation                                                                                                                       |  1.1517166 | 0.6225211 |
| <GO:0050803> | regulation of synapse structure or activity                                                                                                      |  1.1264942 | 0.6225211 |
| <GO:0046902> | regulation of mitochondrial membrane permeability                                                                                                |  1.1770637 | 0.6238212 |
| <GO:0006941> | striated muscle contraction                                                                                                                      |  1.1660203 | 0.6253801 |
| <GO:0034142> | toll-like receptor 4 signaling pathway                                                                                                           |  1.1481432 | 0.6253801 |
| <GO:0050869> | negative regulation of B cell activation                                                                                                         |  1.1737133 | 0.6256678 |
| <GO:1903018> | regulation of glycoprotein metabolic process                                                                                                     |  1.1707823 | 0.6256981 |
| <GO:0007601> | visual perception                                                                                                                                |  1.1706142 | 0.6256981 |
| <GO:0050953> | sensory perception of light stimulus                                                                                                             |  1.1706142 | 0.6256981 |
| <GO:0070231> | T cell apoptotic process                                                                                                                         |  1.1701965 | 0.6269200 |
| <GO:1901136> | carbohydrate derivative catabolic process                                                                                                        |  1.1079718 | 0.6269200 |
| <GO:0001933> | negative regulation of protein phosphorylation                                                                                                   |  1.1485263 | 0.6279647 |
| <GO:0030301> | cholesterol transport                                                                                                                            |  1.1467424 | 0.6300750 |
| <GO:0050708> | regulation of protein secretion                                                                                                                  |  1.1365338 | 0.6300750 |
| <GO:0001894> | tissue homeostasis                                                                                                                               |  1.1228279 | 0.6300750 |
| <GO:0060249> | anatomical structure homeostasis                                                                                                                 |  1.1228279 | 0.6300750 |
| <GO:1990830> | cellular response to leukemia inhibitory factor                                                                                                  |  1.1616044 | 0.6311753 |
| <GO:0006469> | negative regulation of protein kinase activity                                                                                                   |  1.1684751 | 0.6324359 |
| <GO:0006352> | DNA-templated transcription initiation                                                                                                           |  1.1347842 | 0.6324359 |
| <GO:0015807> | L-amino acid transport                                                                                                                           |  1.1740670 | 0.6325915 |
| <GO:0003158> | endothelium development                                                                                                                          |  1.1341289 | 0.6325915 |
| <GO:0031641> | regulation of myelination                                                                                                                        |  1.1662563 | 0.6332240 |
| <GO:0046514> | ceramide catabolic process                                                                                                                       |  1.1592074 | 0.6332240 |
| <GO:0070227> | lymphocyte apoptotic process                                                                                                                     |  1.1589444 | 0.6332240 |
| <GO:0048608> | reproductive structure development                                                                                                               |  1.1276142 | 0.6332240 |
| <GO:0048545> | response to steroid hormone                                                                                                                      |  1.1096420 | 0.6332240 |
| <GO:0001778> | plasma membrane repair                                                                                                                           |  1.1666442 | 0.6336024 |
| <GO:0097300> | programmed necrotic cell death                                                                                                                   |  1.1650691 | 0.6336024 |
| <GO:0048246> | macrophage chemotaxis                                                                                                                            |  1.1645335 | 0.6336024 |
| <GO:0048260> | positive regulation of receptor-mediated endocytosis                                                                                             |  1.1636556 | 0.6347431 |
| <GO:0007369> | gastrulation                                                                                                                                     |  1.1382257 | 0.6347431 |
| <GO:0043270> | positive regulation of monoatomic ion transport                                                                                                  |  1.1363571 | 0.6347431 |
| <GO:0016051> | carbohydrate biosynthetic process                                                                                                                |  1.1224004 | 0.6347431 |
| <GO:0019220> | regulation of phosphate metabolic process                                                                                                        |  1.0797544 | 0.6347431 |
| <GO:0097190> | apoptotic signaling pathway                                                                                                                      |  1.0780655 | 0.6347431 |
| <GO:0170033> | L-amino acid metabolic process                                                                                                                   |  1.1704722 | 0.6355867 |
| <GO:0098586> | cellular response to virus                                                                                                                       |  1.1699164 | 0.6355867 |
| <GO:0042130> | negative regulation of T cell proliferation                                                                                                      |  1.1669312 | 0.6355867 |
| <GO:0046640> | regulation of alpha-beta T cell proliferation                                                                                                    |  1.1608675 | 0.6355867 |
| <GO:0032438> | melanosome organization                                                                                                                          |  1.1582498 | 0.6355867 |
| <GO:0048753> | pigment granule organization                                                                                                                     |  1.1582498 | 0.6355867 |
| <GO:0010565> | regulation of ketone metabolic process                                                                                                           |  1.1556911 | 0.6355867 |
| <GO:0019915> | lipid storage                                                                                                                                    |  1.1550292 | 0.6355867 |
| <GO:0170062> | nutrient storage                                                                                                                                 |  1.1550292 | 0.6355867 |
| <GO:0021915> | neural tube development                                                                                                                          |  1.1423190 | 0.6355867 |
| <GO:0002790> | peptide secretion                                                                                                                                |  1.1324597 | 0.6355867 |
| <GO:0098927> | vesicle-mediated transport between endosomal compartments                                                                                        | -1.1163648 | 0.6355867 |
| <GO:0048010> | vascular endothelial growth factor receptor signaling pathway                                                                                    |  1.1638233 | 0.6365752 |
| <GO:0051403> | stress-activated MAPK cascade                                                                                                                    |  1.1602702 | 0.6372922 |
| <GO:0031529> | ruffle organization                                                                                                                              |  1.1648290 | 0.6395939 |
| <GO:0071384> | cellular response to corticosteroid stimulus                                                                                                     |  1.1623658 | 0.6395939 |
| <GO:0022898> | regulation of transmembrane transporter activity                                                                                                 |  1.1535451 | 0.6395939 |
| <GO:0032412> | regulation of monoatomic ion transmembrane transporter activity                                                                                  |  1.1535451 | 0.6395939 |
| <GO:0090263> | positive regulation of canonical Wnt signaling pathway                                                                                           |  1.1498782 | 0.6395939 |
| <GO:0002520> | immune system development                                                                                                                        |  1.1245602 | 0.6395939 |
| <GO:0044403> | biological process involved in symbiotic interaction                                                                                             |  1.1155003 | 0.6395939 |
| <GO:0006629> | lipid metabolic process                                                                                                                          |  1.0603893 | 0.6395939 |
| <GO:0051338> | regulation of transferase activity                                                                                                               |  1.1042792 | 0.6400651 |
| <GO:0002042> | cell migration involved in sprouting angiogenesis                                                                                                |  1.1467442 | 0.6417657 |
| <GO:0006636> | unsaturated fatty acid biosynthetic process                                                                                                      |  1.1522784 | 0.6429181 |
| <GO:0040008> | regulation of growth                                                                                                                             |  1.0871252 | 0.6429181 |
| <GO:0002262> | myeloid cell homeostasis                                                                                                                         |  1.1090025 | 0.6440666 |
| <GO:0061157> | mRNA destabilization                                                                                                                             | -0.9942560 | 0.6440666 |
| <GO:0007009> | plasma membrane organization                                                                                                                     |  1.1325038 | 0.6454566 |
| <GO:0060761> | negative regulation of response to cytokine stimulus                                                                                             |  1.1583687 | 0.6468974 |
| <GO:0048598> | embryonic morphogenesis                                                                                                                          |  1.0976229 | 0.6468974 |
| <GO:0003014> | renal system process                                                                                                                             |  1.1600932 | 0.6469455 |
| <GO:0071470> | cellular response to osmotic stress                                                                                                              |  1.1560074 | 0.6469455 |
| <GO:0043549> | regulation of kinase activity                                                                                                                    |  1.1112365 | 0.6474889 |
| <GO:0002791> | regulation of peptide secretion                                                                                                                  |  1.1304094 | 0.6479051 |
| <GO:0090087> | regulation of peptide transport                                                                                                                  |  1.1304094 | 0.6479051 |
| <GO:0009266> | response to temperature stimulus                                                                                                                 |  1.1254235 | 0.6479051 |
| <GO:0034138> | toll-like receptor 3 signaling pathway                                                                                                           | -1.0949236 | 0.6479051 |
| <GO:0043330> | response to exogenous dsRNA                                                                                                                      |  1.1493958 | 0.6483562 |
| <GO:2000179> | positive regulation of neural precursor cell proliferation                                                                                       |  1.1472821 | 0.6496653 |
| <GO:0060070> | canonical Wnt signaling pathway                                                                                                                  |  1.1092785 | 0.6598393 |
| <GO:0032607> | interferon-alpha production                                                                                                                      |  1.1450981 | 0.6602316 |
| <GO:0032647> | regulation of interferon-alpha production                                                                                                        |  1.1450981 | 0.6602316 |
| <GO:0030101> | natural killer cell activation                                                                                                                   |  1.1348879 | 0.6602316 |
| <GO:0030177> | positive regulation of Wnt signaling pathway                                                                                                     |  1.1330590 | 0.6602316 |
| <GO:0030258> | lipid modification                                                                                                                               |  1.1236833 | 0.6602316 |
| <GO:2001233> | regulation of apoptotic signaling pathway                                                                                                        |  1.0861141 | 0.6602316 |
| <GO:0051897> | positive regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction                                                        |  1.1113151 | 0.6610604 |
| <GO:0051099> | positive regulation of binding                                                                                                                   |  1.1412566 | 0.6638240 |
| <GO:0001959> | regulation of cytokine-mediated signaling pathway                                                                                                |  1.1118468 | 0.6687618 |
| <GO:0061458> | reproductive system development                                                                                                                  |  1.1114242 | 0.6704048 |
| <GO:0180064> | protein O-linked glycosylation via xylose                                                                                                        |  1.1464417 | 0.6754800 |
| <GO:0046513> | ceramide biosynthetic process                                                                                                                    |  1.1368687 | 0.6754800 |
| <GO:0070059> | intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress                                                                |  1.1354314 | 0.6754800 |
| <GO:0071900> | regulation of protein serine/threonine kinase activity                                                                                           |  1.1245576 | 0.6754800 |
| <GO:0050807> | regulation of synapse organization                                                                                                               |  1.0985562 | 0.6754800 |
| <GO:0042269> | regulation of natural killer cell mediated cytotoxicity                                                                                          | -1.0740093 | 0.6754800 |
| <GO:0001836> | release of cytochrome c from mitochondria                                                                                                        |  1.1399561 | 0.6759035 |
| <GO:0061180> | mammary gland epithelium development                                                                                                             |  1.1389196 | 0.6761027 |
| <GO:0046847> | filopodium assembly                                                                                                                              |  1.1384250 | 0.6761027 |
| <GO:0044409> | symbiont entry into host                                                                                                                         |  1.1176364 | 0.6786462 |
| <GO:0034330> | cell junction organization                                                                                                                       |  1.0696536 | 0.6804035 |
| <GO:0001773> | myeloid dendritic cell activation                                                                                                                |  1.1299079 | 0.6816803 |
| <GO:0006520> | amino acid metabolic process                                                                                                                     |  1.1179699 | 0.6816803 |
| <GO:0002793> | positive regulation of peptide secretion                                                                                                         |  1.1194168 | 0.6821363 |
| <GO:0046434> | organophosphate catabolic process                                                                                                                |  1.0996458 | 0.6831885 |
| <GO:0021879> | forebrain neuron differentiation                                                                                                                 |  1.1304414 | 0.6836762 |
| <GO:0034284> | response to monosaccharide                                                                                                                       |  1.1133369 | 0.6836762 |
| <GO:0003002> | regionalization                                                                                                                                  |  1.1121903 | 0.6836762 |
| <GO:0071398> | cellular response to fatty acid                                                                                                                  |  1.1283169 | 0.6853803 |
| <GO:0007160> | cell-matrix adhesion                                                                                                                             |  1.1016338 | 0.6875676 |
| <GO:0003231> | cardiac ventricle development                                                                                                                    |  1.1339859 | 0.6882307 |
| <GO:0042267> | natural killer cell mediated cytotoxicity                                                                                                        | -1.1483391 | 0.6906835 |
| <GO:0048839> | inner ear development                                                                                                                            |  1.1323918 | 0.6906835 |
| <GO:0072599> | establishment of protein localization to endoplasmic reticulum                                                                                   | -1.1311548 | 0.6906835 |
| <GO:0140115> | export across plasma membrane                                                                                                                    |  1.1308373 | 0.6906835 |
| <GO:0098661> | inorganic anion transmembrane transport                                                                                                          |  1.1307932 | 0.6906835 |
| <GO:0043001> | Golgi to plasma membrane protein transport                                                                                                       |  1.1264414 | 0.6906835 |
| <GO:0021795> | cerebral cortex cell migration                                                                                                                   |  1.1262464 | 0.6906835 |
| <GO:0051101> | regulation of DNA binding                                                                                                                        |  1.1228330 | 0.6906835 |
| <GO:0033059> | cellular pigmentation                                                                                                                            |  1.1223199 | 0.6906835 |
| <GO:0045931> | positive regulation of mitotic cell cycle                                                                                                        |  1.1216680 | 0.6906835 |
| <GO:0061515> | myeloid cell development                                                                                                                         |  1.1113758 | 0.6906835 |
| <GO:0071248> | cellular response to metal ion                                                                                                                   |  1.1045012 | 0.6906835 |
| <GO:1901888> | regulation of cell junction assembly                                                                                                             |  1.1011388 | 0.6906835 |
| <GO:0002753> | cytoplasmic pattern recognition receptor signaling pathway                                                                                       |  1.0808938 | 0.6906835 |
| <GO:0060538> | skeletal muscle organ development                                                                                                                |  1.1104604 | 0.6908940 |
| <GO:0016925> | protein sumoylation                                                                                                                              | -1.0808735 | 0.6908940 |
| <GO:0031398> | positive regulation of protein ubiquitination                                                                                                    |  1.1281863 | 0.6930217 |
| <GO:0051783> | regulation of nuclear division                                                                                                                   |  1.1279713 | 0.6950291 |
| <GO:0010595> | positive regulation of endothelial cell migration                                                                                                |  1.1133260 | 0.6955531 |
| <GO:0001960> | negative regulation of cytokine-mediated signaling pathway                                                                                       |  1.1182849 | 0.6964494 |
| <GO:0071356> | cellular response to tumor necrosis factor                                                                                                       |  1.0879113 | 0.6964494 |
| <GO:0016188> | synaptic vesicle maturation                                                                                                                      |  1.1222231 | 0.6966965 |
| <GO:0060047> | heart contraction                                                                                                                                |  1.1123430 | 0.6966965 |
| <GO:0030183> | B cell differentiation                                                                                                                           |  1.1027415 | 0.6966965 |
| <GO:0010629> | negative regulation of gene expression                                                                                                           |  1.0462642 | 0.6972883 |
| <GO:0003073> | regulation of systemic arterial blood pressure                                                                                                   |  1.1235991 | 0.6979298 |
| <GO:0033028> | myeloid cell apoptotic process                                                                                                                   |  1.1224668 | 0.6997934 |
| <GO:0071260> | cellular response to mechanical stimulus                                                                                                         |  1.1190570 | 0.6997934 |
| <GO:0007265> | Ras protein signal transduction                                                                                                                  |  1.1042880 | 0.7041987 |
| <GO:0043302> | positive regulation of leukocyte degranulation                                                                                                   |  1.1136570 | 0.7045167 |
| <GO:0045123> | cellular extravasation                                                                                                                           |  1.1180077 | 0.7045522 |
| <GO:0046425> | regulation of receptor signaling pathway via JAK-STAT                                                                                            |  1.1234555 | 0.7064259 |
| <GO:1904892> | regulation of receptor signaling pathway via STAT                                                                                                |  1.1234555 | 0.7064259 |
| <GO:0046503> | glycerolipid catabolic process                                                                                                                   |  1.1223298 | 0.7064259 |
| <GO:0035725> | sodium ion transmembrane transport                                                                                                               |  1.1212516 | 0.7064259 |
| <GO:0002823> | negative regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains |  1.1208395 | 0.7064259 |
| <GO:1903078> | positive regulation of protein localization to plasma membrane                                                                                   |  1.1166602 | 0.7064259 |
| <GO:0031670> | cellular response to nutrient                                                                                                                    |  1.1106770 | 0.7064259 |
| <GO:0031401> | positive regulation of protein modification process                                                                                              |  1.0729916 | 0.7064259 |
| <GO:1901655> | cellular response to ketone                                                                                                                      |  1.1096874 | 0.7071747 |
| <GO:0045860> | positive regulation of protein kinase activity                                                                                                   |  1.1080013 | 0.7072930 |
| <GO:0007519> | skeletal muscle tissue development                                                                                                               |  1.1025483 | 0.7074441 |
| <GO:0008637> | apoptotic mitochondrial changes                                                                                                                  |  1.1032958 | 0.7096898 |
| <GO:0031098> | stress-activated protein kinase signaling cascade                                                                                                |  1.1102547 | 0.7112730 |
| <GO:0071495> | cellular response to endogenous stimulus                                                                                                         |  1.0416029 | 0.7112730 |
| <GO:0005975> | carbohydrate metabolic process                                                                                                                   |  1.0529125 | 0.7127919 |
| <GO:0046496> | nicotinamide nucleotide metabolic process                                                                                                        |  1.0874368 | 0.7129442 |
| <GO:0071695> | anatomical structure maturation                                                                                                                  |  1.0906442 | 0.7137478 |
| <GO:1905521> | regulation of macrophage migration                                                                                                               |  1.1135542 | 0.7141245 |
| <GO:0003015> | heart process                                                                                                                                    |  1.0983801 | 0.7141245 |
| <GO:0031294> | lymphocyte costimulation                                                                                                                         |  1.1126905 | 0.7178860 |
| <GO:0031295> | T cell costimulation                                                                                                                             |  1.1126905 | 0.7178860 |
| <GO:0097401> | synaptic vesicle lumen acidification                                                                                                             |  1.1124740 | 0.7178860 |
| <GO:0031644> | regulation of nervous system process                                                                                                             |  1.1089926 | 0.7178860 |
| <GO:0051781> | positive regulation of cell division                                                                                                             |  1.1072171 | 0.7178860 |
| <GO:0048469> | cell maturation                                                                                                                                  |  1.0893140 | 0.7182590 |
| <GO:0035384> | thioester biosynthetic process                                                                                                                   |  1.1161489 | 0.7188599 |
| <GO:0071616> | acyl-CoA biosynthetic process                                                                                                                    |  1.1161489 | 0.7188599 |
| <GO:0009408> | response to heat                                                                                                                                 |  1.1072379 | 0.7188599 |
| <GO:0030048> | actin filament-based movement                                                                                                                    |  1.1124361 | 0.7206622 |
| <GO:0060612> | adipose tissue development                                                                                                                       |  1.1039168 | 0.7235056 |
| <GO:1901186> | positive regulation of ERBB signaling pathway                                                                                                    |  1.1083410 | 0.7241818 |
| <GO:0071322> | cellular response to carbohydrate stimulus                                                                                                       |  1.0950125 | 0.7241818 |
| <GO:0048534> | hematopoietic or lymphoid organ development                                                                                                      |  1.0893499 | 0.7241818 |
| <GO:0001841> | neural tube formation                                                                                                                            |  1.1135242 | 0.7243300 |
| <GO:0097194> | execution phase of apoptosis                                                                                                                     |  1.1107837 | 0.7243300 |
| <GO:0006493> | protein O-linked glycosylation                                                                                                                   |  1.0948122 | 0.7243300 |
| <GO:0010976> | positive regulation of neuron projection development                                                                                             |  1.0943575 | 0.7243300 |
| <GO:0072524> | pyridine-containing compound metabolic process                                                                                                   |  1.0785739 | 0.7243300 |
| <GO:0009725> | response to hormone                                                                                                                              |  1.0422744 | 0.7243300 |
| <GO:0014020> | primary neural tube formation                                                                                                                    |  1.1094598 | 0.7247146 |
| <GO:0046627> | negative regulation of insulin receptor signaling pathway                                                                                        |  1.0938105 | 0.7247647 |
| <GO:1900077> | negative regulation of cellular response to insulin stimulus                                                                                     |  1.0938105 | 0.7247647 |
| <GO:0015850> | organic hydroxy compound transport                                                                                                               |  1.0789118 | 0.7247647 |
| <GO:0030111> | regulation of Wnt signaling pathway                                                                                                              |  1.0780839 | 0.7247647 |
| <GO:0044090> | positive regulation of vacuole organization                                                                                                      |  1.1048402 | 0.7258065 |
| <GO:2000241> | regulation of reproductive process                                                                                                               |  1.1074808 | 0.7300133 |
| <GO:0006111> | regulation of gluconeogenesis                                                                                                                    |  1.1032436 | 0.7300133 |
| <GO:0043303> | mast cell degranulation                                                                                                                          |  1.0955868 | 0.7300133 |
| <GO:0043255> | regulation of carbohydrate biosynthetic process                                                                                                  |  1.1042304 | 0.7335639 |
| <GO:0021700> | developmental maturation                                                                                                                         |  1.0642010 | 0.7341534 |
| <GO:0009101> | glycoprotein biosynthetic process                                                                                                                |  1.0790278 | 0.7363486 |
| <GO:0002448> | mast cell mediated immunity                                                                                                                      |  1.1061607 | 0.7376736 |
| <GO:0009311> | oligosaccharide metabolic process                                                                                                                |  1.0872183 | 0.7394145 |
| <GO:0003006> | developmental process involved in reproduction                                                                                                   |  1.0522751 | 0.7394145 |
| <GO:0008610> | lipid biosynthetic process                                                                                                                       |  1.0452249 | 0.7394145 |
| <GO:0034612> | response to tumor necrosis factor                                                                                                                |  1.0700168 | 0.7395903 |
| <GO:1905953> | negative regulation of lipid localization                                                                                                        |  1.0848169 | 0.7411985 |
| <GO:1902686> | mitochondrial outer membrane permeabilization involved in programmed cell death                                                                  |  1.0900936 | 0.7414362 |
| <GO:0098869> | cellular oxidant detoxification                                                                                                                  |  1.0896444 | 0.7414362 |
| <GO:0007389> | pattern specification process                                                                                                                    |  1.0753083 | 0.7453804 |
| <GO:0016242> | negative regulation of macroautophagy                                                                                                            |  1.0971881 | 0.7461464 |
| <GO:0030073> | insulin secretion                                                                                                                                |  1.0844516 | 0.7464897 |
| <GO:0002562> | somatic diversification of immune receptors via germline recombination within a single locus                                                     |  1.0964390 | 0.7465825 |
| <GO:0016444> | somatic cell DNA recombination                                                                                                                   |  1.0964390 | 0.7465825 |
| <GO:0042157> | lipoprotein metabolic process                                                                                                                    |  1.0884376 | 0.7465825 |
| <GO:0043085> | positive regulation of catalytic activity                                                                                                        |  1.0604256 | 0.7465825 |
| <GO:0043583> | ear development                                                                                                                                  |  1.0950317 | 0.7471153 |
| <GO:0002832> | negative regulation of response to biotic stimulus                                                                                               |  1.0807502 | 0.7471153 |
| <GO:0050853> | B cell receptor signaling pathway                                                                                                                |  1.0941640 | 0.7472682 |
| <GO:0050796> | regulation of insulin secretion                                                                                                                  |  1.0803954 | 0.7479893 |
| <GO:1904019> | epithelial cell apoptotic process                                                                                                                |  1.0908123 | 0.7482096 |
| <GO:0001667> | ameboidal-type cell migration                                                                                                                    |  1.0709694 | 0.7493543 |
| <GO:0046626> | regulation of insulin receptor signaling pathway                                                                                                 |  1.0870218 | 0.7496271 |
| <GO:0032233> | positive regulation of actin filament bundle assembly                                                                                            |  1.0882634 | 0.7547721 |
| <GO:0030178> | negative regulation of Wnt signaling pathway                                                                                                     |  1.0772530 | 0.7572350 |
| <GO:0048738> | cardiac muscle tissue development                                                                                                                |  1.0767907 | 0.7585592 |
| <GO:0002820> | negative regulation of adaptive immune response                                                                                                  |  1.0740907 | 0.7588548 |
| <GO:0051347> | positive regulation of transferase activity                                                                                                      |  1.0717835 | 0.7588548 |
| <GO:0099536> | synaptic signaling                                                                                                                               |  1.0459642 | 0.7637243 |
| <GO:0090207> | regulation of triglyceride metabolic process                                                                                                     |  1.0806333 | 0.7658826 |
| <GO:0033673> | negative regulation of kinase activity                                                                                                           |  1.0801630 | 0.7679356 |
| <GO:1901699> | cellular response to nitrogen compound                                                                                                           |  1.0319125 | 0.7679356 |
| <GO:0034763> | negative regulation of transmembrane transport                                                                                                   |  1.0756991 | 0.7691281 |
| <GO:0055007> | cardiac muscle cell differentiation                                                                                                              |  1.0816677 | 0.7696156 |
| <GO:1900076> | regulation of cellular response to insulin stimulus                                                                                              |  1.0815901 | 0.7696156 |
| <GO:0090114> | COPII-coated vesicle budding                                                                                                                     | -1.0306062 | 0.7696156 |
| <GO:0051310> | metaphase chromosome alignment                                                                                                                   | -1.1801725 | 0.7700048 |
| <GO:0097193> | intrinsic apoptotic signaling pathway                                                                                                            |  1.0521133 | 0.7700048 |
| <GO:0010558> | negative regulation of macromolecule biosynthetic process                                                                                        |  1.0223889 | 0.7700048 |
| <GO:0006692> | prostanoid metabolic process                                                                                                                     |  1.0796048 | 0.7704608 |
| <GO:0060541> | respiratory system development                                                                                                                   |  1.0676890 | 0.7717753 |
| <GO:0043536> | positive regulation of blood vessel endothelial cell migration                                                                                   |  1.0793401 | 0.7727575 |
| <GO:0002931> | response to ischemia                                                                                                                             |  1.0721715 | 0.7727575 |
| <GO:0042180> | ketone metabolic process                                                                                                                         |  1.0685326 | 0.7727575 |
| <GO:0007219> | Notch signaling pathway                                                                                                                          |  1.0679979 | 0.7727575 |
| <GO:0046661> | male sex differentiation                                                                                                                         |  1.0649404 | 0.7727575 |
| <GO:0048706> | embryonic skeletal system development                                                                                                            |  1.0631959 | 0.7727575 |
| <GO:0045124> | regulation of bone resorption                                                                                                                    |  1.0606927 | 0.7727575 |
| <GO:0032689> | negative regulation of type II interferon production                                                                                             |  1.0602457 | 0.7727575 |
| <GO:0048729> | tissue morphogenesis                                                                                                                             |  1.0485331 | 0.7727575 |
| <GO:0051259> | protein complex oligomerization                                                                                                                  |  1.0606723 | 0.7727633 |
| <GO:0060048> | cardiac muscle contraction                                                                                                                       |  1.0716670 | 0.7744281 |
| <GO:0006022> | aminoglycan metabolic process                                                                                                                    |  1.0784620 | 0.7755395 |
| <GO:1902116> | negative regulation of organelle assembly                                                                                                        |  1.0781066 | 0.7755395 |
| <GO:0060411> | cardiac septum morphogenesis                                                                                                                     |  1.0712132 | 0.7755395 |
| <GO:0034605> | cellular response to heat                                                                                                                        |  1.0697859 | 0.7755395 |
| <GO:0099072> | regulation of postsynaptic membrane neurotransmitter receptor levels                                                                             |  1.0628271 | 0.7755395 |
| <GO:0034440> | lipid oxidation                                                                                                                                  |  1.0699248 | 0.7805222 |
| <GO:0007224> | smoothened signaling pathway                                                                                                                     |  1.0615338 | 0.7816136 |
| <GO:0030326> | embryonic limb morphogenesis                                                                                                                     |  1.0711826 | 0.7830400 |
| <GO:0035113> | embryonic appendage morphogenesis                                                                                                                |  1.0711826 | 0.7830400 |
| <GO:0032868> | response to insulin                                                                                                                              |  1.0450036 | 0.7839319 |
| <GO:0070085> | glycosylation                                                                                                                                    |  1.0696866 | 0.7845698 |
| <GO:0034103> | regulation of tissue remodeling                                                                                                                  |  1.0584726 | 0.7845698 |
| <GO:0030072> | peptide hormone secretion                                                                                                                        |  1.0543817 | 0.7845698 |
| <GO:0009132> | nucleoside diphosphate metabolic process                                                                                                         |  1.0491029 | 0.7845698 |
| <GO:0009135> | purine nucleoside diphosphate metabolic process                                                                                                  |  1.0491029 | 0.7845698 |
| <GO:0009179> | purine ribonucleoside diphosphate metabolic process                                                                                              |  1.0491029 | 0.7845698 |
| <GO:0009185> | ribonucleoside diphosphate metabolic process                                                                                                     |  1.0491029 | 0.7845698 |
| <GO:0051496> | positive regulation of stress fiber assembly                                                                                                     |  1.0679513 | 0.7858161 |
| <GO:0019076> | viral release from host cell                                                                                                                     | -1.0559241 | 0.7858161 |
| <GO:0035891> | exit from host cell                                                                                                                              | -1.0559241 | 0.7858161 |
| <GO:0046605> | regulation of centrosome cycle                                                                                                                   | -1.0559241 | 0.7858161 |
| <GO:0010586> | miRNA metabolic process                                                                                                                          |  1.0597365 | 0.7881453 |
| <GO:0099537> | trans-synaptic signaling                                                                                                                         |  1.0402446 | 0.7881453 |
| <GO:0022612> | gland morphogenesis                                                                                                                              |  1.0611400 | 0.7891918 |
| <GO:0008593> | regulation of Notch signaling pathway                                                                                                            |  1.0604153 | 0.7891918 |
| <GO:0001678> | intracellular glucose homeostasis                                                                                                                |  1.0602298 | 0.7891918 |
| <GO:0071385> | cellular response to glucocorticoid stimulus                                                                                                     |  1.0639645 | 0.7892019 |
| <GO:0060041> | retina development in camera-type eye                                                                                                            |  1.0636610 | 0.7892019 |
| <GO:0001892> | embryonic placenta development                                                                                                                   |  1.0632455 | 0.7892019 |
| <GO:0006334> | nucleosome assembly                                                                                                                              |  1.0601916 | 0.7892019 |
| <GO:0007605> | sensory perception of sound                                                                                                                      |  1.0562736 | 0.7892019 |
| <GO:0001913> | T cell mediated cytotoxicity                                                                                                                     |  1.0491282 | 0.7892019 |
| <GO:0002009> | morphogenesis of an epithelium                                                                                                                   |  1.0408985 | 0.7899687 |
| <GO:0050890> | cognition                                                                                                                                        |  1.0447738 | 0.7920905 |
| <GO:1903533> | regulation of protein targeting                                                                                                                  |  1.0544766 | 0.7926488 |
| <GO:0120193> | tight junction organization                                                                                                                      |  1.0629599 | 0.7934685 |
| <GO:1901657> | glycosyl compound metabolic process                                                                                                              |  1.0622175 | 0.7934685 |
| <GO:0048247> | lymphocyte chemotaxis                                                                                                                            |  1.0601292 | 0.7934685 |
| <GO:0006486> | protein glycosylation                                                                                                                            |  1.0566732 | 0.7934685 |
| <GO:0043413> | macromolecule glycosylation                                                                                                                      |  1.0566732 | 0.7934685 |
| <GO:0090382> | phagosome maturation                                                                                                                             |  1.0548770 | 0.7934685 |
| <GO:0034329> | cell junction assembly                                                                                                                           |  1.0384323 | 0.7934685 |
| <GO:0010648> | negative regulation of cell communication                                                                                                        |  1.0209563 | 0.7937612 |
| <GO:0001889> | liver development                                                                                                                                |  1.0584533 | 0.7944861 |
| <GO:1900181> | negative regulation of protein localization to nucleus                                                                                           |  1.0537657 | 0.7952535 |
| <GO:2001242> | regulation of intrinsic apoptotic signaling pathway                                                                                              |  1.0425651 | 0.7955372 |
| <GO:0046887> | positive regulation of hormone secretion                                                                                                         |  1.0531375 | 0.7965661 |
| <GO:0002279> | mast cell activation involved in immune response                                                                                                 |  1.0552786 | 0.7969141 |
| <GO:0051928> | positive regulation of calcium ion transport                                                                                                     |  1.0532779 | 0.7969141 |
| <GO:0090276> | regulation of peptide hormone secretion                                                                                                          |  1.0476684 | 0.7978370 |
| <GO:0030166> | proteoglycan biosynthetic process                                                                                                                |  1.0636213 | 0.7982676 |
| <GO:1903312> | negative regulation of mRNA metabolic process                                                                                                    |  1.0512036 | 0.7993861 |
| <GO:1904427> | positive regulation of calcium ion transmembrane transport                                                                                       | -1.0456988 | 0.8010250 |
| <GO:0010562> | positive regulation of phosphorus metabolic process                                                                                              |  1.0338597 | 0.8019802 |
| <GO:0045937> | positive regulation of phosphate metabolic process                                                                                               |  1.0338597 | 0.8019802 |
| <GO:0006090> | pyruvate metabolic process                                                                                                                       |  1.0392612 | 0.8040848 |
| <GO:1900182> | positive regulation of protein localization to nucleus                                                                                           | -1.0887956 | 0.8041396 |
| <GO:0001954> | positive regulation of cell-matrix adhesion                                                                                                      |  1.0523300 | 0.8047428 |
| <GO:0043648> | dicarboxylic acid metabolic process                                                                                                              |  1.0521042 | 0.8047428 |
| <GO:0023057> | negative regulation of signaling                                                                                                                 |  1.0190820 | 0.8047428 |
| <GO:0014015> | positive regulation of gliogenesis                                                                                                               |  1.0495448 | 0.8051265 |
| <GO:0003044> | regulation of systemic arterial blood pressure mediated by a chemical signal                                                                     |  1.0458864 | 0.8051265 |
| <GO:0050954> | sensory perception of mechanical stimulus                                                                                                        |  1.0443751 | 0.8051265 |
| <GO:0070266> | necroptotic process                                                                                                                              |  1.0422224 | 0.8051265 |
| <GO:0035924> | cellular response to vascular endothelial growth factor stimulus                                                                                 |  1.0421628 | 0.8051265 |
| <GO:0052372> | modulation by symbiont of entry into host                                                                                                        |  1.0411671 | 0.8051265 |
| <GO:0014902> | myotube differentiation                                                                                                                          |  1.0408764 | 0.8051265 |
| <GO:0055002> | striated muscle cell development                                                                                                                 |  1.0406716 | 0.8051265 |
| <GO:0051043> | regulation of membrane protein ectodomain proteolysis                                                                                            |  1.0372493 | 0.8051265 |
| <GO:0009261> | ribonucleotide catabolic process                                                                                                                 |  1.0364151 | 0.8051265 |
| <GO:0071214> | cellular response to abiotic stimulus                                                                                                            |  1.0320383 | 0.8051265 |
| <GO:0104004> | cellular response to environmental stimulus                                                                                                      |  1.0320383 | 0.8051265 |
| <GO:0070848> | response to growth factor                                                                                                                        |  1.0242200 | 0.8051265 |
| <GO:1902532> | negative regulation of intracellular signal transduction                                                                                         |  1.0217898 | 0.8051265 |
| <GO:1901379> | regulation of potassium ion transmembrane transport                                                                                              |  1.0519647 | 0.8062665 |
| <GO:0050806> | positive regulation of synaptic transmission                                                                                                     |  1.0477086 | 0.8062665 |
| <GO:0046364> | monosaccharide biosynthetic process                                                                                                              |  1.0424601 | 0.8062665 |
| <GO:0061008> | hepaticobiliary system development                                                                                                               |  1.0421928 | 0.8062665 |
| <GO:0008584> | male gonad development                                                                                                                           |  1.0400372 | 0.8062665 |
| <GO:0046546> | development of primary male sexual characteristics                                                                                               |  1.0400372 | 0.8062665 |
| <GO:2000108> | positive regulation of leukocyte apoptotic process                                                                                               |  1.0395859 | 0.8062665 |
| <GO:0046718> | symbiont entry into host cell                                                                                                                    |  1.0344142 | 0.8062665 |
| <GO:0022008> | neurogenesis                                                                                                                                     |  1.0150938 | 0.8062665 |
| <GO:0033674> | positive regulation of kinase activity                                                                                                           |  1.0333345 | 0.8107359 |
| <GO:0006936> | muscle contraction                                                                                                                               |  1.0332670 | 0.8107359 |
| <GO:0007626> | locomotory behavior                                                                                                                              |  1.0386624 | 0.8124647 |
| <GO:0031334> | positive regulation of protein-containing complex assembly                                                                                       |  1.0317071 | 0.8124647 |
| <GO:0007088> | regulation of mitotic nuclear division                                                                                                           |  1.0425757 | 0.8127952 |
| <GO:0002088> | lens development in camera-type eye                                                                                                              |  1.0371941 | 0.8127952 |
| <GO:1901264> | carbohydrate derivative transport                                                                                                                |  1.0358600 | 0.8155522 |
| <GO:0048568> | embryonic organ development                                                                                                                      |  1.0277794 | 0.8156979 |
| <GO:0016055> | Wnt signaling pathway                                                                                                                            |  1.0258644 | 0.8161432 |
| <GO:0007411> | axon guidance                                                                                                                                    |  1.0423599 | 0.8162741 |
| <GO:0097485> | neuron projection guidance                                                                                                                       |  1.0423599 | 0.8162741 |
| <GO:0003279> | cardiac septum development                                                                                                                       |  1.0400363 | 0.8163799 |
| <GO:0030323> | respiratory tube development                                                                                                                     |  1.0351087 | 0.8163799 |
| <GO:0006672> | ceramide metabolic process                                                                                                                       |  1.0398662 | 0.8166791 |
| <GO:0030218> | erythrocyte differentiation                                                                                                                      |  1.0285026 | 0.8166791 |
| <GO:0035794> | positive regulation of mitochondrial membrane permeability                                                                                       |  1.0376482 | 0.8169894 |
| <GO:0006749> | glutathione metabolic process                                                                                                                    |  1.0371048 | 0.8169894 |
| <GO:0001837> | epithelial to mesenchymal transition                                                                                                             |  1.0359099 | 0.8169894 |
| <GO:0006275> | regulation of DNA replication                                                                                                                    |  1.0326969 | 0.8169894 |
| <GO:0071363> | cellular response to growth factor stimulus                                                                                                      |  1.0218794 | 0.8172347 |
| <GO:0030902> | hindbrain development                                                                                                                            |  1.0359537 | 0.8178838 |
| <GO:0042326> | negative regulation of phosphorylation                                                                                                           |  1.0352189 | 0.8178838 |
| <GO:0009123> | nucleoside monophosphate metabolic process                                                                                                       |  1.0420889 | 0.8180531 |
| <GO:0051146> | striated muscle cell differentiation                                                                                                             |  1.0273446 | 0.8180531 |
| <GO:0030001> | metal ion transport                                                                                                                              |  1.0201055 | 0.8180531 |
| <GO:0007268> | chemical synaptic transmission                                                                                                                   |  1.0227638 | 0.8205382 |
| <GO:0098916> | anterograde trans-synaptic signaling                                                                                                             |  1.0227638 | 0.8205382 |
| <GO:0014009> | glial cell proliferation                                                                                                                         |  1.0322137 | 0.8206617 |
| <GO:0007517> | muscle organ development                                                                                                                         |  1.0292468 | 0.8206617 |
| <GO:0009100> | glycoprotein metabolic process                                                                                                                   |  1.0196127 | 0.8206617 |
| <GO:0021872> | forebrain generation of neurons                                                                                                                  |  1.0291527 | 0.8208540 |
| <GO:0010605> | negative regulation of macromolecule metabolic process                                                                                           |  1.0079262 | 0.8213415 |
| <GO:0043331> | response to dsRNA                                                                                                                                |  1.0313437 | 0.8214360 |
| <GO:0046031> | ADP metabolic process                                                                                                                            |  1.0278590 | 0.8226942 |
| <GO:0031400> | negative regulation of protein modification process                                                                                              |  1.0236590 | 0.8231075 |
| <GO:0044344> | cellular response to fibroblast growth factor stimulus                                                                                           |  1.0365707 | 0.8269716 |
| <GO:0071277> | cellular response to calcium ion                                                                                                                 |  1.0316077 | 0.8279584 |
| <GO:0001952> | regulation of cell-matrix adhesion                                                                                                               |  1.0278023 | 0.8302597 |
| <GO:0038061> | non-canonical NF-kappaB signal transduction                                                                                                      |  1.0262915 | 0.8306579 |
| <GO:0071229> | cellular response to acid chemical                                                                                                               |  1.0245946 | 0.8312358 |
| <GO:0034109> | homotypic cell-cell adhesion                                                                                                                     |  1.0231123 | 0.8312358 |
| <GO:0097345> | mitochondrial outer membrane permeabilization                                                                                                    |  1.0198727 | 0.8312358 |
| <GO:0008631> | intrinsic apoptotic signaling pathway in response to oxidative stress                                                                            |  1.0233021 | 0.8323689 |
| <GO:0045022> | early endosome to late endosome transport                                                                                                        | -1.0860098 | 0.8352989 |
| <GO:0008333> | endosome to lysosome transport                                                                                                                   |  1.0276372 | 0.8352989 |
| <GO:0032024> | positive regulation of insulin secretion                                                                                                         |  1.0266357 | 0.8352989 |
| <GO:0006821> | chloride transport                                                                                                                               |  1.0247254 | 0.8352989 |
| <GO:0060828> | regulation of canonical Wnt signaling pathway                                                                                                    |  1.0234814 | 0.8352989 |
| <GO:0001649> | osteoblast differentiation                                                                                                                       |  1.0166941 | 0.8352989 |
| <GO:0098926> | postsynaptic signal transduction                                                                                                                 |  1.0165165 | 0.8352989 |
| <GO:0051348> | negative regulation of transferase activity                                                                                                      |  1.0210368 | 0.8369552 |
| <GO:0009154> | purine ribonucleotide catabolic process                                                                                                          |  1.0162191 | 0.8377229 |
| <GO:0045921> | positive regulation of exocytosis                                                                                                                |  1.0245504 | 0.8385755 |
| <GO:0042303> | molting cycle                                                                                                                                    |  1.0203304 | 0.8385755 |
| <GO:0042633> | hair cycle                                                                                                                                       |  1.0203304 | 0.8385755 |
| <GO:0016447> | somatic recombination of immunoglobulin gene segments                                                                                            |  1.0183737 | 0.8385755 |
| <GO:0009134> | nucleoside diphosphate catabolic process                                                                                                         |  1.0161274 | 0.8385755 |
| <GO:0009137> | purine nucleoside diphosphate catabolic process                                                                                                  |  1.0161274 | 0.8385755 |
| <GO:0009181> | purine ribonucleoside diphosphate catabolic process                                                                                              |  1.0161274 | 0.8385755 |
| <GO:0009191> | ribonucleoside diphosphate catabolic process                                                                                                     |  1.0161274 | 0.8385755 |
| <GO:0046032> | ADP catabolic process                                                                                                                            |  1.0161274 | 0.8385755 |
| <GO:0007611> | learning or memory                                                                                                                               |  1.0143353 | 0.8385755 |
| <GO:0015914> | phospholipid transport                                                                                                                           |  1.0139317 | 0.8385755 |
| <GO:0009792> | embryo development ending in birth or egg hatching                                                                                               |  1.0126442 | 0.8385755 |
| <GO:0033108> | mitochondrial respiratory chain complex assembly                                                                                                 | -0.9332500 | 0.8385755 |
| <GO:0090090> | negative regulation of canonical Wnt signaling pathway                                                                                           |  1.0202466 | 0.8390217 |
| <GO:0043524> | negative regulation of neuron apoptotic process                                                                                                  |  1.0158049 | 0.8392125 |
| <GO:0034248> | regulation of amide metabolic process                                                                                                            |  1.0109092 | 0.8392125 |
| <GO:0140894> | endolysosomal toll-like receptor signaling pathway                                                                                               |  1.0188149 | 0.8398800 |
| <GO:0061337> | cardiac conduction                                                                                                                               |  1.0157810 | 0.8398800 |
| <GO:0006516> | glycoprotein catabolic process                                                                                                                   |  1.0146362 | 0.8398800 |
| <GO:1905144> | response to acetylcholine                                                                                                                        |  1.0099954 | 0.8398800 |
| <GO:0022414> | reproductive process                                                                                                                             |  1.0044764 | 0.8398800 |
| <GO:0001919> | regulation of receptor recycling                                                                                                                 | -0.9878918 | 0.8398800 |
| <GO:0006656> | phosphatidylcholine biosynthetic process                                                                                                         | -0.9842010 | 0.8398800 |
| <GO:0050767> | regulation of neurogenesis                                                                                                                       |  1.0177989 | 0.8401660 |
| <GO:0051604> | protein maturation                                                                                                                               |  1.0126013 | 0.8405143 |
| <GO:0090559> | regulation of membrane permeability                                                                                                              |  1.0184269 | 0.8405348 |
| <GO:0001935> | endothelial cell proliferation                                                                                                                   |  1.0166638 | 0.8405348 |
| <GO:0000271> | polysaccharide biosynthetic process                                                                                                              |  1.0128697 | 0.8405348 |
| <GO:0008543> | fibroblast growth factor receptor signaling pathway                                                                                              |  1.0126757 | 0.8405348 |
| <GO:0033209> | tumor necrosis factor-mediated signaling pathway                                                                                                 |  1.0119124 | 0.8405348 |
| <GO:0046850> | regulation of bone remodeling                                                                                                                    |  1.0115673 | 0.8405348 |
| <GO:0007492> | endoderm development                                                                                                                             |  1.0075259 | 0.8405348 |
| <GO:0043434> | response to peptide hormone                                                                                                                      |  1.0094789 | 0.8407535 |
| <GO:0051701> | biological process involved in interaction with host                                                                                             |  1.0087267 | 0.8407535 |
| <GO:1902410> | mitotic cytokinetic process                                                                                                                      | -1.0765936 | 0.8416167 |
| <GO:0019319> | hexose biosynthetic process                                                                                                                      |  1.0188821 | 0.8416167 |
| <GO:0071824> | protein-DNA complex organization                                                                                                                 |  1.0122203 | 0.8417639 |
| <GO:0043500> | muscle adaptation                                                                                                                                |  1.0121732 | 0.8417639 |
| <GO:0048675> | axon extension                                                                                                                                   |  1.0116783 | 0.8417639 |
| <GO:0062014> | negative regulation of small molecule metabolic process                                                                                          |  1.0114815 | 0.8417639 |
| <GO:0015748> | organophosphate ester transport                                                                                                                  |  1.0088034 | 0.8417639 |
| <GO:0002244> | hematopoietic progenitor cell differentiation                                                                                                    |  1.0095090 | 0.8430760 |
| <GO:0002040> | sprouting angiogenesis                                                                                                                           |  1.0136899 | 0.8440541 |
| <GO:0002704> | negative regulation of leukocyte mediated immunity                                                                                               |  1.0121244 | 0.8440541 |
| <GO:0072523> | purine-containing compound catabolic process                                                                                                     |  1.0102723 | 0.8440541 |
| <GO:0010811> | positive regulation of cell-substrate adhesion                                                                                                   |  1.0095957 | 0.8440541 |
| <GO:0001782> | B cell homeostasis                                                                                                                               |  1.0040955 | 0.8440541 |
| <GO:0051653> | spindle localization                                                                                                                             | -0.9927371 | 0.8440541 |
| <GO:0009968> | negative regulation of signal transduction                                                                                                       |  1.0002607 | 0.8445721 |
| <GO:0007612> | learning                                                                                                                                         |  1.0118230 | 0.8452066 |
| <GO:0046825> | regulation of protein export from nucleus                                                                                                        |  1.0068456 | 0.8452066 |
| <GO:0008154> | actin polymerization or depolymerization                                                                                                         |  1.0050309 | 0.8452066 |
| <GO:1903421> | regulation of synaptic vesicle recycling                                                                                                         |  1.0014196 | 0.8452066 |
| <GO:0010594> | regulation of endothelial cell migration                                                                                                         |  1.0074560 | 0.8461916 |
| <GO:0055074> | calcium ion homeostasis                                                                                                                          |  1.0044237 | 0.8461916 |
| <GO:0010631> | epithelial cell migration                                                                                                                        |  1.0032151 | 0.8461916 |
| <GO:0090130> | tissue migration                                                                                                                                 |  1.0032151 | 0.8461916 |
| <GO:0090132> | epithelium migration                                                                                                                             |  1.0032151 | 0.8461916 |
| <GO:0016052> | carbohydrate catabolic process                                                                                                                   |  1.0028735 | 0.8461916 |
| <GO:0006368> | transcription elongation by RNA polymerase II                                                                                                    | -0.8691621 | 0.8461916 |
| <GO:0050435> | amyloid-beta metabolic process                                                                                                                   |  1.0063644 | 0.8464674 |
| <GO:1902774> | late endosome to lysosome transport                                                                                                              |  1.0050298 | 0.8475464 |
| <GO:0060711> | labyrinthine layer development                                                                                                                   |  0.9991981 | 0.8478992 |
| <GO:0042743> | hydrogen peroxide metabolic process                                                                                                              |  1.0084953 | 0.8482727 |
| <GO:0051205> | protein insertion into membrane                                                                                                                  | -1.0834904 | 0.8483837 |
| <GO:0090068> | positive regulation of cell cycle process                                                                                                        |  1.0014547 | 0.8483837 |
| <GO:0035456> | response to interferon-beta                                                                                                                      |  0.9963814 | 0.8483837 |
| <GO:0061572> | actin filament bundle organization                                                                                                               |  1.0013575 | 0.8486884 |
| <GO:0045980> | negative regulation of nucleotide metabolic process                                                                                              |  0.9942461 | 0.8489724 |
| <GO:1900543> | negative regulation of purine nucleotide metabolic process                                                                                       |  0.9942461 | 0.8489724 |
| <GO:1903579> | negative regulation of ATP metabolic process                                                                                                     |  0.9942461 | 0.8489724 |
| <GO:0034101> | erythrocyte homeostasis                                                                                                                          |  1.0000143 | 0.8507473 |
| <GO:0050804> | modulation of chemical synaptic transmission                                                                                                     |  1.0016483 | 0.8513362 |
| <GO:0099177> | regulation of trans-synaptic signaling                                                                                                           |  1.0016483 | 0.8513362 |
| <GO:0016050> | vesicle organization                                                                                                                             |  1.0074965 | 0.8520857 |
| <GO:0010632> | regulation of epithelial cell migration                                                                                                          |  1.0026722 | 0.8520857 |
| <GO:0097178> | ruffle assembly                                                                                                                                  |  1.0018663 | 0.8528430 |
| <GO:0006195> | purine nucleotide catabolic process                                                                                                              |  0.9984489 | 0.8528430 |
| <GO:0001938> | positive regulation of endothelial cell proliferation                                                                                            |  0.9922275 | 0.8528430 |
| <GO:0044272> | sulfur compound biosynthetic process                                                                                                             |  0.9963769 | 0.8532562 |
| <GO:0055001> | muscle cell development                                                                                                                          |  1.0105941 | 0.8534138 |
| <GO:1901224> | positive regulation of non-canonical NF-kappaB signal transduction                                                                               |  0.9996896 | 0.8534138 |
| <GO:1903146> | regulation of autophagy of mitochondrion                                                                                                         |  0.9982214 | 0.8534138 |
| <GO:0006816> | calcium ion transport                                                                                                                            |  1.0017719 | 0.8537757 |
| <GO:0032869> | cellular response to insulin stimulus                                                                                                            |  0.9971854 | 0.8537757 |
| <GO:0050821> | protein stabilization                                                                                                                            |  0.9983279 | 0.8549436 |
| <GO:0051048> | negative regulation of secretion                                                                                                                 |  0.9969826 | 0.8552220 |
| <GO:0002381> | immunoglobulin production involved in immunoglobulin-mediated immune response                                                                    |  0.9932276 | 0.8553515 |
| <GO:0001990> | regulation of systemic arterial blood pressure by hormone                                                                                        |  0.9905508 | 0.8578990 |
| <GO:0003081> | regulation of systemic arterial blood pressure by renin-angiotensin                                                                              |  0.9905508 | 0.8578990 |
| <GO:0043534> | blood vessel endothelial cell migration                                                                                                          |  1.0046227 | 0.8582130 |
| <GO:0030324> | lung development                                                                                                                                 |  0.9957236 | 0.8617846 |
| <GO:0071774> | response to fibroblast growth factor                                                                                                             |  0.9998246 | 0.8618518 |
| <GO:0043502> | regulation of muscle adaptation                                                                                                                  |  0.9992108 | 0.8618518 |
| <GO:0006096> | glycolytic process                                                                                                                               |  0.9946220 | 0.8626726 |
| <GO:0019364> | pyridine nucleotide catabolic process                                                                                                            |  0.9946220 | 0.8626726 |
| <GO:0055085> | transmembrane transport                                                                                                                          |  0.9942789 | 0.8626726 |
| <GO:0046467> | membrane lipid biosynthetic process                                                                                                              |  0.9935577 | 0.8626726 |
| <GO:0062012> | regulation of small molecule metabolic process                                                                                                   |  0.9951073 | 0.8626821 |
| <GO:0007173> | epidermal growth factor receptor signaling pathway                                                                                               |  0.9946514 | 0.8626821 |
| <GO:0030148> | sphingolipid biosynthetic process                                                                                                                |  0.9923165 | 0.8626821 |
| <GO:0070050> | neuron cellular homeostasis                                                                                                                      |  0.9918074 | 0.8626821 |
| <GO:0009166> | nucleotide catabolic process                                                                                                                     |  0.9912603 | 0.8626821 |
| <GO:0007041> | lysosomal transport                                                                                                                              |  0.9911794 | 0.8626821 |
| <GO:1905710> | positive regulation of membrane permeability                                                                                                     |  0.9984638 | 0.8629025 |
| <GO:2000144> | positive regulation of DNA-templated transcription initiation                                                                                    |  0.9918245 | 0.8629025 |
| <GO:0003012> | muscle system process                                                                                                                            |  0.9967533 | 0.8651586 |
| <GO:0005996> | monosaccharide metabolic process                                                                                                                 |  0.9950442 | 0.8655766 |
| <GO:0043087> | regulation of GTPase activity                                                                                                                    |  0.9946131 | 0.8655766 |
| <GO:0140374> | antiviral innate immune response                                                                                                                 |  0.9836087 | 0.8655766 |
| <GO:0002200> | somatic diversification of immune receptors                                                                                                      |  0.9920747 | 0.8660090 |
| <GO:0006970> | response to osmotic stress                                                                                                                       |  0.9915622 | 0.8660090 |
| <GO:0043009> | chordate embryonic development                                                                                                                   |  0.9914831 | 0.8660090 |
| <GO:0097009> | energy homeostasis                                                                                                                               |  0.9901129 | 0.8660090 |
| <GO:1902041> | regulation of extrinsic apoptotic signaling pathway via death domain receptors                                                                   |  0.9897657 | 0.8660090 |
| <GO:0006006> | glucose metabolic process                                                                                                                        |  0.9879743 | 0.8660090 |
| <GO:0050994> | regulation of lipid catabolic process                                                                                                            |  0.9865152 | 0.8660090 |
| <GO:0030901> | midbrain development                                                                                                                             |  0.9841087 | 0.8660090 |
| <GO:0008286> | insulin receptor signaling pathway                                                                                                               |  0.9837311 | 0.8660090 |
| <GO:0150117> | positive regulation of cell-substrate junction organization                                                                                      |  0.9829638 | 0.8660090 |
| <GO:0001942> | hair follicle development                                                                                                                        |  0.9803431 | 0.8660090 |
| <GO:0022404> | molting cycle process                                                                                                                            |  0.9803431 | 0.8660090 |
| <GO:0022405> | hair cycle process                                                                                                                               |  0.9803431 | 0.8660090 |
| <GO:0032526> | response to retinoic acid                                                                                                                        |  0.9799760 | 0.8660090 |
| <GO:0006278> | RNA-templated DNA biosynthetic process                                                                                                           | -0.9620565 | 0.8660090 |
| <GO:0009994> | oocyte differentiation                                                                                                                           | -0.9620565 | 0.8660090 |
| <GO:0032786> | positive regulation of DNA-templated transcription, elongation                                                                                   | -0.9620565 | 0.8660090 |
| <GO:0048599> | oocyte development                                                                                                                               | -0.9620565 | 0.8660090 |
| <GO:0009790> | embryo development                                                                                                                               |  0.9916557 | 0.8670212 |
| <GO:0034766> | negative regulation of monoatomic ion transmembrane transport                                                                                    |  0.9866837 | 0.8670212 |
| <GO:1904063> | negative regulation of cation transmembrane transport                                                                                            |  0.9866837 | 0.8670212 |
| <GO:0006007> | glucose catabolic process                                                                                                                        |  0.9784965 | 0.8670212 |
| <GO:1902882> | regulation of response to oxidative stress                                                                                                       |  0.9754358 | 0.8670212 |
| <GO:0045861> | negative regulation of proteolysis                                                                                                               |  0.9833349 | 0.8670707 |
| <GO:0030308> | negative regulation of cell growth                                                                                                               |  0.9832554 | 0.8670707 |
| <GO:0051495> | positive regulation of cytoskeleton organization                                                                                                 |  0.9788082 | 0.8670707 |
| <GO:0001910> | regulation of leukocyte mediated cytotoxicity                                                                                                    |  0.9890496 | 0.8671008 |
| <GO:0009749> | response to glucose                                                                                                                              |  0.9939738 | 0.8677447 |
| <GO:0071383> | cellular response to steroid hormone stimulus                                                                                                    |  0.9820816 | 0.8677447 |
| <GO:0031341> | regulation of cell killing                                                                                                                       |  0.9802852 | 0.8685546 |
| <GO:1901698> | response to nitrogen compound                                                                                                                    |  0.9903255 | 0.8694161 |
| <GO:0007610> | behavior                                                                                                                                         |  0.9870758 | 0.8694161 |
| <GO:0090277> | positive regulation of peptide hormone secretion                                                                                                 |  0.9826760 | 0.8694161 |
| <GO:1901292> | nucleoside phosphate catabolic process                                                                                                           |  0.9821481 | 0.8694161 |
| <GO:0030201> | heparan sulfate proteoglycan metabolic process                                                                                                   |  0.9771532 | 0.8694161 |
| <GO:1902476> | chloride transmembrane transport                                                                                                                 |  0.9710480 | 0.8720695 |
| <GO:1901989> | positive regulation of cell cycle phase transition                                                                                               |  0.9798236 | 0.8726566 |
| <GO:1901992> | positive regulation of mitotic cell cycle phase transition                                                                                       |  0.9798236 | 0.8726566 |
| <GO:0045661> | regulation of myoblast differentiation                                                                                                           |  0.9747878 | 0.8727489 |
| <GO:0051592> | response to calcium ion                                                                                                                          |  0.9844191 | 0.8728227 |
| <GO:0009952> | anterior/posterior pattern specification                                                                                                         |  0.9766099 | 0.8728227 |
| <GO:1901222> | regulation of non-canonical NF-kappaB signal transduction                                                                                        |  0.9880847 | 0.8730941 |
| <GO:0042692> | muscle cell differentiation                                                                                                                      |  0.9840280 | 0.8730941 |
| <GO:0035282> | segmentation                                                                                                                                     |  0.9703157 | 0.8733951 |
| <GO:0002707> | negative regulation of lymphocyte mediated immunity                                                                                              |  0.9789815 | 0.8738793 |
| <GO:0099003> | vesicle-mediated transport in synapse                                                                                                            |  0.9851278 | 0.8747694 |
| <GO:0007409> | axonogenesis                                                                                                                                     |  0.9826624 | 0.8747694 |
| <GO:0002377> | immunoglobulin production                                                                                                                        |  0.9734804 | 0.8747694 |
| <GO:0060135> | maternal process involved in female pregnancy                                                                                                    |  0.9734045 | 0.8747694 |
| <GO:0051098> | regulation of binding                                                                                                                            |  0.9729048 | 0.8747694 |
| <GO:0006820> | monoatomic anion transport                                                                                                                       |  0.9663590 | 0.8747694 |
| <GO:0051298> | centrosome duplication                                                                                                                           | -0.9544531 | 0.8747694 |
| <GO:0001558> | regulation of cell growth                                                                                                                        |  0.9888719 | 0.8763475 |
| <GO:0045934> | negative regulation of nucleobase-containing compound metabolic process                                                                          |  0.9874454 | 0.8763475 |
| <GO:0007272> | ensheathment of neurons                                                                                                                          |  0.9776429 | 0.8763475 |
| <GO:0008366> | axon ensheathment                                                                                                                                |  0.9776429 | 0.8763475 |
| <GO:0042552> | myelination                                                                                                                                      |  0.9776429 | 0.8763475 |
| <GO:0009755> | hormone-mediated signaling pathway                                                                                                               |  0.9726384 | 0.8763475 |
| <GO:0033143> | regulation of intracellular steroid hormone receptor signaling pathway                                                                           | -1.0529921 | 0.8767539 |
| <GO:0030042> | actin filament depolymerization                                                                                                                  |  0.9640702 | 0.8778399 |
| <GO:1904950> | negative regulation of establishment of protein localization                                                                                     |  0.9777852 | 0.8791714 |
| <GO:0038127> | ERBB signaling pathway                                                                                                                           |  0.9744298 | 0.8791714 |
| <GO:0090596> | sensory organ morphogenesis                                                                                                                      |  0.9713039 | 0.8791714 |
| <GO:0072091> | regulation of stem cell proliferation                                                                                                            |  0.9641996 | 0.8791714 |
| <GO:0006665> | sphingolipid metabolic process                                                                                                                   |  0.9659605 | 0.8814715 |
| <GO:0061061> | muscle structure development                                                                                                                     |  0.9860660 | 0.8818274 |
| <GO:0010563> | negative regulation of phosphorus metabolic process                                                                                              |  0.9715372 | 0.8818274 |
| <GO:0045936> | negative regulation of phosphate metabolic process                                                                                               |  0.9715372 | 0.8818274 |
| <GO:0046822> | regulation of nucleocytoplasmic transport                                                                                                        |  0.9705523 | 0.8818274 |
| <GO:0031399> | regulation of protein modification process                                                                                                       |  0.9851113 | 0.8820000 |
| <GO:0031589> | cell-substrate adhesion                                                                                                                          |  0.9828609 | 0.8820000 |
| <GO:0006790> | sulfur compound metabolic process                                                                                                                |  0.9729359 | 0.8820000 |
| <GO:0050851> | antigen receptor-mediated signaling pathway                                                                                                      |  0.9679480 | 0.8820000 |
| <GO:0051402> | neuron apoptotic process                                                                                                                         |  0.9732810 | 0.8833694 |
| <GO:0001936> | regulation of endothelial cell proliferation                                                                                                     |  0.9696614 | 0.8833694 |
| <GO:0072526> | pyridine-containing compound catabolic process                                                                                                   |  0.9684213 | 0.8833694 |
| <GO:0045185> | maintenance of protein location                                                                                                                  |  0.9558091 | 0.8833694 |
| <GO:0098930> | axonal transport                                                                                                                                 | -1.0786682 | 0.8837815 |
| <GO:0045047> | protein targeting to ER                                                                                                                          | -1.0411248 | 0.8837815 |
| <GO:0062013> | positive regulation of small molecule metabolic process                                                                                          |  0.9656087 | 0.8837815 |
| <GO:0090199> | regulation of release of cytochrome c from mitochondria                                                                                          |  0.9583627 | 0.8837815 |
| <GO:0006110> | regulation of glycolytic process                                                                                                                 |  0.9507782 | 0.8837815 |
| <GO:0030811> | regulation of nucleotide catabolic process                                                                                                       |  0.9507782 | 0.8837815 |
| <GO:0033121> | regulation of purine nucleotide catabolic process                                                                                                |  0.9507782 | 0.8837815 |
| <GO:0071805> | potassium ion transmembrane transport                                                                                                            |  0.9684788 | 0.8838111 |
| <GO:0140112> | extracellular vesicle biogenesis                                                                                                                 |  0.9555518 | 0.8838111 |
| <GO:0019080> | viral gene expression                                                                                                                            | -0.8081989 | 0.8838111 |
| <GO:2000779> | regulation of double-strand break repair                                                                                                         | -1.0492390 | 0.8857079 |
| <GO:0006893> | Golgi to plasma membrane transport                                                                                                               |  0.9683166 | 0.8857538 |
| <GO:0045926> | negative regulation of growth                                                                                                                    |  0.9629132 | 0.8858544 |
| <GO:0009890> | negative regulation of biosynthetic process                                                                                                      |  0.9896632 | 0.8859645 |
| <GO:0048167> | regulation of synaptic plasticity                                                                                                                |  0.9669544 | 0.8859645 |
| <GO:0019377> | glycolipid catabolic process                                                                                                                     |  0.9619144 | 0.8859645 |
| <GO:0045913> | positive regulation of carbohydrate metabolic process                                                                                            |  0.9615356 | 0.8859645 |
| <GO:0045920> | negative regulation of exocytosis                                                                                                                |  0.9553904 | 0.8871800 |
| <GO:0003407> | neural retina development                                                                                                                        |  0.9510485 | 0.8889289 |
| <GO:0006364> | rRNA processing                                                                                                                                  | -1.0287308 | 0.8896797 |
| <GO:0051253> | negative regulation of RNA metabolic process                                                                                                     |  0.9847711 | 0.8896797 |
| <GO:0018193> | peptidyl-amino acid modification                                                                                                                 |  0.9748908 | 0.8896797 |
| <GO:0019318> | hexose metabolic process                                                                                                                         |  0.9678625 | 0.8896797 |
| <GO:0051302> | regulation of cell division                                                                                                                      |  0.9570658 | 0.8896797 |
| <GO:0048562> | embryonic organ morphogenesis                                                                                                                    |  0.9556904 | 0.8896797 |
| <GO:0051646> | mitochondrion localization                                                                                                                       |  0.9504843 | 0.8896797 |
| <GO:0010761> | fibroblast migration                                                                                                                             |  0.9499361 | 0.8896797 |
| <GO:0048261> | negative regulation of receptor-mediated endocytosis                                                                                             |  0.9471011 | 0.8896797 |
| <GO:0071985> | multivesicular body sorting pathway                                                                                                              |  0.9603794 | 0.8898153 |
| <GO:0051128> | regulation of cellular component organization                                                                                                    |  0.9872464 | 0.8913824 |
| <GO:0071375> | cellular response to peptide hormone stimulus                                                                                                    |  0.9664320 | 0.8913824 |
| <GO:0045017> | glycerolipid biosynthetic process                                                                                                                |  0.9643259 | 0.8913824 |
| <GO:0030516> | regulation of axon extension                                                                                                                     |  0.9538456 | 0.8913824 |
| <GO:1904705> | regulation of vascular associated smooth muscle cell proliferation                                                                               |  0.9478030 | 0.8913824 |
| <GO:1990874> | vascular associated smooth muscle cell proliferation                                                                                             |  0.9478030 | 0.8913824 |
| <GO:0006414> | translational elongation                                                                                                                         | -0.9258582 | 0.8913824 |
| <GO:2001234> | negative regulation of apoptotic signaling pathway                                                                                               |  0.9626890 | 0.8917240 |
| <GO:0006367> | transcription initiation at RNA polymerase II promoter                                                                                           |  0.9544530 | 0.8917240 |
| <GO:0048857> | neural nucleus development                                                                                                                       |  0.9519693 | 0.8917687 |
| <GO:0048592> | eye morphogenesis                                                                                                                                |  0.9421384 | 0.8917687 |
| <GO:0009746> | response to hexose                                                                                                                               |  0.9604307 | 0.8919432 |
| <GO:0006937> | regulation of muscle contraction                                                                                                                 |  0.9475814 | 0.8919432 |
| <GO:0090287> | regulation of cellular response to growth factor stimulus                                                                                        |  0.9597948 | 0.8956857 |
| <GO:0010927> | cellular component assembly involved in morphogenesis                                                                                            |  0.9495480 | 0.8970328 |
| <GO:0032989> | cellular anatomical entity morphogenesis                                                                                                         |  0.9495480 | 0.8970328 |
| <GO:0044839> | cell cycle G2/M phase transition                                                                                                                 |  0.9345109 | 0.8996032 |
| <GO:0043523> | regulation of neuron apoptotic process                                                                                                           |  0.9569928 | 0.8996940 |
| <GO:0070527> | platelet aggregation                                                                                                                             |  0.9495909 | 0.8996940 |
| <GO:0051125> | regulation of actin nucleation                                                                                                                   | -0.9382254 | 0.8996940 |
| <GO:0006612> | protein targeting to membrane                                                                                                                    |  0.9581960 | 0.9011644 |
| <GO:0010803> | regulation of tumor necrosis factor-mediated signaling pathway                                                                                   |  0.9415571 | 0.9011967 |
| <GO:0007264> | small GTPase-mediated signal transduction                                                                                                        |  0.9747871 | 0.9024609 |
| <GO:0016239> | positive regulation of macroautophagy                                                                                                            |  0.9472141 | 0.9024609 |
| <GO:0051248> | negative regulation of protein metabolic process                                                                                                 |  0.9757131 | 0.9036642 |
| <GO:0016310> | phosphorylation                                                                                                                                  |  0.9730815 | 0.9036642 |
| <GO:0030522> | intracellular receptor signaling pathway                                                                                                         |  0.9698855 | 0.9036642 |
| <GO:0022411> | cellular component disassembly                                                                                                                   |  0.9621623 | 0.9036642 |
| <GO:0021537> | telencephalon development                                                                                                                        |  0.9521286 | 0.9036642 |
| <GO:0006874> | intracellular calcium ion homeostasis                                                                                                            |  0.9515333 | 0.9036642 |
| <GO:1905477> | positive regulation of protein localization to membrane                                                                                          |  0.9474334 | 0.9036642 |
| <GO:0090398> | cellular senescence                                                                                                                              |  0.9458502 | 0.9036642 |
| <GO:2001243> | negative regulation of intrinsic apoptotic signaling pathway                                                                                     |  0.9452581 | 0.9036642 |
| <GO:0048863> | stem cell differentiation                                                                                                                        |  0.9445779 | 0.9036642 |
| <GO:0071326> | cellular response to monosaccharide stimulus                                                                                                     |  0.9418733 | 0.9036642 |
| <GO:0034764> | positive regulation of transmembrane transport                                                                                                   |  0.9415523 | 0.9036642 |
| <GO:1902003> | regulation of amyloid-beta formation                                                                                                             |  0.9403282 | 0.9036642 |
| <GO:0048821> | erythrocyte development                                                                                                                          |  0.9392745 | 0.9036642 |
| <GO:0071479> | cellular response to ionizing radiation                                                                                                          |  0.9391755 | 0.9036642 |
| <GO:0048278> | vesicle docking                                                                                                                                  |  0.9384569 | 0.9036642 |
| <GO:0070585> | protein localization to mitochondrion                                                                                                            |  0.9349983 | 0.9036642 |
| <GO:0018105> | peptidyl-serine phosphorylation                                                                                                                  |  0.9346574 | 0.9036642 |
| <GO:0018209> | peptidyl-serine modification                                                                                                                     |  0.9346574 | 0.9036642 |
| <GO:0060291> | long-term synaptic potentiation                                                                                                                  |  0.9324337 | 0.9036642 |
| <GO:0048489> | synaptic vesicle transport                                                                                                                       |  0.9282216 | 0.9036642 |
| <GO:0006906> | vesicle fusion                                                                                                                                   |  0.9461179 | 0.9042258 |
| <GO:0034405> | response to fluid shear stress                                                                                                                   |  0.9336131 | 0.9042258 |
| <GO:0007178> | cell surface receptor protein serine/threonine kinase signaling pathway                                                                          |  0.9508304 | 0.9048427 |
| <GO:0061025> | membrane fusion                                                                                                                                  |  0.9463016 | 0.9048427 |
| <GO:0006094> | gluconeogenesis                                                                                                                                  |  0.9450458 | 0.9048427 |
| <GO:0032355> | response to estradiol                                                                                                                            |  0.9404162 | 0.9048427 |
| <GO:0034394> | protein localization to cell surface                                                                                                             |  0.9387703 | 0.9048427 |
| <GO:2000785> | regulation of autophagosome assembly                                                                                                             |  0.9383375 | 0.9048427 |
| <GO:0098656> | monoatomic anion transmembrane transport                                                                                                         |  0.9371450 | 0.9048427 |
| <GO:0030041> | actin filament polymerization                                                                                                                    |  0.9354510 | 0.9048427 |
| <GO:0043535> | regulation of blood vessel endothelial cell migration                                                                                            |  0.9337132 | 0.9048427 |
| <GO:1902991> | regulation of amyloid precursor protein catabolic process                                                                                        |  0.9322986 | 0.9048427 |
| <GO:0007080> | mitotic metaphase chromosome alignment                                                                                                           | -0.9251489 | 0.9048427 |
| <GO:0010824> | regulation of centrosome duplication                                                                                                             | -0.9251489 | 0.9048427 |
| <GO:0045787> | positive regulation of cell cycle                                                                                                                |  0.9420101 | 0.9049645 |
| <GO:0006360> | transcription by RNA polymerase I                                                                                                                | -0.9296297 | 0.9049645 |
| <GO:0032387> | negative regulation of intracellular transport                                                                                                   |  0.9280275 | 0.9049645 |
| <GO:0031113> | regulation of microtubule polymerization                                                                                                         |  0.9235750 | 0.9049645 |
| <GO:0035459> | vesicle cargo loading                                                                                                                            |  0.9229280 | 0.9049645 |
| <GO:0098659> | inorganic cation import across plasma membrane                                                                                                   |  0.9308064 | 0.9056647 |
| <GO:0099587> | inorganic ion import across plasma membrane                                                                                                      |  0.9308064 | 0.9056647 |
| <GO:0010975> | regulation of neuron projection development                                                                                                      |  0.9577491 | 0.9059112 |
| <GO:0032870> | cellular response to hormone stimulus                                                                                                            |  0.9665623 | 0.9063786 |
| <GO:0007007> | inner mitochondrial membrane organization                                                                                                        | -0.8865610 | 0.9063786 |
| <GO:0006811> | monoatomic ion transport                                                                                                                         |  0.9707483 | 0.9068295 |
| <GO:0090066> | regulation of anatomical structure size                                                                                                          |  0.9554782 | 0.9068295 |
| <GO:0006839> | mitochondrial transport                                                                                                                          |  0.9339576 | 0.9068295 |
| <GO:0010821> | regulation of mitochondrion organization                                                                                                         |  0.9324950 | 0.9068295 |
| <GO:0071333> | cellular response to glucose stimulus                                                                                                            |  0.9299208 | 0.9068295 |
| <GO:0032509> | endosome transport via multivesicular body sorting pathway                                                                                       |  0.9241134 | 0.9068295 |
| <GO:2000629> | negative regulation of miRNA metabolic process                                                                                                   |  0.9234117 | 0.9068295 |
| <GO:0048168> | regulation of neuronal synaptic plasticity                                                                                                       |  0.9219657 | 0.9068295 |
| <GO:0022037> | metencephalon development                                                                                                                        |  0.9211396 | 0.9068295 |
| <GO:0032436> | positive regulation of proteasomal ubiquitin-dependent protein catabolic process                                                                 |  0.9202887 | 0.9068295 |
| <GO:0016601> | Rac protein signal transduction                                                                                                                  |  0.9199104 | 0.9068295 |
| <GO:0010810> | regulation of cell-substrate adhesion                                                                                                            |  0.9394971 | 0.9089029 |
| <GO:0032418> | lysosome localization                                                                                                                            |  0.9389052 | 0.9089029 |
| <GO:1990849> | vacuolar localization                                                                                                                            |  0.9389052 | 0.9089029 |
| <GO:0010657> | muscle cell apoptotic process                                                                                                                    |  0.9321058 | 0.9089029 |
| <GO:0061951> | establishment of protein localization to plasma membrane                                                                                         |  0.9171444 | 0.9089029 |
| <GO:0010658> | striated muscle cell apoptotic process                                                                                                           |  0.9167421 | 0.9089029 |
| <GO:0010659> | cardiac muscle cell apoptotic process                                                                                                            |  0.9167421 | 0.9089029 |
| <GO:0043266> | regulation of potassium ion transport                                                                                                            |  0.9301818 | 0.9110750 |
| <GO:0031110> | regulation of microtubule polymerization or depolymerization                                                                                     |  0.9300464 | 0.9110750 |
| <GO:0035023> | regulation of Rho protein signal transduction                                                                                                    |  0.9281354 | 0.9110750 |
| <GO:0060324> | face development                                                                                                                                 |  0.9148158 | 0.9110750 |
| <GO:0030031> | cell projection assembly                                                                                                                         |  0.9507028 | 0.9123445 |
| <GO:1904377> | positive regulation of protein localization to cell periphery                                                                                    |  0.9271850 | 0.9123445 |
| <GO:0042789> | mRNA transcription by RNA polymerase II                                                                                                          |  0.9224782 | 0.9160737 |
| <GO:0045103> | intermediate filament-based process                                                                                                              | -0.8821075 | 0.9160737 |
| <GO:0045104> | intermediate filament cytoskeleton organization                                                                                                  | -0.8821075 | 0.9160737 |
| <GO:1901135> | carbohydrate derivative metabolic process                                                                                                        |  0.9679712 | 0.9173255 |
| <GO:0006399> | tRNA metabolic process                                                                                                                           | -0.9663062 | 0.9173255 |
| <GO:0050852> | T cell receptor signaling pathway                                                                                                                |  0.9288975 | 0.9173255 |
| <GO:1903509> | liposaccharide metabolic process                                                                                                                 |  0.9280880 | 0.9173255 |
| <GO:1903531> | negative regulation of secretion by cell                                                                                                         |  0.9270049 | 0.9173255 |
| <GO:0044088> | regulation of vacuole organization                                                                                                               |  0.9266763 | 0.9173255 |
| <GO:0090257> | regulation of muscle system process                                                                                                              |  0.9237780 | 0.9173255 |
| <GO:0007179> | transforming growth factor beta receptor signaling pathway                                                                                       |  0.9228901 | 0.9173255 |
| <GO:1903322> | positive regulation of protein modification by small protein conjugation or removal                                                              |  0.9227334 | 0.9173255 |
| <GO:0044788> | host-mediated perturbation of viral process                                                                                                      | -0.9225652 | 0.9173255 |
| <GO:0070849> | response to epidermal growth factor                                                                                                              |  0.9209471 | 0.9173255 |
| <GO:0051983> | regulation of chromosome segregation                                                                                                             |  0.9188403 | 0.9173255 |
| <GO:0006611> | protein export from nucleus                                                                                                                      |  0.9163578 | 0.9173255 |
| <GO:0071230> | cellular response to amino acid stimulus                                                                                                         |  0.9160008 | 0.9173255 |
| <GO:2000403> | positive regulation of lymphocyte migration                                                                                                      |  0.9146615 | 0.9173255 |
| <GO:0006584> | catecholamine metabolic process                                                                                                                  |  0.9110348 | 0.9173255 |
| <GO:0009712> | catechol-containing compound metabolic process                                                                                                   |  0.9110348 | 0.9173255 |
| <GO:0042417> | dopamine metabolic process                                                                                                                       |  0.9110348 | 0.9173255 |
| <GO:2001252> | positive regulation of chromosome organization                                                                                                   |  0.9099417 | 0.9173255 |
| <GO:0042446> | hormone biosynthetic process                                                                                                                     |  0.9085079 | 0.9173255 |
| <GO:0001706> | endoderm formation                                                                                                                               |  0.9080626 | 0.9173255 |
| <GO:1903035> | negative regulation of response to wounding                                                                                                      |  0.9071406 | 0.9173255 |
| <GO:0060998> | regulation of dendritic spine development                                                                                                        |  0.9040155 | 0.9173255 |
| <GO:0009116> | nucleoside metabolic process                                                                                                                     |  0.9016933 | 0.9173255 |
| <GO:0034762> | regulation of transmembrane transport                                                                                                            |  0.9459913 | 0.9179015 |
| <GO:0010977> | negative regulation of neuron projection development                                                                                             |  0.9098991 | 0.9179015 |
| <GO:0051896> | regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction                                                                 |  0.9412006 | 0.9195319 |
| <GO:0061512> | protein localization to cilium                                                                                                                   |  0.8976214 | 0.9197904 |
| <GO:0031346> | positive regulation of cell projection organization                                                                                              |  0.9429847 | 0.9204289 |
| <GO:0034315> | regulation of Arp2/3 complex-mediated actin nucleation                                                                                           | -0.9211251 | 0.9204289 |
| <GO:0043470> | regulation of carbohydrate catabolic process                                                                                                     |  0.9124357 | 0.9204289 |
| <GO:0016079> | synaptic vesicle exocytosis                                                                                                                      |  0.9069315 | 0.9204289 |
| <GO:0030071> | regulation of mitotic metaphase/anaphase transition                                                                                              |  0.9054395 | 0.9218992 |
| <GO:1902099> | regulation of metaphase/anaphase transition of cell cycle                                                                                        |  0.9054395 | 0.9218992 |
| <GO:0061387> | regulation of extent of cell growth                                                                                                              |  0.9044738 | 0.9225771 |
| <GO:0000122> | negative regulation of transcription by RNA polymerase II                                                                                        |  0.9591387 | 0.9227352 |
| <GO:0141091> | transforming growth factor beta receptor superfamily signaling pathway                                                                           |  0.9319128 | 0.9231568 |
| <GO:0046470> | phosphatidylcholine metabolic process                                                                                                            |  0.9045227 | 0.9242130 |
| <GO:0090101> | negative regulation of transmembrane receptor protein serine/threonine kinase signaling pathway                                                  |  0.9105466 | 0.9258595 |
| <GO:0048538> | thymus development                                                                                                                               |  0.8928367 | 0.9276536 |
| <GO:0043491> | phosphatidylinositol 3-kinase/protein kinase B signal transduction                                                                               |  0.9394661 | 0.9279461 |
| <GO:0006942> | regulation of striated muscle contraction                                                                                                        |  0.9008597 | 0.9279461 |
| <GO:0016445> | somatic diversification of immunoglobulins                                                                                                       |  0.9008217 | 0.9279461 |
| <GO:0043903> | regulation of biological process involved in symbiotic interaction                                                                               |  0.8959026 | 0.9279461 |
| <GO:0090174> | organelle membrane fusion                                                                                                                        |  0.9162309 | 0.9286272 |
| <GO:0015865> | purine nucleotide transport                                                                                                                      |  0.8888653 | 0.9286272 |
| <GO:0051503> | adenine nucleotide transport                                                                                                                     |  0.8888653 | 0.9286272 |
| <GO:0072089> | stem cell proliferation                                                                                                                          |  0.8897632 | 0.9297174 |
| <GO:0051607> | defense response to virus                                                                                                                        |  0.9385957 | 0.9317068 |
| <GO:0006664> | glycolipid metabolic process                                                                                                                     |  0.9053205 | 0.9317068 |
| <GO:1901184> | regulation of ERBB signaling pathway                                                                                                             |  0.8925310 | 0.9317068 |
| <GO:0055086> | nucleobase-containing small molecule metabolic process                                                                                           |  0.9540410 | 0.9318615 |
| <GO:0051924> | regulation of calcium ion transport                                                                                                              |  0.9146988 | 0.9322916 |
| <GO:0071772> | response to BMP                                                                                                                                  |  0.8957005 | 0.9324757 |
| <GO:0071773> | cellular response to BMP stimulus                                                                                                                |  0.8957005 | 0.9324757 |
| <GO:0042058> | regulation of epidermal growth factor receptor signaling pathway                                                                                 |  0.8953094 | 0.9324757 |
| <GO:0050773> | regulation of dendrite development                                                                                                               |  0.8960301 | 0.9334187 |
| <GO:0006109> | regulation of carbohydrate metabolic process                                                                                                     |  0.9098388 | 0.9335634 |
| <GO:0006487> | protein N-linked glycosylation                                                                                                                   |  0.8902884 | 0.9339096 |
| <GO:1902679> | negative regulation of RNA biosynthetic process                                                                                                  |  0.9602967 | 0.9348103 |
| <GO:1905897> | regulation of response to endoplasmic reticulum stress                                                                                           |  0.9017604 | 0.9348103 |
| <GO:0097061> | dendritic spine organization                                                                                                                     |  0.8930208 | 0.9350675 |
| <GO:0006509> | membrane protein ectodomain proteolysis                                                                                                          |  0.8928238 | 0.9350675 |
| <GO:0045824> | negative regulation of innate immune response                                                                                                    |  0.8853695 | 0.9350675 |
| <GO:0030866> | cortical actin cytoskeleton organization                                                                                                         |  0.8803080 | 0.9350675 |
| <GO:0014013> | regulation of gliogenesis                                                                                                                        |  0.8945899 | 0.9365267 |
| <GO:0042254> | ribosome biogenesis                                                                                                                              | -0.8498107 | 0.9370590 |
| <GO:1902905> | positive regulation of supramolecular fiber organization                                                                                         |  0.9089025 | 0.9374481 |
| <GO:0051345> | positive regulation of hydrolase activity                                                                                                        |  0.8990859 | 0.9374481 |
| <GO:0042987> | amyloid precursor protein catabolic process                                                                                                      |  0.8924124 | 0.9374481 |
| <GO:0045058> | T cell selection                                                                                                                                 |  0.8884229 | 0.9374481 |
| <GO:0007033> | vacuole organization                                                                                                                             |  0.9434753 | 0.9388279 |
| <GO:2001222> | regulation of neuron migration                                                                                                                   |  0.8853938 | 0.9402780 |
| <GO:0120031> | plasma membrane bounded cell projection assembly                                                                                                 |  0.9302745 | 0.9419055 |
| <GO:0006575> | modified amino acid metabolic process                                                                                                            |  0.8933354 | 0.9419055 |
| <GO:0042246> | tissue regeneration                                                                                                                              |  0.8844476 | 0.9419055 |
| <GO:0009566> | fertilization                                                                                                                                    |  0.8814269 | 0.9419055 |
| <GO:1900087> | positive regulation of G1/S transition of mitotic cell cycle                                                                                     |  0.8769008 | 0.9419055 |
| <GO:1902808> | positive regulation of cell cycle G1/S phase transition                                                                                          |  0.8769008 | 0.9419055 |
| <GO:0018958> | phenol-containing compound metabolic process                                                                                                     |  0.8835187 | 0.9420348 |
| <GO:0006140> | regulation of nucleotide metabolic process                                                                                                       |  0.8943468 | 0.9422535 |
| <GO:1900542> | regulation of purine nucleotide metabolic process                                                                                                |  0.8943468 | 0.9422535 |
| <GO:1902749> | regulation of cell cycle G2/M phase transition                                                                                                   |  0.8796099 | 0.9422535 |
| <GO:0009988> | cell-cell recognition                                                                                                                            |  0.8720517 | 0.9422535 |
| <GO:0030182> | neuron differentiation                                                                                                                           |  0.9561739 | 0.9435567 |
| <GO:0099175> | regulation of postsynapse organization                                                                                                           |  0.8956735 | 0.9445834 |
| <GO:0043627> | response to estrogen                                                                                                                             |  0.8792578 | 0.9445834 |
| <GO:0050688> | regulation of defense response to virus                                                                                                          |  0.8727578 | 0.9445834 |
| <GO:0006753> | nucleoside phosphate metabolic process                                                                                                           |  0.9427997 | 0.9458081 |
| <GO:0001701> | in utero embryonic development                                                                                                                   |  0.9276786 | 0.9458081 |
| <GO:0046486> | glycerolipid metabolic process                                                                                                                   |  0.9217179 | 0.9458081 |
| <GO:0048232> | male gamete generation                                                                                                                           |  0.9201892 | 0.9458081 |
| <GO:0051222> | positive regulation of protein transport                                                                                                         |  0.9159846 | 0.9458081 |
| <GO:0051017> | actin filament bundle assembly                                                                                                                   |  0.9007661 | 0.9458081 |
| <GO:0072676> | lymphocyte migration                                                                                                                             |  0.8904669 | 0.9458081 |
| <GO:0140888> | interferon-mediated signaling pathway                                                                                                            |  0.8872214 | 0.9458081 |
| <GO:1903214> | regulation of protein targeting to mitochondrion                                                                                                 |  0.8772401 | 0.9458081 |
| <GO:0010611> | regulation of cardiac muscle hypertrophy                                                                                                         |  0.8750813 | 0.9458081 |
| <GO:0014743> | regulation of muscle hypertrophy                                                                                                                 |  0.8750813 | 0.9458081 |
| <GO:0048709> | oligodendrocyte differentiation                                                                                                                  |  0.8722373 | 0.9458081 |
| <GO:0015868> | purine ribonucleotide transport                                                                                                                  |  0.8678787 | 0.9458081 |
| <GO:1990776> | response to angiotensin                                                                                                                          |  0.8634369 | 0.9458081 |
| <GO:0006405> | RNA export from nucleus                                                                                                                          | -0.7605945 | 0.9458081 |
| <GO:0097479> | synaptic vesicle localization                                                                                                                    |  0.8748263 | 0.9465320 |
| <GO:0071559> | response to transforming growth factor beta                                                                                                      |  0.9030257 | 0.9488973 |
| <GO:0071331> | cellular response to hexose stimulus                                                                                                             |  0.8833090 | 0.9514518 |
| <GO:0050905> | neuromuscular process                                                                                                                            |  0.8829483 | 0.9514518 |
| <GO:1902175> | regulation of oxidative stress-induced intrinsic apoptotic signaling pathway                                                                     |  0.8680460 | 0.9516601 |
| <GO:0032480> | negative regulation of type I interferon production                                                                                              |  0.8670895 | 0.9516601 |
| <GO:0030038> | contractile actin filament bundle assembly                                                                                                       |  0.8856813 | 0.9521901 |
| <GO:0043149> | stress fiber assembly                                                                                                                            |  0.8856813 | 0.9521901 |
| <GO:0055117> | regulation of cardiac muscle contraction                                                                                                         |  0.8573306 | 0.9521901 |
| <GO:0042982> | amyloid precursor protein metabolic process                                                                                                      |  0.8807041 | 0.9538704 |
| <GO:0051881> | regulation of mitochondrial membrane potential                                                                                                   |  0.8715333 | 0.9538704 |
| <GO:0044089> | positive regulation of cellular component biogenesis                                                                                             |  0.9279587 | 0.9549659 |
| <GO:0046677> | response to antibiotic                                                                                                                           |  0.8658727 | 0.9553445 |
| <GO:0045892> | negative regulation of DNA-templated transcription                                                                                               |  0.9502400 | 0.9553473 |
| <GO:0007338> | single fertilization                                                                                                                             |  0.8694957 | 0.9560930 |
| <GO:0051225> | spindle assembly                                                                                                                                 | -0.7541661 | 0.9581954 |
| <GO:0016049> | cell growth                                                                                                                                      |  0.9195655 | 0.9583539 |
| <GO:0070588> | calcium ion transmembrane transport                                                                                                              |  0.8964128 | 0.9583539 |
| <GO:0022412> | cellular process involved in reproduction in multicellular organism                                                                              |  0.8924397 | 0.9583539 |
| <GO:2000142> | regulation of DNA-templated transcription initiation                                                                                             |  0.8622183 | 0.9591770 |
| <GO:0009615> | response to virus                                                                                                                                |  0.9206542 | 0.9596562 |
| <GO:0048667> | cell morphogenesis involved in neuron differentiation                                                                                            |  0.9132159 | 0.9596562 |
| <GO:0099504> | synaptic vesicle cycle                                                                                                                           |  0.8893474 | 0.9596562 |
| <GO:0007062> | sister chromatid cohesion                                                                                                                        | -0.8820902 | 0.9596562 |
| <GO:0043457> | regulation of cellular respiration                                                                                                               | -0.8820902 | 0.9596562 |
| <GO:0007091> | metaphase/anaphase transition of mitotic cell cycle                                                                                              |  0.8565572 | 0.9596562 |
| <GO:0044784> | metaphase/anaphase transition of cell cycle                                                                                                      |  0.8565572 | 0.9596562 |
| <GO:0034205> | amyloid-beta formation                                                                                                                           |  0.8556319 | 0.9596562 |
| <GO:1903747> | regulation of establishment of protein localization to mitochondrion                                                                             |  0.8543532 | 0.9596562 |
| <GO:0051952> | regulation of amine transport                                                                                                                    |  0.8525517 | 0.9596562 |
| <GO:0043647> | inositol phosphate metabolic process                                                                                                             |  0.8520053 | 0.9596562 |
| <GO:0010662> | regulation of striated muscle cell apoptotic process                                                                                             |  0.8513032 | 0.9596562 |
| <GO:0010665> | regulation of cardiac muscle cell apoptotic process                                                                                              |  0.8513032 | 0.9596562 |
| <GO:0008089> | anterograde axonal transport                                                                                                                     | -0.8491078 | 0.9596562 |
| <GO:0051235> | maintenance of location                                                                                                                          |  0.8848505 | 0.9614118 |
| <GO:0006576> | biogenic amine metabolic process                                                                                                                 |  0.8509629 | 0.9614118 |
| <GO:0033619> | membrane protein proteolysis                                                                                                                     |  0.8483866 | 0.9639339 |
| <GO:0030865> | cortical cytoskeleton organization                                                                                                               |  0.8477919 | 0.9670991 |
| <GO:0048814> | regulation of dendrite morphogenesis                                                                                                             |  0.8435077 | 0.9674869 |
| <GO:0019637> | organophosphate metabolic process                                                                                                                |  0.9428056 | 0.9677394 |
| <GO:1902042> | negative regulation of extrinsic apoptotic signaling pathway via death domain receptors                                                          |  0.8439213 | 0.9677394 |
| <GO:0048640> | negative regulation of developmental growth                                                                                                      |  0.8385310 | 0.9677394 |
| <GO:2000406> | positive regulation of T cell migration                                                                                                          |  0.8463579 | 0.9678305 |
| <GO:0071709> | membrane assembly                                                                                                                                |  0.8375617 | 0.9681394 |
| <GO:0045216> | cell-cell junction organization                                                                                                                  |  0.8736107 | 0.9701542 |
| <GO:0006687> | glycosphingolipid metabolic process                                                                                                              |  0.8425574 | 0.9701616 |
| <GO:1903578> | regulation of ATP metabolic process                                                                                                              |  0.8612060 | 0.9707248 |
| <GO:0043066> | negative regulation of apoptotic process                                                                                                         |  0.9330274 | 0.9707885 |
| <GO:0006468> | protein phosphorylation                                                                                                                          |  0.9292918 | 0.9707885 |
| <GO:0043244> | regulation of protein-containing complex disassembly                                                                                             |  0.8572064 | 0.9707885 |
| <GO:0000086> | G2/M transition of mitotic cell cycle                                                                                                            |  0.8478291 | 0.9707885 |
| <GO:0072583> | clathrin-dependent endocytosis                                                                                                                   |  0.8472255 | 0.9707885 |
| <GO:0060338> | regulation of type I interferon-mediated signaling pathway                                                                                       |  0.8338998 | 0.9707885 |
| <GO:1901875> | positive regulation of post-translational protein modification                                                                                   |  0.8565491 | 0.9708328 |
| <GO:0051208> | sequestering of calcium ion                                                                                                                      |  0.8555258 | 0.9708328 |
| <GO:0141188> | nucleic acid catabolic process                                                                                                                   |  0.9088058 | 0.9708670 |
| <GO:0071364> | cellular response to epidermal growth factor stimulus                                                                                            |  0.8357357 | 0.9708670 |
| <GO:0072521> | purine-containing compound metabolic process                                                                                                     |  0.9258225 | 0.9713393 |
| <GO:0051304> | chromosome separation                                                                                                                            |  0.8295754 | 0.9713393 |
| <GO:0010332> | response to gamma radiation                                                                                                                      |  0.8290307 | 0.9713393 |
| <GO:1904951> | positive regulation of establishment of protein localization                                                                                     |  0.9043938 | 0.9716021 |
| <GO:1901524> | regulation of mitophagy                                                                                                                          |  0.8441328 | 0.9716021 |
| <GO:0071902> | positive regulation of protein serine/threonine kinase activity                                                                                  |  0.8402822 | 0.9716021 |
| <GO:0031342> | negative regulation of cell killing                                                                                                              |  0.8380599 | 0.9716021 |
| <GO:0042306> | regulation of protein import into nucleus                                                                                                        |  0.8252100 | 0.9716021 |
| <GO:0035050> | embryonic heart tube development                                                                                                                 |  0.8353803 | 0.9731712 |
| <GO:0048284> | organelle fusion                                                                                                                                 |  0.8692503 | 0.9743709 |
| <GO:0045069> | regulation of viral genome replication                                                                                                           |  0.8414313 | 0.9743709 |
| <GO:0033077> | T cell differentiation in thymus                                                                                                                 |  0.8355598 | 0.9748017 |
| <GO:0071560> | cellular response to transforming growth factor beta stimulus                                                                                    |  0.8797032 | 0.9768818 |
| <GO:0106027> | neuron projection organization                                                                                                                   |  0.8219697 | 0.9768818 |
| <GO:0051209> | release of sequestered calcium ion into cytosol                                                                                                  |  0.8513189 | 0.9775300 |
| <GO:0051282> | regulation of sequestering of calcium ion                                                                                                        |  0.8513189 | 0.9775300 |
| <GO:0051283> | negative regulation of sequestering of calcium ion                                                                                               |  0.8513189 | 0.9775300 |
| <GO:0000018> | regulation of DNA recombination                                                                                                                  |  0.8385962 | 0.9776646 |
| <GO:0031175> | neuron projection development                                                                                                                    |  0.9286281 | 0.9784809 |
| <GO:0030029> | actin filament-based process                                                                                                                     |  0.9225340 | 0.9784809 |
| <GO:0007283> | spermatogenesis                                                                                                                                  |  0.8859012 | 0.9784809 |
| <GO:0009299> | mRNA transcription                                                                                                                               |  0.8402143 | 0.9784809 |
| <GO:0051224> | negative regulation of protein transport                                                                                                         |  0.8399482 | 0.9784809 |
| <GO:0051492> | regulation of stress fiber assembly                                                                                                              |  0.8349072 | 0.9784809 |
| <GO:0097120> | receptor localization to synapse                                                                                                                 |  0.8279272 | 0.9784809 |
| <GO:0042771> | intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator                                                            |  0.8248974 | 0.9784809 |
| <GO:0002028> | regulation of sodium ion transport                                                                                                               |  0.8154964 | 0.9784809 |
| <GO:0031122> | cytoplasmic microtubule organization                                                                                                             | -0.8145074 | 0.9784809 |
| <GO:0007276> | gamete generation                                                                                                                                |  0.8929126 | 0.9818182 |
| <GO:0060260> | regulation of transcription initiation by RNA polymerase II                                                                                      | -0.8390314 | 0.9819639 |
| <GO:0048666> | neuron development                                                                                                                               |  0.9284085 | 0.9833405 |
| <GO:0010959> | regulation of metal ion transport                                                                                                                |  0.8844501 | 0.9833405 |
| <GO:0006606> | protein import into nucleus                                                                                                                      |  0.8603140 | 0.9833405 |
| <GO:0051170> | import into nucleus                                                                                                                              |  0.8603140 | 0.9833405 |
| <GO:0051054> | positive regulation of DNA metabolic process                                                                                                     |  0.8511981 | 0.9833405 |
| <GO:0070534> | protein K63-linked ubiquitination                                                                                                                |  0.8360033 | 0.9833405 |
| <GO:0009247> | glycolipid biosynthetic process                                                                                                                  |  0.8226379 | 0.9833405 |
| <GO:0031503> | protein-containing complex localization                                                                                                          |  0.8198106 | 0.9833405 |
| <GO:0046479> | glycosphingolipid catabolic process                                                                                                              |  0.8175264 | 0.9833405 |
| <GO:0019395> | fatty acid oxidation                                                                                                                             |  0.8116864 | 0.9833405 |
| <GO:0090317> | negative regulation of intracellular protein transport                                                                                           |  0.8113127 | 0.9833405 |
| <GO:0010660> | regulation of muscle cell apoptotic process                                                                                                      |  0.8094950 | 0.9833405 |
| <GO:0051057> | positive regulation of small GTPase mediated signal transduction                                                                                 |  0.8270468 | 0.9841426 |
| <GO:0071482> | cellular response to light stimulus                                                                                                              |  0.8259484 | 0.9849841 |
| <GO:0072659> | protein localization to plasma membrane                                                                                                          |  0.8834520 | 0.9856310 |
| <GO:0061620> | glycolytic process through glucose-6-phosphate                                                                                                   |  0.8055322 | 0.9856310 |
| <GO:0007034> | vacuolar transport                                                                                                                               |  0.8813967 | 0.9880034 |
| <GO:0009799> | specification of symmetry                                                                                                                        |  0.8112915 | 0.9880034 |
| <GO:0009855> | determination of bilateral symmetry                                                                                                              |  0.8112915 | 0.9880034 |
| <GO:0016241> | regulation of macroautophagy                                                                                                                     |  0.8762927 | 0.9919346 |
| <GO:0017148> | negative regulation of translation                                                                                                               |  0.8545154 | 0.9919346 |
| <GO:0021953> | central nervous system neuron differentiation                                                                                                    |  0.8314199 | 0.9919346 |
| <GO:0140029> | exocytic process                                                                                                                                 |  0.8231111 | 0.9919346 |
| <GO:0061615> | glycolytic process through fructose-6-phosphate                                                                                                  |  0.8144362 | 0.9919346 |
| <GO:0003170> | heart valve development                                                                                                                          |  0.8106885 | 0.9919346 |
| <GO:1902683> | regulation of receptor localization to synapse                                                                                                   |  0.8091720 | 0.9919346 |
| <GO:0032784> | regulation of DNA-templated transcription elongation                                                                                             | -0.8050088 | 0.9919346 |
| <GO:0034243> | regulation of transcription elongation by RNA polymerase II                                                                                      | -0.8050088 | 0.9919346 |
| <GO:0019320> | hexose catabolic process                                                                                                                         |  0.8031831 | 0.9949782 |
| <GO:0046365> | monosaccharide catabolic process                                                                                                                 |  0.8031831 | 0.9949782 |
| <GO:0035418> | protein localization to synapse                                                                                                                  |  0.7929309 | 0.9949782 |
| <GO:0050770> | regulation of axonogenesis                                                                                                                       |  0.8323893 | 0.9958900 |
| <GO:0009308> | amine metabolic process                                                                                                                          |  0.8131389 | 0.9970112 |
| <GO:0007043> | cell-cell junction assembly                                                                                                                      |  0.8149661 | 0.9979937 |
| <GO:0006261> | DNA-templated DNA replication                                                                                                                    |  0.8048258 | 0.9981898 |
| <GO:0045666> | positive regulation of neuron differentiation                                                                                                    |  0.7998901 | 0.9981898 |
| <GO:0000423> | mitophagy                                                                                                                                        |  0.8272524 | 0.9986625 |
| <GO:0006997> | nucleus organization                                                                                                                             |  0.8237007 | 0.9994763 |
| <GO:0048699> | generation of neurons                                                                                                                            |  0.9207785 | 1.0000000 |
| <GO:0009628> | response to abiotic stimulus                                                                                                                     |  0.9108628 | 1.0000000 |
| <GO:0000902> | cell morphogenesis                                                                                                                               |  0.9096050 | 1.0000000 |
| <GO:0120036> | plasma membrane bounded cell projection organization                                                                                             |  0.9069480 | 1.0000000 |
| <GO:0043069> | negative regulation of programmed cell death                                                                                                     |  0.8987322 | 1.0000000 |
| <GO:0019953> | sexual reproduction                                                                                                                              |  0.8944521 | 1.0000000 |
| <GO:0030036> | actin cytoskeleton organization                                                                                                                  |  0.8929456 | 1.0000000 |
| <GO:0006793> | phosphorus metabolic process                                                                                                                     |  0.8898851 | 1.0000000 |
| <GO:0006796> | phosphate-containing compound metabolic process                                                                                                  |  0.8898851 | 1.0000000 |
| <GO:0006812> | monoatomic cation transport                                                                                                                      |  0.8891326 | 1.0000000 |
| <GO:0030030> | cell projection organization                                                                                                                     |  0.8886540 | 1.0000000 |
| <GO:0009895> | negative regulation of catabolic process                                                                                                         |  0.8884159 | 1.0000000 |
| <GO:0048609> | multicellular organismal reproductive process                                                                                                    |  0.8818286 | 1.0000000 |
| <GO:0031669> | cellular response to nutrient levels                                                                                                             |  0.8782732 | 1.0000000 |
| <GO:0006644> | phospholipid metabolic process                                                                                                                   |  0.8778050 | 1.0000000 |
| <GO:2000059> | negative regulation of ubiquitin-dependent protein catabolic process                                                                             | -0.8774931 | 1.0000000 |
| <GO:0061024> | membrane organization                                                                                                                            |  0.8753714 | 1.0000000 |
| <GO:0031648> | protein destabilization                                                                                                                          | -0.8741438 | 1.0000000 |
| <GO:0051246> | regulation of protein metabolic process                                                                                                          |  0.8708730 | 1.0000000 |
| <GO:0016072> | rRNA metabolic process                                                                                                                           | -0.8691367 | 1.0000000 |
| <GO:0048812> | neuron projection morphogenesis                                                                                                                  |  0.8650840 | 1.0000000 |
| <GO:0070201> | regulation of establishment of protein localization                                                                                              |  0.8632726 | 1.0000000 |
| <GO:0043254> | regulation of protein-containing complex assembly                                                                                                |  0.8603695 | 1.0000000 |
| <GO:0051052> | regulation of DNA metabolic process                                                                                                              |  0.8595545 | 1.0000000 |
| <GO:0008654> | phospholipid biosynthetic process                                                                                                                |  0.8568339 | 1.0000000 |
| <GO:0043269> | regulation of monoatomic ion transport                                                                                                           |  0.8552593 | 1.0000000 |
| <GO:0051258> | protein polymerization                                                                                                                           |  0.8493278 | 1.0000000 |
| <GO:0046474> | glycerophospholipid biosynthetic process                                                                                                         |  0.8440454 | 1.0000000 |
| <GO:0032984> | protein-containing complex disassembly                                                                                                           |  0.8437362 | 1.0000000 |
| <GO:0009117> | nucleotide metabolic process                                                                                                                     |  0.8429848 | 1.0000000 |
| <GO:0019058> | viral life cycle                                                                                                                                 |  0.8411212 | 1.0000000 |
| <GO:0045944> | positive regulation of transcription by RNA polymerase II                                                                                        |  0.8408551 | 1.0000000 |
| <GO:0051651> | maintenance of location in cell                                                                                                                  |  0.8405351 | 1.0000000 |
| <GO:0034655> | nucleobase-containing compound catabolic process                                                                                                 |  0.8379744 | 1.0000000 |
| <GO:0010508> | positive regulation of autophagy                                                                                                                 |  0.8378397 | 1.0000000 |
| <GO:0006401> | RNA catabolic process                                                                                                                            |  0.8368146 | 1.0000000 |
| <GO:0031032> | actomyosin structure organization                                                                                                                |  0.8359997 | 1.0000000 |
| <GO:0048858> | cell projection morphogenesis                                                                                                                    |  0.8342754 | 1.0000000 |
| <GO:0120039> | plasma membrane bounded cell projection morphogenesis                                                                                            |  0.8342754 | 1.0000000 |
| <GO:0010256> | endomembrane system organization                                                                                                                 |  0.8341500 | 1.0000000 |
| <GO:0010506> | regulation of autophagy                                                                                                                          |  0.8334945 | 1.0000000 |
| <GO:1903828> | negative regulation of protein localization                                                                                                      |  0.8288833 | 1.0000000 |
| <GO:0007015> | actin filament organization                                                                                                                      |  0.8282805 | 1.0000000 |
| <GO:0043933> | protein-containing complex organization                                                                                                          |  0.8280743 | 1.0000000 |
| <GO:0051247> | positive regulation of protein metabolic process                                                                                                 |  0.8253910 | 1.0000000 |
| <GO:0060560> | developmental growth involved in morphogenesis                                                                                                   |  0.8252803 | 1.0000000 |
| <GO:0043467> | regulation of generation of precursor metabolites and energy                                                                                     |  0.8249407 | 1.0000000 |
| <GO:0120035> | regulation of plasma membrane bounded cell projection organization                                                                               |  0.8224074 | 1.0000000 |
| <GO:0007010> | cytoskeleton organization                                                                                                                        |  0.8196802 | 1.0000000 |
| <GO:1990778> | protein localization to cell periphery                                                                                                           |  0.8189571 | 1.0000000 |
| <GO:0071453> | cellular response to oxygen levels                                                                                                               |  0.8186016 | 1.0000000 |
| <GO:0071705> | nitrogen compound transport                                                                                                                      |  0.8176731 | 1.0000000 |
| <GO:0048588> | developmental cell growth                                                                                                                        |  0.8176528 | 1.0000000 |
| <GO:0030162> | regulation of proteolysis                                                                                                                        |  0.8171408 | 1.0000000 |
| <GO:0045664> | regulation of neuron differentiation                                                                                                             |  0.8165102 | 1.0000000 |
| <GO:0051223> | regulation of protein transport                                                                                                                  |  0.8140615 | 1.0000000 |
| <GO:0120034> | positive regulation of plasma membrane bounded cell projection assembly                                                                          |  0.8129528 | 1.0000000 |
| <GO:0048285> | organelle fission                                                                                                                                |  0.8115462 | 1.0000000 |
| <GO:0006813> | potassium ion transport                                                                                                                          |  0.8098062 | 1.0000000 |
| <GO:0045930> | negative regulation of mitotic cell cycle                                                                                                        |  0.8093390 | 1.0000000 |
| <GO:0030509> | BMP signaling pathway                                                                                                                            |  0.8090625 | 1.0000000 |
| <GO:0051261> | protein depolymerization                                                                                                                         |  0.8081570 | 1.0000000 |
| <GO:0007032> | endosome organization                                                                                                                            |  0.8078324 | 1.0000000 |
| <GO:0048511> | rhythmic process                                                                                                                                 |  0.8047873 | 1.0000000 |
| <GO:0003300> | cardiac muscle hypertrophy                                                                                                                       |  0.8044830 | 1.0000000 |
| <GO:0014897> | striated muscle hypertrophy                                                                                                                      |  0.8044830 | 1.0000000 |
| <GO:1903829> | positive regulation of protein localization                                                                                                      |  0.8041508 | 1.0000000 |
| <GO:0048813> | dendrite morphogenesis                                                                                                                           |  0.8028846 | 1.0000000 |
| <GO:0031344> | regulation of cell projection organization                                                                                                       |  0.8027288 | 1.0000000 |
| <GO:0032273> | positive regulation of protein polymerization                                                                                                    |  0.8025716 | 1.0000000 |
| <GO:0045893> | positive regulation of DNA-templated transcription                                                                                               |  0.8025499 | 1.0000000 |
| <GO:1902680> | positive regulation of RNA biosynthetic process                                                                                                  |  0.8025499 | 1.0000000 |
| <GO:0044087> | regulation of cellular component biogenesis                                                                                                      |  0.8023495 | 1.0000000 |
| <GO:1905475> | regulation of protein localization to membrane                                                                                                   |  0.8002128 | 1.0000000 |
| <GO:2000401> | regulation of lymphocyte migration                                                                                                               |  0.7996956 | 1.0000000 |
| <GO:0071478> | cellular response to radiation                                                                                                                   |  0.7996563 | 1.0000000 |
| <GO:0032231> | regulation of actin filament bundle assembly                                                                                                     |  0.7995845 | 1.0000000 |
| <GO:0017015> | regulation of transforming growth factor beta receptor signaling pathway                                                                         |  0.7994577 | 1.0000000 |
| <GO:0006835> | dicarboxylic acid transport                                                                                                                      |  0.7990332 | 1.0000000 |
| <GO:0016358> | dendrite development                                                                                                                             |  0.7983387 | 1.0000000 |
| <GO:1901137> | carbohydrate derivative biosynthetic process                                                                                                     |  0.7977025 | 1.0000000 |
| <GO:0060996> | dendritic spine development                                                                                                                      |  0.7974606 | 1.0000000 |
| <GO:0065003> | protein-containing complex assembly                                                                                                              |  0.7974241 | 1.0000000 |
| <GO:0006626> | protein targeting to mitochondrion                                                                                                               |  0.7958437 | 1.0000000 |
| <GO:0001824> | blastocyst development                                                                                                                           |  0.7954117 | 1.0000000 |
| <GO:0006650> | glycerophospholipid metabolic process                                                                                                            |  0.7934601 | 1.0000000 |
| <GO:0032434> | regulation of proteasomal ubiquitin-dependent protein catabolic process                                                                          |  0.7932609 | 1.0000000 |
| <GO:0034220> | monoatomic ion transmembrane transport                                                                                                           |  0.7920358 | 1.0000000 |
| <GO:0033043> | regulation of organelle organization                                                                                                             |  0.7915641 | 1.0000000 |
| <GO:0060425> | lung morphogenesis                                                                                                                               |  0.7912869 | 1.0000000 |
| <GO:0099173> | postsynapse organization                                                                                                                         |  0.7892689 | 1.0000000 |
| <GO:0047496> | vesicle transport along microtubule                                                                                                              |  0.7881063 | 1.0000000 |
| <GO:0035987> | endodermal cell differentiation                                                                                                                  |  0.7880627 | 1.0000000 |
| <GO:0072332> | intrinsic apoptotic signaling pathway by p53 class mediator                                                                                      |  0.7878753 | 1.0000000 |
| <GO:0032880> | regulation of protein localization                                                                                                               |  0.7878441 | 1.0000000 |
| <GO:0033554> | cellular response to stress                                                                                                                      |  0.7875582 | 1.0000000 |
| <GO:0070482> | response to oxygen levels                                                                                                                        |  0.7873138 | 1.0000000 |
| <GO:0030833> | regulation of actin filament polymerization                                                                                                      |  0.7872243 | 1.0000000 |
| <GO:0043487> | regulation of RNA stability                                                                                                                      |  0.7869099 | 1.0000000 |
| <GO:0000077> | DNA damage checkpoint signaling                                                                                                                  |  0.7864011 | 1.0000000 |
| <GO:0031570> | DNA integrity checkpoint signaling                                                                                                               |  0.7864011 | 1.0000000 |
| <GO:0051894> | positive regulation of focal adhesion assembly                                                                                                   |  0.7848705 | 1.0000000 |
| <GO:1903076> | regulation of protein localization to plasma membrane                                                                                            |  0.7848403 | 1.0000000 |
| <GO:0070830> | bicellular tight junction assembly                                                                                                               |  0.7847799 | 1.0000000 |
| <GO:0120192> | tight junction assembly                                                                                                                          |  0.7847799 | 1.0000000 |
| <GO:0006862> | nucleotide transport                                                                                                                             |  0.7835833 | 1.0000000 |
| <GO:0032507> | maintenance of protein location in cell                                                                                                          |  0.7818644 | 1.0000000 |
| <GO:1903844> | regulation of cellular response to transforming growth factor beta stimulus                                                                      |  0.7809574 | 1.0000000 |
| <GO:1903580> | positive regulation of ATP metabolic process                                                                                                     |  0.7808259 | 1.0000000 |
| <GO:0060416> | response to growth hormone                                                                                                                       |  0.7806579 | 1.0000000 |
| <GO:1990138> | neuron projection extension                                                                                                                      |  0.7801865 | 1.0000000 |
| <GO:0030219> | megakaryocyte differentiation                                                                                                                    |  0.7799849 | 1.0000000 |
| <GO:0110020> | regulation of actomyosin structure organization                                                                                                  |  0.7797612 | 1.0000000 |
| <GO:0051148> | negative regulation of muscle cell differentiation                                                                                               |  0.7783761 | 1.0000000 |
| <GO:0006163> | purine nucleotide metabolic process                                                                                                              |  0.7783224 | 1.0000000 |
| <GO:0008064> | regulation of actin polymerization or depolymerization                                                                                           |  0.7782874 | 1.0000000 |
| <GO:0030832> | regulation of actin filament length                                                                                                              |  0.7782874 | 1.0000000 |
| <GO:0032366> | intracellular sterol transport                                                                                                                   |  0.7780003 | 1.0000000 |
| <GO:0007006> | mitochondrial membrane organization                                                                                                              |  0.7777936 | 1.0000000 |
| <GO:0009267> | cellular response to starvation                                                                                                                  |  0.7776191 | 1.0000000 |
| <GO:0000045> | autophagosome assembly                                                                                                                           |  0.7751244 | 1.0000000 |
| <GO:0045454> | cell redox homeostasis                                                                                                                           |  0.7748748 | 1.0000000 |
| <GO:0044091> | membrane biogenesis                                                                                                                              |  0.7741618 | 1.0000000 |
| <GO:0090407> | organophosphate biosynthetic process                                                                                                             |  0.7737583 | 1.0000000 |
| <GO:0007093> | mitotic cell cycle checkpoint signaling                                                                                                          |  0.7724404 | 1.0000000 |
| <GO:0000422> | autophagy of mitochondrion                                                                                                                       |  0.7724338 | 1.0000000 |
| <GO:0031109> | microtubule polymerization or depolymerization                                                                                                   |  0.7720184 | 1.0000000 |
| <GO:0051851> | host-mediated perturbation of symbiont process                                                                                                   |  0.7704883 | 1.0000000 |
| <GO:0060428> | lung epithelium development                                                                                                                      |  0.7692912 | 1.0000000 |
| <GO:0021987> | cerebral cortex development                                                                                                                      |  0.7686512 | 1.0000000 |
| <GO:0006900> | vesicle budding from membrane                                                                                                                    |  0.7673311 | 1.0000000 |
| <GO:0051129> | negative regulation of cellular component organization                                                                                           |  0.7669388 | 1.0000000 |
| <GO:0051056> | regulation of small GTPase mediated signal transduction                                                                                          |  0.7668602 | 1.0000000 |
| <GO:0031647> | regulation of protein stability                                                                                                                  |  0.7665803 | 1.0000000 |
| <GO:0097435> | supramolecular fiber organization                                                                                                                |  0.7663983 | 1.0000000 |
| <GO:0060236> | regulation of mitotic spindle organization                                                                                                       |  0.7661945 | 1.0000000 |
| <GO:0090224> | regulation of spindle organization                                                                                                               |  0.7661945 | 1.0000000 |
| <GO:0043162> | ubiquitin-dependent protein catabolic process via the multivesicular body sorting pathway                                                        |  0.7660637 | 1.0000000 |
| <GO:0044782> | cilium organization                                                                                                                              |  0.7660536 | 1.0000000 |
| <GO:1901873> | regulation of post-translational protein modification                                                                                            |  0.7657205 | 1.0000000 |
| <GO:0035195> | miRNA-mediated post-transcriptional gene silencing                                                                                               |  0.7651591 | 1.0000000 |
| <GO:1904035> | regulation of epithelial cell apoptotic process                                                                                                  |  0.7647867 | 1.0000000 |
| <GO:0060395> | SMAD protein signal transduction                                                                                                                 |  0.7643614 | 1.0000000 |
| <GO:0009894> | regulation of catabolic process                                                                                                                  |  0.7635463 | 1.0000000 |
| <GO:0010389> | regulation of G2/M transition of mitotic cell cycle                                                                                              |  0.7635401 | 1.0000000 |
| <GO:0007623> | circadian rhythm                                                                                                                                 |  0.7633816 | 1.0000000 |
| <GO:0021549> | cerebellum development                                                                                                                           |  0.7623831 | 1.0000000 |
| <GO:0043271> | negative regulation of monoatomic ion transport                                                                                                  |  0.7617928 | 1.0000000 |
| <GO:0043525> | positive regulation of neuron apoptotic process                                                                                                  |  0.7616062 | 1.0000000 |
| <GO:0005976> | polysaccharide metabolic process                                                                                                                 |  0.7612588 | 1.0000000 |
| <GO:1905037> | autophagosome organization                                                                                                                       |  0.7606332 | 1.0000000 |
| <GO:0034340> | response to type I interferon                                                                                                                    |  0.7578946 | 1.0000000 |
| <GO:0060337> | type I interferon-mediated signaling pathway                                                                                                     |  0.7578946 | 1.0000000 |
| <GO:0071357> | cellular response to type I interferon                                                                                                           |  0.7578946 | 1.0000000 |
| <GO:1900744> | regulation of p38MAPK cascade                                                                                                                    |  0.7577018 | 1.0000000 |
| <GO:0051254> | positive regulation of RNA metabolic process                                                                                                     |  0.7574086 | 1.0000000 |
| <GO:0006366> | transcription by RNA polymerase II                                                                                                               |  0.7555263 | 1.0000000 |
| <GO:0042594> | response to starvation                                                                                                                           |  0.7551200 | 1.0000000 |
| <GO:0098660> | inorganic ion transmembrane transport                                                                                                            |  0.7550780 | 1.0000000 |
| <GO:0072678> | T cell migration                                                                                                                                 |  0.7543084 | 1.0000000 |
| <GO:0045786> | negative regulation of cell cycle                                                                                                                |  0.7532726 | 1.0000000 |
| <GO:0039532> | negative regulation of cytoplasmic pattern recognition receptor signaling pathway                                                                |  0.7531241 | 1.0000000 |
| <GO:0090100> | positive regulation of transmembrane receptor protein serine/threonine kinase signaling pathway                                                  |  0.7528036 | 1.0000000 |
| <GO:0043547> | positive regulation of GTPase activity                                                                                                           |  0.7527502 | 1.0000000 |
| <GO:0033045> | regulation of sister chromatid segregation                                                                                                       |  0.7527096 | 1.0000000 |
| <GO:1903320> | regulation of protein modification by small protein conjugation or removal                                                                       |  0.7511203 | 1.0000000 |
| <GO:1901796> | regulation of signal transduction by p53 class mediator                                                                                          |  0.7510847 | 1.0000000 |
| <GO:0007292> | female gamete generation                                                                                                                         |  0.7507820 | 1.0000000 |
| <GO:0031396> | regulation of protein ubiquitination                                                                                                             |  0.7505739 | 1.0000000 |
| <GO:0097553> | calcium ion transmembrane import into cytosol                                                                                                    |  0.7500055 | 1.0000000 |
| <GO:0032970> | regulation of actin filament-based process                                                                                                       |  0.7493597 | 1.0000000 |
| <GO:0072657> | protein localization to membrane                                                                                                                 |  0.7487837 | 1.0000000 |
| <GO:0051961> | negative regulation of nervous system development                                                                                                |  0.7485652 | 1.0000000 |
| <GO:0035821> | modulation of process of another organism                                                                                                        |  0.7484083 | 1.0000000 |
| <GO:0090150> | establishment of protein localization to membrane                                                                                                |  0.7483685 | 1.0000000 |
| <GO:1901874> | negative regulation of post-translational protein modification                                                                                   |  0.7478958 | 1.0000000 |
| <GO:0080135> | regulation of cellular response to stress                                                                                                        |  0.7475212 | 1.0000000 |
| <GO:0001101> | response to acid chemical                                                                                                                        |  0.7472392 | 1.0000000 |
| <GO:0000288> | nuclear-transcribed mRNA catabolic process, deadenylation-dependent decay                                                                        | -0.7465638 | 1.0000000 |
| <GO:0090092> | regulation of transmembrane receptor protein serine/threonine kinase signaling pathway                                                           |  0.7458473 | 1.0000000 |
| <GO:0045685> | regulation of glial cell differentiation                                                                                                         |  0.7457801 | 1.0000000 |
| <GO:0000082> | G1/S transition of mitotic cell cycle                                                                                                            |  0.7439942 | 1.0000000 |
| <GO:0021543> | pallium development                                                                                                                              |  0.7426861 | 1.0000000 |
| <GO:0006357> | regulation of transcription by RNA polymerase II                                                                                                 |  0.7425944 | 1.0000000 |
| <GO:0042391> | regulation of membrane potential                                                                                                                 |  0.7415839 | 1.0000000 |
| <GO:0060251> | regulation of glial cell proliferation                                                                                                           |  0.7410593 | 1.0000000 |
| <GO:0007368> | determination of left/right symmetry                                                                                                             |  0.7406376 | 1.0000000 |
| <GO:0060972> | left/right pattern formation                                                                                                                     |  0.7406376 | 1.0000000 |
| <GO:0030838> | positive regulation of actin filament polymerization                                                                                             |  0.7402061 | 1.0000000 |
| <GO:0014896> | muscle hypertrophy                                                                                                                               |  0.7401971 | 1.0000000 |
| <GO:0140013> | meiotic nuclear division                                                                                                                         |  0.7401029 | 1.0000000 |
| <GO:0010812> | negative regulation of cell-substrate adhesion                                                                                                   |  0.7395402 | 1.0000000 |
| <GO:0050854> | regulation of antigen receptor-mediated signaling pathway                                                                                        |  0.7392750 | 1.0000000 |
| <GO:0032367> | intracellular cholesterol transport                                                                                                              |  0.7389570 | 1.0000000 |
| <GO:0032956> | regulation of actin cytoskeleton organization                                                                                                    |  0.7386466 | 1.0000000 |
| <GO:0032271> | regulation of protein polymerization                                                                                                             |  0.7385539 | 1.0000000 |
| <GO:0021954> | central nervous system neuron development                                                                                                        |  0.7378957 | 1.0000000 |
| <GO:0042177> | negative regulation of protein catabolic process                                                                                                 |  0.7375219 | 1.0000000 |
| <GO:0140056> | organelle localization by membrane tethering                                                                                                     |  0.7373716 | 1.0000000 |
| <GO:0043297> | apical junction assembly                                                                                                                         |  0.7373671 | 1.0000000 |
| <GO:0065004> | protein-DNA complex assembly                                                                                                                     |  0.7370012 | 1.0000000 |
| <GO:0030834> | regulation of actin filament depolymerization                                                                                                    |  0.7367160 | 1.0000000 |
| <GO:0097581> | lamellipodium organization                                                                                                                       | -0.7366383 | 1.0000000 |
| <GO:0045935> | positive regulation of nucleobase-containing compound metabolic process                                                                          |  0.7362930 | 1.0000000 |
| <GO:1903749> | positive regulation of establishment of protein localization to mitochondrion                                                                    |  0.7359702 | 1.0000000 |
| <GO:1903955> | positive regulation of protein targeting to mitochondrion                                                                                        |  0.7359702 | 1.0000000 |
| <GO:1902750> | negative regulation of cell cycle G2/M phase transition                                                                                          |  0.7349554 | 1.0000000 |
| <GO:0045981> | positive regulation of nucleotide metabolic process                                                                                              |  0.7347317 | 1.0000000 |
| <GO:1900544> | positive regulation of purine nucleotide metabolic process                                                                                       |  0.7347317 | 1.0000000 |
| <GO:0045732> | positive regulation of protein catabolic process                                                                                                 |  0.7339977 | 1.0000000 |
| <GO:0051321> | meiotic cell cycle                                                                                                                               |  0.7330362 | 1.0000000 |
| <GO:0006112> | energy reserve metabolic process                                                                                                                 |  0.7327588 | 1.0000000 |
| <GO:0007617> | mating behavior                                                                                                                                  |  0.7313032 | 1.0000000 |
| <GO:0042176> | regulation of protein catabolic process                                                                                                          |  0.7311123 | 1.0000000 |
| <GO:0006914> | autophagy                                                                                                                                        |  0.7303620 | 1.0000000 |
| <GO:0061919> | process utilizing autophagic mechanism                                                                                                           |  0.7303620 | 1.0000000 |
| <GO:0030512> | negative regulation of transforming growth factor beta receptor signaling pathway                                                                |  0.7294513 | 1.0000000 |
| <GO:0051493> | regulation of cytoskeleton organization                                                                                                          |  0.7291667 | 1.0000000 |
| <GO:1902745> | positive regulation of lamellipodium organization                                                                                                |  0.7288063 | 1.0000000 |
| <GO:0048488> | synaptic vesicle endocytosis                                                                                                                     |  0.7283121 | 1.0000000 |
| <GO:1903077> | negative regulation of protein localization to plasma membrane                                                                                   |  0.7281117 | 1.0000000 |
| <GO:1904376> | negative regulation of protein localization to cell periphery                                                                                    |  0.7281117 | 1.0000000 |
| <GO:1902115> | regulation of organelle assembly                                                                                                                 |  0.7280493 | 1.0000000 |
| <GO:0099068> | postsynapse assembly                                                                                                                             |  0.7277763 | 1.0000000 |
| <GO:0034767> | positive regulation of monoatomic ion transmembrane transport                                                                                    |  0.7257694 | 1.0000000 |
| <GO:0009259> | ribonucleotide metabolic process                                                                                                                 |  0.7253404 | 1.0000000 |
| <GO:0031345> | negative regulation of cell projection organization                                                                                              |  0.7253033 | 1.0000000 |
| <GO:0030837> | negative regulation of actin filament polymerization                                                                                             |  0.7247059 | 1.0000000 |
| <GO:0007422> | peripheral nervous system development                                                                                                            |  0.7246223 | 1.0000000 |
| <GO:0007005> | mitochondrion organization                                                                                                                       |  0.7244199 | 1.0000000 |
| <GO:0019693> | ribose phosphate metabolic process                                                                                                               |  0.7241554 | 1.0000000 |
| <GO:0006402> | mRNA catabolic process                                                                                                                           |  0.7229436 | 1.0000000 |
| <GO:0043488> | regulation of mRNA stability                                                                                                                     |  0.7228166 | 1.0000000 |
| <GO:0032511> | late endosome to vacuole transport via multivesicular body sorting pathway                                                                       |  0.7221673 | 1.0000000 |
| <GO:0045912> | negative regulation of carbohydrate metabolic process                                                                                            |  0.7219039 | 1.0000000 |
| <GO:0016236> | macroautophagy                                                                                                                                   |  0.7217159 | 1.0000000 |
| <GO:0045071> | negative regulation of viral genome replication                                                                                                  |  0.7211161 | 1.0000000 |
| <GO:0140014> | mitotic nuclear division                                                                                                                         |  0.7208448 | 1.0000000 |
| <GO:0006260> | DNA replication                                                                                                                                  |  0.7205175 | 1.0000000 |
| <GO:0061013> | regulation of mRNA catabolic process                                                                                                             |  0.7199577 | 1.0000000 |
| <GO:0007029> | endoplasmic reticulum organization                                                                                                               |  0.7196304 | 1.0000000 |
| <GO:1902017> | regulation of cilium assembly                                                                                                                    |  0.7190307 | 1.0000000 |
| <GO:0030518> | nuclear receptor-mediated steroid hormone signaling pathway                                                                                      |  0.7187332 | 1.0000000 |
| <GO:0043401> | steroid hormone receptor signaling pathway                                                                                                       |  0.7187332 | 1.0000000 |
| <GO:1901185> | negative regulation of ERBB signaling pathway                                                                                                    |  0.7185304 | 1.0000000 |
| <GO:0000280> | nuclear division                                                                                                                                 |  0.7183161 | 1.0000000 |
| <GO:0071867> | response to monoamine                                                                                                                            |  0.7171601 | 1.0000000 |
| <GO:0010507> | negative regulation of autophagy                                                                                                                 |  0.7162821 | 1.0000000 |
| <GO:0019098> | reproductive behavior                                                                                                                            |  0.7161809 | 1.0000000 |
| <GO:0016032> | viral process                                                                                                                                    |  0.7145447 | 1.0000000 |
| <GO:0090288> | negative regulation of cellular response to growth factor stimulus                                                                               |  0.7123098 | 1.0000000 |
| <GO:2000134> | negative regulation of G1/S transition of mitotic cell cycle                                                                                     |  0.7113402 | 1.0000000 |
| <GO:0050691> | regulation of defense response to virus by host                                                                                                  |  0.7111426 | 1.0000000 |
| <GO:1990089> | response to nerve growth factor                                                                                                                  |  0.7107945 | 1.0000000 |
| <GO:1990090> | cellular response to nerve growth factor stimulus                                                                                                |  0.7107945 | 1.0000000 |
| <GO:0051262> | protein tetramerization                                                                                                                          |  0.7103255 | 1.0000000 |
| <GO:0002228> | natural killer cell mediated immunity                                                                                                            |  0.7096928 | 1.0000000 |
| <GO:0009416> | response to light stimulus                                                                                                                       |  0.7091081 | 1.0000000 |
| <GO:0006984> | ER-nucleus signaling pathway                                                                                                                     |  0.7088579 | 1.0000000 |
| <GO:0005977> | glycogen metabolic process                                                                                                                       |  0.7083587 | 1.0000000 |
| <GO:0044042> | glucan metabolic process                                                                                                                         |  0.7083587 | 1.0000000 |
| <GO:0044843> | cell cycle G1/S phase transition                                                                                                                 |  0.7076297 | 1.0000000 |
| <GO:0030520> | estrogen receptor signaling pathway                                                                                                              |  0.7074913 | 1.0000000 |
| <GO:0035270> | endocrine system development                                                                                                                     |  0.7070742 | 1.0000000 |
| <GO:0071539> | protein localization to centrosome                                                                                                               |  0.7066266 | 1.0000000 |
| <GO:0110053> | regulation of actin filament organization                                                                                                        |  0.7066045 | 1.0000000 |
| <GO:0060271> | cilium assembly                                                                                                                                  |  0.7065233 | 1.0000000 |
| <GO:0060491> | regulation of cell projection assembly                                                                                                           |  0.7040656 | 1.0000000 |
| <GO:0120032> | regulation of plasma membrane bounded cell projection assembly                                                                                   |  0.7040656 | 1.0000000 |
| <GO:0046824> | positive regulation of nucleocytoplasmic transport                                                                                               |  0.7037661 | 1.0000000 |
| <GO:1990000> | amyloid fibril formation                                                                                                                         |  0.7036217 | 1.0000000 |
| <GO:0032365> | intracellular lipid transport                                                                                                                    |  0.7032670 | 1.0000000 |
| <GO:0009411> | response to UV                                                                                                                                   |  0.7028560 | 1.0000000 |
| <GO:0035264> | multicellular organism growth                                                                                                                    |  0.7024779 | 1.0000000 |
| <GO:1900271> | regulation of long-term synaptic potentiation                                                                                                    |  0.7024126 | 1.0000000 |
| <GO:0030835> | negative regulation of actin filament depolymerization                                                                                           |  0.7023152 | 1.0000000 |
| <GO:0051693> | actin filament capping                                                                                                                           |  0.7023152 | 1.0000000 |
| <GO:2000404> | regulation of T cell migration                                                                                                                   |  0.7007652 | 1.0000000 |
| <GO:0071869> | response to catecholamine                                                                                                                        |  0.7002178 | 1.0000000 |
| <GO:0031468> | nuclear membrane reassembly                                                                                                                      |  0.6992638 | 1.0000000 |
| <GO:0007346> | regulation of mitotic cell cycle                                                                                                                 |  0.6991636 | 1.0000000 |
| <GO:0043242> | negative regulation of protein-containing complex disassembly                                                                                    |  0.6978541 | 1.0000000 |
| <GO:0048593> | camera-type eye morphogenesis                                                                                                                    |  0.6975587 | 1.0000000 |
| <GO:1902807> | negative regulation of cell cycle G1/S phase transition                                                                                          |  0.6974799 | 1.0000000 |
| <GO:0022406> | membrane docking                                                                                                                                 |  0.6970247 | 1.0000000 |
| <GO:0008306> | associative learning                                                                                                                             |  0.6963253 | 1.0000000 |
| <GO:2000060> | positive regulation of ubiquitin-dependent protein catabolic process                                                                             |  0.6960223 | 1.0000000 |
| <GO:1903900> | regulation of viral life cycle                                                                                                                   |  0.6943523 | 1.0000000 |
| <GO:0140238> | presynaptic endocytosis                                                                                                                          |  0.6940222 | 1.0000000 |
| <GO:1902117> | positive regulation of organelle assembly                                                                                                        |  0.6933852 | 1.0000000 |
| <GO:0000272> | polysaccharide catabolic process                                                                                                                 |  0.6928345 | 1.0000000 |
| <GO:0005980> | glycogen catabolic process                                                                                                                       |  0.6928345 | 1.0000000 |
| <GO:0009251> | glucan catabolic process                                                                                                                         |  0.6928345 | 1.0000000 |
| <GO:0009896> | positive regulation of catabolic process                                                                                                         |  0.6918845 | 1.0000000 |
| <GO:0044772> | mitotic cell cycle phase transition                                                                                                              |  0.6906909 | 1.0000000 |
| <GO:0016441> | post-transcriptional gene silencing                                                                                                              |  0.6901878 | 1.0000000 |
| <GO:0072595> | maintenance of protein localization in organelle                                                                                                 |  0.6896979 | 1.0000000 |
| <GO:0030010> | establishment of cell polarity                                                                                                                   |  0.6881938 | 1.0000000 |
| <GO:1901800> | positive regulation of proteasomal protein catabolic process                                                                                     |  0.6876337 | 1.0000000 |
| <GO:0014823> | response to activity                                                                                                                             |  0.6875872 | 1.0000000 |
| <GO:0007281> | germ cell development                                                                                                                            |  0.6868904 | 1.0000000 |
| <GO:2000765> | regulation of cytoplasmic translation                                                                                                            |  0.6865971 | 1.0000000 |
| <GO:0046488> | phosphatidylinositol metabolic process                                                                                                           |  0.6858618 | 1.0000000 |
| <GO:0097352> | autophagosome maturation                                                                                                                         |  0.6852602 | 1.0000000 |
| <GO:0070972> | protein localization to endoplasmic reticulum                                                                                                    |  0.6849485 | 1.0000000 |
| <GO:0031929> | TOR signaling                                                                                                                                    |  0.6845783 | 1.0000000 |
| <GO:0150116> | regulation of cell-substrate junction organization                                                                                               |  0.6825477 | 1.0000000 |
| <GO:0070925> | organelle assembly                                                                                                                               |  0.6824039 | 1.0000000 |
| <GO:2000045> | regulation of G1/S transition of mitotic cell cycle                                                                                              |  0.6815093 | 1.0000000 |
| <GO:0035773> | insulin secretion involved in cellular response to glucose stimulus                                                                              |  0.6800166 | 1.0000000 |
| <GO:0150115> | cell-substrate junction organization                                                                                                             |  0.6799717 | 1.0000000 |
| <GO:0015031> | protein transport                                                                                                                                |  0.6774521 | 1.0000000 |
| <GO:0045184> | establishment of protein localization                                                                                                            |  0.6764195 | 1.0000000 |
| <GO:0061136> | regulation of proteasomal protein catabolic process                                                                                              |  0.6758199 | 1.0000000 |
| <GO:0006913> | nucleocytoplasmic transport                                                                                                                      |  0.6750065 | 1.0000000 |
| <GO:0051169> | nuclear transport                                                                                                                                |  0.6750065 | 1.0000000 |
| <GO:0034063> | stress granule assembly                                                                                                                          |  0.6745340 | 1.0000000 |
| <GO:0036010> | protein localization to endosome                                                                                                                 |  0.6742769 | 1.0000000 |
| <GO:0009150> | purine ribonucleotide metabolic process                                                                                                          |  0.6737654 | 1.0000000 |
| <GO:0060341> | regulation of cellular localization                                                                                                              |  0.6728546 | 1.0000000 |
| <GO:0098655> | monoatomic cation transmembrane transport                                                                                                        |  0.6711145 | 1.0000000 |
| <GO:0071868> | cellular response to monoamine stimulus                                                                                                          |  0.6705725 | 1.0000000 |
| <GO:1901799> | negative regulation of proteasomal protein catabolic process                                                                                     |  0.6703690 | 1.0000000 |
| <GO:0044770> | cell cycle phase transition                                                                                                                      |  0.6702832 | 1.0000000 |
| <GO:0045862> | positive regulation of proteolysis                                                                                                               |  0.6701706 | 1.0000000 |
| <GO:0048525> | negative regulation of viral process                                                                                                             |  0.6700036 | 1.0000000 |
| <GO:0034644> | cellular response to UV                                                                                                                          |  0.6694255 | 1.0000000 |
| <GO:0001708> | cell fate specification                                                                                                                          |  0.6691278 | 1.0000000 |
| <GO:1905508> | protein localization to microtubule organizing center                                                                                            |  0.6681972 | 1.0000000 |
| <GO:0016311> | dephosphorylation                                                                                                                                |  0.6675518 | 1.0000000 |
| <GO:0034765> | regulation of monoatomic ion transmembrane transport                                                                                             |  0.6674089 | 1.0000000 |
| <GO:0035967> | cellular response to topologically incorrect protein                                                                                             |  0.6673839 | 1.0000000 |
| <GO:0036211> | protein modification process                                                                                                                     |  0.6669403 | 1.0000000 |
| <GO:0031397> | negative regulation of protein ubiquitination                                                                                                    |  0.6668038 | 1.0000000 |
| <GO:1905476> | negative regulation of protein localization to membrane                                                                                          |  0.6661688 | 1.0000000 |
| <GO:0019079> | viral genome replication                                                                                                                         |  0.6648481 | 1.0000000 |
| <GO:0036293> | response to decreased oxygen levels                                                                                                              |  0.6645004 | 1.0000000 |
| <GO:0045324> | late endosome to vacuole transport                                                                                                               |  0.6643561 | 1.0000000 |
| <GO:1903321> | negative regulation of protein modification by small protein conjugation or removal                                                              |  0.6631474 | 1.0000000 |
| <GO:0072329> | monocarboxylic acid catabolic process                                                                                                            |  0.6618666 | 1.0000000 |
| <GO:0051656> | establishment of organelle localization                                                                                                          |  0.6604860 | 1.0000000 |
| <GO:0007044> | cell-substrate junction assembly                                                                                                                 |  0.6592276 | 1.0000000 |
| <GO:0019082> | viral protein processing                                                                                                                         | -0.6581844 | 1.0000000 |
| <GO:0044380> | protein localization to cytoskeleton                                                                                                             |  0.6565214 | 1.0000000 |
| <GO:0009141> | nucleoside triphosphate metabolic process                                                                                                        |  0.6563392 | 1.0000000 |
| <GO:0009144> | purine nucleoside triphosphate metabolic process                                                                                                 |  0.6563392 | 1.0000000 |
| <GO:0045995> | regulation of embryonic development                                                                                                              |  0.6563305 | 1.0000000 |
| <GO:0009205> | purine ribonucleoside triphosphate metabolic process                                                                                             |  0.6554481 | 1.0000000 |
| <GO:0035194> | regulatory ncRNA-mediated post-transcriptional gene silencing                                                                                    |  0.6553702 | 1.0000000 |
| <GO:0050850> | positive regulation of calcium-mediated signaling                                                                                                |  0.6541736 | 1.0000000 |
| <GO:0048009> | insulin-like growth factor receptor signaling pathway                                                                                            |  0.6528812 | 1.0000000 |
| <GO:0150105> | protein localization to cell-cell junction                                                                                                       |  0.6522482 | 1.0000000 |
| <GO:0036257> | multivesicular body organization                                                                                                                 |  0.6512165 | 1.0000000 |
| <GO:0051668> | localization within membrane                                                                                                                     |  0.6507729 | 1.0000000 |
| <GO:0072384> | organelle transport along microtubule                                                                                                            |  0.6496233 | 1.0000000 |
| <GO:0050792> | regulation of viral process                                                                                                                      |  0.6486675 | 1.0000000 |
| <GO:0098662> | inorganic cation transmembrane transport                                                                                                         |  0.6480788 | 1.0000000 |
| <GO:0051279> | regulation of release of sequestered calcium ion into cytosol                                                                                    |  0.6477715 | 1.0000000 |
| <GO:1904375> | regulation of protein localization to cell periphery                                                                                             |  0.6472701 | 1.0000000 |
| <GO:0051649> | establishment of localization in cell                                                                                                            |  0.6472189 | 1.0000000 |
| <GO:0007098> | centrosome cycle                                                                                                                                 |  0.6467611 | 1.0000000 |
| <GO:0031023> | microtubule organizing center organization                                                                                                       |  0.6467611 | 1.0000000 |
| <GO:0045739> | positive regulation of DNA repair                                                                                                                | -0.6455452 | 1.0000000 |
| <GO:0043200> | response to amino acid                                                                                                                           |  0.6454817 | 1.0000000 |
| <GO:0099518> | vesicle cytoskeletal trafficking                                                                                                                 |  0.6444600 | 1.0000000 |
| <GO:0035094> | response to nicotine                                                                                                                             |  0.6440903 | 1.0000000 |
| <GO:0000075> | cell cycle checkpoint signaling                                                                                                                  |  0.6439850 | 1.0000000 |
| <GO:0032259> | methylation                                                                                                                                      |  0.6428561 | 1.0000000 |
| <GO:0010608> | post-transcriptional regulation of gene expression                                                                                               |  0.6423643 | 1.0000000 |
| <GO:0032535> | regulation of cellular component size                                                                                                            |  0.6419855 | 1.0000000 |
| <GO:0015800> | acidic amino acid transport                                                                                                                      |  0.6403203 | 1.0000000 |
| <GO:0051640> | organelle localization                                                                                                                           |  0.6401965 | 1.0000000 |
| <GO:0090148> | membrane fission                                                                                                                                 |  0.6387426 | 1.0000000 |
| <GO:1904064> | positive regulation of cation transmembrane transport                                                                                            |  0.6384063 | 1.0000000 |
| <GO:0098876> | vesicle-mediated transport to the plasma membrane                                                                                                |  0.6374822 | 1.0000000 |
| <GO:0060079> | excitatory postsynaptic potential                                                                                                                |  0.6370578 | 1.0000000 |
| <GO:1901879> | regulation of protein depolymerization                                                                                                           |  0.6368451 | 1.0000000 |
| <GO:0036294> | cellular response to decreased oxygen levels                                                                                                     |  0.6368393 | 1.0000000 |
| <GO:0007286> | spermatid development                                                                                                                            |  0.6365740 | 1.0000000 |
| <GO:0033044> | regulation of chromosome organization                                                                                                            |  0.6358498 | 1.0000000 |
| <GO:0051058> | negative regulation of small GTPase mediated signal transduction                                                                                 |  0.6358010 | 1.0000000 |
| <GO:0051893> | regulation of focal adhesion assembly                                                                                                            |  0.6351717 | 1.0000000 |
| <GO:0090109> | regulation of cell-substrate junction assembly                                                                                                   |  0.6351717 | 1.0000000 |
| <GO:0048041> | focal adhesion assembly                                                                                                                          |  0.6342883 | 1.0000000 |
| <GO:0007163> | establishment or maintenance of cell polarity                                                                                                    |  0.6333686 | 1.0000000 |
| <GO:1904062> | regulation of monoatomic cation transmembrane transport                                                                                          |  0.6329149 | 1.0000000 |
| <GO:0050768> | negative regulation of neurogenesis                                                                                                              |  0.6329131 | 1.0000000 |
| <GO:0044319> | wound healing, spreading of cells                                                                                                                |  0.6320189 | 1.0000000 |
| <GO:0090504> | epiboly                                                                                                                                          |  0.6320189 | 1.0000000 |
| <GO:0090505> | epiboly involved in wound healing                                                                                                                |  0.6320189 | 1.0000000 |
| <GO:0043412> | macromolecule modification                                                                                                                       |  0.6317945 | 1.0000000 |
| <GO:0015931> | nucleobase-containing compound transport                                                                                                         |  0.6300396 | 1.0000000 |
| <GO:0009309> | amine biosynthetic process                                                                                                                       |  0.6295290 | 1.0000000 |
| <GO:0042401> | biogenic amine biosynthetic process                                                                                                              |  0.6295290 | 1.0000000 |
| <GO:0006892> | post-Golgi vesicle-mediated transport                                                                                                            |  0.6293679 | 1.0000000 |
| <GO:1902806> | regulation of cell cycle G1/S phase transition                                                                                                   |  0.6289382 | 1.0000000 |
| <GO:0046034> | ATP metabolic process                                                                                                                            |  0.6263998 | 1.0000000 |
| <GO:0042149> | cellular response to glucose starvation                                                                                                          |  0.6232528 | 1.0000000 |
| <GO:1902903> | regulation of supramolecular fiber organization                                                                                                  |  0.6230553 | 1.0000000 |
| <GO:1903047> | mitotic cell cycle process                                                                                                                       |  0.6192015 | 1.0000000 |
| <GO:0072655> | establishment of protein localization to mitochondrion                                                                                           |  0.6188540 | 1.0000000 |
| <GO:0048477> | oogenesis                                                                                                                                        |  0.6186182 | 1.0000000 |
| <GO:1901988> | negative regulation of cell cycle phase transition                                                                                               |  0.6172158 | 1.0000000 |
| <GO:0032272> | negative regulation of protein polymerization                                                                                                    |  0.6158487 | 1.0000000 |
| <GO:0051726> | regulation of cell cycle                                                                                                                         |  0.6154515 | 1.0000000 |
| <GO:0009451> | RNA modification                                                                                                                                 |  0.6151209 | 1.0000000 |
| <GO:0046785> | microtubule polymerization                                                                                                                       |  0.6137552 | 1.0000000 |
| <GO:1903046> | meiotic cell cycle process                                                                                                                       |  0.6131122 | 1.0000000 |
| <GO:0040029> | epigenetic regulation of gene expression                                                                                                         |  0.6129024 | 1.0000000 |
| <GO:0039702> | viral budding via host ESCRT complex                                                                                                             |  0.6095288 | 1.0000000 |
| <GO:0061952> | midbody abscission                                                                                                                               |  0.6073492 | 1.0000000 |
| <GO:0019068> | virion assembly                                                                                                                                  |  0.6070372 | 1.0000000 |
| <GO:0010948> | negative regulation of cell cycle process                                                                                                        |  0.6068622 | 1.0000000 |
| <GO:0006605> | protein targeting                                                                                                                                |  0.6062209 | 1.0000000 |
| <GO:0051898> | negative regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction                                                        |  0.6061058 | 1.0000000 |
| <GO:0016197> | endosomal transport                                                                                                                              |  0.6056514 | 1.0000000 |
| <GO:0034446> | substrate adhesion-dependent cell spreading                                                                                                      |  0.6043629 | 1.0000000 |
| <GO:0006470> | protein dephosphorylation                                                                                                                        |  0.6037648 | 1.0000000 |
| <GO:0010212> | response to ionizing radiation                                                                                                                   |  0.6030848 | 1.0000000 |
| <GO:0030330> | DNA damage response, signal transduction by p53 class mediator                                                                                   |  0.6028245 | 1.0000000 |
| <GO:1900026> | positive regulation of substrate adhesion-dependent cell spreading                                                                               |  0.6026442 | 1.0000000 |
| <GO:0002715> | regulation of natural killer cell mediated immunity                                                                                              |  0.6022991 | 1.0000000 |
| <GO:0002011> | morphogenesis of an epithelial sheet                                                                                                             |  0.6008718 | 1.0000000 |
| <GO:0006091> | generation of precursor metabolites and energy                                                                                                   |  0.5995439 | 1.0000000 |
| <GO:0034332> | adherens junction organization                                                                                                                   |  0.5994282 | 1.0000000 |
| <GO:1902414> | protein localization to cell junction                                                                                                            |  0.5976454 | 1.0000000 |
| <GO:0044804> | nucleophagy                                                                                                                                      |  0.5960876 | 1.0000000 |
| <GO:0141193> | nuclear receptor-mediated signaling pathway                                                                                                      |  0.5925651 | 1.0000000 |
| <GO:0050856> | regulation of T cell receptor signaling pathway                                                                                                  |  0.5923902 | 1.0000000 |
| <GO:0006508> | proteolysis                                                                                                                                      |  0.5902845 | 1.0000000 |
| <GO:1900180> | regulation of protein localization to nucleus                                                                                                    |  0.5901987 | 1.0000000 |
| <GO:0001825> | blastocyst formation                                                                                                                             |  0.5886085 | 1.0000000 |
| <GO:1901293> | nucleoside phosphate biosynthetic process                                                                                                        |  0.5885176 | 1.0000000 |
| <GO:0032204> | regulation of telomere maintenance                                                                                                               |  0.5879109 | 1.0000000 |
| <GO:0000226> | microtubule cytoskeleton organization                                                                                                            |  0.5874957 | 1.0000000 |
| <GO:0046578> | regulation of Ras protein signal transduction                                                                                                    |  0.5874612 | 1.0000000 |
| <GO:0000184> | nuclear-transcribed mRNA catabolic process, nonsense-mediated decay                                                                              | -0.5870172 | 1.0000000 |
| <GO:0000278> | mitotic cell cycle                                                                                                                               |  0.5865576 | 1.0000000 |
| <GO:0010564> | regulation of cell cycle process                                                                                                                 |  0.5862811 | 1.0000000 |
| <GO:0007052> | mitotic spindle organization                                                                                                                     |  0.5860202 | 1.0000000 |
| <GO:0030534> | adult behavior                                                                                                                                   |  0.5844854 | 1.0000000 |
| <GO:0061462> | protein localization to lysosome                                                                                                                 |  0.5829046 | 1.0000000 |
| <GO:1900024> | regulation of substrate adhesion-dependent cell spreading                                                                                        |  0.5826721 | 1.0000000 |
| <GO:1903169> | regulation of calcium ion transmembrane transport                                                                                                |  0.5824023 | 1.0000000 |
| <GO:0021762> | substantia nigra development                                                                                                                     |  0.5823393 | 1.0000000 |
| <GO:1901990> | regulation of mitotic cell cycle phase transition                                                                                                |  0.5819583 | 1.0000000 |
| <GO:0035329> | hippo signaling                                                                                                                                  |  0.5813090 | 1.0000000 |
| <GO:0060078> | regulation of postsynaptic membrane potential                                                                                                    |  0.5812009 | 1.0000000 |
| <GO:0030968> | endoplasmic reticulum unfolded protein response                                                                                                  |  0.5802075 | 1.0000000 |
| <GO:0001666> | response to hypoxia                                                                                                                              |  0.5795477 | 1.0000000 |
| <GO:0036465> | synaptic vesicle recycling                                                                                                                       |  0.5791475 | 1.0000000 |
| <GO:0071456> | cellular response to hypoxia                                                                                                                     |  0.5776448 | 1.0000000 |
| <GO:0035330> | regulation of hippo signaling                                                                                                                    |  0.5776174 | 1.0000000 |
| <GO:0043112> | receptor metabolic process                                                                                                                       |  0.5761034 | 1.0000000 |
| <GO:0021536> | diencephalon development                                                                                                                         |  0.5759949 | 1.0000000 |
| <GO:0072698> | protein localization to microtubule cytoskeleton                                                                                                 |  0.5756189 | 1.0000000 |
| <GO:0099565> | chemical synaptic transmission, postsynaptic                                                                                                     |  0.5752411 | 1.0000000 |
| <GO:1901880> | negative regulation of protein depolymerization                                                                                                  |  0.5749497 | 1.0000000 |
| <GO:0022402> | cell cycle process                                                                                                                               |  0.5747986 | 1.0000000 |
| <GO:0036258> | multivesicular body assembly                                                                                                                     |  0.5735166 | 1.0000000 |
| <GO:0033555> | multicellular organismal response to stress                                                                                                      |  0.5733310 | 1.0000000 |
| <GO:0043555> | regulation of translation in response to stress                                                                                                  |  0.5731897 | 1.0000000 |
| <GO:0000281> | mitotic cytokinesis                                                                                                                              | -0.5717166 | 1.0000000 |
| <GO:0042770> | signal transduction in response to DNA damage                                                                                                    |  0.5714113 | 1.0000000 |
| <GO:0034504> | protein localization to nucleus                                                                                                                  |  0.5709828 | 1.0000000 |
| <GO:1900151> | regulation of nuclear-transcribed mRNA catabolic process, deadenylation-dependent decay                                                          | -0.5708366 | 1.0000000 |
| <GO:0032922> | circadian regulation of gene expression                                                                                                          |  0.5708152 | 1.0000000 |
| <GO:0038202> | TORC1 signaling                                                                                                                                  |  0.5705796 | 1.0000000 |
| <GO:0006417> | regulation of translation                                                                                                                        |  0.5701036 | 1.0000000 |
| <GO:0010639> | negative regulation of organelle organization                                                                                                    |  0.5676854 | 1.0000000 |
| <GO:0007017> | microtubule-based process                                                                                                                        |  0.5666994 | 1.0000000 |
| <GO:0051028> | mRNA transport                                                                                                                                   |  0.5656250 | 1.0000000 |
| <GO:0009314> | response to radiation                                                                                                                            |  0.5644368 | 1.0000000 |
| <GO:0001881> | receptor recycling                                                                                                                               |  0.5642063 | 1.0000000 |
| <GO:0035020> | regulation of Rac protein signal transduction                                                                                                    |  0.5637039 | 1.0000000 |
| <GO:0010972> | negative regulation of G2/M transition of mitotic cell cycle                                                                                     |  0.5626186 | 1.0000000 |
| <GO:0044818> | mitotic G2/M transition checkpoint                                                                                                               |  0.5626186 | 1.0000000 |
| <GO:1903052> | positive regulation of proteolysis involved in protein catabolic process                                                                         |  0.5626092 | 1.0000000 |
| <GO:0070936> | protein K48-linked ubiquitination                                                                                                                |  0.5622433 | 1.0000000 |
| <GO:0033157> | regulation of intracellular protein transport                                                                                                    |  0.5617762 | 1.0000000 |
| <GO:1902850> | microtubule cytoskeleton organization involved in mitosis                                                                                        |  0.5614348 | 1.0000000 |
| <GO:1901991> | negative regulation of mitotic cell cycle phase transition                                                                                       |  0.5611075 | 1.0000000 |
| <GO:0032386> | regulation of intracellular transport                                                                                                            |  0.5601330 | 1.0000000 |
| <GO:0009062> | fatty acid catabolic process                                                                                                                     |  0.5584254 | 1.0000000 |
| <GO:0021766> | hippocampus development                                                                                                                          |  0.5554158 | 1.0000000 |
| <GO:0003281> | ventricular septum development                                                                                                                   | -0.5536131 | 1.0000000 |
| <GO:0051276> | chromosome organization                                                                                                                          |  0.5531177 | 1.0000000 |
| <GO:0034620> | cellular response to unfolded protein                                                                                                            |  0.5522662 | 1.0000000 |
| <GO:1901987> | regulation of cell cycle phase transition                                                                                                        |  0.5515555 | 1.0000000 |
| <GO:0051926> | negative regulation of calcium ion transport                                                                                                     |  0.5504083 | 1.0000000 |
| <GO:0010591> | regulation of lamellipodium assembly                                                                                                             | -0.5495720 | 1.0000000 |
| <GO:0035966> | response to topologically incorrect protein                                                                                                      |  0.5465001 | 1.0000000 |
| <GO:0043543> | protein acylation                                                                                                                                |  0.5459159 | 1.0000000 |
| <GO:2000736> | regulation of stem cell differentiation                                                                                                          |  0.5454154 | 1.0000000 |
| <GO:0072665> | protein localization to vacuole                                                                                                                  |  0.5448368 | 1.0000000 |
| <GO:0008344> | adult locomotory behavior                                                                                                                        |  0.5421193 | 1.0000000 |
| <GO:0006259> | DNA metabolic process                                                                                                                            |  0.5417030 | 1.0000000 |
| <GO:0007049> | cell cycle                                                                                                                                       |  0.5411638 | 1.0000000 |
| <GO:0048515> | spermatid differentiation                                                                                                                        |  0.5408019 | 1.0000000 |
| <GO:0033365> | protein localization to organelle                                                                                                                |  0.5401632 | 1.0000000 |
| <GO:0006289> | nucleotide-excision repair                                                                                                                       |  0.5388984 | 1.0000000 |
| <GO:0006622> | protein targeting to lysosome                                                                                                                    |  0.5377162 | 1.0000000 |
| <GO:0032465> | regulation of cytokinesis                                                                                                                        |  0.5368891 | 1.0000000 |
| <GO:0072331> | signal transduction by p53 class mediator                                                                                                        |  0.5360313 | 1.0000000 |
| <GO:0046755> | viral budding                                                                                                                                    |  0.5359019 | 1.0000000 |
| <GO:0000381> | regulation of alternative mRNA splicing, via spliceosome                                                                                         |  0.5341258 | 1.0000000 |
| <GO:0032435> | negative regulation of proteasomal ubiquitin-dependent protein catabolic process                                                                 | -0.5330722 | 1.0000000 |
| <GO:0032006> | regulation of TOR signaling                                                                                                                      |  0.5323652 | 1.0000000 |
| <GO:0046907> | intracellular transport                                                                                                                          |  0.5283026 | 1.0000000 |
| <GO:0032008> | positive regulation of TOR signaling                                                                                                             |  0.5266953 | 1.0000000 |
| <GO:0021761> | limbic system development                                                                                                                        |  0.5259115 | 1.0000000 |
| <GO:0006635> | fatty acid beta-oxidation                                                                                                                        |  0.5247977 | 1.0000000 |
| <GO:0001764> | neuron migration                                                                                                                                 |  0.5236681 | 1.0000000 |
| <GO:0050657> | nucleic acid transport                                                                                                                           |  0.5236452 | 1.0000000 |
| <GO:0050658> | RNA transport                                                                                                                                    |  0.5236452 | 1.0000000 |
| <GO:0051236> | establishment of RNA localization                                                                                                                |  0.5236452 | 1.0000000 |
| <GO:0006661> | phosphatidylinositol biosynthetic process                                                                                                        |  0.5232968 | 1.0000000 |
| <GO:2000058> | regulation of ubiquitin-dependent protein catabolic process                                                                                      |  0.5218634 | 1.0000000 |
| <GO:0051650> | establishment of vesicle localization                                                                                                            |  0.5218195 | 1.0000000 |
| <GO:1903050> | regulation of proteolysis involved in protein catabolic process                                                                                  |  0.5216559 | 1.0000000 |
| <GO:0009057> | macromolecule catabolic process                                                                                                                  |  0.5215787 | 1.0000000 |
| <GO:1903051> | negative regulation of proteolysis involved in protein catabolic process                                                                         |  0.5211064 | 1.0000000 |
| <GO:1903432> | regulation of TORC1 signaling                                                                                                                    |  0.5198070 | 1.0000000 |
| <GO:0030705> | cytoskeleton-dependent intracellular transport                                                                                                   |  0.5196291 | 1.0000000 |
| <GO:0045947> | negative regulation of translational initiation                                                                                                  |  0.5193878 | 1.0000000 |
| <GO:0031333> | negative regulation of protein-containing complex assembly                                                                                       |  0.5165651 | 1.0000000 |
| <GO:0006403> | RNA localization                                                                                                                                 |  0.5164693 | 1.0000000 |
| <GO:0000380> | alternative mRNA splicing, via spliceosome                                                                                                       |  0.5161381 | 1.0000000 |
| <GO:0051648> | vesicle localization                                                                                                                             |  0.5124875 | 1.0000000 |
| <GO:0051168> | nuclear export                                                                                                                                   |  0.5089082 | 1.0000000 |
| <GO:0045911> | positive regulation of DNA recombination                                                                                                         |  0.5086906 | 1.0000000 |
| <GO:0000819> | sister chromatid segregation                                                                                                                     |  0.5036860 | 1.0000000 |
| <GO:0015844> | monoamine transport                                                                                                                              |  0.5032280 | 1.0000000 |
| <GO:0071806> | protein transmembrane transport                                                                                                                  |  0.5012464 | 1.0000000 |
| <GO:0030510> | regulation of BMP signaling pathway                                                                                                              |  0.5009164 | 1.0000000 |
| <GO:0000723> | telomere maintenance                                                                                                                             |  0.4980604 | 1.0000000 |
| <GO:0032200> | telomere organization                                                                                                                            |  0.4966469 | 1.0000000 |
| <GO:0090316> | positive regulation of intracellular protein transport                                                                                           |  0.4916517 | 1.0000000 |
| <GO:0006338> | chromatin remodeling                                                                                                                             |  0.4903422 | 1.0000000 |
| <GO:0031146> | SCF-dependent proteasomal ubiquitin-dependent protein catabolic process                                                                          |  0.4895126 | 1.0000000 |
| <GO:0072594> | establishment of protein localization to organelle                                                                                               |  0.4870999 | 1.0000000 |
| <GO:0010970> | transport along microtubule                                                                                                                      |  0.4858060 | 1.0000000 |
| <GO:0008361> | regulation of cell size                                                                                                                          |  0.4813658 | 1.0000000 |
| <GO:1902600> | proton transmembrane transport                                                                                                                   |  0.4783749 | 1.0000000 |
| <GO:0043279> | response to alkaloid                                                                                                                             |  0.4769671 | 1.0000000 |
| <GO:0015980> | energy derivation by oxidation of organic compounds                                                                                              |  0.4761393 | 1.0000000 |
| <GO:0030163> | protein catabolic process                                                                                                                        |  0.4747727 | 1.0000000 |
| <GO:0017004> | cytochrome complex assembly                                                                                                                      |  0.4734757 | 1.0000000 |
| <GO:0070507> | regulation of microtubule cytoskeleton organization                                                                                              |  0.4724253 | 1.0000000 |
| <GO:1905898> | positive regulation of response to endoplasmic reticulum stress                                                                                  |  0.4721356 | 1.0000000 |
| <GO:0006974> | DNA damage response                                                                                                                              |  0.4720274 | 1.0000000 |
| <GO:0032886> | regulation of microtubule-based process                                                                                                          |  0.4699545 | 1.0000000 |
| <GO:0007018> | microtubule-based movement                                                                                                                       |  0.4681565 | 1.0000000 |
| <GO:0006457> | protein folding                                                                                                                                  |  0.4673882 | 1.0000000 |
| <GO:0042752> | regulation of circadian rhythm                                                                                                                   |  0.4603333 | 1.0000000 |
| <GO:0072522> | purine-containing compound biosynthetic process                                                                                                  |  0.4563544 | 1.0000000 |
| <GO:0099111> | microtubule-based transport                                                                                                                      |  0.4547394 | 1.0000000 |
| <GO:0048524> | positive regulation of viral process                                                                                                             |  0.4532433 | 1.0000000 |
| <GO:0043687> | post-translational protein modification                                                                                                          |  0.4527592 | 1.0000000 |
| <GO:0006886> | intracellular protein transport                                                                                                                  |  0.4522011 | 1.0000000 |
| <GO:0006406> | mRNA export from nucleus                                                                                                                         |  0.4520168 | 1.0000000 |
| <GO:0045814> | negative regulation of gene expression, epigenetic                                                                                               |  0.4500851 | 1.0000000 |
| <GO:0032388> | positive regulation of intracellular transport                                                                                                   |  0.4480425 | 1.0000000 |
| <GO:2001251> | negative regulation of chromosome organization                                                                                                   |  0.4477244 | 1.0000000 |
| <GO:1903311> | regulation of mRNA metabolic process                                                                                                             |  0.4457038 | 1.0000000 |
| <GO:2000781> | positive regulation of double-strand break repair                                                                                                |  0.4427713 | 1.0000000 |
| <GO:0046854> | phosphatidylinositol phosphate biosynthetic process                                                                                              |  0.4409149 | 1.0000000 |
| <GO:0031507> | heterochromatin formation                                                                                                                        |  0.4342497 | 1.0000000 |
| <GO:0045333> | cellular respiration                                                                                                                             |  0.4320193 | 1.0000000 |
| <GO:0090307> | mitotic spindle assembly                                                                                                                         |  0.4262584 | 1.0000000 |
| <GO:1902743> | regulation of lamellipodium organization                                                                                                         |  0.4260808 | 1.0000000 |
| <GO:0006325> | chromatin organization                                                                                                                           |  0.4253534 | 1.0000000 |
| <GO:0032007> | negative regulation of TOR signaling                                                                                                             |  0.4215947 | 1.0000000 |
| <GO:0140694> | membraneless organelle assembly                                                                                                                  |  0.4194371 | 1.0000000 |
| <GO:0000209> | protein polyubiquitination                                                                                                                       |  0.4156721 | 1.0000000 |
| <GO:0030032> | lamellipodium assembly                                                                                                                           |  0.4156548 | 1.0000000 |
| <GO:1990928> | response to amino acid starvation                                                                                                                |  0.4139124 | 1.0000000 |
| <GO:0016482> | cytosolic transport                                                                                                                              |  0.4090951 | 1.0000000 |
| <GO:0043161> | proteasome-mediated ubiquitin-dependent protein catabolic process                                                                                |  0.4083856 | 1.0000000 |
| <GO:0034389> | lipid droplet organization                                                                                                                       |  0.4080202 | 1.0000000 |
| <GO:0006986> | response to unfolded protein                                                                                                                     |  0.4075474 | 1.0000000 |
| <GO:0072666> | establishment of protein localization to vacuole                                                                                                 |  0.4012436 | 1.0000000 |
| <GO:0031047> | regulatory ncRNA-mediated gene silencing                                                                                                         |  0.3980510 | 1.0000000 |
| <GO:0048024> | regulation of mRNA splicing, via spliceosome                                                                                                     |  0.3974163 | 1.0000000 |
| <GO:0009060> | aerobic respiration                                                                                                                              |  0.3943663 | 1.0000000 |
| <GO:0098813> | nuclear chromosome segregation                                                                                                                   |  0.3938590 | 1.0000000 |
| <GO:0042147> | retrograde transport, endosome to Golgi                                                                                                          |  0.3926196 | 1.0000000 |
| <GO:0045727> | positive regulation of translation                                                                                                               |  0.3864011 | 1.0000000 |
| <GO:0006623> | protein targeting to vacuole                                                                                                                     |  0.3842814 | 1.0000000 |
| <GO:0007030> | Golgi organization                                                                                                                               |  0.3792003 | 1.0000000 |
| <GO:0016567> | protein ubiquitination                                                                                                                           |  0.3784267 | 1.0000000 |
| <GO:1904262> | negative regulation of TORC1 signaling                                                                                                           |  0.3782513 | 1.0000000 |
| <GO:0006412> | translation                                                                                                                                      |  0.3774343 | 1.0000000 |
| <GO:0051301> | cell division                                                                                                                                    |  0.3759795 | 1.0000000 |
| <GO:0032446> | protein modification by small protein conjugation                                                                                                |  0.3704232 | 1.0000000 |
| <GO:0034976> | response to endoplasmic reticulum stress                                                                                                         |  0.3692195 | 1.0000000 |
| <GO:0034198> | cellular response to amino acid starvation                                                                                                       |  0.3690101 | 1.0000000 |
| <GO:0051603> | proteolysis involved in protein catabolic process                                                                                                |  0.3641209 | 1.0000000 |
| <GO:0006310> | DNA recombination                                                                                                                                |  0.3565692 | 1.0000000 |
| <GO:0006511> | ubiquitin-dependent protein catabolic process                                                                                                    |  0.3533881 | 1.0000000 |
| <GO:0048193> | Golgi vesicle transport                                                                                                                          |  0.3524414 | 1.0000000 |
| <GO:0019941> | modification-dependent protein catabolic process                                                                                                 |  0.3455078 | 1.0000000 |
| <GO:0010498> | proteasomal protein catabolic process                                                                                                            |  0.3429121 | 1.0000000 |
| <GO:0043632> | modification-dependent macromolecule catabolic process                                                                                           |  0.3419158 | 1.0000000 |
| <GO:0016071> | mRNA metabolic process                                                                                                                           |  0.3373852 | 1.0000000 |
| <GO:0007059> | chromosome segregation                                                                                                                           |  0.3355763 | 1.0000000 |
| <GO:0009165> | nucleotide biosynthetic process                                                                                                                  |  0.3322083 | 1.0000000 |
| <GO:0070647> | protein modification by small protein conjugation or removal                                                                                     |  0.3293621 | 1.0000000 |
| <GO:0043484> | regulation of RNA splicing                                                                                                                       |  0.3158158 | 1.0000000 |
| <GO:0050684> | regulation of mRNA processing                                                                                                                    |  0.3084957 | 1.0000000 |
| <GO:0006282> | regulation of DNA repair                                                                                                                         |  0.2776209 | 1.0000000 |
| <GO:0000910> | cytokinesis                                                                                                                                      |  0.2602190 | 1.0000000 |

## Prepare structured evidence

`sn_prepare_*_evidence()` turns stored outputs into compact, reusable
evidence bundles. These are the objects that can later be turned into
prompts or passed into writing helpers.

``` r
knitr::kable(annotation_evidence$cluster_summary)
```

| cluster | n_cells |  fraction | top_markers                                                               |
|:--------|--------:|----------:|:--------------------------------------------------------------------------|
| 0       |     229 | 0.1884774 | BPI, ENSG00000289381, MTARC1, MCEMP1, VNN3P                               |
| 1       |     172 | 0.1415638 | IATPR, TTC39C-AS1, PI16, TNFRSF4, CRIP2                                   |
| 2       |     149 | 0.1226337 | EDAR, ANKRD55, ADTRP, TSHZ2, ENSG00000271774                              |
| 3       |     144 | 0.1185185 | TCL1A, ENSG00000257275, SCN3A, ENSG00000224610, KCNG1                     |
| 4       |     133 | 0.1094650 | LINC02446, CD8B, ENSG00000310107, CRTAM, GZMH                             |
| 5       |     125 | 0.1028807 | CLEC10A, FCER1A, CYP2S1, DTNA, ZBTB46                                     |
| 6       |      75 | 0.0617284 | SLC4A10, ENSG00000228033, IL23R, ADAM12, LINC01644                        |
| 7       |      61 | 0.0502058 | IGHA1, IGHG2, IGHG3, IGHG1, SSPN                                          |
| 8       |      54 | 0.0444444 | ENSG00000288970, ENSG00000291157, MYOM2, ENSG00000288782, ENSG00000294329 |
| 9       |      32 | 0.0263374 | MT-CYB, MT-CO3, MT-ATP6, MT-CO2, MT-ND3                                   |
| 10      |      26 | 0.0213992 | LYPD2, ENSG00000301038, UICLM, PPP1R17, LINC02345                         |
| 11      |      15 | 0.0123457 | CLDN5, PF4V1, CMTM5, PDZK1IP1, PDGFA-DT                                   |

``` r
knitr::kable(enrichment_evidence$top_terms)
```

| ID           | Description                                            |       NES |  p.adjust |
|:-------------|:-------------------------------------------------------|----------:|----------:|
| <GO:0042776> | proton motive force-driven mitochondrial ATP synthesis | -2.515028 | 0.0002216 |
| <GO:0006120> | mitochondrial electron transport, NADH to ubiquinone   | -2.099828 | 0.0427655 |
| <GO:0032543> | mitochondrial translation                              | -1.960416 | 0.0651602 |
| <GO:0098581> | detection of external biotic stimulus                  |  1.921567 | 0.0015675 |
| <GO:0009595> | detection of biotic stimulus                           |  1.912715 | 0.0031257 |

## Build prompts without binding to a specific provider

The first-class package object for the LLM layer is a prompt bundle, not
a network call. This keeps analysis and interpretation cleanly
separated.

``` r
cat(substr(annotation_prompt$user, 1, 900))
#> Task: annotation.
#> 
#> Audience: scientist.
#> 
#> Language: en.
#> 
#> Target style: concise annotation note.
#> 
#> 
#> 
#> Evidence:
#> 
#> task:
#> annotation
#> 
#> cluster_col:
#> seurat_clusters
#> 
#> source_de_name:
#> cluster_markers
#> 
#> analysis_method:
#> wilcox
#> 
#> species:
#> human
#> 
#> cluster_summary:
#> # A tibble: 8 × 4
#>   cluster n_cells fraction top_markers                                          
#>   <fct>     <int>    <dbl> <chr>                                                
#> 1 0           229   0.188  BPI, ENSG00000289381, MTARC1, MCEMP1, VNN3P          
#> 2 1           172   0.142  IATPR, TTC39C-AS1, PI16, TNFRSF4, CRIP2              
#> 3 2           149   0.123  EDAR, ANKRD55, ADTRP, TSHZ2, ENSG00000271774         
#> 4 3           144   0.119  TCL1A, ENSG00000257275
```

## Generate manuscript-style prompts

High-level helpers such as
[`sn_write_results()`](https://songqi.org/shennong/reference/sn_write_results.md)
can return prompts directly.

``` r
cat(substr(results_prompt$user, 1, 900))
#> Task: results.
#> 
#> Audience: scientist.
#> 
#> Language: en.
#> 
#> Target style: manuscript-style Results section.
#> 
#> 
#> 
#> Evidence:
#> 
#> task:
#> results
#> 
#> dataset:
#> n_cells:
#> 1215
#> 
#> n_features:
#> 54872
#> 
#> cluster_col:
#> seurat_clusters
#> 
#> clusters:
#> 12
#> 
#> cluster_summary:
#> # A tibble: 8 × 3
#>   cluster n_cells fraction
#>   <fct>     <int>    <dbl>
#> 1 0           229   0.188 
#> 2 1           172   0.142 
#> 3 2           149   0.123 
#> 4 3           144   0.119 
#> 5 4           133   0.109 
#> 6 5           125   0.103 
#> 7 6            75   0.0617
#> 8 7            61   0.0502
#> 
#> cluster_markers:
#> # A tibble: 8 × 4
#>   cluster n_cells fraction top_markers                                          
#>   
```

## Plug in a provider when needed

If you provide a function that accepts `messages` and returns text, the
same high-level helpers can store the generated interpretation back into
the Seurat object.

``` r
annotation_response
#> [1] "Mock interpretation generated from 2 messages using demo-model"
```

## Interpretation results are stored alongside analysis results

``` r
names(pbmc_interpreted@misc$interpretation_results)
#> [1] "cluster_annotation_note"
```
