# Build an LLM prompt from structured Shennong evidence

Build an LLM prompt from structured Shennong evidence

## Usage

``` r
sn_build_prompt(
  evidence,
  task = c("annotation", "de", "enrichment", "results", "figure_legend",
    "presentation_summary"),
  style = NULL,
  audience = c("scientist", "clinician", "general"),
  language = c("en", "zh"),
  background = NULL,
  output_format = c("llm", "human"),
  include_json_schema = FALSE
)
```

## Arguments

- evidence:

  A structured evidence list created by `sn_prepare_*_evidence()`.

- task:

  Interpretation task type.

- style:

  Optional style instruction, for example `"manuscript"`.

- audience:

  Intended audience such as `"scientist"`.

- language:

  Output language.

- background:

  Optional user-supplied study background or biological context to
  inject into the prompt.

- output_format:

  One of `"llm"` for a model-ready prompt bundle or `"human"` for a
  human-readable markdown brief.

- include_json_schema:

  Whether to request structured JSON output.

## Value

A prompt bundle with `system`, `user`, and `messages`.

## Examples

``` r
evidence <- list(task = "annotation", cluster_summary = data.frame(cluster = "0", top_markers = "CD3D, TRAC"))
prompt <- sn_build_prompt(evidence = evidence, task = "annotation")
names(prompt)
#> [1] "output_format" "task"          "system"        "user"         
#> [5] "messages"      "evidence"     
```
