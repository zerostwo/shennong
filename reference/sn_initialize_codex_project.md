# Initialize Codex-style project guidance for a Shennong analysis

This helper scaffolds a governed analysis project from the packaged
`inst/codex/project-template/` assets.

## Usage

``` r
sn_initialize_codex_project(
  path = ".",
  project_name = NULL,
  objective =
    "Build a reproducible Shennong-based single-cell analysis workflow for this project.",
  overwrite = FALSE
)
```

## Arguments

- path:

  Project directory to initialize. Defaults to the current working
  directory.

- project_name:

  Optional human-readable project name. If `NULL`, the final path
  component of `path` is used.

- objective:

  Short project objective written into the generated project README,
  memory, and config files.

- overwrite:

  Whether to replace existing managed files copied from the template.
  Defaults to `FALSE`.

## Value

Invisibly returns the same structure as
[`sn_initialize_project()`](https://songqi.org/shennong/reference/sn_initialize_project.md).

## Examples

``` r
project_dir <- file.path(tempdir(), "codex-analysis-project")
created <- sn_initialize_codex_project(
  path = project_dir,
  project_name = "PBMC pilot study",
  objective = "Build a reproducible Shennong-based PBMC analysis workflow.",
  overwrite = TRUE
)
file.exists(created$agents)
#> [1] TRUE
file.exists(created$prompt)
#> [1] TRUE
```
