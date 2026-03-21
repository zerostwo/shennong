# Initialize a Shennong analysis project

This helper bootstraps a clean analysis-project structure for end users
who want to analyze their own data with Shennong. When
`with_agent = TRUE`, it materializes the packaged project-governance
template; otherwise it creates the same project skeleton without the
governance layer.

## Usage

``` r
sn_initialize_project(
  path = ".",
  project_name = NULL,
  objective =
    "Build a reproducible Shennong-based single-cell analysis workflow for this project.",
  with_agent = TRUE,
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

- with_agent:

  Logical; if `TRUE`, materialize the full Codex-facing governance
  template. Defaults to `TRUE`.

- overwrite:

  Whether to replace existing managed files copied from the template.
  Defaults to `FALSE`.

## Value

Invisibly returns a named list containing the initialized project
directory plus the main created directory and file paths.

## Examples

``` r
project_dir <- file.path(tempdir(), "analysis-project")
created <- sn_initialize_project(
  path = project_dir,
  project_name = "PBMC pilot study",
  objective = "Build a reproducible Shennong-based PBMC analysis workflow.",
  with_agent = TRUE,
  overwrite = TRUE
)
file.exists(created$readme)
#> [1] TRUE
file.exists(created$agents_md)
#> [1] TRUE
```
