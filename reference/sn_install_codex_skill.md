# Install the bundled Shennong Codex skill for end-user agents

This helper copies the packaged Shennong skill into a user skill
directory, for example `~/.agents/skills`. The installed skill is
intended for analysis environments where the package is already
available.

## Usage

``` r
sn_install_codex_skill(path = "~/.agents/skills", overwrite = TRUE)
```

## Arguments

- path:

  Destination directory that will contain the `shennong` skill folder.
  Defaults to `"~/.agents/skills"`.

- overwrite:

  Whether to replace an existing installed copy of the skill. Defaults
  to `TRUE`.

## Value

Invisibly returns the installed skill directory path.

## Examples

``` r
target_dir <- file.path(tempdir(), "skills")
installed <- sn_install_codex_skill(path = target_dir, overwrite = TRUE)
dir.exists(installed)
#> [1] TRUE
```
