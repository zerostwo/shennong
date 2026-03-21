# Install the bundled Shennong Codex skill for end-user agents

This helper installs the packaged Shennong Codex skills into a user
skill directory, for example `~/.agents/skills`. Package usage skills
are distinct from initialized-project governance skills.

## Usage

``` r
sn_install_codex_skill(
  path = "~/.agents/skills",
  type = c("package_skills", "project_skills", "both"),
  overwrite = TRUE
)
```

## Arguments

- path:

  Destination directory that will contain the `shennong` skill folders.
  Defaults to `"~/.agents/skills"`.

- type:

  One of `"package_skills"`, `"project_skills"`, or `"both"`.

- overwrite:

  Whether to replace an existing installed copy of the target skill
  directories. Defaults to `TRUE`.

## Value

Invisibly returns a named character vector of installed directories.

## Examples

``` r
target_dir <- file.path(tempdir(), "skills")
installed <- sn_install_codex_skill(
  path = target_dir,
  type = "package_skills",
  overwrite = TRUE
)
all(dir.exists(installed))
#> [1] TRUE
```
