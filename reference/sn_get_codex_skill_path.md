# Return installed Shennong Codex asset paths

This helper exposes the Codex-related assets bundled with the installed
Shennong package. It can return the package-level Codex asset root, the
packaged project template, or the package and project skill directories.

## Usage

``` r
sn_get_codex_skill_path(
  component = c("codex_root", "package_skills", "project_template",
    "project_template_skills")
)
```

## Arguments

- component:

  One of `"codex_root"`, `"package_skills"`, `"project_template"`, or
  `"project_template_skills"`.

## Value

A character scalar path to the requested bundled asset directory.

## Examples

``` r
if (requireNamespace("Shennong", quietly = TRUE)) {
  sn_get_codex_skill_path("codex_root")
  sn_get_codex_skill_path("package_skills")
}
#> [1] "/home/runner/work/_temp/Library/Shennong/codex/package-skills"
```
