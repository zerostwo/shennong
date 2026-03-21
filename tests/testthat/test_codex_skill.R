library(testthat)

test_that("sn_get_codex_skill_path exposes the packaged Codex asset paths", {
  codex_root <- sn_get_codex_skill_path("codex_root")
  package_skills <- sn_get_codex_skill_path("package_skills")
  project_template <- sn_get_codex_skill_path("project_template")
  project_skills <- sn_get_codex_skill_path("project_template_skills")

  expect_true(dir.exists(codex_root))
  expect_true(dir.exists(package_skills))
  expect_true(dir.exists(project_template))
  expect_true(dir.exists(project_skills))
  expect_true(file.exists(file.path(project_template, "AGENTS.md")))
  expect_true(file.exists(file.path(package_skills, "use-shennong-project-init", "SKILL.md")))
})

test_that("sn_install_codex_skill installs packaged usage skills into a target directory", {
  target_root <- tempfile("skill-root-")
  dir.create(target_root)

  installed <- sn_install_codex_skill(
    path = target_root,
    type = "package_skills",
    overwrite = TRUE
  )

  expect_true(all(dir.exists(installed)))
  expect_true(file.exists(file.path(target_root, "use-shennong-project-init", "SKILL.md")))
  expect_true(file.exists(file.path(target_root, "manage-shennong-results", "SKILL.md")))
})

test_that("sn_install_codex_skill refuses to overwrite when disabled", {
  target_root <- tempfile("skill-root-")
  dir.create(target_root)
  dir.create(file.path(target_root, "use-shennong-project-init"))

  expect_error(
    sn_install_codex_skill(path = target_root, type = "package_skills", overwrite = FALSE),
    "overwrite = TRUE"
  )
})

test_that("sn_initialize_project writes the base analysis scaffold without governance", {
  project_root <- tempfile("analysis-project-")

  created <- sn_initialize_project(
    path = project_root,
    project_name = "Tumor atlas",
    objective = "Build a reproducible tumor atlas workflow with Shennong.",
    with_agent = FALSE,
    overwrite = FALSE
  )

  expect_true(dir.exists(created$project_dir))
  expect_true(dir.exists(created$data_raw))
  expect_true(dir.exists(created$data_processed))
  expect_true(dir.exists(created$data_metadata))
  expect_true(dir.exists(created$scripts))
  expect_true(dir.exists(created$notebooks))
  expect_true(dir.exists(created$runs))
  expect_true(dir.exists(created$results_figures))
  expect_true(dir.exists(created$results_tables))
  expect_true(dir.exists(created$results_reports))
  expect_true(file.exists(created$readme))
  expect_true(file.exists(created$config_default))
  expect_false("agents_md" %in% names(created))
  expect_false(file.exists(file.path(project_root, "AGENTS.md")))
  expect_false(dir.exists(file.path(project_root, "memory")))
  expect_false(dir.exists(file.path(project_root, "skills")))

  readme_text <- readLines(created$readme, warn = FALSE)
  expect_true(any(grepl("Tumor atlas", readme_text, fixed = TRUE)))
  expect_true(any(grepl("config/default.yaml", readme_text, fixed = TRUE)))
})

test_that("sn_initialize_project writes the agent memory scaffold when enabled", {
  project_root <- tempfile("agent-project-")

  created <- sn_initialize_project(
    path = project_root,
    project_name = "Tumor atlas",
    objective = "Build a reproducible tumor atlas workflow with Shennong.",
    with_agent = TRUE,
    overwrite = FALSE
  )

  expect_true(file.exists(created$agents_md))
  expect_true(file.exists(created$memory_prompt))
  expect_true(file.exists(created$memory_plan))
  expect_true(file.exists(created$memory_status))
  expect_true(file.exists(created$memory_decisions))
  expect_true(file.exists(created$conventions))
  expect_true(file.exists(file.path(created$skills, "update-project-memory", "SKILL.md")))
  expect_true(file.exists(created$config_default))

  agents_text <- readLines(created$agents_md, warn = FALSE)
  prompt_text <- readLines(created$memory_prompt, warn = FALSE)
  config_text <- readLines(created$config_default, warn = FALSE)

  expect_true(any(grepl("Tumor atlas Analysis Governance", agents_text, fixed = TRUE)))
  expect_true(any(grepl("config/default.yaml", agents_text, fixed = TRUE)))
  expect_true(any(grepl("Build a reproducible tumor atlas workflow with Shennong.", prompt_text, fixed = TRUE)))
  expect_true(any(grepl("cellranger_path", config_text, fixed = TRUE)))
})

test_that("sn_initialize_codex_project scaffolds the full packaged project template", {
  project_root <- tempfile("codex-project-")

  created <- sn_initialize_codex_project(path = project_root, overwrite = FALSE)

  expect_true(file.exists(created$agents_md))
  expect_true(file.exists(created$memory_prompt))
  expect_true(file.exists(file.path(created$skills, "create-analysis-run", "SKILL.md")))
  expect_true(file.exists(file.path(project_root, "results", "figures", ".gitkeep")))
})
