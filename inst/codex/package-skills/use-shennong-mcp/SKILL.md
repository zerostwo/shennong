---
name: use-shennong-mcp
description: Use when an agent needs to discover Shennong analysis methods, inspect exported R function help, retrieve workflow recipes, or configure the bundled read-only Shennong MCP server before writing an omics workflow.
---

# Use Shennong MCP

Use the bundled MCP server as the fast discovery layer for the installed
Shennong version. Keep analysis execution in explicit R scripts so that inputs,
parameters, and stored results remain reviewable.

## Procedure

1. Configure a stdio MCP server from `Shennong::sn_mcp_server_config()` or run
   `Rscript -e 'Shennong::sn_mcp_server()'`.
2. Call `package_info` to confirm the installed version and Agent Skill path.
3. Call `list_methods` with an optional task before selecting a backend.
4. Call `method_status` to verify runtime availability, installation action,
   required inputs, and outputs.
5. Call `function_help` for the exact exported function signature.
6. Call `workflow_guide` for `api_map` or `workflow_recipes` when the task spans
   several analysis stages.
7. Write the final workflow with exported `sn_*` functions and retrieve stored
   evidence with `sn_list_results()` and `sn_get_result()`.

## Rules

- Treat all MCP tools as read-only discovery tools.
- Do not invent arguments that are absent from `function_help`.
- Do not replace explicit R analysis code with arbitrary execution hidden
  behind MCP calls.
- Prefer `sn_find_de()` for both single-cell and bulk differential expression;
  use `modality` only when automatic input dispatch is ambiguous.
- For RegVelo, use `sn_run_velocity(method = "regvelo")` and supply a
  regulator-target prior GRN through `backend_control$prior_grn`.
- Check `method_status` before preparing a pixi-backed workflow.

## References

- `../_shared/references/package_api_map.md`
- `../_shared/references/workflow_recipes.md`
