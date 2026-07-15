test_that("MCP server negotiates lifecycle and lists read-only tools", {
  initialize <- Shennong:::.sn_mcp_handle_request(list(
    jsonrpc = "2.0",
    id = 1,
    method = "initialize",
    params = list(protocolVersion = "2025-11-25", capabilities = list())
  ))
  tools <- Shennong:::.sn_mcp_handle_request(list(
    jsonrpc = "2.0", id = 2, method = "tools/list", params = list()
  ))

  expect_identical(initialize$result$protocolVersion, "2025-11-25")
  expect_true("tools" %in% names(initialize$result$capabilities))
  expect_setequal(
    vapply(tools$result$tools, `[[`, character(1), "name"),
    c("list_methods", "method_status", "function_help", "workflow_guide", "package_info")
  )
  expect_true(all(vapply(
    tools$result$tools,
    function(tool) isTRUE(tool$annotations$readOnlyHint) && !isTRUE(tool$annotations$destructiveHint),
    logical(1)
  )))

  fallback <- Shennong:::.sn_mcp_handle_request(list(
    jsonrpc = "2.0",
    id = 3,
    method = "initialize",
    params = list(protocolVersion = "2099-01-01", capabilities = list())
  ))
  expect_identical(fallback$result$protocolVersion, "2025-11-25")
})

test_that("MCP tools return structured Shennong discovery results", {
  response <- Shennong:::.sn_mcp_handle_request(list(
    jsonrpc = "2.0",
    id = 3,
    method = "tools/call",
    params = list(
      name = "method_status",
      arguments = list(method = "regvelo", task = "velocity")
    )
  ))

  expect_false(response$result$isError)
  expect_identical(response$result$structuredContent$status$method, "regvelo")
  expect_true(response$result$structuredContent$status$implemented)

  help_response <- Shennong:::.sn_mcp_handle_request(list(
    jsonrpc = "2.0",
    id = 4,
    method = "tools/call",
    params = list(name = "function_help", arguments = list(function_name = "sn_find_de"))
  ))
  expect_false(help_response$result$isError)
  expect_match(help_response$result$structuredContent$help, "differential expression")

  failure <- Shennong:::.sn_mcp_handle_request(list(
    jsonrpc = "2.0",
    id = 5,
    method = "tools/call",
    params = list(name = "function_help", arguments = list(function_name = "system"))
  ))
  expect_true(failure$result$isError)
  expect_match(failure$result$content[[1]]$text, "sn_\\*")
})

test_that("MCP stdio loop emits one compact JSON-RPC response per request", {
  requests <- c(
    jsonlite::toJSON(list(
      jsonrpc = "2.0", id = 1, method = "initialize",
      params = list(protocolVersion = "2025-11-25", capabilities = list())
    ), auto_unbox = TRUE),
    jsonlite::toJSON(list(jsonrpc = "2.0", method = "notifications/initialized"), auto_unbox = TRUE),
    jsonlite::toJSON(list(jsonrpc = "2.0", id = 2, method = "tools/list"), auto_unbox = TRUE),
    jsonlite::toJSON(list(
      jsonrpc = "2.0", id = 3, method = "tools/call",
      params = list(name = "package_info", arguments = list())
    ), auto_unbox = TRUE)
  )
  input <- textConnection(requests, open = "r")
  output <- textConnection("mcp_output", open = "w", local = TRUE)
  on.exit(close(input), add = TRUE)
  sn_mcp_server(input = input, output = output)
  close(output)

  expect_length(mcp_output, 3L)
  parsed <- lapply(mcp_output, jsonlite::fromJSON, simplifyVector = FALSE)
  expect_identical(vapply(parsed, `[[`, integer(1), "id"), 1:3)
  expect_true(all(vapply(mcp_output, function(line) !grepl("\\n", line), logical(1))))
  expect_identical(parsed[[3]]$result$structuredContent$package, "Shennong")
})

test_that("MCP launcher configuration is directly consumable", {
  config <- sn_mcp_server_config()
  expect_true(file.exists(config$command))
  expect_identical(config$transport, "stdio")
  expect_match(config$args[[2]], "sn_mcp_server")
})
