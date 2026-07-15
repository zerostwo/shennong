.sn_mcp_protocol_versions <- function() {
  c("2025-11-25", "2025-06-18", "2025-03-26", "2024-11-05")
}

.sn_mcp_annotations <- function() {
  list(
    readOnlyHint = TRUE,
    destructiveHint = FALSE,
    idempotentHint = TRUE,
    openWorldHint = FALSE
  )
}

.sn_mcp_tools <- function() {
  annotations <- .sn_mcp_annotations()
  list(
    list(
      name = "list_methods",
      title = "List Shennong Methods",
      description = "List registered Shennong analysis methods and their local availability.",
      inputSchema = list(
        type = "object",
        properties = list(
          task = list(type = "string", description = "Optional task such as velocity, fate, or bulk_de."),
          available = list(type = "boolean", description = "Optionally retain only available or unavailable methods.")
        ),
        additionalProperties = FALSE
      ),
      annotations = annotations
    ),
    list(
      name = "method_status",
      title = "Inspect a Shennong Method",
      description = "Return runtime, installation, input, output, and availability details for one registered method.",
      inputSchema = list(
        type = "object",
        properties = list(
          method = list(type = "string", description = "Registered method name."),
          task = list(type = "string", description = "Optional task used to disambiguate the method.")
        ),
        required = list("method"),
        additionalProperties = FALSE
      ),
      annotations = annotations
    ),
    list(
      name = "function_help",
      title = "Read Shennong Function Help",
      description = "Read the installed R help page for an exported Shennong function without executing it.",
      inputSchema = list(
        type = "object",
        properties = list(
          function_name = list(type = "string", description = "Exported function name beginning with sn_.")
        ),
        required = list("function_name"),
        additionalProperties = FALSE
      ),
      annotations = annotations
    ),
    list(
      name = "workflow_guide",
      title = "Read a Shennong Agent Guide",
      description = "Read the bundled API map or workflow recipes used by Shennong-aware agents.",
      inputSchema = list(
        type = "object",
        properties = list(
          guide = list(
            type = "string",
            enum = list("api_map", "workflow_recipes"),
            description = "Bundled guide to retrieve."
          )
        ),
        required = list("guide"),
        additionalProperties = FALSE
      ),
      annotations = annotations
    ),
    list(
      name = "package_info",
      title = "Inspect the Shennong Agent Surface",
      description = "Return the installed Shennong version plus MCP and Agent Skill entry points.",
      inputSchema = list(type = "object", properties = list(), additionalProperties = FALSE),
      annotations = annotations
    )
  )
}

.sn_mcp_json_native <- function(x) {
  jsonlite::fromJSON(
    jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null", dataframe = "rows"),
    simplifyVector = FALSE
  )
}

.sn_mcp_guide_path <- function(guide) {
  file_name <- switch(
    guide,
    api_map = "package_api_map.md",
    workflow_recipes = "workflow_recipes.md",
    stop("Unknown guide. Choose `api_map` or `workflow_recipes`.", call. = FALSE)
  )
  installed <- system.file(
    "codex", "package-skills", "_shared", "references", file_name,
    package = "Shennong"
  )
  if (nzchar(installed) && file.exists(installed)) return(installed)
  source <- file.path(
    getwd(), "inst", "codex", "package-skills", "_shared", "references", file_name
  )
  if (file.exists(source)) return(source)
  stop("Could not locate the bundled Shennong agent guide.", call. = FALSE)
}

.sn_mcp_function_help <- function(function_name) {
  if (!is.character(function_name) || length(function_name) != 1L || !startsWith(function_name, "sn_")) {
    stop("`function_name` must be one exported Shennong `sn_*` function.", call. = FALSE)
  }
  if (!function_name %in% getNamespaceExports("Shennong")) {
    stop("Unknown exported Shennong function: ", function_name, call. = FALSE)
  }
  package_path <- find.package("Shennong")
  source_man <- file.path(package_path, "man")
  if (dir.exists(source_man)) {
    rd_files <- list.files(source_man, pattern = "[.]Rd$", full.names = TRUE)
    rd_db <- stats::setNames(lapply(rd_files, tools::parse_Rd), basename(rd_files))
  } else {
    rd_db <- tools::Rd_db("Shennong")
  }
  rd <- rd_db[[paste0(function_name, ".Rd")]]
  if (is.null(rd)) {
    matched <- vapply(
      rd_db,
      function(candidate) {
        aliases <- vapply(
          Filter(
            function(element) identical(attr(element, "Rd_tag"), "\\alias"),
            candidate
          ),
          function(element) paste(as.character(element), collapse = ""),
          character(1)
        )
        function_name %in% aliases
      },
      logical(1)
    )
    if (any(matched)) rd <- rd_db[[which(matched)[[1]]]]
  }
  if (is.null(rd)) {
    stop("No installed help page was found for ", function_name, ".", call. = FALSE)
  }
  help_lines <- character()
  output <- textConnection("help_lines", open = "w", local = TRUE)
  on.exit(close(output), add = TRUE)
  tools::Rd2txt(rd, out = output, package = "Shennong")
  help_text <- paste(help_lines, collapse = "\n")
  gsub(paste0(".", intToUtf8(8L)), "", help_text, perl = TRUE)
}

.sn_mcp_call_tool <- function(name, arguments = list()) {
  arguments <- arguments %||% list()
  payload <- switch(
    name,
    list_methods = {
      task <- arguments$task %||% NULL
      available <- arguments$available %||% NULL
      list(methods = .sn_mcp_json_native(sn_list_methods(task = task, available = available)))
    },
    method_status = list(status = .sn_mcp_json_native(sn_method_status(
      method = arguments$method %||% stop("`method` is required.", call. = FALSE),
      task = arguments$task %||% NULL
    ))),
    function_help = list(
      function_name = arguments$function_name %||% stop("`function_name` is required.", call. = FALSE),
      help = .sn_mcp_function_help(arguments$function_name)
    ),
    workflow_guide = {
      guide <- arguments$guide %||% stop("`guide` is required.", call. = FALSE)
      list(guide = guide, markdown = paste(readLines(.sn_mcp_guide_path(guide), warn = FALSE), collapse = "\n"))
    },
    package_info = list(
      package = "Shennong",
      version = as.character(utils::packageVersion("Shennong")),
      mcp = sn_mcp_server_config(),
      package_skills = sn_get_codex_skill_path("package_skills")
    ),
    stop("Unknown tool: ", name, call. = FALSE)
  )
  payload <- .sn_mcp_json_native(payload)
  text <- jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", pretty = TRUE)
  list(
    content = list(list(type = "text", text = as.character(text))),
    structuredContent = payload,
    isError = FALSE
  )
}

.sn_mcp_error <- function(id, code, message, data = NULL) {
  error <- list(code = as.integer(code), message = as.character(message))
  if (!is_null(data)) error$data <- data
  list(jsonrpc = "2.0", id = id, error = error)
}

.sn_mcp_handle_request <- function(request) {
  id <- request$id %||% NULL
  method <- request$method %||% ""
  params <- request$params %||% list()

  if (identical(method, "initialize")) {
    requested <- params$protocolVersion %||% .sn_mcp_protocol_versions()[[1]]
    supported <- .sn_mcp_protocol_versions()
    negotiated <- if (requested %in% supported) requested else supported[[1]]
    return(list(
      jsonrpc = "2.0",
      id = id,
      result = list(
        protocolVersion = negotiated,
        capabilities = list(tools = list(listChanged = FALSE)),
        serverInfo = list(
          name = "shennong",
          title = "Shennong R Package",
          version = as.character(utils::packageVersion("Shennong"))
        ),
        instructions = "Read-only discovery server for Shennong methods, R help, and agent workflow guides."
      )
    ))
  }

  if (identical(method, "ping")) {
    return(list(jsonrpc = "2.0", id = id, result = list()))
  }
  if (identical(method, "tools/list")) {
    return(list(jsonrpc = "2.0", id = id, result = list(tools = .sn_mcp_tools())))
  }
  if (identical(method, "tools/call")) {
    name <- params$name %||% ""
    known <- vapply(.sn_mcp_tools(), `[[`, character(1), "name")
    if (!name %in% known) {
      return(.sn_mcp_error(id, -32602L, paste0("Unknown tool: ", name)))
    }
    result <- tryCatch(
      .sn_mcp_call_tool(name, params$arguments %||% list()),
      error = function(error) {
        list(
          content = list(list(type = "text", text = conditionMessage(error))),
          isError = TRUE
        )
      }
    )
    return(list(jsonrpc = "2.0", id = id, result = result))
  }

  if (is_null(id)) return(NULL)
  .sn_mcp_error(id, -32601L, paste0("Method not found: ", method))
}

#' Return a stdio configuration for the Shennong MCP server
#'
#' @return A list containing the command and arguments needed to launch the
#'   bundled read-only MCP server.
#' @examples
#' sn_mcp_server_config()
#' @export
sn_mcp_server_config <- function() {
  list(
    command = file.path(R.home("bin"), "Rscript"),
    args = c("-e", "Shennong::sn_mcp_server()"),
    transport = "stdio"
  )
}

#' Run the read-only Shennong MCP server over stdio
#'
#' Starts a newline-delimited JSON-RPC 2.0 server implementing MCP lifecycle,
#' tool discovery, and tool calls. The server exposes package metadata and
#' documentation only; it does not execute arbitrary R code or modify analysis
#' files.
#'
#' @param input Input connection. Defaults to standard input.
#' @param output Output connection. Defaults to standard output.
#' @return Invisibly returns \code{NULL} when the input stream closes.
#' @references Model Context Protocol specification:
#'   \url{https://modelcontextprotocol.io/specification/2025-11-25}.
#' @examples
#' config <- sn_mcp_server_config()
#' config$transport
#' \dontrun{
#' sn_mcp_server()
#' }
#' @export
sn_mcp_server <- function(input = stdin(), output = stdout()) {
  repeat {
    line <- readLines(input, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (!nzchar(trimws(line))) next
    response <- tryCatch(
      .sn_mcp_handle_request(jsonlite::fromJSON(line, simplifyVector = FALSE)),
      error = function(error) .sn_mcp_error(NULL, -32700L, conditionMessage(error))
    )
    if (is_null(response)) next
    encoded <- jsonlite::toJSON(
      response,
      auto_unbox = TRUE,
      null = "null",
      na = "null",
      pretty = FALSE
    )
    writeLines(as.character(encoded), con = output, sep = "\n", useBytes = TRUE)
    flush(output)
  }
  invisible(NULL)
}
