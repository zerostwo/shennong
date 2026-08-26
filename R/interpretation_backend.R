# LLM client plumbing for interpretation workflows.
#
# Extracted from interpretation.R: provider defaults and guessing, retryable-error
# hardening, structured-output schemas and message building, response-text extraction,
# and ellmer provider construction plus provider health checks.

.sn_default_llm_api_key <- function() {
  key <- Sys.getenv("OPENAI_API_KEY", unset = "")
  key
}

.sn_default_llm_base_url <- function() {
  base_url <- Sys.getenv("OPENAI_BASE_URL", unset = "")
  if (!nzchar(base_url)) {
    base_url <- "https://api.openai.com/v1"
  }
  base_url
}

.sn_default_llm_model <- function() {
  model <- Sys.getenv("OPENAI_MODEL", unset = "")
  if (!nzchar(model)) {
    model <- "gpt-4.1"
  }
  model
}

.sn_default_reasoning_effort <- function() {
  effort <- Sys.getenv("OPENAI_REASONING_EFFORT", unset = "")
  if (!nzchar(effort)) {
    return(NULL)
  }
  effort
}

.sn_guess_provider_type <- function(base_url) {
  if (grepl("^https://api\\.openai\\.com(/|$)", base_url)) {
    return("openai")
  }
  "openai_compatible"
}

.sn_is_retryable_llm_error <- function(message) {
  if (!is.character(message) || length(message) != 1L || is.na(message)) {
    return(FALSE)
  }
  grepl(
    "lexical error|invalid char in json text|<!DOCTYPE html>|<html|502|503|504|gateway|temporar|timeout|timed out|connection reset",
    message,
    ignore.case = TRUE
  )
}

.sn_harden_llm_error <- function(message) {
  if (!is.character(message) || length(message) != 1L || is.na(message) || !nzchar(message)) {
    return("LLM request failed.")
  }
  if (.sn_is_retryable_llm_error(message)) {
    return(paste(
      "LLM request failed because the upstream provider returned an HTML/error page or another non-JSON proxy response.",
      "This is common with third-party proxy endpoints such as Sub2API when authentication, gateway routing, or temporary upstream health is unstable.",
      "Check the endpoint and credentials, then retry."
    ))
  }
  message
}

.sn_get_default_ellmer_provider <- function() {
  sn_make_ellmer_provider(
    api_key = .sn_default_llm_api_key(),
    base_url = .sn_default_llm_base_url(),
    model = .sn_default_llm_model(),
    provider_type = .sn_guess_provider_type(.sn_default_llm_base_url()),
    reasoning_effort = .sn_default_reasoning_effort()
  )
}

.sn_annotation_json_schema_text <- function() {
  paste(
    "Return valid JSON with keys `cluster_annotations` and `narrative_summary`.",
    "`cluster_annotations` must be an array of objects with:",
    "`cluster`, `primary_label`, `broad_label`, `confidence`, and `status`.",
    "Optional fields may include `alternatives`, `supporting_markers`, `supporting_functions`, `risk_flags`, `note`, and `recommended_checks`.",
    "`confidence` should use one of: high, medium, low.",
    "`status` should use one of: confident, ambiguous, possible_contamination, possible_low_quality, possible_transition.",
    "`risk_flags` should be an array that may include: contamination, low_quality, transitional_state, doublet_like, proliferating_state."
  )
}

.sn_annotation_structured_type <- function() {
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    return(NULL)
  }

  ellmer::type_object(
    "Structured cluster annotation result.",
    cluster_annotations = ellmer::type_array(
      ellmer::type_object(
        "One annotation record per cluster.",
        cluster = ellmer::type_string("Cluster identifier."),
        primary_label = ellmer::type_string("Most plausible cell type or state label."),
        broad_label = ellmer::type_string("Broader lineage or family label."),
        confidence = ellmer::type_enum(
          c("high", "medium", "low"),
          "Calibrated confidence label."
        ),
        status = ellmer::type_enum(
          c("confident", "ambiguous", "possible_contamination", "possible_low_quality", "possible_transition"),
          "Annotation status."
        ),
        alternatives = ellmer::type_array(
          ellmer::type_string("Alternative plausible label."),
          description = "Alternative labels when ambiguity remains.",
          required = FALSE
        ),
        supporting_markers = ellmer::type_array(
          ellmer::type_string("Marker gene supporting the annotation."),
          description = "Directly supportive marker genes.",
          required = FALSE
        ),
        supporting_functions = ellmer::type_array(
          ellmer::type_string("Functional term supporting the annotation."),
          description = "Directly supportive pathways or biological processes.",
          required = FALSE
        ),
        risk_flags = ellmer::type_array(
          ellmer::type_enum(
            c("contamination", "low_quality", "transitional_state", "doublet_like", "proliferating_state")
          ),
          description = "Risk flags or caveats.",
          required = FALSE
        ),
        note = ellmer::type_string(
          "Brief cluster-level rationale.",
          required = FALSE
        ),
        recommended_checks = ellmer::type_array(
          ellmer::type_string("Suggested manual follow-up check."),
          description = "Suggested validation checks.",
          required = FALSE
        )
      ),
      description = "Array of cluster annotation records."
    ),
    narrative_summary = ellmer::type_string(
      "Optional concise narrative summary across all clusters.",
      required = FALSE
    )
  )
}

.sn_build_messages <- function(system_prompt, user_prompt) {
  list(
    list(role = "system", content = system_prompt),
    list(role = "user", content = user_prompt)
  )
}

#' Run an LLM provider on a prepared message list
#'
#' @param messages A message list, typically from \code{sn_build_prompt()}.
#' @param provider A user-supplied function that accepts \code{messages} and
#'   returns text or a list containing \code{text}.
#' @param model Optional model identifier passed through to the provider.
#' @param structured_type Optional structured-output schema or type object
#'   forwarded to providers that support typed responses.
#' @param tools Optional list of tool definitions forwarded to providers that
#'   support tool registration or tool calling.
#' @param ... Additional arguments passed to \code{provider}.
#' @importFrom utils modifyList
#'
#' @return A list containing at least \code{text}.
#'
#' @examples
#' provider <- function(messages, model = NULL, ...) {
#'   list(text = paste("received", length(messages), "messages"))
#' }
#' sn_run_llm(
#'   messages = list(list(role = "user", content = "Summarize this result.")),
#'   provider = provider
#' )
#' @export
sn_run_llm <- function(messages, provider, model = NULL, structured_type = NULL, tools = NULL, ...) {
  if (!is.function(provider)) {
    stop("`provider` must be a function.")
  }

  provider_args <- list(
    messages = messages,
    model = model,
    structured_type = structured_type,
    tools = tools,
    ...
  )
  provider_formals <- names(formals(provider) %||% alist(... = ))
  if (!"..." %in% provider_formals) {
    provider_args <- provider_args[intersect(names(provider_args), provider_formals)]
  }

  response <- do.call(provider, provider_args)
  if (is.character(response) && length(response) == 1) {
    return(list(text = response, model = model))
  }
  if (is.list(response) && any(c("text", "structured") %in% names(response))) {
    if (!"text" %in% names(response) && "structured" %in% names(response)) {
      response$text <- tryCatch(
        jsonlite::toJSON(response$structured, auto_unbox = TRUE, null = "null", dataframe = "rows"),
        error = function(...) NULL
      )
    }
    response$text <- .sn_text_scalar(response$text)
    return(response)
  }

  stop("`provider` must return either a single string or a list containing `text` and/or `structured`.")
}

.sn_text_scalar <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.character(x) && length(x) >= 1L) {
    return(as.character(x[[1]]))
  }

  values <- unlist(x, recursive = TRUE, use.names = FALSE)
  values <- values[!is.na(values)]
  if (length(values) == 0L) {
    return(NULL)
  }
  as.character(values[[1]])
}

.sn_format_openai_base_url <- function(base_url,
                                       require_v1 = TRUE) {
  if (!nzchar(base_url)) {
    stop("`base_url` must be a non-empty URL.", call. = FALSE)
  }

  base_url <- sub("/+$", "", base_url)
  if (grepl("/(responses|chat/completions)$", base_url)) {
    return(base_url)
  }
  if (require_v1 && !grepl("/v1$", base_url)) {
    return(paste0(base_url, "/v1"))
  }
  base_url
}

.sn_extract_openai_response_text <- function(parsed,
                                             wire_api = c("responses", "chat_completions")) {
  wire_api <- match.arg(wire_api)

  if (identical(wire_api, "chat_completions")) {
    choices <- parsed$choices %||% NULL
    choice <- if (is.data.frame(choices)) {
      as.list(choices[1, , drop = FALSE])
    } else if (is.list(choices) && length(choices) > 0) {
      choices[[1]]
    } else {
      NULL
    }
    message <- choice$message %||% NULL
    text <- if (is.list(message)) {
      message$content %||% NULL
    } else if (is.data.frame(message) && "content" %in% colnames(message)) {
      message$content[[1]]
    } else {
      choice$content %||% choice$text %||% NULL
    }
    return(text)
  }

  output_items <- parsed$output %||% list()
  output_text <- NULL
  if (is.list(output_items) && length(output_items) > 0L) {
    for (item in output_items) {
      contents <- item$content %||% list()
      if (is.list(contents) && length(contents) > 0L) {
        for (content_item in contents) {
          if (identical(content_item$type %||% NULL, "output_text") && !is.null(content_item$text)) {
            output_text <- content_item$text
            break
          }
        }
      }
      if (!is.null(output_text)) {
        break
      }
      if (!is.null(item$text)) {
        output_text <- item$text
        break
      }
    }
  }

  parsed$output_text %||%
    output_text %||%
    parsed$text %||%
    NULL
}

#' Create an \pkg{ellmer}-backed provider for Shennong interpretation helpers
#'
#' This adapter is useful when you want to manage transport, streaming, or
#' future structured-output features through \pkg{ellmer} while keeping
#' Shennong's interpretation API unchanged.
#'
#' @param api_key API key. Defaults to \code{OPENAI_API_KEY}.
#' @param base_url Base URL of the API. Defaults to \code{OPENAI_BASE_URL},
#'   then \code{"https://api.openai.com/v1"}.
#' @param model Default model identifier.
#' @param provider_type One of \code{"openai_compatible"} or \code{"openai"}.
#' @param echo Echo mode forwarded to \pkg{ellmer}.
#' @param reasoning_effort Optional reasoning effort forwarded to compatible
#'   GPT-5 chat-completions providers.
#' @param api_args Optional named list appended to each API request.
#' @param retries Number of retry attempts for retryable upstream proxy / HTML
#'   response failures.
#' @param retry_delay_sec Delay between retry attempts in seconds.
#'
#' @return A provider function suitable for \code{sn_run_llm()}.
#'
#' @examples
#' \dontrun{
#' provider <- sn_make_ellmer_provider(
#'   api_key = Sys.getenv("OPENAI_API_KEY"),
#'   base_url = "https://api.catplot.org/v1",
#'   model = "gpt-5.4"
#' )
#' }
#' @export
sn_make_ellmer_provider <- function(api_key = .sn_default_llm_api_key(),
                                    base_url = .sn_default_llm_base_url(),
                                    model = NULL,
                                    provider_type = c("openai_compatible", "openai"),
                                    echo = c("none", "output", "all"),
                                    reasoning_effort = NULL,
                                    api_args = list(),
                                    retries = 2L,
                                    retry_delay_sec = 1) {
  provider_type <- match.arg(provider_type)
  echo <- match.arg(echo)

  check_installed("ellmer")
  if (!nzchar(api_key)) {
    stop("`api_key` is required. Set `OPENAI_API_KEY` or pass it explicitly.", call. = FALSE)
  }

  default_model <- model
  formatted_base <- .sn_format_openai_base_url(base_url)
  provider_api_args <- api_args
  if (!is.null(reasoning_effort) && nzchar(reasoning_effort)) {
    provider_api_args$reasoning_effort <- reasoning_effort
  }

  function(messages, model = NULL, structured_type = NULL, tools = NULL, ...) {
    system_messages <- vapply(
      Filter(function(message) identical(message$role, "system"), messages),
      function(message) as.character(message$content %||% ""),
      character(1)
    )
    system_prompt <- if (length(system_messages) > 0L) {
      paste(system_messages, collapse = "\n\n")
    } else {
      NULL
    }
    non_system_messages <- Filter(function(message) !identical(message$role, "system"), messages)
    prompt_text <- paste(
      vapply(
        non_system_messages,
        function(message) paste0(toupper(message$role), ":\n", as.character(message$content %||% "")),
        character(1)
      ),
      collapse = "\n\n"
    )

    credentials <- function() api_key
    attempts <- max(1L, as.integer(retries %||% 1L))
    last_error <- NULL

    for (attempt in seq_len(attempts)) {
      result <- tryCatch({
        chat <- if (identical(provider_type, "openai")) {
          ellmer::chat_openai(
            system_prompt = system_prompt,
            base_url = formatted_base,
            credentials = credentials,
            model = model %||% default_model,
            api_args = modifyList(provider_api_args, list(...)),
            echo = echo
          )
        } else {
          ellmer::chat_openai_compatible(
            base_url = formatted_base,
            name = "Shennong provider",
            system_prompt = system_prompt,
            credentials = credentials,
            model = model %||% default_model,
            api_args = modifyList(provider_api_args, list(...)),
            echo = echo
          )
        }

        if (length(tools %||% list()) > 0L && is.null(structured_type)) {
          for (current_tool in tools) {
            chat$register_tool(current_tool)
          }
        }

        if (!is.null(structured_type)) {
          structured <- chat$chat_structured(
            prompt_text,
            type = structured_type,
            echo = echo
          )
          return(list(
            text = jsonlite::toJSON(structured, auto_unbox = TRUE, null = "null", dataframe = "rows"),
            structured = structured,
            model = model %||% default_model,
            raw = NULL
          ))
        }

        text <- chat$chat(prompt_text)
        list(
          text = as.character(text),
          model = model %||% default_model,
          raw = NULL
        )
      }, error = identity)

      if (!inherits(result, "error")) {
        return(result)
      }

      last_error <- result
      if (attempt < attempts && .sn_is_retryable_llm_error(conditionMessage(result))) {
        Sys.sleep(retry_delay_sec)
        next
      }
      break
    }

    stop(.sn_harden_llm_error(conditionMessage(last_error)), call. = FALSE)
  }
}

#' Test whether an ellmer-backed LLM provider is reachable and usable
#'
#' @param name Optional label used in the returned summary table.
#' @param model Optional override model used for the test request.
#' @param prompt Prompt used for the connectivity check.
#' @param provider Optional provider function. When supplied, Shennong tests
#'   this function directly; otherwise it builds an ellmer-backed provider from
#'   environment variables.
#'
#' @return A one-row tibble summarizing the result.
#'
#' @examples
#' test_provider <- function(messages, model = NULL, ...) list(text = "OK", model = model %||% "demo")
#' sn_test_llm_provider(provider = test_provider)
#' @export
sn_test_llm_provider <- function(name = NULL,
                                 model = NULL,
                                 prompt = "Reply with exactly OK.",
                                 provider = NULL) {
  provider_name <- name %||% "ellmer"
  provider <- provider %||% .sn_get_default_ellmer_provider()

  started_at <- Sys.time()
  result <- tryCatch(
    sn_run_llm(
      messages = list(list(role = "user", content = prompt)),
      provider = provider,
      model = model
    ),
    error = identity
  )
  elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

  if (inherits(result, "error")) {
    return(tibble::tibble(
      provider = provider_name,
      ok = FALSE,
      model = model %||% NA_character_,
      elapsed_sec = elapsed,
      text = NA_character_,
      error = conditionMessage(result)
    ))
  }

  tibble::tibble(
    provider = provider_name,
    ok = TRUE,
    model = result$model %||% model %||% NA_character_,
    elapsed_sec = elapsed,
    text = result$text %||% NA_character_,
    error = NA_character_
  )
}
