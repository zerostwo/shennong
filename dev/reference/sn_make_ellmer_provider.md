# Create an ellmer-backed provider for Shennong interpretation helpers

This adapter is useful when you want to manage transport, streaming, or
future structured-output features through ellmer while keeping
Shennong's interpretation API unchanged.

## Usage

``` r
sn_make_ellmer_provider(
  api_key = .sn_default_llm_api_key(),
  base_url = .sn_default_llm_base_url(),
  model = NULL,
  provider_type = c("openai_compatible", "openai"),
  echo = c("none", "output", "all"),
  reasoning_effort = NULL,
  api_args = list(),
  retries = 2L,
  retry_delay_sec = 1
)
```

## Arguments

- api_key:

  API key. Defaults to `OPENAI_API_KEY`.

- base_url:

  Base URL of the API. Defaults to `OPENAI_BASE_URL`, then
  `"https://api.openai.com/v1"`.

- model:

  Default model identifier.

- provider_type:

  One of `"openai_compatible"` or `"openai"`.

- echo:

  Echo mode forwarded to ellmer.

- reasoning_effort:

  Optional reasoning effort forwarded to compatible GPT-5
  chat-completions providers.

- api_args:

  Optional named list appended to each API request.

- retries:

  Number of retry attempts for retryable upstream proxy / HTML response
  failures.

- retry_delay_sec:

  Delay between retry attempts in seconds.

## Value

A provider function suitable for
[`sn_run_llm()`](https://songqi.org/shennong/dev/reference/sn_run_llm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
provider <- sn_make_ellmer_provider(
  api_key = Sys.getenv("OPENAI_API_KEY"),
  base_url = "https://api.catplot.org/v1",
  model = "gpt-5.4"
)
} # }
```
