# Test whether an ellmer-backed LLM provider is reachable and usable

Test whether an ellmer-backed LLM provider is reachable and usable

## Usage

``` r
sn_test_llm_provider(
  name = NULL,
  model = NULL,
  prompt = "Reply with exactly OK.",
  provider = NULL
)
```

## Arguments

- name:

  Optional label used in the returned summary table.

- model:

  Optional override model used for the test request.

- prompt:

  Prompt used for the connectivity check.

- provider:

  Optional provider function. When supplied, Shennong tests this
  function directly; otherwise it builds an ellmer-backed provider from
  environment variables.

## Value

A one-row tibble summarizing the result.

## Examples

``` r
test_provider <- function(messages, model = NULL, ...) list(text = "OK", model = model %||% "demo")
sn_test_llm_provider(provider = test_provider)
#> # A tibble: 1 × 6
#>   provider ok    model elapsed_sec text  error
#>   <chr>    <lgl> <chr>       <dbl> <chr> <chr>
#> 1 ellmer   TRUE  demo    0.0000448 OK    NA   
```
