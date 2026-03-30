# Run an LLM provider on a prepared message list

Run an LLM provider on a prepared message list

## Usage

``` r
sn_run_llm(
  messages,
  provider,
  model = NULL,
  structured_type = NULL,
  tools = NULL,
  ...
)
```

## Arguments

- messages:

  A message list, typically from
  [`sn_build_prompt()`](https://songqi.org/shennong/dev/reference/sn_build_prompt.md).

- provider:

  A user-supplied function that accepts `messages` and returns text or a
  list containing `text`.

- model:

  Optional model identifier passed through to the provider.

- structured_type:

  Optional structured-output schema or type object forwarded to
  providers that support typed responses.

- tools:

  Optional list of tool definitions forwarded to providers that support
  tool registration or tool calling.

- ...:

  Additional arguments passed to `provider`.

## Value

A list containing at least `text`.

## Examples

``` r
provider <- function(messages, model = NULL, ...) {
  list(text = paste("received", length(messages), "messages"))
}
sn_run_llm(
  messages = list(list(role = "user", content = "Summarize this result.")),
  provider = provider
)
#> $text
#> [1] "received 1 messages"
#> 
```
