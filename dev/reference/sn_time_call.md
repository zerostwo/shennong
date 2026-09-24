# Time an R expression

Evaluates an expression once, optionally displays elapsed time, and
records it as an ad-hoc workflow when usage tracking is enabled. The
expression's value, visibility, warnings, and errors are preserved.

## Usage

``` r
sn_time_call(expr, label = "ad_hoc", display = TRUE, record = TRUE)
```

## Arguments

- expr:

  Code to evaluate.

- label:

  Short, non-sensitive label used in timing output and the local usage
  database. It must use letters, numbers, dots, underscores, or dashes.

- display:

  Show the elapsed-time message.

- record:

  Record the call when usage tracking is enabled.

## Value

The value of `expr`, with its visibility preserved.

## Examples

``` r
result <- sn_time_call(sum(seq_len(1000)), label = "sum-example")
```
