# Report the status of a registered Shennong method

Report the status of a registered Shennong method

## Usage

``` r
sn_get_method_status(method, task = NULL)
```

## Arguments

- method:

  A registered method name.

- task:

  Optional task used to disambiguate a method registered for more than
  one workflow.

## Value

A named list describing availability, runtime, environment, installation
action, requirements, outputs, and citation.

## Examples

``` r
sn_get_method_status("slingshot", task = "trajectory")
```
