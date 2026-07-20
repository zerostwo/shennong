# Return a stdio configuration for the Shennong MCP server

Return a stdio configuration for the Shennong MCP server

## Usage

``` r
sn_mcp_server_config()
```

## Value

A list containing the command and arguments needed to launch the bundled
read-only MCP server.

## Examples

``` r
sn_mcp_server_config()
#> $command
#> [1] "/opt/R/4.6.1/lib/R/bin/Rscript"
#> 
#> $args
#> [1] "-e"                        "Shennong::sn_mcp_server()"
#> 
#> $transport
#> [1] "stdio"
#> 
```
