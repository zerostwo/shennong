# Run the read-only Shennong MCP server over stdio

Starts a newline-delimited JSON-RPC 2.0 server implementing MCP
lifecycle, tool discovery, and tool calls. The server exposes package
metadata and documentation only; it does not execute arbitrary R code or
modify analysis files.

## Usage

``` r
sn_mcp_server(input = stdin(), output = stdout())
```

## Arguments

- input:

  Input connection. Defaults to standard input.

- output:

  Output connection. Defaults to standard output.

## Value

Invisibly returns `NULL` when the input stream closes.

## References

Model Context Protocol specification:
<https://modelcontextprotocol.io/specification/2025-11-25>.

## Examples

``` r
config <- sn_mcp_server_config()
config$transport
#> [1] "stdio"
if (FALSE) { # \dontrun{
sn_mcp_server()
} # }
```
