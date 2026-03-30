# Resolve a palette into explicit colors

Resolve a palette into explicit colors

## Usage

``` r
sn_get_palette(
  palette = "Paired",
  n = NULL,
  palette_type = c("auto", "discrete", "continuous"),
  direction = 1
)
```

## Arguments

- palette:

  Palette name or explicit character vector of colors.

- n:

  Number of colors to return. When omitted, the palette's native length
  is returned for discrete use and `256` colors are returned for
  continuous use.

- palette_type:

  One of `"auto"`, `"discrete"`, or `"continuous"`. Defaults to
  `"auto"`.

- direction:

  Direction for ordered palettes. Use `1` for the default order and `-1`
  to reverse it.

## Value

A character vector of hex colors.

## Examples

``` r
sn_get_palette("Paired", n = 14)
#>  [1] "#A6CEE3" "#3385BB" "#84BF96" "#6DBD57" "#7F9D55" "#F57C7C"
#>  [7] "#E42622" "#FBB268" "#FE8D19" "#DE9E83" "#9D7BBA" "#926B8F"
#> [13] "#E2C16A" "#B15928"
sn_get_palette("RdBu", palette_type = "continuous", direction = -1)
#>   [1] "#053061" "#063263" "#073466" "#083669" "#09386C" "#0A3A6F"
#>   [7] "#0B3C72" "#0C3E75" "#0D4078" "#0E437B" "#0F457E" "#114781"
#>  [13] "#124984" "#134B87" "#144D8A" "#154F8D" "#165190" "#175493"
#>  [19] "#185695" "#195898" "#1A5A9B" "#1C5C9E" "#1D5EA1" "#1E60A4"
#>  [25] "#1F62A7" "#2064AA" "#2166AC" "#2368AD" "#246AAE" "#256CAF"
#>  [31] "#276DB0" "#286FB0" "#2971B1" "#2B73B2" "#2C75B3" "#2D76B4"
#>  [37] "#2F78B5" "#307AB6" "#317CB7" "#337DB8" "#347FB9" "#3581B9"
#>  [43] "#3783BA" "#3884BB" "#3986BC" "#3B88BD" "#3C8ABE" "#3D8BBF"
#>  [49] "#3F8DC0" "#408FC1" "#4191C2" "#4393C3" "#4694C4" "#4996C5"
#>  [55] "#4C98C6" "#4F9AC7" "#529CC8" "#559EC9" "#58A0CA" "#5BA2CB"
#>  [61] "#5EA4CC" "#61A6CD" "#65A8CE" "#68AACF" "#6BACD0" "#6EAED1"
#>  [67] "#71B0D2" "#74B2D3" "#77B4D5" "#7AB6D6" "#7DB8D7" "#80BAD8"
#>  [73] "#84BCD9" "#87BEDA" "#8AC0DB" "#8DC2DC" "#90C4DD" "#93C5DE"
#>  [79] "#95C6DF" "#98C8DF" "#9AC9E0" "#9DCAE1" "#9FCBE1" "#A2CDE2"
#>  [85] "#A4CEE3" "#A7CFE4" "#A9D0E4" "#ABD2E5" "#AED3E6" "#B0D4E6"
#>  [91] "#B3D5E7" "#B5D7E8" "#B8D8E8" "#BAD9E9" "#BDDAEA" "#BFDCEB"
#>  [97] "#C2DDEB" "#C4DEEC" "#C7DFED" "#C9E1ED" "#CCE2EE" "#CEE3EF"
#> [103] "#D1E5F0" "#D2E5F0" "#D3E6F0" "#D5E7F0" "#D6E7F1" "#D8E8F1"
#> [109] "#D9E9F1" "#DBE9F1" "#DCEAF2" "#DEEBF2" "#DFECF2" "#E1ECF3"
#> [115] "#E2EDF3" "#E4EEF3" "#E5EEF3" "#E7EFF4" "#E8F0F4" "#EAF1F4"
#> [121] "#EBF1F4" "#EDF2F5" "#EEF3F5" "#F0F3F5" "#F1F4F6" "#F3F5F6"
#> [127] "#F4F5F6" "#F6F6F6" "#F7F6F6" "#F7F5F4" "#F7F4F2" "#F7F3F0"
#> [133] "#F8F2EE" "#F8F0EC" "#F8EFEA" "#F8EEE8" "#F9EDE7" "#F9ECE5"
#> [139] "#F9EBE3" "#F9EAE1" "#F9E9DF" "#FAE8DD" "#FAE7DB" "#FAE5D9"
#> [145] "#FAE4D7" "#FBE3D6" "#FBE2D4" "#FBE1D2" "#FBE0D0" "#FCDFCE"
#> [151] "#FCDECC" "#FCDDCA" "#FCDCC8" "#FDDBC7" "#FCD8C4" "#FCD6C1"
#> [157] "#FBD4BE" "#FBD2BC" "#FBD0B9" "#FACEB6" "#FACCB4" "#FACAB1"
#> [163] "#F9C7AE" "#F9C5AB" "#F9C3A9" "#F8C1A6" "#F8BFA3" "#F8BDA1"
#> [169] "#F7BB9E" "#F7B99B" "#F7B698" "#F6B496" "#F6B293" "#F5B090"
#> [175] "#F5AE8E" "#F5AC8B" "#F4AA88" "#F4A886" "#F4A683" "#F3A380"
#> [181] "#F2A07E" "#F19E7C" "#EF9B7A" "#EE9878" "#ED9676" "#EC9374"
#> [187] "#EB9072" "#EA8D70" "#E88B6E" "#E7886C" "#E6856A" "#E58368"
#> [193] "#E48065" "#E27D63" "#E17B61" "#E0785F" "#DF755D" "#DE725B"
#> [199] "#DD7059" "#DB6D57" "#DA6A55" "#D96853" "#D86551" "#D7624F"
#> [205] "#D6604D" "#D45D4B" "#D35A4A" "#D15749" "#D05447" "#CE5146"
#> [211] "#CD4F44" "#CC4C43" "#CA4942" "#C94641" "#C7433F" "#C6403E"
#> [217] "#C53E3C" "#C33B3B" "#C2383A" "#C03538" "#BF3237" "#BE2F36"
#> [223] "#BC2D34" "#BB2A33" "#B92732" "#B82430" "#B6212F" "#B51F2E"
#> [229] "#B41C2D" "#B2192B" "#B0172A" "#AD162A" "#AA1529" "#A71429"
#> [235] "#A41328" "#A11228" "#9E1127" "#9B1027" "#991027" "#960F26"
#> [241] "#930E26" "#900D25" "#8D0C25" "#8A0B24" "#870A24" "#840923"
#> [247] "#810823" "#7E0722" "#7B0622" "#780521" "#750421" "#720320"
#> [253] "#6F0220" "#6C011F" "#69001F" "#67001F"
```
