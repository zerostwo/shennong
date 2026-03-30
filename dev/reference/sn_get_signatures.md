# Retrieve bundled Shennong signature genes by category or tree path

Shennong ships a package-owned signature catalog under `data/`, built
from the upstream SignatuR tree plus any package-maintained additions.
Queries can use either leaf names such as `"mito"` or full tree paths
such as `"Compartments/Mito"`.

## Usage

``` r
sn_get_signatures(species = "human", category = NULL)
```

## Arguments

- species:

  One of `"human"` or `"mouse"`.

- category:

  A character vector of signature names or full tree paths.

## Value

A unique character vector of signature gene symbols.

## Examples

``` r
sn_get_signatures(
  species = "human",
  category = c("mito", "Compartments/Ribo")
)
#>   [1] "MT-ATP6" "MT-ATP8" "MT-CO1"  "MT-CO2"  "MT-CO3"  "MT-CYB" 
#>   [7] "MT-ND1"  "MT-ND2"  "MT-ND3"  "MT-ND4"  "MT-ND4L" "MT-ND5" 
#>  [13] "MT-ND6"  "MT-RNR1" "MT-RNR2" "MT-TA"   "MT-TC"   "MT-TD"  
#>  [19] "MT-TE"   "MT-TF"   "MT-TG"   "MT-TH"   "MT-TI"   "MT-TK"  
#>  [25] "MT-TL1"  "MT-TL2"  "MT-TM"   "MT-TN"   "MT-TP"   "MT-TQ"  
#>  [31] "MT-TR"   "MT-TS1"  "MT-TS2"  "MT-TT"   "MT-TV"   "MT-TW"  
#>  [37] "MT-TY"   "RPL10"   "RPL10A"  "RPL10L"  "RPL11"   "RPL12"  
#>  [43] "RPL13"   "RPL13A"  "RPL14"   "RPL15"   "RPL17"   "RPL18"  
#>  [49] "RPL18A"  "RPL19"   "RPL21"   "RPL22"   "RPL22L1" "RPL23"  
#>  [55] "RPL23A"  "RPL24"   "RPL26"   "RPL26L1" "RPL27"   "RPL27A" 
#>  [61] "RPL28"   "RPL29"   "RPL3"    "RPL30"   "RPL31"   "RPL32"  
#>  [67] "RPL34"   "RPL35"   "RPL35A"  "RPL36"   "RPL36A"  "RPL36AL"
#>  [73] "RPL37"   "RPL37A"  "RPL38"   "RPL39"   "RPL39L"  "RPL3L"  
#>  [79] "RPL4"    "UBA52"   "RPL41"   "RPL5"    "RPL6"    "RPL7"   
#>  [85] "RPL7A"   "RPL7L1"  "RPL8"    "RPL9"    "RPLP0"   "RPLP1"  
#>  [91] "RPLP2"   "MRPL1"   "MRPL10"  "MRPL11"  "MRPL12"  "MRPL13" 
#>  [97] "MRPL14"  "MRPL15"  "MRPL16"  "MRPL17"  "MRPL18"  "MRPL19" 
#> [103] "MRPL2"   "MRPL20"  "MRPL21"  "MRPL22"  "MRPL23"  "MRPL24" 
#> [109] "MRPL27"  "MRPL28"  "MRPL3"   "MRPL30"  "MRPL32"  "MRPL33" 
#> [115] "MRPL34"  "MRPL35"  "MRPL36"  "MRPL37"  "MRPL38"  "MRPL39" 
#> [121] "MRPL4"   "MRPL40"  "MRPL41"  "MRPL42"  "MRPL43"  "MRPL44" 
#> [127] "MRPL45"  "MRPL46"  "MRPL47"  "MRPL48"  "MRPL49"  "MRPL50" 
#> [133] "MRPL51"  "MRPL52"  "MRPL53"  "MRPL54"  "MRPL55"  "MRPL57" 
#> [139] "MRPL58"  "MRPL9"   "MRPS10"  "MRPS11"  "MRPS12"  "MRPS14" 
#> [145] "MRPS15"  "MRPS16"  "MRPS17"  "MRPS18A" "MRPS18B" "MRPS18C"
#> [151] "MRPS2"   "MRPS21"  "MRPS22"  "MRPS23"  "MRPS24"  "MRPS25" 
#> [157] "MRPS26"  "MRPS27"  "MRPS28"  "DAP3"    "MRPS30"  "MRPS31" 
#> [163] "MRPS33"  "MRPS34"  "MRPS35"  "MRPS36"  "MRPS5"   "MRPS6"  
#> [169] "MRPS7"   "MRPS9"   "RPS10"   "RPS11"   "RPS12"   "RPS13"  
#> [175] "RPS14"   "RPS15"   "RPS15A"  "RPS16"   "RPS17"   "RPS18"  
#> [181] "RPS19"   "RPS2"    "RPS20"   "RPS21"   "RPS23"   "RPS24"  
#> [187] "RPS25"   "RPS26"   "RPS27"   "RPS27A"  "RPS27L"  "RPS28"  
#> [193] "RPS29"   "RPS3"    "FAU"     "RPS3A"   "RPS4X"   "RPS4Y1" 
#> [199] "RPS4Y2"  "RPS5"    "RPS6"    "RPS7"    "RPS8"    "RPS9"   
#> [205] "RPSA"   
```
