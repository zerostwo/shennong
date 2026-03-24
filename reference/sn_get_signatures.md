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
#>   [1] "MT-ATP6" "MT-ATP8" "MT-CO1"  "MT-CO2"  "MT-CO3"  "MT-CYB"  "MT-ND1" 
#>   [8] "MT-ND2"  "MT-ND3"  "MT-ND4"  "MT-ND4L" "MT-ND5"  "MT-ND6"  "MT-RNR1"
#>  [15] "MT-RNR2" "MT-TA"   "MT-TC"   "MT-TD"   "MT-TE"   "MT-TF"   "MT-TG"  
#>  [22] "MT-TH"   "MT-TI"   "MT-TK"   "MT-TL1"  "MT-TL2"  "MT-TM"   "MT-TN"  
#>  [29] "MT-TP"   "MT-TQ"   "MT-TR"   "MT-TS1"  "MT-TS2"  "MT-TT"   "MT-TV"  
#>  [36] "MT-TW"   "MT-TY"   "RPL10"   "RPL10A"  "RPL10L"  "RPL11"   "RPL12"  
#>  [43] "RPL13"   "RPL13A"  "RPL14"   "RPL15"   "RPL17"   "RPL18"   "RPL18A" 
#>  [50] "RPL19"   "RPL21"   "RPL22"   "RPL22L1" "RPL23"   "RPL23A"  "RPL24"  
#>  [57] "RPL26"   "RPL26L1" "RPL27"   "RPL27A"  "RPL28"   "RPL29"   "RPL3"   
#>  [64] "RPL30"   "RPL31"   "RPL32"   "RPL34"   "RPL35"   "RPL35A"  "RPL36"  
#>  [71] "RPL36A"  "RPL36AL" "RPL37"   "RPL37A"  "RPL38"   "RPL39"   "RPL39L" 
#>  [78] "RPL3L"   "RPL4"    "UBA52"   "RPL41"   "RPL5"    "RPL6"    "RPL7"   
#>  [85] "RPL7A"   "RPL7L1"  "RPL8"    "RPL9"    "RPLP0"   "RPLP1"   "RPLP2"  
#>  [92] "MRPL1"   "MRPL10"  "MRPL11"  "MRPL12"  "MRPL13"  "MRPL14"  "MRPL15" 
#>  [99] "MRPL16"  "MRPL17"  "MRPL18"  "MRPL19"  "MRPL2"   "MRPL20"  "MRPL21" 
#> [106] "MRPL22"  "MRPL23"  "MRPL24"  "MRPL27"  "MRPL28"  "MRPL3"   "MRPL30" 
#> [113] "MRPL32"  "MRPL33"  "MRPL34"  "MRPL35"  "MRPL36"  "MRPL37"  "MRPL38" 
#> [120] "MRPL39"  "MRPL4"   "MRPL40"  "MRPL41"  "MRPL42"  "MRPL43"  "MRPL44" 
#> [127] "MRPL45"  "MRPL46"  "MRPL47"  "MRPL48"  "MRPL49"  "MRPL50"  "MRPL51" 
#> [134] "MRPL52"  "MRPL53"  "MRPL54"  "MRPL55"  "MRPL57"  "MRPL58"  "MRPL9"  
#> [141] "MRPS10"  "MRPS11"  "MRPS12"  "MRPS14"  "MRPS15"  "MRPS16"  "MRPS17" 
#> [148] "MRPS18A" "MRPS18B" "MRPS18C" "MRPS2"   "MRPS21"  "MRPS22"  "MRPS23" 
#> [155] "MRPS24"  "MRPS25"  "MRPS26"  "MRPS27"  "MRPS28"  "DAP3"    "MRPS30" 
#> [162] "MRPS31"  "MRPS33"  "MRPS34"  "MRPS35"  "MRPS36"  "MRPS5"   "MRPS6"  
#> [169] "MRPS7"   "MRPS9"   "RPS10"   "RPS11"   "RPS12"   "RPS13"   "RPS14"  
#> [176] "RPS15"   "RPS15A"  "RPS16"   "RPS17"   "RPS18"   "RPS19"   "RPS2"   
#> [183] "RPS20"   "RPS21"   "RPS23"   "RPS24"   "RPS25"   "RPS26"   "RPS27"  
#> [190] "RPS27A"  "RPS27L"  "RPS28"   "RPS29"   "RPS3"    "FAU"     "RPS3A"  
#> [197] "RPS4X"   "RPS4Y1"  "RPS4Y2"  "RPS5"    "RPS6"    "RPS7"    "RPS8"   
#> [204] "RPS9"    "RPSA"   
```
