# hwglabr
### Hochwagen lab R package

Compilation of functions used frequently by Hochwagen lab members.
Currently a starting point for the development of a more complete and useful package.

#### Installation

You can install the package directly from the source code on GitHub. For that you will need Hadley Wickham's `devtools` R package:
``` r
install.packages("devtools")
```

Once you have `devtools` you can install and load `hwglabr`:
``` r
devtools::install_github("luisvalesilva/hwglabr")
library(hwglabr)
```

#### Listing included functions

In order to get the names of all included functions directly from within R use the function **`ls()`**:

``` r
library(hwglabr)
ls('package:hwglabr')
```

