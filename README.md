# hwglabr
### Hochwagen lab R package

Compilation of functions written and used frequently by Hochwagen
lab members.

This is a starting point for the development of a more
complete and useful package.

Installation
------------

``` r
# Install Hadley Whickam's `devtools` R package
install.packages("devtools")
library(devtools)

# Install `hwglabr` from the source code on GitHub
install_github("luisvalesilva/hwglabr")

# Load and use :-)
library(hwglabr)
```

Listing included functions
--------------------------

In order to get the names of all included packages directly from within R just use the function **`ls()`**:

``` r
library(hwglabr)
ls('package:hwglabr')
```

