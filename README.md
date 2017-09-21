[![Travis-CI Build Status](https://travis-ci.org/hochwagenlab/hwglabr.svg?branch=master)](https://travis-ci.org/hochwagenlab/hwglabr)

# hwglabr
### Hochwagen lab R package

Compilation of functions used frequently by Hochwagen lab members.

#### Installation

You can install the package directly from the source code on GitHub. For that you
will need Hadley Wickham's `devtools` R package:
``` r
install.packages("devtools")
```

Once you have `devtools` you can install and load `hwglabr`:
``` r
devtools::install_github("hochwagenlab/hwglabr")
library(hwglabr)
```

#### Documentation

Use the package GitHub [repo](https://github.com/hochwagenlab/hwglabr) and the
[documentation website](https://hochwagenlab.github.io/hwglabr/).
Function documentation is accessible within R in the standard way, by typing one of the following:

``` r
help("function_name")

?function_name
```

You can also get the list of included functions directly within R using one of the following:

``` r
ls("package:hwglabr")   # List function names

lsf.str("package:hwglabr")   # List function names and their arguments
```
