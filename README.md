# hwglabr
### Hochwagen lab R package

Compilation of functions used frequently by Hochwagen lab members.
This is a work in progress for the development of a more complete and useful package.

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

#### Documentation

There is a (mostly up to date) [documentation website](http://www.nyu.edu/projects/hochwagen/hwglabr/).

Function documentation is accessible within R in the standard way, by typing `?function_name`. You can
also get the names of all included functions directly from within R using **`ls()`**:

``` r
library(hwglabr)
ls('package:hwglabr')
```