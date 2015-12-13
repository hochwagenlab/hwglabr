# hwglabr
### Hochwagen lab R package

Compilation of functions written and used frequently by Hochwagen
lab members.

This is a starting point for the development of a more
complete and useful package.

Installation
------------
------------

You can install the package directly from the source code on GitHub. For that you will need Hadley Wickham's `devtools` R package:
``` r
install.packages("devtools")
library(devtools)
```

Once you have `devtools` you can install and load `hwglabr`:
``` r
install_github("luisvalesilva/hwglabr")
library(hwglabr)
```

Listing included functions
--------------------------
--------------------------

In order to get the names of all included packages directly from within R just use the function **`ls()`**:

``` r
library(hwglabr)
ls('package:hwglabr')
```

