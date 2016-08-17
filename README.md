# pbdBASE 

* **Version:** 0.4-5
* **License:** [![License](http://img.shields.io/badge/license-MPL%202-orange.svg?style=flat)](https://www.mozilla.org/MPL/2.0/)
* **Author:** See section below.


pbdBASE is a set of bindings to and extensions for the distributed linear algebra libraries BLACS, PBLAS, and ScaLAPACK.  The package is very low-level, and unless you are very familiar with these libraries (or even if you are...), you are instead recommended to see the pbdDMAT and pbdDEMO packages.



## Installation

pbdBASE requires:

* A system installation of MPI
* R version 3.0.0 or higher
* The pbdSLAP and pbdMPI packages, as well as their dependencies.

Assuming you meet the system dependencies, you can install the stable version from CRAN using the usual `install.packages()`:

```r
install.package("pbdBASE")
```

The development version is maintained on GitHub, and can easily be installed by any of the packages that offer installations from GitHub:

```r
### Pick your preference
devtools::install_github("RBigData/pbdBASE")
ghit::install_github("RBigData/pbdBASE")
remotes::install_github("RBigData/pbdBASE")
```

See the vignette for installation troubleshooting.



## Authors

pbdBASE is authored and maintained by the pbdR core team:

* Drew Schmidt
* Wei-Chen Chen
* George Ostrouchov
* Pragneshkumar Patel

With additional contributions from:

* Ewan Higgs
