# pbdBASE 

* **Version:** 0.5-0
* **License:** [License](https://www.mozilla.org/MPL/2.0/)
* **Project home**: https://github.com/RBigData/pbdBASE
* **Bug reports**: https://github.com/RBigData/pbdBASE/issues


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

The development version is maintained on GitHub:

```r
remotes::install_github("RBigData/pbdBASE")
```

See the vignette for installation troubleshooting.
