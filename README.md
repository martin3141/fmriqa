
<!-- README.md is generated from README.Rmd. Please edit that file -->
fmriqa
======

[![Travis build status](https://travis-ci.org/martin3141/fmriqa.svg?branch=master)](https://travis-ci.org/martin3141/fmriqa) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/martin3141/fmriqa?branch=master&svg=true)](https://ci.appveyor.com/project/martin3141/fmriqa) [![](http://cranlogs.r-pkg.org/badges/fmriqa)](http://cran.rstudio.com/web/packages/fmriqa/index.html) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/fmriqa)](https://cran.r-project.org/package=fmriqa) [![Coverage status](https://coveralls.io/repos/github/martin3141/fmriqa/badge.svg)](https://coveralls.io/r/martin3141/fmriqa?branch=master)

Overview
--------

The fmriqa package provides an implemenation of the fMRI quality assurance (QA) analysis protocol detailed by Friedman and Glover (2006) <doi:10.1002/jmri.20583>.

Installation
------------

You can install the stable version of fmriqa from CRAN:

``` r
install.packages("fmriqa", dependencies = TRUE)
```

Or the the development version from GitHub (requires `devtools` package):

``` r
install.packages("devtools")
devtools::install_github("martin3141/fmriqa")
```

Usage
-----

``` r
# load the package
library(fmriqa)

# get help on the options for run_qa
?run_fmriqa

# run the analysis - a file chooser will appear when a data_file argument is not given
run_fmriqa()
```

Real data example
-----------------

``` r
library(fmriqa)
fname <- system.file("extdata", "qa_data.nii.gz", package = "fmriqa")
res <- run_fmriqa(data_file = fname, gen_png = FALSE, gen_res_csv = FALSE, tr = 3)
```

    ## Reading data  : E:/Program Files/R/R-3.4.2/library/fmriqa/extdata/qa_data.nii.gz
    ## 
    ## Basic analysis parameters
    ## -------------------------
    ## X,Y dims      : 80x80
    ## Slices        : 1
    ## TR            : 3s
    ## Slice #       : 1
    ## ROI width     : 21
    ## Total vols    : 200
    ## Analysis vols : 198
    ## 
    ## QA metrics
    ## ----------
    ## Mean signal   : 2561.3
    ## STD           : 6.4
    ## Percent fluc  : 0.25
    ## Drift         : 2.06
    ## Drift fit     : 1.42
    ## SNR           : 141.3
    ## SFNR          : 137
    ## RDC           : 2.84
    ## TC outlier    : 2.66
    ## Spec outlier  : 5.12

Simulation example
------------------

``` r
library(fmriqa)
library(oro.nifti)
```

    ## oro.nifti 0.7.2

``` r
# generate random data
set.seed(1)
sim_data <- array(rnorm(80 * 80 * 1 * 100), dim = c(80, 80, 1, 100))
sim_data[20:60, 20:60, 1, ] <- sim_data[20:60, 20:60, 1, ] + 50
sim_nifti <- oro.nifti::as.nifti(sim_data)
fname <- tempfile()
writeNIfTI(sim_nifti, fname)

# perform qa
res <- run_fmriqa(fname)
```

    ## Reading data  : C:\Users\wilsonmp\AppData\Local\Temp\RtmpW8pDDc\file259c44a83006
    ## 
    ## Basic analysis parameters
    ## -------------------------
    ## X,Y dims      : 80x80
    ## Slices        : 1
    ## TR            : 1s
    ## Slice #       : 1
    ## ROI width     : 21
    ## Total vols    : 100
    ## Analysis vols : 98
    ## 
    ## QA metrics
    ## ----------
    ## Mean signal   : 50
    ## STD           : 0.05
    ## Percent fluc  : 0.09
    ## Drift         : 0.45
    ## Drift fit     : 0.06
    ## SNR           : 51.7
    ## SFNR          : 51.2
    ## RDC           : 21.77
    ## TC outlier    : 2.74
    ## Spec outlier  : 4.04

    ## 
    ## PNG report    : C:\Users\wilsonmp\AppData\Local\Temp\RtmpW8pDDc\file259c44a83006_qa_plot.png
    ## CSV results   : C:\Users\wilsonmp\AppData\Local\Temp\RtmpW8pDDc\file259c44a83006_qa_results.csv

``` r
res$snr
```

    ## [1] 51.74954

Plot output from real data showing RF spiking artifact
------------------------------------------------------

![](SPIKE_EG_qa_plot.png)
