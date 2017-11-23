
<!-- README.md is generated from README.Rmd. Please edit that file -->
fmriqa
======

Overview
--------

The fmriqa package provides an implemenation of the fMRI QA analysis protocol detailed in:

Friedman L and Glover GH. Report on a multicenter fMRI quality assurance protocol. J Magn Reson Imaging. 2006 Jun;23(6):827-39.

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
?run_qa

# run the analysis - a file chooser will appear when a data_file argument is not given
run_qa()
```

Simualated example
------------------

``` r
library(fmriqa)
library(oro.nifti)
```

    ## oro.nifti 0.7.2

``` r
# generate random data
set.seed(1)
sim_data <- array(rnorm(80 * 80 * 11 * 100), dim = c(80, 80, 11, 100))
sim_data[20:60, 20:60, 6, ] <- sim_data[20:60, 20:60, 6, ] + 50
sim_nifti <- oro.nifti::as.nifti(sim_data)
fname <- tempfile()
writeNIfTI(sim_nifti, fname)

# perform qa
res <- run_qa(fname)
```

    ## Reading data  : C:\Users\wilsonmp\AppData\Local\Temp\RtmpELqDND\filec7c66c831c1
    ## 
    ## Basic analysis parameters
    ## -------------------------
    ## X,Y dims      : 80x80
    ## Slices        : 11
    ## TR            : 1s
    ## Slice #       : 6
    ## ROI width     : 21
    ## Total vols    : 100
    ## Analysis vols : 98
    ## 
    ## fBIRN QA metrics
    ## ----------------
    ## Mean signal   : 50
    ## STD           : 0.05
    ## Percent fluc  : 0.1
    ## Drift         : 0.53
    ## Drift fit     : 0.03
    ## SNR           : 51.6
    ## SFNR          : 50.7
    ## RDC           : 20.32
    ## TC outlier    : 2.87
    ## Spec outlier  : 5.99

    ## 
    ## PNG report    : C:\Users\wilsonmp\AppData\Local\Temp\RtmpELqDND\filec7c66c831c1_qa_plot.png
    ## CSV results   : C:\Users\wilsonmp\AppData\Local\Temp\RtmpELqDND\filec7c66c831c1_qa_results.csv

``` r
res$snr
```

    ## [1] 51.59899

Plot output from real data showing RF spiking artifact
------------------------------------------------------

![](SPIKE_EG_qa_plot.png)
