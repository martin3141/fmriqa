Changed doi format in the vignette.

## Test environments

* Windows 7, R 3.4.2
* Linux CentOS 6, R 3.3.3
* OS X 10.9.5, R 3.3.2

## R CMD check results

0 errors | 0 warnings | 1 note

"Package in Depends/Imports which should probably only be in 
LinkingTo: 'RcppEigen'" - this package uses the RcppEigen::fastLmPure function
and therefore RcppEigen is an import.

## Downstream dependancies

There are currently no downstream dependancies for this package.
