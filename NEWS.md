# fmriqa 0.4.0
* Edge detection is now used to improve object position detection.
* Added run_fmriqa_glob function for the analysis of multiple files.
* Added combine_res_glob function to combine fmriqa result files.

# fmriqa 0.3.0
* Improved method for calculating the object center of gravity for automated
ROI placement.
* Added method to calcuate the maximum background (MBG) signal.
* Added an option (pix_dim) to override the image voxel dimensions.
* Mean background image added to output.
* Background TC plot added to output.

# fmriqa 0.2.0
* Added an argument for the polynomial detrending order and changed the default
from 2 to 3.
* Improved k-space spike detection plot limits for very stable data.
* COG calculations are now performed on a thresholded image to reduce bias from
inhomogenity.
* Added coveralls checking.
* Added a basic consistancy test.
* Added an introduction vignette.

# fmriqa 0.1.0
* First public release.
