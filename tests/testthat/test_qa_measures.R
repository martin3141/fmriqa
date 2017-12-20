context("qa_measures")

test_that("default analysis is consistent", {
  fname <- system.file("extdata", "qa_data.nii.gz", package = "fmriqa")
  res <- run_fmriqa(data_file = fname, tr = 3, png_fname = tempfile(),
                    res_fname = tempfile(), spike_detect = TRUE)
  res$data_file <- NULL
  expect_equal_to_reference(res, "ref_results.rds", tolerance = 1e-6)
})
