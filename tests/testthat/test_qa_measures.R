context("qa_measures")

test_that("default analysis is consistent", {
  fname <- system.file("extdata", "qa_data.nii.gz", package = "fmriqa")
  res <- run_fmriqa(data_file = fname, gen_png = FALSE, gen_res_csv = FALSE,
                    tr = 3)
  expect_equal_to_reference(res, "ref_results.rds", tolerance = 1e-6)
})
