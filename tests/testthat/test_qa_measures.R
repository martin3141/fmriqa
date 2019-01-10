context("qa_measures")

test_that("default analysis is consistent", {
  fname <- system.file("extdata", "qa_data.nii.gz", package = "fmriqa")
  res <- run_fmriqa(data_file = fname, tr = 3, pix_dim = c(3,3,3),
                    png_fname = tempfile(), res_fname = tempfile(),
                    spike_detect = TRUE)
  res$data_file <- NULL
  expect_equal_to_reference(res, "ref_results.rds", tolerance = 1e-6)
})

test_that("pdf output runs", {
  fname <- system.file("extdata", "qa_data.nii.gz", package = "fmriqa")
  res <- run_fmriqa(data_file = fname, tr = 3, pix_dim = c(3,3,3),
                    pdf_fname = tempfile(), gen_pdf = TRUE, gen_res_csv = FALSE,
                    gen_png = FALSE)
})

test_that("spec csv output runs", {
  fname <- system.file("extdata", "qa_data.nii.gz", package = "fmriqa")
  res <- run_fmriqa(data_file = fname, tr = 3, pix_dim = c(3,3,3),
                    spec_fname = tempfile(), gen_spec_csv = TRUE,
                    gen_res_csv = FALSE, gen_png = FALSE)
})
