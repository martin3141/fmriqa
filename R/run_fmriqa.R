#' Run fMRI quality assurance procedure on a NIfTI data file
#'
#' @param data_file input data in nifti format, a file chooser will open if not set
#' @param roi_width roi analysis region in pixels (default=21)
#' @param slice_num slice number for analysis (default=middle slice)
#' @param skip number of initial volumes to exclude from the analysis (default=2)
#' @param tr override the TR detected from data (seconds)
#' @param poly_det_ord polynomial order used for detrending (default=3)
#' @param spike_detect generate k-space spike-detection plot (default=FALSE)
#' @param x_pos x position of ROI (default=center of gravity)
#' @param y_pos y position of ROI (default=center of gravity)
#' @param plot_title add a title to the png and pdf plots
#' @param last_vol last volume number to use in the analysis
#' @param gen_png output png plot (default=TRUE)
#' @param gen_res_csv output csv results (default=TRUE)
#' @param gen_pdf output pdf plot (default=FALSE)
#' @param gen_spec_csv output csv of spectral points (default=FALSE)
#' @param png_fname png plot filename
#' @param res_fname csv results filename
#' @param pdf_fname pdf plot filename
#' @param spec_fname csv spectral data filename
#' @param verbose provide text output while running (default=TRUE)
#' @return dataframe of QA metrics
#' @examples
#' fname <- system.file("extdata", "qa_data.nii.gz", package = "fmriqa")
#' res <- run_fmriqa(data_file = fname, gen_png = FALSE, gen_res_csv = FALSE, tr = 3)
#'
#' @import viridisLite
#' @import RNifti
#' @import ggplot2
#' @import reshape2
#' @import gridExtra
#' @import grid
#' @import tidyr
#' @import optparse
#' @import tcltk
#' @import pracma
#' @importFrom grDevices graphics.off pdf png
#' @importFrom stats fft mad poly quantile sd median
#' @importFrom utils write.csv
#' @export
run_fmriqa <- function(data_file = NULL, roi_width = 21, slice_num = NULL,
                   skip = 2, tr = NULL, poly_det_ord = 3, spike_detect = FALSE,
                   x_pos = NULL, y_pos = NULL, plot_title = NULL,
                   last_vol = NULL, gen_png = TRUE, gen_res_csv = TRUE,
                   gen_pdf = FALSE, gen_spec_csv = FALSE, png_fname = NULL,
                   res_fname = NULL, pdf_fname = NULL, spec_fname = NULL,
                   verbose = TRUE) {

  if (is.null(data_file)) {
    filters <- matrix(c("NIfTI", ".nii.gz", "NIfTI", ".nii",
                        "All files", "*"),
                      3, 2, byrow = TRUE)

    data_file <- tk_choose.files(caption = "Select nifti data file for analysis",
                                 multi = FALSE, filters = filters)

    if (length(data_file) == 0) {
      stop("Error : input file not given.")
    }
  }

  basename <- sub(".nii.gz$", "", data_file)
  basename <- sub(".nii$", "", basename)

  if (is.null(res_fname)) {
    csv_file <- paste(basename, "_qa_results.csv", sep = "")
  } else {
    csv_file <- res_fname
  }

  if (is.null(png_fname)) {
    png_file <- paste(basename, "_qa_plot.png", sep = "")
  } else {
    png_file <- png_fname
  }

  if (is.null(pdf_fname)) {
    pdf_file <- paste(basename, "_qa_plot.pdf", sep = "")
  } else {
    pdf_file <- pdf_fname
  }

  if (is.null(spec_fname)) {
    spec_file <- paste(basename, "_qa_spec.csv", sep = "")
  } else {
    spec_file <- spec_fname
  }

  #image_cols <- inferno(64)
  image_cols <- viridis(64)

  if (verbose) cat(paste("Reading data  : ", data_file, "\n\n", sep = ""))
  data <- readNifti(data_file)

  x_dim <- dim(data)[1]
  y_dim <- dim(data)[2]
  z_dim <- dim(data)[3]

  if (is.null(tr)) tr <- pixdim(data)[4]

  if (is.null(slice_num)) slice_num <- ceiling(dim(data)[3] / 2)

  if (is.null(last_vol)) {
    N <- dim(data)[4]
  } else {
    N <- last_vol
  }

  dyns <- N - skip
  t <- seq(from = 0, by = tr, length.out = dyns)
  #t_full <- seq(from = 0, by = tr, length.out = N)

  if (verbose) {
    cat("Basic analysis parameters\n")
    cat("-------------------------\n")
    cat(paste("X,Y dims      : ", x_dim, "x", y_dim, "\n", sep = ""))
    cat(paste("Slices        : ", z_dim, "\n", sep = ""))
    cat(paste("TR            : ", round(tr, 2), "s\n", sep = ""))
    cat(paste("Slice #       : ", slice_num, "\n", sep = ""))
    cat(paste("ROI width     : ", roi_width, "\n", sep = ""))
    cat(paste("Total vols    : ", dim(data)[4], "\n", sep = ""))
    cat(paste("Analysis vols : ", dyns, "\n", sep = ""))
  }

  # scale data
  # scl_slope <- dumpNifti(data)$scl_slope
  # data <- data * scl_slope

  # chop out the slice we will be working with
  data_raw <- data[,,slice_num, (skip + 1):N]

  # detrend data with polynomial
  X <- poly(1:dyns, poly_det_ord)[,]
  X <- cbind(1,X)
  data_detrend <- apply(data_raw, c(1,2), detrend_fast, X)
  data_detrend <- aperm(data_detrend, c(2,3,1))

  # calculate temporal fluctuation noise (TFN)
  TFN <- apply(data_detrend, c(1,2), sd)
  av_image <- apply(data_raw, c(1,2), mean)
  SFNR_full <- av_image / TFN

  # calc diff image
  odd_dynamics <- data_raw[,,c(TRUE, FALSE)]
  even_dynamics <- data_raw[,,c(FALSE, TRUE)]

  if (length(odd_dynamics) > length(even_dynamics)) {
    odd_dynamics <- odd_dynamics[,,-(dim(odd_dynamics)[3])]
    warning("Odd number of dynamic scans, removing last one for the odd even diff calculation.")
  }

  DIFF <- apply(odd_dynamics, c(1, 2), sum) - apply(even_dynamics, c(1, 2), sum)

  # flip lr direction
  # SFNR_full <- flipud(SFNR_full)
  # av_image <- flipud(av_image)
  # DIFF <- flipud(DIFF)
  # TFN <- flipud(TFN)

  # set na values to zero
  SFNR_full[is.na(SFNR_full)] <- 0

  # threshold the image to reduce inhomogenity for cog calc
  cog_image <- imager::threshold(imager::as.cimg(av_image))[,]

  if (is.null(x_pos)) {
    x_pos <- sum(array(1:x_dim, c(x_dim, y_dim)) * cog_image) / sum(cog_image)
    x_pos <- round(x_pos)
  }

  if (is.null(y_pos)) {
    y_pos <- sum(t(array(1:y_dim, c(y_dim, x_dim))) * cog_image) / sum(cog_image)
    y_pos <- round(y_pos)
  }

  # get ROI indices
  ROI_x <- get_pixel_range(x_pos, roi_width)
  ROI_y <- get_pixel_range(y_pos, roi_width)

  SFNR <- SFNR_full[ROI_x, ROI_y]
  av_SFNR <- mean(SFNR)

  DIFF_ROI <- DIFF[ROI_x, ROI_y]

  signal_summary_value <- mean(av_image[ROI_x, ROI_y])

  SNR <- signal_summary_value / sqrt((sd(DIFF_ROI) ^ 2) / dyns)

  slice_data_ROI <- data_raw[ROI_x, ROI_y,]
  mean_sig_intensity_t <- apply(slice_data_ROI, 3, mean)
  mean_sig_intensity <- mean(mean_sig_intensity_t)

  mean_sig_intensity_t_detrend <- detrend_fast(mean_sig_intensity_t, X)
  y_fit <- mean_sig_intensity_t - mean_sig_intensity_t_detrend
  residuals <- mean_sig_intensity_t - y_fit
  sd_roi <- sd(residuals)

  percent_fluc      <- 100.0 * sd_roi / mean_sig_intensity
  percent_drift_fit <- 100.0 * (max(y_fit) - min(y_fit)) / mean_sig_intensity
  percent_drift     <- 100.0 * (max(mean_sig_intensity_t) -
                                  min(mean_sig_intensity_t)) / mean_sig_intensity

  detrend_res <- mean_sig_intensity_t - y_fit
  zp <- 4
  spec <- Mod(fft(c(detrend_res, rep(0,(zp - 1) * dyns))))[1:(dyns * zp / 2)]

  max_spec_outlier <- max(spec) / mad(spec)

  # x <- 1:(zp * N / 2)
  t <- seq(from = 0, by = tr, length.out = dyns)
  vols <- seq(from = skip + 1, by = 1, length.out = dyns)
  freq <- seq(from = 0, to = (1 - 1/(zp * dyns / 2))/(tr * 2),
              length.out = zp * dyns / 2)

  # get a mean time course for each slice
  slice_tc <- apply(data[,,,(skip + 1):N, drop = FALSE], c(3, 4), mean)

  # detrend
  X <- poly(1:dyns, poly_det_ord)[,]
  X <- cbind(1, X)
  slice_tc_dt <- apply(slice_tc, 1, detrend_fast, X)

  max_tc_outlier <- max(abs(slice_tc_dt)) / mad(slice_tc_dt)

  # normalise
  # slice_tc_dt <- scale(slice_tc_dt, center = F)

  # calculate RDC
  CV <- vector(length = roi_width)
  CV_ideal <- vector(length = roi_width)
  for (n in (1:roi_width)) {
    x_range <- get_pixel_range(x_pos, n)
    y_range <- get_pixel_range(y_pos, n)

    slice_data_ROI <- data_raw[x_range, y_range,, drop = F]
    mean_sig_intensity_t <- apply(slice_data_ROI, 3, mean)

    mean_sig_intensity <- mean(mean_sig_intensity_t)

    # detrend
    X <- poly(1:dyns, poly_det_ord)[,]
    X <- cbind(1,X)
    mean_sig_intensity_t_dt <- detrend_fast(y = mean_sig_intensity_t, X = X)

    sd_sig_intensity <- sd(mean_sig_intensity_t_dt)
    CV[n] <- 100 * sd_sig_intensity / mean_sig_intensity
    CV_ideal[n] <- CV[1] / n
  }

  RDC <- CV[1] / CV[length(CV)]

  line1  <- paste("Mean signal   : ", round(mean_sig_intensity, 1), "\n", sep = "")
  line2  <- paste("STD           : ", round(sd_roi, 2), "\n", sep = "")
  line3  <- paste("Percent fluc  : ", round(percent_fluc, 2), "\n", sep = "")
  line4  <- paste("Drift         : ", round(percent_drift, 2), "\n", sep = "")
  line5  <- paste("Drift fit     : ", round(percent_drift_fit, 2), "\n", sep = "")
  line6  <- paste("SNR           : ", round(SNR, 1), "\n", sep = "")
  line7  <- paste("SFNR          : ", round(av_SFNR, 1), "\n", sep = "")
  line8  <- paste("RDC           : ", round(RDC, 2), "\n", sep = "")
  line9  <- paste("TC outlier    : ", round(max_tc_outlier, 2), "\n", sep = "")
  line10 <- paste("Spec outlier  : ", round(max_spec_outlier, 2), "\n", sep = "")

  if (verbose) {
    cat("\nQA metrics\n")
    cat("----------\n")
    cat(line1)
    cat(line2)
    cat(line3)
    cat(line4)
    cat(line5)
    cat(line6)
    cat(line7)
    cat(line8)
    cat(line9)
    cat(line10)
  }

  if (is.null(plot_title)) plot_title <- NA

  results_tab <- data.frame(data_file, title = plot_title,
                            mean_signal = mean_sig_intensity, std = sd_roi,
                            percent_fluc = percent_fluc, drift = percent_drift,
                            drift_fit = percent_drift_fit, snr = SNR,
                            sfnr = av_SFNR, rdc = RDC, tc_outlier = max_tc_outlier,
                            spec_outlier = max_spec_outlier)

  if (gen_res_csv) {
    write.csv(results_tab, csv_file, row.names = FALSE)
  }

  if (gen_spec_csv) {
    spec_out <- data.frame(t(spec))
    colnames(spec_out) <- freq
    spec_out <- cbind(data.frame(data_file, title = plot_title), spec_out)
    write.csv(spec_out, spec_file, row.names = FALSE)
  }

  # plotting stuff below
  if (gen_pdf | gen_png) {

  # spike detection plot
  if (spike_detect) {
    cat("\nCalculating k-space spike detection map...\n")
    # calc diff volumes
    diff_vols <- apply(data[,,,(skip + 1):N, drop = FALSE], c(1,2,3), diff)
    diff_vols <- aperm(diff_vols, c(2,3,4,1))

    # transform all slices into k-space
    diff_vols_fft <- apply(diff_vols, c(3,4), fft)
    dim(diff_vols_fft) <- dim(diff_vols)

    # calc the maximum slice projection in k-space
    max_slice_proj <- apply(abs(diff_vols_fft), c(1,2), max)
    max_slice_proj <- apply(apply(max_slice_proj, 1, fftshift), 1,
                            fftshift)

    #max_z <- max(max_slice_proj) / 4
    max_z <- mad(max_slice_proj) * 8 + median(max_slice_proj)
    max_slice_proj <- ifelse(max_slice_proj > max_z, max_z, max_slice_proj)
    max_slice_proj_plot <- ggplot(melt(max_slice_proj), aes(Var1, Var2, fill = value)) +
      geom_raster(interpolate = TRUE) +
      scale_fill_gradientn(colours = image_cols) +
      coord_fixed(ratio = 1) + labs(x = "",y = "", fill = "Intensity",
                                    title = "Max. proj. of k-space slice differences") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  }


  theme_set(theme_bw())

  marg <- theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

  raw_text <- paste(line1, line2, line3, line4, line5, line6, line7, line8, line9, line10,
                    sep = "")

  text <- textGrob(raw_text, x = 0.2, just = 0, gp = gpar(fontfamily = "mono", fontsize = 14))

  # these are to appease R checks
  Measured <- NULL
  Theoretical <- NULL
  group <- NULL
  roi_width_vec <- NULL
  fit <- NULL
  tc <- NULL
  Var1 <- NULL
  Var2 <- NULL
  value <- NULL

  # RDC plot
  rdc_df <- data.frame(roi_width_vec = 1:roi_width, Theoretical = CV_ideal, Measured = CV)
  rdc_df <- gather(rdc_df, group, CV, c(Measured, Theoretical))

  rdc_plot <- ggplot(rdc_df, aes(x = roi_width_vec, y = CV, colour = group)) + geom_line() +
    geom_point() + scale_x_log10(limits = c(1,100)) +
    scale_y_log10(limits = c(0.01,10), breaks = c(0.01,0.1,1,10)) +
    labs(y = "100*CV", x = "ROI width (pixels)", title = "RDC plot") + marg +
    theme(legend.position = c(0.8, 0.8)) + scale_color_manual(values = c("black","red"))

  tc_fit <- data.frame(t = vols, tc = mean_sig_intensity_t, fit = y_fit)
  tc_plot <- ggplot(tc_fit, aes(t)) + geom_line(aes(y = tc)) +
    geom_line(aes(y = fit), color = "red") +
    theme(legend.position = "none") +
    labs(y = "Intensity (a.u.)", x = "Time (volumes)",
         title = "Intensity drift plot") + marg


  spec_plot <- qplot(freq, spec, xlab = "Frequency (Hz)",
                     ylab = "Intensity (a.u.)", geom = "line",
                     main = "Fluctuation spectrum") + marg

  x_st = ROI_x[1]
  x_end = ROI_x[length(ROI_x)]
  y_st = ROI_y[1]
  y_end = ROI_y[length(ROI_y)]

  lcol <- "white"

  roi_a <- geom_segment(aes(x = x_st, xend = x_st, y = y_st, yend = y_end),
                        colour = lcol)

  roi_b <- geom_segment(aes(x = x_end, xend = x_end, y = y_st, yend = y_end),
                        colour = lcol)

  roi_c <- geom_segment(aes(x = x_st, xend = x_end, y = y_st, yend = y_st),
                        colour = lcol)

  roi_d <- geom_segment(aes(x = x_st, xend = x_end, y = y_end, yend = y_end),
                        colour = lcol)



  top_val <- quantile(SFNR_full,0.999)
  SFNR_full <- ifelse(SFNR_full > top_val, top_val, SFNR_full)

  sfnr_plot <- ggplot(melt(SFNR_full), aes(Var1, Var2, fill = value)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(colours = image_cols) +
    coord_fixed(ratio = 1) + labs(x = "", y = "", fill = "Intensity",
                                  title = "SFNR image") +
    marg + roi_a + roi_b + roi_c + roi_d +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

  # useful for checking where the ROI really is
  # av_image[ROI_x,ROI_y] = 0

  av_plot <- ggplot(melt(av_image), aes(Var1, Var2, fill = value)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(colours = image_cols) +
    coord_fixed(ratio = 1) + labs(x = "",y = "", fill = "Intensity",
                                  title = "Mean image") +
    marg + roi_a + roi_b + roi_c + roi_d +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

  diff_plot <- ggplot(melt(DIFF), aes(Var1, Var2, fill = value)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(colours = image_cols) +
    coord_fixed(ratio = 1) + labs(x = "",y = "", fill = "Intensity",
                                  title = "Odd-even difference") +
    marg + roi_a + roi_b + roi_c + roi_d +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

  tfn_plot <- ggplot(melt(TFN), aes(Var1, Var2, fill = value)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(colours = image_cols) +
    coord_fixed(ratio = 1) +
    labs(x = "", y = "", fill = "Intensity",
         title = "Temporal fluctuation noise") +
    marg + roi_a + roi_b + roi_c + roi_d +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

  slice_tc_plot <- ggplot(melt(slice_tc_dt), aes(x = Var1 + skip, y = value, group = Var2)) +
    geom_line(alpha = 0.5) +
    labs(x = "Time (volumes)", y = "Intensity (a.u.)",
         title = "Mean slice TCs (detrended)")

  if (is.na(plot_title)) {
    title <- NULL
  } else {
    title <- textGrob(plot_title, gp = gpar(fontsize = 25))
  }

  if (gen_pdf) {
    pdf(pdf_file, height = 10, width = 16)
    if (spike_detect) {
      lay <- rbind(c(1,5,6,9),
                   c(4,2,2,3),
                   c(7,8,8,10))
      grid.arrange(text, tc_plot, spec_plot, av_plot, diff_plot, sfnr_plot, tfn_plot,
                   slice_tc_plot, rdc_plot, max_slice_proj_plot, layout_matrix = lay,
                   top = title)
    } else {
      lay <- rbind(c(1,5,6,9),
                   c(4,2,2,3),
                   c(7,8,8,8))
      grid.arrange(text, tc_plot, spec_plot, av_plot, diff_plot, sfnr_plot, tfn_plot,
                   slice_tc_plot, rdc_plot, layout_matrix = lay,
                   top = title)
    }
    graphics.off()
  }

  if (gen_png) {
    png(png_file, height = 800, width = 1200, type = "cairo")
    if (spike_detect) {
      lay <- rbind(c(1,5,6,9),
                   c(4,2,2,3),
                   c(7,8,8,10))
      grid.arrange(text, tc_plot, spec_plot, av_plot, diff_plot, sfnr_plot, tfn_plot,
                   slice_tc_plot, rdc_plot, max_slice_proj_plot, layout_matrix = lay,
                   top = title)
    } else {
      lay <- rbind(c(1,5,6,9),
                   c(4,2,2,3),
                   c(7,8,8,8))
      grid.arrange(text, tc_plot, spec_plot, av_plot, diff_plot, sfnr_plot, tfn_plot,
                   slice_tc_plot, rdc_plot, layout_matrix = lay,
                   top = title)
    }
    graphics.off()
  }
  }

  # end of plotting

  if (verbose) {
    if (gen_pdf)      cat(paste("\nPDF report    : ", pdf_file, sep = ""))
    if (gen_spec_csv) cat(paste("\nCSV spec file : ", spec_file, sep = ""))
    if (gen_png)      cat(paste("\nPNG report    : ", png_file, sep = ""))
    if (gen_res_csv)  cat(paste("\nCSV results   : ", csv_file, "\n\n", sep = ""))
  }

  results_tab
}

#' @import RcppEigen
detrend_fast <- function(y, X) {
  fastLmPure(y = y, X = X)$residual
}

get_pixel_range <- function(center, width) {
  ROI_half_start <- floor(width / 2)
  ROI_half_end <- ceiling(width / 2)
  start <- floor(center - ROI_half_start)
  end <- floor(center + ROI_half_end) - 1
  start:end
}
