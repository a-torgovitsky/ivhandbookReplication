safesavedir <- function(savedir) {
  if (is.null(savedir)) {
    return(FALSE)
  }

  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive = TRUE)
  }
  return(TRUE)
}

timestampdir <- function(topdir) {
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  savedir <- fs::path(topdir, timestamp)
  return(savedir)
}

format_numeric <- function(x) {
  formatted <- if_else(is.na(x), "", sprintf("%.3f", x))
  # Remove leading zero but keep the minus sign if the number is negative
  formatted <- gsub("^(-?)0", "\\1", formatted)
  return(formatted)
}

save_table <- function(tb, colspec, savedir, fnstub) {
  if (safesavedir(savedir)) {
    saveargs <- list(
      tab = tb,
      filename = paste0(fnstub, "-view"),
      positions = colspec,
      output_path = as.character(savedir),
      stand_alone = TRUE,
      compile_tex = TRUE
    )
    do.call(TexSave, saveargs)
    saveargs$filename <- fnstub
    saveargs$stand_alone <- FALSE
    saveargs$compile_tex <- FALSE
    do.call(TexSave, saveargs)
  }
}

#' @title Run everything (replicate)
#'
#' @param resultsdir Directory where results will be saved. This is the
#' top-level directory. A timestamped subdirectory will be created within it,
#' and then results will be saved inside that subdirectory.
#' @export
ivhandbook <- function(resultsdir = fs::path_wd("results")) {
  savedir <- timestampdir(resultsdir)

  message("Saving in ", savedir)
  message(rep("=", 80))

  message("Starting multivalued instrument weighting figure...")
  multi_weighting(savedir)

  message("Starting Angrist and Evans sensitivity figure...")
  ae_sensitivity(savedir)

  message("Starting average monotonicity figure...")
  average_monotonicity(savedir)

  message("Starting unordered treatments group table...")
  unordered_treatments(savedir)

  message("Starting Card data illustration...")
  set.seed(1984) # A fine vintage
  # card_app(savedir, nbs = 500, nsplits = 100)
  card_app(savedir, nbs = 2, nsplits = 2)

  message("Starting Gelbach data illustration...")
  gelbach_bounds(savedir)
  gelbach_comparisons(savedir)

  message(
    ">> Potential benign warnings:\n",
    ">> 1) Namespace overwrites because ddml uses AER and I use ivreg.\n",
    ">> 2) ddml warnings about trimming, which are to be expected.\n",
    ">> 3) ivmte warning about compatible solvers not installed.\n"
  )
}
