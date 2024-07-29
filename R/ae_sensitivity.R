#' @title Sensitivity analysis for Angrist and Evans (1998)
#' @description Evaluate the sensitivity of Angrist and Evans (1998) to the
#' monotonicity condition then create a figure illustrating it.
#'
#' Bias of Wald for LATE is:
#' (Pr(DF)/(Pr(DF) + Pr(CP))) * (LATE - DATE)
#' where
#' DATE = defier average treatment effect
#'
#' @inheritParams shared_savedir
#' @export
ae_sensitivity <- function(savedir = getwd()) {
  wald <- -.133
  fs <- .060
  se <- .026

  expand.grid(
    prdf = seq(from = 0, to = .06, by = .001),
    date = seq(from = -.24, to = .06, by = .06)
  ) %>%
    mutate(
      prcp = prdf + fs,
      ratio = prdf / (prcp - prdf),
      late = (wald + ratio * date) / (1 + ratio),
      bias = (prdf / (prcp - prdf)) * (late - date),
      limit = late + bias,
      bias_prop = bias / late,
      bias_perc = 100 * bias_prop
    ) -> df

  p <- ggplot(data = df, aes(x = prdf, y = late, color = as.factor(date))) +
    geom_hline(yintercept = wald + 2 * se, linetype = "dotted") +
    geom_hline(yintercept = wald - 2 * se, linetype = "dotted") +
    geom_line() +
    geom_point(
      data = data.frame(x = 0, y = wald),
      aes(x = x, y = y),
      color = "black", size = 2
    ) +
    labs(
      x = "Proportion of defiers",
      y = "LATE that would rationalize the Wald estimate",
      color = "Average defier effects"
    ) +
    base_theme()

  if (safesavedir(savedir)) {
    fn <- fs::path(savedir, "ae-sensitivity.pdf")
    ggsave(fn, width = 7, height = 5)
  }
  return(p)
}

globalVariables(c(
  "prdf", "prcp", "ratio", "bias", "bias_prop"
))
