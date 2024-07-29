#' @title Illustrate average monotonicity
#'
#' @inheritParams shared_savedir
#' @export
average_monotonicity <- function(savedir = getwd()) {
  dg <- average_monotonicity_df()
  nj <- log(nrow(dg), base = 2)
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # Prepare annotation
  #
  # Find example groups to highlight in the figure
  # Pull out a dataframe that contains just those groups
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  XLOC_EXGR <- -.0045
  dg %>%
    mutate(dist = abs(abs(gw) - abs(XLOC_EXGR))) %>%
    arrange(dist) %>%
    slice(1:2) %>%
    arrange(gw) %>%
    mutate(gw = round(gw, digits = 4)) -> dg_ann
  # Double check that the weights are equal in magnitude
  stopifnot(abs(dg_ann$gw[1]) == abs(dg_ann$gw[2]))
  YSTARTNEG <- 10
  YENDNEG <- 4.3
  XOFFSET <- 0
  ARROWLENGTH <- .1
  GROUPLABELNEG <- paste(dg_ann[1, 1:nj], collapse = ",")
  GROUPLABELPOS <- paste(dg_ann[2, 1:nj], collapse = ",")
  labelneg <- latex2exp::TeX(paste0("$G_{i} = (", GROUPLABELNEG, ")$"),
    output = "character"
  )
  labelpos <- latex2exp::TeX(paste0("$G_{i} = (", GROUPLABELPOS, ")$"),
    output = "character"
  )

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # Plot
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  XLABEL <- latex2exp::TeX(
    paste0("Distribution of weights across groups $\\omega_{AM}(G_{i})$")
  )
  p <- ggplot(
    data = dg %>% filter(sign != "Zero"),
    aes(x = gw, fill = sign)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_histogram(bins = 40) +
    labs(
      x = XLABEL,
      y = paste0("Number of groups (", nrow(dg), " total)"),
      fill = "Sign"
    ) +
    annotate("segment",
      x = XLOC_EXGR + XOFFSET, y = YSTARTNEG,
      xend = XLOC_EXGR, yend = YENDNEG,
      arrow = arrow(
        type = "closed",
        length = unit(ARROWLENGTH, "inches")
      )
    ) +
    annotate("label",
      x = XLOC_EXGR + XOFFSET, y = YSTARTNEG,
      label = labelneg, parse = TRUE
    ) +
    annotate("segment",
      x = -1 * XLOC_EXGR - XOFFSET, y = YSTARTNEG,
      xend = -1 * XLOC_EXGR, yend = YENDNEG,
      arrow = arrow(
        type = "closed",
        length = unit(ARROWLENGTH, "inches")
      )
    ) +
    annotate("label",
      x = -1 * XLOC_EXGR - XOFFSET, y = YSTARTNEG,
      label = labelpos, parse = TRUE
    ) +
    base_theme()

  if (safesavedir(savedir)) {
    fn <- fs::path(savedir, "average_monotonicity.pdf")
    ggsave(fn, width = 8, height = 4)
  }
  return(p)
}

average_monotonicity_df <- function() {
  stats::lm(
    data = stevenson,
    jail3 ~ 0 + judge_pre_1 + judge_pre_2 + judge_pre_3 + judge_pre_4 +
      judge_pre_5 + judge_pre_6 + judge_pre_7 + judge_pre_8
  ) %>%
    broom::tidy() %>%
    select(term, pscore = estimate) -> fs
  stevenson %>%
    select(starts_with("judge_pre")) %>%
    summarize_all(mean) %>%
    tidyr::pivot_longer(starts_with("judge")) %>%
    select(term = name, pz = value) %>%
    inner_join(fs, by = "term") -> df

  nj <- nrow(df)
  TOL <- 1e-15
  construct_groups(nj) %>%
    compute_group_weights(df$pscore, df$pz) %>%
    mutate(sign = if_else(gw < -1 * TOL, "Negative",
      if_else(gw > TOL, "Positive", "Zero")
    )) -> dg
  return(dg)
}


construct_groups <- function(nj = 3) {
  dg <- expand.grid(rep(list(0:1), nj))
  names(dg) <- paste0("D(", 1:nj, ")")
  dg$num1 <- rowSums(dg)
  return(dg)
}

compute_group_weights <- function(dg, pj, qj, TOL = 1e-6) {
  nj <- log(nrow(dg), base = 2)
  stopifnot(length(pj) == length(qj))
  stopifnot(all(pj >= 0))
  stopifnot(all(pj <= 1))
  stopifnot(nj == length(qj))
  stopifnot(all(qj >= 0))
  stopifnot(all(qj <= 1))
  stopifnot(abs(sum(qj) - 1) < TOL)

  dg$dbar <- NA
  dg$gw <- 0
  epj <- sum(qj * pj)
  for (i in 1:nrow(dg)) {
    dg$dbar[i] <- sum(dg[i, 1:nj] * qj)
    for (j in 1:nj) {
      dg$gw[i] <- dg$gw[i] + qj[j] * (pj[j] - epj) * (dg[i, j] - dg$dbar[i])
    }
  }
  if ("prg" %in% names(dg)) {
    dg$tslsw <- dg$gw * dg$prg / sum(dg$gw * dg$prg)
  }

  return(dg)
}

globalVariables(c(
  "gw", "stevenson", "name", "value", "gw"
))
