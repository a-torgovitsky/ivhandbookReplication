#' @title Bounds for the Gelbach illustration
#'
#' @inheritParams shared_savedir
#' @export
gelbach_bounds <- function(savedir = getwd()) {
  YBDS <- c(.4, .8) # substantive bounds

  # Intersection bounds
  df <- gelbach_moments(gelbach)
  iv_bds_log <- manski_bounds(df)
  iv_bds_sub <- manski_bounds(df, ybds = YBDS)

  # Manski-Robins ("worst-case") bounds as a special case with constant Z
  gelbach %>%
    mutate(const = 1) %>%
    gelbach_moments(instrument = "const") -> df2
  wc_bds_log <- manski_bounds(df2)
  wc_bds_sub <- manski_bounds(df2, ybds = YBDS)

  join_bd_dfs <- function(df1, df2) {
    lapply(list(df1, df2), function(df) select(df, z, d, p, lb, ub)) %>%
      {
        left_join(.[[1]], .[[2]], by = c("z", "d", "p"))
      }
  }

  dfout <- bind_rows(
    join_bd_dfs(wc_bds_log, wc_bds_sub),
    join_bd_dfs(iv_bds_log, iv_bds_sub)
  )

  exp_str <- function(d) {
    paste0("$\\mathbb{E}[Y_{i}(", d, ")]$")
  }
  cexp_str <- function(d, z) {
    s <- exp_str(d)
    s <- stringr::str_sub(s, 1, -3)
    paste0(s, "\\vert Z_{i} = \\text{", z, "}]$")
  }

  format_df <- function(df) {
    df %>%
      mutate(z = dplyr::case_when(
        ((z == 1) | (is.na(z) & !is.na(d))) ~ exp_str(d),
        (is.na(z) & is.na(d)) ~ "$\\mathbb{E}[Y_{i}(1) - Y_{i}(0)]$",
        TRUE ~ cexp_str(d, z)
      )) %>%
      select(!d) %>%
      mutate(across(where(is.numeric), ~ format_numeric(.)))
  }

  join_bd_dfs(wc_bds_log, wc_bds_sub) %>%
    slice(-c(2, 4)) %>%
    format_df() -> wc
  join_bd_dfs(iv_bds_log, iv_bds_sub) %>%
    format_df() -> iv

  # In the table, worst-case will include the un-intersected IV bounds
  wc_tab <- bind_rows(wc, iv[c(1:4, 6:9), ])
  iv_tab <- iv[c(5, 10, 11), ]
  gelbach_table(wc_tab, iv_tab, YBDS, savedir)

  return(list(wc = wc, iv = iv))
}

gelbach_moments <- function(
    datain,
    instrument = "quarter",
    treatment = "public",
    outcome = "work79") {
  datain %>%
    rename(
      z = !!sym(instrument),
      d = !!sym(treatment),
      y = !!sym(outcome)
    ) -> datain
  datain %>%
    group_by(z) %>%
    summarize(
      q = n() / nrow(.),
      p = mean(d),
    ) -> df_pscore

  datain %>%
    group_by(z, d) %>%
    summarize(
      y = mean(y),
      .groups = "drop"
    ) %>%
    left_join(df_pscore, by = "z") %>%
    mutate(z = as.character(haven::as_factor(z))) -> df
  return(df)
}

manski_bounds <- function(df, ybds = c(0, 1)) {
  df %>%
    arrange(d, p) %>%
    mutate(
      lb = d * (y * p + ybds[1] * (1 - p)) +
        (1 - d) * (ybds[1] * p + y * (1 - p)),
      ub = d * (y * p + ybds[2] * (1 - p)) +
        (1 - d) * (ybds[2] * p + y * (1 - p))
    ) -> bds
  bds %>%
    group_by(d) %>%
    summarize(lb = max(lb), ub = min(ub)) %>%
    arrange(d) -> ibds

  bds %>%
    add_row(d = 0, lb = ibds$lb[1], ub = ibds$ub[1], .after = nrow(bds) / 2) %>%
    add_row(d = 1, lb = ibds$lb[2], ub = ibds$ub[2]) %>%
    add_row(lb = ibds$lb[2] - ibds$ub[1], ub = ibds$ub[2] - ibds$lb[1])
}

gelbach_table <- function(wc, iv, ybds, savedir) {
  boundstr <- function(lb, ub) {
    paste0(
      "($y_{\\text{lb}} = ", lb, ", ",
      "y_{\\text{ub}} = ", ub, "$)"
    )
  }

  TexRow(
    c("", "", "\\textbf{Logical}", "\\textbf{Substantive}"),
    cspan = c(1, 1, 2, 2),
    position = "c"
  ) +
    TexRow(
      c("", "", boundstr(0, 1), boundstr(ybds[1], ybds[2])),
      cspan = c(1, 1, 2, 2),
      position = "c"
    ) +
    TexMidrule(list(c(3, 4), c(5, 6))) +
    TexRow(c("", "$\\mathbf{p(z)}$", rep(c("\\textbf{LB}", "\\textbf{UB}"), 2))) +
    TexMidrule() +
    TexRow(
      c("", "\\textit{Manski bounds}"),
      cspan = c(1, 5),
      position = "c"
    ) +
    Reduce("+", wc %>% purrr::pmap(~ TexRow(c(...)))) +
    TexRow(
      c("", "\\textit{Instrumental variable bounds}"),
      cspan = c(1, 5),
      position = "c"
    ) +
    Reduce("+", iv %>% purrr::pmap(~ TexRow(c(...)))) -> tb
  save_table(tb, c("r", rep("c", 5)), savedir, "gelbach-results")
  return(tb)
}

#' @title Comparisons for the Gelbach illustration
#'
#' @inheritParams shared_savedir
#' @param outcome A string for the name of the Y variable
#' @export
gelbach_comparisons <- function(savedir = getwd(), outcome = "work79") {
  gelbach %>%
    mutate(
      age2 = age^2,
      quarter = as.factor(quarter)
    ) -> df

  if (safesavedir(savedir)) {
    fc <- file(fs::path(savedir, "gelbach_comparisons.out"), open = "wt")
    sink(fc)
  }


  postest <- function(fitmodel) {
    lmtest::coeftest(fitmodel, vcov = sandwich::vcovHC, type = "HC0") %>%
      broom::tidy() %>%
      print()
  }

  cat("Standard errors for propensity scores\n")
  stats::lm(
    data = df,
    public ~ 0 + as.factor(quarter)
  ) %>%
    postest()
  cat("\nStandard errors for outcome moments\n")
  stats::lm(
    data = df %>%
      mutate(ind = interaction(as.factor(public), as.factor(quarter))),
    formula(paste(outcome, "~ 0 + ind"))
  ) %>%
    postest()

  controls <- c(
    "num612",
    "num1317",
    "numge18",
    "othrlt18",
    "othrge18",
    "grade",
    "white",
    "centcity",
    "age",
    "age2",
    "as.factor(state)",
    "as.factor(bmodk5)"
  )

  postest <- function(fitmodel) {
    lmtest::coeftest(fitmodel, vcov = sandwich::vcovHC, type = "HC0") %>%
      broom::tidy() %>%
      filter(term == "public") %>%
      print()
  }

  cat("\nGelbach Table 7, Column 1\n")
  stats::lm(
    data = df,
    formula(paste(outcome, "~ public"))
  ) %>% postest()

  cat("\nGelbach Table 7, Column 2\n")
  stats::lm(
    data = df,
    formula(paste(outcome, "~ public +", paste(controls, collapse = "+")))
  ) %>%
    postest()

  cat("\nGelbach Table 7, Column 3\n")
  ivstr <- "| public | quarter"
  ivreg::ivreg(
    data = df,
    formula(paste(outcome, "~", paste(controls, collapse = "+"), ivstr))
  ) %>%
    postest()

  cat("\nGelbach Table 7, Column 3, but no controls\n")
  ivreg::ivreg(
    data = df,
    formula(paste(outcome, "~", substring(ivstr, 3)))
  ) %>%
    postest()

  cat("\nGelbach Table 7, Column 3, but no controls, binarized instrument \n")
  df %>%
    group_by(quarter) %>%
    summarize(mean(public)) # Check order
  df %>%
    mutate(qob_binary = if_else(quarter %in% c(2, 3), 1, 0)) -> df
  ivreg::ivreg(
    data = df,
    formula(paste(outcome, "~", "public | qob_binary"))
  ) %>%
    postest()

  cat("\nLinear MTE estimate of the ATE\n")
  args <- list(
    data = df,
    outcome = outcome,
    propensity = public ~ quarter,
    target = "ate",
    m0 = ~u,
    m1 = ~u,
    point = TRUE,
    bootstraps = 200
  )
  r <- do.call(ivmte::ivmte, args)
  se <- ifelse(exists("point.estimate.se", r), r$point.estimate.se, NA_real_)
  cat(paste(round(r$point.estimate, 3), "(", round(se, 3), ")\n"))

  sink(NULL)
  close(fc)
}

globalVariables(c(
  "gelbach", "z", "d", "p", "lb", "ub", "gelbach", "age", "quarter",
  "public", "y"
))
