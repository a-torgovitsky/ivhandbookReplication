qweightlate_est <- function(df, xlist) {
  stats::glm(
    data = df,
    formula = make_ramsey(xlist),
    family = stats::binomial
  ) %>%
    stats::predict(type = "response") %>%
    tibble(Q = .) %>%
    bind_cols(df, .) -> df
  df %>%
    mutate(
      w1 = Z / Q,
      num1 = Y * w1,
      w0 = (1 - Z) / (1 - Q),
      num0 = Y * w0,
      den1 = D * w1,
      den0 = D * w0
    ) %>%
    summarize(across(c(w0, w1, num0, num1, den0, den1), sum)) -> x
  (x$num1 / x$w1 - x$num0 / x$w0) / (x$den1 / x$w1 - x$den0 / x$w0)
}

qweightlate <- function(df, xlist, bsnum = 0) {
  n <- nrow(df)
  nl <- seq_len(n)
  lapply(seq(from = 0, to = bsnum, by = 1), function(i) {
    if (i == 0) {
      dfbs <- df
    } else {
      idx <- sample(nl, n, replace = TRUE)
      dfbs <- df[idx, ]
    }
    late <- qweightlate_est(dfbs, xlist)
    tibble(late = late, rep = i)
  }) %>%
    bind_rows() -> r
  list(
    coef = r %>% filter(rep == 0) %>% pull(late),
    se = r %>% filter(rep != 0) %>% summarize(se = sd(late)) %>% pull(se)
  )
}

run_qweightlate <- function(df, XList, bsnum = 0) {
  lapply(seq_along(XList), function(i) {
    qweightlate(df, XList[[i]], bsnum = bsnum) %>%
      as_tibble() %>%
      mutate(cov = names(XList[i]))
  }) %>%
    bind_rows()
}

globalVariables(c(
  "Z", "Q", "Y", "w1", "w0", "num0", "num1", "den0", "den1", "late", "se"
))
