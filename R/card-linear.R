run_ols <- function(df, XList) {
  lapply(seq_along(XList), function(i) {
    make_ols(XList[[i]]) %>%
      stats::lm(data = df, formula = .) %>%
      lmtest::coeftest(vcov = sandwich::vcovHC, type = "HC0") %>%
      broom::tidy() %>%
      filter(term == "D") %>%
      select(coef = estimate, se = std.error) %>%
      mutate(cov = names(XList[i]))
  }) %>%
    bind_rows()
}

run_lineariv <- function(df, XList) {
  lapply(seq_along(XList), function(i) {
    make_lineariv(XList[[i]]) %>%
      ivreg::ivreg(data = df, formula = .) %>%
      lmtest::coeftest(vcov = sandwich::vcovHC, type = "HC0") %>%
      broom::tidy() %>%
      filter(term == "D") %>%
      select(coef = estimate, se = std.error) -> dfiv
    make_ramsey(XList[[i]]) %>%
      stats::lm(data = df, formula = .) %>%
      lmtest::resettest() -> rt
    rt[["parameter"]] <- NULL # Otherwise tidy will issue a warning/complaint/notice
    rt %>%
      broom::tidy() %>%
      select(ramsey = p.value) %>%
      bind_cols(dfiv, .) %>%
      mutate(cov = names(XList[i]))
  }) %>%
    bind_rows()
}

globalVariables(c(
  "term", "estimate", "std.error", "p.value", "estimate", "std.error"
))
