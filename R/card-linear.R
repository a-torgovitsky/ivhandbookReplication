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
    # Skip if specification is just a constant
    # RESET test would never reject here
    # And yhat2 yhat3 are collinear, which breaks the code below
    if (length(XList[[i]]) > 1) {
      rf <- make_ramsey(XList[[i]])
      rlm <- stats::lm(data = df, formula = rf)
      df %>%
        mutate(
          yhat = rlm$fitted.values,
          yhat2 = yhat^2,
          yhat3 = yhat^3
        ) -> df
      rf2 <- stats::update.formula(rf, . ~ . + yhat2 + yhat3)
      rlm2 <- stats::lm(data = df, formula = rf2)
      lh <- car::linearHypothesis(
        rlm2, c("yhat2", "yhat3"),
        white.adjust = "hc1",
        singular.ok = TRUE
      )
      pval <- lh$`Pr(>F)`[2]
    } else {
      pval <- 1
    }
    dfiv %>% mutate(ramsey = pval, cov = names(XList[i]))
  }) %>%
    bind_rows()
}

globalVariables(c(
  "term", "estimate", "std.error", "p.value", "estimate", "std.error"
))
