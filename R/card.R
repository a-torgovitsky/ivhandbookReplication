#' @title Empirical results for the Card illustration
#'
#' @inheritParams shared_savedir
#' @param nbs Number of bootstrap replications for p-score weighting estimates
#' @param nsplits Number of split repetitions for DDML
#' @export
card_app <- function(savedir = getwd(), nbs = 2, nsplits = 1) {
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # Setup
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  df <- prepcard()
  XList <- create_covlist(df)

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # Quick estimators
  #   (propensity score weighting takes a little bit with bootstraps,
  #      but still very fast)
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  run_ols(df, XList) %>%
    mutate(est = "OLS") -> dfr
  run_lineariv(df, XList) %>%
    mutate(est = "Linear IV") %>%
    bind_rows(dfr) -> dfr
  run_qweightlate(df, XList[-1], bsnum = nbs) %>%
    mutate(est = "ACR (weighting)") %>%
    bind_rows(dfr) -> dfr

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # DDML (the time consuming part)
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  specs <- expand.grid(
    est = c("pliv", "late"),
    cov = c("XGeo", "XCard"),
    stringsAsFactors = FALSE
  )
  seq_specs <- seq_len(nrow(specs))
  r <- lapply(seq_specs, function(i) {
    run_ddml(
      type = specs$est[i],
      df,
      XList[[specs$cov[i]]],
      nsplits = nsplits,
      testing = FALSE
    ) %>%
      {
        c(., id = paste(specs$est[i], specs$cov[i]))
      }
  })
  lapply(seq_specs, function(i) {
    r[[i]]$results %>%
      # Coding these as being in the with interactions column
      mutate(
        cov = ifelse(specs$cov[i] == "XGeo", "XGeoInt", "XCardInt"),
        est = ifelse(specs$est[i] == "pliv", "PLIV (DDML)", "ACR (DDML)")
      )
  }) %>%
    bind_rows(dfr) -> dfr

  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # Finished; save results and create table
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  if (safesavedir(savedir)) {
    results <- dfr
    ddml <- r
    save(results, ddml, file = fs::path(savedir, "card.RData"))
    fc <- file(fs::path(savedir, "card.out"), open = "wt")
    sink(fc)
    print(results)
    print(paste(
      "Maximum Ramsey p-value is ",
      max(results %>% filter(!is.na(ramsey)) %>% pull(ramsey))
    ))
    sink(NULL)
    close(fc)
  }

  card_table(dfr, savedir)
  return(dfr)
}

make_linearx <- function(xlist) {
  paste(xlist, collapse = " + ")
}

make_ols <- function(xlist) {
  sc <- make_linearx(xlist)
  sc <- if (sc != "") paste(" + ", sc)
  s <- paste0("Y ~ D ", sc)
  stats::as.formula(s)
}

make_lineariv <- function(xlist) {
  sc <- make_linearx(xlist)
  sc <- if (sc != "") paste(sc, " | ")
  s <- paste0("Y ~ ", sc, " D | Z")
  stats::as.formula(s)
}

make_ramsey <- function(xlist) {
  sc <- make_linearx(xlist)
  sc <- if (sc == "") "1" else sc
  s <- paste0("Z ~ ", sc)
  stats::as.formula(s)
}

prepcard <- function() {
  card %>%
    mutate(
      Y = lwage,
      Z = nearc4,
      D = educ,
      south_smsa = south * smsa,
      south_smsa66 = south * smsa66,
      smsa_smsa66 = smsa * smsa66,
      black_exper = black * exper,
      black_expersq = black * expersq,
      black_south = black * south,
      black_smsa = black * smsa,
      black_smsa66 = black * smsa66,
      exper_south = exper * south,
      exper_smsa = exper * smsa,
      exper_smsa66 = exper * smsa66,
      expersq_south = expersq * south,
      expersq_smsa = expersq * smsa,
      expersq_smsa66 = expersq * smsa66
    )
}

create_covlist <- function(df) {
  varnames <- names(df)
  XNull <- c("")
  XGeo <- c(
    "south", "south66", "smsa", "smsa66",
    varnames[stringr::str_detect(varnames, "reg66")]
  )
  XGeoInt <- c(XGeo, "south_smsa", "south_smsa66", "smsa_smsa66")
  XDem <- c("black", "exper", "expersq")
  XCard <- c(XGeo, XDem)
  XCardInt <- c(
    XCard,
    XGeoInt,
    varnames[stringr::str_detect(varnames, "black_")],
    varnames[stringr::str_detect(varnames, "exper_")],
    varnames[stringr::str_detect(varnames, "expersq_")]
  )
  XList <- list(
    NoCov = XNull,
    XGeo = XGeo,
    XGeoInt = XGeoInt,
    XCard = XCard,
    XCardInt = XCardInt
  )
}

card_table <- function(df, savedir = getwd()) {
  df %>%
    mutate(
      estn = case_when( # Order of estimators
        est == "OLS" ~ 1,
        est == "Linear IV" ~ 2,
        stringr::str_detect(est, "PLIV") ~ 3,
        stringr::str_detect(est, "weighting") ~ 4,
        est == "ACR (DDML)" ~ 5,
        TRUE ~ NA_integer_
      ),
      coln = case_when( # Order of covariate specifications
        (cov == "NoCov") ~ 1,
        (cov == "XGeo") ~ 2,
        (cov == "XGeoInt") ~ 3,
        (cov == "XCard") ~ 4,
        (cov == "XCardInt") ~ 5,
        TRUE ~ NA_integer_
      )
    ) -> df

  tb <- list()
  NSPECS <- max(df$coln)
  collist <- tibble(coln = seq_len(NSPECS))
  DIGITS <- 3
  SPRINTF <- paste0("%.", DIGITS, "f")
  for (i in seq_len(max(df$estn))) {
    df %>%
      filter(estn == i) %>%
      select(est, coef, se, coln) %>%
      arrange(coln) -> dfest
    estname <- dfest$est[1]
    dfest <- left_join(collist, dfest, by = "coln")
    dfest %>%
      mutate(
        coef = if_else(
          is.na(coef),
          "---",
          sprintf(SPRINTF, round(coef, DIGITS))
        ),
        se = if_else(
          is.na(se),
          "",
          paste0("(", sprintf(SPRINTF, round(se, DIGITS)), ")")
        )
      ) ->
    dfest
    dfest %>%
      pull(coef) %>%
      c(estname, .) %>%
      TexRow(dec = DIGITS) -> tbest
    dfest %>%
      pull(se) %>%
      c("", .) %>%
      TexRow(dec = DIGITS) -> tbse
    tb[[i]] <- tbest + tbse
  }
  tb <- Reduce("+", tb)

  checkmark_row <- function(title, specs) {
    r <- rep("", NSPECS + 1)
    r[1] <- title
    r[specs + 1] <- "\\checkmark"
    return(r)
  }

  TexRow(c("", paste0("(", seq_len(NSPECS), ")"))) +
    TexMidrule() +
    tb +
    TexMidrule() +
    TexRow(checkmark_row("Geographic controls", c(2:5))) +
    TexRow(checkmark_row("Geographic interactions", c(3, 5))) +
    TexRow(checkmark_row("Demographic controls", c(4, 5))) +
    TexRow(checkmark_row("Demographic interactions", c(5))) -> tb

  save_table(tb, c("l", rep("c", NSPECS)), savedir, "card-results")
}

globalVariables(c(
  "card", "lwage", "nearc4", "educ", "south", "smsa", "black", "exper",
  "expersq", "ramsey", "estn", "est", "se", "coln", "smsa66"
))
