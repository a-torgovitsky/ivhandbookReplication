mdl_nnet <- function(y, X, ...) {
  nn <- nnet::nnet(x = X, y = y, ...)
  mdl <- list(model = nn)
  class(mdl) <- c("mdl_nnet", class(mdl))
  return(mdl)
}

predict.mdl_nnet <- function(object, newdata) {
  if (!inherits(object, "mdl_nnet")) {
    stop("The object is not of class 'mdl_nnet'")
  }
  yhat <- stats::predict(object$model, newdata = newdata)
  return(yhat)
}

define_learners <- function(testing = FALSE) {
  if (testing) { # small ensemble for quick testing
    list(
      list(
        fun = ddml::mdl_ranger,
        args = list(num.trees = 100)
      ),
      list(
        fun = ddml::mdl_ranger,
        args = list(num.trees = 100)
      )
    )
  } else {
    list(
      list(
        fun = ddml::mdl_xgboost,
        args = list(nrounds = 1000, eta = .01, max_depth = 2)
      ),
      list(
        fun = ddml::mdl_xgboost,
        args = list(nrounds = 1000, eta = .01)
      ),
      list(
        fun = ddml::mdl_xgboost,
        args = list(nrounds = 1000, eta = .05, max_depth = 3)
      ),
      list(
        fun = ddml::mdl_xgboost,
        args = list(nrounds = 1000, eta = .05, max_depth = 4)
      ),
      list(
        fun = ddml::mdl_ranger,
        args = list(num.trees = 500 -> NTREES)
      ),
      list(
        fun = ddml::mdl_ranger,
        args = list(num.trees = NTREES, max.depth = 4)
      ),
      list(
        fun = ddml::mdl_ranger,
        args = list(num.trees = NTREES, max.depth = 6)
      ),
      list(fun = mdl_nnet, args = list(
        size = 2, maxit = 1000, decay = 0.02,
        MaxNWts = 10000, linout = T, trace = F
      )),
      list(fun = mdl_nnet, args = list(
        size = 4, maxit = 1000, decay = 0.02,
        MaxNWts = 10000, linout = T, trace = F
      )),
      list(fun = mdl_nnet, args = list(
        size = 6, maxit = 1000, decay = 0.05,
        MaxNWts = 10000, linout = T, trace = F
      ))
    )
  }
}

run_ddml <- function(df, xlist, type = "pliv", nsplits = 1, testing = FALSE) {
  Y <- as.matrix(df$Y)
  D <- as.matrix(df$D)
  Z <- as.matrix(df$Z)
  X <- as.matrix(df[, xlist])

  learners <- define_learners(testing = testing)
  args <- list(
    y = Y, D = D, Z = Z, X = X,
    learners = learners,
    ensemble_type = "nnls1",
    sample_folds = 5,
    shortstack = TRUE,
    silent = TRUE,
    custom_ensemble_weights = diag(length(learners))
  )

  printname <- if_else(type == "pliv", "PLIV", "LATE/ACR")

  clear_line <- function() {
    cat(sprintf("\r"), rep("", 80))
  }

  # Returns a list of estimates, one for each rep
  lddml <- lapply(seq_len(nsplits), function(i) {
    if (type == "pliv") {
      a <- do.call(ddml::ddml_pliv, args)
    } else if (type == "late") {
      a <- do.call(ddml::ddml_late, args)
    } else {
      stop("Type is not one of `pliv` or `late`.")
    }
    clear_line()
    cat(sprintf("\rStarting %s iteration %d / %d", printname, i, nsplits))
    utils::flush.console()
    return(a)
  })

  clear_line()
  cat(sprintf("\rFinished %d %s iterations\n", nsplits, printname))

  r <- process_ddml(lddml)
  diag <- diagnose_ddml(lddml)
  return(list(
    results = r, diagnostics = diag, args = args,
    fulllist = lddml
  ))
}

process_ddml <- function(ddmllist) {
  lapply(ddmllist, function(r) {
    s <- summary(r)
    if (class(r) == "ddml_pliv") {
      df <- tibble(
        est = "PLIV (DDML)",
        coef = s["D_r", "Estimate", "nnls1"],
        se = s["D_r", "Std. Error", "nnls1"]
      )
    } else {
      df <- tibble(
        est = "ACR (DDML)",
        coef = s[1, "Estimate"],
        se = s[1, "Std. Error"]
      )
    }
  }) %>%
    bind_rows() %>%
    group_by(est) %>%
    summarize(
      sdest = sd(coef),
      minest = min(coef),
      maxest = max(coef),
      medcoef = stats::median(coef),
      medse = sqrt(stats::median(se^2 + (coef - stats::median(coef))^2))
    ) %>%
    rename(coef = medcoef, se = medse) %>%
    relocate(coef, se, minest, sdest, maxest)
}

diagnose_ddml <- function(ddmllist) {
  a <- ddmllist

  # The ddml package structures its results differently for pliv and late...
  classes <- sapply(a, class)
  if (all(classes == "ddml_pliv")) {
    type <- classes[1]
  } else if (all(classes == "ddml_late")) {
    type <- classes[1]
  } else {
    stop("ddmllist contains multiple types of classes")
  }

  summarize2 <- function(df, column_name) {
    column_quosure <- ensym(column_name)
    df %>%
      summarize(across(
        {{ column_quosure }},
        list(
          mean = mean,
          median = stats::median,
          sd = sd,
          min = min,
          max = max
        ),
        .names = "{fn}"
      ), .groups = "drop")
  }

  # Pull NNLS1 weights from each iteration and arrange
  lapply(seq_along(a), function(i) {
    rw <- a[[i]]$weights
    lapply(seq_len(length(rw)), function(ii) {
      tibble(fun = names(rw)[ii], weights = rw[[ii]][, "nnls1"]) %>%
        mutate(
          learner = row_number(),
          coef = a[[i]]$late[-1]
        )
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    group_by(learner, fun) %>%
    summarize2(weights) %>%
    arrange(fun, desc(mean)) %>%
    group_by(fun) %>%
    group_split() -> weights

  # Pull MSPE from each iteration and arrange
  lapply(a, function(x) {
    purrr::map_df(x$mspe,
      ~ as.data.frame(t(.)),
      .id = "model"
    )
  }) %>%
    bind_rows() %>%
    tidyr::pivot_longer(
      c(everything(), -model),
      names_to = "estimator",
      values_to = "mspe"
    ) %>%
    group_by(model, estimator) %>%
    summarize2(mspe) %>%
    arrange(model, desc(mean)) %>%
    group_by(model) %>%
    {
      setNames(group_split(.), group_keys(.) %>% pull(model))
    } -> mspe

  # Pull coefficient estimates from each iteration and arrange
  if (type == "ddml_late") {
    m <- sapply(a, function(x) x$late)
  } else {
    m <- sapply(a, function(x) x$coef)
  }
  m %>% t() -> m
  # tidyverse will complain with as_tibble if no column names
  colnames(m) <- colnames(m, do.NULL = FALSE, prefix = "col")
  m %>%
    as_tibble() %>%
    tidyr::pivot_longer(
      everything(),
      names_to = "estimator",
      values_to = "coef"
    ) %>%
    group_by(estimator) %>%
    summarize2(coef) %>%
    arrange(desc(estimator)) -> coefs

  return(list(coefs = coefs, mspe = mspe, weights = weights))
}

globalVariables(c(
  "est", "se", "medcoef", "medse", "minest", "sdest", "maxest",
  "learner", "fun", "model", "estimator"
))
