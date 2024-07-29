#' @title Illustration of weights for IV estimands with multivalued instruments
#' @description Computing the IV estimand for different hypothetical experiments
#' Then create a plot to visualize the difference.
#'
#' @inheritParams shared_savedir
#' @export
multi_weighting <- function(savedir = getwd()) {
  psc <- c(.05, .2, .35)
  late <- c(1, -1)
  comblate <- (psc[2] - psc[1]) * late[1] + (psc[3] - psc[2]) * late[2]

  expand.grid(
    prtreated = c(.25, .5, .75),
    prhigh = seq(from = 0, to = 1, by = .01),
    zeta = c(ZETA_ID, ZETA_PSCORE)
  ) %>%
    mutate(prz = purrr::pmap(
      list(prtreated, prhigh),
      function(prtreated, prhigh, ...) generate_prz(prtreated, prhigh)
    )) %>%
    mutate(weights = purrr::pmap(
      list(prz, zeta),
      function(prz, zeta, ...) compute_tsls_weights(psc, prz, zeta)
    )) %>%
    mutate(
      tsls = purrr::map_dbl(weights, ~ sum(. * late)),
      percinc = 100 * prtreated
    ) -> df

  anndf <- data.frame(
    percinc = rep(min(df$percinc), 3),
    y = c(late, comblate - .15),
    label = c(
      "Average effect for $0 to $10 compliers",
      "Average effect for $10 to $50 compliers",
      "Average effect for all compliers"
    ),
    zeta = rep(NA, 3),
    x = rep(.50, 3)
  )

  p <- ggplot(data = df, aes(x = prhigh, y = tsls, color = zeta)) +
    geom_line() +
    geom_hline(yintercept = c(late, comblate), linetype = "dotted") +
    geom_label(
      data = anndf,
      aes(x = x, y = y, label = label),
      color = "black", fill = "white", size = 2.3
    ) +
    facet_wrap(
      ~percinc,
      labeller = labeller(percinc = function(v) paste0(v, "% incentivized"))
    ) +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(breaks = seq(-1, 1, by = .25)) +
    labs(
      color = "First stage specification",
      x = "Percentage of incentivized assigned to high incentives",
      y = "Value of the IV estimand"
    ) +
    base_theme()

  if (safesavedir(savedir)) {
    fn <- fs::path(savedir, "multiiv-weighting.pdf")
    ggsave(fn, width = 8, height = 4)
  }
  return(p)
}

generate_prz <- function(prtreated, prhigh) {
  prz0 <- 1 - prtreated
  prz1 <- (1 - prhigh) * prtreated
  prz2 <- prhigh * prtreated
  return(c(prz0, prz1, prz2))
}

complier_pr <- function(psc) {
  prcp <- diff(psc) # Pr[D(j+1) = 1, D(j) = 0] (j complier)
  stopifnot(all(prcp >= 0))
  return(prcp)
}

zetaf <- function(zeta, psc) {
  if (zeta == ZETA_ID) {
    function(z) 10 * (z == 1) + 50 * (z == 2)
  } else if (zeta == ZETA_PSCORE) {
    function(z) psc[z + 1]
  } else {
    stop("zeta entry not recognized")
  }
}

compute_tsls_weights <- function(psc, prz, zeta) {
  error_check(psc, prz)

  k <- length(psc) - 1
  zeta <- zetaf(zeta, psc)
  zetaz <- sapply(0:k, zeta)
  ezeta <- sum(zetaz * prz) # E[\zeta(Z)]

  # Cov[D, \zeta(Z)] = E[p(Z)(\zeta(Z) - \zeta(Z))]
  covdzeta <- sum(prz * psc * (zetaz - ezeta))
  prcp <- complier_pr(psc)
  covindzeta <- rep(NA, k - 1)
  w <- rep(NA, k - 1)
  for (j in 1:k) {
    ind <- rep(0, k + 1)
    ind[(j + 1):(k + 1)] <- 1
    covindzeta[j] <- sum(prz * ind * (zetaz - ezeta))
    w[j] <- (prcp[j] * covindzeta[j]) / covdzeta
  }
  names(w) <- paste0(paste0("weight", 0:(k - 1)), 1:k)
  return(w)
}

error_check <- function(psc, prz) {
  k <- length(psc) - 1
  stopifnot(length(prz) == k + 1)
  stopifnot(all(psc >= 0))
  stopifnot(all(psc <= 1))
  stopifnot(all(prz >= 0))
  stopifnot(all(prz <= 1))
  stopifnot(all.equal(sum(prz), 1))
}

#' @title Confirm TSLS weights in simulation
#' @description Just check that I didn't screw up some algebra....
#'
#' @param psc Propensity score vector
#' @param prz Marginal distribution of instrument
#' @param late Vector of LATEs
#' @param zeta Function of instrument used in TSLS
#' @param n Sample size for simulation
#' @param numreps Number of replications for the simulation
tsls_simulation <- function(psc, prz, late, zeta, n = 10, numreps = 2) {
  error_check(psc, prz)
  zeta <- zetaf(zeta, psc)
  k <- length(psc) - 1
  prcp <- complier_pr(psc)
  prat <- psc[1]
  prnt <- 1 - psc[k + 1]
  prg <- c(prat, prcp, prnt)
  stopifnot(all(prg >= 0))
  stopifnot(sum(prg) == 1)
  y0 <- rep(0, length(prg)) # set y0 = 0 for all groups---shouldn't matter
  y1 <- c(0, late, 0) # set y1 = 0 for at/nt and = late for compliers
  tsls <- rep(NA, numreps)
  sapply(seq_len(numreps), function(i) {
    z <- sample(0:k, n, replace = TRUE, prob = prz)
    g <- sample(seq_along(prg) - 1, n, replace = TRUE, prob = prg)
    d <- as.numeric(z >= g)
    y0 <- y0[g + 1]
    y1 <- y1[g + 1]
    y <- y0 + d * (y1 - y0)
    zeta <- zeta(z)
    tsls <- cov(y, zeta) / cov(d, zeta)
    return(tsls)
  }) -> r
  return(list(mean = mean(r), sd = sd(r)))
}

globalVariables(c(
  "prtreated", "prhigh", "prz", "zeta", "tsls", "x", "y", "label"
))
