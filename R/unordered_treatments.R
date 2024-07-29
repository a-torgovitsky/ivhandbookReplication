#' @title Enumerate groups for unordered treatments, compare assumptions,
#' make a table
#'
#' @export
#' @inheritParams shared_savedir
unordered_treatments <- function(savedir = getwd()) {
  dlist <- seq(from = 0, to = 2, by = 1)
  G <- expand.grid(Z0 = dlist, Z1 = dlist, Z2 = dlist)

  cbind(G,
    MON = apply(G, 1, function(g) satisfies_mon(g, dlist)),
    EM = apply(G, 1, satisfies_em),
    IR = apply(G, 1, function(g) satisfies_ir(g, dlist)),
    NB = apply(G, 1, satisfies_nb)
  ) %>%
    as_tibble() %>%
    arrange(Z2, Z1, Z0, MON) %>%
    mutate(
      BCG = MON & EM,
      KLM = MON & IR & NB
    ) -> G

  G %>%
    filter(MON) %>%
    arrange(desc(EM), desc(IR), desc(NB), Z2, Z1, Z0) -> Gtab

  Gtab %>%
    mutate(across(
      !starts_with("Z"),
      function(x) if_else(x, "\\checkmark", "")
    )) %>%
    tibble::rownames_to_column("g") %>%
    mutate(group = paste0(
      "$(d_{", Z0, "}, d_{", Z1, "}, d_{", Z2,
      "}) \\equiv g_{", g, "}$"
    )) %>%
    select(-starts_with("Z"), -g, -BCG) %>%
    relocate(group, .before = 1) %>%
    split(seq(nrow(Gtab))) %>%
    lapply(unlist) %>%
    lapply(TexRow) %>%
    Reduce("+", .) -> tb

  header <- c("$G_{i}$", "Mon.", "EM", "IR", "NB", "KLM")

  TexRow(header) +
    TexMidrule() +
    tb -> tb

  vdots <- rep("", length(header))
  vdots[1] <- "$\\vdots$"
  gnotmono <- vdots
  gnotmono[1] <- "$(d_{1}, d_{0}, d_{2})$"
  tb +
    TexRow(vdots) +
    TexRow(gnotmono) +
    TexRow(vdots) -> tb

  save_table(tb, rep("c", length(header)), savedir, "unordered-groups")
  return(tb)
}

# Would you choose option j = 0,1,2 if you had instrument k = 0,1,2?
ch <- function(g, j, k) {
  return(as.integer(g[k + 1] == j))
}

# Monotonicity: Assumption 4 in KLM and 1(d) in HHKLM
satisfies_mon <- function(g, dlist) {
  for (j in dlist) { # Choices
    for (k in dlist) { # Instruments
      if (!(ch(g, j, j) >= ch(g, j, k))) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

# Condition in Proposition 2(ii) of KLM (``restrictive preferences'')
# Used in BCG (``extended monotonicity'') but not KLM
satisfies_em <- function(g) {
  if ((ch(g, 2, 0) == ch(g, 2, 1)) * (ch(g, 1, 2) == ch(g, 1, 0))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Condition in Proposition 2(iii) of KLM (``irrelevance'')
# or 2(a) in HHKLM
satisfies_ir <- function(g, dlist) {
  for (k in dlist) {
    if (ch(g, k, k) == ch(g, k, 0)) {
      for (j in dlist) {
        if (ch(g, j, k) != ch(g, j, 0)) {
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}

# Condition in Proposition 2(iv) of KLM (``satisfies next best'')
# --> I'm removing the always-takers here as well (KLM/HHKLM)
satisfies_nb <- function(g) {
  if ((ch(g, 1, 0) == 0) * (ch(g, 2, 0) == 0)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

globalVariables(c(
  "Z0", "Z1", "Z2", "MON", "EM", "IR", "NB", "g", "BCG", "group", "."
))
