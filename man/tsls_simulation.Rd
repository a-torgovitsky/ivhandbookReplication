% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_weighting.R
\name{tsls_simulation}
\alias{tsls_simulation}
\title{Confirm TSLS weights in simulation}
\usage{
tsls_simulation(psc, prz, late, zeta, n = 10, numreps = 2)
}
\arguments{
\item{psc}{Propensity score vector}

\item{prz}{Marginal distribution of instrument}

\item{late}{Vector of LATEs}

\item{zeta}{Function of instrument used in TSLS}

\item{n}{Sample size for simulation}

\item{numreps}{Number of replications for the simulation}
}
\description{
Just check that I didn't screw up some algebra....
}
