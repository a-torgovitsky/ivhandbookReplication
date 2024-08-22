# ivhandbookReplication

This repo contains replication code for the empirical results in Mogstad and Torgovitsky ("Instrumental Variables with Heterogeneous Treatment Effects," 2024, _Handbook of Labor Economics_).

If you are just looking for sample code in R and Stata, you probably want the [ivhandbook](https://github.com/a-torgovitsky/ivhandbook) repository instead.

## Installation

```r
devtools::install_github("a-torgovitsky/ivhandbookReplication")
```

Alternatively, if you don't want to install the package you can just clone it:

```shell
git clone git@github.com:a-torgovitsky/ivhandbookReplication.git
```
Then start R with the repo as the current working directory.
This should start `renv` bootstrapping after which you can:

```r
renv::restore() # Make package versions same as lockfile
install.packages("devtools") # Not included in renv
devtools::load_all() # Load the package locally
```

## Replicating

The following code replicates what is reported in the paper.
There may be minor differences in the Card application due to randomness in the bootstrap for the propensity score weighting estimates and randomness in the sample splits for the double/debiased machine learning (DDML) estimates. 

```r
library(ivhandbookReplication)
ivhandbook("path/to/where/you/want/to/save/the/results")
```

The individual components that are run inside `ivhandbook` are:

```r
  multi_weighting() # Figure 2
  ae_sensitivity() # Figure 3
  average_monotonicity() # Figure 4
  unordered_treatments() # Table 2
  card_app(nbs = 500, nsplits = 100) # Table 4
  gelbach_bounds() # Table 7
  gelbach_comparisons() # A few comparison estimates discussed in the text
```

All of these run quickly except for `card_app`, which takes a while (about twelve hours on my machine) because of the DDML estimates.
Reducing `nsplits` will speed this up.
The `nbs` parameter controls the number of bootstrap replications for computing standard errors in the propensity score weighting estimate.
These are pretty quick (a couple of minutes).

## Questions or comments?

Please post an issue.
