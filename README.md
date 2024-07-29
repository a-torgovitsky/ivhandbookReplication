# ivhandbookReplication

This repo contains replication code for the empirical results in Mogstad and Torgovitsky ("Instrumental Variables with Heterogeneous Treatment Effects," 2024, _Handbook of Labor Economics_).

If you are just looking for sample code in R and Stata, you probably want the [ivhandbook](https://github.com/a-torgovitsky/ivhandbook) repository instead.

## Installation

## Replicating

The following code replicates what is reported in the paper.
There may be minor differences in the Card application due to randomness in the bootstrap for the propensity score weighting estimates and randomness in the sample splits for the double/debiased machine learning (DDML) estimates. 

```r
runall("path/to/where/you/want/to/save/the/results")
```

The individual components that are run are:

```r
  multi_weighting() # Figure XXX
  ae_sensitivity() # Figure XXX
  average_monotonicity() # Figure XXX
  unordered_treatments() # Table XXX
  card_app(nbs = 500, nsplits = 100) # Table XXX
  gelbach_bounds(savedir) # Table XXX
  gelbach_comparisons(savedir) # A few comparison estimates discussed in the text
```

All of these run quickly except for `card_app`, which takes a while (about twelve hours on my machine) because of the double/debiased machine learning estimates.
Reducing `nsplits` will speed this up.
The `nbs` (number of bootstraps) parameter refers to the bootstrapped standard errors for the propensity score weighting estimate, which are pretty quick (a couple of minutes).

## Questions or comments?

Please post an issue.
