# GNAR-edge

Utilities for fitting, forecasting, and simulating a grouped GNAR model for
edge-valued time series.

## Usage

Source `R/gnar_edge_groups.R`, create an `igraph` whose edge order matches the
rows of the time-series matrix, and provide a named list assigning edge indices
to groups:

```r
source("R/gnar_edge_groups.R")

groups <- list(route_a = c(1, 3), route_b = c(2, 4))
fit <- gnar_edge_fit(
  ts_data, data_edges, alphaOrder = 2, betaOrder = c(1, 1),
  net, groups = groups
)

forecast <- gnar_edge_predict(
  fit, ts_data, alphaOrder = 2, betaOrder = c(1, 1),
  nedges = nrow(ts_data), npred = 5, wei_mat = fit$wei_mat,
  groups = groups
)
```

The package dependencies are `igraph` and `Matrix`. The functions validate
that groups form a non-overlapping partition of the modelled edges.
