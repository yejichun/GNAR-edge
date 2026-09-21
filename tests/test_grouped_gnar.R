source("R/gnar_edge_groups.R")

edges <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "d"))
network <- igraph::graph_from_data_frame(edges, directed = FALSE)
groups <- list(first = c(1L, 3L), second = 2L)
series <- matrix(seq_len(24), nrow = 3L)

stopifnot(identical(response_vec(series, 2L), as.vector(t(series[, 3:8]))))

design <- design_mat(series, edges, 2L, c(1L, 0L), network, groups)
stopifnot(
  identical(dim(design$dmat), c(18L, 6L)),
  identical(colnames(design$dmat), c(
    "alpha_1_first", "alpha_1_second", "alpha_2_first", "alpha_2_second",
    "beta_1_stage_1_first", "beta_1_stage_1_second"
  )),
  length(design$wei_mat) == 1L
)

fit <- gnar_edge_fit(series, edges, 2L, c(1L, 0L), network, groups = groups)
prediction <- gnar_edge_predict(
  fit, series, 2L, c(1L, 0L), nrow(series), 3L, fit$wei_mat,
  allcoefs = TRUE, groups = groups
)
stopifnot(
  identical(dim(prediction), c(3L, 5L)),
  isTRUE(all.equal(prediction[, 1:2], series[, 7:8], check.attributes = FALSE)),
  all(is.finite(prediction))
)

set.seed(42)
simulation <- gnar_edge_sim(
  n = 10L, net = network,
  alphaParams = list(rep(0.2, 3L)),
  betaParams = list(first = list(0.1), second = list(0.1)),
  sigma = 0.1, nedges = 3L, data_edges = edges, groups = groups
)
stopifnot(identical(dim(simulation), c(3L, 10L)), all(is.finite(simulation)))
