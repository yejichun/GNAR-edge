# Grouped GNAR model for edge-valued time series.

.validate_gnar_inputs <- function(ts_data, alphaOrder, betaOrder, groups = NULL) {
  if (!is.matrix(ts_data) && !is.data.frame(ts_data))
    stop("ts_data must be a matrix or data frame")
  if (length(alphaOrder) != 1L || alphaOrder < 1L || alphaOrder >= ncol(ts_data))
    stop("alphaOrder must be positive and smaller than the number of observations")
  if (length(betaOrder) != alphaOrder || anyNA(betaOrder) || any(betaOrder < 0))
    stop("betaOrder must contain one non-negative value per alpha lag")
  if (!is.null(groups)) {
    if (is.null(names(groups)) || any(names(groups) == "") || anyDuplicated(names(groups)))
      stop("groups must be a uniquely named list")
    members <- unlist(groups, use.names = FALSE)
    if (!identical(sort(as.integer(members)), seq_len(nrow(ts_data))))
      stop("groups must partition the edge indices exactly once")
  }
}

response_vec <- function(ts_data, alphaOrder) {
  .validate_gnar_inputs(ts_data, alphaOrder, rep(0, alphaOrder))
  as.vector(t(as.matrix(ts_data)[, seq.int(alphaOrder + 1L, ncol(ts_data)), drop = FALSE]))
}

edge_neighbors <- function(net, maxstage, stage_nodes, groups) {
  if (maxstage < 1L) return(list())
  if (!igraph::is_connected(net, mode = "weak")) {
    components <- igraph::decompose(net, mode = "weak")
    selected <- Filter(function(x) all(stage_nodes %in% igraph::V(x)$name), components)
    if (!length(selected)) return(NULL)
    net <- selected[[1L]]
  }

  result <- vector("list", maxstage)
  excluded <- igraph::get_edge_ids(net, vp = stage_nodes, directed = FALSE)
  frontier <- as.character(stage_nodes)
  for (stage in seq_len(maxstage)) {
    incident_edges <- unique(unlist(lapply(frontier, function(node) {
      as.integer(igraph::incident(net, node, mode = "all"))
    }), use.names = FALSE))
    current <- setdiff(incident_edges, excluded)
    grouped <- setNames(vector("list", length(groups)), names(groups))
    for (group in names(groups)) grouped[[group]] <- intersect(current, groups[[group]])
    result[[stage]] <- grouped
    if (!length(current)) {
      frontier <- character()
    } else {
      frontier <- setdiff(unique(as.vector(igraph::ends(net, current))),
                          unique(as.vector(igraph::ends(net, excluded))))
    }
    excluded <- unique(c(excluded, current))
  }
  result
}

.neighbor_weights <- function(neighbors, edge, lead_lag_mat, use_lead_lag) {
  if (!length(neighbors)) return(numeric())
  if (!use_lead_lag) return(rep(1 / length(neighbors), length(neighbors)))
  if (is.null(lead_lag_mat)) stop("lead_lag_mat is required when lead_lag_weights is TRUE")
  weights <- as.numeric(lead_lag_mat[neighbors, edge])
  weights[is.na(weights)] <- 0
  weights
}

design_mat <- function(ts_data, data_edges, alphaOrder, betaOrder, net, groups,
                       lead_lag_mat = NULL, globalalpha = TRUE,
                       lead_lag_weights = FALSE) {
  .validate_gnar_inputs(ts_data, alphaOrder, betaOrder, groups)
  if (!globalalpha) stop("edge-specific alpha parameters are not supported by the grouped model")
  ts_data <- as.matrix(ts_data)
  nedges <- nrow(ts_data)
  predt <- ncol(ts_data) - alphaOrder
  group_names <- names(groups)
  alpha_names <- unlist(lapply(seq_len(alphaOrder), function(lag) {
    paste("alpha", lag, group_names, sep = "_")
  }))
  beta_names <- unlist(lapply(seq_len(alphaOrder), function(lag) {
    if (!betaOrder[lag]) return(character())
    unlist(lapply(group_names, function(group) {
      paste("beta", lag, "stage", seq_len(betaOrder[lag]), group, sep = "_")
    }))
  }))
  dmat <- matrix(0, predt * nedges, length(c(alpha_names, beta_names)),
                 dimnames = list(NULL, c(alpha_names, beta_names)))

  for (lag in seq_len(alphaOrder)) {
    columns <- seq.int(alphaOrder + 1L - lag, ncol(ts_data) - lag)
    for (group in group_names) {
      for (edge in groups[[group]]) {
        rows <- seq.int(predt * (edge - 1L) + 1L, predt * edge)
        dmat[rows, paste("alpha", lag, group, sep = "_")] <- ts_data[edge, columns]
      }
    }
  }

  if (!sum(betaOrder)) return(list(dmat = dmat, wei_mat = NULL))
  zero_group_mats <- setNames(lapply(groups, function(x) Matrix::Matrix(0, nedges, nedges,
                                                                       sparse = TRUE)), group_names)
  wei_mat <- lapply(seq_len(max(betaOrder)), function(x) zero_group_mats)
  for (edge in seq_len(nedges)) {
    nodes <- as.character(unlist(data_edges[edge, 1:2], use.names = FALSE))
    neighbors <- edge_neighbors(net, max(betaOrder), nodes, groups)
    if (is.null(neighbors)) next
    rows <- seq.int(predt * (edge - 1L) + 1L, predt * edge)
    for (lag in seq_len(alphaOrder)) {
      columns <- seq.int(alphaOrder + 1L - lag, ncol(ts_data) - lag)
      for (stage in seq_len(betaOrder[lag])) {
        for (group in group_names) {
          adjacent <- neighbors[[stage]][[group]]
          if (!length(adjacent)) next
          weights <- .neighbor_weights(adjacent, edge, lead_lag_mat, lead_lag_weights)
          values <- ts_data[adjacent, columns, drop = FALSE]
          weighted <- vapply(seq_len(ncol(values)), function(i) {
            present <- !is.na(values[, i])
            if (!any(present)) return(0)
            current_weights <- weights[present]
            if (!lead_lag_weights) current_weights <- current_weights / sum(current_weights)
            sum(values[present, i] * current_weights)
          }, numeric(1))
          name <- paste("beta", lag, "stage", stage, group, sep = "_")
          dmat[rows, name] <- weighted
          wei_mat[[stage]][[group]][adjacent, edge] <- weights
        }
      }
    }
  }
  list(dmat = dmat, wei_mat = wei_mat)
}

gnar_edge_fit <- function(ts_data, data_edges, alphaOrder, betaOrder, net,
                          lead_lag_mat = NULL, groups, globalalpha = TRUE,
                          lead_lag_weights = FALSE) {
  design <- design_mat(ts_data, data_edges, alphaOrder, betaOrder, net, groups,
                       lead_lag_mat, globalalpha, lead_lag_weights)
  response <- response_vec(ts_data, alphaOrder)
  keep <- !is.na(response) & stats::complete.cases(design$dmat)
  model_data <- data.frame(response = response[keep],
                           design$dmat[keep, , drop = FALSE], check.names = FALSE)
  model <- stats::lm(stats::reformulate(names(model_data)[-1L], response = "response",
                                        intercept = FALSE), data = model_data)
  list(mod = model, y = response, dd = design$dmat, wei_mat = design$wei_mat)
}

.model_coefficients <- function(fit_mod, allcoefs) {
  estimates <- stats::coef(fit_mod$mod)
  if (allcoefs) {
    estimates[is.na(estimates)] <- 0
    return(estimates)
  }
  table <- summary(fit_mod$mod)$coefficients
  output <- setNames(numeric(length(estimates)), names(estimates))
  available <- intersect(names(estimates), rownames(table))
  significant <- available[!is.na(table[available, 4]) & table[available, 4] < 0.05]
  output[significant] <- estimates[significant]
  output
}

gnar_edge_predict <- function(fit_mod, ts_data, alphaOrder, betaOrder, nedges,
                              npred, wei_mat, set.noise = NULL, allcoefs = FALSE,
                              globalalpha = TRUE, groups = NULL) {
  if (is.null(groups)) stop("groups must be supplied for a grouped GNAR forecast")
  .validate_gnar_inputs(ts_data, alphaOrder, betaOrder, groups)
  if (nedges != nrow(ts_data)) stop("nedges must equal nrow(ts_data)")
  if (npred < 1L) stop("npred must be positive")
  coefficients <- .model_coefficients(fit_mod, allcoefs)
  generated <- matrix(NA_real_, nedges, alphaOrder + npred)
  generated[, seq_len(alphaOrder)] <- as.matrix(ts_data)[,
    seq.int(ncol(ts_data) - alphaOrder + 1L, ncol(ts_data)), drop = FALSE]

  for (time in seq.int(alphaOrder + 1L, alphaOrder + npred)) {
    forecast <- numeric(nedges)
    for (lag in seq_len(alphaOrder)) {
      previous <- generated[, time - lag]
      for (group in names(groups)) {
        alpha_name <- paste("alpha", lag, group, sep = "_")
        forecast[groups[[group]]] <- forecast[groups[[group]]] +
          coefficients[alpha_name] * previous[groups[[group]]]
        for (stage in seq_len(betaOrder[lag])) {
          beta_name <- paste("beta", lag, "stage", stage, group, sep = "_")
          forecast <- forecast + coefficients[beta_name] *
            as.numeric(previous %*% as.matrix(wei_mat[[stage]][[group]]))
        }
      }
    }
    generated[, time] <- forecast
  }
  generated
}

nei_wei_mat <- function(net, data_edges, max.stage, nedges, groups) {
  if (max.stage < 1L) return(list())
  template <- setNames(lapply(groups, function(x) matrix(0, nedges, nedges)), names(groups))
  matrices <- lapply(seq_len(max.stage), function(x) template)
  for (edge in seq_len(nedges)) {
    nodes <- as.character(unlist(data_edges[edge, 1:2], use.names = FALSE))
    neighbors <- edge_neighbors(net, max.stage, nodes, groups)
    if (is.null(neighbors)) next
    for (stage in seq_len(max.stage)) for (group in names(groups)) {
      adjacent <- neighbors[[stage]][[group]]
      if (length(adjacent)) matrices[[stage]][[group]][adjacent, edge] <- 1 / length(adjacent)
    }
  }
  matrices
}

gnar_edge_sim <- function(n = 200, net, alphaParams, betaParams, sigma = 1,
                          meann = 0, nedges, data_edges, groups) {
  if (n < 1L || sigma < 0) stop("n must be positive and sigma must be non-negative")
  lags <- length(alphaParams)
  max_stage <- max(c(0L, unlist(lapply(betaParams, function(group) lengths(group)))))
  weights <- nei_wei_mat(net, data_edges, max_stage, nedges, groups)
  generated <- matrix(0, nedges, n + lags + 50L)
  generated[, seq_len(lags)] <- matrix(stats::rnorm(nedges * lags, sd = sigma), nedges)
  for (time in seq.int(lags + 1L, ncol(generated))) {
    forecast <- numeric(nedges)
    for (lag in seq_len(lags)) for (group in names(groups)) {
      previous <- generated[, time - lag]
      forecast[groups[[group]]] <- forecast[groups[[group]]] +
        alphaParams[[lag]][groups[[group]]] * previous[groups[[group]]]
      for (stage in seq_along(betaParams[[group]][[lag]])) {
        forecast <- forecast + betaParams[[group]][[lag]][stage] *
          as.numeric(previous %*% weights[[stage]][[group]])
      }
    }
    generated[, time] <- forecast + stats::rnorm(nedges, mean = meann, sd = sigma)
  }
  generated[, seq.int(ncol(generated) - n + 1L, ncol(generated)), drop = FALSE]
}
