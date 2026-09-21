# modified GNAR-edge model to include groups

library(igraph)
library(SparseM)
library(Matrix)
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### Modification of GNAR model (Knight et al.) for edge time series ###### 
###### ###### ###### ###### ###### ###### ###### ###### ######  ###### ######

response_vec<-function(ts_data,alphaOrder){
  # ts_data: dataframe with rows edges and columns time
  # alphaOrder: max lag
  # returns: vectorised response (concatenate ts for multiple time series)
  predt <- ncol(ts_data)-alphaOrder
  yvec <- NULL
  for(edge in 1:nrow(ts_data)){
    yvec <- c(yvec, ts_data[edge,((alphaOrder+1):(predt+alphaOrder))])
  }
  return(yvec)
}


design_mat<-function(ts_data,data_edges,alphaOrder,betaOrder,net,groups,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE){
  # ts_data: matrix with rows edges and columns time
  # data_edges: dataframe with edge list, column one node from which edge starts (character type) column two node edge ends (character type)
  # alphaOrder: max lag
  # betaOrder: vector of length as number of lags alphaOrder considered. Each entry give max stage neighbour to be considered for each lag
  # net: network resulting from edge list data_edges
  # lead_lag_mat: matrix with rows and col corresponding to edge time series, entries showing leadingness of each time series over the other
  # globalalpha: if True, non edge-specific alpha parameters, but only lag specific alpha param
  # lead_lag_weights: whether lead-lag matrix used as weights or not
  group_num <- sqrt(length(groups))
  nedges <- nrow(ts_data)
  predt <- ncol(ts_data)-alphaOrder # lenght of time series observations considered when removing max lag
  # CREATE NAMES FOR COLS OF DESIGN MAT ACCORDING TO model PARAMETERS
  parNames <- NULL
  for (ilag in 1:alphaOrder) {
    for (group in names(groups)) {
      parNames <- c(parNames,paste("alpha",ilag,group,sep="_"))
    # if(globalalpha){
    #   parNames <- c(parNames, paste("alpha", ilag, sep = ""))
    # }else{
    #   for(edg in 1:nedges){
    #     parNames <- c(parNames, paste("alpha", ilag, "edge", edg, sep=""))
    #   }
    # }
    }
  }
  for (ilag in 1:alphaOrder) {
    for (group in names(groups)) {
      if (betaOrder[ilag] > 0) {
        for (stag in 1:betaOrder[ilag]) {
          parNames <- c(parNames, paste("beta",ilag,"stage",stag,group,sep="_"))
        }
      }
    }
  }
  
  # CREATE EMPTY DESIGN MATRIX WITH 0s
  if(globalalpha){
    dmat <- matrix(0, nrow = predt * nedges, ncol = length(parNames), dimnames = list(NULL, parNames))
  }
  
  # FILL DESIGN MATRIX COLUMNS CORRESPONDING TO ALPHA PARAMETERS WITH LAGGED TIME SERIES FOR EACH EDGE
  # for (edg in 1:nedges) {
  for (group in names(groups)) {
    for (ilag in 1:alphaOrder) {
    # if(globalalpha){
    #   alphaLoc <- paste("alpha", ilag, sep = "")
    # }else{
    #   alphaLoc <- paste("alpha", ilag, "edge", edg, sep = "")
    # }
      alphaLoc <- paste("alpha",ilag,group,sep="_")
      for (edg in groups[[group]]) {
        dmat[((predt * (edg - 1) + 1):(predt * edg)), alphaLoc] <-  as.vector(ts_data[edg,((alphaOrder +1 - ilag):(predt + (alphaOrder - ilag)))])
      }
    }
  }
  # }
  
  if (sum(betaOrder) > 0) {
    wei_stag <- lapply(1:length(groups),function(x)Matrix(0,nedges,nedges,doDiag = FALSE))
    names(wei_stag) <- names(groups)
    wei_mat <- lapply(1:max(betaOrder),function(x)wei_stag)
    
    for (edg in 1:nedges){
      #print(c("edge",edg))
      
      # GET NEIGHBOUR EDGES 
      nodess <- c(data_edges[edg,1],data_edges[edg,2]) # get nodes involved in edg
      nei <- edge_neighbors(net,max(betaOrder),nodess,groups) # for these nodess get the r-stage neighbours
      
      if(!lead_lag_weights) {
        wei <- lapply(nei, function(x) lapply(x, length))
        wei <- lapply(wei, function(x) lapply(x, function(len) if(len>=1) rep(1/len,len) else 0))
      }

      # if (!lead_lag_weights){# if not using lead-lag matrix, equally weight neighbour edges (mean)
      #   wei <- sapply(nei, length)
      #   wei <- lapply(wei, function(x) rep(1/x,x))
      # }else{ # else use lead-lag matrix
      #   wei <- list()
      #   if ((!is.null(nei)) & (length(nei) > 0)) {
      #     for (ilagg in 1:alphaOrder){
      #       if (betaOrder[ilagg]>=1){
      #         for (stagg in 1:betaOrder[ilagg]){
      #           edg_loc <- edge_matching(data_edges, nei[[stagg]])
      #           wei[[stagg]] <- lead_lag_mat[edg_loc,edg]
      #         }
      #       }
      #     }
      #   }
      # }
      
      
      
      # IF THERE ARE NEIGHBOURS GET IN IF
      # if ((!is.null(nei)) & (length(nei) > 0)) { # check that are non empty lists
      if (!is.null(nei) && length(nei) > 0 && any(lengths(nei) > 0)) {  
        
        
        for (ilag in 1:alphaOrder){
          if (betaOrder[ilag]>=1){
            
            for (stag in 1:betaOrder[ilag]){
              for (g in names(groups)) {
                #print(c("ilag ",ilag, " stage ",stag) )
                betaLoc <- paste("beta", ilag, "stage", stag, g,sep = "_")
                
                if ((!is.null(nei[[stag]][[g]]))){
                  
                  # CASE WHERE stag-STAGE NEIGHBOUR NODES MORE THAN ONE (TO NORMALISE RESPECTIVELY)
                  if (length(nei[[stag]][[g]]) > 1) {
                    # get indexing/locations of neighbouring edges of edg to cut df accordingly
                    #edg_loc <- edge_matching(data_edges, nei[[stag]][[g]])
                    
                    
                    # CUT TIME SERIES MAT FOR COLUMNS CORRESPONDING TO NEIGHBOUR NODES AND ROWS CORRESP TO LAGGED TIME SERIES
                    vts.cut <- ts_data[nei[[stag]][[g]],((alphaOrder + 1 - ilag):(predt + (alphaOrder - ilag)))]
                    
                    # FOR CASE WHERE VTS.CUT CONTAINS NA VALUES -> IMPUTATION
                    # for each col in cut of time series matrix transform the multiple ts to single ts by getting their weighted sum
                    for (t in 1:ncol(vts.cut)) {
                      if (any(is.na(vts.cut[ ,t]))) {
                        if (all(is.na(vts.cut[ ,t]))) {
                          #if there are no neighbours left at any time point, set to zero
                          vts.cut[,t] <- 0
                        }
                        else {
                          new.wei <- wei[[stag]][[g]][which(!is.na(vts.cut[,t]))]
                          if (!lead_lag_weights){ # do the re-normalisation only if not lead-lag weights
                            new.wei <- new.wei/sum(new.wei) # re-normalise after removing weights for NA
                          }
                          sub.val <- vts.cut[which(!is.na(vts.cut[,t])),t] %*% new.wei # create from values of multiple time series without NAs at time t, a single value at time t by the weighted sum
                          vts.cut[which(is.na(vts.cut[,t])),t] <- sub.val # fill in NA with weighted sum of time series at time t without NAs
                        }
                      }
                    }
                    # ASSIGN TO DESIGN MAT THE SUM OF WEIGHTED TIME SERIES FOR LAG
                    dmat[((predt * (edg - 1) + 1):(predt * edg)), betaLoc] <- as.vector(t(vts.cut) %*% wei[[stag]][[g]])
                    wei_mat[[stag]][[g]][nei[[stag]][[g]],edg] <- wei[[stag]][[g]]
                  }
                  
                  # CASE WHERE ONLY ONE stag-STAGE NEIGHBOUR EDGE  (TO NORMALISE RESPECTIVELY)
                  else{
                    # get indexing/locations of neighbouring edges of edg to cut df accordingly
                    #edg_loc <- edge_matching(data_edges, nei[[stag]][[g]])
                    
                    # CASE WHERE stag-STAGE NEIGHBOUR EDGES IS EXACTLY ONE EDGE (NO NEED TO NORMALISE) AND NON NA VALUES
                    if ((length(nei[[stag]][[g]]) == 1) & (!is.na(nei[[stag]][[g]]))) {
                      vts.cut <- as.vector(ts_data[nei[[stag]][[g]],((alphaOrder + 1 - ilag):(predt + (alphaOrder - ilag)))])
                      #and if this is missing at any time point, set to zero (cannot fill in with mean value of neighbours as previously as no other neighbours)
                      vts.cut[is.na(vts.cut)] <- 0
                      dmat[((predt * (edg - 1) + 1):(predt * edg)), betaLoc] <- as.vector(t(vts.cut) * wei[[stag]][[g]])
                    }
                    
                    # CASE WHERE stag-STAGE NEIGHBOUR NAN VALUES
                    else {
                      # set all to 0
                      dmat[((predt * (edg - 1) + 1):(predt * edg)), betaLoc] <- 0
                    }
                    wei_mat[[stag]][[g]][nei[[stag]][[g]],edg] <- wei[[stag]][[g]]
                  }
                } # end if stag, group specific neighbours existing
                
                # case where no neighbours of stag group size
                else{
                  dmat[((predt * (edg - 1) + 1):(predt * edg)), betaLoc] <- 0
                }
              }
            } # end for stage
          } # end if check at least one neighbour considered
        } # end for ilag
        
        
        
      } # end if neighbours existing
      
      
      else{
        for (ilag in 1:alphaOrder){
          for (stag in 1:betaOrder[ilag]){
            betaLoc <- paste("beta", ilag, "stage", stag, group, sep = "_")
            dmat[((predt * (edg - 1) + 1):(predt * edg)),betaLoc] <- 0
          }
        }
      }
      
    } # END FOR LOOP OVER EDGES
  }else{ # END IF
    wei_mat <- NULL
  }
  out <- list(dmat=dmat,wei_mat=wei_mat)
  return(out)
}

gnar_edge_fit <- function(ts_data,data_edges,alphaOrder,betaOrder,net,lead_lag_mat,groups,globalalpha=TRUE,lead_lag_weights=FALSE){
  # ts_data: dataframe with rows edges and columns time
  # data_edges: dataframe with edge list, column one node from which edge starts column two node edge ends
  # alphaOrder: max lag
  # betaOrder: vector of length as number of lags alphaOrder considered. Each entry give max stage neighbour to be considered for each lag
  # net: network resulting from edge list data_edges
  # lead_lag_mat: matrix with rows and col corresponding to edge time series, entries showing leadingness of each time series over the other
  # globalalpha: if True, non edge-specific alpha parameters, but only lag specific alpha param
  # lead_lag_weights: whether lead-lag matrix used as weights or not
  # Returns: list with 1st element fitted model, 2nd element response vec, 3rd element design matrix, 4th element weight matrix used
  des <- design_mat(ts_data=ts_data, data_edges=data_edges,alphaOrder=alphaOrder,betaOrder=betaOrder,net=net,groups=groups,lead_lag_mat=lead_lag_mat,
                    globalalpha=globalalpha,lead_lag_weights = lead_lag_weights)
  dmat <- des$dmat
  yvec <- response_vec(ts_data = ts_data,alphaOrder = alphaOrder)

  # check if NA in response to remove, and respectively remove this rows from design mat
  if(sum(is.na(yvec))>0){
    yvec2 <- yvec[!is.na(yvec)]
    dmat2 <- dmat[!is.na(yvec),]
    modNoIntercept <- lm(yvec2~dmat2+0)
  }else{
    modNoIntercept <- lm(yvec~dmat+0)
  }
  out <- list(mod=modNoIntercept, y=yvec, dd=dmat, wei_mat=des$wei_mat)
  return(out)
}


edge_matching <- function(edge_set, neig_edge_set){
  # edge_set: two column set of edges for data set
  # neig_edge_set: igraph object, neighbouring edges of some specified edge
  # Returns: list with indices/location of edges in original data/data edges of ts
  ids <- as_ids(neig_edge_set)
  split_ids <- strsplit(ids,split='|', fixed=TRUE)
  ids_loc <- sapply(split_ids, function(x) which((edge_set[,1]==x[1]) & (edge_set[,2]==x[2])))
  return(ids_loc)
}


edge_neighbors <- function(net, maxstage, stage_nodes, groups){
  # net: igraph object 
  # maxstage: max stage of neighbour edges to be considered
  # stage_nodes: list of nodes corresponding to edge for which we seek neighbour edges
  # Returns: list of lists, with r element containing neighbour edges for stage r, grouped into the from-to groups
  if (!is_connected(net)){ # if graph is not connected, identify connected component involving stage_nodes to find r-stage neighb
    net_comp <- decompose(net)
    for (i in 1:length(net_comp)){ # find subgraph (component) involving stage_nodes of interest
      if ((stage_nodes[1] %in% vertex_attr(net_comp[[i]])$name) & (stage_nodes[2] %in% vertex_attr(net_comp[[i]])$name)){
        net <- net_comp[[i]]
        break
      }
    }
  }
  neighb_edges <- vector(mode = "list", length = maxstage)
  nei_nodes_acc <- get_edge_ids(net, vp=c(stage_nodes[1],stage_nodes[2]))
  #net <- net - edge(paste(stage_nodes[1],"|",stage_nodes[2],sep = ""))
  if (length(E(net))!=0){
    for (stag in 1:maxstage){
      neighb_edges_loc=list()
      for (node in stage_nodes){
        if (length(neighb_edges_loc)==0){
          neighb_edges_loc<-incident(net, node, mode = "all")#all, out, in
        }else{
          neighb_edges_loc<-igraph::union(neighb_edges_loc,incident(net, node, mode = "all"))#all, out, in
        }
      }
      neighb_edges_loc <- setdiff(as.vector(neighb_edges_loc),nei_nodes_acc)
      #print(neighb_edges_loc)
      #print(nei_nodes_acc)
      neighb_edges_loc_g <- vector(mode='list', length=length(names(groups)))
      names(neighb_edges_loc_g) <- names(groups)
      #print(as.numeric(neighb_edges_loc))
      for (loc in as.numeric(neighb_edges_loc)) {
        # idx <- Position(function(x) as.numeric(loc) %in% x, groups)
        # group <- names(groups)[idx]
        # neighb_edges_loc_g[[group]] <- c(neighb_edges_loc_g[[group]],as.numeric(loc))
        for (group in names(groups)) {
          if (loc %in% groups[[group]]) {
            neighb_edges_loc_g[[group]] <- c(neighb_edges_loc_g[[group]],loc) 
            break
          }
        }
      }
      neighb_edges[[stag]] <- neighb_edges_loc_g
      # create new list of nodes to get their incident edges for stages>1
      stage_nodes <- setdiff(as.vector(ends(net,neighb_edges_loc)),as.vector(ends(net,nei_nodes_acc)))
      nei_nodes_acc <- c(nei_nodes_acc, neighb_edges_loc)
      # remove previous neighbour edge from previous stage from graph, not be included in next set of edges
      #net <- net - neighb_edges_loc
      if (length(E(net))==0) break
    }
    return(neighb_edges)
  }else{
    neighb_edges <- NULL
    return(neighb_edges)
  }
}


gnar_edge_predict <- function(fit_mod, ts_data, alphaOrder, betaOrder, nedges, npred,
                              wei_mat, set.noise=NULL, allcoefs=FALSE, globalalpha=TRUE,
                              groups=NULL){
  # fit_mod    : fitted model from gnar_edge_fit
  # ts_data    : training data (nedges x T matrix)
  # alphaOrder : max lag
  # betaOrder  : vector of length alphaOrder; betaOrder[i] = max stage for lag i
  # nedges     : number of edges
  # npred      : number of steps to forecast
  # wei_mat    : when groups provided: list(stage -> list(group -> nedges x nedges matrix))
  #              otherwise:            list(stage -> nedges x nedges matrix)
  # groups     : named list of edge-index vectors (one entry per group); NULL for non-group model
  if (!is.null(set.noise)) sig <- set.noise else sig <- sigma(fit_mod$mod)

  if (!allcoefs) {
    nas  <- is.na(fit_mod$mod$coefficients)
    pvs  <- summary(fit_mod$mod)$coefficients[, 4] < 0.05
    vals <- rep(0, length(pvs))
    vals[pvs] <- summary(fit_mod$mod)$coefficients[pvs, 1]
    coefvec <- rep(0, length(nas))
    coefvec[!nas] <- vals
  } else {
    coefvec <- fit_mod$mod$coefficients
    coefvec[is.na(coefvec)] <- 0
  }

  xx.init <- ts_data[, (ncol(ts_data) - alphaOrder + 1):ncol(ts_data)]
  xx.gen  <- matrix(NA, ncol = npred + alphaOrder, nrow = nedges)
  xx.gen[, 1:alphaOrder] <- xx.init

  if (!is.null(groups)) {
    # ---- GROUP-STRUCTURED MODEL ----
    # Design matrix column layout (from design_mat in 1_fixed.R):
    #   All alphas first, ordered (ilag, group):
    #     alpha_1_g1, alpha_1_g2, ..., alpha_p_gG   [alphaOrder * n_groups entries]
    #   Then all betas, ordered (ilag, group, stag):
    #     beta_1_stage_1_g1, beta_1_stage_2_g1, ..., beta_1_stage_s1_g1,
    #     beta_1_stage_1_g2, ..., beta_1_stage_s1_g2, ... (all groups for lag 1)
    #     beta_2_stage_1_g1, ...,  (lag 2, etc.)
    n_groups    <- length(groups)
    group_names <- names(groups)

    # Extract alphas: alphaout[[ilag]][[g]] = scalar
    alphaout <- vector("list", alphaOrder)
    for (ilag in 1:alphaOrder) {
      alphaout[[ilag]] <- setNames(
        as.list(coefvec[(ilag - 1) * n_groups + seq_len(n_groups)]),
        group_names)
    }

    # Extract betas: betaout[[ilag]][[stag]][[g]] = scalar
    betaout     <- vector("list", alphaOrder)
    beta_offset <- alphaOrder * n_groups
    for (ilag in 1:alphaOrder) {
      betaout[[ilag]] <- vector("list", betaOrder[ilag])
      if (betaOrder[ilag] > 0) {
        for (gi in seq_along(group_names)) {
          g <- group_names[gi]
          for (stag in 1:betaOrder[ilag]) {
            if (is.null(betaout[[ilag]][[stag]])) betaout[[ilag]][[stag]] <- list()
            betaout[[ilag]][[stag]][[g]] <-
              coefvec[beta_offset + (gi - 1) * betaOrder[ilag] + stag]
          }
        }
        beta_offset <- beta_offset + n_groups * betaOrder[ilag]
      }
    }

    for (tpred in (alphaOrder + 1):(npred + alphaOrder)) {
      # Alpha contribution: each edge uses its group's scalar alpha
      time.forecast <- rep(0, nedges)
      for (ilag in 1:alphaOrder) {
        for (g in group_names) {
          idx <- groups[[g]]
          if (length(idx) > 0)
            time.forecast[idx] <- time.forecast[idx] +
              alphaout[[ilag]][[g]] * xx.gen[idx, tpred - ilag]
        }
      }

      # Neighbor contribution: sum over (lag, stage, group)
      nei.forecast <- rep(0, nedges)
      if (sum(betaOrder) > 0 && !is.null(wei_mat)) {
        for (ilag in 1:alphaOrder) {
          if (betaOrder[ilag] > 0) {
            for (stag in 1:betaOrder[ilag]) {
              for (g in group_names) {
                bval <- betaout[[ilag]][[stag]][[g]]
                wmat <- wei_mat[[stag]][[g]]
                if (!is.null(bval) && !is.null(wmat))
                  nei.forecast <- nei.forecast +
                    bval * as.numeric(xx.gen[, tpred - ilag] %*% as.matrix(wmat))
              }
            }
          }
        }
      }

      xx.gen[, tpred] <- time.forecast + nei.forecast
    }

  } else {
    # ---- ORIGINAL (NON-GROUP) MODEL ----
    if (globalalpha) {
      alphaout <- vector(mode = "list", length = alphaOrder)
      betaout  <- as.list(rep(0, length = alphaOrder))
      count <- 1
      for (ilag in 1:alphaOrder) {
        alphaout[[ilag]] <- coefvec[count]
        if (betaOrder[ilag] > 0)
          betaout[[ilag]] <- coefvec[(count + 1):(count + betaOrder[ilag])]
        count <- count + betaOrder[ilag] + 1
      }
    } else {
      alphaout <- vector(mode = "list", length = alphaOrder)
      betaout  <- as.list(rep(0, length = alphaOrder))
      count <- 1
      for (ilag in 1:alphaOrder) {
        alphaout[[ilag]] <- coefvec[count:(count + nedges - 1)]
        if (betaOrder[ilag] > 0)
          betaout[[ilag]] <- coefvec[(count + nedges):(count + nedges + betaOrder[ilag] - 1)]
        count <- count + nedges + betaOrder[ilag]
      }
    }

    for (tpred in (alphaOrder + 1):(npred + alphaOrder)) {
      time.forecast <- alphaout[[1]] * xx.gen[, tpred - 1]
      for (ilag in 2:alphaOrder)
        time.forecast <- time.forecast + alphaout[[ilag]] * xx.gen[, tpred - ilag]

      nei.forecast <- 0
      if (sum(betaOrder) > 0) {
        for (ilag in 1:alphaOrder) {
          bb <- length(betaout[[ilag]])
          if (bb > 0) {
            for (dd in 1:bb)
              nei.forecast <- nei.forecast +
                betaout[[ilag]][dd] * (xx.gen[, tpred - ilag] %*% wei_mat[[dd]])
          }
        }
      }
      xx.gen[, tpred] <- time.forecast + as.matrix(nei.forecast)
    }
  }

  return(xx.gen)
}

gnar_edge_sim <- function(n=200, net, alphaParams, betaParams, sigma=1, meann=0, nedges,data_edges, groups){
  # n: number of time stamps
  # net: network on whose edges are simulated
  # alphaParams: alpha parameters
  # betaParams: beta parameters
  # nedges: number of edges
  # data_edges: dataframe with edge list, column one node from which edge starts (character type) column two node edge ends (character type)
  
  max.nei <- max(unlist(lapply(betaParams, function(x) lapply(x, length))))
  
  #create weight matrices for neighbours
  nei.mats <- nei_wei_mat(net,data_edges,max.nei,nedges,groups)
  
  #seed the process from normal dist with mean 0 and sigma as given
  #do this for as many time points as needed for alpha
  lags <- length(alphaParams)
  #print(lags)
  # initialise time series to use lags according to max lag size, to produce next time steps
  xx.init <- matrix(rnorm(nedges*lags, mean=0, sd=sigma), nrow=nedges, ncol=lags)
  
  xx.gen <- matrix(0, ncol=n+50, nrow=nedges) # rows edges, columns time series 
  xx.gen[,1:lags] <- xx.init
  print(xx.init)
  
  for(tt in (lags+1):(n+50)){
    for(g in names(groups)){
      tmpi <- 0
      for(ilag in 1:lags){
        #  if(ilag==1){
        #    time.forecast <- alphaParams[[ilag]][groups[[g]]]*xx.gen[groups[[g]],(tt-ilag)]
        #  }else{
        #    tmpi <- alphaParams[[ilag]][groups[[g]]]*xx.gen[groups[[g]],(tt-ilag)]
        #    time.forecast <- time.forecast + tmpi
        #  }
        # }
        # lag_cols <- tt - (1:lags)
        # time.forecast <- xx.gen[groups[[g]], lag_cols,drop=FALSE] %*% alphaParams[[g]]
        
        tmpi <- tmpi + alphaParams[[ilag]][groups[[g]]]*xx.gen[groups[[g]],(tt-ilag)]
      }
      xx.gen[groups[[g]],tt] <- tmpi
    }
    if (tt ==n+50) {
      print(tt)
      print(xx.gen[,tt])
    }
    nei.forecast <- 0
    beta.pos <- NULL
    for(ilag in 1:lags){
      for(g in names(groups)) {
        bb <- length(betaParams[[g]][[ilag]])
        for(dd in seq_len(bb)){
          nei.forecast <- nei.forecast + betaParams[[g]][[ilag]][dd]*(xx.gen[,tt-ilag]%*%nei.mats[[dd]][[g]])
        }
      }
    }
    
    xx.gen[,tt] <- xx.gen[,tt]+nei.forecast+rnorm(nedges, mean=meann, sd=sigma)
  }
  
  return(xx.gen[,51:(n+50)])
}

nei_wei_mat <- function(net,data_edges,max.stage,nedges,groups){
  # net: igraph network
  # data.edges: edge list for net
  # max.stage: (scalar) max r-stage neighbour edges
  # nedges: (scalar) number of edges
  # groups: list with indices as groups(from and to) and values of a vector with edge locations
  # Rerurns: list of length of r, each list containing group # of r-stage matrices, each matrix of size nedges x nedges, 
  # with each column j corresponding to neighbours edges of edge j with equal weights
  wei_stag <- lapply(1:length(groups),function(x)matrix(0,nrow=nedges,ncol=nedges))
  names(wei_stag) <- names(groups)
  wei_mat <- lapply(1:max.stage,function(x) wei_stag)
  
  for (edg in 1:nedges){
    #print(c("edges",edg))
    # GET NEIGHBOUR EDGES 
    nodess <- c(data_edges[edg,1],data_edges[edg,2]) # get nodes involved in edg
    nei <- edge_neighbors(net,max.stage,nodess, groups) # for these nodess get the r-stage neighbours

    # equally weight neighbour edges (mean)
    wei <- lapply(nei, function(x) lapply(x, length))
    wei <- lapply(wei, function(x) lapply(x, function(len) if(len>=1) rep(1/len,len) else 0))

    for (stag in 1:max.stage){
      if (sum(lengths(nei[[stag]]))>0){
        for (g in names(groups)) {
          if (!is.null(nei[[stag]][[g]]) && length(nei[[stag]][[g]]) > 0) wei_mat[[stag]][[g]][nei[[stag]][[g]],edg] <- wei[[stag]][[g]] 
        }
      }else{
        #wei_mat[[stag]] <- NULL
        warning("beta order too large for network, neighbour set ",stag," is empty")
      }
    }
  }
  
  return(wei_mat)
}
