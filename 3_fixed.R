# modified GNAR-edge model to include groups

library(igraph)
library(vars)
library(forecast)
library(reshape2)
library(ggplot2)
library(dplyr)

########################################################
####### SIMULATION EXPERIMENTS FOR SECTION 5.2 #######
########################################################

####### SECTION 5.2.1 #######

# G = 2 groups
group_num <- 2

# Regime 4: GNAR(3,[2,2,2]) — group-structured parameters (matches 2_fixed.R Regime 4)
alpha_par_4 <- list(g1_g1=c(.2,.4,-0.6), g1_g2=c(.3,.1,-0.5), g2_g1=c(.4,.1,-0.5), g2_g2=c(.1,.3,-0.2))
beta_par_4 <- list(g1_g1=list(c(0.3,.1),c(0.1,.1),c(0,-0.2)),
                   g1_g2=list(c(-0.1,-0.1),c(-0.1,0.3),c(-0.3,-0.1)),
                   g2_g1=list(c(0.3,-.1),c(-0.1,-.1),c(0.2,.3)),
                   g2_g2=list(c(-0.3,0.1),c(-0.2,.3),c(0.2,-.3)))

########### PREDICTIVE PERFORMANCE: REPETITIONS OF SBM,ER,RDP GRAPHS WITH 4TH SIM REGIME ###########

set.seed(44)
seeds <- sample(1:500,50,replace=FALSE)
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 12,nrow = 50)
colnames(rmse_mat) <- c("gnarnei_er","gnar_er","var_er","ar_er",
                        "gnarnei_sbm","gnar_sbm","var_sbm","ar_sbm",
                        "gnarnei_grg","gnar_grg","var_grg","ar_grg")

for (i in 1:50){
  print(i)
  set.seed(seeds[i])

  # SBM model
  net_sbm <- sample_sbm(20,probmat,c(10,10),directed = TRUE)
  V(net_sbm)$name <- as.character(seq(1,20,1))
  V(net_sbm)$group <- sample(c('g1','g2'), size=20, replace=TRUE)
  edgelist_sbm <- get.edgelist(net_sbm)
  nedges_sbm <- ecount(net_sbm)
  group_from_sbm <- V(net_sbm)$group[as.numeric(edgelist_sbm[,1])]
  group_to_sbm   <- V(net_sbm)$group[as.numeric(edgelist_sbm[,2])]
  group_from_to_sbm <- paste(group_from_sbm, group_to_sbm, sep='_')
  grouplist_sbm <- split(seq_along(group_from_to_sbm), group_from_to_sbm)
  lst_sbm <- alpha_par_4[group_from_to_sbm]
  alpha_par_4_sbm <- list(unlist(lapply(lst_sbm, function(x) x[1])),
                          unlist(lapply(lst_sbm, function(x) x[2])),
                          unlist(lapply(lst_sbm, function(x) x[3])))

  # ER model
  net_er <- erdos.renyi.game(20,p.or.m = 0.1,type = "gnp",directed = TRUE)
  V(net_er)$name <- as.character(seq(1,20,1))
  V(net_er)$group <- sample(c('g1','g2'), size=20, replace=TRUE)
  edgelist_er <- get.edgelist(net_er)
  nedges_er <- ecount(net_er)
  group_from_er <- V(net_er)$group[as.numeric(edgelist_er[,1])]
  group_to_er   <- V(net_er)$group[as.numeric(edgelist_er[,2])]
  group_from_to_er <- paste(group_from_er, group_to_er, sep='_')
  grouplist_er <- split(seq_along(group_from_to_er), group_from_to_er)
  lst_er <- alpha_par_4[group_from_to_er]
  alpha_par_4_er <- list(unlist(lapply(lst_er, function(x) x[1])),
                         unlist(lapply(lst_er, function(x) x[2])),
                         unlist(lapply(lst_er, function(x) x[3])))

  # Random Geometric Graph
  lpvs2 <- sample_sphere_surface(dim = 2, n = 20, radius = .35)
  net_grg <- sample_dot_product(lpvs2, directed = TRUE)
  V(net_grg)$name <- as.character(seq(1,20,1))
  V(net_grg)$group <- sample(c('g1','g2'), size=20, replace=TRUE)
  edgelist_grg <- get.edgelist(net_grg)
  nedges_grg <- ecount(net_grg)
  group_from_grg <- V(net_grg)$group[as.numeric(edgelist_grg[,1])]
  group_to_grg   <- V(net_grg)$group[as.numeric(edgelist_grg[,2])]
  group_from_to_grg <- paste(group_from_grg, group_to_grg, sep='_')
  grouplist_grg <- split(seq_along(group_from_to_grg), group_from_to_grg)
  lst_grg <- alpha_par_4[group_from_to_grg]
  alpha_par_4_grg <- list(unlist(lapply(lst_grg, function(x) x[1])),
                          unlist(lapply(lst_grg, function(x) x[2])),
                          unlist(lapply(lst_grg, function(x) x[3])))

  # Simulate time series on ER, SBM, GRG
  ts_er <- gnar_edge_sim(n=200,net=net_er,alphaParams=alpha_par_4_er,betaParams=beta_par_4,
                         sigma=1,meann=0,nedges=nedges_er,data_edges=edgelist_er,groups=grouplist_er)
  datasim_train_er <- ts_er[,1:199]

  ts_sbm <- gnar_edge_sim(n=200,net=net_sbm,alphaParams=alpha_par_4_sbm,betaParams=beta_par_4,
                          sigma=1,meann=0,nedges=nedges_sbm,data_edges=edgelist_sbm,groups=grouplist_sbm)
  datasim_train_sbm <- ts_sbm[,1:199]

  ts_grg <- gnar_edge_sim(n=200,net=net_grg,alphaParams=alpha_par_4_grg,betaParams=beta_par_4,
                          sigma=1,meann=0,nedges=nedges_grg,data_edges=edgelist_grg,groups=grouplist_grg)
  datasim_train_grg <- ts_grg[,1:199]

  for(j in c("er","sbm","grg")){
    datasim_train <- get(paste("datasim_train_",j,sep=""))
    edgelist      <- get(paste("edgelist_",j,sep=""))
    ts            <- get(paste("ts_",j,sep=""))
    net           <- get(paste("net_",j,sep=""))
    grouplist     <- get(paste("grouplist_",j,sep=""))

    # Fit GNAR(3,[2,2,2])
    fit_train <- gnar_edge_fit(datasim_train,edgelist,3,rep(2,3),net,lead_lag_mat=NULL,groups=grouplist,globalalpha=TRUE,lead_lag_weights=FALSE)
    fit_pred  <- gnar_edge_predict(fit_train,datasim_train,3,rep(2,3),nrow(edgelist),1,fit_train$wei_mat,groups=grouplist)
    rmse_mat[i,paste("gnarnei_",j,sep="")] <- sqrt(mean((fit_pred[,4]-ts[,200])^2))
    print(c("seed: ",i," gnar nei done"))

    # Fit GNAR(3,[0,0,0])
    simplevar_sim_data_1 <- gnar_edge_fit(datasim_train,edgelist,3,rep(0,3),net,lead_lag_mat=NULL,groups=grouplist,globalalpha=TRUE,lead_lag_weights=FALSE)
    simplevar_pred <- gnar_edge_predict(simplevar_sim_data_1,datasim_train,3,rep(0,3),nrow(edgelist),1,simplevar_sim_data_1$wei_mat,groups=grouplist)
    rmse_mat[i,paste("gnar_",j,sep="")] <- sqrt(mean((simplevar_pred[,4]-ts[,200])^2))
    print(c("seed: ",i," gnar nonei done"))

    # Fit VAR(3)
    simtrainvar <- t(datasim_train)
    varforecast <- predict(restrict(VAR(simtrainvar,p=3,type="none")),n.ahead=1)
    getfcst <- function(x){return(x[,1])}
    varfor <- unlist(lapply(varforecast$fcst, getfcst))
    rmse_mat[i,paste("var_",j,sep="")] <- sqrt(mean((varfor-ts[,200])^2))
    print(c("seed: ",i," var done"))

    # Fit AR max lag 3
    simple_ar <- apply(simtrainvar, 2, function(x){forecast(auto.arima(x,d=0,D=0,max.p=3,max.q=0,max.P=0,max.Q=0,stationary=TRUE,seasonal=FALSE,
                                                                        ic="bic",allowmean=FALSE,allowdrift=FALSE,trace=FALSE),h=1)$mean})
    rmse_mat[i,paste("ar_",j,sep="")] <- sqrt(mean((simple_ar-ts[,200])^2))
    print(c("seed: ",i," ar done"))
  }
}

## PLOTS RESULTS
rmse_df <- as.data.frame(rmse_mat)
colnames(rmse_df) <- rep(c("GNAR(3,[2,2,2])","GNAR(3,[0,0,0])","VAR","AR"),3)
rmse_melt <- melt(rmse_df)
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse")
p1 <- ggplot(rmse_melt, aes(x=m, y=rmse)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.68,1.93)+ggtitle("ER")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(rmse_melt, aes(x=m, y=rmse)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.68,1.93)+ggtitle("SBM")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
p3 <- ggplot(rmse_melt, aes(x=m, y=rmse)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.68,1.93)+ggtitle("RDP")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")
library(patchwork)
p1+p2+p3


#####################################################################################
# Density sweep: different graph densities, comparing neighbour vs no-neighbour
#####################################################################################

# ER model — nedges = 38, 76, 114, 152

set.seed(44)
seeds <- sample(1:500,50,replace=FALSE)
rmse_list_nei   <- vector(mode="list", length=4)
rmse_list_nonei <- vector(mode="list", length=4)
ind <- 1
for (ne in c(38,76,114,152)){
  for (i in 1:50){
    print(c(ne," ",i))
    set.seed(seeds[i])
    net_er <- erdos.renyi.game(20,p.or.m=ne,type="gnm",directed=TRUE)
    V(net_er)$name <- as.character(seq(1,20,1))
    V(net_er)$group <- sample(c('g1','g2'), size=20, replace=TRUE)
    edgelist_er <- get.edgelist(net_er)
    group_from_er    <- V(net_er)$group[as.numeric(edgelist_er[,1])]
    group_to_er      <- V(net_er)$group[as.numeric(edgelist_er[,2])]
    group_from_to_er <- paste(group_from_er, group_to_er, sep='_')
    grouplist_er     <- split(seq_along(group_from_to_er), group_from_to_er)
    lst_er <- alpha_par_4[group_from_to_er]
    alpha_par_4_er <- list(unlist(lapply(lst_er, function(x) x[1])),
                           unlist(lapply(lst_er, function(x) x[2])),
                           unlist(lapply(lst_er, function(x) x[3])))

    simdata <- gnar_edge_sim(n=200,net=net_er,alphaParams=alpha_par_4_er,betaParams=beta_par_4,
                             sigma=1,meann=0,nedges=ne,data_edges=edgelist_er,groups=grouplist_er)
    datasim_train <- simdata[,1:199]

    fit_train <- gnar_edge_fit(datasim_train,edgelist_er,3,rep(2,3),net_er,lead_lag_mat=NULL,groups=grouplist_er,globalalpha=TRUE,lead_lag_weights=FALSE)
    fit_pred  <- gnar_edge_predict(fit_train,datasim_train,3,rep(2,3),nrow(edgelist_er),1,fit_train$wei_mat,groups=grouplist_er)
    rmse_list_nei[[ind]][i] <- sqrt(mean((fit_pred[,4]-simdata[,200])^2))

    fit_train2     <- gnar_edge_fit(datasim_train,edgelist_er,3,rep(0,3),net_er,lead_lag_mat=NULL,groups=grouplist_er,globalalpha=TRUE)
    simplevar_pred <- gnar_edge_predict(fit_train2,datasim_train,3,rep(0,3),nrow(edgelist_er),1,fit_train2$wei_mat,groups=grouplist_er)
    rmse_list_nonei[[ind]][i] <- sqrt(mean((simplevar_pred[,4]-simdata[,200])^2))
  }
  ind <- ind+1
}

rmse_df_nei   <- do.call(cbind, rmse_list_nei)
rmse_df_nonei <- do.call(cbind, rmse_list_nonei)
rmse_df_nonei <- as.data.frame(rmse_df_nonei)
rmse_df_nonei["model"] <- rep("no neighbour",50)
colnames(rmse_df_nonei) <- c("0.1","0.2","0.3","0.4","model")
rmse_df_nei <- as.data.frame(rmse_df_nei)
rmse_df_nei["model"] <- rep("neighbour",50)
colnames(rmse_df_nei) <- c("0.1","0.2","0.3","0.4","model")
rmseall <- rbind(rmse_df_nei,rmse_df_nonei)
rmseall_melt <- melt(rmseall)
colnames(rmseall_melt) <- c("model","density","rmse")
ggplot(rmseall_melt, aes(x=density, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("ER")+theme(plot.title = element_text(hjust = 0.5))


# SBM model

set.seed(44)
seeds <- sample(1:500,50,replace=FALSE)
probmat <- list(matrix(c(.2,.02,.02,.2),nrow=2),matrix(c(.4,.05,.05,.4),nrow=2),
                matrix(c(.5,.1,.1,.5),nrow=2),matrix(c(.7,.2,.1,.7),nrow=2))
rmse_list_nei_sbm   <- vector(mode="list", length=4)
rmse_list_nonei_sbm <- vector(mode="list", length=4)
for (ne in 1:4){
  for (i in 1:50){
    print(c(ne," ",i))
    set.seed(seeds[i])
    net_sbm <- sample_sbm(20,probmat[[ne]],c(10,10),directed=TRUE)
    V(net_sbm)$name <- as.character(seq(1,20,1))
    V(net_sbm)$group <- sample(c('g1','g2'), size=20, replace=TRUE)
    edgelist_sbm <- get.edgelist(net_sbm)
    group_from_sbm    <- V(net_sbm)$group[as.numeric(edgelist_sbm[,1])]
    group_to_sbm      <- V(net_sbm)$group[as.numeric(edgelist_sbm[,2])]
    group_from_to_sbm <- paste(group_from_sbm, group_to_sbm, sep='_')
    grouplist_sbm     <- split(seq_along(group_from_to_sbm), group_from_to_sbm)
    lst_sbm <- alpha_par_4[group_from_to_sbm]
    alpha_par_4_sbm <- list(unlist(lapply(lst_sbm, function(x) x[1])),
                            unlist(lapply(lst_sbm, function(x) x[2])),
                            unlist(lapply(lst_sbm, function(x) x[3])))

    simdata <- gnar_edge_sim(n=200,net=net_sbm,alphaParams=alpha_par_4_sbm,betaParams=beta_par_4,
                             sigma=1,meann=0,nedges=ecount(net_sbm),data_edges=edgelist_sbm,groups=grouplist_sbm)
    datasim_train <- simdata[,1:199]

    fit_train <- gnar_edge_fit(datasim_train,edgelist_sbm,3,rep(2,3),net_sbm,lead_lag_mat=NULL,groups=grouplist_sbm,globalalpha=TRUE,lead_lag_weights=FALSE)
    fit_pred  <- gnar_edge_predict(fit_train,datasim_train,3,rep(2,3),nrow(edgelist_sbm),1,fit_train$wei_mat,groups=grouplist_sbm)
    rmse_list_nei_sbm[[ne]][i] <- sqrt(mean((fit_pred[,4]-simdata[,200])^2))

    fit_train2     <- gnar_edge_fit(datasim_train,edgelist_sbm,3,rep(0,3),net_sbm,lead_lag_mat=NULL,groups=grouplist_sbm,globalalpha=TRUE,lead_lag_weights=FALSE)
    simplevar_pred <- gnar_edge_predict(fit_train2,datasim_train,3,rep(0,3),nrow(edgelist_sbm),1,fit_train2$wei_mat,groups=grouplist_sbm)
    rmse_list_nonei_sbm[[ne]][i] <- sqrt(mean((simplevar_pred[,4]-simdata[,200])^2))
  }
}

rmse_df_nei   <- do.call(cbind, rmse_list_nei_sbm)
rmse_df_nonei <- do.call(cbind, rmse_list_nonei_sbm)
rmse_df_nonei <- as.data.frame(rmse_df_nonei)
rmse_df_nonei["model"] <- rep("no neighbour",50)
colnames(rmse_df_nonei) <- c("0.1","0.2","0.3","0.4","model")
rmse_df_nei <- as.data.frame(rmse_df_nei)
rmse_df_nei["model"] <- rep("neighbour",50)
colnames(rmse_df_nei) <- c("0.1","0.2","0.3","0.4","model")
rmseall <- rbind(rmse_df_nei,rmse_df_nonei)
rmseall_melt <- melt(rmseall)
colnames(rmseall_melt) <- c("model","density","rmse")
ggplot(rmseall_melt, aes(x=density, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("SBM")+theme(plot.title = element_text(hjust = 0.5))


# Random Dot Product model

set.seed(44)
seeds <- sample(1:500,50,replace=FALSE)
radius <- c(.35,0.5,0.61,.7)
rmse_list_nei_rdp   <- vector(mode="list", length=4)
rmse_list_nonei_rdp <- vector(mode="list", length=4)
for (ne in 1:4){
  for (i in 1:50){
    print(c(ne," ",i))
    set.seed(seeds[i])
    lpvs2   <- sample_sphere_surface(dim=2, n=20, radius=radius[ne])
    net_grg <- sample_dot_product(lpvs2, directed=TRUE)
    V(net_grg)$name <- as.character(seq(1,20,1))
    V(net_grg)$group <- sample(c('g1','g2'), size=20, replace=TRUE)
    edgelist_grg <- get.edgelist(net_grg)
    nedges_grg   <- ecount(net_grg)
    group_from_grg    <- V(net_grg)$group[as.numeric(edgelist_grg[,1])]
    group_to_grg      <- V(net_grg)$group[as.numeric(edgelist_grg[,2])]
    group_from_to_grg <- paste(group_from_grg, group_to_grg, sep='_')
    grouplist_grg     <- split(seq_along(group_from_to_grg), group_from_to_grg)
    lst_grg <- alpha_par_4[group_from_to_grg]
    alpha_par_4_grg <- list(unlist(lapply(lst_grg, function(x) x[1])),
                            unlist(lapply(lst_grg, function(x) x[2])),
                            unlist(lapply(lst_grg, function(x) x[3])))

    simdata <- gnar_edge_sim(n=200,net=net_grg,alphaParams=alpha_par_4_grg,betaParams=beta_par_4,
                             sigma=1,meann=0,nedges=nedges_grg,data_edges=edgelist_grg,groups=grouplist_grg)
    datasim_train <- simdata[,1:199]

    fit_train <- gnar_edge_fit(datasim_train,edgelist_grg,3,rep(2,3),net_grg,lead_lag_mat=NULL,groups=grouplist_grg,globalalpha=TRUE,lead_lag_weights=FALSE)
    fit_pred  <- gnar_edge_predict(fit_train,datasim_train,3,rep(2,3),nrow(edgelist_grg),1,fit_train$wei_mat,groups=grouplist_grg)
    rmse_list_nei_rdp[[ne]][i] <- sqrt(mean((fit_pred[,4]-simdata[,200])^2))

    fit_train2     <- gnar_edge_fit(datasim_train,edgelist_grg,3,rep(0,3),net_grg,lead_lag_mat=NULL,groups=grouplist_grg,globalalpha=TRUE,lead_lag_weights=FALSE)
    simplevar_pred <- gnar_edge_predict(fit_train2,datasim_train,3,rep(0,3),nrow(edgelist_grg),1,fit_train2$wei_mat,groups=grouplist_grg)
    rmse_list_nonei_rdp[[ne]][i] <- sqrt(mean((simplevar_pred[,4]-simdata[,200])^2))
  }
}

rmse_df_nei   <- do.call(cbind, rmse_list_nei_rdp)
rmse_df_nonei <- do.call(cbind, rmse_list_nonei_rdp)
rmse_df_nonei <- as.data.frame(rmse_df_nonei)
rmse_df_nonei["model"] <- rep("no neighbour",50)
colnames(rmse_df_nonei) <- c("0.1","0.2","0.3","0.4","model")
rmse_df_nei <- as.data.frame(rmse_df_nei)
rmse_df_nei["model"] <- rep("neighbour",50)
colnames(rmse_df_nei) <- c("0.1","0.2","0.3","0.4","model")
rmseall <- rbind(rmse_df_nei,rmse_df_nonei)
rmseall_melt <- melt(rmseall)
colnames(rmseall_melt) <- c("model","density","rmse")
ggplot(rmseall_melt, aes(x=density, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("RDP")+theme(plot.title = element_text(hjust = 0.5))


####### SECTION 5.2.2 #######

###################################################################################################
################## LARGE NETWORKS ##################
###################################################################################################

###################################################################################################
########### PREDICTIVE PERFORMANCE: REPETITIONS OF SBM, ER, RDP — REGIME 1 ###########
###################################################################################################

# Regime 1: GNAR(4,[1,1,1,1]) — group-structured parameters (matches 2_fixed.R Section 5.1.2)
alpha_par_1 <- list(g1_g1=c(-0.6,-0.4,-0.2,-0.1), g1_g2=c(-0.3,0.4,0.1,0.2),
                    g2_g1=c(0.3,0.1,0.2,0.1),     g2_g2=c(0.1,-0.2,0.2,0.2))
beta_par_1 <- list(g1_g1=list(c(0.2),c(0.1),c(0.3),c(0.05)),
                   g1_g2=list(c(0.1),c(0.2),c(0.1),c(-0.2)),
                   g2_g1=list(c(0.05),c(-0.05),c(0.1),c(-0.1)),
                   g2_g2=list(c(0.1),c(0.2),c(-0.2),c(-0.1)))

set.seed(44)
seeds <- sample(1:500,50,replace=FALSE)
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 9,nrow = 50)
colnames(rmse_mat) <- c("gnarnei_er","gnar_er","ar_er",
                        "gnarnei_sbm","gnar_sbm","ar_sbm",
                        "gnarnei_grg","gnar_grg","ar_grg")

for (i in 1:50){
  print(i)
  set.seed(seeds[i])

  # SBM model
  net_sbm <- sample_sbm(86,probmat,c(43,43),directed=TRUE)
  V(net_sbm)$name <- as.character(seq(1,86,1))
  V(net_sbm)$group <- sample(c('g1','g2'), size=86, replace=TRUE)
  edgelist_sbm <- get.edgelist(net_sbm)
  nedges_sbm   <- ecount(net_sbm)
  group_from_sbm    <- V(net_sbm)$group[as.numeric(edgelist_sbm[,1])]
  group_to_sbm      <- V(net_sbm)$group[as.numeric(edgelist_sbm[,2])]
  group_from_to_sbm <- paste(group_from_sbm, group_to_sbm, sep='_')
  grouplist_sbm <- split(seq_along(group_from_to_sbm), group_from_to_sbm)
  lst_sbm <- alpha_par_1[group_from_to_sbm]
  alpha_par_1_sbm <- list(unlist(lapply(lst_sbm, function(x) x[1])),
                          unlist(lapply(lst_sbm, function(x) x[2])),
                          unlist(lapply(lst_sbm, function(x) x[3])),
                          unlist(lapply(lst_sbm, function(x) x[4])))

  # ER model
  net_er <- erdos.renyi.game(86,p.or.m=0.1,type="gnp",directed=TRUE)
  V(net_er)$name <- as.character(seq(1,86,1))
  V(net_er)$group <- sample(c('g1','g2'), size=86, replace=TRUE)
  edgelist_er <- get.edgelist(net_er)
  nedges_er   <- ecount(net_er)
  group_from_er    <- V(net_er)$group[as.numeric(edgelist_er[,1])]
  group_to_er      <- V(net_er)$group[as.numeric(edgelist_er[,2])]
  group_from_to_er <- paste(group_from_er, group_to_er, sep='_')
  grouplist_er <- split(seq_along(group_from_to_er), group_from_to_er)
  lst_er <- alpha_par_1[group_from_to_er]
  alpha_par_1_er <- list(unlist(lapply(lst_er, function(x) x[1])),
                         unlist(lapply(lst_er, function(x) x[2])),
                         unlist(lapply(lst_er, function(x) x[3])),
                         unlist(lapply(lst_er, function(x) x[4])))

  # Random Geometric Graph
  lpvs2   <- sample_sphere_surface(dim=2, n=86, radius=.35)
  net_grg <- sample_dot_product(lpvs2, directed=TRUE)
  V(net_grg)$name <- as.character(seq(1,86,1))
  V(net_grg)$group <- sample(c('g1','g2'), size=86, replace=TRUE)
  edgelist_grg <- get.edgelist(net_grg)
  nedges_grg   <- ecount(net_grg)
  group_from_grg    <- V(net_grg)$group[as.numeric(edgelist_grg[,1])]
  group_to_grg      <- V(net_grg)$group[as.numeric(edgelist_grg[,2])]
  group_from_to_grg <- paste(group_from_grg, group_to_grg, sep='_')
  grouplist_grg <- split(seq_along(group_from_to_grg), group_from_to_grg)
  lst_grg <- alpha_par_1[group_from_to_grg]
  alpha_par_1_grg <- list(unlist(lapply(lst_grg, function(x) x[1])),
                          unlist(lapply(lst_grg, function(x) x[2])),
                          unlist(lapply(lst_grg, function(x) x[3])),
                          unlist(lapply(lst_grg, function(x) x[4])))

  # Simulate time series
  ts_er <- gnar_edge_sim(n=90,net=net_er,alphaParams=alpha_par_1_er,betaParams=beta_par_1,
                         sigma=1,meann=0,nedges=nedges_er,data_edges=edgelist_er,groups=grouplist_er)
  datasim_train_er <- ts_er[,1:89]

  ts_sbm <- gnar_edge_sim(n=90,net=net_sbm,alphaParams=alpha_par_1_sbm,betaParams=beta_par_1,
                          sigma=1,meann=0,nedges=nedges_sbm,data_edges=edgelist_sbm,groups=grouplist_sbm)
  datasim_train_sbm <- ts_sbm[,1:89]

  ts_grg <- gnar_edge_sim(n=90,net=net_grg,alphaParams=alpha_par_1_grg,betaParams=beta_par_1,
                          sigma=1,meann=0,nedges=nedges_grg,data_edges=edgelist_grg,groups=grouplist_grg)
  datasim_train_grg <- ts_grg[,1:89]

  for(j in c("sbm","er","grg")){
    datasim_train <- get(paste("datasim_train_",j,sep=""))
    edgelist      <- get(paste("edgelist_",j,sep=""))
    ts            <- get(paste("ts_",j,sep=""))
    net           <- get(paste("net_",j,sep=""))
    grouplist     <- get(paste("grouplist_",j,sep=""))

    # Fit GNAR(4,[1,1,1,1])
    fit_train <- gnar_edge_fit(datasim_train,edgelist,4,rep(1,4),net,lead_lag_mat=NULL,groups=grouplist,globalalpha=TRUE,lead_lag_weights=FALSE)
    fit_pred  <- gnar_edge_predict(fit_train,datasim_train,4,rep(1,4),nrow(edgelist),1,fit_train$wei_mat,groups=grouplist)
    rmse_mat[i,paste("gnarnei_",j,sep="")] <- sqrt(mean((fit_pred[,5]-ts[,90])^2))
    print(c("seed: ",i," gnar nei done"))

    # Fit GNAR(4,[0,0,0,0])
    simplevar_sim_data_1 <- gnar_edge_fit(datasim_train,edgelist,4,rep(0,4),net,lead_lag_mat=NULL,groups=grouplist,globalalpha=TRUE,lead_lag_weights=FALSE)
    simplevar_pred <- gnar_edge_predict(simplevar_sim_data_1,datasim_train,4,rep(0,4),nrow(edgelist),1,simplevar_sim_data_1$wei_mat,groups=grouplist)
    rmse_mat[i,paste("gnar_",j,sep="")] <- sqrt(mean((simplevar_pred[,5]-ts[,90])^2))
    print(c("seed: ",i," gnar nonei done"))

    # Fit AR max lag 4 (VAR excluded at this scale)
    simtrainvar <- t(datasim_train)
    simple_ar <- apply(simtrainvar, 2, function(x){forecast(auto.arima(x,d=0,D=0,max.p=4,max.q=0,max.P=0,max.Q=0,stationary=TRUE,seasonal=FALSE,
                                                                        ic="bic",allowmean=FALSE,allowdrift=FALSE,trace=FALSE),h=1)$mean})
    rmse_mat[i,paste("ar_",j,sep="")] <- sqrt(mean((simple_ar-ts[,90])^2))
    print(c("seed: ",i," ar done"))
  }
}

# alias for stage-comparison plots below
rmse_mat_1stag <- rmse_mat

# for ER, SBM, RDP with patchwork
rmse_df <- as.data.frame(rmse_mat[,1:3])
colnames(rmse_df) <- c("GNAR(4,[1,1,1,1])","GNAR(4,[0,0,0,0])","AR")
rmse_melt <- melt(rmse_df)
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse")
p1 <- ggplot(rmse_melt, aes(x=m, y=rmse)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.93,1.1)+ggtitle("ER")+theme(plot.title = element_text(hjust = 0.5))
rmse_df <- as.data.frame(rmse_mat[,4:6])
colnames(rmse_df) <- c("GNAR(4,[1,1,1,1])","GNAR(4,[0,0,0,0])","AR")
rmse_melt <- melt(rmse_df)
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse")
p2 <- ggplot(rmse_melt, aes(x=m, y=rmse)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.93,1.1)+ggtitle("SBM")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
rmse_df <- as.data.frame(rmse_mat[,7:9])
colnames(rmse_df) <- c("GNAR(4,[1,1,1,1])","GNAR(4,[0,0,0,0])","AR")
rmse_melt <- melt(rmse_df)
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse")
p3 <- ggplot(rmse_melt, aes(x=m, y=rmse)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.93,1.1)+ggtitle("RDP")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")
library(patchwork)
p1+p2+p3

# ER large network plot
library(reshape2)
library(ggplot2)
rmse_df <- as.data.frame(rmse_mat[,1:3])
colnames(rmse_df) <- c("GNAR(4,[1,1,1,1])","GNAR(4,[0,0,0,0])","AR")
rmse_melt <- melt(rmse_df)
rmse_melt["model"] <- c(rep("ER",150))
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse","model")
ggplot(rmse_melt, aes(x=m, y=rmse, fill=model)) +
  geom_boxplot()+ theme(axis.text.x = element_text(angle = 60, hjust = 1),legend.position = "none",plot.title = element_text(hjust = 0.5))+xlab(" ")+ggtitle("ER")


###################################################################################################
################## DIFFERENT STAGE NEIGHBOURS ##################
###################################################################################################

#### 2-stage neighbours ####

# Regime 2: GNAR(4,[2,2,2,2]) — group-structured parameters
alpha_par_2 <- list(g1_g1=c(-0.6,-0.4,-0.2,-0.1), g1_g2=c(-0.3,0.4,0.1,0.2),
                    g2_g1=c(0.3,0.1,0.2,0.1),     g2_g2=c(0.1,-0.2,0.2,0.2))
beta_par_2 <- list(g1_g1=list(c(0.4,-0.4),c(0.3,-0.4),c(0.5,-0.3),c(0.05,-0.1)),
                   g1_g2=list(c(0.3,-0.3),c(0.2,-0.3),c(0.4,-0.2),c(0.04,-0.08)),
                   g2_g1=list(c(0.3,-0.3),c(0.2,-0.3),c(0.4,-0.2),c(0.04,-0.08)),
                   g2_g2=list(c(0.4,-0.4),c(0.3,-0.4),c(0.5,-0.3),c(0.05,-0.1)))

set.seed(44)
seeds <- sample(1:500,50,replace=FALSE)
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat_2stag <- matrix(ncol=4,nrow=50)
colnames(rmse_mat_2stag) <- c("gnarnei_er","gnar_er","gnarnei_sbm","gnar_sbm")

for (i in 1:50){
  print(i)
  set.seed(seeds[i])

  # SBM model
  net_sbm <- sample_sbm(86,probmat,c(43,43),directed=TRUE)
  V(net_sbm)$name <- as.character(seq(1,86,1))
  V(net_sbm)$group <- sample(c('g1','g2'), size=86, replace=TRUE)
  edgelist_sbm <- get.edgelist(net_sbm)
  nedges_sbm   <- ecount(net_sbm)
  group_from_sbm    <- V(net_sbm)$group[as.numeric(edgelist_sbm[,1])]
  group_to_sbm      <- V(net_sbm)$group[as.numeric(edgelist_sbm[,2])]
  group_from_to_sbm <- paste(group_from_sbm, group_to_sbm, sep='_')
  grouplist_sbm <- split(seq_along(group_from_to_sbm), group_from_to_sbm)
  lst_sbm <- alpha_par_2[group_from_to_sbm]
  alpha_par_2_sbm <- list(unlist(lapply(lst_sbm, function(x) x[1])),
                          unlist(lapply(lst_sbm, function(x) x[2])),
                          unlist(lapply(lst_sbm, function(x) x[3])),
                          unlist(lapply(lst_sbm, function(x) x[4])))

  # ER model
  net_er <- erdos.renyi.game(86,p.or.m=0.1,type="gnp",directed=TRUE)
  V(net_er)$name <- as.character(seq(1,86,1))
  V(net_er)$group <- sample(c('g1','g2'), size=86, replace=TRUE)
  edgelist_er <- get.edgelist(net_er)
  nedges_er   <- ecount(net_er)
  group_from_er    <- V(net_er)$group[as.numeric(edgelist_er[,1])]
  group_to_er      <- V(net_er)$group[as.numeric(edgelist_er[,2])]
  group_from_to_er <- paste(group_from_er, group_to_er, sep='_')
  grouplist_er <- split(seq_along(group_from_to_er), group_from_to_er)
  lst_er <- alpha_par_2[group_from_to_er]
  alpha_par_2_er <- list(unlist(lapply(lst_er, function(x) x[1])),
                         unlist(lapply(lst_er, function(x) x[2])),
                         unlist(lapply(lst_er, function(x) x[3])),
                         unlist(lapply(lst_er, function(x) x[4])))

  # Simulate
  ts_er <- gnar_edge_sim(n=90,net=net_er,alphaParams=alpha_par_2_er,betaParams=beta_par_2,
                         sigma=1,meann=0,nedges=nedges_er,data_edges=edgelist_er,groups=grouplist_er)
  datasim_train_er <- ts_er[,1:89]

  ts_sbm <- gnar_edge_sim(n=90,net=net_sbm,alphaParams=alpha_par_2_sbm,betaParams=beta_par_2,
                          sigma=1,meann=0,nedges=nedges_sbm,data_edges=edgelist_sbm,groups=grouplist_sbm)
  datasim_train_sbm <- ts_sbm[,1:89]

  for(j in c("er","sbm")){
    datasim_train <- get(paste("datasim_train_",j,sep=""))
    edgelist      <- get(paste("edgelist_",j,sep=""))
    ts            <- get(paste("ts_",j,sep=""))
    net           <- get(paste("net_",j,sep=""))
    grouplist     <- get(paste("grouplist_",j,sep=""))

    # Fit GNAR(4,[2,2,2,2])
    fit_train <- gnar_edge_fit(datasim_train,edgelist,4,rep(2,4),net,lead_lag_mat=NULL,groups=grouplist,globalalpha=TRUE,lead_lag_weights=FALSE)
    fit_pred  <- gnar_edge_predict(fit_train,datasim_train,4,rep(2,4),nrow(edgelist),1,fit_train$wei_mat,groups=grouplist)
    rmse_mat_2stag[i,paste("gnarnei_",j,sep="")] <- sqrt(mean((fit_pred[,5]-ts[,90])^2))
    print(c("seed: ",i," gnar nei done"))

    # Fit GNAR(4,[0,0,0,0])
    simplevar_sim_data_1 <- gnar_edge_fit(datasim_train,edgelist,4,rep(0,4),net,lead_lag_mat=NULL,groups=grouplist,globalalpha=TRUE,lead_lag_weights=FALSE)
    simplevar_pred <- gnar_edge_predict(simplevar_sim_data_1,datasim_train,4,rep(0,4),nrow(edgelist),1,simplevar_sim_data_1$wei_mat,groups=grouplist)
    rmse_mat_2stag[i,paste("gnar_",j,sep="")] <- sqrt(mean((simplevar_pred[,5]-ts[,90])^2))
    print(c("seed: ",i," gnar nonei done"))
  }
}


#### 3-stage neighbours ####

# Regime 3: GNAR(4,[3,3,3,3]) — group-structured parameters
alpha_par_3 <- list(g1_g1=c(-0.6,-0.4,-0.2,-0.1), g1_g2=c(-0.3,0.4,0.1,0.2),
                    g2_g1=c(0.3,0.1,0.2,0.1),     g2_g2=c(0.1,-0.2,0.2,0.2))
beta_par_3 <- list(g1_g1=list(c(0.4,-0.4,-.2),c(0.3,-0.4,.2),c(0.5,-0.3,-0.4),c(0.05,-0.1,0.3)),
                   g1_g2=list(c(0.3,-0.3,-.15),c(0.2,-0.3,.15),c(0.4,-0.2,-0.3),c(0.04,-0.08,0.2)),
                   g2_g1=list(c(0.3,-0.3,-.15),c(0.2,-0.3,.15),c(0.4,-0.2,-0.3),c(0.04,-0.08,0.2)),
                   g2_g2=list(c(0.4,-0.4,-.2),c(0.3,-0.4,.2),c(0.5,-0.3,-0.4),c(0.05,-0.1,0.3)))

set.seed(44)
seeds <- sample(1:500,50,replace=FALSE)
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat_3stag <- matrix(ncol=4,nrow=50)
colnames(rmse_mat_3stag) <- c("gnarnei_er","gnar_er","gnarnei_sbm","gnar_sbm")

for (i in 1:50){
  print(i)
  set.seed(seeds[i])

  # SBM model
  net_sbm <- sample_sbm(86,probmat,c(43,43),directed=TRUE)
  V(net_sbm)$name <- as.character(seq(1,86,1))
  V(net_sbm)$group <- sample(c('g1','g2'), size=86, replace=TRUE)
  edgelist_sbm <- get.edgelist(net_sbm)
  nedges_sbm   <- ecount(net_sbm)
  group_from_sbm    <- V(net_sbm)$group[as.numeric(edgelist_sbm[,1])]
  group_to_sbm      <- V(net_sbm)$group[as.numeric(edgelist_sbm[,2])]
  group_from_to_sbm <- paste(group_from_sbm, group_to_sbm, sep='_')
  grouplist_sbm <- split(seq_along(group_from_to_sbm), group_from_to_sbm)
  lst_sbm <- alpha_par_3[group_from_to_sbm]
  alpha_par_3_sbm <- list(unlist(lapply(lst_sbm, function(x) x[1])),
                          unlist(lapply(lst_sbm, function(x) x[2])),
                          unlist(lapply(lst_sbm, function(x) x[3])),
                          unlist(lapply(lst_sbm, function(x) x[4])))

  # ER model
  net_er <- erdos.renyi.game(86,p.or.m=0.1,type="gnp",directed=TRUE)
  V(net_er)$name <- as.character(seq(1,86,1))
  V(net_er)$group <- sample(c('g1','g2'), size=86, replace=TRUE)
  edgelist_er <- get.edgelist(net_er)
  nedges_er   <- ecount(net_er)
  group_from_er    <- V(net_er)$group[as.numeric(edgelist_er[,1])]
  group_to_er      <- V(net_er)$group[as.numeric(edgelist_er[,2])]
  group_from_to_er <- paste(group_from_er, group_to_er, sep='_')
  grouplist_er <- split(seq_along(group_from_to_er), group_from_to_er)
  lst_er <- alpha_par_3[group_from_to_er]
  alpha_par_3_er <- list(unlist(lapply(lst_er, function(x) x[1])),
                         unlist(lapply(lst_er, function(x) x[2])),
                         unlist(lapply(lst_er, function(x) x[3])),
                         unlist(lapply(lst_er, function(x) x[4])))

  # Simulate
  ts_er <- gnar_edge_sim(n=90,net=net_er,alphaParams=alpha_par_3_er,betaParams=beta_par_3,
                         sigma=1,meann=0,nedges=nedges_er,data_edges=edgelist_er,groups=grouplist_er)
  datasim_train_er <- ts_er[,1:89]

  ts_sbm <- gnar_edge_sim(n=90,net=net_sbm,alphaParams=alpha_par_3_sbm,betaParams=beta_par_3,
                          sigma=1,meann=0,nedges=nedges_sbm,data_edges=edgelist_sbm,groups=grouplist_sbm)
  datasim_train_sbm <- ts_sbm[,1:89]

  for(j in c("er","sbm")){
    datasim_train <- get(paste("datasim_train_",j,sep=""))
    edgelist      <- get(paste("edgelist_",j,sep=""))
    ts            <- get(paste("ts_",j,sep=""))
    net           <- get(paste("net_",j,sep=""))
    grouplist     <- get(paste("grouplist_",j,sep=""))

    # Fit GNAR(4,[3,3,3,3])
    fit_train <- gnar_edge_fit(datasim_train,edgelist,4,rep(3,4),net,lead_lag_mat=NULL,groups=grouplist,globalalpha=TRUE,lead_lag_weights=FALSE)
    fit_pred  <- gnar_edge_predict(fit_train,datasim_train,4,rep(3,4),nrow(edgelist),1,fit_train$wei_mat,groups=grouplist)
    rmse_mat_3stag[i,paste("gnarnei_",j,sep="")] <- sqrt(mean((fit_pred[,5]-ts[,90])^2))
    print(c("seed: ",i," gnar nei done"))

    # Fit GNAR(4,[0,0,0,0])
    simplevar_sim_data_1 <- gnar_edge_fit(datasim_train,edgelist,4,rep(0,4),net,lead_lag_mat=NULL,groups=grouplist,globalalpha=TRUE,lead_lag_weights=FALSE)
    simplevar_pred <- gnar_edge_predict(simplevar_sim_data_1,datasim_train,4,rep(0,4),nrow(edgelist),1,simplevar_sim_data_1$wei_mat,groups=grouplist)
    rmse_mat_3stag[i,paste("gnar_",j,sep="")] <- sqrt(mean((simplevar_pred[,5]-ts[,90])^2))
    print(c("seed: ",i," gnar nonei done"))
  }
}


#### 4-stage neighbours ####

# Regime 4 (large): GNAR(4,[4,4,4,4]) — group-structured parameters
alpha_par_4 <- list(g1_g1=c(-0.6,-0.4,-0.2,-0.1), g1_g2=c(-0.3,0.4,0.1,0.2),
                    g2_g1=c(0.3,0.1,0.2,0.1),     g2_g2=c(0.1,-0.2,0.2,0.2))
beta_par_4 <- list(g1_g1=list(c(0.4,-0.4,-.2,-.1),c(0.3,-0.4,.2,.2),c(0.5,-0.3,-0.4,.3),c(0.05,-0.1,0.3,-.4)),
                   g1_g2=list(c(0.3,-0.3,-.15,-.08),c(0.2,-0.3,.15,.15),c(0.4,-0.2,-0.3,.2),c(0.04,-0.08,0.2,-.3)),
                   g2_g1=list(c(0.3,-0.3,-.15,-.08),c(0.2,-0.3,.15,.15),c(0.4,-0.2,-0.3,.2),c(0.04,-0.08,0.2,-.3)),
                   g2_g2=list(c(0.4,-0.4,-.2,-.1),c(0.3,-0.4,.2,.2),c(0.5,-0.3,-0.4,.3),c(0.05,-0.1,0.3,-.4)))

set.seed(44)
seeds <- sample(1:500,50,replace=FALSE)
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat_4stag <- matrix(ncol=4,nrow=50)
colnames(rmse_mat_4stag) <- c("gnarnei_er","gnar_er","gnarnei_sbm","gnar_sbm")

for (i in 1:50){
  print(i)
  set.seed(seeds[i])

  # SBM model
  net_sbm <- sample_sbm(86,probmat,c(43,43),directed=TRUE)
  V(net_sbm)$name <- as.character(seq(1,86,1))
  V(net_sbm)$group <- sample(c('g1','g2'), size=86, replace=TRUE)
  edgelist_sbm <- get.edgelist(net_sbm)
  nedges_sbm   <- ecount(net_sbm)
  group_from_sbm    <- V(net_sbm)$group[as.numeric(edgelist_sbm[,1])]
  group_to_sbm      <- V(net_sbm)$group[as.numeric(edgelist_sbm[,2])]
  group_from_to_sbm <- paste(group_from_sbm, group_to_sbm, sep='_')
  grouplist_sbm <- split(seq_along(group_from_to_sbm), group_from_to_sbm)
  lst_sbm <- alpha_par_4[group_from_to_sbm]
  alpha_par_4_sbm <- list(unlist(lapply(lst_sbm, function(x) x[1])),
                          unlist(lapply(lst_sbm, function(x) x[2])),
                          unlist(lapply(lst_sbm, function(x) x[3])),
                          unlist(lapply(lst_sbm, function(x) x[4])))

  # ER model
  net_er <- erdos.renyi.game(86,p.or.m=0.1,type="gnp",directed=TRUE)
  V(net_er)$name <- as.character(seq(1,86,1))
  V(net_er)$group <- sample(c('g1','g2'), size=86, replace=TRUE)
  edgelist_er <- get.edgelist(net_er)
  nedges_er   <- ecount(net_er)
  group_from_er    <- V(net_er)$group[as.numeric(edgelist_er[,1])]
  group_to_er      <- V(net_er)$group[as.numeric(edgelist_er[,2])]
  group_from_to_er <- paste(group_from_er, group_to_er, sep='_')
  grouplist_er <- split(seq_along(group_from_to_er), group_from_to_er)
  lst_er <- alpha_par_4[group_from_to_er]
  alpha_par_4_er <- list(unlist(lapply(lst_er, function(x) x[1])),
                         unlist(lapply(lst_er, function(x) x[2])),
                         unlist(lapply(lst_er, function(x) x[3])),
                         unlist(lapply(lst_er, function(x) x[4])))

  # Simulate
  ts_er <- gnar_edge_sim(n=90,net=net_er,alphaParams=alpha_par_4_er,betaParams=beta_par_4,
                         sigma=1,meann=0,nedges=nedges_er,data_edges=edgelist_er,groups=grouplist_er)
  datasim_train_er <- ts_er[,1:89]

  ts_sbm <- gnar_edge_sim(n=90,net=net_sbm,alphaParams=alpha_par_4_sbm,betaParams=beta_par_4,
                          sigma=1,meann=0,nedges=nedges_sbm,data_edges=edgelist_sbm,groups=grouplist_sbm)
  datasim_train_sbm <- ts_sbm[,1:89]

  for(j in c("er","sbm")){
    datasim_train <- get(paste("datasim_train_",j,sep=""))
    edgelist      <- get(paste("edgelist_",j,sep=""))
    ts            <- get(paste("ts_",j,sep=""))
    net           <- get(paste("net_",j,sep=""))
    grouplist     <- get(paste("grouplist_",j,sep=""))

    # Fit GNAR(4,[4,4,4,4])
    fit_train <- gnar_edge_fit(datasim_train,edgelist,4,rep(4,4),net,lead_lag_mat=NULL,groups=grouplist,globalalpha=TRUE,lead_lag_weights=FALSE)
    fit_pred  <- gnar_edge_predict(fit_train,datasim_train,4,rep(4,4),nrow(edgelist),1,fit_train$wei_mat,groups=grouplist)
    rmse_mat_4stag[i,paste("gnarnei_",j,sep="")] <- sqrt(mean((fit_pred[,5]-ts[,90])^2))
    print(c("seed: ",i," gnar nei done"))

    # Fit GNAR(4,[0,0,0,0])
    simplevar_sim_data_1 <- gnar_edge_fit(datasim_train,edgelist,4,rep(0,4),net,lead_lag_mat=NULL,groups=grouplist,globalalpha=TRUE,lead_lag_weights=FALSE)
    simplevar_pred <- gnar_edge_predict(simplevar_sim_data_1,datasim_train,4,rep(0,4),nrow(edgelist),1,simplevar_sim_data_1$wei_mat,groups=grouplist)
    rmse_mat_4stag[i,paste("gnar_",j,sep="")] <- sqrt(mean((simplevar_pred[,5]-ts[,90])^2))
    print(c("seed: ",i," gnar nonei done"))
  }
}


# side by side, various stages, ER
rmse_all_df <- data.frame("1"=c(rmse_mat_1stag[,1],rmse_mat_1stag[,2]),"2"=c(rmse_mat_2stag[,1],rmse_mat_2stag[,2]),
                          "3"=c(rmse_mat_3stag[,1],rmse_mat_3stag[,2]),"4"=c(rmse_mat_4stag[,1],rmse_mat_4stag[,2]),
                          "model"=rep(c("neighbour","no neighbour"),400, each = 50) )
rmseall_melt <- melt(rmse_all_df)
colnames(rmseall_melt) <- c("model","stage","rmse")
rmseall_melt$stage <- as.factor(rmseall_melt$stage)
ggplot(rmseall_melt, aes(x=stage, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("ER")+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(labels=c("1","2","3","4"),name="neighbour stage")

# side by side, various stages, SBM
rmse_all_df <- data.frame("1"=c(rmse_mat_1stag[,4],rmse_mat_1stag[,5]),"2"=c(rmse_mat_2stag[,3],rmse_mat_2stag[,4]),
                          "3"=c(rmse_mat_3stag[,3],rmse_mat_3stag[,4]),"4"=c(rmse_mat_4stag[,3],rmse_mat_4stag[,4]),
                          "model"=rep(c("neighbour","no neighbour"),400, each = 50) )
rmseall_melt <- melt(rmse_all_df)
colnames(rmseall_melt) <- c("model","stage","rmse")
rmseall_melt$stage <- as.factor(rmseall_melt$stage)
ggplot(rmseall_melt, aes(x=stage, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("SBM")+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(labels=c("1","2","3","4"),name="neighbour stage")
