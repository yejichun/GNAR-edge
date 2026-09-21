# modified GNAR-edge model to include groups

library(igraph)
library(vars)
library(forecast)
library(reshape2)
library(ggplot2)
library(dplyr)

########################################################
####### SIMULATION EXPERIMENTS FOR SECTION 5.1 #######
########################################################

####### SECTION 5.1.1 #######

# Global alpha
# Four cases for each network structure: GNAR(1,[1]), GNAR(1,[2]), GNAR(3,[1,1,1]), GNAR(3,[2,2,2])

# Regimes : Parameter specification

# G = 2
# regime_g1_g2
# Regime 1: GNAR(1,[1]) 
alpha_par_1 <- list(g1_g1=.2, g1_g2=.3, g2_g1=.4, g2_g2=.2)
beta_par_1 <- list(g1_g1=list(c(0.3)), g1_g2=list(c(0.2)), g2_g1=list(c(-0.1)), g2_g2=list(c(0.2)))

# Regime 2: GNAR(1,[2])
alpha_par_2 <- list(g1_g1=.2, g1_g2=.3, g2_g1=-0.2, g2_g2=.1)
beta_par_2 <- list(g1_g1=list(c(0.3,0.4)), g1_g2=list(c(0.2,0.1)), g2_g1=list(c(-0.2,-0.4)), g2_g2=list(c(-0.3,0.1)))

# Regime 3: GNAR(3,[1,1,1])
alpha_par_3 <- list(g1_g1=c(.2,.4,-0.6), g1_g2=c(.3,.1,-0.5), g2_g1=c(.4,.1,-0.5), g2_g2=c(.1,.3,-0.2))
beta_par_3 <- list(g1_g1=list(c(0.2),c(0.1),c(-.2)),
                   g1_g2=list(c(-0.1),c(-0.3),c(0.4)),
                   g2_g1=list(c(0.3),c(0.1),c(-0.4)),
                   g2_g2=list(c(0.2),c(0.4),c(-.2))) 

# Regime 4: GNAR(3,[2,2,2])
# alpha_par_4<- list(g1_g1=c(.2,.4,-0.6), g1_g2=c(.3,.1,-0.5), g2_g1=c(.4,.1,-0.5), g2_g2=c(.1,.3,-0.2)) # 0 -0.1 0 0.2
# beta_par_4 <- list(g1_g1=list(c(0.1,-.1),c(0.1,-.1),c(-0.1,0.1)), # 0.1
#                    g1_g2=list(c(0.1,-.1),c(0.1,-.1),c(-0.1,0.1)), # 0
#                    g2_g1=list(c(0.1,-.1),c(0.1,-.1),c(-0.1,0.1)), # -0.1
#                    g2_g2=list(c(0.1,-.1),c(0.1,-.1),c(-0.1,0.1)) # 0.1
#                    )

alpha_par_4<- list(g1_g1=c(.2,.4,-0.6), g1_g2=c(.3,.1,-0.5), g2_g1=c(.4,.1,-0.5), g2_g2=c(.1,.3,-0.2)) # 0 -0.1 0 0.2
beta_par_4 <- list(g1_g1=list(c(0.3,.1),c(0.1,.1),c(0,-0.2)), # 0.1
                   g1_g2=list(c(-0.1,-0.1),c(-0.1,0.3),c(-0.3,-0.1)), # 0
                   g2_g1=list(c(0.3,-.1),c(-0.1,-.1),c(0.2,.3)), # -0.1
                   g2_g2=list(c(-0.3,0.1),c(-0.2,.3),c(0.2,-.3)) # 0.1
)

# beta_par_4 <- list(g1_g1=list(c(0.1,-.1),c(0.1,-.1),c(-0.1,0.1)), # 0.1
#                    g1_g2=list(c(0.1,-.1),c(0.1,-.1),c(-0.1,0.1)), # 0
#                    g2_g1=list(c(0.3,.1),c(0.1,.1),c(-0.2,-0.3)), # -0.1
#                    g2_g2=list(c(-0.3,-0.1),c(-0.1,-.1),c(0.2,.3))# 0.1
# )

# lag x의 neighbor 값들을 더했을 때 1보다 작아야 함

# Regime 5: GNAR(3,[2,0,0])
alpha_par_5 <- list(g1_g1=c(.2,.4,-0.6), g1_g2=c(.3,.1,-0.5), g2_g1=c(.4,.1,-0.5), g2_g2=c(.1,.3,-0.2)) # 0 -0.1 0 0.1
beta_par_5 <- list(g1_g1=list(c(0.3,.4),c(0),c(0)),
                   g1_g2=list(c(0.2,-0.1),c(0),c(0)),
                   g2_g1=list(c(0.2,-0.4),c(0),c(0)),
                   g2_g2=list(c(0.3,-0.4),c(0),c(0))
)


set.seed(44)
seeds <- sample(1:500,50,replace=FALSE)
#SBM parameters
probmat <- matrix(c(.6,.3,.2,.7),nrow = 2)

group_num <- 2
alphaOrder_reg <- c(1,1,3,3,3)
betaOrder_reg <- list(c(1),c(2),rep(1,3),rep(2,3),c(2,0,0))
num_param <- c(2,3,6,9,5) * group_num**2

for (reg in 1:5){
  for (mod in c("grg","er","sbm")){
    print(c("regime ",reg, " model ",mod))
    cov_mat <- matrix(nrow=num_param[reg],ncol = 50)
    rmse_mat <- matrix(nrow=num_param[reg],ncol = 50)
    for (i in 1:50){
      print(c("seed ",i))
      set.seed(seeds[i])
      loc = list()
      if (mod=="sbm"){# SBM model
        net_sbm <- sample_sbm(20,probmat,c(10,10),directed = TRUE)
        V(net_sbm)$name <- as.character(seq(1,20,1)) # name nodes (characters)
        V(net_sbm)$group <- sample(c('g1','g2'), size=20, replace=TRUE)
        edgelist_sbm <- get.edgelist(net_sbm)
        nedges_sbm <- ecount(net_sbm)
        group_from <- V(net_sbm)$group[as.numeric(edgelist_sbm[,1])]
        group_to <- V(net_sbm)$group[as.numeric(edgelist_sbm[,2])]
        group_from_to <- paste(group_from, group_to, sep='_')
        if (reg %in% c(1,2)){
          assign(paste("alpha_par_",reg,"_sbm",sep=""),list(unlist(get(paste("alpha_par_",reg,sep=""))[group_from_to])))
        }else{
          lst <- get(paste('alpha_par_',reg,sep=''))[group_from_to]
          assign(paste("alpha_par_",reg,"_sbm",sep=""),list(unlist(lapply(lst, function(x) x[1])),
                                                            unlist(lapply(lst, function(x) x[2])),
                                                            unlist(lapply(lst, function(x) x[3]))))
        }
        grouplist_sbm <- split(seq_along(group_from_to),group_from_to)
      }
      
      if (mod=="er"){# ER model
        net_er <- erdos.renyi.game(20,p.or.m = 168,type = "gnm",directed = TRUE) # p=.4 fix the number of edges in graph ("gnm") rather than give prob of edge ("gnp")
        V(net_er)$name <- as.character(seq(1,20,1))# name nodes (characters)
        V(net_er)$group <- sample(c('g1','g2'), size=20, replace=TRUE)
        edgelist_er <- get.edgelist(net_er)
        nedges_er <- ecount(net_er)
        group_from <- V(net_er)$group[as.numeric(edgelist_er[,1])]
        group_to <- V(net_er)$group[as.numeric(edgelist_er[,2])]
        group_from_to <- paste(group_from, group_to, sep='_')
        if (reg %in% c(1,2)){
          assign(paste("alpha_par_",reg,"_er",sep=""),list(unlist(get(paste("alpha_par_",reg,sep=""))[group_from_to])))
        }else{
          lst <- get(paste('alpha_par_',reg,sep=''))[group_from_to]
          assign(paste("alpha_par_",reg,"_er",sep=""),list(unlist(lapply(lst, function(x) x[1])),
                                                           unlist(lapply(lst, function(x) x[2])),
                                                           unlist(lapply(lst, function(x) x[3]))))
        }
        grouplist_er <- split(seq_along(group_from_to),group_from_to)
      }
      
      if (mod=="grg"){# Random Geometric Graph
        lpvs2 <- sample_sphere_surface(dim = 2, n = 20,radius = .7)#.35 for 0.1 density, .7 for density 0.4
        net_grg <- sample_dot_product(lpvs2,directed = TRUE)
        #net_grg <- sample_grg(20, 0.45)
        V(net_grg)$name <- as.character(seq(1,20,1))# name nodes (characters)
        V(net_grg)$group <- sample(c('g1','g2'), size=20, replace=TRUE)
        edgelist_grg <- get.edgelist(net_grg)
        nedges_grg <- ecount(net_grg)
        group_from <- V(net_grg)$group[as.numeric(edgelist_grg[,1])]
        group_to <- V(net_grg)$group[as.numeric(edgelist_grg[,2])]
        group_from_to <- paste(group_from, group_to, sep='_')
        if (reg %in% c(1,2)){
          assign(paste("alpha_par_",reg,"_grg",sep=""),list(unlist(get(paste("alpha_par_",reg,sep=""))[group_from_to])))
        }else{
          lst <- get(paste('alpha_par_',reg,sep=''))[group_from_to]
          assign(paste("alpha_par_",reg,"_grg",sep=""),list(unlist(lapply(lst, function(x) x[1])),
                                                            unlist(lapply(lst, function(x) x[2])),
                                                            unlist(lapply(lst, function(x) x[3]))))
        }
        grouplist_grg <- split(seq_along(group_from_to),group_from_to)
      }
      
      # for(group in names(group_list)) print(group)
      # lapply(1:2, function(x)paste('g',x,sep=''))
      # library(Matrix)
      # lapply(names(group_list),function(x)Matrix(0,10,10,doDiag = FALSE))
      # wei_stag <- lapply(1:length(group_list),function(x)Matrix(0,5,5,doDiag = FALSE))
      # names(wei_stag) <- names(group_list)
      # wei_stag
      # wei_mat <- lapply(1:max(4),function(x)wei_stag)
      # lengths(wei_mat)
      # fit$mod
      # length(fit$wei_mat[[1]][['g1_g1']])
      # rowSums(fit$wei_mat[[1]][['g1_g2']]) #[1:10,]
      # simdata
      # all.equal(mat,fit$wei_mat)
      # mat == fit$wei_mat
      # mat[[1]][['g1_g1']]
      # fit$wei_mat[[1]][['g1_g1']]
      # all.equal(mat[[1]][['g1_g1']], fit$wei_mat[[1]][['g1_g1']])
      # vec_A <- as.numeric(mat[[1]][['g1_g1']])
      # vec_B <- as.numeric(fit$wei_mat[[1]][['g1_g1']])
      # Now compare the raw numbers
      #all.equal(vec_A, vec_B)
      
      simdata <- gnar_edge_sim(n=200,net=get(paste("net",mod,sep = "_")),alphaParams=get(paste("alpha_par",reg,mod,sep = "_")),
                               betaParams=get(paste("beta_par",reg,sep = "_")),
                               sigma = 1, meann=0, nedges=get(paste("nedges",mod,sep = "_")),
                               data_edges = get(paste("edgelist_",mod,sep = "")), groups = get(paste('grouplist',mod,sep='_')))
      
      fit <- gnar_edge_fit(simdata,get(paste("edgelist_",mod,sep = "")),alphaOrder_reg[reg],
                           betaOrder_reg[[reg]],get(paste("net_",mod,sep = "")),lead_lag_mat = NULL,groups=get(paste('grouplist',mod,sep='_')),
                           globalalpha = TRUE,lead_lag_weights = FALSE)
      
      comp_truevec <- c()
      for (lag in 1:alphaOrder_reg[reg]){
        # comp_truevec <- c(comp_truevec,get(paste("alpha_par_",reg,sep = ""))[ord],get(paste("beta_par_",reg,sep = ""))[[ord]])
        for (g in names(group_list)) {
          comp_truevec <- c(comp_truevec,get(paste("alpha_par",reg,sep = "_"))[[g]][lag])
          #print(paste('alpha_par',lag,g,sep='_'))
        }
      }
      for (lag in 1:alphaOrder_reg[reg]){
        if (betaOrder_reg[[reg]][lag]>0){
          for (g in names(group_list)) {
            for (stag in 1:betaOrder_reg[[reg]][lag]){
              comp_truevec <- c(comp_truevec,get(paste("beta_par",reg,sep = "_"))[[g]][[lag]][stag])
              #print(paste('beta_par',lag,stag,g,sep='_'))
            }
          }
        }
      }
      
      for(tr in 1:length(comp_truevec)){
        cov_mat[tr,i] <- between(comp_truevec[tr],confint(fit$mod)[tr,1],confint(fit$mod)[tr,2])
      }
      rmse_mat[,i] <- (fit$mod$coefficients-comp_truevec)^2
      
    }
    rownames(cov_mat) <- rownames(confint(fit$mod,level=.95))
    rownames(rmse_mat) <- names(fit$mod$coefficients)
    assign(paste("coverage_reg_",reg,"_mod_",mod,sep = ""),apply(cov_mat,1,sum,na.rm = TRUE)/50) # count NA values as 0
    assign(paste("rmse_reg_",reg,"_mod_",mod,sep = ""),apply(rmse_mat,1,function(x) sqrt(mean(x))))
  }
}

###### simulation scenario relevant to real data ######

# ER model
set.seed(12)
net_er <- erdos.renyi.game(86,p.or.m = 6858,type = "gnm",directed = TRUE) 
# name nodes (characters)
V(net_er)$name <- as.character(seq(1,86,1))
V(net_er)$group <- sample(c('g1','g2'), size=86, replace=TRUE)
edgelist_er <- get.edgelist(net_er)
par(mar = c(1, 1, 1, 1))
plot(net_er,vertex.size=15,edge.color="slategrey", vertex.label.color="black",vertex.color="lightblue",
     vertex.label.font=2, edge.arrow.size=0.18,vertex.label.degree=3.7,cex.sub=15,vertex.label.cex=.8) 
nedges <- ecount(net_er)
group_from <- V(net_er)$group[as.numeric(edgelist_er[,1])]
group_to <- V(net_er)$group[as.numeric(edgelist_er[,2])]
group_from_to <- paste(group_from,group_to,sep='_')
grouplist_er <- split(seq_along(group_from_to),group_from_to)

# Global alpha
# Two cases for ER network: GNAR(4,[1,1,1,1]), GNAR(4,[2,2,2,2]) 

# Regimes : Parameter specification

# Regime 1: GNAR(4,[1,1,1,1])
alpha_par_1 <- list(g1_g1=c(-0.6,-0.4,-0.2,-0.1),g1_g2=c(-0.3,0.4,0.1,0.2),g2_g1=c(0.3,0.1,0.2,0.1),g2_g2=c(0.1,-0.2,0.2,0.2))
beta_par_1 <- list(g1_g1=list(c(0.2),c(0.1),c(0.3),c(0.05)),
                   g1_g2=list(c(0.1),c(0.2),(0.1),(-0.2)),
                   g2_g1=list(c(0.05),c(-0.05),c(0.1),c(-0.1)),
                   g2_g2=list(c(0.1),c(0.2),c(-0.2),c(-0.1)))

# Regime 2: GNAR(4,[2,2,2,2]) 
# alpha_par_2 <- c(-0.6,-0.4,-0.2,-0.1)
# beta_par_2 <- list(c(0.4,-0.4),c(0.3,-0.4),c(0.5,-0.3),c(0.05,-0.1))

alpha_par_2 <- list(
  g1_g1 = c(-0.3, -0.2, 0.1, -0.1), # -0.5
  g1_g2 = c(0.2, -0.1, 0.3, 0.1), # 0.5
  g2_g1 = c(-0.1, 0.4, -0.2, -0.1), # 0
  g2_g2 = c(0.1, -0.3, 0.1, 0.2) # 0.1
)

beta_par_2 <- list(
  g1_g1 = list(c(-0.1, -0.2), c(0.1, -0.1), c(-0.2, 0.1), c(-0.1, -0.1)), # -0.6
  g1_g2 = list(c(-0.1, 0.3), c(0.2, -0.1), c(0.1, 0.1), c(-0.2, 0.1)), # 0.4
  g2_g1 = list(c(-0.2, 0.2), c(0.1, -0.1), c(-0.1, -0.2), c(0.1, -0.1)), # 0.3
  g2_g2 = list(c(-0.2, 0.1), c(0.1, 0.2), c(-0.1, -0.1), c(0.1, 0.2)) # 0.3
)


for (i in 1:2) {
  lst <- get(paste('alpha_par_',i,sep=''))[group_from_to]
  assign(paste("alpha_par_",i,"_real",sep=""),list(unlist(lapply(lst, function(x) x[1])),
                                                   unlist(lapply(lst, function(x) x[2])),
                                                   unlist(lapply(lst, function(x) x[3])),
                                                   unlist(lapply(lst, function(x) x[4]))))  
}

# Simulate data for each regime
seed_data <- 100
for (i in 1:2){
  set.seed(seed_data)
  assign(paste("gnar_edge_sim_data_",i,"_er_real",sep = ""),gnar_edge_sim(n=90,net = net_er,alphaParams = get(paste("alpha_par_",i,'_real',sep = "")),
                                                                          betaParams = get(paste("beta_par_",i,sep = "")),sigma = 1, meann=0,
                                                                          nedges=nedges,data_edges=edgelist_er, groups=grouplist_er))
  
  seed_data <- seed_data+10
}


# Fit the model on Simulated data
alphaOrder_reg <- c(4,4)
betaOrder_reg <- list(rep(1,4),rep(2,4))
for (i in 1:2){
  assign(paste("fit_er_",i,sep = ""),gnar_edge_fit(get(paste("gnar_edge_sim_data_",i,"_er_real",sep = "")),edgelist_er,alphaOrder_reg[i],
                                                   betaOrder_reg[[i]],net_er,lead_lag_mat = NULL,globalalpha = TRUE,lead_lag_weights = FALSE,groups=grouplist_er))
  
}

confint(fit_er_1$mod,  level=0.95)
confint(fit_er_2$mod,  level=0.95)
