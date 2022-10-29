###################
#### Load libraries
###################
library(spdep)
library(sp)
library(MASS)
library(spam)
library(CARBayes)
library(coda)
library(truncnorm)
library(ggplot2)
library(rstan)

set.seed(128)

####################
#### Set up the grid
####################
#### Set up a square lattice region
x.easting <- 1:100
x.northing <- 1:100
Grid <- as.matrix(expand.grid(x.easting, x.northing))
K.grid <- nrow(Grid)

## Number of areal units K
K <- 150 # 100 or 500 or 1000

###########################################
#### Aggregate the grid data to areal units
###########################################
#### Choose the areal unit centroids
area.centroids <- Grid[sample(x=1:K.grid, size=K, replace=FALSE), ]

#### Create the distance matrix between each centroid and each grid square
Dist.mat <- array(NA, c(K, K.grid))
for(i in 1:K)
{
  centroid.temp <- area.centroids[i, ]
  Dist.mat[i, ] <- as.numeric(sqrt((centroid.temp[1] - Grid[ ,1])^2 + (centroid.temp[2] - Grid[ ,2])^2))
}

sp.grid <- SpatialPixelsDataFrame(points=Grid, data=data.frame(region=1:K.grid))

#### Assign each grid square to its nearest centroid
ordering <- apply(Dist.mat, 2, order)
area <- ordering[1, ]
sp.grid@data$area <- area
# spplot(sp.grid, "area")

#### Create the W matrix
W.area.nb.temp <- knn2nb(knearneigh(area.centroids, k=4), row.names=1:K)
W.area.nb <- make.sym.nb(W.area.nb.temp)
W.area <- nb2mat(W.area.nb, style = "B")
W.area.list <- nb2listw(W.area.nb)

####################################################################
### the distance matrix for area.centroids
D <- as.matrix(dist(area.centroids))

n_simulation <- 50


# for stan car model
W.mat<-W.area
W_rowsum <- apply(W.mat,1,sum)

#

N_edges = sum(W.mat)/ 2
node1 =NULL
node2 = NULL
for(i in 1:K)
{
  stuff <- which(W.mat[i,i:ncol(W.mat)]==1)+i-1
  node1 <- c(node1, rep(i, length(stuff)))
  node2 <- c(node2, stuff)
}

#####################################################################
# the above has done the spatial area set up
#####################################################################

fit.code_IAR <- stanc(file = "IAR_conf.stan") # convert to C++ code
fit.model_IAR <- stan_model(stanc_ret = fit.code_IAR) # compile C++ code
# save(fit.model_IAR, file = "fit.model_IAR_conf.RData")
fit.code_IAR_priorMethod <- stanc(file = "priorMethod_IAR_conf.stan") # convert to C++ code
fit.model_IAR_priorMethod <- stan_model(stanc_ret = fit.code_IAR_priorMethod) # compile C++ code
# save(fit.model_IAR_priorMethod, file = "fit.model_IAR_priorMethod_conf.RData")
fit.code_IAR_singleMethod <- stanc(file = "singleMethod_IAR_conf.stan") # convert to C++ code
fit.model_IAR_singleMethod <- stan_model(stanc_ret = fit.code_IAR_singleMethod) # compile C++ code
# save(fit.model_IAR_singleMethod, file = "fit.model_IAR_singleMethod_conf.RData")

fit.code_BYM <- stanc(file = "BYM_conf.stan") # convert to C++ code
fit.model_BYM <- stan_model(stanc_ret = fit.code_BYM) # compile C++ code
# save(fit.model_BYM, file = "fit.model_BYM_conf.RData")
fit.code_BYM_priorMethod <- stanc(file = "priorMethod_BYM_conf.stan") # convert to C++ code
fit.model_BYM_priorMethod <- stan_model(stanc_ret = fit.code_BYM_priorMethod) # compile C++ code
# save(fit.model_BYM_priorMethod, file = "fit.model_BYM_priorMethod_conf.RData")
fit.code_BYM_singleMethod <- stanc(file = "singleMethod_BYM_conf.stan") # convert to C++ code
fit.model_BYM_singleMethod <- stan_model(stanc_ret = fit.code_BYM_singleMethod) # compile C++ code
# save(fit.model_BYM_singleMethod, file = "fit.model_BYM_singleMethod_conf.RData")


for(phi.sd in c(0.04,0.2)) #random effects sd (phi)
{
  for(x.sd in c(1,5)) #sd of the covariate used to generated pollution data (\sigma_z)
  {
    for(xcovariate in c("spatial","nonspatial"))# whether spatial or not for the covariate used to generated pollution data (\sigma_z)
    {
      for(pollution.var in c(1,5))# \sigmasq_x
      {
        ###########################################################
        ######################### IAR model #######################
        set.seed(128)
        for(ii in 1:n_simulation)
        {
          # random effect on each area
          #### Create the quantities of the simulation
          ## Amount of spatial correlation
          nu1 <- 0.5 # Covariate x
          nu2 <- 0.1 # Random effects phi
          
          #### Compute the spatial variance matrices and their cholesky factors
          Sigma1 <- exp(-nu1 * D)
          Sigma2 <- exp(-nu2 * D)
          chol.Sigma1 <-t(chol(Sigma1))
          chol.Sigma2 <-t(chol(Sigma2))
          
          #### Random effects
          phi.area <- chol.Sigma2 %*% rnorm(K)
          phi.area <- phi.sd * (phi.area - mean(phi.area)) / sd(phi.area)
          # ggplot(mapping = aes(x=area.centroids[,1],y=area.centroids[,2], size=phi.area))+geom_point(alpha=0.5)
          
          ###
          #### Covariates for air pollution model
          ## Use the first two lines for an independent x covariate and the
          ## last two for a spatially correlated x covariate
          if (xcovariate=="nonspatial")
          {
            x.area <- rnorm(K)
            x.area <- x.sd * (x.area - mean(x.area)) / sd(x.area)
          }
          if (xcovariate=="spatial")
          {
            x.area <- chol.Sigma1 %*% rnorm(K)
            x.area <- x.sd * (x.area - mean(x.area)) / sd(x.area)
          }
          # ggplot(mapping = aes(x=area.centroids[,1],y=area.centroids[,2], size=x.area))+geom_point(alpha=0.5)
          
          
          ### generate air pollution for each area.
          
          beta_0 <- 0.8
          beta_1 <- 1.5
          sigmasq <- pollution.var # the bigger, the Rsq is higher
          
          x.pollution <- mvrnorm(1,mu = beta_0+beta_1*x.area,Sigma =sigmasq*diag(x=1,length(x.area),length(x.area)))
          # summary(lm(x.pollution~x.area))
          #########################################using fitted values to mimic the unknown x#########
          model <- lm(x.pollution~x.area)
          pred <- predict(model, newdata = data.frame(model$model[,-1]), se.fit = TRUE)
          x.pollution_set <- mvrnorm(n=10, mu=pred$fit, Sigma = diag(pred$se.fit^2))
          x.pollution_nouncertainty <- x.pollution_set[1,]
          ###########################################################
          
          ## add a confounder
          x.conf <- rnorm(length(x.area), mean = 0, sd=1)
          beta_conf <- 0.15
          
          ### generate expected counts
          
          #### Population and offset
          ## Population size P()
          pop.sd <- 10
          pop.min <- 30
          ## Disease prevalence
          e.prop <- 0.05 #
          pop.temp <- pop.sd * chol.Sigma2 %*% rnorm(K) + pop.min
          pop.grid <- pop.temp^2
          e.area <- pop.grid * e.prop
          
          ### generate Y
          ## Pollution - disease effect size
          beta <- 0.1
          theta <- exp(beta*x.pollution+phi.area +x.conf*beta_conf)
          Y.area <- rpois(n=length(theta),lambda = e.area*theta)
          
          # run the stan pollution model --------------------------------------------
          
          mod.data <- list(k_area=K,
                           Y_kt=Y.area,
                           E_kt=c(e.area),
                           area_AP=x.pollution,
                           area_conf=x.conf,
                           N_edges=N_edges,
                           node1=node1,
                           node2=node2
          )
          
          #### NO2, MEAN, EFFECT
          Nsample=5000
          thin <- 10
          ##IAR model
          stan_model <-  rstan::sampling(fit.model_IAR,  data = mod.data,seed=158, control = list(max_treedepth=15),iter = Nsample, chains = 1, thin = thin,verbose=TRUE, par=c("lambda","lambda_conf","nu2", "phi_kt"))
          stan_result <-rstan::extract(stan_model)
          save(stan_result, file=paste0("All_results/IAR","true_phi_",phi.sd,"sigmaZ",x.sd,xcovariate,"sqX",pollution.var,"_n_",ii,".RData"))
          
          ### no uncertainty model
          
          mod.data <- list(k_area=K,
                           Y_kt=Y.area,
                           E_kt=c(e.area),
                           area_AP=x.pollution_nouncertainty,
                           area_conf=x.conf,
                           N_edges=N_edges,
                           node1=node1,
                           node2=node2
          )
          
          #### NO2, MEAN, EFFECT
          Nsample=5000
          thin <- 10
          ##IAR model
          stan_model <-  rstan::sampling(fit.model_IAR,  data = mod.data,seed=158, control = list(max_treedepth=15),iter = Nsample, chains = 1, thin = thin,verbose=TRUE, par=c("lambda","lambda_conf","nu2", "phi_kt"))
          stan_result <-rstan::extract(stan_model)
          save(stan_result, file=paste0("All_results/IAR","noUncertainty_phi_",phi.sd,"sigmaZ",x.sd,xcovariate,"sqX",pollution.var,"_n_",ii,".RData"))    
          
          ### multiset model####
          for(jj in 1:10)
          {
            mod.data <- list(k_area=K,
                             Y_kt=Y.area,
                             E_kt=c(e.area),
                             area_AP=x.pollution_set[jj,],
                             area_conf=x.conf,
                             N_edges=N_edges,
                             node1=node1,
                             node2=node2
            )
            
            #### NO2, MEAN, EFFECT
            Nsample=5000
            thin <- 10
            ##IAR model
            stan_model <-  rstan::sampling(fit.model_IAR,  data = mod.data,seed=158, control = list(max_treedepth=15),iter = Nsample, chains = 1, thin = thin,verbose=TRUE, par=c("lambda","lambda_conf","nu2", "phi_kt"))
            stan_result <-rstan::extract(stan_model)
            save(stan_result, file=paste0("All_results/IAR","multipleSets_phi_",phi.sd,"sigmaZ",x.sd,xcovariate,"sqX",pollution.var,"_set_",jj,"_n_",ii,".RData"))   
          }
          
          ### prior method #######
          mod.data <- list(k_area=K,
                           Y_kt=Y.area,
                           E_kt=c(e.area),
                           n_AP_sample=10,
                           area_AP_sample=t(x.pollution_set),
                           area_conf=x.conf,
                           N_edges=N_edges,
                           node1=node1,
                           node2=node2
          )
          
          #### NO2, MEAN, EFFECT
          Nsample=5000
          thin <- 10
          ##IAR model
          stan_model <-  rstan::sampling(fit.model_IAR_priorMethod,  data = mod.data,seed=158, control = list(max_treedepth=15),iter = Nsample, chains = 1, thin = thin,verbose=TRUE, par=c("lambda","lambda_conf","nu2", "phi_kt"))
          stan_result <-rstan::extract(stan_model)
          save(stan_result, file=paste0("All_results/IAR","priorMethod_phi_",phi.sd,"sigmaZ",x.sd,xcovariate,"sqX",pollution.var,"_n_",ii,".RData"))   
          
          ### single method #######
          mod.data <- list(k_area=K,
                           Y_kt=Y.area,
                           E_kt=c(e.area),
                           area_AP_covariate=c(x.area),
                           area_conf=x.conf,
                           N_edges=N_edges,
                           node1=node1,
                           node2=node2
          )
          
          #### NO2, MEAN, EFFECT
          Nsample=5000
          thin <- 10
          ##IAR model
          stan_model <-  rstan::sampling(fit.model_IAR_singleMethod,  data = mod.data,seed=158, control = list(max_treedepth=15),iter = Nsample, chains = 1, thin = thin,verbose=TRUE, par=c("lambda","lambda_conf","nu2","area_AP","beta_0","beta_1","sigma","phi_kt"))
          stan_result <-rstan::extract(stan_model)
          save(stan_result, file=paste0("All_results/IAR","singleMethod_phi_",phi.sd,"sigmaZ",x.sd,xcovariate,"sqX",pollution.var,"_n_",ii,".RData"))  
          
          
          
          
          
          ###########################################################
          ######################### BYM model #######################
          # need to add bym component to generate Y
          ### generate Y
          ## Pollution - disease effect size
          ## bym component##########
          theta_in <- rnorm(n=K, mean = 0, sd=0.1)
          ##########################
          beta <- 0.1
          ### note that for the bym model, the disease counts need to re-produce####
          theta <- exp(beta*x.pollution+phi.area +x.conf*beta_conf+theta_in)
          Y.area <- rpois(n=length(theta),lambda = e.area*theta)
          ##########################################################################
          
          # run the stan pollution model --------------------------------------------
          
          mod.data <- list(k_area=K,
                           Y_kt=Y.area,
                           E_kt=c(e.area),
                           area_AP=x.pollution,
                           area_conf=x.conf,
                           N_edges=N_edges,
                           node1=node1,
                           node2=node2
          )
          stan_model <-  rstan::sampling(fit.model_BYM,  data = mod.data,seed=158, control = list(max_treedepth=15),
                                         iter = Nsample, chains = 1, thin = thin,verbose=TRUE, par=c("lambda","lambda_conf","nu2", "phi_kt","theta","tau_theta")
          )
          
          stan_result <-rstan::extract(stan_model)
          save(stan_result, file=paste0("All_results/BYM","true_phi_",phi.sd,"sigmaZ",x.sd,xcovariate,"sqX",pollution.var,"_n_",ii,".RData"))
          
          ##### no uncertainty bym ###
          mod.data <- list(k_area=K,
                           Y_kt=Y.area,
                           E_kt=c(e.area),
                           area_AP=x.pollution_nouncertainty,
                           area_conf=x.conf,
                           N_edges=N_edges,
                           node1=node1,
                           node2=node2
          )
          stan_model <-  rstan::sampling(fit.model_BYM,  data = mod.data,seed=158, control = list(max_treedepth=15),
                                         iter = Nsample, chains = 1, thin = thin,verbose=TRUE, par=c("lambda","lambda_conf","nu2", "phi_kt","theta","tau_theta")
          )
          
          stan_result <-rstan::extract(stan_model)
          save(stan_result, file=paste0("All_results/BYM","noUncertainty_phi_",phi.sd,"sigmaZ",x.sd,xcovariate,"sqX",pollution.var,"_n_",ii,".RData"))
          
          ### multiple set model ####
          
          for (jj in 1:10)
          {
            mod.data <- list(k_area=K,
                             Y_kt=Y.area,
                             E_kt=c(e.area),
                             area_AP=x.pollution_set[jj,],
                             area_conf=x.conf,
                             N_edges=N_edges,
                             node1=node1,
                             node2=node2
            )
            stan_model <-  rstan::sampling(fit.model_BYM,  data = mod.data,seed=158, control = list(max_treedepth=15),
                                           iter = Nsample, chains = 1, thin = thin,verbose=TRUE, par=c("lambda","lambda_conf","nu2", "phi_kt","theta","tau_theta")
            )
            
            stan_result <-rstan::extract(stan_model)
            save(stan_result, file=paste0("All_results/BYM","multipleSets_phi_",phi.sd,"sigmaZ",x.sd,xcovariate,"sqX",pollution.var,"_set_",jj,"_n_",ii,".RData"))
          }
          
          ##### prior Method bym ###
          mod.data <- list(k_area=K,
                           Y_kt=Y.area,
                           E_kt=c(e.area),
                           n_AP_sample=10,
                           area_AP_sample=t(x.pollution_set),
                           area_conf=x.conf,
                           N_edges=N_edges,
                           node1=node1,
                           node2=node2
          )
          stan_model <-  rstan::sampling(fit.model_BYM_priorMethod,  data = mod.data,seed=158, control = list(max_treedepth=15),
                                         iter = Nsample, chains = 1, thin = thin,verbose=TRUE, par=c("lambda","lambda_conf","nu2", "phi_kt","theta","tau_theta")
          )
          
          stan_result <-rstan::extract(stan_model)
          save(stan_result, file=paste0("All_results/BYM","priorMethod_phi_",phi.sd,"sigmaZ",x.sd,xcovariate,"sqX",pollution.var,"_n_",ii,".RData"))
          
          ##### single Method bym ###
          mod.data <- list(k_area=K,
                           Y_kt=Y.area,
                           E_kt=c(e.area),
                           area_AP_covariate=c(x.area),
                           area_conf=x.conf,
                           N_edges=N_edges,
                           node1=node1,
                           node2=node2
          )
          stan_model <-  rstan::sampling(fit.model_BYM_singleMethod,  data = mod.data,seed=158, control = list(max_treedepth=15),
                                         iter = Nsample, chains = 1, thin = thin,verbose=TRUE, par=c("lambda","lambda_conf","nu2","area_AP","beta_0","beta_1","sigma","phi_kt","theta","tau_theta")
          )
          
          stan_result <-rstan::extract(stan_model)
          save(stan_result, file=paste0("All_results/BYM","singleMethod_phi_",phi.sd,"sigmaZ",x.sd,xcovariate,"sqX",pollution.var,"_n_",ii,".RData"))
          
          
        }#simulation n end
        
        
      }#pollution.var end
      
    }#xcovariate end
    
  } #x.sd end
}#phi.sd end


