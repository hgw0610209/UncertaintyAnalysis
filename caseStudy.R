
library("shapefiles")
library("sp")
library("tidyverse")
library("spdep") 
library("raster") 


load(file="data.combined.RData")

# Color map
# Layout of the map
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), 
                   offset = c(145000,570000), scale = 30000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), 
                 offset = c(145000, 560000), scale = 50000, fill = c("transparent", "black"))
text1 <- list("sp.text", c(145000, 550000), "0")
text2 <- list("sp.text", c(180000, 550000), "50 km")


spplot(data.combined, "pm25_2019_mean", sp.layout = list(northarrow, scalebar, text1, text2),
       at = seq(min(data.combined$pm25_2019_mean), max(data.combined$pm25_2019_mean), length.out = 8), col="transparent", 
       col.regions=hsv(0.65,seq(0.25,1,length.out=20),1),
       colorkey = list(labels = list( cex = 2)))




##build neibough matrix
# data.combined<-na.omit(data.combined)
W.nb<-poly2nb(data.combined,row.names=rownames(health_Data))
##change into matrix form
W.mat<-nb2mat(W.nb,style="B", zero.policy=TRUE)

# for stan car model
W_rowsum <- apply(W.mat,1,sum)

#

N_edges = sum(W.mat)/ 2
node1 =NULL
node2 = NULL
for(i in 1:1279)
{
  stuff <- which(W.mat[i,i:ncol(W.mat)]==1)+i-1
  node1 <- c(node1, rep(i, length(stuff)))
  node2 <- c(node2, stuff)
}




library(rstan)

fit.code_BYM_caseStudy <- stanc(file = "BYM_caseStudy.stan") # convert to C++ code
fit.model_BYM_caseStudy <- stan_model(stanc_ret = fit.code_BYM_caseStudy) # compile C++ code
# save(fit.model_BYM_caseStudy, file="fit.model_BYM_caseStudy.RData")
# 
fit.code_IAR_caseStudy <- stanc(file = "IAR_caseStudy.stan") # convert to C++ code
fit.model_IAR_caseStudy <- stan_model(stanc_ret = fit.code_IAR_caseStudy) # compile C++ code
# save(fit.model_IAR_caseStudy, file="fit.model_IAR_caseStudy.RData")

# load(file="fit.model_BYM_caseStudy.RData")
# load(file="fit.model_IAR_caseStudy.RData")


# run the stan pollution model --------------------------------------------

mod.data <- list(k_area=1279,
                 n_covariate=4,
                 Y_kt=data.combined$Value,
                 E_kt=data.combined$Expected_count,
                 area_AP=data.combined$pm25_2019_mean,
                 covariates=as.matrix(data.combined@data[,9:12]),
                 N_edges=N_edges,
                 node1=node1,
                 node2=node2
)

#### NO2, MEAN, EFFECT
Nsample=10000
thin=10

for(AP in c("pm25_2019_mean","pm25_2018_mean","pm25_2017_mean","pm25_2016_mean","pm25_2019_max","pm25_2018_max","pm25_2017_max","pm25_2016_max"))
{
  
  mod.data$area_AP <- data.combined[,AP]@data[,1]
  
  stan_model <-  rstan::sampling(fit.model_BYM_caseStudy,  data = mod.data,seed=158, control = list(max_treedepth=15),
                                 iter = Nsample, chains = 2, thin = thin,verbose=TRUE, par=c("alpha","lambda","nu2","phi_kt","theta","tau_theta")
  )
  stan_result <-rstan::extract(stan_model)
  save(stan_result, file=paste0("results/stan_result","BYM_caseStudy_",AP,".RData"))
  
  
  
  stan_model <-  rstan::sampling(fit.model_IAR_caseStudy,  data = mod.data,seed=158, control = list(max_treedepth=15),
                                 iter = Nsample, chains = 2, thin = thin,verbose=TRUE, par=c("alpha","lambda","nu2","phi_kt")
  )
  stan_result <-rstan::extract(stan_model)
  save(stan_result, file=paste0("results/stan_result","IAR_caseStudy_",AP,".RData"))
  
  
  
}

## calculate DIC, MSE, log-likelihood

resultList <- vector(mode = "list",length = 16)
i=1
for(pollution_type in c("pm25_2019_mean","pm25_2018_mean","pm25_2017_mean","pm25_2016_mean","pm25_2019_max","pm25_2018_max","pm25_2017_max","pm25_2016_max"))
{
 load(paste0("results/stan_resultBYM_caseStudy_",pollution_type,".RData"))
  resultList[[i]]= stan_result
i <- i+1
}

for(pollution_type in c("pm25_2019_mean","pm25_2018_mean","pm25_2017_mean","pm25_2016_mean","pm25_2019_max","pm25_2018_max","pm25_2017_max","pm25_2016_max"))
{
 load(paste0("results/stan_resultIAR_caseStudy_",pollution_type,".RData"))
  resultList[[i]]= stan_result
i <- i+1
}

## the counts data
load(file="data.combined.RData")

  #########################
    ## Calculate the deviance
    #########################
offset.mat <- log(data.combined$Expected_count)
covariates <- as.matrix(data.combined@data[,9:12])
exposures <- as.matrix(data.combined@data[,c("pm25_2019_mean","pm25_2018_mean","pm25_2017_mean","pm25_2016_mean","pm25_2019_max","pm25_2018_max","pm25_2017_max","pm25_2016_max")])

DIC <- NULL
likelihood <- NULL
MSE<- NULL

for(i in 1:8)
{
  
BYMmodel <- resultList[[i]]
regression.mat <- covariates%*%t(BYMmodel$alpha) # 1279 (areas) x 1000 (samples)
pollution.mat <- exposures[,i]%*%t(BYMmodel$lambda) # 1279 (areas) x 1000 (samples)
theta.mat <- t(BYMmodel$theta*c(1/sqrt(BYMmodel$tau_theta)))# 1279 (areas) x 1000 (samples)
phi.mat <- t(BYMmodel$phi_kt)# 1279 (areas) x 1000 (samples)

fitted <- as.numeric(exp(offset.mat + regression.mat+ pollution.mat +theta.mat +phi.mat))
deviance.all <- dpois(x=rep(data.combined$Value,1000), lambda=fitted, log=TRUE)
meanofLikelihood <- sum(deviance.all)/1000 # this contains 1000 samples   
    
regression.mat_mean <- covariates%*%(apply(BYMmodel$alpha,2,median)) # 1279 (areas) 
pollution.mat_mean <- exposures[,i]*(median(BYMmodel$lambda)) # 1279 (areas) 
theta.mat_mean <- apply((BYMmodel$theta)*c(1/sqrt(BYMmodel$tau_theta)),2,median)# 1279 (areas) 
phi.mat_mean <- apply(BYMmodel$phi_kt,2,median)# 1279 (areas) 


fitted_mean <- as.numeric(exp(offset.mat + regression.mat_mean+ pollution.mat_mean +theta.mat_mean +phi.mat_mean))
Likelihoodofmean <- sum(dpois(x=data.combined$Value, lambda=fitted_mean, log=TRUE))

likelihood <- c(likelihood, Likelihoodofmean)

DIC <- c(DIC, 2*Likelihoodofmean-4*meanofLikelihood)
MSE <- c(MSE, mean((data.combined$Value-fitted_mean)^2))

}


for(i in 9:16)
{
  
IARmodel <- resultList[[i]]
regression.mat <- covariates%*%t(IARmodel$alpha) # 1279 (areas) x 1000 (samples)
pollution.mat <- exposures[,i-8]%*%t(IARmodel$lambda) # 1279 (areas) x 1000 (samples)
phi.mat <- t(IARmodel$phi_kt)# 1279 (areas) x 1000 (samples)

fitted <- as.numeric(exp(offset.mat + regression.mat+pollution.mat +phi.mat))
deviance.all <- dpois(x=rep(data.combined$Value,1000), lambda=fitted, log=TRUE)
meanofLikelihood <- sum(deviance.all)/1000 # this contains 1000 samples   
    
regression.mat_mean <- covariates%*%(apply(IARmodel$alpha,2,median)) # 1279 (areas) 
pollution.mat_mean <- exposures[,i-8]*(median(IARmodel$lambda)) # 1279 (areas) 
phi.mat_mean <- apply(IARmodel$phi_kt,2,median)# 1279 (areas) 

fitted_mean <- as.numeric(exp(offset.mat + regression.mat_mean+pollution.mat_mean +phi.mat_mean))
Likelihoodofmean <- sum(dpois(x=data.combined$Value, lambda=fitted_mean, log=TRUE))

likelihood <- c(likelihood, Likelihoodofmean)

DIC <- c(DIC, 2*Likelihoodofmean-4*meanofLikelihood)

MSE <- c(MSE, mean((data.combined$Value-fitted_mean)^2))

}

############cumulative exposure model########################

data.combined$pm25meancumul_01 <- 0.5*(data.combined$pm25_2019_mean+data.combined$pm25_2018_mean)
data.combined$pm25meancumul_02 <- (data.combined$pm25_2019_mean+data.combined$pm25_2018_mean+data.combined$pm25_2017_mean)/3
data.combined$pm25meancumul_03 <- (data.combined$pm25_2019_mean+data.combined$pm25_2018_mean+data.combined$pm25_2017_mean+data.combined$pm25_2016_mean)/4

# run the stan pollution model --------------------------------------------

mod.data <- list(k_area=1279,
                 n_covariate=4,
                 Y_kt=data.combined$Value,
                 E_kt=data.combined$Expected_count,
                 area_AP=data.combined$pm25_2019_mean,
                 covariates=as.matrix(data.combined@data[,9:12]),
                 N_edges=N_edges,
                 node1=node1,
                 node2=node2
)

#### NO2, MEAN, EFFECT
Nsample=10000
thin=10

for(AP in c("pm25meancumul_01","pm25meancumul_02","pm25meancumul_03"))
{
  
  mod.data$area_AP <- data.combined[,AP]@data[,1]
  
  stan_model <-  rstan::sampling(fit.model_BYM_caseStudy,  data = mod.data,seed=158, control = list(max_treedepth=15),
                                 iter = Nsample, chains = 2, thin = thin,verbose=TRUE, par=c("alpha","lambda","nu2","phi_kt","theta","tau_theta")
  )
  stan_result <-rstan::extract(stan_model)
  save(stan_result, file=paste0("results/stan_result","BYM_caseStudy_",AP,".RData"))
  
  
  
  stan_model <-  rstan::sampling(fit.model_IAR_caseStudy,  data = mod.data,seed=158, control = list(max_treedepth=15),
                                 iter = Nsample, chains = 2, thin = thin,verbose=TRUE, par=c("alpha","lambda","nu2","phi_kt")
  )
  stan_result <-rstan::extract(stan_model)
  save(stan_result, file=paste0("results/stan_result","IAR_caseStudy_",AP,".RData"))
  
  
  
}
