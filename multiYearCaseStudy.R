library(Epi)
library(sp)
library(spdep)
library(rstan)
# load the ready data 
load(file="data.combined_multiYear.RData")
data.combined <- data.combined_multiYear

##build neibough matrix
data.combined<-na.omit(data.combined)
W.nb<-poly2nb(data.combined,row.names=data.combined$FeatureCode)
##change into matrix form
W.mat<-nb2mat(W.nb,style="B", zero.policy=TRUE)

# for stan car model
W_rowsum <- apply(W.mat,1,sum)

N_edges = sum(W.mat)/ 2
node1 =NULL
node2 = NULL
for(i in 1:1279)
{
  stuff <- which(W.mat[i,i:ncol(W.mat)]==1)+i-1
  node1 <- c(node1, rep(i, length(stuff)))
  node2 <- c(node2, stuff)
}

######################BYM MODEL#############################

fit.code <- stanc(file = "D:/OneDrive - 汕头大学/HGW2021/Project/Shantou_university/STUpaper/inference_from_various_pollutiondata/BYM_caseStudy_multiTimes_4covariates.stan") # convert to C++ code
fit.model <- stan_model(stanc_ret = fit.code) # compile C++ code
# save(fit.model, file="D:/OneDrive - 汕头大学/HGW2021/Project/Shantou_university/STUpaper/inference_from_various_pollutiondata/BYM_caseStudy_multiTimes_4covariates.RData")



# run the stan pollution model --------------------------------------------

mod.data <- list(k_area=1279,
                 t_year=6,
                 n_covariate=4,
                 Y_kt=as.matrix(data.combined@data[,c("Deaths2014","Deaths2015","Deaths2016","Deaths2017","Deaths2018","Deaths2019")]),
                 E_kt=as.matrix(data.combined@data[,c("E_Deaths2014","E_Deaths2015","E_Deaths2016","E_Deaths2017","E_Deaths2018","E_Deaths2019")]),
                 area_AP_one=as.matrix(data.combined@data[,c("pm252013mean","pm252014mean","pm252015mean","pm252016mean","pm252017mean","pm252018mean")]),
                 covariate_one=as.matrix(data.combined@data[,c("Crime","Crime","Crime","Crime","Crime","Crime")]),
                 covariate_two=as.matrix(data.combined@data[,c("EST","EST","EST","EST","EST","EST")]),
                 covariate_three=as.matrix(data.combined@data[,c("AtS","AtS","AtS","AtS","AtS","AtS")]),
                 covariate_four=as.matrix(data.combined@data[,c("Housing","Housing","Housing","Housing","Housing","Housing")]),
                 N_edges=N_edges,
                 node1=node1,
                 node2=node2
)

Nsample=10000
thin=10


stan_model <-  rstan::sampling(fit.model,  data = mod.data,seed=158, control = list(max_treedepth=15),
                               iter = Nsample, chains = 2, thin = thin,verbose=TRUE, par=c("alpha","lambda1","nu2","phi_kt","theta","tau_theta"))
stan_result <-rstan::extract(stan_model)
save(stan_result, file=paste0("results/stan_result","BYM_caseStudy_multiYears_mean",".RData"))


## max 

# the Urban_rural column, urban is 1 and rural is 0
mod.data <- list(k_area=1279,
                 t_year=6,
                 n_covariate=4,
                 Y_kt=as.matrix(data.combined@data[,c("Deaths2014","Deaths2015","Deaths2016","Deaths2017","Deaths2018","Deaths2019")]),
                 E_kt=as.matrix(data.combined@data[,c("E_Deaths2014","E_Deaths2015","E_Deaths2016","E_Deaths2017","E_Deaths2018","E_Deaths2019")]),
                 area_AP_one=as.matrix(data.combined@data[,c("pm252013max","pm252014max","pm252015max","pm252016max","pm252017max","pm252018max")]),
                 covariate_one=as.matrix(data.combined@data[,c("Crime","Crime","Crime","Crime","Crime","Crime")]),
                 covariate_two=as.matrix(data.combined@data[,c("EST","EST","EST","EST","EST","EST")]),
                 covariate_three=as.matrix(data.combined@data[,c("AtS","AtS","AtS","AtS","AtS","AtS")]),
                 covariate_four=as.matrix(data.combined@data[,c("Housing","Housing","Housing","Housing","Housing","Housing")]),
                 N_edges=N_edges,
                 node1=node1,
                 node2=node2
)

Nsample=10000
thin=10


stan_model <-  rstan::sampling(fit.model,  data = mod.data,seed=158, control = list(max_treedepth=15),
                               iter = Nsample, chains = 2, thin = thin,verbose=TRUE, par=c("alpha","lambda1","nu2","phi_kt","theta","tau_theta"))
stan_result <-rstan::extract(stan_model)
save(stan_result, file=paste0("results/stan_result","BYM_caseStudy_multiYears_max",".RData"))


######################IAR MODEL#############################

fit.code <- stanc(file = "D:/OneDrive - 汕头大学/HGW2021/Project/Shantou_university/STUpaper/inference_from_various_pollutiondata/IAR_caseStudy_multiTimes_4covariates.stan") # convert to C++ code
fit.model <- stan_model(stanc_ret = fit.code) # compile C++ code
save(fit.model, file="D:/OneDrive - 汕头大学/HGW2021/Project/Shantou_university/STUpaper/inference_from_various_pollutiondata/IAR_caseStudy_multiTimes_4covariates.RData")



# run the stan pollution model --------------------------------------------

mod.data <- list(k_area=1279,
                 t_year=6,
                 n_covariate=4,
                 Y_kt=as.matrix(data.combined@data[,c("Deaths2014","Deaths2015","Deaths2016","Deaths2017","Deaths2018","Deaths2019")]),
                 E_kt=as.matrix(data.combined@data[,c("E_Deaths2014","E_Deaths2015","E_Deaths2016","E_Deaths2017","E_Deaths2018","E_Deaths2019")]),
                 area_AP_one=as.matrix(data.combined@data[,c("pm252013mean","pm252014mean","pm252015mean","pm252016mean","pm252017mean","pm252018mean")]),
                 covariate_one=as.matrix(data.combined@data[,c("Crime","Crime","Crime","Crime","Crime","Crime")]),
                 covariate_two=as.matrix(data.combined@data[,c("EST","EST","EST","EST","EST","EST")]),
                 covariate_three=as.matrix(data.combined@data[,c("AtS","AtS","AtS","AtS","AtS","AtS")]),
                 covariate_four=as.matrix(data.combined@data[,c("Housing","Housing","Housing","Housing","Housing","Housing")]),
                 N_edges=N_edges,
                 node1=node1,
                 node2=node2
)

Nsample=10000
thin=10


stan_model <-  rstan::sampling(fit.model,  data = mod.data,seed=158, control = list(max_treedepth=15),
                               iter = Nsample, chains = 2, thin = thin,verbose=TRUE, par=c("alpha","lambda1","nu2","phi_kt"))
stan_result <-rstan::extract(stan_model)
save(stan_result, file=paste0("results/stan_result","IAR_caseStudy_multiYears_mean",".RData"))


## max 

mod.data <- list(k_area=1279,
                 t_year=6,
                 n_covariate=4,
                 Y_kt=as.matrix(data.combined@data[,c("Deaths2014","Deaths2015","Deaths2016","Deaths2017","Deaths2018","Deaths2019")]),
                 E_kt=as.matrix(data.combined@data[,c("E_Deaths2014","E_Deaths2015","E_Deaths2016","E_Deaths2017","E_Deaths2018","E_Deaths2019")]),
                 area_AP_one=as.matrix(data.combined@data[,c("pm252013max","pm252014max","pm252015max","pm252016max","pm252017max","pm252018max")]),
                 covariate_one=as.matrix(data.combined@data[,c("Crime","Crime","Crime","Crime","Crime","Crime")]),
                 covariate_two=as.matrix(data.combined@data[,c("EST","EST","EST","EST","EST","EST")]),
                 covariate_three=as.matrix(data.combined@data[,c("AtS","AtS","AtS","AtS","AtS","AtS")]),
                 covariate_four=as.matrix(data.combined@data[,c("Housing","Housing","Housing","Housing","Housing","Housing")]),
                 N_edges=N_edges,
                 node1=node1,
                 node2=node2
)

Nsample=10000
thin=10


stan_model <-  rstan::sampling(fit.model,  data = mod.data,seed=158, control = list(max_treedepth=15),
                               iter = Nsample, chains = 2, thin = thin,verbose=TRUE, par=c("alpha","lambda1","nu2","phi_kt"))
stan_result <-rstan::extract(stan_model)
save(stan_result, file=paste0("results/stan_result","IAR_caseStudy_multiYears_max",".RData"))





