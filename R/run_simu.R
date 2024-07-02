##################################################
# EFHimpacts

#Script: run EFH degrataion/restoration scenarios for a specified species 

# Juliette Champagnat
#################################################

# load packages
library(tidyverse);library(reshape2)
library(rfishbase);library(FishLife)

# 1. Load population parameters ----------------

sp_files <- "Herring" #"Plaice"  #"Cod"# 

if(sp_files=="Herring"){
  genus_name <- "Clupea" #
  sp_name <- "harengus"#
  fec_study = 10}
if(sp_files =="Cod"){
  genus_name <- "Gadus" #
  sp_name <- "morhua" #
  fec_study = 9}
if(sp_files == "Plaice"){
  genus_name <- "Pleuronectes"# 
  sp_name <- "platessa"}

source("R/chap3_git/get_archetype_data.R")
years <- c(1:200);Fval <- Selec*0.2
source("R/chap3_git/make_data.R",local=T)
source("R/chap3_git/Wbar_function.R",local=T)
Wbar_init <- mydata$Wbar


# 2. Run baseline MSY simulations ----
source("R/chap3_git/MSY_wraper.R")

purrr::map_dbl(seq(0,1,0.01),MSY_wraper,value_of_interest = "both" ,
               n_traj=1,years=c(1:500),path_simulation="Simulation/chap3_git",
               simulation.name='tre',
               # comments=paste0("K_",K,"_h_",as.numeric(h)))
)


# 3. Simulate EFH impacts -------
source("R/chap3_git/mult_fun.R")
source("R/chap3_git/Wbar_change_function.R")

path_simulation=paste0("Simulation/chap3_git/",sp_files); #where the output will be save

### 3.1 nurseries impacts ----------

# reset all population parameters
source("R/chap3_git/get_archetype_data.R")
source("R/chap3_git/make_data.R",local=T)
source("R/chap3_git/Wbar_function.R",local=T)
Wbar_init <- mydata$Wbar


# set & apply multiplier 
seq_explo <- c(0.8,1.2)#seq(0.9,1.1,0.1)

marginal_nurs <- mult_fun_surv(surv_wDI = rep(1,length(seq_explo)),
                          surv_rDI= seq_explo,
                          lambda_rDD_type = "rDI" ,
                          lambda_surf = rep(1,length(seq_explo)),
                          h=as.numeric(Predict[[1]]$Mean_pred[13]),
                          Wbar=mydata$Wbar,K=10^6,Mlarvae = Mlarvae,d_l=1,d_r=1) 

# Check viability of pop for all scenarios
marginal_nurs %>% 
  mutate(Mlarvae=Mlarvae,deltaT=1) %>% 
  mutate(Wbar=Wbar_init) %>% 
  mutate(surv_max = 4*h/(Wbar*(1-h))) %>% 
  mutate(alpha = 4*h/(exp(-Mlarvae)*Wbar*(1-h))) %>% 
  mutate(M_DI = -log(alpha)/deltaT) %>% 
  mutate(M_DD = M_DI/(K*(exp(M_DI*deltaT)-1))) %>% summary()
# nothing should be negative

# run simulations
for(i in 1:nrow(marginal_nurs)){
  sc_name_temp=paste0("wDI_",marginal_nurs$surv_wDI[i],"_rDI_",
                      marginal_nurs$surv_rDI[i],"_rDD_",
                      marginal_nurs$surv_rDI[i],
                      "_surf_",marginal_nurs$lambda_surf[i]) #simulation name
  
  
  # run the analysis
  source("R/chap3_git/MSY_wraper.R")
  Fval=0;DYN=F #not using the dynamic options for MSY search
  h <- marginal_nurs$h[i];K <- marginal_nurs$K[i]
  source("R/chap3_git/make_data.R")
  
  purrr::map_dbl(seq(0,1.5,0.005),MSY_wraper,value_of_interest = "both" ,age.plus=F,
                 n_traj=1,years=c(1:500),
                 path_simulation=path_simulation,
                 simulation.name=sc_name_temp,comments=paste0("K_",K,"_h_",h))
  
  
  # seq(0.005,1.495,0.01)
}

### 3.2 spawning ground impacts - eggs hypothesis -----

# reset all population parameters
source("R/chap3_git/get_archetype_data.R")
source("R/chap3_git/make_data.R",local=T)
source("R/chap3_git/Wbar_function.R",local=T)
Wbar_init <- mydata$Wbar

# set & apply multiplier
seq_explo <- c(0.8,1.2)#seq(0.9,1.1,0.1)

marginal_eggs <- mult_fun_surv(surv_wDI = seq_explo,
                          surv_rDI = rep(1,length(seq_explo)),
                          lambda_rDD_type = "rDI" ,
                          lambda_surf = rep(1,length(seq_explo)),
                          h=as.numeric(Predict[[1]]$Mean_pred[13]),
                          Wbar=mydata$Wbar,K=10^6,Mlarvae = Mlarvae,d_l=1,d_r=1) 

# Check viability of pop for all scenarios
marginal_eggs %>% 
  mutate(Mlarvae=Mlarvae,deltaT=1) %>% 
  mutate(Wbar=Wbar_init) %>% 
  mutate(surv_max = 4*h/(Wbar*(1-h))) %>% 
  mutate(alpha = 4*h/(exp(-Mlarvae)*Wbar*(1-h))) %>% 
  mutate(M_DI = -log(alpha)/deltaT) %>% 
  mutate(M_DD = M_DI/(K*(exp(M_DI*deltaT)-1))) %>% summary()
# nothing should be negative

# run simulations
for(i in 1:nrow(marginal_eggs)){
  sc_name_temp=paste0("wDI_",marginal_eggs$surv_wDI[i],
                      "_rDI_",marginal_eggs$surv_rDI[i],
                      "_rDD_",marginal_eggs$surv_rDI[i],
                      "_surf_",marginal_eggs$lambda_surf[i]) #simulation name
  
  
  # run the analysis
  source("R/chap3_git/MSY_wraper.R")
  Fval=0;DYN=F #not using the dynamic options for MSY search
  h <- marginal_eggs$h[i];K <- marginal_eggs$K[i]
  source("R/chap3_git/make_data.R")
  
  purrr::map_dbl(seq(0,1.5,0.005),MSY_wraper,value_of_interest = "both" ,age.plus=F,
                 n_traj=1,years=c(1:500),
                 path_simulation=path_simulation,
                 simulation.name=sc_name_temp,comments=paste0("K_",K,"_h_",h))
  
  
  # seq(0.005,1.495,0.01)
}


### 3.3 spawning ground impacts - spawners hypothesis -----

# reset all population parameters
source("R/chap3_git/get_archetype_data.R")
source("R/chap3_git/make_data.R",local=T)
source("R/chap3_git/Wbar_function.R",local=T)
Wbar_init <- mydata$Wbar

# set & apply multipliers
M_spwgrd_val <- c(0.8,1.2)
M_spwgrd_val <- M_spwgrd_val[which(M_spwgrd_val<exp(M_at_age[1]))] #having something higher than baseline M cannot work

# run simulations
for(i in seq_along(M_spwgrd_val)){
  sc_name_temp=paste0("Mspwgrd_surv",M_spwgrd_val[i]) #simulation name
  
  
  # run the analysis
  source("R/chap3_git/MSY_wraper.R")
  Fval=0;DYN=F #not using the dynamic options for MSY search
  M_spwgrd_surv <- M_spwgrd_val[i]
  source("R/chap3_git/make_data.R")
  source("R/chap3_git/Wbar_function.R")
  
  
  h <- as.numeric(Wbar_change(Wbar_new = mydata$Wbar ,Wbar_0=Wbar_init, h=Predict[[1]]$Mean_pred[13],K=10^6,Mlarvae = Mlarvae,d_l=1,d_r=1))
  
  source("R/chap3_git/make_data.R")
  source("R/chap3_git/Wbar_function.R")
  # mydata$Wbar <- Wbar_init
  
  purrr::map_dbl(seq(0,1.5,0.005),MSY_wraper,value_of_interest = "both" ,age.plus=F,
                 n_traj=1,years=c(1:500),
                 path_simulation=path_simulation,
                 simulation.name=sc_name_temp,
                 comments=paste0("K_",K,"_h_",h,"_Wbar_",mydata$Wbar))
}



### 3.4 Combined nurs x eggs ------------

# reset all population parameters
source("R/chap3_git/get_archetype_data.R")
source("R/chap3_git/make_data.R",local=T)
source("R/chap3_git/Wbar_function.R",local=T)
Wbar_init <- mydata$Wbar


# set & apply multipliers
mult_cross <- data.frame(surv_wDI=c(0.9,1.1),surv_rDI=c(0.9,1.1))#data.frame(surv_wDI=c(0.9,1,1.1),surv_rDI=c(0.9,1,1.1))
combined <- mult_fun_surv(surv_wDI = mult_cross$surv_wDI,
                          surv_rDI= mult_cross$surv_rDI,
                          lambda_rDD_type = "rDI" ,
                          lambda_surf = rep(1,nrow(mult_cross)),
                          h=as.numeric(Predict[[1]]$Mean_pred[13]),
                          Wbar=mydata$Wbar,K=10^6,Mlarvae =Mlarvae,d_l=1,d_r=1) 


# Check viability of pop for all scenarios
combined %>% 
  mutate(Mlarvae=Mlarvae,deltaT=1) %>% 
  mutate(Wbar=Wbar_init) %>% 
  mutate(surv_max = 4*h/(Wbar*(1-h))) %>% 
  mutate(alpha = 4*h/(exp(-Mlarvae)*Wbar*(1-h))) %>% 
  mutate(M_DI = -log(alpha)/deltaT) %>% 
  mutate(M_DD = M_DI/(K*(exp(M_DI*deltaT)-1))) %>% summary()
# nothing should be negative

# run simulations
for(i in 1:nrow(combined)){
  sc_name_temp=paste0("wDI_",combined$surv_wDI[i],
                      "_rDI_",combined$surv_rDI[i],
                      "_rDD_",combined$surv_rDI[i],
                      "_surf_",combined$lambda_surf[i]) #simulation name
  
  
  # run the analysis
  source("R/chap3_git/MSY_wraper.R")
  Fval=0;DYN=F #not using the dynamic options for MSY search
  h <- combined$h[i];K <- combined$K[i]
  source("R/chap3_git/make_data.R")
  
  purrr::map_dbl(seq(0,1.5,0.005),MSY_wraper,value_of_interest = "both" ,age.plus=F,
                 n_traj=1,years=c(1:500),
                 path_simulation=path_simulation,
                 simulation.name=sc_name_temp,
                 comments=paste0("K_",K,"_h_",h,"_Wbar_",mydata$Wbar))
  
}


### 3.5 Combined nurs x spawners ----------

# reset all population parameters
source("R/chap3_git/get_archetype_data.R")
source("R/chap3_git/make_data.R",local=T)
source("R/chap3_git/Wbar_function.R",local=T)
Wbar_init <- mydata$Wbar

# set & apply multiplier
mult_cross_sp_nurs <- data.frame(surv_rDI=c(0.9,1.1),M_spwgrd_val=c(0.9,1.1))


combined_sp_nurs <- mult_fun_surv(surv_wDI = rep(1,nrow(mult_cross_sp_nurs)),
                                  surv_rDI= mult_cross_sp_nurs$surv_rDI,
                                  lambda_rDD_type = "rDI" ,
                                  lambda_surf = rep(1,nrow(mult_cross_sp_nurs)),
                                  h=as.numeric(Predict[[1]]$Mean_pred[13]),
                                  Wbar=mydata$Wbar,K=10^6,Mlarvae =Mlarvae,d_l=1,d_r=1) %>% 
  mutate(Mspw=mult_cross_sp_nurs$M_spwgrd_val)

# run simulations
for(i in 1:nrow(combined_sp_nurs)){
  sc_name_temp=paste0("wDI_",combined_sp_nurs$surv_wDI[i],
                      "_rDI_",combined_sp_nurs$surv_rDI[i],
                      "_rDD_",combined_sp_nurs$surv_rDI[i],
                      "_surf_",combined_sp_nurs$lambda_surf[i]
                      ,"_Mspw_surv_",combined_sp_nurs$Mspw[i]) #simulation name
  
  
  Fval=0;DYN=F #not using the dynamic options for MSY search
  M_spwgrd_surv <- combined_sp_nurs$Mspw[i]
  
  #new Wbar value
  source("R/chap3_git/make_data.R")
  source("R/chap3_git/Wbar_function.R") 
  
  # update h value
  h <- as.numeric(Wbar_change(Wbar_new = mydata$Wbar ,Wbar_0=Wbar_init, h=combined_sp_nurs$h[i],
                              K=combined_sp_nurs$K,Mlarvae = Mlarvae,d_l=1,d_r=1))
  K <- combined_sp_nurs$K[i]
  
  #import those new values in mydata
  source("R/chap3_git/make_data.R")
  
  # add Wbar to mydata (it has been reset by make data)
  source("R/chap3_git/Wbar_function.R")
  
  source("R/chap3_git/MSY_wraper.R")
  purrr::map_dbl(seq(0,1.5,0.005),MSY_wraper,value_of_interest = "both" ,age.plus=F,
                 n_traj=1,years=c(1:500),
                 path_simulation=path_simulation,
                 simulation.name=sc_name_temp,
                 comments=paste0("K_",K,"_h_",h,"_Wbar_",mydata$Wbar))

}
