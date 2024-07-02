###########################################################
# EFHimpacts

# Script: for a given species name grab related data needed for simulation

# Juliette C
###########################################################

DYN=F
#########
# vignette("tutorial","FishLife")

#grab from rfishbase
df_estimates <- estimate(paste0(genus_name," ",sp_name))
df_fecundity <- fecundity(paste0(genus_name," ",sp_name))

#grab from Fishlife
Predict = Plot_taxa( Search_species(Genus=genus_name,Species=sp_name)$match_taxonomy, mfrow=c(2,3),Database=FishLife::FishBase_and_RAM)

#####Age related#####
age_max <- as.numeric(round(exp(Predict[[1]]$Mean_pred[4])))
age_rec <- 1 #age de recrutement ie 1er qui n'est pas englobé dans le h
ages <-  c(age_rec:age_max) # is age_max the age we want to use here? 
age_sa <- age_rec # how to find it?
age_mat <- as.numeric(round(exp(Predict[[1]]$Mean_pred[5]))) #
age_fullselec <- age_mat #how to proceed ? same for all fish? 

#######Spawning related##########

if(genus_name=="Clupea"&sp_name=="harengus"){
  d_sp <- 285 #data of spawning peak 
  ts_sp <- 30#time span of spawning season
}else if(genus_name=="Gadus"&&sp_name=="morhua"){
  d_sp <- 29 #data of spawning peak 
  ts_sp <- 30#time span of spawning season
}else if(genus_name=="Pleuronectes"&&sp_name=="platessa"){
  d_sp <- 74 #data of spawning peak 
  ts_sp <- 30 #time span of spawning season
}else{
  d_sp <- 100 #data of spawning peak 
  ts_sp <- 30 #time span of spawning season
}


######Nurseries related##########
d_sa <- 1 #date (in day) of nurseries leaving for subadult fish 

######Zones###########
spwn_grd <- 1
spwn_grd_names <- 'spw'

nurseries <-1
nurseries_names <- 'nurs'

regions <-1
regions_names <- "region"

#########Female proportion######
Pf <- rep(0.5,length(ages)) #generic value

########## Spawning ground related ###########
# Density dependence in F_spw_grd # 0 = no DD effects
beta_DD <- 0
#surfaces of spawning grounds #with 1 nothing happens
Surf_spwn_grd <- rep(1,length(spwn_grd))
M_spwgrd_mult <- 1 #multiplier of natural mortality on spawning ground
M_spwgrd_surv <- 1 #multiplier of natural mortality on spawning ground SURVIVAL TYPE

##########Process errors########
sd_Sad <- 0
sd_Sjuv <- 0
sd_R <- 0


###### Lenght & Weight at age######
# parameter of Von B relationship 

# for length L(a)=L_inf * (1-exp(-K_growth*a))
K_growth = as.numeric(exp(Predict[[1]]$Mean_pred[2]))
L_inf = as.numeric(exp(Predict[[1]]$Mean_pred[1]))
LatA <- L_inf * (1-exp(-K_growth*(ages-(-0.1))))

# for weight W(a) = w_cst * L(a) ^ w_pow
if(!is.na(df_estimates$a)){w_cst <- df_estimates$a}
if(!is.na(df_estimates$b)){w_pow <- df_estimates$b}
if(!any(!is.na(df_estimates$a))&!any(!is.na(df_estimates$b))){print("You need some length-weight parameters")}

## Sometimes the weight constant has to be rescale by 10^3 sometime no
SWatA <- CWatA <- (w_cst/1000) * LatA ^ w_pow

#I'm using a twisted way to determine is the order of magnitude of Winfinity (from Thorson)
# and the weight of max age fish match

if(round(log10(as.numeric(exp(Predict[[1]]$Mean_pred[3]))))!=round(log10(CWatA[age_max]))){
  SWatA <- CWatA <- (w_cst) * LatA ^ w_pow
}
# if not i'm writing it for user to check
if(round(log10(as.numeric(exp(Predict[[1]]$Mean_pred[3]))))!=round(log10(CWatA[age_max]))){
  print("Seems to have a weight issue")
}

CWatAval <- SWatAval <- CWatA
######## Maturity & Fecondity#############

## maturity 
# from a knife edge function
mat <- ifelse(ages<age_mat,0,1)

## fecundity

# from a fonction of length fec(a) = fec_cst * L(a) ^ fec_pow (estimates from rfishbase)

if(any(!is.na(df_fecundity$a))&any(!is.na(df_fecundity$b))){ # there are fec-L estimates in fish base
  
  if(sum(!is.na(df_fecundity$a))>1|sum(!is.na(df_fecundity$b))>1) #there are more than one L-fec estimate
    
    if(is.null(fec_study)){ #need to choose which fec study will be used
    print("Decide which fec study you want to use and enter fec_study number (line of the study in the df_fecundity)")
      print(df_fecundity %>% select(a,b,r2,Number))
    
    }else{ #take a and b estimates from this study
      fec_cst <- df_fecundity[fec_study,]$a
      fec_pow <-  df_fecundity[fec_study,]$b
    }
}

#for some species i'm gathering fec information from litterature
if(genus_name=="Pleuronectes"&&sp_name=="platessa"){
  fec_cst <- 2.33;fec_pow<- 3.10

}

fec <- ifelse(ages<age_mat,0,fec_cst*(LatA[ages])^fec_pow)

######Eggs survival#######
p_eggs_surv <- rep(1, length(spwn_grd))


#####Selectivity###### 
# for now using a knife edge starting at maturity (c'est aussi ce que fait Kindsvater)
Selec <- ifelse(ages<age_mat,0,1)

########Natural mortality############
M_at_age <- rep(exp(Predict[[1]]$Mean_pred[6]),length(ages))

if(genus_name %in% c("Gadus") &&sp_name%in%c("morhua")){
  Mlarvae <- 5*log(10)
}else {Mlarvae <- 3*log(10)}

###### Nurseries parameters######
alpha <- NULL
Surf <- 1
K <- 10^6
h <- as.numeric(Predict[[1]]$Mean_pred[13])

#######Dispersion matrices############
##useless for now because no multiples zones
D_larvae <- diag(length(spwn_grd)) #Larval key
D_subadult <- diag(length(nurseries)) # Sub adult matrix
D_spawner <- diag(length(spwn_grd)) # Spawner matrix 
D_return <- diag(length(spwn_grd)) # Post spawning matrix 
D_adult <- diag(length(regions)) # Adult connectivity 

####### Inits ##############
# initial numbers as a decreasing function of K and natural mortality
#/!\ Eggs estimated from adult 
NatA_0 <- purrr::map_dbl(ages,function(x){K*exp(-sum(M_at_age[1:x]))})

Egg_0 <- sum(NatA_0*mat*fec*Pf,na.rm = T)


#######################################################################################################################


