##################################################
# EFHimpacts

#Script: for a given F combination (3 values), this script:
# simulates population dynamic for a given number of year 
# then computes # mean Catch and/or SSB for the last half of the simulation period
# saves values in a txt file


# Juliette Champagnat

##########################################

"MSY_wraper" <- function(F_fixed_value,simulation.name = 'test1',
                         what_save = 'stat',   #'all'
                         value_of_interest = "both",#'C',# "SSB" ,  
                         path_simulation=path_simulation,
                         n_traj=1,
                         years=c(1:100),
                         age.plus=F,
                         comments=NULL
                         
){
  #turning non equilibrium flags false
  Flag_Ceq<-Flag_SSBeq<-Flag_Req <- MSY_error <- NULL
  
  #depending if I'm using the script in a Rmd or in a normal script
  if(grepl(pattern="R/",getwd())){wdext <- "../"}else{wdext <- "R/"}
  
  
  Fval <- matrix(data=NA,nrow=length(ages),ncol=length(regions))
  #create Fval from the wanted F combination and selectivity at age
    for (a in 1:length(ages)){
    #for (z in 1:length(regions)){ 
    Fval[a] <- as.numeric(Selec[a] * F_fixed_value)
  }
  
  
  
  #format and arrange data 
  source(paste0(wdext,"chap3_git/make_data.R"),local=T)
  source(paste0(wdext,"chap3_git/Wbar_function.R"),local=T)
  
  
  #### load simulation function #####
  source(paste0(wdext,"chap3_git/sim_function_new.R"),local=T)
  
  #### run simulation ######
  S1 <-  life_cycle_sim(data=mydata,traj=n_traj,param_R = "h",age.plus = age.plus)
  
  # check if equilibrium has been reach #/!\ works only for deterministic i think
  
  #add a condition for if the stock is collapse
  if(S1$C_tot[max(years)]>10^-1){
    if(!round(S1$C_tot[max(years)])==round(S1$C_tot[max(years)-1])){Flag_Ceq <- paste0("for F=",F_fixed_value, " C_tot is not at equilibrium");print(Flag_Ceq)}
    if(!round(S1$SSB[max(years)])==round(S1$SSB[max(years)-1])){Flag_SSBeq<-paste0("for F=",F_fixed_value, " SSB is not at equilibrium");print(Flag_SSBeq)}
    if(!round(S1$N_nurseries[1,max(years),,])==round(S1$N_nurseries[1,max(years)-1,,])){
      Flag_Req <- paste0("for F=",F_fixed_value, " Recruitment is not at equilibrium");print(Flag_Req)}
    # if(any(Flag_Ceq,Flag_SSBeq,Flag_Req)){print}
  }
  ave_y=max(years)-2
  
  ## save outputs
  
  if (value_of_interest %in% c("C","both")){
    
    C_MSY_tot <- S1$C_tot %>%
      melt(value.name = "Catch", varnames = c("Year","Traj")) %>%
      filter(Year>ave_y,Year<max(years)) %>%
      #group_by(Year, Traj) %>%
      #summarize(tot_catch = sum(Catch)) %>%
      #ungroup() %>%
      summarize(median = median(Catch),
                #mode = Find_mode(Catch),
                mean = mean(Catch),
                Q1   = quantile(Catch, probs = 0.025),
                Q2   = quantile(Catch, probs = 0.975)) %>%
      #cbind(F_sim_value) %>%
      mutate(Zone = 'ALL')
    
    C_MSY_temp <- S1$C_spwn_grd %>%
      melt(value.name = "Catch", varnames = c("Age","Year","Zone","Traj")) %>%
      filter(Year>ave_y,Year<max(years))
    
    C_MSY <- S1$C_regions %>%
      melt(value.name = "Catch", varnames = c("Age","Year","Zone","Traj")) %>%
      filter(Year>ave_y,Year<max(years)) %>%  
      bind_rows(C_MSY_temp) %>%
      group_by(Year,Zone,Traj) %>%
      summarize(Catch=sum(Catch,na.rm=T),.groups='drop') %>%
      group_by(Zone) %>%
      summarise(median = median(Catch),
                #mode = Find_mode(Catch),
                mean = mean(Catch),
                Q1   = quantile(Catch, probs = 0.025),
                Q2   = quantile(Catch, probs = 0.975)) %>%
      mutate(Zone = factor(Zone)) %>%
      bind_rows(C_MSY_tot) %>%
      mutate(F =F_fixed_value, type="C", sim=simulation.name)
  }
  
  if(value_of_interest %in% c("SSB", "both")){
    
    SSB_MSY_tot <- S1$SSB %>%
      melt(value.name = "SSB", varnames = c("Year","Zone","Traj")) %>%
      filter(Year>ave_y,Year<max(years)) %>% 
      group_by(Year, Traj) %>%
      summarize(tot_ssb = sum(SSB),.groups='drop') %>%
      ungroup() %>%
      summarize(median = median(tot_ssb),
                #mode = Find_mode(tot_ssb),
                mean = mean(tot_ssb),
                Q1   = quantile(tot_ssb, probs = 0.025),
                Q2   = quantile(tot_ssb, probs = 0.975)) %>%
      #cbind(F_sim_value) %>%
      mutate(Zone = 'ALL')
    
    
    SSB_MSY <- S1$SSB %>%
      melt(value.name = "SSB", varnames = c("Year","Zone","Traj")) %>%
      filter(Year>ave_y,Year<max(years)) %>%
      group_by(Zone) %>%
      summarise(median = median(SSB),
                #mode = Find_mode(SSB),
                mean = mean(SSB),
                Q1   = quantile(SSB, probs = 0.025),
                Q2   = quantile(SSB, probs = 0.975)) %>%
      mutate(Zone = factor(Zone)) %>%
      bind_rows(SSB_MSY_tot) %>%
      mutate(F =F_fixed_value, type="SSB",
             sim=simulation.name)
    
    
  }
  if(value_of_interest == 'both'){
    MSY_values <- bind_rows(C_MSY, SSB_MSY)
    rm(list = c("S1", "C_MSY","SSB_MSY"))
  }
  if(value_of_interest=="C"){
    MSY_values <- C_MSY
    rm(list = c("S1", "C_MSY"))
  }
  if(value_of_interest=="SSB"){
    MSY_values <- SSB_MSY
    rm(list = c("S1","SSB_MSY"))
  }
  
  if(!is.null(comments)){
    MSY_values <- MSY_values %>% 
      mutate(comments=comments)}
  
  
    # file_path <- paste0(path_simulation,'/MSY_',F_dim, 'F_',simulation.name,'_', value_of_interest, '.txt')
  file_path <- ifelse(grepl(pattern="R/",getwd()),paste0("../../",path_simulation,'/MSY_',simulation.name,'.txt'),
                      paste0(path_simulation,'/MSY_',simulation.name,'.txt'))
  
  # errorfile_path <- ifelse(grepl(pattern="R/",getwd()),paste0("../../",path_simulation,'/MSY_',simulation.name,'.txt'),
  #                     paste0(path_simulation,'/MSYerrors_',simulation.name,'.txt'))
  
  write.table(MSY_values, file = file_path, append = TRUE, 
              col.names = !file.exists(file_path),
              row.names = FALSE)
  
  # write.table(MSY_error, file = errorfile_path , append = TRUE, 
  #             col.names = !file.exists(errorfile_path),
  #             row.names = FALSE)
  # 
  
  return(F_fixed_value)
}

