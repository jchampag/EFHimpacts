##################################################
# EFHimpacts

#Script: function which computes Wbar and adds it to mydata() list


# Juliette Champagnat
##################################################

"Wbar_computing" <- function(a,data=mydata,output = "W_a",M_spwgrd_surv=1){
  
  if(grepl(pattern="/R",getwd())){source("../chap3_git/survival_function.R",local=T)
  }else{source("R/chap3_git/survival_function.R",local=T)}

  #compute survival until age a
   S_a <- survival_computing(a=a,data=mydata,Fval=0,first_age = 1,M_spwgrd_surv=M_spwgrd_surv)
   # S_a <- survival_computing_old(a=a,data=mydata,Fval=0,first_age = 1)
    
    #puis le nombre d'oeuf moyen
    W_a = S_a*data$Pf[a-min(ages)+1]*data$mat[a-min(ages)+1]*data$fec[a-min(ages)+1] ##a+1 car ces vecteurs commencent à l'age 0
    
  if(output == "S_a") {return(S_a)}
    else{return(W_a)}
}


mydata$Wbar <- sum(map_dbl(mydata$ages,Wbar_computing,M_spwgrd_surv=M_spwgrd_surv),na.rm=T)

