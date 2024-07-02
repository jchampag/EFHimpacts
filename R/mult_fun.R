##################################################
# EFHimpacts

#Script: this function return the value of h and K for values of spawning ground and/or nursery multipliers


# Juliette Champagnat 

##################################################


"mult_fun_surv" <- function(surv_wDI,surv_rDI,lambda_rDD_type,lambda_surf=1,h,K,Mlarvae,Wbar,d_w,d_l,d_r){
  
  ### BASELINE - computation of M_DI_0 and M_DD_0
  rDI_0 <- -log(4*h/(exp(-Mlarvae*d_l)*Wbar*(1-h)))/d_r
  rDD_0 <- rDI_0/(K*(exp(rDI_0*d_r)-1)) 
  
  ## twist to get lambda_wDI and lambda_rDI from alpha and initial mortality rates
  lambda_wDI <- 1-log(surv_wDI)/(Mlarvae*d_l)
  lambda_rDI <- 1-log(surv_rDI)/(rDI_0*d_r)
  
  ## what lambda_rDD value should be?
  if(lambda_rDD_type == 1){lambda_rDD <- rep(1,length(lambda_rDI))}else{lambda_rDD <- lambda_rDI}

  #computation of vector of h and K associated with the multiplicators
  vec_h <- exp(-lambda_wDI*Mlarvae*d_l)*exp(-lambda_rDI*rDI_0*d_r)*Wbar/(4+(exp(-lambda_wDI*Mlarvae*d_l)*exp(-lambda_rDI*rDI_0*d_r)*Wbar))
  vec_K_temp <- (lambda_rDI*rDI_0)/(lambda_rDD*rDD_0*(exp(lambda_rDI*rDI_0*d_r)-1))
  
  #adding lambda_surf effects
  vec_K <-vec_K_temp/lambda_surf
  
  #defining sim type according to multiplier
  sim_type <- NULL
  for (i in seq_along(lambda_wDI)){
    if(lambda_wDI[i]!=1){ #spw ground modif on
      if(lambda_rDI[i]==1){#nurs Di off
        if(lambda_surf[i]==1){#nurs_surf off
          sim_type_temp="spawning"
        }else{#nurs surf on
          sim_type_temp="spw_surf"
        }
      }else{ #nursery rDI also on
        if(lambda_rDD_type==1){ #nursery DD is off 
          if(lambda_surf[i]==1){sim_type_temp="spw_nurs_DI"
          }else{sim_type_temp="spw_nurs_DI_surf"}
        }else{#nursery DD = nurs DI (called nurs_qual)
          if(lambda_surf[i]==1){sim_type_temp="spw_nurs_qual"
          }else{sim_type_temp="spw_nurs_qual_surf"}
        }
      }
    }else{#spw ground modif off
      if(lambda_rDI[i]==1){ #nurs quality off
        if(lambda_surf[i]==1){ #nursery surface off as well
          sim_type_temp="ref"
        }else{ # surface on
          sim_type_temp="nurs_surf"}
      }else{ # nursery DI on
        if(lambda_rDD_type==1){ # nursery DD is off
          if(lambda_surf[i]==1){ #nursery surface off
            sim_type_temp="nurs_DI"
          }else{ #nursery surface on
            sim_type_temp="nurs_DI_surf" }
        }else{#nursery DD = nurs DI (called nurs_qual)
          if(lambda_surf[i]==1){ #nursery surface off
            sim_type_temp="nurs_qual"
            }else{ #nursery surface on
            sim_type_temp="nurs_qual_surf" }
        }
      }
    }
    
      
    sim_type <- c(sim_type,sim_type_temp)}
  
  return(data.frame(sim_type=sim_type,
                    surv_wDI=surv_wDI, lambda_wDI=lambda_wDI,
                    surv_rDI=surv_rDI,lambda_rDI=lambda_rDI,
                    lambda_rDD=lambda_rDD,lambda_surf=lambda_surf,
                    h=vec_h,K=vec_K))
}



