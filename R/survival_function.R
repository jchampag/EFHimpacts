##################################################
# EFHimpacts

#Script: function which for a given a computes survival until this age


# Juliette Champagnat
##################################################

#arguments
#' @a  age until with survival is computed
#' @Fval fishing mortality value (will be multiplied with a selectivity vector)
#' @data dataset used for grabbing parameters
#' @first_age first_sge considered for survival computation, if 1 the computation starts at an age post recruitment, of 0 the survival density independant from recruitment phase is also considered
                             

"survival_computing" <- function(a,Fval,data=mydata,first_age=1,M_spwgrd_surv=1){

  if(data$age_rec==data$age_sa){ #for this special case no age suffer Mnursery because the fish is directly out of the nursery (so in the transition case)
    S_a <- exp(
      ifelse(a==(data$age_rec-min(data$ages)+1),
             0, #si on est à age rec on ne subit que la mortalité de la phase de recrutement (qui est codée en dessous)
  
      ifelse(a==(data$age_sa-min(data$ages)+2),#si on est = à l'age de sortie des nour +1 on subit les 2 mortalités pondérée par la durée
               -(data$d_sa*(data$M_nurseries[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)])/365)-
                 ((365-data$d_sa)/365*(data$M_regions[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)])),
  
      ifelse(a<=(data$age_mat-min(data$ages)+1), # si on est > à l'age de sortie des nour +1 mais <= a l'age de maturité
                                                # on subit la mortalité transition pondéreé, et la morta adulte
                -(data$d_sa*(data$M_nurseries[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)])/365)-
                ((365-data$d_sa)/365*(data$M_regions[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)]))-
                sum(data$M_regions[(data$age_sa-min(data$ages)+2):(a-min(data$ages)),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+2):(a-min(data$ages))]),
  
               ## si on est > à l'age de maturité +1 on subit la juv, celle de transition pondéreé, et la morta adulte combinée à celle des frayères
                  # cumulated juvenile mortality
                 # ponderated mortality of juv and adult at age_sa
                  -(data$d_sa*(data$M_nurseries[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)])/365)-
                    ((365-data$d_sa)/365*(data$M_regions[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)]))-
                  # unmature adult mortality
                    sum(data$M_regions[(data$age_sa-min(data$ages)+2):(data$age_mat-min(data$ages)),1,1])-
                  # mature adult mortality: ponderated between regions and spawning grounds
                      sum((365-data$ts_sp)/365*data$M_regions[(data$age_mat-min(data$ages)+1):(a-min(data$ages)),1,1]+
                        data$ts_sp/365*data$M_spwn_grd[(data$age_mat-min(data$ages)+1):(a-min(data$ages)),1,1])-
                  # unponderated fishing mortality because no differences between regions and spawning grounds
                         sum(Fval*Selec[(data$age_sa-min(data$ages)+2):(a-min(data$ages)+1)])))
      ) #end of 1st ifelse 
      ) # closing exp
  }else{
    S_a <- exp(
      ifelse(a==(data$age_rec-min(data$ages)+1),
             0, #si on est à age rec on ne subit que la mortalité de la phase de recrutement (qui est codée en dessous)
             
             ifelse(a<=(data$age_sa-min(data$ages)+1), #si on est < à l'age de sortie de nourricerie, on ne subit que la mortalité nat des nourriceries
                    -sum(data$M_nurseries[1:(a-min(data$ages)),1,1]-Fval*Selec[1:(a-min(data$ages))]),
                    
                    ifelse(a==(data$age_sa-min(data$ages)+2),#si on est = à l'age de sortie +1 des nour on subit les 2 mortalités pondérée par la durée
                           -sum(data$M_nurseries[1:(data$age_sa-min(data$ages)),1,1]+Fval*Selec[1:(data$age_sa-min(data$ages))])-
                             (data$d_sa*(data$M_nurseries[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)])/365)-
                             ((365-data$d_sa)/365*(data$M_regions[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)])),
                           
                           ifelse(a<=(data$age_mat-min(data$ages)+1), # si on est > à l'age de sortie des nour +1 mais <= a l'age de maturité
                                  # on subit la juv, celle de transition pondéreé, et la morta adulte
                                  -sum(data$M_nurseries[1:(data$age_sa-min(data$ages)),1,1]+Fval*Selec[1:(data$age_sa-min(data$ages))])-
                                    (data$d_sa*(data$M_nurseries[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)])/365)-
                                    ((365-data$d_sa)/365*(data$M_regions[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)]))-
                                    sum(data$M_regions[(data$age_sa-min(data$ages)+2):(a-min(data$ages)),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+2):(a-min(data$ages))]),
                                  
                                  ## si on est > à l'age de maturité +1 on subit la juv, celle de transition pondéreé, et la morta adulte combinée à celle des frayères
                                  # cumulated juvenile mortality
                                  -sum(data$M_nurseries[1:(data$age_sa-min(data$ages)),1,1]+Fval*Selec[1:(data$age_sa-min(data$ages))])-
                                    # ponderated mortality of juv and adult at age_sa
                                    (data$d_sa*(data$M_nurseries[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)])/365)-
                                    ((365-data$d_sa)/365*(data$M_regions[(data$age_sa-min(data$ages)+1),1,1]+Fval*Selec[(data$age_sa-min(data$ages)+1)]))-
                                    # unmature adult mortality
                                    sum(data$M_regions[(data$age_sa-min(data$ages)+2):(data$age_mat-min(data$ages)),1,1])-
                                    # mature adult mortality: ponderated between regions and spawning grounds
                                    sum((365-data$ts_sp)/365*data$M_regions[(data$age_mat-min(data$ages)+1):(a-min(data$ages)),1,1]+
                                          data$ts_sp/365*data$M_spwn_grd[(data$age_mat-min(data$ages)+1):(a-min(data$ages)),1,1])-
                                    # unponderated fishing mortality because no differences between regions and spawning grounds
                                    sum(Fval*Selec[(data$age_sa-min(data$ages)+2):(a-min(data$ages)+1)]))))
      ) #end of 1st ifelse 
    ) # closing exp
    
  }


  if(first_age==0){S_a <- S_a* (4*h/(data$Wbar*((1-h))))} #si on veut prendre en compte la mortalité DI de la phase de recrutement on l'ajoute a la survie
  if(a>(data$age_mat-min(data$ages)+1)){S_a <- S_a * M_spwgrd_surv^(a-1-data$age_mat+1)}

  return(S_a)
}
  

