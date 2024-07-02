##################################################
# EFHimpacts

#Script: for a new value of Wbar and initial values of other parameters, returns a new h


# Juliette Champagnat
##########################################

Wbar_change <- function(Wbar_new,h,K,Mlarvae,Wbar_0,d_l=1,d_r=1){
  rDI_0 <- -log(4*h/(exp(-Mlarvae*d_l)*Wbar_0*(1-h)))/d_r
  #computation of vector of h and K associated with the new wbar
  vec_h <- exp(-Mlarvae*d_l)*exp(-rDI_0*d_r)*Wbar_new/(4+(exp(-Mlarvae*d_l)*exp(-rDI_0*d_r)*Wbar_new))
  return(vec_h)
}  