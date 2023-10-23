setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#install.packages("emplik")
library(emplik)
#install.packages("moments")
library(moments)

##### Log likelihood for B estimates 
## data: Site frequency spectra that B will be estimated from 

logl<- function(gamma, data, c=0) {
  M=length(data)-1
  y=((1+c):(M-1-c))
  C=length(data[(2+c):(M-c)])
  
  Ly=data[(2+c):(M-c)]
  
  p1<- sum((Ly*gamma*y)/M)
  p2<- sum(Ly*log(M/(y*(M-y))))
  p3<- sum(Ly)*log(sum(exp(gamma*y/M)*(M/(y*(M-y)))))
  logl<- p1+p2-p3
  return(-logl)
}

##### Log likelihood for B estimates with ry
## Symm: GC conservative spectra

logl_wnuissance<- function(gamma, data, c=0, Symm) {
  M=length(data)-1
  y=((1+c):(M-1-c))
  C=length(data[(2+c):(M-c)])
  neutral_exp=1/(1:(M-1))+1/((M-1):1)
  # For focal spectra 
  Ly=data[(2+c):(M-c)]
  
  ### For GC conservartive spectra
  Lyn=Symm[(2+c):(M-c)]
  
  Lyn_f=Lyn+neutral_exp/sum(neutral_exp)*10^-5
  #  ry= (Lyn/sum(Lyn))/((1/y)*(1/(M-y))/sum((1/y)*(1/(M-y))))
  
  ry=(Lyn_f/sum(Lyn_f))/((1/y)*(1/(M-y))/sum((1/y)*(1/(M-y))))
  
  p1<- sum(Ly*log(ry))
  p2<- sum((Ly*gamma*y)/M)
  p3<- sum(Ly*log(M/(y*(M-y))))
  p4<- sum(Ly)*log(sum(ry*exp(gamma*y/M)*(M/(y*(M-y)))))
  logl<- p1+p2+p3-p4
  return(-logl)
}


##### Likelihood confidence intervals ##### ##### ##### ##### ##### #####
## dataf: Site frequency spectra 

WilksLLR_stat<- function(B=0, dataf) {
  ML<- optim(1, logl,data=dataf,c=0, method = "BFGS")$par
  loglML<- -logl(gamma=ML, data=dataf)
  logl0<-  -logl(gamma=B, data=dataf)
  LLR<- 2*(loglML-logl0)
  return(LLR)
}

myfun <- function(B, dataf){
  temp <- WilksLLR_stat(B=B, dataf=dataf)
  list("-2LLR"=temp)
}

##### Likelihood CI with ry parameter

WilksLLR_stat_wry<- function(B=0, dataf, Symmf) {
  ML<- optim(1, logl_wnuissance,data=dataf,c=0,Symm=Symmf ,method = "BFGS")$par
  loglML<- -logl_wnuissance(gamma=ML, data=dataf, Symm = Symmf)
  logl0<-  -logl_wnuissance(gamma=B, data=dataf, Symm=Symmf)
  LLR<- 2*(loglML-logl0)
  return(LLR)
}

myfun_wry <- function(B, dataf, Symmf){
  temp <- WilksLLR_stat_wry(B=B, dataf=dataf, Symmf=Symmf)
  list("-2LLR"=temp)
}


######### Mutation bias conditional on B estimates 
## df: population count data (i.e., Population_data/analyzedSIsites_mel69_autosome)

generalized_harm<- function(gamma, M) {
  y=(1:(M-1))
  res<- sum(exp(gamma*y/M)*(M/(y*(M-y))))
  return(res)
}

corrected_biases_VB15<- function(df,M=69, psig=0.05) {
  results<- matrix(0,1,8)
  colnames(results)<- c("gamma","alpha_VB15","beta_VB15","beta_CI1","beta_CI2","lambda", "eqGC","theta_VB15")
  
  GCspec<- GC_spectra_build(df, M=M)
  
  B_ML<-optim(1, logl,data=GCspec$GC,c=0, method = "BFGS")$par
  results[1,1]<-B_ML
  
  genharm_plus<- generalized_harm(B_ML, M)
  genharm_neg<- generalized_harm(-B_ML, M)*exp(B_ML)
  genharm_total<- genharm_plus + genharm_neg
  
  ll1<- -optim(1, logl,data=GCspec$GC,c=0, method = "BFGS")$value
  ll0<- -logl(gamma=0, data=GCspec$GC,c=0)
  LRTpval<- pchisq(2*(ll1-ll0), df=1, lower.tail = FALSE)
  
  
  likCI<-findUL(fun=myfun,MLE = B_ML, dataf=GCspec$GC)
  gammaLB<- likCI$Low

  genharm_plus_LB<- generalized_harm(gammaLB, M)
  genharm_neg_LB<- generalized_harm(-gammaLB, M)*exp(gammaLB)
  genharm_total_LB<- genharm_plus_LB + genharm_neg_LB
  
  gammaUB<- likCI$Up
  
  genharm_plus_UB<- generalized_harm(gammaUB, M)
  genharm_neg_UB<- generalized_harm(-gammaUB, M)*exp(gammaUB)
  genharm_total_UB<- genharm_plus_UB+genharm_neg_UB
  
  
  if(LRTpval<=psig) {
    ## VB 15
    tmpalpha<- GCspec$GC[(M+1)]/sum(GCspec$GC) + (sum(GCspec$GC[2:M]) * genharm_neg)/(sum(GCspec$GC)*genharm_total)
    tmpalpha_LB<- GCspec$GC[(M+1)]/sum(GCspec$GC) + (sum(GCspec$GC[2:M]) * genharm_neg_LB)/(sum(GCspec$GC)*genharm_total_LB)
    tmpalpha_UB<- GCspec$GC[(M+1)]/sum(GCspec$GC) + (sum(GCspec$GC[2:M]) * genharm_neg_UB)/(sum(GCspec$GC)*genharm_total_UB)
    
    ### 
    results[1,2]<- tmpalpha/(tmpalpha+exp(B_ML)*(1-tmpalpha)) # alpha
    alpha_LB<- tmpalpha_LB/(tmpalpha_LB+exp(gammaLB)*(1-tmpalpha_LB))
    alpha_UB<- tmpalpha_UB/(tmpalpha_UB+exp(gammaUB)*(1-tmpalpha_UB))
    
    results[1,3]<- 1- results[1,2] # beta
    results[1,4]<- 1- alpha_LB
    results[1,5]<- 1- alpha_UB
    
    results[1,6]<- results[1,3]/results[1,2] # lambda
    results[1,7]<- 1/(1+results[1,6]*exp(-B_ML)) # eqGC
    
    tmptheta<- sum(GCspec$GC[2:M])/(sum(GCspec$GC)*genharm_total)
    results[1,8]<- tmptheta*exp(B_ML)/tmpalpha + tmptheta/(1-tmpalpha) # theta
    
  } else {
    ## V14
    tmpalpha<- (GCspec$GC[(M+1)]+1/2*sum(GCspec$GC[2:M]))/sum(GCspec$GC) 
    
    results[1,2]<- tmpalpha
    results[1,3]<- 1- results[1,2]
    results[1,4]<- 1- results[1,2]
    results[1,5]<- 1- results[1,2]
    results[1,6]<- results[1,3]/results[1,2]
    results[1,7]<- mean(df$GCcontent)
    HM <- sum(1/(1:(M-1)))
    results[1,8]<- sum(GCspec$GC[2:M])/(2*sum(GCspec$GC)*HM)
    
  }
  return(results)
}

### Corrected mutation bias with ry parameter
corrected_biases_VB15_wry<- function(df,M=69, psig=0.05) {
  results<- matrix(0,1,8)
  colnames(results)<- c("gamma","alpha_VB15","beta_VB15","beta_CI1","beta_CI2","lambda", "eqGC","theta_VB15")
  
  GCspec<- GC_spectra_build(df, M=M)
  all_spec<- spectra6_build(df,M=M)
  
  B_ML<-optim(1, logl_wnuissance,data=GCspec$GC,c=0,Symm=all_spec$Symm, method = "BFGS")$par
  results[1,1]<-B_ML
  
  genharm_plus<- generalized_harm(B_ML, M)
  genharm_neg<- generalized_harm(-B_ML, M)*exp(B_ML)
  genharm_total<- genharm_plus + genharm_neg
  
  ll1<- -optim(1, logl_wnuissance,data=GCspec$GC,c=0,Symm=all_spec$Symm, method = "BFGS")$value
  ll0<- -logl_wnuissance(gamma=0, data=GCspec$GC,c=0,Symm=all_spec$Symm)
  LRTpval<- pchisq(2*(ll1-ll0), df=1, lower.tail = FALSE)
  
  
  likCI<-findUL(fun=myfun_wry,MLE = B_ML, dataf=GCspec$GC, Symmf=all_spec$Symm)
  gammaLB<- likCI$Low

  genharm_plus_LB<- generalized_harm(gammaLB, M)
  genharm_neg_LB<- generalized_harm(-gammaLB, M)*exp(gammaLB)
  genharm_total_LB<- genharm_plus_LB + genharm_neg_LB
  
  gammaUB<- likCI$Up
  
  genharm_plus_UB<- generalized_harm(gammaUB, M)
  genharm_neg_UB<- generalized_harm(-gammaUB, M)*exp(gammaUB)
  genharm_total_UB<- genharm_plus_UB+genharm_neg_UB
  
  
  if(LRTpval<=psig) {
    ## VB 15
    tmpalpha<- GCspec$GC[(M+1)]/sum(GCspec$GC) + (sum(GCspec$GC[2:M]) * genharm_neg)/(sum(GCspec$GC)*genharm_total)
    tmpalpha_LB<- GCspec$GC[(M+1)]/sum(GCspec$GC) + (sum(GCspec$GC[2:M]) * genharm_neg_LB)/(sum(GCspec$GC)*genharm_total_LB)
    tmpalpha_UB<- GCspec$GC[(M+1)]/sum(GCspec$GC) + (sum(GCspec$GC[2:M]) * genharm_neg_UB)/(sum(GCspec$GC)*genharm_total_UB)
    
    ### 
    results[1,2]<- tmpalpha/(tmpalpha+exp(B_ML)*(1-tmpalpha)) # alpha
    alpha_LB<- tmpalpha_LB/(tmpalpha_LB+exp(gammaLB)*(1-tmpalpha_LB))
    alpha_UB<- tmpalpha_UB/(tmpalpha_UB+exp(gammaUB)*(1-tmpalpha_UB))
    
    results[1,3]<- 1- results[1,2] # beta
    results[1,4]<- 1- alpha_LB
    results[1,5]<- 1- alpha_UB
    
    results[1,6]<- results[1,3]/results[1,2] # lambda
    results[1,7]<- 1/(1+results[1,6]*exp(-B_ML)) # eqGC
    
    tmptheta<- sum(GCspec$GC[2:M])/(sum(GCspec$GC)*genharm_total)
    results[1,8]<- tmptheta*exp(B_ML)/tmpalpha + tmptheta/(1-tmpalpha) # theta
    
  } else {
    ## V14
    tmpalpha<- (GCspec$GC[(M+1)]+1/2*sum(GCspec$GC[2:M]))/sum(GCspec$GC) 
    
    results[1,2]<- tmpalpha
    results[1,3]<- 1- results[1,2]
    results[1,4]<- 1- results[1,2]
    results[1,5]<- 1- results[1,2]
    results[1,6]<- results[1,3]/results[1,2]
    results[1,7]<- mean(df$GCcontent)
    HM <- sum(1/(1:(M-1)))
    results[1,8]<- sum(GCspec$GC[2:M])/(2*sum(GCspec$GC)*HM)
    
  }
  return(results)
}



### Example use;
source("Pol_estimate_functions.R")
mel69_commonreduced5p_spectras_autosome <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_autosome", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)


B_WChr<- optim(1, logl,data=mel69_commonreduced5p_spectras_autosome$Asymm,c=0, method = "BFGS")$par


B_WChr_ll1<- -optim(1, logl,data=mel69_commonreduced5p_spectras_autosome$Asymm,c=0, method = "BFGS")$value
B_WChr_ll0<- -logl(gamma=0, data=mel69_commonreduced5p_spectras_autosome$Asymm,c=0)

likCI_WChr<-findUL(fun=myfun,MLE = B_WChr, dataf=mel69_commonreduced5p_spectras_autosome$Asymm)
likCI_WChr$Low
likCI_WChr$Up
