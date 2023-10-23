setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(readr)
### Site frequency spectra build from population data  
# df: population count data i.e., Population_data/analyzedSIsites_mel69_autosome)
# M: sample size (i.e., M=69 for D.melanogaster)
spectra6_build<- function(df,M) {
  polfor_AC<-df[which(df$A!=0 & df$C!=0),]
  polfor_AG<-df[which(df$A!=0 & df$G!=0),]
  polfor_AT<-df[which(df$A!=0 & df$T!=0),]
  polfor_CG<-df[which(df$C!=0 & df$G!=0),]
  polfor_CT<-df[which(df$C!=0 & df$T!=0),]
  polfor_GT<-df[which(df$G!=0 & df$T!=0),]
  monofor_A<-df[which(df$State=="M" & df$A!=0),]
  monofor_T<-df[which(df$State=="M" & df$T!=0),]
  monofor_G<-df[which(df$State=="M" & df$G!=0),]
  monofor_C<-df[which(df$State=="M" & df$C!=0),]
  spectra<- matrix(0,M+1,6) 
  colnames(spectra)<- c("AC","AG","AT","CG","CT","GT")
  rownames(spectra)<- c(0:M)
  # Place monomorphic
  spectra[1,1]<- nrow(monofor_A)
  spectra[1,2]<- nrow(monofor_A)
  spectra[1,3]<- nrow(monofor_A)
  spectra[1,4]<- nrow(monofor_C)
  spectra[1,5]<- nrow(monofor_C)
  spectra[1,6]<- nrow(monofor_G)
  spectra[M+1,1]<- nrow(monofor_C)
  spectra[M+1,2]<- nrow(monofor_G)
  spectra[M+1,3]<- nrow(monofor_T)
  spectra[M+1,4]<- nrow(monofor_G)
  spectra[M+1,5]<- nrow(monofor_T)
  spectra[M+1,6]<- nrow(monofor_T)
  # Place polymorphic
  for (i in 1:(M-1)) {
    spectra[i+1,1]<- sum(polfor_AC$C==i)
    spectra[i+1,2]<- sum(polfor_AG$G==i)
    spectra[i+1,3]<- sum(polfor_AT$T==i)
    spectra[i+1,4]<- sum(polfor_CG$G==i)
    spectra[i+1,5]<- sum(polfor_CT$T==i)
    spectra[i+1,6]<- sum(polfor_GT$T==i)
  }
  
  spectra<- as.data.frame(spectra)
  spectra$Symm=spectra$AT+spectra$CG
  spectra$Asymm=spectra$AC+spectra$AG+
    spectra$CT[(M+1):1]+spectra$GT[(M+1):1]
  
  
  return(spectra)
}



###### Log likelihood b and e prime
Logl_EandBPrime <- function(ePrimebPrime, data){
  
  M <- dim(data)[1]-1
  y <- 1:(M - 1)
  # Polymorphics
  L_ij_vec <- data[2:M,1:6]
  L_AC_vec <- L_ij_vec$AC
  L_AG_vec <- L_ij_vec$AG
  L_CT_vec <- L_ij_vec$CT
  L_GT_vec <- L_ij_vec$GT
  
  ePrime <- ePrimebPrime[1]
  bPrime <- ePrimebPrime[2]
  logL_eb <- 
    sum((L_AC_vec + L_GT_vec[M-y])*log(ePrime/(M-y) + (1-bPrime)/y) +
          (L_AG_vec + L_CT_vec[M-y])*log((1-ePrime)/(M-y) + bPrime/y))
  return(-logL_eb)
}

##### Function to estimate Mutation rate matrix from site frequency spectrum 
# input df: SFS)


allrate_est<- function(df, region="name", chr="chr") {
  results<- matrix(0,1,10)
  colnames(results)<- c("Data","Chromosome","a","f", "b","d","c","e", "beta","theta")
  results[1,1]<- region
  results[1,2]<- chr
  M <- dim(df)[1]-1		# number of individuals 
  HM <- sum(1/(1:(M-1)))
  y <- 1:(M-1)
  oneOverYSum <- 1/(M-y) + 1/y
  
  # Monomorphics
  L_A <- df$AC[1]
  L_C <- df$CG[1]
  L_G <- df$GT[1]
  L_T <- df$GT[M+1]
  L_i <- c(L_A, L_C, L_G, L_T)
  # Polymorphics
  L_ij_vec <- df[2:M,1:6]
  L_ij <- colSums(L_ij_vec)
  L_AC_vec <- L_ij_vec$AC
  L_AG_vec <- L_ij_vec$AG
  L_AT_vec <- L_ij_vec$AT
  L_CG_vec <- L_ij_vec$CG
  L_CT_vec <- L_ij_vec$CT
  L_GT_vec <- L_ij_vec$GT
  L_AC <- L_CA <- sum(L_AC_vec)
  L_AG <- L_GA <- sum(L_AG_vec)
  L_AT <- L_TA <- sum(L_AT_vec)
  L_CG <- L_GC <- sum(L_CG_vec)
  L_CT <- L_TC <- sum(L_CT_vec)
  L_GT <- L_TG <- sum(L_GT_vec)
  # Totals
  L_m <- sum(L_i)
  L_p <- sum(L_ij)
  L_total <- L_m + L_p
  
  ###
  L_overATCG <- L_AC + L_AG + L_TC + L_TG
  L_overAT <- L_A + L_T + L_AT
  L_overCG <- L_C + L_G + L_CG
  #
  denominatorAT <- L_overAT + 0.5*L_overATCG
  denominatorCG <- L_overCG + 0.5*L_overATCG
  
  eandbPrime<- optim(c(0.5,0.5), Logl_EandBPrime, data=df, method = "Nelder-Mead",control=list(reltol=1.e-16))$par
  
  ePrime <- eandbPrime[1]								
  bPrime <- eandbPrime[2]	
  
  results[1,3]<- L_AT/denominatorAT/HM # a 
  results[1,4]<- L_CG/denominatorCG/HM # f
  
  results[1,5]<- bPrime*L_overATCG/denominatorCG/2/HM  # b: CT and GA (Transition)
  results[1,6]<- (1-bPrime)*L_overATCG/denominatorCG/2/HM # d: GT and CA
  
  results[1,7]<- (1-ePrime)*L_overATCG/denominatorAT/2/HM # c: AG and TC (Transition)
  results[1,8]<- ePrime*L_overATCG/denominatorAT/2/HM # e: AC and TG
  
  results[1,9]<- denominatorAT/L_total
  results[1,10]<- L_overATCG/(L_total*2*HM)
  return(results)
}

####### Intron resample
# input df: population count data i.e., Population_data/analyzedSIsites_mel69_autosome)


sample_5SI<- function(df) {  
  intron_count<-length(unique(df$name))
  names<-unique(df$name)
  samp_intron<- sample(names, intron_count, replace = TRUE)
  idx<-c()
  for (i in 1:length(samp_intron)) {
    idx<- append(idx,which(df$name==samp_intron[i]))
  }
  new_samp<- df[idx,]
  return(new_samp)
}

######### Bootstrap estimates build (input )
# input df: population count data i.e., Population_data/analyzedSIsites_mel69_autosome)

bootstrap_pol_6par<- function(df, region="name", chr="chr", BS=1000, M=69) {
  results<- matrix(0,BS,10)
  colnames(results)<- c("Data","Chromosome","a","f", "b","d","c","e", "beta","theta")
  results[1:BS,1]<- region
  results[1:BS,2]<- chr
  
  for (b in 1:BS) {
    new_samp<- sample_5SI(df)
    spectra<- spectra6_build(new_samp, M=M)
    
    M <- dim(spectra)[1]-1		# number of individuals 
    HM <- sum(1/(1:(M-1)))
    y <- 1:(M-1)
    oneOverYSum <- 1/(M-y) + 1/y
    
    # Monomorphics
    L_A <- spectra$AC[1]
    L_C <- spectra$CG[1]
    L_G <- spectra$GT[1]
    L_T <- spectra$GT[M+1]
    L_i <- c(L_A, L_C, L_G, L_T)
    # Polymorphics
    L_ij_vec <- spectra[2:M,]
    L_ij <- colSums(L_ij_vec)
    L_AC_vec <- L_ij_vec$AC
    L_AG_vec <- L_ij_vec$AG
    L_AT_vec <- L_ij_vec$AT
    L_CG_vec <- L_ij_vec$CG
    L_CT_vec <- L_ij_vec$CT
    L_GT_vec <- L_ij_vec$GT
    L_AC <- L_CA <- sum(L_AC_vec)
    L_AG <- L_GA <- sum(L_AG_vec)
    L_AT <- L_TA <- sum(L_AT_vec)
    L_CG <- L_GC <- sum(L_CG_vec)
    L_CT <- L_TC <- sum(L_CT_vec)
    L_GT <- L_TG <- sum(L_GT_vec)
    # Totals
    L_m <- sum(L_i)
    L_p <- sum(L_ij)
    L_total <- L_m + L_p
    
    ###
    L_overATCG <- L_AC + L_AG + L_TC + L_TG
    L_overAT <- L_A + L_T + L_AT
    L_overCG <- L_C + L_G + L_CG
    #
    denominatorAT <- L_overAT + 0.5*L_overATCG
    denominatorCG <- L_overCG + 0.5*L_overATCG
    
    eandbPrime<- optim(c(0.5,0.5), Logl_EandBPrime, data=spectra, method = "Nelder-Mead",control=list(reltol=1.e-16))$par
    
    ePrime <- eandbPrime[1]								
    bPrime <- eandbPrime[2]	
    
    results[b,3]<- L_AT/denominatorAT/HM # a 
    results[b,4]<- L_CG/denominatorCG/HM # f
    
    results[b,5]<- bPrime*L_overATCG/denominatorCG/2/HM  # b: CT and GA (Transition)
    results[b,6]<- (1-bPrime)*L_overATCG/denominatorCG/2/HM # d: GT and CA
    
    results[b,7]<- (1-ePrime)*L_overATCG/denominatorAT/2/HM # c: AG and TC (Transition)
    results[b,8]<- ePrime*L_overATCG/denominatorAT/2/HM # e: AC and TG
    
    results[b,9]<- denominatorAT/L_total
    results[b,10]<- L_overATCG/(L_total*2*HM)
    
  }
  return(results)
}

### Example use;
analyzedSIsites_mel69_autosome <- read_delim("../../Population_data/Melanogaster/analyzedSIsites_mel69_autosome", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)

spectra_auto<- spectra6_build(analyzedSIsites_mel69_autosome, M=69)
whole_auto<- allrate_est(spectra_auto, region = "WChr",chr="Autosome")
whole_auto<- as.data.frame(whole_auto)
sapply(whole_auto,class)
whole_auto[,3:10]<- lapply(whole_auto[,3:10], function(x) {as.numeric(as.character(x))})


## Mutation matrix build 

df=whole_auto

Q_auto_whole <- t(matrix(c(0 , df$a, df$c, df$e, 
                           df$a,   0 , df$e, df$c, 
                           df$b, df$d,   0 , df$f, 
                           df$d, df$b, df$f,   0 ), nrow=4, ncol=4))
diag(Q_auto_whole) <- -rowSums(Q_auto_whole, na.rm=TRUE)	
dimnames(Q_auto_whole) <- list(c("A", "T", "G", "C"), c("A", "T", "G", "C"))


## Expected heterozygosity  

diagQ<-c(Q_auto_whole[1,1], Q_auto_whole[2,2], Q_auto_whole[3,3], Q_auto_whole[4,4])
pi<-1/2*c(whole_auto$beta,whole_auto$beta, 1-whole_auto$beta,1-whole_auto$beta)

sum(-pi*diagQ) # 0.0196
ExpH<-sum(-pi*diagQ)
