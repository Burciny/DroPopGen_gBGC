
####### Sampling

sample_5SI<- function(df) {  ## input countfile
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

# Spectra
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
  return(spectra)
}

###### Log likelihood b and e prime
Logl_EandBPrime <- function(ePrimebPrime, data){
  
  M <- dim(data)[1]-1
  y <- 1:(M - 1)
  # Polymorphics
  L_ij_vec <- data[2:M,]
  L_AC_vec <- L_ij_vec$AC
  L_AG_vec <- L_ij_vec$AG
  L_CT_vec <- L_ij_vec$CT
  L_GT_vec <- L_ij_vec$GT
  
  ePrime <- ePrimebPrime[1]
  bPrime <- ePrimebPrime[2]
  logL_eb <- 
    sum((L_AC_vec + L_GT_vec[M-y])*log(ePrime/(M-y) + (1 - bPrime)/y) +
          (L_AG_vec + L_CT_vec[M-y])*log((1-ePrime)/(M-y) + bPrime/y))
  return(-logL_eb)
}


#### Matrix 10x10 
mat10by10_reduceddat_create<- function(df) {
  matmelandsim<- matrix(0,10,10)
  colnames(matmelandsim)<- c("A","T","G","C","AT","CG","AC","AG","TC","TG")
  # Column 1 
  
  # 1) sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_Simpopstate=="M") 
  matmelandsim[1,1]<-sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_sim_A==2)
  matmelandsim[2,1]<-sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_sim_T==2)
  matmelandsim[3,1]<-sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_sim_G==2)
  matmelandsim[4,1]<-sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_sim_C==2)
  
  # 2) sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_Simpopstate=="P") 
  matmelandsim[5,1]<-sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_sim_A==1 & df$hap_sim_T==1)
  matmelandsim[6,1]<-sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_sim_C==1 & df$hap_sim_G==1)
  matmelandsim[7,1]<-sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_sim_A==1 & df$hap_sim_C==1)
  matmelandsim[8,1]<-sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_sim_A==1 & df$hap_sim_G==1)
  matmelandsim[9,1]<-sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_sim_T==1 & df$hap_sim_C==1)
  matmelandsim[10,1]<-sum(df$hap_Melstate=="M" & df$hap_A==2 & df$hap_sim_T==1 & df$hap_sim_G==1)
  
  # Column 2
  
  # 1) sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_Simpopstate=="M") 
  matmelandsim[1,2]<-sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_sim_A==2)
  matmelandsim[2,2]<-sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_sim_T==2)
  matmelandsim[3,2]<-sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_sim_G==2)
  matmelandsim[4,2]<-sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_sim_C==2)
  
  # 2) sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_Simpopstate=="P") 
  matmelandsim[5,2]<-sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_sim_A==1 & df$hap_sim_T==1)
  matmelandsim[6,2]<-sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_sim_C==1 & df$hap_sim_G==1)
  matmelandsim[7,2]<-sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_sim_A==1 & df$hap_sim_C==1)
  matmelandsim[8,2]<-sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_sim_A==1 & df$hap_sim_G==1)
  matmelandsim[9,2]<-sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_sim_T==1 & df$hap_sim_C==1)
  matmelandsim[10,2]<-sum(df$hap_Melstate=="M" & df$hap_T==2 & df$hap_sim_T==1 & df$hap_sim_G==1)
  
  # Column 3
  # 1) sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_Simpopstate=="M") 
  matmelandsim[1,3]<-sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_sim_A==2)
  matmelandsim[2,3]<-sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_sim_T==2)
  matmelandsim[3,3]<-sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_sim_G==2)
  matmelandsim[4,3]<-sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_sim_C==2)
  
  # 2) sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_Simpopstate=="P") 
  matmelandsim[5,3]<-sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_sim_A==1 & df$hap_sim_T==1)
  matmelandsim[6,3]<-sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_sim_C==1 & df$hap_sim_G==1)
  matmelandsim[7,3]<-sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_sim_A==1 & df$hap_sim_C==1)
  matmelandsim[8,3]<-sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_sim_A==1 & df$hap_sim_G==1)
  matmelandsim[9,3]<-sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_sim_T==1 & df$hap_sim_C==1)
  matmelandsim[10,3]<-sum(df$hap_Melstate=="M" & df$hap_G==2 & df$hap_sim_T==1 & df$hap_sim_G==1)
  
  # Column 4
  
  # 1) sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_Simpopstate=="M") 
  matmelandsim[1,4]<-sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_sim_A==2)
  matmelandsim[2,4]<-sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_sim_T==2)
  matmelandsim[3,4]<-sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_sim_G==2)
  matmelandsim[4,4]<-sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_sim_C==2)
  
  # 2) sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_Simpopstate=="P") 
  matmelandsim[5,4]<-sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_sim_A==1 & df$hap_sim_T==1)
  matmelandsim[6,4]<-sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_sim_C==1 & df$hap_sim_G==1)
  matmelandsim[7,4]<-sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_sim_A==1 & df$hap_sim_C==1)
  matmelandsim[8,4]<-sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_sim_A==1 & df$hap_sim_G==1)
  matmelandsim[9,4]<-sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_sim_T==1 & df$hap_sim_C==1)
  matmelandsim[10,4]<-sum(df$hap_Melstate=="M" & df$hap_C==2 & df$hap_sim_T==1 & df$hap_sim_G==1)
  
  ## Column 5 
  
  # 1) sum(df$hap_A==1 & df$hap_T==1 & df$hap_Simpopstate=="M")
  
  matmelandsim[1,5]<-sum(df$hap_A==1 & df$hap_T==1 & df$hap_sim_A==2)
  matmelandsim[2,5]<-sum(df$hap_A==1 & df$hap_T==1 & df$hap_sim_T==2)
  matmelandsim[3,5]<-sum(df$hap_A==1 & df$hap_T==1 & df$hap_sim_G==2)
  matmelandsim[4,5]<-sum(df$hap_A==1 & df$hap_T==1 & df$hap_sim_C==2)
  
  # 2) sum(df$hap_A==1 & df$hap_T==1 & df$hap_Simpopstate=="P") 
  
  matmelandsim[5,5]<-sum(df$hap_A==1 & df$hap_T==1 & df$hap_sim_A==1 & df$hap_sim_T==1)
  matmelandsim[6,5]<-sum(df$hap_A==1 & df$hap_T==1 & df$hap_sim_C==1 & df$hap_sim_G==1)
  matmelandsim[7,5]<-sum(df$hap_A==1 & df$hap_T==1 & df$hap_sim_A==1 & df$hap_sim_C==1)
  matmelandsim[8,5]<-sum(df$hap_A==1 & df$hap_T==1 & df$hap_sim_A==1 & df$hap_sim_G==1)
  matmelandsim[9,5]<-sum(df$hap_A==1 & df$hap_T==1 & df$hap_sim_T==1 & df$hap_sim_C==1)
  matmelandsim[10,5]<-sum(df$hap_A==1 & df$hap_T==1 & df$hap_sim_T==1 & df$hap_sim_G==1)
  
  ## Column 6
  
  # 1) sum(df$hap_C==1 & df$hap_G==1 & df$hap_Simpopstate=="M")
  
  matmelandsim[1,6]<-sum(df$hap_C==1 & df$hap_G==1 & df$hap_sim_A==2)
  matmelandsim[2,6]<-sum(df$hap_C==1 & df$hap_G==1 & df$hap_sim_T==2)
  matmelandsim[3,6]<-sum(df$hap_C==1 & df$hap_G==1 & df$hap_sim_G==2)
  matmelandsim[4,6]<-sum(df$hap_C==1 & df$hap_G==1 & df$hap_sim_C==2)
  
  # 2) sum(df$hap_C==1 & df$hap_G==1 & df$hap_Simpopstate=="P") 
  
  matmelandsim[5,6]<-sum(df$hap_C==1 & df$hap_G==1 & df$hap_sim_A==1 & df$hap_sim_T==1)
  matmelandsim[6,6]<-sum(df$hap_C==1 & df$hap_G==1 & df$hap_sim_C==1 & df$hap_sim_G==1)
  matmelandsim[7,6]<-sum(df$hap_C==1 & df$hap_G==1 & df$hap_sim_A==1 & df$hap_sim_C==1)
  matmelandsim[8,6]<-sum(df$hap_C==1 & df$hap_G==1 & df$hap_sim_A==1 & df$hap_sim_G==1)
  matmelandsim[9,6]<-sum(df$hap_C==1 & df$hap_G==1 & df$hap_sim_T==1 & df$hap_sim_C==1)
  matmelandsim[10,6]<-sum(df$hap_C==1 & df$hap_G==1 & df$hap_sim_T==1 & df$hap_sim_G==1)
  
  ## Column 7
  
  # 1) sum(df$hap_A==1 & df$hap_C==1 & df$hap_Simpopstate=="M")
  
  matmelandsim[1,7]<-sum(df$hap_A==1 & df$hap_C==1 & df$hap_sim_A==2)
  matmelandsim[2,7]<-sum(df$hap_A==1 & df$hap_C==1 & df$hap_sim_T==2)
  matmelandsim[3,7]<-sum(df$hap_A==1 & df$hap_C==1 & df$hap_sim_G==2)
  matmelandsim[4,7]<-sum(df$hap_A==1 & df$hap_C==1 & df$hap_sim_C==2)
  
  # 2) sum(df$hap_A==1 & df$hap_C==1 & df$hap_Simpopstate=="P") 
  
  matmelandsim[5,7]<-sum(df$hap_A==1 & df$hap_C==1 & df$hap_sim_A==1 & df$hap_sim_T==1)
  matmelandsim[6,7]<-sum(df$hap_A==1 & df$hap_C==1 & df$hap_sim_C==1 & df$hap_sim_G==1)
  matmelandsim[7,7]<-sum(df$hap_A==1 & df$hap_C==1 & df$hap_sim_A==1 & df$hap_sim_C==1)
  matmelandsim[8,7]<-sum(df$hap_A==1 & df$hap_C==1 & df$hap_sim_A==1 & df$hap_sim_G==1)
  matmelandsim[9,7]<-sum(df$hap_A==1 & df$hap_C==1 & df$hap_sim_T==1 & df$hap_sim_C==1)
  matmelandsim[10,7]<-sum(df$hap_A==1 & df$hap_C==1 & df$hap_sim_T==1 & df$hap_sim_G==1)
  
  ## Column 8
  
  # 1) sum(df$hap_A==1 & df$hap_G==1 & df$hap_Simpopstate=="M")
  
  matmelandsim[1,8]<-sum(df$hap_A==1 & df$hap_G==1 & df$hap_sim_A==2)
  matmelandsim[2,8]<-sum(df$hap_A==1 & df$hap_G==1 & df$hap_sim_T==2)
  matmelandsim[3,8]<-sum(df$hap_A==1 & df$hap_G==1 & df$hap_sim_G==2)
  matmelandsim[4,8]<-sum(df$hap_A==1 & df$hap_G==1 & df$hap_sim_C==2)
  
  # 2) sum(df$hap_A==1 & df$hap_G==1 & df$hap_Simpopstate=="P") 
  
  matmelandsim[5,8]<-sum(df$hap_A==1 & df$hap_G==1 & df$hap_sim_A==1 & df$hap_sim_T==1)
  matmelandsim[6,8]<-sum(df$hap_A==1 & df$hap_G==1 & df$hap_sim_C==1 & df$hap_sim_G==1)
  matmelandsim[7,8]<-sum(df$hap_A==1 & df$hap_G==1 & df$hap_sim_A==1 & df$hap_sim_C==1)
  matmelandsim[8,8]<-sum(df$hap_A==1 & df$hap_G==1 & df$hap_sim_A==1 & df$hap_sim_G==1)
  matmelandsim[9,8]<-sum(df$hap_A==1 & df$hap_G==1 & df$hap_sim_T==1 & df$hap_sim_C==1)
  matmelandsim[10,8]<-sum(df$hap_A==1 & df$hap_G==1 & df$hap_sim_T==1 & df$hap_sim_G==1)
  
  ## Column 9
  
  # 1) sum(df$hap_T==1 & df$hap_C==1 & df$hap_Simpopstate=="M")
  
  matmelandsim[1,9]<-sum(df$hap_T==1 & df$hap_C==1 & df$hap_sim_A==2)
  matmelandsim[2,9]<-sum(df$hap_T==1 & df$hap_C==1 & df$hap_sim_T==2)
  matmelandsim[3,9]<-sum(df$hap_T==1 & df$hap_C==1 & df$hap_sim_G==2)
  matmelandsim[4,9]<-sum(df$hap_T==1 & df$hap_C==1 & df$hap_sim_C==2)
  
  # 2) sum(df$hap_T==1 & df$hap_C==1 & df$hap_Simpopstate=="P") 
  
  matmelandsim[5,9]<-sum(df$hap_T==1 & df$hap_C==1 & df$hap_sim_A==1 & df$hap_sim_T==1)
  matmelandsim[6,9]<-sum(df$hap_T==1 & df$hap_C==1 & df$hap_sim_C==1 & df$hap_sim_G==1)
  matmelandsim[7,9]<-sum(df$hap_T==1 & df$hap_C==1 & df$hap_sim_A==1 & df$hap_sim_C==1)
  matmelandsim[8,9]<-sum(df$hap_T==1 & df$hap_C==1 & df$hap_sim_A==1 & df$hap_sim_G==1)
  matmelandsim[9,9]<-sum(df$hap_T==1 & df$hap_C==1 & df$hap_sim_T==1 & df$hap_sim_C==1)
  matmelandsim[10,9]<-sum(df$hap_T==1 & df$hap_C==1 & df$hap_sim_T==1 & df$hap_sim_G==1)
  
  ## Column 10
  
  # 1) sum(df$hap_T==1 & df$hap_G==1 & df$hap_Simpopstate=="M")
  
  matmelandsim[1,10]<-sum(df$hap_T==1 & df$hap_G==1 & df$hap_sim_A==2)
  matmelandsim[2,10]<-sum(df$hap_T==1 & df$hap_G==1 & df$hap_sim_T==2)
  matmelandsim[3,10]<-sum(df$hap_T==1 & df$hap_G==1 & df$hap_sim_G==2)
  matmelandsim[4,10]<-sum(df$hap_T==1 & df$hap_G==1 & df$hap_sim_C==2)
  
  # 2) sum(df$hap_T==1 & df$hap_G==1 & df$hap_Simpopstate=="P") 
  
  matmelandsim[5,10]<-sum(df$hap_T==1 & df$hap_G==1 & df$hap_sim_A==1 & df$hap_sim_T==1)
  matmelandsim[6,10]<-sum(df$hap_T==1 & df$hap_G==1 & df$hap_sim_C==1 & df$hap_sim_G==1)
  matmelandsim[7,10]<-sum(df$hap_T==1 & df$hap_G==1 & df$hap_sim_A==1 & df$hap_sim_C==1)
  matmelandsim[8,10]<-sum(df$hap_T==1 & df$hap_G==1 & df$hap_sim_A==1 & df$hap_sim_G==1)
  matmelandsim[9,10]<-sum(df$hap_T==1 & df$hap_G==1 & df$hap_sim_T==1 & df$hap_sim_C==1)
  matmelandsim[10,10]<-sum(df$hap_T==1 & df$hap_G==1 & df$hap_sim_T==1 & df$hap_sim_G==1)
  return(matmelandsim)
}

## Prelim data frames
premdf_create<- function(df) {
  TotalAutoPolMel<-colSums(df[,5:10]) 
  TotalAutoPolSim<-rowSums(df[5:10,]) 
  # Absolute Divergence in mat
  matTotalAutoDiv=as.matrix(df[1:4,1:4])
  matTotalAutoDiv
  ## Polymorphism in mat
  matTotalAutoPoly=matrix(rep(0,4*4),nrow=4,ncol=4)
  colnames(matTotalAutoPoly)<- c("A","T","G","C")
  ## AT
  matTotalAutoPoly[1,2]=TotalAutoPolMel[1]
  matTotalAutoPoly[2,1]=TotalAutoPolSim[1]
  ## CG
  matTotalAutoPoly[3,4]=TotalAutoPolMel[2]
  matTotalAutoPoly[4,3]=TotalAutoPolSim[2]
  ## AC
  matTotalAutoPoly[1,4]=TotalAutoPolMel[3]
  matTotalAutoPoly[4,1]=TotalAutoPolSim[3]
  ## AG
  matTotalAutoPoly[1,3]=TotalAutoPolMel[4]
  matTotalAutoPoly[3,1]=TotalAutoPolSim[4]
  ## TC
  matTotalAutoPoly[2,4]=TotalAutoPolMel[5]
  matTotalAutoPoly[4,2]=TotalAutoPolSim[5]
  ## TG
  matTotalAutoPoly[2,3]=TotalAutoPolMel[6]
  matTotalAutoPoly[3,2]=TotalAutoPolSim[6]
  # Fold 
  matTotalAutoDiv_fold=matrix(rep(0,4*4),nrow=4,ncol=4)
  matTotalAutoPoly_fold=matrix(rep(0,4*4),nrow=4,ncol=4)
  for(i in 1:3){
    for(j in (i+1):4){
      matTotalAutoDiv_fold[i,j]=matTotalAutoDiv[i,j]+matTotalAutoDiv[j,i]
      matTotalAutoPoly_fold[i,j]=matTotalAutoPoly[i,j]+matTotalAutoPoly[j,i]
    }
  }
  # Vectorize and combine
  vTotalAutoDiv_fold=rep(0,6)
  vTotalAutoPoly_fold=rep(0,6)
  vTotalAutoDiv_fold[1]=matTotalAutoDiv_fold[1,2] # AT 
  vTotalAutoDiv_fold[2]=matTotalAutoDiv_fold[3,4] # CG
  vTotalAutoDiv_fold[3]=matTotalAutoDiv_fold[1,4] # AC
  vTotalAutoDiv_fold[4]=matTotalAutoDiv_fold[1,3] # TG
  vTotalAutoDiv_fold[5]=matTotalAutoDiv_fold[2,4] # TC
  vTotalAutoDiv_fold[6]=matTotalAutoDiv_fold[2,3] # CG
  
  vTotalAutoPoly_fold[1]=matTotalAutoPoly_fold[1,2]
  vTotalAutoPoly_fold[2]=matTotalAutoPoly_fold[3,4]
  vTotalAutoPoly_fold[3]=matTotalAutoPoly_fold[1,4]
  vTotalAutoPoly_fold[4]=matTotalAutoPoly_fold[1,3]
  vTotalAutoPoly_fold[5]=matTotalAutoPoly_fold[2,4]
  vTotalAutoPoly_fold[6]=matTotalAutoPoly_fold[2,3]
  ## Shared pol
  sharedpol<- c()
  allpolmel<- c()
  allpolsim<- c()
  for (i in 5:10) {
    sharedpol[i-4]<- as.numeric(df[i,i])
    allpolmel[i-4]<- sum(df[5:10,i])
    allpolsim[i-4]<- sum(df[i,5:10])
  }
  
  matAutoDivPoly=rbind(vTotalAutoDiv_fold,vTotalAutoPoly_fold,TotalAutoPolSim, TotalAutoPolMel,sharedpol,allpolsim, allpolmel)
  colnames(matAutoDivPoly)=colnames(df[5:10])
  rownames(matAutoDivPoly)=c("Divergence","Polymorphism", "TotalPolSim","TotalPolMel", "Sharedpol","Allpolsim","Allpolmel")
  return(matAutoDivPoly)
}

denoms_create<- function(df) {
  TotalAutoPolMel<-colSums(df[,5:10]) 
  TotalAutoPolSim<-rowSums(df[5:10,]) 
  vAutoATGCSim=rep(0,4)
  vAutoATGCMel=rep(0,4)
  for(i in 1:4){
    vAutoATGCSim[i]=sum(df[i,])
    vAutoATGCMel[i]=sum(df[,i])
  }
  # A- AT- AC- AG
  vAutoATGCSim[1]=vAutoATGCSim[1]+0.5*(TotalAutoPolSim[1]+TotalAutoPolSim[3]+TotalAutoPolSim[4])
  # T- AT- TC- TG
  vAutoATGCSim[2]=vAutoATGCSim[2]+0.5*(TotalAutoPolSim[1]+TotalAutoPolSim[5]+TotalAutoPolSim[6])
  # G- CG- AG- TG
  vAutoATGCSim[3]=vAutoATGCSim[3]+0.5*(TotalAutoPolSim[2]+TotalAutoPolSim[4]+TotalAutoPolSim[6])
  # C- CG- AC- TC
  vAutoATGCSim[4]=vAutoATGCSim[4]+0.5*(TotalAutoPolSim[2]+TotalAutoPolSim[3]+TotalAutoPolSim[5])
  
  vAutoATGCMel[1]=vAutoATGCMel[1]+0.5*(TotalAutoPolMel[1]+TotalAutoPolMel[3]+TotalAutoPolMel[4])
  vAutoATGCMel[2]=vAutoATGCMel[2]+0.5*(TotalAutoPolMel[1]+TotalAutoPolMel[5]+TotalAutoPolMel[6])
  vAutoATGCMel[3]=vAutoATGCMel[3]+0.5*(TotalAutoPolMel[2]+TotalAutoPolMel[4]+TotalAutoPolMel[6])
  vAutoATGCMel[4]=vAutoATGCMel[4]+0.5*(TotalAutoPolMel[2]+TotalAutoPolMel[3]+TotalAutoPolMel[5])
  denoms<- rbind(vAutoATGCSim, vAutoATGCMel)
  colnames(denoms)=c("A","T","G","C")
  rownames(denoms)=c("Sim","Mel")
  return(denoms)
}

######### Bootstrap build

bootstrap_divoverpol_6par<- function(df, BS=1000) {
  respol<- matrix(0,BS,9)
  colnames(respol)<- c("a","f", "b","d","c","e", "GCneutral", "GCchanging","All")

  resdiv<- matrix(0,BS,9)
  colnames(resdiv)<- c("a","f", "b","d","c","e","GCneutral", "GCchanging","All")
  
  divoverpol<- matrix(0,BS,9)
  colnames(divoverpol)<- c("a","f", "b","d","c","e", "GCneutral", "GCchanging","All")
  
  for (b in 1:BS) {
    new_samp<- sample_5SI(df)
    
    ### Polymorphism 
    spectra<- spectra6_build(new_samp, M=69)
    
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
    
    respol[b,1]<- L_AT/denominatorAT/HM # a 
    respol[b,2]<- L_CG/denominatorCG/HM # f
    
    respol[b,3]<- bPrime*L_overATCG/denominatorCG/2/HM  # b: CT and GA (Transition)
    respol[b,4]<- (1-bPrime)*L_overATCG/denominatorCG/2/HM # d: GT and CA
    
    respol[b,5]<- (1-ePrime)*L_overATCG/denominatorAT/2/HM # c: AG and TC (Transition)
    respol[b,6]<- ePrime*L_overATCG/denominatorAT/2/HM # e: AC and TG

    respol[b,7]<- sum(respol[b,1:2]) # GCneutral
    respol[b,8]<- sum(respol[b,3:6]) # GCchanging
    respol[b,9]<- sum(respol[b,1:6]) # All
    
    ### Divergence 
    
    mat10by10<- mat10by10_reduceddat_create(new_samp)
    
    denoms<- denoms_create(mat10by10)
    summat<- premdf_create(mat10by10)
    
    denominatorAT_mel<- denoms[2,1]+denoms[2,2]
    denominatorCG_mel<- denoms[2,3]+denoms[2,4]
    
    # a
    resdiv[b,1]=summat[1,1]/denominatorAT_mel-summat[7,1]/2/denominatorAT_mel
    # f
    resdiv[b,2]=summat[1,2]/denominatorCG_mel-summat[7,2]/2/denominatorCG_mel
    # b
    resdiv[b,3]=(summat[1,4]+summat[1,5])/denominatorCG_mel-(summat[7,4]+summat[7,5])/2/denominatorCG_mel
    # d 
    resdiv[b,4]=(summat[1,3]+summat[1,6])/denominatorCG_mel-(summat[7,3]+summat[7,6])/2/denominatorCG_mel
    # c
    resdiv[b,5]=(summat[1,4]+summat[1,5])/denominatorAT_mel-(summat[7,4]+summat[7,5])/2/denominatorAT_mel
    # e
    resdiv[b,6]=(summat[1,3]+summat[1,6])/denominatorAT_mel-(summat[7,3]+summat[7,6])/2/denominatorAT_mel
    
    resdiv[b,7]<- sum(resdiv[b,1:2]) # GCneutral
    resdiv[b,8]<- sum(resdiv[b,3:6]) # GCchanging
    resdiv[b,9]<- sum(resdiv[b,1:6]) # All
    
    ## Ratio 
    divoverpol[b,]<- resdiv[b,]/respol[b,]
  }
  return(divoverpol)
}


