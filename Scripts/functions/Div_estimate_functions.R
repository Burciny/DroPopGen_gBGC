setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#### FOR DIVERGENCE 
################################ ################ ################ ################ ################ 
############### Create 10x10 matrix Simulans-Melanogaster data ################
################################################################################################
## input df: Population count data (i.e., Population_data/analyzedSIsites_mel69_autosome)

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


## Preliminary data frames
## input df: 10x10 matrix

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

##### Divergence estimates

### Minus shared polymorphism
## summat: output of premdf_create function ; denom: output of denoms_create function

alldiv_est_sharedpol<- function(summat, denom) {
  denominatorAT_sim<- denom[1,1]+denom[1,2]
  denominatorCG_sim<- denom[1,3]+denom[1,4]
  
  denominatorAT_mel<- denom[2,1]+denom[2,2]
  denominatorCG_mel<- denom[2,3]+denom[2,4]
  
  results<- matrix(0,2,6)
  colnames(results)<- c("a","f", "b","d","c","e")
  rownames(results)<- c("Sim","Mel")
  # Simulans
  # a
  results[1,1]=summat[1,1]/denominatorAT_sim-summat[6,1]/2/denominatorAT_sim
  # f
  results[1,2]=summat[1,2]/denominatorCG_sim-summat[6,2]/2/denominatorCG_sim
  # b
  results[1,3]=(summat[1,4]+summat[1,5])/denominatorCG_sim-(summat[6,4]+summat[6,5])/2/denominatorCG_sim
  # d   
  results[1,4]=(summat[1,3]+summat[1,6])/denominatorCG_sim-(summat[6,3]+summat[6,6])/2/denominatorCG_sim
  # c
  results[1,5]=(summat[1,4]+summat[1,5])/denominatorAT_sim-(summat[6,4]+summat[6,5])/2/denominatorAT_sim
  # e
  results[1,6]=(summat[1,3]+summat[1,6])/denominatorAT_sim-(summat[6,3]+summat[6,6])/2/denominatorAT_sim
  
  # Melanogaster
  results[2,1]=summat[1,1]/denominatorAT_mel-summat[7,1]/2/denominatorAT_mel
  results[2,2]=summat[1,2]/denominatorCG_mel-summat[7,2]/2/denominatorCG_mel
  results[2,3]=(summat[1,4]+summat[1,5])/denominatorCG_mel-(summat[7,4]+summat[7,5])/2/denominatorCG_mel
  results[2,4]=(summat[1,3]+summat[1,6])/denominatorCG_mel-(summat[7,3]+summat[7,6])/2/denominatorCG_mel
  results[2,5]=(summat[1,4]+summat[1,5])/denominatorAT_mel-(summat[7,4]+summat[7,5])/2/denominatorAT_mel
  results[2,6]=(summat[1,3]+summat[1,6])/denominatorAT_mel-(summat[7,3]+summat[7,6])/2/denominatorAT_mel
  
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


## Bootstrap build 

bs_div_6par_sharedpol<- function(df, BS=1000) {
  results<- matrix(0,BS,6)
  colnames(results)<- c("a","f", "b","d","c","e")
  
  for (b in 1:BS) {
    new_samp<- sample_5SI(df)
    mat10by10<- mat10by10_reduceddat_create(new_samp)
    
    denoms<- denoms_create(mat10by10)
    summat<- premdf_create(mat10by10)
    
    denominatorAT_mel<- denoms[2,1]+denoms[2,2]
    denominatorCG_mel<- denoms[2,3]+denoms[2,4]
    
    # Melanogaster
    # a
    results[b,1]=summat[1,1]/denominatorAT_mel-summat[7,1]/2/denominatorAT_mel
    # f
    results[b,2]=summat[1,2]/denominatorCG_mel-summat[7,2]/2/denominatorCG_mel
    # b
    results[b,3]=(summat[1,4]+summat[1,5])/denominatorCG_mel-(summat[7,4]+summat[7,5])/2/denominatorCG_mel
    # d 
    results[b,4]=(summat[1,3]+summat[1,6])/denominatorCG_mel-(summat[7,3]+summat[7,6])/2/denominatorCG_mel
    # c
    results[b,5]=(summat[1,4]+summat[1,5])/denominatorAT_mel-(summat[7,4]+summat[7,5])/2/denominatorAT_mel
    # e
    results[b,6]=(summat[1,3]+summat[1,6])/denominatorAT_mel-(summat[7,3]+summat[7,6])/2/denominatorAT_mel
    
  }
  return(results)
}



### Example use;
analyzedSIsites_mel69_autosome <- read_delim("../../Population_data/Melanogaster/analyzedSIsites_mel69_autosome", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)


Auto5SI_compMelvsSim<- mat10by10_reduceddat_create(analyzedSIsites_mel69_autosome)

denomATGC_auto<- denoms_create(Auto5SI_compMelvsSim)
mats_auto<- premdf_create(Auto5SI_compMelvsSim)

allmuts_auto_div<- alldiv_est_sharedpol(mats_auto, denomATGC_auto)

