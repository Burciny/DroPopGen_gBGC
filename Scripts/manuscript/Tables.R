setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(readr)
################ ################ ################ ################ ################ 
################ Table 1 ################ ################ ################ ################ 


ExpHeterozygosity_all <- read_delim("../../Intermediate_data/ExpHeterozygosity_all", 
                                                  "\t", escape_double = FALSE, trim_ws = TRUE)

ExpHeterozygosity_all


################ ################ ################ ################ ################ 
################ Table 2 ################ ################ ################ ################ 

divergenceoverpolymorphism_ratios <- read_delim("../../Intermediate_data/DivoverPol/divergenceoverpolymorphism_ratios_Mel69", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE)
divergenceoverpolymorphism_ratios[divergenceoverpolymorphism_ratios$Mutationtype=="GCchanging",]
divergenceoverpolymorphism_ratios[divergenceoverpolymorphism_ratios$Mutationtype=="GCneutral",]


################ ################ ################ ################ ################ 
################ Table 3 ################ ################ ################ ################ 
XvsAutosome_pol6par_alongchr <- read_delim("../../Intermediate_data/Pol_estimates/XvsAutosome_pol6par_alongchr", 
                                            "\t", escape_double = FALSE, trim_ws = TRUE)

XvsAutosome_mel_div6par_alongchr <- read_delim("../../Intermediate_data/Div_estimates/XvsAutosome_mel_div6par_alongchr", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)


XvsAutosome_pol6par_alongchr
XvsAutosome_mel_div6par_alongchr



################ ################ ################ ################ ################ 
################ Table 4 ################ ################ ################ ################ 
###
mel69_commonreduced5p_spectras_autosome <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_autosome", 
                                                      "\t", escape_double = FALSE, trim_ws = TRUE)
mel69_commonreduced5p_spectras_X <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_X", 
                                                      "\t", escape_double = FALSE, trim_ws = TRUE)

pool_whole<- mel69_commonreduced5p_spectras_autosome+mel69_commonreduced5p_spectras_X
###
mel69_commonreduced5p_spectras_autosome_middle <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_autosome_middle", 
                                                      "\t", escape_double = FALSE, trim_ws = TRUE)
mel69_commonreduced5p_spectras_X_middle <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_X_middle", 
                                                      "\t", escape_double = FALSE, trim_ws = TRUE)

pool_center<- mel69_commonreduced5p_spectras_autosome_middle+mel69_commonreduced5p_spectras_X_middle

###
mel69_commonreduced5p_spectras_autosome_ends <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_autosome_ends", 
                                                           "\t", escape_double = FALSE, trim_ws = TRUE)
mel69_commonreduced5p_spectras_X_ends <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_X_ends", 
                                                    "\t", escape_double = FALSE, trim_ws = TRUE)

pool_ends<- mel69_commonreduced5p_spectras_autosome_ends+mel69_commonreduced5p_spectras_X_ends



B_Mel69_5SI_likCI <- read_delim("../../Intermediate_data/B_estimates/B_Mel69_5SI_likCI", 
                                          "\t", escape_double = FALSE, trim_ws = TRUE)

### Autosome and X 
B_Mel69_5SI_likCI[B_Mel69_5SI_likCI$Spectra=="GCchanging",]

source("../functions/gBGC_estimate_functions.R")

## Pool 
# WChr 
B_pool_WChr<- optim(1, logl,data=pool_whole$Asymm,c=0, method = "BFGS")$par
B_pool_WChr_ll1<- -optim(1, logl,data=pool_whole$Asymm,c=0, method = "BFGS")$value
B_pool_WChr_ll0<- -logl(gamma=0, data=pool_whole$Asymm,c=0)

likCI_WChr<-findUL(fun=myfun,MLE = B_pool_WChr, dataf=pool_whole$Asymm)
likCI_WChr$Low
likCI_WChr$Up
#### CChr
B_pool_CChr<- optim(1, logl,data=pool_center$Asymm,c=0, method = "BFGS")$par
B_pool_CChr_ll1<- -optim(1, logl,data=pool_center$Asymm,c=0, method = "BFGS")$value
B_pool_CChr_ll0<- -logl(gamma=0, data=pool_center$Asymm,c=0)

likCI_CChr<-findUL(fun=myfun,MLE = B_pool_CChr, dataf=pool_center$Asymm)
likCI_CChr$Low
likCI_CChr$Up

### PChr
B_pool_PChr<- optim(1, logl,data=pool_ends$Asymm,c=0, method = "BFGS")$par
B_pool_PChr_ll1<- -optim(1, logl,data=pool_ends$Asymm,c=0, method = "BFGS")$value
B_pool_PChr_ll0<- -logl(gamma=0, data=pool_ends$Asymm,c=0)

likCI_PChr<-findUL(fun=myfun,MLE = B_pool_PChr, dataf=pool_ends$Asymm)
likCI_PChr$Low
likCI_PChr$Up

################ ################ ################ ################ ################ 
################ Table 5 ################ ################ ################ ################ 

### autosome and X 
B_Mel69_5SI_likCI[B_Mel69_5SI_likCI$Spectra=="Ts",]
B_Mel69_5SI_likCI[B_Mel69_5SI_likCI$Spectra=="Tv",]

### Pool 
# WChr - Transitions
B_pool_WChr_Ts<- optim(1, logl,data=pool_whole$Ts,c=0, method = "BFGS")$par
B_pool_WChr_Ts_ll1<- -optim(1, logl,data=pool_whole$Ts,c=0, method = "BFGS")$value
B_pool_WChr_Ts_ll0<- -logl(gamma=0, data=pool_whole$Ts,c=0)


likCI_WChr_Ts<-findUL(fun=myfun,MLE = B_pool_WChr_Ts, dataf=pool_whole$Ts)
likCI_WChr_Ts$Low
likCI_WChr_Ts$Up

### WChr- Transversions
B_pool_WChr_Tv<- optim(1, logl,data=pool_whole$Tv,c=0, method = "BFGS")$par
B_pool_WChr_Tv_ll1<- -optim(1, logl,data=pool_whole$Tv,c=0, method = "BFGS")$value
B_pool_WChr_Tv_ll0<- -logl(gamma=0, data=pool_whole$Tv,c=0)

likCI_WChr_Tv<-findUL(fun=myfun,MLE = B_pool_WChr_Tv, dataf=pool_whole$Tv)
likCI_WChr_Tv$Low
likCI_WChr_Tv$Up

### CChr- Transitions 

B_pool_CChr_Ts<- optim(1, logl,data=pool_center$Ts,c=0, method = "BFGS")$par
B_pool_CChr_Ts_ll1<- -optim(1, logl,data=pool_center$Ts,c=0, method = "BFGS")$value
B_pool_CChr_Ts_ll0<- -logl(gamma=0, data=pool_center$Ts,c=0)

likCI_CChr_Ts<-findUL(fun=myfun,MLE = B_pool_CChr_Ts, dataf=pool_center$Ts)
likCI_CChr_Ts$Low
likCI_CChr_Ts$Up

### CChr- Transversion
B_pool_CChr_Tv<- optim(1, logl,data=pool_center$Tv,c=0, method = "BFGS")$par
B_pool_CChr_Tv_ll1<- -optim(1, logl,data=pool_center$Tv,c=0, method = "BFGS")$value
B_pool_CChr_Tv_ll0<- -logl(gamma=0, data=pool_center$Tv,c=0)

likCI_CChr_Tv<-findUL(fun=myfun,MLE = B_pool_CChr_Tv, dataf=pool_center$Tv)
likCI_CChr_Tv$Low
likCI_CChr_Tv$Up

#### PChr - Transitions
B_pool_PChr_Ts<- optim(1, logl,data=pool_ends$Ts,c=0, method = "BFGS")$par
B_pool_PChr_Ts_ll1<- -optim(1, logl,data=pool_ends$Ts,c=0, method = "BFGS")$value
B_pool_PChr_Ts_ll0<- -logl(gamma=0, data=pool_ends$Ts,c=0)

likCI_PChr_Ts<-findUL(fun=myfun,MLE = B_pool_PChr_Ts, dataf=pool_ends$Ts)
likCI_PChr_Ts$Low
likCI_PChr_Ts$Up

### PChr- Transversions

B_pool_PChr_Tv<- optim(1, logl,data=pool_ends$Tv,c=0, method = "BFGS")$par
B_pool_ll1<- -optim(1, logl,data=pool_ends$Tv,c=0, method = "BFGS")$value
B_pool_ll0<- -logl(gamma=0, data=pool_ends$Tv,c=0)

likCI_PChr_Tv<-findUL(fun=myfun,MLE = B_pool_PChr_Tv, dataf=pool_ends$Tv)
likCI_PChr_Tv$Low
likCI_PChr_Tv$Up

################ ################ ################ ################ ################ 
################ Table 6 ################ ################ ################ ################ 
###
MDSim21_commonreduced5p_spectras_autosome_middle <- read_delim("../../Intermediate_data/Spectras/Simulans/MDSim21_commonreduced5p_spectras_autosome_middle", 
                                                               "\t", escape_double = FALSE, trim_ws = TRUE)
MDSim21_commonreduced5p_spectras_X_middle <- read_delim("../../Intermediate_data/Spectras/Simulans/MDSim21_commonreduced5p_spectras_X_middle", 
                                                        "\t", escape_double = FALSE, trim_ws = TRUE)

pool_center<- MDSim21_commonreduced5p_spectras_autosome_middle+MDSim21_commonreduced5p_spectras_X_middle

B_MDsim_all <- read_delim("../../Intermediate_data/B_estimates/Simulans/B_MDsim_all", 
                                              "\t", escape_double = FALSE, trim_ws = TRUE)

### Autosome and X 
B_MDsim_all[B_MDsim_all$Part=="center" & B_MDsim_all$Spectra=="GCchanging",]
B_MDsim_all[B_MDsim_all$Part=="center" & B_MDsim_all$Spectra=="Ts",]
B_MDsim_all[B_MDsim_all$Part=="center" & B_MDsim_all$Spectra=="Tv",]

source("../functions/gBGC_estimate_functions.R")

## Pool
# 
B_pool_CChr<- optim(1, logl,data=pool_center$Asymm,c=0, method = "BFGS")$par
B_pool_CChr_ll1<- -optim(1, logl,data=pool_center$Asymm,c=0, method = "BFGS")$value
B_pool_CChr_ll0<- -logl(gamma=0, data=pool_center$Asymm,c=0)

likCI_CChr<-findUL(fun=myfun,MLE = B_pool_CChr, dataf=pool_center$Asymm)
likCI_CChr$Low
likCI_CChr$Up

## Transitions 
B_pool_CChr_Ts<- optim(1, logl,data=pool_center$Ts,c=0, method = "BFGS")$par
B_pool_CChr_Ts_ll1<- -optim(1, logl,data=pool_center$Ts,c=0, method = "BFGS")$value
B_pool_CChr_Ts_ll0<- -logl(gamma=0, data=pool_center$Ts,c=0)

likCI_CChr_Ts<-findUL(fun=myfun,MLE = B_pool_CChr_Ts, dataf=pool_center$Ts)
likCI_CChr_Ts$Low
likCI_CChr_Ts$Up

### Transversion
B_pool_CChr_Tv<- optim(1, logl,data=pool_center$Tv,c=0, method = "BFGS")$par
B_pool_CChr_Tv_ll1<- -optim(1, logl,data=pool_center$Tv,c=0, method = "BFGS")$value
B_pool_CChr_Tv_ll0<- -logl(gamma=0, data=pool_center$Tv,c=0)

likCI_CChr_Tv<-findUL(fun=myfun,MLE = B_pool_CChr_Tv, dataf=pool_center$Tv)
likCI_CChr_Tv$Low
likCI_CChr_Tv$Up

################ ################ ################ ################ ################ 
################ SUPPLLEMENT ################ ################ ################ ################ 

################ ################ ################ ################ ################ 
################ Table S1 ################ ################ ################ ################ 

### WChr 
analyzedSIsites_mel69_autosome_WChr_highr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_autosome_WChr_highr", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_mel69_autosome_WChr_modhighr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_autosome_WChr_modhighr", 
                                                        delim = "\t", escape_double = FALSE, 
                                                        trim_ws = TRUE)
analyzedSIsites_mel69_autosome_WChr_modlowr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_autosome_WChr_modlowr", 
                                                           delim = "\t", escape_double = FALSE, 
                                                           trim_ws = TRUE)
analyzedSIsites_mel69_autosome_WChr_lowr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_autosome_WChr_lowr", 
                                                           delim = "\t", escape_double = FALSE, 
                                                           trim_ws = TRUE)

analyzedSIsites_mel69_X_WChr_highr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_X_WChr_highr", 
                                                        delim = "\t", escape_double = FALSE, 
                                                        trim_ws = TRUE)
analyzedSIsites_mel69_X_WChr_modhighr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_X_WChr_modhighr", 
                                                           delim = "\t", escape_double = FALSE, 
                                                           trim_ws = TRUE)
analyzedSIsites_mel69_X_WChr_modlowr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_X_WChr_modlowr", 
                                                          delim = "\t", escape_double = FALSE, 
                                                          trim_ws = TRUE)
analyzedSIsites_mel69_X_WChr_lowr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_X_WChr_lowr", 
                                                       delim = "\t", escape_double = FALSE, 
                                                       trim_ws = TRUE)


###### CChr

analyzedSIsites_mel69_autosome_CChr_highr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_autosome_CChr_highr", 
                                                        delim = "\t", escape_double = FALSE, 
                                                        trim_ws = TRUE)
analyzedSIsites_mel69_autosome_CChr_modhighr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_autosome_CChr_modhighr", 
                                                           delim = "\t", escape_double = FALSE, 
                                                           trim_ws = TRUE)
analyzedSIsites_mel69_autosome_CChr_modlowr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_autosome_CChr_modlowr", 
                                                          delim = "\t", escape_double = FALSE, 
                                                          trim_ws = TRUE)
analyzedSIsites_mel69_autosome_CChr_lowr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_autosome_CChr_lowr", 
                                                       delim = "\t", escape_double = FALSE, 
                                                       trim_ws = TRUE)

analyzedSIsites_mel69_X_CChr_highr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_X_CChr_highr", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_mel69_X_CChr_modhighr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_X_CChr_modhighr", 
                                                    delim = "\t", escape_double = FALSE, 
                                                    trim_ws = TRUE)
analyzedSIsites_mel69_X_CChr_modlowr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_X_CChr_modlowr", 
                                                   delim = "\t", escape_double = FALSE, 
                                                   trim_ws = TRUE)
analyzedSIsites_mel69_X_CChr_lowr <- read_delim("../../Population_data/Melanogaster/Rbins/WChr/analyzedSIsites_mel69_X_CChr_lowr", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE)


################ ################ ################ ################ ################ 
################ Table S2 ################ ################ ################ ################ 

analyzedSIsites_mel69_autosome_GC1 <- read_delim("../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_autosome_GC1", 
                                             delim = "\t", escape_double = FALSE, 
                                             trim_ws = TRUE)
analyzedSIsites_mel69_autosome_GC2 <- read_delim("../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_autosome_GC2", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_mel69_autosome_GC3 <- read_delim("../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_autosome_GC3", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_mel69_autosome_GC4 <- read_delim("../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_autosome_GC4", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_mel69_autosome_GC5 <- read_delim("../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_autosome_GC5", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)

analyzedSIsites_mel69_autosome_GC1$bins[1]
analyzedSIsites_mel69_autosome_GC2$bins[1]
analyzedSIsites_mel69_autosome_GC3$bins[1]
analyzedSIsites_mel69_autosome_GC4$bins[1]
analyzedSIsites_mel69_autosome_GC5$bins[1]

analyzedSIsites_mel69_X_GC1 <- read_delim("../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_X_GC1", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_mel69_X_GC2 <- read_delim("../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_X_GC2", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
analyzedSIsites_mel69_X_GC3 <- read_delim("../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_X_GC3", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
analyzedSIsites_mel69_X_GC4 <- read_delim("../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_X_GC4", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
analyzedSIsites_mel69_X_GC5 <- read_delim("../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_X_GC5", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

analyzedSIsites_mel69_X_GC1$bins[1]
analyzedSIsites_mel69_X_GC2$bins[1]
analyzedSIsites_mel69_X_GC3$bins[1]
analyzedSIsites_mel69_X_GC4$bins[1]
analyzedSIsites_mel69_X_GC5$bins[1]

################ ################ ################ ################ ################ 
################ Table S3 ################ ################ ################ ################ 

analyzedSIsites_MDsim_autosome_GC1 <- read_delim("../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_autosome_GC1", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_MDsim_autosome_GC2 <- read_delim("../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_autosome_GC2", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_MDsim_autosome_GC3 <- read_delim("../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_autosome_GC3", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_MDsim_autosome_GC4 <- read_delim("../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_autosome_GC4", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_MDsim_autosome_GC5 <- read_delim("../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_autosome_GC5", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_MDsim_autosome_GC1$bins[1]
analyzedSIsites_MDsim_autosome_GC2$bins[1]
analyzedSIsites_MDsim_autosome_GC3$bins[1]
analyzedSIsites_MDsim_autosome_GC4$bins[1]
analyzedSIsites_MDsim_autosome_GC5$bins[1]

analyzedSIsites_MDsim_X_GC1 <- read_delim("../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_X_GC1", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
analyzedSIsites_MDsim_X_GC2 <- read_delim("../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_X_GC2", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
analyzedSIsites_MDsim_X_GC3 <- read_delim("../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_X_GC3", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
analyzedSIsites_MDsim_X_GC4 <- read_delim("../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_X_GC4", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
analyzedSIsites_MDsim_X_GC5 <- read_delim("../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_X_GC5", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

analyzedSIsites_MDsim_X_GC1$bins[1]
analyzedSIsites_MDsim_X_GC2$bins[1]
analyzedSIsites_MDsim_X_GC3$bins[1]
analyzedSIsites_MDsim_X_GC4$bins[1]
analyzedSIsites_MDsim_X_GC5$bins[1]

################ ################ ################ ################ ################ 
################ Table 4- Juraj  ################ ################ ################ ################ 
###
###
mel69_commonreduced5p_spectras_autosome <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_autosome", 
                                                      "\t", escape_double = FALSE, trim_ws = TRUE)
mel69_commonreduced5p_spectras_X <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_X", 
                                               "\t", escape_double = FALSE, trim_ws = TRUE)

pool_whole<- mel69_commonreduced5p_spectras_autosome+mel69_commonreduced5p_spectras_X
###
mel69_commonreduced5p_spectras_autosome_middle <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_autosome_middle", 
                                                             "\t", escape_double = FALSE, trim_ws = TRUE)
mel69_commonreduced5p_spectras_X_middle <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_X_middle", 
                                                      "\t", escape_double = FALSE, trim_ws = TRUE)

pool_center<- mel69_commonreduced5p_spectras_autosome_middle+mel69_commonreduced5p_spectras_X_middle

###
mel69_commonreduced5p_spectras_autosome_ends <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_autosome_ends", 
                                                           "\t", escape_double = FALSE, trim_ws = TRUE)
mel69_commonreduced5p_spectras_X_ends <- read_delim("../../Intermediate_data/Spectras/Melanogaster/mel69_commonreduced5p_spectras_X_ends", 
                                                    "\t", escape_double = FALSE, trim_ws = TRUE)

pool_ends<- mel69_commonreduced5p_spectras_autosome_ends+mel69_commonreduced5p_spectras_X_ends



B_Mel69_5SI_wry <- read_delim("../../Intermediate_data/B_estimates/with_ry/B_Mel69_5SI_wry", 
                                              "\t", escape_double = FALSE, trim_ws = TRUE)

### Autosome and X 
B_Mel69_5SI_wry[B_Mel69_5SI_wry$Spectra=="GCchanging",]

source("../functions/gBGC_estimate_functions.R")


## Pool 
# WChr 
B_pool_WChr<- optim(1, logl_wnuissance,data=pool_whole$Asymm,c=0,Symm=pool_whole$Symm, method = "BFGS")$par
B_pool_WChr_ll1<- -optim(1, logl_wnuissance,data=pool_whole$Asymm,c=0,Symm=pool_whole$Symm, method = "BFGS")$value
B_pool_WChr_ll0<- -logl_wnuissance(gamma=0, data=pool_whole$Asymm,c=0, Symm=pool_whole$Symm)

pchisq(2*(B_pool_WChr_ll1-B_pool_WChr_ll0), df=1, lower.tail = FALSE)

likCI_WChr<-findUL(fun=myfun_wry,MLE = B_pool_WChr, dataf=pool_whole$Asymm,Symmf=pool_whole$Symm)
likCI_WChr$Low
likCI_WChr$Up

#### CChr
B_pool_CChr<- optim(1, logl_wnuissance,data=pool_center$Asymm,c=0, Symm=pool_center$Symm, method = "BFGS")$par
B_pool_CChr_ll1<- -optim(1, logl_wnuissance,data=pool_center$Asymm,c=0,Symm=pool_center$Symm, method = "BFGS")$value
B_pool_CChr_ll0<- -logl_wnuissance(gamma=0, data=pool_center$Asymm,c=0, Symm=pool_center$Symm)

pchisq(2*(B_pool_CChr_ll1-B_pool_CChr_ll0), df=1, lower.tail = FALSE)

likCI_CChr<-findUL(fun=myfun_wry,MLE = B_pool_CChr, dataf=pool_center$Asymm, Symmf=pool_center$Symm)
likCI_CChr$Low
likCI_CChr$Up

### PChr
B_pool_PChr<- optim(1, logl_wnuissance,data=pool_ends$Asymm,c=0, Symm=pool_ends$Symm,method = "BFGS")$par
B_pool_PChr_ll1<- -optim(1, logl_wnuissance,data=pool_ends$Asymm,c=0,Symm=pool_ends$Symm, method = "BFGS")$value
B_pool_PChr_ll0<- -logl_wnuissance(gamma=0, data=pool_ends$Asymm,c=0, Symm=pool_ends$Symm)

pchisq(2*(B_pool_PChr_ll1-B_pool_PChr_ll0), df=1, lower.tail = FALSE)

likCI_CChr<-findUL(fun=myfun_wry,MLE = B_pool_PChr, dataf=pool_ends$Asymm, Symmf=pool_ends$Symm)
likCI_CChr$Low
likCI_CChr$Up

################ ################ ################ ################ ################ 
################ Table 5- Juraj  ################ ################ ################ ################ 
###


### autosome and X 
B_Mel69_5SI_wry[B_Mel69_5SI_wry$Spectra=="Ts",]
B_Mel69_5SI_wry[B_Mel69_5SI_wry$Spectra=="Tv",]

### Pool 
# WChr - Transitions
B_pool_WChr_Ts<- optim(1, logl_wnuissance,data=pool_whole$Ts,c=0, Symm=pool_whole$Symm, method = "BFGS")$par
B_pool_WChr_Ts_ll1<- -optim(1, logl_wnuissance,data=pool_whole$Ts,c=0, Symm=pool_whole$Symm, method = "BFGS")$value
B_pool_WChr_Ts_ll0<- -logl_wnuissance(gamma=0, data=pool_whole$Ts,c=0, Symm=pool_whole$Symm)

likCI_WChr_Ts<-findUL(fun=myfun_wry,MLE = B_pool_WChr_Ts, dataf=pool_whole$Ts, Symmf=pool_whole$Symm)
likCI_WChr_Ts$Low
likCI_WChr_Ts$Up

### WChr- Transversions
B_pool_WChr_Tv<- optim(1, logl_wnuissance,data=pool_whole$Tv,c=0, Symm=pool_whole$Symm, method = "BFGS")$par
B_pool_WChr_Tv_ll1<- -optim(1, logl_wnuissance,data=pool_whole$Tv,c=0, Symm=pool_whole$Symm, method = "BFGS")$value
B_pool_WChr_Tv_ll0<- -logl_wnuissance(gamma=0, data=pool_whole$Tv,c=0, Symm=pool_whole$Symm)

likCI_WChr_Tv<-findUL(fun=myfun_wry,MLE = B_pool_WChr_Tv, dataf=pool_whole$Tv, Symmf=pool_whole$Symm)
likCI_WChr_Tv$Low
likCI_WChr_Tv$Up

### CChr- Transitions 

B_pool_CChr_Ts<- optim(1, logl_wnuissance,data=pool_center$Ts,c=0, Symm=pool_center$Symm, method = "BFGS")$par
B_pool_CChr_Ts_ll1<- -optim(1, logl_wnuissance,data=pool_center$Ts,c=0, Symm=pool_center$Symm, method = "BFGS")$value
B_pool_CChr_Ts_ll0<- -logl_wnuissance(gamma=0, data=pool_center$Ts,c=0, Symm=pool_center$Symm)

likCI_CChr_Ts<-findUL(fun=myfun_wry,MLE = B_pool_CChr_Ts, dataf=pool_center$Ts, Symmf=pool_center$Symm)
likCI_CChr_Ts$Low
likCI_CChr_Ts$Up


### CChr- Transversion
B_pool_CChr_Tv<- optim(1, logl_wnuissance,data=pool_center$Tv,c=0, Symm=pool_center$Symm, method = "BFGS")$par
B_pool_CChr_Tv_ll1<- -optim(1, logl_wnuissance,data=pool_center$Tv,c=0, Symm=pool_center$Symm, method = "BFGS")$value
B_pool_CChr_Tv_ll0<- -logl_wnuissance(gamma=0, data=pool_center$Tv,c=0, Symm=pool_center$Symm)

likCI_CChr_Tv<-findUL(fun=myfun_wry,MLE = B_pool_CChr_Tv, dataf=pool_center$Tv, Symmf=pool_center$Symm)
likCI_CChr_Tv$Low
likCI_CChr_Tv$Up


#### PChr - Transitions
B_pool_PChr_Ts<- optim(1, logl_wnuissance,data=pool_ends$Ts,c=0, Symm=pool_ends$Symm, method = "BFGS")$par
B_pool_PChr_Ts_ll1<- -optim(1, logl_wnuissance,data=pool_ends$Ts,c=0, Symm=pool_ends$Symm, method = "BFGS")$value
B_pool_PChr_Ts_ll0<- -logl_wnuissance(gamma=0, data=pool_ends$Ts,c=0, Symm=pool_ends$Symm)

likCI_PChr_Ts<-findUL(fun=myfun_wry,MLE = B_pool_PChr_Ts, dataf=pool_ends$Ts, Symmf=pool_ends$Symm)
likCI_PChr_Ts$Low
likCI_PChr_Ts$Up

### PChr- Transversions

B_pool_PChr_Tv<- optim(1, logl_wnuissance,data=pool_ends$Tv,c=0, Symm=pool_ends$Symm, method = "BFGS")$par
B_pool_ll1<- -optim(1, logl_wnuissance,data=pool_ends$Tv,c=0, Symm=pool_ends$Symm, method = "BFGS")$value
B_pool_ll0<- -logl_wnuissance(gamma=0, data=pool_ends$Tv,c=0, Symm=pool_ends$Symm)

likCI_PChr_Tv<-findUL(fun=myfun_wry,MLE = B_pool_PChr_Tv, dataf=pool_ends$Tv, Symmf=pool_ends$Symm)
likCI_PChr_Tv$Low
likCI_PChr_Tv$Up


################ ################ ################ ################ ################ 
################ Table 6- Juraj ################ ################ ################ ################ 
###
MDSim21_commonreduced5p_spectras_autosome_middle <- read_delim("../../Intermediate_data/Spectras/Simulans/MDSim21_commonreduced5p_spectras_autosome_middle", 
                                                               "\t", escape_double = FALSE, trim_ws = TRUE)
MDSim21_commonreduced5p_spectras_X_middle <- read_delim("../../Intermediate_data/Spectras/Simulans//MDSim21_commonreduced5p_spectras_X_middle", 
                                                        "\t", escape_double = FALSE, trim_ws = TRUE)

pool_center<- MDSim21_commonreduced5p_spectras_autosome_middle+MDSim21_commonreduced5p_spectras_X_middle

B_MDsim_all_wry <- read_delim("../../Intermediate_data/B_estimates/Simulans/with_ry/B_MDsim_all_wry", 
                                                        "\t", escape_double = FALSE, trim_ws = TRUE)

### Autosome and X 
B_MDsim_all_wry[B_MDsim_all_wry$Part=="CChr" & B_MDsim_all_wry$Spectra=="GCchanging",]
B_MDsim_all_wry[B_MDsim_all_wry$Part=="CChr" & B_MDsim_all_wry$Spectra=="Ts",]
B_MDsim_all_wry[B_MDsim_all_wry$Part=="CChr" & B_MDsim_all_wry$Spectra=="Tv",]

source("../functions/gBGC_estimate_functions.R")

## Pool
# 
B_pool_CChr<- optim(1, logl_wnuissance,data=pool_center$Asymm,c=0,Symm=pool_center$Symm, method = "BFGS")$par
B_pool_CChr_ll1<- -optim(1, logl_wnuissance,data=pool_center$Asymm,c=0,Symm=pool_center$Symm, method = "BFGS")$value
B_pool_CChr_ll0<- -logl_wnuissance(gamma=0, data=pool_center$Asymm,c=0,Symm=pool_center$Symm)

likCI_CChr<-findUL(fun=myfun_wry,MLE = B_pool_CChr, dataf=pool_center$Asymm, Symmf=pool_center$Symm)
likCI_CChr$Low
likCI_CChr$Up

## Transitions 
B_pool_CChr_Ts<- optim(1, logl_wnuissance,data=pool_center$Ts,c=0,Symm=pool_center$Symm, method = "BFGS")$par
B_pool_CChr_Ts_ll1<- -optim(1, logl_wnuissance,data=pool_center$Ts,c=0,Symm=pool_center$Symm, method = "BFGS")$value
B_pool_CChr_Ts_ll0<- -logl_wnuissance(gamma=0, data=pool_center$Ts,c=0,Symm=pool_center$Symm)

likCI_CChr_Ts<-findUL(fun=myfun_wry,MLE = B_pool_CChr_Ts, dataf=pool_center$Ts, Symmf=pool_center$Symm)
likCI_CChr_Ts$Low
likCI_CChr_Ts$Up

### Transversion
B_pool_CChr_Tv<- optim(1, logl_wnuissance,data=pool_center$Tv,c=0, Symm=pool_center$Symm, method = "BFGS")$par
B_pool_CChr_Tv_ll1<- -optim(1, logl_wnuissance,data=pool_center$Tv,c=0, Symm=pool_center$Symm,method = "BFGS")$value
B_pool_CChr_Tv_ll0<- -logl_wnuissance(gamma=0, data=pool_center$Tv,c=0,Symm=pool_center$Symm)

likCI_CChr_Tv<-findUL(fun=myfun_wry,MLE = B_pool_CChr_Tv, dataf=pool_center$Tv, Symmf=pool_center$Symm)
likCI_CChr_Tv$Low
likCI_CChr_Tv$Up


