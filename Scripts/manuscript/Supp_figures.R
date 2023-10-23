setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(readr)
library(gtools)
library(tidyr)
library(dplyr)
library(ggplot2)

###################################################################################################
########################### FIGURE S1 ###############################################################

divergenceoverpolymorphism_ratios_Mel69 <- read_delim("../../Intermediate_data/DivoverPol/divergenceoverpolymorphism_ratios_Mel69", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE)

divergenceoverpolymorphism_ratios_Mel69$Part<- factor(divergenceoverpolymorphism_ratios_Mel69$Part, levels = c("whole","center","ends"))
divergenceoverpolymorphism_ratios_Mel69$Mutationtype<- factor(divergenceoverpolymorphism_ratios_Mel69$Mutationtype, levels = c("a","f","b","c","d","e","GCneutral","GCchanging", "All"))


center<- divergenceoverpolymorphism_ratios_Mel69[divergenceoverpolymorphism_ratios_Mel69$Part=="center",]

postscript("FigS1.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

ggplot(center[1:12,], aes(x=Chromosome, y=log(v)))+geom_point()+
  geom_errorbar(aes(ymax=log(Lowerbound), ymin=log(Upperbound)))+
  facet_grid(.~Mutationtype, scales="free")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),strip.text.x = element_text(size = 10))+xlab("Chromosome")+ylab("Divergence/Polymorphism")

dev.off()

######################## ############ ############ ############ ############ ############  
############ ############ FIGURE S2  ############ ############ ############ 
XoverA <- read_delim("../../Intermediate_data/XoverA_PolandDiv_together", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)


XoverA$Part<- factor(XoverA$Part, levels = c("whole","center","ends"))
XoverA$Mutationtype<- factor(XoverA$Mutationtype, levels = c("a","f","b","c","d","e", "GCneutral","GCchanging"))
label_names<- c("GCneutral"="GC-conservative", "GCchanging"="GC-changing")
XoverA_reduced<- XoverA[XoverA$Data!="Divergence:X/A",]

postscript("FigS2.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 9, height = 4, colormodel = "cmyk")

ggplot(XoverA_reduced[c(43:48,19:24),], aes(x=Part, y=ratio, colour=Data))+geom_point()+
  facet_grid(.~Mutationtype, scales="free",labeller = as_labeller(label_names))+geom_errorbar(aes(ymax=ratioUB, ymin=ratioLB))+
  geom_hline(yintercept = c(0.75,1), linetype="dashed")+ylab("X/A")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),strip.text.x = element_text(size = 10))+xlab("Chromosome Part")+
  scale_x_discrete(labels=c("WChr","CChr","PChr"))

dev.off()

###################################################################################################
########################### FIGURE S3 ###############################################################

XoverA_r4bin <- read_delim("../../Intermediate_data/XoverA_PolandDiv_together_r4bin", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)


XoverA_r4bin_WChr<- XoverA_r4bin[XoverA_r4bin$binning=="WChr",]

XoverA_r4bin_WChr$Recombination<- factor(XoverA_r4bin_WChr$Recombination, levels = c("high","modhigh","modlow","low"))
XoverA_r4bin_WChr$Mutationtype<- factor(XoverA_r4bin_WChr$Mutationtype, levels = c("a","f","b","c","d","e", "GCneutral","GCchanging"))
label_names<- c("a"="a","f"="f","b"="b","c"="c","d"="d","e"="e",
                "GCneutral"="GC-conservative", "GCchanging"="GC-changing")



postscript("FigS3.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 9, height = 4, colormodel = "cmyk")

ggplot(XoverA_r4bin_WChr[-c(25:32,57:64),], aes(x=Recombination, y=ratio, colour=Data))+geom_point()+
  facet_grid(.~Mutationtype, scales="free")+geom_errorbar(aes(ymax=ratioUB, ymin=ratioLB))+
  geom_hline(yintercept = c(0.75,1), linetype="dashed")+ylab("X/A")+theme(axis.text.x = element_text(size=10, angle=90),
                                                                          panel.background = element_rect(fill = "white", colour="black"),
                                                                          panel.grid.minor = element_line(colour = "grey90"),
                                                                          panel.grid.major = element_line(colour = "grey90"),
                                                                          strip.text.x = element_text(size = 10))+
  xlab("Recombination rate bins")

dev.off()

###################################################################################################
########################### FIGURE S4 ###############################################################

XoverA_r4bin_CChr<- XoverA_r4bin[XoverA_r4bin$binning=="CChr",]

XoverA_r4bin_CChr$Recombination<- factor(XoverA_r4bin_CChr$Recombination, levels = c("high","modhigh","modlow","low"))
XoverA_r4bin_CChr$Mutationtype<- factor(XoverA_r4bin_CChr$Mutationtype, levels = c("a","f","b","c","d","e", "GCneutral","GCchanging"))
label_names<- c("a"="a","f"="f","b"="b","c"="c","d"="d","e"="e",
                "GCneutral"="GC-conservative", "GCchanging"="GC-changing")


postscript("FigS4.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 9, height = 4, colormodel = "cmyk")

ggplot(XoverA_r4bin_CChr[-c(25:32,57:64),], aes(x=Recombination, y=ratio, colour=Data))+geom_point()+
  facet_grid(.~Mutationtype, scales="free",labeller = as_labeller(label_names))+geom_errorbar(aes(ymax=ratioUB, ymin=ratioLB))+
  geom_hline(yintercept = c(0.75,1), linetype="dashed")+ylab("X/A")+theme(axis.text.x = element_text(size=10, angle=90),
                                                                          panel.background = element_rect(fill = "white", colour="black"),
                                                                          panel.grid.minor = element_line(colour = "grey90"),
                                                                          panel.grid.major = element_line(colour = "grey90"),
                                                                          strip.text.x = element_text(size = 10))+
  xlab("Recombination rate bins")


dev.off()


################ ################ ################ ################ ################ 
################ Fig S5 (Table 4 and 5 - Juraj as figure)  ################ ################ ################ ################ 
###

B_Mel69_5SI_wry <- read_delim("../../Intermediate_data/B_estimates/with_ry/B_Mel69_5SI_wry", 
                                              "\t", escape_double = FALSE, trim_ws = TRUE)

some<- c("GCchanging", "Ts","Tv")
reduced_ry<-filter(B_Mel69_5SI_wry, Spectra %in% some)
reduced_ry$Part<-factor(reduced_ry$Part, levels=c("WChr", "CChr","PChr"))

### Without ry 

B_Mel69_5SI_likCI <- read_delim("../../Intermediate_data/B_estimates/B_Mel69_5SI_likCI", 
                                          "\t", escape_double = FALSE, trim_ws = TRUE)


some<- c("GCchanging", "Ts","Tv")
reduced<-filter(B_Mel69_5SI_likCI, Spectra %in% some)
reduced$Part<- c(rep("WChr",6),rep("CChr",6), rep("PChr",6))
reduced$Part<-factor(reduced$Part, levels=c("WChr", "CChr","PChr"))


reduced$Method<- rep("without ry", nrow(reduced))
reduced_ry$Method<- rep("with ry", nrow(reduced_ry))

colnames(reduced)
colnames(reduced_ry)

gamma_together<- rbind(reduced, reduced_ry) 

label_names<- c("Autosome"="Autosome", "X"="X",
                "GCchanging"="GC-changing","Ts"="Transitions", "Tv"="Transversions")


postscript("FigS5.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 8, height = 6, colormodel = "cmyk")

ggplot(gamma_together, aes(x=Part, y=gamma, colour=Method))+geom_point(aes(size=LRTpval<0.051))+
  facet_grid(Spectra~Chromosome, scales = "free",labeller = as_labeller(label_names))+theme(axis.text.x = element_text(size=10),
                                                                                            panel.background = element_rect(fill = "white", colour="black"),
                                                                                            panel.grid.minor = element_line(colour = "grey90"),
                                                                                            panel.grid.major = element_line(colour = "grey90"),
                                                                                            strip.text.x = element_text(size = 10),strip.text.y = element_text(size = 10))+
  scale_size_manual(values=c(2,5))+guides(size = "none")+
  geom_errorbar(aes(ymax=gammaUB, ymin=gammaLB))+xlab("Chromosome Part")+ylab("B")+
  scale_color_manual(values = c("#CCCC33", "#6600CC"))

dev.off()

###################################################################################################
########################### FIGURE S6 ###############################################################

B_Mel69_requal_WChr <- read_delim("../../Intermediate_data/B_estimates/B_Mel69_requal_WChr", 
                                                           "\t", escape_double = FALSE, trim_ws = TRUE)

GCchanging_requal<- B_Mel69_requal_WChr[B_Mel69_requal_WChr$Spectra=="GCchanging",]
GCchanging_requal$Recombination<- factor(GCchanging_requal$Recombination, levels = c("high","modhigh","modlow","low"))



B_Mel69_requal_CChr <- read_delim("../../Intermediate_data/B_estimates/B_Mel69_requal_CChr", 
                                                                  "\t", escape_double = FALSE, trim_ws = TRUE)


GCchanging_requalmiddle<- B_Mel69_requal_CChr[B_Mel69_requal_CChr$Spectra=="GCchanging",]
GCchanging_requalmiddle$Recombination<- factor(GCchanging_requalmiddle$Recombination, levels = c("high","modhigh","modlow","low"))



GCchanging_requal$data<- rep("WChr", nrow(GCchanging_requal))
GCchanging_requalmiddle$data<- rep("CChr", nrow(GCchanging_requalmiddle))

GCchanging_all<- rbind(GCchanging_requal, GCchanging_requalmiddle)
GCchanging_all$Recombination<- factor(GCchanging_all$Recombination, levels = c("high","modhigh","modlow","low"))
GCchanging_all$data<- factor(GCchanging_all$data, levels = c("WChr","CChr"))


postscript("FigS6.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 8, height = 4, colormodel = "cmyk")

ggplot(GCchanging_all, aes(x=Recombination, y=gamma, colour=data))+geom_point(size=3)+xlab("Recombination rate bins")+ylab("B")+
  facet_grid(.~Chromosome, scales = "free")+theme(axis.text.x = element_text(size=10),
                                                  panel.background = element_rect(fill = "white", colour="black"),
                                                  panel.grid.minor = element_line(colour = "grey90"),
                                                  panel.grid.major = element_line(colour = "grey90"))+
  geom_errorbar(aes(ymax=gammaUB, ymin=gammaLB))+
  scale_colour_manual(name=NULL, values = c("#3A6B35","#E3B448"))

dev.off()

###################################################################################################
########################### FIGURE S7- (S6 Juraj method) ###############################################################

B_Mel69_requal_WChr_wry <- read_delim("../../Intermediate_data/B_estimates/with_ry/B_Mel69_requal_WChr_wry", 
                                                           "\t", escape_double = FALSE, trim_ws = TRUE)

GCchanging_requal<- B_Mel69_requal_WChr_wry[B_Mel69_requal_WChr_wry$Spectra=="GCchanging",]
GCchanging_requal$Recombination<- factor(GCchanging_requal$Recombination, levels = c("high","modhigh","modlow","low"))



B_Mel69_requal_CChr_wry <- read_delim("../../Intermediate_data/B_estimates/with_ry/B_Mel69_requal_CChr_wry", 
                                                                  "\t", escape_double = FALSE, trim_ws = TRUE)


GCchanging_requalmiddle<- B_Mel69_requal_CChr_wry[B_Mel69_requal_CChr_wry$Spectra=="GCchanging",]
GCchanging_requalmiddle$Recombination<- factor(GCchanging_requalmiddle$Recombination, levels = c("high","modhigh","modlow","low"))



GCchanging_requal$data<- rep("WChr", nrow(GCchanging_requal))
GCchanging_requalmiddle$data<- rep("CChr", nrow(GCchanging_requalmiddle))

GCchanging_all<- rbind(GCchanging_requal, GCchanging_requalmiddle)
GCchanging_all$Recombination<- factor(GCchanging_all$Recombination, levels = c("high","modhigh","modlow","low"))
GCchanging_all$data<- factor(GCchanging_all$data, levels = c("WChr","CChr"))


postscript("FigS7.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 8, height = 4, colormodel = "cmyk")

ggplot(GCchanging_all, aes(x=Recombination, y=gamma, colour=data))+geom_point(size=3)+xlab("Recombination rate bins")+ylab("B")+
  facet_grid(.~Chromosome, scales = "free")+theme(axis.text.x = element_text(size=10),
                                                  panel.background = element_rect(fill = "white", colour="black"),
                                                  panel.grid.minor = element_line(colour = "grey90"),
                                                  panel.grid.major = element_line(colour = "grey90"))+
  geom_errorbar(aes(ymax=gammaUB, ymin=gammaLB))+
  scale_colour_manual(name=NULL, values = c("#3A6B35","#E3B448"))

dev.off()

######################## ############ ############ ############ ############ ############  
############ ############ FIGURE S8  ############ ############ ############ 

analyzedSIsites_mel69_autosome <- read_delim("../../Population_data/Melanogaster/analyzedSIsites_mel69_autosome", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)
analyzedSIsites_mel69_X <- read_delim("../../Population_data/Melanogaster/analyzedSIsites_mel69_X", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)


geneRGC<- analyzedSIsites_mel69_autosome %>% group_by(name) %>% summarise(meanGCcontent=mean(GCcontent), meanR=mean(R),
                                                                      meanFF_GC=mean(FFDS_GCcontent), meanRRC=mean(RRC))


geneRGC_X<- analyzedSIsites_mel69_X %>% group_by(name) %>% summarise(meanGCcontent=mean(GCcontent), meanR=mean(R),
                                                                 meanFF_GC=mean(FFDS_GCcontent),meanRRC=mean(RRC))

geneRGC$chromosome=rep("Autosome",nrow(geneRGC))
geneRGC_X$chromosome=rep("X",nrow(geneRGC_X))

geneRGC_long<- gather(geneRGC, key="region", value="meanGCcont",c(2,4))
geneRGC_X_long<- gather(geneRGC_X, key="region", value="meanGCcont",c(2,4))

geneRGC_all<-rbind(geneRGC_long, geneRGC_X_long)


postscript("FigS8.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")

ggplot(geneRGC_all, aes(x=(meanGCcont), fill=region, colour=region))+
  geom_density(alpha=0.3)+facet_grid(.~chromosome)+
  xlab("GC content")+ylab("Density")+ 
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+scale_colour_manual(labels=c("FFDS","5SI"), values=c("#999999", "#E69F00"))+
  scale_fill_manual(labels=c("FFDS","5SI"), values = c("#999999", "#E69F00"))+
  geom_vline(data=filter(geneRGC_all, chromosome=="Autosome"), aes(xintercept=0.3501911), colour="#E69F00", linetype="dashed") +
  geom_vline(data=filter(geneRGC_all, chromosome=="X"), aes(xintercept=0.3958597), colour="#E69F00",linetype="dashed") +
  geom_vline(data=filter(geneRGC_all, chromosome=="Autosome"), aes(xintercept=0.6315019), colour="#999999", linetype="dashed") +
  geom_vline(data=filter(geneRGC_all, chromosome=="X"), aes(xintercept=0.6971818), colour="#999999",linetype="dashed")

dev.off()

######################## ############ ############ ############ ############ ############  
############ ############ FIGURE S9  ############ ############ ############ 

####### Dmel GC bins ########
### Without ry 
B_Mel69_backgroundFFGCbins <- read_delim("../../Intermediate_data/B_estimates/B_Mel69_backgroundFFGCbins", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)


reduced<- B_Mel69_backgroundFFGCbins[B_Mel69_backgroundFFGCbins$Spectra=="GCchanging",]

B_Mel69_5SI_likCI <- read_delim("../../Intermediate_data/B_estimates/B_Mel69_5SI_likCI", 
                                          "\t", escape_double = FALSE, trim_ws = TRUE)

whole<- B_Mel69_5SI_likCI[B_Mel69_5SI_likCI$Part=="whole" & B_Mel69_5SI_likCI$Spectra=="GCchanging",]


GCcontent<- c("all","all")
whole$GCcontent<- GCcontent
whole[,c(1,10,2,3,4,5,6,8,9)]
reduced
reduced_2<- rbind(reduced,whole[,c(1,10,2,3,4,5,6,8,9)])


### With ry 

### Gamma
B_Mel69_backgroundFFGCbins_wry <- read_delim("../../Intermediate_data/B_estimates/with_ry/B_Mel69_backgroundFFGCbins_wry", 
                                                     delim = "\t", escape_double = FALSE, 
                                                     trim_ws = TRUE)


reduced_ry<- B_Mel69_backgroundFFGCbins_wry[B_Mel69_backgroundFFGCbins_wry$Spectra=="GCchanging",]

B_Mel69_5SI_wry <- read_delim("../../Intermediate_data/B_estimates/with_ry/B_Mel69_5SI_wry", 
                                              "\t", escape_double = FALSE, trim_ws = TRUE)

whole_ry<- B_Mel69_5SI_wry[B_Mel69_5SI_wry$Part=="WChr" & B_Mel69_5SI_wry$Spectra=="GCchanging",]


GCcontent<- c("all","all")
whole_ry$GCcontent<- GCcontent
whole_ry[,c(1,10,2,3,4,5,6,8,9)]
reduced_ry
reduced_2_ry<- rbind(reduced_ry,whole_ry[,c(1,10,2,3,4,5,6,8,9)])


reduced_2$Method<- rep("without ry", nrow(reduced_2))
reduced_2_ry$Method<- rep("with ry", nrow(reduced_2_ry))

gamma_together<- rbind(reduced_2, reduced_2_ry) 

gamma_together$Method<- factor(gamma_together$Method, levels = c("without ry","with ry"))


postscript("FigS9A_dodge.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")

ggplot(gamma_together, aes(x=GCcontent, y=gamma,group=interaction(Chromosome,Method) ,colour=Chromosome, shape=Method))+
  geom_point(aes(size=LRTpval<0.051),position = position_dodge(0.5))+
  xlab("GC content bins")+ylab("B")+  scale_size_manual(values=c(2,5))+guides(size = "none")+
  theme(axis.text.x = element_text(size=10),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  geom_pointrange(aes(ymax=ifelse(LRTpval<0.051, gammaUB,gamma), ymin=ifelse(LRTpval<0.051, gammaLB,gamma)),position = position_dodge(0.5))+
  scale_color_manual(values = c("#101820FF","chocolate"))

dev.off()

##### Dmel GC bins mutation bias #######

### without ry
correctedbias_GCbins <- read_delim("../../Intermediate_data/B_estimates/correctedbias_GCbins", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

correctedbias_alongchr <- read_delim("../../Intermediate_data/B_estimates/correctedbias_alongchr", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)

whole<- correctedbias_alongchr[correctedbias_alongchr$Part=="whole" & correctedbias_alongchr$Spectra=="GCchanging",]


GCcontent<- c("all","all")
whole$GCcontent<- GCcontent
colnames(whole[,c(1,18,2:6,8:17)])

colnames(correctedbias_GCbins)

mutbias<- rbind(correctedbias_GCbins,whole[,c(1,18,2:6,8:17)])

### with ry

### Bias
correctedbias_GCbins_wry <- read_delim("../../Intermediate_data/B_estimates/with_ry/correctedbias_GCbins_wry", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

correctedbias_alongchr_wry <- read_delim("../../Intermediate_data/B_estimates/with_ry/correctedbias_alongchr_wry", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)

whole_ry<- correctedbias_alongchr_wry[correctedbias_alongchr_wry$Part=="WChr" & correctedbias_alongchr_wry$Spectra=="GCchanging",]


GCcontent<- c("all","all")
whole_ry$GCcontent<- GCcontent
colnames(whole_ry[,c(1,18,2:6,8:17)])

colnames(correctedbias_GCbins_wry)

mutbias_ry<- rbind(correctedbias_GCbins_wry,whole_ry[,c(1,18,2:6,8:17)])

mutbias$Method<- rep("without ry", nrow(mutbias))
mutbias_ry$Method<- rep("with ry", nrow(mutbias_ry))

colnames(mutbias)
colnames(mutbias_ry)

bias_together<- rbind(mutbias, mutbias_ry) 

bias_together$Method<- factor(bias_together$Method, levels = c("without ry","with ry"))

postscript("FigS9B_dodge.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")

ggplot(bias_together, aes(x=GCcontent, y=beta_VB15, group=interaction(Chromosome,Method), colour=Chromosome, shape=Method))+
  geom_point(aes(size=LRTpval<0.051),position = position_dodge(0.5))+
  geom_pointrange(aes(ymax=beta_CI2, ymin=beta_CI1),position = position_dodge(0.5))+
  scale_size_manual(values=c(2,5))+guides(size = "none")+
  theme(axis.text.x = element_text(size=10),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+ylab(expression(hat(beta)))+
  xlab("GC content bins")+ylim(0.5,0.9)+
  scale_color_manual(values = c("#101820FF","chocolate"))

dev.off()

######################## ############ ############ ############ ############ ############  
############ ############ FIGURE S10  ############ ############ ############ 

divergenceoverpolymorphism_ratios_MDsim <- read_delim("../../Intermediate_data/DivoverPol/divergenceoverpolymorphism_ratios_MDsim", 
                                                    delim = "\t", escape_double = FALSE, 
                                                    trim_ws = TRUE)

divergenceoverpolymorphism_ratios_MDsim$Part<- factor(divergenceoverpolymorphism_ratios_MDsim$Part, levels = c("whole","center","ends"))
divergenceoverpolymorphism_ratios_MDsim$Mutationtype<- factor(divergenceoverpolymorphism_ratios_MDsim$Mutationtype, levels = c("a","f","b","c","d","e","GCneutral","GCchanging", "All"))

reduced=divergenceoverpolymorphism_ratios_MDsim[1:36,]

postscript("FigS10.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

ggplot(reduced[reduced$Part=="center",], aes(x=Chromosome, y=v))+geom_point()+
  geom_errorbar(aes(ymax=ratioLB, ymin=ratioUB))+
  facet_grid(.~Mutationtype, scales="free")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.x = element_text(size = 10))+xlab("Chromosome")+ylab("Divergence/Polymorphism")

dev.off()

######################## ############ ############ ############ ############ ############  
############ ############ FIGURE S11  ############ ############ ############ 

B_MDsim_all_wry <- read_delim("../../Intermediate_data/B_estimates/Simulans/with_ry/B_MDsim_all_wry", 
                                                        "\t", escape_double = FALSE, trim_ws = TRUE)

center_ry<- B_MDsim_all_wry[B_MDsim_all_wry$Part=="CChr",]

some<- c("GCchanging", "Ts","Tv")
center_ry<-filter(center_ry, Spectra %in% some)

### Without ry 

B_MDsim_all <- read_delim("../../Intermediate_data/B_estimates/Simulans/B_MDsim_all", 
                                              "\t", escape_double = FALSE, trim_ws = TRUE)


center<- B_MDsim_all[B_MDsim_all$Part=="center",]
some<- c("GCchanging", "Ts","Tv")
center<-filter(center, Spectra %in% some)

center$Part<- rep("CChr",6)


center$Method<- rep("without ry", nrow(center))
center_ry$Method<- rep("with ry", nrow(center_ry))

colnames(center)<-colnames(center_ry)

gamma_together<- rbind(center, center_ry) 
gamma_together$LRTpval<0.05

postscript("FigS11.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 8, height = 4, colormodel = "cmyk")

ggplot(gamma_together, aes(x=Spectra, y=gamma, colour=Method))+geom_point(aes(size=LRTpval<0.051))+
  facet_grid(.~Chromosome, scales = "free")+theme(axis.text.x = element_text(size=9),
                                                  panel.background = element_rect(fill = "white", colour="black"),
                                                  panel.grid.minor = element_line(colour = "grey90"),
                                                  panel.grid.major = element_line(colour = "grey90"),
                                                  strip.text.x = element_text(size = 10),strip.text.y = element_text(size = 10))+
  scale_size_manual(values=c(5,5))+guides(size = "none")+
  geom_errorbar(aes(ymax=gammaUB, ymin=gammaLB))+xlab("GC changing mutations")+ylab("B")+
  scale_color_manual(values = c("#CCCC33", "#6600CC"))+ scale_x_discrete(breaks=c("GCchanging", "Ts", "Tv"),
                                                                         labels=c("All", "Transitions", "Transversions"))

dev.off()

######################## ############ ############ ############ ############ ############  
############ ############ FIGURE S12  ############ ############ ############ 

####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 
##############  Dsim GC bins ############### ####### ####### ####### ####### 

### Without ry 
B_MDsim_backgroundFFGCbins <- read_delim("../../Intermediate_data/B_estimates/Simulans/B_MDsim_backgroundFFGCbins", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)

reduced<- B_MDsim_backgroundFFGCbins[B_MDsim_backgroundFFGCbins$Spectra=="GCchanging",]

B_MDsim_all <- read_delim("../../Intermediate_data/B_estimates/Simulans/B_MDsim_all", 
                                              "\t", escape_double = FALSE, trim_ws = TRUE)


whole<- B_MDsim_all[B_MDsim_all$Part=="whole" & B_MDsim_all$Spectra=="GCchanging",]


GCcontent<- c("all","all")
whole$GCcontent<- GCcontent

colnames(reduced)<- colnames(whole[,c(1,10,2,3,4,5,6,8,9)])
reduced_2<- rbind(reduced,whole[,c(1,10,2,3,4,5,6,8,9)])


### With ry 

B_MDsim_backgroundFFGCbins_wry <- read_delim("../../Intermediate_data/B_estimates/Simulans/with_ry/B_MDsim_backgroundFFGCbins_wry", 
                                                   delim = "\t", escape_double = FALSE, 
                                                   trim_ws = TRUE)


reduced_ry<- B_MDsim_backgroundFFGCbins_wry[B_MDsim_backgroundFFGCbins_wry$Spectra=="GCchanging",]


B_MDsim_all_wry <- read_delim("../../Intermediate_data/B_estimates/Simulans/with_ry/B_MDsim_all_wry", 
                                                        "\t", escape_double = FALSE, trim_ws = TRUE)


whole_ry<- B_MDsim_all_wry[B_MDsim_all_wry$Part=="WChr" & B_MDsim_all_wry$Spectra=="GCchanging",]


GCcontent<- c("all","all")
whole_ry$GCcontent<- GCcontent

#colnames(reduced_ry)<- colnames(whole_ry[,c(1,10,2,3,4,5,6,8,9)])
reduced_2_ry<- rbind(reduced_ry,whole_ry[,c(1,10,2,3,4,5,6,8,9)])

reduced_2$Method<- rep("without ry", nrow(reduced_2))
reduced_2_ry$Method<- rep("with ry", nrow(reduced_2_ry))

colnames(reduced_2)<- colnames(reduced_2_ry)

gamma_together<- rbind(reduced_2, reduced_2_ry) 

gamma_together$Method<- factor(gamma_together$Method, levels = c("without ry","with ry"))

postscript("FigS12A_dodge.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")

ggplot(gamma_together, aes(x=GCcontent, y=gamma,group=interaction(Chromosome,Method), colour=Chromosome, shape=Method))+
  geom_point(aes(size=LRTpval<0.01),position = position_dodge(0.5))+
  xlab("GC content bins")+ylab("B")+  scale_size_manual(values=c(2,5))+guides(size = "none")+
  theme(axis.text.x = element_text(size=10),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  geom_pointrange(aes(ymax=ifelse(LRTpval<0.01, gammaUB,gamma), ymin=ifelse(LRTpval<0.01, gammaLB,gamma)),position = position_dodge(0.5))+
  scale_color_manual(values = c("#101820FF","chocolate"))

dev.off()


##### Dsim GC bins mutation bias #######

### without ry

correctedbias_sim_GCbins <- read_delim("../../Intermediate_data/B_estimates/Simulans/correctedbias_sim_GCbins", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)


correctedbias_sim_alongchr <- read_delim("../../Intermediate_data/B_estimates/Simulans/correctedbias_sim_alongchr", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)


whole<- correctedbias_sim_alongchr[correctedbias_sim_alongchr$Part=="whole" & correctedbias_sim_alongchr$Spectra=="GCchanging",]


GCcontent<- c("all","all")
whole$GCcontent<- GCcontent
colnames(whole[,c(1,18,2:6,8:17)])

mutbias<- rbind(correctedbias_sim_GCbins,whole[,c(1,18,2:6,8:17)])

### with ry

correctedbias_sim_GCbins_wry <- read_delim("../../Intermediate_data/B_estimates/Simulans/with_ry/correctedbias_sim_GCbins_wry", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)


correctedbias_sim_alongchr_wry <- read_delim("../../Intermediate_data/B_estimates/Simulans/with_ry/correctedbias_sim_alongchr_wry", 
                                                   delim = "\t", escape_double = FALSE, 
                                                   trim_ws = TRUE)


whole_ry<- correctedbias_sim_alongchr_wry[correctedbias_sim_alongchr_wry$Part=="WChr" & correctedbias_sim_alongchr_wry$Spectra=="GCchanging",]


GCcontent<- c("all","all")
whole_ry$GCcontent<- GCcontent
colnames(whole_ry[,c(1,18,2:6,8:17)])
colnames(correctedbias_sim_GCbins_wry)

mutbias_ry<- rbind(correctedbias_sim_GCbins_wry,whole_ry[,c(1,18,2:6,8:17)])

mutbias$Method<- rep("without ry", nrow(mutbias))
mutbias_ry$Method<- rep("with ry", nrow(mutbias_ry))

colnames(mutbias)
colnames(mutbias_ry)

bias_together<- rbind(mutbias, mutbias_ry) 

bias_together$Method<- factor(bias_together$Method, levels = c("without ry","with ry"))


postscript("FigS12B_dodge.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")

ggplot(bias_together, aes(x=GCcontent, y=beta_VB15,group=interaction(Chromosome,Method), colour=Chromosome, shape=Method))+
  geom_point(aes(size=LRTpval<0.01),position = position_dodge(0.5))+
  geom_pointrange(aes(ymax=beta_CI2, ymin=beta_CI1),position = position_dodge(0.5))+
  scale_size_manual(values=c(2,5))+guides(size = "none")+
  theme(axis.text.x = element_text(size=10),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+ylab(expression(hat(beta)))+
  xlab("GC content bins")+ylim(0.5,0.9)+
  scale_color_manual(values = c("#101820FF","chocolate"))

dev.off()
