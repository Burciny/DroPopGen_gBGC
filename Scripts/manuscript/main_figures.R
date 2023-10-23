setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(readr)
################ ################ ################ ################ ################ 
################ FIGURE 2 ################ ################ ################ ################ 

#### Autosome

PolandDiv_together_autosome <- read_delim("../../Intermediate_data/PolandDiv_together_autosome", 
                                                  "\t", escape_double = FALSE, trim_ws = TRUE)


PolandDiv_together_autosome$Part<- factor(PolandDiv_together_autosome$Part, levels = c("whole", "center","ends"))
PolandDiv_together_autosome$Data<- factor(PolandDiv_together_autosome$Data, levels = c("Polymorphism", "Divergence"))
PolandDiv_together_autosome$Mutationtype<- factor(PolandDiv_together_autosome$Mutationtype, levels = c("a", "f", "b","c","d","e"))

postscript("Fig2A.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

ggplot(PolandDiv_together_autosome[1:36,], aes(x=Part, y=Value))+geom_point()+
  geom_errorbar(aes(ymax=Upperbound, ymin=Lowerbound))+
  facet_grid(Data~Mutationtype, scales="free")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.y = element_text(size = 10),strip.text.x = element_text(size = 10))+xlab("Chromosome Part")+ylab("Estimate")+
  scale_x_discrete(labels=c("WChr","CChr","PChr"))

dev.off()

######## X 
PolandDiv_together_X <- read_delim("../../Intermediate_data/PolandDiv_together_X", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)


PolandDiv_together_X$Part<- factor(PolandDiv_together_X$Part, levels = c("whole", "center","ends"))
PolandDiv_together_X$Data<- factor(PolandDiv_together_X$Data, levels = c("Polymorphism", "Divergence"))
PolandDiv_together_X$Mutationtype<- factor(PolandDiv_together_X$Mutationtype, levels = c("a", "f", "b","c","d","e"))

postscript("Fig2B.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

ggplot(PolandDiv_together_X[1:36,], aes(x=Part, y=Value))+geom_point()+
  geom_errorbar(aes(ymax=Upperbound, ymin=Lowerbound))+
  facet_grid(Data~Mutationtype, scales="free")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.y = element_text(size = 10),strip.text.x = element_text(size = 10))+xlab("Chromosome Part")+ylab("Estimate")+
  scale_x_discrete(labels=c("WChr","CChr","PChr"))
dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 3 ################ ################ ################ ################ 
XoverA <- read_delim("../../Intermediate_data/XoverA_PolandDiv_together", delim = "\t", 
                     escape_double = FALSE, trim_ws = TRUE)

XoverA$Part<- factor(XoverA$Part, levels = c("whole", "center","ends"))
XoverA$Mutationtype<- factor(XoverA$Mutationtype, levels = c("a", "f", "b","c","d","e", "GCchanging","GCneutral"))

XoverA_reduced<- XoverA[XoverA$Data!="Divergence:X/A",]

postscript("Fig3.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 9, height = 4, colormodel = "cmyk")

ggplot(XoverA_reduced[-c(43:48,19:24),], aes(x=Part, y=ratio, colour=Data))+geom_point()+
  facet_grid(.~Mutationtype, scales="free")+geom_errorbar(aes(ymax=ratioUB, ymin=ratioLB))+
  geom_hline(yintercept = c(0.75,1), linetype="dashed")+ylab("X/A")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.x = element_text(size = 10))+xlab("Chromosome Part")+
  scale_x_discrete(labels=c("WChr","CChr","PChr"))

dev.off()

postscript("Fig3_no1line.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 9, height = 4, colormodel = "cmyk")


ggplot(XoverA_reduced[-c(43:48,19:24),], aes(x=Part, y=ratio, colour=Data))+geom_point()+
  facet_grid(.~Mutationtype, scales="free")+geom_errorbar(aes(ymax=ratioUB, ymin=ratioLB))+
  geom_hline(yintercept = c(0.75), linetype="dashed")+ylab("X/A")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.x = element_text(size = 10))+xlab("Chromosome Part")+
  scale_x_discrete(labels=c("WChr","CChr","PChr"))

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 4 ################ ################ ################ ################ 
PolandDiv_together_autosome_r4bin <- read_delim("../../Intermediate_data/PolandDiv_together_autosome_r4bin", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

PolandDiv_together_autosome_r4bin$Data<- factor(PolandDiv_together_autosome_r4bin$Data, levels = c("Polymorphism","Divergence"))
PolandDiv_together_autosome_r4bin$binning<- factor(PolandDiv_together_autosome_r4bin$binning, levels = c("WChr","CChr"))
PolandDiv_together_autosome_r4bin$Recombination<- factor(PolandDiv_together_autosome_r4bin$Recombination, levels = c("high","modhigh","modlow","low"))
PolandDiv_together_autosome_r4bin$Mutationtype<- factor(PolandDiv_together_autosome_r4bin$Mutationtype, levels = c("a", "f", "b","c","d","e","GCchanging","GCneutral"))

postscript("Fig4A.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

ggplot(PolandDiv_together_autosome_r4bin[c(1:48,65:112),], aes(x=Recombination, y=Value, colour=binning))+geom_point()+
  geom_errorbar(aes(ymax=Upperbound, ymin=Lowerbound))+
  facet_grid(Data~Mutationtype, scales="free")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.y = element_text(size = 10),strip.text.x = element_text(size = 10))+xlab("Recombination rate bins")+ylab("Estimate")+
  scale_colour_manual(name=NULL, values = c("#3A6B35","#E3B448"))

dev.off()

PolandDiv_together_X_r4bin <- read_delim("../../Intermediate_data/PolandDiv_together_X_r4bin", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)

PolandDiv_together_X_r4bin$Data<- factor(PolandDiv_together_X_r4bin$Data, levels = c("Polymorphism","Divergence"))
PolandDiv_together_X_r4bin$binning<- factor(PolandDiv_together_X_r4bin$binning, levels = c("WChr","CChr"))
PolandDiv_together_X_r4bin$Recombination<- factor(PolandDiv_together_X_r4bin$Recombination, levels = c("high","modhigh","modlow","low"))
PolandDiv_together_X_r4bin$Mutationtype<- factor(PolandDiv_together_X_r4bin$Mutationtype, levels = c("a", "f", "b","c","d","e", "GCneutral","GCchanging"))

postscript("Fig4B.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

ggplot(PolandDiv_together_X_r4bin[c(1:48,65:112),], aes(x=Recombination, y=Value,colour=binning))+geom_point()+
  geom_errorbar(aes(ymax=Upperbound, ymin=Lowerbound))+
  facet_grid(Data~Mutationtype, scales="free")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.y = element_text(size = 10),strip.text.x = element_text(size = 10))+xlab("Recombination rate bins")+ylab("Estimate")+
  scale_colour_manual(name=NULL, values = c("#3A6B35","#E3B448"))

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 5 ################ ################ ################ ################ 
Divoverpol_AandX_r4bin_together <- read_delim("../../Intermediate_data/DivoverPol/Divoverpol_AandX_r4bin_together", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

Divoverpol_AandX_r4bin_together$Chromosome<- factor(Divoverpol_AandX_r4bin_together$Chromosome, levels = c("Autosome","X"))
Divoverpol_AandX_r4bin_together$binning<- factor(Divoverpol_AandX_r4bin_together$binning, levels = c("WChr","CChr"))
Divoverpol_AandX_r4bin_together$Recombination<- factor(Divoverpol_AandX_r4bin_together$Recombination, levels = c("high","modhigh","modlow","low"))
Divoverpol_AandX_r4bin_together$Mutationtype<- factor(Divoverpol_AandX_r4bin_together$Mutationtype, levels = c("a", "f", "b","c","d","e"))


postscript("Fig5.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

ggplot(Divoverpol_AandX_r4bin_together, aes(x=Recombination, y=v, colour=binning))+geom_point()+
  geom_errorbar(aes(ymax=Lowerbound, ymin=Upperbound))+
  facet_grid(Chromosome~Mutationtype, scales="free")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.y = element_text(size = 10),strip.text.x = element_text(size = 10))+xlab("Recombination rate bins")+ylab("Divergence/Polymorphism")+
  scale_colour_manual(name=NULL, values = c("#3A6B35","#E3B448"))

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 6 ################ ################ ################ ################ 
XoverA_together_r4bin <- read_delim("../../Intermediate_data/XoverA_PolandDiv_together_r4bin", delim = "\t", 
                           escape_double = FALSE, trim_ws = TRUE)


some<- c("GCchanging","GCneutral")
gBGCmuts<- filter(XoverA_together_r4bin, Mutationtype %in% some)

gBGCmuts$Mutationtype<- factor(gBGCmuts$Mutationtype, levels = c("GCneutral","GCchanging"))
gBGCmuts$Recombination<- factor(gBGCmuts$Recombination, levels = c("high","modhigh","modlow","low"))
gBGCmuts$binning<- factor(gBGCmuts$binning, levels = c("WChr","CChr"))

label_names<- c("GCneutral"="GC-conservative", "GCchanging"="GC-changing",
                "WChr"="WChr","CChr"="CChr")


postscript("Fig6.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")

ggplot(gBGCmuts, aes(x=Recombination, y=ratio, colour=Data))+geom_point()+
  facet_grid(binning~Mutationtype, scales="free",labeller = as_labeller(label_names))+geom_errorbar(aes(ymax=ratioUB, ymin=ratioLB))+
  geom_hline(yintercept = c(0.75,1), linetype="dashed")+ylab("X/A")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.y = element_text(size = 10),strip.text.x = element_text(size = 10))+xlab("Recombination rate bins")

dev.off()


postscript("Fig6_no1line.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")

ggplot(gBGCmuts, aes(x=Recombination, y=ratio, colour=Data))+geom_point()+
  facet_grid(binning~Mutationtype, scales="free",labeller = as_labeller(label_names))+geom_errorbar(aes(ymax=ratioUB, ymin=ratioLB))+
  geom_hline(yintercept = c(0.75), linetype="dashed")+ylab("X/A")+
  theme(axis.text.x = element_text(size=10, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.y = element_text(size = 10),strip.text.x = element_text(size = 10))+xlab("Recombination rate bins")

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 7 ################ ################ ################ ################ 

### Gamma
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


postscript("Fig7A_dodge.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")

ggplot(reduced_2, aes(x=GCcontent, y=gamma, group=Chromosome, colour=Chromosome))+
  geom_point(aes(size=LRTpval<0.051),position = position_dodge(0.5))+
  xlab("GC content bins")+ylab("B")+  scale_size_manual(values=c(2,5))+guides(size = "none")+
  theme(axis.text.x = element_text(size=10),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  geom_pointrange(aes(ymax=ifelse(LRTpval<0.051, gammaUB,gamma), ymin=ifelse(LRTpval<0.051, gammaLB,gamma)),position = position_dodge(0.5))+
  scale_color_manual(values = c("#101820FF","chocolate"))


dev.off()

### Bias
correctedbias_GCbins <- read_delim("../../Intermediate_data/B_estimates/correctedbias_GCbins", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

correctedbias_alongchr <- read_delim("../../Intermediate_data/B_estimates/correctedbias_alongchr", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)

whole<- correctedbias_alongchr[correctedbias_alongchr$Part=="whole" & correctedbias_alongchr$Spectra=="GCchanging",]


GCcontent<- c("all","all")
whole$GCcontent<- GCcontent

Mutbiases<- rbind(correctedbias_GCbins,whole[,c(1,18,2:6,8:17)])


postscript("Fig7B_dodge.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")

ggplot(Mutbiases, aes(x=GCcontent, y=beta_VB15,group=Chromosome ,colour=Chromosome))+
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

################ ################ ################ ################ ################ 
################ FIGURE 8 ################ ################ ################ ################ 

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


postscript("Fig8A_dodge.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")

ggplot(reduced_2, aes(x=GCcontent, y=gamma,group=Chromosome, colour=Chromosome))+
  geom_point(aes(size=LRTpval<0.01),position = position_dodge(0.5))+
  xlab("GC content bins")+ylab("B")+  scale_size_manual(values=c(2,5))+guides(size = "none")+
  theme(axis.text.x = element_text(size=10),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  geom_pointrange(aes(ymax=ifelse(LRTpval<0.01, UB,gamma), ymin=ifelse(LRTpval<0.01, LB,gamma)),position = position_dodge(0.5))+
  scale_color_manual(values = c("#101820FF","chocolate"))+ylim(-0.07,1)

dev.off()

### Bias 

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

Mutbiases_sim<- rbind(correctedbias_sim_GCbins,whole[,c(1,18,2:6,8:17)])

postscript("Fig8B_dodge.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")

ggplot(Mutbiases_sim, aes(x=GCcontent, y=beta_VB15,group=Chromosome, colour=Chromosome))+
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

