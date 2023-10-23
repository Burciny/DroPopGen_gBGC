setwd(dirname(rstudioapi::getSourceEditorContext()$path))

####### Melanogaster  ####### ####### ####### #######
####### ####### ####### ####### ####### #######
##### R bins - WChr ##### 
analyzedSIsites_mel69_autosome <- read_delim("../../Population_data/Melanogaster/analyzedSIsites_mel69_autosome", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)

bins<-cut_number(analyzedSIsites_mel69_autosome$R,4)

levels(bins)
bins_df<- as.data.frame(bins)
ggplot(bins_df,aes(x=bins))+geom_bar()
bins

analyzedSIsites_mel69_autosome$bins<-bins

Lowr<- analyzedSIsites_mel69_autosome[which(analyzedSIsites_mel69_autosome$bins=="[0,0.976]"),]
ModLowr<- analyzedSIsites_mel69_autosome[which(analyzedSIsites_mel69_autosome$bins=="(0.976,1.96]"),]
ModHighr<- analyzedSIsites_mel69_autosome[which(analyzedSIsites_mel69_autosome$bins=="(1.96,3.25]"),]
Highr<- analyzedSIsites_mel69_autosome[which(analyzedSIsites_mel69_autosome$bins=="(3.25,10.2]"),]

# write.table(Highr, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_autosome_WChr_highr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(ModHighr, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_autosome_WChr_modhighr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(ModLowr, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_autosome_WChr_modlowr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(Lowr, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_autosome_WChr_lowr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


# X
analyzedSIsites_mel69_X <- read_delim("../../Population_data/Melanogaster/analyzedSIsites_mel69_X", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)

bins<-cut_number(analyzedSIsites_mel69_X$R,4)

levels(bins)
bins_df<- as.data.frame(bins)
ggplot(bins_df,aes(x=bins))+geom_bar()
bins

analyzedSIsites_mel69_X$bins<-bins

LowrX<- analyzedSIsites_mel69_X[which(analyzedSIsites_mel69_X$bins=="[0,0.723]"),]
ModLowrX<- analyzedSIsites_mel69_X[which(analyzedSIsites_mel69_X$bins=="(0.723,2.26]"),]
ModHighrX<- analyzedSIsites_mel69_X[which(analyzedSIsites_mel69_X$bins=="(2.26,3.25]"),]
HighrX<- analyzedSIsites_mel69_X[which(analyzedSIsites_mel69_X$bins=="(3.25,12.2]"),]

# write.table(HighrX, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_X_WChr_highr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(ModHighrX, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_X_WChr_modhighr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(ModLowrX, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_X_WChr_modlowr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(LowrX, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_X_WChr_lowr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

####### ####### ####### ####### ####### #######
##### R bins - CChr ##### 

analyzedSIsites_mel69_autosome_CChr <- read_delim("../../Population_data/Melanogaster/analyzedSIsites_mel69_autosome_CChr", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)

bins<-cut_number(analyzedSIsites_mel69_autosome_CChr$R,4)

levels(bins)
bins_df<- as.data.frame(bins)
ggplot(bins_df,aes(x=bins))+geom_bar()
bins

analyzedSIsites_mel69_autosome_CChr$bins<-bins

Lowr<- analyzedSIsites_mel69_autosome_CChr[which(analyzedSIsites_mel69_autosome_CChr$bins=="[0.117,1.64]"),]
ModLowr<- analyzedSIsites_mel69_autosome_CChr[which(analyzedSIsites_mel69_autosome_CChr$bins=="(1.64,2.49]"),]
ModHighr<- analyzedSIsites_mel69_autosome_CChr[which(analyzedSIsites_mel69_autosome_CChr$bins=="(2.49,3.69]"),]
Highr<- analyzedSIsites_mel69_autosome_CChr[which(analyzedSIsites_mel69_autosome_CChr$bins=="(3.69,10.2]"),]

# write.table(Highr, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_autosome_CChr_highr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(ModHighr, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_autosome_CChr_modhighr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(ModLowr, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_autosome_CChr_modlowr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(Lowr, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_autosome_CChr_lowr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


# X
analyzedSIsites_mel69_X_CChr <- read_delim("../../Population_data/Melanogaster/analyzedSIsites_mel69_X_CChr", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)

bins<-cut_number(analyzedSIsites_mel69_X_CChr$R,4)

levels(bins)
bins_df<- as.data.frame(bins)
ggplot(bins_df,aes(x=bins))+geom_bar()
bins

analyzedSIsites_mel69_X_CChr$bins<-bins

LowrX<- analyzedSIsites_mel69_X_CChr[which(analyzedSIsites_mel69_X_CChr$bins=="[0.362,1.9]"),]
ModLowrX<- analyzedSIsites_mel69_X_CChr[which(analyzedSIsites_mel69_X_CChr$bins=="(1.9,2.62]"),]
ModHighrX<- analyzedSIsites_mel69_X_CChr[which(analyzedSIsites_mel69_X_CChr$bins=="(2.62,3.98]"),]
HighrX<- analyzedSIsites_mel69_X_CChr[which(analyzedSIsites_mel69_X_CChr$bins=="(3.98,12.2]"),]

# write.table(HighrX, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_X_CChr_highr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(ModHighrX, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_X_CChr_modhighr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(ModLowrX, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_X_CChr_modlowr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# write.table(LowrX, "../../Population_data/Melanogaster/Rbins/analyzedSIsites_mel69_X_CChr_lowr", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

####### ####### ####### ####### ####### ####### ####### #######
##### Background GC bins ####### ####### ####### #######

# Autosome
analyzedSIsites_mel69_autosome <- read_delim("../../Population_data/Melanogaster/analyzedSIsites_mel69_autosome", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)

bins<-cut_number(analyzedSIsites_mel69_autosome$FFDS_GCcontent,5)

levels(bins)
bins_df<- as.data.frame(bins)
ggplot(bins_df,aes(x=bins))+geom_bar()

analyzedSIsites_mel69_autosome$bins<-bins

GC1<- analyzedSIsites_mel69_autosome[which(analyzedSIsites_mel69_autosome$bins=="[0.198,0.556]"),]
GC2<- analyzedSIsites_mel69_autosome[which(analyzedSIsites_mel69_autosome$bins=="(0.556,0.613]"),]
GC3<- analyzedSIsites_mel69_autosome[which(analyzedSIsites_mel69_autosome$bins=="(0.613,0.657]"),]
GC4<- analyzedSIsites_mel69_autosome[which(analyzedSIsites_mel69_autosome$bins=="(0.657,0.709]"),]
GC5<- analyzedSIsites_mel69_autosome[which(analyzedSIsites_mel69_autosome$bins=="(0.709,0.963]"),]

mean(GC1$R) # 2.057051
mean(GC2$R) # 2.273045
mean(GC3$R) # 2.276601
mean(GC4$R) # 2.356648
mean(GC5$R) # 2.537986

# write.table(GC1,"../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_autosome_GC1", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC2,"../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_autosome_GC2", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC3,"../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_autosome_GC3", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC4,"../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_autosome_GC4", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC5,"../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_autosome_GC5", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)

### X 
analyzedSIsites_mel69_X <- read_delim("../../Population_data/Melanogaster/analyzedSIsites_mel69_X", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)

bins<-cut_number(analyzedSIsites_mel69_X$FFDS_GCcontent,5)

levels(bins)
bins_df<- as.data.frame(bins)
ggplot(bins_df,aes(x=bins))+geom_bar()


analyzedSIsites_mel69_X$bins<-bins

GC1<- analyzedSIsites_mel69_X[which(analyzedSIsites_mel69_X$bins=="[0.364,0.615]"),]
GC2<- analyzedSIsites_mel69_X[which(analyzedSIsites_mel69_X$bins=="(0.615,0.68]"),]
GC3<- analyzedSIsites_mel69_X[which(analyzedSIsites_mel69_X$bins=="(0.68,0.731]"),]
GC4<- analyzedSIsites_mel69_X[which(analyzedSIsites_mel69_X$bins=="(0.731,0.776]"),]
GC5<- analyzedSIsites_mel69_X[which(analyzedSIsites_mel69_X$bins=="(0.776,0.944]"),]


mean(GC1$R) # 2.304338
mean(GC2$R) # 2.353341
mean(GC3$R) # 2.564982
mean(GC4$R) # 2.339466
mean(GC5$R) # 2.679921

# write.table(GC1,"../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_X_GC1", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC2,"../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_X_GC2", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC3,"../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_X_GC3", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC4,"../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_X_GC4", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC5,"../../Population_data/Melanogaster/GCbins/analyzedSIsites_mel69_X_GC5", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)

####### Simulans  ####### ####### ####### #######
####### ####### ####### ####### ####### #######

#### GC bins ######
## Autosome 
analyzedSIsites_MDsim_autosome <- read_delim("../../Population_data/Simulans/analyzedSIsites_MDsim_autosome", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)


bins<-cut_number(analyzedSIsites_MDsim_autosome$FFDS_GCcontent,5)

levels(bins)
bins_df<- as.data.frame(bins)
ggplot(bins_df,aes(x=bins))+geom_bar()


analyzedSIsites_MDsim_autosome$bins<-bins
split(analyzedSIsites_MDsim_autosome, bins)

GC1<- analyzedSIsites_MDsim_autosome[which(analyzedSIsites_MDsim_autosome$bins=="[0,0.564]"),]
GC2<- analyzedSIsites_MDsim_autosome[which(analyzedSIsites_MDsim_autosome$bins=="(0.564,0.628]"),]
GC3<- analyzedSIsites_MDsim_autosome[which(analyzedSIsites_MDsim_autosome$bins=="(0.628,0.674]"),]
GC4<- analyzedSIsites_MDsim_autosome[which(analyzedSIsites_MDsim_autosome$bins=="(0.674,0.726]"),]
GC5<- analyzedSIsites_MDsim_autosome[which(analyzedSIsites_MDsim_autosome$bins=="(0.726,0.951]"),]

mean(GC1$RR) # 2.975123
mean(GC2$RR) # 3.0853
mean(GC3$RR) # 3.180112
mean(GC4$RR) # 3.166913
mean(GC5$RR) # 3.073402

# write.table(GC1,"../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_autosome_GC1", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC2,"../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_autosome_GC2", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC3,"../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_autosome_GC3", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC4,"../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_autosome_GC4", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC5,"../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_autosome_GC5", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)

### X 

analyzedSIsites_MDsim_X <- read_delim("../../Population_data/Simulans/analyzedSIsites_MDsim_X", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)


bins<-cut_number(analyzedSIsites_MDsim_X$FFDS_GCcontent,5)

levels(bins)
bins_df<- as.data.frame(bins)
ggplot(bins_df,aes(x=bins))+geom_bar()



analyzedSIsites_MDsim_X$bins<-bins
split(analyzedSIsites_MDsim_X, bins)

GC1<- analyzedSIsites_MDsim_X[which(analyzedSIsites_MDsim_X$bins=="[0,0.636]"),]
GC2<- analyzedSIsites_MDsim_X[which(analyzedSIsites_MDsim_X$bins=="(0.636,0.698]"),]
GC3<- analyzedSIsites_MDsim_X[which(analyzedSIsites_MDsim_X$bins=="(0.698,0.751]"),]
GC4<- analyzedSIsites_MDsim_X[which(analyzedSIsites_MDsim_X$bins=="(0.751,0.799]"),]
GC5<- analyzedSIsites_MDsim_X[which(analyzedSIsites_MDsim_X$bins=="(0.799,1]"),]


mean(GC1$RR) # 3.362282
mean(GC2$RR) # 2.545624
mean(GC3$RR) # 2.997811
mean(GC4$RR) # 3.021675
mean(GC5$RR) # 3.076495

# write.table(GC1,"../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_X_GC1", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC2,"../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_X_GC2", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC3,"../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_X_GC3", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC4,"../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_X_GC4", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
# write.table(GC5,"../../Population_data/Simulans/GCbins/analyzedSIsites_MDsim_X_GC5", col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
