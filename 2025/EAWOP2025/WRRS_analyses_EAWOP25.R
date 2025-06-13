# WRRS DATA ANALYSIS FOR EAWOP 2025
# Towards the Italian version of the Work-related Rumination Scale: Development and Psychometrics
# Authors: Luca Menghini, Merylin Monaro, Luciano Gamberini, University of Padova

# emptying worskpace and loading used packages and functions
rm(list=ls())
library(ggplot2); library(gridExtra); library(reshape2); library(plyr); 
library(lavaan); library(MuMIn); library(htmlTable)
source("functions/items.hist.R")
source("functions/black.cor.R")

##################################################################################
# Sample S1 #####################################################################
s1 <- read.csv("data/DATASETS/t1_clean.csv") # reading data
s1$wrrs_pd_2 <- 6 - s1$wrrs_pd_2 # reversing item pd_2
s2 <- s1[s1$version=="v5",] # version v5 in s2
s1 <- s1[s1$version=="v4",] # version v4 in s1

# ...............................................................................
## sample stats ####
nrow(s1) # No. participants

# descriptives
100*table(s1$gender)/nrow(s1)
round(mean(s1$age),1); round(sd(s1$age),1)
100*table(s1$job.type.isco1)/nrow(s1)
round(mean(s1$job.tenure),1); round(sd(s1$job.tenure),1)
round(mean(s1$job.hours),1); round(sd(s1$job.hours),1)
100*table(s1$job.position)/nrow(s1)
s1$remote <- FALSE
s1$remote[s1$job.hours_remote>0] <- TRUE
100*table(s1$remote)/nrow(s1)

# online vs tablet questionnaire administration
100*table(s1$type)/nrow(s1)

# ...............................................................................
## item distribution and correlations ####

# plotting response count
p <- list()
wrrs <- colnames(s1)[grep("wrrs",colnames(s1))]
wrrs_i <- c(1,5,7,9,15,2,4,8,11,13,3,6,10,12,14)
items.hist(s1,wrrs,wrrs_i)

# polycoric correlation matrix
cors <- s1
colnames(cors)[grep("wrrs",colnames(cors))] <- wrrs_i
cors <- reshape2::melt(polychoric(cors[,as.character(wrrs_i)])[["rho"]])
cors$Var1 <- factor(cors$Var1,levels=wrrs_i)
cors$Var2 <- factor(cors$Var2,levels=wrrs_i)
black.corr(cors)

# ...............................................................................
## CFA #####################################################################

# 3-factor model
m3_15 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4 + wrrs_ar_5
      PP =~ wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_4 + wrrs_pp_5
      PD =~ wrrs_pd_1 + wrrs_pd_2 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
'

# 2-factor model with AR and PP as single factor
m2a_15 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4 + wrrs_ar_5 + wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_4 + wrrs_pp_5
      PD =~ wrrs_pd_1 + wrrs_pd_2 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
'

# 2-factor model with PP and PD as single factor
m2b_15 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4 + wrrs_ar_5
      PP =~ wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_4 + wrrs_pp_5 + wrrs_pd_1 + wrrs_pd_2 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
'

# 2-factor model with AR and PD as single factor
m2c_15 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4 + wrrs_ar_5 + wrrs_pd_1 + wrrs_pd_2 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
      PP =~ wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_4 + wrrs_pp_5
'
m1_15 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4 + wrrs_ar_5 + wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_4 + wrrs_pp_5 + wrrs_pd_1 + wrrs_pd_2 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
'

# Fitting MCFA models
orderedVars <- colnames(s1)[grepl("wrrs",colnames(s1))]
fit3_15 <- cfa(m3_15,data=s1,ordered=orderedVars,std.lv=TRUE)
fit2a_15 <- cfa(m2a_15,data=s1,ordered=orderedVars,std.lv=TRUE)
fit2b_15 <- cfa(m2b_15,data=s1,ordered=orderedVars,std.lv=TRUE)
fit2c_15 <- cfa(m2c_15,data=s1,ordered=orderedVars,std.lv=TRUE)
fit1_15 <- cfa(m1_15,data=s1,ordered=orderedVars,std.lv=TRUE)

# model fit
fitind <- c("npar","rmsea","cfi","tli","srmr")
for(fit in list(fit2a_15,fit2b_15,fit2c_15)){ # firstly comparing alternative 2-factor models
  print(lavInspect(fit,"fit")[fitind]) } # best one is fit2b (PP and PD as unique) 
htmlTable(cbind(data.frame(model=c("3","2b","1")), # fit indices table
                round(rbind(lavInspect(fit3_15,"fit")[fitind],
                            lavInspect(fit2b_15,"fit")[fitind],
                            lavInspect(fit1_15,"fit")[fitind]),3)))
standardizedsolution(fit3_15) # standardized loadings
modificationindices(fit3_15)[order(modificationindices(fit3_15)$mi, # modification indices
                                decreasing=TRUE),]

# ...............................................................................
## Item reduction #####################################################################

# fitting models without items pp4, ar5, and pd2

# 3-factor
m3_14 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4 + wrrs_ar_5
      PP =~ wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_5
      PD =~ wrrs_pd_1 + wrrs_pd_2 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
'
m3_13 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4
      PP =~ wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_5
      PD =~ wrrs_pd_1 + wrrs_pd_2 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
'
m3_12 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4
      PP =~ wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_5
      PD =~ wrrs_pd_1 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
'

# refitting model 3 without items pp4, ar5, and pd2
fit3_14 <- cfa(m3_14,ordered=orderedVars,data=s1,std.lv=TRUE)
fit3_13 <- cfa(m3_13,ordered=orderedVars,data=s1,std.lv=TRUE)
fit3_12 <- cfa(m3_12,ordered=orderedVars,data=s1,std.lv=TRUE)

# model fit
htmlTable(cbind(data.frame(model=c("3f without pp4","3f without pp4 & ar5",
                                   "3f without pp4, ar5, & pd2")),
                round(rbind(lavInspect(fit3_14,"fit")[fitind],
                            lavInspect(fit3_13,"fit")[fitind],
                            lavInspect(fit3_12,"fit")[fitind]),3)))

# ...............................................................................
## Robustness checks ###########################################################

# replicating with ML
for(fit in list(cfa(m3_15,data=s1,std.lv=TRUE),
                cfa(m2b_15,data=s1,std.lv=TRUE),
                cfa(m1_15,data=s1,std.lv=TRUE),
                cfa(m3_14,data=s1,std.lv=TRUE),
                cfa(m3_13,data=s1,std.lv=TRUE),
                cfa(m3_12,data=s1,std.lv=TRUE))){
  print(lavInspect(fit,"fit")[fitind]) }

# replicating with MLR
for(fit in list(cfa(m3_15,data=s1,std.lv=TRUE,estimator="MLR"),
                cfa(m2b_15,data=s1,std.lv=TRUE,estimator="MLR"),
                cfa(m1_15,data=s1,std.lv=TRUE,estimator="MLR"),
                cfa(m3_14,data=s1,std.lv=TRUE,estimator="MLR"),
                cfa(m3_13,data=s1,std.lv=TRUE,estimator="MLR"),
                cfa(m3_12,data=s1,std.lv=TRUE,estimator="MLR"))){
  print(lavInspect(fit,"fit")[c("npar","rmsea.robust","cfi.robust","tli.robust",
                                "srmr")]) }

#################################################################################
# Sample S2 #####################################################################
#################################################################################
s2[,c("wrrs_ar_5","wrrs_pp_4","wrrs_pd_2")] <- NULL
s2$project <- "other"
s2.bil <- read.csv("data/DATASETS/t1_Bilendi_clean.csv") # reading data
colnames(s2.bil)[grep("wrrs",colnames(s2.bil))] <-
  colnames(s2)[grep("wrrs",colnames(s2))]
s2.bil$project <- "Bilendi"
s2 <- rbind(s2[,colnames(s2)%in%colnames(s2.bil)],
            s2.bil[,colnames(s2.bil)%in%colnames(s2)])

# descriptives
nrow(s2)
100*table(s2$gender)/nrow(s2)
round(mean(s2$age),1); round(sd(s2$age),1)
100*table(s2$job.type.isco1)/nrow(s2)
round(mean(s2$job.tenure),1); round(sd(s2$job.tenure),1)
round(mean(s2$job.hours),1); round(sd(s2$job.hours),1)
100*table(s2$job.position)/nrow(s2)
s2$remote <- FALSE
s2$remote[s2$job.hours_remote>0] <- TRUE
100*table(s2$remote)/nrow(s2)

# ...............................................................................
## item distribution and correlations ####

# plotting response count
p <- list()
wrrs <- colnames(s2)[grep("wrrs",colnames(s2))]
wrrs_i <- c(1,5,7,9,2,4,8,13,3,10,12,14)
items.hist(s2,wrrs,wrrs_i)

# polycoric correlation matrix
cors <- s2
colnames(cors)[grep("wrrs",colnames(cors))] <- wrrs_i
cors <- reshape2::melt(polychoric(cors[,as.character(wrrs_i)])[["rho"]])
cors$Var1 <- factor(cors$Var1,levels=wrrs_i)
cors$Var2 <- factor(cors$Var2,levels=wrrs_i)
black.corr(cors)

# ...............................................................................
## CFA #####################################################################

# specifying 2-F and 3-F models without items pp4, ar5, and pd2
m2a_12 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4 + wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_5
      PD =~ wrrs_pd_1 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
'
m2b_12 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4
      PP =~ wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_5 + wrrs_pd_1 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
'
m2c_12 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4 + wrrs_pd_1 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
      PP =~ wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_5
'
m1_12 <- '
      AR =~ wrrs_ar_1 + wrrs_ar_2 + wrrs_ar_3 + wrrs_ar_4 + wrrs_pp_1 + wrrs_pp_2 + wrrs_pp_3 + wrrs_pp_5 + wrrs_pd_1 + wrrs_pd_3 + wrrs_pd_4 + wrrs_pd_5
'

# Fitting MCFA models
orderedVars <- colnames(s2)[grepl("wrrs",colnames(s2))]
fit3_12 <- cfa(m3_12,data=s2,ordered=orderedVars,std.lv=TRUE)
fit2a_12 <- cfa(m2a_12,data=s2,ordered=orderedVars,std.lv=TRUE)
fit2b_12 <- cfa(m2b_12,data=s2,ordered=orderedVars,std.lv=TRUE)
fit2c_12 <- cfa(m2c_12,data=s2,ordered=orderedVars,std.lv=TRUE)
fit1_12 <- cfa(m1_12,data=s2,ordered=orderedVars,std.lv=TRUE)

# model fit
fitind <- c("npar","rmsea","cfi","tli","srmr")
for(fit in list(fit2a_12,fit2b_12,fit2c_12)){ # firstly comparing alternative 2-factor models
  print(lavInspect(fit,"fit")[fitind]) } # best one is fit2a (AR and PP as unique) 
htmlTable(cbind(data.frame(model=c("3","2a","1")), # fit indices table
                round(rbind(lavInspect(fit3_12,"fit")[fitind],
                            lavInspect(fit2a_12,"fit")[fitind],
                            lavInspect(fit1_12,"fit")[fitind]),3)))
standardizedsolution(fit3_12) # standardized loadings
modificationindices(fit3_12)[order(modificationindices(fit3_12)$mi, # modification indices
                                decreasing=TRUE),][1:10,]


# ...............................................................................
## Robustness checks ###########################################################

# replicating with ML
for(fit in list(cfa(m3_12,data=s2,std.lv=TRUE),
                cfa(m2b_12,data=s2,std.lv=TRUE),
                cfa(m1_12,data=s2,std.lv=TRUE))){
  print(lavInspect(fit,"fit")[fitind]) }

# replicating with MLR
for(fit in list(cfa(m3_12,data=s2,std.lv=TRUE,estimator="MLR"),
                cfa(m2b_12,data=s2,std.lv=TRUE,estimator="MLR"),
                cfa(m1_12,data=s2,std.lv=TRUE,estimator="MLR"))){
  print(lavInspect(fit,"fit")[c("npar","rmsea.robust","cfi.robust","tli.robust",
                                "srmr")]) }

# replicating with WLSMV
for(fit in list(cfa(m3_12,data=s2,std.lv=TRUE,ordered=orderedVars,estimator="WLSMV"),
                cfa(m2b_12,data=s2,std.lv=TRUE,ordered=orderedVars,estimator="WLSMV"),
                cfa(m1_12,data=s2,std.lv=TRUE,ordered=orderedVars,estimator="WLSMV"))){
  print(lavInspect(fit,"fit")[c("npar","rmsea.robust","cfi.robust","tli.robust",
                                "srmr")]) }

# ...............................................................................
## Invariance ###########################################################

# Invariance by gender
s2.gender <- s2[s2$gender!="Other",]
nrow(s2.gender)
s2.gender$gender <- as.factor(s2.gender$gender)
summary(s2.gender$gender)
fit3.inv.conf <- cfa(m3_12,data=s2.gender,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="gender")
fit3.inv.metr <- cfa(m3_12,data=s2.gender,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="gender",group.equal="loadings")
fit3.inv.scal <- cfa(m3_12,data=s2.gender,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="gender",group.equal=c("loadings","intercepts"))
fits <- cbind(data.frame(model=c("conf","metr","scal")),
              round(rbind(lavInspect(fit3.inv.conf,"fit")[fitind],
                          lavInspect(fit3.inv.metr,"fit")[fitind],
                          lavInspect(fit3.inv.scal,"fit")[fitind]),3))
fits$deltaCFI[2:nrow(fits)] <- diff(fits$cfi)
htmlTable::htmlTable(fits)

# plotting
factor.scores <- lavPredict(fit3.inv.scal)
s2.gender <- rbind(cbind(s2.gender[s2.gender$gender=="F",],factor.scores[["F"]]),
                   cbind(s2.gender[s2.gender$gender=="M",],factor.scores[["M"]]))
p1 <- ggplot(s2.gender,aes(x=gender,y=AR)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white")) 
p2 <- ggplot(s2.gender,aes(x=gender,y=PP)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white")) 
p3 <- ggplot(s2.gender,aes(x=gender,y=PD)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white")) 
grid.arrange(p1,p2,p3,nrow=1)

# Invariance by age
s2[s2$age>=18 & s2$age<=24,"age.group"] <- "a18_24"
s2[s2$age>=25 & s2$age<=34,"age.group"] <- "a25_34"
s2[s2$age>=35 & s2$age<=44,"age.group"] <- "a35_44"
s2[s2$age>=45 & s2$age<=54,"age.group"] <- "a45_54"
s2[s2$age>=55 & s2$age<=64,"age.group"] <- "a55_64"
s2[s2$age>=65 & s2$age<=80,"age.group"] <- "a65_80"
s2$age.group <- as.factor(s2$age.group)
summary(s2$age.group)
s2.age <- s2[s2$age>=25 & s2$age<=64,]
nrow(s2.age)
s2.age$age.group <- as.factor(as.character(s2.age$age.group))
summary(s2.age$age.group)

fit3.inv.conf <- cfa(m3_12,data=s2.age,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="age.group")
fit3.inv.metr <- cfa(m3_12,data=s2.age,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="age.group",group.equal="loadings")
fit3.inv.scal <- cfa(m3_12,data=s2.age,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="age.group",group.equal=c("loadings","intercepts"))
fits <- cbind(data.frame(model=c("conf","metr","scal")),
              round(rbind(lavInspect(fit3.inv.conf,"fit")[fitind],
                          lavInspect(fit3.inv.metr,"fit")[fitind],
                          lavInspect(fit3.inv.scal,"fit")[fitind]),3))
fits$deltaCFI[2:nrow(fits)] <- diff(fits$cfi)
htmlTable::htmlTable(fits)

# plotting
factor.scores <- lavPredict(fit3.inv.scal)
s2.age <- rbind(cbind(s2.age[s2.age$age.group=="a25_34",],factor.scores[["a25_34"]]),
                cbind(s2.age[s2.age$age.group=="a35_44",],factor.scores[["a35_44"]]),
                cbind(s2.age[s2.age$age.group=="a45_54",],factor.scores[["a45_54"]]),
                cbind(s2.age[s2.age$age.group=="a55_64",],factor.scores[["a55_64"]]))
p1 <- ggplot(s2.age,aes(x=gsub("a","",gsub("_","-",age.group)),y=AR)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1)) 
p2 <- ggplot(s2.age,aes(x=gsub("a","",gsub("_","-",age.group)),y=PP)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1)) 
p3 <- ggplot(s2.age,aes(x=gsub("a","",gsub("_","-",age.group)),y=PD)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1))
grid.arrange(p1,p2,p3,nrow=1)

# Invariance by occupation
s2.isco <- s2[s2$job.type.isco1 %in%
                c("Professionals","Technicians and associate professionals",
                  "Managers","Clerical support workers"),]
nrow(s2.isco)
s2.isco$job.type.isco1 <- as.factor(s2.isco$job.type.isco1)
summary(s2.isco$job.type.isco1)
fit3.inv.conf <- cfa(m3_12,data=s2.isco,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="job.type.isco1")
fit3.inv.metr <- cfa(m3_12,data=s2.isco,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="job.type.isco1",group.equal="loadings")
fit3.inv.scal <- cfa(m3_12,data=s2.isco,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="job.type.isco1",group.equal=c("loadings","intercepts"))
fits <- cbind(data.frame(model=c("conf","metr","scal")),
              round(rbind(lavInspect(fit3.inv.conf,"fit")[fitind],
                          lavInspect(fit3.inv.metr,"fit")[fitind],
                          lavInspect(fit3.inv.scal,"fit")[fitind]),3))
fits$deltaCFI[2:nrow(fits)] <- diff(fits$cfi)
htmlTable::htmlTable(fits)

# plotting
factor.scores <- lavPredict(fit3.inv.scal)
s2.isco <- rbind(cbind(s2.isco[s2.isco$job.type.isco1=="Clerical support workers",],factor.scores[["Clerical support workers"]]),
                cbind(s2.isco[s2.isco$job.type.isco1=="Technicians and associate professionals",],factor.scores[["Technicians and associate professionals"]]),
                cbind(s2.isco[s2.isco$job.type.isco1=="Managers",],factor.scores[["Managers"]]),
                cbind(s2.isco[s2.isco$job.type.isco1=="Professionals",],factor.scores[["Professionals"]]))
s2.isco$job.type.isco1 <-
  as.factor(gsub("Clerical support workers","Clerical",
                 gsub("Technicians and associate professionals","Tech/Ass",
                      gsub("Professionals","Profess.",s2.isco$job.type.isco1))))
p1 <- ggplot(s2.isco,aes(x=gsub("a","",gsub("_","-",job.type.isco1)),y=AR)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1)) 
p2 <- ggplot(s2.isco,aes(x=gsub("a","",gsub("_","-",job.type.isco1)),y=PP)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1)) 
p3 <- ggplot(s2.isco,aes(x=gsub("a","",gsub("_","-",job.type.isco1)),y=PD)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1))
grid.arrange(p1,p2,p3,nrow=1)

# Invariance by area
s2$residence <- as.factor(s2$residence)
summary(s2$residence)
fit3.inv.conf <- cfa(m3_12,data=s2,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="residence")
fit3.inv.metr <- cfa(m3_12,data=s2,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="residence",group.equal="loadings")
fit3.inv.scal <- cfa(m3_12,data=s2,ordered=orderedVars,std.lv=TRUE,estimator="WLSMV",
                     group="residence",group.equal=c("loadings","intercepts"))
fits <- cbind(data.frame(model=c("conf","metr","scal")),
              round(rbind(lavInspect(fit3.inv.conf,"fit")[fitind],
                          lavInspect(fit3.inv.metr,"fit")[fitind],
                          lavInspect(fit3.inv.scal,"fit")[fitind]),3))
fits$deltaCFI[2:nrow(fits)] <- diff(fits$cfi)
htmlTable::htmlTable(fits)

# plotting
factor.scores <- lavPredict(fit3.inv.scal)
s2 <- rbind(cbind(s2[s2$residence=="Center",],factor.scores[["Center"]]),
                 cbind(s2[s2$residence=="Insular",],factor.scores[["Insular"]]),
                 cbind(s2[s2$residence=="Northeast",],factor.scores[["Northeast"]]),
                 cbind(s2[s2$residence=="Northwest",],factor.scores[["Northwest"]]),
                 cbind(s2[s2$residence=="South",],factor.scores[["South"]]))
s2$residence <-
  as.factor(gsub("Northeast","NE",
                 gsub("Northwest","NO",s2$residence)))
p1 <- ggplot(s2,aes(x=residence,y=AR)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1)) 
p2 <- ggplot(s2,aes(x=residence,y=PP)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1)) 
p3 <- ggplot(s2,aes(x=residence,y=PD)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1))
grid.arrange(p1,p2,p3,nrow=1)

# ...............................................................................
## Predictive validity ###########################################################


# composite scores
mh <- paste0("phq.gad_",1:4)
for(item in mh){ s2[,item] <- s2[,item] - 1 }
s2$depr.anx <- apply(s2[,c("phq.gad_1","phq.gad_2",
                           "phq.gad_3","phq.gad_4")],1,sum,na.rm=TRUE)
round(omega(s2[,c("phq.gad_1","phq.gad_2",
                  "phq.gad_3","phq.gad_4")],nfactor=1,poly=TRUE)$omega.tot,2)
s2$burnout <- apply(s2[,paste0("cbi_",1:5)],1,mean,na.rm=TRUE)
round(omega(s2[,paste0("cbi_",1:5)],nfactor=1,poly=TRUE)$omega.tot,2)
s2$wl <- apply(s2[,paste0("qws_",1:5)],1,mean,na.rm=TRUE)
round(omega(s2[,paste0("qws_",1:5)],nfactor=1,poly=TRUE)$omega.tot,2)
s2$ts.ti <- apply(s2[,paste0("ti_",1:3)],1,mean,na.rm=TRUE)
round(omega(s2[,paste0("ti_",1:3)],nfactor=1,poly=TRUE)$omega.tot,2)
s2$tasw <- apply(s2[,paste0("tasw_",1:4)],1,mean,na.rm=TRUE)
round(omega(s2[,paste0("tasw_",1:4)],nfactor=1,poly=TRUE)$omega.tot,2)
s2$wfc <- apply(s2[,paste0("wfc_",1:5)],1,mean,na.rm=TRUE)
round(omega(s2[,paste0("wfc_",1:5)],nfactor=1,poly=TRUE)$omega.tot,2)
s2$wrrs.ar <- apply(s2[,grepl("wrrs_ar_",colnames(s2))],1,mean,na.rm=TRUE)
s2$wrrs.pp <- apply(s2[,grepl("wrrs_pp_",colnames(s2))],1,mean,na.rm=TRUE)
s2$wrrs.pd <- apply(s2[,grepl("wrrs_pd_",colnames(s2))],1,mean,na.rm=TRUE)
s2$sd <- apply(s2[,paste0("msq_",1:4)],1,mean,na.rm=TRUE)
round(omega(s2[,paste0("msq_",1:4)],nfactor=1,poly=TRUE)$omega.tot,2)

# correlation matrices
cors <- reshape2::melt(cor(s2[,c("wrrs.ar","wrrs.pp","wrrs.pd",
                                 "wl","ts.ti","tasw","wfc",
                                 "sd","depr.anx","burnout")]))
cors <- cors[grepl("wrrs",cors$Var2) & !grepl("wrrs",cors$Var1),]
cors.stressors <- cors[cors$Var1%in%c("wl","ts.ti","tasw","wfc"),]
ggplot(data = cors.stressors,
       aes(x=Var1,y=Var2, fill=value)) + geom_tile() + labs(x="",y="") +
  # geom_text(aes(x=Var1,y=Var2,label=round(value,2)),color="black",size=3.5)+
  scale_fill_gradient2(low="darkblue",high="#f03b20",mid="white",midpoint=0,
                       limit = c(-1,1),space="Lab",name="Polychoric\nCorrelation",
                       guide="legend",breaks=round(seq(1,-1,length.out = 11),2),
                       minor_breaks=round(seq(1,-1,length.out = 11),2))+
  geom_text(aes(label=round(value,2))) +
  theme(axis.text=element_text(face="bold",size=10,color="white"),
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        legend.position = "none",
        panel.grid = element_blank())

cors.outcomes <- cors[cors$Var1%in%c("burnout","sd","depr.anx"),]
cors.outcomes$Var1 <- factor(cors.outcomes$Var1,levels=c("burnout","sd","depr.anx"))
ggplot(data = cors.outcomes,
       aes(x=Var1,y=Var2, fill=value)) + geom_tile() + labs(x="",y="") +
  # geom_text(aes(x=Var1,y=Var2,label=round(value,2)),color="black",size=3.5)+
  scale_fill_gradient2(low="darkblue",high="#f03b20",mid="white",midpoint=0,
                       limit = c(-1,1),space="Lab",name="Pearson\nCorrelation",
                       guide="legend",breaks=round(seq(1,-1,length.out = 11),2),
                       minor_breaks=round(seq(1,-1,length.out = 11),2))+
  geom_text(aes(label=round(value,2))) +
  theme(axis.text=element_text(face="bold",size=10,color="white"),
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        legend.position = "none",
        # legend.background = element_rect(fill="black"),
        # legend.text = element_text(color="white"),
        # legend.title = element_text(color="white",face="bold"),
        panel.grid = element_blank())

s2$mh.risk <- "absent"
s2[s2$depr.anx >= 3 & s2$depr.anx <= 5,"mh.risk"] <- "fair"
s2[s2$depr.anx >= 6 & s2$depr.anx <= 8,"mh.risk"] <- "moderate"
s2[s2$depr.anx >= 9,"mh.risk"] <- "severe"
s2$mh.risk <- as.factor(s2$mh.risk)

ggplot(s2,aes(x=mh.risk,y=wrrs.ar)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1))
ggplot(s2,aes(x=mh.risk,y=wrrs.pp)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1))
ggplot(s2,aes(x=mh.risk,y=wrrs.pd)) + 
  geom_point(size=1,position=position_jitter(width=0.1),color="darkgray") +
  geom_violin(aes(),alpha=0.5,color="black") +
  geom_boxplot(alpha=0.5,width=0.1) +
  theme(legend.position = "none",
        plot.background = element_rect(fill="black",color="black"),
        panel.background = element_rect(fill="black",color="black"),
        panel.grid = element_line(color="darkgray"),
        axis.text=element_text(face="bold",size=10,color="white"),
        axis.text.x=element_text(angle=45,hjust=1))
