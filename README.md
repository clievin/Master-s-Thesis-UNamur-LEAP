# Master-s-Thesis-UNamur-LEAP
Delayed effects of permethrin exposure on personality traits and (epi)genetic associated mechanisms in the mangrove rivulus, Kryptolebias marmoratus

######Working directory######

setwd("~/Master/Mémoire/Excel/Excelfinal") #definition of the working directory

######packages installation and loading######


library(tidyverse)
library(data.table)
library(lubridate)
library(lme4)
library(lmerTest)
library(sjPlot)
library(rptR)#calculer la repetabilite
library(ggpubr)
library(car)
library(bestNormalize)
library(kableExtra)
library(psych)
library(dunn.test)
library(ggplot2)
install.packages("survminer") #install from CRAN for survival curves
library(survminer)

###### PM impact on life history traits######
##### --- MORTALITY --- #####
##1 Mortality data loading ##
mort<-read.table("mortalite_cassandra.txt", header=T) #loading data into a dataframe named "mort"

head(mort)#6 first lines of mort

##2 Data preparation : transformation of the non-continuous variables ##
str(mort)#structure of the variables
mort$treatment <- as.factor(mort$treatment) #transformation of this non-continuous into a factor
mort$age_death=as.numeric(mort$age_death) #transformation of this non-continuous into a numeric variable
mort$death=as.logical(mort$death) #transformation of this non-continuous into a logical variable

##3 Exploration of the data in order to know in which direction to go for data treatment ##

summary(mort) #result summaries of the results of my dataframe mort, for each variable
describeBy(mort, group="treatment") #Report basic summary statistics by a grouping variable (treatment here), since we are mainly looking for a potential effect of the PM treatment

##4 Model creation ##
mod1 <-glm(death~treatment,data=mort,family="binomial") #Generalized Linear Model (GLM)  this model is mostly used when the response variable is a count or a proportion 

library(car)#loading car package for Anova function
Anova(mod1) # p-value=0.001886 (significative effect of the treatment on the mortality)

library(dunn.test) #test post-hoc for GLM
mort$death=as.numeric(mort$death)#x (death in this dataframe), has to be numeric for using dunn test
dunn.test(mort$death,mort$treatment) # significative difference (p-value < 0.025) of mortality between the treatments control-high, and low-high

##5 Model validation ##
#nb: GLM binomial: no model verification needed since we assume that the distribution of the variable is normal. Next lines is just for exploring.

##6 Model vizualization ##
#barplot of the percentage of survivals at the end of the experiment for each treatment#

prop_alive_controle<-100 #% of alive individuals at the end of the experiment for the treatment "control"
prop_alive_low<-100 #% of alive individuals at the end of the experiment for the treatment "low"
prop_alive_high<-60 #% of alive individuals at the end of the experiment for the treatment "high"

meanlc1<-mean((mort$death[mort$treatment=="controle"])) #get mean for the treatment "control"
errorlc1<-(sd((mort$death[mort$treatment=="controle"])))/(sqrt(length((mort$death[mort$treatment=="controle"])))) #get errors for the treatment "control"

meanlc2<-mean((mort$death[mort$treatment=="low"])) #get mean for the treatment "low"
errorlc2<-(sd((mort$death[mort$treatment=="low"])))/(sqrt(length((mort$death[mort$treatment=="low"])))) #get errors for the treatment "low"

meanlc3<-mean((mort$death[mort$treatment=="high"])) #get mean for the treatment "high"
errorlc3<-(sd((mort$death[mort$treatment=="high"])))/(sqrt(length((mort$death[mort$treatment=="high"])))) #get errors for the treatment "high"



matobs1 <-matrix(c(prop_alive_controle,prop_alive_low,prop_alive_high),nrow=1,dimnames=list(c("")))#matrix creation named matobs1

barplot(matobs1
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (?g/L)"
        ,ylab=" Survival (%)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,120)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("white","grey","grey30"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1)
axis(side=2,at=c(0,20,40,60,80,100),cex.axis=1,las=2)
abline(h= 0, col = "black")

text(0.5,110,"a") #addition of letters to indicate significative differences between treatments (same letter= no significative differences, and conversely
text(2.5,110,"a")
text(4.5,70,"b")




##### --- SURVIVAL (just for survival curves) --- #####

##1 Survival data loading ##

alive <-read.table("survival_cassandra.txt", header=T) 
head(alive)#col age_death.1 contains correction for a better visualization of the graphe. 138 in high refers to 140,136 in low refers to 140, it just helps avoiding overlapping graphs.

##2 Data preparation : transformation of the non-continuous variables into factors ##

alive$age_death=as.numeric(alive$age_death)#into a numeric variable
alive$death=as.logical(alive$death)#into a logical variable (0= alive and 1= dead). 

##3 Model vizualization : Survival curves ##

library(survival) #? voir si vraiment utile
surv <- survdiff(Surv(age_death,death)~treatment,data=alive)
surv #p= 8e-04 

library(survminer)

alivegraph <- survfit(Surv(age_death.1, death) ~ treatment, data =alive)#survfit() function creates survival curves using the Kaplan-Meier method based on a formula. 

ggsurvplot(
  alivegraph,                       #  fit object with calculated statistics
  legend="top",           # change the legend for x axis is written on the bottom of the curves
  legend.labs = c("control","low","high"), #change the labels for the x axis
  #risk.table = TRUE,         #show risk table.
  #risk.table.col = "treatment",# Change risk table color by treatment
  pval = TRUE,               #show p-value
  #conf.int = TRUE,           #show confidence intervals for point estimates of survival curves
  xlim = c(0,150),           #range of values on the x axis
  xlab = "Time in days",     # x axis title
  break.time.by = 20,        #break X axis in time intervals by 20
  ggtheme = theme_light(),   # customize plot and risk table with a theme.
  #risk.table.y.text.col = T, # colour risk table text annotations
  #surv.median.line = "hv"    #add the median survival pointer
)

##### --- GROWTH --- #####

##1 Growth data loading ##

growth<-read.table("growth.csv", header=T, sep =";", dec = "." )  #growth rate= lenght_140 - lenght_7
head(growth)
str(growth)

##2 Data preparation : transformation of the non-continuous variables into factors ##

str(mort)#structure of the variables
growth$treatment <- as.factor(growth$treatment) #transformation of this non-continuous


##3 Exploration of the data in order to know in which direction to go for data treatment ##

describeBy(growth, group="treatment") #we are mainly looking for a potential effect of the PM treatment here


##4 Model creation ##

hist(growth$growth_rate)#+- normal distribution
shapiro.test(growth$growth_rate)#normal distribution
mod2 <-aov(growth_rate ~ treatment, data=growth) #anova because there is one explicative variable (treatment) with more than 2 factors (need confirmation)

summary(mod2) #p-value= 0.101 : no significative effect of the treatment on the growth rate has been highlited


##5 Model validation ##
#for linear or mixed models, need to check mainly the residuals'distribution. We can have a normal variable distribution but not for residuals.

plot (mod2, which=1)#general analysis of residuals : ok
plot (mod2, which=2)#distribution of residuals: not sure it is ok
plot (mod2, which=4)#Cooke's distance for outlayers. Not exceeds 1 so ok
plot (mod2, which=5) #prob if leverage > 0.5. So here: not ok I suppose


# Homogeneity of variance: residuals vs predicted values
plot(resid(mod2)~fitted(mod2))# i think ok but not sure

#Independence of residuals versus individual explanatory variables   (bug)

op <- par(mfrow=c(2,2), mar=c(4,4,.5,.5))
plot(resid(mod2)~mod2$treatment), xlab="Treatment", ylab="Normalized residuals")+
  abline(h=0, lty=2)
par(op)

#Check for normality of residuals
hist(resid(mod2)) # bof
qqnorm(resid(mod2)) # bof
qqline(resid(mod2)) # bof

#maybe check if there is outlayers ?
### search for extreme values

dotchart(growth$growth_rate, group=growth$treatment)#extreme values

##6?Model vizualization : Graph barplot (no boxplot because not a normal distribution i think so)

meanlc4<-mean((growth$growth_rate[growth$treatment=="Controle"])) #get mean for the treatment "controle"
errorlc4<-(sd((growth$growth_rate[growth$treatment=="Controle"])))/(sqrt(length((growth$growth_rate[growth$treatment=="Controle"])))) #get errors for the treatment "controle"

meanlc5<-mean((growth$growth_rate[growth$treatment=="Low"])) 
errorlc5<-(sd((growth$growth_rate[growth$treatment=="Low"])))/(sqrt(length((growth$growth_rate[growth$treatment=="Low"]))))

meanlc6<-mean((growth$growth_rate[growth$treatment=="High"])) 
errorlc6<-(sd((growth$growth_rate[growth$treatment=="High"])))/(sqrt(length((growth$growth_rate[growth$treatment=="High"]))))


matgrowth <-matrix(c(meanlc4,meanlc5,meanlc6),nrow=1,dimnames=list(c("")))

barplot(matgrowth
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (?g/L)"
        ,ylab=" Growth rate mean from J7 dph to J140 dph (cm)"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,20)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("white","grey","grey30"))

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1)
axis(side=2,at=c(0,20,40,60,80,100),cex.axis=1,las=2)
abline(h= 0, col = "black")

text(0.5,110,"a")
text(2.5,110,"a")
text(4.5,70,"a")

#errors bars
arrows(0.5, meanlc4 - errorlc4, 0.5,meanlc4 + errorlc4, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc5 - errorlc5, 2.5,meanlc5 + errorlc5, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc6 - errorlc6, 4.5,meanlc6 + errorlc6, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

#faire un histogramme avec un barplot avec direct moyenne et erreur standards: je pense que c'est ok mtn du coup

##### --- FECUNDITY --- #####

##1? fecundity data loading ##

fecundity <-read.table("fecundity.csv", header=T, sep =";", dec = "." )
head(fecundity)

##2? Data preparation : transformation of the non-continuous variables into factors ##

str(fecundity)
fecundity$hatching_date=parse_date_time(fecundity$hatching_date,"dmy")
fecundity$date_first_egg=parse_date_time(fecundity$date_first_egg,"dmy")

fecundity$sexual_maturity <- as.logical(fecundity$sexual_maturity) 

##3? Model creation ## : analysis of the variable "sexual_maturity" (did the individual hatch or not?)

mod3 <-glm(sexual_maturity~treatment,data=fecundity,family="binomial") #GLM binomial

anova(mod3,test = "Chisq") #p-value= 0.08651 so no significant effect of the treatment on the probability of hatching has been highlited



##4? Model validation ##

# Check for homogeneity of variance: residuals vs predicted values
plot(resid(mod3)~fitted(mod3), xlab="Predicted values", ylab="Normalized residuals")+
  abline(h=0, lty=2) #5 points a bit lost but I think it is ok ?

#Check for independence of residuals versus individual explanatory variables (I don't know how to use these commands)

op <- par(mfrow=c(1,1), mar=c(4,4,.5,.5))
plot(resid(mod3)~test1$Replicate, xlab="Replicate", ylab="Normalized residuals")+
  abline(h=0, lty=2)
par(op)

summary(mod3) #intercept et high sont signif

#Dispersion of residuals (GLM)
mod3$deviance/mod3$df.resid #0.3637901 : mild-overdispersion (1.227392) so more variability than expected. It is closed to 1 so ok. Is it like hist and shapiro de resid ?

# These commands are important ?

plot (mod1, which=1)#analyse g?n?rale des r?sidus
plot (mod1, which=2)#distribution  des r?sidus
plot (mod1, which=4)#distance de cooke pour voir outlayers. je sais plus comment interpr?ter ?
plot (mod1, which=5)#


##5?Model vizualization ##

# Barplot : probability of hatching for an individual in each treatment 

mature <- fecundity[fecundity$sexual_maturity=="1",] #subset of the "mature" individuals
view(mature)

propmaturecontrole <- 3/17 #proportion of individuals which have hatched within the treatment controle on the total number of individuals in the treatment controle
propmatureld <- 5/14 #proportion of individuals which have hatched within the treatment low dose
propmaturehd <- 4/6 #proportion of individuals which have hatched within the treatment high dose


# no mean and sd and no error bars as the response variable is a logical. The R response is : NA_real for errorlc and meanlc
matobsmature <-matrix(c(propmaturecontrole,propmatureld,propmaturehd),nrow=1,dimnames=list(c("")))


barplot(matobsmature
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (?g/L)"
        ,ylab="Proportion of individuals that have laid eggs"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(0,1)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red")) 

axis(side=1.5,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(0,0.4,0.8,1.2, 1.6),cex.axis=1,las=2)
abline(h= 0, col = "black")


text(0.5,1,"a")  
text(2.5,1.6,"a")
text(4.5,1.7,"a")

##6? Descriptive statistics on the "mature" individuals (N too low) ##

#a) Total number of eggs hatched

summary(mature$total_eggs) #mean, median, min, max

hist(mature$total_eggs) # distribution on the left

summary(mature$total_eggs[mature$treatment== "Controle"])
summary(mature$total_eggs[mature$treatment== "Low"])
summary(mature$total_eggs[mature$treatment== "High"])

sd(mature$total_eggs[mature$treatment=="Controle"]) #6.806859
sd(mature$total_eggs[mature$treatment=="Low"])#8.774964
sd(mature$total_eggs[mature$treatment=="High"])#3.559026

#b) Age of the first hatching

summary(mature$age_first_egg)

hist(mature$age_first_egg) # distribution on the left

summary(mature$age_first_egg[mature$treatment== "Controle"])
summary(mature$age_first_egg[mature$treatment== "Low"])
summary(mature$age_first_egg[mature$treatment== "High"])

sd(mature$age_first_egg[mature$treatment=="Controle"]) #3.21455
sd(mature$age_first_egg[mature$treatment=="Low"])#5.549775
sd(mature$age_first_egg[mature$treatment=="High"])#1.414214


#c): the occurence of hatching ?
###### PM impact on personality traits (boldness - shyness) ######

##1? Shelter data loading ##

m <-read.table("memoire_cassandra_txt.txt", header=T)
head(m,n=3)

##2? Data preparation : transformation of the non-continuous variables into factors ##

str(m)

m$treatment<-as.factor(m$treatment)
m$repetition<-as.factor(m$repetition)
m$out_shelter<-as.factor(m$out_shelter)


##3? Exploration of the data in order to know in which direction to go for data treatment ##

summary(m)
describeBy(m, group="treatment") #we are mainly looking for a potential effect of the PM treatment here


####4? --- How time spent in shelter varied with repetition and treatment length of the individual? --- ####


##4?1 Model selection ##

#Distribution of the variable

hist(m$cumulative_time_shelter) # distribution does not follow a normal distribution
library(bestNormalize) # bestNormalize helps determining the best transformation
bestNormalize(m$cumulative_time_shelter) # orderNorm 

cumulative_shelter <- orderNorm(m$cumulative_time_shelter) #new vector containing the values of cumulative time in the shelter on which orderNorm transformation has been applied
x1 <- predict(cumulative_shelter) #extraction of the predicted values given by orderNorm
m$predict_cumul_shelter <- x1 #addition of a new column in the df "m" containing the new transformed variable 

hist(m$predict_cumul_shelter)# distribution does follow a normal distribution
shapiro.test(m$predict_cumul_shelter) #p-value = 0.9633 so confirmation that the distribution of the variable follows a normal distribution

# Model (selected by using ANOVA)

#Linear mixed-effects model 

M1 <- lmer(predict_cumul_shelter~treatment + repetition + length +length:repetition + treatment:repetition + treatment:length+ (1|fish), data=m, REML=F)#M1 contains all the potential explicative variables and interactions
summary(M1) 
Anova(M1) #delete first the non-significative interaction (with highest p-value): interaction treatment - repetition


M2 <- lmer(predict_cumul_shelter~treatment + repetition + length + length:repetition  + treatment:length+ (1|fish), data=m, REML=F)
summary(M2)
Anova(M2)# delete the interaction treatment - lenght

M3 <- lmer(predict_cumul_shelter~treatment + repetition + length + length:repetition + (1|fish), data=m, REML=F)
summary(M3)
Anova(M3) #delete the variable treatment

M4 <- lmer(predict_cumul_shelter~repetition + length + length:repetition +(1|fish), data=m, REML=F)
summary(M4)
Anova(M4) #keep it that way. Even though length is described by a p-value >0.05, keep it because the interaction repetition-lenght is significative.


#Comparison of the models
anova <- anova(M4, M3, M2, M1) #compare the models with an anova and add this value in a new vector called anova
anova # the best model is M4 (lower AIC value)

library(kableExtra) 
kbl <- kbl(anova)
kbl
table1 <- anova %>% kbl(caption="AOV") %>% kable_classic("striped",full_width=F) %>%
  column_spec(1, bold=T)
table1 # better presentation of models' comparison, in a table

# Run the best model with REML=TRUE (for having the right p-value)

M4bis <- lmer(predict_cumul_shelter~repetition + length + length:repetition +(1|fish), data=m, REML=T)

ranova(M4bis) # check the random effect with the likelihood ratio test ; if the random variable is significant, it means that the model is worse without the random effect > keep the random effect.
#p-value = 0.006211 : Keep the random effect

tab <- tab_model(M4bis, p.val="kr", show.df=T, show.reflvl=T, p.style="stars")
tab
save_kable(tab, file="aovCDIS.pdf") # ICC = 0.28. ICC < 0.3 so less probability to have a personality 


Anova(M4bis)#at least one of the means of the time spent in the shelter in repetition 1,2 or 3 is different from the others

library(emmeans) #post hoc comparisons
emmeans(M4bis,pairwise~repetition, adjust= "tukey")# pareil que Tukey HSD ? (car je gal?re ? l'utiliser)
#no significant p-value ... ?
#au cas o?
TukeyHSD(m, lmer(predict_cumul_shelter~repetition + length + length:repetition +(1|fish), conf.level=.95))#post hoc comparisons: pairwise comparisons between the means of each group
plot(m$predict_cumul_shelter~m$repetition)



##4?2 Model validation ##

# Check for colinearity with the Variance Inflation Factor
vif(M4bis) # GVIF <5 so no colinearity

# Check for homogeneity of variance: residuals vs predicted values

plot(resid(M4bis)~fitted(M4bis), xlab="Predicted values", ylab="Normalized residuals")+
  abline(h=0, lty=2) #integer (0) donc veut pas me l'afficher ...

#Check for independence of residuals versus individual explanatory variables (marche pas non plus)
 
par(mfrow=c(1,2), mar=c(4,4,.5,.5))

plot(resid(M4bis)~m$repetition, xlab="Repetition", ylab="Normalized residuals")+
  abline(h=0, lty=2)

plot(resid(M4bis)~m$length, xlab="Length", ylab="Normalized residuals")+
  abline(h=0, lty=2)

dev.off()

#Check for normality of residuals (rien ne s'affiche...)

hist(resid(M4bis))

qqnorm(resid(M4bis))
qqline(resid(M4bis))


##4?3 Model vizualization ##
#Effect of repetition

plot1bis<- ggplot(m) +
  aes(x = treatment, y = predict_cumul_shelter, colour = repetition) +
  geom_boxplot() +
  scale_color_hue(direction = 1) +
  xlab("Permethrin treatment")+
  ylab("Time spent in the shelter (sec)")+
  theme_minimal()+ 
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "black")
plot1bis

m$treatment <- factor(m$treatment, levels= levels (m$treatment[c("Controle", "Low", "High")]))

matrixconcentration <- c("0, 5, 200")   
plot1 <- ggplot(m) +
  aes(x = treatment, y = predict_cumul_shelter, colour = repetition) +
  geom_boxplot() +
  xlab("Permethrin concentration during exposure (?g/L)")+
  ylab("Time spent in the shelter (sec)")+
  scale_color_hue(direction = 1) +
  scale_x_discrete(breaks=c(0.5, 2 , 5), labels=c("Controle", "Low","High"))+
  theme_minimal()+
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "black") #line= mean and asterisk = sd
plot1 #pas de label sur l'axe x ??

##4?4 Repeatability ##

### Calculate conditional (adjusted) repeatability for each treatment

# The repeatability after controlling for the fixed effects.
#https://cran.r-project.org/web/packages/rptR/vignettes/rptR.html 

tab1 <- m[m$treatment=="controle",] #subset of our data set "m" with only the lines concerning the treatment controle
rpt1 <- rpt(predict_cumul_shelter~ repetition +(1|fish),grname=c("fish", "Fixed"), data= tab1, datatype="Gaussian", nboot=1000, npermut=0)#calculation of the repeatability for individuals of the controle treatment
print(rpt1) # R  = 0.02:  the "control" fish do not differ from each other for the time spent in the shelter, over the repetitions
#R for fixed effects: R  = 0.111 so it means that fixed effects are important in the model
plot(rpt1, cex.main=1) # ?a sert ? voir les valeurs de r?p?tabilit? ? ?a prouve qd m?me que qql indivs ont s?rement une personnalit? ?


tab2 <- m[m$treatment=="low",] 
rpt2 <- rpt(predict_first_latency_arena~repetition+ (1|fish),grname=c("fish", "Fixed"), data= tab2, datatype="Gaussian", nboot=1000, npermut=0)
print(rpt2) #R for fish = 0.182 : the "low dose" fish do not differ from each other for the time spent in the shelter, over the repetitions
plot(rpt2, cex.main=1)


tab3 <- m[m$treatment=="high",] 
rpt3 <- rpt(predict_first_latency_arena~repetition+ (1|fish),grname=c("fish", "Fixed"), data= tab3, datatype="Gaussian", nboot=1000, npermut=0)
print(rpt3)#R for Fish =  0.444 : the "high dose" fish do not differ from each other for the time spent in the shelter, over the repetitions

plot(rpt3, cex.main=1)#




####5? --- How first latency in the arena varied with repetition and treatment length of the individual? --- ####

##5?1 Model selection ##

#Distribution of the variable

hist(m$latency_to_first_extern) # distribution does not follow a normal distribution
library(bestNormalize) 
bestNormalize(m$latency_to_first_extern) # Standardized Yeo-Johnson Transformation

firstlatency_arena <- yeojohnson(m$latency_to_first_extern) 
x2 <- predict(firstlatency_arena) 
m$predict_first_latency_arena <- x2
hist(m$predict_first_latency_arena) # normal distribution
shapiro.test(m$predict_first_latency_arena)#p-value= 0.08501 so normal distribution

#Model 

M5 <- lmer(predict_first_latency_arena~treatment+repetition+length+length:repetition + treatment:repetition+ (1|fish), data=m, REML=F)
summary(M5)
Anova(M5) #delete the interaction treatment:repetition


M6 <- lmer(predict_first_latency_arena~treatment+repetition+length+length:repetition + (1|fish), data=m, REML=F)
summary(M6)
Anova(M6)#delete "treatment"


M7 <- lmer(predict_first_latency_arena~ repetition+length+length:repetition +  (1|fish), data=m, REML=F)
summary(M7)
Anova(M7)#delete repetition:lenght 

M8 <- lmer(predict_first_latency_arena~repetition+length+ (1|fish), data=m, REML=F)
summary(M8)
Anova(M8)#delete "length"

M9 <- lmer(predict_first_latency_arena~repetition+ (1|fish), data=m, REML=F)
summary(M9)
Anova(M9) #at least one mean of the first latency to the arena is different between the other repetitions 


#Comparison of the models
anova2 <- anova(M5,M6, M7, M8, M9)
anova2 # best model is M9 (lower AIC)


library(kableExtra) 
kbl2 <- kbl(anova2)
kbl2
table2 <- anova2 %>% kbl(caption="AOV") %>% kable_classic("striped",full_width=F) %>%
  column_spec(1, bold=T)
table2 # better presentation of models' comparison, in a table

# Run the best model with REML=TRUE (for having the right p-value)

M9bis <- lmer(predict_first_latency_arena~repetition+ (1|fish), data=m, REML=T)

ranova(M9bis) # check the random effect with the likelihood ratio test ; if the random variable is significant, it means that the model is worse without the random effect > keep the random effect.
#p-value = 0.006802 : Keep the random effect

tab2 <- tab_model(M9bis, p.val="kr", show.df=T, show.reflvl=T, p.style="stars")
tab2
save_kable(tab2, file="aovCDIS.pdf") # ICC = 0.28. ICC < 0.3 so less probability to have a personality 


Anova(M9bis)#at least one of the means of the time spent in the shelter in repetition 1,2 or 3 is different from the others

library(emmeans) #post hoc comparisons
emmeans(M9bis,pairwise~repetition, adjust= "tukey")#significative difference between repetition 1 and 3


##5?2 Model validation ##

# Check for colinearity with the Variance Inflation Factor : not applicable because there is only one explanatory variable

# Check for homogeneity of variance: residuals vs predicted values

plot(resid(M9bis)~fitted(M9bis), xlab="Predicted values", ylab="Normalized residuals")+
  abline(h=0, lty=2) #integer (0) donc veut pas me l'afficher ...

#Check for independence of residuals versus individual explanatory variables (marche pas non plus)

par(mfrow=c(1,1), mar=c(4,4,.5,.5))

plot(resid(M9bis)~m$repetition, xlab="Repetition", ylab="Normalized residuals")+
  abline(h=0, lty=2) #integer(0) ...

dev.off()

#Check for normality of residuals (rien ne s'affiche...)

hist(resid(M9bis))#distribution of residuals does follow a normal distribution

qqnorm(resid(M9bis))# ok
qqline(resid(M9bis))# ok


##5?3 Model vizualization ##

#Effect of repetition

plot3 <- ggplot(m) +
  aes(x = treatment, y = predict_first_latency_arena, colour = repetition) +
  geom_boxplot() +
  scale_color_hue(direction = 1) +
  xlab("Permethrin treatment")+
  ylab("First latency to the arena (sec)")+
  theme_minimal()+ 
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "black")
plot3

m$treatment <- factor(m$treatment, levels= levels (m$treatment[c("Controle", "Low", "High")]))

matrixconcentration <- c("0, 5, 200")   
plot1 <- ggplot(m) +
  aes(x = treatment, y = predict_cumul_shelter, colour = repetition) +
  geom_boxplot() +
  xlab("Permethrin concentration during exposure (?g/L)")+
  ylab("Time spent in the shelter (sec)")+
  scale_color_hue(direction = 1) +
  scale_x_discrete(breaks=c(0.5, 2 , 5), labels=c("Controle", "Low","High"))+
  theme_minimal()+
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "black") #line= mean and asterisk = sd
plot1 #pas de label sur l'axe x ??


##4?4 Repeatability ##

### Calculate conditional (adjusted) repeatability for each treatment

tab1 <- m[m$treatment=="controle",] #subset of our data set "m" with only the lines concerning the treatment controle
rpt4 <- rpt(predict_first_latency_arena~ repetition +(1|fish),grname=c("fish", "Fixed"), data= tab1, datatype="Gaussian", nboot=1000, npermut=0)#calculation of the repeatability for individuals of the controle treatment
print(rpt4) # R  = 0.267:  the "control" fish do not differ from each other for the time spent in the shelter, over the repetitions
#R for fixed effects: R  = 0.137 so it means that fixed effects are important in the model
plot(rpt4, cex.main=1) # ?a sert ? voir les valeurs de r?p?tabilit? ? ?a prouve qd m?me que qql indivs ont s?rement une personnalit? ?


tab2 <- m[m$treatment=="low",] 
rpt5 <- rpt(predict_first_latency_arena~repetition+ (1|fish),grname=c("fish", "Fixed"), data= tab2, datatype="Gaussian", nboot=1000, npermut=0)
print(rpt5) #R for fish = 0.182 : the "low dose" fish do not differ from each other for the time spent in the shelter, over the repetitions
#R for fixed effects = 0.04
plot(rpt5, cex.main=1)


tab3 <- m[m$treatment=="high",] 
rpt6 <- rpt(predict_first_latency_arena~repetition+ (1|fish),grname=c("fish", "Fixed"), data= tab3, datatype="Gaussian", nboot=1000, npermut=0)
print(rpt6)#R for Fish =  0.444 : the "high dose" fish do not differ from each other for the time spent in the shelter, over the repetitions
#R for fixed effects = 0.063
plot(rpt6, cex.main=1)




####6? --- How the ratio of  time in  the intern zone on the time spent in the arena varied with repetition, treatment and lenght of the individual? --- ####

##6?1 Model selection ##

#Distribution of the variable
m$cumularena <- (m$cumulative_duration_intern+ m$cumulative_duration_extern) #new column with the cumulated time in the arena (intern + extern parts)
m$ratio <- m$cumulative_duration_intern/m$cumularena #new column with the new variable needed for analysis

hist(m$ratio)# distribution does not follow a normal distribution
shapiro.test(m$ratio) # p-value = 2.407e-08 
bestNormalize(m$ratio) # Standardized Yeo-Johnson Transformation

newratio <- yeojohnson(m$ratio) 
x3 <- predict(newratio) 
m$predict_ratio <- x3
hist(m$predict_ratio)
shapiro.test(m$predict_ratio) # p-value = 0.03392. Distribution is still not following a normal distribution
#Model 
library(dplyr)
library(MASS)
library(nlme)
#glmm : Generalized Linear Mixed Model): not a GLM because here the data are dependent. Mixed because there is a random effect.
glmmPQL1 <- glmmPQL(ratio~treatment + repetition + length + repetition:treatment + repetition:length, random=~1|fish, family = quasibinomial, data=m)
glmm(PQL1)# bug: il me marque iteration 1, 2, 3 et 4.

##6?2 Model validation ##
##6?3 Model vizualization ##
#Effect of ?





####7?--- How the total distance moved in the arena varied with repetition, treatment and length of the individual ? --- ####

##7?1 Model selection ##
#Distribtution of the variable

m$newreltdm <- m$total_relative_distance/m$length #creation of the variable total distance moved divided by the lenght  
head(m)

hist(m$newreltdm) # does not follow a normal distribution
library(bestNormalize) 
bestNormalize(m$newreltdm) # Standardized Box Cox Transformation. lambda = 0.237862. mean (before standardization) = 11.2422. sd (before standardization) = 5.496728

Boxcox #bug

newrelativetdm <- boxcox(m$newreltdm)
x4 <- predict(newrelativetdm)
m$predict_rtdm <- x4

hist(m$predict_rtdm)
shapiro.test(m$predict_rtdm) 

hist(m$newreltdm) # distribution does not follow a normal distribution

# test avec Yeo Johnson

newrelativetdm <- yeojohnson(m$newreltdm) 
x4 <- predict(newrelativetdm) 
m$predict_rtdm <- x4 

hist(m$predict_rtdm) #distribution follows a normal distribution
shapiro.test(m$predict_rtdm) #p-value = 0.07074 so normal distribution

#Model

M10 <- lmer(predict_rtdm~treatment + repetition + treatment:repetition + (1|fish), data=m, REML=F)
summary(M10)
Anova(M10)#delete interaction treatment:repetition

M11 <- lmer(predict_rtdm~treatment + repetition + (1|fish), data=m, REML=F)
summary(M11)
Anova(M11)#delete treatment

M12 <- lmer(predict_rtdm~repetition + (1|fish), data=m, REML=F)
summary(M12)
Anova(M12)

#Comparison of the models
anova4 <- anova(M12, M11, M10)
anova4 #best model: M12 (lower AIC)

library(kableExtra) 
kbl4 <- kbl(anova4)
kbl4

table4 <- anova4 %>% kbl(caption="AOV") %>% kable_classic("striped",full_width=F) %>%
  column_spec(1, bold=T)
table4 # better presentation of models' comparison, in a table

# Run the best model with REML=TRUE (for having the right p-value)

M12bis <- lmer(predict_rtdm~repetition + (1|fish), data=m, REML=T)

ranova(M12bis) # check the random effect with the likelihood ratio test ; if the random variable is significant, it means that the model is worse without the random effect > keep the random effect.
#p-value = 0.0003932 : Keep the random effect

tab4 <- tab_model(M12bis, p.val="kr", show.df=T, show.reflvl=T, p.style="stars")
tab4
save_kable(tab4, file="aovCDIS.pdf") # ICC = 0.0.36. ICC > 0.3 so probability to have a personality 


Anova(M12bis)#at least one of the means of the total distance moved in the arena in repetition 1,2 or 3 is different from the others

library(emmeans) #post hoc comparisons
emmeans(M12bis,pairwise~repetition, adjust= "tukey")#significative difference between repetitions 1-3, and repetitions 1-2



##7?2 Model validation ##

# Check for colinearity with the Variance Inflation Factor : not applicable because there is only one explanatory variable

# Check for homogeneity of variance: residuals vs predicted values

plot(resid(M12bis)~fitted(M12bis), xlab="Predicted values", ylab="Normalized residuals")+
  abline(h=0, lty=2) #integer (0) donc veut pas me l'afficher ...

#Check for independence of residuals versus individual explanatory variables (marche pas non plus)

par(mfrow=c(1,1), mar=c(4,4,.5,.5))

plot(resid(M12bis)~m$repetition, xlab="Repetition", ylab="Normalized residuals")+
  abline(h=0, lty=2) #integer(0) ...

dev.off()

#Check for normality of residuals (rien ne s'affiche...)

hist(resid(M12bis))#distribution of residuals does follow a normal distribution
shapiro.test(resid(M12bis))# p-value = 0.7628 : normal distribution
qqnorm(resid(M12bis))# ok
qqline(resid(M12bis))# ok


##7?3 Model vizualization ##
#Effect of the repetition

plot4 <- ggplot(m) +
  aes(x = treatment, y = predict_rtdm, colour = repetition) +
  geom_boxplot() +
  scale_color_hue(direction = 1) +
  xlab("Permethrin treatment")+
  ylab("Total distance moved in the arena (cm)")+
  theme_minimal()+ 
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "black")
plot4



####8?--- ACP --- ####


###### PM impact on relative genetic expression ######
## Data loading and preparation ##

expression<-read.table("relative_gene_expression_cassandra_csv.csv", header=T, sep =";", dec = "." )
head(expression) 

str(expression)

expression$treatment<-as.factor(expression$treatment)
expression$gene<-as.factor(expression$gene)
expression$relative_gene_expression<-as.numeric(expression$relative_gene_expression)


####1? --- NipBL --- #####

## NipBL data preparation ##
nipbl <- expression[expression$gene=="nipbl",] #subset of the dataframe expression with only lines concerning nipbl


hist(nipbl$relative_gene_expression)#does not follow a normal distribution
shapiro.test(nipbl$relative_gene_expression)#does not follow a normal distribution

library(bestNormalize) 
bestNormalize(nipbl$relative_gene_expression)#OrderNorm transformation


REG <- orderNorm(nipbl$relative_gene_expression)
y1 <- predict(REG) 
nipbl$predict_REG <- y1

hist(nipbl$predict_REG)#does follow a normal distribution

## Model preparation ##

mod1 <- aov(predict_REG ~ treatment, data= nipbl)#jsp si je dois bien faire comme ?a
summary(mod1)# p-value=0.374 : no nipbl genetic expression's mean is significantly different from a treatment to another
shapiro.test(resid(mod3)) #normalit? des r?sidus suivie (p-value = 0.998)


## Model validation ##

# Check for homogeneity of variance: residuals vs predicted values

plot(resid(mod1)~fitted(mod1)), xlab="Predicted values", ylab="Normalized residuals")+
  abline(h=0, lty=2) # ? d?passe 1 en y et en x, ? -0.3 et ? 0.2.

#Check for independence of residuals versus individual explanatory variables (marche pas non plus)

par(mfrow=c(1,1), mar=c(4,4,.5,.5))

plot(resid(mod1)~expression$treatment, xlab="Treatment", ylab="Normalized residuals")+
  abline(h=0, lty=2) #va pas non plus

dev.off()

#Check for normality of residuals

hist(resid(mod1))#distribution of residuals does follow a normal distribution
qqnorm(resid(mod1))# ok
qqline(resid(mod1))# ok m?me si au debut un peu sous la ligne

shapiro.test(nipbl$predict_REG)#p-value=1 donc normalit? des r?sidus suivie


## Model visualization ##

meanlc7<-mean((nipbl$predict_REG[nipbl$treatment=="controle"]))
errorlc7<-(sd((nipbl$predict_REG[nipbl$treatment=="controle"])))/(sqrt(length((nipbl$predict_REG[nipbl$treatment=="controle"]))))

meanlc8<-mean((nipbl$predict_REG[nipbl$treatment=="low"]))
errorlc8<-(sd((nipbl$predict_REG[nipbl$treatment=="low"])))/(sqrt(length((nipbl$predict_REG[nipbl$treatment=="low"]))))

meanlc9<-mean((nipbl$predict_REG[nipbl$treatment=="high"]))
errorlc9<-(sd((nipbl$predict_REG[nipbl$treatment=="high"])))/(sqrt(length((nipbl$predict_REG[nipbl$treatment=="high"]))))


matnipbl <-matrix(c(meanlc7,meanlc8,meanlc9),nrow=1,dimnames=list(c("")))

barplot(matnipbl
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (?g/L)"
        ,ylab="Nipbl relative gene expression"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(-1,1)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(-0.8, -0.4, 0,0.4,0.8),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc7 - errorlc7, 0.5,meanlc7 + errorlc7, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc8 - errorlc8, 2.5,meanlc8 + errorlc8, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc9 - errorlc9, 4.5,meanlc9 + errorlc9, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,1,"a")  
text(2.5,1.6,"a")
text(4.5,1.7,"a")


####2?--- DNMT3a1 --- #####
## DNMT3a1 data preparation ##
dnmt3a <- expression[expression$gene=="dnmt3a1_",] #subset of the dataframe expression with only lines concerning dnmt3a


hist(dnmt3a$relative_gene_expression)#does not follow a normal distribution
shapiro.test(dnmt3a$relative_gene_expression)#does not follow a normal distribution

library(bestNormalize) 
bestNormalize(dnmt3a$relative_gene_expression)#Standardized Yeo-Johnson Transformation


REG2 <- yeojohnson(dnmt3a$relative_gene_expression)
y2 <- predict(REG2) 
dnmt3a$predict_REG <- y2

hist(dnmt3a$predict_REG)#not sure it does follow a normal distribution
shapiro.test(dnmt3a$predict_REG)#p-value = 0.1757 : does follow a normal distribution
## Model preparation ##


mod2 <- aov(predict_REG ~ treatment, data= dnmt3a)#jsp si je dois bien faire comme ?a
summary(mod2)# p-value=0.649 : no dnmt3a genetic expression's mean is significantly different from a treatment to another
shapiro.test(resid(mod2))#  p-value: 0.3885 so it does


## Model validation ##

# Check for homogeneity of variance: residuals vs predicted values

plot(resid(mod2)~fitted(mod2), xlab="Predicted values", ylab="Normalized residuals")+
  abline(h=0, lty=2) # rep: integer(0)

#Check for independence of residuals versus individual explanatory variables (marche pas non plus)

par(mfrow=c(1,1), mar=c(4,4,.5,.5))

plot(resid(mod2)~expression$treatment, xlab="Treatment", ylab="Normalized residuals")+
  abline(h=0, lty=2) #va pas non plus

dev.off()

#Check for normality of residuals

hist(resid(mod2))#distribution of residuals does follow a normal distribution
qqnorm(resid(mod2))# ok
qqline(resid(mod2))# ok 


## Model visualization ##

meanlc11<-mean((dnmt3a$predict_REG[dnmt3a$treatment=="controle"]))
errorlc11<-(sd((dnmt3a$predict_REG[dnmt3a$treatment=="controle"])))/(sqrt(length((dnmt3a$predict_REG[dnmt3a$treatment=="controle"]))))

meanlc12<-mean((dnmt3a$predict_REG[dnmt3a$treatment=="low"]))
errorlc12<-(sd((dnmt3a$predict_REG[dnmt3a$treatment=="low"])))/(sqrt(length((dnmt3a$predict_REG[dnmt3a$treatment=="low"]))))

meanlc13<-mean((dnmt3a$predict_REG[dnmt3a$treatment=="high"]))
errorlc13<-(sd((dnmt3a$predict_REG[dnmt3a$treatment=="high"])))/(sqrt(length((dnmt3a$predict_REG[dnmt3a$treatment=="high"]))))


matdnmt <-matrix(c(meanlc11,meanlc12,meanlc13),nrow=1,dimnames=list(c("")))

barplot(matdnmt
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (?g/L)"
        ,ylab="DNMT3a relative gene expression"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(-1,1)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(-0.8, -0.4, 0,0.4,0.8),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc11 - errorlc11, 0.5,meanlc11 + errorlc11, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc12 - errorlc12, 2.5,meanlc12 + errorlc12, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc13 - errorlc13, 4.5,meanlc13 + errorlc13, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,1,"a")  
text(2.5,1.6,"a")
text(4.5,1.7,"a")


####3? --- Mecp2 --- #####


## Mecp2 data preparation ##

mecp2 <- expression[expression$gene=="mecp2_",] #subset of the dataframe expression with only lines concerning mecp2


hist(mecp2$relative_gene_expression)#does not follow a normal distribution
shapiro.test(mecp2$relative_gene_expression)#does not follow a normal distribution

library(bestNormalize) 
bestNormalize(mecp2$relative_gene_expression)#Standardized Yeo-Johnson Transformation


REG3 <- yeojohnson(mecp2$relative_gene_expression)
y3 <- predict(REG3) 
mecp2$predict_REG <- y3

hist(mecp2$predict_REG)#does follow a normal distribution
shapiro.test(mecp2$predict_REG)# p-value = 0.1327 : does follow a normal distribution


## Model preparation ##

mod3 <- aov(predict_REG ~ treatment, data= mecp2)#jsp si je dois bien faire comme ?a
summary(mod3)# p-value=0.543 : no mecp2 genetic expression's mean is significantly different from a treatment to another


## Model validation ##

# Check for homogeneity of variance: residuals vs predicted values

plot(resid(mod3)~fitted(mod3), xlab="Predicted values", ylab="Normalized residuals")+
  abline(h=0, lty=2) # rep: integer(0) mais graphic there. ok.

#Check for independence of residuals versus individual explanatory variables (marche pas non plus)

par(mfrow=c(1,1), mar=c(4,4,.5,.5))

plot(resid(mod3)~expression$treatment, xlab="Treatment", ylab="Normalized residuals")+
  abline(h=0, lty=2) #va pas non plus

dev.off()

#Check for normality of residuals

hist(resid(mod3))#distribution of residuals does follow a normal distribution
qqnorm(resid(mod3))
qqline(resid(mod3))# ok 


## Model visualization ##

meanlc14<-mean((mecp2$predict_REG[mecp2$treatment=="controle"]))
errorlc14<-(sd((mecp2$predict_REG[mecp2$treatment=="controle"])))/(sqrt(length((mecp2$predict_REG[mecp2$treatment=="controle"]))))

meanlc15<-mean((mecp2$predict_REG[mecp2$treatment=="low"]))
errorlc15<-(sd((mecp2$predict_REG[mecp2$treatment=="low"])))/(sqrt(length((mecp2$predict_REG[mecp2$treatment=="low"]))))

meanlc16<-mean((mecp2$predict_REG[mecp2$treatment=="high"]))
errorlc16<-(sd((mecp2$predict_REG[mecp2$treatment=="high"])))/(sqrt(length((mecp2$predict_REG[mecp2$treatment=="high"]))))


matmecp2 <-matrix(c(meanlc14,meanlc15,meanlc16),nrow=1,dimnames=list(c("")))

barplot(matmecp2
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (?g/L)"
        ,ylab="Mecp2 relative gene expression"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(-1,1)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(-0.8, -0.4, 0,0.4,0.8),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc14 - errorlc14, 0.5,meanlc14 + errorlc14, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc15 - errorlc15, 2.5,meanlc15 + errorlc15, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc16 - errorlc16, 4.5,meanlc16 + errorlc16, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,1,"a")  
text(2.5,1.6,"a")
text(4.5,1.7,"a")



####4? --- sep15 --- #####
## Sep15 data preparation ##


sep <- expression[expression$gene=="sep_15_",] #subset of the dataframe expression with only lines concerning sep15


hist(sep$relative_gene_expression)#does not follow a normal distribution
shapiro.test(sep$relative_gene_expression)#does not follow a normal distribution

library(bestNormalize) 
bestNormalize(sep$relative_gene_expression)#Box Cox Transformation
#faut denouveau faire boxcox mais je gal?re

#Standardized Yeo-Johnson Transformation

REG4 <- yeojohnson(sep$relative_gene_expression)
y4 <- predict(REG4) 
sep$predict_REG <- y4

hist(sep$predict_REG)#does follow a normal distribution
shapiro.test(sep$predict_REG)# p-value = p-value = 0.947 : does follow a normal distribution


## Model preparation ##

mod4 <- aov(predict_REG ~ treatment, data= sep)
summary(mod4)# p-value=0.749 : no sep15 genetic expression's mean is significantly different from a treatment to another


## Model validation ##

# Check for homogeneity of variance: residuals vs predicted values

plot(resid(mod4)~fitted(mod4), xlab="Predicted values", ylab="Normalized residuals")+
  abline(h=0, lty=2) # rep: integer(0) mais graphic there. jsp si ok.

#Check for independence of residuals versus individual explanatory variables (marche pas non plus)

par(mfrow=c(1,1), mar=c(4,4,.5,.5))

plot(resid(mod4)~expression$treatment, xlab="Treatment", ylab="Normalized residuals")+
  abline(h=0, lty=2) #va pas non plus

dev.off()

#Check for normality of residuals

hist(resid(mod4))#distribution of residuals does follow a normal distribution
qqnorm(resid(mod4))
qqline(resid(mod4))# ok mais bouge bcp au d?but qd m?me


## Model visualization ##

meanlc17<-mean((sep$predict_REG[sep$treatment=="controle"]))
errorlc17<-(sd((sep$predict_REG[sep$treatment=="controle"])))/(sqrt(length((sep$predict_REG[sep$treatment=="controle"]))))

meanlc18<-mean((sep$predict_REG[sep$treatment=="low"]))
errorlc18<-(sd((sep$predict_REG[sep$treatment=="low"])))/(sqrt(length((sep$predict_REG[sep$treatment=="low"]))))

meanlc19<-mean((sep$predict_REG[sep$treatment=="high"]))
errorlc19<-(sd((sep$predict_REG[sep$treatment=="high"])))/(sqrt(length((sep$predict_REG[sep$treatment=="high"]))))


matsep <-matrix(c(meanlc17,meanlc18,meanlc19),nrow=1,dimnames=list(c("")))

barplot(matsep
        ,beside = TRUE
        , horiz = FALSE
        , legend.text = FALSE
        ,xlab="Permethrin concentration during exposure (?g/L)"
        ,ylab="Sep15 relative gene expression"
        ,cex.lab=1.2
        ,xlim=c(0,5.5)
        ,ylim=c(-1,1)
        ,lwd = 2
        ,pch=16
        ,axes=FALSE
        ,space=c(0,1,1)
        ,col=c("blue","green","red"))

axis(side=1.5,,at=c(0.5,2.5,4.5),labels=c("0","5","200"),tick=FALSE,cex.axis=1) 
axis(side=2,at=c(-0.8, -0.4, 0,0.4,0.8),cex.axis=1,las=2)
abline(h= 0, col = "black")

#errors bars
arrows(0.5, meanlc17 - errorlc17, 0.5,meanlc17 + errorlc17, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(2.5, meanlc18 - errorlc18, 2.5,meanlc18 + errorlc18, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)
arrows(4.5, meanlc19 - errorlc19, 4.5,meanlc19 + errorlc19, col = "black", angle = 90, code = 3, length = 0.1,lwd = 2)

text(0.5,1,"a")  
text(2.5,1.6,"a")
text(4.5,1.7,"a")


####5? --- ACP --- ####
###### PM impact on methylation level######

mort<-read.table("mortalite_cassandra.txt", header=T) #loading data into a dataframe named "mort"

head(mort)#6 first lines of mort





#acp pour voir si changement d'expression dans le m?me sens
#hit map pour voir visuellement quels sont les g?nes qui augmentent et/ou descendent: faire ?a avec prism.
#ou faire avec cluster (voir app)
#anova : suivre la normalit? et homog?n?it? des variances: si ok, faire ?a. si pas ok, 

