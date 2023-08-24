#############################################################################
#### This file undertakes the meta-regression analysis based on the database 
#### from the Orellano et al (2020) paper
#############################################################################

library(readxl)
library(nloptr)
library(meta)
library(metafor)


##################################################
#### Read in the data and subset to fixed criteria
##################################################
#### Read in the data
db0 <- read_excel("Orellano data.xlsx", sheet="Database_for_R")


#### Subset to studies based on all ages and both sexes
db1 <- db0[which(db0$Age=='All' & db0$Sex=='Both'), ]


#### Subset to studies not using 1 hour exposure windows
db2 <- db1[which(db1$Time != "1 hr"), ]


##### Subset to studies having just all-cause, cardio, resp outcomes
db3 <- db2[which(db2$Outcome!="Cerebrovascular"), ]


#### Subset to studies in America, Asia, Europe only
table(db3$Continent)
dbfinal <- db3[which(db3$Continent %in% c("America", "Asia", "Europe")), ]


#### Create a table of pollutant-outcome pairs by continent
table(dbfinal$Outcome, dbfinal$Pollutant, dbfinal$Continent)




######################################################################
#### Run the meta-analytic methods for a single pollutant-outcome pair
######################################################################
#### Subset the data to the required outcome and pollutant
dbnew <- dbfinal[ which(dbfinal$Pollutant=='O3' & dbfinal$Outcome=='All-cause'), ]
dbnew$LnEE.spm.unadjusted <- as.numeric(dbnew$LnEE.spm.unadjusted)
dbnew$SELnEE.spm.unadjusted <- as.numeric(dbnew$SELnEE.spm.unadjusted)
dim(dbnew)


#### Random effects meta-analysis of all studies
model <- metagen(TE=dbnew$LnEE.spm.unadjusted, seTE=dbnew$SELnEE.spm.unadjusted, studlab=dbnew$Article, sm="RR", comb.random=TRUE, method.tau="DL", hakn=TRUE, backtransf=TRUE, prediction=TRUE, level.predict=0.80)
summary(model)

#### Forest plot
png(file="Forest_plot_o3_allcause.png",width=4800,height=5000,res=300)
forest(model, comb.fixed=FALSE, hetstat=TRUE, prediction=TRUE, plotwidth="8cm", backtransf=TRUE, digits=4, leftcols=c("studlab"), leftlabs="Study Name", xlab="Relative risk (RR)", smlab="Relative risk (RR)",
         fontsize = 17)       
dev.off()


#### Meta-regression - estimates and confidence intervals
model2 <- metareg(x=model, formula=~factor(dbnew$Continent), intercept=FALSE, hakn=TRUE)
summary(model2)
coef <- model2$beta
ci <- cbind(model2$ci.lb, model2$ci.ub)
round(exp(cbind(coef, ci)),4)


#### Meta-regression - Test for heterogeneity
model3 <- metareg(x=model, formula=~factor(dbnew$Continent), intercept=TRUE, hakn=TRUE)
summary(model3)


#### % of variation between studies due to continent
round(100 * (model$tau2 - model2$tau2) / model$tau2,1)


#### Sub-group analysis
model4 <- update.meta(model, data=dbnew, subgroup = Continent, tau.common = FALSE)
summary(model4)