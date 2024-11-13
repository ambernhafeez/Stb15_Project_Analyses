###############################
# Step 1 - read and format data
###############################

# read data
dat <- read.csv("Septoria_AUDPC_data.csv", head=T, na.strings = c("", "-","NA"), stringsAsFactors = F)

# set factors and numerical columns
dat[,1:8]<-lapply(dat[,1:8], factor);
dat[,9:ncol(dat)] <- sapply((sapply(dat[,9:ncol(dat)], as.character)), as.numeric);

# remove missing data
dat <- dat[!(is.na(dat$dAUDPC) & is.na(dat$pAUDPC)), ];

###############################
# Step 2 - calculate % max dAUDPC 
###############################

#calculate maximum possible AUDPC value for each assay (isolate and batch)
max_score_323_b1 <- 100*(32-14)
max_score_323_b2 <- 100*(32-10)
max_score_88004 <- 100*(31-17)
max_score_90012_b1 <- 100*(31-17)
max_score_90012_b2 = 100*(28-14)
#divide dAUDPC by the maximum possible AUDPC and multiply by 100
dat[dat$Isolate == "IPO323" & dat$Batch == 1, "pmax_dAUDPC"] <- (dat[dat$Isolate == "IPO323" & dat$Batch==1, "dAUDPC"]/max_score_323_b1)*100
dat[dat$Isolate == "IPO323" & dat$Batch == 2, "pmax_dAUDPC"] <- (dat[dat$Isolate == "IPO323" & dat$Batch==2, "dAUDPC"]/max_score_323_b2)*100
dat[dat$Isolate == "IPO88004", "pmax_dAUDPC"] <- (dat[dat$Isolate == "IPO88004", "dAUDPC"]/max_score_88004)*100
dat[dat$Isolate == "IPO90012" & dat$Batch == 1, "pmax_dAUDPC"] <- (dat[dat$Isolate == "IPO90012" & dat$Batch == 1, "dAUDPC"]/max_score_90012_b1)*100
dat[dat$Isolate == "IPO90012" & dat$Batch == 2, "pmax_dAUDPC"] <- (dat[dat$Isolate == "IPO90012" & dat$Batch == 2, "dAUDPC"]/max_score_90012_b2)*100

###############################
# Step 3 - fit linear mixed models
###############################

if(!require(lmerTest)){install.packages(lmerTest)}
library(lmerTest)

# pycnidia 
modp <- lmer(lgt_pAUDPC ~ Isolate*Line + Scorer +
               (1| Isolate:Batch) + (1|Isolate:Batch:Rep) + (1|Isolate:Batch:Rep:Box),
             data=dat)
anova(modp)
summary(modp)

# fit the model to damage data
modd <- lmer(pmax_dAUDPC ~ Isolate*Line + Scorer + 
               (1| Isolate:Batch) + (1|Isolate:Batch:Rep) + (1|Isolate:Batch:Rep:Box) + (1|Isolate:Batch:Rep:Box:Tray), data=dat2);
anova(modd);
summary(modd)

######################################
# Step 4 - generate plots of residuals
######################################
res.lev1 <- resid(modp,scaled=T) 
fit.fix<-predict(modp, re.form=~0) 

par(mfrow=c(2,2))
qqnorm(res.lev1)
plot(res.lev1~fit.fix, main="Fitted-Value Plot") 
hist(res.lev1, main = "Histogram of Residuals")
mtext("lgt_pAUDPC", side = 3, line = -1, outer = TRUE)

######################################
# Step 5 - calculate estimated marginal means
######################################

if(!require(emmeans)){install.packages(emmeans)}
library(emmeans)

emm.modp <- emmeans(modp, specs="Isolate", by="Line", adjust="bonferroni")
est.modp <- as.data.frame(emm.modp) 

emm.modd <- emmeans(modd, specs="Isolate", by="Line", adjust="bonferroni")
est.modd <- as.data.frame(emm.modd) 

######################################
# Step 6 - format and combine the datasets
######################################

trimp <- est.modp[,1:3]
colnames(trimp)[3] <- "emmean_lgt_pAUDPC"
trimd <- ests.modd[,1:3]
colnames(trimd)[3] <- "emmean_dAUDPC"
means <- merge(trimp, trimd)
