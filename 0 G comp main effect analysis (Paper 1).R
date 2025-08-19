################################################################################
#### AUTHOR: Siva Balakrishnan
#### UPDATED: Oct 21 2024
################################################################################
######################## SORTED EFFECTS ANALYSIS ############################### 
################################################################################
# Libraries
rm(list = ls())

# Installing (if needed) and setting libraries 
vec.pac= c("GenericML", "labelled", "hdm", "readxl", "openxlsx", 
           'haven', 'writexl',"here",  'caret', "gbm", "randomForest", "rpart", "grf",
           "dplyr","fastDummies","multcomp",'lfe',"doParallel", 'parallel','matrixStats', 
           'MLeval', "MLmetrics",'tibble', 'tidyr','rmarkdown', 'flextable', 'table1', 'psych', 
           'lsr', 'ragg', "pbapply", "boot", "stats", "labelled", "graphics", "SortedEffects", 
           "Hmisc", "xtable", "knitr", "gt", "rsample", "ggplot2") 

lapply(vec.pac,
       FUN = function(x) {
         if (!require(x, character.only = TRUE)) {
           install.packages(x, dependencies = TRUE)
           library(x, character.only = TRUE)
         }
       }
)
rm(vec.pac)

# Get edited functions
source("/Users/sivabalakrishnan/Library/CloudStorage/Dropbox/Siva's dissertation/Dissertation paper 1/Analysis/wtd.quantile hmisc function.R")
source("/Users/sivabalakrishnan/Library/CloudStorage/Dropbox/Siva's dissertation/Dissertation paper 1/Analysis/peestimate function.R")
source("/Users/sivabalakrishnan/Library/CloudStorage/Dropbox/Siva's dissertation/Dissertation paper 1/Analysis/New SPE function.R")

################################################################################
################################################################################
# Dataset
setwd("/Users/sivabalakrishnan/Library/CloudStorage/Dropbox/Siva's dissertation/Dissertation paper 1")
data <- read_dta("/Users/sivabalakrishnan/Library/CloudStorage/Dropbox/Siva's dissertation/0 Master dataset/MASTER.dta")
data <- as.data.frame(lapply(data, as.numeric))

# All Outcomes: "FLVAR_FEED_EBFCOMB", "FLVAR_FEED_DDS", "FLVAR_CH_HEMO"
# All Exposures: data$TREAT
# Heterogeneity variables
control_vars  <- c("BLVAR_BORN","BB_FEMALE","FLVAR_AGE_MON","FLVAR_AGE_MON_SQ",
                   "BB_CSEC","BBVAR_BIRTH_GESWK","BB_BIRTH_WGHT_KG",
                   "BL_BFKNOW_SCORE_ST","BL_CONFLICT_ST","BL_BF_ATTITUDE_SCORE_ST",
                   "BL_SS_ST","BL_DECISION_ST","BLVAR_DASS_PRI_TOT_ST","BL_INFO_HOSP",
                   "BL_INFO_ONLINE","BL_HH_MOM_AGE","BLVAR_HH_MOM_HIEDU","HHVAR_MOM_MIGHIS",
                   "HH_MOM_BIRTHHIS","BLVAR_HYGN_PRI_NUM_ST","BL_HH_SIZE","BLVAR_HH_TOWN",
                   "HHVAR_ASSET","BLVAR_HH_DAD_AGE","BLVAR_HH_DAD_HIEDU","BL_HH_DAD_LIVE",
                   "BLVAR_HH_GRND","BLVAR_HH_SCARE_ANY","COUNTY2","COUNTY3","COUNTY4" 
)

clan_vars     <- c("BLVAR_BORN","BB_FEMALE","FLVAR_AGE_MON",
                   "BB_CSEC","BBVAR_BIRTH_GESWK","BB_BIRTH_WGHT_KG",
                   "BL_BFKNOW_SCORE_ST","BL_CONFLICT_ST","BL_BF_ATTITUDE_SCORE_ST",
                   "BL_SS_ST","BL_DECISION_ST","BLVAR_DASS_PRI_TOT_ST","BL_INFO_HOSP",
                   "BL_INFO_ONLINE","BL_HH_MOM_AGE","BLVAR_HH_MOM_HIEDU","HHVAR_MOM_MIGHIS",
                   "HH_MOM_BIRTHHIS","BLVAR_HYGN_PRI_NUM_ST","BL_HH_SIZE","BLVAR_HH_TOWN",
                   "HHVAR_ASSET","BLVAR_HH_DAD_AGE","BLVAR_HH_DAD_HIEDU","BL_HH_DAD_LIVE",
                   "BLVAR_HH_GRND","BLVAR_HH_SCARE_ANY","COUNTY1","COUNTY2","COUNTY3","COUNTY4"
)     

# Formulas (some missing vars cause convergence issues so remove them)
## HEMO
step1_hemo          <- paste(control_vars, collapse=" + ")
step2_hemo          <- paste("TREAT", step1_hemo,  sep=" + ")
step3_hemo          <- paste("TREAT*(", step1_hemo, ")", sep="")
step4_hemo          <- paste(step2_hemo, step3_hemo, sep=" + "  )
fm1                 <- paste("FLVAR_CH_HEMO", step4_hemo, sep=" ~ ")
rm(step1_hemo, step2_hemo, step3_hemo, step4_hemo)

## EBF
step1_ebf           <- paste(control_vars, collapse=" + ")
step2_ebf           <- paste("TREAT", step1_ebf,  sep=" + ")
step3_ebf           <- paste("TREAT*(", step1_ebf, ")", sep="")
step4_ebf           <- paste(step2_ebf, step3_ebf, sep=" + "  )
fm2                 <- paste("FLVAR_FEED_EBFCOMB", step4_ebf, sep=" ~ ")
rm(step1_ebf, step2_ebf, step3_ebf, step4_ebf)

## DDS
step1_dds           <- paste(control_vars, collapse=" + ")
step2_dds           <- paste("TREAT", step1_dds,  sep=" + ")
step3_dds           <- paste("TREAT*(", step1_dds, ")", sep="")
step4_dds           <- paste(step2_dds, step3_dds, sep=" + "  )
fm3                 <- paste("FLVAR_FEED_DDS", step4_dds, sep=" ~ ")
rm(step1_dds, step2_dds, step3_dds, step4_dds)

# Drop missing 
data_ebf  <- data[!is.na(data$FLVAR_FEED_EBFCOMB), ]
data_dds  <- data[!is.na(data$FLVAR_FEED_DDS), ]
data_hemo <- data[!is.na(data$FLVAR_CH_HEMO), ]

# Test Models
set.seed(2024)
modela <- glm(fm1, data=data_hemo)
modelb <- glm(fm2, data=data_ebf, family=binomial(link="logit"))
modelc <- glm(fm3, data=data_dds)
rm(modela, modelb, modelc)
################################################################################
################################################################################
# G-comp + bootstrapping

## Funtion 
bootgcomp <- function(seed, fm, data, family, fm_rhs, bootnum, alpha, us){
  set.seed(seed)
  # Gcomp function and main estimates
  gcomp <- function(fm, data_model, family, data_org, fm_rhs){
    # Model
    model       <- glm(formula=fm, data=data_model, family=family)
    data1       <- data2 <- data_org
    data1$TREAT <- 1
    data2$TREAT <- 0
    y1          <- predict(model, newdata=data1[,fm_rhs])
    y0          <- predict(model, newdata=data2[,fm_rhs])
    # Estimates
    return(y1 - y0)
  }
  est <- gcomp(fm=fm, data_model=data, family=family, data_org=data, fm_rhs=fm_rhs) 
  # Boostrapping dataset 500 times
  D        <- data %>% nest(-CLUSTER)
  bs       <- bootstraps(D, times = bootnum)  
  bs_final <- list()
  for(i in 1:bootnum){
    as_tibble(bs$splits[[i]]) %>% arrange(CLUSTER)
    bs_final[[i]] <- unnest(as_tibble(bs$splits[[i]]))
  }
  # Running gcomp on bootstrapped datasets
  bs_est <- data.frame(matrix(nrow=nrow(data), ncol=bootnum))
  for(i in 1:bootnum){
    bs_est[,i] <- gcomp(fm=fm, data_model=bs_final[[i]], family=family, data_org=data, fm_rhs=fm_rhs)
  }
  ci_low1 <- ci_high1 <- se <- sigma <- t_hat <- ci_low3 <- ci_high3 <- c()
  # Getting std errors for each row
  for(a in 1:nrow(data)){
    # 1st method CI
    ci_low1[a]  <- quantile(as.numeric(bs_est[a,]), probs=alpha/2)
    ci_high1[a] <- quantile(as.numeric(bs_est[a,]), probs=1-(alpha/2))
    # 2nd method CI
    se[a]       <- sqrt(sum((as.numeric(bs_est[a,])-mean(as.numeric(bs_est[a,])))^2)/(bootnum-1))
    # 3rd method CI (Bundle function)
    z_bs        <- as.numeric(bs_est[a,]) - est[a]
    sigma[a]    <- (quantile(z_bs, 0.75) - quantile(z_bs, 0.25))/(qnorm(0.75)-qnorm(.25))
    t_hat[a]    <- max(abs(z_bs /sigma[a]))
    crt         <- quantile(t_hat, 1 - alpha)
    ci_low3[a]  <- sort(est[a] - crt * sigma[a])
    ci_high3[a] <- sort(est[a] + crt * sigma[a])
  }
  # Using std errors to calculate confidence intervals (2nd method)
  ci_low2  <- est - qnorm(1-(alpha/2))*se
  ci_high2 <- est + qnorm(1-(alpha/2))*se
  result   <- data.frame(est, se, ci_low1, ci_high1, ci_low2, ci_high2, ci_low3, ci_high3)
  # APE
  ape_se   <- sqrt(mean(se^2)) 
  ape      <- c(mean(est), mean(est)-qnorm(1-(alpha/2))*ape_se, mean(est)+qnorm(1-(alpha/2))*ape_se)
  # Quantiles
  result_quantiles <- data.frame(quantiles=us,
                                 est =quantile(est, us, names=FALSE),
                                 low =quantile(ci_low2, us, names=FALSE),
                                 high=quantile(ci_high2, us, names=FALSE)
                                 )
  # Top and bottom
  ttest  <- t.test(est[est>=quantile(est, 0.8)], est[est<=quantile(est, 0.2)])
  topbot  <- c(mean(est[est>=quantile(est, 0.8)]),
               mean(est[est<=quantile(est, 0.2)]),
               mean(est[est>=quantile(est, 0.8)])-mean(est[est<=quantile(est, 0.2)]),
               ttest$p.value
              )
return(list(result=result, ape=ape, result_quantiles=result_quantiles, alpha=alpha, topbot=topbot))
}

## Results
hemo <- bootgcomp(seed=2024, fm=fm1, data=data_hemo, family=gaussian, fm_rhs=c("TREAT", control_vars), 
                  bootnum=500, alpha=0.1, us=c(1:99)/100)
ebf  <- bootgcomp(seed=2024, fm=fm2, data=data_ebf, family=binomial(link="logit"), 
                  fm_rhs=c("TREAT", control_vars), bootnum=500, alpha=0.1, us=c(1:99)/100)
dds  <- bootgcomp(seed=2024, fm=fm3, data=data_dds, family=gaussian, fm_rhs=c("TREAT", control_vars), 
                  bootnum=500, alpha=0.1, us=c(1:99)/100)

# SPE plots
## Function
plot.sorted <- function(x, ylim = NULL, main = NULL, sub=NULL, xlab = "Percentile Index", 
                        ylab = "Sorted Effects", ...){
SE     <- x[[3]]
AE     <- x[[2]]
us     <- unlist(x[[3]][1])
alpha  <- x[[4]]
xlim   <- range(x[[3]][1])
#print(SE, AE, us, alpha, xlim)
# Plot function
plot(us, SE[,2] , type = "l", xlim, ylim, log = "", main, sub, xlab, ylab,
     col = 4, lwd = 2)
polygon(c(us, rev(us)),c(SE[,3], rev(SE[,4])), density = 60, border = F,
        col = 'light blue', lty = 1, lwd = 1)
lines(us, SE[,2], lwd = 2, col = 4 )
# APE and CI
abline(h = AE[[1]], col = 1)
abline(h = AE[[2]], col = 1, lty = 2)
abline(h = AE[[3]], col = 1, lty = 2)
points(c(min(us),.2,.4,.6,.8,max(us)), rep(AE[[1]], 6), col = 1, pch = 15)
# Legend
legend(x = "topleft", col = c(4, 1, "light blue", 1), lwd = c(1, 1, 5, 1),
       lty = c(1, 1, 1, 2), pch = c(NA, 15, NA, NA), pt.cex = c(2, 1),
       bty = 'n', legend = c("SPE","APE", paste0((1 - alpha)*100,"% CB(SPE)"),
                             paste0((1 - alpha)*100,"% CB(APE)")))

}

# Table 2 (GATEs)
# Top and bottom 20% ATEs with propagated errors (root of the mean of the square of st error)
table2           <- data.frame(matrix(nrow=3, ncol=5))
colnames(table2) <- c("Variable", "ATE in top 20%", "ATE in bottom 20%", "Difference", "p-value")
table2[1,]       <- c("Hemoglobin (g/L)", round(hemo$topbot,2)) 
table2[2,]       <- c("EBF (1/0)", round(ebf$topbot,2)) 
table2[3,]       <- c("DDS (0-12)", round(dds$topbot,2)) 


















# TEST area
hemo <- hemo[order(hemo$est),]
ggplot(data=hemoq, aes(x=c(1:nrow(hemoq)), y=est)) +
  geom_ribbon(aes(ymin=low, ymax=high, pattern="stripe"), color="light blue", fill="light blue") +
  geom_line(color="deepskyblue3") + 
  xlab("Percentile") + 
  ylab("Treatment effect on hemoglobin") + 
  geom_line(y=mean(hemo$est)) + 
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

