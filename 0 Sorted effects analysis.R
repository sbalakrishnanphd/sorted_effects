################################################################################
#### AUTHOR: Siva Balakrishnan
#### UPDATED: Nov 19 2024
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
           "Hmisc", "xtable", "knitr", "gt", "rsample") 

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
# Table 3: CLAN
# Classification analyses
samp_weight <- rep(1, nrow(data_hemo))
ca_hemo     <- ca(fm     = fm1, 
                  data   = data_hemo,
                  method = "ols",
                  var    = "TREAT",
                  u      = 0.2,
                  alpha  = 0.05,
                  seed   = 2024,
                  t      = clan_vars,
                  cl     = "diff"
)

samp_weight <- rep(1, nrow(data_ebf))
ca_ebf      <- ca(fm     = fm2, 
                  data   = data_ebf,
                  method = "logit",
                  var    = "TREAT",
                  u      = 0.2,
                  alpha  = 0.05,
                  seed   = 2024,
                  t      = clan_vars,
                  cl     = "diff"
)

samp_weight <- rep(1, nrow(data_dds))
ca_dds      <- ca(fm     = fm3, 
                  data   = data_dds,
                  method = "ols",
                  var    = "TREAT",
                  u      = 0.2,
                  alpha  = 0.05,
                  seed   = 2024,
                  t      = clan_vars,
                  cl     = "diff"
)

# CLAN Results table
## Stars
hemo_stars <- cut(data.frame(summary(ca_hemo))[,3], 
                  breaks = c(0, 0.001, 0.01, 0.05, 2), 
                  labels = c("***", "**", "*", ""), 
                  right = FALSE)
ebf_stars  <- cut(data.frame(summary(ca_ebf))[,3], 
                  breaks = c(0, 0.001, 0.01, 0.05, 2), 
                  labels = c("***", "**", "*", ""), 
                  right = FALSE)
dds_stars  <- cut(data.frame(summary(ca_dds))[,3], 
                  breaks = c(0, 0.001, 0.01, 0.05, 2), 
                  labels = c("***", "**", "*", ""), 
                  right = FALSE)
## Template
table3a <- read_xlsx("/Users/sivabalakrishnan/Library/CloudStorage/Dropbox/Siva's dissertation/Dissertation paper 1/Analysis/0 datavars_final.xlsx")
## Differences + stars
table3b <- data.frame(Varnames=rownames(summary(ca_hemo)),
                      Hemoglobin=paste(round(data.frame(summary(ca_hemo))$Estimate, 2),
                                       hemo_stars), 
                      EBF=paste(round(data.frame(summary(ca_ebf))$Estimate, 2),
                                ebf_stars),
                      DDS=paste(round(data.frame(summary(ca_dds))$Estimate, 2),
                                dds_stars)
)
## Join to make table 3
table3  <- left_join(table3a, table3b, by=join_by(Varnames == Varnames))              

## Save table 3
table3$Varnames = NULL
gtsave(gt(table3), filename="table3.docx")
rm(ca_hemo, ca_ebf, ca_dds, hemo_stars, ebf_stars, dds_stars, table3a, table3b)
################################################################################
################################################################################
# Sorted partial effects (Figures 1-3)

# Hemoglobin
samp_weight <- rep(1, nrow(data_hemo))
spe_hemo    <- spe1(fm     = fm1,
                    data   = data_hemo,
                    method = "ols",
                    var    = "TREAT",
                    us     = c(1:99)/100,
                    alpha  = 0.1,
                    seed   = 2024,
                    bc     = FALSE,
                    b      = 500
)

# EBF
samp_weight <- rep(1, nrow(data_ebf))
spe_ebf     <- spe1(fm     = fm2,
                    data   = data_ebf,
                    method = "logit",
                    var    = "TREAT",
                    us     = c(1:99)/100,
                    alpha  = 0.1,
                    seed   = 2024,
                    bc     = FALSE,
                    b      = 500
)

# DDS
samp_weight <- rep(1, nrow(data_dds))
spe_dds     <- spe1(fm     = fm3,
                    data   = data_dds,
                    method = "ols",
                    var    = "TREAT",
                    us     = c(1:99)/100,
                    alpha  = 0.1,
                    seed   = 2024,
                    bc     = FALSE,
                    b      = 500
)

# Plots
par(mfrow=c(3,1))
  plot.spe(spe_hemo, ylim=c(-12, 12), main="Hemoglobin", cex=1.5)
  plot.spe(spe_ebf, ylim=c(-0.6, 0.6), main="Exclusive breastfeeding", cex=1.5)
  plot.spe(spe_dds, ylim=c(-2, 2), main="Dietary diversity score", cex=1.5)
#rm(spe_hemo, spe_ebf, spe_dds)
################################################################################
################################################################################
# Table 2

################################################################################
################################################################################





