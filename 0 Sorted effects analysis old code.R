################################################################################
#### AUTHOR: Siva Balakrishnan
#### UPDATED: Oct 11 2024
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
           "Hmisc", "xtable", "knitr", "gt") 

package.check <- lapply(
  vec.pac,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

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
missing_vars  <- c("mi_BB_BIRTH_WGHT_KG","mi_BB_CSEC","mi_BBVAR_BIRTH_GESWK",      
                  "mi_BL_BF_ATTITUDE_SCORE_ST", "mi_BL_BFKNOW_SCORE_ST",
                  "mi_BL_CONFLICT_ST", "mi_BL_DECISION_ST","mi_BL_HH_MOM_AGE", 
                  "mi_BL_INFO_HOSP", "mi_BL_INFO_ONLINE", "mi_BL_SS_ST",
                  "mi_BLVAR_DASS_PRI_TOT_ST", "mi_BLVAR_HH_DAD_AGE", "mi_BLVAR_HH_DAD_HIEDU",
                  "mi_BLVAR_HH_MOM_HIEDU", "mi_BLVAR_HYGN_PRI_NUM_ST", "mi_HH_MOM_BIRTHHIS",
                  "mi_HHVAR_ASSET", "mi_HHVAR_MOM_MIGHIS"   
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

# Setting all CLAN table 3
data_vars           <- read_xlsx("/Users/sivabalakrishnan/Library/CloudStorage/Dropbox/Siva's dissertation/Dissertation paper 1/Analysis/0 datavars_final.xlsx")
data_vars$Group     <- NULL
data_vars$Subgroup  <- NULL

# Formulas
# Some missing vars cause convergence issues
## EBF
missing_issues_ebf  <- c("mi_BL_INFO_ONLINE", "mi_BL_SS_ST", "mi_BLVAR_DASS_PRI_TOT_ST",
                        "mi_BLVAR_HYGN_PRI_NUM_ST", "mi_BL_CONFLICT_ST", "mi_BL_INFO_HOSP", 
                        "mi_BLVAR_HH_DAD_AGE", "mi_HHVAR_ASSET", "mi_HHVAR_MOM_MIGHIS")
step1_ebf           <- paste(c(control_vars, missing_vars[!missing_vars %in% missing_issues_ebf]), collapse=" + ")
step2_ebf           <- paste("TREAT", step1_ebf,  sep=" + ")
step3_ebf           <- paste("TREAT*(", step1_ebf, ")", sep="")
step4_ebf           <- paste(step2_ebf, step3_ebf, sep=" + "  )
fm1                 <- paste("FLVAR_FEED_EBFCOMB", step4_ebf, sep=" ~ ")

## DDS
missing_issues_dds  <- c("mi_BL_INFO_ONLINE", "mi_BL_SS_ST", "mi_BLVAR_DASS_PRI_TOT_ST",
                         "mi_BLVAR_HYGN_PRI_NUM_ST", "mi_BB_CSEC", "mi_BL_CONFLICT_ST",
                         "mi_BL_DECISION_ST", "mi_BL_INFO_HOSP", "mi_BLVAR_HH_DAD_AGE", 
                         "mi_HHVAR_ASSET", "mi_HHVAR_MOM_MIGHIS")
step1_dds           <- paste(c(control_vars, missing_vars[!missing_vars %in% missing_issues_dds]), collapse=" + ")
step2_dds           <- paste("TREAT", step1_dds,  sep=" + ")
step3_dds           <- paste("TREAT*(", step1_dds, ")", sep="")
step4_dds           <- paste(step2_dds, step3_dds, sep=" + "  )
fm2                 <- paste("FLVAR_FEED_DDS", step4_dds, sep=" ~ ")

## HEMO
missing_issues_hemo <- c("mi_BL_INFO_ONLINE", "mi_BL_SS_ST", "mi_BLVAR_DASS_PRI_TOT_ST",
                         "mi_BLVAR_HYGN_PRI_NUM_ST", "mi_HHVAR_MOM_MIGHIS", "mi_BB_CSEC",
                         "mi_BL_CONFLICT_ST", "mi_BL_INFO_HOSP", "mi_BLVAR_HH_DAD_AGE",
                         "mi_HHVAR_ASSET")
step1_hemo          <- paste(c(control_vars, missing_vars[!missing_vars %in% missing_issues_hemo]), collapse=" + ")
step2_hemo          <- paste("TREAT", step1_hemo,  sep=" + ")
step3_hemo          <- paste("TREAT*(", step1_hemo, ")", sep="")
step4_hemo          <- paste(step2_hemo, step3_hemo, sep=" + "  )
fm3                 <- paste("FLVAR_CH_HEMO", step4_hemo, sep=" ~ ")

################################################################################
###### EBF
################################################################################
# Drop missing EBF
data_ebf <- data[!is.na(data$FLVAR_FEED_EBFCOMB), ]

# GLM check for EBF
modela <- glm(fm1, data=data_ebf, family=binomial(link="logit"))
summary(modela)

# Some parameters
samp_weight <- rep(1, nrow(data_ebf))

# SPE
spe_ebf <- spe1( fm  = fm1,
                data   = data_ebf,
                method = "logit",
                var    = "TREAT",
                us     = c(1:99)/100,
                alpha  = 0.1,
                seed   = 2024,
                bc     = TRUE,
                b      = 3000
                )

# Classification analysis
ca_ebf1 <- ca(fm     = fm1, 
              data   = data_ebf,
              method = "logit",
              var    = "TREAT",
              u      = 0.2,
              alpha  = 0.05,
              seed   = 2024,
              t      = clan_vars,
              cl     = "both"
)

ca_ebf2 <- ca(fm     = fm1, 
              data   = data_ebf,
              method = "logit",
              var    = "TREAT",
              u      = 0.2,
              alpha  = 0.05,
              seed   = 2024,
              t      = clan_vars,
              cl     = "diff"
)

# EBF CLAN Results table
table_ebf                    <- data.frame(matrix(ncol=5, nrow=length(clan_vars)))
colnames(table_ebf)          <- c("Varnames", "Most", "Least", "Difference", "Adjusted p-value")
table_ebf$Varnames           <- rownames(summary(ca_ebf1))
table_ebf$Most               <- round(data.frame(summary(ca_ebf1))$Most, 2)
table_ebf$Least              <- round(data.frame(summary(ca_ebf1))$Least, 2)
table_ebf$Difference         <- round(data.frame(summary(ca_ebf2))$Estimate, 2)
table_ebf$`Adjusted p-value` <- data.frame(summary(ca_ebf2))[,3]

################################################################################
###### DDS
################################################################################
# Drop missing DDS
data_dds <- data[!is.na(data$FLVAR_FEED_DDS), ]

# GLM check for DDS
modelb <- glm(fm2, data=data_dds)
summary(modelb)

# Some parameters
samp_weight <- rep(1, nrow(data_dds))

# SPE
spe_dds <- spe1( fm     = fm2,
                 data   = data_dds,
                 method = "ols",
                 var    = "TREAT",
                 us     = c(1:99)/100,
                 alpha  = 0.1,
                 seed   = 2024,
                 bc     = TRUE,
                 b      = 3000
)

# Classification analysis
ca_dds1 <- ca(fm     = fm2, 
              data   = data_dds,
              method = "ols",
              var    = "TREAT",
              u      = 0.2,
              alpha  = 0.05,
              seed   = 2024,
              t      = clan_vars,
              cl     = "both"
)

ca_dds2 <- ca(fm     = fm2, 
              data   = data_dds,
              method = "ols",
              var    = "TREAT",
              u      = 0.2,
              alpha  = 0.05,
              seed   = 2024,
              t      = clan_vars,
              cl     = "diff"
)

# DDS CLAN Results table
table_dds                    <- data.frame(matrix(ncol=5, nrow=length(clan_vars)))
colnames(table_dds)          <- c("Varnames", "Most", "Least", "Difference", "Adjusted p-value")
table_dds$Varnames           <- rownames(summary(ca_dds1))
table_dds$Most               <- round(data.frame(summary(ca_dds1))$Most, 2)
table_dds$Least              <- round(data.frame(summary(ca_dds1))$Least, 2)
table_dds$Difference         <- round(data.frame(summary(ca_dds2))$Estimate, 2)
table_dds$`Adjusted p-value` <- data.frame(summary(ca_dds2))[,3]

################################################################################
###### Hemoglobin
################################################################################
# Drop missing Hemo
data_hemo <- data[!is.na(data$FLVAR_CH_HEMO), ]

# GLM check for EBF
modelc <- glm(fm3, data=data_hemo)
summary(modelc)

# Some parameters
samp_weight <- rep(1, nrow(data_hemo))

# SPE
spe_hemo <- spe1(fm     = fm3,
                 data   = data_hemo,
                 method = "ols",
                 var    = "TREAT",
                 us     = c(1:99)/100,
                 alpha  = 0.1,
                 seed   = 2024,
                 bc     = TRUE,
                 b      = 3000
)

# Classification analysis
ca_hemo1 <- ca(fm     = fm3, 
              data   = data_hemo,
              method = "ols",
              var    = "TREAT",
              u      = 0.2,
              alpha  = 0.05,
              seed   = 2024,
              t      = clan_vars,
              cl     = "both"
)

ca_hemo2 <- ca(fm     = fm3, 
              data   = data_hemo,
              method = "ols",
              var    = "TREAT",
              u      = 0.2,
              alpha  = 0.05,
              seed   = 2024,
              t      = clan_vars,
              cl     = "diff"
)


# Result table
table_hemo <- data.frame(matrix(ncol=5, nrow=length(clan_vars)))
colnames(table_hemo) <- c("Varnames", "Most", "Least", "Difference", "Adjusted p-value")
table_hemo$Varnames           <- rownames(summary(ca_hemo1))
table_hemo$Most               <- round(data.frame(summary(ca_hemo1))$Most, 2)
table_hemo$Least              <- round(data.frame(summary(ca_hemo1))$Least, 2)
table_hemo$Difference         <- round(data.frame(summary(ca_hemo2))$Estimate, 2)
table_hemo$`Adjusted p-value` <- data.frame(summary(ca_hemo2))[,3]

################################################################################
##### CLAN tables
################################################################################
# Adding names
table_ebf_fin <- left_join(data_vars, table_ebf, by = join_by(Varnames == Varnames))
table_ebf_fin$Varnames <- NULL
table_ebf_fin$Variables.y <- NULL
colnames(table_ebf_fin) <- c("Level", "Construct", "Variable", "Most", "Least", "Difference", "Adjusted p-value")

table_dds_fin <- left_join(data_vars, table_dds, by = join_by(Varnames == Varnames))
table_dds_fin$Varnames <- NULL
table_dds_fin$Variables.y <- NULL
colnames(table_dds_fin) <- c("Level", "Construct", "Variable", "Most", "Least", "Difference", "Adjusted p-value")

table_hemo_fin <- left_join(data_vars, table_hemo, by = join_by(Varnames == Varnames))
table_hemo_fin$Varnames <- NULL
table_hemo_fin$Variables.y <- NULL
colnames(table_hemo_fin) <- c("Level", "Construct", "Variable", "Most", "Least", "Difference", "Adjusted p-value")

# Save CLAN tables
gtsave(gt(table_ebf_fin), filename="table3_ebf.docx")
gtsave(gt(table_dds_fin), filename="table4_dds.docx")
gtsave(gt(table_hemo_fin), filename="table5_hemo.docx")

################################################################################
##### Sorted partial effects (FIGURES 1-3)
################################################################################

### Plots
#plot.spe(spe_ebf, ylim=c(-0.6, 0.6), main="Exclusive breastfeeding", cex=1.5)
#plot.spe(spe_dds, ylim=c(-2, 2), main="Dietary diversity score", cex=1.5)
#plot.spe(spe_hemo, ylim=c(-12, 12), main="Hemoglobin", cex=1.5)

################################################################################
##### Gates analysis (TABLE 2)
################################################################################

#("FLVAR_FEED_EBFCOMB", "FLVAR_FEED_DDS", "FLVAR_CH_HEMO")

# EBF
ebf           <- data.frame(spe_ebf$individual)
ebf$outcome   <- data_ebf$FLVAR_FEED_EBFCOMB
ebf$treat     <- data_ebf$TREAT
ebf$quintile  <- ntile(ebf$spe_ebf.individual, 5)
ebf_top       <- ebf[ebf$quintile==5,]
ebf_bot       <- ebf[ebf$quintile==1,]
ebf_t         <- t.test(ebf_top$spe_ebf.individual, ebf_bot$spe_ebf.individual)

table_gates           <- data.frame(matrix(ncol=8, nrow=3))
colnames(table_gates) <- c("Outcome",  "Mean in treated", "Mean in control", "Top ATE", "Mean in treated2", "Mean in control2", "Bot ATE", "P-value for difference in ATEs")
table_gates[1,]       <- c("EBF (Yes/No)", 
                           round(mean(ebf_top[ebf_top$treat==1,]$outcome), 2),
                           round(mean(ebf_top[ebf_top$treat==0,]$outcome), 2),
                           round(exp(as.numeric(ebf_t$estimate[1])), 2),
                           round(mean(ebf_bot[ebf_bot$treat==1,]$outcome), 2),
                           round(mean(ebf_bot[ebf_bot$treat==0,]$outcome), 2),  
                           round(exp(as.numeric(ebf_t$estimate[2])), 2),
                           round(as.numeric(ebf_t$p.value), 3)
                          )

# DDS
dds             <- data.frame(spe_dds$individual)
dds$outcome     <- data_dds$FLVAR_FEED_DDS      
dds$quintile    <- ntile(dds$spe_dds.individual, 5)
dds$treat       <- data_dds$TREAT
dds_top         <- dds[dds$quintile==5,]
dds_bot         <- dds[dds$quintile==1,]
dds_t           <- t.test(dds_top$spe_dds.individual, dds_bot$spe_dds.individual)
table_gates[2,]       <- c("DDS (0-12)", 
                           round(mean(dds_top[dds_top$treat==1,]$outcome), 2),
                           round(mean(dds_top[dds_top$treat==0,]$outcome), 2),
                           round(as.numeric(dds_t$estimate[1]), 2),
                           round(mean(dds_bot[dds_bot$treat==1,]$outcome), 2),
                           round(mean(dds_bot[dds_bot$treat==0,]$outcome), 2),  
                           round(as.numeric(dds_t$estimate[2]), 2),
                           round(as.numeric(dds_t$p.value), 3)
)

# HEMO
hemo            <- data.frame(spe_hemo$individual)
hemo$outcome    <- data_hemo$FLVAR_CH_HEMO
hemo$quintile   <- ntile(hemo$spe_hemo.individual, 5)
hemo$treat      <- data_hemo$TREAT
hemo_top        <- hemo[hemo$quintile==5,]
hemo_bot        <- hemo[hemo$quintile==1,]
hemo_t          <- t.test(hemo_top$spe_hemo.individual, hemo_bot$spe_hemo.individual)
table_gates[3,]       <- c("Hemoglobin (g/L)", 
                           round(mean(hemo_top[hemo_top$treat==1,]$outcome), 2),
                           round(mean(hemo_top[hemo_top$treat==0,]$outcome), 2),
                           round(as.numeric(hemo_t$estimate[1]), 2),
                           round(mean(hemo_bot[hemo_bot$treat==1,]$outcome), 2),
                           round(mean(hemo_bot[hemo_bot$treat==0,]$outcome), 2),  
                           round(as.numeric(hemo_t$estimate[2]), 2),
                           round(as.numeric(hemo_t$p.value), 3)
)

# Table 2
table_gates$`P-value for difference in ATEs` <- c("<0.001", "<0.001", "<0.001")

gtsave(gt(table_gates), filename="table2_gates.docx")


# 90% CI
####
ebf_top <- c(exp(mean(spe_ebf$spe$est_bc [80:99])), 
             exp(mean(spe_ebf$spe$lbound_est_bc [80:99])),
             exp(mean(spe_ebf$spe$ubound_est_bc [80:99]))
             )

ebf_bot <- c(exp(mean(spe_ebf$spe$est_bc [1:20])), 
             exp(mean(spe_ebf$spe$lbound_est_bc [1:20])),
             exp(mean(spe_ebf$spe$ubound_est_bc [1:20]))
)
####
dds_top <- c(mean(spe_dds$spe$est_bc [80:99]), 
             mean(spe_dds$spe$lbound_est_bc [80:99]),
             mean(spe_dds$spe$ubound_est_bc [80:99])
)

dds_bot <- c(mean(spe_dds$spe$est_bc [1:20]), 
             mean(spe_dds$spe$lbound_est_bc [1:20]),
             mean(spe_dds$spe$ubound_est_bc [1:20])
)
####
hemo_top <- c(mean(spe_hemo$spe$est_bc [80:99]), 
             mean(spe_hemo$spe$lbound_est_bc [80:99]),
             mean(spe_hemo$spe$ubound_est_bc [80:99])
)

hemo_bot <- c(mean(spe_hemo$spe$est_bc [1:20]), 
             mean(spe_hemo$spe$lbound_est_bc [1:20]),
             mean(spe_hemo$spe$ubound_est_bc [1:20])
)
################################################################################ 
################################################################################
################################################################################ 
################################################################################

# top and bottom 20%
# EBF
samp_weight <- rep(1, nrow(data_ebf))
spe_ebf_5 <- spe1( fm  = fm1,
                 data   = data_ebf,
                 method = "logit",
                 var    = "TREAT",
                 us     = c(1:4)/5,
                 alpha  = 0.1,
                 seed   = 2024,
                 bc     = TRUE,
                 b      = 3000
)
# DDS
samp_weight <- rep(1, nrow(data_dds))
spe_dds_5 <- spe1( fm     = fm2,
                 data   = data_dds,
                 method = "ols",
                 var    = "TREAT",
                 us     = c(1:4)/5,
                 alpha  = 0.1,
                 seed   = 2024,
                 bc     = TRUE,
                 b      = 3000
)
# HEMO
samp_weight <- rep(1, nrow(data_hemo))
spe_hemo_5 <- spe1(fm     = fm3,
                 data   = data_hemo,
                 method = "ols",
                 var    = "TREAT",
                 us     = c(1:4)/5,
                 alpha  = 0.1,
                 seed   = 2024,
                 bc     = TRUE,
                 b      = 3000
)

##################

# top and bottom 20%
# EBF
samp_weight <- rep(1, nrow(data_ebf))
spe_ebf_10 <- spe1( fm  = fm1,
                   data   = data_ebf,
                   method = "logit",
                   var    = "TREAT",
                   us     = c(1:9)/10,
                   alpha  = 0.1,
                   seed   = 2024,
                   bc     = TRUE,
                   b      = 3000
)
# DDS
samp_weight <- rep(1, nrow(data_dds))
spe_dds_10 <- spe1( fm     = fm2,
                   data   = data_dds,
                   method = "ols",
                   var    = "TREAT",
                   us     = c(1:9)/10,
                   alpha  = 0.1,
                   seed   = 2024,
                   bc     = TRUE,
                   b      = 3000
)
# HEMO
samp_weight <- rep(1, nrow(data_hemo))
spe_hemo_10 <- spe1(fm     = fm3,
                   data   = data_hemo,
                   method = "ols",
                   var    = "TREAT",
                   us     = c(1:9)/10,
                   alpha  = 0.1,
                   seed   = 2024,
                   bc     = TRUE,
                   b      = 3000
)

################################################################################
################################################################################
# Gates analysis (Table 2)
## HEMO
samp_weight   <- rep(1, nrow(data_hemo))
spe_hemo_1    <- spe1(fm     = fm1,
                      data   = data_hemo,
                      method = "ols",
                      var    = "TREAT",
                      us     = c(1:99)/100,
                      alpha  = 0.1,
                      seed   = 2024,
                      bc     = TRUE,
                      b      = 500
)
## EBF
samp_weight   <- rep(1, nrow(data_ebf))
spe_ebf_1     <- spe1(fm     = fm2,
                      data   = data_ebf,
                      method = "logit",
                      var    = "TREAT",
                      us     = c(1:99)/100,
                      alpha  = 0.1,
                      seed   = 2024,
                      bc     = TRUE,
                      b      = 500
)
## DDS
samp_weight   <- rep(1, nrow(data_dds))
spe_dds_1     <- spe1(fm     = fm3,
                      data   = data_dds,
                      method = "ols",
                      var    = "TREAT",
                      us     = c(1:99)/100,
                      alpha  = 0.1,
                      seed   = 2024,
                      bc     = TRUE,
                      b      = 500
)
# us     = c(1:(nrow(data_hemo)-1))/nrow(data_hemo),

table2a <- data.frame(var=c("Hemo", "EBF", "DDS"),
                      #Bottom
                      bot=c(mean(spe_hemo_1$spe$est_bc[names(spe_hemo_1$spe$est_bc)<=20]),
                            mean(spe_ebf_1$spe$est_bc[names(spe_ebf_1$spe$est_bc)<=20]),
                            mean(spe_dds_1$spe$est_bc[names(spe_dds_1$spe$est_bc)<=20])),
                      botl=c(mean(spe_hemo_1$spe$lbound_est_bc[names(spe_hemo_1$spe$lbound_est_bc)<=20]),
                             mean(spe_ebf_1$spe$lbound_est_bc[names(spe_ebf_1$spe$lbound_est_bc)<=20]),
                             mean(spe_dds_1$spe$lbound_est_bc[names(spe_dds_1$spe$lbound_est_bc)<=20])),
                      botu=c(mean(spe_hemo_1$spe$ubound_est_bc[names(spe_hemo_1$spe$ubound_est_bc)<=20]),
                             mean(spe_ebf_1$spe$ubound_est_bc[names(spe_ebf_1$spe$ubound_est_bc)<=20]),
                             mean(spe_dds_1$spe$ubound_est_bc[names(spe_dds_1$spe$ubound_est_bc)<=20])),
                      #Top
                      top=c(mean(spe_hemo_1$spe$est_bc[names(spe_hemo_1$spe$est_bc)>=80]),
                            mean(spe_ebf_1$spe$est_bc[names(spe_ebf_1$spe$est_bc)>=80]),
                            mean(spe_dds_1$spe$est_bc[names(spe_dds_1$spe$est_bc)>=80])),
                      topl=c(mean(spe_hemo_1$spe$lbound_est_bc[names(spe_hemo_1$spe$lbound_est_bc)>=80]),
                             mean(spe_ebf_1$spe$lbound_est_bc[names(spe_ebf_1$spe$lbound_est_bc)>=80]),
                             mean(spe_dds_1$spe$lbound_est_bc[names(spe_dds_1$spe$lbound_est_bc)>=80])),
                      topu=c(mean(spe_hemo_1$spe$ubound_est_bc[names(spe_hemo_1$spe$ubound_est_bc)>=80]),
                             mean(spe_ebf_1$spe$ubound_est_bc[names(spe_ebf_1$spe$ubound_est_bc)>=80]),
                             mean(spe_dds_1$spe$ubound_est_bc[names(spe_dds_1$spe$ubound_est_bc)>=80]))
)

################################################################################
################################################################################

# Table 2
top1 <- t.test(spe_hemo$individual[spe_hemo$individual>= quantile(spe_hemo$individual, 0.8)])
top2 <- t.test(spe_ebf$individual[spe_ebf$individual>= quantile(spe_ebf$individual, 0.8)])
top3 <- t.test(spe_dds$individual[spe_dds$individual>= quantile(spe_dds$individual, 0.8)])

bot1 <- t.test(spe_hemo$individual[spe_hemo$individual<= quantile(spe_hemo$individual, 0.2)])
bot2 <- t.test(spe_ebf$individual[spe_ebf$individual<= quantile(spe_ebf$individual, 0.2)])
bot3 <- t.test(spe_dds$individual[spe_dds$individual<= quantile(spe_dds$individual, 0.2)])


diff1 <- t.test(spe_hemo$individual[spe_hemo$individual>= quantile(spe_hemo$individual, 0.8)],
                spe_hemo$individual[spe_hemo$individual<= quantile(spe_hemo$individual, 0.2)])
diff2 <- t.test(spe_ebf$individual[spe_ebf$individual>= quantile(spe_ebf$individual, 0.8)],
                spe_ebf$individual[spe_ebf$individual<= quantile(spe_ebf$individual, 0.2)])
diff3 <- t.test(spe_dds$individual[spe_dds$individual>= quantile(spe_dds$individual, 0.8)],
                spe_dds$individual[spe_dds$individual<= quantile(spe_dds$individual, 0.2)])

table2 <- data.frame(Variables=c("Hemoglobin", "EBF", "DDS"),
                     
                     Most=as.numeric(round(c(top1$estimate, exp(top2$estimate), top3$estimate),2)),
                     Most_l=as.numeric(round(c(top1$conf.int[1], exp(top2$conf.int[1]), top3$conf.int[1]),2)),
                     Most_u=as.numeric(round(c(top1$conf.int[2], exp(top2$conf.int[2]), top3$conf.int[2]),2)),
                     
                     Least=as.numeric(round(c(bot1$estimate, exp(bot2$estimate), bot3$estimate),2)),
                     Least_l=as.numeric(round(c(bot1$conf.int[1], exp(bot2$conf.int[1]), bot3$conf.int[1]),2)),
                     Least_u=as.numeric(round(c(bot1$conf.int[2], exp(bot2$conf.int[2]), bot3$conf.int[2]),2)),
                     
                     Diff=as.numeric(round(c(diff1$estimate[1]-diff1$estimate[2], exp(diff2$estimate[1]-diff2$estimate[2]), diff3$estimate[1]-diff3$estimate[2]),2)),
                     Diff_l=as.numeric(round(c(diff1$conf.int[1], exp(diff2$conf.int[1]), diff3$conf.int[1]),2)),
                     Diff_u=as.numeric(round(c(diff1$conf.int[2], exp(diff2$conf.int[2]), diff3$conf.int[2]),2))
)
