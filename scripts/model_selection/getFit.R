library(dplyr)
library(lme4)
library(plotly)

# Load mwMod class definition
source("~/Google Drive/Project :: ipmwam/Papers/SivansPaper/DataAnalysis/mwMod.R")

metadata=tbl_df(read.csv("~/Google Drive/Project :: ipmwam/Papers/SivansPaper/Milkweed/DATA SETS/Data Frames/metadata13to15.csv"))

# Where there was a plant in June and NA for height in Sept, make surv 0
metadata <- metadata %>% mutate(surv = ifelse(is.na(surv) & is.na(h_apical.next) & !is.na(h_apical), 0, surv), 
                                munched = ifelse(herb_avg > 0, 1, 0),
                                log_herb_avg = log(0.1+herb_avg))

# Set numerical options for glmer
glmerCtrl <- glmerControl(optimizer = c("bobyqa"), optCtrl = list(maxfun=50000))

# =========
# Flowering
# =========

metadata_usc <- metadata %>% filter(!is.na(h_apical),
                                    !is.na(h_apical.next),
                                    !is.na(herb_avg),
                                    !is.na(fec.flower),
                                    !is.na(surv),
                                    h_apical.next > 50,
                                    !(site == "YTB" & h_apical <= 50))

metadata_sc <- metadata_usc %>% mutate_each(funs(sc = as.numeric(scale(.))), h_apical, log_herb_avg)

cat("Computing flowering fit...")
flower.mdl <- glmer(fec.flower ~ h_apical + h_apical:log_herb_avg - 1 + (1|site/transect)+(h_apical+log_herb_avg|year), data=metadata_sc, nAGQ=1, family=binomial(), control=glmerCtrl)
cat("done!\n")

# Create model class
flower.fit <- mwMod(list(mdl = flower.mdl, 
                         vars = c("h_apical", "log_herb_avg"), 
                         scaled = list(h_apical = scale(metadata_usc$h_apical), 
                                       log_herb_avg = scale(metadata_usc$log_herb_avg))))
# Check parameters
cat("Checking parameters:\n")
checkPars(flower.fit)

# ========
# Survival
# ========

metadata_usc <- metadata %>% filter(!is.na(h_apical),
                                    !is.na(h_apical.next),
                                    !is.na(herb_avg),
                                    fec.flower == 1,
                                    !is.na(surv),
                                    h_apical.next > 50,
                                    !(site == "YTB" & h_apical <= 50))

metadata_sc <- metadata_usc %>% mutate_each(funs(sc = as.numeric(scale(.))), h_apical, log_herb_avg)

cat("Computing survival fit...")
surv.mdl <- glmer(surv ~ 1 + (log_herb_avg|site/transect)+(1|year), data=metadata_sc, family=binomial(), nAGQ=1, control=glmerCtrl)
cat("done!\n")

surv.fit <- mwMod(list(mdl = surv.mdl, 
                       vars = c("h_apical", "log_herb_avg"), 
                       scaled = list(h_apical = scale(metadata_usc$h_apical), 
                                     log_herb_avg = scale(metadata_usc$log_herb_avg))))
# Check parameters
cat("Checking parameters:\n")
checkPars(surv.fit)

# ======
# Growth
# ======

metadata_usc <- metadata %>% filter(!is.na(h_apical),
                                    !is.na(h_apical.next),
                                    !is.na(herb_avg),
                                    fec.flower == 1,
                                    surv == 1,
                                    h_apical.next > 50,
                                    !(site == "YTB" & h_apical <= 50))

metadata_sc <- metadata_usc %>% mutate_each(funs(sc = as.numeric(scale(.))), h_apical, h_apical.next, log_herb_avg)

cat("Computing growth fit...")
growth.mdl <- lmer(h_apical.next ~ h_apical + h_apical:log_herb_avg + (h_apical+log_herb_avg|site/transect)+(h_apical|year), data=metadata_sc, REML=T)
cat("done!\n")

growth.fit <- mwMod(list(mdl = growth.mdl, 
                         vars = c("h_apical", "log_herb_avg"), 
                         scaled = list(h_apical = scale(metadata_usc$h_apical),
                                       h_apical.next = scale(metadata_usc$h_apical.next),
                                       log_herb_avg = scale(metadata_usc$log_herb_avg))))
# Check parameters
cat("Checking parameters:\n")
checkPars(growth.fit)

# ====
# Pods
# ====

metadata_usc <- metadata %>% filter(!is.na(h_apical),
                                    !is.na(h_apical.next),
                                    !is.na(herb_avg),
                                    fec.flower == 1,
                                    surv == 1,
                                    !is.na(N_pods),
                                    h_apical.next > 50,
                                    !(site == "YTB" & h_apical <= 50))

metadata_sc <- metadata_usc %>% mutate_each(funs(sc = as.numeric(scale(.))), h_apical, h_apical.next, log_herb_avg)

cat("Computing pods fit...")
pods.mdl <- glmer(N_pods ~ h_apical.next + log_herb_avg - 1 + (h_apical.next+log_herb_avg|site/transect)+(h_apical.next+log_herb_avg|year), data=metadata_sc, nAGQ=1, family=poisson(), control=glmerCtrl)
cat("done!\n")

pods.fit <- mwMod(list(mdl = pods.mdl, 
                       vars = c("h_apical.next", "log_herb_avg"), 
                       scaled = list(h_apical = scale(metadata_usc$h_apical),
                                     h_apical.next = scale(metadata_usc$h_apical.next),
                                     log_herb_avg = scale(metadata_usc$log_herb_avg))))
# Check parameters
cat("Checking parameters:\n")
checkPars(pods.fit)

save(flower.fit, surv.fit, growth.fit, pods.fit, file = "vitalFits.RData")