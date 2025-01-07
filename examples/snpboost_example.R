rm(list=ls())

library(data.table)
library(tidyverse)
library(parallel)
library(bigsnpr)
library(ggplot2)
library(cowplot)
library(readr)

### load snpboost functions
source("/daten/epi/PRS/snpnet/SNP_boost/functions/Clean/lm_boost.R")
source("/daten/epi/PRS/snpnet/SNP_boost/functions/Clean/snpboost.R")
source("/daten/epi/PRS/snpnet/SNP_boost/functions/Clean/functions_snpboost.R")
source("/daten/epi/PRS/snpnet/SNP_boost/functions/Clean/functions_snpnet.R")

genotype.pfile <- "/daten/epi/PRS/UKBB_genotyped/ukb_cal_allChrs_XY" ### genotype file, pgen format, .pvar.zst
phenotype.file <- "/daten/epi/PRS/snpnet/non_imputed_data/phenotypes/phenotypes_unrelated.phe" ### phenotype file, .phe, tabular data
phenotype <- "LDL_adj" ### name of phenotype column in phenotype file
covariates <- NULL ### vector of covariates name if any (used during the fitting)


configs <- list(plink2.path = "/daten/epi/PRS/GWAS/Tools/plink2", ## path to plink file
                zstdcat.path = "/home/klinkhammer/anaconda3/bin/zstdcat", ## path to zstdcat
                results.dir = paste0("/daten/epi/PRS/snpnet/SNP_boost/functions/test_",phenotype), ## results folder
                save = FALSE, ## if true, results are saved after each batch
                prevIter = 0,
                missing.rate = 0.1,
                MAF.thresh = 0.001,
                num.snps.batch = 1000,  # not important here
                early.stopping = TRUE, # must be true to stop with validation set
                stopping.lag = 2, # stop after two batches
                verbose = TRUE,
                mem=16000,
                niter=30000, # not important here
                nCores=16,
                standardize.variant=TRUE
)


###fit fit_snpboost
time_snpboost_start <- Sys.time()

fit_snpboost <- snpboost(genotype.pfile = genotype.pfile,
                     phenotype.file = phenotype.file,
                     phenotype = phenotype,
                     covariates = covariates,
                     configs = configs,
                     split.col = "train_test", ### name of column in phenotype file indication train, val, test set (named "train", "val", "test")
                     p_batch = 1000, # batch size
                     m_batch = 1000, # maximum number of boosting steps per batch
                     b_max = 20000, # maximum number of batches
                     b_stop = 2, # stop if performance not improving for more than 2 batches
                     sl= 0.1, # learning rate (default=0.1)
                     coeff_path_save = TRUE, # should be true, coeff will be saved anyways
                     give_residuals = FALSE, # TRUE if residuals should be saved for each boosting step, attention: large files!
                     family = "gaussian", # family argument: gaussian, binomial, cox or count
                     metric = "MSEP" # loss function/metric: gaussian: MSEP or quantilereg, binomial: log_loss, count: Poisson or negative_binomial, cox: weightedL2, AFT-Weibull, AFT-logistic, AFT-normal, Cox, C (last two slow + memory issues)
                     ) 


time_snpboost_end <- Sys.time()

#predict PRS for subset "test"
pred_test_snpboost <- predict_snpboost(fit_snpboost,genotype.pfile,phenotype.file,phenotype,subset="test")

#predict PRS for all samples
pred_all_snpboost <- predict_snpboost(fit_snpboost,genotype.pfile,phenotype.file,phenotype)

# extract coefficients and number of chosen SNPs
idx <- which.min(fit_snpboost$metric.val)
betas <- get_coefficients(fit_snpboost$beta,idx,covariates=fit_snpboost$configs[['covariates']])
nmb_snps <- length(betas)-1-length(fit_snpboost$configs[['covariates']])

## fit glm with PRS and covariates on train and val data and predict on test data

data <- fread(phenotype.file) %>% mutate(IID = as.character(IID))
pred <- data.table(IID=rownames(pred_all_snpboost$prediction), PRS =pred_all_snpboost$prediction)

data <- full_join(data,pred)%>%rename(PRS=PRS.V1) %>% filter(complete.cases(.))

formula_model <- as.formula(paste0(phenotype,"~ sex + age_at_2015 +",paste("PC",1:10,collapse="+",sep=""),"+ PRS"))

full_model <- glm(formula_model,data=data %>% filter(train_test %in% c("train","val")),family="gaussian")

pred_full_model_test <- predict.glm(full_model, newdata = data %>% filter(train_test %in% c("test")))

MSEP_test_full_model=mean((data[data$train_test %in% c("test"),]%>%pull(all_of(phenotype))-pred_full_model_test)^2)

cor_squared_test_full_model=cor(data[data$train_test %in% c("test"),]%>%pull(all_of(phenotype)),pred_full_model_test)^2

MSEP_only_PRS <- pred_test_snpboost$metric

cor_squared_only_PRS <- cor(data[data$train_test=="test",]%>%pull(all_of(phenotype)),pred_test_snpboost$prediction)^2

# model only covariates

formula_model_cov <- as.formula(paste0(phenotype,"~ sex + age_at_2015 +",paste("PC",1:10,collapse="+",sep="")))

cov_model <- glm(formula_model_cov,data=data %>% filter(train_test %in% c("train","val")),family="gaussian")

pred_cov_model_test <- predict.glm(cov_model, newdata = data %>% filter(train_test=="test"))

MSEP_cov_model=mean((data[data$train_test=="test",]%>%pull(all_of(phenotype))-pred_cov_model_test)^2)

cor_squared_test_cov_model=cor(data[data$train_test %in% c("test"),]%>%pull(all_of(phenotype)),pred_cov_model_test)^2


save.image(file=paste0("/daten/epi/PRS/snpnet/SNP_boost/functions/test.RData"))