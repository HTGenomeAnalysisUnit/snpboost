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

plot_performance <- FALSE


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

#predict PRS for all samples
pred_all_snpboost  <- predict_snpboost(fit_snpboost,genotype.pfile,phenotype.file,phenotype)

# extract coefficients and number of chosen SNPs
idx <- which.min(fit_snpboost$metric.val)
betas <- get_coefficients(fit_snpboost$beta,idx,covariates=fit_snpboost$configs[['covariates']])

intercept <- betas[1]

#option 1: save intercept in a seperate file
fwrite(list(intercept=intercept), file=paste0(fit_snpboost$configs[['results.dir']],"/intercept.txt"))

#option 2: save complete beta vector including intercept
fwrite(list(rsID=str_split(names(betas),"_",simplify=T)[,1],A1=str_split(names(betas),"_",simplify=T)[,2],beta=betas), file = paste0(fit_snpboost$configs[['results.dir']],"/betas.txt"),sep="\t",row.names = F)

#additional: save predictions
fwrite(list(IID = rownames(pred_all_snpboost$prediction), pred = pred_all_snpboost$prediction), file = paste0(fit_snpboost$configs[['results.dir']],"/pred.txt"),sep="\t",row.names = F)

if(plot_performance){
  sparsity <- sapply(1:nrow(fit_snpboost$beta),function(x)length(unique(fit_snpboost$beta$name[1:x])))
  
  df_plot <- data.frame(iteration=1:length(sparsity), sparsity = sparsity-1,
                        metric.train = fit_snpboost$metric.train, metric.val = fit_snpboost$metric.val)
  
  #r2 defined as 1-MSEP/MSEP(null model)
  df_plot$r2.train <- 1-df_plot$metric.train/df_plot$metric.train[1]
  
  df_plot$r2.val <- 1-df_plot$metric.val/df_plot$metric.val[1]
  
  fwrite(df_plot, file=file=paste0(fit_snpboost$configs[['results.dir']],"/df_plot.txt"),sep="\t",row.names = F)
  
  ggplot(df_plot, aes(x = iteration)) +
    geom_line(aes(y = sparsity, color = "sparsity")) +
    geom_line(aes(y = r2.val* max(sparsity)/max(r2.val), color = "r2")) +
    scale_color_manual(name = "", values = c("sparsity" = "blue", "r2" = "red"), labels= c(expression(r^2), "sparsity")) +
    scale_y_continuous(name = "sparsity", sec.axis = sec_axis(~ . * max(df_plot$r2.val)/max(df_plot$sparsity), name = expression(r^2)))+
    theme_minimal()+
    theme(text=element_text(size=18))
  
  ggsave(paste0(fit_snpboost$configs[['results.dir']],"/performance_plot.png"))
  
}