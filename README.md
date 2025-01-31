# Polygenic Risk Scoring Protocol for PGC SUI
This repository details a Standard Operating Procedure for the data analysts that will be performing polygenic risk score analyses in individual PGC SUI cohorts. We provide a fully automated analysis pipeline that will analyze and format the results to be returned.

The following analyses should be repeated for each phenotype available in your sample. Analyses should also be conducted separately for each genetic ancestry available in your sample. 

If you encounter any issues or have any questions please reach out to Sarah Colbert (sarah.colbert@icahn.mssm.edu)

## Pipeline Prerequisites and Requirements

#### Software Requirements

* [Plink 1.9](https://www.cog-genomics.org/plink/)
* R >= 3.3
  * [tidyverse](https://github.com/tidyverse/tidyverse)
  * [lme4](https://cran.r-project.org/web/packages/lme4/index.html)
  * [fmsb](https://cran.r-project.org/web/packages/fmsb/index.html)
  * [pROC](https://cran.r-project.org/web/packages/pROC/index.html)
  * [boot](https://cran.r-project.org/web/packages/boot/index.html)

#### Data Requirements

1) SNP weights for each phenotype-ancestry combination available in your cohort.

2) Imputed allele dosages or best guess genotypes (imputed allele dosages converted to hard calls) in plink binary format are both fine. Please only include individuals who were included in the GWAS. 

3) Phenotype data for all available phenotypes (these should exactly match the cases and controls used for the GWAS). Phenotypes should be coded as binary variables such that cases = 1 and controls = 0 (so keep this in mind if you are using a plink .pheno file which uses 2,1,-9 coding!!)

4) Complete covariate data (these should be the exact covariates that were used when running the GWAS). Please do not include age or sex. 


## Step 1: Score individuals using PLINK

To calculate polygenic risk scores in your target dataset using the provided SNP weights file on google drive, you can use the [score function](https://www.cog-genomics.org/plink/1.9/score) from plink. Please note that after downloading the weights files, the code below decompresses them, so you don't need to do so prior. Example code: 

```
## set variables
target=namex ## replace namex your 5 character cohort code (same as what was used in the weights file sent to you)
ancestry=anc ## replace anc the 3 character ancestry code (same as what was used in the weights file sent to you: afr, eas, eur, or lat)
phenotype=phenox ## replace phenox with the phenotype code (either si, sa_B1, sa_B2, or sd)

## decompress weights file for use in plink
gunzip /my/path/weights/${target}_${ancestry}_${phenotype}_META_pst_eff_a1_b0.5_phi1e-02_all.txt.gz

## run scoring in plink (you will need to replace with your paths)
/my/path/plink/plink \
--bfile /my/path/plink_bfiles \ 
--score /my/path/weights/${target}_${ancestry}_${phenotype}_META_pst_eff_a1_b0.5_phi1e-02_all.txt 2 4 6 \
--out /my/path/${target}_${ancestry}_${phenotype}_prscsx
```

## Step 2: Test PRS associations

Below we provide code for how to test the association between a single PRS and phenotype of interest. You will need to make some edits/replacements in the code chunks which start with **"## !!EDIT:"**

The code will first merge together prs, phenotype and covariate data. **These should be the exact phenotypes and covariates used in the GWAS.** Then it will run two logistic regressions, one which includes the PRS as a predictor (PRS+PCs) and one which does not (PCs only). It will use the regression output to calculate Nagelkerke's R^2^ (using the fmsb package) and liability R^2^. To calculate R^2^ on the liability scale, it uses the following population prevalences (K): 

| Phenotype | K             |
| :-------- |--------------:|
| SI        | 0.09          |
| SA        | 0.02          |
| SD        | 0.001         |

```
# Author: Sarah Colbert
# Title: Run PRS associations for PGC SUI cohorts
# Date: 20250131

## -----------------------------------------------------------
## load packages ---------------------------------------------
## -----------------------------------------------------------

## load packages
library(dplyr)
library(lme4)
library(fmsb)
library(pROC)
library(boot)

## -----------------------------------------------------------
## !!EDIT: set up variables ----------------------------------
## -----------------------------------------------------------

## please create the following variables
target_name <- "target" ## this is your 5 character cohort code (same as what was used in the weights file sent to you)
ancestry <- "eur" ## this is the 3 character ancestry code (same as what was used in the weights file sent to you: afr, eas, eur, or lat)
phenotype <- "phenox" ## this is the phenotype code (either si, sa_B1, sa_B2, or sd)
phe_col <- "MYPHENO" ## this is the name of the phenotype column in your phenotype file (incase it is not the same as the phenotype above)
analyst <- "initials" ## this is the analyst's initials
out_dir <- "/my/path/output/" ## this is the output directory for the results file

## -----------------------------------------------------------
## !!EDIT: load target files ---------------------------------
## -----------------------------------------------------------

## load scores (replace scores file name)
scores_df <- read.table("/path/to/target_ancestry_phenotype_prscsx.profile", h = T) %>% ## replace with your score file from plink
 select(FID, IID, SCORE) %>% mutate(FID=as.character(FID), IID=as.character(IID)) ## read ID cols as characters

## load phenotypes (replace phenotype file name)
## (assumes that phenotype file has an FID and IID column like in PLINK)
## please remove the mutate() line if your phenotype is already correctly coded with 0,1,NA
phenos_df <- read.table("target_phenos.txt", h = T) %>%
 mutate(!!phe_col := ifelse(.data[[phe_col]] == 1, 0, ifelse(.data[[phe_col]] == 2, 1, NA))) %>% ## (remove this line if already coded correctly!!) 
 filter(!is.na(.data[[phe_col]])) %>% ## remove individuals missing phenotype data
 mutate(FID=as.character(FID), IID=as.character(IID)) ## read ID cols as characters

## load covariates
## (assumes that covariate file has an FID and IID column like in PLINK)
covs_df <- read.table("target_covs.txt", h = T) %>% mutate(FID=as.character(FID), IID=as.character(IID)) ## read ID cols as characters

## -----------------------------------------------------------
## merge target data -----------------------------------------
## -----------------------------------------------------------

## merge all datasets together (assuming all individuals have a unique IID)
joint_df <- inner_join(phenos_df, scores_df, by = c("FID", "IID")) %>%
  inner_join(covs_df, by = c("FID", "IID"))

## -----------------------------------------------------------
## !!EDIT: set up regression models --------------------------
## -----------------------------------------------------------

## for this code chunk please edit the below line to set the regression covariates to include the proper PCs for your cohort
covariates <- "C1 + C2 + C3 + C4 + C6 + C8 + C14 + C16"

## don't need to edit anything else in this chunk
## make model with PRS + covars
prs_model_text <- paste0(phe_col, " ~ scale(SCORE) + ", covariates)
## make model without PRS
base_model_text <- paste0(phe_col, " ~ ", covariates)

## make sure that going forward for all other calcs we only use individuals that were included in full regression (i.e. have all vars available)
complete_rows <- which(complete.cases(joint_df[, all.vars(as.formula(prs_model_text))]))
full_df <- joint_df[complete_rows,] 

## -----------------------------------------------------------
## R2 calculations -------------------------------------------
## -----------------------------------------------------------

## make function to calculate liability R2
lR2 <- function(df, type) {
  
  ## run model with PRS + covars
  prs_model <- glm(prs_model_text, family = binomial(link = 'logit'), data = df)
  ## run model without PRS (i.e. covars only)
  base_model <- glm(base_model_text, family = binomial(link = 'logit'), data = df) 

  ## calculate R2N
  R2N <- (NagelkerkeR2(prs_model)$R2)-(NagelkerkeR2(base_model)$R2)

  ## set some variables that will be used to calc liability R2
  K <- ifelse(phenotype=="si", 0.09, ifelse(phenotype=="sa", 0.02, ifelse(phenotype=="sd", 0.001, NA))) ## sets population prevalence
  N <- nrow(df)
  N_cases <- length(which(df[[phe_col]] == 1))
  N_controls <- N-N_cases
  N_eff <- (4*N_cases*N_controls)/(N_cases+N_controls)
  P <- N_cases/N

  ## calculate liability R2
  thd = -qnorm(K,0,1) #threshold on normal distribution which truncates the proportion of #disease prevalence
  zv = dnorm(thd) #z (normal density)
  mv = zv/K #mean liability for case
  theta = mv*(P-K)/(1-K)*(mv*(P-K)/(1-K)-thd) #theta in equation (15)
  cv = K*(1-K)/zv^2*K*(1-K)/(P*(1-P)) #C in equation (15) 
  max=1-P^(2*P)*(1-P)^(2*(1-P))
  R2O=R2N*max #Convert NKR2 back to Cox & Snell R2, equivalent to R2 on observed scale #from a linear model
  R2 = R2O*cv/(1+R2O*theta*cv)

  ## if main analysis, then return these results
  if(type=="main"){
    return(list(R2 = R2, N = N, N_cases = N_cases, N_controls = N_controls, N_eff = N_eff, R2N = R2N, model_beta = (c(summary(prs_model)$coefficients[2,1])), model_se = (c(summary(prs_model)$coefficients[2,2])), model_p = (c(summary(prs_model)$coefficients[2,4]))))
  }

  # if doing bootstrapping, then only return R2
  if(type=="boot"){
    return(R2)
  }
}

## calculate main PRS results using the function
main_prs_results <- lR2(full_df, "main")

## -----------------------------------------------------------
## bootstrapping to get R2 CIs -------------------------------
## -----------------------------------------------------------

## create a wrapper function to pass to boot
wrapper_function <- function(data, indices) {
  ## resample the data
  resampled_data <- data[indices, ]
  ## call lR2 with the "boot" type argument
  return(lR2(resampled_data, type = "boot"))
}

## do bootstrapping
n_resamp <- 1000
boot_results <- boot(data = full_df, statistic = wrapper_function, R = n_resamp)
R2_se <- sd(boot_results$t)

## calculate confidence interval
R2_95CI_low <- boot_results$t0-(1.96*R2_se)
R2_95CI_high <- boot_results$t0+(1.96*R2_se)

## -----------------------------------------------------------
## AUC calculation -------------------------------------------
## -----------------------------------------------------------

## logistic model with score only (copying formulas from ricopili danscore_3)
tstS_text <- paste0(phe_col, " ~ scale(SCORE)")
tstS <- glm(tstS_text, family = binomial(link = 'logit'), data = full_df)

## calculate AUC and it's SE
roc_obj <- roc(full_df[[phe_col]],tstS$linear.predictors)
aucvS <- auc(roc_obj)
aucvS_se <- sqrt(var(roc_obj)) 

## -----------------------------------------------------------
## OR calculations -------------------------------------------
## -----------------------------------------------------------

## residualize
resids_text <- paste0("SCORE ~ ", covariates)
resids_model <- lm(resids_text, data = full_df)

## add standardized residuals to dataframe
resids_df <- full_df %>% mutate(prs_resid=rstandard(resids_model))

## get different quantiles
resids_df$quint <- ntile(resids_df$prs_resid, 5)
resids_df$dec <- ntile(resids_df$prs_resid, 10)

## quintile groups
top_q <- resids_df %>% filter(quint==5) %>% mutate(top=1)
mid_q <- resids_df %>% filter(quint==3) %>% mutate(top=0) ## middle quintile
bot_q <- resids_df %>% filter(quint==1) %>% mutate(top=0) ## bottom quintile
## decile groups
top_d <- resids_df %>% filter(dec==10) %>% mutate(top=1)
mid_d <- resids_df %>% filter(dec==5 | dec==6) %>% mutate(quart=ntile(prs_resid, 4)) %>% filter(quart==2 | quart==3) %>% select(-quart) %>% mutate(top=0) ## middle quintile
bot_d <- resids_df %>% filter(dec==1) %>% mutate(top=0) ## bottom quintile

## make different combos of quartiles to compare
combos <- list(
  top_mid_q = bind_rows(top_q, mid_q),
  top_bot_q = bind_rows(top_q, bot_q),
  top_mid_d = bind_rows(top_d, mid_d),
  top_bot_d = bind_rows(top_d, bot_d)
)

## OR functions
or_model_text <- paste0(phe_col, " ~ top")
ors <- function(combo){
    df <- combos[[combo]]
    or_model <- glm(or_model_text, family = binomial(link = 'logit'), data = df) 
    OR <- exp(or_model$coefficients[2])
    ORL <- exp(or_model$coefficients[2]-1.96*summary(or_model)$coefficients[2,2])
    ORH <- exp(or_model$coefficients[2]+1.96*summary(or_model)$coefficients[2,2])
    as.data.frame(cbind(OR, ORL, ORH)) %>% 
        rename(!!paste0("OR_", combo):="OR") %>% 
        rename(!!paste0("ORL_", combo):="ORL") %>% 
        rename(!!paste0("ORH_", combo):="ORH")  %>% `rownames<-`( NULL )
}

## calc ORs and 95% CIs for all combos
or_results <- bind_cols(lapply(names(combos), ors))

## save some of the Ns from the quantiles to report
quant_df_names <- c("top_q", "mid_q", "bot_q", "top_d", "mid_d", "bot_d")
## create list with counts of cases and controls in each quantile group
counts_list <- lapply(quant_df_names, function(df) {
  tab <- table(get(df)[[phe_col]]) 
  data.frame(id = 1, quant_name = df,
    Ncases = ifelse("1" %in% names(tab), tab["1"], 0), Ncontrols = ifelse("0" %in% names(tab), tab["0"], 0)
  )
})
## make list of counts into a df
OR_Ns_df <- reshape(as.data.frame(do.call(rbind, counts_list)), direction = "wide", idvar="id", timevar="quant_name") %>% select(-id)

## -----------------------------------------------------------
## Compile results -------------------------------------------
## -----------------------------------------------------------

## compile results into a table
prs_results <- cbind("cohort" = target_name,  
                    "phenotype" = phenotype,
                    "ancestry" = ancestry, 
                    "beta" = main_prs_results$model_beta, ## this is the effect size of the PRS so long as the PRS is the first predictor in the model 
                    "se" = main_prs_results$model_se, ## this is the standard error of the effect size of the PRS so long as the PRS is the first predictor in the model 
                    "p" = main_prs_results$model_p, ## this is the p value of the association with the PRS so long as the PRS is the first predictor in the model
                    "Nagelkerke_R2" = main_prs_results$R2N, ## this is Nagelkerke's R2
                    "liability_R2" = main_prs_results$R2, ## this is the liability R2
                    "liability_R2_se" = R2_se, ## this is se of liability R2 from bootstrapping
                    "liability_R2_95CI_low" = R2_95CI_low, ## lower range of 95% CI for liability R2
                    "liability_R2_95CI_high" = R2_95CI_high, ## upper range of 95% CI for liability R2
                    "AUC" = aucvS, ## this is what we think is the most appropriate estimate of AUC attributed to the score (even tho covars are ignored)
                    "AUC_se" = aucvS_se, ## this is standard error of AUC
                    "N_cases" = main_prs_results$N_cases, 
                    "N_controls" = main_prs_results$N_controls, 
                    "N"= main_prs_results$N, 
                    "N_eff" = main_prs_results$N_eff, 
                    or_results, ## OR calc results
                    OR_Ns_df) ## OR quantile group counts

## -----------------------------------------------------------
## Generate results file -----------------------------
## -----------------------------------------------------------

## save the results in specified output directory
write.csv(prs_results, paste0(out_dir, "/", target_name, "_", ancestry, "_", phenotype, "_", analyst, "_prs_results.csv"), row.names = F) 

```

## Step 3: Share PRS results

Please email the results/csv file(s) to Sarah Colbert (sarah.colbert@icahn.mssm.edu). 
