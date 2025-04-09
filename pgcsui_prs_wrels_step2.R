# Author: Sarah Colbert
# Title: Run PRS associations for PGC SUI cohorts
# Date: 20250409

## -----------------------------------------------------------
## load packages ---------------------------------------------
## -----------------------------------------------------------

## load packages
library(dplyr)
library(lme4)
library(fmsb)
library(pROC)
library(boot)
library(performance)

## -----------------------------------------------------------
## !!EDIT: set up variables ----------------------------------
## -----------------------------------------------------------

## please create the following variables
target_name <- "target" ## this is your 5 character cohort code (same as what was used in the weights file sent to you)
ancestry <- "eur" ## this is the 3 character ancestry code (same as what was used in the weights file sent to you: afr, eas, eur, or lat)
phenotype <- "phenox" ## this is the phenotype code (either si, sa, or sd)
phe_col <- "MYPHENO" ## this is the name of the phenotype column in your phenotype file (incase it is not the same as the phenotype above)
analyst <- "initials" ## this is the analyst's initials
out_dir <- "/my/path/output/" ## this is the output directory for the results file

## -----------------------------------------------------------
## !!EDIT: load target files ---------------------------------
## -----------------------------------------------------------

## load scores (replace scores file name)
scores_df <- read.table("./path/to/target_ancestry_phenotype_prscsx.profile", h = T) %>% ## replace with your score file from plink
 select(IID, SCORE) %>% mutate(IID=as.character(IID)) ## read ID col as characters

## load phenotypes (replace phenotype file name with the new phenotype file generated in the pre-step)
phenos_df <- read.table("./path/to/target_phenos.txt", h = T) %>%
 filter(!is.na(.data[[phe_col]])) %>% ## remove individuals missing phenotype data
 mutate(FID=as.character(FID), IID=as.character(IID)) ## read ID cols as characters

## load covariates
## (IMPORTANT: if there is an FID col, remove it since it may not match the corrected ones in the pheno file!)
covs_df <- read.table("./path/to/target_covs.txt", h = T) %>% 
    select(-FID) %>% ## can remove this line if no FID col in your covariate file
    mutate(IID=as.character(IID)) ## read ID cols as characters

## -----------------------------------------------------------
## merge target data -----------------------------------------
## -----------------------------------------------------------

## merge all datasets together (assuming all individuals have a unique IID)
joint_df <- inner_join(phenos_df, scores_df, by = "IID") %>%
  inner_join(covs_df, by = "IID")

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
  
  ## set optimizer for main model vs. boot
  if(type == "boot"){
    ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e4))
  } else {
    ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  }

  ## run model with PRS + covars
  prs_model <- glmer(prs_model_text, family = binomial(link = 'probit'), data = df, control = ctrl)
  ## run model without PRS (i.e. covars only)
  base_model <- glmer(base_model_text, family = binomial(link = 'probit'), data = df, control = ctrl) 

  ## calculate p-value
  prs_pval <- anova(prs_model, base_model, test="Chi")[2,8] 
  ## calculate R2N (Nakagawa's)
  R2N <- as.numeric(r2_nakagawa(prs_model)$R2_marginal-r2_nakagawa(base_model)$R2_marginal[[1]])

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
    return(list(R2 = R2, N = N, N_cases = N_cases, N_controls = N_controls, N_eff = N_eff, R2N = R2N, model_beta = (c(summary(prs_model)$coefficients[2,1])), model_se = (c(summary(prs_model)$coefficients[2,2])), model_p = prs_pval))
  }

  # if doing bootstrapping, then only return R2
  if(type=="boot"){
    return(R2)
  }
}

## calculate main PRS results using the function
main_prs_results <- lR2(full_df, "main")

