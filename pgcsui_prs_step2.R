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

  ## calculate p-value
  prs_pval <- anova(prs_model, base_model, test="Chi")[2,5] 
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
    return(list(R2 = R2, N = N, N_cases = N_cases, N_controls = N_controls, N_eff = N_eff, R2N = R2N, model_beta = (c(summary(prs_model)$coefficients[2,1])), model_se = (c(summary(prs_model)$coefficients[2,2])), model_p = prs_pval))
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
n_resamp <- 10000
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
write.csv(prs_results, paste0(out_dir, "/", target_name, "_", ancestry, "_", phenotype, "_", analyst, format(Sys.Date(),"_%Y%m%d"), "_prs_results.csv"), row.names = F) 
