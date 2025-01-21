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

#### Data Requirements

1) SNP weights for each phenotype-ancestry combination available in your cohort.

2) Imputed allele dosages or best guess genotypes (imputed allele dosages converted to hard calls) in plink binary format are both fine. Please only include individuals who were included in the GWAS. 

3) Phenotype data for all available phenotypes (these should exactly match the cases and controls used for the GWAS). Phenotypes should be coded as binary variables such that cases = 1 and controls = 0 (so keep this in mind if you are using a plink .pheno file which uses 2,1,-9 coding!!)

4) Complete covariate data (these should be the exact covariates that were used when running the GWAS). Please do not include age or sex. 


## Step 1: Score individuals using PLINK

To calculate polygenic risk scores in your target dataset using the provided SNP weights file, you can use the [score function](https://www.cog-genomics.org/plink/1.9/score) from plink. Example code: 

```
/my/software/plink/plink \
--bfile /my/target_plink \
--score ${weights_dir}/si_weights_META_all 2 4 6 \
--out /my/scores/target_si_prs
```

## Step 2: Test PRS associations

Below we provide an example for how to test the association between a single PRS and phenotype of interest. Some cohorts may have several phenotypes or genetic ancestries to test and in those cases it may be easier to loop through those groups, so we also provide an example script in this repostiory called multiple_PRS_associations_example.R for those cohorts to refer to. 

The code below will first merge together prs, phenotype and covariate data. Then it will run two logistic regressions, one which includes the PRS as a predicor and one which does not. It is very important that the PRS is scaled (i.e., mean 0, sd 1) so that the effect sizes are standardized and comparable across target datasets. 

The code also shows how to calculate Nagelkerke's R^2^ and liability R^2^. To calculate R^2^ on the liability scale, you should use the following population prevalences (K): 

| Phenotype | K             |
| :-------- |--------------:|
| SI        | 0.09          |
| SA        | 0.02          |
| SD        | 0.001         |

```
## -----------------------------------------------------------
## set up ----------------------------------------------------
## -----------------------------------------------------------

## load packages
library(dplyr)
library(lme4)
library(fmsb)

## load scores
scores_df <- read.table("/my/scores/target_PHENO_prs.profile")

## load phenotypes
phenos_df <- read.table("target_phenos.txt") %>%
 mutate(PHENO=ifelse(PHENO=1, 0, ifelse(PHENO=2, 1, NA)) ## recode the phenotype properly

## load covariates
covs_df <- read.table("target_covs.txt")

## merge all datasets together (assuming all individuals have a unique IID)
full_df <- left_join(scores_df, phenos_df, by = "IID") %>%
  left_join(covs_df, by = "IID")

## -----------------------------------------------------------
## test associations -----------------------------------------
## -----------------------------------------------------------

## run model with PRS
prs_model <- glm("PHENO ~ scale(SCORE) + PC1 + PC2 + PC3 + PC4 + PC5", family = binomial(link = 'logit'), data = full_df)

## run model without PRS
base_model <- glm("PHENO ~ PC1 + PC2 + PC3 + PC4 + PC5", family = binomial(link = 'logit'), data = full_df)

## -----------------------------------------------------------
## R2 calculations -------------------------------------------
## -----------------------------------------------------------

## calculate R2N
R2N <- (NagelkerkeR2(prs_model)$R2)-(NagelkerkeR2(base_model)$R2)

## set some variables that will be used to calc liability R2
K <- 0.09
N <- (c(nobs(prs_model))) ## set the N (grabbing this from the model)
N_cases <- length(which(full_df[[PHENO]]==1))
P <- N_cases/N

## calculate liability R2
## do not need to change anything in this code chunk
thd = -qnorm(K,0,1) #threshold on normal distribution which truncates the proportion of #disease prevalence
zv = dnorm(thd) #z (normal density)
mv = zv/K #mean liability for case
theta = mv*(P-K)/(1-K)*(mv*(P-K)/(1-K)-thd) #theta in equation (15)
cv = K*(1-K)/zv^2*K*(1-K)/(P*(1-P)) #C in equation (15) 
max=1-P^(2*P)*(1-P)^(2*(1-P))
R2O=R2N*max #Convert NKR2 back to Cox & Snell R2, equivalent to R2 on observed scale #from a linear model
R2 = R2O*cv/(1+R2O*theta*cv)

## compile results into a table
## please make sure to change the lines that are indicated
prs_results <- cbind("cohort" = "target", ## please use your 5 character cohort code 
                    "phenotype" = "PHENO", ## please indicate one: "si", "sa", "sd"
                    "ancestry" = "eur", ## please indicate one: "afr", "csa", "eas", "eur", "lat" 
                    "estimate" = (c(summary(prs_model)$coefficients[2,1])), ## this is the effect size of the PRS so long as the PRS is the first predictor in the model 
                    "se" = (c(summary(prs_model)$coefficients[2,2])), ## this is the standard error of the effect size of the PRS so long as the PRS is the first predictor in the model 
                    "p" = (c(summary(prs_model)$coefficients[2,4])), ## this is the p value of the association with the PRS so long as the PRS is the first predictor in the model
                    "Nagelkerke_R2" = R2N, ## this is Nagelkerke's R2
                    "liability_R2" = R2, ## this is the liability R2
                    "N_cases" = N_cases, 
                    "N"= N))

## save results
write.csv(prs_results, "/my/target_PHENO_eur_prs_results.csv", row.names = F)

```

## Step 3: Share PRS results

Please email the results/csv file(s) to the PGC SUI DRC rep, Sarah Colbert (sarah.colbert@icahn.mssm.edu). 
