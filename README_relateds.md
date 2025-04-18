# Polygenic Risk Scoring Protocol for PGC SUI (in cohort with relatives)
This repository details a Standard Operating Procedure for the data analysts that will be performing polygenic risk score analyses in individual PGC SUI cohorts **that include related individuals**. We provide a fully automated analysis pipeline that will analyze and format the results to be returned.

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
  * [boot](https://cran.r-project.org/web/packages/boot/index.html) (this may already be installed in base R)
  * [performance](https://cran.r-project.org/web/packages/performance/index.html) 

#### Data Requirements

1) SNP weights for each phenotype-ancestry combination available in your cohort.

2) Imputed allele dosages or best guess genotypes (imputed allele dosages converted to hard calls) in plink binary format are both fine. Please only include individuals who were included in the GWAS. Furthermore, please filter your sample to include at max two first-degree relatives per family, while prioritizing cases. The process for doing so is described in [Pre-Step: Relatives Identification and Filtering](#pre-step-relatives-identification-and-filtering).

3) Phenotype data for all available phenotypes (these should exactly match the cases and controls used for the GWAS). As phenotypes are likely to be coded in PLINK .pheno file format (2,1,-9 coding), we provide code in the pre-step which will recode the phenotype such that cases = 1, controls = 0 and those with missing phenotypes are NA. 

4) Complete covariate data (these should be the exact covariates that were used when running the GWAS). Please do not include age or sex. 

## Pre-Step: Relatives Identification and Filtering

The PRS tests should only be performed in a subset of your sample that has no more than two first-degree relatives per family, with a preference for retaining cases. Relatives more distant than first-degree will be removed, unless they are the only case in a family, then they are retained and all others in the family are removed. 

To determine relatedness, please use [KING's kinship inference](https://www.kingrelatedness.com/manual.shtml#WITHIN) to get kinship values amongst individuals in your cohort who are third-degree relatives or closer. The KING output file (suffix .kin0) should then be used to ID relatives and perform filtering with the script [id_rels_and_filter.R](https://github.com/sarahcolbert/PGC_SUI_PRS/blob/main/id_rels_and_filter.R). **Every individual must have a unique IID in every file for all of these script to work.** 

This script will update your phenotype file to the correct FIDs and remove the appropriate related individuals. **The FIDs and set of individuals in your PLINK and covariate files may therefore not match, but this is okay and will be dealt with in step 2.** 

## Step 1: Score individuals using PLINK

To calculate polygenic risk scores in your target dataset using the provided SNP weights file on google drive, you can use the [score function](https://www.cog-genomics.org/plink/1.9/score) from plink. Please note that after downloading the weights files, the code below decompresses them, so you don't need to do so prior. Example code: 

```
## set variables
target=namex ## replace namex with your 5 character cohort code (same as in the weights file sent to you)
ancestry=anc ## replace anc with the 3 character ancestry code (same as in the weights file sent to you: afr, eas, eur, or lat)
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

The script [pgcsui_prs_wrels_step2.R](https://github.com/sarahcolbert/PGC_SUI_PRS/blob/main/pgcsui_prs_wrels_step2.R) provides code for how to test the association between a single PRS and phenotype of interest. You will need to make some edits/replacements in the code chunks which start with **"## !!EDIT:"**

The code will first merge together prs, phenotype and covariate data. **These should be the exact phenotypes and covariates used in the GWAS.** 

Then it will run two mixed-effects models, one which includes the PRS as a predictor (PRS+PCs+FID) and one which does not (PCs+FID only). It will use the regression output to calculate Nakagawa's R<sup>2</sup> (using the performance package) and liability R<sup>2</sup>. The code will also calculate odds ratios and the AUC. It will then compile all the results into a single csv file that will be shared. 

The script uses a bootstrapping approach that performs 10k resamples to calculate the SE and CIs for the liability R<sup>2</sup>. This step will take ~10 minutes so you may want to run the Rscript inside of a job or interactive session.


## Step 3: Share PRS results

Please email the results/csv file(s) to Sarah Colbert (sarah.colbert@icahn.mssm.edu). 
