# Polygenic Risk Scoring Protocol for PGC SUI
This repository details a Standard Operating Procedure for the data analysts that will be performing polygenic risk score analyses in individual PGC SUI cohorts. We provide a fully automated analysis pipeline that will analyze and format the results to be returned.

The following analyses should be repeated for each phenotype available in your sample. Analyses should also be conducted separately for each genetic ancestry available in your sample. If you have related individuals in your sample, there are two options. The first would be to remove second degree relatives (or whatever degree of relatedness is optimal for your sample), while preferentially retaining cases, which can be done using the [prune_relateds.R script](https://github.com/sarahcolbert/PGC_SUI_PRS/blob/main/prune_relateds.R). The second option is to retain some related individuals and account for FID in a mixed effects model. To implement option 2, please refer to the [README_relateds.md](https://github.com/sarahcolbert/PGC_SUI_PRS/blob/main/README_relateds.md) instead.

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

#### Data Requirements

1) SNP weights for each phenotype-ancestry combination available in your cohort.

2) Imputed allele dosages or best guess genotypes (imputed allele dosages converted to hard calls) in plink binary format are both fine. Please only include individuals who were included in the GWAS. 

3) Phenotype data for all available phenotypes (these should exactly match the cases and controls used for the GWAS). As phenotypes are likely to be coded in PLINK .pheno file format (2,1,-9 coding), we provide code in step 2 which will recode the phenotype such that cases = 1, controls = 0 and those with missing phenotypes are NA. **IF THE PHENOTYPE IS ALREADY CODED CORRECTLY (1,0,NA) THEN YOU MUST COMMENT OUT LINE 40 IN THE R SCRIPT IN STEP 2.** For SA PRS testing, use SA B2 phenotypes if available; otherwise, use SA B1 phenotypes.
 
4) Complete covariate data (these should be the exact covariates that were used when running the GWAS). Please do not include age or sex. 


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

The script [pgcsui_prs_step2.R](https://github.com/sarahcolbert/PGC_SUI_PRS/blob/main/pgcsui_prs_step2.R) provides code for how to test the association between a single PRS and phenotype of interest. You will need to make some edits/replacements in the code chunks which start with **"## !!EDIT:"**

The code will first merge together prs, phenotype and covariate data. **These should be the exact phenotypes and covariates used in the GWAS.** 

Then it will run two logistic regressions, one which includes the PRS as a predictor (PRS+PCs) and one which does not (PCs only). It will use the regression output to calculate Nagelkerke's R<sup>2</sup> (using the fmsb package) and liability R<sup>2</sup>. The code will also calculate odds ratios and the AUC. It will then compile all the results into a single csv file that will be shared. 

The script uses a bootstrapping approach that performs 10k resamples to calculate the SE and CIs for the liability R<sup>2</sup>. This step will take ~10 minutes so you may want to run the Rscript inside of a job or interactive session.


## Step 3: Share PRS results

Please email the results/csv file(s) to Sarah Colbert (sarah.colbert@icahn.mssm.edu). 
