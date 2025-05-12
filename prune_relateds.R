# Author: Sarah Colbert
# Title: Optimal retention of most cases/individuals when pruning for relatedness
# Date: 20250512

## -----------------------------------------------------------
## load packages ---------------------------------------------
## -----------------------------------------------------------

library(dplyr)
library(data.table)
library(igraph)


## -----------------------------------------------------------
## !!EDIT: load your files -----------------------------------
## -----------------------------------------------------------

## load KING output (replace with your KING file)
king_raw <- read.table("./path/to/my_cohort_relatedness.kin0", header = TRUE) %>% 
  select(ID1, ID2, Kinship) %>% ## select columns
  filter(Kinship > 0.0884) %>% ## retains only pairs >= 2nd degree, can swap this out to be more conservative though
  mutate(ID1=as.character(ID1), ID2=as.character(ID2)) ## make sure IDs are read as characters

## replace MYPHENO with the name of the phenotype column in your phenotype file 
phe_col <- "MYPHENO" 

## load phenotypes (replace phenotype file name) 
## assumes pheno file has only 3 cols: FID, IID, PHENO
## also assumes that all IIDs are unique
## please remove the mutate() line if your phenotype is already correctly coded with 0,1,NA
meta <- read.table("./path/to/original_phenos.txt", h = T) %>% 
  rename("PHENO"=all_of(phe_col)) %>% 
  #mutate(PHENO = ifelse(PHENO == 1, 0, ifelse(PHENO == 2, 1, NA))) %>% ## (remove this line if already coded correctly!!) 
  filter(!is.na(PHENO)) %>% 
  select(FID, IID, PHENO) %>% 
  mutate(FID=IID) ## reassign FID to each unique IID

## set output file name for new phenotype file
outname <- "./path/to/target_phenos.txt"

# ----------------------------------------------------------------------
# Make FIDs ------------------------------------------------------------
# ----------------------------------------------------------------------

## use igraph to find groups where ALL members appear to be relatives
## an individual (A) with two relatives (B and C) that aren't relatives will get both FIDs
## but individuals B and C will each only have their own FID, they won't be in the same family as each other
g <- graph.data.frame(king_raw[,c("ID1","ID2")],directed = F)
cliques <- max_cliques(g, min=2)

## make a data table where FIDs are recoded for each ID
df_family <- data.table(id=names(unlist(cliques)),fid=rep(paste0("f", seq_along(cliques)), sapply(cliques, length)))


# ----------------------------------------------------------------------
# Prune the related individuals ----------------------------------------
# ----------------------------------------------------------------------

## remove individuals that are in more than one family (aka "relateds to the most" that show up in df more than once)
df_family2 <- df_family[!(duplicated(df_family$id) | duplicated(df_family$id, fromLast = TRUE)), ]
rels_to_the_most_remove <- df_family %>% filter(!(id %in% df_family2$id)) %>% distinct(id) %>% pull(id)

## make list of cases from the phenotype file
case_ids <- meta %>% filter(PHENO==1) %>% pull(IID)

## prune rels to max out cases 
## need to exclude related individuals, but we want to make sure the ones we exclude don't have the icds code we want, so we can max our # of cases
## if one or more id in family is a case, then randomly keep one case
## if no cases in the family, then make random list of one ID to keep

## create a new dataframe to store the selected IDs
df_keep <- data.frame(id=integer(), fid=integer())

## get list of family ids
all_fids <- unique(df_family2$fid)
# loop through each unique family ID
for(i in 1:length(all_fids)) {
  # subset the dataframe to include only members of the current family
  fid_i <- all_fids[i]
  df_family_i <- df_family2 %>% filter(fid==fid_i)
  # check if any of the family members are in the "case_ids" vector
  if(any(df_family_i$id %in% case_ids)) {
    # if so, randomly select one of the cases and add it to the new dataframe
    df_family_i_case <- df_family_i %>% filter(id %in% case_ids)
    set.seed(1)
    df_keep <- rbind(df_keep, df_family_i_case[sample(nrow(df_family_i_case), 1),])
  } else {
    # otherwise, randomly select one member of the family to keep
    set.seed(1)
    df_keep <- rbind(df_keep, df_family_i[sample(nrow(df_family_i), 1),])
  }
}


# ----------------------------------------------------------------------
# Combine singletons and pruned relateds -------------------------------
# ----------------------------------------------------------------------

## get pheno df of singletons (e.g. those not in king file)
singletons <- meta %>% filter(!(IID %in% king_raw$ID1 | IID %in% king_raw$ID2))

## get pheno df of rels we are keeping 
rels_keep <- meta %>% filter(IID %in% df_keep$id)

## merge 
final_set <- rbind(singletons, rels_keep)

## rename phenotype column with original name of the phenotype
colnames(final_set)[3] <- phe_col

## save
write.table(final_set, outname, row.names = F, quote = F)
