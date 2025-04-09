# Author: Sarah Colbert
# Title: Making a family dataset in BioMe for testing PRS analyses
# Last Update: 20250409

## -----------------------------------------------------------
## load packages ---------------------------------------------
## -----------------------------------------------------------

## load packages
library(igraph)
library(dplyr)


## -----------------------------------------------------------
## !!EDIT: load your files -----------------------------------
## -----------------------------------------------------------

## load KING output (replace with your KING file)
king_raw <- read.table("./path/to/my_cohort_relatedness.kin0", header = TRUE) %>% select(ID1, ID2, Kinship)

## replace MYPHENO with the name of the phenotype column in your phenotype file 
phe_col <- "MYPHENO" 

## load phenotypes (replace phenotype file name) 
## assumes pheno file has only 3 cols: FID, IID, PHENO
## please remove the mutate() line if your phenotype is already correctly coded with 0,1,NA
meta <- read.table("./path/to/original_phenos.txt", h = T) %>% 
  rename("PHENO"=all_of(phe_col)) %>% 
  mutate(PHENO = ifelse(PHENO == 1, 0, ifelse(PHENO == 2, 1, NA))) %>% ## (remove this line if already coded correctly!!) 
  filter(!is.na(PHENO)) %>% 
  select(FID, IID, PHENO)

## set output file name for new phenotype file
outname <- "./path/to/target_phenos.txt"


# ----------------------------------------------------------------------
# Cleaning -------------------------------------------------------------
# ----------------------------------------------------------------------

## add relationship type according to kinship 
## from KING: dup/MZtwin >0.354, 1st deg = 0.177-0.354, 2nd deg = 0.0884-0.177, 3rd deg  = 0.0442-0.0884
king <- king_raw %>%
  mutate(Relationship = case_when(
    Kinship >= 0.177 ~ "1st",
    TRUE ~ "2nd+"
  )) 

## force to characters
king$ID1 <- as.character(king$ID1)
king$ID2 <- as.character(king$ID2)


## filter to just related individuals
## get unique individuals in related pairs
related_ids <- unique(c(king$ID1, king$ID2))
## filter meta df to just these relateds
meta_related <- meta %>% filter(IID %in% related_ids) %>% rename(name = IID)
meta_related$name <- as.character(meta_related$name)

# ----------------------------------------------------------------------
# !!STOP: Check for TRUE statements ------------------------------------
# ----------------------------------------------------------------------

## confirm all IDs match before building graph
all(king$ID1 %in% meta_related$name)
all(king$ID2 %in% meta_related$name)
## OK to continue if both TRUE -----------------------------------------

## build graph of families
g <- graph_from_data_frame(king, directed = FALSE, vertices = meta_related)

## set node weights based on case/control status
V(g)$weight <- ifelse(V(g)$PHENO == 1, 10, 1)  ## this will prefer cases
## add components
components <- components(g)


# ----------------------------------------------------------------------
# Make decisions on which to keep in families --------------------------
# ----------------------------------------------------------------------

## intialize empty vector to hold IDs to keep 
keep_ids <- c()

for (comp_id in unique(components$membership)) {
  subgraph <- induced_subgraph(g, which(components$membership == comp_id))
  
  ## remove individuals to eliminate 2nd+ relationships
  repeat {
    ## check for 2nd+ edges
    second_deg_edges <- E(subgraph)[E(subgraph)$Relationship == "2nd+"]
    if (length(second_deg_edges) == 0) break
    
    ## find node with lowest weight (this prefers removing controls)
    edge_ends <- ends(subgraph, second_deg_edges)

    ## compute combined weights for each edge pair
    weights <- apply(edge_ends, 1, function(pair) {
        V(subgraph)[pair[1]]$weight + V(subgraph)[pair[2]]$weight
    })
    worst_edge_idx <- which.min(weights)
    worst_pair <- edge_ends[worst_edge_idx, ]
    
    ## remove one node (prefer controls but if both then randomly sample)
    to_remove <- ifelse(V(subgraph)[worst_pair[1]]$weight < V(subgraph)[worst_pair[2]]$weight,
                        worst_pair[1], worst_pair[2])
    subgraph <- delete_vertices(subgraph, to_remove)
  }

  ## now prune if more than 2 individuals remain (limit to ≤2 first-degree relatives)
  if (vcount(subgraph) > 2) {
    ## cluster by 1st deg edges
    first_deg <- E(subgraph)[E(subgraph)$Relationship == "1st"]
    g_first <- subgraph.edges(subgraph, first_deg)
    components_first <- components(g_first)
    
    for (sub_comp in unique(components_first$membership)) {
      subsubgraph <- induced_subgraph(g_first, which(components_first$membership == sub_comp))
      if (vcount(subsubgraph) > 2) {
        ## keep two 1st deg rels (preference for cases according to weight)
        top2 <- names(sort(V(subsubgraph)$weight, decreasing = TRUE))[1:2]
        subsubgraph <- induced_subgraph(subsubgraph, top2)
      }
      keep_ids <- c(keep_ids, V(subsubgraph)$name)
    }
  } else {
    keep_ids <- c(keep_ids, V(subgraph)$name)
  }
}


# ----------------------------------------------------------------------
# Merge IDs we are keeping from families with unrelateds ---------------
# ----------------------------------------------------------------------

## get unrelateds
singletons <- meta %>% filter(!(IID %in% related_ids)) %>% pull(IID)

## now make final analytic dataset
meta_final <- meta %>% filter(IID %in% c(singletons, keep_ids))


# ----------------------------------------------------------------------
# !!STOP: SANITY CHECKS! (all must be TRUE) ----------------------------
# ----------------------------------------------------------------------

## check that individuals actually were removed
nrow(meta_final) < nrow(meta) ## (check for TRUE)

## check for reasonable case/control distribution in new set vs old
table(meta_final$PHENO)
table(meta$PHENO)

## check that no 2nd+ relationship pairs remain
king_final_check <- king %>% filter(ID1 %in% meta_final$IID & ID2 %in% meta_final$IID)
all(king_final_check$Relationship=="1st") ## this should only be 1st (check for TRUE)

## check that no more than two 1st rels per fam
king_1st_final <- king_final_check %>% filter(Relationship == "1st")
g_check <- graph_from_data_frame(king_1st_final, directed = FALSE)
all(components(g_check)$csize<=2) ## should all be <= 2 (check for TRUE)


# ----------------------------------------------------------------------
# Assign family IDs ----------------------------------------------------
# ----------------------------------------------------------------------

## get components for all first degree pairs that were retained 
fid_components <- components(g_check)

## make vector with component assignments
fid_map <- data.frame(IID = names(fid_components$membership),FID = as.character(fid_components$membership))

## make the FIDs by collapsing all IDs in the family
fid_map <- fid_map %>%
  group_by(FID) %>%
  mutate(FID = paste0(sort(IID), collapse = "-")) %>%
  ungroup()
## for singletons they just get ID as FID 
singleton_ids <- meta_final$IID[!(meta_final$IID %in% fid_map$IID)]
singleton_fids <- data.frame(IID = as.character(singleton_ids), FID = as.character(singleton_ids))
## all FIDs
fid_all <- bind_rows(fid_map, singleton_fids) 

## add on to dataframe
meta_final_with_fid <- meta_final %>% select(-FID) %>% mutate(IID=as.character(IID)) %>% left_join(fid_all, by = "IID") %>% select(FID, IID, PHENO)
## rename phenotype column with original name of the phenotype
colnames(meta_final_with_fid)[3] <- phe_col


# ----------------------------------------------------------------------
# Save updated phenotype file ------------------------------------------
# ----------------------------------------------------------------------

## save pheno file with filters and new FIDs (replace path and filename)
write.table(meta_final_with_fid, outname, row.names = F, quote = F)
