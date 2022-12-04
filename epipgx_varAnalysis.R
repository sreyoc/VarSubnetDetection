#This script includes the nonsynonymous and lof functions separately for EpiPGX data

## Load libraries
library(tidyverse)
library(here)
library(kableExtra)
library(org.Hs.eg.db)
library(annotate)
library(clusterProfiler)
library(ssh)
library(here)
# not required : library(mygene)

## Load the tsv files
epipgx_lof <- read_tsv(here::here("pvalues","epipgx_lof.tsv"))
epipgx_nonsyn <- read_tsv(here::here("pvalues","epipgx_nonsyn.tsv"))

##  Select the minimum of the pvalues as heat scores
epipgx_scores  <- epipgx_lof %>% full_join(epipgx_nonsyn, by ="Range") %>% 
  mutate(scores = ifelse(PvalueTwoSide.x <= PvalueTwoSide.y & is.na(PvalueTwoSide.x) == FALSE | 
                           is.na(PvalueTwoSide.y) == TRUE, PvalueTwoSide.x, PvalueTwoSide.y)) %>% 
  dplyr::select(Range, scores)

#### Functions
make_scores <- function(file, gene, col1) {
  select_score <- read_tsv(here::here("pvalues", file)) %>% 
    dplyr::select(gene, col1)
}

writetoFile <- function(var, name){
  write_tsv(var, path
            =paste0("/Users/sreyoshi.chatterjee/Documents/workspace/projects/VarSubnetDetection/",name), 
            col_names = FALSE)
}

#####

# EpiPGX_lof
epipgx_lof_pvalue <- make_scores("epipgx_lof.tsv", "Range", "PvalueTwoSide") %>%
  mutate(hscores = round(-log10(PvalueTwoSide), 1)) %>% dplyr::select(Range, hscores)
epipgx_lof_odds <-  make_scores("epipgx_lof.tsv", "Range", "odds_ratio")

## EpiPGX_nonsyn
epipgx_nonsyn_pvalue <- make_scores("epipgx_nonsyn.tsv", "Range", "PvalueTwoSide") %>% 
  mutate(hscores = round(-log10(PvalueTwoSide), 1)) %>% dplyr::select(Range, hscores)
epipgx_nonsyn_odds <-  make_scores("epipgx_nonsyn.tsv","Range", "odds_ratio")

## String data
string_hs_interaction <- read.table(here("9606.protein.links.v11.0.txt"), header = TRUE)
string_hs_geneNames <- read_delim(here("9606.protein.info.v11.0.txt"), delim = "\t")

# Select string data with filtering combined score
string_hs_interaction_geneNames <- 
  string_hs_interaction %>% 
  dplyr::inner_join(string_hs_geneNames, by = c("protein1" = "protein_external_id")) %>%
  filter(combined_score > 900) %>% 
  dplyr::select(-protein_size, -annotation, -combined_score) %>% 
  dplyr::rename(GeneA= preferred_name) %>% 
  inner_join(string_hs_geneNames, by = c("protein2" = "protein_external_id")) %>%
  dplyr::rename(GeneB = preferred_name) %>% dplyr::select(GeneA, GeneB)

## HumanBase topedges data

hb_topedges_raw <-  read_delim(here::here("humanbase", "brain_top.txt"), delim = "\t", 
                               col_names = FALSE) %>% dplyr::filter(X3 >= 0.5)

# Convert the entrez gene ids to gene symbol
hb_topedges <- hb_topedges_raw %>% 
  mutate(GeneA = getSYMBOL(as.character(hb_topedges_raw$X1), data='org.Hs.eg'),
         GeneB = getSYMBOL(as.character(hb_topedges_raw$X2),data = 'org.Hs.eg')) %>% 
  dplyr::select(GeneA, GeneB) 

# Function to create hhotnet input files 
generate_networkFiles <- function(int_data, epipgx_data, flag, network) {
  
  if (flag == 0) { 
    index_gene_id <- int_data %>% 
      dplyr::select(GeneA) %>% 
      distinct(GeneA, .keep_all = TRUE) %>% 
      rowid_to_column(var = "ID") 
  }
  else {
    BdiffA <- setdiff(int_data$GeneB, 
                      int_data$GeneA) 
    
    index_gene_id <- int_data %>% dplyr::select(GeneA)   %>% 
      add_row(GeneA = BdiffA) %>% distinct(GeneA, .keep_all = TRUE) %>% 
      rowid_to_column(var = "ID")
  }
  
  common_index <- intersect(index_gene_id$GeneA, epipgx_data$Range)
  
  index_gene_sub <- index_gene_id %>% 
    filter(GeneA %in% common_index) %>% dplyr::select(GeneA) %>% 
    distinct(GeneA, .keep_all = T)  %>% rowid_to_column(var = "ID") 
  
  edge_list_sub <- index_gene_sub %>% inner_join(int_data, by = "GeneA") %>% 
    inner_join(index_gene_sub, by = c("GeneB" = "GeneA")) %>% 
    dplyr::select(ID.x, ID.y) %>% 
    drop_na() %>% distinct()
  
  network_sub <- index_gene_sub %>% inner_join(int_data, by = "GeneA") %>% 
    inner_join(index_gene_sub, by = c("GeneB" = "GeneA")) %>% 
    dplyr::select(GeneA, GeneB) %>% 
    drop_na() %>% distinct()
  
  writetoFile(index_gene_sub, paste0(network,"_index_gene.tsv"))
  writetoFile(edge_list_sub, paste0(network,"_edge_list.tsv"))
  writetoFile(network_sub, paste0(network,".tsv"))
}

## Checking if the function works correctly
d <- generate_networkFiles(string_hs_interaction_geneNames, epipgx_scores, 0, "network_1")
d2 <- generate_networkFiles(hb_topedges, epipgx_scores, 1, "network_2")

## EPIPGX LOF
## Calling for a single run of hierarchical hotnet
writetoFile(epipgx_lof_pvalue, "scores_1.tsv")
writetoFile(epipgx_lof_odds, "scores_2.tsv")
generate_networkFiles(string_hs_interaction_geneNames, epipgx_lof_odds, 0, "network_1")
generate_networkFiles(hb_topedges, epipgx_lof_odds, 1, "network_2")

#Upload to iris-cluster via scp
session <- ssh_connect("schatterjee@access-iris.uni.lu:8022")
file <- c("network_1_index_gene.tsv", "network_1_edge_list.tsv", "network_1.tsv", "scores_1.tsv",
          "network_2_index_gene.tsv", "network_2_edge_list.tsv","network_2.tsv", "scores_2.tsv")
scp_upload(session, file, to = "~/mechEPI_variantAnalysis/hhotnet/epipgx/2nw_2sc_lof/data")
ssh_disconnect(session)


## EPIPGX NONSYN
## Calling for a single run of hierarchical hotnet
writetoFile(epipgx_nonsyn_pvalue, "scores_1.tsv")
writetoFile(epipgx_nonsyn_odds, "scores_2.tsv")
## Only for this part, the combined protein score has been kept as 900
## as the number of edges were exceeding 6,00,000 
generate_networkFiles(string_hs_interaction_geneNames, epipgx_nonsyn_odds, 0, "network_1")
generate_networkFiles(hb_topedges, epipgx_nonsyn_odds, 1, "network_2")

#Upload to iris-cluster via scp
session <- ssh_connect("schatterjee@access-iris.uni.lu:8022")
file <- c("network_1_index_gene.tsv", "network_1_edge_list.tsv", "network_1.tsv", "scores_1.tsv",
          "network_2_index_gene.tsv", "network_2_edge_list.tsv","network_2.tsv", "scores_2.tsv")
scp_upload(session, file, to = "~/mechEPI_variantAnalysis/hhotnet/epipgx/2nw_2sc_nonsyn/data")
ssh_disconnect(session)


## Download from cluster
session <- ssh_connect("schatterjee@access-iris.uni.lu:8022")
file <- c( "consensus_edges.tsv", "consensus_nodes.tsv")
scp_download(session, "~/mechEPI_variantAnalysis/hhotnet/epipgx/2nw_2sc_nonsyn/tmp/*", 
             to = "~/Documents/workspace/projects/VarSubnetDetection/results/epipgx/2nw_2sc_nonsyn/")
ssh_disconnect(session)


### Over-representation analysis

#### Important functions
# Calling enrichGO
call_enrichGO <- function(enrich_nodes){
  enrichGO(
    gene = enrich_nodes,
    keyType = "SYMBOL",
    OrgDb   = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    pAdjustMethod = "BH"
  )
}

# creating dotplot
create_dotplot <- function( enrich_out,num, name) {
  dotplot(enrich_out, showCategory = num, font.size = 8, title = name)
}


## @knitr lofplots
# LOF
# total 31 genes have been identified

#Preparing  the consensus nodes
consensus_lof <-
  read_delim(
    here::here("results", "epipgx", "2nw_2sc_lof", "consensus_nodes.tsv"),
    delim = "\t",
    col_names = FALSE) 
consensus_lof_C <-paste(unlist(consensus_lof)) %>% 
  as_tibble() %>% 
  filter(!value == "NA")
colnames(consensus_lof_C) <- NULL
consensus_lof_C <- as.matrix(consensus_lof_C)

lof_enrich_epi <- call_enrichGO(consensus_lof_C)
# 
head(lof_enrich_epi) %>% kable %>% 
  kable_styling(bootstrap_options = "striped", font_size= 10, position = "left")

# Creating the plots
create_dotplot(lof_enrich_epi, 10, "")
cnetplot(lof_enrich_epi, 5)

## @knitr nonsynplots
# NON-SYN
# Total genes identified : 175 (for threshold = 2)
consensus_nonsyn <-
  read_delim(
    here::here("results", "epipgx", "2nw_2sc_nonsyn", "consensus_nodes.tsv"),
    delim = "\t",
    col_names = FALSE) 
consensus_nonsyn_C <-paste(unlist(consensus_nonsyn)) %>% 
  as_tibble() %>% 
  filter(!value == "NA")
colnames(consensus_nonsyn_C) <- NULL
consensus_nonsyn_C <- as.matrix(consensus_nonsyn_C)

nonsyn_enrich_epi <- call_enrichGO(consensus_nonsyn_C)

#t <- gofilter(nonsyn_enrich_epi, level = 4)
head(nonsyn_enrich_epi, 10) %>% kable() %>% 
  kable_styling(bootstrap_options = "striped", font_size= 10, position = "left") 

create_dotplot(nonsyn_enrich_epi, 10, " ")
cnetplot(nonsyn_enrich_epi, showCategory = 7) +
  ggtitle("Gene Network relationship of top 5 enriched terms")
#cnetplot(t, showCategory = 5) +   ggtitle("\n\n Level 4 Gene Network relationship")
