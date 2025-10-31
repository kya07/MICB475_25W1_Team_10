#!/usr/bin/env RScript 

#### Load packages ####
# Load in the 'phyloseq', 'tidyverse', and 'ape' package
library(phyloseq)
library(tidyverse)
library(ape) 

#### Load data ####
# Load in the cn_metadata, OTU table, taxonomy file, and phylogenetic tree
cn_meta <- read_delim(file = "cn_metadata_new.tsv", delim = "\t")
cn_otu <- read_delim(file="cn-feature-table.txt", delim = "\t", skip=1)
cn_taxonomy <- read_delim(file = "cn-taxonomy.tsv", delim="\t")
cn_tree <- read.tree("cn-tree.nwk")

#### Format data ####
# Format cn_otu table and read into a phyloseq object
cn_otu_mat <- as.matrix(cn_otu[,-1])
rownames(cn_otu_mat) <- cn_otu$`#OTU ID`
cn_OTU <- otu_table(cn_otu_mat, taxa_are_rows = TRUE)
class(cn_OTU) 

# Format cn metadata and read into a phyloseq object
cn_samp_df <- as.data.frame(cn_meta[,-1])
rownames(cn_samp_df)<- cn_meta$'sample-id'
cn_SAMP <- sample_data(cn_samp_df)
class(cn_SAMP)
    
# Format taxonomy and read into a phyloseq object
cn_tax_mat <- cn_taxonomy |> 
  select(-Confidence) |> 
  separate(col=Taxon, sep="; "   
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) |> 
  as.matrix() 
cn_tax_mat <- cn_tax_mat[,-1]
rownames(cn_tax_mat) <- cn_taxonomy$`Feature ID`
cn_TAX <- tax_table(cn_tax_mat)
class(cn_TAX)

# No adjustments for cn-tree is needed. 
#### Create phyloseq object ####
cn_phyloseq_obj <- phyloseq(cn_OTU, cn_SAMP, cn_TAX, cn_tree)

# Verify by viewing components of phyloseq object
otu_table(cn_phyloseq_obj)  
sample_data(cn_phyloseq_obj)  
tax_table(cn_phyloseq_obj)
phy_tree(cn_phyloseq_obj)