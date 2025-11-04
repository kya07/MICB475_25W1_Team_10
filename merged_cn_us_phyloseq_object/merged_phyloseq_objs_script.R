#!/usr/bin/env RScript 
# Script for Merging US and CN phyloseq objects together (without phylo trees as 
                                      # one tree has a different number of tips) 

# Note: The U.S. dataset targets 16S V4, while the Chinese dataset targets 16S V4–V5.
# QIIME 2 builds trees from sequence alignments (rep-seqs), not from known taxonomy.
# Because the regions differ, the algorithm would treat the datasets as having distinct phylogenies
# even if similar taxa are present. 
# We will therefore omit the tree from the merged phyloseq object
# and not generate tree-based diversity metrics,  
# as it’s not possible to run the algorithm based on known taxonomies instead of read sequences. 
# This limitation will be noted in the manuscript.

#### Load packages ####
# Load in the 'phyloseq', 'tidyverse', and 'ape' package
library(phyloseq)
library(tidyverse)
library(ape) 

#### Load data #####
cn_meta <- read_delim(file = "cn_metadata_new.tsv", delim = "\t")
cn_otu <- read_delim(file="cn-feature-table.txt", delim = "\t", skip=1)
cn_taxonomy <- read_delim(file = "cn-taxonomy.tsv", delim="\t")
cn_tree <- read.tree("cn-tree.nwk")

meta <- read_delim(file= "us_metadata_new.tsv", delim = "\t")
otu <- read_delim(file= "feature-table.txt", delim = "\t", skip=1)
tax <- read_delim("taxonomy.tsv", delim="\t")
phylotree <- read.tree("tree.nwk")

#### U.S. Phyloseq Object #####
##### Format sample metadata #####
# Save everything except sample-id as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sample-id the rownames
rownames(samp_df) <- meta$ "sample-id"
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

##### Format OTU table #####
# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
class(OTU)

##### Formatting taxonomy #####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax |>  select(-Confidence)|>
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) |>
  as.matrix() # Saving as a matrix
# Save everything except feature IDs
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX) 

##### Create us_phyloseq object without tree #####
# Merge all into a phyloseq object
# mpt <- phyloseq(OTU, SAMP, TAX, phylotree) # with phylo tree
mpt <- phyloseq(OTU, SAMP, TAX) # phyloseq object WITHOUT phylotree


#### Chinese Phyloseq Object #####
##### Format OTU table #####
# Format cn_otu table and read into a phyloseq object
cn_otu_mat <- as.matrix(cn_otu[,-1])
rownames(cn_otu_mat) <- cn_otu$`#OTU ID`
cn_OTU <- otu_table(cn_otu_mat, taxa_are_rows = TRUE)
class(cn_OTU) 

##### Format metadata #####
# Format cn metadata and read into a phyloseq object
cn_samp_df <- as.data.frame(cn_meta[,-1])
rownames(cn_samp_df)<- cn_meta$'sample-id'
cn_SAMP <- sample_data(cn_samp_df)
class(cn_SAMP)

##### Formatting taxonomy #####
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
##### Create cn_phyloseq object without tree #####
# cn_phyloseq_obj <- phyloseq(cn_OTU, cn_SAMP, cn_TAX, cn_tree) # with phylo tree
cn_phyloseq_obj <- phyloseq(cn_OTU, cn_SAMP, cn_TAX) # WITHOUT phylo tree

#### Merging us_and_cn phyloseq objects ######
ps_merged <- merge_phyloseq(mpt, cn_phyloseq_obj)
# if with phylo_tree:
# ERROR: Error in FUN(X[[i]], ...) : one tree has a different number of tips


#### Integrity checks #### 
ps_merged #merged WITHOUT phylo trees

# Expect these:
nsamples(ps_merged)                  # total samples = n_mpt + n_cn (90+33=123)



