#### Load in the 'phyloseq', 'ape', 'tidyverse', and 'vegan' packages ####
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

## The taxonomy.qza and rooted-tree.qza files were already exported from the server into a tsv and nwk file respectively ##
## The feature-table was converted from biom to txt file and exported from the server ##
## All these files are found in the us_exports folder which is found in the project_2_team_10 folder in the qiime_pipeline_uk folder ##

#### Load in the us_metadata_new.tsv file,  OTU table, taxonomy file, and phylogenetic tree ####
metaFP <- file.path("metadata", "us_metadata_new.tsv")
meta <- read_delim(file=metaFP, delim = "\t")

otuFP <- file.path("qiime_pipeline_uk", "project_2_team_10", "us_exports", "table_export", "feature-table.txt")
otu <- read_delim(file=otuFP, delim = "\t", skip=1)

taxfp <- file.path("qiime_pipeline_uk", "project_2_team_10", "us_exports", "taxonomy_export", "taxonomy.tsv")
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- file.path("qiime_pipeline_uk", "project_2_team_10", "us_exports","rooted-tree_export", "tree.nwk")
phylotree <- read_tree(phylotreefp)
class(phylotree)

#### Format sample metadata ####
# Save everything except sample-id as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sample-id the rownames
rownames(samp_df) <- meta$ "sample-id"
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Format OTU table ####
# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
class(OTU)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX) 

#### Create phyloseq object ####
# Merge all into a phyloseq object
mpt <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(mpt)
sample_data(mpt)
tax_table(mpt)
phy_tree(mpt)
