##Functional analysis comparing e-cigarettes vs control in us samples -->  
#diffferential abundance analysis to start -> errobar plots and heatmaps on ko terms

## Loading packages

library(tidyverse)
library(phyloseq)
library(ggpicrust2)
library(dplyr)
library(MicrobiomeStat)
library(ggplot2)
library(stringr)
library(pheatmap)


## Load Data
### Metadata can be extracted from phyloseq object -> copy is in this directory
###Merged and rarefied phyloseq was added

load("rarefied_phyloseq.RData")
meta<- rarefied_phyloseq@sam_data |>
  data.frame() |>
  rownames_to_column ('sample_name')

#Using KO results for this set of analysis 

ko <- read.delim('usa/KO_metagenome_out/pred_metagenome_unstrat.tsv', row.names = 1)

###Differential Analysis###

#Filter only to compare us control vs ecigs 

ecig_v_cont_meta <- meta |>
  filter(Combined %in% c("us_ecigarettes", "us_control"))

ecig_v_cont_ko <- ko |>
  select (all_of(ecig_v_cont_meta$sample_name))

####Run Differential abundance analysis (DAA)
ecig_v_cont_DAA <- pathway_daa(abundance = ecig_v_cont_ko,
                               metadata = ecig_v_cont_meta,
                               group = "Combined",
                               daa_method = "LinDA",
                               select = NULL, reference = NULL)

#Annotate pathway results using KO to KEGG conversion 
#Only need to run pathway_annotate once, then save it on folder (hence why its commented out)
ecig_v_cont_daa_annotated_results_df <- pathway_annotation(pathway = "KO", 
                                                           daa_results_df = ecig_v_cont_DAA, 
                                                           ko_to_kegg = TRUE)

#saveRDS(ecig_v_cont_daa_annotated_results_df, 'us_ecigs_v_cont_results/ecig_v_cont_daa_annotated_results_df.rds')

ecig_v_cont_daa_annotated_results_df = readRDS('us_ecigs_v_cont_results/ecig_v_cont_daa_annotated_results_df.rds')

##Generate pathway error bar plot
#Loading R script provided by Avril withthe fixed version for plotting errorbars
source('ggpicrust2_errorbar_function_fixed.R')

#first check how many hits u get, dont over annotate, change p-value and log2FoldChange if needed 
ecig_v_cont_daa_signif <- ecig_v_cont_daa_annotated_results_df|>
  dplyr::filter (
    p_adjust < 0.0005, 
    abs(log2FoldChange) > 5 )

##if I do >4 is 81, >5 is 13 
#nrow(ecig_v_cont_daa_signif)

ecig_v_cont_peb <- pathway_errorbar_fixed(abundance = ecig_v_cont_ko, 
                                          daa_results_df = ecig_v_cont_daa_signif, 
                                          Group = ecig_v_cont_meta$Combined, 
                                          wrap_label = T, wraplength=60,
                                          fc_cutoff = 5, order_by_log = F,
                                          p_values_threshold = 0.0005, 
                                          order = "pathway_class", 
                                          ko_to_kegg = TRUE, 
                                          p_value_bar = TRUE, 
                                          x_lab = "pathway_name")

ecig_v_cont_peb

ggsave("us_ecigs_v_cont_results/Plots/KEGG Error Bar Fixed.png",
       plot = ecig_v_cont_peb, width = 20, height = 10)

###Making Heat Map 
# pulling out the significant features that we worked with oreviously 
sig_features <- ecig_v_cont_daa_annotated_results_df |>
  filter(p_adjust < 0.0005, abs(log2FoldChange) > 5) |>
  pull(feature) |> unique()

ko_relab = read.delim('usa/KO_metagenome_out/pred_metagenome_unstrat.tsv',row.names = 1) |> 
  apply(2,function(x) x/sum(x)) |>
  as.data.frame()
colSums(ko_relab) # Should be all 1

stats <- ko_relab |>
  t() |>  # switch rows and columns
  as.data.frame() |> # Lets us use select()
  select(all_of(sig_features)) |> 
  cor(method = 'pearson') #Pearson tests between every feature pair

# Define the color palette
color_palette <- colorRampPalette(c("blue", "white", "red"))(40)

# Define breaks for the color scale
breaks <- seq(-1, 1, length.out = 41)  # 40 colors + 1 for the endpoint

# Create the clustered heatmap with centered color at zero
pheatmap(stats, 
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",          
         color = color_palette, 
         breaks = breaks,
         main = '',
         fontsize_row = 10, 
         fontsize_col = 10,
         filename = "us_ecigs_v_cont_results/Plots/KEGG Heatmap.png",
         height= 9, width = 9)

# Generate pathway PCA plot
p <- pathway_pca(abundance = ecig_v_cont_ko,
                 metadata = ecig_v_cont_meta, 
                 group = "Combined")

ggsave("us_ecigs_v_cont_results/Plots/Pathway_PCA.png", plot = p, width = 8, height = 6)

