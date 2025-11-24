##Functional analysis comparing e-cigarettes vs control in chinese samples --> do 
#diffferential abundance analysis to start -> errobar plots and heatmaps  

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

ko <- read.delim('china/KO_metagenome_out/pred_metagenome_unstrat.tsv', row.names = 1)

###Differential Analysis###

#Filter only to compare china control vs ecigs 

ecig_v_cont_meta <- meta |>
  filter(Combined %in% c("cn_ecigarettes", "cn_control"))

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
##ecig_v_cont_daa_annotated_results_df <- pathway_annotation(pathway = "KO", 
                                               daa_results_df = ecig_v_cont_DAA, 
                                               ko_to_kegg = TRUE)

##saveRDS(ecig_v_cont_daa_annotated_results_df, 'cn_ecigs_v_cont_results/ecig_v_cont_daa_annotated_results_df.rds')

ecig_v_cont_daa_annotated_results_df = readRDS('cn_ecigs_v_cont_results/ecig_v_cont_daa_annotated_results_df.rds')

##Generate pathway error bar plot

#first check how many hits u get, dont over annotate, change p-value and log2FoldChange if needed 
ecig_v_cont_daa_signif <- ecig_v_cont_daa_annotated_results_df|>
  dplyr::filter (
    p_adjust < 0.0005, 
    abs(log2FoldChange) > 5 )

##if I do >4 is 81, >5 is 13 
  #nrow(ecig_v_cont_daa_signif)

ecig_v_cont_paths <- pathway_errorbar(abundance = ecig_v_cont_ko,
                                 daa_results_df = ecig_v_cont_daa_signif, 
                                 Group = ecig_v_cont_meta$Combined,
                                 p_values_threshold = 0.0005,
                                 order = "pathway_class",
                                 ko_to_kegg = TRUE,
                                 p_value_bar = TRUE,
                                 x_lab = "pathway_name") 

# Making it look better --> more to be added later on
paths <- ecig_v_cont_paths +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    strip.text = element_text(face = "bold"))

#ggsave("cn_ecigs_v_cont_results/KEGG Error Bar Fixed.png")

###Making Heat Map 
# pulling out the significant features that we worked with oreviously 
sig_features <- ecig_v_cont_daa_annotated_results_df |>
  filter(p_adjust < 0.0005, abs(log2FoldChange) > 5) |>
  pull(feature)

#feeding only the important KO that we have been using into the pipeline
stats <- ecig_v_cont_ko |>
  t() |> 
  as.data.frame ()|> 
  select(all_of(sig_features))|>
  cor (method = 'spearman')

map <- pheatmap(stats, 
                clustering_distance_rows = "euclidean",
                clustering_distance_cols = "euclidean",
                clustering_method = "complete",
                fontsize_row = 10, fontsize_col = 10)

  
  
