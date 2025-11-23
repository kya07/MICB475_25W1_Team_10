#!/usr/bin/env RScript 
# MICB 475 Relative Abundance analysis/ Abundance Modeling using DESeq2 RScript
# Determine the # of points abs log2FoldChange above/below 2 and
# Determine which phylum most abundant
# 18 November, 2025

# DESeq is a parametric statistical tool that uses a negative binomial distribution to model read counts
# It is powerful, but does not handle zero-inflated data well.
# DESeq2 visualized with volcano plot and bar plots

#### Load Packages ####
library(tidyverse)
library(phyloseq)
library(DESeq2) 

#### Load Data (non-rarefied) ####
filtered_phyloseq <- get(load("filtered_phyloseq.RData")) # get() around load() gives the actual phyloseq object

#### DESeq ####
##### Count Normalization ##### 
## Above returned zeros error, thus have to add '1' count to all reads
# Transform counts to add a pseudocount of 1
merged_plus1 <- transform_sample_counts(filtered_phyloseq, function(x) x+1)
merged_deseq <- phyloseq_to_deseq2(merged_plus1, ~`Combined`)
DESEQ_merged <- DESeq(merged_deseq)

##### WithIN dataset DESeq2 #####
###### US dataset(1-3) ###### 
us_ecigarettes_vs_us_control <- results(DESEQ_merged, tidy=TRUE, 
                                        #c(condition, comparison/treatment group, control/reference)
                                        contrast = c("Combined","us_ecigarettes","us_control"))

us_tobacco_vs_us_control <- results(DESEQ_merged, tidy=TRUE, 
                                    contrast = c("Combined","us_tobacco","us_control"))

us_tobacco_vs_us_ecigarettes <- results(DESEQ_merged, tidy=TRUE, 
                                        contrast = c("Combined","us_tobacco","us_ecigarettes"))

###### CN dataset (4-6) ###### 
cn_ecigarettes_vs_cn_control <- results(DESEQ_merged, tidy=TRUE, 
                                        contrast = c("Combined","cn_ecigarettes","cn_control"))

cn_tobacco_vs_cn_control <- results(DESEQ_merged, tidy=TRUE, 
                                    contrast = c("Combined","cn_tobacco","cn_control"))

cn_tobacco_vs_cn_ecigarettes <- results(DESEQ_merged, tidy=TRUE, 
                                        contrast = c("Combined","cn_tobacco","cn_ecigarettes"))

##### Between datasets DESeq2 (7-9)##### 
cn_control_vs_us_control <- results(DESEQ_merged, tidy=TRUE, 
                                    contrast = c("Combined","cn_control","us_control"))

cn_ecigarettes_vs_us_ecigarettes <- results(DESEQ_merged, tidy=TRUE, 
                                            contrast = c("Combined","cn_ecigarettes","us_ecigarettes"))

cn_tobacco_vs_us_tobacco <- results(DESEQ_merged, tidy=TRUE, 
                                    contrast = c("Combined","cn_tobacco","us_tobacco"))


#### Bar Plot Visualization (Phylum level)####
##### WithIN dataset #####
###### U.S. Bar Plots with ASV on x-axis (1-3) ######
# To get table of results
sigASVs_us_ecigarettes_vs_us_control <- us_ecigarettes_vs_us_control  |>  
  filter(padj<0.01 & abs(log2FoldChange)>2) |>    
  dplyr::rename(ASV=row) # rename column "row" to be called ASV

# Vector of significant ASVs
sigASVs_us_ecigarettes_vs_us_control_vec <- sigASVs_us_ecigarettes_vs_us_control |> 
  pull(ASV)

# Prune phyloseq file
us_ecigarettes_vs_us_control_merged_DESeq <- prune_taxa(sigASVs_us_ecigarettes_vs_us_control_vec,filtered_phyloseq)

# Add taxonomy onto DESeq results table
merged_results_us_ecigarettes_vs_us_control <- tax_table(us_ecigarettes_vs_us_control_merged_DESeq) |> 
  as.data.frame() |> 
  rownames_to_column(var="ASV")  |> 
  right_join(sigASVs_us_ecigarettes_vs_us_control, by = "ASV")  |> 
  arrange(log2FoldChange) |> # in Ascending order 
  mutate(Phylum = make.unique(Phylum)) |>  #only list unique Phylum to avoid repeating Genera
  mutate(Phylum = str_remove(Phylum, "^p__"))  |>  # Clean Phylum column
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

# Generate Horizontal Bar Plot
us_ecigarettes_vs_us_control_bar_plot <- ggplot(merged_results_us_ecigarettes_vs_us_control) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+ 
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + # error bar included 
  ggtitle("U.S. E-Cigarette Users vs. U.S. Controls (padj <0.01) ") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ 
  coord_flip() +
  theme_classic()
us_ecigarettes_vs_us_control_bar_plot

# ggsave("us_ecigarettes_vs_us_control_bar_plot.png",
#        plot = us_ecigarettes_vs_us_control_bar_plot,
#        width = 10, height = 7)

# Plot 2 us_tobacco_vs_us_control
sigASVs_us_tobacco_vs_us_control <- us_tobacco_vs_us_control  |>  
  filter(padj<0.01 & abs(log2FoldChange)>2) |>   
  dplyr::rename(ASV=row) 

# Vector of significant ASVs
sigASVs_us_tobacco_vs_us_control_vec <- sigASVs_us_tobacco_vs_us_control |> 
  pull(ASV)

# Prune phyloseq file
us_tobacco_vs_us_control_merged_DESeq <- prune_taxa(sigASVs_us_tobacco_vs_us_control_vec,filtered_phyloseq)

# Add taxonomy onto DESeq results table
merged_results_us_tobacco_vs_us_control <- tax_table(us_tobacco_vs_us_control_merged_DESeq) |>  
  as.data.frame() |> 
  rownames_to_column(var="ASV")  |> 
  right_join(sigASVs_us_tobacco_vs_us_control, by = "ASV")  |> 
  arrange(log2FoldChange) |>   
  mutate(Phylum = make.unique(Phylum)) |>  
  mutate(Phylum = str_remove(Phylum, "^p__"))  |> 
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

# Generate Horizontal Bar Plot
us_tobacco_vs_us_control_bar_plot <- ggplot(merged_results_us_tobacco_vs_us_control) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+ 
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +  
  ggtitle("U.S. Tobacco Users vs. U.S. Controls (padj <0.01) ") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ 
  coord_flip() +
  theme_classic()
us_tobacco_vs_us_control_bar_plot

# ggsave("us_tobacco_vs_us_control_bar_plot.png",
#        plot = us_tobacco_vs_us_control_bar_plot,
#        width = 10, height = 7)

# Plot 3 us_tobacco_vs_us_ecigarettes
sigASVs_us_tobacco_vs_us_ecigarettes <- us_tobacco_vs_us_ecigarettes  |>  
  filter(padj<0.01 & abs(log2FoldChange)>2) |>  
  dplyr::rename(ASV=row) 

# Vector of significant ASVs
sigASVs_us_tobacco_vs_us_ecigarettes_vec <- sigASVs_us_tobacco_vs_us_ecigarettes |> 
  pull(ASV)

# Prune phyloseq file
us_tobacco_vs_us_ecigarettes_merged_DESeq <- prune_taxa(sigASVs_us_tobacco_vs_us_ecigarettes_vec,filtered_phyloseq)

# Add taxonomy onto DESeq results table
merged_results_us_tobacco_vs_us_ecigarettes <- tax_table(us_tobacco_vs_us_ecigarettes_merged_DESeq) |>  
  as.data.frame() |> 
  rownames_to_column(var="ASV")  |> 
  right_join(sigASVs_us_tobacco_vs_us_ecigarettes, by = "ASV")  |> 
  arrange(log2FoldChange) |> # in Ascending order 
  mutate(Phylum = make.unique(Phylum)) |> 
  mutate(Phylum = str_remove(Phylum, "^p__"))  |> 
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

# Generate Horizontal Bar Plot
us_tobacco_vs_us_ecigarettes_bar_plot <- ggplot(merged_results_us_tobacco_vs_us_ecigarettes) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+ 
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + # error bar included 
  ggtitle("U.S. Tobacco Users vs. U.S. E-Cigarette Users (padj <0.01) ") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ 
  coord_flip() +
  theme_classic()
us_tobacco_vs_us_ecigarettes_bar_plot

# ggsave("us_tobacco_vs_us_ecigarettes_bar_plot.png",
#        plot = us_tobacco_vs_us_ecigarettes_bar_plot,
#        width = 10, height = 7)


###### CN Bar Plots with ASV on x-axis (4-6) ######
# Plot 4 cn_ecigarettes_vs_cn_control
sigASVs_cn_ecigarettes_vs_cn_control <- cn_ecigarettes_vs_cn_control  |>  
  filter(padj<0.01 & abs(log2FoldChange)>2) |>  
  dplyr::rename(ASV=row) 

# Vector of significant ASVs
sigASVs_cn_ecigarettes_vs_cn_control_vec <- sigASVs_cn_ecigarettes_vs_cn_control |> 
  pull(ASV)

# Prune phyloseq file
cn_ecigarettes_vs_cn_control_merged_DESeq <- prune_taxa(sigASVs_cn_ecigarettes_vs_cn_control_vec,filtered_phyloseq)

# Add taxonomy onto DESeq results table
merged_results_cn_ecigarettes_vs_cn_control <- tax_table(cn_ecigarettes_vs_cn_control_merged_DESeq) |>  
  as.data.frame() |> 
  rownames_to_column(var="ASV")  |> 
  right_join(sigASVs_cn_ecigarettes_vs_cn_control, by = "ASV")  |> 
  arrange(log2FoldChange) |>  
  mutate(Phylum = make.unique(Phylum)) |>  
  mutate(Phylum = str_remove(Phylum, "^p__"))  |> 
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

# Generate Horizontal Bar Plot
cn_ecigarettes_vs_cn_control_bar_plot <- ggplot(merged_results_cn_ecigarettes_vs_cn_control) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+ 
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + # error bar included 
  ggtitle("CN E-Cigarette Users vs. CN Controls (padj <0.01) ") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ 
  coord_flip() +
  theme_classic()
cn_ecigarettes_vs_cn_control_bar_plot

# ggsave("cn_ecigarettes_vs_cn_control_bar_plot.png",
#        plot = cn_ecigarettes_vs_cn_control_bar_plot,
#        width = 10, height = 7)

# Plot 5 cn_tobacco_vs_cn_control
sigASVs_cn_tobacco_vs_cn_control <- cn_tobacco_vs_cn_control  |>  
  filter(padj<0.01 & abs(log2FoldChange)>2) |>  
  dplyr::rename(ASV=row) 

# Vector of significant ASVs
sigASVs_cn_tobacco_vs_cn_control_vec <- sigASVs_cn_tobacco_vs_cn_control |> 
  pull(ASV)

# Prune phyloseq file
cn_tobacco_vs_cn_control_merged_DESeq <- prune_taxa(sigASVs_cn_tobacco_vs_cn_control_vec,filtered_phyloseq)

# Add taxonomy onto DESeq results table
merged_results_cn_tobacco_vs_cn_control <- tax_table(cn_tobacco_vs_cn_control_merged_DESeq) |>  
  as.data.frame() |> 
  rownames_to_column(var="ASV")  |> 
  right_join(sigASVs_cn_tobacco_vs_cn_control, by = "ASV")  |> 
  arrange(log2FoldChange) |> # in Ascending order 
  mutate(Phylum = make.unique(Phylum)) |>  
  mutate(Phylum = str_remove(Phylum, "^p__"))  |> 
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

# Generate Horizontal Bar Plot
cn_tobacco_vs_cn_control_bar_plot <- ggplot(merged_results_cn_tobacco_vs_cn_control) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+ 
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + # error bar included 
  ggtitle("CN Tobacco Users vs. CN Controls (padj <0.01) ") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ 
  coord_flip() +
  theme_classic()
cn_tobacco_vs_cn_control_bar_plot

# ggsave("cn_tobacco_vs_cn_control_bar_plot.png",
#        plot = cn_tobacco_vs_cn_control_bar_plot,
#        width = 12, height = 8)

# Plot 6 cn_tobacco_vs_cn_ecigarettes
sigASVs_cn_tobacco_vs_cn_ecigarettes <- cn_tobacco_vs_cn_ecigarettes  |>  
  filter(padj<0.01 & abs(log2FoldChange)>2) |>  
  dplyr::rename(ASV=row) 

# Vector of significant ASVs
sigASVs_cn_tobacco_vs_cn_ecigarettes_vec <- sigASVs_cn_tobacco_vs_cn_ecigarettes |> 
  pull(ASV)

# Prune phyloseq file
cn_tobacco_vs_cn_ecigarettes_merged_DESeq <- prune_taxa(sigASVs_cn_tobacco_vs_cn_ecigarettes_vec,filtered_phyloseq)

# Add taxonomy onto DESeq results table
merged_results_cn_tobacco_vs_cn_ecigarettes <- tax_table(cn_tobacco_vs_cn_ecigarettes_merged_DESeq) |>  
  as.data.frame() |> 
  rownames_to_column(var="ASV")  |> 
  right_join(sigASVs_cn_tobacco_vs_cn_ecigarettes, by = "ASV")  |>  # right_join(sigASVs_us_tobacco_vs_us_control, by = "ASV")  |> 
  arrange(log2FoldChange) |> # in Ascending order 
  mutate(Phylum = make.unique(Phylum)) |> 
  mutate(Phylum = str_remove(Phylum, "^p__"))  |> 
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

# Generate Horizontal Bar Plot
cn_tobacco_vs_cn_ecigarettes_bar_plot <- ggplot(merged_results_cn_tobacco_vs_cn_ecigarettes) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+ 
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + # error bar included 
  ggtitle("CN Tobacco Users vs. CN E-Cigarette Users (padj <0.01) ") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ 
  coord_flip() +
  theme_classic()
cn_tobacco_vs_cn_ecigarettes_bar_plot

# ggsave("cn_tobacco_vs_cn_ecigarettes_bar_plot.png",
#        plot = cn_tobacco_vs_cn_ecigarettes_bar_plot,
#        width = 12, height = 8)


##### Between datasets (7-9) #####
# Plot 7 cn_control_vs_us_control
sigASVs_cn_control_vs_us_control <- cn_control_vs_us_control  |>  
  filter(padj<0.01 & abs(log2FoldChange)>2) |>  
  dplyr::rename(ASV=row) 

# Vector of significant ASVs
sigASVs_cn_control_vs_us_control_vec <- sigASVs_cn_control_vs_us_control |> 
  pull(ASV)

# Prune phyloseq file
cn_control_vs_us_control_merged_DESeq <- prune_taxa(sigASVs_cn_control_vs_us_control_vec,filtered_phyloseq)

# Add taxonomy onto DESeq results table
merged_results_cn_control_vs_us_control <- tax_table(cn_control_vs_us_control_merged_DESeq) |>  
  as.data.frame() |> 
  rownames_to_column(var="ASV")  |> 
  right_join(sigASVs_cn_control_vs_us_control, by = "ASV")  |> 
  arrange(log2FoldChange) |>
  mutate(Phylum = make.unique(Phylum)) |> 
  mutate(Phylum = str_remove(Phylum, "^p__"))  |> 
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

# Generate Horizontal Bar Plot
cn_control_vs_us_control_bar_plot <- ggplot(merged_results_cn_control_vs_us_control) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+ 
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + # error bar included 
  ggtitle("CN Controls vs. U.S. Controls (padj <0.01) ") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ 
  coord_flip() +
  theme_classic()
cn_control_vs_us_control_bar_plot

# ggsave("cn_control_vs_us_control_bar_plot.png",
#        plot = cn_control_vs_us_control_bar_plot,
#        width = 20, height = 22)

# Plot 8 cn_ecigarettes_vs_us_ecigarettes
sigASVs_cn_ecigarettes_vs_us_ecigarettes <- cn_ecigarettes_vs_us_ecigarettes  |>  
  filter(padj<0.01 & abs(log2FoldChange)>2) |>  
  dplyr::rename(ASV=row) 

# Vector of significant ASVs
sigASVs_cn_ecigarettes_vs_us_ecigarettes_vec <- sigASVs_cn_ecigarettes_vs_us_ecigarettes |> 
  pull(ASV)

# Prune phyloseq file
cn_ecigarettes_vs_us_ecigarettes_merged_DESeq <- prune_taxa(sigASVs_cn_ecigarettes_vs_us_ecigarettes_vec,filtered_phyloseq)

# Add taxonomy onto DESeq results table
merged_results_cn_ecigarettes_vs_us_ecigarettes <- tax_table(cn_ecigarettes_vs_us_ecigarettes_merged_DESeq) |>  
  as.data.frame() |> 
  rownames_to_column(var="ASV")  |> 
  right_join(sigASVs_cn_ecigarettes_vs_us_ecigarettes, by = "ASV")  |> 
  arrange(log2FoldChange) |>  
  mutate(Phylum = make.unique(Phylum)) |> 
  mutate(Phylum = str_remove(Phylum, "^p__"))  |> 
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

# Generate Horizontal Bar Plot
cn_ecigarettes_vs_us_ecigarettes_bar_plot <- ggplot(merged_results_cn_ecigarettes_vs_us_ecigarettes) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+ 
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + # error bar included 
  ggtitle("CN E-Cigarette Users vs. U.S. E-Cigarette Users (padj <0.01) ") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ 
  coord_flip() +
  theme_classic()
cn_ecigarettes_vs_us_ecigarettes_bar_plot

# ggsave("cn_ecigarettes_vs_us_ecigarettes_bar_plot.png",
#        plot = cn_ecigarettes_vs_us_ecigarettes_bar_plot,
#        width = 20, height = 22)


# Plot 9 cn_tobacco_vs_us_tobacco
sigASVs_cn_tobacco_vs_us_tobacco <- cn_tobacco_vs_us_tobacco  |>  
  filter(padj<0.01 & abs(log2FoldChange)>2) |>  
  dplyr::rename(ASV=row) 

# Vector of significant ASVs
sigASVs_cn_tobacco_vs_us_tobacco_vec <- sigASVs_cn_tobacco_vs_us_tobacco |> 
  pull(ASV)

# Prune phyloseq file
cn_tobacco_vs_us_tobacco_merged_DESeq <- prune_taxa(sigASVs_cn_tobacco_vs_us_tobacco_vec,filtered_phyloseq)

# Add taxonomy onto DESeq results table
merged_results_cn_tobacco_vs_us_tobacco <- tax_table(cn_tobacco_vs_us_tobacco_merged_DESeq) |>  
  as.data.frame() |> 
  rownames_to_column(var="ASV")  |> 
  right_join(sigASVs_cn_tobacco_vs_us_tobacco, by = "ASV")  |> 
  arrange(log2FoldChange) |> # in Ascending order 
  mutate(Phylum = make.unique(Phylum)) |>  
  mutate(Phylum = str_remove(Phylum, "^p__"))  |> 
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

# Generate Horizontal Bar Plot
cn_tobacco_vs_us_tobacco_bar_plot <- ggplot(merged_results_cn_tobacco_vs_us_tobacco) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+ 
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + # error bar included 
  ggtitle("CN Tobacco Users vs. U.S. Tobacco Users (padj <0.01) ") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ 
  coord_flip() +
  theme_classic()
cn_tobacco_vs_us_tobacco_bar_plot

# ggsave("cn_tobacco_vs_us_tobacco_bar_plot.png",
#        plot = cn_tobacco_vs_us_tobacco_bar_plot,
#        width = 30, height = 35)


#### Determining the # of points abs log2FoldChange above/below 2#### 
##### Within U.S. or China and Between U.S. vs. China##### 
# Add a new column to each merged results table to indicate Smoking Description and Cohort
merged_results_us_ecigarettes_vs_us_control <- merged_results_us_ecigarettes_vs_us_control |> 
  mutate(Type = "us_ecigarettes_vs_us_control")

merged_results_us_tobacco_vs_us_control <- merged_results_us_tobacco_vs_us_control |> 
  mutate(Type = "us_tobacco_vs_us_control")

merged_results_us_tobacco_vs_us_ecigarettes <- merged_results_us_tobacco_vs_us_ecigarettes |> 
  mutate(Type = "us_tobacco_vs_us_ecigarettes")

merged_results_cn_ecigarettes_vs_cn_control <- merged_results_cn_ecigarettes_vs_cn_control |> 
  mutate(Type = "cn_ecigarettes_vs_cn_control")

merged_results_cn_tobacco_vs_cn_control <- merged_results_cn_tobacco_vs_cn_control |> 
  mutate(Type = "cn_tobacco_vs_cn_control")

merged_results_cn_tobacco_vs_cn_ecigarettes <- merged_results_cn_tobacco_vs_cn_ecigarettes |> 
  mutate(Type = "cn_tobacco_vs_cn_ecigarettes")

merged_results_cn_control_vs_us_control <- merged_results_cn_control_vs_us_control |> 
  mutate(Type = "cn_control_vs_us_control")
  
merged_results_cn_ecigarettes_vs_us_ecigarettes <- merged_results_cn_ecigarettes_vs_us_ecigarettes |> 
  mutate(Type = "cn_ecigarettes_vs_us_ecigarettes")

merged_results_cn_tobacco_vs_us_tobacco <- merged_results_cn_tobacco_vs_us_tobacco |> 
  mutate(Type = "cn_tobacco_vs_us_tobacco")

# Combine all three data frames into one
combined_results_type <- bind_rows(merged_results_us_ecigarettes_vs_us_control, 
                                   merged_results_us_tobacco_vs_us_control, 
                                   merged_results_us_tobacco_vs_us_ecigarettes, 
                                   merged_results_cn_ecigarettes_vs_cn_control, 
                                   merged_results_cn_tobacco_vs_cn_control, 
                                   merged_results_cn_tobacco_vs_cn_ecigarettes,
                                   merged_results_cn_control_vs_us_control, 
                                   merged_results_cn_ecigarettes_vs_us_ecigarettes,
                                   merged_results_cn_tobacco_vs_us_tobacco)

# Calculate mean and standard deviation for each Phylum and Type (cohort and smoking description)
count_values_type <- combined_results_type  |> 
  group_by(Type)  |> 
  summarise(
    count_above_2 = sum(log2FoldChange > 2, na.rm = TRUE),
    count_below_2 = sum(log2FoldChange < 2, na.rm = TRUE))
count_values_type

# Save the merged results as a CSV file
# write.csv(count_values_type, "DESeq_graphvalues_combined.csv", row.names = FALSE)


#### Determining which phylum most abundant #### 
# Clean Phylum by removing the ".number" part
phylum_abundance_type <- combined_results_type  |> 
  mutate(Phylum = str_remove(Phylum, "\\..*")) |>   # Remove ".number" part
  group_by(Type, Phylum) |> 
  summarise(
    count_above_2 = sum(log2FoldChange > 2, na.rm = TRUE),
    count_below_2 = sum(log2FoldChange < 2, na.rm = TRUE))
phylum_abundance_type

# Save the merged results as a CSV file
# write.csv(phylum_abundance_type, "DESeq_abundvalues_cohort_status.csv", row.names = FALSE)
