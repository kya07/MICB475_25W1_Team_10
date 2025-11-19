#!/usr/bin/env RScript 
# MICB 475 Relative Abundance analysis/ Abundance Modeling using DESeq2 RScript
# Bar plots resolved to Phylum level
# 07 November, 2025

# DESeq is a parametric statistical tool that uses a negative binomial distribution to model read counts
# It is powerful, but does not handle zero-inflated data well.
# DESeq2 visualized with volcano plot and bar plots

#### Load Packages ####
library(tidyverse)
library(phyloseq)
library(DESeq2) 

#### Load Data (non-rarefied) ####
filtered_phyloseq <- get(load("filtered_phyloseq.RData")) # IMPORTANT: get() around load() gives the actual phyloseq object
print(filtered_phyloseq)

#### DESeq ####
# # Covert phyloseq object to DESeq object
# merged_deseq <- phyloseq_to_deseq2(filtered_phyloseq, ~`Combined`)
# DESEQ_merged<- DESeq(merged_deseq)

##### Count Normalization ##### 
## Above returned zeros error, thus have to add '1' count to all reads
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


#### Volcano Plot Visualization ####
## Volcano plot: effect size VS significance
## Make variable to color by whether it is significant + large change

##### WithIN dataset #####
###### US_volcano_plots (1-3) ######
us_ecigarettes_vs_us_control_vol_plot <- us_ecigarettes_vs_us_control |> 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |>  # new significant column created with set thresholds 
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) + # using the -log10() function to transform padj values
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + # set boundary lines on x- and y-axes to plot data into quadrants
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  xlab("Log2FoldChange") + # change axes names
  ylab("-Log10(Adjusted P-Value)") +
  ggtitle("U.S. E-Cigarette Users vs. U.S. Controls (padj < 0.01, |log2FC| > 2)") +
  theme_test() + # apply minimalist theme, grey background removed
  theme(
    plot.title = element_text(face = "bold", size = 15), 
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), # bold the axis titles
    panel.border = element_rect(linewidth = 2)) # border the graph 
us_ecigarettes_vs_us_control_vol_plot

us_tobacco_vs_us_control_vol_plot <- us_tobacco_vs_us_control |> 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |> 
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) + 
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + 
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  xlab("Log2FoldChange") + 
  ylab("-Log10(Adjusted P-Value)") +
  ggtitle("U.S. Tobacco Users vs. U.S. Controls (padj < 0.01, |log2FC| > 2)") +
  theme_test() + 
  theme(
    plot.title = element_text(face = "bold", size = 15), 
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), 
    panel.border = element_rect(linewidth = 2)) 
us_tobacco_vs_us_control_vol_plot


us_tobacco_vs_us_ecigarettes_vol_plot <- us_tobacco_vs_us_ecigarettes |> 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |> 
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) + 
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + 
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  xlab("Log2FoldChange") + 
  ylab("-Log10(Adjusted P-Value)") +
  ggtitle("U.S. Tobacco Users vs. U.S. E-Cigarette Users (padj < 0.01, |log2FC| > 2)") +
  theme_test() + 
  theme(
    plot.title = element_text(face = "bold", size = 15), 
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), 
    panel.border = element_rect(linewidth = 2)) 
us_tobacco_vs_us_ecigarettes_vol_plot

###### CN_volcano_plots (4-6) ######
cn_ecigarettes_vs_cn_control_vol_plot <- cn_ecigarettes_vs_cn_control |> 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |>   
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) +  
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +  
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  xlab("Log2FoldChange") + # change axes names
  ylab("-Log10(Adjusted P-Value)") +
  ggtitle("CN E-Cigarette Users vs. CN Controls (padj < 0.01, |log2FC| > 2)") +
  theme_test() +  
  theme(
    plot.title = element_text(face = "bold", size = 15), 
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),  
    panel.border = element_rect(linewidth = 2))  
cn_ecigarettes_vs_cn_control_vol_plot

cn_tobacco_vs_cn_control_vol_plot <- cn_tobacco_vs_cn_control |> 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |> 
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) + 
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + 
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  xlab("Log2FoldChange") + 
  ylab("-Log10(Adjusted P-Value)") +
  ggtitle("CN Tobacco Users vs. CN Controls (padj < 0.01, |log2FC| > 2)") +
  theme_test() + 
  theme(
    plot.title = element_text(face = "bold", size = 15), 
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), 
    panel.border = element_rect(linewidth = 2)) 
cn_tobacco_vs_cn_control_vol_plot

cn_tobacco_vs_cn_ecigarettes_vol_plot <- cn_tobacco_vs_cn_ecigarettes |> 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |> 
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) + 
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + 
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  xlab("Log2FoldChange") + 
  ylab("-Log10(Adjusted P-Value)") +
  ggtitle("CN Tobacco Users vs. CN E-Cigarette Users (padj < 0.01, |log2FC| > 2)") +
  theme_test() + 
  theme(
    plot.title = element_text(face = "bold", size = 15), 
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), 
    panel.border = element_rect(linewidth = 2)) 
cn_tobacco_vs_cn_ecigarettes_vol_plot

##### Between datasets visualization (7-9)##### 
cn_control_vs_us_control_vol_plot <- cn_control_vs_us_control |> 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |>  # new significant column created with set thresholds 
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) + # using the -log10() function to transform padj values
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + # set boundary lines on x- and y-axes to plot data into quadrants
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  xlab("Log2FoldChange") + # change axes names
  ylab("-Log10(Adjusted P-Value)") +
  ggtitle("CN Controls vs. U.S. Controls (padj < 0.01, |log2FC| > 2)") +
  theme_test() + # apply minimalist theme, grey background removed
  theme(
    plot.title = element_text(face = "bold", size = 15), 
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), # bold the axis titles
    panel.border = element_rect(linewidth = 2)) # border the graph 
cn_control_vs_us_control_vol_plot

cn_ecigarettes_vs_us_ecigarettes_vol_plot <- cn_ecigarettes_vs_us_ecigarettes |> 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |>   
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) +  
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +  
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  xlab("Log2FoldChange") + # change axes names
  ylab("-Log10(Adjusted P-Value)") +
  ggtitle("CN E-Cigarette Users vs. U.S. E-Cigarette Users (padj < 0.01, |log2FC| > 2)") +
  theme_test() +  
  theme(
    plot.title = element_text(face = "bold", size = 15), 
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),  
    panel.border = element_rect(linewidth = 2))  
cn_ecigarettes_vs_us_ecigarettes_vol_plot

cn_tobacco_vs_us_tobacco_vol_plot <- cn_tobacco_vs_us_tobacco |> 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |> 
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) + 
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + 
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  xlab("Log2FoldChange") + 
  ylab("-Log10(Adjusted P-Value)") +
  ggtitle("CN Tobacco Users vs. U.S. Tobacco Users (padj < 0.01, |log2FC| > 2)") +
  theme_test() + 
  theme(
    plot.title = element_text(face = "bold", size = 15), 
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), 
    panel.border = element_rect(linewidth = 2)) 
cn_tobacco_vs_us_tobacco_vol_plot


##### Combine into a panel grid ##### 
# 1) Bind your 9 DESeq2 result tables and tag them
volc_all <- bind_rows(
  cn_ecigarettes_vs_cn_control |> mutate(dataset="China",      comparison="E-Cigarettes vs Controls"),
  cn_tobacco_vs_cn_control |> mutate(dataset="China",      comparison="Tobacco vs Controls"),
  cn_tobacco_vs_cn_ecigarettes |> mutate(dataset="China",      comparison="Tobacco vs E-Cigarettes"),
  
  us_ecigarettes_vs_us_control |> mutate(dataset="U.S.",      comparison="E-Cigarettes vs Controls"),
  us_tobacco_vs_us_control |> mutate(dataset="U.S.",      comparison="Tobacco vs Controls"),
  us_tobacco_vs_us_ecigarettes |> mutate(dataset="U.S.",      comparison="Tobacco vs E-Cigarettes"),
  
  cn_control_vs_us_control |> mutate(dataset="China vs. U.S.", comparison="Controls"),
  cn_ecigarettes_vs_us_ecigarettes |> mutate(dataset="China vs. U.S.", comparison="E-Cigarettes"),
  cn_tobacco_vs_us_tobacco |> mutate(dataset="China vs. U.S.", comparison="Tobacco")
) |> 
mutate(
  dataset    = factor(dataset, levels = c("China", "U.S.", "China vs. U.S.")),
  comparison = factor(comparison,
                      levels = c("E-Cigarettes vs Controls",
                                 "Tobacco vs Controls",
                                 "Tobacco vs E-Cigarettes",
                                 "Controls", "E-Cigarettes", "Tobacco"))) 
# above: 1.5) Lock facet order to 3x3 grid

# 2) Faceted volcano
facted_volcano_plot <- volc_all |> 
        mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |>  
        ggplot(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) + 
        geom_point(alpha = 0.6, size = 1) +  
        geom_hline(yintercept = -log10(0.01), linetype = "dashed") +  
        geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
        facet_grid(dataset ~ comparison) +
        xlab("Log2FoldChange") +  
        ylab("-Log10(Adjusted P-Value)") +
        ggtitle("Within-China, Within-U.S., and Between-Dataset DESeq2 Analysis (padj < 0.01, |log2FC| > 2)") +
        theme_test() +  
        theme(
          plot.title = element_text(face = "bold", size = 15, hjust = 0.5), 
          axis.text = element_text(colour = "black", size = 11),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),  
          panel.border = element_rect(linewidth = 2)) 
facted_volcano_plot

##### Save Volcano Plots #####
# put 9 individual volcano plots and 1 faceted plot in a named list
plots_list <- list(
  "us_ecigarettes_vs_us_control"     = us_ecigarettes_vs_us_control_vol_plot,
  "us_tobacco_vs_us_control"         = us_tobacco_vs_us_control_vol_plot,
  "us_tobacco_vs_us_ecigarettes"     = us_tobacco_vs_us_ecigarettes_vol_plot,
  "cn_ecigarettes_vs_cn_control"     = cn_ecigarettes_vs_cn_control_vol_plot,
  "cn_tobacco_vs_cn_control"         = cn_tobacco_vs_cn_control_vol_plot,
  "cn_tobacco_vs_cn_ecigarettes"     = cn_tobacco_vs_cn_ecigarettes_vol_plot,
  "cn_control_vs_us_control"         = cn_control_vs_us_control_vol_plot,
  "cn_ecigarettes_vs_us_ecigarettes" = cn_ecigarettes_vs_us_ecigarettes_vol_plot,
  "cn_tobacco_vs_us_tobacco"         = cn_tobacco_vs_us_tobacco_vol_plot, 
  "faceted_volcano_plot"             = facted_volcano_plot   # add the combined grid
)

# # create output folder for PNGs
dir.create("volcano_pngs_padj0.01", showWarnings = FALSE)

# Save each plot automatically
for (nm in names(plots_list)) {
  ggsave(
    filename = file.path("volcano_pngs_padj0.01", paste0(nm, ".png")),
    plot = plots_list[[nm]],
    width = 10.5, height = 5
  )
}


#### Bar Plot Visualization (Phylum level)####
# Plots not only the fold change but also the taxonomic information (on x-axis)

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

ggsave("us_ecigarettes_vs_us_control_bar_plot.png",
       plot = us_ecigarettes_vs_us_control_bar_plot,
       width = 10, height = 7)

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

ggsave("us_tobacco_vs_us_control_bar_plot.png",
       plot = us_tobacco_vs_us_control_bar_plot,
       width = 10, height = 7)

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

ggsave("us_tobacco_vs_us_ecigarettes_bar_plot.png",
       plot = us_tobacco_vs_us_ecigarettes_bar_plot,
       width = 10, height = 7)


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

ggsave("cn_ecigarettes_vs_cn_control_bar_plot.png",
       plot = cn_ecigarettes_vs_cn_control_bar_plot,
       width = 10, height = 7)


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

ggsave("cn_tobacco_vs_cn_control_bar_plot.png",
       plot = cn_tobacco_vs_cn_control_bar_plot,
       width = 12, height = 8)


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

ggsave("cn_tobacco_vs_cn_ecigarettes_bar_plot.png",
       plot = cn_tobacco_vs_cn_ecigarettes_bar_plot,
       width = 12, height = 8)


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

ggsave("cn_control_vs_us_control_bar_plot.png",
       plot = cn_control_vs_us_control_bar_plot,
       width = 20, height = 22)

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

ggsave("cn_ecigarettes_vs_us_ecigarettes_bar_plot.png",
       plot = cn_ecigarettes_vs_us_ecigarettes_bar_plot,
       width = 20, height = 22)


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

ggsave("cn_tobacco_vs_us_tobacco_bar_plot.png",
       plot = cn_tobacco_vs_us_tobacco_bar_plot,
       width = 30, height = 35)


##### Combine into a panel grid ##### 
all_res_bar_plot <- bind_rows(
  cn_ecigarettes_vs_cn_control |> mutate(dataset="China",      comparison="E-Cigarettes vs Controls"),
  cn_tobacco_vs_cn_control |> mutate(dataset="China",      comparison="Tobacco vs Controls"),
  cn_tobacco_vs_cn_ecigarettes |> mutate(dataset="China",      comparison="Tobacco vs E-Cigarettes"),
  
  us_ecigarettes_vs_us_control |> mutate(dataset="U.S.",      comparison="E-Cigarettes vs Controls"),
  us_tobacco_vs_us_control |> mutate(dataset="U.S.",      comparison="Tobacco vs Controls"),
  us_tobacco_vs_us_ecigarettes |> mutate(dataset="U.S.",      comparison="Tobacco vs E-Cigarettes"),
  
  cn_control_vs_us_control |> mutate(dataset="China vs. U.S.", comparison="Controls"),
  cn_ecigarettes_vs_us_ecigarettes |> mutate(dataset="China vs. U.S.", comparison="E-Cigarettes"),
  cn_tobacco_vs_us_tobacco |> mutate(dataset="China vs. U.S.", comparison="Tobacco")
) |> 
  mutate(
    dataset    = factor(dataset, levels = c("China", "U.S.", "China vs. U.S.")),
    comparison = factor(comparison,
                        levels = c("E-Cigarettes vs Controls",
                                   "Tobacco vs Controls",
                                   "Tobacco vs E-Cigarettes",
                                   "Controls", "E-Cigarettes", "Tobacco"))) 
# above: Lock facet order to 3x3 grid

# To get table of results
sigASVs_all <- all_res_bar_plot |> 
  filter(padj<0.01 & abs(log2FoldChange)>2) |>  
  dplyr::rename(ASV=row) 

# Get only asv names into a vector 
sigASVs_vec_all <- sigASVs_all |> 
  pull(ASV)

# Prune phyloseq file
all_merged_DESeq <- prune_taxa(sigASVs_vec_all,filtered_phyloseq)
sigASVs_all <- tax_table(all_merged_DESeq) |>  as.data.frame() |> 
  rownames_to_column(var="ASV")  |> 
  right_join(sigASVs_all, by = "ASV")  |> 
  arrange(log2FoldChange) |> # in Ascending order 
  mutate(Phylum = make.unique(Phylum)) |>  #only list unique Phylum to avoid repeating Genera
  mutate(Phylum = str_remove(Phylum, "^p__"))  |> 
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))



###### Facet Within US ONLY ######
p_us <- sigASVs_all |> 
  dplyr::filter(dataset == "U.S.") |> 
  ggplot(aes(x = Phylum, y = log2FoldChange)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x = Phylum,
                    ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE),
                width = 0.3) +
  coord_flip() +
  facet_grid(~ comparison, scales = "free_y") +
  ggtitle("Significant ASVs within U.S. Dataset (padj < 0.01, |log2FC| > 2)") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 5.5),
    axis.text.x = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    legend.title = element_text(face = "bold"),  
    panel.border = element_rect(linewidth = 2)
  )
p_us

ggsave("facet_US_bar_plot_asv.png",
       plot = p_us,
       width = 16, height = 12)


###### Facet Within CN only ######
p_cn <- sigASVs_all |> 
  dplyr::filter(dataset == "China") |> 
  ggplot(aes(x = Phylum, y = log2FoldChange)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE),
                width = 0.3) +
  coord_flip() +
  facet_grid(~ comparison, scales = "free_y") +
  ggtitle("Significant ASVs within Chinese Dataset (padj < 0.01, |log2FC| > 2)") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 5.5),
    axis.text.x = element_text(size = 9),
    strip.text  = element_text(face = "bold", size = 10),
    plot.title  = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(linewidth = 2, fill = NA)
  )
p_cn

ggsave("facet_CN_bar_plot_asv.png",
       plot = p_cn,
       width = 16, height = 14)


###### Facet Beweetn Datasets ######
p_between <- sigASVs_all |> 
  dplyr::filter(dataset == "China vs. U.S.") |> 
  ggplot(aes(x = Phylum, y = log2FoldChange)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE),
                width = 0.3) +
  coord_flip() +
  facet_grid(~ comparison, scales = "free_y") +
  ggtitle("Significant ASVs between China vs. U.S. (padj < 0.01, |log2FC| > 2)") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 5.5),
    axis.text.x = element_text(size = 9),
    strip.text  = element_text(face = "bold", size = 10),
    plot.title  = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(linewidth = 2, fill = NA)
  )
p_between

ggsave("facet_between_bar_plot_asv.png",
       plot = p_between,
       width = 25, height = 45)
