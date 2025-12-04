#!/usr/bin/env RScript 
# MICB 475 Relative Abundance analysis/ Abundance Modeling using DESeq2 RScript
# Create visualizations for Manuscript 
# 25 November, 2025

#### Load Packages ####
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(DESeq2))

#### Load Data (non-rarefied) ####
filtered_phyloseq <- get(load("filtered_phyloseq.RData")) # IMPORTANT: get() around load() gives the actual phyloseq object
print(filtered_phyloseq)

#### DESeq ####
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


#### Faceted Within-Group Volcano Plot Visualization ####
# 1) Bind 6 within-group DESeq2 result tables and tag them
volc_within <- bind_rows(
  us_ecigarettes_vs_us_control |> mutate(dataset="U.S.A.",      comparison="E-Cigarettes vs. Control"),
  us_tobacco_vs_us_control |> mutate(dataset="U.S.A.",      comparison="Tobacco vs. Control"),
  us_tobacco_vs_us_ecigarettes |> mutate(dataset="U.S.A.",      comparison="Tobacco vs. E-Cigarettes"),
  
  cn_ecigarettes_vs_cn_control |> mutate(dataset="China",      comparison="E-Cigarettes vs. Control"),
  cn_tobacco_vs_cn_control |> mutate(dataset="China",      comparison="Tobacco vs. Control"),
  cn_tobacco_vs_cn_ecigarettes |> mutate(dataset="China",      comparison="Tobacco vs. E-Cigarettes"),
) |> 
  mutate(
    dataset    = factor(dataset, levels = c("U.S.A.", "China")),
    comparison = factor(comparison,
                        levels = c("E-Cigarettes vs. Control",
                                   "Tobacco vs. Control",
                                   "Tobacco vs. E-Cigarettes"))) 

# 2) Faceted volcano (i.e., combined into a 2 row × 3 columns table)
facted_within_volcano_plot <- volc_within |> 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |>  
  ggplot(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) + 
  geom_point(alpha = 0.5, size = 1.95) +  
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +  
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  facet_grid(dataset ~ comparison) +
  scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) + 
  xlab("Log2FoldChange") +  
  ylab("-Log10(Adjusted P-Value)") +
  ggtitle("Within-U.S.A. and Within-China DESeq2 Analysis (padj < 0.01, |log2FC| > 2)") +
  theme_test() +  
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5), 
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(face = "bold",  size = 13),
    panel.border = element_rect(linewidth = 2), 
    strip.text = element_text(size = 12, face = "bold", colour = "black"),
    legend.position = "none")   # remove legend  
facted_within_volcano_plot

# Save plot 
# ggsave("manuscript_within_cohort_volcano.png",
#        plot = facted_within_volcano_plot,
#        width = 12, height = 8)

#### Supplementary Between-group Volcano Plot ####
# 1) Bind the 3 between-group DESeq2 result tables and tag them
volc_between <- bind_rows(
  cn_control_vs_us_control |> mutate(dataset="China vs. U.S.A.", comparison="Control"),
  cn_ecigarettes_vs_us_ecigarettes |> mutate(dataset="China vs. U.S.A.", comparison="E-Cigarettes"),
  cn_tobacco_vs_us_tobacco |> mutate(dataset="China vs. U.S.A.", comparison="Tobacco")
) |> 
  mutate(
    dataset    = factor(dataset, levels = c("China vs. U.S.A.")),
    comparison = factor(comparison,
                        levels = c("Control", "E-Cigarettes", "Tobacco")))

# 2) Generate supplementary between-group comparison (1 row × 3 columns) volcano plot 
facted_between_volcano_plot <- volc_between |> 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2)  |>  
  ggplot(aes(x = log2FoldChange, y = -log10(padj),  col=significant)) + 
  geom_point(alpha = 0.5, size = 1.95) +  
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +  
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  facet_grid(dataset ~ comparison) +
  scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) + 
  xlab("Log2FoldChange") +  
  ylab("-Log10(Adjusted P-Value)") +
  ggtitle("Between-Dataset (China vs U.S.A.) DESeq2 Analysis (padj < 0.01, |log2FC| > 2)") +
  theme_test() +  
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5), 
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(face = "bold",  size = 13),
    panel.border = element_rect(linewidth = 2), 
    strip.text = element_text(size = 12, face = "bold", colour = "black"),
    legend.position = "none")    
facted_between_volcano_plot

# Save plot 
# ggsave("supplementary_between_cohort_volcano.png",
#        plot = facted_between_volcano_plot,
#        width = 12, height = 8)



#### Bar graph based on the absolute “phylum count” table #### 
# 1. Load in Phylum count table 
phylum_counts <- read.csv("DESeq_abundvalues_cohort_status.csv")

# 2. Keep only the 4 comparisons of interest. Ignored E-Cig vs Tobacco for both cohorts
comp_levels <- c(
  "us_ecigarettes_vs_us_control",
  "us_tobacco_vs_us_control",
  "cn_ecigarettes_vs_cn_control",
  "cn_tobacco_vs_cn_control"
)

comp_labels <- c(
  "U.S.A. E-Cigarettes vs. Control",
  "U.S.A. Tobacco vs. Control",
  "China E-Cigarettes vs. Control",
  "China Tobacco vs. Control"
)

# 3. Convert the wide table into long format for plotting.
# wide to long
# turn columns "count_above_2" and "count_below_2"
# into two rows per phylum: one enriched, one depleted
phylum_long_all <- phylum_counts |>
  filter(Type %in% comp_levels) |>
  mutate(
    Type = factor(Type, levels = comp_levels, labels = comp_labels)
  ) |>
  pivot_longer(
    cols = c(count_below_2, count_above_2),
    names_to = "direction_raw",
    values_to = "Count"
  ) |>
  mutate(
    Count = as.numeric(Count),
    direction = case_when(
      direction_raw == "count_below_2" ~ "Depleted",
      direction_raw == "count_above_2" ~ "Enriched"
    ),
    direction = factor(direction, levels = c("Depleted", "Enriched"))
  )

# 4. Global phylum ordering 
# Alphabetical order per phylum across ALL conditions
global_totals <- phylum_long_all |>
  distinct(Phylum) |>
  arrange(Phylum)      # alphabetical A to Z

# use this order for ALL phyla
phylum_order <- global_totals$Phylum

phylum_long_ordered <- phylum_long_all |>
  mutate(
    Phylum = factor(Phylum, levels = phylum_order)  # all phyla, ordered alphabetically
  )

# 5. Faceted bar plot. Phylum on x-axis, Count on y-axis. 
# Two bars side-by-side per phylum: red  = depleted and blue = enriched
dodge_width <- 0.9
bar_width   <- 0.9

phylum_bar_facet <- ggplot(
  phylum_long_ordered,
  aes(x = Phylum, y = Count, fill = direction)
) +
  geom_col(
    position = position_dodge(width = dodge_width),
    width = bar_width
  ) +
  facet_wrap(~ Type, ncol = 1) + #vertical facets and free y-scales (scales = "free_y")
  scale_fill_manual(
    name   = "Direction",
    values = c(
      "Depleted" = "#C0392B",   # "#000000", black
      "Enriched" = "#2980B9"     # "#7F7F7F" medium gray
    )
  ) +
  xlab("Phylum") +
  ylab("Count") +
  ggtitle("Phylum-Level Enrichment and Depletion Across Tobacco and E-Cigarette Use in China and U.S.A. Cohorts") +
  theme_test() +
  theme(
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title   = element_text(face = "bold", size = 13),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 10),
    strip.text   = element_text(size = 12, face = "bold"),
    panel.border = element_rect(linewidth = 1.1),
    legend.title = element_text(face = "bold")
  )

phylum_bar_facet

# Save plot
# ggsave("barplot_phylum_count_fix_y_scale.png",
#        plot = phylum_bar_facet,
#        width = 12, height = 8)
