#### Load required library packages ####
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(ggpubr)

#### Load data ####
## Load the merged phyloseq object that have been filtered, but not rarefied
load("filtered_phyloseq.RData")
# metadata <- sample_data(merged_phyloseq_filtered)

#### Core microbiome ####
## Convert to relative abundance
phyloseq_RA <- transform_sample_counts(merged_phyloseq_filtered, fun=function(x) x/sum(x))

## Subset dataset into each cohort
us_non_smokers <- subset_samples(phyloseq_RA, `Combined` == "us_control")
us_ecig <- subset_samples(phyloseq_RA, `Combined` == "us_ecigarettes")
us_tobacco <- subset_samples(phyloseq_RA, `Combined` == "us_tobacco")

cn_non_smokers <- subset_samples(phyloseq_RA, `Combined` == "cn_control")
cn_ecig <- subset_samples(phyloseq_RA, `Combined` == "cn_ecigarettes")
cn_tobacco <- subset_samples(phyloseq_RA, `Combined` == "cn_tobacco")

## Identify core members using an abundance threshold of 0.001 and prevalence threshold of 0.5 
us_non_smokers_core <- core_members(us_non_smokers, detection=0.001, prevalence = 0.5)
us_ecig_core <- core_members(us_ecig, detection=0.001, prevalence = 0.5)
us_tobacco_core <- core_members(us_tobacco, detection=0.001, prevalence = 0.5)

cn_non_smokers_core <- core_members(cn_non_smokers, detection=0.001, prevalence = 0.5)
cn_ecig_core <- core_members(cn_ecig, detection=0.001, prevalence = 0.5)
cn_tobacco_core <- core_members(cn_tobacco, detection=0.001, prevalence = 0.5)

## Visualize core members as a Venn diagram
# Venn diagram for all cohorts
venn_all <- ggVennDiagram(x=list(US_Non_Smokers = us_non_smokers_core, 
                                      US_Ecigarettes = us_ecig_core,
                                      US_Tobacco = us_tobacco_core,
                                      CN_Non_Smokers = cn_non_smokers_core, 
                                      CN_Ecigarettes = cn_ecig_core,
                                      CN_Tobacco = cn_tobacco_core), 
                          label_color = "black",
                          label_alpha=0,
                          set_color = c("red","brown","green","blue","purple","magenta" )) + 
            scale_x_continuous(expand = expansion(mult = .2)) +
            scale_fill_distiller(palette = "Paired") +
            labs(fill = "Number of ASVs") +
            theme(legend.position="right")
venn_all

ggsave("venn_all.PNG", 
       venn_all, 
       bg = "white",
       width = 15, 
       height = 10, 
       units = "in")

## Venn diagrams grouped by geographical location
all_cohorts <- list(US_NS = us_non_smokers_core, 
                    US_EC = us_ecig_core, 
                    US_TOB = us_tobacco_core, 
                    CN_NS = cn_non_smokers_core, 
                    CN_EC = cn_ecig_core, 
                    CN_TOB = cn_tobacco_core)

# Venn diagram for U.S.A. cohorts
US_venn <-ggVennDiagram(x=all_cohorts[1:3], category.names = c("Non-Smokers", "E-cigarettes", "Tobacco"),
                        label_color = "white", label_size = 10, label_alpha = 0.2, set_size = 10) + 
          scale_x_continuous(expand = expansion(mult = .2)) +
          scale_fill_distiller(palette = "Purples", direction = 1, limits = c(2, 22)) +
          labs(fill = "Number of \nASVs \n", title = "U.S.A.") +
          theme(legend.position="right", legend.text = element_text(size = 25), legend.title = element_text(size = 30), 
                title = element_text(size = 30), plot.title = element_text(hjust = 0.5))
US_venn

ggsave("US_venn.PNG", 
       US_venn, 
       bg = "white",
       width = 12, 
       height = 8, 
       units = "in")

# Venn diagram for China cohorts
CN_venn <- ggVennDiagram(x=all_cohorts[4:6], category.names = c("Non-Smokers", "E-cigarettes", "Tobacco"),
                         label_color = "white", label_size = 10, label_alpha = 0.2, set_size = 10) + 
            scale_x_continuous(expand = expansion(mult = .2)) +
            scale_fill_distiller(palette = "Oranges", direction = 1, limits = c(2, 22)) +
            labs(fill = "Number of \nASVs \n", title = "China") +
            theme(legend.position="right", legend.text = element_text(size = 25), legend.title = element_text(size = 30), 
                  title = element_text(size = 30), plot.title = element_text(hjust = 0.5))
CN_venn

ggsave("CN_venn.PNG", 
       CN_venn, 
       bg = "white",
       width = 12, 
       height = 8, 
       units = "in")

# Combine the U.S.A. and China venn diagrams
venn_us_cn <- ggarrange(US_venn,
                        CN_venn, 
                        legend = "right")
venn_us_cn 
ggsave("venn_us_cn.PNG", venn_us_cn, width = 22, height = 12, bg = "white")