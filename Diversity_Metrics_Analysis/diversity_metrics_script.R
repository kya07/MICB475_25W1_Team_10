library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(vegan)
library(ggsignif)


##### Load in RData #####
load("merged_phyloseq.RData")


##### Filter Merged Phyloseq Object ####
# 1. Remove non-bacterial sequences (After filtering: 123 samples)
phyloseq_bacteria <- merged_phyloseq %>% 
                      subset_taxa(Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# 2. Include only saliva samples (After filtering: 63 samples)
phyloseq_saliva <- phyloseq_bacteria %>% 
                    subset_samples(SampleType == "Saliva" | Isolation.source == "saliva")

# 3. Remove NA from Description column (After filtering: 54 samples)
phyloseq_no_NA <- phyloseq_saliva %>% 
                    subset_samples(!is.na(Description))

# 4. Remove marijuana use (After filtering: 50 samples)
phyloseq_no_marijuana <- phyloseq_no_NA %>% 
                          subset_samples(Marijuana == "No" | is.na(Marijuana))

# 5. Remove BMI < 18.5 (After filtering: 48 samples)
phyloseq_BMI <- phyloseq_no_marijuana %>% 
                  subset_samples(BMI>=18.5 | is.na(BMI))

# 6. Include e-cigarette smokers whose age fall within 24-48 years old, 
#    tobacco smokers whose age fall within 28-47 years old,
#    non-smokers whose age fall within 29-54 years old (After filtering: 44 samples)
phyloseq_age <- phyloseq_BMI %>% 
                subset_samples((Description == "ecigarettes" & Age >= 24 & Age <= 48) | 
                               (Description == "tobacco" & Age >= 28 & Age <= 47) | 
                               (Description == "control" & Age >= 29 & Age <= 54) |
                               is.na(Age))
# Save filtered phyloseq object
merged_phyloseq_filtered <- phyloseq_age
save(merged_phyloseq_filtered, 
     file="filtered_phyloseq.RData")


#### Generating Rarefaction Curve #### 
# Before filtering the merged phyloseq object
rarecurve(t(as.data.frame(otu_table(merged_phyloseq))), cex=0.1)

# After filtering the merged phyloseq object
rarecurve(t(as.data.frame(otu_table(merged_phyloseq_filtered))), cex=0.1)


#### Rarefy Samples #### 
rarefied_phyloseq <- rarefy_even_depth(merged_phyloseq_filtered, rngseed = 1, sample.size = 10870)

# Save rarefied phyloseq object
save(rarefied_phyloseq, 
     file="rarefied_phyloseq.RData")


##### Diversity metrics ##### 
#### Shannon Diversity ####
## Generate Shannon diversity boxplot
boxplot_Shannon <- plot_richness(rarefied_phyloseq, x = "Combined", measures = "Shannon") +
                xlab("Cohorts \n(Geographic Location_Smoking Status)") +
                geom_boxplot()
boxplot_Shannon

# Save the boxplot
ggsave(filename = "boxplot_Shannon.png", 
       boxplot_Shannon)

#### Kruskal Wallis Test ####
alphadiv <- estimate_richness(rarefied_phyloseq)
samp_dat <- sample_data(rarefied_phyloseq)
samp_dat_alpha <- data.frame(samp_dat, alphadiv)
View(samp_dat_alpha)

kruskal_Shannon <- kruskal.test(Shannon ~ `Combined`, data = samp_dat_alpha)
View(kruskal_Shannon) #Results: p-value = 8.912e-05 (suggests significant difference as p-value < 0.05)

## Finding out which comparisons have significant difference 
lm_shannon_vs_combined_log <- lm(log(Shannon) ~ `Combined`, data=samp_dat_alpha)
anova_shannon_vs_combined_log <- aov(lm_shannon_vs_combined_log)
summary(anova_shannon_vs_combined_log)
TukeyHSD(anova_shannon_vs_combined_log)

## Mapping the significance to the box plot [INCOMPLETE]
# ggplot(samp_dat_alpha, aes(x=`Combined`, y=Shannon)) +
#        geom_boxplot() +
#        geom_signif()


#### Bray Curtis #### 
dm_bc <- distance(rarefied_phyloseq, method= "bray")

pcoa_bc <- ordinate(rarefied_phyloseq, method= "PCoA", distance= dm_bc)

pcoa_plot_bc <- plot_ordination(rarefied_phyloseq, pcoa_bc, color = "Combined") +
                labs(col = "Cohorts \n(Geographic Location_Smoking Status)")
pcoa_plot_bc

ggsave("PCOA_plot_BrayCurtis.png", 
       pcoa_plot_bc)

#### PERMANOVA ####
samp_dat_wdiv <- data.frame(sample_data(rarefied_phyloseq), estimate_richness(rarefied_phyloseq))

dm_bray <- vegdist(t(otu_table(rarefied_phyloseq)), method="bray")
adonis2(dm_bray ~ `Combined`, data=samp_dat_wdiv)


## Re-plot the above PCoA with ellipses to show a significant difference between cohorts
ord.bray <- ordinate(rarefied_phyloseq, method="PCoA", distance="bray")
pcoa_plot_bc_annotated <- plot_ordination(rarefied_phyloseq, ord.bray, color = "Combined") +
                          stat_ellipse(type = "norm") +
                          labs(col = "Cohorts \n(Geographic Location_Smoking Status)")
ggsave("PCOA_plot_BrayCurtis_annotated.png", 
       pcoa_plot_bc_annotated)