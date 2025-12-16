library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(vegan)
library(ggsignif)
library(dunn.test)
library(ggpubr)
library(RColorBrewer)
library(cowplot)

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


###### Alpha Diversity ###### 
##### Shannon Diversity #####
# Visualize Shannon Diversity results as box plots
boxplot_Shannon <- plot_richness(rarefied_phyloseq, x = "Combined", 
                                 measures = "Shannon", color = "Combined") +
                    geom_boxplot() +
                    scale_color_manual(values = c("#A63603","#D94801","#FD8D3C","#4A1486","#6A51A3","#807DBA")) +
                    scale_x_discrete(labels=c('China, Non-Smokers', 
                                              'China, E-cigarettes', 
                                              'China, Tobacco', 
                                              'USA, Non-Smokers',
                                              'USA, E-cigarettes',
                                              'USA, Tobacco')) +
                    theme_test(base_size = 18) +
                    theme(legend.position="none", 
                          axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
                          axis.title.x = element_blank())
boxplot_Shannon

ggsave(filename = "boxplot_Shannon.png",
       width = 5, height = 6, 
       boxplot_Shannon)

#### Statistical analysis for Shannon Diversity ####
# Perform Kruskall Wallis test on Shannon Diversity results
alphadiv <- estimate_richness(rarefied_phyloseq)
samp_dat <- sample_data(rarefied_phyloseq)
samp_dat_alpha <- data.frame(samp_dat, alphadiv)
View(samp_dat_alpha)

kruskal_Shannon <- kruskal.test(Shannon ~ `Combined`, data = samp_dat_alpha)
View(kruskal_Shannon) #Results: p-value = 8.912e-05 (suggests significant difference as p-value < 0.05)

# Perform Dunn's test on Shannon Diversity results
## Aim: Determine which comparisons have significant difference
dunn_Shannon <- dunn.test(samp_dat_alpha$Shannon, samp_dat_alpha$Combined, kw = TRUE)
# Found significance in
## cn_control vs cn_ecig = 0.0086
## cn_control vs us_control = 0.0002
## cn_tobacco vs us_control = 0.0029
## cn_control vs us_ecig = 0.0227
## cn_control vs us_tobacco = 0.0000
## cn_tobacco vs us_tobacco = 0.0002
## us_ecig vs us_tobacco = 0.0226


#### Mapping significance to Shannon Diversity box plots ####
# Re-plot the Shannon Diversity box plots to include significance
Shannon_signif_comparisons = list(c("cn_control","cn_ecigarettes"),
                                  c("cn_control","us_control"),
                                  c("cn_tobacco","us_control"),
                                  c("cn_control","us_ecigarettes"),
                                  c("cn_control","us_tobacco"),
                                  c("cn_tobacco","us_tobacco"),
                                  c("us_ecigarettes","us_tobacco"))
boxplot_Shannon_signif <- boxplot_Shannon +
                            theme_test(base_size = 18) +
                            theme(legend.position="none", 
                                  axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
                                  axis.title.x = element_blank()) +
                            geom_signif(comparisons = Shannon_signif_comparisons,
                                        annotations = c("0.0086","0.0002","0.0029","0.0227","0.0000","0.0002","0.0226"),
                                        step_increase = 0.1,
                                        textsize = 8) 
boxplot_Shannon_signif

ggsave("boxplot_Shannon_signif.PNG", boxplot_Shannon_signif, width = 10, height = 12)


##### Observed features #####
# Visualize Observed Features results as box plots
boxplot_Observed <- plot_richness(rarefied_phyloseq, x = "Combined", measures = "Observed", color = "Combined") +
                    geom_boxplot() +
                    scale_color_manual(values = c("#A63603","#D94801","#FD8D3C","#4A1486","#6A51A3","#807DBA")) +
                    scale_x_discrete(labels=c('China, Non-Smokers', 
                                              'China, E-cigarettes', 
                                              'China, Tobacco', 
                                              'USA, Non-Smokers',
                                              'USA, E-cigarettes',
                                              'USA, Tobacco')) +
                    theme_test(base_size = 18) +
                    theme(legend.position="none", 
                          axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
                          axis.title.x = element_blank())
boxplot_Observed 

# Save the boxplot
ggsave(filename = "boxplot_Observed.png", 
       boxplot_Observed)

#### Statistical analysis for Observed Features ####
# Perform Kruskall Wallis test on Observed Features results
alphadiv_obs <- estimate_richness(rarefied_phyloseq)
samp_dat_obs <- sample_data(rarefied_phyloseq)
samp_dat_alpha_obs <- data.frame(samp_dat_obs, alphadiv_obs)
View(samp_dat_alpha_obs)

kruskal_Observed <- kruskal.test(Observed ~ `Combined`, data = samp_dat_alpha_obs)
View(kruskal_Observed) #Results: p-value = 0.00113 (suggests significant difference as p-value < 0.05)

# Perform Dunn's test on Observed Features results
## Aim: Determine which comparisons have significant difference
dunn_Obs <-dunn.test(samp_dat_alpha_obs$Observed, samp_dat_alpha_obs$Combined, kw = TRUE)
# Found significance in
## cn_control vs cn_ecig = 0.0360 
## cn_control vs cn_tobacco = 0.0066
## cn_control vs us_control = 0.0004
## cn_control vs us_ecig = 0.0494
## cn_control vs us_tobacco = 0.0000
## cn_ecig vs us_tobacco = 0.0169
## cn_tobacco vs us_tobacco = 0.0149
## us_ecig vs us_tobacco = 0.0111

#### Mapping significance to Observed Features box plots ####
# Re-plot the Observed Features box plots to include significance
Observed_signif_comparisons = list(c("cn_control","cn_ecigarettes"),
                                   c("cn_control","cn_tobacco"),
                                   c("cn_control","us_control"),
                                   c("cn_control","us_ecigarettes"),
                                   c("cn_control","us_tobacco"),
                                   c("cn_ecigarettes","us_tobacco"),
                                   c("cn_tobacco","us_tobacco"),
                                   c("us_ecigarettes","us_tobacco"))

boxplot_Observed_signif <- boxplot_Observed +
                            theme_test(base_size = 18) +
                            theme(legend.position="none", 
                                  axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
                                  axis.title.x = element_blank()) +
                            geom_signif(comparisons = Observed_signif_comparisons,
                                        annotations = c("0.0360","0.0066","0.0004","0.0494","0.0000",
                                                        "0.0169","0.0149","0.0111"),
                                        step_increase = 0.1,
                                        textsize = 8) 
boxplot_Observed_signif

ggsave("boxplot_Observed_signif.PNG", boxplot_Observed_signif, width = 10, height = 12)

##### Combined Alpha Diversity plots #####
# Combine Shannon Diversity and Observed Features box plots
alphadiv_plots <- ggarrange(boxplot_Shannon_signif + theme(axis.title.x = element_blank()), 
                            boxplot_Observed_signif + theme(axis.title.y = element_blank(),
                                                            axis.title.x = element_blank()), 
                            common.legend = TRUE, legend = "none")
alphadiv_plots

# Save the alpha diversity plot
ggsave("alpha_diversity_boxplots.PNG", alphadiv_plots, width = 14, height = 12)


###### Beta Diversity ###### 
##### Bray Curtis #####
# Visualize Bray Curtis results as PCoA plot
dm_bc <- distance(rarefied_phyloseq, method= "bray")
pcoa_bc <- ordinate(rarefied_phyloseq, method= "PCoA", distance= dm_bc)

pcoa_plot_bc <- plot_ordination(rarefied_phyloseq, pcoa_bc, color = "Combined") +
                labs(col = "Cohorts \n(Country, Smoking Status)")
pcoa_plot_bc

ggsave("PCOA_plot_BrayCurtis.png", 
       pcoa_plot_bc)

#### Statistical analysis for Bray Curtis ####
# Perform PERMANOVA on Bray Curtis results
samp_dat_wdiv <- data.frame(sample_data(rarefied_phyloseq), estimate_richness(rarefied_phyloseq))

dm_bray <- vegdist(t(otu_table(rarefied_phyloseq)), method="bray")
adonis2(dm_bray ~ `Combined`, data = samp_dat_wdiv)
adonis2(dm_bray ~ `Description`*Cohort, data = samp_dat_wdiv)

#### Mapping significance to Bray Curtis PCoA plot ####
# Re-plot the Bray Curtis PCoA plot to include ellipses that show significant difference
ord.bray <- ordinate(rarefied_phyloseq, method="PCoA", distance="bray")
pcoa_bc_signif <- plot_ordination(rarefied_phyloseq, ord.bray, color =  "Description", shape = "Cohort") +
                  stat_ellipse(type = "norm") +
                  labs(title = "Bray Curtis", col = "Smoking Status", shape = "Country") +
                  scale_shape_manual(labels = c("China", "U.S.A."), values = c(17,15)) +
                  scale_color_manual (values = c("#000000","#56B4E9", "#CC79A7"), 
                                      labels = c("Control", "E-cigarettes", "Tobacco")) +
                  theme_classic(base_size = 18) +
                  theme(plot.title = element_text(hjust = 0.5))
pcoa_bc_signif

ggsave("PCOA_BrayCurtis_signif.png",
       pcoa_bc_signif, width = 8, height = 6, units = "in")

##### Jaccard #####
# Visualize Jaccard results as PCoA plot
dm_jac <- distance(rarefied_phyloseq, method= "jaccard")
pcoa_jac <- ordinate(rarefied_phyloseq, method= "PCoA", distance= dm_jac)

pcoa_plot_jac <- plot_ordination(rarefied_phyloseq, pcoa_jac, color = "Combined") +
                  labs(col = "Cohorts \n(Country, Smoking Status)")
pcoa_plot_jac

ggsave("PCOA_plot_Jaccard.png", 
       pcoa_plot_jac)

#### Statistical analysis for Bray Curtis ####
# Perform PERMANOVA on Bray Curtis results
samp_dat_wdiv_jac <- data.frame(sample_data(rarefied_phyloseq), estimate_richness(rarefied_phyloseq))

dm_jaccard <- vegdist(t(otu_table(rarefied_phyloseq)), method="jaccard")
adonis2(dm_jaccard ~ `Combined`, data= samp_dat_wdiv_jac)

#### Mapping significance to Jaccard PCoA plot ####
# Re-plot the Jaccard PCoA plot to include ellipses that show significant difference
ord.jaccard <- ordinate(rarefied_phyloseq, method="PCoA", distance="jaccard")
pcoa_jac_signif <- plot_ordination(rarefied_phyloseq, ord.jaccard, color =  "Description", shape = "Cohort") +
                    stat_ellipse(type = "norm") +
                    labs(title = "Jaccard", col = "Smoking Status", shape = "Country") +
                    scale_shape_manual(labels = c("China", "U.S.A."), values = c(17,15)) +
                    scale_color_manual (values = c("#000000","#56B4E9", "#CC79A7"), 
                                        labels = c("Control", "E-cigarettes", "Tobacco")) +
                    theme_classic(base_size = 18) +
                    theme(plot.title = element_text(hjust = 0.5))
pcoa_jac_signif
ggsave("PCOA_Jaccard_signif.png",
       pcoa_jac_signif, width = 8, height = 6, units = "in")

##### Combined Beta Diversity plots #####
# Combine Bray Curtis and Jaccard PCoA plots
betadiv_plots <- ggarrange(pcoa_bc_signif,
                            pcoa_jac_signif, 
                            common.legend = TRUE, legend = "right")
betadiv_plots

# Save the beta diversity plot
ggsave("beta_diversity_plots.PNG", betadiv_plots, width = 12, height = 6)