#### Aim 2: Identify the microbial taxa that are strongly associated with each category (Indicator Species Analysis) ####

#### Load in the 'tidyverse', 'phyloseq', and 'indicscpecies' packages ####
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load the filtered, non-rarefied and merged phyloseq object ####
load("../Diversity_Metrics_Analysis/filtered_phyloseq.RData")

#### Perform Indicator Species Analysis on each of the 6 categories ####
## The six categories are: us_control, us_tobacco, us_e-cigarettes, cn_control, cn_tobacco and cn_e-cigarettes

# Step 1: glom to Phylum level while keeping the NAs
merged_phyloseq_phylum <- tax_glom(merged_phyloseq_filtered, "Phylum", NArm = FALSE)
# Step 2: Convert the counts in the OTU table to relative abundance
merged_phyloseq_phylum_RA <- transform_sample_counts(merged_phyloseq_phylum, fun=function(x) x/sum(x))

# Step 3: Run ISA using the Combined column
isa <- multipatt(t(otu_table(merged_phyloseq_phylum_RA)), cluster = sample_data(merged_phyloseq_phylum_RA)$`Combined`)
summary(isa) #view results

# Step 4: Convert taxonomy to table with ASV IDs
taxtable <- tax_table(merged_phyloseq_filtered) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Step 5: Make final table by joining the indicator ASVs back to their taxonomy to interpret the results biologically.
final_table_phylum <- isa$sign %>%
  rownames_to_column(var="ASV") %>% #Convert the ASV from row names into a separate column
  left_join(taxtable) %>% #Join with taxa table based on the ASVs ID
  filter(p.value<0.05) #Filter for anything that has a p value of <0.05 

# Step 6: Save the table in the working directory
write.csv(final_table_phylum, file = "final_table_phylum.csv", row.names = FALSE)


