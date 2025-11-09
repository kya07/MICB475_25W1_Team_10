#### Aim 2: Identify the microbial taxa that are strongly associated with each category (Indicator Species Analysis) ####

#### Load in the 'tidyverse', 'phyloseq', and 'indicscpecies' packages ####
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load the filtered, non-rarefied and merged phyloseq object ####
load("../Diversity_Metrics_Analysis/filtered_phyloseq.RData")

#### Perform Indicator Species Analysis on each of the 9 groups selected####
## The nine groups selected were: us_control v.s us_tobacco, us_control v.s us_e-cigarette, cn_control v.s cn _tobacco, 
## cn_control v.s cn_e-cigarette, us_control v.s cn_control, us_tobacco v.s cn_tobacco, us_e-cigarette v.s cn_e-cigarette, 
## us_e-cigarette v.s us_tobacco, cn_e-cigarette v.s cn_tobacco ##

# Step 1: glom to Genus level while keeping the NAs
merged_phyloseq_genus <- tax_glom(merged_phyloseq_filtered, "Genus", NArm = FALSE)
# Step 2: Convert the counts in the OTU table to relative abundance
merged_phyloseq_genus_RA <- transform_sample_counts(merged_phyloseq_genus, fun=function(x) x/sum(x))

# Step 3: Make a list of our comparison groups
comparisons <- list(
  c("us_control","us_tobacco"),
  c("us_control","us_ecigarettes"),
  c("cn_control","cn_tobacco"),
  c("cn_control","cn_ecigarettes"),
  c("us_control","cn_control"),
  c("us_tobacco","cn_tobacco"),
  c("us_ecigarettes","cn_ecigarettes"),
  c("us_ecigarettes","us_tobacco"),
  c("cn_ecigarettes","cn_tobacco")
)

# Step 4: Run Automated ISA for All Pairwise Comparisons
results <- list()   # initialize the list to store output

for(i in seq_along(comparisons)){
  
  pair <- comparisons[[i]]
  
  # subset to the two groups (use Combined column)
  pair_ps <- subset_samples(merged_phyloseq_genus_RA, Combined %in% pair)
  
  # extract OTU table
  otu <- as.data.frame(t(otu_table(pair_ps)))
  # correct grouping variable (MUST match subset)
  grp <- sample_data(pair_ps)$Combined
  
  # Run ISA
  res <- multipatt(otu, grp, func = "IndVal.g", control = how(nperm = 999))
  
  # store summary in results list
  results[[paste(pair, collapse = "_vs_")]] <- res
}

# View results
names(results)


# Step 5: Convert taxonomy to table with ASV IDs
taxtable <- tax_table(merged_phyloseq_filtered) %>% as.data.frame() %>% rownames_to_column(var="ASV") #convert the ASV from row names into a separate column

# Step 6: Combine all res$sign tables from the list
all_results <- map_df(names(results), function(name){
  results[[name]]$sign %>%
    as.data.frame() %>%
    rownames_to_column("ASV") %>% #convert the ASV from row names into a separate column
    mutate(Comparison = name)
})

#Step 7: Make a final table by joining the indicator ASVs back to their taxonomy to interpret the results biologically
final_table <- all_results %>%
  left_join(taxtable, by = "ASV") %>% #join with taxa table based on the ASVs ID 
  filter(p.value < 0.05) ##Filter for anything that has a p value of <0.05 

#Reorder the columns in the existing final_table for better visualization
final_table <- final_table %>%
  relocate(s.us_ecigarettes,s.cn_control, s.cn_tobacco, s.cn_ecigarettes,
           .before = index)

# View final result
View(final_table)

