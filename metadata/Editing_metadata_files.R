library(tidyverse)

#install "realxl" package to be able to import excel file 
install.packages("readxl")   # only once
library(readxl)


##Editing metadata tables so that they each contain a column called "Description"
#wher eis determines if the subject consumes only "Ecigarretes", only "Tobacco", 
#or "None", any other possibility is NA 

##Oct 16th -> adding two new columns, one for cohort and another one for 'combined'
# which will be 'cohort_description'




#### Editing the vaping_uk_meatadata ####

meta_us_fp <- "metadata/vaping_uk_metadata.xlsx"
meta_us <- read_excel(meta_us_fp)

## Creating a new column, looking at the "ecig" and "tobacco" conditons
meta_us_description_column <- meta_us %>%
  mutate(Description = case_when(
    Ecig == "No" & Tobacco == "No" ~ "control",
    Ecig == "Yes" & Tobacco == "No" ~ "ecigarettes",
    Ecig == "No" & Tobacco == "Yes" ~ "tobacco",
    TRUE ~ NA_character_  
  ))

## Creating a new column, for cohort--> all will be us 
meta_us_description_cohort_columns <- meta_us_description_column %>%
  mutate(Cohort = "us")

## Creating a new column combining the two new columns
meta_us_description_cohort_combined_columns <- meta_us_description_cohort_columns %>%
  mutate (Combined = paste(Cohort, Description, sep = "_"))

##Renaming the final metadata file
meta_us_new <- meta_us_description_cohort_combined_columns
  
  


#### Editing the vaping_chinese_meatadata ####

meta_cn_fp <- "metadata/vaping_chinese_metadata.xlsx"
meta_cn <- read_excel(meta_cn_fp)

## Creating a new column, whith the same contidions as above
meta_cn_description_column <- meta_cn  %>%
  mutate(Description = case_when(
    str_detect(`Public description`, regex("E-cig", ignore_case = TRUE)) ~ "ecigarettes",
    str_detect(`Public description`, regex("Non", ignore_case = TRUE)) ~ "control",
    str_detect(`Public description`, regex("Common", ignore_case = TRUE)) ~ "tobacco",
    TRUE ~ NA_character_
  ))

## Creating a new column, for cohort--> all will be us 
meta_cn_description_cohort_columns <- meta_cn_description_column %>%
  mutate(Cohort = "cn")

## Creating a new column combining the two new columns
meta_cn_description_cohort_combined_columns <- meta_cn_description_cohort_columns %>%
  mutate (Combined = paste(Cohort, Description, sep = "_"))

##Renaming the final metadata file
meta_cn_new <- meta_cn_description_cohort_combined_columns

## Save both files as cn_metadata_new and us_metadata_new (the "uk" is actually from us volunteers)

write_tsv(meta_cn_new, "metadata/cn_metadata_new.tsv")

write_tsv(meta_us_new, "metadata/us_metadata_new.tsv")

