library(tidyverse)

#install "realxl" package to be able to import excel file 
install.packages("readxl")   # only once
library(readxl)


##Editing metadata tables so that they each contain a column called "Description"
#wher eis determines if the subject consumes only "Ecigarretes", only "Tobacco", 
#or "None", any other possibility is NA 


## Editing the vaping_uk_meatadata

meta_uk_fp <- "vaping_uk_metadata.xlsx"
meta_uk <- read_excel(meta_uk_fp)

## Creating a new column, looking at the "Ecig" and "Tobacco" conditons
meta_uk_new <- meta_uk %>%
  mutate(Description = case_when(
    Ecig == "No" & Tobacco == "No" ~ "None",
    Ecig == "Yes" & Tobacco == "No" ~ "Ecigarettes",
    Ecig == "No" & Tobacco == "Yes" ~ "Tobacco",
    TRUE ~ NA_character_  
  ))

## Editing the vaping_chinese_meatadata

meta_cn_fp <- "vaping_chinese_metadata.xlsx"
meta_cn <- read_excel(meta_cn_fp)

## Creating a new column, whith the same contidions as above
meta_cn_new <- meta_cn  %>%
  mutate(Description = case_when(
    str_detect(`Public description`, regex("E-cig", ignore_case = TRUE)) ~ "Ecigarettes",
    str_detect(`Public description`, regex("Non", ignore_case = TRUE)) ~ "None",
    str_detect(`Public description`, regex("Common", ignore_case = TRUE)) ~ "Tobacco",
    TRUE ~ NA_character_
  ))


## Save both files as cn_metadata_new and uk_metadata_new

write_tsv(meta_cn_new, "cn_metadata_new.tsv")

write_tsv(meta_uk_new, "uk_metadata_new.tsv")

