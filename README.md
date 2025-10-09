# MICB475_25W1_Team_10
This repository contains weekly **meeting agendas and minutes** for **Team 10 (Aimee, Kaiya, Kayla and Sophie)** in UJEMI Project 2 (MICB475).
Agendas are posted before each meeting for mentor review, and minutes are added to record key discussions, decisions, and next steps. Team members will alternate responsibility for maintaining these records to ensure consistent tracking of project progress.

## Weekly Records
### Week 2 (09-Oct-2025)
#### Agenda
* Proposal is due on the 26th
* We have to choose

#### Meeting Minutes
 Ideas (✅):
- **Smoking and vaping w sample locations** dataset (#18): Ethnic differences in the USA between smoking (=vaping) populations (Caucasian vs. non-Caucasian) 
  - **What should we look at:**
      - We are looking at the American one
      - Multiple oral samples were taken, and it was said that there were significant differences between them, but idk....
      - May not be a good idea
  - **Backup (✅):**
      - Comparing Chinese and US-based (analysis done in the UK) cohorts --> only e-cigs 
      - UK but US-based --> saliva, fecal and oral
      - China --> saliva
      - Three groups --> only e-cigs, only tobacco, and none(control)
      - Have to merge, first filter all the UK metadata, and then make a new column with the description of the three groups 
      - We can do all the taxonomic analyses possible
      - We are only exploring one column, so we gotta do all the possible analyses later on
      - Merging
          - We have to check if they got the same variable regions
            - if they did --> talk to Ritu on how to merge it from the very beginning
            -  If they didn't --> there is an easy 5-line code to merge at the end 
          - have to merge the etadata too 
          - Have an origin, a sample ID, and a group version
       
  - For the proposal
      - We should try to do all the merging before the deadline (they are flexible, but it's  best to do it) --> ideally get to the denoising step 

  - How do we do it? Prepare the data to merge!
      - Variables are different, so we will merge at the end before we go
      - Start doing the Qiime pipeline up to the rarefaction for both
          - one person per metadata (2 people)
          - 2 people working on the metadata
          - Do not delete anything from the metadata
              - Create a new column in both metadata tables
                   - for the ones that don't fall into either of those, put NA 
              - Do everything separately until we build the phyloseq object before putting it to R
              - Do all filtering stuff after merging
           
      - Every file that we modify, keep a copy and work on it --> make copies for everything
      - check with Ritu before we trim anything
      - To make the new metadata 
 
#### Action items
- Have the new column in both metadata files done and ready--> they have to be exactly the same 
- Read the proposal docs to know what we have to do 
- Make action items for next week
- Start the pipeline!!!
- We have to upload all of our files on GitHub --> include all qza files too

### Week 1 (02-Oct-2025)
#### Agenda
* Potential datasets of interest: #11 (Anemia), #14 (Diabetes), #18 (Vaping and Smoking), #24 (Gastric cancer)
* Present possible research questions and discuss with TA
* Ensure all members have access to the team’s GitHub repository<br/><br/> 

#### Meeting Minutes
Approved Ideas (✅):
- **Smoking and vaping** dataset (#18): Ethnic differences in the USA between smoking (=vaping) populations (Caucasian vs. non-Caucasian) 
  - **Main factors of interest (3):**  
      1. Sample location  
      2. Ethnicity  
      3. Smoker vs. non-smoker  
   - **Proposed workflow:**
     - Categorize by sample location and compare between locations
     - Within each location, analyze ethnicity
     - Conduct functional analysis (i.e., gender comparisons) only *after narrowing locations where clear differences are observed* <br/><br/>
    
Not Approved / Rejected Ideas (❌):
- **Anemia** dataset (#11):
  - gender comparison 
  - diet comparison  
  - comparison of different months 

- **Smoking and vaping** dataset (#18): 
  - comparison between geographical locations, but some with a limited sample size and fewer annotations
  - Different sample sources (fecal vs. oral) 
  - Population study (cohort vs. cohort comparison) by combining all datasets, focusing only on e-cigarette ✅ (backup plan)
  - compare ethnicity in the UK (limited representation) 

- **Diabetes** dataset (#14):
  - both datasets are Texas-based. The "Mexican" dataset was conducted in Mexico. 
  - BMI and diabetes 
  - Look more in depth into the dataset: Does diabetes influence bacterial pathway enrichment versus depletion?
  - Functional analyses (e.g., comparison between genders)
    
- **Gastric cancer** dataset (#24):
  - Most significant difference is between gastric cancer vs. healthy
  - Comparison of different subtype histopathologies   
  
- **Alcohol consumption** database (#22):
  - alcohol consumption vs BMI<br/><br/>    
 
#### Action items
- add TA Ritu onto the repository and the team #
- explore the USA_Smoking_Datasets
- look into literature to search for additional smoking/microbiome datasets with additional ethnicity (by tomorrow 03-Oct-24) 



