# MICB475_25W1_Team_10
This repository contains weekly **meeting agendas and minutes** for **Team 10 (Aimee, Kaiya, Kayla and Sophie)** in UJEMI Project 2 (MICB475).
Agendas are posted before each meeting for mentor review, and minutes are added to record key discussions, decisions, and next steps. Team members will alternate responsibility for maintaining these records to ensure consistent tracking of project progress.

## Weekly Records
### Week 4 (23-Oct-2025)
#### Agenda
Questions about the Project workflow/execution: 
* Discuss what the ideal sampling depth range is for each dataset.
* Need more clarification and directions on metadata processing. Merging of metadata? Merging of phyloseq objects? Filter before merging or the other way around?
* Regarding data processing:
  * Is alpha rarefaction curve created differently for paired-end sequences?
  * How to merge phyloseq objects?
  * Should we merge the phyloseq objects before determining the sampling depth for alpha diversity?
  * What variables should we control for?
* Regarding diversity metrics:
  * What alpha and beta diversity metrics should we use? How to do this in R?
  * Do we need to subset the combined dataset before performing diversity metrics?
  * Should we use Kruskall-Wallis test for alpha diversity and PERMANOVA for beta diversity?
* Regarding Indicator Taxa Analysis:
  * Should we use the taxonomic rank "genus" as the basis for grouping?
* Regarding core microbiome:
  * How should we specify the detection abundance and prevalence thresholds?
  * Should we create both venn diagram and bar plot or just one of them?
* Regarding functional analysis:
  * Should we create a heatmap, volcano plot, and bar plot, or is one of them enough?
  * Do we need to perform statistical analysis?
  
Questions about the Proposal: 
* After citing qimme2 article, do we need to cite DADA2 denoising separately?
* Confirm "Aims and rationale" and "Proposed approach" sections

---
### Week 3 (16-Oct-2025)
#### Agenda
* Q: How do we choose trimming/truncation depth for paired-end reads in QIIME2? Do the forward and reverse reads need to be trimmed to the same length?

#### Meeting Minutes
<ins> Overview of the progress that has been made </ins>
* **Metadata file**
  - Made a code on R to add a new ‘Description’ column to both the UK and China metadata files
  - For the UK metadata file
    * If the **‘Ecig’** column is **Yes** and the **‘Tobacco’** column is **No**, the **‘Description’** column says **ecigarretes**
    * If the **‘Ecig’** column is **No** and the **‘Tobacco’** column is **Yes**, the **‘Description’** column says **tobacco**
    * If both the **‘Ecig’** column and the **‘Tobacco’** column is **No**, the **‘Description’** column says **none**
    * If both the **‘Ecig’** column and the **‘Tobacco’** column is **Yes**, the **‘Description’** column says **NA**
    * Saved the file as **“us_metadata_new.tsv”**
  - **For the China metadata file**
    * It has a column called **‘Public description’** which states if they do **“E-cigarretes smoking #“**, **“Common tobacco smoking #”**, **“Non-smoking #”**, or **”Quitting smoking with tobacco #”** where **# is a number**, so each one of the entries is different
    * If the **‘Public Description’** column starts with **“E-cig”**, the **‘Description’** column says **ecigarretes**
    * If the **‘Public Description’** column starts with **“Common”**, the **‘Description’** column says **tobacco**
    * If the **‘Public Description’** column starts with **“Non”**, the **‘Description’** column says **none**
    * If the **‘Public Description’** column starts with **“Quitting”**, the **‘Description’** column says **NA**
    * Saved the file as **“cn_metadata_new.tsv”**
      
* **UK Dataset**
  - Have performed the Qiime2 pipeline until taxonomy analysis and generated a taxa-bar-plots.qzv and tree files
    * Found out that the dataset was paired ends
    * When performing manifest, received an error message that was resolved by renaming all the  uk_seqs and uk_manifest.tsv files from “.fastq” to “.fast.gz”
    * Performed the taxonomic analysis using the “silva-138-99-515-806-nb-classifier.qza” classifier since the UK dataset uses the V4 region
    * See details on the Code script qiime_pipeline_uk file
   - All the files are stored in the qiime_pipeline_uk folder
     
* **China Dataset**
  - In progress

<ins> Q&A </ins>
1. For paired ends sequences, do we need to trim both the reverse and forward sequences at the same length?
   - No, it doesn't really matter
   - During discussion, Dr. Sun was added to the GitHub and after taking a look at the quality graphs from the demux files, she recommended to trim both forward and reverse at 200. However, since the qiime2 pipeline took hours to run, she looked at the table.qzv file and found that there were minimal sample loss when the trimming parameters (243 for the forward and 221 for the reverse) were used. Therefore, no changes were need to be made.
2. How to directly move files into GitHub?
   - There is no direct way. To move files into GitHub, export the files from the server into local computer, then upload the files into GitHub
     
<ins> Discussion on the proposal and the next steps </ins>
* Before writing the proposal, we need to reach the step right after denoising in the Qiime2 pipeline for both of the datasets
* For experimental aims, the first aim should be filtering and wrangling data to prepare for analysis
* The steps below are considered as experimental approaches
1. ***Edit the metadata***
   * In the 'Description' column of both datasets, change "none" to "control"
   * Add a column called "cohort" which specify where the dataset is from, i.e. UK or China
   * Add a column called "combined" which sepcify where the dataset is from and the main category which the data belongs to, e.g. UK_tobacco, UK_ecig, UK_control, etc.
2. ***Merge the UK and China phyloseq objects***
   * Ask Ritu how to merge the phyloseq objects
3. ***Filter the merged phyloseq object***
   * Filter out non-bacteria sequences, e.g. mitochondria, eukarya, chloroplast, etc.
   * Filter out the NAs
   * Control the confounding variables
     - Keep in mind that
       - Trying to control all the confounding variables (i.e. more filtering) will result in fewer samples being retained
       - Controlling some of the confounding variables (i.e. less filtering) will result in more samples being retained, but with limitations
       - There are some variables that just can't be controlled
     - Recommendations from Dr. Sun
       - Filter in such a way that both the datasets are comparable
       - Need to do more filtering in the UK dataset compared to the China dateset due to the bigger size and more metadata
       - No need to control for antibiotics (Antibiotics_Last6mo column) in the UK dataset since all of them are "No"
       - Controlling for marijuana in the UK dataset may be a good idea, but check the literature on the impact of marijuana on the microbiome before filtering
       - Controlling for alcohol in the UK dataset will result in losing half of the dataset
       - China dataset does not provide information for each participant, but the standard deviation (SD) mentioned in Table 1 of the paper can be used as a basis to filter
4. ***Performed analysis***
   * Make sure that both the UK and China datasets are comparable before performing the analysis
   1. **Diversity metrics**
      - This is done to look at the difference between the 6 main categories, which are UK_ecig, UK_tobacco, UK_control, China_ecig, China_tobacco, and China_control
        - Use the 'combined' column
        - Cross-comparing the 6 individuals may result in a more interesting data, but there's also an option to compare the UK and China datasets more generally
      - The results from this analysis will drive the rest of the research questions
      - Aftr running all of the diversity metrics:
        - We will start to see interesting patterns or not and this will help us decide what to focus on
          - We will not be able to focus on all of the interesting observations, instead we need to focus on the interesting observation in which we can collect enough data on
        - We can also know whether the control are the same or not. IF it's the same, the controls should be binned together
   2. **Indicator taxa analysis**
      - This is done to know whether there are specific species of bacteria that are indicative of the different conditions
        - It will answer questions such as: "which bacteria are always present in the UK_tobacco versus the China_tobacco and which bacteria are shared between the two cohorts?" 
   3. **Core Microbiome**
      - This is done to know how many bacteria are share between UK and China versus how many of them are different
   4. **Differential abundance**
      - If the results from Bray Curtis or Shannon Diveristy are interesting, this analysis need to be done
      - This is done to determine the abundance of bacteria, i.e. whether a bacteria is more enriched in one cohort versus another if the cohorts have shared taxa
   5. **Functional analysis**
      - This is done to know what metabolic pathways are present and which metabolic pathways are likely to be upregulated or downregulated based on the taxa that exists on those conditions
        - For example: In the mouth, there are some metabolic pathways that bacteria relied on. By doing functional analysis, we can determine how tobaccor or e-cigarettes disrupt the pathways and whether it differs between UK and China cohorts
      - This analysis will take a while to run, so make sure to not perform the analysis on all of the potential comparisons, just do the analysis on the most interesting comparisons (this is based on the results from the previous analysis)
      - Refer to module optional Module 19 (picrust2)

#### Action items
1. Finish Qiime2 pipeline for China dataset
   - Running the Qiime2 pipeline on the China dataset should be faster since the dataset invovles single end sequences
   - If encountered any problem with the Qiime2 pipeline, email the teaching team and ask for alternative rubric of the proposal
2. Ask Ritu how to merge the phyloseq objects

---
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

---
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



