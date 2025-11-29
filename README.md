# MICB475_25W1_Team_10
This repository contains weekly **meeting agendas and minutes** for **Team 10 (Aimee, Kaiya, Kayla and Sophie)** in UJEMI Project 2 (MICB475).
Agendas are posted before each meeting for mentor review, and minutes are added to record key discussions, decisions, and next steps. Team members will alternate responsibility for maintaining these records to ensure consistent tracking of project progress.

---
### Weekly Agenda, Meeting Minutes, and Action Items
- [Week 1 (02-Oct-2025)](#week-1-02-oct-2025)
- [Week 2 (09-Oct-2025)](#week-2-09-oct-2025)
- [Week 3 (16-Oct-2025)](#week-3-16-oct-2025)
- [Week 4 (23-Oct-2025)](#week-4-23-oct-2025)
- [Week 5 (30-Oct-2025)](#week-5-30-oct-2025)
- [Week 6 (06-Nov-2025)](#week-6-06-nov-2025)
- [Week 7 (20-Nov-2025)](#week-7-20-nov-2025)
- [Week 8 (20-Nov-2025)](#week-8-27-nov-2025)

### Folder Clarification/ Notes
* Add in notes if applicable. 
---
### Week 8 (27-Nov-2025)
#### Agenda
* Obtain feedback for the presentation slide draft
* Confirm interpretation of data
* Feedback for 'conclusive statement' titles
* Discuss which information doesn't add to our presentation 

#### Meeting Minutes
* **Presentation slides notes for Sophie**
   - Introduction of topic and Importance of the topic: 
      - 3 columns:
        - Brief introduction of e-cig and tobacco
        - Link microbiome and smoking
        - Importance of study: Not much research across cohorts
   - Dataset overview and Methods: 
      - Comparison of the dataset (2-column comparison) 
      - Boxes showcasing what matrices we ran
   - Others:
        - Only show the ecig slide 
        - Ignore smoking one 
        - No difference in US
        - Differences in China, and that the changes are BAD. (functionally, they don’t look) 


* **Presentation slides notes for Kayla**
   - Overall: 
        - Consistency in group graphs by removing grey line 
        - Colourblind-saved colour?  
        - Alpha diversity title is good as a conclusive statement 
        - China box plot 1 shape and US plot another shape 
        - Need to simplify text 
   - Must mention for ALPHA diversity: ← Can be cut for presentation
        - Baselines between populations are different 
        - US_control baseline is different from CN_control
        - Less changes in microbiome in US cohort 
        - Double check the p-values → observed might be missing some significant ones 
        - Above 3 points → highlight each one with a box as an animation as points appear 
        - **Baseline Difference**
        - Cut out alpha as show similar to Beta
        - Less change in US cohort’s microbiome 
   - Must mention for Beta diversity - Bray-Curtis 
        - NO difference in US 
        - Only changes in the CN cohort → a different title is needed 
        - No text needed → easily understandable 
   - Venn diagrams 
        - Conclusions/title change → reflect how much of the profile is shared (i.e., More significant changes in CN cohort than US) 
        - Us ⅓ 
        - CHINA only 15% 
        - Supports the overall conclusion

* **Presentation slides notes for Aimee**
   - Remove the first statement for each of the points for each cohort 
   - Make the font larger
   - 2nd slide 
        - Have the table on the side and only add the bolded statements 
        - No full sentences needed 
   - Title → for indicator taxamaybe mention that it's mainly pathogenic or non-pathogenic → describe the profile a bit  
        - US seems to be mainly shaped with commensal non-pathogenic bacteria 
        - CN seems to be mainly shaped with pathogenic bacteria 
   - Maybe remove phylum to make it bigger 

* **Presentation slides notes for Kaiya**
   - Box plot, revise: 
        - Orange - China
        - Blue - US 
        - Not much of a concern though, since present phylum-count 
   - Use over and underrepresented instead 
   - Don't use changes as it is not a time-course analysis, but differences 
   - Microbes in the blue are 5x greater relative to the control
   - Remove tobacco vs ecigarettes for the presentation only
   - DIFFERENCES (not word “changes”)

#### Action items
* Make your title (conclusive statement) shorter → put into conclusions!
* For whole story (will be your conclusion slide): 
   - Differences between cohorts suggest that each cohort’s microbiome responds differently to tobacco and e-cigarette use.
   - Not much changes in US 
   - More changes in China cohort 
   - Is this profile change (ie increase) good or bad tho? 
      - Accordingly to functional analysis, profile changes in Chinese cohort is NOT positive 
* Change in diversity is not always good, can be unhealthy / negative profile 
* Cut alpha? Confirm. Cut alpha diversity for presentation. 

---
### Week 7 (20-Nov-2025)
#### Agenda
* Present Results and find interesting findings
* Get feedback on ways to improve visualizations

#### Meeting Minutes
* **Diversity Metrics**
  * Alpha Diversity:
      1. Shannon Diversity:
         - cn_control have lower baseline compared to cn_ecig and cn_tobacco
         - cn_control lower than us_control
         - us_control has very high baseline
         - us_control has the same diversity as us_tobacco
         - Similar to what was observed in the Columbia study where participants are Columbians but eat westernized diet
         - Suggests that diet may play a part in the difference in baseline between US and CN
      2. Observed features: 
         - The observed features of US and CN tobacco are very similar to their Shannon Diversity and that the US cohort has higher baseline than CN cohort
         - Implies that the effect of tobacco on diversity may depend on how much bacteria you have in the first place
  * Beta Diversity:
    - Beta diversity of each smoking status overlaps in the US cohort meaning that the beta diversity between the groups doesn’t change much
    - Meanwhile in the CN cohort, the non-smokers cluster is distinct from the tobacco and ecig smokers clusters, which suggests that beta diversity between smoking status in the CN cohort change significantly
    - In addition, beta diversity of CN ecig and tobacco overlaps
  * Statistical Analysis
    - For alpha diversity, use Kruskal-Wallis and Dunn test
    - It will provide a more accurate p-value compared to the anova
* **Diversity Metrics - *Tips on visualizations***:
  * Alpha Diversity:
    - Delete the non-significance, only show the significance ones
    - In the figure legend, mentioned that “p-value above 0.05 is not shown”
    - Re-label x-axis using the format: (Cohort, Smoking Status). E.g. (China, Control), etc.
    - Re-format figure legend accordingly
    - Boxplot should be colored based on the color used in beta diversity PCOA plot
    - Alternative option: when saving the plot from R, use a bigger width so that the significance can be manually added in powerpoint
  * Beta diversity: 
    - PERMANOVA results should be mentioned in the figure legend
    - Re-label the legends; keep it consistent with alpha diversity
  * Overall: 
    - All figures throughout the paper should have a consistent format
    - Same theme, same color coding, etc.
    - All the figures should have the same plot theme
    - Currently, the plot theme (i.e. background) is grey, recommended to change to a different theme
    - Recommended to use text size 18
    - When making plots, keep in mind that the plot will be small when put into the manuscript
    - So any text or labels should be visible even if the plot size is small
    - In the manuscript, usually alpha diversity plots will be above beta diversity plots

* **Core Microbiome**
  * Takeaways: 
    - Different pattern between us and china cohorts
    - Use DESeq or Indicator Taxa Analysis to identify the bacteria unique for each category
    - Initially Dr. Sun recommended to use a table generated from the core microbiome analysis to identify the unique bacteria for each category, but Dr. Avril disagree since there will be too much data
  * TIPS ON VISUALIZATION: 
    - Use only one color
    - Use the venn diagram separated based on cohort, i.e. US, CN
    - No need to include the large venn diagram
    - For publication, it is advised to not display one result in multiple formats
    - Use bigger fonts for the label

* **Indicator Taxa Analysis**
  * TIPS:
    - Look at past UJEMI papers to see how they format the indicator taxa analysis table
    - The indicator taxa analysis table should be a big table which includes the unique taxa for each categories
    - Just need to remove the ASV column since the name of the ASV will change every time you run the code
    - Columns required: Stat =use 3 d.p. since the p-value uses 3 d.p.; P-value; Taxa = includes at least the phylum, family, genus
    - Treat “NA” or “Incertae_sedis” as  “unclassified”
  *Takeaways: 
    - Taxa in cn_control is pathogenic
    - There’s some weird results in us_ecig since there’s lactobacillus and commensal bacteria, which are considered “good” bacteria
  *TIPS ON ANALYSIS: 
    - Oral health can affect the microbiome
    - Look for literature explaining this relationship
    - Also, look for periodontitis
    - Need to see functional analysis results to make sense of it

* **DESeq2**
  * TIPS ON VISUALIZATION: 
    1. Volcano Plot: 
        - Remove the “significant” legend in the volcano plot. Instead will be added to the Manuscript’s figure legend
        - In the figure legend, mentioned that red is not significant and blue is significant
        - Remove volcano plots that compare between US vs China
        - Diversity metrics (see Part I) is enough to show that there’s a significant difference between cohorts
        - Differences between cohorts will not be the main focus for our manuscript; we will be focusing differences WITHIN cohorts
        - Alternative: Can include the plots that compare between US vs China in the supplemental figures section
        - Only keep WITHIN cohorts (2 by 3 table)
        - Volcano plot information is very general → then connect to Bar Plot’s Count table

    2. Bar Plot & associated Count Table: 
        - Since there’s a lot of phylum, using the DESeq bar plot won’t be meaningful; instead, make a “phylum count” table that contains a list of phyla and their count (basically what we have done is correct)
        - Resolving to the phylum level will be the best choice
        - Using p-adj <0.01 and log 2-fold are okay
        - The phylum count table won’t be included in the manuscript, but the bar plot that was created based on this table will be included
        - As the previous volcano plot and the present bar plot showcasing counting results can be presented together in manuscript
        - For the “phylum count” table, no need to include
        - the ecig vs tobacco comparison for both US and CN
        - All comparisons between US and CN
        - Create a **bar graph based on the “phylum count” table** (i.e., graphing the absolute counts of the microbe that are enriched or depleted):
           - x-axis = Phylum
           - y-axis = Count
           - Each phylum will have 2 bars: count below 2 = color: red, label: “depleted”  & count above 2 = color: blue, label: “enriched”
           - The coloring of the bars should be consistent with the coloring used in the volcano plots   (Enriched vs Depleted)
           - In total there will be 4 bar graphs; one graph for each comparison us_ecig vs us_control, us_tobacco vs us_control, cn_ecig vs cn_control, cn_tobacco vs cn_control
           - ignore the ecig vs tobacco comparison for US and CN
  * TAKEAWAYS: 
    - The phylum that changed the most in each cohort
    - US = Bacteriodota
    - CN = Bacillota (previous name: Firmicutes)
    - The phyla are similar within each cohort
    - DESeq results correlate with alpha diversity results

* **Functional analysis**
    - In progress
    - Sophie had issues with running the functional analysis in the server and that’s because the way the server was set up for this term is that it has limited RAM storage and that the storage is not big enough to run the functional analysis
    - Dr. Avril Metcalfe-Roache promised to send Sophie the functional analysis later today
    - UPDATE: Sophie received Dr. Avril’s email and the functional analysis results at 2AM on Nov 21

* **Overall project**
    - The main storyline is that US and CN cohorts have very distinct microbiome and the effect of smoking status on those two cohorts are very different
    - I.e. Depending on where you’re from, you’ll have different baseline and the smoking status will affect the salivary microbiome differently

#### Action items
  *  Watch Module 17 regarding presentation:
  * Create presentation slides, including the edited plots
    - So that we can receive feedback before our presentation submission
  * TIPS FOR PRESENTATION SLIDES: 
    - Keep in mind that the presentation is only 10 minutes
    - Aim for 10 slides
    - Don’t take too much time explaining the background
    - Aim for 2 slides of background
    - Slides should be a balance between text and figures
    - Don’t write paragraphs; use point form instead
    - The text should be the main points of the results that you want other group to address
    - For functional analysis, include the plots only if the result is not complicated
    - It’s not allowed to add speaker notes when making our presentation
    - But, when we’re presenting other group’s slides, we can bring and read our notes during the presentation
    - If the slides that we are presenting is very bad, it will affect the other group’s grade, not ours


### Week 6 (06-Nov-2025)
#### Agenda
* Recap merging of US and CN phyloseq objects (tree is dropped for downstream analyses)
* Discuss results obtained: 
  * merged phyloseq object's rarefaction curve before vs after filtering (the number of samples went from 123 to 44)
  * α-diversity: Shannon Index with Kruskal-Wallis results
  * β-diversity: Bray-Curtis with PERMANOVA
  * The filtering that was done is very logical according to our TA

#### Meeting Minutes
* Discussion about the results obtained from the diversity metrics
  * 9 groups were selected to proceed with the next phase of our aims: us_control v.s us_tobacco, us_control v.s us_e-cigarrete, cn_control v.s cn _tobacco, cn_control v.s cn_e-cigarrete, us_control v.s cn_control, us_tobacco v.s cn_tobacco, us_e-cigarrete v.s cn_e-cigarrete, us_e-cigarrete v.s us_tobacco, cn_e-cigarrete v.s cn_tobacco
  * Significant differences were noted in: us_control v.s cn_control, us_tobacco v.s cn_tobacco, cn_control v.s cn_e-cigarrete and us_e-cigarrete v.s us_tobacco
  * Curiosly no differences were noted in: cn_control v.s cn _tobacco, cn_e-cigarrete v.s cn_tobacco, us_e-cigarrete v.s cn_e-cigarrete, us_control v.s us_e-cigarrete, us_control v.s us_tobacco
  * Differences between the control groups from each cohort may be influenced by external factors
  * us_control and cn_e-cigarretes don’t have a significance difference, which is interesting but might not add to the story
  * PCoA plots showed that cn control is different from cn tobacco and cn e-cigarettes groups which overlap with one another. On the other hand, for the US groups, all of them overlap with one another so no difference between them is noted.
 
#### Action items
*	Make another table with only the significant ones to ensure clarity
*	Indicator Taxa, DEseq and Core Microbiome analyses are expected to be completed by next week on each of 9 groups selected
*	Functional analysis will be performed after Aims 2-4 are completed
*	Tasks were splited between each member: Kaiya-DEseq, Kayla-Core Microbiome, Aimee-Indicator Taxa Analysis, Sophie- Functional Analysis


---
### Week 5 (30-Oct-2025)
#### Agenda
* Address questions both within and outside the scope of the proposal. 
* Attempt and troubleshoot merging of the two phyloseq objects.

#### Meeting Minutes
* Generation of alpha rarefaction curve 
   * For the alpha rarefaction curve, there are two options: either you can make one curve per dataset or make one curve representing both datasets.
   * If decided to go for one curve per dataset, use the sequencing depth of the dataset that has the lower sampling depth.
   *  Sequencing depth for CN is 10870 and US is 16,432; thus, CN has the lower sequenicng depth. 
   * Regardless of which option the team decides to go with, remember to mention that in the dataset overview section clearly.
   * May generate another rarefaction curve after phyloseq objects merging and compare with before merging (for your interest only)
* Filtering parameters _after merging_ check
    * Filter out: NAs, Weed column (from us_dataset), bmi for e-cigarette < 18.5, bmi for non-smoker(cntl) < 18.5
    * Restrict age for each group to: e-cigarette (age 24-48), tabacco (age 28-47) and non-smoking control (29-54)
 
#### Action items
* Continue proposal writing
* Merge us and cn phyloseq objects into one single phyloseq object (Merging command: merge_phyloseq(phyloseq1, phyloseq2))
* Update R scripts to ensure consistency in GitHub organization

---
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

#### Meeting Minutes
* Logistics / Clarifications: 
     * Build separate phyloseq objects per dataset, then merge.
     * Required components for each phyloseq object: Rooted tree (phy_tree); OTU/ASV table (otu_table); Representative sequences (refseq); Taxonomy table (tax_table)
* Sequencing depth for each dataset should be reported separately before merging (document per-dataset depth stats)
     * Sequencing depth for CN is 10870 and US is 16,432
     * Filter accordingly after merging based on dataset limitations and to control for confounding variables
     * Proceed to diversity and differential analyses
* Poposal's introduction
   * Provide concise scientific background from prior studies directly tied to our research question.
   * Keep narrative connected; every citation should motivate or justify the question.
   * End with a clear gap statement → leads to our hypothesis.
* Methods:
   * Specify exact statistical tests (e.g., Wilcoxon rank-sum, Kruskal–Wallis, PERMANOVA/ADONIS).
   * Specify diversity metrics: Alpha: Shannon, Observed ASVs, Faith’s PD (if tree available). Beta: Bray–Curtis, UniFrac (unweighted/weighted).
   * Document rarefaction depth decision and filtering criteria.  

#### Action items
* Discuss on filtering parameters that come after phyloseq object merging
* Come up with a feasible project timeline 
  
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



