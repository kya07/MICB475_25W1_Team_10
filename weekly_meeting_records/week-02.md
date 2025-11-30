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
