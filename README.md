# Files Overview

- `preprocess_16S.py`  
  Input: from Yamas- otu.csv + taxonomy.csv  
  Output: formatted csv for the pipeline: first column SampleId, rest taxa
       
  Note: *need to add group to output name*    
  Note: suits to this structure:  
  Yoram_Omri_Preg/datasets_after_yamas/16S/  
│
├── PRJNA669650/  
│   ├── otu_PRJNA669650.csv  
│   └── taxonomy_PRJNA669650.csv  
│  
├── PRJNA1254708/  
│   ├── otu_PRJNA1254708.csv  
│   └── taxonomy_PRJNA1254708.csv  

- `1_set_cols_per_ds.py`  
Input: csv from Yamas  (SG)  
Output: formatted csv for the pipeline: first column   SampleId, second group (Pregnant/Control)  
Note:

- `fig_1B_full_SG.py`  
  Input: formatted csv  
  Output: stackbar for each csv + violin plot for all  
  Note:  

- `fig_1B_16S.py`  
  Input: formatted csv after MIPMLP  
  Output: stackbar for each csv + violin plot for all  
  Note: 

- `fig_1C_16S.py`    
  Input: formatted csv   
  Output: shannon plot   
  Note:  


- `fig_1C_SG.py`    
  Input: formatted csv   
  Output: shannon plot   
  Note: 


- `fig_1D_SG.py`  
  Input: formatted csv (going over SG folder)  
  Output: GIMIC plot   
  Note: 

- `fig_1D_16S.py`  
  Input: formatted csv (going over 16S folder)  
  Output: GIMIC plot   
  Note: 

- `fig_1E_16S.py`  
  Input: formatted csv (going over 16S folder)  
  Output: miMic / SAMBA outputs and plots 
  Note: taxonomy level 2

- `fig_1E_16S_family.py`  
  Input: formatted csv (going over 16S folder)  
  Output: miMic / SAMBA outputs and plots 
  Note: unite by family level

- `fig_1E_16S_genus.py`  
  Input: formatted csv (going over 16S folder)  
  Output: miMic / SAMBA outputs and plots 
  Note: unite by genus level

- `fig_2A.py`  
  Input:
  1. Microbiome table:
    - rows = samples
    - columns = microbial taxa / bacteria
    - first column = sample IDs / index

  2. Metadata table:
    - rows = samples
    - columns = metadata variables
    - first column = sample IDs / index

  Output:  
  For each project, the script saves:
  - all_microbe_metadata_correlations.csv
  - significant_correlation_counts.csv
  - number_of_significant_correlations.png
  - number_of_significant_correlations.pdf

  It also saves one combined table across all projects.  
  


---

# Projects Overview

| Project Name | Group    | Type | Country | # Samples | V   | metadata |
|--------------|----------|------|---------|--------|------| -------- | 
| PRJNA1247940 | Pregnant | WGS  | USA   |  74    | V     | BMI, Consent_Age, Delivery_EGA, country, Preeclampsia_During_Labor, PreTerm_Labor_EGA, tobacco_use, Birth Weight, delivery_type, sex_of_baby, PretermLabor, GDM, hypertension, Preeclampsia_At_PP, StillBirth |  
| ERP020710    | Pregnant | WGS  | China   |  55    | V     | Age, Body-mass index, Disease status, Gestational weeks (fecal sample collection), Gestational weeks (OGTT testing) , Height, Total mass, (manually - https://gigadb.org/dataset/100326 )  | 
| PRJNA1254708 | Pregnant | 16S  | China   |  63    | V     | none - need to contact authors |
| omri_stool   | Pregnant | 16S  | Israel? |  382   | V     | BMI, Age, glucose, number of pregnancies, number of deliveries, way of conception, Antibiotics, Aspirin, Delivery_week, chronic medications, Smoking, Newborn_weight, stress test, education years, calories per day, carbs per day, food preferces |
|--------------|----------|------|---------|--------| ----  | ------- |
| PRJEB37731   | Control  | WGS  | Denmark |  160   | V     | host_body_mass_index, travel_outside_the_country_in_last_six_months, host_disease_status, presence_of_pets_or_farm_animals, drug_usage, urine/urogenital_tract_disorder |
| PRJNA48479   | Control  | WGS  | USA     |  100   | V     | nothing relevant
| PRJNA1067170 | Control  | WGS  | USA     |  124   | need permision to download  | 
| PRJNA669650  | Control  | 16S  | Finland |  208   | V     | Host_age, Host_disease (null/control)  |
| PRJNA388263  | Control  | 16S  | USA     |  400   |problem exporting |  
 


---
# fig 2A status  


| Project Name | Group    | Type | Status |
|--------------|----------|------|--------|
| PRJNA1247940 | Pregnant | WGS  | 55 significant correlations (to batch) |
| ERP020710    | Pregnant | WGS  | 0 significant correlations | 
| PRJNA1254708 | Pregnant | 16S  | skip - no metadata | 
| omri_stool   | Pregnant | 16S  | 165 significant correlations  |
|--------------|----------|------|--------|
| PRJEB37731   | Control  | WGS  | 20 significant correlations       | 
| PRJNA48479   | Control  | WGS  | skip - no metadata | 
| PRJNA669650  | Control  | 16S  | no usable metadata  | 




---

# Reminders

-  fig 1B - first building stackbar for each DS and then create the violin plot - the code for only stackbars is at archive folder.
-  
-  

