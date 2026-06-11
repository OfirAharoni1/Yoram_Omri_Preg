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

- `fig_2A_16S.py`  
  Input: 2 datasets- Preg + control + metadata file for each  
  
  Output:    
  pregnancy_vs_control_significant_correlations/   
├── all_microbe_metadata_correlations_by_group.csv  
├── significant_correlation_counts_by_group.csv  
├── pregnancy_vs_control_number_of_significant_correlations.png
└── pregnancy_vs_control_number_of_significant_correlations.pdf       
  
  Note: THIS IS AN INITIAL DRAFT - NEED TO MAKE CHANGES - JUST A BASE FOR NOW


---

# Projects Overview

| Project Name | Group    | Type | Country | # Samples | V   | metadata |
|--------------|----------|------|---------|--------|------| -------- | 
| PRJNA1247940 | Pregnant | WGS  | USA   |  74    | V     | BMI, Consent_Age, Delivery_EGA, country, Preeclampsia_During_Labor, PreTerm_Labor_EGA, tobacco_use, Birth Weight, delivery_type, sex_of_baby, PretermLabor, GDM, hypertension, Preeclampsia_At_PP, StillBirth |  
| ERP020710    | Pregnant | WGS  | China   |  55    | V     | Age, Body-mass index, Disease status, Gestational weeks (fecal sample collection), Gestational weeks (OGTT testing) , Height, Total mass, (manually - https://gigadb.org/dataset/100326 )  | 
| PRJNA1254708 | Pregnant | 16S  | China   |  63    | V     | none - need to contact authors |
| omri_stool   | Pregnant | 16S  | Israel? |  382   | V     | BMI, Age, glucose, number of pregnancies, number of deliveries, way of conception, Antibiotics, Aspirin, Delivery_week, chronic medications, Smoking, Newborn_weight, stress test, education years, calories per day, carbs per day, food preferces |
|--------------|----------|------|---------|--------| ----  | ------- |
| PRJEB37731   | Control  | WGS  | Denmark |  160   | V     |
| PRJNA48479   | Control  | WGS  | USA     |  100   | V     |
| PRJNA1067170 | Control  | WGS  | USA     |  124   |       |
| PRJNA669650  | Control  | 16S  | Finland |  208   | V     |
| PRJNA388263  | Control  | 16S  | USA     |  400   |problrm exporting |  



---

# Reminders

-  fig 1B - first building stackbar for each DS and then create the violin plot - the code for only stackbars is at archive folder.
-  
-  

