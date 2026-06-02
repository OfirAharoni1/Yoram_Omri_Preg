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
---

# Projects Overview

| Project Name | Group    | Type | Country | # Samples | V   |
|--------------|----------|------|---------|--------|------|
| PRJNA1247940 | Pregnant | WGS  | China   |  74    | V
| ERP020710    | Pregnant | WGS  | China   |  55    | V
| PRJNA1254708 | Pregnant | 16S  | China   |  63    | V
| omri_stool   | Pregnant | 16S  | Israel? |  382   | V
|--------------|----------|------|---------|--------| --
| PRJEB37731   | Control  | WGS  | Denmark |  160   | V
| PRJNA48479   | Control  | WGS  | USA     |  100   | V  
| PRJNA1067170 | Control  | WGS  | USA     |  124   |    
| PRJNA669650  | Control  | 16S  | Finland |  208   | V  
| PRJNA388263  | Control  | 16S  | USA     |  400   | exporting   



---

# Reminders

-  fig 1B - first building stackbar for each DS and then create the violin plot - the code for only stackbars is at archive folder.
-  
-  

