"""
RUN ME FIRST
on the datasets from yamas

we get formatted csv where:
1) first column named: SampleID
2) second column is group: Pregnant or Control

"""


import pandas as pd
import os

# ===== SETTINGS =====
file_path = r"Yoram_Omri_Preg/datasets_after_yamas/ERP020710_Pregnant_SG.csv"
group_value = "Pregnant"   # or "Control"

# ===== LOAD =====
df = pd.read_csv(file_path, skiprows=1, sep=",")

# ===== RENAME FIRST COLUMN (ID -> SampleID) =====
df.rename(columns={df.columns[0]: "SampleID"}, inplace=True)

# ===== ADD Group COLUMN AS SECOND COLUMN =====
df.insert(1, "Group", group_value)

# ===== SAVE AS NEW FILE =====
base, ext = os.path.splitext(file_path)
out_path = f"{base}_formatted{ext}"
df.to_csv(out_path, index=False)

print("Done!")