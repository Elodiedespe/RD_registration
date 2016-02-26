# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 12:01:19 2016

@author: edogerde
"""

import os 
import pandas as pd
import shutil 

Path= "/home/edogerde/Bureau/extraction_roi"
Path_label= "/media/edogerde/MY PASSPORT/data_dosimetry/results_from_label_to_rd_after_ants"
Path_rd= "/media/edogerde/MY PASSPORT/these/sujets/subject_24_rd_resample"
# Create new empty directories with the name of the patients
Suj = pd.read_csv("/media/edogerde/MY PASSPORT/clinical_dataFinal.csv")
sub_name= Suj["anonym_nom"]
for name in sub_name:
    print(name)
    os.mkdir(os.path.join(Path, name))

#Create a list of patients directories
valid_subject_dirs = [os.path.join(Path, dir_name)
                      for dir_name in os.listdir(Path)
                      if os.path.isdir(os.path.join(Path, dir_name))]


#copy all the file in the new directories
"""For each patient, check if the registration is done or not. If the 
reagistration is done, take the register label and the RD of the patient and 
copy it in the new directories"""

for sujb_path in valid_subject_dirs:
    suj_id= os.path.basename(sujb_path)
    Register= Suj[["Recalage","anonym_nom"]]
    Subj_Recalage = Register.Recalage[Register.anonym_nom == suj_id].values[0]
    if Subj_Recalage == 1.0:
        label = os.path.join(Path_label, suj_id,"labels_to_rd.nii.gz")
        rd = os.path.join(Path_rd, suj_id, "rd", suj_id + ".nii")
        NewDirlabel = shutil.copy(label, os.path.join(Path, suj_id, "labels_to_rd.nii.gz"))
        NewDirRD = shutil.copy(rd, os.path.join(Path, suj_id, suj_id +".nii"))    
    else:
        continue
    
   
 

