# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 10:13:35 2015

@author: edogerde
"""
import os
import subprocess
import pandas as pd
import shutil

def reorient_image(input_axes, in_file, typeFile, output_dir):
    """ Rectify the orientation of an image.
    """
    reoriented_file = os.path.join(output_dir, '%s_%s.nii.gz' % (subj_id, typeFile))
    print("ok1")
    cmd = ["AimsFlip", "-i", in_file, "-o", reoriented_file, "-m", input_axes]
    #print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)
    
    return reoriented_file

if __name__ == "__main__":

    path_t1 = '/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/bundle_rd/sujet_18_rt'
    path_RDCT = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/ct_labels_processing"
    output_path = os.path.join("/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/bundle_rd/WM_Atlas")
    subjects_csv = os.path.join(path_t1, "clinical_data.csv")

    # Read the dataframe clinical data: age, orientation of the ct
    df_subjects = pd.read_csv(subjects_csv)


    valid_subject_dirs = [os.path.join(path_t1, dir_name)
                          for dir_name in os.listdir(path_t1)
                          if os.path.isdir(os.path.join(path_t1, dir_name))]
                          
    for subject_path in valid_subject_dirs:
        
            subj_id = os.path.basename(subject_path)
            print subj_id

            t1_nii = os.path.join(subject_path, "mri", subj_id + '_T1.nii')
            ct_nii = os.path.join(path_RDCT, subj_id, "ct", "ct_cut_brain.nii.gz")
           
            if os.path.exists(ct_nii):
                print(subj_id + " have good ct")
                
            else:
                print (subj_id + " no ct") 
                continue 
        
            rd_nii = os.path.join(path_RDCT, subj_id, "ct", "rd_rescale_brain.nii.gz")
        
        
            output_dir = os.path.join(output_path, subj_id)
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)
            
            # Flip the rd and the ct
            flip = df_subjects.check_flip[df_subjects.anonym_nom == subj_id].values[0]
            shutil.copy(t1_nii, os.path.join(output_dir, subj_id +"_T1.nii.gz" ))
        
            if flip == 1:
                reorient_image("XXYY", ct_nii, 'ct', output_dir)
                reorient_image("XXYY", rd_nii, 'rd', output_dir)
             
                print("FLIP")
            else:
                shutil.copy(ct_nii,os.path.join(output_dir, subj_id +"_ct.nii.gz" ))
                shutil.copy(rd_nii,os.path.join(output_dir, subj_id +"_rd.nii.gz" ))