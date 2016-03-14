# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 14:23:16 2015

@author: edogerde
"""
# System import
import numpy
import os
import scipy.signal
import glob
import pandas as pd


# IO import
import nibabel

def mri_to_ct(t1_nii, ct_nii, min_thr, output_dir, flip, verbose= 0):
    
    # Output autocompletion
    ct_modify_nii = os.path.join(output_dir, "ct_modify.nii.gz")
    cut_brain_index_fileName = os.path.join(output_dir, "ct_brain_index.txt")
    
    # Load ct and modify the data for brain extraction
    print 'ct_nii: ', ct_nii
    ct_im = nibabel.load(ct_nii)
    print "ok load"
    ct_data = ct_im.get_data()
    ct_data[numpy.where(ct_data < 0)] = 0
    nibabel.save(ct_im, ct_modify_nii)
    # Detect the neck
    ct_im = nibabel.load(ct_modify_nii)
    ct_data = ct_im.get_data()
    power = numpy.sum(numpy.sum(ct_data, axis=0), axis=0)
    powerfilter = scipy.signal.savgol_filter(power, window_length=11, polyorder=1)
    mins = (numpy.diff(numpy.sign(numpy.diff(powerfilter))) > 0).nonzero()[0] + 1
    global_min = numpy.inf
    global_min_index = -1
    for index in mins:
        if powerfilter[index] > min_thr and global_min > powerfilter[index]:
            global_min = powerfilter[index]
            global_min_index = index
    cut_brain_index_file = open(cut_brain_index_fileName, "w")
    cut_brain_index_file.write(str(global_min_index))
    cut_brain_index_file.close()

    return global_min_index
    
if __name__ == "__main__":


    # Global parameters
    BASE_PATH = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence"
    ATLAS_PATH = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/atlas"
    nii_path = os.path.join(BASE_PATH, "sujet_18_rt")
    output_path = os.path.join(BASE_PATH, "results_from_label_to_rd_after_ants")
    subjects_csv = os.path.join(nii_path, "clinical_data.csv")
    
    # get the appropriate labels
    
    labels_2 = os.path.join(ATLAS_PATH, "atlas_0-2/ANTS2-0Years_label_regis_head_trans.nii.gz")
    labels_5 = os.path.join(ATLAS_PATH, "atlas_5_9/ANTS9-5Years3T_label_regis_head_trans.nii.gz")
    labels_2_5 = os.path.join(ATLAS_PATH, "atlas_2_5/ANTS2-5Years_label_regis_head_trans.nii.gz")


    # Keep the valid subject
    valid_subject_dirs = [os.path.join(nii_path, dir_name)
                      for dir_name in os.listdir(nii_path)
                      if os.path.isdir(os.path.join(nii_path, dir_name))]
    
    #valid_subject_dirs.sort()
    
    # Read the dataframe clinical data: age, orientation of the ct
    df_subjects = pd.read_csv(subjects_csv)  

                         

    # Go through all subjects
    for subject_path in valid_subject_dirs:
        #subject_path = os.path.join(nii_path, 'sujet_024_VM')
        print "Processing: '{0}'...".format(subject_path)
        
        # Get subject id
        if not nii_path.endswith(os.path.sep):
            nii_path = nii_path + os.path.sep
        #subj_id = subject_path.replace(nii_path, "").split(os.path.sep)[0]
        subj_id = os.path.basename(subject_path)
        print subj_id
        # Select the correct atlas according to age
        subj_age = df_subjects.ART[df_subjects.anonym_nom == subj_id].values[0]
        print subj_age
        
        if subj_age < 2:
            template_labels = labels_2
            print " under 2 years old"
        elif subj_age > 5:
            template_labels = labels_5
            print " over 5 years old"
        else:
            template_labels = labels_2_5
            print " between 2 and 5 years old"
        
        # Create output directory and skip processing if already
        output_dir = os.path.join(output_path, subj_id)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        

        #Check whether ct image has to be flipped
        flip = df_subjects.check_flip[df_subjects.anonym_nom == subj_id].values[0]
        
        # Get the t1, the inv_trans (from atlas to t1), the rd and the ct of the patient
        t1_nii = glob.glob(os.path.join(subject_path, "mri", subj_id + '_T1.nii'))[0]
        ct_nii = glob.glob(os.path.join(subject_path, "ct", subj_id + "_ct.nii.gz"))[0]
        rd_nii = glob.glob(os.path.join(subject_path, "rd", "*.nii"))[0]
        transform_warped = os.path.join(subject_path, "transform_warped", "transformInverseComposite.h5")
        
        print "Executing: %s" % (t1_nii)
        print "Executing: %s" % (ct_nii)
        print "Executing: %s" % (rd_nii)

        print "Executing: %s" % (transform_warped)
        
        
        
        global_min_index = mri_to_ct(t1_nii, ct_nii, 50000, output_dir, flip, verbose=0)