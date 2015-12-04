# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 16:22:28 2015

@author: edogerde
"""
# System import

import os
import sys

def transformation(transfo_input, direction, source, output_dir):
    """
    Convert .txt file to .trm file
    """
    transfo_output = os.path.join(output_dir, "t1_To_ct.trm" )
    print ("ok")
    cmd = "cartoFSLmatToTrm.py " + \
            " -i " + transfo_input + \
            " -o " + transfo_output + \
            " -d " + direction + \
            " -s " + source 
    print("ok2")               
    #print "Executing: ", command
    os.system(cmd)
    print("ok3")
    
    return transfo_output
    
if __name__ == "__main__":
    
    BASE_PATH="/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence"
    input_path = os.path.join(BASE_PATH, "results_from_label_to_rd_after_ants")
    im_path =  os.path.join(BASE_PATH, "elodie_dosimetry", "subjects_data")
    output_path= os.path.join(BASE_PATH, "elodie_dosimetry", "subjects_data")
    
    
    image_dirs = [os.path.join(im_path, dir_name)
                  for dir_name in os.listdir(im_path)
                  if os.path.isdir(os.path.join(im_path, dir_name))]
     
    for image_path in image_dirs:
        subj_id = os.path.basename(image_path)
        transfo_input= os.path.join(input_path, subj_id, "t1_to_cut_ct.txt")
        direction= os.path.join(image_path, "%s_ct.nii.gz" % subj_id)
        source = os.path.join(image_path, "%s_T1.nii.gz" % subj_id)
        output_dir = os.path.join(output_path, subj_id)
        print subj_id
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
    
        
        if not os.path.exists(direction):
            print (subj_id + " ct found")
            continue    
        
        transfo_output = transformation(transfo_input, direction, source, output_dir)
                          