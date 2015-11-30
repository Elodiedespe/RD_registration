# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 12:20:12 2015

@author: edogerde
"""

# System import
import numpy
import os
import subprocess
import pandas as pd
import glob

# IO import
import nibabel

def to_id(in_file):
    """ Convert file to itself (get .minf file).
    """
    cmd = ["AimsFileConvert", "-i", in_file, "-o", in_file]
    #print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)
    return in_file

def reorient_image(input_axes, in_file, output_dir):
    """ Rectify the orientation of an image.
    """
    reoriented_file = in_file.split('.')[0] + "_reorient.nii.gz"
    reoriented_file = os.path.join(output_dir, os.path.basename(reoriented_file))
    cmd = ["AimsFlip", "-i", in_file, "-o", reoriented_file, "-m", input_axes]
    #print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)
    return reoriented_file


def labels_to_ct(global_min_index, t1_nii, ct_cut_reoriented_nii, transformation, labels_nii, transform_warped, output_dir):
    """ Register the labels to the CT.
    """

    # Output autocompletion
    registered_labels_nii = os.path.join(output_dir, "labels_to_cut_ct.nii.gz")
    labels_ct_nii = os.path.join(output_dir, "labels_to_ct.nii.gz")
    convert_trans_itk = os.path.join(output_dir, "convert_trans_itk.txt")
    labels_ct_native_nii = os.path.join(output_dir, "labels_to_ct_native.nii.gz")

    # Convert affine transformation from fsl2Ras (ants)
    cmd = "c3d_affine_tool " + \
          " -ref " + ct_cut_reoriented_nii + \
          " -src " + t1_nii + " " + transformation + \
          " -fsl2ras " + \
          " -oitk " + convert_trans_itk
    print ("Executing: " + cmd)
    os.system(cmd)

    # Apply the affine and warp transformation
    print (labels_nii)
    print (registered_labels_nii)
    print (ct_cut_reoriented_nii)
    print (transform_warped)
    print (convert_trans_itk)
    cmd = "antsApplyTransforms " + \
          " --float " + \
          " --default-value 0 " + \
          " --input %s " % (labels_nii) + \
          " --input-image-type 3 " + \
          " --interpolation NearestNeighbor " + \
          " --output %s " % (registered_labels_nii) + \
          " --reference-image %s " % (ct_cut_reoriented_nii) + \
          " --transform [ %s , 0] [%s , 0] " % (convert_trans_itk, transform_warped)
    print ( "executing: " + cmd)
    os.system(cmd)

    if flip == 1:
        labels_ct_nii = reorient_image("XXYY", registered_labels_nii, output_dir)
    else:
        labels_ct_nii = to_id(registered_labels_nii)

    # Send the label to the native ct

   # Load the image
    ct_im = nibabel.load(ct_nii)
    ct_data = ct_im.get_data()
    labels_data = nibabel.load(labels_ct_nii).get_data()
    print 'label shape: ', labels_data.shape
    print 'ct shape: ', ct_data.shape
    print 'global index = ', global_min_index
    labels_to_ct_native = numpy.zeros(ct_data.shape)
    labels_to_ct_native[:, :, :global_min_index] = 0
    labels_to_ct_native[:, :, global_min_index:] = labels_data
    labels_ct_native_im = nibabel.Nifti1Image(labels_to_ct_native, ct_im.get_affine())
    nibabel.save(labels_ct_native_im, labels_ct_native_nii)

    return labels_ct_nii, labels_ct_native_nii



if __name__ == "__main__":


    # Global parameters
    BASE_PATH = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence"
    ATLAS_PATH = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/atlas"
    nii_path = os.path.join(BASE_PATH, "sujet_18_rt")
    output_path = os.path.join(BASE_PATH, "results_from_label_to_rd_after_ants")
    subjects_csv = os.path.join(nii_path, "clinical_data.csv")

    # get the appropriate labels

    labels_2 = os.path.join(ATLAS_PATH, "atlas_0-2/ANTS2-0Years_label_regis_head.nii.gz")
    labels_5 = os.path.join(ATLAS_PATH, "atlas_5_9/ANTS9-5Years3T_label_regis_head.nii.gz")
    labels_2_5 = os.path.join(ATLAS_PATH, "atlas_2_5/ANTS2-5Years_label_regis_head.nii.gz")


    # Keep the valid subject
    valid_subject_dirs = [os.path.join(nii_path, dir_name)
                      for dir_name in os.listdir(nii_path)
                      if os.path.isdir(os.path.join(nii_path, dir_name))]

    #valid_subject_dirs.sort()

    # Read the dataframe clinical data: age, orientation of the ct
    df_subjects = pd.read_csv(subjects_csv)



    # Go through all subjects
    for subject_path in valid_subject_dirs[1:9]:
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

        if flip == 1:
             ct_cut_reoriented_nii = os.path.join(output_dir, "ct_cut_brain_reorient.nii.gz")
        else:
             ct_cut_reoriented_nii = os.path.join(output_dir, "ct_cut_brain.nii.gz")
        print "NO flip"

        transformation= os.path.join(output_dir, "t1_to_cut_ct.txt")

        print "Executing: %s" % (t1_nii)
        print "Executing: %s" % (ct_nii)
        print "Executing: %s" % (rd_nii)

        print "Executing: %s" % (transform_warped)

        ct_cut_brain = os.path.join(output_dir, 'ct_cut_brain.nii.gz')
        cut_brain_index_fileName = os.path.join(output_dir, "ct_brain_index.txt")
        cut_brain_index_file = open(cut_brain_index_fileName, "r")
        cut_brain_index = int(cut_brain_index_file.read())
        cut_brain_index_file.close()
        print "cut_brain_index: %s" % (cut_brain_index)
        labels_ct_nii = labels_to_ct(cut_brain_index, t1_nii, ct_cut_reoriented_nii, transformation, template_labels, transform_warped, output_dir)
