# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:15:42 2015

@author: cp243490
"""

import os
import pandas as pd
import glob
import shutil
import pandas as pd
import numpy as np
import nibabel as nib


def registerVolume(reference, Volume, mode, outputTrm, outputITrm):
    command = 'AimsMIRegister ' + \
              ' -r ' + reference + \
              ' -t ' + Volume + \
              ' --dir ' + outputTrm + \
              ' --inv ' + outputITrm + \
              ' --error ' + str(0.0005) + \
              ' --interpolation 3 '
    #print "Executing: ", command
    os.system(command)


def composeTransformation(trmFile1, trmFile2, outputTrm):
    command = 'AimsComposeTransformation ' + \
              ' -i ' + trmFile2 + \
              ' -j ' + trmFile1 + \
              ' -o ' + outputTrm
    #print "Executing: ", command
    os.system(command)


def resampleVolume(scalarVolumeFile, trmFile, order, output_file):
    '''
    Scalar parameter map normalisation in Children Template referential
    '''
    # Dimensionof the image in template referential
    xTemplate = 88
    yTemplate = 81
    zTemplate = 71
    # Resolution of the Volume
    xRes = 2
    yRes = 2
    zRes = 2
    # Volume Resampling
    command = "AimsResample " + \
              " -i " + scalarVolumeFile + \
              " --dx %i --dy %i --dz %i " % (xTemplate, yTemplate, zTemplate) + \
              " --sx %i --sy %i --sz %i " % (xRes, yRes, zRes) + \
              " -m " + trmFile + \
              " -o " + output_file + \
              " -t %i" % (order)
    #print "Executing: ", command
    os.system(command)
    return output_file


def template_to_ct(templateFile, t1File, ctFile, T1ToCT_trm, output_path):
    # templateT1 to MRI
    TemplateToT1_trm = os.path.join(output_path, 'Template_To_T1.trm')
    T1ToTemplate_trm = os.path.join(output_path, 'T1_To_template.trm')
    registerVolume(t1File, templateT1,
                   'mutual-information',
                   TemplateToT1_trm, T1ToTemplate_trm)

    # Template to CT
    TemplateToCT_trm = os.path.join(output_path, 'Template_To_CT.trm')
    composeTransformation(TemplateToT1_trm, T1ToCT_trm, TemplateToCT_trm)


if __name__ == "__main__":

    ROOT_PATH = '/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/elodie_dosimetry'

    SUBJECTS_DATA_PATH = os.path.join(ROOT_PATH, 'subjects_data')
    ATLAS_DIFF_PATH = os.path.join(ROOT_PATH, 'Atlas_Pediatric_2Hemispheres')
    subjects_csv = os.path.join(ROOT_PATH, "clinical_data.csv")

    templateT1 = os.path.join(ATLAS_DIFF_PATH, 'template_children_final.ima')

    TO_CT_PATH = os.path.join(ROOT_PATH, '01-TO_CT')
    if not os.path.exists(TO_CT_PATH):
        os.mkdir(TO_CT_PATH)

    # Read subjects
    df_subjects = pd.read_csv(subjects_csv)
    subjects = df_subjects['anonym_nom'].values

    for subject in subjects[4:5]:
        subject_path = os.path.join(SUBJECTS_DATA_PATH, subject)
        age = df_subjects.ART[df_subjects['anonym_nom'] == subject].values[0]

        if os.path.exists(subject_path):
            ctFile = os.path.join(subject_path, subject +"_ct.nii.gz")
            t1File = os.path.join(subject_path, subject +"_T1.nii.gz")
            rdFile = os.path.join(subject_path, subject +"_rd.nii.gz")
            t1TOct_trm = os.path.join(subject_path, 't1_To_ct.trm')

            # Register atlas template into ct referential
            #################################################
            toCT_subject_path = os.path.join(TO_CT_PATH, subject)
            if not os.path.exists(toCT_subject_path):
                os.mkdir(toCT_subject_path)
    
            template_to_ct(templateT1, t1File, ctFile, t1TOct_trm,
                           toCT_subject_path)

        else:
            print "path %s does not exist" % subject_path
        
