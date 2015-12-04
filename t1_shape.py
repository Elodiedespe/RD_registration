# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 19:09:05 2015

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
t1_nii = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/sujet_18_rt/sujet_015_NI/mri/sujet_015_NI_T1.nii"
index = -0.009
t1_im = nibabel.load(t1_nii)
t1_data = t1_im.get_data()
print 't1 shape: ', 
t1_data.shape

t1_modified = numpy.zeros(512,512, 143)
labels_to_ct_native[:, :, :global_min_index] = 0
labels_to_ct_native[:, :, global_min_index:] = labels_data
labels_ct_native_im = nibabel.Nifti1Image(labels_to_ct_native, ct_im.get_affine())
nibabel.save(labels_ct_native_im, labels_ct_native_nii)