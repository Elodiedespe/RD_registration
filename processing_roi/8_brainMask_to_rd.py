# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 17:47:48 2016

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

def Extrac_brain(in_file, output_dir):
    """ Extract the brain.
    """
    BrainMask_nii = os.path.join(subject_path, subj_id + "brain_mask.nii.gz")  
    Brainextrac_nii = os.path.join(subject_path, subj_id + "T1_brain.nii.gz") 
    
    cmd = ["Bet", in_file, Brainextrac_nii, "-f", str(0.5), "-g", str(0), "-m",BrainMask_nii ]
    print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)
    
    return BrainMask_nii
    
def BrainMask_to_ct(global_min_index, ct_cut_reoriented_nii, transformation, BrainMask_nii, output_dir):
    """ Register the brainMask to the CT.
    """

    # Output autocompletion
    registered_BrainMask_nii = os.path.join(output_dir, "BrainMask_to_cut_ct.nii.gz")
    BrainMask_ct_nii = os.path.join(output_dir, "BrainMask_to_ct.nii.gz")
    BrainMask_ct_native_nii = os.path.join(output_dir, "BrainMask_to_ct_native.nii.gz")

    # Apply the affine and warp transformation
    registered_BrainMask_nii = os.path.join(output_dir, "BrainMask_to_cut_ct.nii.gz")

    # Warp the labels
    cmd = ["flirt", "-in", BrainMask_nii,
           "-ref", ct_cut_reoriented_nii, "-applyxfm", "-init", transformation, "-out", registered_BrainMask_nii, "-interp", "nearestneighbour"]
    print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)
    
    if flip == 1:
        BrainMask_ct_nii = reorient_image("XXYY", registered_BrainMask_nii, output_dir)
    else:
        BrainMask_ct_nii = to_id(registered_BrainMask_nii)

    # Send the label to the native ct

   # Load the image
    ct_im = nibabel.load(ct_nii)
    ct_data = ct_im.get_data()
    BrainMask_data = nibabel.load(BrainMask_ct_nii).get_data()
    print 'BrainMask shape: ', BrainMask_data.shape
    print 'ct shape: ', ct_data.shape
    print 'global index = ', global_min_index
    BrainMask_to_ct_native = numpy.zeros(ct_data.shape)
    BrainMask_to_ct_native[:, :, :global_min_index] = 0
    BrainMask_to_ct_native[:, :, global_min_index:] = BrainMask_data
    BrainMask_ct_native_im = nibabel.Nifti1Image(BrainMask_to_ct_native, ct_im.get_affine())
    nibabel.save(BrainMask_ct_native_im, BrainMask_ct_native_nii)

    return BrainMask_ct_nii, BrainMask_ct_native_nii

def threed_dot(matrice, vector):
    """ Dot product between a 3d matrix and an image of 3d vectors.
    """
    res = numpy.zeros(vector.shape)
    for i in range(3):
      res[..., i] = (matrice[i, 0] * vector[..., 0] +
                       matrice[i, 1] * vector[..., 1] +
                       matrice[i, 2] * vector[..., 2] +
                       matrice[i, 3])
    return res

def inverse_affine(affine):
    """ Invert an affine transformation.
    """
    invr = numpy.linalg.inv(affine[:3, :3])
    inv_affine = numpy.zeros((4, 4))
    inv_affine[3, 3] = 1
    inv_affine[:3, :3] = invr
    inv_affine[:3, 3] =  - numpy.dot(invr, affine[:3, 3])
    return inv_affine


def BrainMask_to_rd(BrainMask_ct_native_nii, rd_nii, correct_cta, output_dir):
    """ Register the Brain%ask to rd space.
    """

    # Output autocompletion
    BrainMask_rescale_file = os.path.join(output_dir, "BrainMask_to_rd.nii.gz")

    # Load images
    BrainMask_ct_native_im = nibabel.load(BrainMask_ct_native_nii)
    BrainMask_data = BrainMask_ct_native_im.get_data()
    rd_im = nibabel.load(rd_nii)
    rd_data = rd_im.get_data()
    cta = BrainMask_ct_native_im.get_affine()
    rda = rd_im.get_affine()

    # Correct the rda affine matrix
   # cta[1, 3] += 21
    cta[2, 2] = correct_cta
    print cta


    # Inverse affine transformation
    icta = inverse_affine(cta)
    t = numpy.dot(icta, rda)
    numpy.savetxt(os.path.join(output_dir,"t_icta_rda.txt"), t)
    # Matricial dot product
    BrainMask_rescale = numpy.zeros(rd_data.shape)
    dot_image = numpy.zeros(rd_data.shape + (3, ))
    x = numpy.linspace(0, rd_data.shape[0] - 1, rd_data.shape[0])
    y = numpy.linspace(0, rd_data.shape[1] - 1, rd_data.shape[1])
    z = numpy.linspace(0, rd_data.shape[2] - 1, rd_data.shape[2])
    xg, yg, zg = numpy.meshgrid(y, x, z)
    print 'dot image shape: ', dot_image.shape
    print 'yg shape: ', yg.shape
    print 'xg shape: ', xg.shape
    print 'zg shape: ', zg.shape
    dot_image[..., 0] = yg
    dot_image[..., 1] = xg
    dot_image[..., 2] = zg
    dot_image = threed_dot(t, dot_image)

    cnt = 0
    print rd_data.size
    for x in range(rd_data.shape[0]):
        for y in range(rd_data.shape[1]):
            for z in range(rd_data.shape[2]):
            #for z in range(cut_brain_index, rd_data.shape[2]):
                if cnt % 100000 == 0:
                    print cnt
                cnt += 1
                voxel_BrainMask = dot_image[x, y, z]
                if (voxel_BrainMask > 0).all() and (voxel_BrainMask < (numpy.asarray(BrainMask_data.shape) - 1)).all():
                    voxel_BrainMask = numpy.round(voxel_BrainMask)
                    BrainMask_rescale[x, y, z] = BrainMask_data[voxel_BrainMask[0], voxel_BrainMask[1], voxel_BrainMask[2]]

    BrainMask_rescale_im = nibabel.Nifti1Image(BrainMask_rescale, rda)
    nibabel.save(BrainMask_rescale_im, BrainMask_rescale_file)

    return BrainMask_rescale_file    




if __name__ == "__main__":


    # Global parameters
    BASE_PATH = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence"
    ATLAS_PATH = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/atlas"
    nii_path = os.path.join(BASE_PATH, "sujet_18_rt")
    output_path = os.path.join(BASE_PATH, "results_from_label_to_rd_after_ants")
    subjects_csv = os.path.join(nii_path, "clinical_data.csv")


    # Keep the valid subject
    valid_subject_dirs = [os.path.join(nii_path, dir_name)
                      for dir_name in os.listdir(nii_path)
                      if os.path.isdir(os.path.join(nii_path, dir_name))]

    #valid_subject_dirs.sort()

    # Read the dataframe clinical data: age, orientation of the ct
    df_subjects = pd.read_csv(subjects_csv)



    # Go through all subjects
    for subject_path in valid_subject_dirs[1:2]:
        #subject_path = os.path.join(nii_path, 'sujet_024_VM')
        print "Processing: '{0}'...".format(subject_path)

        # Get subject id
        if not nii_path.endswith(os.path.sep):
            nii_path = nii_path + os.path.sep
        #subj_id = subject_path.replace(nii_path, "").split(os.path.sep)[0]
        subj_id = os.path.basename(subject_path)

        # Create output directory and skip processing if already
        output_dir = os.path.join(output_path, subj_id)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)


        #Check whether ct image has to be flipped
        flip = df_subjects.check_flip[df_subjects.anonym_nom == subj_id].values[0]

        # Get the t1, the inv_trans (from atlas to t1), the rd and the ct of the patient
        #t1_nii = glob.glob(os.path.join(subject_path, "mri", subj_id + '_T1.nii'))[0]
        #BrainMask_nii = glob.glob(os.path.join(subject_path, "brainMask", subj_id + "_BM.nii.gz"))[0]
        t1_nii = glob.glob(os.path.join(subject_path, "mri", subj_id + '_T1.nii'))[0]
        ct_nii = glob.glob(os.path.join(subject_path, "ct", subj_id + "_ct.nii.gz"))[0]
        rd_nii = glob.glob(os.path.join(subject_path, "rd", "*.nii"))[0]
        
        
        if flip == 1:
             ct_cut_reoriented_nii = os.path.join(output_dir, "ct_cut_brain_reorient.nii.gz")
        else:
             ct_cut_reoriented_nii = os.path.join(output_dir, "ct_cut_brain.nii.gz")
        print "NO flip"
        
        if os.path.isfile(os.path.join(output_dir, 'ct_cut_brain.nii.gz')):
            print("ok")
        elif os.path.isfile(os.path.join(output_dir, 'ct_cut_brain_reorient.nii.gz')):
            print("ok2")
        else:
            continue

        transformation= os.path.join(output_dir, "t1_to_cut_ct.txt")
        
        

        print "Executing: %s" % (ct_nii)
        print "Executing: %s" % (rd_nii)


        cut_brain_index_fileName = os.path.join(output_dir, "ct_brain_index.txt")
        cut_brain_index_file = open(cut_brain_index_fileName, "r")
        cut_brain_index = int(cut_brain_index_file.read())
        cut_brain_index_file.close()
        print "cut_brain_index: %s" % (cut_brain_index)
        
        # Correct the Rzz in the ct affine to find the correct correspondance in the physical coordonnates
        correct = {'sujet_005_BZ': 1.21, 'sujet_007_MM': 1.21,
                   'sujet_010_SA': 1.21, 'sujet_011_PA': 1,
                   'sujet_014_WS': 1, 'sujet_015_NI': 1.21,
                   'sujet_016_DG': 0.997, 'sujet_017_AG': 1.396,
                   'sujet_022_SK': 1, 'sujet_024_VM': 1.334,
                   'sujet_027_BL': 1.292, 'sujet_028_CH': 1,
                   'sujet_029_CT': 1.21,
                   'sujet_032_HB': 1.053, 'sujet_033_HL': 1.174,
                   'sujet_034_HI': 1, 'sujet_038_ZH': 1,
                   'sujet_012_OY':1.65}        
        
        BrainMask_nii = Extrac_brain(t1_nii, output_dir)        
        BrainMask_ct_nii = BrainMask_to_ct(cut_brain_index, ct_cut_reoriented_nii, transformation, BrainMask_nii, output_dir)
        
        BrainMask_ct_native_nii = os.path.join(output_dir,"BrainMask_to_ct_native.nii.gz" )
        if not os.path.isfile(os.path.join(output_dir,"BrainMask_to_ct_native.nii.gz" )):
            continue
        else:
		BrainMask_rd_nii = BrainMask_to_rd(BrainMask_ct_native_nii, rd_nii, correct[subj_id], output_dir)

                         