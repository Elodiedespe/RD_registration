# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 11:57:49 2015

@author: edogerde
"""

""" Register and resample the RD over the CT. Rd rescale will allow us to 
get the dose profile and mean of the WM tract with bundle properties.py and 
01_TO_CT.py """

# System import
import numpy
import os

# IO import
import nibabel

def inverse_affine(affine):
    """ Invert an affine transformation.
    """
    invr = numpy.linalg.inv(affine[:3, :3])
    inv_affine = numpy.zeros((4, 4))
    inv_affine[3, 3] = 1
    inv_affine[:3, :3] = invr
    inv_affine[:3, 3] =  - numpy.dot(invr, affine[:3, 3])
    return inv_affine


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

def rd_to_ct(ct_nii, rd_nii, cut_brain_index, output_dir):
    """ Register the rd to the ct space.
    """
    # Output autocompletion
    rd_rescale_file = os.path.join(output_dir, "rd_rescale.nii.gz")
    rd_rescale_brain_file = os.path.join(output_dir, "rd_rescale_brain.nii.gz")
    # Load images
    ct_im = nibabel.load(ct_nii)
    ct_data = ct_im.get_data()
    rd_im = nibabel.load(rd_nii)
    rd_data = rd_im.get_data()
    cta = ct_im.get_affine()
    rda = rd_im.get_affine()
    print "ok 1"
    # Correct the rda affine matrix
    rda[2, 2] = 3
    # Inverse affine transformation
    irda = inverse_affine(rda)
    t = numpy.dot(irda, cta)
    print "ok2"
    # Matricial dot product
    rd_rescale = numpy.zeros(ct_data.shape)
    dot_image = numpy.zeros(ct_data.shape + (3, ))
    x = numpy.linspace(0, ct_data.shape[0] - 1, ct_data.shape[0])
    y = numpy.linspace(0, ct_data.shape[1] - 1, ct_data.shape[1])
    z = numpy.linspace(0, ct_data.shape[2] - 1, ct_data.shape[2])
    xg, yg, zg = numpy.meshgrid(y, x, z)
    dot_image[..., 0] = yg
    dot_image[..., 1] = xg
    dot_image[..., 2] = zg
    dot_image = threed_dot(t, dot_image)  
    cnt = 0
    print ct_data.size
    for x in range(ct_data.shape[0]):
        for y in range(ct_data.shape[1]):
            for z in range(cut_brain_index, ct_data.shape[2]):
                if cnt % 100000 == 0:
                    print cnt  
                cnt += 1          
                voxel_rd = dot_image[x, y, z]
                if (voxel_rd > 0).all() and (voxel_rd < (numpy.asarray(rd_data.shape) - 1)).all():
                    rd_voxel = numpy.round(voxel_rd)
                    rd_rescale[x, y, z] = rd_data[rd_voxel[0], rd_voxel[1], rd_voxel[2]]
    #Create a nifti object 
    rd_rescale_im = nibabel.Nifti1Image(rd_rescale, cta)
    nibabel.save(rd_rescale_im, rd_rescale_file)
    rd_rescale_brain = numpy.zeros((ct_data.shape[0], ct_data.shape[1],
                                   ct_data.shape[2] - cut_brain_index))
    rd_rescale_brain = rd_rescale[:, :, cut_brain_index: ct_data.shape[2]]
    rd_rescale_brain_im = nibabel.Nifti1Image(rd_rescale_brain, cta)
    nibabel.save(rd_rescale_brain_im, rd_rescale_brain_file)
    
    return rd_rescale_file
   
if __name__ == "__main__": 
   
    # Global parameters
    BASE_PATH = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence"
    nii_path = os.path.join(BASE_PATH, "sujet_18_rt")
    SUBJ_ID = ["sujet_010_SA", "sujet_011_PA", "sujet_012_OY", "sujet_017_AG",
               "sujet_038_ZH", "sujet_014_WS" ]
    Index = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/results_from_label_to_rd_after_ants"
   
    # Find the appropriate ct, rd and ct_cut_brain_index data                      
    for subj_id in SUBJ_ID:             
    
        ct_nii = os.path.join(nii_path, subj_id, "ct", subj_id + "_ct.nii.gz")
        rd_nii = os.path.join(nii_path, subj_id, "rd", subj_id + ".nii")
        cut_brain_index_fileName = os.path.join(Index ,subj_id, "ct_brain_index.txt")
        cut_brain_index_file = open(cut_brain_index_fileName, "r")
        cut_brain_index = int(cut_brain_index_file.read())
        cut_brain_index_file.close()
        print "cut_brain_index: %s" % (cut_brain_index)
        print "Executing: %s" % (ct_nii)
        print "Executing: %s" % (rd_nii)
        
        # Create a output directory- create the output file where the output 
        #of labels_to_ct.py where
        output_dir = os.path.join(BASE_PATH,"ct_labels_processing", subj_id)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        
        # Correct the Rzz in the ct affine to find the correct correspondance 
        #in the physical coordonnates
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
        
        # Call the function
        rd_to_ct(ct_nii, rd_nii, cut_brain_index, output_dir)