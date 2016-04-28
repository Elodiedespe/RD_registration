# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 12:38:34 2015

@author: edogerde
"""

import numpy
import os

# IO import
import nibabel


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

    
def ct_to_rd(ct_nii, rd_nii, correct_cta, output_dir):
    """ Register the rd to the ct space.
    """
    
    # Output autocompletion
    ct_rescale_file = os.path.join(output_dir, "ct_rescale.nii.gz")

    # Load images
    ct_im = nibabel.load(ct_nii)
    ct_data = ct_im.get_data()
    rd_im = nibabel.load(rd_nii)
    rd_data = rd_im.get_data()
    cta = ct_im.get_affine()
    rda = rd_im.get_affine()

    # Correct the rda affine matrix
    
    cta[2, 2] = correct_cta
    #rda[2, 2] = 3

    # Inverse affine transformation
    icta = inverse_affine(cta)
    t = numpy.dot(icta, rda)

    # Matricial dot product
    ct_rescale = numpy.zeros(rd_data.shape)
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
                voxel_ct = dot_image[x, y, z]
                if (voxel_ct > 0).all() and (voxel_ct < (numpy.asarray(ct_data.shape) - 1)).all():
                    ct_voxel = numpy.round(voxel_ct)
                    ct_rescale[x, y, z] = ct_data[ct_voxel[0], ct_voxel[1], ct_voxel[2]]

    ct_rescale_im = nibabel.Nifti1Image(ct_rescale, rda)
    nibabel.save(ct_rescale_im, ct_rescale_file)
    
  
    

    return ct_rescale_file

if __name__ == "__main__":

    # Global parameters
    BASE_PATH = "/home/edogerde/Bureau/MB_10"
    nii_path = os.path.join(BASE_PATH, "sujet_10_rt")
    output_path = os.path.join(BASE_PATH, "results")

    # Keep the valid subject
    valid_subject_dirs = [os.path.join(nii_path, dir_name)
                      for dir_name in os.listdir(nii_path)
                      if os.path.isdir(os.path.join(nii_path, dir_name))]
    
    #valid_subject_dirs.sort()

    # Go through all subjects
    for subject_path in valid_subject_dirs[1:5]:
        #subject_path = os.path.join(nii_path, 'sujet_024_VM')
        print "Processing: '{0}'...".format(subject_path)

        # Get subject id
        if not nii_path.endswith(os.path.sep):
            nii_path = nii_path + os.path.sep
        #subj_id = subject_path.replace(nii_path, "").split(os.path.sep)[0]
        subj_id = os.path.basename(subject_path)
        print subj_id
        

        # Create output directory and skip processing if already
        output_dir = os.path.join(output_path, subj_id)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        # Correct the Rzz in the ct affine to find the correct correspondance in the physical coordonnates
        correct = {'sujet_027_BL': 1.274, 'sujet_040_DK': 1.293,
                    'sujet_041_CA': 1.25, 'sujet_042_CO': 1,
                    'sujet_043_MI': 1.25, 'sujet_044_DS': 1,
                    'sujet_046_EN': 1.299, 'sujet_047_RS': 1,}
        # Get the t1, the inv_trans (from atlas to t1), the rd and the ct of the patient
        ct_nii = os.path.join(subject_path, "ct", subj_id + "_ct.nii.gz")
        rd_nii = os.path.join(subject_path, "rd", "rd_resample.nii")

        ct_rescale_file = ct_to_rd(ct_nii, rd_nii, correct[subj_id], output_dir)