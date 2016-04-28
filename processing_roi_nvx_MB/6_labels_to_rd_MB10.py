

# System import
import numpy
import os
import pandas as pd
import glob
import subprocess
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
    
    
def label_to_ctspace(ct_nii, label_to_ct_nii, global_min_index, output_dir):
    # Apply the affine and warp transformation

    labels_ct_reo_nii =  os.path.join(output_dir, "labels_to_ct_reor.nii.gz")
    labels_ct_native_nii = os.path.join(output_dir, "labels_to_ct_native.nii.gz")    
    
    labels_ct_reo_nii = reorient_image("XXYY", label_to_ct_nii, output_dir)    
    ct_im = nibabel.load(ct_nii)
    ct_data = ct_im.get_data()
    labels_data = nibabel.load(labels_ct_reo_nii).get_data()
    print'label shape: ', labels_data.shape
    print 'ct shape: ', ct_data.shape
    print 'global index = ', global_min_index
    labels_to_ct_native = numpy.zeros(ct_data.shape)
    labels_to_ct_native[:, :, :global_min_index] = 0
    labels_to_ct_native[:, :, global_min_index:] = labels_data
    labels_ct_native_im = nibabel.Nifti1Image(labels_to_ct_native, ct_im.get_affine())
    nibabel.save(labels_ct_native_im, labels_ct_native_nii)
    

    return labels_ct_native_nii    

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


def labels_to_rd(labels_ct_native_nii, rd_nii, correct_cta, output_dir):
    """ Register the rd to the ct space.
    """

    # Output autocompletion
    labels_rescale_file = os.path.join(output_dir, "labels_to_rd.nii.gz")

    # Load images
    labels_ct_native_im = nibabel.load(labels_ct_native_nii)
    labels_data = labels_ct_native_im.get_data()
    rd_im = nibabel.load(rd_nii)
    rd_data = rd_im.get_data()
    cta = labels_ct_native_im.get_affine()
    rda = rd_im.get_affine()

    # Correct the rda affine matrix
    print cta	
    #cta[1, 3] += 21
    cta[2, 2] = correct_cta
    print cta


    # Inverse affine transformation
    icta = inverse_affine(cta)
    t = numpy.dot(icta, rda)
    numpy.savetxt(os.path.join(output_dir,"t_icta_rda.txt"), t)
    # Matricial dot product
    labels_rescale = numpy.zeros(rd_data.shape)
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
                voxel_labels = dot_image[x, y, z]
                if (voxel_labels > 0).all() and (voxel_labels < (numpy.asarray(labels_data.shape) - 1)).all():
                    voxel_labels = numpy.round(voxel_labels)
                    labels_rescale[x, y, z] = labels_data[voxel_labels[0], voxel_labels[1], voxel_labels[2]]

    labels_rescale_im = nibabel.Nifti1Image(labels_rescale, rda)
    nibabel.save(labels_rescale_im, labels_rescale_file)

    return labels_rescale_file


if __name__ == "__main__":


    # Global parameters
    BASE_PATH = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/MB_10"
    nii_path = os.path.join(BASE_PATH, "sujet_10_rt")
    output_path = os.path.join(BASE_PATH, "results")
    subjects_csv = os.path.join(BASE_PATH, "flipfileNMB.csv")


    # Keep the valid subject
    valid_subject_dirs = [os.path.join(nii_path, dir_name)
                      for dir_name in os.listdir(nii_path)
                      if os.path.isdir(os.path.join(nii_path, dir_name))]


    # Read the dataframe clinical data: age, orientation of the ct
    df_subjects = pd.read_csv(subjects_csv)

    

    # Go through all subjects
    for subject_path in valid_subject_dirs[4:5]:
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


        #Check whether ct image has to be flipped
        flip = df_subjects.flip[df_subjects.anonym_nom == subj_id].values[0]

        # Get the ct, the label_to_ct, the rd and brain_cut_index the patient
        label_to_ct_nii = glob.glob(os.path.join(output_dir, "atlas", "*.nii.gz"))[0]
        # reorient the registered label 
    
            
        ct_nii = os.path.join(subject_path, "ct", subj_id + "_ct.nii.gz")
        rd_nii = os.path.join(subject_path, "rd", "rd_resample.nii")

        print "Executing: %s" % (ct_nii)
        print "Executing: %s" % (rd_nii)
        print "Executing: %s" % (label_to_ct_nii)
        
        cut_brain_index_fileName = os.path.join(output_dir, "ct_brain_index.txt")
        cut_brain_index_file = open(cut_brain_index_fileName, "r")
        cut_brain_index = int(cut_brain_index_file.read())
        cut_brain_index_file.close()

        # Correct the Rzz in the ct affine to find the correct correspondance in the physical coordonnates
        correct = {'sujet_027_BL': 1.274, 'sujet_040_DK': 1.293,
                   'sujet_041_CA': 1.25, 'sujet_042_CO': 1,
                   'sujet_043_MI': 1.25, 'sujet_044_DS': 1,
                   'sujet_046_EN': 1.299, 'sujet_047_RS': 1,}
                   
        
        labels_ct_native_nii = label_to_ctspace(ct_nii, label_to_ct_nii, cut_brain_index, output_dir)
        labels_rd_nii = labels_to_rd(labels_ct_native_nii, rd_nii, correct[subj_id], output_dir)

        
