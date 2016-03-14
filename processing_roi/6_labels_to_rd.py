

# System import
import numpy
import os
import subprocess
import pandas as pd
import glob

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
    cta[1, 3] += 21
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
    for subject_path in valid_subject_dirs[7:8]:
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
        t1_nii = glob.glob(os.path.join(subject_path, "mri", "*", subj_id + '_T1.nii'))
        ct_nii = glob.glob(os.path.join(subject_path, "ct", "*", subj_id + "_ct.nii.gz"))
        rd_nii = glob.glob(os.path.join(subject_path, "rd", "*.nii"))[0]
        transform_warped = os.path.join(subject_path, "transform_warped", "transformInverseComposite.h5")

        print "Executing: %s" % (t1_nii)
        print "Executing: %s" % (ct_nii)
        print "Executing: %s" % (rd_nii)

        print "Executing: %s" % (transform_warped)


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

        labels_ct_native_nii = os.path.join(output_dir,"labels_to_ct_native.nii.gz" )
        if not os.path.isfile(os.path.join(output_dir,"labels_to_ct_native.nii.gz" )):
            continue
        else:
		labels_rd_nii = labels_to_rd(labels_ct_native_nii, rd_nii, correct[subj_id], output_dir)
