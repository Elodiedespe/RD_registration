
# System import
import numpy
import os
import subprocess
import scipy.signal
import glob
import shutil

# Plot import
import matplotlib.pyplot as plt

# IO import
import nibabel

POSSIBLE_AXES_ORIENTATIONS = [
    "LAI", "LIA", "ALI", "AIL", "ILA", "IAL",
    "LAS", "LSA", "ALS", "ASL", "SLA", "SAL",
    "LPI", "LIP", "PLI", "PIL", "ILP", "IPL",
    "LPS", "LSP", "PLS", "PSL", "SLP", "SPL",
    "RAI", "RIA", "ARI", "AIR", "IRA", "IAR",
    "RAS", "RSA", "ARS", "ASR", "SRA", "SAR",
    "RPI", "RIP", "PRI", "PIR", "IRP", "IPR",
    "RPS", "RSP", "PRS", "PSR", "SRP", "SPR"]

CORRECTION_MATRIX_COLUMNS = {
    "R": (1, 0, 0),
    "L": (-1, 0, 0),
    "A": (0, 1, 0),
    "P": (0, -1, 0),
    "S": (0, 0, 1),
    "I": (0, 0, -1)
}


def mri_to_template(t1_nii, atlas_nii, output_dir):
    """ Register the template to the t1 subject reference.
    """
    # Output autocompletion
    register_t1_aff_nii = os.path.join(output_dir, "t1_to_atlas_aff.nii.gz")
    trans_aff = os.path.join(output_dir, "t1_to_atlas_aff.txt")
    register_t1_nl_nii = os.path.join(output_dir, "t1_to_atlas_nl.nii.gz")
    trans_nl = os.path.join(output_dir, "t1_to_atlas_nl_field.nii.gz")

    # Affine registration
    cmd = ["flirt", "-cost", "normmi", "-omat", trans_aff, "-in", t1_nii,
           "-ref", atlas_nii, "-out", register_t1_aff_nii]
    print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)

    # NL registration
    cmd = ["fnirt", "--ref={0}".format(atlas_nii), "--in={0}".format(t1_nii),
           "--iout={0}".format(register_t1_nl_nii),
           "--fout={0}".format(trans_nl), "--aff={0}".format(trans_aff),
           "--config=T1_2_MNI152_2mm.cnf"]
    print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)

    return trans_aff, trans_nl


def swap_affine(axes):
    """ Build a correction matrix, from the given orientation
    of axes to RAS.
    """
    rotation = numpy.eye(4)
    rotation[:3, 0] = CORRECTION_MATRIX_COLUMNS[axes[0]]
    rotation[:3, 1] = CORRECTION_MATRIX_COLUMNS[axes[1]]
    rotation[:3, 2] = CORRECTION_MATRIX_COLUMNS[axes[2]]
    return rotation


def reorient_image(input_axes, in_file, output_dir):
    """ Rectify the orientation of an image.
    """
    # get the transformation to the RAS space
    rotation = swap_affine(input_axes)
    det = numpy.linalg.det(rotation)
    if det != 1:
        raise Exception("Determinant must be equal to "
                        "one got: {0}.".format(det))

    # load image
    image = nibabel.load(in_file)

    # get affine transform (qform or sform)
    affine = image.get_affine()

    # apply transformation
    transformation = numpy.dot(rotation, affine)
    image.set_qform(transformation)
    image.set_sform(transformation)

    # save result
    reoriented_file = os.path.join(output_dir, "im_reorient.nii.gz")
    nibabel.save(image, reoriented_file)

    return reoriented_file


def mri_to_ct(t1_nii, ct_nii, min_thr, output_dir, verbose=0):
    """ Register the mri t1 scan to the ct image.
    """
    # Output autocompletion
    ct_modify_nii = os.path.join(output_dir, "ct_modify.nii.gz")
    ct_brain_nii = os.path.join(output_dir, "ct_cut_brain.nii.gz")
    register_t1_nii = os.path.join(output_dir, "t1_to_cut_ct.nii.gz")
    transformation = os.path.join(output_dir, "t1_to_cut_ct.txt")
    t1_ct_nii = os.path.join(output_dir, "t1_to_ct.nii.gz")

    # Load ct and modify the data for brain extraction
    ct_im = nibabel.load(ct_nii)
    ct_data = ct_im.get_data()
    ct_shape = ct_data.shape
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

    # Diplay if verbose mode
    if verbose == 1:
        x = range(power.shape[0])
        plt.plot(x, power, '.', linewidth=1)
        plt.plot(x, powerfilter, '--', linewidth=1)
        plt.plot(x[global_min_index], powerfilter[global_min_index], "o")
        plt.show()

    # Cut the image
    ct_cut_data = ct_data[:, :, range(global_min_index, ct_data.shape[2])]
    brain_im = nibabel.Nifti1Image(ct_cut_data, ct_im.get_affine())
    nibabel.save(brain_im, ct_brain_nii)

    # Reorient ct brain image
    reorient_image("LPS", ct_brain_nii, output_dir)

    # Register
    cmd = ["flirt", "-cost", "normmi", "-omat", transformation, "-in", t1_nii,
           "-ref", ct_brain_nii, "-out", register_t1_nii]
    print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)

    # Send the t1 to the ct original space
    t1_data = nibabel.load(register_t1_nii).get_data()
    ct_t1_data = numpy.zeros(ct_data.shape)
    ct_t1_data[:, :, range(global_min_index, ct_data.shape[2])] = t1_data
    t1_ct_im = nibabel.Nifti1Image(ct_t1_data, ct_im.get_affine())
    nibabel.save(t1_ct_im, t1_ct_nii)

    return transformation, ct_brain_nii, global_min_index


def labels_to_ct(t1_trans, t2_trans, t1_nii, labels_nii, ct_brain_nii, ct_nii,
                 output_dir):
    """ Register the labels to the CT.
    """
    # Output autocompletion
    invt1_trans = os.path.join(output_dir, "inv_t1_to_atlas_nl_field.nii.gz")
    combined_trans = os.path.join(output_dir, "combined_atlas_to_ct_field.nii.gz")
    registered_labels_nii = os.path.join(output_dir, "labels_to_cut_ct.nii.gz")
    labels_ct_nii = os.path.join(output_dir, "labels_to_ct.nii.gz")

    # Invert the nl warp
    cmd = ["invwarp", "-w", t1_trans, "-o", invt1_trans, "-r", t1_nii]
    print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)

    # Combine transformation
    cmd = ["convertwarp", "--ref={0}".format(ct_brain_nii),
           "--postmat={0}".format(t2_trans), "--warp1={0}".format(invt1_trans),
           "--out={0}".format(combined_trans)]
    print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)

    # Warp the labels
    cmd = ["applywarp", "-i", labels_nii, "-o", registered_labels_nii,
           "-r", ct_brain_nii, "-w", combined_trans, "--interp=nn"]
    print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)

    # Send the labels to the ct original space
    ct_im = nibabel.load(ct_nii)
    ct_data = ct_im.get_data()
    labels_data = nibabel.load(registered_labels_nii).get_data()
    ct_labels_data = numpy.zeros(ct_data.shape)
    ct_labels_data[:, :, -labels_data.shape[2]:] = labels_data
    labels_ct_im = nibabel.Nifti1Image(ct_labels_data, ct_im.get_affine())
    nibabel.save(labels_ct_im, labels_ct_nii)

    return labels_ct_nii


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

    # Load images
    ct_im = nibabel.load(ct_nii)
    ct_data = ct_im.get_data()
    rd_im = nibabel.load(rd_nii)
    rd_data = rd_im.get_data()
    cta = ct_im.get_affine()
    rda = rd_im.get_affine()

    # Correct the rda affine matrix
    rda[2, 2] = 3

    # Inverse affine transformation
    irda = inverse_affine(rda)
    t = numpy.dot(irda, cta)

    # Matricial dot product
    rd_rescale = numpy.zeros(ct_data.shape)
    dot_image = numpy.zeros(ct_data.shape + (3, ))
    x = numpy.linspace(0, ct_data.shape[0] - 1, ct_data.shape[0])
    y = numpy.linspace(0, ct_data.shape[1] - 1, ct_data.shape[1])
    z = numpy.linspace(0, ct_data.shape[2] - 1, ct_data.shape[2])
    xg, yg, zg = numpy.meshgrid(x, y, z)
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

    rd_rescale_im = nibabel.Nifti1Image(rd_rescale, cta)
    nibabel.save(rd_rescale_im, rd_rescale_file)

    return rd_rescale_file




if __name__ == "__main__":

    # Global parameters
    nii_path = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/data_set_n_7_nifti"
    output_path = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/ct_labels_processing"
    atlas_nii = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/atlas/atlas_t1.nii.gz"
    labels_nii = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/atlas/atlas_labels.nii.gz"

    # Keep only valid folder
    valid_subject_dirs = [os.path.join(nii_path, dir_name)
                          for dir_name in os.listdir(nii_path)
                          if os.path.isdir(os.path.join(nii_path, dir_name))]

    # Go through all subjects
    for subject_path in valid_subject_dirs:

        print "Processing: '{0}'...".format(subject_path)

        # Get subject id
        if not nii_path.endswith(os.path.sep):
            nii_path = nii_path + os.path.sep
        subj_id = subject_path.replace(nii_path, "").split(os.path.sep)[0]

        # Get the t1
        t1_nii = glob.glob(os.path.join(subject_path, "mri", "mri_bravo", "*", "*.nii*"))
        if len(t1_nii) > 0:
            t1_nii = t1_nii[0]
        else:
            t1_nii = glob.glob(os.path.join(subject_path, "mri", "mri_sagcube", "*", "*.nii*"))
            if len(t1_nii) > 0:
                t1_nii = t1_nii[0]
            else:
                print "Can't process subject '{0}', no t1 found.".format(subject_path)
                continue

        # Get all the ct
        ct_niis = glob.glob(os.path.join(subject_path, "ct", "*", "*.nii.gz"))

        # Go though all ct
        for ct_nii in ct_niis:

            print "    Processing: '{0}'...".format(ct_nii)

            # Find the associated rd
            rd_niis = glob.glob(
                os.path.join(os.path.dirname(ct_nii).replace("ct", "rd"),
                             "*.nii.gz"))
            if len(rd_niis) != 1:
                raise ValueError(
                    "Can't associate a rd file with ct: '{0}'.".format(rd_niis))
            rd_nii = rd_niis[0]

            # Get the ct identifier
            ct_id = ct_nii.split(os.path.sep)[-2]

            # Build output directory
            output_dir = os.path.join(output_path, subj_id, ct_id)

            # Create output directory and skip processing if already
            # created
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)
            else:
                print "-" * 20
                print "already processed."
                print "-" * 20
                continue

            # Execute ina try except in order to clean folder when a processing
            # error occured
            print "-" * 20
            try:
                t1_trans_aff, t1_trans_nl = mri_to_template(
                    t1_nii, atlas_nii, output_dir)
                t2_trans, ct_brain_nii, cut_brain_index = mri_to_ct(
                    t1_nii, ct_nii, 50000, output_dir, verbose=0)
                labels_ct_nii = labels_to_ct(
                    t1_trans_nl, t2_trans, t1_nii, labels_nii, ct_brain_nii,
                    ct_nii, output_dir)
                rd_rescale_file = rd_to_ct(
                    ct_nii, rd_nii, cut_brain_index, output_dir)
                print "Result in: '{0}' - '{1}'.".format(labels_ct_nii,
                                                         rd_rescale_file)
            except:
                shutil.rmtree(output_dir)
                raise
            print "-" * 20

