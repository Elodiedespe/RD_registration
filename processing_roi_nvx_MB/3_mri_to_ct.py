"""Cette partie dy script est dans mri_to_ct_registration10032015.py"""
# System import
import numpy
import os
import subprocess
import scipy.signal
import glob


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


def mri_to_ct( ct_nii, min_thr, output_dir, verbose=0):
    """ Register the mri t1 scan to the ct image.
    """
    # Output autocompletion
    ct_modify_nii = os.path.join(output_dir, "ct_modify.nii.gz")
    ct_brain_nii = os.path.join(output_dir, "ct_cut_brain.nii.gz")
    ct_cut_reoriented_nii = os.path.join(output_dir, "ct_cut_reoriented.nii.gz")
    cut_brain_index_fileName = os.path.join(output_dir, "ct_brain_index.txt")

    # Load ct and modify the data for brain extraction
    ct_im = nibabel.load(ct_nii)
    ct_data = ct_im.get_data()
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

    cut_brain_index_file = open(cut_brain_index_fileName, "w")
    cut_brain_index_file.write(str(global_min_index))
    cut_brain_index_file.close()

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
    
    ct_cut_reoriented_nii = reorient_image("LPS", ct_brain_nii, output_dir)

    return ct_cut_reoriented_nii, global_min_index

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
    for subject_path in valid_subject_dirs[5:6]:
        #subject_path = os.path.join(nii_path, 'sujet_024_VM')
        print "Processing: '{0}'...".format(subject_path)

        # Get subject id
        if not nii_path.endswith(os.path.sep):
            nii_path = nii_path + os.path.sep
        #subj_id = subject_path.replace(nii_path, "").split(os.path.sep)[0]
        subj_id = os.path.basename(subject_path)
        print subj_id
        
        
        # Select the correct atlas according to age
        if subj_id == "sujet_044_DS":
        	continue

        else:

            # Create output directory and skip processing if already
            output_dir = os.path.join(output_path, subj_id)
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)

            # Get the t1, the inv_trans (from atlas to t1), the rd and the ct of the patient
            ct_nii = glob.glob(os.path.join(subject_path, "ct", subj_id + "_ct.nii.gz"))[0]

        
            ct_cut_reoriented_nii, global_min_index = mri_to_ct(ct_nii, 50000, output_dir, verbose=1)

