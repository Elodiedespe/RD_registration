# -*- coding: utf-8 -*-

from nipype.interfaces.ants.registration import Registration
from nipype.interfaces.ants.resampling import ApplyTransforms
from nipype.caching import Memory                  # caching mechanism
import os
import glob
import pandas as pd

n_proc = 7


def save_file(file_to_save, save_to, target_name=False):
    print 'Processing ' + str(file_to_save)
    # If file_to_save is a list then call save_file recursively
    if type(file_to_save) is list:
        [save_file(f, save_to, target_name) for f in file_to_save]
        return
    # Determine new file name and create required directories
    if target_name:
        new_file_name = save_to
    else:
        new_file_name = os.path.join(save_to,
                                     os.path.split(file_to_save)[1])
    if not os.path.exists(os.path.dirname(new_file_name)):
        os.makedirs(os.path.dirname(new_file_name))

    if os.path.isfile(new_file_name):
        os.unlink(new_file_name)
    os.link(file_to_save, new_file_name)


def anat_preproc(file_to_register, register_to, warp_back, pipeline_dir):
    # DATA CONFIGURATION. FOLLOWING OPENFMRI STANDARD.
    save_to = os.path.join(pipeline_dir,
                           file_to_register.split('/')[-1].split('.')[0])
    # Run pipeline imperatively with caching (without workflow object)
    mem = Memory(pipeline_dir)
    antsreg = mem.cache(Registration)
    transform = mem.cache(ApplyTransforms)
    save_list = []
    # nodes manual parameter configuration and run
    reg = antsreg(args='--float',
                  collapse_output_transforms=True,
                  moving_image=file_to_register,
                  fixed_image=register_to,
                  initial_moving_transform_com=True,
                  num_threads=n_proc,
                  output_inverse_warped_image=True,
                  output_warped_image=True,
                  sigma_units=['vox']*3,
                  transforms=['Rigid', 'Affine', 'SyN'],
                  terminal_output='file',
                  winsorize_lower_quantile=0.005,
                  winsorize_upper_quantile=0.995,
                  convergence_threshold=[1e-06],
                  convergence_window_size=[10],
                  metric=['MI', 'MI', 'CC'],
                  metric_weight=[1.0]*3,
                  number_of_iterations=[[1000, 500, 250, 100],
                                        [1000, 500, 250, 100],
                                        [100, 70, 50, 20]],
                  radius_or_number_of_bins=[32, 32, 4],
                  sampling_percentage=[0.25, 0.25, 1],
                  sampling_strategy=['Regular',
                                     'Regular',
                                     'None'],
                  shrink_factors=[[8, 4, 2, 1]]*3,
                  smoothing_sigmas=[[3, 2, 1, 0]]*3,
                  transform_parameters=[(0.1,),
                                        (0.1,),
                                        (0.1, 3.0, 0.0)],
                  use_histogram_matching=True,
                  write_composite_transform=True)
    save_list.append([reg.outputs.composite_transform, save_to])
    save_list.append([reg.outputs.warped_image, save_to])
    save_list.append([reg.outputs.inverse_composite_transform, save_to])
    save_list.append([reg.outputs.inverse_warped_image, save_to])
    transformed = transform(args='--float',
                            input_image_type=3,
                            interpolation='NearestNeighbor',
                            invert_transform_flags=[False],
                            num_threads=n_proc,
                            reference_image=file_to_register,
                            terminal_output='file',
                            transforms=reg.outputs.inverse_composite_transform,
                            input_image=warp_back)
    save_list.append([transformed.outputs.output_image, save_to])
    return save_list


if __name__ == "__main__":

    # MAIN PATH
    ROOT_PATH = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas"
    ATLAS_PATH = os.path.join(ROOT_PATH,"atlas")

    # get the appropriate richard_template
    template_2 = os.path.join(ATLAS_PATH, "atlas_0-2/ANTS2-0Years_head.nii.gz")
    labels_2 = os.path.join(ATLAS_PATH, "atlas_0-2/ANTS2-0Years_label_regis_head.nii.gz")
    template_5 = os.path.join(ATLAS_PATH, "atlas_5_9/ANTS9-5Years3T_head.nii.gz")
    labels_5 = os.path.join(ATLAS_PATH, "atlas_5_9/ANTS9-5Years3T_label_regis_head.nii.gz")
    template_2_5 = os.path.join(ATLAS_PATH, "atlas_2_5/ANTS2-5Years_head.nii.gz")
    labels_2_5 = os.path.join(ATLAS_PATH, "atlas_2_5/ANTS2-5Years_label_regis_head.nii.gz")

    # get the t1
    files = glob.glob(os.path.join(ROOT_PATH, 'children_mri', 'sujet*',
                                   'mri', '*', '*.nii'))
    files.sort()
    results = []
    pipeline_dir = os.path.join(ROOT_PATH, 'registered')
    if not os.path.exists(pipeline_dir):
        os.makedirs(pipeline_dir)

    for f in files:
        # Find the age of the subject
        subjects_csv = os.path.join(ROOT_PATH, "clinical_data.csv")
        df_subjects = pd.read_csv(subjects_csv)

        subj_id = f.split('/')[-4]
        subj_age = df_subjects.ART[df_subjects.anonym_nom == subj_id].values[0]
        print "subject: ", subj_id
        print "age: ", subj_age
        print '\n'

        if subj_age < 2:
            template = template_2
            template_labels = labels_2
            print " under 2 years old"
        elif subj_age > 5:
            template = template_5
            template_labels = labels_5
            print " over 5 years old"
        else:
            template = template_2_5
            template_labels = labels_2_5
            print " between 2 and 5 years old"

        results.append(anat_preproc(f, template, template_labels,
                       pipeline_dir))

    # Save the results
    for r in results:
        for s in r:
            save_file(*s)
