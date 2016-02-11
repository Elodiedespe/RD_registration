from nipype.caching import Memory
from nipype.interfaces.ants import ApplyTransforms
import os
from labels import DKTprotocol
import pandas as pd
import glob
import nibabel as nib
import numpy as np
from nilearn.image import resample_img
from nilearn._utils.numpy_conversions import as_ndarray
from nilearn.image.resampling import coord_transform
import math
from scipy import ndimage
import numpy as np
from scipy.ndimage import label
from scipy.stats import norm
from nilearn.input_data import NiftiMasker

n_proc = 7


def from_native_to_mni(img, sub_id, include_trans=[True, True, True],
                       interpolation='Linear'):
    '''Maps image from native space to mni.

    WARNING THERE IS A CLEAR PROBLEM IN THE UNDERSTANDING OF TRANSFORM ORDER
    WHEN ONLY USING THE LAST TWO TRANSFORMS THE ORDER SHOULD BE INVERTED

    We assume that the transformation files already exist for the mappings
    between:
    1) mean bold and anatomy
    2) anatomy and oasis template
    3) oasis template and mni template

    The transforms to include are:
    1) From bold to anat
    2) From anat to oasis
    3) From oasis to mni

    The include transforms should be sequential to have meaninful output,
    which means that transformations sequence [True, False, True] is invalid.
    '''
    check = (include_trans == [True, False, True])
    if check:
        raise Exception('Invalid transformation sequence')
    pipeline_dir = 'pipelines/transformations'
    if not os.path.exists(pipeline_dir):
        os.makedirs(pipeline_dir)
    mem = Memory(pipeline_dir)
    transform = mem.cache(ApplyTransforms)

    anat = os.path.join('pipelines',
                        'preprocessing',
                        'sub{0}'.format(sub_id),
                        'highres001.nii')
    oasis_template = os.path.join('pipelines',
                                  'OASIS-30_Atropos_template',
                                  'T_template0.nii.gz')
    mni_template = os.path.join('pipelines',
                                'mni_icbm152_nlin_asym_09a_nifti',
                                'mni_icbm152_nlin_asym_09a',
                                'mni_icbm152_t1_tal_nlin_asym_09a.nii')
    bold_to_anat = os.path.join('pipelines', 'preprocessing',
                                'sub{0}'.format(sub_id),
                                'bold_to_anat.txt')
    anat_to_oasis = os.path.join('pipelines', 'preprocessing',
                                 'sub{0}'.format(sub_id),
                                 'anat_to_oasis.h5')
    oasis_to_mni = os.path.join('pipelines', 'preprocessing',
                                'registered_templates', 'oasis_to_mni.h5')
    all_references = [anat, oasis_template, mni_template]
    all_trans = [bold_to_anat, anat_to_oasis, oasis_to_mni]
    all_inv_trans = [False, False, False]
    transforms = []
    inv_trans_flags = []
    reference = None
    for idx, flag in enumerate(include_trans):
        if flag:
            transforms.append(all_trans[idx])
            inv_trans_flags.append(all_inv_trans[idx])
            # Use latest transformation as reference
            reference = all_references[idx]

    trans = transform(args='--float',
                      input_image_type=3,
                      interpolation=interpolation,
                      invert_transform_flags=inv_trans_flags[::-1],
                      num_threads=n_proc,
                      reference_image=reference,
                      terminal_output='file',
                      transforms=transforms[::-1],
                      input_image=img)

    return trans.outputs.output_image


def from_mni_to_native(img, sub_id, include_trans=[True, True, True],
                       interpolation='Linear'):
    '''Maps image from native space to mni.

    We assume that the transformation files already exist for the mappings
    between:
    1) mean bold and anatomy
    2) anatomy and oasis template
    3) oasis template and mni template

    The transforms to include are:
    1) From mni to oasis
    2) From oasis to anat
    3) From anat to bold

    The include transforms should be sequential to have meaninful output,
    which means that transformations sequence [True, False, True] is invalid.
    '''
    check = (include_trans == [True, False, True])
    if check:
        raise Exception('Invalid transformation sequence')
    pipeline_dir = 'pipelines/transformations'
    if not os.path.exists(pipeline_dir):
        os.makedirs(pipeline_dir)
    mem = Memory(pipeline_dir)
    transform = mem.cache(ApplyTransforms)

    oasis_template = os.path.join('pipelines',
                                  'OASIS-30_Atropos_template',
                                  'T_template0.nii.gz')
    anat = os.path.join('pipelines',
                        'preprocessing',
                        'sub{0}'.format(sub_id),
                        'highres001.nii')
    mean_bold = os.path.join('pipelines', 'preprocessing',
                             'sub{0}'.format(sub_id),
                             'mean_bold.nii')
    mni_to_oasis = os.path.join('pipelines', 'preprocessing',
                                'registered_templates', 'mni_to_oasis.h5')
    oasis_to_anat = os.path.join('pipelines', 'preprocessing',
                                 'sub{0}'.format(sub_id),
                                 'oasis_to_anat.h5')
    bold_to_anat = os.path.join('pipelines', 'preprocessing',
                                'sub{0}'.format(sub_id),
                                'bold_to_anat.txt')

    all_references = [oasis_template, anat, mean_bold]
    all_trans = [mni_to_oasis, oasis_to_anat, bold_to_anat]
    all_inv_trans = [False, False, True]
    transforms = []
    inv_trans_flags = []
    reference = None
    for idx, flag in enumerate(include_trans):
        if flag:
            transforms.append(all_trans[idx])
            inv_trans_flags.append(all_inv_trans[idx])
            # Use latest transformation as reference
            reference = all_references[idx]

    trans = transform(args='--float',
                      input_image_type=3,
                      interpolation=interpolation,
                      invert_transform_flags=inv_trans_flags[::-1],
                      num_threads=n_proc,
                      reference_image=reference,
                      terminal_output='file',
                      transforms=transforms[::-1],
                      input_image=img)

    return trans.outputs.output_image


def resample_roi_to_mni_template(img_path):
    pipeline_dir = 'pipelines/transformations/resampled'
    if not os.path.exists(pipeline_dir):
        os.makedirs(pipeline_dir)
    path_parts = os.path.split(img_path)
    new_img_path = os.path.join(pipeline_dir, path_parts[1])
    if os.path.isfile(new_img_path):
        return new_img_path

    mni_template = os.path.join('pipelines',
                                'mni_icbm152_nlin_asym_09a_nifti',
                                'mni_icbm152_nlin_asym_09a',
                                'mni_icbm152_t1_tal_nlin_asym_09a.nii')
    template_img = nib.load(mni_template)
    img = nib.load(img_path)
    r_img = resample_img(img, template_img.affine, template_img.shape,
                         'nearest')
    nib.save(r_img, new_img_path)
    return new_img_path


def resample_img_to_template(img_path, template_path, output_path,
                             interpolation='continuous'):
    pipeline_dir = 'pipelines/transformations/resampled'
    if not os.path.exists(os.path.join(pipeline_dir, output_path)):
        os.makedirs(os.path.join(pipeline_dir, output_path))
    path_parts = os.path.split(img_path)
    new_img_path = os.path.join(pipeline_dir, output_path, path_parts[1])
    if os.path.isfile(new_img_path):
        return new_img_path

    template_img = nib.load(template_path)
    img = nib.load(img_path)
    r_img = resample_img(img, template_img.affine, template_img.shape,
                         interpolation)
    nib.save(r_img, new_img_path)
    return new_img_path


def get_roi_mask(roi_filename, label_number=1):
    mask = nib.load(roi_filename).get_data()
    mask = mask.reshape(mask.shape[0:3])
    mask = mask == label_number
    return mask


def get_roi_niftiimage(roi_filename, label_number=1):
    mask = nib.load(roi_filename)
    affine = mask.get_affine()
    mask = mask.get_data()
    mask = mask.reshape(mask.shape[0:3])
    mask[mask != label_number] = 0
    mask[mask == label_number] = 1
    return nib.Nifti1Image(mask, affine)


def extract_data_in_mask(imgs, mask):
    """ Intersect the imgs with the specified ROI

    The current code assumes that mask and the imgs have the same shape!

    Parameters
    ----------
    imgs : is a list of 3D img files,
    roi_filename : is the name of an img file containing the ROI mask

    Returns
    -------
    data : array with imgs information in the given ROI
    """
    if type(imgs) is not list:
        imgs = [imgs]

    data = np.zeros([mask.shape[0], mask.shape[1], mask.shape[2], len(imgs)])
    for i in range(len(imgs)):
        assert os.path.isfile(imgs[i])
        img = nib.load(imgs[i]).get_data()
        img = img.reshape(img.shape[:3])
        assert img.shape == mask.shape
        data[:, :, :, i] = img

    return data[mask, :]

    # Alternative way of extracting data in ROI that is binary mask
    # roi = nib.load(roi_filename).get_data().astype(np.bool)
    # nvoxels = roi.sum()
    # X = np.zeros(shape=(nvoxels, ncon))
    # for icon, con in enumerate(conds_name):
    #     tmp = nib.load(imgs[icon]).get_data()
    #     this_mask = np.logical_and(roi, np.isfinite(tmp))
    #     X[:, icon] = tmp[this_mask]


def get_rois_mask_in_native_space(sub_id):
    """Gives the roi boolean numpy array and label.
    Also relevant dir and whole img dir in native space"""
    roi_dir = os.path.join('pipelines', 'ROIs')
    # Get all ROI file paths
    anat_roi_images = glob.glob(os.path.join(roi_dir, 'anatomical',
                                             '*.nii*'))
    func_roi_images = glob.glob(os.path.join(roi_dir, 'functional',
                                             '*', '*.nii*'))
    subj_roi_images = glob.glob(os.path.join(roi_dir, 'subject_specific',
                                             'sub{0}'.format(sub_id),
                                             '*', '*.nii*'))
    # Yield ROI for anat and func
    for img in (anat_roi_images + func_roi_images):
        path_parts = os.path.split(img)
        csv_file = os.path.join(path_parts[0],
                                path_parts[1].split('.')[0]+'.csv')
        img = resample_roi_to_mni_template(img)
        img_in_native_space = from_mni_to_native(img, sub_id,
                                                 interpolation='NearestNeighbor')
        if os.path.isfile(csv_file):
            labels = pd.read_csv(csv_file, index_col=False)
            for idx, label in labels.iterrows():
                # yield get_roi_mask(img, label['number']), label['name']
                yield (get_roi_mask(img_in_native_space,
                                    label['number']),
                       label['name'],
                       path_parts[0].split('/')[2:],
                       get_roi_niftiimage(img_in_native_space, label['number']),
                       get_roi_niftiimage(img, label['number']))
        else:
            # yield get_roi_mask(img), path_parts[1].split('.')[0]
            yield (get_roi_mask(img_in_native_space),
                   path_parts[1].split('.')[0],
                   path_parts[0].split('/')[2:],
                   img_in_native_space,
                   img)
    # Yield ROI for subject specific contrasts
    for img_in_native_space in subj_roi_images:
        path_parts = os.path.split(img_in_native_space)
        yield (get_roi_mask(img_in_native_space),
               path_parts[1].split('.')[0],
               ['subject_specific']+path_parts[0].split('/')[4:],
               img_in_native_space,
               from_native_to_mni(img_in_native_space, sub_id,
                                  interpolation='NearestNeighbor'))


def get_roi_center(roi_native_path, roi_mni_path):
    """Get ROI center of mass.
    Get back coordinate in img space and in coordinate space.
    Also actual center of mass.
    """
    # computations in native space
    if type(roi_native_path) is str:
        img = nib.load(roi_native_path)
    else:
        img = roi_native_path
    data = img.get_data()
    data = as_ndarray(data)
    my_map = data.copy()
    center_coords = ndimage.center_of_mass(np.abs(my_map))

    x_map, y_map, z_map = center_coords[:3]
    native_coords = np.asarray(coord_transform(x_map, y_map, z_map,
                                               img.get_affine())).tolist()
    voxel = [round(x) for x in center_coords]
    # computations in mni space
    if type(roi_mni_path) is str:
        img = nib.load(roi_mni_path)
    else:
        img = roi_mni_path
    data = img.get_data()
    data = as_ndarray(data)
    my_map = data.copy()
    mni_center_coords = ndimage.center_of_mass(np.abs(my_map))
    x_map, y_map, z_map = mni_center_coords[:3]
    mni_coords = np.asarray(coord_transform(x_map, y_map, z_map,
                                            img.get_affine())).tolist()
    # returns voxel and true center mass coords
    # returns also native and mni space coords
    return (voxel[:3], center_coords[:3], [round(x) for x in native_coords],
            [round(x) for x in mni_coords])


def fdr_threshold(z_vals, alpha):
    """ return the BH fdr for the input z_vals"""
    z_vals_ = - np.sort(- z_vals)
    p_vals = norm.sf(z_vals_)
    n_samples = len(p_vals)
    pos = p_vals < alpha * np.linspace(
        .5 / n_samples, 1 - .5 / n_samples, n_samples)
    if pos.any():
        return (z_vals_[pos][-1] - 1.e-8)
    else:
        return np.infty


def map_threshold(stat_img, mask_img, threshold, height_control='fpr',
                  cluster_threshold=0):
    """ Threshold the provvided map

    Parameters
    ----------
    stat_img : Niimg-like object,
       statistical image (presumably in z scale)

    mask_img : Niimg-like object,
        mask image

    threshold: float,
        cluster forming threshold (either a p-value or z-scale value)

    height_control: string
        false positive control meaning of cluster forming
        threshold: 'fpr'|'fdr'|'bonferroni'|'none'

    cluster_threshold : float, optional
        cluster size threshold

    Returns
    -------
    thresholded_map : Nifti1Image,
        the stat_map theresholded at the prescribed voxel- and cluster-level
    """
    # Masking
    masker = NiftiMasker(mask_img=mask_img)
    stats = np.ravel(masker.fit_transform(stat_img))
    n_voxels = np.size(stats)

    # Thresholding
    if height_control == 'fpr':
        z_th = norm.isf(threshold)
    elif height_control == 'fdr':
        z_th = fdr_threshold(stats, threshold)
    elif height_control == 'bonferroni':
        z_th = norm.isf(threshold / n_voxels)
    else:  # Brute-force thresholding
        z_th = threshold
    stats *= (stats > z_th)

    stat_map = masker.inverse_transform(stats).get_data()

    # Extract connected components above threshold
    label_map, n_labels = label(stat_map > z_th)
    labels = label_map[(masker.mask_img_.get_data() > 0)]
    for label_ in range(1, n_labels + 1):
        if np.sum(labels == label_) < cluster_threshold:
            stats[labels == label_] = 0

    return masker.inverse_transform(stats)


def create_rois_from_clusters(contrast_tmap, mask, threshold=3.09,
                              height_control='brute', cluster_threshold=10,
                              save_path=None):
    if save_path is not None:
        if not os.path.exists(save_path):
            os.makedirs(save_path)

    thresholded = map_threshold(contrast_tmap, mask, threshold,
                                height_control, cluster_threshold)
    cluster_map, n_cluster = label(thresholded.get_data() > 0)

    clusters = []
    masker = NiftiMasker(mask_img=mask)
    masker.fit()
    mask_affine = nib.load(mask).get_affine()
    for label_ in range(1, n_cluster + 1):
        cluster = cluster_map.copy()
        cluster[cluster_map != label_] = 0
        cluster[cluster_map == label_] = 1
        cluster = nib.Nifti1Image(cluster, mask_affine)
        clusters.append(cluster)
        if save_path is not None:
            nib.save(cluster, os.path.join(save_path,
                     'cluster_{0}.nii'.format(label_)))

    return clusters


if __name__ == "__main__":
    # Map Oasis anatomical labels to mni space and put file in ROIs
    oasis_roi_name = 'OASIS-TRT-20_jointfusion_DKT31_CMA_labels'
    # Add selection of labels to automatically extract and analyze
    roi_dir = os.path.join('pipelines', 'ROIs')
    dkt = DKTprotocol()
    roi_info = {}
    roi_info['number'] = (dkt.left_cerebrum_cortex_DKT31_numbers +
                          dkt.right_cerebrum_cortex_DKT31_numbers)
    roi_info['name'] = (dkt.left_cerebrum_cortex_DKT31_names +
                        dkt.right_cerebrum_cortex_DKT31_names)
    oasis_rois = pd.DataFrame(roi_info)
    oasis_rois.to_csv(os.path.join(roi_dir, 'anatomical',
                                   oasis_roi_name + '.csv'),
                      index=False)
