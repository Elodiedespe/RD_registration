# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 18:13:18 2015

@author: cp243490
"""


import os
import pandas as pd
import glob
import shutil
import pandas as pd
import numpy as np
import nibabel as nib


def to_nii(image, outdir):
    """
    Convert ima file to nii
    """
    image_nii = image.split('.')[0] + '.nii'
    image_nii = os.path.join(outdir, os.path.basename(image_nii))
    os.system('AimsFileConvert -i %s -o %s'
              % (image, image_nii))
    return image_nii


def create_centroid(inputBundle, outputCentroid):
    '''
    Compute teh centroid's bundle
    Centroid : average of fibers
    '''
    command = ' PtkDwiBundleOperator ' + \
              ' -i ' + inputBundle + \
              ' -o ' + outputCentroid + \
              ' -op compute-centroids ' + \
              ' -stringParameters ' + ' average-fiber ' + \
                                      ' all ' + \
                                      'symmetric-maximum-of-mean-closest-point-distance'
    #print "Executing: ", command
    os.system(command)


def get_bundleSubjectInfo(segment, bundleName):
    if segment == 'interhemispheric':
        hemisphere = 'interhemispheric'
        bundleShortName = bundleName
        segmentShortName = segment
    else:
        hemisphere = segment.split('-')[0]
        bundleShortName = bundleName[:- (len(hemisphere) + 1)]
        segmentShortName = segment[len(hemisphere) + 1:]

    return hemisphere, bundleShortName, segmentShortName


def weighted_mean(scalarVolume, weightVolume):
    '''
    Weighted mean
    map to be averaged: scalar parameter map
    weight in each voxel: bundle density in the voxel
    '''
    if not scalarVolume.split('.')[1] == 'nii':
        scalarVolume_nii = to_nii(scalarVolume, os.path.dirname(scalarVolume))
    else:
        scalarVolume_nii = scalarVolume
    if not weightVolume.split('.')[1] == 'nii':
        weightVolume_nii = to_nii(weightVolume, os.path.dirname(weightVolume))
    else:
        weightVolume_nii = weightVolume

    scalarVolume_im = nib.load(scalarVolume_nii).get_data()
    weightVolume_im = nib.load(weightVolume_nii).get_data()
    intersectScalarWeight = scalarVolume_im[weightVolume_im > 0][scalarVolume_im[weightVolume_im > 0] > 0]
    if np.max(weightVolume_im > 0):
        mean = np.average(scalarVolume_im, weights=weightVolume_im)
    else:
        mean = None
    if len(intersectScalarWeight) != 0:
        minVal = np.min(intersectScalarWeight)
        maxVal = np.max(intersectScalarWeight)
        stdVal = np.std(intersectScalarWeight)
    else:
        minVal = None
        maxVal = None
        stdVal = None

    stats = [minVal, maxVal, mean, stdVal]

    if not scalarVolume == scalarVolume_nii:
        os.remove(scalarVolume_nii)
        os.remove(scalarVolume_nii + '.minf')
    if not weightVolume == weightVolume_nii:
        os.remove(weightVolume_nii)
        os.remove(weightVolume_nii + '.minf')

    return stats


def createDensityMask(bundleFile, scalarVolume, tmpdir, output_fileName):
    """
    create density mask of a bundle
    """
    scalarVolume_S16 = '%s_S16.%s' % (scalarVolume.split('.')[0],
                                      scalarVolume.split('.')[1])
    scalarVolume_S16 = os.path.join(tmpdir,
                                    os.path.basename(scalarVolume_S16))
    os.system('AimsFileConvert -i %s -o %s -t S16'
              % (scalarVolume, scalarVolume_S16))
    os.system('PtkDwiBundleMapDensityMask -b %s -t2 %s -o %s'
              % (bundleFile, scalarVolume_S16, output_fileName))
    os.remove(scalarVolume_S16)
    os.remove(scalarVolume_S16 + '.minf')
    return output_fileName


def createSectionMask(bundleFile, centroidFile, width, output_sectionMask):
    """
    create section mask of a bundle
    """
    command = "PtkDwiBundleMapSectionMask " + \
              " -b " + bundleFile + \
              " -c " + centroidFile + \
              " -os " + output_sectionMask + \
              " -od None " + \
              " -w " + str(width)
    #print "Executing: ", command
    os.system(command)


def resampleBundle(bundleFile, trmFile, output_dir):
    '''
    Bundle resampling in a volume referential
    '''
    bundleResampled = "%s_resampled.%s" % (bundleFile.split('.')[0],
                                           bundleFile.split('.')[1])
    bundleResampled = os.path.join(output_dir,
                                   os.path.basename(bundleResampled))
    # bundle Resampling
    command = "ConnectBundleAnalysis " + \
              " -i " + bundleFile + \
              " -t " + trmFile + \
              " -b " + bundleResampled
    #print "Executing: ", command
    os.system(command)
    return bundleResampled


def maskVolume(volumeFile, maskFile, outdir):
    """
    Mask a volume
    """
    maskedVolumeFile = "%s_threshold.%s" % (volumeFile.split('.')[0],
                                            volumeFile.split('.')[1])
    maskedVolumeFile = os.path.join(outdir,
                                    os.path.basename(maskedVolumeFile))
    command = "AimsMask " + \
              " -i " + volumeFile + \
              " -o " + maskedVolumeFile + \
              " -m " + maskFile
    #print "Executing: ", command
    os.system(command)
    return maskedVolumeFile


def binarizeVolume(volumeFile, mode, thresh1, thresh2, outFile):
    """
    Create a mask of a thresholded volume
    """
    if thresh1 is None:
        outFile = volumeFile
    elif thresh2 is not None:
        command = "PtkBinarizer " + \
                  " -i " + volumeFile + \
                  " -o " + outFile + \
                  " -t int16_t " + \
                  " -m " + mode + \
                  " -t1 " + str(thresh1) + \
                  " -t2 " + str(thresh2)
        #print "Executing: ", command
        os.system(command)
    else:
        command = "PtkBinarizer " + \
                  " -i " + volumeFile + \
                  " -o " + outFile + \
                  " -t int16_t " + \
                  " -m " + mode + \
                  " -t1 " + str(thresh1)
        #print "Executing: ", command
        os.system(command)
    return outFile


def thresholdVolume(volumeFile, threshold, outdir):
    """
    threshold a volume
    """
    thresholdVolumeFile = "%s_threshold.%s" % (volumeFile.split('.')[0],
                                               volumeFile.split('.')[1])
    thresholdVolumeFile = os.path.join(outdir,
                                       os.path.basename(thresholdVolumeFile))
    command = "PtkThresholder " + \
              " -i " + volumeFile + \
              " -o " + thresholdVolumeFile + \
              " -m bt " + \
              " -t1 0 " + \
              " -t2 " + str(threshold)
    #print "Executing: ", command
    os.system(command)
    return thresholdVolumeFile


def fiber_statistics(bundle, output_fileName):
    """
    For the bundle, create a file that gives
    its number of fibers, its length ans its tortuosity
    """
    command = "PtkDwiBundleMeasure " + \
              " -b " + bundle + \
              " -o " + output_fileName + \
              " -m " + " bundle_fiber_count " + \
                       " bundle_fiber_length_statistics " + \
                       " bundle_fiber_tortuosity_statistics " + \
              " -scalarParameters 1 1 1 "
    #print "Executing: ", command
    os.system(command)


def statistics(scalarVolume, densityMask, sectionMask, bundleName, output_dir):
    '''
    Mean value over a bundle of a scalar parameter map
    diffusion maps and relaxometry maps
    '''
    stats = {}
    if os.path.isfile(scalarVolume):
        if sectionMask is not None:
            # Get labels of the sectionsMask
            sectionMask_nii = to_nii(sectionMask, output_dir)
        else :
            sectionMask_nii = os.path.join(output_dir,
                                           '%_sectionMask.nii' % bundleName)
            sectionMask_nii = binarizeVolume(densityMask, 'gt',
                                             0, None, sectionMask_nii)
        sectionMask_im = nib.load(sectionMask_nii)
        sectionMask_arr = sectionMask_im.get_data()
        sections = np.unique(sectionMask_arr[sectionMask_arr > 0])

        for section in sections:
            sectionMask_label = os.path.join(output_dir, "%s_label%s.%s"
                                             % (os.path.basename(sectionMask.split('.')[0]),
                                                str(section),
                                                sectionMask.split('.')[1]))
            sectionMask_label = binarizeVolume(sectionMask, 'eq',
                                               section, None, sectionMask_label)
            densityMask_label = os.path.join(output_dir, "%s_label%s.%s"
                                             % (os.path.basename(densityMask.split('.')[0]),
                                                str(section),
                                                densityMask.split('.')[1]))
            densityMask_label = maskVolume(densityMask, sectionMask_label,
                                           output_dir)
            stats[section] = weighted_mean(scalarVolume, densityMask_label)

            fileTrashList =  glob.glob(sectionMask_label.split('.')[0] + '.*') + \
                             glob.glob(densityMask_label.split('.')[0] + '.*')
            for fileTrash in fileTrashList:
                os.remove(fileTrash)

        if sectionMask_nii != sectionMask:
            os.remove(sectionMask_nii)
            os.remove(sectionMask_nii + '.minf')
    else:
        print "no file " + scalarVolume

    return stats


def stats_rd(scalarVolume, bundleFile, bundleName,
                    densityMask, sectionMask, centroidFile, sectionWidth,
                    output_dir):
    '''
    Mean value over a bundle of a scalar parameter map
    diffusion maps and relaxometry maps
    '''
    if os.path.isfile(scalarVolume):
        # create density mask of the bundle
        createDensityMask(bundleFile, scalarVolume, output_dir, densityMask)
        if not os.path.exists(sectionMask):
            createSectionMask(bundleFile, centroidFile, sectionWidth,
                              sectionMask)
        # Mean over bundle
        mean = statistics(scalarVolume, densityMask, sectionMask,
                          bundleName, output_dir)

        # Mean over bundle
        profile = statistics(scalarVolume, densityMask, sectionMask,
                             bundleName, output_dir)

    else:
        print "no file " + scalarVolume
        mean = {}
        profile = {}

    return mean, profile


if __name__ == "__main__":

    ROOT_PATH = '/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/elodie_dosimetry'
    TRANSFO_PATH = os.path.join(ROOT_PATH, '01-TO_CT')
    DATA_PATH = os.path.join(ROOT_PATH, 'subjects_data')
    ATLAS_PATH = os.path.join(ROOT_PATH, 'Atlas_Pediatric_2Hemispheres')

    OUTPUT_PATH = os.path.join(ROOT_PATH, '02-Bundles-RD-Properties')
    if not os.path.exists(OUTPUT_PATH):
        os.mkdir(OUTPUT_PATH)

    ####################
    ## Read Subjects csv
    INPUT_SUBJECTS_CSV = os.path.join(ROOT_PATH, 'clinical_data.csv')
    df_subjects = pd.read_csv(INPUT_SUBJECTS_CSV)
    subjects = df_subjects['anonym_nom'].values

    ####################
    ##  OUtPUT csv file
    columns = ['Segment', 'Bundle', 'Subject', 'Age', 'Traitement',
               'Hemisphere', 'section'] + ['min', 'max', 'mean', 'std']

    # When density map is sectionned (profile along bundle):
    # sectionSize= length of 1 section
    sectionSize = 4  # mm

    ## MRI parameters Profile along bundle
    ################################################
    print "* MRI parameters Profile along bundle"
    print "*" * 20

    DENSITYMASK_PATH = os.path.join(OUTPUT_PATH, 'DensityMask')

    OUTPUT_PROFILE_PATH = os.path.join(OUTPUT_PATH, 'Profile-along-bundle')
    OUTPUT_MEAN_PATH = os.path.join(OUTPUT_PATH, 'Mean-over-bundle')

    if not os.path.exists(DENSITYMASK_PATH):
        os.mkdir(DENSITYMASK_PATH)
    if not os.path.exists(OUTPUT_PROFILE_PATH):
        os.mkdir(OUTPUT_PROFILE_PATH)
    if not os.path.exists(OUTPUT_MEAN_PATH):
        os.mkdir(OUTPUT_MEAN_PATH)

    for idx in df_subjects[df_subjects.RT > 0].index:
        cur = df_subjects.loc[idx]
        subject = cur.anonym_nom
        age = cur.age_today
        chimio = cur.Traitement
        print "subject: ", subject

        SUBJECT_DATA_PATH = os.path.join(DATA_PATH, subject)
        if not os.path.exists(SUBJECT_DATA_PATH):
            continue
        rdFile = os.path.join(SUBJECT_DATA_PATH,
                              '%s_rd.nii.gz' % subject)
        template_to_ct_trm = os.path.join(TRANSFO_PATH, subject,
                                          'Template_To_CT.trm')

        bundles = glob.glob(os.path.join(ATLAS_PATH, '*bundles', 'bundles', '*',
                                         '*.bundles'))

        for bundle in bundles:
            bundleName = os.path.basename(bundle.split('.')[0])[len('atlas') + 1:]
            segment = bundle.split('/')[-2]
            bundleType = bundle.split('/')[-4][len('atlas') + 1:]
            print bundleName

            bundleSubjectInfo = get_bundleSubjectInfo(segment, bundleName)
            hemisphere = bundleSubjectInfo[0]
            bundleShortName = bundleSubjectInfo[1]
            segmentName = bundleSubjectInfo[2]

            DENSITYMASK_SUBJ_PATH = os.path.join(DENSITYMASK_PATH, subject,
                                                 bundleType, segment)
            for i, dir in enumerate(DENSITYMASK_SUBJ_PATH.split('/')[-4:-1]):
                if not os.path.exists('/'.join(DENSITYMASK_SUBJ_PATH.split('/')[:-3+i])):
                    os.mkdir('/'.join(DENSITYMASK_SUBJ_PATH.split('/')[:-3+i]))
            if not os.path.exists(DENSITYMASK_SUBJ_PATH):
                os.mkdir(DENSITYMASK_SUBJ_PATH)


            OUTPUT_SUBJ_MEAN_PATH = os.path.join(OUTPUT_MEAN_PATH, subject,
                                                 bundleType, segment) 
            for i, dir in enumerate( OUTPUT_SUBJ_MEAN_PATH.split('/')[-4:-1]):
                if not os.path.exists('/'.join( OUTPUT_SUBJ_MEAN_PATH.split('/')[:-3+i])):
                    os.mkdir('/'.join( OUTPUT_SUBJ_MEAN_PATH.split('/')[:-3+i]))
            if not os.path.exists( OUTPUT_SUBJ_MEAN_PATH):
                os.mkdir( OUTPUT_SUBJ_MEAN_PATH)
                
            OUTPUT_SUBJ_PROFILE_PATH = os.path.join(OUTPUT_PROFILE_PATH, subject,
                                                 bundleType, segment)
          
            for i, dir in enumerate(OUTPUT_SUBJ_PROFILE_PATH.split('/')[-4:-1]):
                if not os.path.exists('/'.join(OUTPUT_SUBJ_PROFILE_PATH.split('/')[:-3+i])):
                    os.mkdir('/'.join(OUTPUT_SUBJ_PROFILE_PATH.split('/')[:-3+i]))
            if not os.path.exists(OUTPUT_SUBJ_PROFILE_PATH):
                os.mkdir(OUTPUT_SUBJ_PROFILE_PATH)

            TMP_PATH = os.path.join(OUTPUT_MEAN_PATH, subject,
                                    bundleType, 'tmp')
            if not os.path.exists(TMP_PATH):
                os.mkdir(TMP_PATH)

            # Resample bundle into native referential
            bundleFile_subj = resampleBundle(bundle, template_to_ct_trm,
                                             TMP_PATH)
            bundleCentroidFile_subj = os.path.join(TMP_PATH,
                                                   '%s_centroid.bundles'
                                                   % bundleName)
            if not os.path.exists(bundleCentroidFile_subj):
                create_centroid(bundleFile_subj, bundleCentroidFile_subj)

            # density_mask
            densityMask_subj = os.path.join(DENSITYMASK_SUBJ_PATH,
                                            bundleName + '_densityMask.ima')
            sectionMask_subj = os.path.join(DENSITYMASK_SUBJ_PATH,
                                            bundleName + '_sectionMask.ima')

            if os.path.isfile(rdFile):
                # mean value over bundle
                mean, profile = stats_rd(rdFile, bundleFile_subj, bundleName,
                                         densityMask_subj, sectionMask_subj,
                                         bundleCentroidFile_subj, sectionSize,
                                         TMP_PATH)
                df_mean = pd.DataFrame(columns=columns)
                df_mean.loc[0] = [segmentName, bundleShortName,
                                  subject, age, chimio, hemisphere, 1] + \
                                 mean[1]
                df_mean.to_csv(os.path.join(OUTPUT_SUBJ_MEAN_PATH, '%s_RD.csv'
                                            % bundleName),
                            index=False)

                df_profile = pd.DataFrame(columns=columns)
                for i, item in enumerate(profile.iteritems()):
                    section = item[0]
                    stats = item[1]
                    distance = section * sectionSize
                    df_profile.loc[i] = [segmentName, bundleShortName, subject,
                                         age, chimio, hemisphere, distance] + \
                                        stats
                df_profile.to_csv(os.path.join(OUTPUT_SUBJ_PROFILE_PATH, '%s_RD.csv'
                                            % bundleName),
                                  index=False)

            else:
                print "no file " + rdFile

        trashFiles = [os.path.join(TMP_PATH, f) for f in os.listdir(TMP_PATH)
                      if bundleName in f] + \
                     [os.path.join(TMP_PATH, f) for f in os.listdir(TMP_PATH)
                      if bundleName in f]
        for trashFile in trashFiles:
            if os.path.isfile(trashFile):
                os.remove(trashFile)

    shutil.rmtree(TMP_PATH)
