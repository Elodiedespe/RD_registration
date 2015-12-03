# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 16:19:16 2015

@author: cp243490
"""

import pandas as pd
import numpy as np
import os
import glob


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


if __name__ == "__main__":

    ROOT_PATH = '/volatile/pinaud/elodie_dosimetry'

    INPUT_PATH = os.path.join(ROOT_PATH, '02-Bundles-RD-Properties')

    ANALYSIS_PATH = os.path.join(ROOT_PATH, '03-Analysis')
    if not os.path.exists(ANALYSIS_PATH):
        os.mkdir(ANALYSIS_PATH)

    ####################
    ## Read Subjects csv
    INPUT_SUBJECTS_CSV = os.path.join(ROOT_PATH, 'subjects_data.csv')
    df_subjects = pd.read_csv(INPUT_SUBJECTS_CSV)
    subjects = df_subjects['anonym_nom'].values

    #########################
    ## extract and summarize statistics information from each subject and each
    ## bundle in a csv file:
    # Parameters properties associated to bundles
    # CSV summary
    columns = ['Segment', 'Bundle', 'Subject', 'Age', 'Traitement',
               'Hemisphere', 'section']
    measures = ['min', 'max', 'mean', 'std']
    statList = ['Mean-over-bundle', 'Profile-along-bundle']

    output_mean_fileName = os.path.join(ANALYSIS_PATH, 'bundles-mean.csv')
    output_profile_fileName = os.path.join(ANALYSIS_PATH,
                                           'bundles-profile.csv')

    #########################
    ## extract and summarize statistics information from each subject and each
    ## bundle in a csv file
    for stat in statList:
        print stat

        if stat == statList[0]:
            output_filename = output_mean_fileName
        elif stat == statList[1]:
            output_filename = output_profile_fileName

        df_stat = pd.DataFrame(columns=columns + measures)

        compt = 0
        for idx in df_subjects.index:
            cur = df_subjects.loc[idx]
            subject = cur.Subject
            chimio = cur.Traitement
            age = cur.Age
            print "subject: ", subject

            bundles = glob.glob(os.path.join(INPUT_PATH, stat, subject,
                                             '*', '*', '*.csv'))
            for bundle in bundles:
                bundleName = os.path.basename(bundle.split('.')[0])
                param = bundleName.split('_')[-1]
                bundleName = '_'.join(bundleName.split('_')[:-1])
                segment = bundle.split('/')[-2]
                bundleType = bundle.split('/')[-3]

                bundleSubjectInfo = get_bundleSubjectInfo(segment, bundleName)
                hemisphere = bundleSubjectInfo[0]
                bundleShortName = bundleSubjectInfo[1]
                segmentName = bundleSubjectInfo[2]

                df_bundle = pd.read_csv(bundle)
                for row in df_bundle.index:
                    cur = df_bundle.loc[row]
                    section = cur.section
                    stats = cur[measures]
                    df_stat.loc[compt] = [segmentName, bundleShortName,
                                          subject, age, chimio, hemisphere,
                                          section] + \
                                         stats.values.tolist()
                    compt += 1

        df_stat.to_csv(output_filename, index=False)
