import roi_manager as rm
import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import matplotlib
import statsmodels.formula.api as sm
import argparse
from joblib import Parallel, delayed
import nitime.timeseries as ts
import nitime.analysis as nta
import nitime.viz as viz
from nilearn import input_data
import glob
import nibabel as nib
from nipype.interfaces.base import Bunch
import itertools
from scipy import signal
import math
import scipy.stats as sst

n_proc = 1
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
font = {#'family': 'arial',
        'weight': 'bold',
        'size': 12}
matplotlib.rc('font', **font)

dataset_name = 'openfmri_dataset'
OPENFMRI_SEP = '\t'


def subject_model(sub_id, task_id, model_id, runs, concat):
    '''Extracts model from openfmri dataset for given task and subject'''
    onsets_template = os.path.join(dataset_name, 'sub{0}', 'model', 'model{1}',
                                   'onsets', 'task{2}_run{3}', '{4}.txt')
    param_template = os.path.join(dataset_name, 'sub{0}', 'model', 'model{1}',
                                  'parameters', 'task{2}_run{3}', '{4}.txt')
    condkey_template = os.path.join(dataset_name, 'models', 'model{0}',
                                    'condition_key.txt')
    conditions_df = pd.read_csv(condkey_template.format(model_id),
                                sep=OPENFMRI_SEP,
                                header=None)
    task_conditions = conditions_df[conditions_df[0] == 'task'+task_id]
    conds = task_conditions[1].tolist()
    conds_name = task_conditions[2].tolist()
    subject_info = []
    all_event_files = []
    ec = 1
    for run in runs:
        run_id = run[run.index('run')+3:run.index('run')+6]
        names = []
        allonsets = []
        alldurations = []
        event_files = []
        # ASSUMING SAME NUMBER OF REPETITION FOR ALL RUNS (should change)
        # This will be used to add dummy events.
        vol = nib.load(os.path.join(dataset_name, 'sub{0}', 'BOLD',
                                    'task{1}_run001', 'bold.nii.gz')
                       .format(sub_id, task_id)).shape[3]
        TR = float(open(os.path.join(dataset_name,
                                     'scan_key.txt')).read().split('\t')[1])
        pmod = None
        # configure parametric modulation as necessary
        have_param = len(glob.glob(param_template.format(sub_id, model_id,
                                                         task_id, run_id,
                                                         '*'))) > 0
        if have_param:
            pmod = []
        # Build the subject_info object
        for i in range(len(conds)):
            names.append(conds_name[i])
            cond_file = onsets_template.format(sub_id, model_id,
                                               task_id, run_id, conds[i])
            event_files.append(cond_file)
            onsets = []
            durations = []
            # ULTRA IMPORTANT: WE NEED AT LEAST TWO EVENTS FOR ANY CONDITION
            # THERE IS A WEIRD BUG IN RELEASE VERSION OF SPECIFYSPMMODEL
            # THIS BUG IS NOW ADDRESSED IN THE NIPYPE DEVELOPER RELEASE
            if os.stat(cond_file).st_size > 0:
                cond_info = pd.read_csv(cond_file,
                                        sep=OPENFMRI_SEP,
                                        header=None)
                onsets = cond_info[0].tolist()
                durations = cond_info[1].tolist()
                # Uncomment these lines if working with nipype Release 0.10.0
                # if len(onsets) < 2:
                #     onsets.append(vol*TR - 0.1)
                #     durations.append(0.0)
            elif not concat or len(runs) == 1:
                print 'empty file found: ' + cond_file
                onsets = [vol*(TR-1) - 0.06*ec]
                durations = [0.05]
                ec += 1
                # Uncomment these lines if working with nipype Release 0.10.0
                # onsets = [vol*TR - 0.1, vol*TR - 0.2]
                # durations = [0.0, 0.0]
            allonsets.append(onsets)
            alldurations.append(durations)
            n_on = len(list(itertools.chain.from_iterable(allonsets)))
            n_dur = len(list(itertools.chain.from_iterable(alldurations)))
            assert(n_on == n_dur)
            # check existence of parametric regressors
            if have_param:
                pmod_name = []
                pmod_poly = []
                pmod_param = []
                p_files = glob.glob(param_template.format(sub_id, model_id,
                                                          task_id, run_id,
                                                          conds[i]+'_*'))
                p_files.sort()
                for p_file in p_files:
                    name = os.path.basename(p_file).split('.')[0].split('_')[2]
                    param_info = pd.read_csv(p_file,
                                             sep=OPENFMRI_SEP,
                                             header=None)
                    # if param_info[0].tolist():
                    pmod_name.append(name)
                    pmod_poly.append(1)
                    pmod_param.append(param_info[0].tolist())
                if p_files:
                    pmod.append(Bunch(name=pmod_name, poly=pmod_poly,
                                      param=pmod_param))
                else:
                    pmod.append(None)

        subject_info.append(Bunch(conditions=names,
                                  onsets=allonsets,
                                  durations=alldurations,
                                  amplitudes=None,
                                  tmod=None,
                                  pmod=pmod,
                                  regressor_names=None,
                                  regressors=None))
        all_event_files.append(event_files)
    # print sub_id, task_id

    return subject_info, all_event_files


def roi_histogram(mask, mask_name, imgs, results_path, title_extra=''):
    roi_data = []
    for idx, img in enumerate(imgs):
        roi_data.append(rm.extract_data_in_mask(img, mask))
    roi_data = np.concatenate(roi_data)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    fig.suptitle(mask_name+'\n'+title_extra)
    sns.distplot(roi_data, ax=ax)
    fig.savefig(os.path.join(results_path, mask_name+'.png'))
    plt.close()


def roi_barplot(mask, mask_name, contrasts, labels, results_path, title_extra=''):
    roi_data = pd.DataFrame(np.zeros(shape=(mask.sum(), len(contrasts))),
                            columns=labels)
    for idx, img in enumerate([contrast_path.format(c) for c in contrasts]):
        roi_data[labels[idx]] = rm.extract_data_in_mask(img, mask)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    fig.suptitle(mask_name+'\n'+title_extra)
    sns.barplot(data=roi_data, ax=ax)
    fig.savefig(os.path.join(results_path, mask_name+'.png'))
    plt.close()


def roi_statistical_model(mask, mask_name, contrasts, labels, factors,
                          results_path, model_desc):
    # factors should be a pandas dataframe with named columns
    # factors correspond to labels by idx
    roi_data = pd.DataFrame(np.zeros(shape=(mask.sum()*len(contrasts),
                                            len(factors.columns)+3)),
                            columns=['condition']+factors.columns.tolist() + \
                                    ['voxel', 'beta']
                           )
    counter = 0
    print 'preparing data in ' + mask_name
    for idx, img in enumerate([contrast_path.format(c) for c in contrasts]):
        data_cat = [labels[idx]] + factors.loc[idx].tolist()
        for v_idx, beta in enumerate(rm.extract_data_in_mask(img, mask)):
            row = data_cat + [v_idx, beta[0]]
            roi_data.loc[counter] = row
            counter += 1
    print 'computing model: ' + model_desc + ' in ' + mask_name
    regress = sm.ols(model_desc, data=roi_data).fit()
    f = open(os.path.join(results_path, mask_name+'.txt'), 'w')
    f.write(regress.summary().as_text())
    f.close()


def coord_event_related_analysis(coord, raw_bold, clean_bold, onset_events, TR,
                                 len_et, results_path, mask_name):
    # Prepare the event related analyzer from nitime
    t0 = ts.TimeSeries(raw_bold, sampling_interval=TR, time_unit='s')
    t1 = ts.TimeSeries(clean_bold, sampling_interval=TR, time_unit='s')
    # t2 = ts.Events(onset_events, time_unit='s')
    t2 = ts.TimeSeries(onset_events, sampling_interval=TR, time_unit='s')
    E = nta.EventRelatedAnalyzer(t1, t2, len_et)

    # Create figures
    fig = viz.plot_tseries(t0, ylabel='BOLD (raw) '+str(coord))
    fig.savefig(os.path.join(results_path, mask_name + '_rawbold.png'))
    plt.close()
    fig = viz.plot_tseries(t1, ylabel='BOLD (clean) '+str(coord))
    fig.savefig(os.path.join(results_path, mask_name + '_cleanbold.png'))
    plt.close()

    fig = viz.plot_tseries(E.eta, ylabel='BOLD (% signal change) '+str(coord),
                           yerror=E.ets)
    fig.savefig(os.path.join(results_path, mask_name + '_eta.png'))
    plt.close()

    fig = viz.plot_tseries(E.FIR, ylabel='BOLD (% signal change) '+str(coord))
    fig.savefig(os.path.join(results_path, mask_name + '_fir.png'))
    plt.close()

    # fig = viz.plot_tseries(E.xcorr_eta, ylabel='BOLD (% signal change) '+str(coord))
    # fig.savefig(os.path.join(results_path, mask_name + '_xcorr.png'))
    # plt.close()


def cluster_extraction(subject, model, task, contrasts, tz=3.09, tc=10):
    sub_roi_path = os.path.join('pipelines', 'ROIs', 'subject_specific',
                                subject)
    brain_mask = os.path.join('pipelines', 'preprocessing', subject,
                              'highres001_brain_mask_bold_native.nii')
    glm_path = os.path.join('pipelines', 'glm_sub_analysis', subject,
                            model, task, 'basic',
                            'allruns_concat_noderivatives')
    contrast_path = os.path.join(glm_path, 'con_{0:04d}.nii')
    contrast_desc_file = os.path.join(glm_path, 'contrasts.csv')
    contrast_df = pd.read_csv(contrast_desc_file)
    for c in contrasts:
        lookup = contrast_df['nii'] == 'con_{0:04d}.nii'.format(c)
        c_name = contrast_df[lookup]['name'].tolist()[0]
        c_folder = task + '_' + c_name
        rm.create_rois_from_clusters(contrast_path.format(c),
                                     brain_mask,
                                     save_path=os.path.join(sub_roi_path,
                                                            c_folder),
                                     threshold=tz,
                                     cluster_threshold=tc)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run ROI analysis')
    parser.add_argument('-sub',
                        type=str,
                        nargs=1,
                        default='001',
                        help='subject to consider')
    args = parser.parse_args()
    subject = 'sub'+args.sub[0]
    sub_id = args.sub[0]

    preproc_pipeline_dir = 'pipelines/preprocessing'
    run_file_template = os.path.join(preproc_pipeline_dir,
                                     'sub{0}',
                                     'task{1}',
                                     'run{2}_sr_bold.nii')

    # # Extract clusters
    # model = 'model001'
    # task = 'task002'
    # contrasts = [3]
    # cluster_extraction(subject, model, task, contrasts)

    # task = 'task003'
    # contrasts = [4, 6]
    # cluster_extraction(subject, model, task, contrasts)

    # task = 'task004'
    # contrasts = [18, 19, 22, 23, 24, 25, 30, 31, 32, 33]
    # cluster_extraction(subject, model, task, contrasts)

    # TSNR distribution in ROIS
    # template_results_path = os.path.join('pipelines', 'roi_analysis', subject,
    #                                      'tsnr', 'task{0}')
    # tsnr_path = 'pipelines/preprocessing/sub{0}/task{1}/tsnr_s_run*.nii'
    # for roi in rm.get_rois_mask_in_native_space(sub_id):
    #     for task in ['001', '002', '003', '004']:
    #         imgs = glob.glob(tsnr_path.format(sub_id, task))
    #         results_path = template_results_path.format(task)
    #         if not os.path.exists(os.path.join(results_path, *roi[2])):
    #             os.makedirs(os.path.join(results_path, *roi[2]))
    #         roi_histogram(roi[0], roi[1], imgs, os.path.join(results_path, *roi[2]))

    # # Searchlight distributions
    # # variance test
    # temp_results_path = os.path.join('pipelines', 'roi_analysis', subject,
    #                                  'var_ratio_test')
    # img_path = 'pipelines/searchlight/sub{0}/ratio_b_over_nb.nii'.format(sub_id)
    # for roi in rm.get_rois_mask_in_native_space(sub_id):
    #     imgs = [img_path]
    #     results_path = temp_results_path
    #     if not os.path.exists(os.path.join(results_path, *roi[2])):
    #         os.makedirs(os.path.join(results_path, *roi[2]))
    #     roi_histogram(roi[0], roi[1], imgs, os.path.join(results_path, *roi[2]))

    # # category prediction
    # temp_results_path = os.path.join('pipelines', 'roi_analysis', subject,
    #                                  'cat_pred')
    # img_path = 'pipelines/searchlight/sub{0}/number_cat_pred.nii'.format(sub_id)
    # for roi in rm.get_rois_mask_in_native_space(sub_id):
    #     imgs = [img_path]
    #     results_path = temp_results_path
    #     if not os.path.exists(os.path.join(results_path, *roi[2])):
    #         os.makedirs(os.path.join(results_path, *roi[2]))
    #     roi_histogram(roi[0], roi[1], imgs, os.path.join(results_path, *roi[2]))

    # # quantity prediction
    # temp_results_path = os.path.join('pipelines', 'roi_analysis', subject,
    #                                  'quant_pred')
    # img_path = 'pipelines/searchlight/sub{0}/number_quant_pred.nii'.format(sub_id)
    # for roi in rm.get_rois_mask_in_native_space(sub_id):
    #     imgs = [img_path]
    #     results_path = temp_results_path
    #     if not os.path.exists(os.path.join(results_path, *roi[2])):
    #         os.makedirs(os.path.join(results_path, *roi[2]))
    #     roi_histogram(roi[0], roi[1], imgs, os.path.join(results_path, *roi[2]))

    # # Language localizer
    # model = 'model001'
    # task = 'task002'
    # results_path = os.path.join('pipelines', 'roi_analysis', subject, model,
    #                             task, 'basic')
    # if not os.path.exists(results_path):
    #     os.makedirs(results_path)

    # labels = ['Sentence', 'Consonants']
    # glm_path = os.path.join('pipelines', 'glm_sub_analysis', subject,
    #                         model, task, 'basic',
    #                         'allruns_concat_noderivatives')
    # contrast_path = os.path.join(glm_path, 'con_{0:04d}.nii')
    # contrasts = [1, 2]

    # assert(len(contrasts) == len(labels))
    # for roi in rm.get_rois_mask_in_native_space(sub_id):
    #     if not os.path.exists(os.path.join(results_path, *roi[2])):
    #         os.makedirs(os.path.join(results_path, *roi[2]))
    # #     roi_barplot(roi[0], roi[1], contrasts, labels,
    # #                 os.path.join(results_path, *roi[2]),
    # #                 str([math.ceil(x) for x in rm.get_roi_center(roi[3], roi[4])[3]]))
    #     model_id = '001'
    #     task_id = '002'
    #     TR = 1.5
    #     runs = glob.glob(run_file_template.format(sub_id, task_id, '*'))
    #     runs.sort()
    #     subject_info, _ = subject_model(sub_id, task_id, model_id, runs, True)
    #     time_series = []
    #     time_series_events = []

    #     allcoords = rm.get_roi_center(roi[3], roi[4])
    #     xc, yc, zc = allcoords[0]
    #     #coords = allcoords[2]
    #     #print allcoords
    #     #masker = input_data.NiftiSpheresMasker([tuple(coords)], radius=1.)
    #     for idx, run in enumerate(runs):
    #         # time_series.append(masker.fit_transform(run))
    #         runts = nib.load(run).get_data()[xc][yc][zc]
    #         time_series.append(runts)
    #         events = np.zeros(len(runts))
    #         for oidx, onsets in enumerate(subject_info[idx].onsets):
    #             for ons in onsets:
    #                 events[int(round(ons/TR))] = oidx + 1
    #         time_series_events.append(events)

    #     session_points = [len(runts)*i for i in range(len(runs))]
    #     processed_signal = signal.detrend(np.concatenate(time_series),
    #                                       bp=session_points)
    #     processed_signal = processed_signal/float(np.mean(processed_signal))*100

    #     raw_bold = np.concatenate(time_series)
    #     onset_ts = np.concatenate(time_series_events).tolist()
    #     assert(len(processed_signal) == len(onset_ts))

    #     coord_event_related_analysis(allcoords[3], raw_bold, processed_signal,
    #                                  onset_ts, TR, 14,
    #                                  os.path.join(results_path, *roi[2]),
    #                                  roi[1])

    # # Number localizer
    # model = 'model001'
    # task = 'task003'
    # results_path = os.path.join('pipelines', 'roi_analysis', subject, model,
    #                             task, 'basic')
    # if not os.path.exists(results_path):
    #     os.makedirs(results_path)

    # labels = ['Numbers', 'Pseudo', 'Consonants']
    # glm_path = os.path.join('pipelines', 'glm_sub_analysis', subject,
    #                         model, task, 'basic',
    #                         'allruns_concat_noderivatives')
    # contrast_path = os.path.join(glm_path, 'con_{0:04d}.nii')
    # contrasts = [1, 2, 3]

    # assert(len(contrasts) == len(labels))
    # for roi in rm.get_rois_mask_in_native_space(sub_id):
    #     if not os.path.exists(os.path.join(results_path, *roi[2])):
    #         os.makedirs(os.path.join(results_path, *roi[2]))
    #     # roi_barplot(roi[0], roi[1], contrasts, labels,
    #     #             os.path.join(results_path, *roi[2]),
    #     #             str([math.ceil(x) for x in rm.get_roi_center(roi[3], roi[4])[3]]))
    #     model_id = '001'
    #     task_id = '003'
    #     TR = 1.5
    #     runs = glob.glob(run_file_template.format(sub_id, task_id, '*'))
    #     runs.sort()
    #     subject_info, _ = subject_model(sub_id, task_id, model_id, runs, True)
    #     time_series = []
    #     time_series_events = []

    #     allcoords = rm.get_roi_center(roi[3], roi[4])
    #     xc, yc, zc = allcoords[0]
    #     #coords = allcoords[2]
    #     #print allcoords
    #     #masker = input_data.NiftiSpheresMasker([tuple(coords)], radius=1.)
    #     for idx, run in enumerate(runs):
    #         # time_series.append(masker.fit_transform(run))
    #         runts = nib.load(run).get_data()[xc][yc][zc]
    #         time_series.append(runts)
    #         events = np.zeros(len(runts))
    #         for oidx, onsets in enumerate(subject_info[idx].onsets):
    #             for ons in onsets:
    #                 events[int(round(ons/TR))] = oidx + 1
    #         time_series_events.append(events)

    #     session_points = [len(runts)*i for i in range(len(runs))]

    #     # processed_signal = signal.detrend(np.concatenate(time_series),
    #     #                                   bp=session_points)
    #     processed_signal = np.concatenate(time_series)
    #     processed_signal = (processed_signal/float(np.mean(processed_signal))*100)-100
    #     processed_signal = signal.detrend(processed_signal)#, bp=session_points)
    #     raw_bold = np.concatenate(time_series)
    #     onset_ts = np.concatenate(time_series_events).tolist()
    #     assert(len(processed_signal) == len(onset_ts))

    #     coord_event_related_analysis(allcoords[3], raw_bold, processed_signal,
    #                                  onset_ts, TR, 14,
    #                                  os.path.join(results_path, *roi[2]),
    #                                  roi[1])


    # # Number query
    # model = 'model002'
    # task = 'task001'
    # results_path = os.path.join('pipelines', 'roi_analysis', subject, model,
    #                             task, 'basic', 'concat')
    # if not os.path.exists(results_path):
    #     os.makedirs(results_path)

    # labels = ['22', '23', '25', '28', '32', '33', '35', '38', '52', '53', '55',
    #           '58', '82', '83', '85', '88']
    # glm_path = os.path.join('pipelines', 'glm_sub_analysis', subject,
    #                         model, task, 'basic',
    #                         'allruns_concat_noderivatives')
    # contrast_path = os.path.join(glm_path, 'con_{0:04d}.nii')
    # contrasts = range(1, 17)

    # assert(len(contrasts) == len(labels))
    # factors = pd.DataFrame(columns=['decade', 'unit'])
    # factors['decade'] = [x[0] for x in labels]
    # factors['unit'] = [x[1] for x in labels]

    # for roi in rm.get_rois_mask_in_native_space(sub_id):
    #     if not os.path.exists(os.path.join(results_path, *roi[2])):
    #         os.makedirs(os.path.join(results_path, *roi[2]))
    #     roi_barplot(roi[0], roi[1], contrasts, labels,
    #                 os.path.join(results_path, *roi[2]),
    #                 str([math.ceil(x) for x in rm.get_roi_center(roi[3], roi[4])[3]]))
        # roi_statistical_model(mask, mask_name, contrasts, labels, factors,
        #                       results_path, 'beta ~ C(decade)*C(unit)')
        # model_id = '002'
        # task_id = '001'
        # TR = 1.5
        # runs = glob.glob(run_file_template.format(sub_id, task_id, '*'))
        # runs.sort()
        # subject_info, _ = subject_model(sub_id, task_id, model_id, runs, True)
        # time_series = []
        # time_series_events = []

        # allcoords = rm.get_roi_center(roi[3], roi[4])
        # xc, yc, zc = allcoords[0]
        # #coords = allcoords[2]
        # #print allcoords
        # #masker = input_data.NiftiSpheresMasker([tuple(coords)], radius=1.)
        # for idx, run in enumerate(runs):
        #     # time_series.append(masker.fit_transform(run))
        #     runts = nib.load(run).get_data()[xc][yc][zc]
        #     time_series.append(runts)
        #     events = np.zeros(len(runts))
        #     for oidx, onsets in enumerate(subject_info[idx].onsets):
        #         for ons in onsets:
        #             events[int(round(ons/TR))] = oidx + 1
        #     time_series_events.append(events)

        # session_points = [len(runts)*i for i in range(len(runs))]
        # processed_signal = signal.detrend(np.concatenate(time_series),
        #                                   bp=session_points)
        # processed_signal = processed_signal/float(np.mean(processed_signal))*100

        # raw_bold = np.concatenate(time_series)
        # onset_ts = np.concatenate(time_series_events).tolist()
        # assert(len(processed_signal) == len(onset_ts))

        # coord_event_related_analysis(allcoords[3], raw_bold, processed_signal,
        #                              onset_ts, TR, 14,
        #                              os.path.join(results_path, *roi[2]),
        #                              roi[1])

    # Parallel(n_jobs=n_proc, verbose=5)(delayed(roi_statistical_model)(mask, mask_name, contrasts, labels, factors,
    #          results_path, 'beta ~ C(decade)*C(unit)') for mask, mask_name, _ in rm.get_rois_mask_in_native_space(sub_id))

    # # Number query F tests
    # model = 'model002'
    # task = 'task001'
    # results_path = os.path.join('pipelines', 'roi_analysis', subject, model,
    #                             task, 'f_tests', 'concat')
    # if not os.path.exists(results_path):
    #     os.makedirs(results_path)

    # labels = ['decade_effect', 'unit_effect', 'interaction',
    #           'f_allnumbereffects']
    # glm_path = os.path.join('pipelines', 'glm_sub_analysis', subject,
    #                         model, task, 'f_tests',
    #                         'allruns_concat_noderivatives')
    # contrast_path = os.path.join(glm_path, 'spmF_{0:04d}.nii')
    # contrasts = [52, 24, 20, 40]

    # assert(len(contrasts) == len(labels))
    # for roi in rm.get_rois_mask_in_native_space(sub_id):
    #     if not os.path.exists(os.path.join(results_path, *roi[2])):
    #         os.makedirs(os.path.join(results_path, *roi[2]))
    #     roi_barplot(roi[0], roi[1], contrasts, labels,
    #                 os.path.join(results_path, *roi[2]),
    #                 str([math.ceil(x) for x in rm.get_roi_center(roi[3], roi[4])[3]]))


    # # Priming
    # model = 'model005'
    # task = 'task001'
    # results_path = os.path.join('pipelines', 'roi_analysis', subject, model,
    #                             task, 'basic')
    # if not os.path.exists(results_path):
    #     os.makedirs(results_path)

    # labels = ['allP', '10sP', '1sP', 'nP', 'fS']
    # glm_path = os.path.join('pipelines', 'glm_sub_analysis', subject,
    #                         model, task, 'basic',
    #                         'allruns_noconcat_noderivatives')
    # contrast_path = os.path.join(glm_path, 'con_{0:04d}.nii')
    # contrasts = [1, 2, 3, 4, 5]

    # assert(len(contrasts) == len(labels))
    # for roi in rm.get_rois_mask_in_native_space(sub_id):
    #     if not os.path.exists(os.path.join(results_path, *roi[2])):
    #         os.makedirs(os.path.join(results_path, *roi[2]))
    #     # roi_barplot(roi[0], roi[1], contrasts, labels,
    #     #             os.path.join(results_path, *roi[2]),
    #     #             str([math.ceil(x) for x in rm.get_roi_center(roi[3], roi[4])[3]]))
    #     model_id = '005'
    #     task_id = '001'
    #     TR = 1.5
    #     runs = glob.glob(run_file_template.format(sub_id, task_id, '*'))
    #     runs.sort()
    #     subject_info, _ = subject_model(sub_id, task_id, model_id, runs, True)
    #     time_series = []
    #     time_series_events = []

    #     allcoords = rm.get_roi_center(roi[3], roi[4])
    #     xc, yc, zc = allcoords[0]
    #     #coords = allcoords[2]
    #     #print allcoords
    #     #masker = input_data.NiftiSpheresMasker([tuple(coords)], radius=1.)
    #     for idx, run in enumerate(runs):
    #         # time_series.append(masker.fit_transform(run))
    #         runts = nib.load(run).get_data()[xc][yc][zc]
    #         time_series.append(runts)
    #         events = np.zeros(len(runts))
    #         for oidx, onsets in enumerate(subject_info[idx].onsets):
    #             for ons in onsets:
    #                 events[int(round(ons/TR))] = oidx + 1
    #         time_series_events.append(events)

    #     session_points = [len(runts)*i for i in range(len(runs))]
    #     processed_signal = signal.detrend(np.concatenate(time_series),
    #                                       bp=session_points)
    #     processed_signal = processed_signal/float(np.mean(processed_signal))*100

    #     raw_bold = np.concatenate(time_series)
    #     onset_ts = np.concatenate(time_series_events).tolist()
    #     assert(len(processed_signal) == len(onset_ts))

    #     coord_event_related_analysis(allcoords[3], raw_bold, processed_signal,
    #                                  onset_ts, TR, 14,
    #                                  os.path.join(results_path, *roi[2]),
    #                                  roi[1])

    # # Priming parametric
    # model = 'model009'
    # task = 'task001'
    # results_path = os.path.join('pipelines', 'roi_analysis', subject, model,
    #                             task, 'basic')
    # if not os.path.exists(results_path):
    #     os.makedirs(results_path)

    # labels = ['number', '10s', '1s']
    # glm_path = os.path.join('pipelines', 'glm_sub_analysis', subject,
    #                         model, task, 'basic',
    #                         'allruns_noconcat_noderivatives')
    # contrast_path = os.path.join(glm_path, 'con_{0:04d}.nii')
    # contrasts = [1, 10, 11]

    # assert(len(contrasts) == len(labels))
    # for roi in rm.get_rois_mask_in_native_space(sub_id):
    #     if not os.path.exists(os.path.join(results_path, *roi[2])):
    #         os.makedirs(os.path.join(results_path, *roi[2]))
    #     # roi_barplot(roi[0], roi[1], contrasts, labels,
    #     #             os.path.join(results_path, *roi[2]),
    #     #             str([math.ceil(x) for x in rm.get_roi_center(roi[3], roi[4])[3]]))
    #     model_id = '009'
    #     task_id = '001'
    #     TR = 1.5
    #     runs = glob.glob(run_file_template.format(sub_id, task_id, '*'))
    #     runs.sort()
    #     subject_info, _ = subject_model(sub_id, task_id, model_id, runs, True)
    #     time_series = []
    #     time_series_events = []

    #     allcoords = rm.get_roi_center(roi[3], roi[4])
    #     xc, yc, zc = allcoords[0]
    #     #coords = allcoords[2]
    #     #print allcoords
    #     #masker = input_data.NiftiSpheresMasker([tuple(coords)], radius=1.)
    #     for idx, run in enumerate(runs):
    #         # time_series.append(masker.fit_transform(run))
    #         runts = nib.load(run).get_data()[xc][yc][zc]
    #         time_series.append(runts)
    #         events = np.zeros(len(runts))
    #         for oidx, onsets in enumerate(subject_info[idx].onsets):
    #             for ons in onsets:
    #                 events[int(round(ons/TR))] = oidx + 1
    #         time_series_events.append(events)

    #     session_points = [len(runts)*i for i in range(len(runs))]
    #     processed_signal = signal.detrend(np.concatenate(time_series),
    #                                       bp=session_points)
    #     processed_signal = processed_signal/float(np.mean(processed_signal))*100

    #     raw_bold = np.concatenate(time_series)
    #     onset_ts = np.concatenate(time_series_events).tolist()
    #     assert(len(processed_signal) == len(onset_ts))

    #     coord_event_related_analysis(allcoords[3], raw_bold, processed_signal,
    #                                  onset_ts, TR, 14,
    #                                  os.path.join(results_path, *roi[2]),
    #                                  roi[1])




    # mask = os.path.join('/volatile/martin_local/Share_Martin_Christophe/Martin_PHD_files/QA_checks/reading_lefthemisphere.nii')
    # mask = rm.resample_roi_to_mni_template(mask)
    # mask = rm.from_mni_to_native(mask, '001', '002')
    # mask = nib.load(mask).get_data()
    # mask = mask.reshape(mask.shape[0:3])
    # mask = mask == 1
    # mask_name = 'aha'


    # # Number query F tests individual effect and interactions
    # model = 'model002'
    # task = 'task001'
    # results_path = os.path.join('pipelines', 'roi_analysis', subject, model,
    #                             task, 'anova_effects')
    # if not os.path.exists(results_path):
    #     os.makedirs(results_path)

    # labels = ['d1', 'd2', 'd3', 'u1', 'u2', 'u3',
    #           'i1', 'i2', 'i3', 'i4', 'i5', 'i6', 'i7', 'i8', 'i9']
    # glm_path = os.path.join('pipelines', 'glm_sub_analysis', subject,
    #                         model, task, 'f_tests',
    #                         'allruns_noconcat_noderivatives')
    # contrast_path = os.path.join(glm_path, 'con_{0:04d}.nii')
    # contrasts = [17, 18, 19, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31, 32, 33]

    # assert(len(contrasts) == len(labels))
    # for mask, mask_name in rm.get_rois_mask_in_native_space(sub_id):
    #     roi_barplot(mask, mask_name, contrasts, labels, results_path)

   # Density plots p values
    model = 'model002'
    task = 'task001'
    results_path = os.path.join('pipelines', 'roi_analysis', subject, model,
                                task, 'f_tests', 'density')
    if not os.path.exists(results_path):
        os.makedirs(results_path)

    labels = ['interaction','f_allnumbereffects']
    glm_path = os.path.join('pipelines', 'glm_sub_analysis', subject,
                            model, task, 'f_tests',
                            'allruns_concat_noderivatives')
    contrast_path = os.path.join(glm_path, 'spmF_{0:04d}.nii')
    contrasts = [20, 40]

    assert(len(contrasts) == len(labels))
    gm_img = os.path.join('pipelines', 'preprocessing', subject,
                          'highres001_graymatter_mask_bold_native.nii')
    gm_mask = rm.get_roi_mask(gm_img)
    diff_map = contrast_path.format(40)
    inter_map = contrast_path.format(20)
    df2 = 219*8 - 32
    pedges = [1, 0.01, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001, 0]
    pedges.reverse()
    for roi in rm.get_rois_mask_in_native_space(sub_id):
        if not os.path.exists(os.path.join(results_path, *roi[2])):
            os.makedirs(os.path.join(results_path, *roi[2]))
        gm_adj_roi = gm_mask * roi[0]
        diff = np.ravel(rm.extract_data_in_mask(diff_map, gm_adj_roi))
        inter = np.ravel(rm.extract_data_in_mask(inter_map, gm_adj_roi))
        diffp = [-np.log10(sst.f.sf(x, 15, df2)) for x in diff]
        interp = [-np.log10(sst.f.sf(x, 9, df2)) for x in inter]
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        extra_text = str([math.ceil(x) for x in rm.get_roi_center(roi[3], roi[4])[3]])
        g = sns.jointplot(x=np.array(diffp), y=np.array(interp), kind='scatter', stat_func=None, xlim=(0,np.max(diffp)+1), ylim=(0,np.max(interp)+1)).set_axis_labels("-log10(p_value) number difference", "-log10(p_value) interaction")
        # , xlim=(0,1), ylim=(0,1)
        # plt.show()

        # H, xedges, yedges = np.histogram2d(interp, diffp, bins=(pedges, pedges),
        #                                    normed=False)
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        extra_text = str([math.ceil(x) for x in rm.get_roi_center(roi[3], roi[4])[3]])
        # ax.set_title(extra_text + ' '+ str(np.sum(gm_adj_roi)))
        # im = plt.imshow(H, interpolation='none', origin='lower',
        #                 extent=None)


        g.savefig(os.path.join(os.path.join(results_path, *roi[2]),
                    roi[1]+'.png'))

        plt.close()

        # roi_barplot(roi[0], roi[1], contrasts, labels,
        #             os.path.join(results_path, *roi[2]),
        #             str([math.ceil(x) for x in rm.get_roi_center(roi[3], roi[4])[3]]))
