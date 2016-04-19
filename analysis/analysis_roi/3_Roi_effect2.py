
# Modules
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import scipy
from scipy import stats
import statsmodels.formula.api as smfrmla
import statsmodels.api as sm
import xml.etree.ElementTree as ET
import statsmodels.sandbox.stats.multicomp as multicomp

def read_roiLabel(xml_path):

    tree = ET.parse(xml_path)
    root = tree.getroot()
    roi_label = []

    for label in root.findall('label'):
        rank = label.get('id')
        name = label.get('fullname')
        # print rank, name
        roi_label.append(rank)

    roi_label = [int(r) for r in roi_label]

    return name, rank, roi_label

def get_roi_mask(roi_filename, label_number):
    mask1 = nib.load(roi_filename)
    roi_affine = mask1.get_affine()
    mask= mask1.get_data()  
    mask = mask.reshape(mask.shape[0:3])

    mask = mask == label_number
    roi_nii = nib.Nifti1Image(mask.astype(np.int8), roi_affine)

    nib.save(roi_nii,os.path.join(maskfile,"%s.nii.gz"%(label_number)))
    return mask , roi_nii



if __name__ == '__main__':

    df_data = pd.read_csv("/home/edogerde/Bureau/data_analysis/df_FINAL.csv")

    xml_path = "/home/edogerde/Bureau/delineation_space/lpba40.label.xml"

    resSGRT = []
    
    # Path for plotting.plot_roi
    atlas_nii = "/home/edogerde/Bureau/Atlas_plt_roi/atlas.nii.gz"
    anat_nii = "/home/edogerde/Bureau/Atlas_plt_roi/anat.nii"
    maskfile = "/home/edogerde/Bureau/Atlas_plt_roi"
    
    df_final= pd.DataFrame(columns=['Score','roi', 'coeff', 't value', 'pvalue' , 'pval_bonferroni', 'signi_bonferonni', 'Rsquare', 'std'])
    df_final_GR= pd.DataFrame(columns=['Score','roi', 'coeff', 'pvalue' , 'pval_bonferroni', 'signi_bonferonni', 'Rsquare'])
    f_regression = open(os.path.join('regressionGlobalmeanScoresQi.txt'), 'w')
    df_ScoreGlobalmean = pd.DataFrame(columns=['variable','test', 'coeff', 'r_squared' , 'p_value', 'std_err'])
    df_corr = pd.DataFrame(columns=['variable','test', 'coeff', 'p_value'])
    df_corrAgDvsRt = pd.DataFrame(columns=['variable','test', 'coeff', 'p_value'])

# List of VD
    #VDs = ["QIT", "QIP", "QIV", "IMT", "IVT", "Cubes", "Vocabulaire"]
    VDs = ["dQIT", "dQIP", "dQIV", "dIMT", "dIVT", "dCubes", "dVocabulaire"]
    
    #for VD in VDs:
        #load the dataframe and set the seaborn
        #df_behavioural= df_data[["Patients", VD, "Traitement", "AgeAtDiagnosis"]]
        #sns.set(style="whitegrid", color_codes=True)
        
        #sns.barplot(x="Traitement", y=VD, data=df_behavioural)       
        #plt.show()

    
    for VD in VDs[0:7]:
        # Slect the dataframe
        #df = df_data[(df_data.last_delta==1) & (df_data.RT=="YES")][["Patients", "mean", "Global_mean", "label", "AgeAtDiagnosis", "CSP"]]
        #df = df_data[(df_data.Rang==1) & (df_data.RT=="YES")][["Patients", VD, "mean", "Global_mean", "label", "AgeAtDiagnosis", "CSP"]]
        df= df_data[["Patients", "mean", "Global_mean", VD, "label", "AgeAtDiagnosis"]]
        df_drop = df.dropna(subset=[VD])

        # Select x and y 
        y, x = df_drop[VD], df_drop["Global_mean"]
        
        # Run the regression
        beta, beta0, r_value, p_value, std_err = stats.linregress(x, y)
        print("y=%f x + %f, r:%f, r-squared:%f, \np-value:%f, std_err:%f"
        % (beta, beta0, r_value, r_value**2, p_value, std_err))
        
        f_regression.write(str("Global_mean %s, y=%f x + %f, r:%f, r-squared:%f, \np-value:%f, std_err:%f"
        % (VD, beta, beta0, r_value, r_value**2, p_value, std_err)))
        df_ScoreGlobalmean.loc[len(df_ScoreGlobalmean)] = ["Global_mean", VD, beta, r_value**2, p_value, std_err ]


     
        # Plot the line
        yhat = beta * x + beta0  # regression line
        plt.plot(x, yhat, 'r-', x, y, 'o')
        plt.xlabel('Global mean')
        plt.ylabel(VD)
        #plt.show()


        # Create list for roi
        name, rank, roi_label = read_roiLabel(xml_path)

        # Run the multiple regression models for each Roi

        for cnt, r in enumerate(roi_label):

            X1 = np.array(df_drop[df_drop.label==r][['mean','Global_mean']])

            #X1 = np.array(df_drop[df_drop.label==r]["mean"])
            Y1 = np.array(df_drop[df_drop.label==r][VD])
            ## Fit and summary:
            rho, pval = scipy.stats.spearmanr(X1,Y1)
            model_SGRT = sm.OLS(Y1, X1).fit()
            #print(model_SGRT.summary())
            yAgeDia = np.array(df_drop[df_drop.label==r][['AgeAtDiagnosis']])
            rho2, pval2 = scipy.stats.spearmanr(X1,yAgeDia)


            # Correction of multiple comparaison with Bonferroninipyth
            pvals = model_SGRT.pvalues
            pvals_fwer = multicomp.multipletests(pvals, alpha = 0.05, method = 'fdr_bh')

            # concatenate all the results into a pd dataframe
            df_corrAgDvsRt.loc[len(df_corrAgDvsRt)] = [VD, r, rho2, pval2]
            df_corr.loc[len(df_corr)] = [VD, r, rho, pval]
            df_final.loc[len(df_final)] = [VD, r, model_SGRT.params, model_SGRT.tvalues, model_SGRT.pvalues, pvals_fwer[1], pvals_fwer[0], model_SGRT.rsquared, model_SGRT.bse]
            #df_final_GR.loc[len(df_final)] = [VD, r, model_GRT.params, model_GRT.pvalues, pvals_GR_fwer[1], pvals_GR_fwer[0], model_SGRT.rsquared]




    f_regression.close()
