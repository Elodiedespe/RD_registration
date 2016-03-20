
# Modules
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
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


if __name__ == '__main__':

    df_data = pd.read_csv("/home/edogerde/Bureau/data_analysis/df_F.csv")

    xml_path = "/home/edogerde/Bureau/delineation_space/lpba40.label.xml"

    resSGRT = []
    
    df_final= pd.DataFrame(columns=['Score','roi', 'coeff', 'pvalue' , 'pval_bonferroni', 'signi_bonferonni'])
    
# List of VD
    VDs = ["dQIT", "dQIP", "dQIV", "dIMT", "dIVT", "dCubes", "dVocabulaire"]

    for VD in VDs:
        # Slect the dataframe
        df = df_data[(df_data.last_delta==1) & (df_data.RT=="YES")][["Patients", VD, "mean", "Global_mean", "label", "AgeAtDiagnosis", "CSP"]]
        df_drop = df.dropna(subset=[VD])

        # Select x and y
        y, x = df_drop[VD], df_drop["Global_mean"]

        # Run the regression
        beta, beta0, r_value, p_value, std_err = stats.linregress(x, y)
        print("y=%f x + %f, r:%f, r-squared:%f, \np-value:%f, std_err:%f"
        % (beta, beta0, r_value, r_value**2, p_value, std_err))

        # Plot the line
        yhat = beta * x + beta0  # regression line
        plt.plot(x, yhat, 'r-', x, y, 'o')
        plt.xlabel('Global mean')
        plt.ylabel(VD)
        # plt.show()

        # New dataframe for VD without global effect
        df_roi = pd.DataFrame(columns = ['Roi', VD, 'AgeAuDiagnostic', "meanRoi", "CSP"])
        df_roi['Roi'] = df_drop["label"]
        df_roi['AgeAuDiagnostic'] = df_drop["AgeAtDiagnosis"]
        df_roi['meanRoi'] = df_drop["mean"]
        df_roi['CSP'] = df_drop["CSP"]

        # Detrend global mean effect
        Roi_effect = y - beta*x
        df_roi[VD] = Roi_effect

        # Create list for roi
        name, rank, roi_label = read_roiLabel(xml_path)

        # Run the multiple regression models for each Roi
        f_wthgl = open(os.path.join('elusion.txt'), 'w')
        f_wgl = open(os.path.join('aloa.txt'), 'w')

        for cnt, r in enumerate(roi_label):
            # print "print handling roi %i" % r
            X = np.array(df_roi[df_roi.Roi==r][['AgeAuDiagnostic', 'meanRoi', "CSP"]])
            Y = np.array(df_roi[df_roi.
                Roi==r][VD])
            # Fit and summary:
            model_GRT = sm.OLS(Y, X).fit()
            print (" %s Regression sans effet global RT" % (r))
            print (" %s Regression sans effet global RT" % (VD))
            print (" %s Regression sans effet global RT" % (r))
            print(model_GRT.summary())
            f_wthgl.write(VD)
            f_wthgl.write(str(model_GRT.summary()))
            # concatenate all the results into a list
            """for resGRT in model_GRT:
                low, upp = res.confint().T   # unpack columns 
                res_all.append(numpy.concatenate(([resGRT.llf], resGRT.params, resGRT.tvalues, resGRT.pvalues, 
                                   low, upp)))"""

            X1 = np.array(df_drop[df_drop.label==r][['AgeAtDiagnosis', 'mean',"CSP"]])
            Y1 = np.array(df_drop[df_drop.label==r][VD])
            ## Fit and summary:
            model_SGRT = sm.OLS(Y1, X1).fit()
            print (" %s Regression avec effet global RT" % (VD))  
            print (" %s Regression avec effet global RT" % (r))
            print(model_SGRT.summary())
            pvals = model_SGRT.pvalues
            #f_wgl.write(VD)
            #f_wgl.write(str(model_SGRT.summary()))

            # concatenate all the results into a list
            """for resSGRT in model_SGRT:
                low, upp = resSGRT.confint().T   # unpack columns 
                res_all.append(numpy.concatenate(([resSGRT.llf], resSGRT.params, resSGRT.tvalues, resSGRT.pvalues, 
                                   low, upp)))"""

            pvals_fwer = multicomp.multipletests(pvals, alpha = 0.05, method = 'bonferroni')
            
            df_final.loc[len(df_final)] = [VD, r, model_SGRT.params, model_SGRT.pvalues, pvals_fwer[1], pvals_fwer[0]]

    f_wgl.write(str(resSGRT))
    f_wthgl.close()
    f_wgl.close()