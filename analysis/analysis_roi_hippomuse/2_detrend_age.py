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
    roi_label=[]

    for label in root.findall('label'):
        rank = label.get('id')
        name = label.get('fullname')
        #print rank, name
        roi_label.append(rank)

    roi_label=map(int,roi_label)

    return name, rank, roi_label

if __name__ == '__main__':

    df_data = pd.read_csv("/home/edogerde/Bureau/data_analysis/hippomuseDataBaseFinal.csv")

    TEST = ["reconnaissanceImage", "reconnaissanceOdeur","musique", "what", "repereJournalier","heure","when", "episodicite"]
    
    df_final= pd.DataFrame(columns=['Score','roi', 'coeff', 'pvalue' , 'pval_bonferroni', 'signi_bonferonni'])
    
    df_ScoreAge = pd.DataFrame(columns=['variable','test', 'coeff', 'r_squared' , 'p_value', 'std_err'])
    
    df_ScoreGlobalmean = pd.DataFrame(columns=['variable','test', 'coeff', 'r_squared' , 'p_value', 'std_err'])

    #df_detrendAge = pd.DataFrame(columns=[['Patients',"test",'Result']])
    df_detrendAge = dict()
    xml_path = "/home/edogerde/Bureau/delineation_space/lpba40.label.xml"

    # Write the result into a txt file
    f_regression = open(os.path.join('regressionAgeScores.txt'), 'w')

    for t in TEST[3:8]:

        if t == 'what' or t == 'repereJournalier' or t == 'heure' or t == 'when' or t == "episodicite":
            df = df_data[(df_data.Test == t) & (df_data.Phase =="P2") & (df_data.Item =="delta")][["Patients","Age", "Result", "Global_mean", "label",'AgeAtDiagnosis', 'mean']]
            df_drop = df.dropna(subset=["Result"])
        else:
            # Slect the dataframe
            df = df_data[(df_data.Test == t) & (df_data.Item =="totale")][["Patients","Age","Result", "Global_mean", "label", 'AgeAtDiagnosis', 'mean']]
            df_drop = df.dropna(subset=["Result"])

        # Select x and y
        y, x = df_drop["Result"], df_drop["Age"]

        # Run the regression
        beta, beta0, r_value, p_value, std_err = stats.linregress(x.astype(float).values, y.astype(float).values)
        print("y=%f x + %f, r:%f, r-squared:%f, \np-value:%f, std_err:%f"
        % (beta, beta0, r_value, r_value**2, p_value, std_err))

        df_ScoreAge.loc[len(df_ScoreAge)] = ["Age", t, beta, r_value**2, p_value, std_err ]
        df_ScoreAge.to_csv("regressionScoreAge.csv")        

        ScoreWotAge = y.values - beta*x.values
        
        # Detrend Age effect
        df_detrendAge.setdefault('Patients', []).extend(df_drop["Patients"].values.tolist())
        df_detrendAge.setdefault('Test', []).extend([t] * len(ScoreWotAge))
        df_detrendAge.setdefault('Result', []).extend(ScoreWotAge.tolist())
                
        """
        df_detrendAge["Patients"] = df_drop["Patients"]
        df_detrendAge["Result"] = ScoreWotAge 
        #df_detrendAge["test"] =  df_drop[t]
        
        import pdb; pdb.set_trace()  # breakpoint d4b8c8ab //
        
        
        # Plot the line
        yhat = beta * x + beta0  # regression line
        plt.plot(x, yhat, 'r-', x, y, 'o')
        #plt.ylim((-1,4))
        #plt.xlim((5,13))
        plt.xlabel('Age')
        plt.ylabel('score de %s'%(t))
        #plt.show()
        
        x1 = df_drop["Global_mean"]
        beta, beta0, r_value, p_value, std_err = stats.linregress(x1.astype(float).values, y.astype(float).values)
        print("y=%f x + %f, r:%f, r-squared:%f, \np-value:%f, std_err:%f"
        % (beta, beta0, r_value, r_value**2, p_value, std_err))

        df_ScoreGlobalmean.loc[len(df_ScoreGlobalmean)] = ["Global_mean", t, beta, r_value**2, p_value, std_err ]
        df_ScoreGlobalmean.to_csv("regressionScoreGlobalmean.csv")
        
        # Plot the line
        yhat = beta * x1 + beta0  # regression line
        plt.plot(x1, yhat, 'r-', x1, y, 'o')
        #plt.ylim((-1,4))
        #plt.xlim((5,13))
        plt.xlabel('Global_mean')
        plt.ylabel('score de %s'%(t))
        #plt.show()

        f_regression.write(str(" Global mean,  %s, y=%f x + %f, r:%f, r-squared:%f, \np-value:%f, std_err:%f"
        % (t, beta, beta0, r_value, r_value**2, p_value, std_err)))


        # Create list for roi
        name, rank, roi_label = read_roiLabel(xml_path)


        for cnt, r in enumerate(roi_label):

            X1 = np.array(df_drop[df_drop.label==r][['mean', "AgeAtDiagnosis"]])
            Y1 = np.array(df_drop[(df_drop.label==r)]["Result"])
            
            ## Fit and summary:
            model_SGRT = sm.OLS(Y1, X1).fit()

            # Correction of multiple comparaison with Bonferronin
            pvals = model_SGRT.pvalues
            pvals_fwer = multicomp.multipletests(pvals, alpha = 0.05, method = 'bonferroni')

            # concatenate all the results into a pd dataframe
            df_final.loc[len(df_final)] = [t, r, model_SGRT.params, model_SGRT.pvalues, pvals_fwer[1], pvals_fwer[0]]
    
    # Close the txt file    
    f_regression.close()
"""
df_detrendAge = pd.DataFrame(df_detrendAge)