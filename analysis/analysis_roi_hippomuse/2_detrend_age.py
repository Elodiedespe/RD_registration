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
        print rank, name
        roi_label.append(rank)

    roi_label=map(int,roi_label)

    return name, rank, roi_label

if __name__ == '__main__':

    df_data = pd.read_csv("/home/edogerde/Bureau/data_analysis/hippomuseDataBaseFinal.csv")

    TEST = ["reconnaissanceImage", "reconnaissanceOdeur","hedonicite","musique","episodicite", "what", "repereJournalier","heure","when"] 

    df_GlobalRt= pd.read_csv ('/home/edogerde/Bureau/data_analysis/Global_mean.csv')
    df_Add_GlobalRT = pd.merge(df_GlobalRt,df_data, on= "Patients") 

    # Write the result into a txt file
    f_regression = open(os.path.join('regressionAgeScores.txt'), 'w')

    for t in TEST:

        if t == 'what' or t == 'repereJournalier' or t == 'heure' or t == 'when':
            df = df_data[(df_data.Test == t) & (df_data.Phase =="P2")& (df_data.label==0)][["Age", "Result"]]
            df_drop = df.dropna(subset=["Result"])
        else:
            # Slect the dataframe
            df = df_data[(df_data.Test == t) & (df_data.Item =="totale") & (df_data.label==0)][["Age","Result"]]
            df_drop = df.dropna(subset=["Result"])

        # Select x and y
        y, x = df_drop["Result"], df_drop["Age"]

        # Run the regression
        beta, beta0, r_value, p_value, std_err = stats.linregress(x.astype(float).values, y.astype(float).values)
        print("y=%f x + %f, r:%f, r-squared:%f, \np-value:%f, std_err:%f"
        % (beta, beta0, r_value, r_value**2, p_value, std_err))

        f_regression.write(str("Age , %s, y=%f x + %f, r:%f, r-squared:%f, \np-value:%f, std_err:%f"
        % (t, beta, beta0, r_value, r_value**2, p_value, std_err)))

        # Plot the line
        yhat = beta * x + beta0  # regression line
        plt.plot(x, yhat, 'r-', x, y, 'o')
        #plt.ylim((-1,4))
        #plt.xlim((5,13))
        plt.xlabel('Age')
        plt.ylabel('score de %s'%(t))
        plt.show()
        
        z = df_Add_GlobalRT["Global_mean"]
        beta, beta0, r_value, p_value, std_err = stats.linregress(x.astype(float).values, z.astype(float).values)
        print("y=%f x + %f, r:%f, r-squared:%f, \np-value:%f, std_err:%f"
        % (beta, beta0, r_value, r_value**2, p_value, std_err))
        
        # Plot the line
        yhat = beta * x + beta0  # regression line
        plt.plot(x, yhat, 'r-', x, z, 'o')
        #plt.ylim((-1,4))
        #plt.xlim((5,13))
        plt.xlabel('Global_mean')
        plt.ylabel('score de %s'%(t))
        plt.show()
        
        f_regression.write(str(" Global mean,  %s, y=%f x + %f, r:%f, r-squared:%f, \np-value:%f, std_err:%f"
        % (t, beta, beta0, r_value, r_value**2, p_value, std_err)))

    # Close the txt file    
    f_regression.close()

