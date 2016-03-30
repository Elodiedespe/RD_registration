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



if __name__ == '__main__':

    df_data = pd.read_csv("/home/edogerde/Bureau/data_analysis/df_F_Bundles.csv")
	
	# Add a columns RT to dataframes 
    RT=[]
    for row in df_data["Traitement"]:
    	if row == "C":
    		RT.append('NO')
        else:
            RT.append("YES")

    df_data["RT"]= RT

    Lat_hemisphere = ["right", "left", "interhemispheric"]

	# Create New dataframe
    df_final_bundle = pd.DataFrame(columns=['Score','bundle', 'lateralization','coeff', 'pvalue' , 'pval_bonferroni', 'signi_bonferonni'])
    

# List of VD
    VDs = ["dQIT", "dQIP", "dQIV", "dIMT", "dIVT", "dCubes", "dVocabulaire"]

    for VD in VDs:

        for lat in Lat_hemisphere:


        # Slect the dataframe
            df = df_data[(df_data.last_delta==1) & (df_data.RT=="YES") & (df_data.Hemisphere==lat)][["Patients", VD, "mean", "Bundle", "AgeAtDiagnosis", "CSP"]]
            df_drop = df.dropna(subset=[VD])


            # Run the multiple regression models for each Roi
        
            for cnt, r in enumerate(df["Bundle"]):
            # print "print handling roi %i" % r
                X = np.array(df_drop[df_drop.Bundle==r][['AgeAtDiagnosis', 'mean', "CSP"]])
                Y = np.array(df_drop[df_drop.Bundle==r][VD])
                
                # Fit and summary:
                model_Bundles = sm.OLS(Y, X).fit()
                print(model_Bundles.summary())

                # Correction of multiple comparaison with Bonferronin
                pvals = model_Bundles.pvalues
                pvals_fwer = multicomp.multipletests(pvals, alpha = 0.05, method = 'bonferroni')

                # concatenate all the results into a pd dataframe
                df_final_bundle.loc[len(df_final_bundle)] = [VD, r, lat, model_Bundles.params, model_Bundles.pvalues, pvals_fwer[1], pvals_fwer[0]]

 