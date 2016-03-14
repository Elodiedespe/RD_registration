
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

	
	
	# Path
	Path = "/home/edogerde/Bureau/extraction_roi"
	df_data = pd.read_csv("/home/edogerde/Bureau/data_analysis/df_F.csv")
	df = df_data[df_data.last_delta==1][["Patients", "label", "mean"]]
	df_2= df_data[["Patients","AgeAtRT", "Traitement"]]
	
	xml_path = "/home/edogerde/Bureau/delineation_space/lpba40.label.xml"


 	# Get subject id
 	
	valid_subject_dirs = [os.path.join(Path, dir_name)
                        for dir_name in os.listdir(Path)
                        if os.path.isdir(os.path.join(Path, dir_name))]

	############# Global mean #####################################################

	#Create new dataframe and merge it to the main one 
	new_Gm_df= pd.DataFrame(columns=['Patients','Global_mean'])

	Patients=[]
	Global_means=[]
	import pdb; pdb.set_trace()  # breakpoint 247c0009 //

	# Get the global mean for each subject
	for sujb_path in valid_subject_dirs:
		suj_id = os.path.basename(sujb_path)   
   	 	global_mean = df[df.Patients == suj_id].mean(axis=0).iloc[1]
    	print global_mean    
    	Patients.append(suj_id)
    	Global_means.append(global_mean)

	new_Gm_df["Patients"]=Patients
	new_Gm_df["Global_mean"]=Global_means

	# Add the results to the dataframe
	df_Gm_data = pd.merge(df_data, new_Gm_df, on='Patients')
	import pdb; pdb.set_trace()  # breakpoint 98000e93 //


	# Add a columns RT to dataframes 
	"""RT=[]
	for row in df_data["Traitement"]:
		if row == "C":
			RT.append('NO')
    	else:
        	RT.append("YES")

	df_data["RT"]= RT
	df_Gm_data["RT"]= RT"""
	import pdb; pdb.set_trace()  # breakpoint 59b7055b //

	# List of VD
	VDs= ["dQIP", "dQIV", "dIMT", "dIVT","dCubes", "dVocabulaire"]

	for VD in VDs:
		
		# Slect the dataframe
		df = df_Gm_data[(df_Gm_data.last_delta==1) & (df_Gm_data.RT=='YES')][["Patients", VD, "mean", "label","AgeAtDiagnosis","CSP"]]
		df_drop=df.dropna(subset=[VD])
		
		import pdb; pdb.set_trace()  # breakpoint 0acdc1f0 //
		
		# Select x and y
		y, x = df_drop[VD], df_drop["mean"]
		import pdb; pdb.set_trace()  # breakpoint 9de9f010 //
		
		df_roi= pd.DataFrame(columns=['Roi', VD, 'AgeAuDiagnostic', "meanRoi", "CSP"])
		df_roi['Roi'] = df_drop["label"]
		df_roi['AgeAuDiagnostic']= df_drop["AgeAtDiagnosis"]
		df_roi['meanRoi']= df_drop["mean"]
		df_roi['CSP']= df_drop["CSP"]

		import pdb; pdb.set_trace()  # breakpoint 075a7bab //

		# Run the regression
		beta, beta0, r_value, p_value, std_err = stats.linregress(x,y)
		print("y=%f x + %f, r:%f, r-squared:%f, \np-value:%f, std_err:%f"
		% (beta, beta0, r_value, r_value**2, p_value, std_err))
		import pdb; pdb.set_trace()  # breakpoint d298a084 //
		
		# Plot the line
		yhat = beta * x + beta0 # regression line
		plt.plot(x, yhat, 'r-', x, y,'o')
		plt.xlabel('Global mean')
		plt.ylabel(VD)
		plt.show()
		import pdb; pdb.set_trace()  # breakpoint 9e52a840 //

		# Detrend global mean effect
		Roi_effect = y - beta*x
		df_roi[VD] = Roi_effect
		
		# Create list for roi
		name, rank, roi_label = read_roiLabel(xml_path)

		# Run the multiple regression models for each Roi
		for r in roi_label:

			X = np.array(df_roi[df_roi.Roi==r][['AgeAuDiagnostic', 'meanRoi',"CSP"]])
			Y = np.array(df_roi[df_roi.Roi==r][VD])
			## Fit and summary:
    		model_GRT = sm.OLS(Y, X).fit()
    		print (VD + "Regression sans effet global RT") 
    		print (r + "Regression sans effet global RT")
    		print(model_GRT.summary())
    		
    		# concatenate all the results into a list
    		resGRT_all = []
    		for resGRT in model_GRT:
    			low, upp = res.confint().T   # unpack columns 
    			res_all.append(numpy.concatenate(([resGRT.llf], resGRT.params, resGRT.tvalues, resGRT.pvalues, 
                   				low, upp)))
			import pdb; pdb.set_trace()  # breakpoint 79fc44c4 //

    		X1 = np.array(df[df.label==r][['AgeAtDiagnostic', 'mean',"CSP"]])
    		Y2 = np.array(df[df.label==r][VD])
			## Fit and summary:
    		model_SGRT = sm.OLS(Y1, X2).fit()
    		print (VD + "Regression avec effet global RT")  
    		print (r + "Regression avec effet global RT") 
    		print(model_SGRT.summary())

    		# concatenate all the results into a list
    		res_all = []
    		for resSGRT in model_SGRT:
    			low, upp = resSGRT.confint().T   # unpack columns 
    			res_all.append(numpy.concatenate(([resSGRT.llf], resSGRT.params, resSGRT.tvalues, resSGRT.pvalues, 
                   				low, upp)))
