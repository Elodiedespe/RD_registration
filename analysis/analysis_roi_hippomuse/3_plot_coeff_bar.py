import pandas as pd
import seaborn as sns
#import matplotlib libary
import matplotlib.pyplot as plt
import numpy as np

#df = pd.read_csv('/home/edogerde/Bureau/betaGlobalmean.csv' )
#df = pd.read_csv('/home/edogerde/Bureau/betaQIT.csv' )
df = pd.read_csv('/home/edogerde/Bureau/for_git/RD_registration/analysis/analysis_roi/lalala.csv' )
#labels = ["Item","Repere journalier","heure","Quand"]
labels=["dQIT","dQIP","dQIV","dIMT","dIVT","dCubes","dVocabulaire"]
Y = np.array(df[df.variable == "dQIP"]['coeff'])

#Y = np.array(df['coeff'].astype(float).values)
X = np.array(range(56))
#Y_error = np.array(df["std"].astype(float).values)
#Y = np.array(df['coeff'])
#X = np.array([1,2,3,4,5,6,7])
#Y_error = np.array(df["std_err"])


#plot data
plt.plot(X,Y,linestyle="dashed", marker="o", color="green")
#plt.errorbar(X, Y, yerr=Y_error, linestyle="None", marker="None", color="green")

#labels
plt.xticks(X)
plt.xlabel("Roi")
plt.ylabel("coeff")


#show plot
plt.show()
