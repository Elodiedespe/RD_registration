
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

df_data = pd.read_csv("/home/edogerde/Bureau/data_analysis/df_F_Bundles.csv")
df = df_data[(df_data.last_delta==1) & (df_data.RT=="YES") & (df_data.Hemisphere=="left") & (df_data.Bundle=="Uncinate")][["Patients", "dVocabulaire", "mean", "Bundle", "AgeAtDiagnosis", "CSP"]]
df_drop = df.dropna(subset=["dVocabulaire"])

g = sns.jointplot(x="mean", y="dVocabulaire", data=df_drop, kind="reg")
g = sns.jointplot(x="AgeAtDiagnosis", y="dVocabulaire", data=df_drop, kind="reg")
g = sns.jointplot(x="CSP", y="dVocabulaire", data=df_drop, kind="reg")
