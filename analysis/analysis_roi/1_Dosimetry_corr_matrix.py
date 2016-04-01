
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
"""    
Plotting a diagonal correlation matrix
======================================

"""

# import the dataframe

df = pd.read_csv("/home/edogerde/Bureau/ROI.csv")

# seaborn
sns.set(style="white")

# Generate a large random dataset
d = df[["label", "Patients","mean" ]]
dt= d.pivot(columns="label", index = "Patients" , values= "mean")

# Compute the correlation matrix
corr = dt.corr()

### Create the matrix figure with seaborn
# Generate a mask for the upper triangle
mask = np.zeros_like(corr, dtype=np.bool)
mask[np.triu_indices_from(mask)] = True

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(20, 30))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(220, 10, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.7,
            square=True, xticklabels=5, yticklabels=5,
            linewidths=.8, cbar_kws={"shrink": .5}, ax=ax)
plt.show()

### Create the matrix figure with mathplot
plt.imshow(corr, cmap=plt.cm.hot)
plt.colorbar()
plt.show()

### Pattern amount subject
sns.set(style="whitegrid")
sns.barplot(y="mean", x="Patients", hue = "label" , data= d,)   
plt.show()