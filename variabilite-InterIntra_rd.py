import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

df = pd.read_csv('/home/edogerde/Desktop/variation_interIntra_RD.csv') 
plt.plot(df.Results[df.Dose=='Dose 1'])
plt.plot(df.Results[df.Dose=='Dose 2'])
plt.plot(df.Results[df.Dose=='Dose 3'])
plt.plot(df.Results[df.Dose=='Dose 4'])
plt.plot(df.Results[df.Dose=='Dose 5'])
plt.xlim(0, 3)
plt.ylim(20, 750)
plt.xlabel("Dose en fonction de plusieurs points dans le cerveau")
plt.ylabel("Dose en u.a") 

