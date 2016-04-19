import numpy as np
import pandas as pd


dfInput = pd.read_csv("/home/edogerde/Bureau/Hippomuse/patients_hippomuse/results_381.csv")

dfOutput = pd.DataFrame(columns=['Sujet', 'Age', 'GroupeAge', 'AgeRT', 'RT', 'TypeRT', 'Phase', 'Item', 'Test', 'Result'])

compt = 0

for idx in dfInput.index:
    cur = dfInput.loc[idx]
    for j in range(1,110):
        colName = dfInput.columns[j]
        itemName = colName.split('_')[0]
        testName = colName.split('_')[1]
        phaseName = colName.split('_')[2]
        result = cur[colName]
        sujet = cur['sujet']
        age = cur['ageToday']
        groupeAge = cur['groupeAge']
        ageRT = cur['ageRT']
        rt = cur['presenceRT']
        typeRT = cur['typeRT']
        dfOutput.loc[compt] = [sujet, age, groupeAge, ageRT, rt, typeRT, phaseName, itemName, testName, result]
        compt += 1

dfOutput.to_csv('/home/edogerde/Bureau/hippomuseDataBasetest.csv', index=False)

""" I merge the Df with Roi.csv and Clinical_data.csv"""
