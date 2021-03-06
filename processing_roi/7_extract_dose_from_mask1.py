# import module
import os
import numpy as np
import nibabel as nib
import xml.etree.ElementTree as ET 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import shutil 



###################################################################

"""
Path= "/home/edogerde/Bureau/extraction_roi"
Path_label= "/media/edogerde/MY PASSPORT/data_dosimetry/results_from_label_to_rd_after_ants"
Path_rd= "/media/edogerde/MY PASSPORT/these/sujets/subject_24_rd_resample"
# Create new empty directories with the name of the patients
Suj = pd.read_csv("/media/edogerde/MY PASSPORT/clinical_dataFinal.csv")
sub_name= Suj["anonym_nom"]
for name in sub_name:
    print(name)
    os.mkdir(os.path.join(Path, name))

#Create a list of patients directories
valid_subject_dirs = [os.path.join(Path, dir_name)
                      for dir_name in os.listdir(Path)
                      if os.path.isdir(os.path.join(Path, dir_name))]


#copy all the file in the new directories
For each patient, check if the registration is done or not. If the 
reagistration is done, take the register label and the RD of the patient and 
copy it in the new directories

for sujb_path in valid_subject_dirs:
    suj_id= os.path.basename(sujb_path)
    Register= Suj[["Recalage","anonym_nom"]]
    Subj_Recalage = Register.Recalage[Register.anonym_nom == suj_id].values[0]
    if Subj_Recalage == 1.0:
        label = os.path.join(Path_label, suj_id,"labels_to_rd.nii.gz")
        rd = os.path.join(Path_rd, suj_id, "rd", suj_id + ".nii")
        NewDirlabel = shutil.copy(label, os.path.join(Path, suj_id, "labels_to_rd.nii.gz"))
        NewDirRD = shutil.copy(rd, os.path.join(Path, suj_id, suj_id +".nii"))    
    else:
        continue

###############################################################
# Plot the dose 

def factorplot (df, X, Y, Hue, Col, Data, path):
    
    sns.set(style="whitegrid")
    g = sns.factorplot(x=X , y=Y, hue=Hue, col=Col, data=Data,
                       size=6, aspect=.75)
    g.despine(left=True)
    plt.show()
    plt.save(path)

"""

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

   
def get_roi_mask(roi_filename, label_number=1):
    mask = nib.load(roi_filename).get_data()
    mask = mask.reshape(mask.shape[0:3])
    mask = mask == label_number
    return mask


def extract_data_in_mask(imgs, mask):
    """ Intersect the imgs with the specified ROI

    The current code assumes that mask and the imgs have the same shape!

    Parameters
    ----------
    imgs : is a list of 3D img files,
    roi_filename : is the name of an img file containing the ROI mask

    Returns
    -------
    data : array with imgs information in the given ROI
    """
    if type(imgs) is not list:
        imgs = [imgs]

    data = np.zeros([mask.shape[0], mask.shape[1], mask.shape[2], len(imgs)])
    for i in range(len(imgs)):
        assert os.path.isfile(imgs[i])
        img = nib.load(imgs[i]).get_data()
        img = img.reshape(img.shape[:3])
        assert img.shape == mask.shape
        data[:, :, :, i] = img

    return data[mask, :]

if __name__ == "__main__":


    # Global parameters
    Path = "/home/edogerde/Bureau/extraction_roi"
    xml_path = "/home/edogerde/Bureau/delineation_space/lpba40.label.xml"
    
    # Create a CSV file

    columns = ["Patients","label"]
    measures = ["min", "max", "mean", "median", "std"]

    df = pd.DataFrame(columns= columns + measures)

    Patients = []
    Label=[]
    Min=[]
    Max=[]
    Mean=[]
    Median=[]
    std=[]
     
   # get all the subject in the directory
   
    
    valid_subject_dirs = [os.path.join(Path, dir_name)
                          for dir_name in os.listdir(Path)
                          if os.path.isdir(os.path.join(Path, dir_name))]

    # Extract the dose from the roi and put it in a list
    for sujb_path in valid_subject_dirs:
        suj_id= os.path.basename(sujb_path)
        print(suj_id)
        rd = os.path.join(Path, suj_id, suj_id +".nii")
        print("rd OK")
        if os.path.isfile(rd):
            
            
            label = os.path.join(Path, suj_id, "labels_to_rd.nii.gz")

            name, rank, roi_label = read_roiLabel(xml_path)
            
            for labelNumber in roi_label:
                mask = get_roi_mask(label, label_number = labelNumber)
                data = extract_data_in_mask(rd, mask)
                print(suj_id,labelNumber)
                a = np.mean(data)
                b = np.median(data)
                c = np.max(data)
                d = np.min(data)
                e = np.std(data)
                Patients.append(suj_id)
                Label.append(labelNumber)
                Min.append(d)
                Max.append(c)
                Mean.append(a)
                Median.append(b)
                std.append(e)
                print("done")
        else:
            name, rank, roi_label = read_roiLabel(xml_path)
            for labelNumber in roi_label:
                Patients.append(suj_id)
                Label.append(labelNumber)
                Min.append(0)
                Max.append(0)
                Mean.append(0)
                Median.append(0)
                std.append(0)

            print("no file for"  +  suj_id)
            continue
    
    # Organize the data in a dataframe 
    df["Patients"]=Patients
    df["min"]=Min
    df["max"]=Max
    df["median"]=Median
    df["std"]=std
    df["median"]=Median
    df["mean"]=Mean
    df["label"]=Label
    df.to_csv("/home/edogerde/Bureau/Roi_RC_C_R.csv")

   # dfAllRoi= pd.read_csv("/home/edogerde/Bureau/Roi_RC_C_R.csv")
"""    
    # TEST : Plot the dose 
    df_hippocampeR = dfAllRoi[(df.label == 89) &(df.label == 90)][["mean"]]
    sns.set(style="whitegrid")
    g = sns.factorplot(y="mean", x= "Patients", hue="label", data=df_hippocampeR,
                       size=6, aspect=.75)
    g.despine(left=True)
    plt.show()

    sns.barplot(y="mean", x="label", data= df_hippocampeR)   
    plt.show()
"""
"""    
Plotting a diagonal correlation matrix
======================================

"""
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(style="white")

# Generate a large random dataset
d = dfAllRoi[["label", "Patients","mean" ]]
dt=d.pivot(columns="label", index = "Patients" , values= "mean")

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
"""

"""
# custom way
mask = get_roi_mask(path, label)
data = extract_data_in_mask(radio_img, mask)
np.max(data)
np.mean(data)

# nilearn way
masker = NiftiMasker(mask)
data = masker.fit_transform()
result_img = masker.inverse_transform(results)
nib.save(result_img, path)
"""
