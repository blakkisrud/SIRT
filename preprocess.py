import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import nrrd
from scipy.io import savemat

# Script to perform pre-processing of PET/CT and voxel-dumps from PMOD to be used 
# for voxel-dosimetry

# Input values:
# 

pet_img_path = "" # Path to the pet-img
voxel_data = "" # Path to the voxel dump
output_matlab_path = "output_to_matlab.mat" # Path to the output

# First load the PET-image. This have been saved from before as an nrrd-file

pet_img, pet_header = nrrd.read(pet_img_path, index_order="F")

pet_img = pet_img*1e-3 # Convert to kBq/ml

# The excel-reader from the pandas library is used. It might not be the most
# intended use of the library, but it works
pd_data = pd.read_excel(voxel_data)

print("Dataset contains " + str(len(np.unique(pd_data["VoiName(Region) [string]"]))) + " number of VOIs")

voi_names = np.unique(pd_data["VoiName(Region) [string]"])
labels = np.arange(len(voi_names))+1

name_label_dict = dict(zip(voi_names, labels))

label_map = np.zeros(pet_img.shape) # This matrix contains the labels

print(voi_names)

# Get the values and place them in the labe_map matrix

for row in pd_data.iterrows():
    label = name_label_dict[row[1]['VoiName(Region) [string]']]
    x = row[1]["X [pixel]"]-1
    y = row[1]["Y [pixel]"]-1
    z = row[1]["Z [pixel]"]-1
    label_map[x,y,z] = label

# Save the resulting matrices as matlab-files to use further

# Construct dictionary

output_dict = {"pet_img": pet_img, "label_matrix" : label_map}

savemat(output_matlab_path , output_dict)



