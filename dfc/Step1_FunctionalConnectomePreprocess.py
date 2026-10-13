####################################################### Functional Connectome Preprocessing #######################################

# A.Zarkali | 4/11/2024
# Script to derive timeseries from functional data
# Assumes that data have been preprocessed using fmriprep

# INPUTS: 
## -i data_folder (absolute path prefered) 
## -o: output_folder (absolute path prefered)
## -atlas: atlas path (nifti image in MNIspace)
## -n: number of regions in the parcellation
## -subs: subject_list (txt file of included subjects, assumes BIDS compatible format and no "sub" pre-fix)
### EXAMPLE HARD-CODED SHARED-DRIVE PATH (documentation only): replace for your site;
### actual source folder is provided with -i.
### For windows paths to shared drives format is: '//rdp.arc.ucl.ac.uk/ritd-ag-project-rd00hy-rswei70/DLB/FMRI/fmriprep' or include r before path

# OUTPUTS in the output_folder:
## 1: An np.array with shape=(subjects, volumes, ROIs)
## 2: A connectome matrix, raw and normalised (fisher transformed for each participant)

# Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats as stats
from scipy.stats import norm
from scipy.stats import shapiro
from scipy.stats import percentileofscore
import statsmodels.api as sm
from statsmodels.formula.api import ols
import seaborn as sns
from pathlib import Path
import bct as bct
import networkx as nx
from sklearn.covariance import EmpiricalCovariance as EmpiricalCovariance
from nilearn.connectome import ConnectivityMeasure
import nibabel
import nilearn
from nilearn import datasets, plotting, input_data, signal, image # datasets: for fetching atlas
from nilearn.input_data import NiftiLabelsMasker
import argparse
import os

# Set up argument parsing
parser=argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="path to data folder")
parser.add_argument("-o", type=str, help="path to output folder")
parser.add_argument("-atlas", type=str, help="path to atlas nifti file")
parser.add_argument("-n", type=int, help="number of regions in parcellation")
parser.add_argument("-subs", type=str, help="path to list of participants")
args=parser.parse_args()

# # Set parameters
data_folder = Path(args.i)
output_folder = Path(args.o)
atlas_filename = Path(args.atlas)
num_regions = args.n
path_list = Path(args.subs)
subject_list = np.loadtxt(path_list, dtype=str)
atlas_filename = image.load_img(atlas_filename, wildcards=True)

# Extract Motion confounds
# ## Specify the confounds I want to extract: WM and CSF signal + 6 motion parameters
# Cols=["trans_x", "trans_x_derivative1", "trans_x_derivative1_power2", "trans_x_power2", "trans_y", "trans_y_derivative1", "trans_y_power2", "trans_y_derivative1_power2", "trans_z", "trans_z_derivative1", "trans_z_derivative1_power2", "trans_z_power2", "rot_x", "rot_x_derivative1", "rot_x_derivative1_power2", "rot_x_power2", "rot_y", "rot_y_derivative1", "rot_y_derivative1_power2", "rot_y_power2", "rot_z", "rot_z_derivative1", "rot_z_derivative1_power2", "rot_z_power2", "csf", "csf_derivative1", "csf_derivative1_power2", "csf_power2", "white_matter", "white_matter_derivative1", "white_matter_derivative1_power2", "white_matter_power2", "global_signal", "global_signal_derivative1", "global_signal_power2", "global_signal_derivative1_power2"]
# ## Loop through all subjects and generate new confound in csv format
# for i in range(len(subject_list)):
#     # tsv = "sub-%s_task-rest_desc-confounds_regressors.tsv" %(subject_list[i]) # name of each confound tsv file
#     tsv = "sub-%s_task-rest_desc-confounds_timeseries.tsv" %(subject_list[i])
#     path = data_folder / str("sub-" + str(subject_list[i])) / "func" / tsv # path to the tsv file 
#     df = pd.read_csv(path, sep='\t') # read tsv file
#     df = df[Cols] # select the subset of dataset
#     df = df.replace([np.inf, np.nan], 0)
#     export_path = data_folder / str("sub-" + subject_list[i]) / "func" / str("sub-" + subject_list[i] + "_motion.csv") # clarify export path
#     df.to_csv(export_path, header=True, index=False) # save csv file with extracted confounds

# Load lists of fmris and confounds
fmri_list = []
confounds_list = []
## Loop through list of subjects and append lists
for i in range(len(subject_list)):
    mri = "sub-%s_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz" %(subject_list[i]) # name of the preprocessed fmri file
    # mri = "sub-%s_task-rest_bold.nii" %(subject_list[i])
    path_mri = data_folder / str("sub-" + str(subject_list[i])) / "func" / mri # path to the fmri file (this need to be inputed as a windows path)
    fmri_list.append(path_mri) # append path to the fmri list
    path_tsv = str(data_folder / str("sub-" + str(subject_list[i])) / "func" / str("sub-" + str(subject_list[i]) + "_motion.csv")) # path to the confound file but as a string not a path (!)
    confounds_list.append(path_tsv) # append path to list

# Load Parcellation 
num = [range(1,num_regions+1)] # its 1 rather than 0 numbered
labelsNum = np.array(num)
time_series_All = []
masker = NiftiLabelsMasker(labels_img=atlas_filename, standardize=True, low_pass=0.08, high_pass=0.01, t_r= 0.7, resampling_target="data") # include bandpass filter 0.01 to 0.08
for sub in range(len(subject_list)):
    print(subject_list[sub])
    ## load confounds file
    confounds = np.loadtxt(confounds_list[sub], delimiter=",",skiprows=1)
    ## fit timeseries
    time_series = masker.fit_transform(fmri_list[sub], confounds = confounds)
    ## discard dummy volumes (first 4 volumes)
    time_series = time_series[3:,]
    ## append timeseries to list
    time_series_All.append(time_series)
# Stack all timeseries to an Np array
allts = np.stack(time_series_All)
# Save array
array_path = output_folder / str("Timeseries_array_" + str(num_regions) + ".npy")
np.save(array_path, allts)

# Save for melodic 
S, W, R = allts.shape
ts2d = allts.reshape(S * W, R)  # shape (timepoints, ROIs)
ts_path = output_folder / str("Timeseries_array_" + str(num_regions) + "melocid.txt")
np.savetxt('melodic_input.txt', ts2d)

correlation_measure = ConnectivityMeasure(cov_estimator=EmpiricalCovariance(), kind="correlation", standardize="zscore_sample")
## Compute correlation matrix
for sub in range(len(subject_list)):
    correlation_matrix = correlation_measure.fit_transform([allts[sub]])[0]
    correlation_matrix = np.nan_to_num(correlation_matrix)
    correlation_matrix_norm = np.arctanh(correlation_matrix) ## Fisher transform
    # Save correlation matrix
    out_filename = str("sub-" + subject_list[sub] + "_shaefer232.csv")
    out_path = output_folder/ out_filename
    out_filename_norm = str("sub-" + subject_list[sub] + "_shaefer232_norm.csv")
    out_path_norm = output_folder / out_filename_norm
    np.savetxt(out_path, correlation_matrix, delimiter=",")
    np.savetxt(out_path_norm, correlation_matrix_norm, delimiter=",")
