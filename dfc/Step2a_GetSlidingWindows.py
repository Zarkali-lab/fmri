####################################################### Functional Connectome Preprocessing #######################################

# A.Zarkali | 4/11/2024
# Script to derive metrics to find optimal k means clustering number from a timeseries
# INPUTS: 
## -i: array of timeseries (absolute path prefered) 
## -subs: subject_list (txt file of included subjects, assumes BIDS compatible format and no "sub" pre-fix)
## -o: output_folder (absolute path prefered)

# OUTPUTS in the output_folder:
## 1: An np.array of sliding windows FC matrices with shape=(subjects, number of sliding windows, volumes*ROIs)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats as stats
from scipy.stats import norm, shapiro, percentileofscore
from pathlib import Path
from sklearn.metrics import silhouette_score
from sklearn.covariance import EmpiricalCovariance
from sklearn.cluster import MiniBatchKMeans
from nilearn.connectome import ConnectivityMeasure
from dfc import *
import argparse

# Set paths and constants
window = 44
parser=argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="path to timeseries as npy stacked array")
parser.add_argument("-o", type=str, help="path to output folder")
parser.add_argument("-subs", type=str, help="path to list of participants")
args=parser.parse_args()

fc_path = Path(args.i)
output_folder = Path(args.o)
list_path = Path(args.subs)

# Load data
fc = np.load(fc_path)
subject_list = np.loadtxt(list_path, dtype=str)
dFC,_ = dfc_slid_window(fc,44)
print("Calculing sliding windows")

# Compute FC matrices in vector form
estimator = ConnectivityMeasure(kind='covariance',vectorize=True,discard_diagonal=True)
dFC_vector = []
for sl in dFC:
    dFC_vector.append(estimator.fit_transform(sl))    
dFC_vector = np.stack(dFC_vector)
dFC_vector_fisher = dFC_vector
dFC_vector_fisher = np.nan_to_num(dFC_vector_fisher)
for i in range(dFC_vector.shape[0]):
    dFC_vector_fisher[i] = np.arctanh(dFC_vector_fisher[i])
datapath = output_folder / "FC_matrices_vectors.npz"
np.savez_compressed(datapath,FCmatrices=dFC_vector,ids=subject_list)
print(str("Saved vectors at: " + str(datapath))) ### Message of progress