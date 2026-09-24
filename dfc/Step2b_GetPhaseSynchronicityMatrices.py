###################################################### Phase & Synchronicity #######################################

# A.Zarkali | 8/11/2024
# Script to get phase and synchronicity matrices for subsequent Kuramoto and Leidan Analysis
# INPUTS: 
## -i: An np.array with shape=(subjects, volumes, ROIs)
## -subs: Optimal number of clusters
## -o: Output folder (absolute path prefered)
## -tr: the TR of the acquisition

# OUTPUTS in the output_folder:
### Stacked subjects, array
## 1: For LEIDA
### a. Phases array for all parcels/voxels
### b. Synchronicity matrix for all parcels/voxels
### c. Leading eigenvector of synchornicity matrix
## 2: For Kuramoto model
### a. Phases array for all parcels/voxels
### b. Synchronicity matrix for all parcels/voxels
### c. Leading eigenvector of synchornicity matrix

import scipy as sc
import scipy.io
import numpy as np
import dynfc as dyn
import pickle
import pandas as pd
from pathlib import Path
import argparse

# Set paths and constants
parser=argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="np array of timeseries, shape=(subjects, volumes, number of ROIs in parcellation")
parser.add_argument("-subs", type=str, help="path to list of participants")
parser.add_argument("-tr", type=float, help="TR of the acquisition")
parser.add_argument("-o", type=str, help="path output_folder")
args=parser.parse_args()

TR=args.tr
timeseries_path = Path(args.i)
timeseries= np.load(timeseries_path)
timeseries = timeseries.transpose(1, 2, 0)
path_list = Path(args.subs)
subject_list = np.loadtxt(path_list, dtype=str)
output_folder = Path(args.o)

lthr = 0.02
uthr = 0.1
k=2

# Get Phase and Synchronicity Matrices for LEIDA
Phases, SyncConn, Leida = dyn.run_multiPatLEiDA(timeseries, TR, lthr, uthr, k)
# Transpose all results to make sure outputs have shape: subjects, volumes, ROIs
Phases = Phases.transpose(2,1,0)
SyncConn = SyncConn.transpose(3,2,0,1)
Leida = Leida.transpose(2,0,1)
np.savez_compressed((output_folder / "LEIDA_Phases.npz"),Phases=Phases,ids=subject_list)
np.savez_compressed((output_folder / "LEIDA_SynchronicityMatrices.npz"),FCmatrices=SyncConn,ids=subject_list)
np.savez_compressed((output_folder / "LEIDA_LeadingEigenvectors.npz"),Eigenvectors=Leida,ids=subject_list)

# Get Phase and Synchronicity Matrices for Kuramoto
Metastability, KOPSyncConn, Entropy = dyn.run_multiPatKOP(timeseries, k, 8, TR, lthr, uthr, k)
## Collate
output = pd.DataFrame(columns=["id", "metastability", "entropy"], data=np.zeros((len(subject_list),3)))
output.id = subject_list
output.metastability = Metastability
output.entropy = Entropy
out_path = output_folder / "KOP.csv"
output.to_csv(out_path)

KOPSyncConn = KOPSyncConn.transpose(1,0)
np.savez_compressed((output_folder / "KOP_SynchronicityMatrices.npz"),FCmatrices=KOPSyncConn,ids=subject_list)
