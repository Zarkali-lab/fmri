####################################################### Functional Connectome Preprocessing #######################################

# A.Zarkali | 4/11/2024
# Script to derive metrics to find optimal k means clustering number from a timeseries
# INPUTS: 
## -i: An np.array of things to cluster:
#### if sliding windows this should FC matrices with shape=(subjects, number of sliding windows, volumes*ROIs)
#### if leida this should the leading eigenvectors with shape=(subject, number of windows, eigenvectors)
## -subs: subject_list (txt file of included subjects, assumes BIDS compatible format and no "sub" pre-fix)
## -pipeline: either 'window' for sliding windows, and 'leida' for LEIDA
## -o: output_folder (absolute path prefered)

# OUTPUTS in the output_folder:
## 1: An np.array of sliding windows FC matrices with shape=(subjects, number of sliding windows, volumes*ROIs)
## 2: A csv file with clustering statistics saved in the output folder
## 3: Images for each clustering statistic saved in the output folder

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
from dfc import *
import argparse

# Set paths and constants
parser=argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="np.array of stacked inputs for clustering of shape: (subjects, number of windows, volumes*ROIS) for sliding window, (subjects, leading eigenvector) for LEIDA ")
parser.add_argument("-pipeline", type=str, help="'window' for sliding windows, and 'leida' for LEIDA")
parser.add_argument("-o", type=str, help="path to output folder")
parser.add_argument("-subs", type=str, help="path to list of participants")
args=parser.parse_args()

filename = Path(args.i)
output_folder = Path(args.o)
list_path = Path(args.subs)
subject_list = np.loadtxt(list_path, dtype=str)
pipeline = args.pipeline

plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

if pipeline == "window": 
    dFC = np.load(filename)
    dFC_vector_fisher = dFC["FCmatrices"]
    # Reshape for analysis
    cluster = dFC_vector_fisher.reshape((dFC_vector_fisher.shape[0]*dFC_vector_fisher.shape[1],dFC_vector_fisher.shape[2]))
    # Calculate metrics
elif pipeline == "leida":
    leida = np.load(filename)
    leida_vector = leida["Eigenvectors"]
    cluster = leida_vector.reshape((leida_vector.shape[0]*leida_vector.shape[1],leida_vector.shape[2]))
else:
    print("incorrect pipeline input")
    quit()

## Elbow criterion
print("Calculating elbow criterion")
elbow = elbow(cluster)
K = range(1,10)
fig=plt.figure(figsize=(10,10))
plt.plot(K, elbow, linewidth=2, c="black", marker="o")
plt.xlabel('Number of clusters')
plt.ylabel('Sum of Squared Distances')
plt.title('Elbow Method For Optimal k')
fig_path = output_folder / str(pipeline + "Elbow.png")
plt.savefig(fig_path)
plt.close()

## Silhouette criterion
print("Calculating silhouette score")
silhouette = silhouette(cluster)
K = range(2, 10)
fig=plt.figure(figsize=(10,10))
plt.plot(K, silhouette, linewidth=2, c="black", marker="o")
plt.ylabel("Silhouette", fontsize=15)
plt.xlabel("k", fontsize=15)
plt.title("Silhouette for K-means cell's behaviour")
fig_path = output_folder / str(pipeline + "Silhouette.png")
plt.savefig(fig_path)
plt.close()

## Davies Bouldin
print("Calculating Davies Bouldin score")
davies_bouldin = davies_bouldin(cluster)
fig=plt.figure(figsize=(10,10))
K = range(2,10)
plt.plot(K, davies_bouldin, linewidth=2, c="black", marker="o")
plt.grid(False)
plt.xlabel('Number of States')
plt.ylabel('Davies Bouldin Score')
plt.title('Davies Bouldin Score by Cluster Count')
fig_path = output_folder /  str(pipeline + "DaviesBouldin.png")
plt.savefig(fig_path)
plt.close()

## Calinski Harabasz
print("Calculating Calinski Harabasz score")
calinski = calinski(cluster)
fig=plt.figure(figsize=(10,10))
K = range(2,10)
plt.plot(K, calinski, linewidth=2, c="black", marker="o")
plt.grid(False)
plt.xlabel('Number of States')
plt.ylabel('Calinski Harabasz  Score')
plt.title('Calinski Harabasz  Score by Cluster Count')
fig_path = output_folder /  str(pipeline + "Calinski.png")
plt.savefig(fig_path)
plt.close()

## Collate all scores
output = pd.DataFrame(columns=["elbow", "silhouette", "davies_bouldin", "calinski"], data=np.zeros((9,4)))
output.elbow = elbow
output.silhouette = np.concatenate([["NA"], silhouette], axis=0)
output.davies_bouldin = np.concatenate([["NA"], davies_bouldin], axis=0)
output.calinski = np.concatenate([["NA"], calinski], axis=0)
out_path = output_folder /  str(pipeline + "OptimalClusteringScores.csv") 
output.to_csv(out_path)