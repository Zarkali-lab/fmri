####################################################### Perform Clustering #######################################

# A.Zarkali | 5/11/2024
# Script to perform kmeans clustering on timeseries
# INPUTS: 
## i: An np.array of sliding windows FC matrices with shape=(subjects, number of sliding windows, volumes*ROIs)
## n: Optimal number of clusters
## o: Output folder (absolute path prefered)

# OUTPUTS in the output_folder:
## 1: An np.array of the resulting metastates and the allocation of each timewindow to each state for each participant
## 2: A csv with dwell, frequency and total transitions per participant

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats as stats
from scipy.stats import norm, shapiro, percentileofscore
from pathlib import Path
from sklearn.metrics import silhouette_score
from sklearn.covariance import EmpiricalCovariance
from sklearn.cluster import KMeans
from nilearn.connectome import vec_to_sym_matrix, sym_matrix_to_vec
from nilearn.plotting import plot_matrix
import seaborn as sns
from dfc import *
import argparse

# Set paths and constants
parser=argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="sliding windows FC matrices with shape=(subjects, number of sliding windows, volumes*ROIs")
parser.add_argument("-n", type=int, help="number of clusters to derive")
parser.add_argument("-pipeline", type=str, help="'window' for sliding windows, and 'leida' for LEIDA")
parser.add_argument("-o", type=str, help="path output_folder")
parser.add_argument("-rois", type=int, help="number of rois in the parcellation")
args=parser.parse_args()

filename = Path(args.i)
output_folder = Path(args.o)
pipeline = args.pipeline
optimalK= args.n
rois=args.rois
output_folder = Path(args.o)

plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
cmap = "jet"

# Load data
if pipeline == "window": 
    dFC = np.load(filename)
    dFC_vector = dFC["FCmatrices"]
    # Reshape for analysis
    cluster = dFC_vector.reshape((dFC_vector.shape[0]*dFC_vector.shape[1],dFC_vector.shape[2]))
    # Calculate metrics
elif pipeline == "leida":
    dFC = np.load(filename)
    dFC_vector = dFC["Eigenvectors"]
    cluster = dFC_vector.reshape((dFC_vector.shape[0]*dFC_vector.shape[1],dFC_vector.shape[2]))
else:
    print("incorrect pipeline input")
    quit()

# Perform clustering
reassignement = 0.01
batch_size = 50
kclust = KMeans(n_clusters=optimalK,init='k-means++')
kclust.fit(cluster)

# Save clustering results
if pipeline == "window":
    centroids = np.ndarray((kclust.n_clusters,rois,rois))
    for i in range(centroids.shape[0]):
        centroids[i] = vec_to_sym_matrix(kclust.cluster_centers_[i],diagonal=np.zeros(rois))
    clustlabels = kclust.labels_.reshape((dFC_vector.shape[0],dFC_vector.shape[1]))
    filename = output_folder/str("kmeans" + str(optimalK) + "_results.npz")
    np.savez_compressed(filename,labels = clustlabels,centroids = centroids)
elif pipeline == "leida":
    centroids = np.ndarray((kclust.n_clusters,rois))
    for i in range(centroids.shape[0]):
        centroids[i] = kclust.cluster_centers_[i]
    clustlabels = kclust.labels_.reshape((dFC_vector.shape[0], dFC_vector.shape[1]))
    filename = output_folder/str("kmeans" + str(optimalK) + "_results.npz")
    np.savez_compressed(filename,labels = clustlabels,centroids = centroids)
else:
    print("incorrect pipeline input")
    quit()

# Get metrics
kinds= ['correlation','partial correlation','tangent','covariance']
### returns dwel, frequency, mean dwell, mean freq, transitions
def freq_dwell(labels,n_states=optimalK):
    n_subj=labels.shape[0]
    counter=np.zeros((n_subj,n_states))
    trans=np.ones((n_subj,n_states))
    dwell=np.zeros((n_subj,n_states))
    freq=np.zeros((n_subj,n_states))
    for i in range(n_subj):
        val_old=0
        for  val in labels[i]:
            for j in np.arange(0,n_states):
                if val == j:
                    counter[i,j] = counter[i,j] +1          
                if val != val_old and val==j:
                    trans[i,j] += 1
            val_old=val
    dwell = counter/trans
    freq = counter*100/labels[0].shape
    return dwell,freq,trans

# Save output
d,f,t = freq_dwell(clustlabels, optimalK)
df = pd.DataFrame(data=d)
df.columns = [f"Dwell{i+1}" for i in range(optimalK)]
for n in range(optimalK):
	df[str("Freq" + str(n))] = f[:,n] # Frequency of State 1
df["Transitions"] = np.sum(t,axis=1)
out_path = output_folder / str(str(pipeline) + "_" + str(optimalK)  +  "_FrequencyDwell.csv")
df.to_csv(out_path)

# Save plot
if pipeline == "window":
    plt.figure(figsize=(20,20))
    fig, axes = plt.subplots(1, optimalK, figsize=(12, 12))
    for i in range(centroids.shape[0]):
        ax = plt.subplot(3,2,i+1)
        sns.heatmap(centroids[i], cmap=cmap)
        ax.set_title(str("Metastate" + str(i+1)))
    plt.savefig(output_folder / str(str(pipeline) + "_" + str(optimalK)  +  "_Metastates.png"))
else:
    quit()