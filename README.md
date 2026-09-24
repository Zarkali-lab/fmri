# fmri
Contains fmri pipelines and tools commonly used in the lab.

## dfc
The folder contains the analysis pipeline to perform different dynamic functional connectivity analyses at whole connectome level (not GICA). 

#### Dependencies:
- Conda environment with installed packages (incl versions) are listed on environment.yml 
- dfc.py (contains functions required in the steps)

#### For all pipelines:
Each step is run in the command line
List of inputs in each script but also can be seen by running: 
$ python <script>.py --help

**STEP 1:**
- Preprocessing: extract timeseries from each participant that has

##### ANALYSIS A) Sliding window and K-means clustering
**STEP 2a:**
- Derive sliding windows (default 44sec - editable in script)
- Calculate elbow, silhouette and Davies Bouldin score in order to find the optimal K for future clustering
**STEP 3a:**
- Perform k-means clustering with optimal K
- Outputs a csv with dwell, frequency and transitions

##### ANALYSIS B) LEIDA
**STEP 2b:**
- Get Phase and Synchronicity Matrices
