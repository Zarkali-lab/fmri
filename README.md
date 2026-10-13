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

## Participant lists and data paths

The `dfc` Python steps accepting `-subs` already load participant identifiers from an external text file. Create a local `participants.txt` with one participant pseudocode per line (without `sub-` prefix), then pass `-subs participants.txt` to each relevant step. Do not commit this file.

Input and output file paths for `dfc` scripts are supplied via command-line arguments such as `-i`, `-o` and `-atlas`, not embedded site-specific directories.

The `preprocessing/*myriad.sh` scripts use **hard-coded Myriad paths** marked `HARD-CODED MYRIAD PATH`. Replace them for your cluster account before running. Existing path strings have deliberately been left unchanged. The fMRIPrep script's quoted `~` log paths also need correction before use.
