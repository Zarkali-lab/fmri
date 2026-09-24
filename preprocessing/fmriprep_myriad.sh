#$ -l h_rt=12:00:00

#$ -l mem=3.8G

# Job name
#$ -N fmriprep

#$ -cwd 

#$ -o /dev/null

#$ -e /dev/null

#$ -pe smp 4

subjects=($@)

idx=$((${SGE_TASK_ID}-1))

subj=${subjects[${idx}]}

bidsDir=~/Scratch/MMPD/MRI_DATA/raw_data_fmriprep/3T
derivDir=~/Scratch/MMPD/MRI_DATA/derivatives/3T
FSlicense=~/Scratch/FS_license.txt

FMRIPREP=~/Scratch/fmriprep-23.2.0.simg


mkdir -p ~/Scratch/MMPD/MRI_DATA/derivatives/3T/sub-$subj/ses-01/func/

exec >"~Scratch/MMPD/MRI_DATA/derivatives/3T/sub-$subj/ses-01/func/log_fmriprep_hpc_${JOB_ID}.${SGE_TASK_ID}.out" \
	2>"~/Scratch/MMPD/MRI_DATA/derivatives/3T/sub-$subj/ses-01/func/log_fmriprep_hpc_${JOB_ID}.${SGE_TASK_ID}.err"


mkdir -p $derivDir/fmriprep_work

apptainer run --cleanenv \
	-B $bidsDir:/data:ro \
	-B $derivDir:/out \
	-B $derivDir/fmriprep_work:/work \
	-B $(dirname $FSlicense):$(dirname $FSlicense):ro \
	$FMRIPREP \
	/data /out participant \
	--participant-label $subj \
	--fs-license-file $FSlicense \
	--work-dir /work \
	--nprocs $NSLOTS

