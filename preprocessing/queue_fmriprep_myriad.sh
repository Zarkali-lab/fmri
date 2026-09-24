#!/bin/bash

d0=$PWD

##### EDIT FROM HERE 

# e.g. /home/user/Scratch/MMPD/MRI_DATA/raw_data/3T/
rawDataDir=/home/zcqsgth/Scratch/MMPD/MRI_DATA/raw_data_fmriprep/3T/

#### TO HERE


# hold for job id

hold_jid=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -hold_jid)
            hold_jid="$2"
            shift 2
            ;;
        *)
            break
            ;;
    esac
done

# subject inputs

if [ $# -eq 0 ]; then
        subjects=($(find $rawDataDir -maxdepth 1 -type d -name "sub-*" | \
		xargs -n 1 basename | \
		sed 's/^sub-//' | \
		sort))
else
        subjects=($@)
fi

numSubs=${#subjects[@]}

echo $numSubs

# qsub cmd

qsub_cmd=(qsub -wd "$d0" -t 1-"$numSubs")

if [ -n "$hold_jid" ]; then
    qsub_cmd+=(-hold_jid "$hold_jid")
fi

qsub_cmd+=(fmriprep_myriad.sh "${subjects[@]}")

echo "${qsub_cmd[@]}"
"${qsub_cmd[@]}"

