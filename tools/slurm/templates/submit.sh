#!/bin/bash

for job in $@
do
    LINK_COUNT=1
    # Check if job is a chain job
    echo $job | gsed '/chain/!{q1}' > /dev/null
    if [[ $? == 0 ]]; then
        LINK_COUNT=$(echo $job | gsed -r 's/.*chain_([0-9]*)\.sh/\1/g')
    fi

    CURRENT_LINK=0
    while [ $CURRENT_LINK -lt $LINK_COUNT ] ; do

    if [ $CURRENT_LINK -eq 0 ] ; then
        SLURM_OPT=""
    else
        # if job is a chain job wait for successful completion of last chain link
        SLURM_OPT="-d afterok:$ID"
    fi
    echo "sbatch --export=LINK=$CURRENT_LINK $SLURM_OPT $job"|tee -a manifest.txt
    ID=$CURRENT_LINK
    # ID=$(sbatch --export=LINK=$CURRENT_LINK $SLURM_OPT $job 2>&1|sed 's/[S,a-z]* //g')

    ## Check if ERROR occured
    if [[ "$jobID" =~ "ERROR" ]] ; then
        echo "---> submission failed!" |tee -a manifest.txt; exit 1
    else
        echo "---> job number = $ID"|tee -a manifest.txt
    fi

   let CURRENT_LINK+=1
done
done
