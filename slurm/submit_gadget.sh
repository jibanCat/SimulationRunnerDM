#!/bin/bash

# loop over folders
for f in test*
do  
    # move into the folder and submit
    cd $f

    output_size=`ls -lh output | wc -l`
    if [ $output_size -eq 1 ]
    then
        echo "checking output folder size: $f/output ...$output_size"
        echo " ... submitting gadget files: $f"
        sbatch mpi_submit
    fi

    # move back to base folder
    cd ..
done
