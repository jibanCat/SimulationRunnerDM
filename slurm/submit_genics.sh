#!/bin/bash

# loop over folders
for f in test*
do
    # move into the folder and submit
    cd $f
    num_files=`ls ICS | wc -l`
    if [ $num_files -eq 0 ]
    then
        echo "submitting genic IC files: $f"
        sbatch mpi_submit_genic
    fi
    # move back to base folder
    cd ..
done
