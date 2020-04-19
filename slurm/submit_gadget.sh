#!/bin/bash

# loop over folders
for f in test*
do
    echo "submitting gadget files: $f ..."
    
    # move into the folder and submit
    cd $f
    sbatch mpi_submit

    # move back to base folder
    cd ..
done
