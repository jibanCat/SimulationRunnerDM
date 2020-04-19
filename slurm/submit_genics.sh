#!/bin/bash

# loop over folders
for f in test*
do
    echo "submitting genic files: $f ..."

    # move into the folder and submit
    cd $f
    sbatch mpi_submit_genic

    # move back to base folder
    cd ..
done
