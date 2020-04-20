#!/bin/bash

# loop over folders
for f in test*
do
    # move into the folder and check the ICs folder
    cd $f
    
    num_files=`ls ICS | wc -l`
    echo "checking genic IC files: $f ...$num_files"

    # move back to base folder
    cd ..
done
