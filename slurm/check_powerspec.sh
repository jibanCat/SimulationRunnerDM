#!/bin/bash

# loop over folders
for f in test*
do
    # move into the folder and check the ICs folder
    cd $f
    
    last_powerspec=`ls output | grep power | tail -n 1`
    echo "checking last powerspec file: $f ...$last_powerspec"

    # move back to base folder
    cd ..
done
