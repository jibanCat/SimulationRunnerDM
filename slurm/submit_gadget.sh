#!/bin/bash

# loop over folders
for f in test*
do  
    # move into the folder and submit
    cd $f

    last_powerspec=`ls output | grep power | tail -n 1`
    if [ $last_powerspec == 'powerspectrum-1.0000.txt' ]
    then
        echo "checking last powerspec file: $f ...$last_powerspec"
        echo " ... submitting gadget files: $f"
        sbatch mpi_submit
    fi

    # move back to base folder
    cd ..
done
