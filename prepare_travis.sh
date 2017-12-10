#!/bin/bash
#Script to fetch and build dependencies for the tests, for travis.

wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b -p $HOME/miniconda
export PATH=$HOME/miniconda/bin:$PATH
conda update --yes conda
conda create --yes -n test python=3.4
source $HOME/miniconda/bin/activate test
conda install --yes -c gfortran_linux-64 bccp nbodykit matplotlib numpy scipy configobj bigfile nose

mkdir tests
mkdir depends
cd depends
##Get and make CAMB.
git clone https://github.com/sbird/camb.git
cd camb
cd pycamb
python3 setup.py install --user

#Get and make N-GenIC.
cd ../
git clone https://github.com/sbird/S-GenIC.git
cd S-GenIC
git submodule update --init --recursive
make

#Get and make GenPK.
cd ../
git clone https://github.com/sbird/GenPK.git
cd GenPK
git submodule update --init --recursive
make
