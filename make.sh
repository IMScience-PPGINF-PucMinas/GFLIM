#! /bin/bash

# install pyIFT (required for build GFLIM)
pipenv install
pipenv shell
cd PyIFT/ ; chmod +x install-pyift.sh; cd ..
export NEWIFT_DIR=${PWD}/
cd PyIFT/ ; ./install-pyift.sh; cd ..; exit

# compile GFLIM
mkdir build; cd build; cmake ..; make; cd ..;


