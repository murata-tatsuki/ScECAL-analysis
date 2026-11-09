#!/bin/bash

# Geant4 10.7.2  ROOT 6.24.06
#source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

# Geant4 11.0.3  ROOT 6.26.08
# source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/setup.sh

# Geant4 11.1.1  ROOT 6.28.00
#source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc11-opt/setup.sh

source /meg/home/meg/pro/cern/setenv.sh
export CMAKE_INSTALL_PREFIX=/megraid01/users/data_beamtest/simulation/SimCalModule:$CMAKE_INSTALL_PREFIX
