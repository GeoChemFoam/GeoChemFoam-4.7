#!/bin/bash

###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation
TotalTime=2000

# Define flow rate
flowrate=3.41e-8

## Define the kinematic viscocity of the fluid (m^2/s)
##(e.g for water this is 1e-6, for air this would be 1.478e-5)
Visc=1e-06

## Number of processor
NP=24

#### END OF USER INPUT

cp constant/transportPropertiesSP constant/transportProperties
cp system/fvSolutionSP system/fvSolution
cp system/fvSchemesSP system/fvSchemes
cp system/controlDictSP system/controlDict
sed -i "s/Visc/$Visc/g" constant/transportProperties
sed -i "s/TotalTime/$TotalTime/g" system/controlDict

# Load user environment variables 
source ~/.bashrc

source $HOME/works/GeoChemFoam-dev/etc/bashrc

set -e

MPIRUN=mpirun 

cp -r 0_org_SP 0

sed -i "s/flowrate/$flowrate/g" 0/U

# Decompose
echo -e "DecomposePar"
decomposePar > decomposeParSP.out

# Run simpleFoam in parallel
echo -e "Run simpleFoam in parallel"
$MPIRUN -np $NP simpleFoam -parallel  > simpleFoamSP.out

# ReconstructPar
echo -e "reconstructPar"
reconstructPar -latestTime > reconstructParSP.out

rm -rf 0
rm -rf processor*
echo -e "processPoroPerm"
processPoroPerm > poroPermSP.out

