#!/bin/bash

###### USERS INPUT ############################################################

## Define the total time of the simulation and how often to output concentration field

TotalTime=900
WriteTimestep=150
runTimestep=1

## Define the diffusion coefficient of the species as a solute (m^2/s)
## Eg for water vapor in air this would be 2.42e-5, for CO2(aq) in water this would be 3e-9
Diff=1e-8

## Number of processor
NP=24

#### END OF USER INPUT #######################################################

cp system/controlDictTransport system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict
sed -i "s/runTimestep/$runTimestep/g" system/controlDict


cp constant/thermoPhysicalProperties1 constant/thermoPhysicalProperties
sed -i "s/Diff/$Diff/g" constant/thermoPhysicalProperties


# Load user environment variables 
source ~/.bashrc

source $HOME/works/GeoChemFoam-4.7/etc/bashrc

set -e

MPIRUN=mpirun 

#rm -rf *.out

cp 0_orig/Species 0/.
rm -rf 0/uniform

cp system/decomposeParDict1 system/decomposeParDict
sed -i "s/NP/$NP/g" system/decomposeParDict

# Decompose
echo -e "DecomposePar"
decomposePar > decomposeParTransport.out

# Run multiSpeciesTransportFoam in parallel
echo -e "Run multiSpeciesTransportFoam in parallel"
$MPIRUN -np $NP multiSpeciesTransportFoam -parallel  > multiSpeciesTransport.out

# ReconstructPar
echo -e "reconstructPar"
reconstructPar > reconstructParTransport.out

rm -rf processor*

echo "process concentration"
processConcentration > processConcentrationTransport.out
