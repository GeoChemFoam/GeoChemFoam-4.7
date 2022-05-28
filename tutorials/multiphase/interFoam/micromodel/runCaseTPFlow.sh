#!/bin/bash

###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation and how often to output
TotalTime=0.5
WriteTimestep=0.025
runTimestep=2e-5

#Define flow rate
flowrate=3.41e-8

#fluid viscosities (m2/s) 
Visc1=1e-06
Visc2=1.65e-5

#fluid densities (kg/m3)
rho1=1000
rho2=864

#interfacial tension (N/m)
ift=0.03

#### END OF USER INPUT #######################################################

source $HOME/works/GeoChemFoam-4.7/etc/bashrc

set -e

nsatpoints=$(expr $TotalTime/$WriteTimestep | bc)

cp constant/transportPropertiesTP constant/transportProperties
cp system/fvSolutionTP system/fvSolution
cp system/fvSchemesTP system/fvSchemes
cp system/controlDictTP system/controlDict
sed -i "s/Visc1/$Visc1/g" constant/transportProperties
sed -i "s/Visc2/$Visc2/g" constant/transportProperties
sed -i "s/rho1/$rho1/g" constant/transportProperties
sed -i "s/rho2/$rho2/g" constant/transportProperties
sed -i "s/ift/$ift/g" constant/transportProperties
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict
sed -i "s/runTimestep/$runTimestep/g" system/controlDict
sed -i "s/nsatpoints/$nsatpoints/g" system/postProcessDict

rm -rf 0 [1-9]*

cp -r 0_org 0

sed -i "s/flowrate/$flowrate/g" 0/U

echo -e "set alpha field"
setFields > setFieldsTP.out
echo -e "Decompose parallel mesh"
decomposePar > decomposeParTP.out
echo -e "Run interGCFoam"
mpiexec -np 24 interGCFoam  -parallel > interGCFoamTP.out
echo -e "reconstruct parallel mesh"
reconstructPar > reconstructParMeshTP.out
rm -rf proc*

echo -e "process relative permeability"
processRelPerm > processRelPermTP.out
