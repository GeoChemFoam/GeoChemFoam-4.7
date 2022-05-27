#!/bin/bash

###### USERS INPUT ############################################################

#Define image name
Image_name="HM12_2_4"

# Define image dimensions
x_dim=2000
y_dim=2000
z_dim=10

# Define cropping parameters
# NOTE: The image MUST be cropped by at least 1 voxel on every face
x_min=501
x_max=1501
y_min=701
y_max=1301
z_min=1
z_max=6

# Percent of image cut for permeability calculation (boundary and capillary end effect)
cut_x=0.25
cut_y=0.05

# number of cells of initial mesh
# n*(nlevel+1) should be equal to image dimension when not binning 
n_x=500
n_y=300
#In 2D images, z has to be the empty direction
n_z=1

# flow direction 0 or 1, 2 is empty
direction=0

# define resolution (m)
res=0.00003

# define value of the pores and solid in the image 
pores_value=255
solid_value=0

# Number of processors
NP=24

#### END OF USER INPUT #######################################################

# Load user environment variables 
source ~/.bashrc

source $HOME/works/GeoChemFoam-dev/etc/bashrc

source $OF4X_DIR/OpenFOAM-4.x/etc/bashrc WM_LABEL_SIZE=64 FOAMY_HEX_MESH=yes

MPIRUN=mpirun 
echo -e "make stl"
cd constant/triSurface
tar -xf $Image_name\.raw.tar.gz
python raw2stl.py --x_min=$x_min --x_max=$x_max --y_min=$y_min --y_max=$y_max --z_min=$z_min --z_max=$z_max --pores_value=$pores_value --solid_value=$solid_value  --image_name=$Image_name --x_dim=$x_dim --y_dim=$y_dim --z_dim=$z_dim
rm $Image_name\.raw
cd ../..
# ./runSmooth.sh

export pore_index_0="$(cat constant/triSurface/pore_indx)"
export pore_index_1="$(cat constant/triSurface/pore_indy)"
export pore_index_2="$(cat constant/triSurface/pore_indz)"
 
echo -e "Coordinates at center of a pore = ($pore_index_0,$pore_index_1,$pore_index_2)" 
 
rm constant/triSurface/pore_ind*

# Create background mesh
echo -e "Create background mesh"
cp system/blockMeshDict2D$direction system/blockMeshDict
dx=1
let dx+=$x_max
let dx-=$x_min
dy=1
let dy+=$y_max
let dy-=$y_min
dz=1
let dz+=$z_max
let dz-=$z_min

sed -i "s/dx/$dx/g" system/blockMeshDict
sed -i "s/dy/$dy/g" system/blockMeshDict
sed -i "s/dz/$dz/g" system/blockMeshDict

sed -i "s/nx/$n_x/g" system/blockMeshDict
sed -i "s/ny/$n_y/g" system/blockMeshDict
sed -i "s/nz/$n_z/g" system/blockMeshDict

cp system/controlDictInit system/controlDict

blockMesh  > blockMesh.out

cp system/decomposeParDictRun system/decomposeParDict
sed -i "s/NP/$NP/g" system/decomposeParDict


# Decompose background mesh
echo -e "Decompose background mesh"
decomposePar > decomposeBlockMesh.out
rm -rf processor*/0/*

# Run snappyHexMesh in parallel
echo -e "Run snappyHexMesh in parallel"
cp system/snappyHexMeshDict1 system/snappyHexMeshDict

sed -i "s/poreIndex0/$pore_index_0/g" system/snappyHexMeshDict
sed -i "s/poreIndex1/$pore_index_1/g" system/snappyHexMeshDict
sed -i "s/poreIndex2/$pore_index_2/g" system/snappyHexMeshDict


x_1=$(expr $cut_x*$dx*$res | bc)
y_1=$(expr $cut_y*$dy*$res | bc)
z_1=$res

x_2=$(expr $dx*$res-$x_1 | bc)
y_2=$(expr $dy*$res-$y_1 | bc)
z_2=$(expr $dz*$res | bc)


cp system/postProcessDict1 system/postProcessDict
sed -i "s/x_1/$x_1/g" system/postProcessDict
sed -i "s/y_1/$y_1/g" system/postProcessDict
sed -i "s/z_1/$z_1/g" system/postProcessDict

sed -i "s/x_2/$x_2/g" system/postProcessDict
sed -i "s/y_2/$y_2/g" system/postProcessDict
sed -i "s/z_2/$z_2/g" system/postProcessDict

sed -i "s/flowdir/$direction/g" system/postProcessDict


$MPIRUN -np $NP snappyHexMesh -overwrite -parallel  > snappyHexMesh.out

echo -e "reconstruct parallel mesh"
reconstructParMesh -constant > reconstructParMesh.out

echo -e "transformPoints" 
vector="($res $res $res)"
transformPoints -scale "$vector" > transformPoints.out

rm -rf processor*

echo -e "Image Initialised. It is advised to check in paraview to confirm mesh of porespace is reasonable before running flow" 
