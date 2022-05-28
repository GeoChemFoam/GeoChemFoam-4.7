#!/bin/bash

###### USERS INPUT ############################################################

#Define image name
Image_name="calcitePost"

# Define image dimensions
x_dim=536
y_dim=300
z_dim=40

# Define cropping parameters
# NOTE: The image MUST be cropped by at least 1 voxel on every face
x_min=0
x_max=536
y_min=0
y_max=300
z_min=0
z_max=40

# number of cells of initial mesh
n_x=134
n_y=75
n_z=10

#Mesh refinement level
nlevel=1
nRef=200
lowRef=0.01
upRef=0.99


# flow direction 0 or 1, 2 is empty
direction=0

# Number of processors
NP=24

# define resolution (m)
res=0.000005

#### END OF USER INPUT #######################################################

# Load user environment variables 
source ~/.bashrc

source $HOME/works/GeoChemFoam-4.7/etc/bashrc

set -e

#Insert dimensions in postProcessDict
x_1=$(expr $x_min*$res | bc)
y_1=$(expr $y_min*$res | bc)
z_1=$(expr $z_min*$res | bc)

x_2=$(expr $x_max*$res | bc)
y_2=$(expr $y_max*$res | bc)
z_2=$(expr $z_max*$res | bc)

cp system/postProcessDict1 system/postProcessDict
sed -i "s/x_1/$x_min/g" system/postProcessDict
sed -i "s/y_1/$y_min/g" system/postProcessDict
sed -i "s/z_1/$z_min/g" system/postProcessDict

sed -i "s/x_2/$x_max/g" system/postProcessDict
sed -i "s/y_2/$y_max/g" system/postProcessDict
sed -i "s/z_2/$z_max/g" system/postProcessDict

sed -i "s/flowdir/$direction/g" system/postProcessDict

cp constant/dynamicMeshDictInit constant/dynamicMeshDict
cp constant/dynamicMeshDictInit constant/dynamicMeshDict0
cp constant/dynamicMeshDictInit constant/dynamicMeshDict1

cp constant/dynamicMeshDictAMR constant/dynamicMeshDictAMRInit
cp constant/dynamicMeshDictAMR constant/dynamicMeshDictAMR0
cp constant/dynamicMeshDictAMR constant/dynamicMeshDictAMR1

cp system/fvSolutionAMRInit system/fvSolution

sed -i "s/nRef/1/g" constant/dynamicMeshDictAMRInit
sed -i "s/lowRef/0.01/g" constant/dynamicMeshDictAMRInit
sed -i "s/upRef/1/g" constant/dynamicMeshDictAMRInit
sed -i "s/refLevel/$nlevel/g" constant/dynamicMeshDictAMRInit

sed -i "s/nRef/1/g" constant/dynamicMeshDictAMR0
sed -i "s/lowRef/$lowRef/g" constant/dynamicMeshDictAMR0
sed -i "s/upRef/$upRef/g" constant/dynamicMeshDictAMR0
sed -i "s/refLevel/$nlevel/g" constant/dynamicMeshDictAMR0

sed -i "s/nRef/$nRef/g" constant/dynamicMeshDictAMR1
sed -i "s/lowRef/$lowRef/g" constant/dynamicMeshDictAMR1
sed -i "s/upRef/$upRef/g" constant/dynamicMeshDictAMR1
sed -i "s/refLevel/$nlevel/g" constant/dynamicMeshDictAMR1


cp system/decomposeParDict0 system/decomposeParDict
sed -i "s/NP/$NP/g" system/decomposeParDict

# Create background mesh
echo -e "Create background mesh"
cp system/blockMeshDict$direction constant/polyMesh/blockMeshDict

sed -i "s/x_min/$x_1/g" constant/polyMesh/blockMeshDict
sed -i "s/y_min/$y_1/g" constant/polyMesh/blockMeshDict
sed -i "s/z_min/$z_1/g" constant/polyMesh/blockMeshDict

sed -i "s/x_max/$x_2/g" constant/polyMesh/blockMeshDict
sed -i "s/y_max/$y_2/g" constant/polyMesh/blockMeshDict
sed -i "s/z_max/$z_2/g" constant/polyMesh/blockMeshDict

sed -i "s/nx/$n_x/g" constant/polyMesh/blockMeshDict
sed -i "s/ny/$n_y/g" constant/polyMesh/blockMeshDict
sed -i "s/nz/$n_z/g" constant/polyMesh/blockMeshDict

cp system/controlDictInit system/controlDict
blockMesh > blockMesh.out


rm -rf 0
cp -r 0_org 0

cd constant/triSurface
tar -zxvf $Image_name.raw.tgz 
cd ../..

#Dummy fluid properties
flowRate=1e-30
Visc=1e-6
Diff=1e-30
kreac=0
scoeff=1
rhos=2710
Mws=100
cinlet=0
kf=0 

cp constant/transportProperties0 constant/transportProperties
sed -i "s/Visc/$Visc/g" constant/transportProperties
sed -i "s/rho_s/$rhos/g" constant/transportProperties
sed -i "s/Mw_s/$Mws/g" constant/transportProperties
sed -i "s/k_f/$kf/g" constant/transportProperties

cp constant/thermoPhysicalProperties0 constant/thermoPhysicalProperties
sed -i "s/Diff/$Diff/g" constant/thermoPhysicalProperties
sed -i "s/s_coeff/$scoeff/g" constant/thermoPhysicalProperties
sed -i "s/k_reac/$kreac/g" constant/thermoPhysicalProperties

sed -i "s/flow_rate/$flowRate/g" 0/U
sed -i "s/c_inlet/$cinlet/g" 0/C

echo -e "calculate eps for initial image"
python << END
import numpy as np
import h5py
import array
import os

f = open('constant/triSurface/$Image_name.raw','rb')
img = np.fromfile(f, dtype=np.uint8)
img = np.reshape(img,($z_dim,$y_dim,$x_dim))
my_array=img


ncell = $n_x*$n_y*$n_z 
p=($x_max-$x_min)/$n_x
q=($y_max-$y_min)/$n_y
r=($z_max-$z_min)/$n_z
f=open('0/eps','a')
f.seek(0) #get to the first position
f.write("FoamFile"+'\n')
f.write("{"+'\n')
f.write("    version     2.0;"+'\n')
f.write("    format      ascii;"+'\n')
f.write("    class       volScalarField;"+'\n')
f.write("    object      eps;"+'\n')
f.write("}"+'\n')
f.write(""+'\n')
f.write("dimensions      [0 0 0 0 0 0 0];"+'\n')
f.write("internalField   nonuniform List<scalar>"+'\n')
f.write(str(ncell)+'\n')
f.write("("+'\n')
eps=np.zeros(($n_x,$n_y,$n_z),dtype=float)
for i in range (0,$n_x):
    for ii in range (0,p):
        for j in range (0,$n_y):
            for jj in range (0,q):
                for k in range (0,$n_z):
                    for kk in range (0,r):
                        eps[i,j,k]+=(my_array[k*r+kk,j*q+jj,i*p+ii]/255.0+ (1-my_array[k*r+kk,j*q+jj,i*p+ii]/255.0)*0.001)/p/q/r
for k in range (0,$n_z):
    for j in range (0, $n_y):
        for i in range (0, $n_x):
            f.write(str(eps[i,j,k])+'\n')
f.write(")"+'\n')
f.write(";"+'\n')
f.write(""+'\n')
f.write("boundaryField"+'\n')
f.write("{"+'\n')
f.write("    walls"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    inlet"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    outlet"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("}"+'\n')
f.close()

os.system('echo "decompose parallel mesh"')
os.system('decomposePar > decomposePar.out')

if ($nlevel>0):
    os.system('cp constant/dynamicMeshDictAMRInit constant/dynamicMeshDict')
    os.system('cp constant/dynamicMeshDictAMR0 constant/dynamicMeshDict0')
    os.system('cp constant/dynamicMeshDictAMR1 constant/dynamicMeshDict1')

    os.system('sed -i "s/nRef/1/g" constant/dynamicMeshDict')
    os.system('sed -i "s/lowRef/0.01/g" constant/dynamicMeshDict')
    os.system('sed -i "s/upRef/1/g" constant/dynamicMeshDict')
    os.system('sed -i "s/refLevel/$nlevel/g" constant/dynamicMeshDict')

    os.system('sed -i "s/nRef/1/g" constant/dynamicMeshDictAMR0')
    os.system('sed -i "s/lowRef/$lowRef/g" constant/dynamicMeshDictAMR0')
    os.system('sed -i "s/upRef/$upRef/g" constant/dynamicMeshDictAMR0')
    os.system('sed -i "s/refLevel/$nlevel/g" constant/dynamicMeshDictAMR0')

    os.system('sed -i "s/nRef/$nref/g" constant/dynamicMeshDictAMR1')
    os.system('sed -i "s/lowRef/$lowRef/g" constant/dynamicMeshDictAMR1')
    os.system('sed -i "s/upRef/$upRef/g" constant/dynamicMeshDictAMR1')
    os.system('sed -i "s/refLevel/$nlevel/g" constant/dynamicMeshDictAMR1')

    dt=1e-6/$nlevel

    os.system('cp system/controlDictAMRInit system/controlDict')
    os.system('sed -i "s/dt/'+str(dt)+'/g" system/controlDict')

    os.system('echo "refine mesh at interface by running reactiveTransportDBSFoam for small time"')
    os.system("mpiexec -np $NP reactiveTransportDBSFoam -parallel > reactiveTransportDBSFoamInit.out")

    os.system('echo "reconstruct parallel mesh"')
    os.system("reconstructParMesh -latestTime > reconstructParMeshInit.out")

    os.system("rm -rf processor*/0")
    for i in range(0,$NP):
        os.system('mv processor'+str(i)+'/1e-06 processor'+str(i)+'/0')

    os.system("rm processor*/0/phi")
    os.system("rm processor*/0/ddt*")
    os.system("rm processor*/0/R")
    os.system("rm -rf processor*/0/uniform")

    os.system("rm -rf 0")
    os.system("mv 1e-06 0")

    os.system("rm 0/phi")
    os.system("rm 0/ddt*")
    os.system("rm 0/R")
    os.system("rm -rf 0/uniform")
    os.system("cp 0_org/U 0/.")
    os.system("cp 0_org/p 0/.")
    os.system("cp 0_org/C 0/.")
    os.system('sed -i "s/flow_rate/$flowRate/g" 0/U')
    os.system('sed -i "s/c_inlet/$cinlet/g" 0/C')

    os.system('echo "process mesh centers"')
    os.system("processMeshCellCenters > processMeshCellCenters.out")

    nz=2*$n_z*$nlevel
    ny=2*$n_y*$nlevel
    nx=2*$n_x*$nlevel

    p=($x_max-$x_min)/nx
    q=($y_max-$y_min)/ny
    r=($z_max-$z_min)/nz

    os.system('echo "calculate eps for fine mesh"')
    eps=np.zeros((nx,ny,nz),dtype=float)
    for i in range (0,nx):
        for ii in range (0,p):
            for j in range (0,ny):
                for jj in range (0,q):
                    for k in range (0,nz):
                        for kk in range (0,r):
                            eps[i,j,k]+=(my_array[k*r+kk,j*q+jj,i*p+ii]/255.0+ (1-my_array[k*r+kk,j*q+jj,i*p+ii]/255.0)*0.001)/p/q/r

    dx = ($x_2-$x_1)/nx
    dy = ($y_2-$y_1)/ny
    dz = ($z_2-$z_1)/nz
    
    ix = np.zeros(nx*ny*nz)
    iy = np.zeros(nx*ny*nz)
    iz = np.zeros(nx*ny*nz)

    file = open("0/cellCenters","r")
    Lines = file.readlines()
    count =0
    wbool=0
    for line in Lines:
      ls = line.strip()
      if (ls==")"):
          break
      if (wbool==1):
          x=float(ls.split("(")[1].split(")")[0].split()[0])
          ix[count] = np.floor(x/dx);
          y=float(ls.split("(")[1].split(")")[0].split()[1])
          iy[count] = np.floor(y/dy);
          z=float(ls.split("(")[1].split(")")[0].split()[2])
          iz[count] = np.floor(z/dz);
          count +=1
      if (ls=="("):
          wbool=1

    ncell = count

    newEps = np.zeros(ncell)

    os.system('echo "modify eps according to finer mesh for cells at interface using cell center"')

    file = open("0/eps","r")
    Lines = file.readlines()
    count =0
    wbool=0
    for line in Lines:
      ls = line.strip()
      if (ls==")"):
          break
      if (wbool==1):
          epsVal=float(ls)
          if (epsVal<1) and (epsVal>0.01):
            newEps[count] = eps[ix[count].astype(int),iy[count].astype(int),iz[count].astype(int)]
          else:
            newEps[count] = epsVal
          count +=1
      if (ls=="("):
          wbool=1

    os.system('rm 0/eps')
    f=open('0/eps','a')
    f.seek(0) #get to the first position
    f.write("FoamFile"+'\n')
    f.write("{"+'\n')
    f.write("    version     2.0;"+'\n')
    f.write("    format      ascii;"+'\n')
    f.write("    class       volScalarField;"+'\n')
    f.write("    object      eps;"+'\n')
    f.write("}"+'\n')
    f.write(""+'\n')
    f.write("dimensions      [0 0 0 0 0 0 0];"+'\n')
    f.write("internalField   nonuniform List<scalar>"+'\n')
    f.write(str(ncell)+'\n')
    f.write("("+'\n')

    for k in range(0,ncell):
            f.write(str(newEps[k])+'\n')
    f.write(")"+'\n')
    f.write(";"+'\n')
    f.write(""+'\n')
    f.write("boundaryField"+'\n')
    f.write("{"+'\n')
    f.write("    walls"+'\n')
    f.write("    {"+'\n')
    f.write("        type zeroGradient;"+'\n')
    f.write("    }"+'\n')
    f.write("    inlet"+'\n')
    f.write("    {"+'\n')
    f.write("        type zeroGradient;"+'\n')
    f.write("    }"+'\n')
    f.write("    outlet"+'\n')
    f.write("    {"+'\n')
    f.write("        type zeroGradient;"+'\n')
    f.write("    }"+'\n')
    f.write("}"+'\n')
    f.close()

    os.system('rm -f 0/cellCenters')


END

rm constant/triSurface/$Image_name.raw
echo -e "Case initialised"





