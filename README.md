# IEC/IEEE 62704-1: spatial-average SAR
Python-C implementation of the spatial-average SAR alogrithm compliant to the requirements of the IEC/IEEE 62704-1 standard. The core parts of the algorithm are implemented in C in order to improve the code efficiency. Input and output data are managed in Python.

## Requirements
**Python**
Python (>=3.7)
numpy
ctypes

**For C code compilation (optional)**
[cmake](https://cmake.org/) (>= 3.24)
[LAPACK](https://www.netlib.org/lapack/) (optional)

## The IEC/IEEE 62704-1 spatial-average SAR procedure
The patial-average SAR procedure described in the standard is based on two steps. The procedure is reported in Figure 2 of the standard and briefly described below. 
Before starting, all the voxels in the computation domain are flagged as UNUSED.

### Step 1
For each body voxel, an averaging cube is evenly expanded around it until the target averaging mass is reached. If the obtained cube contain less than 10 % of background and all it faces touch or intersect tissues, the voxel is flagged as VALID and is assigned the SAR averaged over the identified cube. All the voxels entirely included in the averaging volume are temporarily flagged as USED. At the end of the first iteration through all the voxels, the voxels still flagged as USED are assigned the largest spatial-average SAR of the averaging cubes in which they were enclosed.

### Step 2
For each UNUSED body voxel, six averaging cubes are constructed having the voxel centered on one of their faces. The other five faces are evenly extended until the target mass is reached, irrespective of the volume of background included in the averaging cube. The volume of the smallest of these six cubes is computed and the SAR is averaged on the cubes whose volume is not more than 5 % larger than the smallest cube. The maximum among the computed SAR values is assigned to the UNUSED voxel.

The following figure shows the different averaging cubes and voxel flags during the two steps of the procedure

![alt text](https://github.com/umbertozanovello/IEC-IEEE-62704-1-spatial-average-SAR/blob/main/images/AveragingCubes.jpg?raw=true)

## The spatial-average SAR script

### Implementation pipeline
The implementation uses Python for pre- and post-process the data. The core part of the algorithm is instead implemented in C and interfaced with Python through the *ctypes* library. The following figure clarifies the pipeline.

![alt text](https://github.com/umbertozanovello/IEC-IEEE-62704-1-spatial-average-SAR/blob/main/images/pipeline.jpg?raw=true)

### Usage
The spatial-average SAR can be computed calling the Python function `spatialAverageSAR` defined in the spatialAverageSAR.py file. The function returns two numpy ndarrays, one containing the spatial-average SAR values (`avgSARArray`) and the other containing the flags assigned to the voxel (`voxStatusArray`). The flags are associated to integers according to the following definition:
- 0 : INVALID (background)
- 1 : UNUSED
- 2 : USED
- 3 : VALID

Both `avgSARArray` and `voxStatusArray` dimensions are *nx* x *ny* x *nz*, that is size of the domain along the *x*-, *y*- and *z*-Cartesian direction.

The `spatialAverageSAR` Python function takes the following arguments:
- `massArray` (numpy ndarray): *nx* x *ny* x *nz* points array containing the voxel masses in kg. The background is identified by `nan` values
- `localSARArray` (numpy ndarray): *nx* x *ny* x *nz* points array containing the local SAR distribution in W/kg
- `targetMass` (float): target mass in kg over which averaging the local SAR
- `step1_libPath` (string): path to the "Step 1" shared library
- `step2_libPath` (string): path to the "Step 2" shared library
- `additionalBackground` (list): three element list indicating the additional background voxels to be added for computation along the x,y,z directions

The C shared libraries can be either compiled from the source C code (see section below) or downloaded from the release assets. precompiled libraries are indeed provided both for Linux (.so extension) and Windows (.dll extension) OS. The algorithm requires to repeatedly compute the roots of a third degree polynomial. In the provided C code this is achieved either with a faster but approximated algorithm or by a slower but more accurate algorithm based on the computation of the eigenvalues of the companion matrix. Libraries named as "libAvgSARStep*_cmroots.*" compute the polynomial roots with the latter strategy.

### C libraries compilation
Whereas it is possible to use the precompiled libraries provieded along with the current release, it is also possible to compile them from the source code. This has the advantage that the libraries are optimised for the system on which they are compiled. Compilation can be easily performed with [cmake](https://cmake.org/) and the provided "CMakeLists.txt" file provided. Linux user can run the following commands in a terminal opened in the "CMakeLists.txt" directory:
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```
According to the value assigned to the `CMROOTS_ENABLE` option in the "CMakeLists.txt" file, the libraries rely on the companion matrix algorithm (`CMROOTS_ENABLE=ON`) or approximated algorithm (`CMROOTS_ENABLE=OFF`) to compute the third degree polynomial roots. 
If the cmake option `REPORT_ENABLE` is set to `ON`, a report is generated with the name and folder specified in the `REPORT_PATH` variable. To be compliant with the report format of the IEC/IEEE 62704-1 standard, the generated report should be reorganized with the Python script "sortReport.py" available in the "tools" folder.

## Testing
The IEC/IEEE 62704-1 standard procides the so-called “SAR Star” for testing different algorithm implementations. The “SAR Star” is a bi-compartmental solid designed to involve all the averaging volume types for all directions of space during the spatial-averaging algorithm execution.

![alt text](https://github.com/umbertozanovello/IEC-IEEE-62704-1-spatial-average-SAR/blob/main/images/SARStarExploded.png?raw=true)

The standard provides an *stl* model of the “SAR Star” together with local SAR values and reference results in terms of:
- spatial-average SAR
- averaging cube masses
- averaging cube volumes
- voxel flags
- averaging cube expansion directions (relevant only to Step 2)

These data are collected in text files and can be downloaded [here](https://assets.iec.ch/public/tc106/62704-1_supplemental_files.zip?2024111817)

To improve the usability of the test data, binary *mat* files are provided as assets with the release. Two files are provided where the “SAR Star” is discretized on a uniform grid and reference results are given for 1 g and 10 g averaging masses. Each file contains the following attributes:
- *x_offset*, *y_offset*, *z_offset*: *x*, *y*, and *z* coordinates;
- *local_SAR*: local SAR values in W/kg to be averaged;
- *star*: material code of each voxel from 0 (background) to 2;
- *densities*: Mass density in kg/m^3 of the materials. First element is background, second element is material 1 and last element is material 2;
- *average_SAR*: spatial-average SAR values;
- *mass*: masses of the averaging cubes in g;
- *volume*: volumes of the averaging cubes in mm^3;
- *status*: flags associated to the voxels according to the following definition:
    - 0 : INVALID (background)
    - 1 : UNUSED
    - 2 : USED
    - 3 : VALID
- *orientation*: orientation of the averaging cubes with respect to the reference location:
    - face-centred cube (Step 2): -x=1, +x=2, -y=3, +y=4, -z=5, +z=6
    - volume centred cube (Step 1): =7

Finally, the IEC/IEEE 62704-1 standard reports a code (see *sarstar_evaluation_script_V1.1.m* provided along with the [text reference data](https://assets.iec.ch/public/tc106/62704-1_supplemental_files.zip?2024111817)) for generating a log file comparing the reference data against those obtained with the developed algorithm. Log files relevant to the provided algorithm are uploaded [here](https://github.com/umbertozanovello/IEC-IEEE-62704-1-spatial-average-SAR/tree/main/others)
These results are obtained using the companion matrix algorithm to compute the polynomial roots (see above).

## Acknowledgement
This script has been developed in the context of the 21NRM05 STASIS project. The project (21NRM05 STASIS) has received funding from the European Partnership on Metrology, co-financed from the European Union's Horizon Europe Research and Innovation Programme and by the Participating States.