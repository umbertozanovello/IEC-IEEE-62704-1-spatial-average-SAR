# IEC/IEEE 62704-1: spatial-average SAR
Python-C implementation of the spatial-average SAR alogrithm compliant to the requirements of the IEC/IEEE 62704-1 standard. The core parts of the algorithm are implemented in C in order to improve the code efficiency. Input and output data are managed in Python.

## Requirements
**Python**
Python (>=3.7)
numpy
ctypes

**For C code compiation (optional)**
[cmake](https://cmake.org/) (>= 3.24)
[LAPACK](https://www.netlib.org/lapack/) (optional)

## The IEC/IEEE 62704-1 spatial-average SAR procedure
The patial-average SAR procedure described in the standard is based on two steps. The procedure is reported in Figure 2 of the standard and briefly described below. 
Before starting, all the voxels in the computation domain are flagged as UNUSED.

### Step 1
For each body voxel, an averaging cube is evenly expanded around it until the target averaging mass is reached. If the obtained cube contain less than 10 % of background and all it faces touch or intersect tissues, the voxel is flagged as VALID and is assigned the SAR averaged over the identified cube. All the voxels entirely included in the averaging volume are temporarily flagged as USED. At the end of the first iteration through all the voxels, the voxels still flagged as USED are assigned the largest spatial-average SAR of the averaging cubes in which they were enclosed.

### Step 2
For each UNUSED body voxel, six averaging cubes are constructed having the voxel centered on one of their faces. The other five faces are evenly extended until the target mass is reached, irrespective of the volume of background included in the averaging cube. The volume of the smallest of these six cubes is computed and the SAR is averaged on the cubes whose volume is not more than 5 % larger than the smallest cube. The maximum among the computed SAR values is assigned to the UNUSED voxel.

## Implementation pipeline
The implementation uses Python for pre- and post-process the data. The core part of the algorithm is instead implemented in C and interfaced with Python through the *ctypes* library. The following figure carifies the pipeline.




## Usage
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

## C libraries compilation
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
???

## Acknowledgement
This script has been developed in the context of the 21NRM05 STASIS project. The project (21NRM05 STASIS) has received funding from the European Partnership on Metrology, co-financed from the European Union’s Horizon Europe Research and Innovation Programme and by the Participating States.