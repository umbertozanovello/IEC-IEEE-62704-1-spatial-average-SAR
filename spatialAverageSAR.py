import numpy as np
from ctypes import *

def spatialAverageSAR(massArray, localSARArray, targetMass, step1_libPath, step2_libPath, additionalBackground=[0,0,0]):
    """It computes the spatial average SAR averaged over the targetMass

    Args:
        massArray (numpy ndarray): nx x ny x nz points array containing the voxel masses in kg. The background is identified by nan values
        localSARArray (numpy ndarray): nx x ny x nz points array containing the local SAR distribution in W/kg
        targetMass (float): target mass in kg over which averaging the local SAR
        step1_libPath (string): path to the "Step 1" shared library
        step2_libPath (string): path to the "Step 2" shared library
        additionalBackground (list): 
    """

    # Loading the libraries
    step1_libPath = CDLL(step1_libPath)
    step2_libPath = CDLL(step2_libPath)

    avgSARStep1 = step1_libPath.main
    avgSARStep1.restype = c_int
    avgSARStep2 = step2_libPath.main
    avgSARStep2.restype = c_int

    # Preparing the numpy data
    original_n_points = massArray.shape
    additionalBackground = np.array(additionalBackground)
    n_points = original_n_points + additionalBackground*2

    voxStatusArray = np.ones_like(massArray, dtype=int) # INVALID=0, UNUSED=1, USED=2, VALID=3

    original_slices = []
    for i in range(3):
        original_slices.append(slice(additionalBackground[i],additionalBackground[i]+original_n_points[i]))

    if ((additionalBackground)>0).any():
        
        
        box = np.full(n_points, np.nan)
        box[*original_slices] = massArray
        massArray = np.copy(box)
        
        box[*original_slices] = localSARArray
        localSARArray = np.copy(box)
        localSARArray[np.isnan(localSARArray)] = 0

        box[*original_slices] = voxStatusArray
        voxStatusArray = np.copy(box)
    
    voxStatusArray[np.isnan(massArray)] = 0
    voxStatusArray = voxStatusArray.astype(int)

    # ctypes definitions

    voxStatusArray_c = (c_int * voxStatusArray.size)(*voxStatusArray.flatten())
    massArray_c = (c_double * massArray.size)(*massArray.flatten())
    avgSARArray_c = (c_double * massArray.size)(*(np.zeros(massArray.size)))
    localSARArray_c = (c_double * localSARArray.size)(*localSARArray.flatten())
    n_points_c = (c_int * 3)(*n_points)
    targetMass_c = c_double(targetMass)

    # C functions execution
    
    print("Step one in progress ...\n\n")
    ret = avgSARStep1(massArray_c, localSARArray_c, avgSARArray_c, voxStatusArray_c, targetMass_c, n_points_c)

    print("Step two in progress ...\n\n")
    ret = avgSARStep2(massArray_c, localSARArray_c, avgSARArray_c, voxStatusArray_c, targetMass_c, n_points_c)

    avgSARArray = np.array(avgSARArray_c).reshape(n_points)[*original_slices]
    voxStatusArray = np.array(voxStatusArray_c).reshape(n_points)[*original_slices].astype(float)
    localSARArray = localSARArray[*original_slices]

    return avgSARArray, voxStatusArray