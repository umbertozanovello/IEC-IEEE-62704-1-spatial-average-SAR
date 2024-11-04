#include "globals.h"

#ifndef AUXFUNCS_H
#define AUXFUNCS_H

int from3DIdxTo1DIdx(int* idxs, int* n_points);
void setValue_f(double* array, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double value);
void setValue_i(int* array, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, int value);
void setValue_ifPositive_ff(double* array, double* condArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double value);
void setWithOther_ifPositive_ff(double* array, double* otherArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points);
double nansum_f(double* array, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points);
double nansum_ifPositive_ff(double* array, double* testArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points);
double nansum_ifPositive_fi(double* array, int* testArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points);
double nansum_product_ff(double* array1, double* array2, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points);
void incrementByOther(double* array1, double* array2, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double multFactor);
double* copyArray_f(double* array, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points);
int checkNanFaces(double* array, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points);
double backgroundRatio(double* massArray, double* massUsageMask, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points);
double targetMassAverage(double* massArray, double* massUsageMask, double* localSARArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double targetMass);
void markAndSetUSEDvox(double* avgSARArray, Status* voxStatusArray, double* massUsageMask, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double avgSARValueVALID);
double solveForK(double a, double b, double c, double d);
void getAnnularMasses_CC(double* massArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double* masses);
void setPartialMassUsageMask_CC(double* massUsageMask, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double k);
void getAnnularMasses_SC(double* massArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double* masses, int delta, int face);
void setPartialMassUsageMask_SC(double* massUsageMask, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double k, int delta, int face);
double nanmax_f(double* array, int length);

#endif
