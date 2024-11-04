# include "include/globals.h"

int from3DIdxTo1DIdx(int* idxs, int* n_points)
{ /* It returns the index of the flatten array with c-order corresponding to the indices contained in the idxs array relevant to a 3D array with shape equal to n_points*/
    int delta_j = n_points[1];
    int delta_k = n_points[2];

    return idxs[0]*delta_j*delta_k + idxs[1]*delta_k + idxs[2];
}

void setValue_f(double* array, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double value)
{ /* array[i_min:i_max, j_min:j_max, k_min:k_max] = value*/   

    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = 0;
    int j_range = 0;
    int k_range = 0;

    i_range = (i_max>=0) ? (i_max-i_min) : 1;
    j_range = (j_max>=0) ? (j_max-j_min) : 1;
    k_range = (k_max>=0) ? (k_max-k_min) : 1;
    
    for (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                array[idx] = value;
            }
        }
    }
    return;
}

void setValue_i(int* array, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, int value)
{ /* array[i_min:i_max, j_min:j_max, k_min:k_max] = value*/   

    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = 0;
    int j_range = 0;
    int k_range = 0;

    i_range = (i_max>=0) ? (i_max-i_min) : 1;
    j_range = (j_max>=0) ? (j_max-j_min) : 1;
    k_range = (k_max>=0) ? (k_max-k_min) : 1;
    
    for (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                array[idx] = value;
            }
        }
    }
    return;
}

void setValue_ifPositive_ff(double* array, double* condArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double value)
{ /* array[i_min:i_max, j_min:j_max, k_min:k_max] = value*/   

    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = 0;
    int j_range = 0;
    int k_range = 0;

    i_range = (i_max>=0) ? (i_max-i_min) : 1;
    j_range = (j_max>=0) ? (j_max-j_min) : 1;
    k_range = (k_max>=0) ? (k_max-k_min) : 1;
    
    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                if(condArray[idx]>0)
                {
                    array[idx] = value;
                }
            }
        }
    }
    return;
}

void setWithOther_ifPositive_ff(double* array, double* otherArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points)
{ /* array[i_min:i_max, j_min:j_max, k_min:k_max] = value*/   

    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = 0;
    int j_range = 0;
    int k_range = 0;

    i_range = (i_max>=0) ? (i_max-i_min) : 1;
    j_range = (j_max>=0) ? (j_max-j_min) : 1;
    k_range = (k_max>=0) ? (k_max-k_min) : 1;
    
    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                if(otherArray[idx]>0)
                {
                    array[idx] = otherArray[idx];
                }
            }
        }
    }
    return;
}

double nansum_f(double* array, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points)
{/*nansum of array within the indices*/

    double sum = 0;
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = 0;
    int j_range = 0;
    int k_range = 0;

    i_range = (i_max>=0) ? (i_max-i_min) : 1;
    j_range = (j_max>=0) ? (j_max-j_min) : 1;
    k_range = (k_max>=0) ? (k_max-k_min) : 1;
    
    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                if (array[idx] == array[idx]) // if array[idx] is not nan...
                {
                    sum += array[idx];
                }
            }
        }
    }
    return sum;
}

double nansum_ifPositive_ff(double* array, double* testArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points)
{ /*nansum within the pecified indices if the testArray value is higher than zero*/
    double sum = 0;
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = 0;
    int j_range = 0;
    int k_range = 0;

    i_range = (i_max>=0) ? (i_max-i_min) : 1;
    j_range = (j_max>=0) ? (j_max-j_min) : 1;
    k_range = (k_max>=0) ? (k_max-k_min) : 1;
    
    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                if (array[idx] == array[idx] &&  testArray[idx] > 0) // if array[idx] is not nan...
                {
                    sum += array[idx];
                }
            }
        }
    }
    return sum;
}

double nansum_ifPositive_fi(double* array, int* testArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points)
{ /*nansum within the pecified indices if the testArray value is higher than zero*/
    double sum = 0;
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = 0;
    int j_range = 0;
    int k_range = 0;

    i_range = (i_max>=0) ? (i_max-i_min) : 1;
    j_range = (j_max>=0) ? (j_max-j_min) : 1;
    k_range = (k_max>=0) ? (k_max-k_min) : 1;
    
    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                if (array[idx] == array[idx] &&  testArray[idx] > 0) // if array[idx] is not nan...
                {
                    sum += array[idx];
                }
            }
        }
    }
    return sum;
}

double nansum_product_ff(double* array1, double* array2, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points)
{ /*nansum within the pecified indices if the testArray value is higher than zero*/
    double sum = 0;
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = 0;
    int j_range = 0;
    int k_range = 0;

    i_range = (i_max>=0) ? (i_max-i_min) : 1;
    j_range = (j_max>=0) ? (j_max-j_min) : 1;
    k_range = (k_max>=0) ? (k_max-k_min) : 1;
    
    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                if (array1[idx] == array1[idx] && array2[idx] == array2[idx]) // if neither array1[idx] nor array2[idx] is nan
                {
                    sum += array1[idx] * array2[idx];
                }
            }
        }
    }
    return sum;
}

void incrementByOther(double* array1, double* array2, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double multFactor)
{ /* array1[i_min:i_max, j_min:j_max, k_min:k_max] += multFactor*array2[i_min:i_max, j_min:j_max, k_min:k_max]*/

    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = 0;
    int j_range = 0;
    int k_range = 0;

    i_range = (i_max>=0) ? (i_max-i_min) : 1;
    j_range = (j_max>=0) ? (j_max-j_min) : 1;
    k_range = (k_max>=0) ? (k_max-k_min) : 1;
    
    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                array1[idx] += multFactor*array2[idx];
            }
        }
    }
    return;
}

double* copyArray_f(double* array, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points)
{
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx, copIdx=0;
    int length;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = 0;
    int j_range = 0;
    int k_range = 0;

    i_range = (i_max>=0) ? (i_max-i_min) : 1;
    j_range = (j_max>=0) ? (j_max-j_min) : 1;
    k_range = (k_max>=0) ? (k_max-k_min) : 1;

    length = i_range * j_range * k_range;

    double *copiedArray = (double*)malloc(length * sizeof(double));

    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                copiedArray[copIdx++] = array[idx];
                
            }
        }
    }
    return copiedArray;
}

int checkNanFaces(double* array, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points)
{ /*Returns 1 if for all faces of the cube within the n_min and n_max indices of array there's at least one value different from nan*/
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    char flag_im=0, flag_ip=0, flag_jm=0, flag_jp=0, flag_km=0, flag_kp=0;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = 0;
    int j_range = 0;
    int k_range = 0;

    i_range = (i_max-i_min);
    j_range = (j_max-j_min);
    k_range = (k_max-k_min);
    
    // i-
    idx_t1 = (i_min)*delta_j*delta_k;
    for  (int idx_j=0; idx_j<j_range; idx_j++)
    {
        idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

        for  (int idx_k=0; idx_k<k_range; idx_k++)
        {
            idx = idx_t2 + (idx_k+k_min);

            if (array[idx]==array[idx]) // At least one element on the cube face is not nan
            {   
                flag_im = 1;
                break;
            }
        }
    }

    // i+
    idx_t1 = (i_max-1)*delta_j*delta_k;
    for  (int idx_j=0; idx_j<j_range; idx_j++)
    {
        idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

        for  (int idx_k=0; idx_k<k_range; idx_k++)
        {
            idx = idx_t2 + (idx_k+k_min);

            if (array[idx]==array[idx]) // At least one element on the cube face is not nan
            {   
                flag_ip = 1;
                break;
            }
        }
    }

    // j-
    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k + j_min*delta_k;

        for  (int idx_k=0; idx_k<k_range; idx_k++)
        {
            idx = idx_t1 + (idx_k+k_min);
            
            if (array[idx]==array[idx]) // At least one element on the cube face is not nan
            {   
                flag_jm = 1;
                break;
            }
        }
    }

    // j+
    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k + (j_max-1)*delta_k;

        for  (int idx_k=0; idx_k<k_range; idx_k++)
        {
            idx = idx_t1 + (idx_k+k_min);
            
            if (array[idx]==array[idx]) // At least one element on the cube face is not nan
            {   
                flag_jp = 1;
                break;
            }
        }
    }

    // k-
    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx = idx_t1 + (idx_j+j_min)*delta_k + k_min;
            
            if (array[idx]==array[idx]) // At least one element on the cube face is not nan
            {   
                flag_km = 1;
                break;
            }
        }
    }

    // k+
    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx = idx_t1 + (idx_j+j_min)*delta_k + (k_max-1);
            
            if (array[idx]==array[idx]) // At least one element on the cube face is not nan
            {   
                flag_kp = 1;
                break;
            }
        }
    }

    return flag_im && flag_ip && flag_jm && flag_jp && flag_km && flag_kp;
}

double backgroundRatio(double* massArray, double* massUsageMask, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points)
{
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = (i_max-i_min);
    int j_range = (j_max-j_min);
    int k_range = (k_max-k_min);

    double totVolume = 0, backgroundVolume = 0;

    for  (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);

                if (massArray[idx] != massArray[idx])
                {
                    backgroundVolume += massUsageMask[idx];
                }

                totVolume += massUsageMask[idx];
            }
        }
    }
    return backgroundVolume/totVolume;
}

double targetMassAverage(double* massArray, double* massUsageMask, double* localSARArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double targetMass)
{ /* Given a mask (massUsageMask) of double 0<=value<=1 representing the fraction of used voxel and sum(massUsageMask * massArray = targetMass), it computes the mass average SAR according to the targetMass*/
    double avgSAR = 0;
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = (i_max-i_min);
    int j_range = (j_max-j_min);
    int k_range = (k_max-k_min);

    for (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                if (massArray[idx] == massArray[idx])
                    avgSAR += (massArray[idx] * massUsageMask[idx] * localSARArray[idx]);
            }
        }
    }

    return avgSAR/targetMass;
}

void markAndSetUSEDvox(double* avgSARArray, Status* voxStatusArray, double* massUsageMask, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double avgSARValueVALID)
{ /* If a voxel within the i,j,k limits is USED and its associated average SAR is lower than the avgSARValueVALID, the latter value is associated to that voxel*/
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = (i_max-i_min);
    int j_range = (j_max-j_min);
    int k_range = (k_max-k_min);

    for (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);
                if ((voxStatusArray[idx] == UNUSED || voxStatusArray[idx] == USED) && massUsageMask[idx]>TH_FOR_INC)
                {
                    if (voxStatusArray[idx] == UNUSED) {voxStatusArray[idx] = USED;}
                    if (avgSARArray[idx] < avgSARValueVALID)
                        avgSARArray[idx] = avgSARValueVALID;
                }
            }
        }
    }
    return;
}

double solveForK(double a, double b, double c, double d)
{ /* a*k^3 + b*k^2 + c*k + d = 0 
    It returns 0 if no solution exists between 0 and 1*/
    

    double p, q, delta, k, z_re, z_imm, cbrt_absz2, arg_z;

    if (a == 0)
    { 
        if (b != 0) // Second degree equation
        {
            delta = c*c-4*b*d;
            k = (-c+sqrt(delta))/(2*b);
            if (k>0 && k<=1)
                return k;
            else
            {
                k = (-c-sqrt(delta))/(2*b);
                if (k>0 && k<=1)
                    return k;
                else
                    #ifdef WARNINGS
                        printf("WARNING: No solution for at least one voxel. The SAR is averaged on the voxel mass");
                    #endif
                    return 0;
            }
        }
        else // b=0: First degree equation
        {
            k = -d/c; // Always between 0 and 1
            return k;
        }
    }

    #ifdef CM_ROOTS
        double m[9] = {0,1,0,0,0,1,-d/a,-c/a,-b/a};
        double wr[3], wi[3];
        int info;
        int n=3, lwork=18;
        double* work=(double *)malloc(lwork*sizeof(double *));
        char flag = 'N';

        dgeev_( &flag, &flag, &n, m,
                        &n, wr, wi, NULL,&n,
                        NULL, &n, work, &lwork, &info );

        free(work);
        
        if (!info)
        {
            for (int i=0; i<3; i++)
            {
                if (wi[i]==0 && wr[i]>0)
                    return wr[i];
            }
            #ifdef WARNINGS
                    printf("WARNING: No solution for at least one voxel. The SAR is averaged on the voxel mass\nk = %f, delta = %f\n", k, delta);
                #endif
            return 0;
        }
        else
        {
            #ifdef WARNINGS
                    printf("WARNING: dgeev failed on at least one voxel. The SAR is averaged on the voxel mass\nk = %f, delta = %f\n", k, delta);
            #endif
            return 0;
        }
        
    #else
        b = b/a;
        c = c/a;
        d = d/a;

        p = c-b*b/3;
        q = d+2*b*b*b/27-b*c/3;

        delta = q*q/4 + p*p*p/27;

        if (delta >= 0)
        {
            delta = sqrt(delta);
            k = cbrt(-q/2 + delta) + cbrt(-q/2 - delta) - b/3;
            if (k>0 && k<=1)
                return k;
            else
            {
                #ifdef WARNINGS
                    printf("WARNING: No solution for at least one voxel. The SAR is averaged on the voxel mass\nk = %f, delta = %f\n", k, delta);
                #endif
                return 0;
            }
        }
        else
        {
            z_re = -q/2;
            z_imm = sqrt(-delta);

            cbrt_absz2 = 2 * cbrt(sqrt(z_re*z_re + z_imm*z_imm));
            arg_z = atan2(z_imm, z_re);

            for (int i=0; i<3; i++)
            {
                k = cbrt_absz2 * cos(arg_z/3 + 2*i*PI/3) - b/3;
                if (k>0 && k<=1)
                    return k;
            }
        }
        #ifdef WARNINGS
            printf("WARNING: No solution for at least one voxel. The SAR is averaged on the voxel mass. delta: %f\n", delta);
        #endif
        return 0;
    #endif
}

void getAnnularMasses_CC(double* massArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double* masses)
{   /*

    It returns the masses of the 
        - vertex voxels: ma
        - edge voxels: mb
        - face voxels: mc
    for the centred averaging cube (step 1)

    */
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = (i_max-i_min);
    int j_range = (j_max-j_min);
    int k_range = (k_max-k_min);

    masses[0] = 0; // ma
    masses[1] = 0; // mb
    masses[2] = 0; // mc


    for (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);

                if (((idx_i == 0 || idx_i == i_range-1) && (idx_j == 0 || idx_j == j_range-1) && (idx_k == 0 || idx_k == k_range-1)) && (massArray[idx]==massArray[idx])) // Cube vertex
                {
                    masses[0] += massArray[idx];
                }
                else if ((((idx_i == 0 || idx_i == i_range-1) && idx_j >0 && idx_j < j_range-1 && idx_k >0 && idx_k < k_range-1) || // x faces
                ((idx_j == 0 || idx_j == j_range-1) && idx_i >0 && idx_i < i_range-1 && idx_k >0 && idx_k < k_range-1) || // j faces
                ((idx_k == 0 || idx_k == k_range-1) && idx_j >0 && idx_j < j_range-1 && idx_i >0 && idx_i < i_range-1)) // k faces
                 && (massArray[idx]==massArray[idx]))
                {
                    masses[2] += massArray[idx];
                }
                else if (!(idx_i > 0 && idx_i < i_range-1 && idx_j > 0 && idx_j < j_range-1 && idx_k > 0 && idx_k < k_range-1) && (massArray[idx]==massArray[idx])) // Is not a point inside the cube therefore it's an edge
                {
                    masses[1] += massArray[idx];
                }
            }
        }
    }
}

void getAnnularMasses_SC(double* massArray, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double* masses, int delta, int face)
{   /*

    It returns the masses of the 
        - vertex voxels: ma
        - edge voxels: mb
        - face voxels: mc

    for the side averaging cube (step 1)
    */

    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = (i_max-i_min);
    int j_range = (j_max-j_min);
    int k_range = (k_max-k_min);

    double mf = 0, me = 0, mv = 0, meo = 0, mfo = 0;

    masses[0] = 0; // ma -> k**3 contribution
    masses[1] = 0; // mb -> k**2 contribution
    masses[2] = 0; // mc -> k contribution


    for (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);

                switch (face)
                {
                    case 0: // x- --> x+
                    {
                        if (((idx_i == i_range-1) && (idx_j == 0 || idx_j == j_range-1) && (idx_k == 0 || idx_k == k_range-1)) && (massArray[idx]==massArray[idx])) // Cube vertex
                        {
                            mv += massArray[idx];
                        }
                        else if ((((idx_j == 0 || idx_j == j_range-1) && idx_i >=0 && idx_i < i_range-1 && idx_k >0 && idx_k < k_range-1) || // j faces
                        ((idx_k == 0 || idx_k == k_range-1) && idx_j >0 && idx_j < j_range-1 && idx_i >=0 && idx_i < i_range-1))             // k faces
                        && (massArray[idx]==massArray[idx]))
                        {
                            mf += massArray[idx];
                        }
                        else if ((( idx_i == i_range-1) && idx_j >0 && idx_j < j_range-1 && idx_k >0 && idx_k < k_range-1) && (massArray[idx]==massArray[idx])) // x+ face
                        {
                            mfo += massArray[idx];
                        }
                        else if (!(idx_i >= 0 && idx_i < i_range-1 && idx_j > 0 && idx_j < j_range-1 && idx_k > 0 && idx_k < k_range-1) && (massArray[idx]==massArray[idx])) // Is not a point inside the cube or inside the starting face therefore it's an edge
                        {   
                            if (idx_i == i_range-1) // Opposite edge
                            {
                                meo += massArray[idx];
                            }
                            else
                            {
                                me += massArray[idx];
                            }
                        }
                        break;
                    }

                    case 1: // x+ --> x-
                    {
                        if (((idx_i == 0) && (idx_j == 0 || idx_j == j_range-1) && (idx_k == 0 || idx_k == k_range-1)) && (massArray[idx]==massArray[idx])) // Cube vertex
                        {
                            mv += massArray[idx];
                        }
                        else if ((((idx_j == 0 || idx_j == j_range-1) && idx_i >0 && idx_i <= i_range-1 && idx_k >0 && idx_k < k_range-1) || // j faces
                                    ((idx_k == 0 || idx_k == k_range-1) && idx_j >0 && idx_j < j_range-1 && idx_i >0 && idx_i <= i_range-1)) // k faces
                                    && (massArray[idx]==massArray[idx]))
                        {
                            mf += massArray[idx];
                        }
                        else if (((idx_i == 0) && idx_j >0 && idx_j < j_range-1 && idx_k >0 && idx_k < k_range-1) && massArray[idx]==massArray[idx]) // x- face
                        {
                            mfo += massArray[idx];
                        }
                        else if (!(idx_i > 0 && idx_i <= i_range-1 && idx_j > 0 && idx_j < j_range-1 && idx_k > 0 && idx_k < k_range-1) && (massArray[idx]==massArray[idx])) // Is not a point inside the cube or inside the starting face therefore it's an edge
                        {
                            if (idx_i == 0) // Opposite edge
                            {
                                meo += massArray[idx];
                            }
                            else
                            {
                                me += massArray[idx];
                            }
                        }
                        break;
                    }

                    case 2: // y- --> y+
                    {
                        if (((idx_i == 0 || idx_i == i_range-1) && (idx_j == j_range-1) && (idx_k == 0 || idx_k == k_range-1)) && (massArray[idx]==massArray[idx])) // Cube vertex
                        {
                            mv += massArray[idx];
                        }
                        else if ((((idx_i == 0 || idx_i == i_range-1) && idx_j >= 0 && idx_j < j_range-1 && idx_k >0 && idx_k < k_range-1) || // x faces
                                    ((idx_k == 0 || idx_k == k_range-1) && idx_j >= 0 && idx_j < j_range-1 && idx_i >0 && idx_i < i_range-1)) // k faces
                         && (massArray[idx]==massArray[idx]))
                        {
                            mf += massArray[idx];
                        }
                        else if (((idx_j == j_range-1) && idx_i >0 && idx_i < i_range-1 && idx_k >0 && idx_k < k_range-1) && massArray[idx]==massArray[idx])  // j+ face
                        {
                            mfo += massArray[idx];
                        }
                        else if (!(idx_i > 0 && idx_i < i_range-1 && idx_j >= 0 && idx_j < j_range-1 && idx_k > 0 && idx_k < k_range-1) && (massArray[idx]==massArray[idx])) // Is not a point inside the cube or on inside starting face therefore it's an edge
                        {
                            if (idx_j == j_range-1) // Opposite edge
                            {
                                meo += massArray[idx];
                            }
                            else
                            {
                                me += massArray[idx];
                            }
                        }
                        break;
                    }

                    case 3: // y+ --> y-
                    {
                        if (((idx_i == 0 || idx_i == i_range-1) && (idx_j == 0) && (idx_k == 0 || idx_k == k_range-1)) && (massArray[idx]==massArray[idx])) // Cube vertex
                        {
                            mv += massArray[idx];
                        }
                        else if ((((idx_i == 0 || idx_i == i_range-1) && idx_j >0 && idx_j <= j_range-1 && idx_k >0 && idx_k < k_range-1) || // x faces
                        ((idx_k == 0 || idx_k == k_range-1) && idx_j >0 && idx_j <= j_range-1 && idx_i >0 && idx_i < i_range-1)) // k faces
                         && (massArray[idx]==massArray[idx]))
                        {
                             mf += massArray[idx];
                        }
                        else if (((idx_j == 0) && idx_i >0 && idx_i < i_range-1 && idx_k >0 && idx_k < k_range-1) && (massArray[idx]==massArray[idx])) // j- face
                        {
                            mfo += massArray[idx];
                        }
                        else if (!(idx_i > 0 && idx_i < i_range-1 && idx_j > 0 && idx_j <= j_range-1 && idx_k > 0 && idx_k < k_range-1) && (massArray[idx]==massArray[idx])) // Is not a point inside the cube or insed the starting face therefore it's an edge
                        {
                            if (idx_j == 0) // Opposite edge
                            {
                                meo += massArray[idx];
                            }
                            else
                            {
                                me += massArray[idx];
                            }
                        }
                        break;
                    }

                    case 4: // z- --> z+
                    {
                        if (((idx_i == 0 || idx_i == i_range-1) && (idx_j == 0 || idx_j == j_range-1) && (idx_k == k_range-1)) && (massArray[idx]==massArray[idx])) // Cube vertex
                        {
                            mv += massArray[idx];
                        }
                        else if ((((idx_i == 0 || idx_i == i_range-1) && idx_j >0 && idx_j < j_range-1 && idx_k >= 0 && idx_k < k_range-1) || // x faces
                                ((idx_j == 0 || idx_j == j_range-1) && idx_i >0 && idx_i < i_range-1 && idx_k >= 0 && idx_k < k_range-1)) // j faces
                                && (massArray[idx]==massArray[idx]))
                        {
                            mf += massArray[idx];
                        }
                        else if (((idx_k == k_range-1) && idx_j >0 && idx_j < j_range-1 && idx_i >0 && idx_i < i_range-1) && (massArray[idx]==massArray[idx])) // k+ face
                        {
                            mfo += massArray[idx];
                        }
                        else if (!(idx_i > 0 && idx_i < i_range-1 && idx_j > 0 && idx_j < j_range-1 && idx_k >= 0 && idx_k < k_range-1) && (massArray[idx]==massArray[idx])) // Is not a point inside the cube or inside the starting face therefore it's an edge
                        {
                            if (idx_k == k_range-1) // Opposite edge
                            {
                                meo += massArray[idx];
                            }
                            else
                            {
                                me += massArray[idx];
                            }
                        }
                        break;
                    }

                    case 5: // z+ --> z-
                    {
                        if (((idx_i == 0 || idx_i == i_range-1) && (idx_j == 0 || idx_j == j_range-1) && (idx_k == 0)) && (massArray[idx]==massArray[idx])) // Cube vertex
                        {
                            mv += massArray[idx];
                        }
                        else if ((((idx_i == 0 || idx_i == i_range-1) && idx_j >0 && idx_j < j_range-1 && idx_k > 0 && idx_k <= k_range-1) || // x faces
                                ((idx_j == 0 || idx_j == j_range-1) && idx_i >0 && idx_i < i_range-1 && idx_k > 0 && idx_k <= k_range-1)) // j faces
                                && (massArray[idx]==massArray[idx]))
                        {
                            mf += massArray[idx];
                        }
                        else if (((idx_k == 0) && idx_j >0 && idx_j < j_range-1 && idx_i >0 && idx_i < i_range-1) && (massArray[idx]==massArray[idx])) // k- face
                        {
                            mfo += massArray[idx];
                        }
                        else if (!(idx_i > 0 && idx_i < i_range-1 && idx_j > 0 && idx_j < j_range-1 && idx_k > 0 && idx_k <= k_range-1) && (massArray[idx]==massArray[idx])) // Is not a point inside the cube or inside the starting face therefore it's an edge
                        {
                            if (idx_k == 0) // Opposite edge
                            {
                                meo += massArray[idx];
                            }
                            else
                            {
                                me += massArray[idx];
                            }
                        }
                        break;
                    }
                }
            }
        }
    }

    if (delta%2==1) // delta is odd
    {
        masses[0] = 0.25 * mv;
        masses[1] = 0.25*me + 0.5*meo;
        masses[2] = 0.5*mf + mfo;
    }
    else //delta is even
    {
        masses[0] = 0.25 * mv;
        masses[1] = me*0.25 + meo*0.5 + mv*0.5;
        masses[2] = 0.5*mf + me*0.5 + mv*0.25 + meo*0.5 + mfo;
    }
}

void setPartialMassUsageMask_CC(double* massUsageMask, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double k)
{
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = (i_max-i_min);
    int j_range = (j_max-j_min);
    int k_range = (k_max-k_min);

    for (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);

                if ((idx_i == 0 || idx_i == i_range-1) && (idx_j == 0 || idx_j == j_range-1) && (idx_k == 0 || idx_k == k_range-1)) // Cube vertex
                {
                    massUsageMask[idx] = k*k*k;
                }
                else if (((idx_i == 0 || idx_i == i_range-1) && idx_j >0 && idx_j < j_range-1 && idx_k >0 && idx_k < k_range-1) || // x faces
                ((idx_j == 0 || idx_j == j_range-1) && idx_i >0 && idx_i < i_range-1 && idx_k >0 && idx_k < k_range-1) || // j faces
                ((idx_k == 0 || idx_k == k_range-1) && idx_j >0 && idx_j < j_range-1 && idx_i >0 && idx_i < i_range-1)) // k faces
                {
                    massUsageMask[idx] = k;
                }
                else if (!(idx_i >0 && idx_i < i_range-1 && idx_j >0 && idx_j < j_range-1 && idx_k >0 && idx_k < k_range-1)) // Is not a point inside the cube therefore it's an edge
                {
                   massUsageMask[idx] = k*k;
                }
            }
        }
    }
}


void setPartialMassUsageMask_SC(double* massUsageMask, int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int* n_points, double k, int delta, int face)
{
    int idx_j;
    int idx_k;
    int idx_t1, idx_t2, idx;

    int delta_j = n_points[1];
    int delta_k = n_points[2];

    int i_range = (i_max-i_min);
    int j_range = (j_max-j_min);
    int k_range = (k_max-k_min);

    for (int idx_i=0; idx_i<i_range; idx_i++)
    {
        idx_t1 = (idx_i+i_min)*delta_j*delta_k;

        for  (int idx_j=0; idx_j<j_range; idx_j++)
        {
            idx_t2 = idx_t1 + (idx_j+j_min)*delta_k;

            for  (int idx_k=0; idx_k<k_range; idx_k++)
            {
                idx = idx_t2 + (idx_k+k_min);

                switch (face)
                {
                    case 0: // x- --> x+
                    {
                        if ((idx_i == i_range-1) && (idx_j == 0 || idx_j == j_range-1) && (idx_k == 0 || idx_k == k_range-1)) // Cube vertex
                        {
                            massUsageMask[idx] = (0.5*k)*(0.5*k)*k * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k)*k * ((delta+1)%2);
                        }
                        else if (((idx_j == 0 || idx_j == j_range-1) && idx_i >=0 && idx_i < i_range-1 && idx_k >0 && idx_k < k_range-1) || // j faces
                        ((idx_k == 0 || idx_k == k_range-1) && idx_j >0 && idx_j < j_range-1 && idx_i >=0 && idx_i < i_range-1))            // k faces
                        {
                            massUsageMask[idx] = 0.5 * ((delta+1)%2) + 0.5*k;
                        }
                        else if (( idx_i == i_range-1) && idx_j >0 && idx_j < j_range-1 && idx_k >0 && idx_k < k_range-1) // x+ face
                        {
                            massUsageMask[idx] = k;
                        }
                        else if (!(idx_i >=0 && idx_i < i_range-1 && idx_j >0 && idx_j < j_range-1 && idx_k >0 && idx_k < k_range-1)) // Is not a point inside the cube or on the starting face therefore it's an edge
                        {
                            if (idx_i == i_range-1) // Opposite edge
                            {
                                massUsageMask[idx] = 0.5*k*k * (delta%2) + (0.5+0.5*k)*k * ((delta+1)%2);
                            }
                            else
                            {
                                massUsageMask[idx] = (0.5*k)*(0.5*k) * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k) * ((delta+1)%2);
                            }
                        }
                        break;
                    }

                    case 1: // x+ --> x-
                    {
                        if ((idx_i == 0) && (idx_j == 0 || idx_j == j_range-1) && (idx_k == 0 || idx_k == k_range-1)) // Cube vertex
                        {
                            massUsageMask[idx] = (0.5*k)*(0.5*k)*k * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k)*k * ((delta+1)%2);
                        }
                        else if (((idx_j == 0 || idx_j == j_range-1) && idx_i >0 && idx_i <= i_range-1 && idx_k >0 && idx_k < k_range-1) || // j faces
                                ((idx_k == 0 || idx_k == k_range-1) && idx_j >0 && idx_j < j_range-1 && idx_i >0 && idx_i <= i_range-1)) // k faces
                        {
                            massUsageMask[idx] = 0.5 * ((delta+1)%2) + 0.5*k;
                        }
                        else if ((idx_i == 0) && idx_j >0 && idx_j < j_range-1 && idx_k >0 && idx_k < k_range-1) // x- face
                        {
                            massUsageMask[idx] = k;
                        }
                        else if (!(idx_i >0 && idx_i <= i_range-1 && idx_j >0 && idx_j < j_range-1 && idx_k >0 && idx_k < k_range-1)) // Is not a point inside the cube or on the starting face therefore it's an edge
                        {
                            if (idx_i == 0) // Opposite edge
                            {
                                massUsageMask[idx] = 0.5*k*k * (delta%2) + (0.5+0.5*k)*k * ((delta+1)%2);
                            }
                            else
                            {
                                massUsageMask[idx] = (0.5*k)*(0.5*k) * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k) * ((delta+1)%2);
                            }
                        }
                        break;
                    }

                    case 2: // y- --> y+
                    {
                        if (((idx_i == 0 || idx_i == i_range-1) && (idx_j == j_range-1) && (idx_k == 0 || idx_k == k_range-1))) // Cube vertex
                        {
                            massUsageMask[idx] = (0.5*k)*(0.5*k)*k * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k)*k * ((delta+1)%2);
                        }
                        else if ((((idx_i == 0 || idx_i == i_range-1) && idx_j >= 0 && idx_j < j_range-1 && idx_k >0 && idx_k < k_range-1) || // x faces
                                    ((idx_k == 0 || idx_k == k_range-1) && idx_j >= 0 && idx_j < j_range-1 && idx_i >0 && idx_i < i_range-1))) // k faces
                        {
                            massUsageMask[idx] = 0.5 * ((delta+1)%2) + 0.5*k;
                        }
                        else if (((idx_j == j_range-1) && idx_i >0 && idx_i < i_range-1 && idx_k >0 && idx_k < k_range-1))  // j+ face
                        {
                            massUsageMask[idx] = k;
                        }
                        else if (!(idx_i > 0 && idx_i < i_range-1 && idx_j >= 0 && idx_j < j_range-1 && idx_k > 0 && idx_k < k_range-1)) // Is not a point inside the cube or on inside starting face therefore it's an edge
                        {
                             if (idx_j == j_range-1) // Opposite edge
                            {
                                massUsageMask[idx] = 0.5*k*k * (delta%2) + (0.5+0.5*k)*k * ((delta+1)%2);
                            }
                            else
                            {
                                massUsageMask[idx] = (0.5*k)*(0.5*k) * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k) * ((delta+1)%2);
                            }
                        }
                        break;
                    }

                    case 3: // y+ --> y-
                    {
                        if (((idx_i == 0 || idx_i == i_range-1) && (idx_j == 0) && (idx_k == 0 || idx_k == k_range-1))) // Cube vertex
                        {
                            massUsageMask[idx] = (0.5*k)*(0.5*k)*k * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k)*k * ((delta+1)%2);
                        }
                        else if ((((idx_i == 0 || idx_i == i_range-1) && idx_j >0 && idx_j <= j_range-1 && idx_k >0 && idx_k < k_range-1) || // x faces
                        ((idx_k == 0 || idx_k == k_range-1) && idx_j >0 && idx_j <= j_range-1 && idx_i >0 && idx_i < i_range-1))) // k faces
                        {
                            massUsageMask[idx] = 0.5 * ((delta+1)%2) + 0.5*k;
                        }
                        else if (((idx_j == 0) && idx_i >0 && idx_i < i_range-1 && idx_k >0 && idx_k < k_range-1)) // j- face
                        {
                            massUsageMask[idx] = k;
                        }
                        else if (!(idx_i > 0 && idx_i < i_range-1 && idx_j > 0 && idx_j <= j_range-1 && idx_k > 0 && idx_k < k_range-1)) // Is not a point inside the cube or inside the starting face therefore it's an edge
                        {
                            if (idx_j == 0) // Opposite edge
                            {
                                massUsageMask[idx] = 0.5*k*k * (delta%2) + (0.5+0.5*k)*k * ((delta+1)%2);
                            }
                            else
                            {
                                massUsageMask[idx] = (0.5*k)*(0.5*k) * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k) * ((delta+1)%2);
                            }
                        }
                        break;
                    }

                    case 4: // z- --> z+
                    {
                        if (((idx_i == 0 || idx_i == i_range-1) && (idx_j == 0 || idx_j == j_range-1) && (idx_k == k_range-1))) // Cube vertex
                        {
                            massUsageMask[idx] = (0.5*k)*(0.5*k)*k * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k)*k * ((delta+1)%2);
                        }
                        else if (((idx_i == 0 || idx_i == i_range-1) && idx_j >0 && idx_j < j_range-1 && idx_k >= 0 && idx_k < k_range-1) || // x faces
                                ((idx_j == 0 || idx_j == j_range-1) && idx_i >0 && idx_i < i_range-1 && idx_k >= 0 && idx_k < k_range-1)) // j faces
                        {
                            massUsageMask[idx] = 0.5 * ((delta+1)%2) + 0.5*k;
                        }
                        else if (((idx_k == k_range-1) && idx_j >0 && idx_j < j_range-1 && idx_i >0 && idx_i < i_range-1)) // k+ face
                        {
                            massUsageMask[idx] = k;
                        }
                        else if (!(idx_i > 0 && idx_i < i_range-1 && idx_j > 0 && idx_j < j_range-1 && idx_k >= 0 && idx_k < k_range-1)) // Is not a point inside the cube or inside the starting face therefore it's an edge
                        {
                            if (idx_k == k_range-1) // Opposite edge
                            {
                                massUsageMask[idx] = 0.5*k*k * (delta%2) + (0.5+0.5*k)*k * ((delta+1)%2);
                            }
                            else
                            {
                                massUsageMask[idx] = (0.5*k)*(0.5*k) * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k) * ((delta+1)%2);
                            }
                        }
                        break;
                    }

                    case 5: // z+ --> z-
                    {
                        if (((idx_i == 0 || idx_i == i_range-1) && (idx_j == 0 || idx_j == j_range-1) && (idx_k == 0))) // Cube vertex
                        {
                            massUsageMask[idx] = (0.5*k)*(0.5*k)*k * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k)*k * ((delta+1)%2);
                        }
                        else if (((idx_i == 0 || idx_i == i_range-1) && idx_j >0 && idx_j < j_range-1 && idx_k > 0 && idx_k <= k_range-1) || // x faces
                                ((idx_j == 0 || idx_j == j_range-1) && idx_i >0 && idx_i < i_range-1 && idx_k > 0 && idx_k <= k_range-1)) // j faces
                        {
                            massUsageMask[idx] = 0.5 * ((delta+1)%2) + 0.5*k;
                        }
                        else if (((idx_k == 0) && idx_j >0 && idx_j < j_range-1 && idx_i >0 && idx_i < i_range-1)) // k- face
                        {
                            massUsageMask[idx] = k;
                        }
                        else if (!(idx_i > 0 && idx_i < i_range-1 && idx_j > 0 && idx_j < j_range-1 && idx_k > 0 && idx_k <= k_range-1)) // Is not a point inside the cube or inside the starting face therefore it's an edge
                        {
                            if (idx_k == 0) // Opposite edge
                            {
                                massUsageMask[idx] = 0.5*k*k * (delta%2) + (0.5+0.5*k)*k * ((delta+1)%2);
                            }
                            else
                            {
                                massUsageMask[idx] = (0.5*k)*(0.5*k) * (delta%2) + (0.5+0.5*k)*(0.5+0.5*k) * ((delta+1)%2);
                            }
                        }
                        break;
                    }
                }
            }
        }
    }
}

double nanmax_f(double* array, int length)
{
    double maxValue = 0;

    for (int i=0; i<length; i++)
    {
        if (array[i]==array[i] && array[i]>maxValue)
            maxValue = array[i];
    }

    return maxValue;
}