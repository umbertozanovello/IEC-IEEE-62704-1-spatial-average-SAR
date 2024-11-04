# include "include/auxFuncs.h"
# include "include/globals.h"

int averagingMaskCenter3D(int* ic, int* n_points, double* massArray, int delta_min, double targetMass, Status* voxStatusArray, double* massUsageMask, double* annularMask, int* ss)
{
    uint8_t stopFlag = 0;
    double cumMass, prevCumMass, k;
    uint8_t delta = delta_min;
    double masses[3] = {0,0,0};

    /* Preset voxels within delta_min*/
    if (delta>0) // If delta == 0, there's nothing to preset
    {
        for (int i=0; i<DIMS; i++)  
        {
            if ((ic[i] - delta<0) || (ic[i] + delta>=n_points[i])) // The averaging cube goes outside the calculation domain. There's therefore at least face of the averaging cube does not intersect or touch tissues
            {
                return 0;
            }
        }

        for (int i=0; i<DIMS; i++)
        {
            ss[2*i] = ic[i] - delta;
            ss[2*i+1] = ic[i] + delta + 1;
        }

        setValue_f(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 1); // I set to full used all the voxel I'm sure will sum to a mass lower than the target mass
        prevCumMass = nansum_f(massArray, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points);
        cumMass = prevCumMass;

        if (cumMass == targetMass) {stopFlag = 1;} // I already identified the correct averaging cube

        delta++;
    }

    while(!stopFlag)
    {
        for (int i=0; i<DIMS; i++)
            {
                if ((ic[i] - delta<0) || (ic[i] + delta>=n_points[i])) // The averaging cube goes outside the calculation domain. There's therefore at least face of the averaging cube does not intersect or touch tissues
                {
                    return 0;
                }
            }

        for (int i=0; i<DIMS; i++)
        {
            ss[2*i] = ic[i] - delta;
            ss[2*i+1] = ic[i] + delta + 1;
        }

        if (!checkNanFaces(massArray, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points)) // The averaging cube has at least one face intersecting only background
            return 0;

        setValue_f(annularMask, ic[0]-delta, -1, ic[1]-delta, ic[1]+delta+1, ic[2]-delta, ic[2]+delta+1, n_points, 1); // x-
        setValue_f(annularMask, ic[0]+delta, -1, ic[1]-delta, ic[1]+delta+1, ic[2]-delta, ic[2]+delta+1, n_points, 1); // x+
        setValue_f(annularMask, ic[0]-delta, ic[0]+delta+1, ic[1]-delta, -1, ic[2]-delta, ic[2]+delta+1, n_points, 1); // y-
        setValue_f(annularMask, ic[0]-delta, ic[0]+delta+1, ic[1]+delta, -1, ic[2]-delta, ic[2]+delta+1, n_points, 1); // y+
        setValue_f(annularMask, ic[0]-delta, ic[0]+delta+1, ic[1]-delta, ic[1]+delta+1, ic[2]-delta, -1, n_points, 1); // z-
        setValue_f(annularMask, ic[0]-delta, ic[0]+delta+1, ic[1]-delta, ic[1]+delta+1, ic[2]+delta, -1, n_points, 1); // z+

        cumMass += nansum_ifPositive_ff(massArray, annularMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points);

        if ((cumMass-targetMass)/targetMass < -(AVG_VOL_TOLL)) // cumMass is lower than target mass and outside the accepted tollerance
        {
            setValue_ifPositive_ff(massUsageMask, annularMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 1);
        }
        else if ((cumMass-targetMass)/targetMass > abs(AVG_VOL_TOLL)) // cumMass is higher than target mass and outside the accepted tollerance
        {
            if (delta == 0)
                massUsageMask[from3DIdxTo1DIdx(ic, n_points)] = targetMass / massArray[from3DIdxTo1DIdx(ic, n_points)];
            else
            {
                getAnnularMasses_CC(massArray, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, masses);
                k = solveForK(masses[0], masses[1], masses[2], prevCumMass-targetMass);
                
                #ifdef DEBUG
                    printf("ic: %d,%d,%d\n",ic[0],ic[1],ic[2]);
                    printf("ma: %f, mb: %f, mc: %f, d: %f\n",masses[0], masses[1], masses[2], prevCumMass-targetMass);
                    printf("k: %.10f\n",k);
                #endif
                
                setPartialMassUsageMask_CC(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, k);
            }
            stopFlag = 1;
        }
        else
        {
            setValue_ifPositive_ff(massUsageMask, annularMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 1);
            stopFlag = 1;
        }
        
        // Set again annularMask to all zeros
        setValue_f(annularMask, ic[0]-delta, -1, ic[1]-delta, ic[1]+delta+1, ic[2]-delta, ic[2]+delta+1, n_points, 0); // x-
        setValue_f(annularMask, ic[0]+delta, -1, ic[1]-delta, ic[1]+delta+1, ic[2]-delta, ic[2]+delta+1, n_points, 0); // x+
        setValue_f(annularMask, ic[0]-delta, ic[0]+delta+1, ic[1]-delta, -1, ic[2]-delta, ic[2]+delta+1, n_points, 0); // y-
        setValue_f(annularMask, ic[0]-delta, ic[0]+delta+1, ic[1]+delta, -1, ic[2]-delta, ic[2]+delta+1, n_points, 0); // y+
        setValue_f(annularMask, ic[0]-delta, ic[0]+delta+1, ic[1]-delta, ic[1]+delta+1, ic[2]-delta, -1, n_points, 0); // z-
        setValue_f(annularMask, ic[0]-delta, ic[0]+delta+1, ic[1]-delta, ic[1]+delta+1, ic[2]+delta, -1, n_points, 0); // z+

        prevCumMass = cumMass;
        delta++;
    }

    // Voxel status check

    if (backgroundRatio(massArray, massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points)<=MAX_BACKGROUND_RATIO)
    {
        // markUsedVoxels(voxStatusArray, massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points);
        if (massUsageMask[from3DIdxTo1DIdx(ic, n_points)]>TH_FOR_INC)
        {
            voxStatusArray[from3DIdxTo1DIdx(ic, n_points)] = VALID;
        }
        return 1;
    }
    else
    {
        return 0;
    }
}

int main(double* massArray, double* localSARArray, double* avgSARArray, Status* voxStatusArray, double targetMass, int* n_points)
{
    int tot_points = n_points[0] * n_points[1] * n_points[2];
    int idx = 0;
    int ss[2*DIMS] = {0, 0, 0, 0, 0, 0}; // slice of array as [i_min,i_max,j_min,j_max,k_min,k_max];
    int ic[3] = {0,0,0};
    int delta_min, retValue;
    double* massUsageMask = calloc(tot_points, sizeof(double));
    double* annularMask = calloc(tot_points, sizeof(double));

    #ifdef REPORT
        FILE* report_ptr;
        report_ptr = fopen(REPORT, "w"); 
        double avgMass=0, avgVolume=0;
    #endif

    delta_min = floor((cbrt(targetMass / nanmax_f(massArray, tot_points))-1)/2);
    
    for (int i=0; i<n_points[0]; i++)
    {
        ic[0] = i;
        for (int j=0; j<n_points[1]; j++)
        {
            ic[1] = j;
            for (int k=0; k<n_points[2]; k++)
            {
                ic[2] = k;

                if (massArray[idx] == massArray[idx]) // massArray[i] is not nan (i.e., background)
                {
                    retValue = averagingMaskCenter3D(ic, n_points, massArray, delta_min, targetMass, voxStatusArray, massUsageMask, annularMask, ss);

                    if (retValue==1) // The averaging cube has been defined
                    {
                        avgSARArray[idx] = targetMassAverage(massArray, massUsageMask, localSARArray, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, targetMass);
                        markAndSetUSEDvox(avgSARArray, voxStatusArray, massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, avgSARArray[idx]);

                        #ifdef REPORT
                            // VALID voxels
                            avgVolume = nansum_f(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points)*VOX_VOL;
                            avgMass = nansum_product_ff(massArray, massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points); // kg
                            fprintf(report_ptr, "%d %d %d %d %d %.6e %.6e %d %.6e %.6e\n", idx, i, j, k, voxStatusArray[idx], avgMass*1e3, avgVolume, 7, localSARArray[idx], avgSARArray[idx]);
                        #endif
                    }
                    setValue_f(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 0); // I set to zero the values I already know are not zero

                }
                #ifdef DEBUG
                //printf("\rStep one: %.2f %%", 100.0*idx/tot_points);
                #endif
                idx++;
            }
        }
    }

    free(massUsageMask);
    free(annularMask);
    #ifdef REPORT
        // USED voxels
        avgVolume = 0;
        avgMass = 0;
        idx = 0;
        for (int i=0; i<n_points[0]; i++)
        {
            ic[0] = i;
            for (int j=0; j<n_points[1]; j++)
            {
                ic[1] = j;
                for (int k=0; k<n_points[2]; k++)
                {
                    ic[2] = k;
                    if (massArray[idx]==massArray[idx] && voxStatusArray[idx] == USED)
                        {
                        fprintf(report_ptr, "%d %d %d %d %d %.6e %.6e %d %.6e %.6e\n", idx, i, j, k, voxStatusArray[idx], avgMass, avgVolume, 0, localSARArray[idx], avgSARArray[idx]);
                        }
                    idx++;
                }
            }
        }
        fclose(report_ptr);
    #endif
    return 0;
}