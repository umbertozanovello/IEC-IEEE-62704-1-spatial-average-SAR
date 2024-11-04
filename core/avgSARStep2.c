# include "include/auxFuncs.h"
# include "include/globals.h"

int averagingMaskSide3D(int* ic, int* n_points, double* massArray, int delta_min, double targetMass, int face, double* massUsageMask, double* annularMask, int* ss)
{
    int delta = delta_min, stopFlag = 0;
    double cumMass = 0, prevCumMass = 0, k;
    double masses[3] = {0,0,0};

    /* I preset the massUsageMask within delta_min */
    if (delta>0) // If delta == 0, there's nothing to preset
    {
        switch (face)
        {
            case 0: // side x- -> x+
                if (ic[0]+delta >= n_points[0] || ic[1]-ceil(delta/2.) < 0 ||  ic[1]+ceil(delta/2.)>=n_points[1] || ic[2]-ceil(delta/2.) < 0 ||  ic[2]+ceil(delta/2.)>=n_points[2]) // Outside the domain
                {
                    return 0;
                }

                ss[0] = (int)(ic[0]);
                ss[1] = (int)(ic[0]+delta+1);
                ss[2] = (int)(ic[1]-(int)(ceil(delta/2.)));
                ss[3] = (int)(ic[1]+(int)(ceil(delta/2.))+1);
                ss[4] = (int)(ic[2]-(int)(ceil(delta/2.)));
                ss[5] = (int)(ic[2]+(int)(ceil(delta/2.))+1);

                setValue_f(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 1); // x+ face
                setValue_f(massUsageMask, ic[0], ic[0]+delta+delta%2, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // y- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0], ic[0]+delta+delta%2, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // y+ face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0], ic[0]+delta+delta%2, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // z- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0], ic[0]+delta+delta%2, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // z+ face, 1 if delta is even, 0.5 if delta is odd

                setValue_f(massUsageMask, ic[0], ic[0]+delta+delta%2, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // y-,z- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0], ic[0]+delta+delta%2, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // y-,z+ edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0], ic[0]+delta+delta%2, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // y+,z- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0], ic[0]+delta+delta%2, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // y+,z+ edge, 1 if delta is even, 0.25 if delta is odd

                break;
            
            case 1: // side x+ -> x-
                if (ic[0]-delta < 0 || ic[1]-ceil(delta/2.) < 0 ||  ic[1]+ceil(delta/2.)>=n_points[1] || ic[2]-ceil(delta/2.) < 0 ||  ic[2]+ceil(delta/2.)>=n_points[2]) // Outside the domain
                {
                    return 0;
                }

                ss[0] = (int)(ic[0]-delta);
                ss[1] = (int)(ic[0]+1);
                ss[2] = (int)(ic[1]-(int)(ceil(delta/2.)));
                ss[3] = (int)(ic[1]+(int)(ceil(delta/2.))+1);
                ss[4] = (int)(ic[2]-(int)(ceil(delta/2.)));
                ss[5] = (int)(ic[2]+(int)(ceil(delta/2.))+1);

                setValue_f(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 1); // x- face
                setValue_f(massUsageMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // y- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // y+ face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // z- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // z+ face, 1 if delta is even, 0.5 if delta is odd

                setValue_f(massUsageMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // y-,z- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // y-,z+ edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // y+,z- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // y+,z+ edge, 1 if delta is even, 0.25 if delta is odd

                break;

            case 2: // side y- -> y+
                if (ic[1]+delta >= n_points[1] || ic[0]-ceil(delta/2.) < 0 ||  ic[0]+ceil(delta/2.)>=n_points[0] || ic[2]-ceil(delta/2.) < 0 ||  ic[2]+ceil(delta/2.)>=n_points[2]) // Outside the domain
                {
                    return 0;
                }

                ss[0] = (int)(ic[0]-(int)(ceil(delta/2.)));
                ss[1] = (int)(ic[0]+(int)(ceil(delta/2.))+1);
                ss[2] = (int)(ic[1]);
                ss[3] = (int)(ic[1]+delta+1);
                ss[4] = (int)(ic[2]-(int)(ceil(delta/2.)));
                ss[5] = (int)(ic[2]+(int)(ceil(delta/2.))+1);

                setValue_f(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 1); // y- face
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // x- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // x+ face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1], ic[1]+delta+delta%2, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // z- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1], ic[1]+delta+delta%2, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // z+ face, 1 if delta is even, 0.5 if delta is odd

                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x-,z- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x-,z+ edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x+,z- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x+,z+ edge, 1 if delta is even, 0.25 if delta is odd

                break;

            case 3: // side y+ -> y-
                if (ic[1]-delta < 0 || ic[0]-ceil(delta/2.) < 0 ||  ic[0]+ceil(delta/2.)>=n_points[0] || ic[2]-ceil(delta/2.) < 0 ||  ic[2]+ceil(delta/2.)>=n_points[2]) // Outside the domain
                {
                    return 0;
                }

                ss[0] = (int)(ic[0]-(int)(ceil(delta/2.)));
                ss[1] = (int)(ic[0]+(int)(ceil(delta/2.))+1);
                ss[2] = (int)(ic[1]-delta);
                ss[3] = (int)(ic[1]+1);
                ss[4] = (int)(ic[2]-(int)(ceil(delta/2.)));
                ss[5] = (int)(ic[2]+(int)(ceil(delta/2.))+1);

                setValue_f(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 1); // y- face
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // x- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // x+ face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // z- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // z+ face, 1 if delta is even, 0.5 if delta is odd

                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x-,z- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x-,z+ edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x+,z- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x+,z+ edge, 1 if delta is even, 0.25 if delta is odd

                break;
            
            case 4: // side z- -> z+
                if (ic[2]+delta >= n_points[2] || ic[0]-ceil(delta/2.) < 0 ||  ic[0]+ceil(delta/2.)>=n_points[0] || ic[1]-ceil(delta/2.) < 0 ||  ic[1]+ceil(delta/2.)>=n_points[1]) // Outside the domain
                {
                    return 0;
                }

                ss[0] = (int)(ic[0]-(int)(ceil(delta/2.)));
                ss[1] = (int)(ic[0]+(int)(ceil(delta/2.))+1);
                ss[2] = (int)(ic[1]-(int)(ceil(delta/2.)));
                ss[3] = (int)(ic[1]+(int)(ceil(delta/2.))+1);
                ss[4] = (int)(ic[2]);
                ss[5] = (int)(ic[2]+delta+1);

                setValue_f(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 1); // z+ face
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2], ic[2]+delta+delta%2, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // # x- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2], ic[2]+delta+delta%2, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // # x+ face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // # y- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // # y+ face, 1 if delta is even, 0.5 if delta is odd

                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x-,y- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x-,y+ edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x+,y- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x+,y+ edge, 1 if delta is even, 0.25 if delta is odd

                break;

            case 5: // side z+ -> z-
                if (ic[2]-delta < 0 || ic[0]-ceil(delta/2.) < 0 ||  ic[0]+ceil(delta/2.)>=n_points[0] || ic[1]-ceil(delta/2.) < 0 ||  ic[1]+ceil(delta/2.)>=n_points[1]) // Outside the domain
                {
                    return 0;
                }

                ss[0] = (int)(ic[0]-(int)(ceil(delta/2.)));
                ss[1] = (int)(ic[0]+(int)(ceil(delta/2.))+1);
                ss[2] = (int)(ic[1]-(int)(ceil(delta/2.)));
                ss[3] = (int)(ic[1]+(int)(ceil(delta/2.))+1);
                ss[4] = (int)(ic[2]-delta);
                ss[5] = (int)(ic[2]+1);

                setValue_f(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 1); // z- face
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // # x- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // # x+ face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // # y- face, 1 if delta is even, 0.5 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.5*(delta%2) + 1*(abs((delta-1)%2))); // # y+ face, 1 if delta is even, 0.5 if delta is odd

                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x-,y- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x-,y+ edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x+,y- edge, 1 if delta is even, 0.25 if delta is odd
                setValue_f(massUsageMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.25*(delta%2) + 1*(abs((delta-1)%2))); // x+,y+ edge, 1 if delta is even, 0.25 if delta is odd

                break;

        }

        prevCumMass = nansum_product_ff(massArray, massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points);
        cumMass = prevCumMass;

        if (cumMass == targetMass) {return 1;} // I already identified the correct averaging cube

        delta++;

    }

    while (!stopFlag)
    {
        switch (face)
        {
            case 0: // side x- -> x+
                if (ic[0]+delta >= n_points[0] || ic[1]-ceil(delta/2.) < 0 ||  ic[1]+ceil(delta/2.)>=n_points[1] || ic[2]-ceil(delta/2.) < 0 ||  ic[2]+ceil(delta/2.)>=n_points[2]) // Outside the domain
                {
                    return 0;
                }

                setValue_f(annularMask, ic[0]+delta, -1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 1); // x+ face
                setValue_f(annularMask, ic[0], ic[0]+delta+delta%2, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5); // y- face
                setValue_f(annularMask, ic[0], ic[0]+delta+delta%2, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5); // y+ face
                setValue_f(annularMask, ic[0], ic[0]+delta+delta%2, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.5); // z- face
                setValue_f(annularMask, ic[0], ic[0]+delta+delta%2, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.5); // z+ face

                setValue_f(annularMask, ic[0], ic[0]+delta+delta%2, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // y-,z- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0], ic[0]+delta+delta%2, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // y-,z+ edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0], ic[0]+delta+delta%2, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // y+,z- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0], ic[0]+delta+delta%2, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // y+,z+ edge, 0.75 if delta is even, 0.25 if delta is odd

                ss[0] = (int)(ic[0]);
                ss[1] = (int)(ic[0]+delta+1);
                ss[2] = (int)(ic[1]-(int)(ceil(delta/2.)));
                ss[3] = (int)(ic[1]+(int)(ceil(delta/2.))+1);
                ss[4] = (int)(ic[2]-(int)(ceil(delta/2.)));
                ss[5] = (int)(ic[2]+(int)(ceil(delta/2.))+1);

                break;
            
            case 1: // side x+ -> x-
                if (ic[0]-delta < 0 || ic[1]-ceil(delta/2.) < 0 ||  ic[1]+ceil(delta/2.)>=n_points[1] || ic[2]-ceil(delta/2.) < 0 ||  ic[2]+ceil(delta/2.)>=n_points[2]) // Outside the domain
                {
                    return 0;
                }

                setValue_f(annularMask, ic[0]-delta, -1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 1); // x- face
                setValue_f(annularMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5); // y- face
                setValue_f(annularMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5); // y+ face
                setValue_f(annularMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.5); // z- face
                setValue_f(annularMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.5); // z+ face

                setValue_f(annularMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // y-,z- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // y-,z+ edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // y+,z- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]-delta+abs((delta-1)%2), ic[0]+1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // y+,z+ edge, 0.75 if delta is even, 0.25 if delta is odd
                
                ss[0] = (int)(ic[0]-delta);
                ss[1] = (int)(ic[0]+1);
                ss[2] = (int)(ic[1]-(int)(ceil(delta/2.)));
                ss[3] = (int)(ic[1]+(int)(ceil(delta/2.))+1);
                ss[4] = (int)(ic[2]-(int)(ceil(delta/2.)));
                ss[5] = (int)(ic[2]+(int)(ceil(delta/2.))+1);

                break;

            case 2: // side y- -> y+
                if (ic[1]+delta >= n_points[1] || ic[0]-ceil(delta/2.) < 0 ||  ic[0]+ceil(delta/2.)>=n_points[0] || ic[2]-ceil(delta/2.) < 0 ||  ic[2]+ceil(delta/2.)>=n_points[2]) // Outside the domain
                {
                    return 0;
                }

                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]+delta, -1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 1); // y- face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5); // x- face
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5); // x+ face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1], ic[1]+delta+delta%2, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.5); // z- face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1], ic[1]+delta+delta%2, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.5); // z+ face

                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x-,z- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x-,z+ edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x+,z- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1], ic[1]+delta+delta%2, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x+,z+ edge, 0.75 if delta is even, 0.25 if delta is odd

                ss[0] = (int)(ic[0]-(int)(ceil(delta/2.)));
                ss[1] = (int)(ic[0]+(int)(ceil(delta/2.))+1);
                ss[2] = (int)(ic[1]);
                ss[3] = (int)(ic[1]+delta+1);
                ss[4] = (int)(ic[2]-(int)(ceil(delta/2.)));
                ss[5] = (int)(ic[2]+(int)(ceil(delta/2.))+1);

                break;

            case 3: // side y+ -> y-
                if (ic[1]-delta < 0 || ic[0]-ceil(delta/2.) < 0 ||  ic[0]+ceil(delta/2.)>=n_points[0] || ic[2]-ceil(delta/2.) < 0 ||  ic[2]+ceil(delta/2.)>=n_points[2]) // Outside the domain
                {
                    return 0;
                }

                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]-delta, -1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 1); // y- face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5); // x- face
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]-(int)(ceil(delta/2.)), ic[2]+(int)(ceil(delta/2.))+1, n_points, 0.5); // x+ face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.5); // z- face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.5); // z+ face

                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x-,z- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x-,z+ edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]-(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x+,z- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-delta+abs((delta-1)%2), ic[1]+1, ic[2]+(int)(ceil(delta/2.)), -1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x+,z+ edge, 0.75 if delta is even, 0.25 if delta is odd

                ss[0] = (int)(ic[0]-(int)(ceil(delta/2.)));
                ss[1] = (int)(ic[0]+(int)(ceil(delta/2.))+1);
                ss[2] = (int)(ic[1]-delta);
                ss[3] = (int)(ic[1]+1);
                ss[4] = (int)(ic[2]-(int)(ceil(delta/2.)));
                ss[5] = (int)(ic[2]+(int)(ceil(delta/2.))+1);

                break;
            
            case 4: // side z- -> z+
                if (ic[2]+delta >= n_points[2] || ic[0]-ceil(delta/2.) < 0 ||  ic[0]+ceil(delta/2.)>=n_points[0] || ic[1]-ceil(delta/2.) < 0 ||  ic[1]+ceil(delta/2.)>=n_points[1]) // Outside the domain
                {
                    return 0;
                }

                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]+delta, -1, n_points, 1); // z+ face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2], ic[2]+delta+delta%2, n_points, 0.5); // # x- face
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2], ic[2]+delta+delta%2, n_points, 0.5); // # x+ face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.5); // # y- face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.5); // # y+ face

                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x-,y- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x-,y+ edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x+,y- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2], ic[2]+delta+delta%2, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x+,y+ edge, 0.75 if delta is even, 0.25 if delta is odd

                ss[0] = (int)(ic[0]-(int)(ceil(delta/2.)));
                ss[1] = (int)(ic[0]+(int)(ceil(delta/2.))+1);
                ss[2] = (int)(ic[1]-(int)(ceil(delta/2.)));
                ss[3] = (int)(ic[1]+(int)(ceil(delta/2.))+1);
                ss[4] = (int)(ic[2]);
                ss[5] = (int)(ic[2]+delta+1);

                break;

            case 5: // side z+ -> z-
                if (ic[2]-delta < 0 || ic[0]-ceil(delta/2.) < 0 ||  ic[0]+ceil(delta/2.)>=n_points[0] || ic[1]-ceil(delta/2.) < 0 ||  ic[1]+ceil(delta/2.)>=n_points[1]) // Outside the domain
                {
                    return 0;
                }

                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]-delta, -1, n_points, 1); // z- face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.5); // # x- face
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), ic[1]+(int)(ceil(delta/2.))+1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.5); // # x+ face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.5); // # y- face
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), ic[0]+(int)(ceil(delta/2.))+1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.5); // # y+ face

                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x-,y- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]-(int)(ceil(delta/2.)), -1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x-,y+ edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]-(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x+,y- edge, 0.75 if delta is even, 0.25 if delta is odd
                setValue_f(annularMask, ic[0]+(int)(ceil(delta/2.)), -1, ic[1]+(int)(ceil(delta/2.)), -1, ic[2]-delta+abs((delta-1)%2), ic[2]+1, n_points, 0.25*(delta%2) + 0.75*(abs((delta-1)%2))); // x+,y+ edge, 0.75 if delta is even, 0.25 if delta is odd

                ss[0] = (int)(ic[0]-(int)(ceil(delta/2.)));
                ss[1] = (int)(ic[0]+(int)(ceil(delta/2.))+1);
                ss[2] = (int)(ic[1]-(int)(ceil(delta/2.)));
                ss[3] = (int)(ic[1]+(int)(ceil(delta/2.))+1);
                ss[4] = (int)(ic[2]-delta);
                ss[5] = (int)(ic[2]+1);

                break;
        }
        cumMass += nansum_product_ff(massArray, annularMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points);

        if ((cumMass-targetMass)/targetMass < -(AVG_VOL_TOLL)) // cumMass is lower than target mass and outside the accepted tollerance
        {
            incrementByOther(massUsageMask, annularMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 1);
        }
        else if ((cumMass-targetMass)/targetMass > AVG_VOL_TOLL) // cumMass is higher than target mass and outside the accepted tollerance)
        {
            if (delta == 0)
                massUsageMask[from3DIdxTo1DIdx(ic, n_points)] = targetMass / massArray[from3DIdxTo1DIdx(ic, n_points)];
            else
            {
                getAnnularMasses_SC(massArray, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, masses, delta, face);
                k = solveForK(masses[0], masses[1], masses[2], prevCumMass-targetMass);
                #ifdef WARNINGS
                    if (k==0 || k >= 1)
                        printf("WARNING on voxel %d, %d, %d\n", ic[0], ic[1], ic[2]);
                #endif
                setPartialMassUsageMask_SC(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, k, delta, face);
            }
            stopFlag = 1;
        }
        else // cumMass is within the target mass accepted tollerance
        {
            incrementByOther(massUsageMask, annularMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 1);
            stopFlag = 1;
        }
        prevCumMass = cumMass;
        delta++;
        setValue_f(annularMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 0); // I reset the annularMask
        #ifdef DEBUG
            if (stopFlag == 1)
            {
            printf("side: %d\n", face);
            printf("delta: %d\n", delta);
            printf("prevCumMass: %f\n", prevCumMass);
            printf("CumMass: %f\n", cumMass);
            printf("targetMass: %f\n", targetMass);
            printf("Relative Ratio: %e\n", (cumMass-targetMass)/targetMass);
            printf("ic: %d,%d,%d\n",ic[0],ic[1],ic[2]);
            printf("ma: %f, mb: %f, mc: %f, d: %f\n",masses[0], masses[1], masses[2], prevCumMass-targetMass);
            printf("k: %.10f\n",k);
            printf("slices: %d, %d, %d, %d, %d, %d\n", ss[0], ss[1], ss[2], ss[3], ss[4], ss[5]);
            printf("\n\n");
            }
        #endif
    }
    return 1;
}

int main(double* massArray, double* localSARArray, double* avgSARArray, Status* voxStatusArray, double targetMass, int* n_points)
{
    int tot_points = n_points[0] * n_points[1] * n_points[2];
    int idx = 0;
    int ss[2*DIMS] = {0, 0, 0, 0, 0, 0}; // slice of array as [i_min,i_max,j_min,j_max,k_min,k_max];
    int ic[3] = {0,0,0};
    int delta_min, retValue, face=0;
    double* massUsageMask = calloc(tot_points, sizeof(double));
    double* annularMask = calloc(tot_points, sizeof(double));
    double volumes[6] = {0, 0, 0, 0, 0, 0};
    double avgSARs[6] = {0, 0, 0, 0, 0, 0};
    double minVolume = 0;
    #ifdef REPORT
        FILE* report_ptr;
        report_ptr = fopen(REPORT, "a"); 
        double avgMasses[6] = {0, 0, 0, 0, 0, 0};
        double avgMass = 0, avgVolume = 0;
        int  selected_face = 0;
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

                if (massArray[idx] == massArray[idx] && voxStatusArray[idx] == UNUSED) // massArray[i] is not nan and the voxel is UNUSED
                {
                    minVolume = 0;
                    for (face = 0; face < 2*DIMS; face++)
                    {
                        retValue = averagingMaskSide3D(ic, n_points, massArray, delta_min, targetMass, face, massUsageMask, annularMask, ss);

                        if (retValue == 0)
                        {
                            volumes[face] = 0;
                            avgSARs[face] = 0;
                            #ifdef REPORT
                                avgMasses[face] = 0;
                            #endif
                        }
                        else
                        {
                            avgSARs[face] = targetMassAverage(massArray, massUsageMask, localSARArray, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, targetMass);
                            volumes[face] = nansum_f(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points);

                            #ifdef REPORT
                                avgMasses[face] = nansum_product_ff(massArray, massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points); // kg
                            #endif

                            if ((minVolume == 0 && volumes[face]>0)||(minVolume > volumes[face])){minVolume = volumes[face];}
                        }
                        setValue_f(massUsageMask, ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], n_points, 0); // I set to zero the values I already know are not zero
                    }
                    for (face = 0; face < 2*DIMS; face++)
                    {
                        if (volumes[face] <= minVolume*MAX_VOL_RATIO && avgSARs[face] > avgSARArray[idx])
                        {
                            avgSARArray[idx] = avgSARs[face];

                            #ifdef REPORT
                                avgMass = avgMasses[face];
                                avgVolume = volumes[face];
                                selected_face = face;
                            #endif
                        }
                    }
                    #ifdef REPORT
                        // UNUSED voxels
                        fprintf(report_ptr, "%d %d %d %d %d %.6e %.6e %d %.6e %.6e\n", idx, i, j, k, voxStatusArray[idx], avgMass*1e3, avgVolume*VOX_VOL, selected_face+1, localSARArray[idx], avgSARArray[idx]);
                    #endif
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
        fclose(report_ptr);
    #endif

    return 0;
}