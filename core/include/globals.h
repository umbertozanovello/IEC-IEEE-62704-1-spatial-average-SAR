#ifndef GLOBALS_H
#define GLOBALS_H

# define CM_ROOTS
/* #undef DEBUG */
/* #undef REPORT */
/* #undef WARNINGS */

# include <stdint.h>
# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>


# define VOX_VOL 1 // mm^3 Used only for IEC report
# define TH_FOR_INC 0.999
# define AVG_VOL_TOLL 0.000001
# define MAX_VOL_RATIO 1.05
# define DIMS 3
# define MAX_BACKGROUND_RATIO 0.1
# define PI 3.141592653589793

typedef enum {
  INVALID,
  UNUSED,
  USED,
  VALID,
} Status; 

extern void dgeev_( char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info );
                
#endif
