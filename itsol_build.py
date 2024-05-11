from cffi import FFI
import os
ffibuilder = FFI()
ffibuilder.cdef("""typedef struct ITS_SparMat_
{
    int n;
    int *nzcount;  /* length of each row */
    int **ja;      /* pointer-to-pointer to store column indices  */
    double **ma;   /* pointer-to-pointer to store nonzero entries */

} ITS_SparMat;
typedef struct ITS_ILUSpar_
{
    int n;
    ITS_SparMat *L;   /* L part elements                            */
    double *D;        /* diagonal elements                          */
    ITS_SparMat *U;   /* U part elements                            */
    int *work;        /* working buffer */

} ITS_ILUSpar;
typedef struct ITS_CSRnumpy_
{
    int n;
    int *indices;
    int *indptr; 
    double *data;

} ITS_CSRnumpy;
int itsol_pc_lofC(int lofM, ITS_SparMat *csmat, ITS_ILUSpar *lu, FILE *fp); 
int itsol_pc_ilukC(int lofM, ITS_SparMat *csmat, ITS_ILUSpar *lu, FILE *fp);
int itsol_setupILU(ITS_ILUSpar *lu, int n);
int itsol_cleanILU(ITS_ILUSpar *lu);
ITS_ILUSpar *itsol_from_numpy(int n, int lfil, double *data, int *indices, int *indptr);
ITS_CSRnumpy * itsol_SpaFmtNumpy(ITS_SparMat *a);
int itsol_lusolC(double *y, double *x, ITS_ILUSpar *lu);
""")
ffibuilder.set_source("_itsol_cffi",
"""
     #include "include/pc-iluk.h"   // the C header of the library
     #include "include/data-types.h"
     #include "include/utils.h"
     #include "include/mat-utils.h"
""",
     libraries=['itsol','gfortran'], library_dirs=['./src/'])   # library name, for the linker
if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

