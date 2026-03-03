#include <R.h>
#include <Rinternals.h>

#include <R_ext/BLAS.h>
#include <R_ext/libextern.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

void F77_NAME(hc)(int *n, int *len, int *iopt, int *ia, int *ib,
                 double *crit, double *membr, int *nn,
                 double *disnn, double *diss);

void F77_NAME(hcass2)(int *n, int *ia, int *ib, int *iorder, int *iia, int *iib);

/* =============================================================================================
 *
 *  Register native routines here.
 *
 * =============================================================================================*/
  
void attribute_visible R_init_flashClust(DllInfo * info)
{ 
  static R_NativePrimitiveArgType 
     hc_t[] = { INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,
                REALSXP, REALSXP, // CRIT and MEMBR
                INTSXP, REALSXP, // NN, DISNN
                REALSXP }; // DISS
  
  static R_NativePrimitiveArgType 
     hcass2_t[] = { INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP };
      
  static const R_FortranMethodDef FMethods[]  = {
    {"hc", (DL_FUNC) &F77_NAME(hc), 10, hc_t},
    {"hcass2", (DL_FUNC) &F77_NAME(hcass2), 6, hcass2_t},
    {NULL, NULL, 0}
  };

  R_registerRoutines(info, NULL, NULL, FMethods, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
} 

