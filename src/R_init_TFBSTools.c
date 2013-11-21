#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {
  /* matrixAlignerDynamic.c */
  CALLMETHOD_DEF(matrixAligner, 4),

  {NULL, NULL, 0}
};

void R_init_TFBSTools(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  return;
}
