#ifndef ZHEEVD_H
#define ZHEEVD_H

#ifdef __cplusplus
extern "C" {
#endif
  typedef struct { double r, i; } doublecomplex;

  /* Subroutine */ int zheevd_(char *jobz, char *uplo, int *n,
                               doublecomplex *a, int *lda, double *w, doublecomplex *work,
                               int *lwork, double *rwork, int *lrwork, int *iwork,
                               int *liwork, int *info);
#ifdef __cplusplus
}
#endif
#endif
