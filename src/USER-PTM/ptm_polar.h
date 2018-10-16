#ifndef PTM_POLAR_H
#define PTM_POLAR_H

#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>

namespace ptm {

int polar_decomposition_3x3(double* _A, bool right_sided, double* U, double* P);
void InnerProduct(double *A, int num, const double (*coords1)[3], double (*coords2)[3], int8_t* permutation);
int FastCalcRMSDAndRotation(double *A, double E0, double *p_nrmsdsq, double *q, double* U);

}

#endif

