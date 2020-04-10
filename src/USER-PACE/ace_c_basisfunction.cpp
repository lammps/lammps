//
// Created by Yury Lysogorskiy on 26.02.20.
//
#include <cstdio>

#include "ace_c_basisfunction.h"

void print_C_tilde_B_basis_function(const ACECTildeBasisFunction &func) {
    printf("C_tilde_B_basis_function: ndensity= %d, mu0 = %d ", func.ndensity, func.mu0);
    printf(" XS=(");
    for (RANK_TYPE r = 0; r < func.rank; r++)
        printf("%d ", int(func.mus[r]));
    printf("), ns=(");
    for (RANK_TYPE r = 0; r < func.rank; r++)
        printf("%d ", func.ns[r]);
    printf("), ls=(");
    for (RANK_TYPE r = 0; r < func.rank; r++)
        printf("%d ", func.ls[r]);

    printf("), %d m_s combinations: {\n", func.num_ms_combs);
    for (int i = 0; i < func.num_ms_combs; i++) {
        printf("\t< ");
        for (RANK_TYPE r = 0; r < func.rank; r++)
            printf("%d ", func.ms_combs[i * func.rank + r]);
        printf(" >: c_tilde=");
        for (DENSITY_TYPE p = 0; p < func.ndensity; ++p)
            printf(" %f ", func.ctildes[i * func.ndensity + p]);
        printf("\n");
    }
    printf("}\n");
}
