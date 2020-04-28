//
// Created by Yury Lysogorskiy on 31.01.20.
//

#include "ace_evaluator.h"
#include "ace_types.h"
#include "ace_abstract_basis.h"

void ACEEvaluator::init(ACEAbstractBasisSet *basis_set) {
    A.init(basis_set->nelements, basis_set->nradmax + 1, basis_set->lmax + 1, "A");
    A_rank1.init(basis_set->nelements, basis_set->nradbase, "A_rank1");

    rhos.init(basis_set->ndensitymax + 1, "rhos"); // +1 density for core repulsion
    dF_drho.init(basis_set->ndensitymax + 1, "dF_drho"); // +1 density for core repulsion
}

void ACEEvaluator::init_timers() {
    loop_over_neighbour_timer.init();

    forces_calc_loop_timer.init();
    forces_calc_neighbour_timer.init();
    phi_calc_timer.init();
    phi_recalc_timer.init();
    energy_calc_timer.init();
    bond_calc_timer.init();
    A_calc_timer.init();

    per_atom_calc_timer.init();


    basis_func_calc_timer.init();
    total_time_calc_timer.init();

}

//================================================================================================================

void ACECTildeEvaluator::set_basis(ACECTildeBasisSet &bas) {
    basis_set = &bas;
    init(basis_set);
}

void ACECTildeEvaluator::init(ACECTildeBasisSet *basis_set) {

    ACEEvaluator::init(basis_set);


    weights.init(basis_set->nelements, basis_set->nradmax + 1, basis_set->lmax + 1,
                 "weights");

    weights_rank1.init(basis_set->nelements, basis_set->nradbase, "weights_rank1");


    DG_cache.init(1, basis_set->nradbase, "DG_cache");
    DG_cache.fill(0);

    R_cache.init(1, basis_set->nradmax, basis_set->lmax + 1, "R_cache");
    R_cache.fill(0);

    DR_cache.init(1, basis_set->nradmax, basis_set->lmax + 1, "DR_cache");
    DR_cache.fill(0);

    Y_cache.init(1, basis_set->lmax + 1, "Y_cache");
    Y_cache.fill({0, 0});

    DY_cache.init(1, basis_set->lmax + 1, "dY_dense_cache");
    DY_cache.fill({0.});

    //hard-core repulsion
    DCR_cache.init(1, "DCR_cache");
    DCR_cache.fill(0);
    dB_flatten.init(basis_set->max_dB_array_size, "dB_flatten");


}

void ACECTildeEvaluator::resize_neighbours_cache(int max_jnum) {
    if (R_cache.get_dim(0) < max_jnum) {

        //TODO: implement grow
        R_cache.resize(max_jnum, basis_set->nradmax, basis_set->lmax + 1);
        R_cache.fill(0);

        DR_cache.resize(max_jnum, basis_set->nradmax, basis_set->lmax + 1);
        DR_cache.fill(0);

        DG_cache.resize(max_jnum, basis_set->nradbase);
        DG_cache.fill(0);

        Y_cache.resize(max_jnum, basis_set->lmax + 1);
        Y_cache.fill({0, 0});

        DY_cache.resize(max_jnum, basis_set->lmax + 1);
        DY_cache.fill({0});

        //hard-core repulsion
        DCR_cache.init(max_jnum, "DCR_cache");
        DCR_cache.fill(0);
    }
}



// double** r - atomic coordinates of atom I
// int* types - atomic types if atom I
// int **firstneigh -  ptr to 1st J int value of each I atom. Usage: jlist = firstneigh[i];
// Usage: j = jlist_of_i[jj];
// jnum - number of J neighbors for each I atom.  jnum = numneigh[i];

void
ACECTildeEvaluator::compute_atom(int i, DOUBLE_TYPE **x, const SPECIES_TYPE *type, const int jnum, const int *jlist) {
    per_atom_calc_timer.start();
#ifdef PRINT_MAIN_STEPS
    printf("\n ATOM: ind = %d r_norm=(%f, %f, %f)\n",i, x[i][0], x[i][1], x[i][2]);
#endif
    DOUBLE_TYPE evdwl = 0, evdwl_cut = 0, rho_core = 0;
    DOUBLE_TYPE r_norm;
    DOUBLE_TYPE xn, yn, zn, r_xyz;
    DOUBLE_TYPE R, GR, DGR, R_over_r, DR, DCR;
    DOUBLE_TYPE *r_hat;

    SPECIES_TYPE mu_j;
    DENSITY_TYPE ndensity; //TODO: extract from basis set, as it is equal to all functions
    RANK_TYPE r, rank, t;
    NS_TYPE n;
    LS_TYPE l;
    MS_TYPE m, m_t;

    SPECIES_TYPE *mus;
    NS_TYPE *ns;
    LS_TYPE *ls;
    MS_TYPE *ms;

    int j, jj, func_ind, ms_ind;
    SHORT_INT_TYPE factor;

    ACEComplex Y{0}, Y_DR{0.};
    ACEComplex B{0.};
    ACEComplex dB{0};
    ACEComplex A_cache[basis_set->rankmax];

    dB_flatten.fill({0.});

    Dycomponent grad_phi_nlm{0}, DY{0.};

    //size is +1 of max to avoid out-of-boundary array access in double-triangular scheme
    ACEComplex A_forward_prod[basis_set->rankmax + 1];
    ACEComplex A_backward_prod[basis_set->rankmax + 1];

    DOUBLE_TYPE inv_r_norm;
    DOUBLE_TYPE r_norms[jnum];
    DOUBLE_TYPE inv_r_norms[jnum];
    DOUBLE_TYPE rhats[jnum][3];//normalized vector
    SPECIES_TYPE elements[jnum];
    const DOUBLE_TYPE xtmp = x[i][0];
    const DOUBLE_TYPE ytmp = x[i][1];
    const DOUBLE_TYPE ztmp = x[i][2];
    DOUBLE_TYPE f_ji[3];

    bool is_element_mapping = element_type_mapping.get_size() > 0;
    SPECIES_TYPE mu_i;
    if (is_element_mapping)
        mu_i = element_type_mapping(type[i]);
    else
        mu_i = type[i];

    const SHORT_INT_TYPE total_basis_size_rank1 = basis_set->total_basis_size_rank1[mu_i];
    const SHORT_INT_TYPE total_basis_size = basis_set->total_basis_size[mu_i];

    ACECTildeBasisFunction *basis_rank1 = basis_set->basis_rank1[mu_i];
    ACECTildeBasisFunction *basis = basis_set->basis[mu_i];

    DOUBLE_TYPE rho_cut, drho_cut, fcut, dfcut;
    DOUBLE_TYPE dF_drho_core;

    //TODO: lmax -> lmaxi
    const LS_TYPE lmaxi = basis_set->lmax;

    //TODO: nradmax -> nradiali
    const NS_TYPE nradiali = basis_set->nradmax;

    //TODO: nradbase -> nradbasei
    const NS_TYPE nradbasei = basis_set->nradbase;

    neighbours_forces.resize(jnum, 3);
    neighbours_forces.fill(0);

    //TODO: shift nullifications to place where arrays are used
    weights.fill({0});
    weights_rank1.fill(0);
    A.fill({0});
    A_rank1.fill(0);
    rhos.fill(0);
    dF_drho.fill(0);

    //proxy references to spherical harmonics and radial functions arrays
    const Array2DLM<ACEComplex> &ylm = basis_set->spherical_harmonics.ylm;
    const Array2DLM<Dycomponent> &dylm = basis_set->spherical_harmonics.dylm;

    const Array2D<DOUBLE_TYPE> &fr = basis_set->radial_functions.fr;
    const Array2D<DOUBLE_TYPE> &dfr = basis_set->radial_functions.dfr;

    const Array1D<DOUBLE_TYPE> &gr = basis_set->radial_functions.gr;
    const Array1D<DOUBLE_TYPE> &dgr = basis_set->radial_functions.dgr;

    loop_over_neighbour_timer.start();

    int jj_actual = 0;
    SPECIES_TYPE type_j = 0;
    int neighbour_index_mapping[jnum]; // jj_actual -> jj
    //loop over neighbours, compute distance, consider only atoms within with r<cutoff(mu_i, mu_j)
    for (jj = 0; jj < jnum; ++jj) {

        j = jlist[jj];
        xn = x[j][0] - xtmp;
        yn = x[j][1] - ytmp;
        zn = x[j][2] - ztmp;
        type_j = type[j];
        if (is_element_mapping)
            mu_j = element_type_mapping(type_j);
        else
            mu_j = type_j;

        DOUBLE_TYPE current_cutoff = basis_set->radial_functions.cut(mu_i, mu_j);
        r_xyz = sqrt(xn * xn + yn * yn + zn * zn);

        if (r_xyz >= current_cutoff)
            continue;

        inv_r_norm = 1 / r_xyz;

        r_norms[jj_actual] = r_xyz;
        inv_r_norms[jj_actual] = inv_r_norm;
        rhats[jj_actual][0] = xn * inv_r_norm;
        rhats[jj_actual][1] = yn * inv_r_norm;
        rhats[jj_actual][2] = zn * inv_r_norm;
        elements[jj_actual] = mu_j;
        neighbour_index_mapping[jj_actual] = jj;
        jj_actual++;
    }

    int jnum_actual = jj_actual;

    //ALGORITHM 1: Atomic base A
    for (jj = 0; jj < jnum_actual; ++jj) {
        r_norm = r_norms[jj];
        mu_j = elements[jj];
        r_hat = rhats[jj];

        //proxies
        Array2DLM<ACEComplex> &Y_jj = Y_cache(jj);
        Array2DLM<Dycomponent> &DY_jj = DY_cache(jj);


        basis_set->radial_functions.lookupRadspline(r_norm, basis_set->nradbase, nradiali, mu_i, mu_j);
        basis_set->spherical_harmonics.compute_ylm(r_hat[0], r_hat[1], r_hat[2], lmaxi);
        //loop for computing A's
        //rank = 1
        for (n = 0; n < basis_set->nradbase; n++) {
            GR = gr(n);
#ifdef DEBUG_ENERGY_CALCULATIONS
            printf("-neigh atom %d\n", jj);
            printf("gr(n=%d)(r=%f) = %f\n", n, r_norm, gr(n));
            printf("dgr(n=%d)(r=%f) = %f\n", n, r_norm, dgr(n));
#endif
            DG_cache(jj, n) = dgr(n);
            A_rank1(mu_j, n) += GR * Y00;
        }
        //loop for computing A's
        // for rank > 1
        for (n = 0; n < nradiali; n++) {
            auto &A_lm = A(mu_j, n);
            for (l = 0; l <= lmaxi; l++) {
                R = fr(n, l);
#ifdef DEBUG_ENERGY_CALCULATIONS
                printf("R(nl=%d,%d)(r=%f)=%f\n", n + 1, l, r_norm, R);
#endif

                DR_cache(jj, n, l) = dfr(n, l);
                R_cache(jj, n, l) = R;

                for (m = 0; m <= l; m++) {
                    Y = ylm(l, m);
#ifdef DEBUG_ENERGY_CALCULATIONS
                    printf("Y(lm=%d,%d)=(%f, %f)\n", l, m, Y.real, Y.img);
#endif
                    A_lm(l, m) += R * Y; //accumulation sum over neighbours
                    Y_jj(l, m) = Y;
                    DY_jj(l, m) = dylm(l, m);
                }
            }
        }

        //hard-core repulsion
        rho_core += basis_set->radial_functions.cr;
        DCR_cache(jj) = basis_set->radial_functions.dcr;

    } //end loop over neighbours

    //complex conjugate A's (for NEGATIVE (-m) terms)
    // for rank > 1
    for (mu_j = 0; mu_j < basis_set->nelements; mu_j++) {
        for (n = 0; n < nradiali; n++) {
            auto &A_lm = A(mu_j, n);
            for (l = 0; l <= lmaxi; l++) {
                //fill in -m part in the outer loop using the same m <-> -m symmetry as for Ylm
                for (m = 1; m <= l; m++) {
                    factor = m % 2 == 0 ? 1 : -1;
                    A_lm(l, -m) = A_lm(l, m).conjugated() * factor;
                }
            }
        }
    }    //now A's are constructed
    loop_over_neighbour_timer.stop();

    // ==================== ENERGY ====================

    energy_calc_timer.start();

    //ALGORITHM 2: Basis functions B with iterative product and density rho(p) calculation
    //rank=1
    for (int f_ind = 0; f_ind < total_basis_size_rank1; ++f_ind) {
        ACECTildeBasisFunction *func = &basis_rank1[f_ind];
        ndensity = func->ndensity;
#ifdef PRINT_LOOPS_INDICES
        printf("Num density = %d r = 0\n",(int) ndensity );
        print_C_tilde_B_basis_function(*func);
#endif
        double A_cur = A_rank1(func->mus[0], func->ns[0] - 1);
#ifdef DEBUG_ENERGY_CALCULATIONS
        printf("A_r=1(x=%d, n=%d)=(%f)\n", func->mus[0], func->ns[0], A_cur);
#endif
        for (DENSITY_TYPE p = 0; p < ndensity; ++p) {
            //for rank=1 (r=0) only 1 ms-combination exists (ms_ind=0), so index of func.ctildes is 0..ndensity-1
            rhos(p) += func->ctildes[p] * A_cur;
        }
    } // end loop for rank=1

    //rank>1
    int func_ms_ind = 0;
    int func_ms_t_ind = 0;// index for dB

    for (func_ind = 0; func_ind < total_basis_size; ++func_ind) {
        ACECTildeBasisFunction *func = &basis[func_ind];
        ndensity = func->ndensity;
        rank = func->rank;
        r = rank - 1;
#ifdef PRINT_LOOPS_INDICES
        printf("Num density = %d r = %d\n",(int) ndensity, (int)r );
        print_C_tilde_B_basis_function(*func);
#endif
        mus = func->mus;
        ns = func->ns;
        ls = func->ls;

        //loop over {ms} combinations in sum
        for (ms_ind = 0; ms_ind < func->num_ms_combs; ++ms_ind, ++func_ms_ind) {
            ms = &func->ms_combs[ms_ind * rank]; // current ms-combination (of length = rank)

            //loop over m, collect B  = product of A with given ms
            A_forward_prod[0] = 1;
            A_backward_prod[r] = 1;

            //fill forward A-product triangle
            for (t = 0; t < rank; t++) {
                //TODO: optimize ns[t]-1 -> ns[t] during functions construction
                A_cache[t] = A(mus[t], ns[t] - 1, ls[t], ms[t]);
#ifdef DEBUG_ENERGY_CALCULATIONS
                printf("A(x=%d, n=%d, l=%d, m=%d)=(%f,%f)\n", mus[t], ns[t], ls[t], ms[t], A_cache[t].real,
                       A_cache[t].img);
#endif
                A_forward_prod[t + 1] = A_forward_prod[t] * A_cache[t];
            }

            B = A_forward_prod[t];
#ifdef DEBUG_FORCES_CALCULATIONS
            printf("B = (%f, %f)\n", (B).real, (B).img);
#endif
            //fill backward A-product triangle
            for (t = r; t >= 1; t--) {
                A_backward_prod[t - 1] =
                        A_backward_prod[t] * A_cache[t];
            }

            for (t = 0; t < rank; ++t, ++func_ms_t_ind) {
                dB = A_forward_prod[t] * A_backward_prod[t]; //dB - product of all A's except t-th
                dB_flatten(func_ms_t_ind) = dB;
#ifdef DEBUG_FORCES_CALCULATIONS
                m_t = ms[t];
                printf("dB(n,l,m)(%d,%d,%d) = (%f, %f)\n", ns[t], ls[t], m_t, (dB).real, (dB).img);
#endif
            }

            for (DENSITY_TYPE p = 0; p < ndensity; ++p) {
                //real-part only multiplication
                rhos(p) += B.real_part_product(func->ctildes[ms_ind * ndensity + p]);
#ifdef PRINT_INTERMEDIATE_VALUES
                printf("rhos(%d) += %f\n", p, B.real_part_product(func->ctildes[ms_ind * ndensity + p]));
                printf("Rho[i = %d][p = %d] = %f\n",  i , p , rhos(p));
#endif
            }
        }//end of loop over {ms} combinations in sum
    }// end loop for rank>1

#ifdef DEBUG_FORCES_CALCULATIONS
    printf("rhos = ");
    for(DENSITY_TYPE p =0; p<ndensity; ++p) printf(" %.20f ",rhos(p));
    printf("\n");
#endif


    // energy cutoff
    rho_cut = basis_set->rho_core_cutoffs(mu_i);
    drho_cut = basis_set->drho_core_cutoffs(mu_i);

    basis_set->inner_cutoff(rho_core, rho_cut, drho_cut, fcut, dfcut);
    basis_set->FS_values_and_derivatives(rhos, evdwl, dF_drho, ndensity);

    dF_drho_core = evdwl * dfcut + 1;
    for (DENSITY_TYPE p = 0; p < ndensity; ++p)
        dF_drho(p) *= fcut;
    evdwl_cut = evdwl * fcut + rho_core;

#ifdef DEBUG_FORCES_CALCULATIONS
    printf("dFrhos = ");
    for(DENSITY_TYPE p =0; p<ndensity; ++p) printf(" %f ",dF_drho(p));
    printf("\n");
#endif

    //ALGORITHM 3: Weights and theta calculation
    // rank = 1
    for (int f_ind = 0; f_ind < total_basis_size_rank1; ++f_ind) {
        ACECTildeBasisFunction *func = &basis_rank1[f_ind];
        ndensity = func->ndensity;
        for (DENSITY_TYPE p = 0; p < ndensity; ++p) {
            //for rank=1 (r=0) only 1 ms-combination exists (ms_ind=0), so index of func.ctildes is 0..ndensity-1
            weights_rank1(func->mus[0], func->ns[0] - 1) += dF_drho(p) * func->ctildes[p];
        }
    }

    // rank>1
    func_ms_ind = 0;
    func_ms_t_ind = 0;// index for dB
    DOUBLE_TYPE theta = 0;
    for (func_ind = 0; func_ind < total_basis_size; ++func_ind) {
        ACECTildeBasisFunction *func = &basis[func_ind];
        ndensity = func->ndensity;
        rank = func->rank;
        mus = func->mus;
        ns = func->ns;
        ls = func->ls;
        for (ms_ind = 0; ms_ind < func->num_ms_combs; ++ms_ind, ++func_ms_ind) {
            ms = &func->ms_combs[ms_ind * rank];
            theta = 0;
            for (DENSITY_TYPE p = 0; p < ndensity; ++p) {
                theta += dF_drho(p) * func->ctildes[ms_ind * ndensity + p];
#ifdef DEBUG_FORCES_CALCULATIONS
                printf("(p=%d) theta += dF_drho[p] * func.ctildes[ms_ind * ndensity + p] = %f * %f = %f\n",p, dF_drho(p), func->ctildes[ms_ind * ndensity + p],dF_drho(p)*func->ctildes[ms_ind * ndensity + p]);
                printf("theta=%f\n",theta);
#endif
            }

            theta *= 0.5; // 0.5 factor due to possible double counting ???
            for (t = 0; t < rank; ++t, ++func_ms_t_ind) {
                m_t = ms[t];
                factor = (m_t % 2 == 0 ? 1 : -1);
                dB = dB_flatten(func_ms_t_ind);
                weights(mus[t], ns[t] - 1, ls[t], m_t) += theta * dB; //Theta_array(func_ms_ind);
                // update -m_t (that could also be positive), because the basis is half_basis
                weights(mus[t], ns[t] - 1, ls[t], -m_t) +=
                        theta * (dB).conjugated() * factor;// Theta_array(func_ms_ind);
#ifdef DEBUG_FORCES_CALCULATIONS
                printf("dB(n,l,m)(%d,%d,%d) = (%f, %f)\n", ns[t], ls[t], m_t, (dB).real, (dB).img);
                printf("theta = %f\n",theta);
                printf("weights(n,l,m)(%d,%d,%d) += (%f, %f)\n", ns[t], ls[t], m_t, (theta * dB * 0.5).real,
                       (theta * dB * 0.5).img);
                printf("weights(n,l,-m)(%d,%d,%d) += (%f, %f)\n", ns[t], ls[t], -m_t,
                       ( theta * (dB).conjugated() * factor * 0.5).real,
                       ( theta * (dB).conjugated() * factor * 0.5).img);
#endif
            }
        }
    }
    energy_calc_timer.stop();

// ==================== FORCES ====================
#ifdef PRINT_MAIN_STEPS
    printf("\nFORCE CALCULATION\n");
    printf("loop over neighbours\n");
#endif

    forces_calc_loop_timer.start();
// loop over neighbour atoms for force calculations
    for (jj = 0; jj < jnum_actual; ++jj) {
        mu_j = elements[jj];
        r_hat = rhats[jj];
        inv_r_norm = inv_r_norms[jj];

        Array2DLM<ACEComplex> &Y_cache_jj = Y_cache(jj);
        Array2DLM<Dycomponent> &DY_cache_jj = DY_cache(jj);

#ifdef PRINT_LOOPS_INDICES
        printf("\nneighbour atom #%d\n", jj);
        printf("rhat = (%f, %f, %f)\n", r_hat[0], r_hat[1], r_hat[2]);
#endif

        forces_calc_neighbour_timer.start();

        f_ji[0] = f_ji[1] = f_ji[2] = 0;

//for rank = 1
        for (n = 0; n < nradbasei; ++n) {
            if (weights_rank1(mu_j, n) == 0)
                continue;
            auto &DG = DG_cache(jj, n);
            DGR = DG * Y00;
            DGR *= weights_rank1(mu_j, n);
#ifdef DEBUG_FORCES_CALCULATIONS
            printf("r=1: (n,l,m)=(%d, 0, 0)\n",n+1);
            printf("\tGR(n=%d, r=%f)=%f\n",n+1,r_norm, gr(n));
            printf("\tDGR(n=%d, r=%f)=%f\n",n+1,r_norm, dgr(n));
            printf("\tdF+=(%f, %f, %f)\n",DGR * r_hat[0], DGR * r_hat[1], DGR * r_hat[2]);
#endif
            f_ji[0] += DGR * r_hat[0];
            f_ji[1] += DGR * r_hat[1];
            f_ji[2] += DGR * r_hat[2];
        }

//for rank > 1
        for (n = 0; n < nradiali; n++) {
            for (l = 0; l <= lmaxi; l++) {
                R_over_r = R_cache(jj, n, l) * inv_r_norm;
                DR = DR_cache(jj, n, l);

                // for m>=0
                for (m = 0; m <= l; m++) {
                    ACEComplex w = weights(mu_j, n, l, m);
                    if (w == 0)
                        continue;
                    //counting for -m cases if m>0
                    if (m > 0) w *= 2;
                    DY = DY_cache_jj(l, m);
                    Y_DR = Y_cache_jj(l, m) * DR;

                    grad_phi_nlm.a[0] = Y_DR * r_hat[0] + DY.a[0] * R_over_r;
                    grad_phi_nlm.a[1] = Y_DR * r_hat[1] + DY.a[1] * R_over_r;
                    grad_phi_nlm.a[2] = Y_DR * r_hat[2] + DY.a[2] * R_over_r;
#ifdef DEBUG_FORCES_CALCULATIONS
                    printf("d_phi(n=%d, l=%d, m=%d) = ((%f,%f), (%f,%f), (%f,%f))\n",n+1,l,m,
                           grad_phi_nlm.a[0].real, grad_phi_nlm.a[0].img,
                           grad_phi_nlm.a[1].real, grad_phi_nlm.a[1].img,
                           grad_phi_nlm.a[2].real, grad_phi_nlm.a[2].img);

                    printf("weights(n,l,m)(%d,%d,%d) = (%f,%f)\n", n+1, l, m,w.real, w.img);
                    //if (m>0) w*=2;
                    printf("dF(n,l,m)(%d, %d, %d) += (%f, %f, %f)\n", n + 1, l, m,
                           w.real_part_product(grad_phi_nlm.a[0]),
                           w.real_part_product(grad_phi_nlm.a[1]),
                           w.real_part_product(grad_phi_nlm.a[2])
                    );
#endif
// real-part multiplication only
                    f_ji[0] += w.real_part_product(grad_phi_nlm.a[0]);
                    f_ji[1] += w.real_part_product(grad_phi_nlm.a[1]);
                    f_ji[2] += w.real_part_product(grad_phi_nlm.a[2]);
                }
            }
        }


#ifdef PRINT_INTERMEDIATE_VALUES
        printf("f_ji(jj=%d, i=%d)=(%f, %f, %f)\n", jj, i,
               f_ji[0], f_ji[1], f_ji[2]
        );
#endif

        //hard-core repulsion
        DCR = DCR_cache(jj);
#ifdef   DEBUG_FORCES_CALCULATIONS
        printf("DCR = %f\n",DCR);
#endif
        f_ji[0] += dF_drho_core * DCR * r_hat[0];
        f_ji[1] += dF_drho_core * DCR * r_hat[1];
        f_ji[2] += dF_drho_core * DCR * r_hat[2];
#ifdef PRINT_INTERMEDIATE_VALUES
        printf("with core-repulsion\n");
        printf("f_ji(jj=%d, i=%d)=(%f, %f, %f)\n", jj, i,
               f_ji[0], f_ji[1], f_ji[2]
        );
        printf("neighbour_index_mapping[jj=%d]=%d\n",jj,neighbour_index_mapping[jj]);
#endif

        neighbours_forces(neighbour_index_mapping[jj], 0) = f_ji[0];
        neighbours_forces(neighbour_index_mapping[jj], 1) = f_ji[1];
        neighbours_forces(neighbour_index_mapping[jj], 2) = f_ji[2];

        forces_calc_neighbour_timer.stop();
    }// end loop over neighbour atoms for forces

    forces_calc_loop_timer.stop();

    //now, energies and forces are ready
    //energies(i) = evdwl + rho_core;
    e_atom = evdwl_cut;

#ifdef PRINT_INTERMEDIATE_VALUES
    printf("energies(i) = FS(...rho_p_accum...) = %f\n", evdwl);
#endif
    per_atom_calc_timer.stop();
}