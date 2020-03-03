//
// Created by lysogy36 on 31.01.20.
//

#include "ace.h"

#include "ace_types.h"
#include "ace_atoms.h"

void ACE::set_basis(ACEBasisSet &bas) {
    basis_set = &bas;
    init();
}

void ACE::init() {
    A.init(basis_set->nelements, basis_set->nradmax + 1, basis_set->lmax + 1, "A");
    A_rank1.init(basis_set->nelements, basis_set->nradbase, "A_rank1");

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


    rhos.init(basis_set->ndensitymax, "rhos");
    dF_drho.init(basis_set->ndensitymax, "dF_drho");

    //Theta_array.init(basis_set->max_B_array_size, "Theta_array");
    dB_flatten.init(basis_set->max_dB_array_size, "dB_flatten");

}

void ACE::compute(ACEAtomicEnvironment &atomic_environment, bool verbose) {
    ACE_TIMER_INIT(loop_over_neighbour)
    ACE_TIMER_INIT(forces_calc_loop)
    ACE_TIMER_INIT(forces_calc_neighbour)
    ACE_TIMER_INIT(phi_calc)
    ACE_TIMER_INIT(phi_recalc)
    ACE_TIMER_INIT(energy_calc)
    ACE_TIMER_INIT(bond_calc)
    ACE_TIMER_INIT(A_calc)
    ACE_TIMER_INIT(per_atom_calc)
    ACE_TIMER_INIT(basis_func_calc)
    ACE_TIMER_INIT(total_time_calc)
    ACE_TIMER_START(total_time_calc)

    int i, j, jj;
    energy = 0;
    energies.resize(atomic_environment.n_atoms_tot);
    energies.fill(0);
    forces.resize(atomic_environment.n_atoms_tot, 3);// per-atom forces
    forces.fill(0);


    //loop over atoms
#ifdef PRINT_MAIN_STEPS
    cout << "=====LOOP OVER ATOMS=====" << endl;
#endif
    //derermine the maximum numer of neighbours
    int max_jnum = -1;
    for (i = 0; i < atomic_environment.n_atoms_tot; ++i)
        if (atomic_environment.num_neighbours[i] > max_jnum)
            max_jnum = atomic_environment.num_neighbours[i];

    resize_neighbours_cache(max_jnum);

    for (i = 0; i < atomic_environment.n_atoms_tot; ++i) {
        compute_atom(i,
                     atomic_environment.x,
                     atomic_environment.species_type,
                     atomic_environment.num_neighbours[i],
                     atomic_environment.neighbour_list[i]);
        //this will also update the energies(i) and neighbours_forces(jj, alpha) arrays
        //update global energies and forces accumulators
        energy += energies(i);
        for (jj = 0; jj < atomic_environment.num_neighbours[i]; jj++) {
            j = atomic_environment.neighbour_list[i][jj];

            forces(i, 0) += neighbours_forces(jj, 0);
            forces(i, 1) += neighbours_forces(jj, 1);
            forces(i, 2) += neighbours_forces(jj, 2);

            forces(j, 0) -= neighbours_forces(jj, 0);
            forces(j, 1) -= neighbours_forces(jj, 1);
            forces(j, 2) -= neighbours_forces(jj, 2);

#ifdef DEBUG_FORCES_CALCULATIONS
            printf("accumulated forces: F(i=%d)=(%f,%f,%f)\n", i, forces(i, 0), forces(i, 1), forces(i, 2));
            printf("accumulated forces: F(j=%d)=(%f,%f,%f)\n", j, forces(j, 0), forces(j, 1), forces(j, 2));
#endif
        }
    } // loop over atoms (i_at)
    ACE_TIMER_STOP(total_time_calc)

#ifdef FINE_TIMING
    if (verbose) {
        printf("   Total time: %ld microseconds\n", ACE_TIMER_MICROSECONDS(total_time_calc));
        printf("Per atom time:    %ld microseconds\n",
               ACE_TIMER_MICROSECONDS(per_atom_calc) / atomic_environment.n_atoms_tot);


        printf("Loop_over_nei/atom: %ld microseconds\n",
               ACE_TIMER_MICROSECONDS(loop_over_neighbour) / atomic_environment.n_atoms_tot);

        printf("       Energy/atom: %ld microseconds\n",
               ACE_TIMER_MICROSECONDS(energy_calc) / atomic_environment.n_atoms_tot);

        printf("       Forces/atom: %ld microseconds\n",
               ACE_TIMER_MICROSECONDS(forces_calc_loop) / atomic_environment.n_atoms_tot);

        printf("phi_recalcs/atom: %ld microseconds\n",
               ACE_TIMER_MICROSECONDS(phi_recalc) / atomic_environment.n_atoms_tot);

        printf("     forces_neig: %ld microseconds\n",
               ACE_TIMER_MICROSECONDS(forces_calc_neighbour) / atomic_environment.n_atoms_tot);

    }
#endif


}

void ACE::resize_neighbours_cache(int max_jnum) {
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
    }
}



// double** r - atomic coordinates of atom I
// int* types - atomic types if atom I
// int **firstneigh -  ptr to 1st J int value of each I atom. Usage: jlist = firstneigh[i];
// Usage: j = jlist_of_i[jj];
// jnum - number of J neighbors for each I atom.  jnum = numneigh[i];

void ACE::compute_atom(int i, DOUBLE_TYPE **x, const SPECIES_TYPE *type, int jnum, const int *jlist) {
    ACE_TIMER_START(per_atom_calc)
#ifdef PRINT_MAIN_STEPS
    cout << endl << "ATOM: ind = " << i << " r_norm=(" << x[i][0] << ","
         << x[i][1] << ","
         << x[i][2] << ")" << endl;
#endif
    DOUBLE_TYPE evdwl = 0;
    DOUBLE_TYPE r_norm;
    DOUBLE_TYPE xn, yn, zn, r_xyz;
    DOUBLE_TYPE R, GR, DGR, R_over_r, DR;
    DOUBLE_TYPE *r_hat;
    int j;
    SPECIES_TYPE elej;
    DENSITY_TYPE ndensity; //TODO: extract from basis set, as it is equal to all functions
    RANK_TYPE r, rank;
    NS_TYPE n;
    LS_TYPE l;
    MS_TYPE m, m_t;

    SPECIES_TYPE *mus;
    NS_TYPE *ns;
    LS_TYPE *ls;
    MS_TYPE *ms;

    int jj, func_ind, ms_ind;
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

    const int itype = type[i];

    const SPECIES_TYPE mu = map_lammps_at_type_to_element[itype];
    const SHORT_INT_TYPE total_basis_size_rank1 = basis_set->total_basis_size_rank1[mu];
    const SHORT_INT_TYPE total_basis_size = basis_set->total_basis_size[mu];

    ACECTildeBasisFunction *basis_rank1 = basis_set->basis_rank1[mu];
    ACECTildeBasisFunction *basis = basis_set->basis[mu];

    //TODO: lmax -> lmaxi
    const LS_TYPE lmaxi = basis_set->lmax;

    //TODO: nradmax -> nradiali
    const NS_TYPE nradiali = basis_set->nradmax;

    //TODO: nradbase -> nradbasei
    const NS_TYPE nradbasei = basis_set->nradbase;

    neighbours_forces.resize(jnum, 3);
    neighbours_forces.fill(0);

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

    ACE_TIMER_START(loop_over_neighbour)

    //loop over neighbours
    for (jj = 0; jj < jnum; ++jj) {

        j = jlist[jj];
        xn = x[j][0] - xtmp;
        yn = x[j][1] - ytmp;
        zn = x[j][2] - ztmp;

        r_xyz = sqrt(xn * xn + yn * yn + zn * zn);

        rhats[jj][0] = xn / r_xyz;
        rhats[jj][1] = yn / r_xyz;
        rhats[jj][2] = zn / r_xyz;

        r_norms[jj] = r_xyz;

        elements[jj] = map_lammps_at_type_to_element[type[j]];
    }


    //ALGORITHM 1: Atomic base construction
    for (jj = 0; jj < jnum; ++jj) {
        r_norm = r_norms[jj];
        elej = elements[jj];
        r_hat = rhats[jj];


        inv_r_norm = 1. / r_norm;
        inv_r_norms[jj] = inv_r_norm;
        //proxies
        Array2DLM<ACEComplex> &Y_jj = Y_cache(jj);
        Array2DLM<Dycomponent> &DY_jj = DY_cache(jj);


        basis_set->radial_functions.lookupRadspline(r_norm, basis_set->nradbase, nradiali, mu, elej);
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
            A_rank1(elej, n) += GR * Y00;
        }
        //loop for computing A's
        // for rank > 1
        for (n = 0; n < nradiali; n++) {
            auto &A_lm = A(elej, n);
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
                    A_lm(l, m) += Y * R; //accumulation sum over neighbours
                    Y_jj(l, m) = Y;
                    DY_jj(l, m) = dylm(l, m);
                }
            }
        }
    } //end loop over neighbours

    //complex conjugate A's (for NEGATIVE (-m) terms)
    // for rank > 1
    for (elej = 0; elej < basis_set->nelements; elej++) {
        for (n = 0; n < nradiali; n++) {
            auto &A_lm = A(elej, n);
            for (l = 0; l <= lmaxi; l++) {
                //fill in -m part in the outer loop using the same m <-> -m symmetry as for Ylm
                for (m = 1; m <= l; m++) {
                    factor = m % 2 == 0 ? 1 : -1;
                    A_lm(l, -m) = A_lm(l, m).conjugated() * factor;
                }
            }
        }
    }    //now A's are constructed
    ACE_TIMER_STOP(loop_over_neighbour)

    // ==================== ENERGY ====================

    ACE_TIMER_START(energy_calc)

    //ALGORITHM 2+3+4: B-basis functions with iterative product and density rho(p) calculation
    //rank=1
    for (int f_ind = 0; f_ind < total_basis_size_rank1; ++f_ind) {
        ACECTildeBasisFunction *func = &basis_rank1[f_ind];
        ndensity = func->ndensity;
#ifdef PRINT_LOOPS_INDICES
        cout << "Num density = " << (int) ndensity << " r = 0 " << endl;
        print_C_tilde_B_basis_function(*func);
#endif
        double A_cur = A_rank1(func->mus[0], func->ns[0] - 1);
#ifdef DEBUG_ENERGY_CALCULATIONS
        printf("A_r=1(x=%d, n=%d)=(%f)\n", func->mus[0], func->ns[0], A_cur);
#endif
        for (DENSITY_TYPE p = 0; p < ndensity; ++p) {
            //for rank=1 (r=0) only 1 ms-combination exists (ms_ind=0), so index of func.ctildes is 0..ndensity-1
            rhos(p) += A_cur * func->ctildes[p];
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
        cout << "Num density = " << (int) ndensity << " r = " << (int) r << endl;
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
            for (j = 0; j < rank; j++) {
                //TODO: optimize ns[j]-1 -> ns[j] during functions construction
                A_cache[j] = A(mus[j], ns[j] - 1, ls[j], ms[j]);
#ifdef DEBUG_ENERGY_CALCULATIONS
                printf("A(x=%d, n=%d, l=%d, m=%d)=(%f,%f)\n", mus[j], ns[j], ls[j], ms[j], A_cache[j].real,
                       A_cache[j].img);
#endif
                A_forward_prod[j + 1] = A_forward_prod[j] * A_cache[j];
            }

            B = A_forward_prod[j];
#ifdef DEBUG_FORCES_CALCULATIONS
            printf("B = (%f, %f)\n", (B).real, (B).img);
#endif
            //fill backward A-product triangle
            for (j = r; j >= 1; j--) {
                A_backward_prod[j - 1] =
                        A_backward_prod[j] * A_cache[j];
            }

            for (j = 0; j < rank; ++j, ++func_ms_t_ind) {
                dB = A_forward_prod[j] * A_backward_prod[j]; //dB - product of all A's except j-th
                dB_flatten(func_ms_t_ind) = dB;
#ifdef DEBUG_FORCES_CALCULATIONS
                m_t = ms[j];
                printf("dB(n,l,m)(%d,%d,%d) = (%f, %f)\n", ns[j], ls[j], m_t, (dB).real, (dB).img);
#endif
            }

            for (DENSITY_TYPE p = 0; p < ndensity; ++p) {
                //real-part only multiplication
                rhos(p) += B.real_part_product(func->ctildes[ms_ind * ndensity + p]);
#ifdef PRINT_INTERMEDIATE_VALUES
                printf("rhos(%d) += %f\n", p, B.real_part_product(func->ctildes[ms_ind * ndensity + p]));
                cout << "Rho[i = " << i << "][p = " << p << "] = " << rhos(p) << endl;
#endif
            }
        }//end of loop over {ms} combinations in sum
    }// end loop for rank>1

#ifdef DEBUG_FORCES_CALCULATIONS
    printf("rhos = ");
    for(DENSITY_TYPE p =0; p<ndensity; ++p) printf(" %f ",rhos(p));
    printf("\n");
#endif

    basis_set->FS_values_and_derivatives(rhos, evdwl, dF_drho, ndensity);

#ifdef DEBUG_FORCES_CALCULATIONS
    printf("dFrhos = ");
    for(DENSITY_TYPE p =0; p<ndensity; ++p) printf(" %f ",dF_drho(p));
    printf("\n");
#endif

    //Algorithm 5 + 6 - weights and theta calculation
    // rank = 1
    for (int f_ind = 0; f_ind < total_basis_size_rank1; ++f_ind) {
        ACECTildeBasisFunction *func = &basis_rank1[f_ind];
        ndensity = func->ndensity;
        for (DENSITY_TYPE p = 0; p < ndensity; ++p) {
            //TODO: multiply by dF/dRho (on store in different densities-dimensions)
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
                theta += dF_drho(p) * func->ctildes[ms_ind * ndensity + p]; //*0.5 ?
#ifdef DEBUG_FORCES_CALCULATIONS
                printf("(p=%d) theta += dF_drho[p] * func.ctildes[ms_ind * ndensity + p] = %f * %f = %f\n",p, dF_drho(p), func->ctildes[ms_ind * ndensity + p],dF_drho(p)*func->ctildes[ms_ind * ndensity + p]);
                printf("theta=%f\n",theta);
#endif
            }

            for (j = 0; j < rank; ++j, ++func_ms_t_ind) {
                m_t = ms[j];
                factor = (m_t % 2 == 0 ? 1 : -1);
                dB = dB_flatten(func_ms_t_ind);
                weights(mus[j], ns[j] - 1, ls[j], m_t) += dB * theta * 0.5; // Theta_array(func_ms_ind);
                // update -m_t (that could also be positive), because the basis is half_basis
                weights(mus[j], ns[j] - 1, ls[j], -m_t) +=
                        (dB).conjugated() * factor * theta * 0.5;// Theta_array(func_ms_ind);

#ifdef DEBUG_FORCES_CALCULATIONS
                printf("dB(n,l,m)(%d,%d,%d) = (%f, %f)\n", ns[j], ls[j], m_t, (dB).real, (dB).img);
                printf("theta = %f\n",theta);
                printf("weights(n,l,m)(%d,%d,%d) += (%f, %f)\n", ns[j], ls[j], m_t, (dB * theta * 0.5).real,
                       (dB * theta * 0.5).img);
                printf("weights(n,l,-m)(%d,%d,%d) += (%f, %f)\n", ns[j], ls[j], -m_t,
                       ( (dB).conjugated() * factor * theta * 0.5).real,
                       ( (dB).conjugated() * factor * theta * 0.5).img);
#endif
            }
        }
    }
    ACE_TIMER_STOP(energy_calc)

// ==================== FORCES ====================
#ifdef PRINT_MAIN_STEPS
    cout << endl << "FORCE CALCULATION" <<
         endl;
    cout << "loop over neighbours" <<
         endl;
#endif
    ACE_TIMER_START(forces_calc_loop)
// loop over neighbour atoms for force calculations
    for (jj = 0; jj < jnum; ++jj) {
        elej = elements[jj];
        r_hat = rhats[jj];
        inv_r_norm = inv_r_norms[jj];

        Array2DLM<ACEComplex> &Y_cache_jj = Y_cache(jj);
        Array2DLM<Dycomponent> &DY_cache_jj = DY_cache(jj);

#ifdef PRINT_LOOPS_INDICES
        cout << endl << "neighbour atom #" << jj <<
             endl;
        printf("rhat = (%f, %f, %f)\n", r_hat[0], r_hat[1], r_hat[2]);
#endif

        ACE_TIMER_START(forces_calc_neighbour);

        f_ji[0] = f_ji[1] = f_ji[2] = 0;

//for rank = 1
        for (n = 0; n < nradbasei; ++n) {
            if (weights_rank1(elej, n) == 0)
                continue;
            auto &DG = DG_cache(jj, n);
            DGR = DG * Y00;
            DGR *= weights_rank1(elej, n);
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
                    auto w = weights(elej, n, l, m);
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
        printf("f_ji(k=%d, i=%d)=(%f, %f, %f)\n", jj, i,
               f_ji[0], f_ji[1], f_ji[2]
        );
#endif

        neighbours_forces(jj, 0) = f_ji[0];
        neighbours_forces(jj, 1) = f_ji[1];
        neighbours_forces(jj, 2) = f_ji[2];
        ACE_TIMER_STOP(forces_calc_neighbour);
    }// end loop over neighbour atoms for forces
    ACE_TIMER_STOP(forces_calc_loop)

    //now, energies and forces are ready
    energies(i) = evdwl;
#ifdef PRINT_INTERMEDIATE_VALUES
    cout << "energies(i) = FS(...rho_p_accum...) = " << evdwl <<
         endl;
#endif
    ACE_TIMER_STOP(per_atom_calc)
}
