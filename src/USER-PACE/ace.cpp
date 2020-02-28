//
// Created by lysogy36 on 31.01.20.
//

#include "ace.h"
#include "ace_types.h"
#include "atomic_environment.h"

void ACE::set_basis(ACEBasisSet &bas) {
    basis_set = &bas;
    init();
}

void ACE::init() {
    A_dense.init(basis_set->nelements, basis_set->nradmax + 1, basis_set->lmax + 1, "A_dense");
    A_rank1.init(basis_set->nelements, basis_set->nradbase, "A_rank1");

    weights_dense.init(basis_set->nelements, basis_set->nradmax + 1, basis_set->lmax + 1,
                       "weights_dense");

    weights_rank1.init(basis_set->nelements, basis_set->nradbase, "weights_rank1");


    DG_dense_cache.init(1, basis_set->nradbase, "DG_dense_cache");
    DG_dense_cache.fill(0);

    R_dense_cache.init(1, basis_set->nradmax, basis_set->lmax + 1, "R_dense_cache");
    R_dense_cache.fill(0);

    DR_dense_cache.init(1, basis_set->nradmax, basis_set->lmax + 1, "DR_dense_cache");
    DR_dense_cache.fill(0);

    Y_dense_cache.init(1, basis_set->lmax + 1, "Y_dense_cache");
    Y_dense_cache.fill({0, 0});

    DY_dense_cache.init(1, basis_set->lmax + 1, "dY_dense_cache");
    DY_dense_cache.fill({0.});


    rhos.init(basis_set->ndensitymax, "rhos");
    dF_dRho.init(basis_set->ndensitymax, "dF_dRho");

    //Theta_array.init(basis_set->max_B_array_size, "Theta_array");
    dB_array.init(basis_set->max_dB_array_size, "dB_array");

}

void ACE::compute(AtomicEnvironment &atomic_environment, bool verbose) {
    TIMER_INIT(loop_over_neighbour)
    TIMER_INIT(forces_calc_loop)
    TIMER_INIT(forces_calc_neighbour)
    TIMER_INIT(phi_calc)
    TIMER_INIT(phi_recalc)
    TIMER_INIT(energy_calc)
    TIMER_INIT(bond_calc)
    TIMER_INIT(A_calc)
    TIMER_INIT(per_atom_calc)
    TIMER_INIT(basis_func_calc)
    TIMER_INIT(total_time_calc)
    TIMER_START(total_time_calc)

    int i, j, jj;
    total_energy = 0;
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

    if (R_dense_cache.get_dim(0) < max_jnum) {

        //TODO: implement grow
        R_dense_cache.resize(max_jnum, basis_set->nradmax, basis_set->lmax + 1);
        R_dense_cache.fill(0);

        DR_dense_cache.resize(max_jnum, basis_set->nradmax, basis_set->lmax + 1);
        DR_dense_cache.fill(0);

        DG_dense_cache.resize(max_jnum, basis_set->nradbase);
        DG_dense_cache.fill(0);

        Y_dense_cache.resize(max_jnum, basis_set->lmax + 1);
        Y_dense_cache.fill({0, 0});

        DY_dense_cache.resize(max_jnum, basis_set->lmax + 1);
        DY_dense_cache.fill({0});
    }

    for (i = 0; i < atomic_environment.n_atoms_tot; ++i) {
        compute_atom(i,
                     atomic_environment.x,
                     atomic_environment.species_type,
                     atomic_environment.num_neighbours[i],
                     atomic_environment.neighbour_list[i]);
        //this will also update the energies(i) and pair_forces_ji(jj, alpha) arrays
        //update global energies and forces accumulators
        total_energy += energies(i);
        for (jj = 0; jj < atomic_environment.num_neighbours[i]; jj++) {
            j = atomic_environment.neighbour_list[i][jj];

            forces(i, 0) += pair_forces_ji(jj, 0);
            forces(i, 1) += pair_forces_ji(jj, 1);
            forces(i, 2) += pair_forces_ji(jj, 2);

            forces(j, 0) -= pair_forces_ji(jj, 0);
            forces(j, 1) -= pair_forces_ji(jj, 1);
            forces(j, 2) -= pair_forces_ji(jj, 2);

#ifdef DEBUG_FORCES_CALCULATIONS
            printf("accumulated forces: F(i=%d)=(%f,%f,%f)\n", i, forces(i, 0), forces(i, 1), forces(i, 2));
            printf("accumulated forces: F(j=%d)=(%f,%f,%f)\n", j, forces(j, 0), forces(j, 1), forces(j, 2));
#endif
        }
    } // loop over atoms (i_at)
    TIMER_STOP(total_time_calc)

#ifdef FINE_TIMING
    if (verbose) {
        std::cout << "   Total time: " << TIMER_MICROSECONDS(total_time_calc) << " microseconds" << endl;
        std::cout << "Per atom time:    "
                  << TIMER_MICROSECONDS(per_atom_calc) / atomic_environment.n_atoms_tot
                  << " microseconds" << endl;

        std::cout << "Loop_over_nei/atom:"
                  << TIMER_MICROSECONDS(loop_over_neighbour) / atomic_environment.n_atoms_tot << " microseconds"
                  << endl;
        std::cout << "Energy/atom:      "
                  << TIMER_MICROSECONDS(energy_calc) /
                     atomic_environment.n_atoms_tot << " microseconds" << endl;
//        std::cout << "         Basis_func/atom: "
//                  << TIMER_NANOSECONDS(basis_func_calc) / 1000. /
//                     atomic_environment.n_atoms_tot << " microseconds" << endl;

        std::cout << "Forces/atom:      "
                  << TIMER_MICROSECONDS(forces_calc_loop) / atomic_environment.n_atoms_tot << " microseconds" << endl;
        std::cout << "        phi_recalcs/atom: " << TIMER_MICROSECONDS(phi_recalc) / atomic_environment.n_atoms_tot
                  << " microseconds" << endl;
        std::cout << "             forces_neig: "
                  << TIMER_MICROSECONDS(forces_calc_neighbour) /
                     atomic_environment.n_atoms_tot
                  << " microseconds" << endl;
    }
#endif


}



// double** r - atomic coordinates of atom I
// int* types - atomic types if atom I
// int **firstneigh -  ptr to 1st J int value of each I atom. Usage: jlist = firstneigh[i];
// Usage: j = jlist_of_i[jj];
// jnum - number of J neighbors for each I atom.  jnum = numneigh[i];

void ACE::compute_atom(int i, DOUBLE_TYPE **x, const SPECIES_TYPE *type, int jnum, const int *jlist) {
    TIMER_START(per_atom_calc)
#ifdef PRINT_MAIN_STEPS
    cout << endl << "ATOM: ind = " << i << " r_norm=(" << x[i][0] << ","
         << x[i][1] << ","
         << x[i][2] << ")" << endl;
#endif
    DOUBLE_TYPE evdwl = 0;
    DOUBLE_TYPE r_norm;
    DOUBLE_TYPE xn, yn, zn, r_xyz;
    DOUBLE_TYPE R, GR, DGR, R_over_r, DR_cache;
    DOUBLE_TYPE *r_hat;

    SPECIES_TYPE elej;
    DENSITY_TYPE ndensity; //TODO: extract from basis set, as it is equal to all functions
    RANK_TYPE r, rank, t;
    NS_TYPE n;
    LS_TYPE l;
    MS_TYPE m, m_t;

    SPECIES_TYPE *mus;
    NS_TYPE *ns;
    LS_TYPE *ls;
    MS_TYPE *ms;

    int jj, func_ind, ms_ind;
    SHORT_INT_TYPE factor;

    Complex Y{0}, Y_DR{0.};
    Complex B{0.};
    Complex dB{0};
    Complex A_cache[basis_set->rankmax];

    //Theta_array.fill({0.});
    dB_array.fill({0.});

    Dycomponent grad_phi_nlm{0}, DY_cache{0.};

    //size is +1 of max to avoid out-of-boundary array access in double-triangular scheme
    Complex A_forward_prod[basis_set->rankmax + 1];
    Complex A_backward_prod[basis_set->rankmax + 1];

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

    C_tilde_B_basis_function *basis_rank1 = basis_set->basis_rank1[mu];
    C_tilde_B_basis_function *basis = basis_set->basis[mu];

    //TODO: lmax -> lmaxi
    const LS_TYPE lmaxi = basis_set->lmax;

    //TODO: nradmax -> nradiali
    const NS_TYPE nradiali = basis_set->nradmax;

    //TODO: nradbase -> nradbasei
    const NS_TYPE nradbasei = basis_set->nradbase;

    pair_forces_ji.resize(jnum, 3);
    pair_forces_ji.fill(0);

    weights_dense.fill({0});
    weights_rank1.fill(0);
    A_dense.fill({0});
    A_rank1.fill(0);
    rhos.fill(0);
    dF_dRho.fill(0);

    //proxy references to spherical harmonics and radial functions arrays
    const Array2DLM <Complex> &ylm = basis_set->spherical_harmonics.ylm;
    const Array2DLM <Dycomponent> &dylm = basis_set->spherical_harmonics.dylm;

    const Array2D<DOUBLE_TYPE> &fr = basis_set->radial_functions.fr;
    const Array2D<DOUBLE_TYPE> &dfr = basis_set->radial_functions.dfr;

    const Array1D<DOUBLE_TYPE> &gr = basis_set->radial_functions.gr;
    const Array1D<DOUBLE_TYPE> &dgr = basis_set->radial_functions.dgr;

    TIMER_START(loop_over_neighbour)

    //loop over neighbours
    for (jj = 0; jj < jnum; ++jj) {

        t = jlist[jj];

        xn = x[t][0] - xtmp;
        yn = x[t][1] - ytmp;
        zn = x[t][2] - ztmp;

        r_xyz = sqrt(xn * xn + yn * yn + zn * zn);

        rhats[jj][0] = xn / r_xyz;
        rhats[jj][1] = yn / r_xyz;
        rhats[jj][2] = zn / r_xyz;

        r_norms[jj] = r_xyz;

        elements[jj] = map_lammps_at_type_to_element[type[t]];
    }


    //ALGORITHM 1: Atomic base construction
    for (jj = 0; jj < jnum; ++jj) {
        r_norm = r_norms[jj];
        elej = elements[jj];
        r_hat = rhats[jj];


        inv_r_norm = 1. / r_norm;
        inv_r_norms[jj] = inv_r_norm;
        //proxies
        Array2DLM <Complex> &Y_dense_cache_jj = Y_dense_cache(jj);
        Array2DLM <Dycomponent> &DY_dense_cache_jj = DY_dense_cache(jj);


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
            DG_dense_cache(jj, n) = dgr(n);
            A_rank1(elej, n) += GR * Y00;
        }
        //loop for computing A's
        // for rank > 1
        for (n = 0; n < nradiali; n++) {
            auto &A_lm = A_dense(elej, n);
            for (l = 0; l <= lmaxi; l++) {
                R = fr(l, n);
#ifdef DEBUG_ENERGY_CALCULATIONS
                printf("R(nl=%d,%d)(r=%f)=%f\n", n + 1, l, r_norm, R);
#endif

                DR_dense_cache(jj, n, l) = dfr(l, n);
                R_dense_cache(jj, n, l) = R;

                for (m = 0; m <= l; m++) {
                    Y = ylm(l, m);
#ifdef DEBUG_ENERGY_CALCULATIONS
                    printf("Y(lm=%d,%d)=(%f, %f)\n", l, m, Y.real, Y.img);
#endif
                    A_lm(l, m) += Y * R; //accumulation sum over neighbours
                    Y_dense_cache_jj(l, m) = Y;
                    DY_dense_cache_jj(l, m) = dylm(l, m);
                }
            }
        }
    } //end loop over neighbours

    //complex conjugate A's (for NEGATIVE (-m) terms)
    // for rank > 1
    for (elej = 0; elej < basis_set->nelements; elej++) {
        for (n = 0; n < nradiali; n++) {
            auto &A_lm = A_dense(elej, n);
            for (l = 0; l <= lmaxi; l++) {
                //fill in -m part in the outer loop using the same m <-> -m symmetry as for Ylm
                for (m = 1; m <= l; m++) {
                    factor = m % 2 == 0 ? 1 : -1;
                    A_lm(l, -m) = A_lm(l, m).conjugated() * factor;
                }
            }
        }
    }    //now A's are constructed
    TIMER_STOP(loop_over_neighbour)

    // ==================== ENERGY ====================

    TIMER_START(energy_calc)

    //ALGORITHM 2+3+4: B-basis functions with iterative product and density rho(p) calculation
    //r=0, rank=1
    for (int f_ind = 0; f_ind < total_basis_size_rank1; ++f_ind) {
        C_tilde_B_basis_function *func = &basis_rank1[f_ind];
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
    }

    //r>0, rank>1
    int func_ms_ind = 0;
    int func_ms_t_ind = 0;// index for dB

    for (func_ind = 0; func_ind < total_basis_size; ++func_ind) {
        C_tilde_B_basis_function *func = &basis[func_ind];
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
        for (ms_ind = 0; ms_ind < func->num_of_ms_combinations; ++ms_ind, ++func_ms_ind) {
            ms = &func->ms[ms_ind * rank];

            //loop over m, collect B  = product of A with given ms
            A_forward_prod[0] = 1;
            A_backward_prod[r] = 1;

            //fill forward A-product triangle
            for (t = 0; t < rank; t++) {
                //TODO: optimize ns[t]-1 -> ns[t] during functions construction
                A_cache[t] = A_dense(mus[t], ns[t] - 1, ls[t], ms[t]);
#ifdef DEBUG_ENERGY_CALCULATIONS
                printf("A_dense(x=%d, n=%d, l=%d, m=%d)=(%f,%f)\n", mus[t], ns[t], ls[t], ms[t], A_cache[t].real,
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
                dB_array(func_ms_t_ind) = dB;
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
                cout << "Rho[i = " << i << "][p = " << p << "] = " << rhos(p) << endl;
#endif
            }
        }//end of loop over {ms} combinations in sum
    }// loop over current basis function

#ifdef DEBUG_FORCES_CALCULATIONS
    printf("rhos = ");
    for(DENSITY_TYPE p =0; p<ndensity; ++p) printf(" %f ",rhos(p));
    printf("\n");
#endif

    basis_set->FS_values_and_derivatives(rhos, evdwl, dF_dRho, ndensity);

#ifdef DEBUG_FORCES_CALCULATIONS
    printf("dFrhos = ");
    for(DENSITY_TYPE p =0; p<ndensity; ++p) printf(" %f ",dF_dRho(p));
    printf("\n");
#endif

    //Algorithm 5 + 6 - weights and theta calculation
    //r=0, rank = 1
    for (int f_ind = 0; f_ind < total_basis_size_rank1; ++f_ind) {
        C_tilde_B_basis_function *func = &basis_rank1[f_ind];
        ndensity = func->ndensity;
        for (DENSITY_TYPE p = 0; p < ndensity; ++p) {
            //TODO: multiply by dF/dRho (on store in different densities-dimensions)
            //for rank=1 (r=0) only 1 ms-combination exists (ms_ind=0), so index of func.ctildes is 0..ndensity-1
            weights_rank1(func->mus[0], func->ns[0] - 1) += dF_dRho(p) * func->ctildes[p];
        }
    }
    //r>0, rank>1
    func_ms_ind = 0;
    func_ms_t_ind = 0;// index for dB
    DOUBLE_TYPE theta = 0;
    for (func_ind = 0; func_ind < total_basis_size; ++func_ind) {
        C_tilde_B_basis_function *func = &basis[func_ind];
        ndensity = func->ndensity;
        rank = func->rank;
        mus = func->mus;
        ns = func->ns;
        ls = func->ls;
        for (ms_ind = 0; ms_ind < func->num_of_ms_combinations; ++ms_ind, ++func_ms_ind) {
            ms = &func->ms[ms_ind * rank];
            theta = 0;
            for (DENSITY_TYPE p = 0; p < ndensity; ++p) {
                theta += dF_dRho(p) * func->ctildes[ms_ind * ndensity + p]; //*0.5 ?
#ifdef DEBUG_FORCES_CALCULATIONS
                printf("(p=%d) theta += dF_dRho[p] * func.ctildes[ms_ind * ndensity + p] = %f * %f = %f\n",p, dF_dRho(p), func->ctildes[ms_ind * ndensity + p],dF_dRho(p)*func->ctildes[ms_ind * ndensity + p]);
                printf("theta=%f\n",theta);
#endif
            }

            for (t = 0; t < rank; ++t, ++func_ms_t_ind) {
                m_t = ms[t];
                factor = (m_t % 2 == 0 ? 1 : -1);
                dB = dB_array(func_ms_t_ind);
                weights_dense(mus[t], ns[t] - 1, ls[t], m_t) += dB * theta * 0.5;// Theta_array(func_ms_ind);
                // update -m_t (that could also be positive), because the basis is half_basis
                weights_dense(mus[t], ns[t] - 1, ls[t], -m_t) +=
                        (dB).conjugated() * factor * theta * 0.5;// Theta_array(func_ms_ind);

#ifdef DEBUG_FORCES_CALCULATIONS
                printf("dB(n,l,m)(%d,%d,%d) = (%f, %f)\n", ns[t], ls[t], m_t, (dB).real, (dB).img);
                printf("theta = %f\n",theta);
                printf("weights_dense(n,l,m)(%d,%d,%d) += (%f, %f)\n", ns[t], ls[t], m_t, (dB * theta * 0.5).real,
                       (dB * theta * 0.5).img);
                printf("weights_dense(n,l,-m)(%d,%d,%d) += (%f, %f)\n", ns[t], ls[t], -m_t,
                       ( (dB).conjugated() * factor * theta * 0.5).real,
                       ( (dB).conjugated() * factor * theta * 0.5).img);
#endif
            }
        }
    }
    TIMER_STOP(energy_calc)

// ==================== FORCES ====================
#ifdef PRINT_MAIN_STEPS
    cout << endl << "FORCE CALCULATION" <<
         endl;
    cout << "loop over neighbours" <<
         endl;
#endif
    TIMER_START(forces_calc_loop)
// loop over neighbour atoms for force calculations
    for (jj = 0; jj < jnum; ++jj) {
        elej = elements[jj];
        r_hat = rhats[jj];
        inv_r_norm = inv_r_norms[jj];

        Array2DLM <Complex> &Y_dense_cache_jj = Y_dense_cache(jj);
        Array2DLM <Dycomponent> &DY_dense_cache_jj = DY_dense_cache(jj);

#ifdef PRINT_LOOPS_INDICES
        cout << endl << "neighbour atom #" << jj <<
             endl;
        printf("rhat = (%f, %f, %f)\n", r_hat[0], r_hat[1], r_hat[2]);
#endif

        TIMER_START(forces_calc_neighbour);

        f_ji[0] = f_ji[1] = f_ji[2] = 0;

//for rank = 1
        for (n = 0; n < nradbasei; ++n) {
            if (weights_rank1(elej, n) == 0)
                continue;
            auto &DG_cache = DG_dense_cache(jj, n);
            DGR = DG_cache * Y00;
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
                R_over_r = R_dense_cache(jj, n, l) * inv_r_norm;
                DR_cache = DR_dense_cache(jj, n, l);
//m<=-1
                for (m = 0; m <= l; m++) {
                    auto w = weights_dense(elej, n, l, m);
                    if (w == 0)
                        continue;
//counting for -m cases if m>0
                    if (m > 0) w *= 2;
                    DY_cache = DY_dense_cache_jj(l, m);
                    Y_DR = Y_dense_cache_jj(l, m) * DR_cache;

                    grad_phi_nlm.a[0] = Y_DR * r_hat[0] + DY_cache.a[0] * R_over_r;
                    grad_phi_nlm.a[1] = Y_DR * r_hat[1] + DY_cache.a[1] * R_over_r;
                    grad_phi_nlm.a[2] = Y_DR * r_hat[2] + DY_cache.a[2] * R_over_r;
#ifdef DEBUG_FORCES_CALCULATIONS
                    printf("d_phi(n=%d, l=%d, m=%d) = ((%f,%f), (%f,%f), (%f,%f))\n",n+1,l,m,
                           grad_phi_nlm.a[0].real, grad_phi_nlm.a[0].img,
                           grad_phi_nlm.a[1].real, grad_phi_nlm.a[1].img,
                           grad_phi_nlm.a[2].real, grad_phi_nlm.a[2].img);

                    printf("weights_dense(n,l,m)(%d,%d,%d) = (%f,%f)\n", n+1, l, m,w.real, w.img);
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
#ifdef DEBUG_FORCES_CALCULATIONS
                    // printf("accumulated f_ji = (%f, %f, %f)\n", f_ji[0], f_ji[1], f_ji[2]);
#endif
                }
            }
        }

#ifdef PRINT_INTERMEDIATE_VALUES
        printf("f_ji(k=%d, i=%d)=(%f, %f, %f)\n", jj, i,
               f_ji[0], f_ji[1], f_ji[2]
        );
#endif

        pair_forces_ji(jj, 0) = f_ji[0];
        pair_forces_ji(jj, 1) = f_ji[1];
        pair_forces_ji(jj, 2) = f_ji[2];
        TIMER_STOP(forces_calc_neighbour);
    }// end loop over neighbour atoms for forces
    TIMER_STOP(forces_calc_loop)

    //now, energies and forces are ready
    energies(i) = evdwl;
#ifdef PRINT_INTERMEDIATE_VALUES
    cout << "energies(i) = FS(...rho_p_accum...) = " << evdwl <<
         endl;
#endif
    TIMER_STOP(per_atom_calc)
}
