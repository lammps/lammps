#include "ace_c_basis.h"

using namespace std;

ACEBasisSet::ACEBasisSet(const ACEBasisSet &other) {
    _copy_scalar_memory(other);
    _copy_dynamic_memory(other);
    pack_flatten_basis();
}

void ACEBasisSet::_clean() {
    delete[] full_ns;
    delete[] full_ls;
    delete[] full_Xs;
    delete[] full_ms;
    delete[] full_c_tildes;

    delete[] full_ns_rank1;
    delete[] full_ls_rank1;
    delete[] full_Xs_rank1;
    delete[] full_ms_rank1;
    delete[] full_c_tildes_rank1;

    if (basis_rank1 != nullptr)
        for (SPECIES_TYPE mu = 0; mu < this->nelements; ++mu) {
            delete[] basis_rank1[mu];
        }

    if (basis != nullptr)
        for (SPECIES_TYPE mu = 0; mu < this->nelements; ++mu) {
            delete[] basis[mu];
        }
    delete[] basis;
    delete[] basis_rank1;

    delete[] total_basis_size;
    delete[] total_basis_size_rank1;

    delete[] elements_name;
}

ACEBasisSet::~ACEBasisSet() {
    _clean();
}

ACEBasisSet &ACEBasisSet::operator=(const ACEBasisSet &other) {
    if (this != &other) {
        // deallocate old memory
        _clean();
        //copy scalar values
        _copy_scalar_memory(other);
        //copy dynamic memory
        _copy_dynamic_memory(other);

        //pack basis functions
        pack_flatten_basis();
    }
    return *this;
}

void ACEBasisSet::_copy_scalar_memory(const ACEBasisSet &other) {
    radial_functions = other.radial_functions;
    spherical_harmonics = other.spherical_harmonics;
    lmax = other.lmax;
    nradbase = other.nradbase;
    nradmax = other.nradmax;
    nelements = other.nelements;
    rankmax = other.rankmax;
    ndensitymax = other.ndensitymax;
    num_ctilde_max = other.num_ctilde_max;
    num_ms_combinations_max = other.num_ms_combinations_max;
    ntot = other.ntot;
    rank_array_total_size = other.rank_array_total_size;
    ms_array_total_size = other.ms_array_total_size;
    coeff_array_total_size = other.coeff_array_total_size;

    rank_array_total_size_rank1 = other.rank_array_total_size_rank1;
    coeff_array_total_size_rank1 = other.coeff_array_total_size_rank1;
}

void ACEBasisSet::_copy_dynamic_memory(const ACEBasisSet &other) {//allocate new memory
    basis = new ACECTildeBasisFunction *[nelements];
    basis_rank1 = new ACECTildeBasisFunction *[nelements];

    total_basis_size = new SHORT_INT_TYPE[nelements];
    total_basis_size_rank1 = new SHORT_INT_TYPE[nelements];

    elements_name = new string[nelements];

    //copy
    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        total_basis_size_rank1[mu] = other.total_basis_size_rank1[mu];
        basis_rank1[mu] = new ACECTildeBasisFunction[total_basis_size_rank1[mu]];

        for (size_t i = 0; i < total_basis_size_rank1[mu]; i++) {
            basis_rank1[mu][i] = other.basis_rank1[mu][i];
        }

        total_basis_size[mu] = other.total_basis_size[mu];
        basis[mu] = new ACECTildeBasisFunction[total_basis_size[mu]];
        for (size_t i = 0; i < total_basis_size[mu]; i++) {
            basis[mu][i] = other.basis[mu][i];
        }

        elements_name[mu] = other.elements_name[mu];
    }
}

//embedding function
//case nemb = 1 only implementation
//F = sign(x)*(  ( ( 1 - exp(-w*x**2) )*abs(x) )^m +  m*exp(-w*x**2)*abs(x) )
// !! no prefactor wpre
void Fexp(DOUBLE_TYPE rho, DOUBLE_TYPE mexp, DOUBLE_TYPE &F, DOUBLE_TYPE &DF) {
    DOUBLE_TYPE w = 10.0;
    DOUBLE_TYPE eps = 1e-10;

    if (abs(rho) > eps) {
        DOUBLE_TYPE g, a, omg, y2, y1 = w * rho * rho;
        DOUBLE_TYPE sign_factor = (signbit(rho) ? -1 : 1);
        if (y1 > 30.0) g = 0;
        else g = exp(-y1);

        omg = 1. - g;
        a = abs(rho);
        y1 = pow(omg * a, mexp);
        y2 = mexp * g * a;
        F = sign_factor * (y1 + y2);

        DOUBLE_TYPE dg, da, dy, dy1, dy2;
        dg = -2.0 * w * rho * g;
        da = sign_factor;
        if (abs(y1) < eps) dy = 0.;
        else dy = mexp * y1 / (omg * a);


        dy1 = dy * (-dg * a + omg * da);
        dy2 = mexp * (dg * a + g * da);
        DF = sign_factor * (dy1 + dy2);

    } else {
        F = mexp * rho;
        DF = mexp;
    }
}

void ACEBasisSet::FS_values_and_derivatives(Array1D<DOUBLE_TYPE> &rhos, DOUBLE_TYPE &value,
                                            Array1D<DOUBLE_TYPE> &derivatives, DENSITY_TYPE ndensity) {
    DOUBLE_TYPE F, DF = 0;
    for (int p = 0; p < ndensity; p++) {
        Fexp(rhos(p), parameters.at(p * ndensity + 1), F, DF);
        value += F * parameters.at(p * ndensity + 0); // * weight
        derivatives(p) = DF * parameters.at(p * ndensity + 0);// * weight
    }
}


void ACEBasisSet::inner_cutoff(DOUBLE_TYPE rho_core, DOUBLE_TYPE rho_cut, DOUBLE_TYPE drho_cut,
                               DOUBLE_TYPE &fcut, DOUBLE_TYPE &dfcut) {

    DOUBLE_TYPE rho_low = rho_cut - drho_cut;
    if (rho_core >= rho_cut) {
        fcut = 0;
        dfcut = 0;
    } else if (rho_core <= rho_low) {
        fcut = 1;
        dfcut = 0;
    } else {
        fcut = 0.5 * (1 + cos(M_PI * (rho_core - rho_low) / drho_cut));
        dfcut = -0.5 * sin(M_PI * (rho_core - rho_low) / drho_cut) * M_PI / drho_cut;
    }
}

//re-pack the constinuent dynamic arrays of all basis functions in contiguous arrays
void ACEBasisSet::pack_flatten_basis() {

    //1. compute arrays sizes
    rank_array_total_size_rank1 = 0;
    //ms_array_total_size_rank1 = rank_array_total_size_rank1;
    coeff_array_total_size_rank1 = 0;

    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        if (total_basis_size_rank1[mu] > 0) {
            rank_array_total_size_rank1 += total_basis_size_rank1[mu];

            ACECTildeBasisFunction &func = basis_rank1[mu][0];
            coeff_array_total_size_rank1 += total_basis_size_rank1[mu] * func.ndensity;
        }
    }

    rank_array_total_size = 0;
    coeff_array_total_size = 0;

    ms_array_total_size = 0;
    max_dB_array_size = 0;


    max_B_array_size = 0;

    size_t cur_ms_size = 0;
    size_t cur_ms_rank_size = 0;

    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {

        cur_ms_size = 0;
        cur_ms_rank_size = 0;
        for (int func_ind = 0; func_ind < total_basis_size[mu]; ++func_ind) {
            ACECTildeBasisFunction &func = basis[mu][func_ind];
            rank_array_total_size += func.rank;
            ms_array_total_size += func.rank * func.num_ms_combs;
            coeff_array_total_size += func.ndensity * func.num_ms_combs;

            cur_ms_size += func.num_ms_combs;
            cur_ms_rank_size += func.rank * func.num_ms_combs;
        }

        if (cur_ms_size > max_B_array_size)
            max_B_array_size = cur_ms_size;

        if (cur_ms_rank_size > max_dB_array_size)
            max_dB_array_size = cur_ms_rank_size;
    }

    //2. allocate contiguous arrays
    full_ns_rank1 = new NS_TYPE[rank_array_total_size_rank1];
    full_ls_rank1 = new NS_TYPE[rank_array_total_size_rank1];
    full_Xs_rank1 = new SPECIES_TYPE[rank_array_total_size_rank1];
    full_ms_rank1 = new MS_TYPE[rank_array_total_size_rank1];
    full_c_tildes_rank1 = new DOUBLE_TYPE[coeff_array_total_size_rank1];


    full_ns = new NS_TYPE[rank_array_total_size];
    full_ls = new LS_TYPE[rank_array_total_size];
    full_Xs = new SPECIES_TYPE[rank_array_total_size];
    full_c_tildes = new DOUBLE_TYPE[coeff_array_total_size];
    full_ms = new MS_TYPE[ms_array_total_size];

    //TODO: allocate dB, Theta, B
    //ACEComplex** dB = new ACEComplex*[ms_array_total_size];

    //3. copy the values from private C_tilde_B_basis_function arrays to new contigous space
    //4. clean private memory
    //5. reassign private array pointers

    //r = 0, rank = 1
    size_t rank_array_ind_rank1 = 0;
    size_t coeff_array_ind_rank1 = 0;
    size_t ms_array_ind_rank1 = 0;

    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        for (int func_ind_r1 = 0; func_ind_r1 < total_basis_size_rank1[mu]; ++func_ind_r1) {
            ACECTildeBasisFunction &func = basis_rank1[mu][func_ind_r1];

            //copy values ns from c_tilde_basis_function private memory to contigous memory part
            full_ns_rank1[rank_array_ind_rank1] = func.ns[0];

            //copy values ls from c_tilde_basis_function private memory to contigous memory part
            full_ls_rank1[rank_array_ind_rank1] = func.ls[0];

            //copy values mus from c_tilde_basis_function private memory to contigous memory part
            full_Xs_rank1[rank_array_ind_rank1] = func.mus[0];

            //copy values ctildes from c_tilde_basis_function private memory to contigous memory part
            memcpy(&full_c_tildes_rank1[coeff_array_ind_rank1], func.ctildes,
                   func.ndensity * sizeof(DOUBLE_TYPE));


            //copy values mus from c_tilde_basis_function private memory to contigous memory part
            memcpy(&full_ms_rank1[ms_array_ind_rank1], func.ms_combs,
                   func.num_ms_combs *
                   func.rank * sizeof(MS_TYPE));

            //release memory of each ACECTildeBasisFunction if it is not proxy
            func._clean();

            func.ns = &full_ns_rank1[rank_array_ind_rank1];
            func.ls = &full_ls_rank1[rank_array_ind_rank1];
            func.mus = &full_Xs_rank1[rank_array_ind_rank1];
            func.ctildes = &full_c_tildes_rank1[coeff_array_ind_rank1];
            func.ms_combs = &full_ms_rank1[ms_array_ind_rank1];
            func.is_proxy = true;

            rank_array_ind_rank1 += func.rank;
            ms_array_ind_rank1 += func.rank *
                                  func.num_ms_combs;
            coeff_array_ind_rank1 += func.num_ms_combs * func.ndensity;

            func_ind_r1++;
        }
    }


    //rank>1, r>0
    size_t rank_array_ind = 0;
    size_t coeff_array_ind = 0;
    size_t ms_array_ind = 0;

    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        for (int func_ind = 0; func_ind < total_basis_size[mu]; ++func_ind) {
            ACECTildeBasisFunction &func = basis[mu][func_ind];

            //copy values ns from c_tilde_basis_function private memory to contigous memory part
            memcpy(&full_ns[rank_array_ind], func.ns,
                   func.rank * sizeof(NS_TYPE));
            //copy values ls from c_tilde_basis_function private memory to contigous memory part
            memcpy(&full_ls[rank_array_ind], func.ls,
                   func.rank * sizeof(LS_TYPE));
            //copy values mus from c_tilde_basis_function private memory to contigous memory part
            memcpy(&full_Xs[rank_array_ind], func.mus,
                   func.rank * sizeof(SPECIES_TYPE));

            //copy values ctildes from c_tilde_basis_function private memory to contigous memory part
            memcpy(&full_c_tildes[coeff_array_ind], func.ctildes,
                   func.num_ms_combs * func.ndensity * sizeof(DOUBLE_TYPE));


            //copy values mus from c_tilde_basis_function private memory to contigous memory part
            memcpy(&full_ms[ms_array_ind], func.ms_combs,
                   func.num_ms_combs *
                   func.rank * sizeof(MS_TYPE));

            //release memory of each ACECTildeBasisFunction if it is not proxy
            func._clean();

            func.ns = &full_ns[rank_array_ind];
            func.ls = &full_ls[rank_array_ind];
            func.mus = &full_Xs[rank_array_ind];
            func.ctildes = &full_c_tildes[coeff_array_ind];
            func.ms_combs = &full_ms[ms_array_ind];
            func.is_proxy = true;

            rank_array_ind += func.rank;
            ms_array_ind += func.rank *
                            func.num_ms_combs;
            coeff_array_ind += func.num_ms_combs * func.ndensity;

            func_ind++;
        }
    }
}

void fwrite_c_tilde_b_basis_func(FILE *fptr, ACECTildeBasisFunction &func) {
    RANK_TYPE r;
    fprintf(fptr, "ctilde_basis_func: ");
    fprintf(fptr, "rank=%d ndens=%d mu0=%d ", func.rank, func.ndensity, func.mu0);

    fprintf(fptr, "mu=(");
    for (r = 0; r < func.rank; ++r)
        fprintf(fptr, " %d ", func.mus[r]);
    fprintf(fptr, ")\n");

    fprintf(fptr, "n=(");
    for (r = 0; r < func.rank; ++r)
        fprintf(fptr, " %d ", func.ns[r]);
    fprintf(fptr, ")\n");

    fprintf(fptr, "l=(");
    for (r = 0; r < func.rank; ++r)
        fprintf(fptr, " %d ", func.ls[r]);
    fprintf(fptr, ")\n");

    fprintf(fptr, "num_ms=%d\n", func.num_ms_combs);

    for (int m_ind = 0; m_ind < func.num_ms_combs; m_ind++) {
        fprintf(fptr, "<");
        for (r = 0; r < func.rank; ++r)
            fprintf(fptr, " %d ", func.ms_combs[m_ind * func.rank + r]);
        fprintf(fptr, ">: ");
        for (DENSITY_TYPE p = 0; p < func.ndensity; p++)
            fprintf(fptr, " %.18f ", func.ctildes[m_ind * func.ndensity + p]);
        fprintf(fptr, "\n");
    }

}

void ACEBasisSet::save(const string &filename) {
    FILE *fptr;
    fptr = fopen(filename.c_str(), "w");
    fprintf(fptr, "lmax=%d\n", lmax);
    fprintf(fptr, "nradbase=%d\n", nradbase);
    fprintf(fptr, "nradmax=%d\n", nradmax);
    fprintf(fptr, "nelements=%d\n", nelements);
    fprintf(fptr, "rankmax=%d\n", rankmax);
    fprintf(fptr, "ndensitymax=%d\n", ndensitymax);
    fprintf(fptr, "cutoffmax=%f\n", cutoffmax);

    fprintf(fptr, "ntot=%d\n", ntot);

    fprintf(fptr, "cutoff=%f\n", cutoff);

    fprintf(fptr, "%ld parameters: ", parameters.size());
    for (int i = 0; i < parameters.size(); ++i) {
        fprintf(fptr, " %f", parameters.at(i));
    }
    fprintf(fptr, "\n");

    //hard-core repulsion
    fprintf(fptr, "core repulsion parameters: ");
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i)
        for (SPECIES_TYPE mu_j = 0; mu_j < nelements; ++mu_j)
            fprintf(fptr, "%.18f %.18f\n", radial_functions.prehc(mu_i, mu_j), radial_functions.lambdahc(mu_j, mu_j));

    //hard-core energy cutoff repulsion
    fprintf(fptr, "core energy-cutoff parameters: ");
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i)
        fprintf(fptr, "%.18f %.18f\n", rho_core_cutoffs(mu_i), drho_core_cutoffs(mu_i));

    //elements mapping
    fprintf(fptr, "elements:");
    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu)
        fprintf(fptr, " %s", elements_name[mu].c_str());
    fprintf(fptr, "\n");

    //TODO: radial functions
    //radparameter
    fprintf(fptr, "radparameter=");
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i)
        for (SPECIES_TYPE mu_j = 0; mu_j < nelements; ++mu_j)
            fprintf(fptr, " %.18f", radial_functions.lambda(mu_i, mu_j));
    fprintf(fptr, "\n");

    fprintf(fptr, "cutoff=");
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i)
        for (SPECIES_TYPE mu_j = 0; mu_j < nelements; ++mu_j)
            fprintf(fptr, " %.18f", radial_functions.cut(mu_i, mu_j));
    fprintf(fptr, "\n");

    fprintf(fptr, "dcut=");
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i)
        for (SPECIES_TYPE mu_j = 0; mu_j < nelements; ++mu_j)
            fprintf(fptr, " %.18f", radial_functions.dcut(mu_i, mu_j));
    fprintf(fptr, "\n");

    fprintf(fptr, "crad=");
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i)
        for (SPECIES_TYPE mu_j = 0; mu_j < nelements; ++mu_j) {
            for (NS_TYPE idx = 1; idx <= nradbase; idx++) {
                for (NS_TYPE nr = 1; nr <= nradmax; nr++) {
                    for (LS_TYPE l = 0; l <= lmax; l++) {
                        fprintf(fptr, " %.18f", radial_functions.crad(mu_i, mu_j, l, nr - 1, idx - 1));
                    }
                    fprintf(fptr, "\n");
                }
            }
        }

    fprintf(fptr, "\n");

    //num_c_tilde_max
    fprintf(fptr, "num_c_tilde_max=%d\n", num_ctilde_max);
    fprintf(fptr, "num_ms_combinations_max=%d\n", num_ms_combinations_max);


    //write total_basis_size and total_basis_size_rank1
    fprintf(fptr, "total_basis_size_rank1: ");
    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        fprintf(fptr, "%d ", total_basis_size_rank1[mu]);
    }
    fprintf(fptr, "\n");

    for (SPECIES_TYPE mu = 0; mu < nelements; mu++)
        for (SHORT_INT_TYPE func_ind = 0; func_ind < total_basis_size_rank1[mu]; ++func_ind)
            fwrite_c_tilde_b_basis_func(fptr, basis_rank1[mu][func_ind]);

    fprintf(fptr, "total_basis_size: ");
    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        fprintf(fptr, "%d ", total_basis_size[mu]);
    }
    fprintf(fptr, "\n");

    for (SPECIES_TYPE mu = 0; mu < nelements; mu++)
        for (SHORT_INT_TYPE func_ind = 0; func_ind < total_basis_size[mu]; ++func_ind)
            fwrite_c_tilde_b_basis_func(fptr, basis[mu][func_ind]);


    fclose(fptr);
}

void fread_c_tilde_b_basis_func(FILE *fptr, ACECTildeBasisFunction &func) {
    RANK_TYPE r;
    int res;
    char buf[3][128];

    res = fscanf(fptr, " ctilde_basis_func: ");

    res = fscanf(fptr, "rank=%s ndens=%s mu0=%s ", buf[0], buf[1], buf[2]);
    if (res != 3) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }

    func.rank = (RANK_TYPE) stol(buf[0]);
    func.ndensity = (DENSITY_TYPE) stol(buf[1]);
    func.mu0 = (SPECIES_TYPE) stol(buf[2]);

    func.mus = new SPECIES_TYPE[func.rank];
    func.ns = new NS_TYPE[func.rank];
    func.ls = new LS_TYPE[func.rank];

    res = fscanf(fptr, " mu=(");
    for (r = 0; r < func.rank; ++r) {
        res = fscanf(fptr, "%s", buf[0]);
        if (res != 1) {
            printf("Error while reading file");
            exit(EXIT_FAILURE);
        }
        func.mus[r] = (SPECIES_TYPE) stol(buf[0]);
    }
    res = fscanf(fptr, " )"); // ")"

    res = fscanf(fptr, " n=("); // "n="
    for (r = 0; r < func.rank; ++r) {
        res = fscanf(fptr, "%s", buf[0]);
        if (res != 1) {
            printf("Error while reading file");
            exit(EXIT_FAILURE);
        }

        func.ns[r] = (NS_TYPE) stol(buf[0]);
    }
    res = fscanf(fptr, " )");

    res = fscanf(fptr, " l=(");
    for (r = 0; r < func.rank; ++r) {
        res = fscanf(fptr, "%s", buf[0]);
        if (res != 1) {
            printf("Error while reading file");
            exit(EXIT_FAILURE);
        }
        func.ls[r] = (NS_TYPE) stol(buf[0]);
    }
    res = fscanf(fptr, " )");

    res = fscanf(fptr, " num_ms=%s\n", buf[0]);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    func.num_ms_combs = (SHORT_INT_TYPE) stoi(buf[0]);

    func.ms_combs = new MS_TYPE[func.rank * func.num_ms_combs];
    func.ctildes = new DOUBLE_TYPE[func.ndensity * func.num_ms_combs];

    for (int m_ind = 0; m_ind < func.num_ms_combs; m_ind++) {
        res = fscanf(fptr, " <");
        for (r = 0; r < func.rank; ++r) {
            res = fscanf(fptr, "%s", buf[0]);
            if (res != 1) {
                printf("Error while reading file");
                exit(EXIT_FAILURE);
            }
            func.ms_combs[m_ind * func.rank + r] = stoi(buf[0]);
        }
        res = fscanf(fptr, " >:");
        for (DENSITY_TYPE p = 0; p < func.ndensity; p++) {
            res = fscanf(fptr, "%s", buf[0]);
            if (res != 1) {
                printf("Error while reading file");
                exit(EXIT_FAILURE);
            }
            func.ctildes[m_ind * func.ndensity + p] = stod(buf[0]);
        }
    }
}

void ACEBasisSet::load(string filename) {
    int res;
    FILE *fptr;
    char buffer[1024], buffer2[1024];
    fptr = fopen(filename.c_str(), "r");
    if (fptr == NULL) {
        printf("Could not open file %s.\n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    res = fscanf(fptr, "lmax=%s", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    lmax = stoi(buffer);

    res = fscanf(fptr, " nradbase=");
    res = fscanf(fptr, "%s", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    nradbase = stoi(buffer);

    res = fscanf(fptr, " nradmax=");
    res = fscanf(fptr, "%s", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    nradmax = stoi(buffer);

    res = fscanf(fptr, " nelements=");
    res = fscanf(fptr, "%s", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    nelements = stoi(buffer);

    res = fscanf(fptr, " rankmax=");
    res = fscanf(fptr, "%s", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    rankmax = stoi(buffer);

    res = fscanf(fptr, " ndensitymax=");
    res = fscanf(fptr, "%s", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    ndensitymax = stoi(buffer);


    res = fscanf(fptr, " cutoffmax=");
    res = fscanf(fptr, "%s", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    cutoffmax = stod(buffer);


    res = fscanf(fptr, " ntot=");
    res = fscanf(fptr, "%s", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    ntot = stoi(buffer);

    res = fscanf(fptr, " cutoff=");
    res = fscanf(fptr, "%s", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    //cutoff = stod(buffer);
    //TODO: write-read it properly
    cutoff = cutoffmax;

    int parameters_size;
    res = fscanf(fptr, "%s parameters:", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    parameters_size = stoi(buffer);
    parameters.resize(parameters_size);

    spherical_harmonics.init(lmax);
    radial_functions.init(nradbase, lmax, nradmax,
                          ntot,
                          nelements,
                          cutoff);
    rho_core_cutoffs.init(nelements, "rho_core_cutoffs");
    drho_core_cutoffs.init(nelements, "drho_core_cutoffs");

    for (int i = 0; i < parameters.size(); ++i) {
        res = fscanf(fptr, "%s", buffer);
        if (res != 1) {
            printf("Error while reading file");
            exit(EXIT_FAILURE);
        }
        if (!res) exit(EXIT_FAILURE);
        parameters[i] = stof(buffer);
    }

    //hard-core repulsion
    res = fscanf(fptr, " core repulsion parameters:");
    if (res != 0) {
        printf("Error while reading core repulsion parameters\n");
        exit(EXIT_FAILURE);
    }
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i)
        for (SPECIES_TYPE mu_j = 0; mu_j < nelements; ++mu_j) {
            res = fscanf(fptr, "%s %s", buffer, buffer2);
            if (res != 2) {
                printf("Error while reading file core repulsion parameters (values)\n");
                exit(EXIT_FAILURE);
            }
            radial_functions.prehc(mu_i, mu_j) = stod(buffer);
            radial_functions.lambdahc(mu_i, mu_j) = stod(buffer2);

//            printf("Read: prehc(mu_i, mu_j)=%f\n",radial_functions.prehc(mu_i, mu_j));
//            printf("Read: lambdahc(mu_i, mu_j)=%f\n",radial_functions.lambdahc(mu_i, mu_j));
        }

    //hard-core energy cutoff repulsion
    res = fscanf(fptr, " core energy-cutoff parameters:");
    if (res != 0) {
        printf("Error while reading core energy-cutoff parameters\n");
        exit(EXIT_FAILURE);
    }
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i) {
        res = fscanf(fptr, "%s %s", buffer, buffer2);
        if (res != 2) {
            printf("Error while reading file core energy-cutoff parameters (values)\n");
            exit(EXIT_FAILURE);
        }
        rho_core_cutoffs(mu_i) = stod(buffer);
        drho_core_cutoffs(mu_i) = stod(buffer2);

//        printf("Read: rho_core_cutoffs(mu_i)=%f\n", rho_core_cutoffs(mu_i));
//        printf("Read: drho_core_cutoffs(mu_i)=%f\n", drho_core_cutoffs(mu_i));
    }


    //elements mapping
    elements_name = new string[nelements];
    res = fscanf(fptr, " elements:");
    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        res = fscanf(fptr, "%s", buffer);
        if (res != 1) {
            printf("Error while reading file");
            exit(EXIT_FAILURE);
        }
        elements_name[mu] = buffer;
    }

    //read radial functions parameter

    res = fscanf(fptr, " radparameter=");
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i)
        for (SPECIES_TYPE mu_j = 0; mu_j < nelements; ++mu_j) {
            res = fscanf(fptr, "%s", buffer);
            if (res != 1) {
                printf("Error while reading file");
                exit(EXIT_FAILURE);
            }
            radial_functions.lambda(mu_i, mu_j) = stod(buffer);
        }


    res = fscanf(fptr, " cutoff=");
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i)
        for (SPECIES_TYPE mu_j = 0; mu_j < nelements; ++mu_j) {
            res = fscanf(fptr, "%s", buffer);
            if (res != 1) {
                printf("Error while reading file");
                exit(EXIT_FAILURE);
            }
            radial_functions.cut(mu_i, mu_j) = stod(buffer);
        }


    res = fscanf(fptr, " dcut=");
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i)
        for (SPECIES_TYPE mu_j = 0; mu_j < nelements; ++mu_j) {
            res = fscanf(fptr, " %s", buffer);
            if (res != 1) {
                printf("Error while reading file");
                exit(EXIT_FAILURE);
            }
            radial_functions.dcut(mu_i, mu_j) = stod(buffer);
        }


    res = fscanf(fptr, " crad=");
    for (SPECIES_TYPE mu_i = 0; mu_i < nelements; ++mu_i)
        for (SPECIES_TYPE mu_j = 0; mu_j < nelements; ++mu_j)
            for (NS_TYPE idx = 1; idx <= nradbase; idx++)
                for (NS_TYPE nr = 1; nr <= nradmax; nr++)
                    for (LS_TYPE l = 0; l <= lmax; l++) {
                        res = fscanf(fptr, "%s", buffer);
                        if (res != 1) {
                            printf("Error while reading file");
                            exit(EXIT_FAILURE);
                        }
                        radial_functions.crad(mu_i, mu_j, l, nr - 1, idx - 1) = stod(buffer);
                    }

    radial_functions.setuplookupRadspline();

    //num_c_tilde_max
    res = fscanf(fptr, " num_c_tilde_max=");
    res = fscanf(fptr, "%s\n", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    num_ctilde_max = stol(buffer);

    res = fscanf(fptr, " num_ms_combinations_max=");
    res = fscanf(fptr, "%s", buffer);
    if (res != 1) {
        printf("Error while reading file");
        exit(EXIT_FAILURE);
    }
    num_ms_combinations_max = stol(buffer);

    //read total_basis_size_rank1
    total_basis_size_rank1 = new SHORT_INT_TYPE[nelements];
    basis_rank1 = new ACECTildeBasisFunction *[nelements];
    res = fscanf(fptr, " total_basis_size_rank1: ");


    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        res = fscanf(fptr, "%s", buffer);
        if (res != 1) {
            printf("Error while reading file");
            exit(EXIT_FAILURE);
        }
        total_basis_size_rank1[mu] = stoi(buffer);
        //printf("total_basis_size[%d] = %d\n", mu, total_basis_size_rank1[mu]);
        basis_rank1[mu] = new ACECTildeBasisFunction[total_basis_size_rank1[mu]];
    }

    for (SPECIES_TYPE mu = 0; mu < nelements; mu++)
        for (SHORT_INT_TYPE func_ind = 0; func_ind < total_basis_size_rank1[mu]; ++func_ind) {
            fread_c_tilde_b_basis_func(fptr, basis_rank1[mu][func_ind]);
            //print_C_tilde_B_basis_function(basis_rank1[mu][func_ind]);
        }

    //read total_basis_size
    res = fscanf(fptr, " total_basis_size: ");
    total_basis_size = new SHORT_INT_TYPE[nelements];
    basis = new ACECTildeBasisFunction *[nelements];

    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        res = fscanf(fptr, "%s", buffer);
        if (res != 1) {
            printf("Error while reading file");
            exit(EXIT_FAILURE);
        }
        total_basis_size[mu] = stoi(buffer);
        basis[mu] = new ACECTildeBasisFunction[total_basis_size[mu]];
    }

    for (SPECIES_TYPE mu = 0; mu < nelements; mu++)
        for (SHORT_INT_TYPE func_ind = 0; func_ind < total_basis_size[mu]; ++func_ind) {
            fread_c_tilde_b_basis_func(fptr, basis[mu][func_ind]);
            //print_C_tilde_B_basis_function(basis[mu][func_ind]);
        }

    fclose(fptr);

    pack_flatten_basis();
}

