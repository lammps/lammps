/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
Copyright 2022 Yury Lysogorskiy^1, Anton Bochkarev^1, Matous Mrovec^1, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
 */


//
// Created by Lysogorskiy Yury on 2.01.22.
//

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include "pair_pace_al.h"
#include "atom.h"
#include "dump_custom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "update.h"


#include "math_const.h"

#include "ace_b_evaluator.h"
#include "ace_b_basis.h"
#include "ace_recursive.h"
#include "ace_version.h"

namespace LAMMPS_NS {
    struct ACEALImpl {
        ACEALImpl() : basis_set(nullptr), ace(nullptr), ctilde_basis_set(nullptr), rec_ace(nullptr) {}
        ~ACEALImpl()
        {
            delete basis_set;
            delete ace;
        }
        ACEBBasisSet *basis_set;
        ACEBEvaluator* ace;
        ACECTildeBasisSet *ctilde_basis_set;
        ACERecursiveEvaluator* rec_ace;
    };
}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define DELTA 4
#define PACE_AL_EXTRAPOLATION_GRADE_FNAME "grade.dat"
//added YL



int elements_num_pace_al = 104;
char const *const elements_pace_al[104] = {"X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
                                           "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
                                           "Mn",
                                           "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
                                           "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
                                           "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
                                           "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
                                           "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                                           "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
};

int AtomicNumberByName_pace_al(char *elname) {
    for (int i = 1; i < elements_num_pace_al; i++)
        if (strcmp(elname, elements_pace_al[i]) == 0)
            return i;
    return -1;
}

/* ----------------------------------------------------------------------
 * Append extrapolation grade to file
 ---------------------------------------------------------------------- */
void dump_extrapolation_grade_header() {
    FILE *gamma_file = fopen(PACE_AL_EXTRAPOLATION_GRADE_FNAME, "w");
    fprintf(gamma_file, "Step\tgamma\n");
    fclose(gamma_file);
}

/* ----------------------------------------------------------------------
 * Append extrapolation grade to file
 ---------------------------------------------------------------------- */
void dump_extrapolation_grade(int timestep, double gamma) {
    FILE *gamma_file = fopen(PACE_AL_EXTRAPOLATION_GRADE_FNAME, "a");
    fprintf(gamma_file, "%d\t%f\n", timestep, gamma);
    fclose(gamma_file);
}


/* ---------------------------------------------------------------------- */
PairPACEActiveLearning::PairPACEActiveLearning(LAMMPS *lmp) : Pair(lmp) {
    single_enable = 0;
    restartinfo = 0;
    one_coeff = 1;
    manybody_flag = 1;


    aceimpl = new ACEALImpl;
    scale= nullptr;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairPACEActiveLearning::~PairPACEActiveLearning() {
    if (copymode) return;

    delete aceimpl;

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(scale);
    }

    if (dump)
        delete dump;
}

/* ---------------------------------------------------------------------- */

void PairPACEActiveLearning::compute(int eflag, int vflag) {
    int i, j, ii, jj, inum, jnum;
    double delx, dely, delz, evdwl;
    double fij[3];
    int *ilist, *jlist, *numneigh, **firstneigh;
    double gamma_grade = 0;
    double global_gamma_grade = 0;
    ev_init(eflag, vflag);

    // downwards modified by YL

    double **x = atom->x;
    double **f = atom->f;
    tagint *tag = atom->tag;
    int *type = atom->type;
    // number of atoms in cell
    int nlocal = atom->nlocal;

    int newton_pair = force->newton_pair;

    // number of atoms including ghost atoms
    int nall = nlocal + atom->nghost;

    // inum: length of the neighborlists list
    inum = list->inum;

    // ilist: list of "i" atoms for which neighbor lists exist
    ilist = list->ilist;

    //numneigh: the length of each these neigbor list
    numneigh = list->numneigh;

    // the pointer to the list of neighbors of "i"
    firstneigh = list->firstneigh;

    if (inum != nlocal) {
        char str[128];
        snprintf(str, 128, "inum: %d nlocal: %d are different", inum, nlocal);
        error->all(FLERR, str);
    }


    // Aidan Thompson told RD (26 July 2019) that practically always holds:
    // inum = nlocal
    // i = ilist(ii) < inum
    // j = jlist(jj) < nall
    // neighborlist contains neighbor atoms plus skin atoms,
    //       skin atoms can be removed by setting skin to zero but here
    //       they are disregarded anyway


    //determine the maximum number of neighbours
    int max_jnum = 0;
    int nei = 0;
    for (ii = 0; ii < list->inum; ii++) {
        i = ilist[ii];
        jnum = numneigh[i];
        nei = nei + jnum;
        if (jnum > max_jnum)
            max_jnum = jnum;
    }

    int current_timestep = update->ntimestep;
    bool is_bevaluator = current_timestep % gamma_grade_eval_freq == 0;

    if (is_bevaluator)
        aceimpl->ace->resize_neighbours_cache(max_jnum);
    else
        aceimpl->rec_ace->resize_neighbours_cache(max_jnum);

    //loop over atoms
    for (ii = 0; ii < list->inum; ii++) {
        i = list->ilist[ii];
        const int itype = type[i];

        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];

        jlist = firstneigh[i];
        jnum = numneigh[i];

        // checking if neighbours are actually within cutoff range is done inside compute_atom
        // mapping from LAMMPS atom types ('type' array) to ACE species is done inside compute_atom
        //      by using 'ace->element_type_mapping' array
        // x: [r0 ,r1, r2, ..., r100]
        // i = 0 ,1
        // jnum(0) = 50
        // jlist(neigh ind of 0-atom) = [1,2,10,7,99,25, .. 50 element in total]
        try {
            if (is_bevaluator)
                aceimpl->ace->compute_atom(i, x, type, jnum, jlist);
            else
                aceimpl->rec_ace->compute_atom(i, x, type, jnum, jlist);
        } catch (std::exception &e) {
            error->all(FLERR, e.what());
            exit(EXIT_FAILURE);
        }
        // 'compute_atom' will update the `ace->e_atom` and `ace->neighbours_forces(jj, alpha)` arrays and max_gamma_grade

        if (is_bevaluator) {
            if (gamma_grade < aceimpl->ace->max_gamma_grade)
                gamma_grade = aceimpl->ace->max_gamma_grade;
        }

        Array2D<DOUBLE_TYPE> &neighbours_forces = (is_bevaluator ? aceimpl->ace->neighbours_forces : aceimpl->rec_ace->neighbours_forces);
        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            const int jtype = type[j];
            j &= NEIGHMASK;
            delx = x[j][0] - xtmp;
            dely = x[j][1] - ytmp;
            delz = x[j][2] - ztmp;

            fij[0] = scale[itype][jtype] * neighbours_forces(jj, 0);
            fij[1] = scale[itype][jtype] * neighbours_forces(jj, 1);
            fij[2] = scale[itype][jtype] * neighbours_forces(jj, 2);


            f[i][0] += fij[0];
            f[i][1] += fij[1];
            f[i][2] += fij[2];
            f[j][0] -= fij[0];
            f[j][1] -= fij[1];
            f[j][2] -= fij[2];

            // tally per-atom virial contribution
            if (vflag)
                ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0,
                             fij[0], fij[1], fij[2],
                             -delx, -dely, -delz);
        }

        // tally energy contribution
        if (eflag) {
            // evdwl = energy of atom I
            DOUBLE_TYPE e_atom;
            if (is_bevaluator)
                e_atom= aceimpl->ace->e_atom;
            else
                e_atom= aceimpl->rec_ace->e_atom;
            evdwl = scale[itype][itype] * e_atom;
            ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
    }


    if (vflag_fdotr) virial_fdotr_compute();

    //TODO: check correctness of MPI usage, maybe use comm->me==0 instead ?
    if (is_bevaluator) {
        //gather together global_gamma_grade
        MPI_Allreduce(&gamma_grade, &global_gamma_grade, 1, MPI_DOUBLE, MPI_MAX, world);
        int mpi_rank;
        MPI_Comm_rank(world, &mpi_rank);

        if (global_gamma_grade > gamma_upper_bound) {
            if (mpi_rank == 0) dump_extrapolation_grade(update->ntimestep, global_gamma_grade);
            dump->write();
            MPI_Barrier(world);

            if (mpi_rank == 0) {
                error->all(FLERR, "Extrapolation grade is too large, stopping...\n");
            }

            MPI_Abort(world, 1); //abort properly with error code '1' if not using 4 processes
            exit(EXIT_FAILURE);
        } else if (global_gamma_grade > gamma_lower_bound) {
            if (mpi_rank == 0) dump_extrapolation_grade(update->ntimestep, global_gamma_grade);
            dump->write();
        }
    }

// end modifications YL
}

/* ---------------------------------------------------------------------- */

void PairPACEActiveLearning::allocate() {
    allocated = 1;
    int n = atom->ntypes;

    memory->create(setflag, n + 1, n + 1, "pair:setflag");
    memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
    //TODO: remove ?
    memory->create(map, n + 1, "pair:map");
    memory->create(scale, n + 1, n + 1, "pair:scale");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPACEActiveLearning::settings(int narg, char **arg) {
    if (narg > 3) {
        error->all(FLERR,
                   "Illegal pair_style command. Correct form:\n\tpair_style pace/al [gamma_lower_bound] [gamma_upper_bound]");
    }

    if (narg > 0) {
        double glb = atof(arg[0]); // gamma lower bound
        if (glb < 1.0)
            error->all(FLERR,
                       "Illegal gamma_lower_bound value: it should be real number >= 1.0");
        else
            gamma_lower_bound = glb;
    }

    if (narg > 1) {
        double gub = atof(arg[1]); // gamma upper bound
        if (gub < gamma_lower_bound)
            error->all(FLERR,
                       "Illegal gamma_upper_bound value: it should be real number >= gamma_lower_bound >= 1.0");
        else
            gamma_upper_bound = gub;
    }

    if (narg > 2) {
        gamma_grade_eval_freq = atoi(arg[2]);
        if (gamma_grade_eval_freq < 1)
            error->all(FLERR,
                       "Illegal gamma_grade_eval_freq value: it should be integer number >= 1");
    }

    if (comm->me == 0) {
        if (screen) {
            utils::logmesg(lmp, "ACE/AL version: {}.{}.{}\n", VERSION_YEAR, VERSION_MONTH, VERSION_DAY);
            utils::logmesg(lmp, "Extrapolation grade thresholds (lower/upper): {}/{}\n", gamma_lower_bound,
                    gamma_upper_bound);
            utils::logmesg(lmp, "Extrapolation grade evaluation frequency: {}\n", gamma_grade_eval_freq);
        }
    }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPACEActiveLearning::coeff(int narg, char **arg) {

    if (narg < 5)
        error->all(FLERR,
                   "Incorrect args for pair coefficients. Correct form:\npair_coeff * * <potential.yaml> <potential.asi> elem1 elem2 ...");

    if (!allocated) allocate();

    map_element2type(narg - 4, arg + 4);

    auto potential_file_name = utils::get_potential_file_path(arg[2]);
    auto active_set_inv_filename = utils::get_potential_file_path(arg[3]);
    char **elemtypes = &arg[4];


    delete aceimpl->basis_set;
    delete aceimpl->ctilde_basis_set;

    //load potential file
    aceimpl->basis_set = new ACEBBasisSet();
    if (comm->me == 0) utils::logmesg(lmp, "Loading {}\n", potential_file_name);
    aceimpl->basis_set->load(potential_file_name);

    //convert the basis set to CTilde format
    aceimpl->ctilde_basis_set = new ACECTildeBasisSet();
    *aceimpl->ctilde_basis_set = aceimpl->basis_set->to_ACECTildeBasisSet();

    if (comm->me == 0) {
        utils::logmesg(lmp, "Total number of basis functions\n");

        for (SPECIES_TYPE mu = 0; mu < aceimpl->basis_set->nelements; mu++) {
            int n_r1 = aceimpl->basis_set->total_basis_size_rank1[mu];
            int n = aceimpl->basis_set->total_basis_size[mu];
            utils::logmesg(lmp, "\t{}: {} (r=1) {} (r>1)\n", aceimpl->basis_set->elements_name[mu], n_r1,
                           n);
        }
    }


    // read args that map atom types to PACE elements
    // map[i] = which element the Ith atom type is, -1 if not mapped
    // map[0] is not used
    delete aceimpl->ace;
    delete aceimpl->rec_ace;

    aceimpl->ace = new ACEBEvaluator();
    aceimpl->ace->element_type_mapping.init(atom->ntypes + 1);

    aceimpl->rec_ace = new ACERecursiveEvaluator();
    aceimpl->rec_ace->set_recursive(true);
    aceimpl->rec_ace->element_type_mapping.init(atom->ntypes + 1);
    aceimpl->rec_ace->element_type_mapping.fill(-1); //-1 means atom not included into potential

    FILE *species_type_file = fopen("species_types.dat", "w");
    const int n = atom->ntypes;
    for (int i = 1; i <= n; i++) {
        char *elemname = elemtypes[i - 1];
        if (strcmp(elemname, "NULL") == 0) {
            // species_type=-1 value will not reach ACE Evaluator::compute_atom,
            // but if it will ,then error will be thrown there
            aceimpl->ace->element_type_mapping(i) = -1;
            map[i] = -1;
            if (comm->me == 0) utils::logmesg(lmp, "Skipping LAMMPS atom type #{}(NULL)\n", i);
        } else {

            // dump species types for reconstruction of atomic configurations
            fprintf(species_type_file, "%s ", elemname);
            int atomic_number = AtomicNumberByName_pace_al(elemname);
            if (atomic_number == -1) error->all(FLERR, "'{}' is not a valid element\n", elemname);
            SPECIES_TYPE mu = aceimpl->basis_set->get_species_index_by_name(elemname);
            if (mu != -1) {
                if (comm->me == 0)
                    utils::logmesg(lmp, "Mapping LAMMPS atom type #{}({}) -> ACE species type #{}\n", i,
                                   elemname, mu);
                map[i] = mu;
                // set up LAMMPS atom type to ACE species  mapping for ace evaluators
                aceimpl->ace->element_type_mapping(i) = mu;
                aceimpl->rec_ace->element_type_mapping(i) = mu;
            } else {
                error->all(FLERR, "Element {} is not supported by ACE-potential from file {}", elemname,
                           potential_file_name);
            }
        }
    }
    fclose(species_type_file);
    aceimpl->ace->set_basis(*aceimpl->basis_set);
    aceimpl->rec_ace->set_basis(*aceimpl->ctilde_basis_set);

    if (comm->me == 0) utils::logmesg(lmp, "Loading ASI {}\n", active_set_inv_filename);
    aceimpl->ace->load_active_set(active_set_inv_filename);

    // clear setflag since coeff() called once with I,J = * *
    for (int i = 1; i <= n; i++) 
        for (int j = i; j <= n; j++) 
            scale[i][j] = 1.0;
    
    // prepare dump class
    if (!dump) {
        // dump WRITE_DUMP all cfg 10 dump.snap.*.cfg mass type xs ys zs
        // dump WRITE_DUMP all custom freq extrapolation.dat id type mass x y z
        char **dumpargs = new char *[11];
        dumpargs[0] = (char *) "WRITE_DUMP"; // dump id
        dumpargs[1] = (char *) "all";                // group
        dumpargs[2] = (char *) "custom";                // dump style
        dumpargs[3] = (char *) "1";          // dump frequency
        dumpargs[4] = (char *) "extrapolative_structures.dat";          // fname
        dumpargs[5] = (char *) "id";
        dumpargs[6] = (char *) "type";
        dumpargs[7] = (char *) "mass";
        dumpargs[8] = (char *) "x";
        dumpargs[9] = (char *) "y";
        dumpargs[10] = (char *) "z";
        dump = new DumpCustom(lmp, 11, dumpargs);
        dump->init();

        // dump_modify WRITE_DUMP element X Y Z
        char **dumpargs3 = new char *[atom->ntypes + 1];
        dumpargs3[0] = (char *) "element";
        for (int k = 0; k < atom->ntypes; k++)
            dumpargs3[k + 1] = elemtypes[k];
        dump->modify_params(atom->ntypes + 1, dumpargs3);
    }

    // write extrapolation_grade.dat file header
    if (comm->me == 0)
        dump_extrapolation_grade_header();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPACEActiveLearning::init_style() {
    if (atom->tag_enable == 0)
        error->all(FLERR, "Pair style PACE requires atom IDs");
    if (force->newton_pair == 0)
        error->all(FLERR, "Pair style PACE requires newton pair on");

    // request a full neighbor list
    neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPACEActiveLearning::init_one(int i, int j) {
    if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
    //cutoff from the basis set's radial functions settings
    scale[j][i] = scale[i][j];
    return aceimpl->basis_set->radial_functions->cut(map[i], map[j]);
}

/* ---------------------------------------------------------------------- 
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairPACEActiveLearning::extract(const char *str, int &dim) {
    dim = 2;
    if (strcmp(str, "scale") == 0) return (void *) scale;
    return nullptr;
}

