/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
Copyright 2021 Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
^2: University of Cambridge, Cambridge, United Kingdom
^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
^4: University of British Columbia, Vancouver, BC, Canada
*/

//
// Created by Lysogorskiy Yury on 27.02.20.
//

#include "pair_grace_fs.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cstring>
#include <exception>

#include "ace-evaluator/ace_version.h"
#include "ace/grace_fs_evaluator.h"

#include "utils_pace.h"

namespace LAMMPS_NS {
    struct ACEImpl {
        ACEImpl() : basis_set(nullptr), ace(nullptr) {}

        ~ACEImpl() {
            delete basis_set;
            delete ace;
        }

        GRACEFSBasisSet *basis_set;
        GRACEFSBEvaluator *ace;
    };
}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;


/* ---------------------------------------------------------------------- */
PairGRACEFS::PairGRACEFS(LAMMPS *lmp) : Pair(lmp) {
    single_enable = 0;
    restartinfo = 0;
    one_coeff = 1;
    manybody_flag = 1;

    aceimpl = new ACEImpl;
    flag_compute_extrapolation_grade = 0;
    extrapolation_grade_gamma = nullptr;
    scale = nullptr;

    chunksize = 4096;

    centroidstressflag = CENTROID_AVAIL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairGRACEFS::~PairGRACEFS() {
    if (copymode) return;

    delete aceimpl;

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(scale);
        memory->destroy(extrapolation_grade_gamma);
    }
}

/* ---------------------------------------------------------------------- */

void PairGRACEFS::compute(int eflag, int vflag) {
    int i, j, ii, jj, inum, jnum;
    double delx, dely, delz, evdwl;
    double fij[3];
    int *ilist, *jlist, *numneigh, **firstneigh;

    ev_init(eflag, vflag);

    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;

    // number of atoms in cell
    int nlocal = atom->nlocal;
    int newton_pair = force->newton_pair;

    // inum: length of the neighborlists list
    inum = list->inum;

    // ilist: list of "i" atoms for which neighbor lists exist
    ilist = list->ilist;

    //numneigh: the length of each these neigbor list
    numneigh = list->numneigh;

    // the pointer to the list of neighbors of "i"
    firstneigh = list->firstneigh;

    if (flag_compute_extrapolation_grade && atom->nlocal > nmax) {
        memory->destroy(extrapolation_grade_gamma);
        nmax = atom->nlocal;
        memory->create(extrapolation_grade_gamma, nmax, "pace/atom:gamma");
        //zeroify array
        memset(extrapolation_grade_gamma, 0, nmax * sizeof(*extrapolation_grade_gamma));
    }
    //determine the maximum number of neighbours
    int max_jnum = 0;
    int nei = 0;
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        jnum = numneigh[i];
        nei = nei + jnum;
        if (jnum > max_jnum) max_jnum = jnum;
    }

    aceimpl->ace->resize_neighbours_cache(max_jnum);
    int *my_neigh_jlist = new int[max_jnum];

    //loop over atoms
    for (ii = 0; ii < inum; ii++) {
        i = list->ilist[ii];
        const int itype = type[i];

        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];

        jlist = firstneigh[i];
        jnum = numneigh[i];

        // checking if neighbours are actually within cutoff range is done inside compute_atom
        // mapping from LAMMPS atom types ('type' array) to ACE species is done inside compute_atom
        //      by using 'aceimpl->ace->element_type_mapping' array
        // x: [r0 ,r1, r2, ..., r100]
        // i = 0 ,1
        // jnum(0) = 50
        // jlist(neigh ind of 0-atom) = [1,2,10,7,99,25, .. 50 element in total]

        // apply NEIGHMASK
        for (jj = 0; jj < jnum; ++jj) {
            my_neigh_jlist[jj] = jlist[jj] & NEIGHMASK;
        }
        try {
            aceimpl->ace->compute_projections = flag_compute_extrapolation_grade;
            aceimpl->ace->compute_atom(i, x, type, jnum, my_neigh_jlist);
        } catch (std::exception &e) {
            error->one(FLERR, e.what());
        }
        // 'compute_atom' will update the `aceimpl->ace->e_atom` and `aceimpl->ace->neighbours_forces(jj, alpha)` arrays

        if (flag_compute_extrapolation_grade)
            extrapolation_grade_gamma[i] = aceimpl->ace->max_gamma_grade;

        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];

            fij[0] = scale[itype][itype] * aceimpl->ace->neighbours_forces(jj, 0);
            fij[1] = scale[itype][itype] * aceimpl->ace->neighbours_forces(jj, 1);
            fij[2] = scale[itype][itype] * aceimpl->ace->neighbours_forces(jj, 2);

            f[i][0] += fij[0];
            f[i][1] += fij[1];
            f[i][2] += fij[2];
            f[j][0] -= fij[0];
            f[j][1] -= fij[1];
            f[j][2] -= fij[2];

            // tally per-atom virial contribution
            if (vflag_either)
                ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fij[0], fij[1], fij[2], delx, dely,
                             delz);

                // Centroid Stress
                if (cvflag_atom) {
                    double fx=fij[0], fy=fij[1], fz=fij[2];

                    cvatom[i][0] += 0.5*delx * fx; // xx
                    cvatom[i][1] += 0.5*dely * fy; // yy
                    cvatom[i][2] += 0.5*delz * fz; // zz
                    cvatom[i][3] += 0.5*delx * fy;  // xy
                    cvatom[i][4] += 0.5*delx * fz; // xz
                    cvatom[i][5] += 0.5*dely * fz; // yz
                    cvatom[i][6] += 0.5*dely * fx; // yx
                    cvatom[i][7] += 0.5*delz * fx; // zx
                    cvatom[i][8] += 0.5*delz * fy; // zy


                    cvatom[j][0] += 0.5*delx * fx; // xx
                    cvatom[j][1] += 0.5*dely * fy; // yy
                    cvatom[j][2] += 0.5*delz * fz; // zz
                    cvatom[j][3] += 0.5*delx * fy;  // xy
                    cvatom[j][4] += 0.5*delx * fz; // xz
                    cvatom[j][5] += 0.5*dely * fz; // yz
                    cvatom[j][6] += 0.5*dely * fx; // yx
                    cvatom[j][7] += 0.5*delz * fx; // zx
                    cvatom[j][8] += 0.5*delz * fy; // zy
                }
        }

        // tally energy contribution
        if (eflag_either) {
            // evdwl = energy of atom I
            evdwl = scale[itype][itype] * aceimpl->ace->e_atom;
            ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
    }

    if (vflag_fdotr) virial_fdotr_compute();
    delete[] my_neigh_jlist;

    // end modifications YL
}

/* ---------------------------------------------------------------------- */

void PairGRACEFS::allocate() {
    allocated = 1;
    int n = atom->ntypes + 1;

    memory->create(setflag, n, n, "pair:setflag");
    memory->create(cutsq, n, n, "pair:cutsq");
    memory->create(scale, n, n, "pair:scale");
    map = new int[n];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGRACEFS::settings(int narg, char **arg) {
    if (narg > 3) utils::missing_cmd_args(FLERR, "pair_style grace/fs", error);

    // ACE potentials are parameterized in metal units
    if (strcmp("metal", update->unit_style) != 0)
        error->all(FLERR, "GRACE/FS potentials require 'metal' units");


    if (comm->me == 0) {
        utils::logmesg(lmp, "[GRACE/FS]\n");
    }

    int iarg = 0;
    while (iarg < narg) {
        if (strcmp(arg[iarg], "chunksize") == 0) {
            chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;
        } else if (strcmp(arg[iarg], "extrapolation") == 0) {
            request_extrapolation = true;
            iarg += 1;
        } else
            error->all(FLERR, "Unknown pair_style grace_fs keyword: {}", arg[iarg]);
    }

    if (comm->me == 0) {
        utils::logmesg(lmp, "[GRACE-FS] Version: {}.{}.{}\n", VERSION_YEAR, VERSION_MONTH, VERSION_DAY);
//    if (recursive)
//      utils::logmesg(lmp, "Recursive evaluator is used\n");
//    else
        utils::logmesg(lmp, "[GRACE-FS] Product evaluator is used\n");
    }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGRACEFS::coeff(int narg, char **arg) {

    if (!allocated) allocate();

    int arg_shift = 2;
    auto potential_file_name = utils::get_potential_file_path(arg[2]);

    //load potential file
    delete aceimpl->basis_set;
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE-FS] Loading {}\n", potential_file_name);
    aceimpl->basis_set = new GRACEFSBasisSet(potential_file_name);

    std::string asi_file_name;
    if (request_extrapolation) {
        asi_file_name = utils::get_potential_file_path(arg[3]);
        arg_shift += 1;
    }

    map_element2type(narg - (arg_shift + 1), arg + (arg_shift + 1));

    if (comm->me == 0) {
        utils::logmesg(lmp, "[GRACE-FS] Supporting {} elements: ", aceimpl->basis_set->nelements);
        for (SPECIES_TYPE mu = 0; mu < aceimpl->basis_set->nelements; mu++) {
            utils::logmesg(lmp, "{} ", aceimpl->basis_set->elements_name[mu]);
        }
        utils::logmesg(lmp, "\n");
        utils::logmesg(lmp, "[GRACE-FS] Number of functions per element: {} \n", aceimpl->basis_set->basis[0].size());
    }

    // read args that map atom types to PACE elements
    // map[i] = which element the Ith atom type is, -1 if not mapped
    // map[0] is not used

    delete aceimpl->ace;
    aceimpl->ace = new GRACEFSBEvaluator();
    aceimpl->ace->element_type_mapping.init(atom->ntypes + 1);

    const int n = atom->ntypes;
    for (int i = 1; i <= n; i++) {
        char *elemname = arg[arg_shift + i];
        if (strcmp(elemname, "NULL") == 0) {
            // species_type=-1 value will not reach ACE Evaluator::compute_atom,
            // but if it will ,then error will be thrown there
            aceimpl->ace->element_type_mapping(i) = -1;
            map[i] = -1;
            if (comm->me == 0) utils::logmesg(lmp, "[GRACE-FS] Skipping LAMMPS atom type #{}(NULL)\n", i);
        } else {
            int atomic_number = PACE::AtomicNumberByName(elemname);
            if (atomic_number == -1) error->all(FLERR, "'{}' is not a valid element\n", elemname);
            SPECIES_TYPE mu = aceimpl->basis_set->get_species_index_by_name(elemname);
            if (mu != -1) {
                if (comm->me == 0)
                    utils::logmesg(lmp, "[GRACE-FS] Mapping LAMMPS atom type #{}({}) -> ACE species type #{}\n", i,
                                   elemname, mu);
                map[i] = mu;
                // set up LAMMPS atom type to ACE species  mapping for ace evaluator
                aceimpl->ace->element_type_mapping(i) = mu;
            } else {
                error->all(FLERR, "[GRACE-FS] Element {} is not supported by ACE-potential from file {}", elemname,
                           potential_file_name);
            }
        }
    }

    // initialize scale factor
    for (int i = 1; i <= n; i++) {
        for (int j = i; j <= n; j++) scale[i][j] = 1.0;
    }

    aceimpl->ace->set_basis(*aceimpl->basis_set);

    if (request_extrapolation) {
        if (comm->me == 0) utils::logmesg(lmp, "[GRACE-FS] Loading ASI {}\n", asi_file_name);
        aceimpl->ace->load_active_set(asi_file_name);
    }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGRACEFS::init_style() {
    if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace/fs requires atom IDs");
    if (force->newton_pair == 0) error->all(FLERR, "Pair style grace/fs requires newton pair on");

    // request a full neighbor list
    neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGRACEFS::init_one(int i, int j) {
    if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
    //cutoff from the basis set's radial functions settings
    scale[j][i] = scale[i][j];
    return aceimpl->basis_set->radial_functions.rcut; // the same for all
}

/* ----------------------------------------------------------------------
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairGRACEFS::extract(const char *str, int &dim) {
    dim = 0;
    //check if str=="gamma_flag" then compute extrapolation grades on this iteration
    if (strcmp(str, "gamma_flag") == 0) return (void *) &flag_compute_extrapolation_grade;

    dim = 2;
    if (strcmp(str, "scale") == 0) return (void *) scale;
    return nullptr;
}

/* ----------------------------------------------------------------------
   peratom requests from FixPair
   return ptr to requested data
   also return ncol = # of quantites per atom
     0 = per-atom vector
     1 or more = # of columns in per-atom array
   return NULL if str is not recognized
---------------------------------------------------------------------- */
void *PairGRACEFS::extract_peratom(const char *str, int &ncol) {
    if (strcmp(str, "gamma") == 0) {
        ncol = 0;
        return (void *) extrapolation_grade_gamma;
    }

//  if (strcmp(str, "corerep") == 0) {
//    ncol = 0;
//    return (void *) corerep_factor;
//  }

    return nullptr;
}
