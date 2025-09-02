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
This file is a modified version of src/ML-PACE/pair_pace.cpp.

Original copyright:
Copyright 2021 Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
^2: University of Cambridge, Cambridge, United Kingdom
^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
^4: University of British Columbia, Vancouver, BC, Canada
*/

//
// Originally created by Lysogorskiy Yury on 27.02.20.
// Adaptive precision added by David Immel in 2025
// (d.immel@fz-juelich.de, FZJ, Germany).
//

#include "pair_pace_apip.h"

#include "atom.h"
#include "atom_vec_apip.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <exception>
#include <cstring>

#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_recursive.h"
#include "ace-evaluator/ace_version.h"
#include "ace/ace_b_basis.h"

namespace LAMMPS_NS {
struct ACEImpl {
  ACEImpl() : basis_set(nullptr), ace(nullptr) {}
  ~ACEImpl()
  {
    delete basis_set;
    delete ace;
  }
  ACECTildeBasisSet *basis_set;
  ACERecursiveEvaluator *ace;
};
}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

static char const *const elements_pace[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};
static constexpr int elements_num_pace = sizeof(elements_pace) / sizeof(const char *);

static int AtomicNumberByName_pace(char *elname)
{
  for (int i = 1; i < elements_num_pace; i++)
    if (strcmp(elname, elements_pace[i]) == 0) return i;
  return -1;
}

/* ---------------------------------------------------------------------- */
PairPACEAPIP::PairPACEAPIP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nmax_corerep = 0;
  flag_corerep_factor = 0;
  corerep_factor = nullptr;

  aceimpl = new ACEImpl;
  recursive = false;

  scale = nullptr;

  chunksize = 4096;

  // start of adaptive-precision modifications by DI
  lambda_thermostat = true;

  n_computations_accumulated = 0;
  time_wall_accumulated = 0;
  time_per_atom = -1;
  // end of adaptive-precision modifications by DI
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairPACEAPIP::~PairPACEAPIP()
{
  if (copymode) return;

  delete aceimpl;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
    memory->destroy(corerep_factor);
  }
}

/**
  * Set lambda_required based on lambda and lambda_const
  * @return true if this calculation is not required
  */

// written by DI. This function is required for the adaptive-precision.
int PairPACEAPIP::check_abort_condition(double *lambda, double *lambda_const, int *lambda_required,
                                        int i)
{
  if ((lambda[i] == 1) && ((!lambda_thermostat) || (lambda_thermostat && lambda_const[i] == 1))) {
    lambda_required[i] |= ApipLambdaRequired::NO_COMPLEX;
    return 1;
  }
  lambda_required[i] |= ApipLambdaRequired::COMPLEX;
  return 0;
}

/**
  * @return prefactor 1-lambda which is used for a precise ACE potential
  */

// written by DI. This function is required for the adaptive-precision.
double PairPACEAPIP::compute_factor_lambda(double lambda)
{
  return 1 - lambda;
}

/**
  * @return atom->apip_e_precise which is used for a precise ACE potential
  */

// written by DI. This function is required for the adaptive-precision.
double *PairPACEAPIP::get_e_ref_ptr()
{
  return atom->apip_e_precise;
}

/* ---------------------------------------------------------------------- */

void PairPACEAPIP::compute(int eflag, int vflag)
{

  // start of adaptive-precision modifications by DI
  // start timers
  double time_wall_start = platform::walltime();
  double *lambda = atom->apip_lambda;
  int *lambda_required = atom->apip_lambda_required;

  double **f_const_lambda = nullptr;
  double **f_dyn_lambda = nullptr;
  double *e_ref = nullptr;
  double *lambda_const = nullptr;
  if (lambda_thermostat) {
    f_const_lambda = atom->apip_f_const_lambda;
    f_dyn_lambda = atom->apip_f_dyn_lambda;
    e_ref = get_e_ref_ptr();
    lambda_const = atom->apip_lambda_const;
  }
  int n_computations = 0;
  // end of adaptive-precision modifications by DI

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

  if (flag_corerep_factor && atom->nlocal > nmax_corerep) {
    memory->destroy(corerep_factor);
    nmax_corerep = atom->nlocal;
    memory->create(corerep_factor, nmax_corerep, "pace/atom:corerep_factor");
    //zeroify array
    memset(corerep_factor, 0, nmax_corerep * sizeof(*corerep_factor));
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

  //loop over atoms
  for (ii = 0; ii < inum; ii++) {
    i = list->ilist[ii];
    const int itype = type[i];

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // start of adaptive-precision modifications by DI
    // Abort calculation when ace is not required for this atom.
    // All force and energy contributions calculated in the following are weighted with
    // 1-lambda when ace is used as precise potential (and with lambda if ace is used as simple potential).
    // As this weighting factor can be 0, i.e. there is no contribution, one can abort the calculation
    // of this atom in this case.
    if (check_abort_condition(lambda, lambda_const, lambda_required, i)) { continue; }
    n_computations++;
    // set factor required for forces and for energy summation
    // fast potential: lambda
    // precise potential: 1-lambda
    const double factor_lambda_i = compute_factor_lambda(lambda[i]);
    const double factor_lambdaconst_i =
        (lambda_thermostat ? compute_factor_lambda(lambda_const[i]) : 0);
    // end of adaptive-precision modifications by DI

    // checking if neighbours are actually within cutoff range is done inside compute_atom
    // mapping from LAMMPS atom types ('type' array) to ACE species is done inside compute_atom
    //      by using 'aceimpl->ace->element_type_mapping' array
    // x: [r0 ,r1, r2, ..., r100]
    // i = 0 ,1
    // jnum(0) = 50
    // jlist(neigh ind of 0-atom) = [1,2,10,7,99,25, .. 50 element in total]

    try {
      aceimpl->ace->compute_atom(i, x, type, jnum, jlist);
    } catch (std::exception &e) {
      error->one(FLERR, e.what());
    }

    if (flag_corerep_factor) corerep_factor[i] = 1 - aceimpl->ace->ace_fcut;

    // 'compute_atom' will update the `aceimpl->ace->e_atom` and `aceimpl->ace->neighbours_forces(jj, alpha)` arrays

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;

      fij[0] = scale[itype][itype] * aceimpl->ace->neighbours_forces(jj, 0);
      fij[1] = scale[itype][itype] * aceimpl->ace->neighbours_forces(jj, 1);
      fij[2] = scale[itype][itype] * aceimpl->ace->neighbours_forces(jj, 2);

      // start of adaptive-precision modifications by DI
      // The force contributions fij need to be weighted with 1-lambda[i]
      // (or lambda[i] in case of a fast ace potential).
      f[i][0] += factor_lambda_i * fij[0];
      f[i][1] += factor_lambda_i * fij[1];
      f[i][2] += factor_lambda_i * fij[2];
      f[j][0] -= factor_lambda_i * fij[0];
      f[j][1] -= factor_lambda_i * fij[1];
      f[j][2] -= factor_lambda_i * fij[2];

      if (lambda_thermostat) {
        f_dyn_lambda[i][0] += factor_lambda_i * fij[0];
        f_dyn_lambda[i][1] += factor_lambda_i * fij[1];
        f_dyn_lambda[i][2] += factor_lambda_i * fij[2];
        f_dyn_lambda[j][0] -= factor_lambda_i * fij[0];
        f_dyn_lambda[j][1] -= factor_lambda_i * fij[1];
        f_dyn_lambda[j][2] -= factor_lambda_i * fij[2];
        f_const_lambda[i][0] += factor_lambdaconst_i * fij[0];
        f_const_lambda[i][1] += factor_lambdaconst_i * fij[1];
        f_const_lambda[i][2] += factor_lambdaconst_i * fij[2];
        f_const_lambda[j][0] -= factor_lambdaconst_i * fij[0];
        f_const_lambda[j][1] -= factor_lambdaconst_i * fij[1];
        f_const_lambda[j][2] -= factor_lambdaconst_i * fij[2];
      }
      // end of adaptive-precision modifications by DI

      // tally per-atom virial contribution
      if (vflag_either)
        // following line of code modified by DI
        ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, factor_lambda_i * fij[0],
                     factor_lambda_i * fij[1], factor_lambda_i * fij[2], -delx, -dely, -delz);
    }

    // tally energy contribution
    // start of adaptive-precision modifications by DI
    if (eflag_either || lambda_thermostat) {
      // The potential energy needs to be stored to apply the
      // energy correction with the local thermostat.
      if (e_ref) e_ref[i] = scale[itype][itype] * aceimpl->ace->e_atom;
      // evdwl = energy of atom I
      // The potential energy is weighted with lambda[i] as well.
      evdwl = factor_lambda_i * scale[itype][itype] * aceimpl->ace->e_atom;
      // end of adaptive-precision modifications by DI
      ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

  // end modifications YL

  // start of adaptive-precision modifications by DI
  // stop timers
  time_wall_accumulated += platform::walltime() - time_wall_start;
  n_computations_accumulated += n_computations;
  // end of adaptive-precision modifications by DI
}

/* ---------------------------------------------------------------------- */

void PairPACEAPIP::allocate()
{
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

void PairPACEAPIP::settings(int narg, char **arg)
{
  if (narg > 3) utils::missing_cmd_args(FLERR, "pair_style pace", error);

  // ACE potentials are parameterized in metal units
  if (strcmp("metal", update->unit_style) != 0)
    error->all(FLERR, "ACE potentials require 'metal' units");

  recursive = true;    // default evaluator style: RECURSIVE

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "recursive") == 0) {
      recursive = true;
      iarg += 1;
    } else if (strcmp(arg[iarg], "product") == 0) {
      recursive = false;
      iarg += 1;
    } else if (strcmp(arg[iarg], "chunksize") == 0) {
      chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Unknown pair_style pace keyword: {}", arg[iarg]);
  }

  if (comm->me == 0) {
    utils::logmesg(lmp, "ACE version: {}.{}.{}\n", VERSION_YEAR, VERSION_MONTH, VERSION_DAY);
    if (recursive)
      utils::logmesg(lmp, "Recursive evaluator is used by ACE\n");
    else
      utils::logmesg(lmp, "Product evaluator is used by ACE\n");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPACEAPIP::coeff(int narg, char **arg)
{

  if (!allocated) allocate();

  map_element2type(narg - 3, arg + 3);

  auto potential_file_name = utils::get_potential_file_path(arg[2]);

  //load potential file
  delete aceimpl->basis_set;
  if (comm->me == 0) utils::logmesg(lmp, "Loading {}\n", potential_file_name);
  // if potential is in ACEBBasisSet (YAML) format, then convert to ACECTildeBasisSet automatically
  if (utils::strmatch(potential_file_name, ".*\\.yaml$")) {
    ACEBBasisSet bBasisSet = ACEBBasisSet(potential_file_name);
    ACECTildeBasisSet cTildeBasisSet = bBasisSet.to_ACECTildeBasisSet();
    aceimpl->basis_set = new ACECTildeBasisSet(cTildeBasisSet);
  } else {
    aceimpl->basis_set = new ACECTildeBasisSet(potential_file_name);
  }

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
  aceimpl->ace = new ACERecursiveEvaluator();
  aceimpl->ace->set_recursive(recursive);
  aceimpl->ace->element_type_mapping.init(atom->ntypes + 1);

  const int n = atom->ntypes;
  for (int i = 1; i <= n; i++) {
    char *elemname = arg[2 + i];
    if (strcmp(elemname, "NULL") == 0) {
      // species_type=-1 value will not reach ACE Evaluator::compute_atom,
      // but if it will ,then error will be thrown there
      aceimpl->ace->element_type_mapping(i) = -1;
      map[i] = -1;
      if (comm->me == 0) utils::logmesg(lmp, "Skipping LAMMPS atom type #{}(NULL)\n", i);
    } else {
      int atomic_number = AtomicNumberByName_pace(elemname);
      if (atomic_number == -1) error->all(FLERR, "'{}' is not a valid element\n", elemname);
      SPECIES_TYPE mu = aceimpl->basis_set->get_species_index_by_name(elemname);
      if (mu != -1) {
        if (comm->me == 0)
          utils::logmesg(lmp, "Mapping LAMMPS atom type #{}({}) -> ACE species type #{}\n", i,
                         elemname, mu);
        map[i] = mu;
        // set up LAMMPS atom type to ACE species  mapping for ace evaluator
        aceimpl->ace->element_type_mapping(i) = mu;
      } else {
        error->all(FLERR, "Element {} is not supported by ACE-potential from file {}", elemname,
                   potential_file_name);
      }
    }
  }

  // initialize scale factor
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) scale[i][j] = 1.0;
  }

  aceimpl->ace->set_basis(*aceimpl->basis_set, 1);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPACEAPIP::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style pace requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style pace requires newton pair on");

  // start of adaptive-precision modifications by DI
  if (!atom->apip_lambda_required_flag)
    error->all(FLERR, "pair style pace/apip requires an atom style with lambda_required.");
  if (!atom->apip_lambda_flag)
    error->all(FLERR, "Pair style pace/apip requires an atom style with lambda");
  // end of adaptive-precision modifications by DI

  // request a full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPACEAPIP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
  //cutoff from the basis set's radial functions settings
  scale[j][i] = scale[i][j];
  return aceimpl->basis_set->radial_functions->cut(map[i], map[j]);
}

/**
  * setup specific to this pair style
  * Determine whether there is a fix lambda_thermostat/apip or not and set
  * lambda_thermostat.
  */

// written by DI. This function is required for the adaptive-precision.
void PairPACEAPIP::setup()
{
  if (modify->get_fix_by_style("^lambda_thermostat/apip$").size() == 0) {
    lambda_thermostat = false;
  } else {
    lambda_thermostat = true;
    if (!atom->apip_lambda_const_flag)
      error->all(
          FLERR,
          "Pair style pace/apip requires an atom style with lambda_const for a local thermostat.");
    if (!atom->apip_e_fast_flag)
      error->all(
          FLERR,
          "Pair style pace/apip requires an atom style with e_simple for a local thermostat.");
    if (!atom->apip_e_precise_flag)
      error->all(
          FLERR,
          "Pair style pace/apip requires an atom style with e_complex for a local thermostat.");
    if (!atom->apip_f_const_lambda_flag)
      error->all(FLERR,
                 "Pair style pace/apip requires an atom style with f_const_lambda for a local "
                 "thermostat.");
    if (!atom->apip_f_dyn_lambda_flag)
      error->all(FLERR,
                 "Pair style pace/apip requires an atom style with f_const_lambda for a local "
                 "thermostat.");
  }
}

/**
  * set return values for timers and number of computed particles
  */

// written by DI. This function is required for the adaptive-precision.
void PairPACEAPIP::calculate_time_per_atom()
{
  if (n_computations_accumulated > 0)
    time_per_atom = time_wall_accumulated / n_computations_accumulated;
  else
    time_per_atom = -1;

  // reset
  time_wall_accumulated = 0;
  n_computations_accumulated = 0;
}

/* ----------------------------------------------------------------------
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairPACEAPIP::extract(const char *str, int &dim)
{
  dim = 0;
  //check if str=="corerep_flag" then compute extrapolation grades on this iteration
  if (strcmp(str, "corerep_flag") == 0) return (void *) &flag_corerep_factor;
  // DI: The following option is required for the adaptive precision.
  if (strcmp(str, "pace/apip:time_per_atom") == 0) {
    calculate_time_per_atom();
    return (void *) &time_per_atom;
  }

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
void *PairPACEAPIP::extract_peratom(const char *str, int &ncol)
{
  if (strcmp(str, "corerep") == 0) {
    ncol = 0;
    return (void *) corerep_factor;
  }

  return nullptr;
}
