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
Copyright 2022 Yury Lysogorskiy^1, Anton Bochkarev^1, Matous Mrovec^1, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
 */

//
// Created by Lysogorskiy Yury on 2.01.22.
//

#include "pair_pace_extrapolation.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "ace/ace_b_basis.h"
#include "ace/ace_b_evaluator.h"
#include "ace-evaluator/ace_recursive.h"
#include "ace-evaluator/ace_version.h"

namespace LAMMPS_NS {
struct ACEALImpl {
  ACEALImpl() : basis_set(nullptr), ace(nullptr), ctilde_basis_set(nullptr), rec_ace(nullptr) {}

  ~ACEALImpl()
  {
    delete basis_set;
    delete ace;

    delete ctilde_basis_set;
    delete rec_ace;
  }

  ACEBBasisSet *basis_set;
  ACEBEvaluator *ace;
  ACECTildeBasisSet *ctilde_basis_set;
  ACERecursiveEvaluator *rec_ace;
};
}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

static char const *const elements_pace_al[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};
static constexpr int elements_num_pace_al = sizeof(elements_pace_al) / sizeof(const char *);

int AtomicNumberByName_pace_al(char *elname)
{
  for (int i = 1; i < elements_num_pace_al; i++)
    if (strcmp(elname, elements_pace_al[i]) == 0) return i;
  return -1;
}

/* ---------------------------------------------------------------------- */
PairPACEExtrapolation::PairPACEExtrapolation(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nmax = 0;
  nmax_corerep = 0;

  aceimpl = new ACEALImpl;
  scale = nullptr;
  flag_compute_extrapolation_grade = 0;
  extrapolation_grade_gamma = nullptr;
  flag_corerep_factor = 0;
  corerep_factor = nullptr;

  chunksize = 4096;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairPACEExtrapolation::~PairPACEExtrapolation()
{
  if (copymode) return;

  delete aceimpl;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
    memory->destroy(map);
    memory->destroy(extrapolation_grade_gamma);
    memory->destroy(corerep_factor);
  }
}

/* ---------------------------------------------------------------------- */

void PairPACEExtrapolation::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum;
  double delx, dely, delz, evdwl;
  double fij[3];
  int *ilist, *jlist, *numneigh, **firstneigh;
  ev_init(eflag, vflag);

  // downwards modified by YL

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

  // this happens when used as substyle in pair style hybrid.
  // So this check and error effectively disallows use with pair style hybrid.
  if (inum != nlocal) { error->all(FLERR, "inum: {} nlocal: {} are different", inum, nlocal); }

  //if flag_compute_extrapolation_grade at this iteration then
  // grow extrapolation_grade_gamma array, that store per-atom extrapolation grades
  if (flag_compute_extrapolation_grade && atom->nlocal > nmax) {
    memory->destroy(extrapolation_grade_gamma);
    nmax = atom->nlocal;
    memory->create(extrapolation_grade_gamma, nmax, "pace/atom:gamma");
    //zeroify array
    memset(extrapolation_grade_gamma, 0, nmax * sizeof(*extrapolation_grade_gamma));
  }
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

  if (flag_compute_extrapolation_grade)
    aceimpl->ace->resize_neighbours_cache(max_jnum);
  else
    aceimpl->rec_ace->resize_neighbours_cache(max_jnum);

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
    //      by using 'ace->element_type_mapping' array
    // x: [r0 ,r1, r2, ..., r100]
    // i = 0 ,1
    // jnum(0) = 50
    // jlist(neigh ind of 0-atom) = [1,2,10,7,99,25, .. 50 element in total]
    try {
      if (flag_compute_extrapolation_grade) {
        aceimpl->ace->compute_projections = true;
        aceimpl->ace->compute_atom(i, x, type, jnum, jlist);
      }
      else
        aceimpl->rec_ace->compute_atom(i, x, type, jnum, jlist);
    } catch (std::exception &e) {
      error->one(FLERR, e.what());
    }
    // 'compute_atom' will update the `ace->e_atom` and `ace->neighbours_forces(jj, alpha)` arrays and max_gamma_grade

    if (flag_compute_extrapolation_grade)
      extrapolation_grade_gamma[i] = aceimpl->ace->max_gamma_grade;

    if (flag_corerep_factor) {
      corerep_factor[i] = 1 - (flag_compute_extrapolation_grade ? aceimpl->ace->ace_fcut
                              : aceimpl->rec_ace->ace_fcut);
    }

    Array2D<DOUBLE_TYPE> &neighbours_forces =
        (flag_compute_extrapolation_grade ? aceimpl->ace->neighbours_forces
                                          : aceimpl->rec_ace->neighbours_forces);
    //optionally assign global forces arrays

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
        ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fij[0], fij[1], fij[2], -delx, -dely,
                     -delz);
    }

    // tally energy contribution
    if (eflag) {
      // evdwl = energy of atom I
      DOUBLE_TYPE e_atom;
      if (flag_compute_extrapolation_grade)
        e_atom = aceimpl->ace->e_atom;
      else
        e_atom = aceimpl->rec_ace->e_atom;
      evdwl = scale[itype][itype] * e_atom;
      ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

  // end modifications YL
}

/* ---------------------------------------------------------------------- */

void PairPACEExtrapolation::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(map, np1, "pair:map");
  memory->create(scale, np1, np1, "pair:scale");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPACEExtrapolation::settings(int narg, char **arg)
{
//  if (narg > 2) error->all(FLERR, "Pair style pace/extrapolation supports no keywords");
  if (narg > 2) utils::missing_cmd_args(FLERR, "pair_style pace/extrapolation", error);
  // ACE potentials are parameterized in metal units
  if (strcmp("metal", update->unit_style) != 0)
    error->all(FLERR, "ACE potentials require 'metal' units");

  int iarg = 0;
  while (iarg < narg) {
      if (strcmp(arg[iarg], "chunksize") == 0) {
          chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
          iarg += 2;
      } else
          error->all(FLERR, "Unknown pair_style pace keyword: {}", arg[iarg]);
  }

  if (comm->me == 0)
    utils::logmesg(lmp, "ACE/AL version: {}.{}.{}\n", VERSION_YEAR, VERSION_MONTH, VERSION_DAY);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPACEExtrapolation::coeff(int narg, char **arg)
{

  if (narg < 5) utils::missing_cmd_args(FLERR, "pair_coeff", error);

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
  aceimpl->rec_ace->element_type_mapping.fill(-1);    //-1 means atom not included into potential


  const int n = atom->ntypes;
  element_names.resize(n);
  for (int i = 1; i <= n; i++) {
    char *elemname = elemtypes[i - 1];
    element_names[i - 1] = elemname;
    if (strcmp(elemname, "NULL") == 0) {
      // species_type=-1 value will not reach ACE Evaluator::compute_atom,
      // but if it will ,then error will be thrown there
      aceimpl->ace->element_type_mapping(i) = -1;
      map[i] = -1;
      if (comm->me == 0) utils::logmesg(lmp, "Skipping LAMMPS atom type #{}(NULL)\n", i);
    } else {
      // dump species types for reconstruction of atomic configurations
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

  aceimpl->ace->set_basis(*aceimpl->basis_set);
  aceimpl->rec_ace->set_basis(*aceimpl->ctilde_basis_set);

  if (comm->me == 0) utils::logmesg(lmp, "Loading ASI {}\n", active_set_inv_filename);
  aceimpl->ace->load_active_set(active_set_inv_filename);
  bool is_linear_extrapolation_grade = aceimpl->ace->get_is_linear_extrapolation_grade();
  if (comm->me == 0) {
    if (is_linear_extrapolation_grade)
      utils::logmesg(lmp, "LINEAR ASI is loaded\n");
    else
      utils::logmesg(lmp, "FULL ASI is loaded\n");
  }

  // clear setflag since coeff() called once with I,J = * *
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) scale[i][j] = 1.0;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPACEExtrapolation::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style pace/extrapolation requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR, "Pair style pace/extrapolation requires newton pair on");

  // request a full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPACEExtrapolation::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
  //cutoff from the basis set's radial functions settings
  scale[j][i] = scale[i][j];
  return aceimpl->basis_set->radial_functions->cut(map[i], map[j]);
}

/* ----------------------------------------------------------------------
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairPACEExtrapolation::extract(const char *str, int &dim)
{
  dim = 0;
  //check if str=="gamma_flag" then compute extrapolation grades on this iteration
  if (strcmp(str, "gamma_flag") == 0) return (void *) &flag_compute_extrapolation_grade;
  //check if str=="corerep_flag" then compute extrapolation grades on this iteration
  if (strcmp(str, "corerep_flag") == 0) return (void *) &flag_corerep_factor;

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
void *PairPACEExtrapolation::extract_peratom(const char *str, int &ncol)
{
  if (strcmp(str, "gamma") == 0) {
    ncol = 0;
    return (void *) extrapolation_grade_gamma;
  }

  if (strcmp(str, "corerep") == 0) {
    ncol = 0;
    return (void *) corerep_factor;
  }

  return nullptr;
}
