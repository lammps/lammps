/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_pace.h"

#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_types.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <exception>

namespace LAMMPS_NS {
struct ACECimpl {
  ACECimpl() : basis_set(nullptr), ace(nullptr) {}
  ~ACECimpl()
  {
    delete basis_set;
    delete ace;
  }
  ACECTildeBasisSet *basis_set;
  ACECTildeEvaluator *ace;
};
}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;

template <int PERATOM>
ComputePACE<PERATOM>::ComputePACE(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), list(nullptr), pace(nullptr), paceall(nullptr), pace_peratom(nullptr),
    pace_atom(nullptr), c_pe(nullptr), c_virial(nullptr), acecimpl(nullptr)
{
  bikflag = 0;
  dgradflag = 0;

  // compute pace needs: file bikflag dgradflag, compute pace/atom needs: file

  int nargmin = PERATOM ? 4 : 6;
  if (narg < nargmin) utils::missing_cmd_args(FLERR, fmt::format("compute {}", style), error);
  if (narg > nargmin)
    error->all(FLERR, nargmin, "Unknown compute {} keyword: {}", style, arg[nargmin]);

  if (!PERATOM) {
    bikflag = utils::inumeric(FLERR, arg[4], false, lmp);
    dgradflag = utils::inumeric(FLERR, arg[5], false, lmp);
    if (bikflag != 0 && bikflag != 1)
      error->all(FLERR, 4, "Compute {} bikflag must be 0 or 1", style);
    if (dgradflag != 0 && dgradflag != 1)
      error->all(FLERR, 5, "Compute {} dgradflag must be 0 or 1", style);
    if (dgradflag && !bikflag)
      error->all(FLERR, 5, "Compute {} dgradflag = 1 requires bikflag = 1", style);
  }

  acecimpl = new ACECimpl;

  //read in file with CG coefficients or c_tilde coefficients

  auto potential_file_name = utils::get_potential_file_path(arg[3]);
  acecimpl->basis_set = new ACECTildeBasisSet(potential_file_name);
  cutmax = acecimpl->basis_set->cutoffmax;

  //# of rank 1, rank > 1 functions

  int n_r1 = acecimpl->basis_set->total_basis_size_rank1[0];
  int n_rp = acecimpl->basis_set->total_basis_size[0];
  nvalues = n_r1 + n_rp;

  // the descriptors of all elements must have the same size

  for (int mu = 1; mu < (int) acecimpl->basis_set->nelements; mu++)
    if (acecimpl->basis_set->total_basis_size_rank1[mu] != n_r1 ||
        acecimpl->basis_set->total_basis_size[mu] != n_rp)
      error->all(FLERR,
                 "Compute {}: the number of basis functions must be the same for all "
                 "elements",
                 style);

  ndims_force = 3;
  ndims_virial = 6;
  ndims_peratom = ndims_force;
  size_peratom = ndims_peratom * nvalues * atom->ntypes;
  yoffset = nvalues;
  zoffset = 2 * nvalues;
  nmax = 0;
  natoms = 0;
  bik_rows = 1;
  lastcol = 0;

  if (PERATOM) {

    // per-atom array with one row per atom and one column per descriptor

    peratom_flag = 1;
    size_peratom_cols = nvalues;
    return;
  }

  // global array of descriptors, forces, and virial.
  // rows are indexed by atom ID, so the size must fit in an int.

  array_flag = 1;
  extarray = 0;

  const bigint n = atom->natoms;
  if (bikflag) bik_rows = n;
  bigint nrows, ncols;
  if (dgradflag) {
    nrows = bik_rows + 3 * n * n + 1;
    ncols = nvalues + 3;
    if (comm->me == 0)
      error->warning(FLERR, "dgradflag=1 creates a N^2 array, beware of large systems.");
  } else {
    nrows = bik_rows + ndims_force * n + ndims_virial;
    ncols = (bigint) nvalues * atom->ntypes + 1;
  }
  if (nrows * ncols > MAXSMALLINT)
    error->all(FLERR, "Compute {} global array is too large for {} atoms", style, n);
  natoms = n;
  size_array_rows = nrows;
  size_array_cols = ncols;
  lastcol = size_array_cols - 1;
}

/* ---------------------------------------------------------------------- */

template <int PERATOM> ComputePACE<PERATOM>::~ComputePACE()
{
  // the compute for the reference virial tensor exists only after init()

  if (modify->get_compute_by_id(id_virial)) modify->delete_compute(id_virial);

  delete acecimpl;
  memory->destroy(pace);
  memory->destroy(paceall);
  memory->destroy(pace_peratom);
  memory->destroy(pace_atom);
}

/* ---------------------------------------------------------------------- */

template <int PERATOM> void ComputePACE<PERATOM>::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Compute {} requires a pair style be defined", style);

  if (cutmax > force->pair->cutforce)
    error->all(FLERR, "Compute {} cutoff is longer than pairwise cutoff", style);

  // the global array was sized at definition and is indexed by atom ID

  if (!PERATOM && atom->natoms != natoms)
    error->all(FLERR, "Compute {} was defined for {} atoms but there are now {}", style, natoms,
               atom->natoms);
  if (!PERATOM && (atom->natoms > 0) && !atom->tag_consecutive())
    error->all(FLERR, "Compute {} requires consecutive atom IDs", style);

  // atom types are mapped to elements by order

  const int ntypes = atom->ntypes;
  if (ntypes > (int) acecimpl->basis_set->nelements)
    error->all(FLERR,
               "Compute {}: number of atom types ({}) exceeds number of elements in "
               "potential file ({})",
               style, ntypes, acecimpl->basis_set->nelements);

  // need an occasional full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style(fmt::format("^{}$", style)).size() > 1 && comm->me == 0)
    error->warning(FLERR, "More than one compute {}", style);

  // create the evaluator, it is reused for all atoms

  delete acecimpl->ace;
  acecimpl->ace = new ACECTildeEvaluator(*acecimpl->basis_set);
  acecimpl->ace->compute_projections = true;
  acecimpl->ace->compute_b_grad = !PERATOM;
  acecimpl->ace->element_type_mapping.init(ntypes + 1);
  for (int ik = 1; ik <= ntypes; ik++) acecimpl->ace->element_type_mapping(ik) = ik - 1;

  if (PERATOM) return;

  // allocate memory for global array
  memory->destroy(pace);
  memory->destroy(paceall);
  memory->create(pace, size_array_rows, size_array_cols, "pace:pace");
  memory->create(paceall, size_array_rows, size_array_cols, "pace:paceall");
  array = paceall;

  // find compute for reference energy

  c_pe = modify->get_compute_by_id("thermo_pe");
  if (!c_pe) error->all(FLERR, "Compute thermo_pe does not exist.");

  // add compute for reference virial tensor, unless it exists from a previous run

  id_virial = id + std::string("_press") + std::to_string(instance_me);
  c_virial = modify->get_compute_by_id(id_virial);
  if (!c_virial) c_virial = modify->add_compute(id_virial + " all pressure NULL virial");
}

/* ---------------------------------------------------------------------- */

template <int PERATOM> void ComputePACE<PERATOM>::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

template <int PERATOM> void ComputePACE<PERATOM>::compute_array()
{
  int ntotal = atom->nlocal + atom->nghost;
  invoked_array = update->ntimestep;

  // atoms may have been created or deleted since init(), e.g. by fix deposit

  if (atom->natoms != natoms)
    error->all(FLERR, "Compute {} was defined for {} atoms but there are now {}", style, natoms,
               atom->natoms);

  // grow pace_peratom array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(pace_peratom);
    nmax = atom->nmax;
    memory->create(pace_peratom, nmax, size_peratom, "pace:pace_peratom");
  }

  // clear global array

  for (int irow = 0; irow < size_array_rows; irow++) {
    for (int icoeff = 0; icoeff < size_array_cols; icoeff++) { pace[irow][icoeff] = 0.0; }
  }

  // clear local peratom array

  for (int i = 0; i < ntotal; i++) {
    for (int icoeff = 0; icoeff < size_peratom; icoeff++) { pace_peratom[i][icoeff] = 0.0; }
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  int **const firstneigh = list->firstneigh;
  int *const type = atom->type;

  //determine the maximum number of neighbours
  int max_jnum = 0;
  for (int iitmp = 0; iitmp < inum; iitmp++) max_jnum = MAX(max_jnum, numneigh[ilist[iitmp]]);

  acecimpl->ace->resize_neighbours_cache(max_jnum);

  // compute pace derivatives for each atom in group
  // use full neighbor list to count atoms less than cutoff

  const int *const mask = atom->mask;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      const int ti = atom->tag[i] - 1;
      const int irow = bikflag ? ti : 0;
      const int itype = type[i];
      const int *const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      const int typeoffset_local = ndims_peratom * nvalues * (itype - 1);
      const int typeoffset_global = nvalues * (itype - 1);

      if (dgradflag) {

        // index tags (i, j, component) of the dBi/dRj rows

        for (int j = 0; j < natoms; j++) {
          for (int d = 0; d < 3; d++) {
            pace[bik_rows + (j * 3 * natoms) + 3 * ti + d][0] = ti;
            pace[bik_rows + (j * 3 * natoms) + 3 * ti + d][1] = j;
            pace[bik_rows + (j * 3 * natoms) + 3 * ti + d][2] = d;
          }
        }
      }

      try {
        acecimpl->ace->compute_atom(i, atom->x, atom->type, numneigh[i], firstneigh[i]);
      } catch (std::exception &e) {
        error->one(FLERR, e.what());
      }
      const auto &Bs = acecimpl->ace->projections;

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        if (!dgradflag) {
          double *pacedi = pace_peratom[i] + typeoffset_local;
          double *pacedj = pace_peratom[j] + typeoffset_local;

          // neighbours_dB(func_ind, jj, xyz) = dB_{i,func_ind}/dr_j for neighbor jj of atom i

          for (int func_ind = 0; func_ind < nvalues; func_ind++) {
            DOUBLE_TYPE fx_dB = acecimpl->ace->neighbours_dB(func_ind, jj, 0);
            DOUBLE_TYPE fy_dB = acecimpl->ace->neighbours_dB(func_ind, jj, 1);
            DOUBLE_TYPE fz_dB = acecimpl->ace->neighbours_dB(func_ind, jj, 2);
            pacedi[func_ind] += fx_dB;
            pacedi[func_ind + yoffset] += fy_dB;
            pacedi[func_ind + zoffset] += fz_dB;
            pacedj[func_ind] -= fx_dB;
            pacedj[func_ind + yoffset] -= fy_dB;
            pacedj[func_ind + zoffset] -= fz_dB;
          }
        } else {
          const int tj = atom->tag[j] - 1;
          for (int iicoeff = 0; iicoeff < nvalues; iicoeff++) {

            // add to pace array for this proc
            // dBi/dRj
            DOUBLE_TYPE fx_dB = acecimpl->ace->neighbours_dB(iicoeff, jj, 0);
            DOUBLE_TYPE fy_dB = acecimpl->ace->neighbours_dB(iicoeff, jj, 1);
            DOUBLE_TYPE fz_dB = acecimpl->ace->neighbours_dB(iicoeff, jj, 2);
            pace[bik_rows + (tj * 3 * natoms) + 3 * ti + 0][iicoeff + 3] -= fx_dB;
            pace[bik_rows + (tj * 3 * natoms) + 3 * ti + 1][iicoeff + 3] -= fy_dB;
            pace[bik_rows + (tj * 3 * natoms) + 3 * ti + 2][iicoeff + 3] -= fz_dB;

            // dBi/dRi
            pace[bik_rows + (ti * 3 * natoms) + 3 * ti + 0][iicoeff + 3] += fx_dB;
            pace[bik_rows + (ti * 3 * natoms) + 3 * ti + 1][iicoeff + 3] += fy_dB;
            pace[bik_rows + (ti * 3 * natoms) + 3 * ti + 2][iicoeff + 3] += fz_dB;
          }
        }
      }    // loop over jj inside
      if (!dgradflag) {

        int k = typeoffset_global;

        for (int icoeff = 0; icoeff < nvalues; icoeff++) { pace[irow][k++] += Bs(icoeff); }
      } else {
        int k = 3;
        for (int icoeff = 0; icoeff < nvalues; icoeff++) { pace[irow][k++] += Bs(icoeff); }
      }
    }    //group bit
  }      // for ii loop
  // accumulate force contributions to global array
  if (!dgradflag) {
    for (int itype = 0; itype < atom->ntypes; itype++) {
      const int typeoffset_local = ndims_peratom * nvalues * itype;
      const int typeoffset_global = nvalues * itype;
      for (int icoeff = 0; icoeff < nvalues; icoeff++) {
        for (int i = 0; i < ntotal; i++) {
          double *pacedi = pace_peratom[i] + typeoffset_local;
          int iglobal = atom->tag[i];
          int irow = 3 * (iglobal - 1) + bik_rows;
          pace[irow++][icoeff + typeoffset_global] += pacedi[icoeff];
          pace[irow++][icoeff + typeoffset_global] += pacedi[icoeff + yoffset];
          pace[irow][icoeff + typeoffset_global] += pacedi[icoeff + zoffset];
        }
      }
    }
  }

  if (!dgradflag) {
    // accumulate forces to global array
    for (int i = 0; i < atom->nlocal; i++) {
      int iglobal = atom->tag[i];
      int irow = 3 * (iglobal - 1) + bik_rows;
      pace[irow++][lastcol] = atom->f[i][0];
      pace[irow++][lastcol] = atom->f[i][1];
      pace[irow][lastcol] = atom->f[i][2];
    }
  } else {

    // for dgradflag=1, put forces at first 3 columns of bik rows

    for (int i = 0; i < atom->nlocal; i++) {
      int iglobal = atom->tag[i];
      pace[iglobal - 1][0 + 0] = atom->f[i][0];
      pace[iglobal - 1][0 + 1] = atom->f[i][1];
      pace[iglobal - 1][0 + 2] = atom->f[i][2];
    }
  }

  dbdotr_compute();

  // sum up over all processes
  MPI_Allreduce(&pace[0][0], &paceall[0][0], size_array_rows * size_array_cols, MPI_DOUBLE, MPI_SUM,
                world);

  // assign energy to last column

  if (!dgradflag) {
    for (int i = 0; i < bik_rows; i++) paceall[i][lastcol] = 0;
    int irow = 0;
    double reference_energy = c_pe->compute_scalar();
    paceall[irow][lastcol] = reference_energy;
  } else {

    // assign reference energy right after the dgrad rows, first column

    int irow = bik_rows + 3 * natoms * natoms;
    double reference_energy = c_pe->compute_scalar();
    paceall[irow][0] = reference_energy;
  }

  // assign virial stress to last column
  // switch to Voigt notation

  if (!dgradflag) {
    c_virial->compute_vector();
    int irow = 3 * natoms + bik_rows;
    paceall[irow++][lastcol] = c_virial->vector[0];
    paceall[irow++][lastcol] = c_virial->vector[1];
    paceall[irow++][lastcol] = c_virial->vector[2];
    paceall[irow++][lastcol] = c_virial->vector[5];
    paceall[irow++][lastcol] = c_virial->vector[4];
    paceall[irow++][lastcol] = c_virial->vector[3];
  }
}

/* ----------------------------------------------------------------------
   compute the ACE descriptors of each atom in the group
   rows of atoms outside the group are zero
------------------------------------------------------------------------- */

template <int PERATOM> void ComputePACE<PERATOM>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow pace_atom array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(pace_atom);
    nmax = atom->nmax;
    memory->create(pace_atom, nmax, size_peratom_cols, "pace/atom:pace_atom");
    array_atom = pace_atom;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  int **const firstneigh = list->firstneigh;

  //determine the maximum number of neighbours
  int max_jnum = 0;
  for (int iitmp = 0; iitmp < inum; iitmp++) max_jnum = MAX(max_jnum, numneigh[ilist[iitmp]]);

  acecimpl->ace->resize_neighbours_cache(max_jnum);

  // zero all local atoms: atoms outside the group or not in the neighbor list
  // (e.g. excluded with neigh_modify) must report zeros

  for (int i = 0; i < atom->nlocal; i++)
    for (int icoeff = 0; icoeff < nvalues; icoeff++) pace_atom[i][icoeff] = 0.0;

  // compute descriptors for each atom in group
  // use full neighbor list to count atoms less than cutoff

  const int *const mask = atom->mask;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    try {
      acecimpl->ace->compute_atom(i, atom->x, atom->type, numneigh[i], firstneigh[i]);
    } catch (std::exception &e) {
      error->one(FLERR, e.what());
    }
    const auto &Bs = acecimpl->ace->projections;
    for (int icoeff = 0; icoeff < nvalues; icoeff++) pace_atom[i][icoeff] = Bs(icoeff);
  }
}

/* ----------------------------------------------------------------------
   compute global virial contributions via summing r_i.dB^j/dr_i over
   own & ghost atoms
------------------------------------------------------------------------- */
template <int PERATOM> void ComputePACE<PERATOM>::dbdotr_compute()
{

  if (dgradflag) return;

  double **x = atom->x;
  int irow0 = bik_rows + ndims_force * natoms;

  // sum over ace contributions to forces
  // on all particles including ghosts

  int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; i++)
    for (int itype = 0; itype < atom->ntypes; itype++) {
      const int typeoffset_local = ndims_peratom * nvalues * itype;
      const int typeoffset_global = nvalues * itype;
      double *pacedi = pace_peratom[i] + typeoffset_local;
      for (int icoeff = 0; icoeff < nvalues; icoeff++) {
        double dbdx = pacedi[icoeff];
        double dbdy = pacedi[icoeff + yoffset];
        double dbdz = pacedi[icoeff + zoffset];
        int irow = irow0;
        pace[irow++][icoeff + typeoffset_global] += dbdx * x[i][0];
        pace[irow++][icoeff + typeoffset_global] += dbdy * x[i][1];
        pace[irow++][icoeff + typeoffset_global] += dbdz * x[i][2];
        pace[irow++][icoeff + typeoffset_global] += dbdz * x[i][1];
        pace[irow++][icoeff + typeoffset_global] += dbdz * x[i][0];
        pace[irow++][icoeff + typeoffset_global] += dbdy * x[i][0];
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

template <int PERATOM> double ComputePACE<PERATOM>::memory_usage()
{
  double bytes = 0.0;
  if (PERATOM) {
    bytes += (double) nmax * nvalues * sizeof(double);    // pace_atom
  } else {
    bytes += (double) size_array_rows * size_array_cols * sizeof(double);    // pace
    bytes += (double) size_array_rows * size_array_cols * sizeof(double);    // paceall
    bytes += (double) nmax * size_peratom * sizeof(double);                  // pace_peratom
  }
  return bytes;
}

namespace LAMMPS_NS {
template class ComputePACE<0>;
template class ComputePACE<1>;
}    // namespace LAMMPS_NS
