// clang-format off
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

#include <cstdio>


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

enum { SCALAR, VECTOR, ARRAY };
ComputePACE::ComputePACE(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), list(nullptr), pace(nullptr), paceall(nullptr),
    c_pe(nullptr), c_virial(nullptr), acecimpl(nullptr)
{
  array_flag = 1;
  extarray = 0;
  bikflag = 0;
  dgradflag = 0;

  const int ntypes = atom->ntypes;
  const int nargmin = 4;

  acecimpl = new ACECimpl;
  if (narg < nargmin) error->all(FLERR,"Illegal compute pace command");

  bikflag = utils::inumeric(FLERR, arg[4], false, lmp);
  dgradflag = utils::inumeric(FLERR, arg[5], false, lmp);
  if (dgradflag && !bikflag)
    error->all(FLERR,"Illegal compute pace command: dgradflag=1 requires bikflag=1");

  //read in file with CG coefficients or c_tilde coefficients
  auto potential_file_name = utils::get_potential_file_path(arg[3]);
  delete acecimpl->basis_set;
  acecimpl->basis_set = new ACECTildeBasisSet(potential_file_name);
  cutmax = acecimpl->basis_set->cutoffmax;
  delete acecimpl->ace;
  acecimpl->ace = new ACECTildeEvaluator(*acecimpl->basis_set);
  acecimpl->ace->compute_projections = true;
  acecimpl->ace->compute_b_grad = true;
  acecimpl->ace->element_type_mapping.init(ntypes+1);
  for (int ik = 1; ik <= ntypes; ik++) {
    for(int mu = 0; mu < acecimpl->basis_set->nelements; mu++){
      if (mu == ik - 1) acecimpl->ace->element_type_mapping(ik) = mu;
    }
  }

  // ACECTildeEvaluator::get_func_ind_shift() has a bug so instead lets do it manually here
  ncoeff = 0;
  number_of_functions = std::vector<int>(ntypes+1, 0);
  type_offsets = std::vector<int>(ntypes+1, 0);
  for( int i=1 ; i<=ntypes ; i++ ) {
    type_offsets.at(i) = ncoeff;
    const SPECIES_TYPE mu = acecimpl->ace->element_type_mapping(i);
    number_of_functions.at(i) = acecimpl->basis_set->total_basis_size_rank1[mu];
    number_of_functions.at(i) += acecimpl->basis_set->total_basis_size[mu];
    ncoeff += number_of_functions.at(i);
  }

  ndims_force = 3;
  ndims_virial = 6;
  bik_rows = 1;
  natoms = atom->natoms;
  if (bikflag) bik_rows = natoms;
  dgrad_rows = ndims_force*natoms;
  size_array_rows = bik_rows + dgrad_rows + ndims_virial;
  if (dgradflag) {
    size_array_rows = bik_rows + 3*natoms*natoms + 1;
    size_array_cols = ncoeff + 3;
    if (comm->me == 0) error->warning(FLERR,"dgradflag=1 creates a N^2 array, beware of large systems.");
  } else size_array_cols = ncoeff + 1;
  lastcol = size_array_cols-1;

}

/* ---------------------------------------------------------------------- */

ComputePACE::~ComputePACE()
{
  if( modify->find_compute(id_virial) != -1 ) modify->delete_compute(id_virial);

  delete acecimpl;
  memory->destroy(pace);
  memory->destroy(paceall);
}

/* ---------------------------------------------------------------------- */

void ComputePACE::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute pace requires a pair style be defined");

  if (cutmax > force->pair->cutforce)
    error->all(FLERR,"Compute pace cutoff {} is longer than pairwise cutoff {}", cutmax, force->pair->cutforce);

  // need an occasional full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style("pace").size() > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute pace");

  // allocate memory for global array
  memory->create(pace,size_array_rows,size_array_cols, "pace:pace");
  memory->create(paceall,size_array_rows,size_array_cols, "pace:paceall");
  array = paceall;

  // find compute for reference energy

  c_pe = modify->get_compute_by_id("thermo_pe");
  if (!c_pe) error->all(FLERR,"Compute thermo_pe does not exist.");

  // add compute for reference virial tensor

  id_virial = id + std::string("_press");
  c_virial = modify->add_compute(id_virial + " all pressure NULL virial");
}

/* ---------------------------------------------------------------------- */

void ComputePACE::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePACE::compute_array()
{
  int ntotal = atom->nlocal + atom->nghost;
  invoked_array = update->ntimestep;

  // clear global array
  for (int irow = 0; irow < size_array_rows; irow++){
    for (int icoeff = 0; icoeff < size_array_cols; icoeff++) pace[irow][icoeff] = 0.0;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int* const ilist = list->ilist;
  const int* const numneigh = list->numneigh;
  int** const firstneigh = list->firstneigh;
  int * const type = atom->type;
  double **x = atom->x;

  //determine the maximum number of neighbours
  int max_jnum = -1;
  int nei = 0;
  int jtmp =0;
  for (int iitmp = 0; iitmp < list->inum; iitmp++) {
    int itmp = ilist[iitmp];
    jtmp = numneigh[itmp];
    nei = nei + jtmp;
    if (jtmp > max_jnum) max_jnum = jtmp;
  }

  // compute pace derivatives for each atom in group
  // use full neighbor list to count atoms less than cutoff

  const int* const mask = atom->mask;
  const int ntypes = atom->ntypes;

  for (int ii = 0; ii < inum; ii++) {
    int irow = 0;
    if (bikflag) irow = atom->tag[ilist[ii] & NEIGHMASK]-1;
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      const int itype = type[i];
      const int* const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      const int row_offset_i = bik_rows + 3*(atom->tag[i]-1);
      const int type_offset = type_offsets.at(itype);

      if (dgradflag) {
        // dBi/dRi tags
        pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 0][0] = atom->tag[i]-1;
        pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 0][1] = atom->tag[i]-1;
        pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 0][2] = 0;
        pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 1][0] = atom->tag[i]-1;
        pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 1][1] = atom->tag[i]-1;
        pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 1][2] = 1;
        pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 2][0] = atom->tag[i]-1;
        pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 2][1] = atom->tag[i]-1;
        pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 2][2] = 2;
        // dBi/dRj tags
        for (int j=0; j<natoms; j++) {
          pace[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 0][0] = atom->tag[i]-1;
          pace[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 0][1] = j;
          pace[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 0][2] = 0;
          pace[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 1][0] = atom->tag[i]-1;
          pace[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 1][1] = j;
          pace[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 1][2] = 1;
          pace[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 2][0] = atom->tag[i]-1;
          pace[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 2][1] = j;
          pace[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 2][2] = 2;
        }
      }

      // resize the neighbor cache after setting the basis
      acecimpl->ace->resize_neighbours_cache(max_jnum);
      acecimpl->ace->compute_atom(i, atom->x, atom->type, list->numneigh[i], list->firstneigh[i]);
      Array1D<DOUBLE_TYPE> Bs = acecimpl->ace->projections;

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj];
        const int row_offset_j = bik_rows + 3*(atom->tag[j]-1);
        for (int func_ind=0; func_ind < number_of_functions.at(itype); func_ind++){
          DOUBLE_TYPE fx_dB = acecimpl->ace->neighbours_dB(func_ind,jj,0);
          DOUBLE_TYPE fy_dB = acecimpl->ace->neighbours_dB(func_ind,jj,1);
          DOUBLE_TYPE fz_dB = acecimpl->ace->neighbours_dB(func_ind,jj,2);
          if (!dgradflag) {
            // forces
            pace[row_offset_i    ][type_offset + func_ind] += fx_dB;
            pace[row_offset_i + 1][type_offset + func_ind] += fy_dB;
            pace[row_offset_i + 2][type_offset + func_ind] += fz_dB;
            pace[row_offset_j    ][type_offset + func_ind] -= fx_dB;
            pace[row_offset_j + 1][type_offset + func_ind] -= fy_dB;
            pace[row_offset_j + 2][type_offset + func_ind] -= fz_dB;
            // virial
            pace[size_array_rows-6][type_offset + func_ind] += (fx_dB*x[i][0] - fx_dB*x[j][0]);
            pace[size_array_rows-5][type_offset + func_ind] += (fy_dB*x[i][1] - fy_dB*x[j][1]);
            pace[size_array_rows-4][type_offset + func_ind] += (fz_dB*x[i][2] - fz_dB*x[j][2]);
            pace[size_array_rows-3][type_offset + func_ind] += (fz_dB*x[i][1] - fz_dB*x[j][1]);
            pace[size_array_rows-2][type_offset + func_ind] += (fz_dB*x[i][0] - fz_dB*x[j][0]);
            pace[size_array_rows-1][type_offset + func_ind] += (fy_dB*x[i][0] - fy_dB*x[j][0]);
          } else {
              const int column_offset = 3 + type_offset + func_ind;
              // dBi/dRj
              pace[bik_rows + ((atom->tag[j]-1)*3*natoms) + 3*(atom->tag[i]-1) + 0][column_offset] -= fx_dB;
              pace[bik_rows + ((atom->tag[j]-1)*3*natoms) + 3*(atom->tag[i]-1) + 1][column_offset] -= fy_dB;
              pace[bik_rows + ((atom->tag[j]-1)*3*natoms) + 3*(atom->tag[i]-1) + 2][column_offset] -= fz_dB;
              // dBi/dRi
              pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 0][column_offset] += fx_dB;
              pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 1][column_offset] += fy_dB;
              pace[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 2][column_offset] += fz_dB;
          }
        }
      } // loop over jj inside
      for (int func_ind=0; func_ind < number_of_functions.at(itype); func_ind++) {
        if (!dgradflag)
          pace[irow][type_offset+func_ind] += Bs(func_ind);
        else
          pace[irow][3+type_offset+func_ind] += Bs(func_ind);
      }
    } //group bit
  } // for ii loop

  if (!dgradflag) {
    // accumulate forces to global array
    for (int i = 0; i < atom->nlocal; i++) {
      int iglobal = atom->tag[i];
      int irow = 3*(iglobal-1)+1;
      pace[irow++][lastcol] = atom->f[i][0];
      pace[irow++][lastcol] = atom->f[i][1];
      pace[irow][lastcol] = atom->f[i][2];
    }
  } else {
    // for dgradflag=1, put forces at first 3 columns of bik rows
    for (int i=0; i<atom->nlocal; i++) {
      int iglobal = atom->tag[i];
      pace[iglobal-1][0] = atom->f[i][0];
      pace[iglobal-1][1] = atom->f[i][1];
      pace[iglobal-1][2] = atom->f[i][2];
    }
  }

  // sum up over all processes
  MPI_Allreduce(&pace[0][0],&paceall[0][0],size_array_rows*size_array_cols,MPI_DOUBLE,MPI_SUM,world);

  // assign energy to last column
  if (!dgradflag) {
    for (int i = 0; i < bik_rows; i++) paceall[i][lastcol] = 0;
    int irow = 0;
    double reference_energy = c_pe->compute_scalar();
    paceall[irow][lastcol] = reference_energy;
  } else {
    // assign reference energy right after the dgrad rows, first column
    int irow = bik_rows + 3*natoms*natoms;
    double reference_energy = c_pe->compute_scalar();
    paceall[irow][0] = reference_energy;
  }

  // assign virial stress to last column
  // switch to Voigt notation
  if (!dgradflag) {
    c_virial->compute_vector();
    int irow = 3*natoms+bik_rows;
    paceall[irow++][lastcol] = c_virial->vector[0];
    paceall[irow++][lastcol] = c_virial->vector[1];
    paceall[irow++][lastcol] = c_virial->vector[2];
    paceall[irow++][lastcol] = c_virial->vector[5];
    paceall[irow++][lastcol] = c_virial->vector[4];
    paceall[irow++][lastcol] = c_virial->vector[3];
  }

}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputePACE::memory_usage()
{
  double bytes = (double)size_array_rows*size_array_cols*sizeof(double); // pace
  bytes += (double)size_array_rows*size_array_cols*sizeof(double);       // paceall
  return bytes;
}
