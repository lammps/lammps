// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------
   BOCS written by: Nicholas J. H. Dunn and Michael R. DeLyser
   from The Pennsylvania State University
------------------------------------------------------------------------- */

#include "compute_pressure_bocs.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePressureBocs::ComputePressureBocs(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  vptr(nullptr), id_temp(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR,"compute pressure/bocs", error);
  if (igroup) error->all(FLERR,"Compute pressure/bocs must use group all");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 0;
  pressflag = 1;
  timeflag = 1;

  p_match_flag = 0;
  phi_coeff = nullptr;

  // store temperature ID used by pressure computation
  // ensure it is valid for temperature computation

  if (strcmp(arg[3],"NULL") == 0) id_temp = nullptr;
  else {
    id_temp = utils::strdup(arg[3]);

    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature)
      error->all(FLERR,"Could not find compute pressure/bocs temperature compute {}", id_temp);
    if (temperature->tempflag == 0)
      error->all(FLERR,"Compute pressure/bocs temperature compute {} does not compute "
                 "temperature", id_temp);
  }

  // process optional args

  if (narg == 4) {
    keflag = 1;
    pairflag = 1;
    bondflag = angleflag = dihedralflag = improperflag = 1;
    kspaceflag = fixflag = 1;
  } else {
    keflag = 0;
    pairflag = 0;
    bondflag = angleflag = dihedralflag = improperflag = 0;
    kspaceflag = fixflag = 0;
    int iarg = 4;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"ke") == 0) keflag = 1;
      else if (strcmp(arg[iarg],"pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg],"bond") == 0) bondflag = 1;
      else if (strcmp(arg[iarg],"angle") == 0) angleflag = 1;
      else if (strcmp(arg[iarg],"dihedral") == 0) dihedralflag = 1;
      else if (strcmp(arg[iarg],"improper") == 0) improperflag = 1;
      else if (strcmp(arg[iarg],"kspace") == 0) kspaceflag = 1;
      else if (strcmp(arg[iarg],"fix") == 0) fixflag = 1;
      else if (strcmp(arg[iarg],"virial") == 0) {
        pairflag = 1;
        bondflag = angleflag = dihedralflag = improperflag = 1;
        kspaceflag = fixflag = 1;
      } else error->all(FLERR,"Illegal compute pressure/bocs command");
      iarg++;
    }
  }

  // error check

  if (keflag && id_temp == nullptr)
    error->all(FLERR,"Compute pressure/bocs requires temperature ID "
               "to include kinetic energy");

  vector = new double[size_vector];
  nvirial = 0;
  vptr = nullptr;

  splines = nullptr;
  spline_length = 0;
}

/* ---------------------------------------------------------------------- */

ComputePressureBocs::~ComputePressureBocs()
{
  delete [] id_temp;
  delete [] vector;
  delete [] vptr;
  if (phi_coeff) free(phi_coeff);
}

/* ---------------------------------------------------------------------- */

void ComputePressureBocs::init()
{
  boltz = force->boltz;
  nktv2p = force->nktv2p;
  dimension = domain->dimension;

  // set temperature compute, must be done in init()
  // fixes could have changed or compute_modify could have changed it

  if (keflag) {
    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature)
      error->all(FLERR,"Could not find compute pressure/bocs temperature compute {}", id_temp);
  }

  // detect contributions to virial
  // vptr points to all virial[6] contributions

  delete [] vptr;
  nvirial = 0;
  vptr = nullptr;

  if (pairflag && force->pair) nvirial++;
  if (atom->molecular != Atom::ATOMIC) {
    if (bondflag && force->bond) nvirial++;
    if (angleflag && force->angle) nvirial++;
    if (dihedralflag && force->dihedral) nvirial++;
    if (improperflag && force->improper) nvirial++;
  }
  if (fixflag) {
    for (const auto &ifix : modify->get_fix_list())
      if (ifix->thermo_virial) nvirial++;
  }

  if (nvirial) {
    vptr = new double*[nvirial];
    nvirial = 0;
    if (pairflag && force->pair) vptr[nvirial++] = force->pair->virial;
    if (bondflag && force->bond) vptr[nvirial++] = force->bond->virial;
    if (angleflag && force->angle) vptr[nvirial++] = force->angle->virial;
    if (dihedralflag && force->dihedral)
      vptr[nvirial++] = force->dihedral->virial;
    if (improperflag && force->improper)
      vptr[nvirial++] = force->improper->virial;
    if (fixflag) {
      for (const auto &ifix : modify->get_fix_list())
        if (ifix->virial_global_flag && ifix->thermo_virial)
          vptr[nvirial++] = ifix->virial;
    }
  }

  // flag Kspace contribution separately, since not summed across procs

  if (kspaceflag && force->kspace) kspace_virial = force->kspace->virial;
  else kspace_virial = nullptr;
}

/* ----------------------------------------------------------------------
   Compute the pressure correction for the analytical basis set
------------------------------------------------------------------------- */

double ComputePressureBocs::get_cg_p_corr(int N_basis, double *phi_coeff,
                                      int N_mol, double vavg, double vCG)
{
  double correction = 0.0;
  for (int i = 1; i <= N_basis; ++i)
    correction -= phi_coeff[i-1] * ( N_mol * i / vavg ) *
      pow( ( 1 / vavg ) * ( vCG - vavg ),i-1);
  return correction;
}

/* ----------------------------------------------------------------------
   Find the relevant index position if using a spline basis set
------------------------------------------------------------------------- */

double ComputePressureBocs::find_index(double * grid, double value)
{
  int i;
  double spacing = fabs(grid[1]-grid[0]);
  int gridsize = spline_length;
  for (i = 0; i < (gridsize-1); ++i)
  {
    if (value >= grid[i] && value <= grid[i+1]) { return i; }
  }

  if (value >= grid[i] && value <= (grid[i] + spacing)) { return i; }

  error->all(FLERR,"find_index could not find value in grid for value: {}", value);
  for (int i = 0; i < gridsize; ++i)
  {
    fprintf(stderr, "grid %d: %f\n",i,grid[i]);
  }

  exit(1);
}

/* ----------------------------------------------------------------------
   Compute the pressure correction for a spline basis set
------------------------------------------------------------------------- */

double ComputePressureBocs::get_cg_p_corr(double ** grid, int basis_type,
                                          double vCG)
{
  int i = find_index(grid[0],vCG);
  double deltax = vCG - grid[0][i];

  if (basis_type == BASIS_LINEAR_SPLINE)
    return grid[1][i] + (deltax) * ( grid[1][i+1] - grid[1][i] ) / ( grid[0][i+1] - grid[0][i] );
  else if (basis_type == BASIS_CUBIC_SPLINE)
    return grid[1][i] + (grid[2][i] * deltax) + (grid[3][i] * pow(deltax,2)) + (grid[4][i] * pow(deltax,3));
  else error->all(FLERR,"bad spline type passed to get_cg_p_corr()\n");
  return 0.0;
}

/* ----------------------------------------------------------------------
   send cg info from fix_bocs to compute_pressure_bocs for the analytical
   basis set
------------------------------------------------------------------------- */

void ComputePressureBocs::send_cg_info(int basis_type, int sent_N_basis,
                                       double *sent_phi_coeff, int sent_N_mol,
                                       double sent_vavg)
{
  if (basis_type == BASIS_ANALYTIC) p_basis_type = BASIS_ANALYTIC;
  else error->all(FLERR,"Incorrect basis type passed to ComputePressureBocs\n");

  p_match_flag = 1;

  N_basis = sent_N_basis;
  if (phi_coeff) free(phi_coeff);
  phi_coeff = ((double *) calloc(N_basis, sizeof(double)) );
  for (int i=0; i<N_basis; i++) { phi_coeff[i] = sent_phi_coeff[i]; }

  N_mol = sent_N_mol;
  vavg = sent_vavg;
}

/* ----------------------------------------------------------------------
   send cg info from fix_bocs to compute_pressure_bocs for a spline basis
   set
------------------------------------------------------------------------- */

void ComputePressureBocs::send_cg_info(int basis_type,
                                       double ** in_splines, int gridsize)
{
  if (basis_type == BASIS_LINEAR_SPLINE) { p_basis_type = BASIS_LINEAR_SPLINE; }
  else if (basis_type == BASIS_CUBIC_SPLINE) { p_basis_type = BASIS_CUBIC_SPLINE; }
  else
  {
    error->all(FLERR,"Incorrect basis type passed to ComputePressureBocs\n");
  }
  splines = in_splines;
  spline_length = gridsize;
  p_match_flag = 1;
}

/* End of new functions for BOCS */

/* ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
------------------------------------------------------------------------- */

double ComputePressureBocs::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->vflag_global != invoked_scalar)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  // invoke temperature if it hasn't been already

  double t;
  double volume, correction = 0;
  if (keflag) {
    if (temperature->invoked_scalar != update->ntimestep)
      t = temperature->compute_scalar();
    else t = temperature->scalar;
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    volume = (domain->xprd * domain->yprd * domain->zprd);

    /* MRD NJD if block */
    if (p_basis_type == BASIS_ANALYTIC)
    {
      correction = get_cg_p_corr(N_basis,phi_coeff,N_mol,vavg,volume);
    }
    else if (p_basis_type == BASIS_LINEAR_SPLINE || p_basis_type == BASIS_CUBIC_SPLINE)
    {
      correction = get_cg_p_corr(splines, p_basis_type, volume);
    }

    virial_compute(3,3);
    if (keflag)
      scalar = (temperature->dof * boltz * t +
                virial[0] + virial[1] + virial[2]) / 3.0 *
                inv_volume * nktv2p + (correction);
    else
      scalar = (virial[0] + virial[1] + virial[2]) / 3.0 *
               inv_volume * nktv2p + (correction);
  } else {
    if (p_match_flag)
    {
      error->all(FLERR,"Pressure matching not implemented in 2-d.\n");
      exit(1);
    } // The rest of this can probably be deleted.
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(2,2);
    if (keflag)
      scalar = (temperature->dof * boltz * t +
                virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
    else
      scalar = (virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
  }

  return scalar;
}

/* ----------------------------------------------------------------------
   compute pressure tensor
   assume KE tensor has already been computed
------------------------------------------------------------------------- */

void ComputePressureBocs::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->vflag_global != invoked_vector)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  if (force->kspace && kspace_virial && force->kspace->scalar_pressure_flag)
    error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' for "
               "tensor components with kspace_style msm");

  // invoke temperature if it hasn't been already

  double *ke_tensor;
  if (keflag) {
    if (temperature->invoked_vector != update->ntimestep)
      temperature->compute_vector();
    ke_tensor = temperature->vector;
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(6,3);
    if (keflag) {
      for (int i = 0; i < 6; i++)
        vector[i] = (ke_tensor[i] + virial[i]) * inv_volume * nktv2p;
    } else
      for (int i = 0; i < 6; i++)
        vector[i] = virial[i] * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(4,2);
    if (keflag) {
      vector[0] = (ke_tensor[0] + virial[0]) * inv_volume * nktv2p;
      vector[1] = (ke_tensor[1] + virial[1]) * inv_volume * nktv2p;
      vector[3] = (ke_tensor[3] + virial[3]) * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    } else {
      vector[0] = virial[0] * inv_volume * nktv2p;
      vector[1] = virial[1] * inv_volume * nktv2p;
      vector[3] = virial[3] * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePressureBocs::virial_compute(int n, int ndiag)
{
  int i,j;
  double v[6],*vcomponent;

  for (i = 0; i < n; i++) v[i] = 0.0;

  // sum contributions to virial from forces and fixes

  for (j = 0; j < nvirial; j++) {
    vcomponent = vptr[j];
    for (i = 0; i < n; i++) v[i] += vcomponent[i];
  }

  // sum virial across procs

  MPI_Allreduce(v,virial,n,MPI_DOUBLE,MPI_SUM,world);

  // KSpace virial contribution is already summed across procs

  if (kspace_virial)
    for (i = 0; i < n; i++) virial[i] += kspace_virial[i];

  // LJ long-range tail correction, only if pair contributions are included

  if (force->pair && pairflag && force->pair->tail_flag)
    for (i = 0; i < ndiag; i++) virial[i] += force->pair->ptail * inv_volume;
}

/* ---------------------------------------------------------------------- */

void ComputePressureBocs::reset_extra_compute_fix(const char *id_new)
{
  delete [] id_temp;
  id_temp = utils::strdup(id_new);
}
