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
------------------------------------------------------------------------- */

#include "integrate.h"

#include "citeme.h"
#include "compute.h"
#include "force.h"
#include "kspace.h"
#include "modify.h"
#include "pair.h"
#include "output.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Integrate::Integrate(LAMMPS *lmp, int /*narg*/, char ** /*arg*/) : Pointers(lmp)
{
  elist_global = elist_atom = nullptr;
  vlist_global = vlist_atom = cvlist_atom = nullptr;
  external_force_clear = 0;
}

/* ---------------------------------------------------------------------- */

Integrate::~Integrate()
{
  delete [] elist_global;
  delete [] elist_atom;
  delete [] vlist_global;
  delete [] vlist_atom;
  delete [] cvlist_atom;
}

/* ---------------------------------------------------------------------- */

void Integrate::init()
{
  if (lmp->citeme) lmp->citeme->flush();
  update->atimestep = update->ntimestep;

  // allow pair and Kspace compute() to be turned off via modify flags

  if (force->pair && force->pair->compute_flag) pair_compute_flag = 1;
  else pair_compute_flag = 0;
  if (force->kspace && force->kspace->compute_flag) kspace_compute_flag = 1;
  else kspace_compute_flag = 0;

  // should add checks:
  // for any acceleration package that has its own integrate/minimize
  // in case input script has reset the run or minimize style explicitly
  // e.g. invalid to have kokkos pair style with non-kokkos verlet
  // but OK to have kokkos verlet with non kokkos pair style (just warn)
  // making these checks would require all the pair, fix, etc styles have
  //   kokkos, intel flags
}

/* ----------------------------------------------------------------------
   setup lists of computes for global and per-atom PE and pressure
------------------------------------------------------------------------- */

void Integrate::ev_setup()
{
  delete [] elist_global;
  delete [] elist_atom;
  delete [] vlist_global;
  delete [] vlist_atom;
  delete [] cvlist_atom;
  elist_global = elist_atom = nullptr;
  vlist_global = vlist_atom = cvlist_atom = nullptr;

  nelist_global = nelist_atom = 0;
  nvlist_global = nvlist_atom = ncvlist_atom = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->peflag) nelist_global++;
    if (modify->compute[i]->peatomflag) nelist_atom++;
    if (modify->compute[i]->pressflag) nvlist_global++;
    if (modify->compute[i]->pressatomflag & 1) nvlist_atom++;
    if (modify->compute[i]->pressatomflag & 2) ncvlist_atom++;
  }

  if (nelist_global) elist_global = new Compute*[nelist_global];
  if (nelist_atom) elist_atom = new Compute*[nelist_atom];
  if (nvlist_global) vlist_global = new Compute*[nvlist_global];
  if (nvlist_atom) vlist_atom = new Compute*[nvlist_atom];
  if (ncvlist_atom) cvlist_atom = new Compute*[ncvlist_atom];

  nelist_global = nelist_atom = 0;
  nvlist_global = nvlist_atom = ncvlist_atom = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->peflag)
      elist_global[nelist_global++] = modify->compute[i];
    if (modify->compute[i]->peatomflag)
      elist_atom[nelist_atom++] = modify->compute[i];
    if (modify->compute[i]->pressflag)
      vlist_global[nvlist_global++] = modify->compute[i];
    if (modify->compute[i]->pressatomflag & 1)
      vlist_atom[nvlist_atom++] = modify->compute[i];
    if (modify->compute[i]->pressatomflag & 2)
      cvlist_atom[ncvlist_atom++] = modify->compute[i];
  }
}

/* ----------------------------------------------------------------------
   set eflag,vflag for current iteration
   based on
     (1) computes that need energy/virial info on this timestep
     (2) time dumps that may need per-atom compute info on this timestep
     NOTE: inefficient to add all per-atom eng/virial computes
             but don't know which ones the dump needs
           see NOTE in output.cpp
   invoke matchstep() on all timestep-dependent computes to clear their arrays
   eflag: set any or no bits
     ENERGY_GLOBAL bit for global energy
     ENERGY_ATOM   bit for per-atom energy
   vflag: set any or no bits, but PAIR/FDOTR bits cannot both be set
     VIRIAL_PAIR     bit for global virial as sum of pairwise terms
     VIRIAL_FDOTR    bit for global virial via F dot r
     VIRIAL_ATOM     bit for per-atom virial
     VIRIAL_CENTROID bit for per-atom centroid virial
   all force components (pair,bond,angle,...,kspace) use eflag/vflag
     in their ev_setup() method to set local energy/virial flags
------------------------------------------------------------------------- */

void Integrate::ev_set(bigint ntimestep)
{
  int i,flag;

  int tdflag = 0;
  if (output->any_time_dumps &&
      output->next_time_dump_any == ntimestep) tdflag = 1;

  flag = 0;
  int eflag_global = 0;
  for (i = 0; i < nelist_global; i++)
    if (elist_global[i]->matchstep(ntimestep)) flag = 1;
  if (flag) eflag_global = ENERGY_GLOBAL;

  flag = 0;
  int eflag_atom = 0;
  for (i = 0; i < nelist_atom; i++)
    if (elist_atom[i]->matchstep(ntimestep)) flag = 1;
  if (flag || (tdflag && nelist_atom)) eflag_atom = ENERGY_ATOM;

  if (eflag_global) update->eflag_global = ntimestep;
  if (eflag_atom) update->eflag_atom = ntimestep;
  eflag = eflag_global + eflag_atom;

  flag = 0;
  int vflag_global = 0;
  for (i = 0; i < nvlist_global; i++)
    if (vlist_global[i]->matchstep(ntimestep)) flag = 1;
  if (flag) vflag_global = virial_style;

  flag = 0;
  int vflag_atom = 0;
  for (i = 0; i < nvlist_atom; i++)
    if (vlist_atom[i]->matchstep(ntimestep)) flag = 1;
  if (flag || (tdflag && nvlist_atom)) vflag_atom = VIRIAL_ATOM;

  flag = 0;
  int cvflag_atom = 0;
  for (i = 0; i < ncvlist_atom; i++)
    if (cvlist_atom[i]->matchstep(ntimestep)) flag = 1;
  if (flag || (tdflag && ncvlist_atom)) cvflag_atom = VIRIAL_CENTROID;

  if (vflag_global) update->vflag_global = ntimestep;
  if (vflag_atom || cvflag_atom) update->vflag_atom = ntimestep;
  vflag = vflag_global + vflag_atom + cvflag_atom;
}
