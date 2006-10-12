/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "pressure.h"
#include "output.h"
#include "thermo.h"
#include "temperature.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "fix_npt.h"
#include "fix_nph.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "pair.h"
#include "kspace.h"
#include "atom.h"
#include "domain.h"

/* ---------------------------------------------------------------------- */

void Pressure::init()
{
  boltz = force->boltz;
  nktv2p = force->nktv2p;

  // set tensor flag if p_tensor will ever be used
  // possibly used by thermo and fix_npt, fix_nph

  tensorflag = 0;
  if (output->thermo->tensorflag) tensorflag = 1;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"npt") == 0)
      if (((FixNPT *) modify->fix[i])->press_couple > 0) tensorflag = 1;
    if (strcmp(modify->fix[i]->style,"nph") == 0)
      if (((FixNPH *) modify->fix[i])->press_couple > 0) tensorflag = 1;
  }
    
  // set flags/ptrs for all contributions to virial
    
  pairflag = bondflag = angleflag = dihedralflag = improperflag = 0;
  kspaceflag = 0;
  shakeflag = bodyflag = rigidflag = poemsflag = 0;

  if (force->pair) {
    pairflag = 1;
    pair_virial = force->pair->virial;
  }
  if (atom->molecular) {
    if (force->bond) {
      bondflag = 1;
      bond_virial = force->bond->virial;
    }
    if (force->angle) {
      angleflag = 1;
      angle_virial = force->angle->virial;
    }
    if (force->dihedral) {
      dihedralflag = 1;
      dihedral_virial = force->dihedral->virial;
    }
    if (force->improper) {
      improperflag = 1;
      improper_virial = force->improper->virial;
    }
  }

  if (force->kspace) {
    kspaceflag = 1;
    kspace_virial = force->kspace->virial;
  }

  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"shake") == 0) {
      shakeflag = 1;
      shake_virial = modify->fix[i]->virial;
    }
    if (strcmp(modify->fix[i]->style,"rigid") == 0 ||
	strcmp(modify->fix[i]->style,"poems") == 0) {
      bodyflag = 1;
      if (strcmp(modify->fix[i]->style,"rigid") == 0) {
	rigidflag = 1;
	rigid_virial = modify->fix[i]->virial;
      } else {
	poemsflag = 1;
	poems_virial = modify->fix[i]->virial;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void Pressure::compute(Temperature *temperature)
{
  int i,n;
  double v[6];

  if (tensorflag) n = 6;
  else n = 3;
  for (i = 0; i < n; i++) v[i] = 0.0;

  // sum contributions to virial from various forces and fixes

  if (pairflag)
    for (i = 0; i < n; i++) v[i] += pair_virial[i];

  if (atom->molecular) {
    if (bondflag)
      for (i = 0; i < n; i++) v[i] += bond_virial[i];
    if (angleflag)
      for (i = 0; i < n; i++) v[i] += angle_virial[i];
    if (dihedralflag)
      for (i = 0; i < n; i++) v[i] += dihedral_virial[i];
    if (improperflag)
      for (i = 0; i < n; i++) v[i] += improper_virial[i];
    if (shakeflag)
      for (i = 0; i < n; i++) v[i] += shake_virial[i];
  }

  if (bodyflag) {
    if (rigidflag) for (i = 0; i < n; i++) v[i] += rigid_virial[i];
    if (poemsflag) for (i = 0; i < n; i++) v[i] += poems_virial[i];
  }

  // sum virial across procs

  MPI_Allreduce(v,virial,n,MPI_DOUBLE,MPI_SUM,world);

  // KSpace virial contribution is already summed across procs

  if (force->kspace)
    for (i = 0; i < n; i++) virial[i] += kspace_virial[i];

  // LJ long-range tail correction

  double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);

  if (force->pair && force->pair->tail_flag)
    for (i = 0; i < n; i++) virial[i] += force->pair->ptail * inv_volume;

  // compute just total average pressure or entire pressure tensor
  // p_total = total pressure = average of trace of tensor
  // for tensor, must first compute 6 components of kinetic energy tensor

  if (tensorflag) {
    temperature->tensor();
    double *ke_tensor = temperature->ke_tensor;
    for (i = 0; i < n; i++)
      p_tensor[i] = (ke_tensor[i] + virial[i]) * inv_volume * nktv2p;
    p_total = (p_tensor[0] + p_tensor[1] + p_tensor[2]) / 3.0;
  } else
    p_total = (temperature->dof * boltz * temperature->t_total + 
	       virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
}
