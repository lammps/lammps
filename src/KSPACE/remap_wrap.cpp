// clang-format off
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

#include "remap_wrap.h"

#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Remap::Remap(LAMMPS *lmp, MPI_Comm comm,
             int in_ilo, int in_ihi, int in_jlo, int in_jhi,
             int in_klo, int in_khi,
             int out_ilo, int out_ihi, int out_jlo, int out_jhi,
             int out_klo, int out_khi,
             int nqty, int permute, int memory,
             int precision, int usecollective) : Pointers(lmp)
{
  plan = remap_3d_create_plan(comm,
                              in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
                              out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
                              nqty,permute,memory,precision,usecollective);
  if (plan == nullptr) error->one(FLERR,"Could not create 3d remap plan");
}

/* ---------------------------------------------------------------------- */

Remap::~Remap()
{
  remap_3d_destroy_plan(plan);
}

/* ---------------------------------------------------------------------- */

void Remap::perform(FFT_SCALAR *in, FFT_SCALAR *out, FFT_SCALAR *buf)
{
  remap_3d(in,out,buf,plan);
}
