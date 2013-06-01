/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Rob Hoy
------------------------------------------------------------------------- */

#include "string.h"
#include "compute_msd_nongauss.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "fix_store.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMSDNonGauss::ComputeMSDNonGauss(LAMMPS *lmp, int narg, char **arg) :
  ComputeMSD(lmp, narg, arg)
{
  size_vector = 3;
}

/* ---------------------------------------------------------------------- */

void ComputeMSDNonGauss::compute_vector()
{
  invoked_vector = update->ntimestep;

  // cm = current center of mass

  double cm[3];
  if (comflag) group->xcm(igroup,masstotal,cm);
  else cm[0] = cm[1] = cm[2] = 0.0;

  // dx,dy,dz = displacement of atom from original position
  // original unwrapped position is stored by fix
  // relative to center of mass if comflag is set
  // for triclinic, need to unwrap current atom coord via h matrix

  double **xoriginal = fix->astore;

  double **x = atom->x;
  int *mask = atom->mask;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double dx,dy,dz;
  int xbox,ybox,zbox;

  double msd[2];
  msd[0] = msd[1] = 0.0;

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + xbox*xprd - cm[0] - xoriginal[i][0];
        dy = x[i][1] + ybox*yprd - cm[1] - xoriginal[i][1];
        dz = x[i][2] + zbox*zprd - cm[2] - xoriginal[i][2];
        msd[0] += dx*dx + dy*dy + dz*dz;
        msd[1] += (dx*dx + dy*dy + dz*dz)*(dx*dx + dy*dy + dz*dz);
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox -
          cm[0] - xoriginal[i][0];
        dy = x[i][1] + h[1]*ybox + h[3]*zbox - cm[1] - xoriginal[i][1];
        dz = x[i][2] + h[2]*zbox - cm[2] - xoriginal[i][2];
        msd[0] += dx*dx + dy*dy + dz*dz;
        msd[1] += (dx*dx + dy*dy + dz*dz)*(dx*dx + dy*dy + dz*dz);
      }
  }

  MPI_Allreduce(msd,vector,2,MPI_DOUBLE,MPI_SUM,world);
  if (nmsd) {
    vector[0] /= nmsd;
    vector[1] /= nmsd;
    vector[2] = (3*vector[1])/(5*vector[0]*vector[0]) - 1;
  }
}
