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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "compute_event_displace.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "fix_event.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;

#define INVOKED_SCALAR 1

/* ---------------------------------------------------------------------- */

ComputeEventDisplace::ComputeEventDisplace(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute event/displace command");

  scalar_flag = 1;
  extscalar = 0;

  double displace_dist = force->numeric(FLERR,arg[3]);
  if (displace_dist <= 0.0)
    error->all(FLERR,"Distance must be > 0 for compute event/displace");
  displace_distsq = displace_dist * displace_dist;

  // fix event ID will be set later by PRD

  id_event = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEventDisplace::~ComputeEventDisplace()
{
  delete [] id_event;
}

/* ---------------------------------------------------------------------- */

void ComputeEventDisplace::init()
{
  // if id_event is not set, this compute is not active
  // if set by PRD, then find fix which stores original atom coords
  // check if it is correct style

  if (id_event != NULL) {
    int ifix = modify->find_fix(id_event);
    if (ifix < 0) error->all(FLERR,
                             "Could not find compute event/displace fix ID");
    fix_event = (FixEvent*) modify->fix[ifix];

    if (strcmp(fix_event->style,"EVENT/PRD") != 0 &&
        strcmp(fix_event->style,"EVENT/TAD") != 0)
      error->all(FLERR,"Compute event/displace has invalid fix event assigned");
  }

  triclinic = domain->triclinic;
}

/* ----------------------------------------------------------------------
   return non-zero if an atom has moved > displace_dist since last event
------------------------------------------------------------------------- */

double ComputeEventDisplace::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (id_event == NULL) return 0.0;

  double event = 0.0;
  double **xevent = fix_event->array_atom;

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  int xbox,ybox,zbox;
  double dx,dy,dz,rsq;

  if (triclinic == 0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + xbox*xprd - xevent[i][0];
        dy = x[i][1] + ybox*yprd - xevent[i][1];
        dz = x[i][2] + zbox*zprd - xevent[i][2];
        rsq = dx*dx + dy*dy + dz*dz;
        if (rsq >= displace_distsq) {
          event = 1.0;
          break;
        }
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox - xevent[i][0];
        dy = x[i][1] + h[1]*ybox + h[3]*zbox - xevent[i][1];
        dz = x[i][2] + h[2]*zbox - xevent[i][2];
        rsq = dx*dx + dy*dy + dz*dz;
        if (rsq >= displace_distsq) {
          event = 1.0;
          break;
        }
      }
  }

  MPI_Allreduce(&event,&scalar,1,MPI_DOUBLE,MPI_SUM,world);

  return scalar;
}


/* ---------------------------------------------------------------------- */

void ComputeEventDisplace::reset_extra_compute_fix(const char *id_new)
{
  delete [] id_event;
  id_event = NULL;
  if (id_new == NULL) return;

  int n = strlen(id_new) + 1;
  id_event = new char[n];
  strcpy(id_event,id_new);
}
