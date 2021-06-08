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

/*------------------------------------------------------------------------
  Contributing Authors : Romain Vermorel (LFCR), Laurent Joly (ULyon)
  --------------------------------------------------------------------------*/

#include "compute_stress_mop_profile.h"

#include <cmath>
#include <cstring>

#include "atom.h"
#include "update.h"
#include "domain.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

enum{X,Y,Z};
enum{LOWER,CENTER,UPPER,COORD};
enum{TOTAL,CONF,KIN};

#define BIG 1000000000

/* ---------------------------------------------------------------------- */

ComputeStressMopProfile::ComputeStressMopProfile(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal compute stress/mop/profile command");

  MPI_Comm_rank(world,&me);

  // set compute mode and direction of plane(s) for pressure calculation

  if (strcmp(arg[3],"x")==0) {
    dir = X;
  } else if (strcmp(arg[3],"y")==0) {
    dir = Y;
  } else if (strcmp(arg[3],"z")==0) {
    dir = Z;
  } else error->all(FLERR,"Illegal compute stress/mop/profile command");

  // bin parameters

  if (strcmp(arg[4],"lower") == 0) originflag = LOWER;
  else if (strcmp(arg[4],"center") == 0) originflag = CENTER;
  else if (strcmp(arg[4],"upper") == 0) originflag = UPPER;
  else originflag = COORD;
  if (originflag == COORD)
    origin = utils::numeric(FLERR,arg[4],false,lmp);
  delta = utils::numeric(FLERR,arg[5],false,lmp);
  invdelta = 1.0/delta;

  // parse values until one isn't recognized

  which = new int[3*(narg-6)];
  nvalues = 0;
  int i;

  int iarg=6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"conf") == 0) {
      for (i=0; i<3; i++) {
        which[nvalues] = CONF;
        nvalues++;
      }
    } else if (strcmp(arg[iarg],"kin") == 0) {
      for (i=0; i<3; i++) {
        which[nvalues] = KIN;
        nvalues++;
      }
    } else if (strcmp(arg[iarg],"total") == 0) {
      for (i=0; i<3; i++) {
        which[nvalues] = TOTAL;
        nvalues++;
      }
    } else error->all(FLERR, "Illegal compute stress/mop/profile command"); //break;

    iarg++;
  }

  // check domain related errors

  // 3D only

  if (domain->dimension < 3)
    error->all(FLERR, "Compute stress/mop/profile incompatible with simulation dimension");

  // orthogonal simulation box

  if (domain->triclinic != 0)
    error->all(FLERR, "Compute stress/mop/profile incompatible with triclinic simulation box");

  // initialize some variables

  nbins = 0;
  coord = coordp = nullptr;
  values_local = values_global = array = nullptr;

  // bin setup

  setup_bins();

  // this fix produces a global array

  memory->create(array,nbins,1+nvalues,"stress/mop/profile:array");
  size_array_rows = nbins;
  size_array_cols = 1 + nvalues;

  array_flag = 1;
  extarray = 0;
}

/* ---------------------------------------------------------------------- */

ComputeStressMopProfile::~ComputeStressMopProfile()
{

  delete [] which;

  memory->destroy(coord);
  memory->destroy(coordp);
  memory->destroy(values_local);
  memory->destroy(values_global);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeStressMopProfile::init()
{

  // conversion constants

  nktv2p = force->nktv2p;
  ftm2v = force->ftm2v;

  // plane area

  area = 1;
  int i;
  for (i=0; i<3; i++) {
    if (i!=dir) area = area*domain->prd[i];
  }

  // timestep Value

  dt = update->dt;

  // Error check
  // Compute stress/mop/profile requires fixed simulation box

  if (domain->box_change_size || domain->box_change_shape || domain->deform_flag)
    error->all(FLERR, "Compute stress/mop/profile requires a fixed simulation box");

  //This compute requires a pair style with pair_single method implemented

  if (force->pair == nullptr)
    error->all(FLERR,"No pair style is defined for compute stress/mop/profile");
  if (force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support compute stress/mop/profile");

  // Warnings

  if (me==0) {

    //Compute stress/mop/profile only accounts for pair interactions.
    // issue a warning if any intramolecular potential or Kspace is defined.

    if (force->bond!=nullptr)
      error->warning(FLERR,"compute stress/mop/profile does not account for bond potentials");
    if (force->angle!=nullptr)
      error->warning(FLERR,"compute stress/mop/profile does not account for angle potentials");
    if (force->dihedral!=nullptr)
      error->warning(FLERR,"compute stress/mop/profile does not account for dihedral potentials");
    if (force->improper!=nullptr)
      error->warning(FLERR,"compute stress/mop/profile does not account for improper potentials");
    if (force->kspace!=nullptr)
      error->warning(FLERR,"compute stress/mop/profile does not account for kspace contributions");
  }

  // need an occasional half neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeStressMopProfile::init_list(int /* id */, NeighList *ptr)
{
  list = ptr;
}


/* ----------------------------------------------------------------------
   compute output array
   ------------------------------------------------------------------------- */

void ComputeStressMopProfile::compute_array()
{
  invoked_array = update->ntimestep;

  //Compute pressures on separate procs
  compute_pairs();

  // sum pressure contributions over all procs
  MPI_Allreduce(&values_local[0][0],&values_global[0][0],nbins*nvalues,
                MPI_DOUBLE,MPI_SUM,world);

  int ibin,m,mo;
  for (ibin=0; ibin<nbins; ibin++) {
    array[ibin][0] = coord[ibin][0];
    mo=1;

    m = 0;
    while (m<nvalues) {
      array[ibin][m+mo] = values_global[ibin][m];
      m++;
    }
  }
}


/*------------------------------------------------------------------------
  compute pressure contribution of local proc
  -------------------------------------------------------------------------*/

void ComputeStressMopProfile::compute_pairs()

{
  int i,j,m,ii,jj,inum,jnum,itype,jtype,ibin;
  double delx,dely,delz;
  double rsq,fpair,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double pos,pos1;

  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;


  // zero out arrays for one sample
  for (m = 0; m < nbins; m++) {
    for (i = 0; i < nvalues; i++) values_local[m][i] = 0.0;
  }

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  // parse values

  double xi[3];
  double vi[3];
  double fi[3];
  double xj[3];

  m = 0;
  while (m<nvalues) {
    if (which[m] == CONF || which[m] == TOTAL) {

      // Compute configurational contribution to pressure

      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];

        xi[0] = atom->x[i][0];
        xi[1] = atom->x[i][1];
        xi[2] = atom->x[i][2];
        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          factor_lj = special_lj[sbmask(j)];
          factor_coul = special_coul[sbmask(j)];
          j &= NEIGHMASK;

          // skip if neither I nor J are in group

          if (!(mask[i] & groupbit || mask[j] & groupbit)) continue;

          xj[0] = atom->x[j][0];
          xj[1] = atom->x[j][1];
          xj[2] = atom->x[j][2];
          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          delz = xi[2] - xj[2];
          rsq = delx*delx + dely*dely + delz*delz;
          jtype = type[j];
          if (rsq >= cutsq[itype][jtype]) continue;

          if (newton_pair || j < nlocal) {

            for (ibin=0;ibin<nbins;ibin++) {
              pos = coord[ibin][0];
              pos1 = coordp[ibin][0];

              //check if ij pair is across plane, add contribution to pressure

              if ( ((xi[dir]>pos) && (xj[dir]<pos))
                   || ((xi[dir]>pos1) && (xj[dir]<pos1))) {

                pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

                values_local[ibin][m] += fpair*(xi[0]-xj[0])/area*nktv2p;
                values_local[ibin][m+1] += fpair*(xi[1]-xj[1])/area*nktv2p;
                values_local[ibin][m+2] += fpair*(xi[2]-xj[2])/area*nktv2p;

              } else if ( ((xi[dir]<pos) && (xj[dir]>pos))
                          || ((xi[dir]<pos1) && (xj[dir]>pos1))) {

                pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

                values_local[ibin][m] -= fpair*(xi[0]-xj[0])/area*nktv2p;
                values_local[ibin][m+1] -= fpair*(xi[1]-xj[1])/area*nktv2p;
                values_local[ibin][m+2] -= fpair*(xi[2]-xj[2])/area*nktv2p;
              }
            }
          } else {

            for (ibin=0;ibin<nbins;ibin++) {
              pos = coord[ibin][0];
              pos1 = coordp[ibin][0];

              //check if ij pair is across plane, add contribution to pressure

              if ( ((xi[dir]>pos) && (xj[dir]<pos))
                   || ((xi[dir]>pos1) && (xj[dir]<pos1))) {

                pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

                values_local[ibin][m] += fpair*(xi[0]-xj[0])/area*nktv2p;
                values_local[ibin][m+1] += fpair*(xi[1]-xj[1])/area*nktv2p;
                values_local[ibin][m+2] += fpair*(xi[2]-xj[2])/area*nktv2p;
              }
            }
          }
        }
      }
    }

    // compute kinetic contribution to pressure
    // counts local particles transfers across the plane

    if (which[m] == KIN || which[m] == TOTAL) {

      double sgn;

      for (int i = 0; i < nlocal; i++) {

        // skip if I is not in group

        if (mask[i] & groupbit) {

          itype = type[i];

          //coordinates at t
          xi[0] = atom->x[i][0];
          xi[1] = atom->x[i][1];
          xi[2] = atom->x[i][2];

          //velocities at t
          vi[0] = atom->v[i][0];
          vi[1] = atom->v[i][1];
          vi[2] = atom->v[i][2];

          //forces at t
          fi[0] = atom->f[i][0];
          fi[1] = atom->f[i][1];
          fi[2] = atom->f[i][2];

          //coordinates at t-dt (based on Velocity-Verlet alg.)
          if (rmass) {
            xj[0] = xi[0]-vi[0]*dt+fi[0]/2/rmass[i]*dt*dt*ftm2v;
            xj[1] = xi[1]-vi[1]*dt+fi[1]/2/rmass[i]*dt*dt*ftm2v;
            xj[2] = xi[2]-vi[2]*dt+fi[2]/2/rmass[i]*dt*dt*ftm2v;
          } else {
            xj[0] = xi[0]-vi[0]*dt+fi[0]/2/mass[itype]*dt*dt*ftm2v;
            xj[1] = xi[1]-vi[1]*dt+fi[1]/2/mass[itype]*dt*dt*ftm2v;
            xj[2] = xi[2]-vi[2]*dt+fi[2]/2/mass[itype]*dt*dt*ftm2v;
          }

          for (ibin=0;ibin<nbins;ibin++) {
            pos = coord[ibin][0];
            pos1 = coordp[ibin][0];

            if (((xi[dir]-pos)*(xj[dir]-pos)*(xi[dir]-pos1)*(xj[dir]-pos1)<0)) {

              sgn = copysign(1.0,vi[dir]);

              //approximate crossing velocity by v(t-dt/2) (based on Velocity-Verlet alg.)
              if (rmass) {
                double vcross[3];
                vcross[0] = vi[0]-fi[0]/rmass[i]/2*ftm2v*dt;
                vcross[1] = vi[1]-fi[1]/rmass[i]/2*ftm2v*dt;
                vcross[2] = vi[2]-fi[2]/rmass[i]/2*ftm2v*dt;

                values_local[ibin][m] += rmass[i]*vcross[0]*sgn/dt/area*nktv2p/ftm2v;
                values_local[ibin][m+1] += rmass[i]*vcross[1]*sgn/dt/area*nktv2p/ftm2v;
                values_local[ibin][m+2] += rmass[i]*vcross[2]*sgn/dt/area*nktv2p/ftm2v;
              } else {
                double vcross[3];
                vcross[0] = vi[0]-fi[0]/mass[itype]/2*ftm2v*dt;
                vcross[1] = vi[1]-fi[1]/mass[itype]/2*ftm2v*dt;
                vcross[2] = vi[2]-fi[2]/mass[itype]/2*ftm2v*dt;

                values_local[ibin][m] += mass[itype]*vcross[0]*sgn/dt/area*nktv2p/ftm2v;
                values_local[ibin][m+1] += mass[itype]*vcross[1]*sgn/dt/area*nktv2p/ftm2v;
                values_local[ibin][m+2] += mass[itype]*vcross[2]*sgn/dt/area*nktv2p/ftm2v;
              }
            }
          }
        }
      }
    }
    m+=3;
  }
}

/* ----------------------------------------------------------------------
   setup 1d bins and their extent and coordinates
   called at init()
   ------------------------------------------------------------------------- */

void ComputeStressMopProfile::setup_bins()
{
  int i,n;
  double lo = 0.0, hi = 0.0;

  double *boxlo,*boxhi;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;

  if (originflag == LOWER) origin = boxlo[dir];
  else if (originflag == UPPER) origin = boxhi[dir];
  else if (originflag == CENTER)
    origin = 0.5 * (boxlo[dir] + boxhi[dir]);

  if (origin < boxlo[dir]) {
    error->all(FLERR,"Origin of bins for compute stress/mop/profile is out of bounds" );
  } else {
    n = static_cast<int> ((origin - boxlo[dir]) * invdelta);
    lo = origin - n*delta;
  }
  if (origin < boxhi[dir]) {
    n = static_cast<int> ((boxhi[dir] - origin) * invdelta);
    hi = origin + n*delta;
  } else {
    error->all(FLERR,"Origin of bins for compute stress/mop/profile is out of bounds" );
  }

  offset = lo;
  nbins = static_cast<int> ((hi-lo) * invdelta + 1.5);

  //allocate bin arrays
  memory->create(coord,nbins,1,"stress/mop/profile:coord");
  memory->create(coordp,nbins,1,"stress/mop/profile:coordp");
  memory->create(values_local,nbins,nvalues,"stress/mop/profile:values_local");
  memory->create(values_global,nbins,nvalues,"stress/mop/profile:values_global");

  // set bin coordinates
  for (i = 0; i < nbins; i++) {
    coord[i][0] = offset + i*delta;
    if (coord[i][0] < (domain->boxlo[dir]+domain->prd_half[dir])) {
      coordp[i][0] = coord[i][0] + domain->prd[dir];
    } else {
      coordp[i][0] = coord[i][0] - domain->prd[dir];
    }
  }
}
