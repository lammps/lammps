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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_pour.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "fix_gravity.h"
#include "domain.h"
#include "region.h"
#include "region_block.h"
#include "region_cylinder.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

FixPour::FixPour(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix pour command");

  time_depend = 1;

  if (!atom->radius_flag || !atom->rmass_flag)
    error->all(FLERR,"Fix pour requires atom attributes radius, rmass");

  // required args

  ninsert = atoi(arg[3]);
  ntype = atoi(arg[4]);
  seed = atoi(arg[5]);

  if (seed <= 0) error->all(FLERR,"Illegal fix pour command");

  // option defaults

  int iregion = -1;
  radius_lo = radius_hi = 0.5;
  density_lo = density_hi = 1.0;
  volfrac = 0.25;
  maxattempt = 50;
  rate = 0.0;
  vxlo = vxhi = vylo = vyhi = vy = vz = 0.0;

  // optional args

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pour command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->all(FLERR,"Fix pour region ID does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"diam") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix pour command");
      radius_lo = 0.5 * atof(arg[iarg+1]);
      radius_hi = 0.5 * atof(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"dens") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix pour command");
      density_lo = atof(arg[iarg+1]);
      density_hi = atof(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"vol") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix pour command");
      volfrac = atof(arg[iarg+1]);
      maxattempt = atoi(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"rate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pour command");
      rate = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (domain->dimension == 3) {
	if (iarg+6 > narg) error->all(FLERR,"Illegal fix pour command");
	vxlo = atof(arg[iarg+1]);
	vxhi = atof(arg[iarg+2]);
	vylo = atof(arg[iarg+3]);
	vyhi = atof(arg[iarg+4]);
	vz = atof(arg[iarg+5]);
	iarg += 6;
      } else {
	if (iarg+4 > narg) error->all(FLERR,"Illegal fix pour command");
	vxlo = atof(arg[iarg+1]);
	vxhi = atof(arg[iarg+2]);
	vy = atof(arg[iarg+3]);
	vz = 0.0;
	iarg += 4;
      }
    } else error->all(FLERR,"Illegal fix pour command");
  }

  // error checks on region and its extent being inside simulation box

  if (iregion == -1) error->all(FLERR,"Must specify a region in fix pour");
  if (domain->regions[iregion]->bboxflag == 0)
    error->all(FLERR,"Fix pour region does not support a bounding box");
  if (domain->regions[iregion]->dynamic_check())
    error->all(FLERR,"Fix pour region cannot be dynamic");

  if (strcmp(domain->regions[iregion]->style,"block") == 0) {
    region_style = 1;
    xlo = ((RegBlock *) domain->regions[iregion])->xlo;
    xhi = ((RegBlock *) domain->regions[iregion])->xhi;
    ylo = ((RegBlock *) domain->regions[iregion])->ylo;
    yhi = ((RegBlock *) domain->regions[iregion])->yhi;
    zlo = ((RegBlock *) domain->regions[iregion])->zlo;
    zhi = ((RegBlock *) domain->regions[iregion])->zhi;
    if (xlo < domain->boxlo[0] || xhi > domain->boxhi[0] || 
	ylo < domain->boxlo[1] || yhi > domain->boxhi[1] || 
	zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all(FLERR,"Insertion region extends outside simulation box");
  } else if (strcmp(domain->regions[iregion]->style,"cylinder") == 0) {
    region_style = 2;
    char axis = ((RegCylinder *) domain->regions[iregion])->axis;
    xc = ((RegCylinder *) domain->regions[iregion])->c1;
    yc = ((RegCylinder *) domain->regions[iregion])->c2;
    rc = ((RegCylinder *) domain->regions[iregion])->radius;
    zlo = ((RegCylinder *) domain->regions[iregion])->lo;
    zhi = ((RegCylinder *) domain->regions[iregion])->hi;
    if (axis != 'z')
      error->all(FLERR,"Must use a z-axis cylinder with fix pour");
    if (xc-rc < domain->boxlo[0] || xc+rc > domain->boxhi[0] || 
	yc-rc < domain->boxlo[1] || yc+rc > domain->boxhi[1] || 
	zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all(FLERR,"Insertion region extends outside simulation box");
  } else error->all(FLERR,"Must use a block or cylinder region with fix pour");

  if (region_style == 2 && domain->dimension == 2)
    error->all(FLERR,"Must use a block region with fix pour for 2d simulations");

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);

  // allgather arrays

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  recvcounts = new int[nprocs];
  displs = new int[nprocs];

  // grav = gravity in distance/time^2 units
  // assume grav = -magnitude at this point, enforce in init()

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++) {
    if (strcmp(modify->fix[ifix]->style,"gravity") == 0) break;
    if (strcmp(modify->fix[ifix]->style,"gravity/omp") == 0) break;
  }
  if (ifix == modify->nfix) 
    error->all(FLERR,"No fix gravity defined for fix pour");
  grav = - ((FixGravity *) modify->fix[ifix])->magnitude * force->ftm2v;

  // nfreq = timesteps between insertions
  // should be time for a particle to fall from top of insertion region
  //   to bottom, taking into account that the region may be moving
  // set these 2 eqs equal to each other, solve for smallest positive t
  //   x = zhi + vz*t + 1/2 grav t^2
  //   x = zlo + rate*t
  //   gives t = [-(vz-rate) - sqrt((vz-rate)^2 - 2*grav*(zhi-zlo))] / grav
  //   where zhi-zlo > 0, grav < 0, and vz & rate can be either > 0 or < 0

  double v_relative,delta;
  if (domain->dimension == 3) {
    v_relative = vz - rate;
    delta = zhi - zlo;
  } else {
    v_relative = vy - rate;
    delta = yhi - ylo;
  }
  double t = 
    (-v_relative - sqrt(v_relative*v_relative - 2.0*grav*delta)) / grav;
  nfreq = static_cast<int> (t/update->dt + 0.5);

  // 1st insertion on next timestep

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  nfirst = next_reneighbor;
  ninserted = 0;

  // nper = # to insert each time
  // depends on specified volume fraction
  // volume = volume of insertion region
  // volume_one = volume of inserted particle (with max possible radius)
  // in 3d, insure dy >= 1, for quasi-2d simulations

  double volume,volume_one;
  if (domain->dimension == 3) {
    if (region_style == 1) {
      double dy = yhi - ylo;
      if (dy < 1.0) dy = 1.0;
      volume = (xhi-xlo) * dy * (zhi-zlo);
    } else volume = MY_PI*rc*rc * (zhi-zlo);
    volume_one = 4.0/3.0 * MY_PI * radius_hi*radius_hi*radius_hi;
  } else {
    volume = (xhi-xlo) * (yhi-ylo);
    volume_one = MY_PI * radius_hi*radius_hi;
  }

  nper = static_cast<int> (volfrac*volume/volume_one);
  int nfinal = update->ntimestep + 1 + (ninsert-1)/nper * nfreq;

  // print stats

  if (me == 0) {
    if (screen)
      fprintf(screen,
	      "Particle insertion: %d every %d steps, %d by step %d\n",
	      nper,nfreq,ninsert,nfinal);
    if (logfile)
      fprintf(logfile,
	      "Particle insertion: %d every %d steps, %d by step %d\n",
	      nper,nfreq,ninsert,nfinal);
  }
}

/* ---------------------------------------------------------------------- */

FixPour::~FixPour()
{
  delete random;
  delete [] recvcounts;
  delete [] displs;
}

/* ---------------------------------------------------------------------- */

int FixPour::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPour::init()
{
  if (domain->triclinic) error->all(FLERR,"Cannot use fix pour with triclinic box");

  // insure gravity fix exists
  // for 3d must point in -z, for 2d must point in -y
  // else insertion cannot work

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++) {
    if (strcmp(modify->fix[ifix]->style,"gravity") == 0) break;
    if (strcmp(modify->fix[ifix]->style,"gravity/omp") == 0) break;
  }
  if (ifix == modify->nfix) 
    error->all(FLERR,"No fix gravity defined for fix pour");

  double xgrav = ((FixGravity *) modify->fix[ifix])->xgrav;
  double ygrav = ((FixGravity *) modify->fix[ifix])->ygrav;
  double zgrav = ((FixGravity *) modify->fix[ifix])->zgrav;

  if (domain->dimension == 3) {
    if (fabs(xgrav) > EPSILON || fabs(ygrav) > EPSILON ||
	fabs(zgrav+1.0) > EPSILON)
      error->all(FLERR,"Gravity must point in -z to use with fix pour in 3d");
  } else {
    if (fabs(xgrav) > EPSILON || fabs(ygrav+1.0) > EPSILON ||
	fabs(zgrav) > EPSILON)
      error->all(FLERR,"Gravity must point in -y to use with fix pour in 2d");
  }

  double gnew = - ((FixGravity *) modify->fix[ifix])->magnitude * force->ftm2v;
  if (gnew != grav)
    error->all(FLERR,"Gravity changed since fix pour was created");
}

/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixPour::pre_exchange()
{
  int i;

  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  // nnew = # to insert this timestep

  int nnew = nper;
  if (ninserted + nnew > ninsert) nnew = ninsert - ninserted;

  // lo/hi current = z (or y) bounds of insertion region this timestep

  if (domain->dimension == 3) {
    lo_current = zlo + (update->ntimestep - nfirst) * update->dt * rate;
    hi_current = zhi + (update->ntimestep - nfirst) * update->dt * rate;
  } else {
    lo_current = ylo + (update->ntimestep - nfirst) * update->dt * rate;
    hi_current = yhi + (update->ntimestep - nfirst) * update->dt * rate;
  }

  // ncount = # of my atoms that overlap the insertion region
  // nprevious = total of ncount across all procs
  
  int ncount = 0;
  for (i = 0; i < atom->nlocal; i++)
    if (overlap(i)) ncount++;

  int nprevious;
  MPI_Allreduce(&ncount,&nprevious,1,MPI_INT,MPI_SUM,world);

  // xmine is for my atoms
  // xnear is for atoms from all procs + atoms to be inserted

  double **xmine,**xnear;
  memory->create(xmine,ncount,4,"fix_pour:xmine");
  memory->create(xnear,nprevious+nnew,4,"fix_pour:xnear");
  int nnear = nprevious;

  // setup for allgatherv

  int n = 4*ncount;
  MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,world);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++)
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];

  // load up xmine array
  
  double **x = atom->x;
  double *radius = atom->radius;

  ncount = 0;
  for (i = 0; i < atom->nlocal; i++)
    if (overlap(i)) {
      xmine[ncount][0] = x[i][0];
      xmine[ncount][1] = x[i][1];
      xmine[ncount][2] = x[i][2];
      xmine[ncount][3] = radius[i];
      ncount++;
    }

  // perform allgatherv to acquire list of nearby particles on all procs

  double *ptr = NULL;
  if (ncount) ptr = xmine[0];
  MPI_Allgatherv(ptr,4*ncount,MPI_DOUBLE,
		 xnear[0],recvcounts,displs,MPI_DOUBLE,world);

  // insert new atoms into xnear list, one by one
  // check against all nearby atoms and previously inserted ones
  // if there is an overlap then try again at same z (3d) or y (2d) coord
  // else insert by adding to xnear list
  // max = maximum # of insertion attempts for all particles
  // h = height, biased to give uniform distribution in time of insertion

  int success;
  double coord[3],radtmp,delx,dely,delz,rsq,radsum,rn,h;

  int attempt = 0;
  int max = nnew * maxattempt;
  int ntotal = nprevious+nnew;

  while (nnear < ntotal) {
    rn = random->uniform();
    h = hi_current - rn*rn * (hi_current-lo_current);
    radtmp = radius_lo + random->uniform() * (radius_hi-radius_lo);
    success = 0;
    while (attempt < max) {
      attempt++;
      xyz_random(h,coord);
      for (i = 0; i < nnear; i++) {
	delx = coord[0] - xnear[i][0];
	dely = coord[1] - xnear[i][1];
	delz = coord[2] - xnear[i][2];
	rsq = delx*delx + dely*dely + delz*delz;
	radsum = radtmp + xnear[i][3];
	if (rsq <= radsum*radsum) break;
      }
      if (i == nnear) {
	success = 1;
	break;
      }
    }
    if (success) {
      xnear[nnear][0] = coord[0];
      xnear[nnear][1] = coord[1];
      xnear[nnear][2] = coord[2];
      xnear[nnear][3] = radtmp;
      nnear++;
    } else break;
  }

  // warn if not all insertions were performed

  ninserted += nnear-nprevious;
  if (nnear - nprevious < nnew && me == 0)
    error->warning(FLERR,"Less insertions than requested",0);

  // check if new atom is in my sub-box or above it if I'm highest proc
  // if so, add to my list via create_atom()
  // initialize info about the atom
  // type, diameter, density set from fix parameters
  // group mask set to "all" plus fix group
  // z velocity set to what velocity would be if particle
  //   had fallen from top of insertion region
  //   this gives continuous stream of atoms
  //   solution for v from these 2 eqs, after eliminate t:
  //     v = vz + grav*t
  //     coord[2] = hi_current + vz*t + 1/2 grav t^2
  // set npartner for new atom to 0 (assume not touching any others)

  AtomVec *avec = atom->avec;
  int j,m,flag;
  double denstmp,vxtmp,vytmp,vztmp;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  int nfix = modify->nfix;
  Fix **fix = modify->fix;

  for (i = nprevious; i < nnear; i++) {
    coord[0] = xnear[i][0];
    coord[1] = xnear[i][1];
    coord[2] = xnear[i][2];
    radtmp = xnear[i][3];
    denstmp = density_lo + random->uniform() * (density_hi-density_lo);
    if (domain->dimension == 3) {
      vxtmp = vxlo + random->uniform() * (vxhi-vxlo);
      vytmp = vylo + random->uniform() * (vyhi-vylo);
      vztmp = -sqrt(vz*vz + 2.0*grav*(coord[2]-hi_current));
    } else {
      vxtmp = vxlo + random->uniform() * (vxhi-vxlo);
      vytmp = -sqrt(vy*vy + 2.0*grav*(coord[1]-hi_current));
      vztmp = 0.0;
    }

    flag = 0;
    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
	coord[1] >= sublo[1] && coord[1] < subhi[1] &&
	coord[2] >= sublo[2] && coord[2] < subhi[2]) flag = 1;
    else if (domain->dimension == 3 && coord[2] >= domain->boxhi[2] &&
	     comm->myloc[2] == comm->procgrid[2]-1 &&
	     coord[0] >= sublo[0] && coord[0] < subhi[0] &&
	     coord[1] >= sublo[1] && coord[1] < subhi[1]) flag = 1;
    else if (domain->dimension == 2 && coord[1] >= domain->boxhi[1] &&
	     comm->myloc[1] == comm->procgrid[1]-1 &&
	     coord[0] >= sublo[0] && coord[0] < subhi[0]) flag = 1;

    if (flag) {
      avec->create_atom(ntype,coord);
      m = atom->nlocal - 1;
      atom->type[m] = ntype;
      atom->radius[m] = radtmp;
      atom->rmass[m] = 4.0*MY_PI/3.0 * radtmp*radtmp*radtmp * denstmp;
      atom->mask[m] = 1 | groupbit;
      atom->v[m][0] = vxtmp;
      atom->v[m][1] = vytmp;
      atom->v[m][2] = vztmp;
      for (j = 0; j < nfix; j++)
	if (fix[j]->create_attribute) fix[j]->set_arrays(m);
    }
  }

  // reset global natoms
  // set tag # of new particles beyond all previous atoms
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts

  if (nnear - nprevious > 0) {
    atom->natoms += nnear - nprevious;
    if (atom->tag_enable) {
      atom->tag_extend();
      if (atom->map_style) {
	atom->nghost = 0;
	atom->map_init();
	atom->map_set();
      }
    }
  }

  // free local memory

  memory->destroy(xmine);
  memory->destroy(xnear);

  // next timestep to insert

  if (ninserted < ninsert) next_reneighbor += nfreq;
  else next_reneighbor = 0;
}

/* ----------------------------------------------------------------------
   check if particle i could overlap with a particle inserted into region
   return 1 if yes, 0 if no
   use maximum diameter for inserted particle
------------------------------------------------------------------------- */

int FixPour::overlap(int i)
{
  double delta = radius_hi + atom->radius[i];
  double **x = atom->x;

  if (domain->dimension == 3) {
    if (region_style == 1) {
      if (x[i][0] < xlo-delta || x[i][0] > xhi+delta ||
	  x[i][1] < ylo-delta || x[i][1] > yhi+delta ||
	  x[i][2] < lo_current-delta || x[i][2] > hi_current+delta) return 0;
    } else {
      if (x[i][2] < lo_current-delta || x[i][2] > hi_current+delta) return 0;
      double delx = x[i][0] - xc;
      double dely = x[i][1] - yc;
      double rsq = delx*delx + dely*dely;
      double r = rc + delta;
      if (rsq > r*r) return 0;
    }
  } else {
      if (x[i][0] < xlo-delta || x[i][0] > xhi+delta ||
	  x[i][1] < lo_current-delta || x[i][1] > hi_current+delta) return 0;
  }

  return 1;
}

/* ---------------------------------------------------------------------- */

void FixPour::xyz_random(double h, double *coord)
{
  if (domain->dimension == 3) {
    if (region_style == 1) {
      coord[0] = xlo + random->uniform() * (xhi-xlo);
      coord[1] = ylo + random->uniform() * (yhi-ylo);
      coord[2] = h;
    } else {
      double r1,r2;
      while (1) {
	r1 = random->uniform() - 0.5;
	r2 = random->uniform() - 0.5;
	if (r1*r1 + r2*r2 < 0.25) break;
      }
      coord[0] = xc + 2.0*r1*rc;
      coord[1] = yc + 2.0*r2*rc;
      coord[2] = h;
    }
  } else {
    coord[0] = xlo + random->uniform() * (xhi-xlo);
    coord[1] = h;
    coord[2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixPour::reset_dt()
{
  error->all(FLERR,"Cannot change timestep with fix pour");
}
