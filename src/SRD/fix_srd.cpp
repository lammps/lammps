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
   Contributing authors: Jeremy Lechman (SNL), Pieter in 't Veld (BASF)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_srd.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "group.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "comm.h"
#include "modify.h"
#include "fix_deform.h"
#include "fix_wall_srd.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{SLIP,NOSLIP};
enum{SPHERE,ELLIPSOID,LINE,TRIANGLE,WALL};
enum{INSIDE_ERROR,INSIDE_WARN,INSIDE_IGNORE};
enum{BIG_MOVE,SRD_MOVE,SRD_ROTATE};
enum{CUBIC_ERROR,CUBIC_WARN};
enum{SHIFT_NO,SHIFT_YES,SHIFT_POSSIBLE};

enum{NO_REMAP,X_REMAP,V_REMAP};                   // same as fix_deform.cpp

#define EINERTIA 0.2          // moment of inertia prefactor for ellipsoid

#define ATOMPERBIN 30
#define BIG 1.0e20
#define VBINSIZE 5
#define TOLERANCE 0.00001
#define MAXITER 20

//#define SRD_DEBUG 1
//#define SRD_DEBUG_ATOMID 58
//#define SRD_DEBUG_TIMESTEP 449

/* ---------------------------------------------------------------------- */

FixSRD::FixSRD(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 8) error->all(FLERR,"Illegal fix srd command");

  restart_pbc = 1;
  vector_flag = 1;
  size_vector = 12;
  global_freq = 1;
  extvector = 0;

  nevery = atoi(arg[3]);

  bigexist = 1;
  if (strcmp(arg[4],"NULL") == 0) bigexist = 0;
  else biggroup = group->find(arg[4]);

  temperature_srd = atof(arg[5]);
  gridsrd = atof(arg[6]);
  int seed = atoi(arg[7]);

  // parse options

  lamdaflag = 0;
  collidestyle = NOSLIP;
  overlap = 0;
  insideflag = INSIDE_ERROR;
  exactflag = 1;
  radfactor = 1.0;
  maxbounceallow = 0;
  gridsearch = gridsrd;
  cubicflag = CUBIC_ERROR;
  cubictol = 0.01;
  shiftuser = SHIFT_NO;
  shiftseed = 0;
  tstat = 0;

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"lamda") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix srd command");
      lamda = atof(arg[iarg+1]);
      lamdaflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"collision") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix srd command");
      if (strcmp(arg[iarg+1],"slip") == 0) collidestyle = SLIP;
      else if (strcmp(arg[iarg+1],"noslip") == 0) collidestyle = NOSLIP;
      else error->all(FLERR,"Illegal fix srd command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"overlap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix srd command");
      if (strcmp(arg[iarg+1],"yes") == 0) overlap = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) overlap = 0;
      else error->all(FLERR,"Illegal fix srd command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"inside") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix srd command");
      if (strcmp(arg[iarg+1],"error") == 0) insideflag = INSIDE_ERROR;
      else if (strcmp(arg[iarg+1],"warn") == 0) insideflag = INSIDE_WARN;
      else if (strcmp(arg[iarg+1],"ignore") == 0) insideflag = INSIDE_IGNORE;
      else error->all(FLERR,"Illegal fix srd command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"exact") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix srd command");
      if (strcmp(arg[iarg+1],"yes") == 0) exactflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) exactflag = 0;
      else error->all(FLERR,"Illegal fix srd command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"radius") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix srd command");
      radfactor = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"bounce") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix srd command");
      maxbounceallow = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"search") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix srd command");
      gridsearch = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"cubic") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix srd command");
      if (strcmp(arg[iarg+1],"error") == 0) cubicflag = CUBIC_ERROR;
      else if (strcmp(arg[iarg+1],"warn") == 0) cubicflag = CUBIC_WARN;
      else error->all(FLERR,"Illegal fix srd command");
      cubictol = atof(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"shift") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix srd command");
      else if (strcmp(arg[iarg+1],"no") == 0) shiftuser = SHIFT_NO;
      else if (strcmp(arg[iarg+1],"yes") == 0) shiftuser = SHIFT_YES;
      else if (strcmp(arg[iarg+1],"possible") == 0) shiftuser = SHIFT_POSSIBLE;
      else error->all(FLERR,"Illegal fix srd command");
      shiftseed = atoi(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"tstat") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix srd command");
      if (strcmp(arg[iarg+1],"no") == 0) tstat = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) tstat = 1;
      else error->all(FLERR,"Illegal fix srd command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix srd command");
  }

  // error check

  if (nevery <= 0) error->all(FLERR,"Illegal fix srd command");
  if (bigexist && biggroup < 0) 
    error->all(FLERR,"Could not find fix srd group ID");
  if (gridsrd <= 0.0) error->all(FLERR,"Illegal fix srd command");
  if (temperature_srd <= 0.0) error->all(FLERR,"Illegal fix srd command");
  if (seed <= 0) error->all(FLERR,"Illegal fix srd command");
  if (radfactor <= 0.0) error->all(FLERR,"Illegal fix srd command");
  if (maxbounceallow < 0) error->all(FLERR,"Illegal fix srd command");
  if (lamdaflag && lamda <= 0.0) error->all(FLERR,"Illegal fix srd command");
  if (gridsearch <= 0.0) error->all(FLERR,"Illegal fix srd command");
  if (cubictol < 0.0 || cubictol > 1.0) 
    error->all(FLERR,"Illegal fix srd command");
  if ((shiftuser == SHIFT_YES || shiftuser == SHIFT_POSSIBLE) && 
      shiftseed <= 0) error->all(FLERR,"Illegal fix srd command");

  // initialize Marsaglia RNG with processor-unique seed

  me = comm->me;
  nprocs = comm->nprocs;

  random = new RanMars(lmp,seed + me);

  // if requested, initialize shift RNG, same on every proc

  if (shiftuser == SHIFT_YES || shiftuser == SHIFT_POSSIBLE)
    randomshift = new RanPark(lmp,shiftseed);
  else randomshift = NULL;

  // initialize data structs and flags

  if (bigexist) biggroupbit = group->bitmask[biggroup];
  else biggroupbit = 0;

  nmax = 0;
  binhead = NULL;
  maxbin1 = 0;
  binnext = NULL;
  maxbuf = 0;
  sbuf1 = sbuf2 = rbuf1 = rbuf2 = NULL;

  shifts[0].maxvbin = shifts[1].maxvbin = 0;
  shifts[0].vbin = shifts[1].vbin = NULL;

  shifts[0].maxbinsq = shifts[1].maxbinsq = 0;
  for (int ishift = 0; ishift < 2; ishift++)
    for (int iswap = 0; iswap < 6; iswap++)
      shifts[ishift].bcomm[iswap].sendlist = 
	shifts[ishift].bcomm[iswap].recvlist = NULL;

  maxbin2 = 0;
  nbinbig = NULL;
  binbig = NULL;
  binsrd = NULL;

  nstencil = maxstencil = 0;
  stencil = NULL;

  maxbig = 0;
  biglist = NULL;

  stats_flag = 1;
  for (int i = 0; i < size_vector; i++) stats_all[i] = 0.0;

  initflag = 0;

  srd_bin_temp = 0.0;
  srd_bin_count = 0;

  // atom style pointers to particles that store bonus info

  avec_ellipsoid = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  avec_line = (AtomVecLine *) atom->style_match("line");
  avec_tri = (AtomVecTri *) atom->style_match("tri");

  // fix parameters

  if (collidestyle == SLIP) comm_reverse = 3;
  else comm_reverse = 6;
  force_reneighbor = 1;
}

/* ---------------------------------------------------------------------- */

FixSRD::~FixSRD()
{
  delete random;
  delete randomshift;

  memory->destroy(binhead);
  memory->destroy(binnext);
  memory->destroy(sbuf1);
  memory->destroy(sbuf2);
  memory->destroy(rbuf1);
  memory->destroy(rbuf2);

  memory->sfree(shifts[0].vbin);
  memory->sfree(shifts[1].vbin);
  for (int ishift = 0; ishift < 2; ishift++)
    for (int iswap = 0; iswap < 6; iswap++) {
      memory->destroy(shifts[ishift].bcomm[iswap].sendlist);
      memory->destroy(shifts[ishift].bcomm[iswap].recvlist);
    }

  memory->destroy(nbinbig);
  memory->destroy(binbig);
  memory->destroy(binsrd);
  memory->destroy(stencil);
  memory->sfree(biglist);
}

/* ---------------------------------------------------------------------- */

int FixSRD::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSRD::init()
{
  // error checks

  if (force->newton_pair == 0) 
    error->all(FLERR,"Fix srd requires newton pair on");
  if (bigexist && comm->ghost_velocity == 0)
    error->all(FLERR,"Fix srd requires ghost atoms store velocity");
  if (bigexist && collidestyle == NOSLIP && !atom->torque_flag)
    error->all(FLERR,"Fix SRD no-slip requires atom attribute torque");
  if (initflag && update->dt != dt_big)
    error->all(FLERR,"Cannot change timestep once fix srd is setup");

  // orthogonal vs triclinic simulation box
  // could be static or shearing box

  triclinic = domain->triclinic;

  // wallexist = 1 if SRD wall(s) are defined

  wallexist = 0;
  for (int m = 0; m < modify->nfix; m++) {
    if (strcmp(modify->fix[m]->style,"wall/srd") == 0) {
      if (wallexist) error->all(FLERR,"Cannot use fix wall/srd more than once");
      wallexist = 1;
      wallfix = (FixWallSRD *) modify->fix[m];
      nwall = wallfix->nwall;
      wallvarflag = wallfix->varflag;
      wallwhich = wallfix->wallwhich;
      xwall = wallfix->xwall;
      xwallhold = wallfix->xwallhold;
      vwall = wallfix->vwall;
      fwall = wallfix->fwall;
      walltrigger = 0.5 * neighbor->skin;
      if (wallfix->overlap && overlap == 0 && me == 0)
	error->warning(FLERR,
		       "Fix SRD walls overlap but fix srd overlap not set");
    }
  }

  // set change_flags if box size or shape changes

  change_size = change_shape = deformflag = 0;
  if (domain->nonperiodic == 2) change_size = 1;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->box_change) {
      if (modify->fix[i]->box_change_size) change_size = 1;
      if (modify->fix[i]->box_change_shape) change_shape = 1;
      if (strcmp(modify->fix[i]->style,"deform") == 0) {
	deformflag = 1;
	FixDeform *deform = (FixDeform *) modify->fix[i];
	if (deform->box_change_shape && deform->remapflag != V_REMAP)
	  error->all(FLERR,"Using fix srd with inconsistent "
		     "fix deform remap option");
      }
    }

  if (deformflag && tstat == 0 && me == 0)
    error->warning(FLERR,
		   "Using fix srd with box deformation but no SRD thermostat");

  // parameterize based on current box volume

  dimension = domain->dimension;
  parameterize();

  // limit initial SRD velocities if necessary

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double vsq;
  nrescale = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > vmaxsq) {
	nrescale++;
	MathExtra::scale3(vmax/sqrt(vsq),v[i]);
      }
    }

  int all;
  MPI_Allreduce(&nrescale,&all,1,MPI_INT,MPI_SUM,world);
  if (me == 0) {
    if (screen)
      fprintf(screen,"  # of rescaled SRD velocities = %d\n",all);
    if (logfile)
      fprintf(logfile,"  # of rescaled SRD velocities = %d\n",all);
  }

  velocity_stats(igroup);
  if (bigexist) velocity_stats(biggroup);

  // zero per-run stats

  bouncemaxnum = 0;
  bouncemax = 0;
  reneighcount = 0;
  initflag = 1;

  next_reneighbor = -1;
}

/* ---------------------------------------------------------------------- */

void FixSRD::setup(int vflag)
{
  setup_bounds();

  if (dist_srd_reneigh < nevery*dt_big*vmax && me == 0)
    error->warning(FLERR,
		   "Fix srd SRD moves may trigger frequent reneighboring");

  // setup search bins and search stencil based on these distances

  if (bigexist || wallexist) {
    setup_search_bins();
    setup_search_stencil();
  } else nbins2 = 0;

  // perform first bining of SRD and big particles and walls
  // set reneighflag to turn off SRD rotation
  // don't do SRD rotation in setup, only during timestepping

  reneighflag = BIG_MOVE;
  pre_neighbor();
}

/* ----------------------------------------------------------------------
   assign SRD particles to bins
   assign big particles to all bins they overlap
------------------------------------------------------------------------- */

void FixSRD::pre_neighbor()
{
  int i,j,m,ix,iy,iz,jx,jy,jz,ibin,jbin,lo,hi;
  double rsq,cutbinsq;
  double xlamda[3];

  // grow SRD per-atom bin arrays if necessary

  if (atom->nlocal > nmax) {
    nmax = atom->nmax;
    memory->destroy(binsrd);
    memory->destroy(binnext);
    memory->create(binsrd,nmax,"fix/srd:binsrd");
    memory->create(binnext,nmax,"fix/srd:binnext");
  }

  // setup and grow BIG info list if necessary
  // set index ptrs to BIG particles and to WALLS
  // big_static() adds static properties to info list

  if (bigexist || wallexist) {
    if (bigexist) {
      if (biggroup == atom->firstgroup) nbig = atom->nfirst + atom->nghost;
      else {
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	nbig = atom->nghost;
	for (i = 0; i < nlocal; i++)
	  if (mask[i] & biggroupbit) nbig++;
      }
    } else nbig = 0;

    int ninfo = nbig;
    if (wallexist) ninfo += nwall;

    if (ninfo > maxbig) {
      maxbig = ninfo;
      memory->destroy(biglist);
      biglist = (Big *) memory->smalloc(maxbig*sizeof(Big),"fix/srd:biglist");
    }

    if (bigexist) {
      int *mask = atom->mask;
      int nlocal = atom->nlocal;
      if (biggroup == atom->firstgroup) nlocal = atom->nfirst;
      nbig = 0;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & biggroupbit) biglist[nbig++].index = i;
      int nall = atom->nlocal + atom->nghost;
      for (i = atom->nlocal; i < nall; i++)
	if (mask[i] & biggroupbit) biglist[nbig++].index = i;
      big_static();
    }

    if (wallexist) {
      for (m = 0; m < nwall; m++) {
	biglist[nbig+m].index = m;
	biglist[nbig+m].type = WALL;
      }
      wallfix->wall_params(1);
    }
  }

  // if simulation box size changes, reset velocity bins
  // if big particles exist, reset search bins if box size or shape changes,
  //   b/c finite-size particles will overlap different bins as the box tilts

  if (change_size) setup_bounds();
  if (change_size) setup_velocity_bins();
  if ((change_size || change_shape) && (bigexist || wallexist)) {
    setup_search_bins();
    setup_search_stencil();
  }

  // map each owned & ghost big particle to search bins it overlaps
  // zero out bin counters for big particles
  // if firstgroup is defined, only loop over first and ghost particles
  // for each big particle: loop over stencil to find overlap bins

  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int nfirst = nlocal;
  if (bigexist && biggroup == atom->firstgroup) nfirst = atom->nfirst;

  if (bigexist || wallexist)
    for (i = 0; i < nbins2; i++)
      nbinbig[i] = 0;

  if (bigexist) {
    i = nbig = 0;
    while (i < nall) {
      if (mask[i] & biggroupbit) {
	ix = static_cast<int> ((x[i][0]-xblo2)*bininv2x);
	iy = static_cast<int> ((x[i][1]-yblo2)*bininv2y);
	iz = static_cast<int> ((x[i][2]-zblo2)*bininv2z);
	ibin = iz*nbin2y*nbin2x + iy*nbin2x + ix;

	if (ix < 0 || ix >= nbin2x || iy < 0 || iy >= nbin2y || 
	    iz < 0 || iz >= nbin2z)
	  error->one(FLERR,"Fix SRD: bad search bin assignment");
      
	cutbinsq = biglist[nbig].cutbinsq;
	for (j = 0; j < nstencil; j++) {
	  jx = ix + stencil[j][0];
	  jy = iy + stencil[j][1];
	  jz = iz + stencil[j][2];
	  
	  if (jx < 0 || jx >= nbin2x || jy < 0 || jy >= nbin2y || 
	      jz < 0 || jz >= nbin2z) {
	    printf("Big particle %d %d %g %g %g\n",
		   atom->tag[i],i,x[i][0],x[i][1],x[i][2]);
	    printf("Bin indices: %d %d %d, %d %d %d, %d %d %d\n",
		   ix,iy,iz,jx,jy,jz,nbin2x,nbin2y,nbin2z);
	    error->one(FLERR,"Fix SRD: bad stencil bin for big particle");
	  }
	  rsq = point_bin_distance(x[i],jx,jy,jz);
	  if (rsq < cutbinsq) {
	    jbin = ibin + stencil[j][3];
	    if (nbinbig[jbin] == ATOMPERBIN)
	      error->one(FLERR,"Fix SRD: too many big particles in bin");
	    binbig[jbin][nbinbig[jbin]++] = nbig;
	  }
	}
	nbig++;
      }

      i++;
      if (i == nfirst) i = nlocal;
    }
  }

  // map each wall to search bins it covers, up to non-periodic boundary
  // if wall moves, add walltrigger to its position
  // this insures it is added to all search bins it may move into
  // may not overlap any of my search bins

  if (wallexist) {
    double delta = 0.0;
    if (wallvarflag) delta = walltrigger;

    for (m = 0; m < nwall; m++) {
      int dim = wallwhich[m] / 2;
      int side = wallwhich[m] % 2;

      if (dim == 0) {
	if (side == 0) {
	  hi = static_cast<int> ((xwall[m]+delta-xblo2)*bininv2x);
	  if (hi < 0) continue;
	  if (hi >= nbin2x) error->all(FLERR,
				       "Fix SRD: bad search bin assignment");
	  lo = 0;
	} else {
	  lo = static_cast<int> ((xwall[m]-delta-xblo2)*bininv2x);
	  if (lo >= nbin2x) continue;
	  if (lo < 0) error->all(FLERR,"Fix SRD: bad search bin assignment");
	  hi = nbin2x-1;
	}

	for (ix = lo; ix <= hi; ix++)
	  for (iy = 0; iy < nbin2y; iy++)
	    for (iz = 0; iz < nbin2z; iz++) {
	      ibin = iz*nbin2y*nbin2x + iy*nbin2x + ix;
	      if (nbinbig[ibin] == ATOMPERBIN)
		error->all(FLERR,"Fix SRD: too many walls in bin");
	      binbig[ibin][nbinbig[ibin]++] = nbig+m;
	    }

      } else if (dim == 1) {
	if (side == 0) {
	  hi = static_cast<int> ((xwall[m]+delta-yblo2)*bininv2y);
	  if (hi < 0) continue;
	  if (hi >= nbin2y) error->all(FLERR,
				       "Fix SRD: bad search bin assignment");
	  lo = 0;
	} else {
	  lo = static_cast<int> ((xwall[m]-delta-yblo2)*bininv2y);
	  if (lo >= nbin2y) continue;
	  if (lo < 0) error->all(FLERR,"Fix SRD: bad search bin assignment");
	  hi = nbin2y-1;
	}

	for (iy = lo; iy <= hi; iy++)
	  for (ix = 0; ix < nbin2x; ix++)
	    for (iz = 0; iz < nbin2z; iz++) {
	      ibin = iz*nbin2y*nbin2x + iy*nbin2x + ix;
	      if (nbinbig[ibin] == ATOMPERBIN)
		error->all(FLERR,"Fix SRD: too many walls in bin");
	      binbig[ibin][nbinbig[ibin]++] = nbig+m;
	    }

      } else if (dim == 2) {
	if (side == 0) {
	  hi = static_cast<int> ((xwall[m]+delta-zblo2)*bininv2z);
	  if (hi < 0) continue;
	  if (hi >= nbin2z) error->all(FLERR,
				       "Fix SRD: bad search bin assignment");
	  lo = 0;
	} else {
	  lo = static_cast<int> ((xwall[m]-delta-zblo2)*bininv2z);
	  if (lo >= nbin2z) continue;
	  if (lo < 0) error->all(FLERR,"Fix SRD: bad search bin assignment");
	  hi = nbin2z-1;
	}

	for (iz = lo; iz < hi; iz++)
	  for (ix = 0; ix < nbin2x; ix++)
	    for (iy = 0; iy < nbin2y; iy++) {
	      ibin = iz*nbin2y*nbin2x + iy*nbin2x + ix;
	      if (nbinbig[ibin] == ATOMPERBIN)
		error->all(FLERR,"Fix SRD: too many walls in bin");
	      binbig[ibin][nbinbig[ibin]++] = nbig+m;
	    }
      }
    }
  }

  // rotate SRD velocities on SRD timestep
  // done now since all SRDs are currently inside my sub-domain

  if (reneighflag == SRD_ROTATE) reset_velocities();

  // log stats if reneighboring occurred b/c SRDs moved too far

  if (reneighflag == SRD_MOVE) reneighcount++;
  reneighflag = BIG_MOVE;
}

/* ----------------------------------------------------------------------
   advect SRD particles and detect collisions between SRD and BIG particles
   when collision occurs, change x,v of SRD, force,torque of BIG particle
------------------------------------------------------------------------- */

void FixSRD::post_force(int vflag)
{
  int i,m,ix,iy,iz;
  double xlamda[3];

  // zero per-timestep stats

  stats_flag = 0;
  ncheck = ncollide = nbounce = ninside = nrescale = 0;

  // zero ghost forces & torques on BIG particles

  double **f = atom->f;
  double **torque = atom->torque;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  if (bigexist == 0) nall = 0;

  for (i = nlocal; i < nall; i++)
    f[i][0] = f[i][1] = f[i][2] = 0.0;

  if (collidestyle == NOSLIP)
    for (i = nlocal; i < nall; i++)
      torque[i][0] = torque[i][1] = torque[i][2] = 0.0;

  // advect SRD particles
  // assign to search bins if big particles or walls exist

  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;

  if (bigexist || wallexist) {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	x[i][0] += dt_big*v[i][0];
	x[i][1] += dt_big*v[i][1];
	x[i][2] += dt_big*v[i][2];
	
	ix = static_cast<int> ((x[i][0]-xblo2)*bininv2x);
	iy = static_cast<int> ((x[i][1]-yblo2)*bininv2y);
	iz = static_cast<int> ((x[i][2]-zblo2)*bininv2z);
	binsrd[i] = iz*nbin2y*nbin2x + iy*nbin2x + ix;
	
	if (ix < 0 || ix >= nbin2x || iy < 0 || iy >= nbin2y || 
	    iz < 0 || iz >= nbin2z) {
	  if (screen) {
	    fprintf(screen,"SRD particle %d on step " BIGINT_FORMAT "\n",
		    atom->tag[i],update->ntimestep);
	    fprintf(screen,"v = %g %g %g\n",v[i][0],v[i][1],v[i][2]);
	    fprintf(screen,"x = %g %g %g\n",x[i][0],x[i][1],x[i][2]);
	    fprintf(screen,"ix,iy,iz nx,ny,nz = %d %d %d %d %d %d\n",
		    ix,iy,iz,nbin2x,nbin2y,nbin2z);
	  }
	  error->one(FLERR,"Fix SRD: bad bin assignment for SRD advection");
	}
      }

  } else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	x[i][0] += dt_big*v[i][0];
	x[i][1] += dt_big*v[i][1];
	x[i][2] += dt_big*v[i][2];
      }
  }

  // detect collision of SRDs with BIG particles or walls

  if (bigexist || wallexist) {
    if (bigexist) big_dynamic();
    if (wallexist) wallfix->wall_params(0);
    if (overlap) collisions_multi();
    else collisions_single();
  }

  // reverse communicate forces & torques on BIG particles

  if (bigexist) {
    flocal = f;
    tlocal = torque;
    comm->reverse_comm_fix(this);
  }

  // if any SRD particle has moved too far, trigger reneigh on next step
  // for triclinic, perform check in lamda units

  int flag = 0;

  if (triclinic) domain->x2lamda(nlocal);
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (x[i][0] < srdlo_reneigh[0] || x[i][0] > srdhi_reneigh[0] ||
	  x[i][1] < srdlo_reneigh[1] || x[i][1] > srdhi_reneigh[1] ||
	  x[i][2] < srdlo_reneigh[2] || x[i][2] > srdhi_reneigh[2]) flag = 1;
    }
  if (triclinic) domain->lamda2x(nlocal);

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) {
    next_reneighbor = update->ntimestep + 1;
    reneighflag = SRD_MOVE;
  }

  // if wall has moved too far, trigger reneigh on next step
  // analagous to neighbor check for big particle moving 1/2 of skin distance

  if (wallexist) {
    for (m = 0; m < nwall; m++)
      if (fabs(xwall[m]-xwallhold[m]) > walltrigger)
	next_reneighbor = update->ntimestep + 1;
  }

  // if next timestep is SRD timestep, trigger reneigh

  if ((update->ntimestep+1) % nevery == 0) {
    next_reneighbor = update->ntimestep + 1;
    reneighflag = SRD_ROTATE;
  }
}

/* ----------------------------------------------------------------------
   reset SRD velocities
   may perform random shifting by up to 1/2 bin in each dimension
   called at pre-neighbor stage when all SRDs are now inside my sub-domain
   if tstat, then thermostat SRD particles as well, including streaming effects
------------------------------------------------------------------------- */

void FixSRD::reset_velocities()
{
  int i,j,n,ix,iy,iz,ibin,axis,sign,irandom;
  double u[3],vsum[3];
  double vx,vy,vz,vsq,tbin,scale;
  double *vave,*vnew,*xlamda;
  double vstream[3];

  // if requested, perform a dynamic shift of bin positions

  if (shiftflag) {
    double *boxlo;
    if (triclinic == 0) boxlo = domain->boxlo;
    else boxlo = domain->boxlo_lamda;
    shifts[1].corner[0] = boxlo[0] - binsize1x*randomshift->uniform();
    shifts[1].corner[1] = boxlo[1] - binsize1y*randomshift->uniform();
    if (dimension == 3)
      shifts[1].corner[2] = boxlo[2] - binsize1z*randomshift->uniform();
    else shifts[1].corner[2] = boxlo[2];
    setup_velocity_shift(1,1);
  }

  double *corner = shifts[shiftflag].corner;
  int *binlo = shifts[shiftflag].binlo;
  int *binhi = shifts[shiftflag].binhi;
  int nbins = shifts[shiftflag].nbins;
  int nbinx = shifts[shiftflag].nbinx;
  int nbiny = shifts[shiftflag].nbiny;
  BinAve *vbin = shifts[shiftflag].vbin;

  // binhead = 1st SRD particle in each bin
  // binnext = index of next particle in bin
  // bin assignment is done in lamda units for triclinic

  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  int nlocal = atom->nlocal;

  if (triclinic) domain->x2lamda(nlocal);

  for (i = 0; i < nbins; i++) binhead[i] = -1;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ix = static_cast<int> ((x[i][0]-corner[0])*bininv1x);
      ix = MAX(ix,binlo[0]);
      ix = MIN(ix,binhi[0]);
      iy = static_cast<int> ((x[i][1]-corner[1])*bininv1y);
      iy = MAX(iy,binlo[1]);
      iy = MIN(iy,binhi[1]);
      iz = static_cast<int> ((x[i][2]-corner[2])*bininv1z);
      iz = MAX(iz,binlo[2]);
      iz = MIN(iz,binhi[2]);

      ibin = (iz-binlo[2])*nbiny*nbinx + (iy-binlo[1])*nbinx + (ix-binlo[0]);
      binnext[i] = binhead[ibin];
      binhead[ibin] = i;
    }

  if (triclinic) domain->lamda2x(nlocal);

  // for each bin I have particles contributing to:
  // compute summed v of particles in that bin
  // if I own the bin, set its random value, else set to 0.0

  for (i = 0; i < nbins; i++) {
    n = 0;
    vsum[0] = vsum[1] = vsum[2] = 0.0;
    for (j = binhead[i]; j >= 0; j = binnext[j]) {
      vsum[0] += v[j][0];
      vsum[1] += v[j][1];
      vsum[2] += v[j][2];
      n++;
    }

    vbin[i].vsum[0] = vsum[0];
    vbin[i].vsum[1] = vsum[1];
    vbin[i].vsum[2] = vsum[2];
    vbin[i].n = n;
    if (vbin[i].owner) vbin[i].random = random->uniform();
    else vbin[i].random = 0.0;
  }

  // communicate bin info for bins which more than 1 proc contribute to

  if (shifts[shiftflag].commflag) vbin_comm(shiftflag);

  // for each bin I have particles contributing to:
  // compute vave over particles in bin
  // thermal velocity of each particle = v - vave
  // rotate thermal vel of each particle around one of 6 random axes
  // add vave back to each particle
  // thermostat if requested:
  //   if no deformation, rescale thermal vel to temperature
  //   if deformation, rescale thermal vel and change vave to vstream
  //   these are settings for extra dof_temp, dof_tstat to subtract
  //     (not sure why these settings work best)
  //     no deformation, no tstat: dof_temp = 1
  //     yes deformation, no tstat: doesn't matter, system will not equilibrate
  //     no deformation, yes tstat: dof_temp = dof_tstat = 1
  //     yes deformation, yes tstat: dof_temp = dof_tstat = 0
  // accumulate final T_srd for each bin I own

  double tfactor = force->mvv2e * mass_srd / (dimension * force->boltz);
  int dof_temp = 1;
  int dof_tstat;
  if (tstat) {
    if (deformflag) dof_tstat = dof_temp = 0;
    else dof_tstat = 1;
  }

  srd_bin_temp = 0.0;
  srd_bin_count = 0;

  if (dimension == 2) axis = 2;
  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  for (i = 0; i < nbins; i++) {
    n = vbin[i].n;
    if (n == 0) continue;
    vave = vbin[i].vsum;
    vave[0] /= n;
    vave[1] /= n;
    vave[2] /= n;

    irandom = static_cast<int> (6.0*vbin[i].random);
    sign = irandom % 2;
    if (dimension == 3) axis = irandom / 2;

    vsq = 0.0;
    for (j = binhead[i]; j >= 0; j = binnext[j]) {
      if (axis == 0) {
	u[0] = v[j][0]-vave[0];
	u[1] = sign ? v[j][2]-vave[2] : vave[2]-v[j][2];
	u[2] = sign ? vave[1]-v[j][1] : v[j][1]-vave[1];
      } else if (axis == 1) {
	u[1] = v[j][1]-vave[1];
	u[0] = sign ? v[j][2]-vave[2] : vave[2]-v[j][2];
	u[2] = sign ? vave[0]-v[j][0] : v[j][0]-vave[0];
      } else {
	u[2] = v[j][2]-vave[2];
	u[1] = sign ? v[j][0]-vave[0] : vave[0]-v[j][0];
	u[0] = sign ? vave[1]-v[j][1] : v[j][1]-vave[1];
      }
      vsq += u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
      v[j][0] = u[0] + vave[0];
      v[j][1] = u[1] + vave[1];
      v[j][2] = u[2] + vave[2];
    }

    if (tstat && n > 1) {
      if (deformflag) {
	xlamda = vbin[i].xctr;
	vstream[0] = h_rate[0]*xlamda[0] + h_rate[5]*xlamda[1] + 
	  h_rate[4]*xlamda[2] + h_ratelo[0];
	vstream[1] = h_rate[1]*xlamda[1] + h_rate[3]*xlamda[2] + h_ratelo[1];
	vstream[2] = h_rate[2]*xlamda[2] + h_ratelo[2];
      } else {
	vstream[0] = vave[0];
	vstream[1] = vave[1];
	vstream[2] = vave[2];
      }

      // tbin = thermal temperature of particles in bin
      // scale = scale factor for thermal velocity

      tbin = vsq/(n-dof_tstat) * tfactor;
      scale = sqrt(temperature_srd/tbin);

      vsq = 0.0;
      for (j = binhead[i]; j >= 0; j = binnext[j]) {
	u[0] = (v[j][0] - vave[0]) * scale;
	u[1] = (v[j][1] - vave[1]) * scale;
	u[2] = (v[j][2] - vave[2]) * scale;
	vsq += u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
	v[j][0] = u[0] + vstream[0];
	v[j][1] = u[1] + vstream[1];
	v[j][2] = u[2] + vstream[2];
      }
    }

    // sum partial contribution of my particles to T even if I don't own bin
    // but only count bin if I own it, so each bin is counted exactly once

    if (n > 1) srd_bin_temp += vsq/(n-dof_temp);
    if (vbin[i].owner) srd_bin_count++;
  }

  srd_bin_temp *= tfactor;

  // rescale any too-large velocities 

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > vmaxsq) {
	nrescale++;
	MathExtra::scale3(vmax/sqrt(vsq),v[i]);
      }
    }
}

/* ----------------------------------------------------------------------
   communicate summed particle info for bins that overlap 1 or more procs
------------------------------------------------------------------------- */

void FixSRD::vbin_comm(int ishift)
{
  BinComm *bcomm1,*bcomm2;
  MPI_Request request1,request2;
  MPI_Status status;
  
  // send/recv bins in both directions in each dimension
  // don't send if nsend = 0
  //   due to static bins aliging with proc boundary
  //   due to dynamic bins across non-periodic global boundary
  // copy to self if sendproc = me
  // MPI send to another proc if sendproc != me
  // don't recv if nrecv = 0
  // copy from self if recvproc = me
  // MPI recv from another proc if recvproc != me
  
  BinAve *vbin = shifts[ishift].vbin;
  int *procgrid = comm->procgrid;

  int iswap = 0;
  for (int idim = 0; idim < dimension; idim++) {
    bcomm1 = &shifts[ishift].bcomm[iswap++];
    bcomm2 = &shifts[ishift].bcomm[iswap++];
    
    if (procgrid[idim] == 1) {
      if (bcomm1->nsend)
	vbin_pack(vbin,bcomm1->nsend,bcomm1->sendlist,sbuf1);
      if (bcomm2->nsend)
	vbin_pack(vbin,bcomm2->nsend,bcomm2->sendlist,sbuf2);
      if (bcomm1->nrecv)
	vbin_unpack(sbuf1,vbin,bcomm1->nrecv,bcomm1->recvlist);
      if (bcomm2->nrecv)
	vbin_unpack(sbuf2,vbin,bcomm2->nrecv,bcomm2->recvlist);

    } else {
      if (bcomm1->nrecv)
	MPI_Irecv(rbuf1,bcomm1->nrecv*VBINSIZE,MPI_DOUBLE,bcomm1->recvproc,0,
		  world,&request1);
      if (bcomm2->nrecv)
	MPI_Irecv(rbuf2,bcomm2->nrecv*VBINSIZE,MPI_DOUBLE,bcomm2->recvproc,0,
		  world,&request2);
      if (bcomm1->nsend) {
	vbin_pack(vbin,bcomm1->nsend,bcomm1->sendlist,sbuf1);
	MPI_Send(sbuf1,bcomm1->nsend*VBINSIZE,MPI_DOUBLE,
		 bcomm1->sendproc,0,world);
      }
      if (bcomm2->nsend) {
	vbin_pack(vbin,bcomm2->nsend,bcomm2->sendlist,sbuf2);
	MPI_Send(sbuf2,bcomm2->nsend*VBINSIZE,MPI_DOUBLE,
		 bcomm2->sendproc,0,world);
      }
      if (bcomm1->nrecv) {
	MPI_Wait(&request1,&status);
	vbin_unpack(rbuf1,vbin,bcomm1->nrecv,bcomm1->recvlist);
      }
      if (bcomm2->nrecv) {
	MPI_Wait(&request2,&status);
	vbin_unpack(rbuf2,vbin,bcomm2->nrecv,bcomm2->recvlist);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   pack velocity bin data into a message buffer for sending
------------------------------------------------------------------------- */

void FixSRD::vbin_pack(BinAve *vbin, int n, int *list, double *buf)
{
  int j;
  int m = 0;
  for (int i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = vbin[j].n;
    buf[m++] = vbin[j].vsum[0];
    buf[m++] = vbin[j].vsum[1];
    buf[m++] = vbin[j].vsum[2];
    buf[m++] = vbin[j].random;
  }
}

/* ----------------------------------------------------------------------
   unpack velocity bin data from a message buffer and sum values to my bins
------------------------------------------------------------------------- */

void FixSRD::vbin_unpack(double *buf, BinAve *vbin, int n, int *list)
{
  int j;
  int m = 0;
  for (int i = 0; i < n; i++) {
    j = list[i];
    vbin[j].n += static_cast<int> (buf[m++]);
    vbin[j].vsum[0] += buf[m++];
    vbin[j].vsum[1] += buf[m++];
    vbin[j].vsum[2] += buf[m++];
    vbin[j].random += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   detect all collisions between SRD and BIG particles or WALLS
   assume SRD can be inside at most one BIG particle or WALL at a time
   unoverlap SRDs for each collision
------------------------------------------------------------------------- */

void FixSRD::collisions_single()
{
  int i,j,k,m,type,nbig,ibin,ibounce,inside,collide_flag,lineside;
  double dt,t_remain;
  double norm[3],xscoll[3],xbcoll[3],vsnew[3];
  Big *big;

  // outer loop over SRD particles
  // inner loop over BIG particles or WALLS that overlap SRD particle bin
  // if overlap between SRD and BIG particle or wall:
  //   for exact, compute collision pt in time
  //   for inexact, push SRD to surf of BIG particle or WALL
  // update x,v of SRD and f,torque on BIG particle
  // re-bin SRD particle after collision
  // iterate until the SRD particle has no overlaps with BIG particles or WALLS

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    ibin = binsrd[i];
    if (nbinbig[ibin] == 0) continue;

    ibounce = 0;
    collide_flag = 1;
    dt = dt_big;

    while (collide_flag) {
      nbig = nbinbig[ibin];
      if (ibounce == 0) ncheck += nbig;

      collide_flag = 0;
      for (m = 0; m < nbig; m++) {
	k = binbig[ibin][m];
	big = &biglist[k];
	j = big->index;
	type = big->type;

	if (type == SPHERE) inside = inside_sphere(x[i],x[j],big);
	else if (type == ELLIPSOID) inside = inside_ellipsoid(x[i],x[j],big);
	else inside = inside_wall(x[i],j);

	if (inside) {
	  if (exactflag) {
	    if (type == SPHERE)
	      t_remain = collision_sphere_exact(x[i],x[j],v[i],v[j],big,
						xscoll,xbcoll,norm);
	    else if (type == ELLIPSOID)
	      t_remain = collision_ellipsoid_exact(x[i],x[j],v[i],v[j],big,
						   xscoll,xbcoll,norm);
	    else 
	      t_remain = collision_wall_exact(x[i],j,v[i],xscoll,xbcoll,norm);

	  } else {
	    t_remain = 0.5*dt;
	    if (type == SPHERE)
	      collision_sphere_inexact(x[i],x[j],big,xscoll,xbcoll,norm);
	    else if (type == ELLIPSOID)
	      collision_ellipsoid_inexact(x[i],x[j],big,xscoll,xbcoll,norm);
	    else
	      collision_wall_inexact(x[i],j,xscoll,xbcoll,norm);
	  }

#ifdef SRD_DEBUG
	  if (update->ntimestep == SRD_DEBUG_TIMESTEP &&
	      atom->tag[i] == SRD_DEBUG_ATOMID)
	    print_collision(i,j,ibounce,t_remain,dt,xscoll,xbcoll,norm,type);
#endif

	  if (t_remain > dt) {
	    ninside++;
	    if (insideflag == INSIDE_ERROR || insideflag == INSIDE_WARN) {
	      char str[128];
	      if (type != WALL)
		sprintf(str,
			"SRD particle %d started "
			"inside big particle %d on step " BIGINT_FORMAT 
			" bounce %d",
			atom->tag[i],atom->tag[j],update->ntimestep,ibounce+1);
	      else
		sprintf(str,
			"SRD particle %d started "
			"inside big particle %d on step " BIGINT_FORMAT 
			" bounce %d",
			atom->tag[i],atom->tag[j],update->ntimestep,ibounce+1);
	      if (insideflag == INSIDE_ERROR) error->one(FLERR,str);
	      error->warning(FLERR,str);
	    }
	    break;
	  }

	  if (collidestyle == SLIP) {
	    if (type != WALL) slip(v[i],v[j],x[j],big,xscoll,norm,vsnew);
	    else slip_wall(v[i],j,norm,vsnew);
	  } else {
	    if (type != WALL) noslip(v[i],v[j],x[j],big,-1, xscoll,norm,vsnew);
	    else noslip(v[i],NULL,x[j],big,j,xscoll,norm,vsnew);
	  }

	  if (dimension == 2) vsnew[2] = 0.0;

	  // check on rescaling of vsnew

	  double vsq = vsnew[0]*vsnew[0] + vsnew[1]*vsnew[1] + 
	    vsnew[2]*vsnew[2];
	  if (vsq > vmaxsq) {
	    nrescale++;
	    MathExtra::scale3(vmax/sqrt(vsq),vsnew);
	  }

	  // update BIG particle and WALL and SRD
	  // BIG particle is not torqued if sphere and SLIP collision

	  if (collidestyle == SLIP && type == SPHERE)
	    force_torque(v[i],vsnew,xscoll,xbcoll,f[j],NULL);
	  else if (type != WALL)
	    force_torque(v[i],vsnew,xscoll,xbcoll,f[j],torque[j]);
	  else if (type == WALL)
	    force_wall(v[i],vsnew,j);

	  ibin = binsrd[i] = update_srd(i,t_remain,xscoll,vsnew,x[i],v[i]);

	  if (ibounce == 0) ncollide++;
	  ibounce++;
	  if (ibounce < maxbounceallow || maxbounceallow == 0)
	    collide_flag = 1;
	  dt = t_remain;
	  break;
	}
      }
    }

    nbounce += ibounce;
    if (maxbounceallow && ibounce >= maxbounceallow) bouncemaxnum++;
    if (ibounce > bouncemax) bouncemax = ibounce;
  }
}

/* ----------------------------------------------------------------------
   detect all collisions between SRD and big particles
   an SRD can be inside more than one big particle at a time
   requires finding which big particle SRD collided with first
   unoverlap SRDs for each collision
------------------------------------------------------------------------- */

void FixSRD::collisions_multi()
{
  int i,j,k,m,type,nbig,ibin,ibounce,inside,jfirst,typefirst,jlast;
  double dt,t_remain,t_first;
  double norm[3],xscoll[3],xbcoll[3],vsnew[3];
  double normfirst[3],xscollfirst[3],xbcollfirst[3];
  Big *big;

  // outer loop over SRD particles
  // inner loop over BIG particles or WALLS that overlap SRD particle bin
  // loop over all BIG and WALLS to find which one SRD collided with first
  // if overlap between SRD and BIG particle or wall:
  //   compute collision pt in time
  // update x,v of SRD and f,torque on BIG particle
  // re-bin SRD particle after collision
  // iterate until the SRD particle has no overlaps with BIG particles or WALLS

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    ibin = binsrd[i];
    if (nbinbig[ibin] == 0) continue;

    ibounce = 0;
    jlast = -1;
    dt = dt_big;

    while (1) {
      nbig = nbinbig[ibin];
      if (ibounce == 0) ncheck += nbig;

      t_first = 0.0;
      for (m = 0; m < nbig; m++) {
	k = binbig[ibin][m];
	big = &biglist[k];
	j = big->index;
	if (j == jlast) continue;
	type = big->type;

	if (type == SPHERE)
	  inside = inside_sphere(x[i],x[j],big);
	else if (type == ELLIPSOID)
	  inside = inside_ellipsoid(x[i],x[j],big);
	else if (type == LINE)
	  inside = inside_line(x[i],x[j],v[i],v[j],big,dt);
	else if (type == TRIANGLE) 
	  inside = inside_tri(x[i],x[j],v[i],v[j],big,dt);
	else
	  inside = inside_wall(x[i],j);

	if (inside) {
	  if (type == SPHERE)
	    t_remain = collision_sphere_exact(x[i],x[j],v[i],v[j],big,
					      xscoll,xbcoll,norm);
	  else if (type == ELLIPSOID)
	    t_remain = collision_ellipsoid_exact(x[i],x[j],v[i],v[j],big,
						 xscoll,xbcoll,norm);
	  else if (type == LINE)
	    t_remain = collision_line_exact(x[i],x[j],v[i],v[j],big,dt,
					    xscoll,xbcoll,norm);
	  else if (type == TRIANGLE)
	    t_remain = collision_tri_exact(x[i],x[j],v[i],v[j],big,dt,
					   xscoll,xbcoll,norm);
	  else 
	    t_remain = collision_wall_exact(x[i],j,v[i],xscoll,xbcoll,norm);

#ifdef SRD_DEBUG
	  if (update->ntimestep == SRD_DEBUG_TIMESTEP &&
	      atom->tag[i] == SRD_DEBUG_ATOMID)
	    print_collision(i,j,ibounce,t_remain,dt,xscoll,xbcoll,norm,type);
#endif

	  if (t_remain > dt || t_remain < 0.0) {
	    ninside++;
	    if (insideflag == INSIDE_ERROR || insideflag == INSIDE_WARN) {
	      char str[128];
	      sprintf(str,
		      "SRD particle %d started "
		      "inside big particle %d on step " BIGINT_FORMAT 
		      " bounce %d",
		      atom->tag[i],atom->tag[j],update->ntimestep,ibounce+1);
	      if (insideflag == INSIDE_ERROR) error->one(FLERR,str);
	      error->warning(FLERR,str);
	    }
	    t_first = 0.0;
	    break;
	  }

	  if (t_remain > t_first) {
	    t_first = t_remain;
	    jfirst = j;
	    typefirst = type;
	    xscollfirst[0] = xscoll[0];
	    xscollfirst[1] = xscoll[1];
	    xscollfirst[2] = xscoll[2];
	    xbcollfirst[0] = xbcoll[0];
	    xbcollfirst[1] = xbcoll[1];
	    xbcollfirst[2] = xbcoll[2];
	    normfirst[0] = norm[0];
	    normfirst[1] = norm[1];
	    normfirst[2] = norm[2];
	  }
	}
      }

      if (t_first == 0.0) break;
      j = jlast = jfirst;
      type = typefirst;
      xscoll[0] = xscollfirst[0];
      xscoll[1] = xscollfirst[1];
      xscoll[2] = xscollfirst[2];
      xbcoll[0] = xbcollfirst[0];
      xbcoll[1] = xbcollfirst[1];
      xbcoll[2] = xbcollfirst[2];
      norm[0] = normfirst[0];
      norm[1] = normfirst[1];
      norm[2] = normfirst[2];

      if (collidestyle == SLIP) {
	if (type != WALL) slip(v[i],v[j],x[j],big,xscoll,norm,vsnew);
	else slip_wall(v[i],j,norm,vsnew);
      } else {
	if (type != WALL) noslip(v[i],v[j],x[j],big,-1,xscoll,norm,vsnew);
	else noslip(v[i],NULL,x[j],big,j,xscoll,norm,vsnew);
      }

      if (dimension == 2) vsnew[2] = 0.0;

      // check on rescaling of vsnew

      double vsq = vsnew[0]*vsnew[0] + vsnew[1]*vsnew[1] + vsnew[2]*vsnew[2];
      if (vsq > vmaxsq) {
	nrescale++;
	MathExtra::scale3(vmax/sqrt(vsq),vsnew);
      }

      // update BIG particle and WALL and SRD
      // BIG particle is not torqued if sphere and SLIP collision

      if (collidestyle == SLIP && type == SPHERE)
	force_torque(v[i],vsnew,xscoll,xbcoll,f[j],NULL);
      else if (type != WALL)
	force_torque(v[i],vsnew,xscoll,xbcoll,f[j],torque[j]);
      else if (type == WALL)
	force_wall(v[i],vsnew,j);

      ibin = binsrd[i] = update_srd(i,t_first,xscoll,vsnew,x[i],v[i]);

      if (ibounce == 0) ncollide++;
      ibounce++;
      if (ibounce == maxbounceallow) break;
      dt = t_first;
    }

    nbounce += ibounce;
    if (maxbounceallow && ibounce >= maxbounceallow) bouncemaxnum++;
    if (ibounce > bouncemax) bouncemax = ibounce;
  }
}

/* ----------------------------------------------------------------------
   check if SRD particle S is inside spherical big particle B
------------------------------------------------------------------------- */

int FixSRD::inside_sphere(double *xs, double *xb, Big *big)
{
  double dx,dy,dz;

  dx = xs[0] - xb[0];
  dy = xs[1] - xb[1];
  dz = xs[2] - xb[2];

  if (dx*dx + dy*dy + dz*dz <= big->radsq) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if SRD particle S is inside ellipsoidal big particle B
------------------------------------------------------------------------- */

int FixSRD::inside_ellipsoid(double *xs, double *xb, Big *big)
{
  double x,y,z;

  double *ex = big->ex;
  double *ey = big->ey;
  double *ez = big->ez;

  double xs_xb[3];
  xs_xb[0] = xs[0] - xb[0];
  xs_xb[1] = xs[1] - xb[1];
  xs_xb[2] = xs[2] - xb[2];

  x = xs_xb[0]*ex[0] + xs_xb[1]*ex[1] + xs_xb[2]*ex[2];
  y = xs_xb[0]*ey[0] + xs_xb[1]*ey[1] + xs_xb[2]*ey[2];
  z = xs_xb[0]*ez[0] + xs_xb[1]*ez[1] + xs_xb[2]*ez[2];

  if (x*x*big->aradsqinv + y*y*big->bradsqinv + z*z*big->cradsqinv <= 1.0)
    return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if SRD particle S is inside line big particle B
   collision only possible if:
     S starts on positive side of infinite line,
       which means it will collide with outside of rigid body made of lines
       since line segments have outward normals,
       when vector from first to last point is crossed with +z axis
     S ends on negative side of infinite line
   unlike most other inside() routines, then calculate exact collision:
     solve for collision pt along infinite line
     collision if pt is within endpoints of B
------------------------------------------------------------------------- */

int FixSRD::inside_line(double *xs, double *xb, double *vs, double *vb,
			Big *big, double dt_step)
{
  double pmc0[2],pmc1[2],n0[2],n1[2];
  double n1_n0[2],pmc1_pmc0[2];

  // 1 and 2 = start and end of timestep
  // pmc = P - C, where P = position of S, C = position of B
  // n = normal to line = [-sin(theta),cos(theta)], theta = orientation of B
  // (P-C) dot N = side of line that S is on
  // side0 = -1,0,1 for which side of line B that S is on at start of step
  // side1 = -1,0,1 for which side of line B that S is on at end of step

  xs1[0] = xs[0];
  xs1[1] = xs[1];
  xb1[0] = xb[0];
  xb1[1] = xb[1];

  xs0[0] = xs1[0] - dt_step*vs[0];
  xs0[1] = xs1[1] - dt_step*vs[1];
  xb0[0] = xb1[0] - dt_step*vb[0];
  xb0[1] = xb1[1] - dt_step*vb[1];

  theta1 = big->theta;
  theta0 = theta1 - dt_step*big->omega[2];

  pmc0[0] = xs0[0] - xb0[0];
  pmc0[1] = xs0[1] - xb0[1];
  n0[0] = sin(theta0);
  n0[1] = -cos(theta0);

  pmc1[0] = xs1[0] - xb1[0];
  pmc1[1] = xs1[1] - xb1[1];
  n1[0] = sin(theta1);
  n1[1] = -cos(theta1);

  double side0 = pmc0[0]*n0[0] + pmc0[1]*n0[1];
  double side1 = pmc1[0]*n1[0] + pmc1[1]*n1[1];

  if (side0 <= 0.0 || side1 >= 0.0) return 0;

  // solve for time t (0 to 1) at which moving particle
  //   crosses infinite moving/rotating line

  // Newton-Raphson solve of full non-linear parametric equation

  tfraction = newton_raphson(0.0,1.0);

  // quadratic equation solve of approximate parametric equation

  /*
  n1_n0[0] = n1[0]-n0[0]; n1_n0[1] = n1[1]-n0[1];
  pmc1_pmc0[0] = pmc1[0]-pmc0[0]; pmc1_pmc0[1] = pmc1[1]-pmc0[1];

  double a = pmc1_pmc0[0]*n1_n0[0] + pmc1_pmc0[1]*n1_n0[1];
  double b = pmc1_pmc0[0]*n0[0] + pmc1_pmc0[1]*n0[1] + 
  n1_n0[0]*pmc0[0] + n1_n0[1]*pmc0[1];
  double c = pmc0[0]*n0[0] + pmc0[1]*n0[1];
  
  if (a == 0.0) {
    double dot0 = pmc0[0]*n0[0] + pmc0[1]*n0[1];
    double dot1 = pmc1[0]*n0[0] + pmc1[1]*n0[1];
    double root = -dot0 / (dot1 - dot0);
    //printf("Linear root: %g %g\n",root,tfraction);
    tfraction = root;

  } else {

    double term = sqrt(b*b - 4.0*a*c);
    double root1 = (-b + term) / (2.0*a);
    double root2 = (-b - term) / (2.0*a);

    //printf("ABC vecs: %g %g: %g %g\n",
    //	   pmc1_pmc0[0],pmc1_pmc0[1],n1_n0[0],n1_n0[1]);
    //printf("ABC vecs: %g %g: %g %g: %g %g %g\n",
    //	   n0[0],n0[1],n1[0],n1[1],theta0,theta1,big->omega[2]);
    //printf("ABC root: %g %g %g: %g %g %g\n",a,b,c,root1,root2,tfraction);
    
    if (0.0 <= root1 && root1 <= 1.0) tfraction = root1;
    else if (0.0 <= root2 && root2 <= 1.0) tfraction = root2;
    else error->one(FLERR,"Bad quadratic solve for particle/line collision");
  }
  */

  // check if collision pt is within line segment at collision time

  xsc[0] = xs0[0] + tfraction*(xs1[0]-xs0[0]);
  xsc[1] = xs0[1] + tfraction*(xs1[1]-xs0[1]);
  xbc[0] = xb0[0] + tfraction*(xb1[0]-xb0[0]);
  xbc[1] = xb0[1] + tfraction*(xb1[1]-xb0[1]);
  double delx = xsc[0] - xbc[0];
  double dely = xsc[1] - xbc[1];
  double rsq = delx*delx + dely*dely;
  if (rsq > 0.25*big->length*big->length) return 0;

  //nbc[0] = n0[0] + tfraction*(n1[0]-n0[0]);
  //nbc[1] = n0[1] + tfraction*(n1[1]-n0[1]);

  nbc[0] = sin(theta0 + tfraction*(theta1-theta0));
  nbc[1] = -cos(theta0 + tfraction*(theta1-theta0));

  return 1;
}

/* ----------------------------------------------------------------------
   check if SRD particle S is inside triangle big particle B
   collision only possible if:
     S starts on positive side of triangle plane,
       which means it will collide with outside of rigid body made of tris
       since triangles have outward normals,
     S ends on negative side of triangle plane
   unlike most other inside() routines, then calculate exact collision:
     solve for collision pt on triangle plane
     collision if pt is inside triangle B
------------------------------------------------------------------------- */

int FixSRD::inside_tri(double *xs, double *xb, double *vs, double *vb,
		       Big *big, double dt_step)
{
  double pmc0[3],pmc1[3],n0[3];
  double n1_n0[3],pmc1_pmc0[3];
  double excoll[3],eycoll[3],ezcoll[3];
  double dc1[3],dc2[3],dc3[3];
  double c1[3],c2[3],c3[3];
  double c2mc1[3],c3mc2[3],c1mc3[3];
  double pvec[3],xproduct[3];

  // 1 and 2 = start and end of timestep
  // pmc = P - C, where P = position of S, C = position of B
  // n = normal to triangle
  // (P-C) dot N = side of tri that S is on
  // side0 = -1,0,1 for which side of tri B that S is on at start of step
  // side1 = -1,0,1 for which side of tri B that S is on at end of step

  double *omega = big->omega;
  double *n1 = big->norm;

  n0[0] = n1[0] - dt_step*(omega[1]*n1[2] - omega[2]*n1[1]);
  n0[1] = n1[1] - dt_step*(omega[2]*n1[0] - omega[0]*n1[2]);
  n0[2] = n1[2] - dt_step*(omega[0]*n1[1] - omega[1]*n1[0]);

  pmc0[0] = xs[0] - dt_step*vs[0] - xb[0] + dt_step*vb[0];
  pmc0[1] = xs[1] - dt_step*vs[1] - xb[1] + dt_step*vb[1];
  pmc0[2] = xs[2] - dt_step*vs[2] - xb[2] + dt_step*vb[2];
  pmc1[0] = xs[0] - xb[0];
  pmc1[1] = xs[1] - xb[1];
  pmc1[2] = xs[2] - xb[2];

  double side0 = MathExtra::dot3(pmc0,n0);
  double side1 = MathExtra::dot3(pmc1,n1);

  if (side0 <= 0.0 || side1 >= 0.0) return 0;

  // solve for time t (0 to 1) at which moving particle
  //   crosses moving/rotating tri
  // quadratic equation solve of approximate parametric equation

  n1_n0[0] = n1[0]-n0[0];
  n1_n0[1] = n1[1]-n0[1];
  n1_n0[2] = n1[2]-n0[2];
  pmc1_pmc0[0] = pmc1[0]-pmc0[0];
  pmc1_pmc0[1] = pmc1[1]-pmc0[1];
  pmc1_pmc0[2] = pmc1[2]-pmc0[2];
  
  double a = MathExtra::dot3(pmc1_pmc0,n1_n0);
  double b = MathExtra::dot3(pmc1_pmc0,n0) + MathExtra::dot3(n1_n0,pmc0);
  double c = MathExtra::dot3(pmc0,n0);

  if (a == 0.0) {
    double dot0 = MathExtra::dot3(pmc0,n0);
    double dot1 = MathExtra::dot3(pmc1,n0);
    double root = -dot0 / (dot1 - dot0);
    tfraction = root;
  } else {
    double term = sqrt(b*b - 4.0*a*c);
    double root1 = (-b + term) / (2.0*a);
    double root2 = (-b - term) / (2.0*a);
    if (0.0 <= root1 && root1 <= 1.0) tfraction = root1;
    else if (0.0 <= root2 && root2 <= 1.0) tfraction = root2;
    else error->one(FLERR,"Bad quadratic solve for particle/tri collision");
  }

  // calculate position/orientation of S and B at collision time
  // dt = time previous to now at which collision occurs
  // point = S position in plane of triangle at collision time
  // Excoll,Eycoll,Ezcoll = orientation of tri at collision time
  // c1,c2,c3 = corner points of triangle at collision time
  // nbc = normal to plane of triangle at collision time

  AtomVecTri::Bonus *tbonus;
  tbonus = avec_tri->bonus;

  double *ex = big->ex;
  double *ey = big->ey;
  double *ez = big->ez;
  int index = atom->tri[big->index];
  double *c1body = tbonus[index].c1;
  double *c2body = tbonus[index].c2;
  double *c3body = tbonus[index].c3;

  double dt = (1.0-tfraction)*dt_step;

  xsc[0] = xs[0] - dt*vs[0];
  xsc[1] = xs[1] - dt*vs[1];
  xsc[2] = xs[2] - dt*vs[2];
  xbc[0] = xb[0] - dt*vb[0];
  xbc[1] = xb[1] - dt*vb[1];
  xbc[2] = xb[2] - dt*vb[2];

  excoll[0] = ex[0] - dt*(omega[1]*ex[2] - omega[2]*ex[1]);
  excoll[1] = ex[1] - dt*(omega[2]*ex[0] - omega[0]*ex[2]);
  excoll[2] = ex[2] - dt*(omega[0]*ex[1] - omega[1]*ex[0]);

  eycoll[0] = ey[0] - dt*(omega[1]*ey[2] - omega[2]*ey[1]);
  eycoll[1] = ey[1] - dt*(omega[2]*ey[0] - omega[0]*ey[2]);
  eycoll[2] = ey[2] - dt*(omega[0]*ey[1] - omega[1]*ey[0]);

  ezcoll[0] = ez[0] - dt*(omega[1]*ez[2] - omega[2]*ez[1]);
  ezcoll[1] = ez[1] - dt*(omega[2]*ez[0] - omega[0]*ez[2]);
  ezcoll[2] = ez[2] - dt*(omega[0]*ez[1] - omega[1]*ez[0]);

  MathExtra::matvec(excoll,eycoll,ezcoll,c1body,dc1);
  MathExtra::matvec(excoll,eycoll,ezcoll,c2body,dc2);
  MathExtra::matvec(excoll,eycoll,ezcoll,c3body,dc3);

  MathExtra::add3(xbc,dc1,c1);
  MathExtra::add3(xbc,dc2,c2);
  MathExtra::add3(xbc,dc3,c3);

  MathExtra::sub3(c2,c1,c2mc1);
  MathExtra::sub3(c3,c2,c3mc2);
  MathExtra::sub3(c1,c3,c1mc3);

  MathExtra::cross3(c2mc1,c3mc2,nbc);
  MathExtra::norm3(nbc);

  // check if collision pt is within triangle
  // pvec = vector from tri vertex to intersection point
  // xproduct = cross product of edge vec with pvec
  // if dot product of xproduct with nbc < 0.0 for any of 3 edges,
  //   intersection point is outside tri

  MathExtra::sub3(xsc,c1,pvec);
  MathExtra::cross3(c2mc1,pvec,xproduct);
  if (MathExtra::dot3(xproduct,nbc) < 0.0) return 0;

  MathExtra::sub3(xsc,c2,pvec);
  MathExtra::cross3(c3mc2,pvec,xproduct);
  if (MathExtra::dot3(xproduct,nbc) < 0.0) return 0;

  MathExtra::sub3(xsc,c3,pvec);
  MathExtra::cross3(c1mc3,pvec,xproduct);
  if (MathExtra::dot3(xproduct,nbc) < 0.0) return 0;

  return 1;
}

/* ----------------------------------------------------------------------
   check if SRD particle S is inside wall IWALL
------------------------------------------------------------------------- */

int FixSRD::inside_wall(double *xs, int iwall)
{
  int dim = wallwhich[iwall] / 2;
  int side = wallwhich[iwall] % 2;

  if (side == 0 && xs[dim] < xwall[iwall]) return 1;
  if (side && xs[dim] > xwall[iwall]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   collision of SRD particle S with surface of spherical big particle B
   exact because compute time of collision
   dt = time previous to now at which collision occurs
   xscoll = collision pt = position of SRD at time of collision
   xbcoll = position of big particle at time of collision
   norm = surface normal of collision pt at time of collision
------------------------------------------------------------------------- */

double FixSRD::collision_sphere_exact(double *xs, double *xb,
				      double *vs, double *vb, Big *big,
				      double *xscoll, double *xbcoll,
				      double *norm)
{
  double vs_dot_vs,vb_dot_vb,vs_dot_vb;
  double vs_dot_xb,vb_dot_xs,vs_dot_xs,vb_dot_xb;
  double xs_dot_xs,xb_dot_xb,xs_dot_xb;
  double a,b,c,scale;

  vs_dot_vs = vs[0]*vs[0] + vs[1]*vs[1] + vs[2]*vs[2];
  vb_dot_vb = vb[0]*vb[0] + vb[1]*vb[1] + vb[2]*vb[2];
  vs_dot_vb = vs[0]*vb[0] + vs[1]*vb[1] + vs[2]*vb[2];

  vs_dot_xb = vs[0]*xb[0] + vs[1]*xb[1] + vs[2]*xb[2];
  vb_dot_xs = vb[0]*xs[0] + vb[1]*xs[1] + vb[2]*xs[2];
  vs_dot_xs = vs[0]*xs[0] + vs[1]*xs[1] + vs[2]*xs[2];
  vb_dot_xb = vb[0]*xb[0] + vb[1]*xb[1] + vb[2]*xb[2];

  xs_dot_xs = xs[0]*xs[0] + xs[1]*xs[1] + xs[2]*xs[2];
  xb_dot_xb = xb[0]*xb[0] + xb[1]*xb[1] + xb[2]*xb[2];
  xs_dot_xb = xs[0]*xb[0] + xs[1]*xb[1] + xs[2]*xb[2];

  a = vs_dot_vs + vb_dot_vb - 2.0*vs_dot_vb;
  b = 2.0 * (vs_dot_xb + vb_dot_xs - vs_dot_xs - vb_dot_xb);
  c = xs_dot_xs + xb_dot_xb - 2.0*xs_dot_xb - big->radsq;

  double dt = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);

  xscoll[0] = xs[0] - dt*vs[0];
  xscoll[1] = xs[1] - dt*vs[1];
  xscoll[2] = xs[2] - dt*vs[2];

  xbcoll[0] = xb[0] - dt*vb[0];
  xbcoll[1] = xb[1] - dt*vb[1];
  xbcoll[2] = xb[2] - dt*vb[2];

  norm[0] = xscoll[0] - xbcoll[0];
  norm[1] = xscoll[1] - xbcoll[1];
  norm[2] = xscoll[2] - xbcoll[2];
  scale = 1.0/sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
  norm[0] *= scale;
  norm[1] *= scale;
  norm[2] *= scale;

  return dt;
}

/* ----------------------------------------------------------------------
   collision of SRD particle S with surface of spherical big particle B
   inexact because just push SRD to surface of big particle at end of step
   time of collision = end of step
   xscoll = collision pt = position of SRD at time of collision
   xbcoll = xb = position of big particle at time of collision
   norm = surface normal of collision pt at time of collision
------------------------------------------------------------------------- */

void FixSRD::collision_sphere_inexact(double *xs, double *xb,
				      Big *big,
				      double *xscoll, double *xbcoll,
				      double *norm)
{
  double scale;

  norm[0] = xs[0] - xb[0];
  norm[1] = xs[1] - xb[1];
  norm[2] = xs[2] - xb[2];
  scale = 1.0/sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
  norm[0] *= scale;
  norm[1] *= scale;
  norm[2] *= scale;

  xscoll[0] = xb[0] + big->radius*norm[0];
  xscoll[1] = xb[1] + big->radius*norm[1];
  xscoll[2] = xb[2] + big->radius*norm[2];

  xbcoll[0] = xb[0];
  xbcoll[1] = xb[1];
  xbcoll[2] = xb[2];
}

/* ----------------------------------------------------------------------
   collision of SRD particle S with surface of ellipsoidal big particle B
   exact because compute time of collision
   dt = time previous to now at which collision occurs
   xscoll = collision pt = position of SRD at time of collision
   xbcoll = position of big particle at time of collision
   norm = surface normal of collision pt at time of collision
------------------------------------------------------------------------- */

double FixSRD::collision_ellipsoid_exact(double *xs, double *xb,
					 double *vs, double *vb, Big *big,
					 double *xscoll, double *xbcoll,
					 double *norm)
{
  double vs_vb[3],xs_xb[3],omega_ex[3],omega_ey[3],omega_ez[3];
  double excoll[3],eycoll[3],ezcoll[3],delta[3],xbody[3],nbody[3];
  double ax,bx,cx,ay,by,cy,az,bz,cz;
  double a,b,c,scale;

  double *omega = big->omega;
  double *ex = big->ex;
  double *ey = big->ey;
  double *ez = big->ez;

  vs_vb[0] = vs[0]-vb[0]; vs_vb[1] = vs[1]-vb[1]; vs_vb[2] = vs[2]-vb[2];
  xs_xb[0] = xs[0]-xb[0]; xs_xb[1] = xs[1]-xb[1]; xs_xb[2] = xs[2]-xb[2];

  MathExtra::cross3(omega,ex,omega_ex);
  MathExtra::cross3(omega,ey,omega_ey);
  MathExtra::cross3(omega,ez,omega_ez);

  ax = vs_vb[0]*omega_ex[0] + vs_vb[1]*omega_ex[1] + vs_vb[2]*omega_ex[2];
  bx = -(vs_vb[0]*ex[0] + vs_vb[1]*ex[1] + vs_vb[2]*ex[2]);
  bx -= xs_xb[0]*omega_ex[0] + xs_xb[1]*omega_ex[1] + xs_xb[2]*omega_ex[2];
  cx = xs_xb[0]*ex[0] + xs_xb[1]*ex[1] + xs_xb[2]*ex[2];

  ay = vs_vb[0]*omega_ey[0] + vs_vb[1]*omega_ey[1] + vs_vb[2]*omega_ey[2];
  by = -(vs_vb[0]*ey[0] + vs_vb[1]*ey[1] + vs_vb[2]*ey[2]);
  by -= xs_xb[0]*omega_ey[0] + xs_xb[1]*omega_ey[1] + xs_xb[2]*omega_ey[2];
  cy = xs_xb[0]*ey[0] + xs_xb[1]*ey[1] + xs_xb[2]*ey[2];

  az = vs_vb[0]*omega_ez[0] + vs_vb[1]*omega_ez[1] + vs_vb[2]*omega_ez[2];
  bz = -(vs_vb[0]*ez[0] + vs_vb[1]*ez[1] + vs_vb[2]*ez[2]);
  bz -= xs_xb[0]*omega_ez[0] + xs_xb[1]*omega_ez[1] + xs_xb[2]*omega_ez[2];
  cz = xs_xb[0]*ez[0] + xs_xb[1]*ez[1] + xs_xb[2]*ez[2];

  a = (bx*bx + 2.0*ax*cx)*big->aradsqinv +
    (by*by + 2.0*ay*cy)*big->bradsqinv + 
    (bz*bz + 2.0*az*cz)*big->cradsqinv;
  b = 2.0 * (bx*cx*big->aradsqinv + by*cy*big->bradsqinv +
	     bz*cz*big->cradsqinv);
  c = cx*cx*big->aradsqinv + cy*cy*big->bradsqinv +
    cz*cz*big->cradsqinv - 1.0;

  double dt = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);

  xscoll[0] = xs[0] - dt*vs[0];
  xscoll[1] = xs[1] - dt*vs[1];
  xscoll[2] = xs[2] - dt*vs[2];

  xbcoll[0] = xb[0] - dt*vb[0];
  xbcoll[1] = xb[1] - dt*vb[1];
  xbcoll[2] = xb[2] - dt*vb[2];

  // calculate normal to ellipsoid at collision pt
  // Excoll,Eycoll,Ezcoll = orientation of ellipsoid at collision time
  // nbody = normal in body frame of ellipsoid (Excoll,Eycoll,Ezcoll)
  // norm = normal in space frame
  // only worry about normalizing final norm vector

  excoll[0] = ex[0] - dt*(omega[1]*ex[2] - omega[2]*ex[1]);
  excoll[1] = ex[1] - dt*(omega[2]*ex[0] - omega[0]*ex[2]);
  excoll[2] = ex[2] - dt*(omega[0]*ex[1] - omega[1]*ex[0]);

  eycoll[0] = ey[0] - dt*(omega[1]*ey[2] - omega[2]*ey[1]);
  eycoll[1] = ey[1] - dt*(omega[2]*ey[0] - omega[0]*ey[2]);
  eycoll[2] = ey[2] - dt*(omega[0]*ey[1] - omega[1]*ey[0]);

  ezcoll[0] = ez[0] - dt*(omega[1]*ez[2] - omega[2]*ez[1]);
  ezcoll[1] = ez[1] - dt*(omega[2]*ez[0] - omega[0]*ez[2]);
  ezcoll[2] = ez[2] - dt*(omega[0]*ez[1] - omega[1]*ez[0]);

  MathExtra::sub3(xscoll,xbcoll,delta);
  MathExtra::transpose_matvec(excoll,eycoll,ezcoll,delta,xbody);

  nbody[0] = xbody[0]*big->aradsqinv;
  nbody[1] = xbody[1]*big->bradsqinv;
  nbody[2] = xbody[2]*big->cradsqinv;

  MathExtra::matvec(excoll,eycoll,ezcoll,nbody,norm);
  MathExtra::norm3(norm);

  return dt;
}

/* ----------------------------------------------------------------------
   collision of SRD particle S with surface of ellipsoidal big particle B
   inexact because just push SRD to surface of big particle at end of step
   time of collision = end of step
   xscoll = collision pt = position of SRD at time of collision
   xbcoll = xb = position of big particle at time of collision
   norm = surface normal of collision pt at time of collision
------------------------------------------------------------------------- */

void FixSRD::collision_ellipsoid_inexact(double *xs, double *xb,
					 Big *big,
					 double *xscoll, double *xbcoll,
					 double *norm)
{
  double xs_xb[3],delta[3],xbody[3],nbody[3];

  double *ex = big->ex;
  double *ey = big->ey;
  double *ez = big->ez;

  MathExtra::sub3(xs,xb,xs_xb);
  double x = MathExtra::dot3(xs_xb,ex);
  double y = MathExtra::dot3(xs_xb,ey);
  double z = MathExtra::dot3(xs_xb,ez);

  double scale = 1.0/sqrt(x*x*big->aradsqinv + y*y*big->bradsqinv +
			  z*z*big->cradsqinv);
  x *= scale;
  y *= scale;
  z *= scale;

  xscoll[0] = x*ex[0] + y*ey[0] + z*ez[0] + xb[0];
  xscoll[1] = x*ex[1] + y*ey[1] + z*ez[1] + xb[1];
  xscoll[2] = x*ex[2] + y*ey[2] + z*ez[2] + xb[2];

  xbcoll[0] = xb[0];
  xbcoll[1] = xb[1];
  xbcoll[2] = xb[2];

  // calculate normal to ellipsoid at collision pt
  // nbody = normal in body frame of ellipsoid
  // norm = normal in space frame
  // only worry about normalizing final norm vector

  MathExtra::sub3(xscoll,xbcoll,delta);
  MathExtra::transpose_matvec(ex,ey,ez,delta,xbody);

  nbody[0] = xbody[0]*big->aradsqinv;
  nbody[1] = xbody[1]*big->bradsqinv;
  nbody[2] = xbody[2]*big->cradsqinv;

  MathExtra::matvec(ex,ey,ez,nbody,norm);
  MathExtra::norm3(norm);
}

/* ----------------------------------------------------------------------
   collision of SRD particle S with surface of line big particle B
   exact because compute time of collision
   dt = time previous to now at which collision occurs
   xscoll = collision pt = position of SRD at time of collision
   xbcoll = position of big particle at time of collision
   norm = surface normal of collision pt at time of collision
------------------------------------------------------------------------- */

double FixSRD::collision_line_exact(double *xs, double *xb,
				    double *vs, double *vb, Big *big,
				    double dt_step,
				    double *xscoll, double *xbcoll,
				    double *norm)
{
  xscoll[0] = xsc[0];
  xscoll[1] = xsc[1];
  xscoll[2] = 0.0;
  xbcoll[0] = xbc[0];
  xbcoll[1] = xbc[1];
  xbcoll[2] = 0.0;

  norm[0] = nbc[0];
  norm[1] = nbc[1];
  norm[2] = 0.0;

  return (1.0-tfraction)*dt_step;
}

/* ----------------------------------------------------------------------
   collision of SRD particle S with surface of triangle big particle B
   exact because compute time of collision
   dt = time previous to now at which collision occurs
   xscoll = collision pt = position of SRD at time of collision
   xbcoll = position of big particle at time of collision
   norm = surface normal of collision pt at time of collision
------------------------------------------------------------------------- */

double FixSRD::collision_tri_exact(double *xs, double *xb,
				   double *vs, double *vb, Big *big,
				   double dt_step,
				   double *xscoll, double *xbcoll,
				   double *norm)
{
  xscoll[0] = xsc[0];
  xscoll[1] = xsc[1];
  xscoll[2] = xsc[2];
  xbcoll[0] = xbc[0];
  xbcoll[1] = xbc[1];
  xbcoll[2] = xbc[2];

  norm[0] = nbc[0];
  norm[1] = nbc[1];
  norm[2] = nbc[2];

  return (1.0-tfraction)*dt_step;
}

/* ----------------------------------------------------------------------
   collision of SRD particle S with wall IWALL
   exact because compute time of collision
   dt = time previous to now at which collision occurs
   xscoll = collision pt = position of SRD at time of collision
   xbcoll = position of wall at time of collision
   norm = surface normal of collision pt at time of collision
------------------------------------------------------------------------- */

double FixSRD::collision_wall_exact(double *xs, int iwall, double *vs,
				     double *xscoll, double *xbcoll,
				     double *norm)
{
  int dim = wallwhich[iwall] / 2;

  double dt = (xs[dim] - xwall[iwall]) / (vs[dim] - vwall[iwall]);
  xscoll[0] = xs[0] - dt*vs[0];
  xscoll[1] = xs[1] - dt*vs[1];
  xscoll[2] = xs[2] - dt*vs[2];

  xbcoll[0] = xbcoll[1] = xbcoll[2] = 0.0;
  xbcoll[dim] = xwall[iwall] - dt*vwall[iwall];

  int side = wallwhich[iwall] % 2;
  norm[0] = norm[1] = norm[2] = 0.0;
  if (side == 0) norm[dim] = 1.0;
  else norm[dim] = -1.0;

  return dt;
}

/* ----------------------------------------------------------------------
   collision of SRD particle S with wall IWALL
   inexact because just push SRD to surface of wall at end of step
   time of collision = end of step
   xscoll = collision pt = position of SRD at time of collision
   xbcoll = position of wall at time of collision
   norm = surface normal of collision pt at time of collision
------------------------------------------------------------------------- */

void FixSRD::collision_wall_inexact(double *xs, int iwall, double *xscoll, 
				     double *xbcoll, double *norm)
{
  int dim = wallwhich[iwall] / 2;

  xscoll[0] = xs[0];
  xscoll[1] = xs[1];
  xscoll[2] = xs[2];
  xscoll[dim] = xwall[iwall];

  xbcoll[0] = xbcoll[1] = xbcoll[2] = 0.0;
  xbcoll[dim] = xwall[iwall];

  int side = wallwhich[iwall] % 2;
  norm[0] = norm[1] = norm[2] = 0.0;
  if (side == 0) norm[dim] = 1.0;
  else norm[dim] = -1.0;
}

/* ----------------------------------------------------------------------
   SLIP collision with BIG particle with omega
   vs = velocity of SRD, vb = velocity of BIG
   xb = position of BIG, omega = rotation of BIG
   xsurf = collision pt on surf of BIG
   norm = unit normal from surface of BIG at collision pt
   v of BIG particle in direction of surf normal is added to v of SRD
   includes component due to rotation of BIG
   return vsnew of SRD
------------------------------------------------------------------------- */

void FixSRD::slip(double *vs, double *vb, double *xb, Big *big,
		  double *xsurf, double *norm, double *vsnew)
{
  double r1,r2,vnmag,vs_dot_n,vsurf_dot_n;
  double tangent[3],vsurf[3];
  double *omega = big->omega;

  while (1) {
    r1 = sigma * random->gaussian();
    r2 = sigma * random->gaussian();
    vnmag = sqrt(r1*r1 + r2*r2);
    if (vnmag*vnmag <= vmaxsq) break;
  }

  vs_dot_n = vs[0]*norm[0] + vs[1]*norm[1] + vs[2]*norm[2];

  tangent[0] = vs[0] - vs_dot_n*norm[0];
  tangent[1] = vs[1] - vs_dot_n*norm[1];
  tangent[2] = vs[2] - vs_dot_n*norm[2];

  // vsurf = velocity of collision pt = translation/rotation of BIG particle
  // NOTE: for sphere could just use vsurf = vb, since w x (xsurf-xb)
  //       is orthogonal to norm and thus doesn't contribute to vsurf_dot_n

  vsurf[0] = vb[0] + omega[1]*(xsurf[2]-xb[2]) - omega[2]*(xsurf[1]-xb[1]);
  vsurf[1] = vb[1] + omega[2]*(xsurf[0]-xb[0]) - omega[0]*(xsurf[2]-xb[2]);
  vsurf[2] = vb[2] + omega[0]*(xsurf[1]-xb[1]) - omega[1]*(xsurf[0]-xb[0]);

  vsurf_dot_n = vsurf[0]*norm[0] + vsurf[1]*norm[1] + vsurf[2]*norm[2];

  vsnew[0] = (vnmag+vsurf_dot_n)*norm[0] + tangent[0];
  vsnew[1] = (vnmag+vsurf_dot_n)*norm[1] + tangent[1];
  vsnew[2] = (vnmag+vsurf_dot_n)*norm[2] + tangent[2];
}

/* ----------------------------------------------------------------------
   SLIP collision with wall IWALL
   vs = velocity of SRD
   norm = unit normal from WALL at collision pt
   v of WALL in direction of surf normal is added to v of SRD
   return vsnew of SRD
------------------------------------------------------------------------- */

void FixSRD::slip_wall(double *vs, int iwall, double *norm, double *vsnew)
{
  double vs_dot_n,scale,r1,r2,vnmag,vtmag1,vtmag2;
  double tangent1[3],tangent2[3];

  vs_dot_n = vs[0]*norm[0] + vs[1]*norm[1] + vs[2]*norm[2];

  tangent1[0] = vs[0] - vs_dot_n*norm[0];
  tangent1[1] = vs[1] - vs_dot_n*norm[1];
  tangent1[2] = vs[2] - vs_dot_n*norm[2];
  scale = 1.0/sqrt(tangent1[0]*tangent1[0] + tangent1[1]*tangent1[1] +
		   tangent1[2]*tangent1[2]);
  tangent1[0] *= scale;
  tangent1[1] *= scale;
  tangent1[2] *= scale;

  tangent2[0] = norm[1]*tangent1[2] - norm[2]*tangent1[1];
  tangent2[1] = norm[2]*tangent1[0] - norm[0]*tangent1[2];
  tangent2[2] = norm[0]*tangent1[1] - norm[1]*tangent1[0];

  while (1) {
    r1 = sigma * random->gaussian();
    r2 = sigma * random->gaussian();
    vnmag = sqrt(r1*r1 + r2*r2);
    vtmag1 = sigma * random->gaussian();
    vtmag2 = sigma * random->gaussian();
    if (vnmag*vnmag + vtmag1*vtmag1 + vtmag2*vtmag2 <= vmaxsq) break;
  }

  vsnew[0] = vnmag*norm[0] + vtmag1*tangent1[0] + vtmag2*tangent2[0];
  vsnew[1] = vnmag*norm[1] + vtmag1*tangent1[1] + vtmag2*tangent2[1];
  vsnew[2] = vnmag*norm[2] + vtmag1*tangent1[2] + vtmag2*tangent2[2];

  // add in velocity of collision pt = velocity of wall

  int dim = wallwhich[iwall] / 2;
  vsnew[dim] += vwall[iwall];
}

/* ----------------------------------------------------------------------
   NO-SLIP collision with BIG particle including WALL
   vs = velocity of SRD, vb = velocity of BIG
   xb = position of BIG, omega = rotation of BIG
   xsurf = collision pt on surf of BIG
   norm = unit normal from surface of BIG at collision pt
   v of collision pt is added to v of SRD
   includes component due to rotation of BIG
   return vsnew of SRD
------------------------------------------------------------------------- */

void FixSRD::noslip(double *vs, double *vb, double *xb, Big *big, int iwall,
		    double *xsurf, double *norm, double *vsnew)
{
  double vs_dot_n,scale,r1,r2,vnmag,vtmag1,vtmag2;
  double tangent1[3],tangent2[3];

  vs_dot_n = vs[0]*norm[0] + vs[1]*norm[1] + vs[2]*norm[2];

  tangent1[0] = vs[0] - vs_dot_n*norm[0];
  tangent1[1] = vs[1] - vs_dot_n*norm[1];
  tangent1[2] = vs[2] - vs_dot_n*norm[2];
  scale = 1.0/sqrt(tangent1[0]*tangent1[0] + tangent1[1]*tangent1[1] +
		   tangent1[2]*tangent1[2]);
  tangent1[0] *= scale;
  tangent1[1] *= scale;
  tangent1[2] *= scale;

  tangent2[0] = norm[1]*tangent1[2] - norm[2]*tangent1[1];
  tangent2[1] = norm[2]*tangent1[0] - norm[0]*tangent1[2];
  tangent2[2] = norm[0]*tangent1[1] - norm[1]*tangent1[0];

  while (1) {
    r1 = sigma * random->gaussian();
    r2 = sigma * random->gaussian();
    vnmag = sqrt(r1*r1 + r2*r2);
    vtmag1 = sigma * random->gaussian();
    vtmag2 = sigma * random->gaussian();
    if (vnmag*vnmag + vtmag1*vtmag1 + vtmag2*vtmag2 <= vmaxsq) break;
  }

  vsnew[0] = vnmag*norm[0] + vtmag1*tangent1[0] + vtmag2*tangent2[0];
  vsnew[1] = vnmag*norm[1] + vtmag1*tangent1[1] + vtmag2*tangent2[1];
  vsnew[2] = vnmag*norm[2] + vtmag1*tangent1[2] + vtmag2*tangent2[2];

  // add in velocity of collision pt
  // for WALL: velocity of wall in one dim
  // else: translation/rotation of BIG particle

  if (big->type == WALL) {
    int dim = wallwhich[iwall] / 2;
    vsnew[dim] += vwall[iwall];

  } else {
    double *omega = big->omega;
    vsnew[0] += vb[0] + omega[1]*(xsurf[2]-xb[2]) - omega[2]*(xsurf[1]-xb[1]);
    vsnew[1] += vb[1] + omega[2]*(xsurf[0]-xb[0]) - omega[0]*(xsurf[2]-xb[2]);
    vsnew[2] += vb[2] + omega[0]*(xsurf[1]-xb[1]) - omega[1]*(xsurf[0]-xb[0]);
  }
}

/* ----------------------------------------------------------------------
   impart force and torque to BIG particle
   force on BIG particle = -dp/dt of SRD particle
   torque on BIG particle = r cross (-dp/dt)
------------------------------------------------------------------------- */

void FixSRD::force_torque(double *vsold, double *vsnew,
			  double *xs, double *xb,
			  double *fb, double *tb)
{
  double dpdt[3],xs_xb[3];

  double factor = mass_srd / dt_big / force->ftm2v;
  dpdt[0] = factor * (vsnew[0] - vsold[0]);
  dpdt[1] = factor * (vsnew[1] - vsold[1]);
  dpdt[2] = factor * (vsnew[2] - vsold[2]);

  fb[0] -= dpdt[0];
  fb[1] -= dpdt[1];
  fb[2] -= dpdt[2];

  // no torque if SLIP collision and BIG is a sphere

  if (tb) {
    xs_xb[0] = xs[0]-xb[0];
    xs_xb[1] = xs[1]-xb[1];
    xs_xb[2] = xs[2]-xb[2];
    tb[0] -= xs_xb[1]*dpdt[2] - xs_xb[2]*dpdt[1];
    tb[1] -= xs_xb[2]*dpdt[0] - xs_xb[0]*dpdt[2];
    tb[2] -= xs_xb[0]*dpdt[1] - xs_xb[1]*dpdt[0];
  }
}

/* ----------------------------------------------------------------------
   impart force to WALL
   force on WALL = -dp/dt of SRD particle
------------------------------------------------------------------------- */

void FixSRD::force_wall(double *vsold, double *vsnew, int iwall)

{
  double dpdt[3];

  double factor = mass_srd / dt_big / force->ftm2v;
  dpdt[0] = factor * (vsnew[0] - vsold[0]);
  dpdt[1] = factor * (vsnew[1] - vsold[1]);
  dpdt[2] = factor * (vsnew[2] - vsold[2]);

  fwall[iwall][0] -= dpdt[0];
  fwall[iwall][1] -= dpdt[1];
  fwall[iwall][2] -= dpdt[2];
}

/* ----------------------------------------------------------------------
   update SRD particle position & velocity & search bin assignment
   check if SRD moved outside of valid region
   if so, may overlap off-processor BIG particle
------------------------------------------------------------------------- */

int FixSRD::update_srd(int i, double dt, double *xscoll, double *vsnew,
		       double *xs, double *vs)
{
  int ix,iy,iz;

  vs[0] = vsnew[0];
  vs[1] = vsnew[1];
  vs[2] = vsnew[2];

  xs[0] = xscoll[0] + dt*vsnew[0];
  xs[1] = xscoll[1] + dt*vsnew[1];
  xs[2] = xscoll[2] + dt*vsnew[2];

  if (triclinic) domain->x2lamda(xs,xs);

  if (xs[0] < srdlo[0] || xs[0] > srdhi[0] || 
      xs[1] < srdlo[1] || xs[1] > srdhi[1] || 
      xs[2] < srdlo[2] || xs[2] > srdhi[2]) {
    if (screen) {
      error->warning(FLERR,"Fix srd particle moved outside valid domain");
      fprintf(screen,"  particle %d on proc %d at timestep " BIGINT_FORMAT,
	      atom->tag[i],me,update->ntimestep);
      fprintf(screen,"  xnew %g %g %g\n",xs[0],xs[1],xs[2]);
      fprintf(screen,"  srdlo/hi x %g %g\n",srdlo[0],srdhi[0]);
      fprintf(screen,"  srdlo/hi y %g %g\n",srdlo[1],srdhi[1]);
      fprintf(screen,"  srdlo/hi z %g %g\n",srdlo[2],srdhi[2]);
    }
  }

  if (triclinic) domain->lamda2x(xs,xs);

  ix = static_cast<int> ((xs[0]-xblo2)*bininv2x);
  iy = static_cast<int> ((xs[1]-yblo2)*bininv2y);
  iz = static_cast<int> ((xs[2]-zblo2)*bininv2z);
  return iz*nbin2y*nbin2x + iy*nbin2x + ix;
}

/* ----------------------------------------------------------------------
   setup all SRD parameters with big particles
------------------------------------------------------------------------- */

void FixSRD::parameterize()
{
  // timesteps

  dt_big = update->dt;
  dt_srd = nevery * update->dt;

  // maxbigdiam,minbigdiam = max/min diameter of any big particle
  // big particle must either have radius > 0 or shape > 0 defined
  // apply radfactor at end

  AtomVecEllipsoid::Bonus *ebonus;
  if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
  AtomVecLine::Bonus *lbonus;
  if (avec_line) lbonus = avec_line->bonus;
  AtomVecTri::Bonus *tbonus;
  if (avec_tri) tbonus = avec_tri->bonus;
  double *radius = atom->radius;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int any_ellipsoids = 0;
  int any_lines = 0;
  int any_tris = 0;
  maxbigdiam = 0.0;
  minbigdiam = BIG;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & biggroupbit) {
      if (radius && radius[i] > 0.0) {
	maxbigdiam = MAX(maxbigdiam,2.0*radius[i]);
	minbigdiam = MIN(minbigdiam,2.0*radius[i]);
      } else if (ellipsoid && ellipsoid[i] >= 0) {
	any_ellipsoids = 1;
	double *shape = ebonus[ellipsoid[i]].shape;
	maxbigdiam = MAX(maxbigdiam,2.0*shape[0]);
	maxbigdiam = MAX(maxbigdiam,2.0*shape[1]);
	maxbigdiam = MAX(maxbigdiam,2.0*shape[2]);
	minbigdiam = MIN(minbigdiam,2.0*shape[0]);
	minbigdiam = MIN(minbigdiam,2.0*shape[1]);
	minbigdiam = MIN(minbigdiam,2.0*shape[2]);
      } else if (line && line[i] >= 0) {
	any_lines = 1;
	double length = lbonus[line[i]].length;
	maxbigdiam = MAX(maxbigdiam,length);
	minbigdiam = MIN(minbigdiam,length);
      } else if (tri && tri[i] >= 0) {
	any_tris = 1;
	double length1 = MathExtra::len3(tbonus[tri[i]].c1);
	double length2 = MathExtra::len3(tbonus[tri[i]].c2);
	double length3 = MathExtra::len3(tbonus[tri[i]].c3);
	double length = MAX(length1,length2);
	length = MAX(length,length3);
	maxbigdiam = MAX(maxbigdiam,length);
	minbigdiam = MIN(minbigdiam,length);
      } else 
	error->one(FLERR,"Big particle in fix srd cannot be point particle");
    }

  double tmp = maxbigdiam;
  MPI_Allreduce(&tmp,&maxbigdiam,1,MPI_DOUBLE,MPI_MAX,world);
  tmp = minbigdiam;
  MPI_Allreduce(&tmp,&minbigdiam,1,MPI_DOUBLE,MPI_MIN,world);

  maxbigdiam *= radfactor;
  minbigdiam *= radfactor;

  int itmp = any_ellipsoids;
  MPI_Allreduce(&itmp,&any_ellipsoids,1,MPI_INT,MPI_MAX,world);
  itmp = any_lines;
  MPI_Allreduce(&itmp,&any_lines,1,MPI_INT,MPI_MAX,world);
  itmp = any_tris;
  MPI_Allreduce(&itmp,&any_tris,1,MPI_INT,MPI_MAX,world);

  if (any_lines && overlap == 0)
    error->all(FLERR,"Cannot use lines with fix srd unless overlap is set");
  if (any_tris && overlap == 0)
    error->all(FLERR,"Cannot use tris with fix srd unless overlap is set");

  // big particles are only torqued if ellipsoids/lines/tris or NOSLIP

  if (any_ellipsoids == 0 && any_lines == 0 && any_tris == 0 &&
      collidestyle == SLIP) torqueflag = 0;
  else torqueflag = 1;

  // mass of SRD particles, require monodispersity

  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;

  int flag = 0;
  mass_srd = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) {
	if (mass_srd == 0.0) mass_srd = rmass[i];
	else if (rmass[i] != mass_srd) flag = 1;
      } else {
	if (mass_srd == 0.0) mass_srd = mass[type[i]];
	else if (mass[type[i]] != mass_srd) flag = 1;
      }
    }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall) 
    error->all(FLERR,"Fix srd requires SRD particles all have same mass");

  // set temperature and lamda of SRD particles from each other
  // lamda = dt_srd * sqrt(boltz * temperature_srd / mass_srd);

  if (lamdaflag == 0)
    lamda = dt_srd * sqrt(force->boltz*temperature_srd/mass_srd/force->mvv2e);
  else
    temperature_srd = force->mvv2e * 
      (lamda/dt_srd)*(lamda/dt_srd) * mass_srd/force->boltz;

  // vmax = maximum velocity of an SRD particle
  // dmax = maximum distance an SRD can move = 4*lamda = vmax * dt_srd

  sigma = lamda/dt_srd;
  dmax = 4.0*lamda;
  vmax = dmax/dt_srd;
  vmaxsq = vmax*vmax;

  // volbig = total volume of all big particles
  // LINE/TRI particles have no volume
  //   incorrect if part of rigid particles, so add fudge factor with WIDTH
  // apply radfactor to reduce final volume

  double volbig = 0.0;
  double WIDTH = 1.0;

  if (dimension == 3) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & biggroupbit) {
	if (radius && radius[i] > 0.0) {
	  double r = radfactor * radius[i];
	  volbig += 4.0/3.0*MY_PI * r*r*r;;
	} else if (ellipsoid && ellipsoid[i] >= 0) {
	  double *shape = ebonus[ellipsoid[i]].shape;
	  volbig += 4.0/3.0*MY_PI * shape[0]*shape[1]*shape[2] *
	    radfactor*radfactor*radfactor;
	} else if (tri && tri[i] >= 0) {
	  double *c1 = tbonus[tri[i]].c1;
	  double *c2 = tbonus[tri[i]].c2;
	  double *c3 = tbonus[tri[i]].c3;
	  double c2mc1[3],c3mc1[3],cross[3];
	  MathExtra::sub3(c2,c1,c2mc1);
	  MathExtra::sub3(c3,c1,c3mc1);
	  MathExtra::cross3(c2mc1,c3mc1,cross);
	  volbig += 0.5 * MathExtra::len3(cross);
	}
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & biggroupbit) {
	if (radius && radius[i] > 0.0) {
	  double r = radfactor * radius[i];
	  volbig += MY_PI * r*r;
	} else if (ellipsoid && ellipsoid[i] >= 0) {
	  double *shape = ebonus[ellipsoid[i]].shape;
	  volbig += MY_PI * shape[0]*shape[1] * radfactor*radfactor;
	} else if (line && line[i] >= 0) {
	  double length = lbonus[line[i]].length;
	  volbig += length * WIDTH;
	}
      }
  }

  tmp = volbig;
  MPI_Allreduce(&tmp,&volbig,1,MPI_DOUBLE,MPI_SUM,world);

  // particle counts

  bigint mbig = 0;
  if (bigexist) mbig = group->count(biggroup);
  bigint nsrd = group->count(igroup);

  // mass_big = total mass of all big particles

  mass_big = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & biggroupbit) {
      if (rmass) mass_big += rmass[i];
      else mass_big += mass[type[i]];
    }

  tmp = mass_big;
  MPI_Allreduce(&tmp,&mass_big,1,MPI_DOUBLE,MPI_SUM,world);

  // mass density ratio = big / SRD

  double density_big = 0.0;
  if (bigexist) density_big = mass_big / volbig;

  double volsrd,density_srd;
  if (dimension == 3) {
    volsrd = (domain->xprd * domain->yprd * domain->zprd) - volbig;
    density_srd = nsrd * mass_srd / 
      (domain->xprd*domain->yprd*domain->zprd - volbig);
  } else {
    volsrd = (domain->xprd * domain->yprd) - volbig;
    density_srd = nsrd * mass_srd /
      (domain->xprd*domain->yprd - volbig);
  }

  double mdratio = density_big/density_srd;

  // create grid for binning/rotating SRD particles from gridsrd

  setup_velocity_bins();

  // binsize3 = binsize1 in box units (not lamda units for triclinic)

  double binsize3x = binsize1x;
  double binsize3y = binsize1y;
  double binsize3z = binsize1z;
  if (triclinic) {
    binsize3x = binsize1x * domain->xprd;
    binsize3y = binsize1y * domain->yprd;
    binsize3z = binsize1z * domain->zprd;
  }

  // srd_per_grid = # of SRD particles per SRD grid cell

  double ncell;
  if (dimension == 3)
    ncell = volsrd / (binsize3x*binsize3y*binsize3z);
  else
    ncell = volsrd / (binsize3x*binsize3y);

  srd_per_cell = (double) nsrd / ncell;

  // kinematic viscosity of SRD fluid
  // output in cm^2/sec units, converted by xxt2kmu

  double viscosity;
  if (dimension == 3)
    viscosity = gridsrd*gridsrd/(18.0*dt_srd) * 
      (1.0-(1.0-exp(-srd_per_cell))/srd_per_cell) + 
      (force->boltz*temperature_srd*dt_srd/(4.0*mass_srd*force->mvv2e)) *
      ((srd_per_cell+2.0)/(srd_per_cell-1.0));
  else
    viscosity = 
      (force->boltz*temperature_srd*dt_srd/(2.0*mass_srd*force->mvv2e)) *
      (srd_per_cell/(srd_per_cell-1.0 + exp(-srd_per_cell)) - 1.0) +
      (gridsrd*gridsrd)/(12.0*dt_srd) * 
      ((srd_per_cell-1.0 + exp(-srd_per_cell))/srd_per_cell);
  viscosity *= force->xxt2kmu;

  // print SRD parameters

  if (me == 0) {
    if (screen) {
      fprintf(screen,"SRD info:\n");
      fprintf(screen,
	      "  SRD/big particles = " BIGINT_FORMAT " " BIGINT_FORMAT "\n",
	      nsrd,mbig);
      fprintf(screen,"  big particle diameter max/min = %g %g\n",
	      maxbigdiam,minbigdiam);
      fprintf(screen,"  SRD temperature & lamda = %g %g\n",
	      temperature_srd,lamda);
      fprintf(screen,"  SRD max distance & max velocity = %g %g\n",dmax,vmax);
      fprintf(screen,"  SRD grid counts: %d %d %d\n",nbin1x,nbin1y,nbin1z);
      fprintf(screen,"  SRD grid size: request, actual (xyz) = %g, %g %g %g\n",
	      gridsrd,binsize3x,binsize3y,binsize3z);
      fprintf(screen,"  SRD per actual grid cell = %g\n",srd_per_cell);
      fprintf(screen,"  SRD viscosity = %g\n",viscosity);
      fprintf(screen,"  big/SRD mass density ratio = %g\n",mdratio);
    }
    if (logfile) {
      fprintf(logfile,"SRD info:\n");
      fprintf(logfile,
	      "  SRD/big particles = " BIGINT_FORMAT " " BIGINT_FORMAT "\n",
	      nsrd,mbig);
      fprintf(logfile,"  big particle diameter max/min = %g %g\n",
	      maxbigdiam,minbigdiam);
      fprintf(logfile,"  SRD temperature & lamda = %g %g\n",
	      temperature_srd,lamda);
      fprintf(logfile,"  SRD max distance & max velocity = %g %g\n",dmax,vmax);
      fprintf(logfile,"  SRD grid counts: %d %d %d\n",nbin1x,nbin1y,nbin1z);
      fprintf(logfile,"  SRD grid size: request, actual (xyz) = %g, %g %g %g\n",
	      gridsrd,binsize3x,binsize3y,binsize3z);
      fprintf(logfile,"  SRD per actual grid cell = %g\n",srd_per_cell);
      fprintf(logfile,"  SRD viscosity = %g\n",viscosity);
      fprintf(logfile,"  big/SRD mass density ratio = %g\n",mdratio);
    }
  }

  // error if less than 1 SRD bin per processor in some dim

  if (nbin1x < comm->procgrid[0] || nbin1y < comm->procgrid[1] || 
      nbin1z < comm->procgrid[2]) 
    error->all(FLERR,"Fewer SRD bins than processors in some dimension");

  // check if SRD bins are within tolerance for shape and size

  int tolflag = 0;
  if (binsize3y/binsize3x > 1.0+cubictol ||
      binsize3x/binsize3y > 1.0+cubictol) tolflag = 1;
  if (dimension == 3) {
    if (binsize3z/binsize3x > 1.0+cubictol ||
	binsize3x/binsize3z > 1.0+cubictol) tolflag = 1;
  }

  if (tolflag) {
    if (cubicflag == CUBIC_ERROR)
      error->all(FLERR,"SRD bins for fix srd are not cubic enough");
    if (me == 0)
      error->warning(FLERR,"SRD bins for fix srd are not cubic enough");
  }

  tolflag = 0;
  if (binsize3x/gridsrd > 1.0+cubictol || gridsrd/binsize3x > 1.0+cubictol)
    tolflag = 1;
  if (binsize3y/gridsrd > 1.0+cubictol || gridsrd/binsize3y > 1.0+cubictol)
    tolflag = 1;
  if (dimension == 3) {
    if (binsize3z/gridsrd > 1.0+cubictol || gridsrd/binsize3z > 1.0+cubictol)
      tolflag = 1;
  }

  if (tolflag) {
    if (cubicflag == CUBIC_ERROR)
      error->all(FLERR,"SRD bin size for fix srd differs from user request");
    if (me == 0)
      error->warning(FLERR,
		     "SRD bin size for fix srd differs from user request");
  }

  // error if lamda < 0.6 of SRD grid size and no shifting allowed
  // turn on shifting in this case if allowed

  double maxgridsrd = MAX(binsize3x,binsize3y);
  if (dimension == 3) maxgridsrd = MAX(maxgridsrd,binsize3z);

  shiftflag = 0;
  if (lamda < 0.6*maxgridsrd && shiftuser == SHIFT_NO)
    error->all(FLERR,"Fix srd lamda must be >= 0.6 of SRD grid size");
  else if (lamda < 0.6*maxgridsrd && shiftuser == SHIFT_POSSIBLE) {
    shiftflag = 1;
    if (me == 0) 
      error->warning(FLERR,"SRD bin shifting turned on due to small lamda");
  } else if (shiftuser == SHIFT_YES) shiftflag = 1;

  // warnings

  if (bigexist && maxgridsrd > 0.25 * minbigdiam && me == 0)
    error->warning(FLERR,"Fix srd grid size > 1/4 of big particle diameter");
  if (viscosity < 0.0 && me == 0)
    error->warning(FLERR,"Fix srd viscosity < 0.0 due to low SRD density");
  if (bigexist && dt_big*vmax > minbigdiam && me == 0)
    error->warning(FLERR,"Fix srd particles may move > big particle diameter");
}

/* ----------------------------------------------------------------------
   set static parameters of each big particle, owned and ghost
   called each reneighboring
   use radfactor in distance parameters as appropriate
------------------------------------------------------------------------- */

void FixSRD::big_static()
{
  int i;
  double rad,arad,brad,crad,length,length1,length2,length3;
  double *shape,*c1,*c2,*c3;
  double c2mc1[3],c3mc1[3];

  AtomVecEllipsoid::Bonus *ebonus;
  if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
  AtomVecLine::Bonus *lbonus;
  if (avec_line) lbonus = avec_line->bonus;
  AtomVecTri::Bonus *tbonus;
  if (avec_tri) tbonus = avec_tri->bonus;
  double *radius = atom->radius;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;
  int *type = atom->type;

  double skinhalf = 0.5 * neighbor->skin;

  for (int k = 0; k < nbig; k++) {
    i = biglist[k].index;

    // sphere
    // set radius and radsq and cutoff based on radius

    if (radius && radius[i] > 0.0) {
      biglist[k].type = SPHERE;
      rad = radfactor*radius[i];
      biglist[k].radius = rad;
      biglist[k].radsq = rad*rad;
      biglist[k].cutbinsq = (rad+skinhalf) * (rad+skinhalf);

    // ellipsoid
    // set abc radsqinv and cutoff based on max radius

    } else if (ellipsoid && ellipsoid[i] >= 0) {
      shape = ebonus[ellipsoid[i]].shape;
      biglist[k].type = ELLIPSOID;
      arad = radfactor*shape[0];
      brad = radfactor*shape[1];
      crad = radfactor*shape[2];
      biglist[k].aradsqinv = 1.0/(arad*arad);
      biglist[k].bradsqinv = 1.0/(brad*brad);
      biglist[k].cradsqinv = 1.0/(crad*crad);
      rad = MAX(arad,brad);
      rad = MAX(rad,crad);
      biglist[k].cutbinsq = (rad+skinhalf) * (rad+skinhalf);

    // line
    // set length and cutoff based on 1/2 length

    } else if (line && line[i] >= 0) {
      length = lbonus[line[i]].length;
      biglist[k].type = LINE;
      biglist[k].length = length;
      rad = 0.5*length;
      biglist[k].cutbinsq = (rad+skinhalf) * (rad+skinhalf);

    // tri
    // set normbody based on c1,c2,c3
    // set cutoff based on point furthest from centroid

    } else if (tri && tri[i] >= 0) {
      biglist[k].type = TRIANGLE;
      c1 = tbonus[tri[i]].c1;
      c2 = tbonus[tri[i]].c2;
      c3 = tbonus[tri[i]].c3;
      MathExtra::sub3(c2,c1,c2mc1);
      MathExtra::sub3(c3,c1,c3mc1);
      MathExtra::cross3(c2mc1,c3mc1,biglist[k].normbody);
      length1 = MathExtra::len3(c1);
      length2 = MathExtra::len3(c2);
      length3 = MathExtra::len3(c3);
      rad = MAX(length1,length2);
      rad = MAX(rad,length3);
      biglist[k].cutbinsq = (rad+skinhalf) * (rad+skinhalf);
    }
  }
}

/* ----------------------------------------------------------------------
   set dynamic parameters of each big particle, owned and ghost
   called each timestep
------------------------------------------------------------------------- */

void FixSRD::big_dynamic()
{
  int i;
  double *shape,*quat,*inertia;
  double inertiaone[3];

  AtomVecEllipsoid::Bonus *ebonus;
  if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
  AtomVecLine::Bonus *lbonus;
  if (avec_line) lbonus = avec_line->bonus;
  AtomVecTri::Bonus *tbonus;
  if (avec_tri) tbonus = avec_tri->bonus;
  double **omega = atom->omega;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;

  for (int k = 0; k < nbig; k++) {
    i = biglist[k].index;

    // sphere
    // set omega from atom->omega directly

    if (biglist[k].type == SPHERE) {
      biglist[k].omega[0] = omega[i][0];
      biglist[k].omega[1] = omega[i][1];
      biglist[k].omega[2] = omega[i][2];

    // ellipsoid
    // set ex,ey,ez from quaternion
    // set omega from angmom & ex,ey,ez

    } else if (biglist[k].type == ELLIPSOID) {
      quat = ebonus[ellipsoid[i]].quat;
      MathExtra::q_to_exyz(quat,biglist[k].ex,biglist[k].ey,biglist[k].ez);
      shape = ebonus[ellipsoid[i]].shape;
      inertiaone[0] = EINERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
      inertiaone[1] = EINERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
      inertiaone[2] = EINERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);
      MathExtra::angmom_to_omega(angmom[i],
				 biglist[k].ex,biglist[k].ey,biglist[k].ez,
				 inertiaone,biglist[k].omega);

    // line
    // set omega from atom->omega directly

    } else if (biglist[k].type == LINE) {
      biglist[k].theta = lbonus[line[i]].theta;
      biglist[k].omega[0] = omega[i][0];
      biglist[k].omega[1] = omega[i][1];
      biglist[k].omega[2] = omega[i][2];

    // tri
    // set ex,ey,ez from quaternion
    // set omega from angmom & ex,ey,ez
    // set unit space-frame norm from body-frame norm

    } else if (biglist[k].type == TRIANGLE) {
      quat = tbonus[tri[i]].quat;
      MathExtra::q_to_exyz(quat,biglist[k].ex,biglist[k].ey,biglist[k].ez);
      inertia = tbonus[tri[i]].inertia;
      MathExtra::angmom_to_omega(angmom[i],
				 biglist[k].ex,biglist[k].ey,biglist[k].ez,
				 inertia,biglist[k].omega);
      MathExtra::matvec(biglist[k].ex,biglist[k].ey,biglist[k].ez,
		        biglist[k].normbody,biglist[k].norm);
      MathExtra::norm3(biglist[k].norm);
    }
  }
}

/* ----------------------------------------------------------------------
   set bounds for big and SRD particle movement
   called at setup() and when box size changes (but not shape)
------------------------------------------------------------------------- */

void FixSRD::setup_bounds()
{
  // triclinic scale factors
  // convert a real distance (perpendicular to box face) to a lamda distance

  double length0,length1,length2;
  if (triclinic) {
    double *h_inv = domain->h_inv;
    length0 = sqrt(h_inv[0]*h_inv[0] + h_inv[5]*h_inv[5] + h_inv[4]*h_inv[4]);
    length1 = sqrt(h_inv[1]*h_inv[1] + h_inv[3]*h_inv[3]);
    length2 = h_inv[2];
  }

  // collision object = CO = big particle or wall
  // big particles can be owned or ghost or unknown, walls are all owned
  // dist_ghost = distance from sub-domain (SD) that
  //   owned/ghost CO may move to before reneigh,
  //   used to bound search bins in setup_search_bins()
  // dist_srd = distance from SD at which SRD could collide with unknown CO,
  //   used to error check bounds of SRD movement after collisions via srdlo/hi
  // dist_srd_reneigh = distance from SD at which an SRD should trigger
  //   a reneigh, b/c next SRD move might overlap with unknown CO,
  //   used for SRD triggering of reneighboring via srdlo/hi_reneigh
  // onemove = max distance an SRD can move in one step
  // if big_particles (and possibly walls):
  //   dist_ghost = cut + 1/2 skin due to moving away before reneigh
  //   dist_srd = cut - 1/2 skin - 1/2 diam due to ghost CO moving towards
  //   dist_reneigh = dist_srd - onemove
  // if walls and no big particles:
  //   dist_ghost = 0.0, since not used
  // if no big particles or walls:
  //   dist_ghost and dist_srd = 0.0, since not used since no search bins
  //   dist_srd_reneigh = subsize - onemove = 
  //     max distance to move without being lost during comm->exchange()
  //   subsize = perp distance between sub-domain faces (orthog or triclinic)

  double cut = MAX(neighbor->cutneighmax,comm->cutghostuser);
  double onemove = dt_big*vmax;

  if (bigexist) {
    dist_ghost = cut + 0.5*neighbor->skin;
    dist_srd = cut - 0.5*neighbor->skin - 0.5*maxbigdiam;
    dist_srd_reneigh = dist_srd - onemove;
  } else if (wallexist) {
    dist_ghost = 4*onemove;
    dist_srd = 4*onemove;
    dist_srd_reneigh = 4*onemove - onemove;
  } else {
    dist_ghost = dist_srd = 0.0;
    double subsize;
    if (triclinic == 0) {
      subsize = domain->prd[0]/comm->procgrid[0];
      subsize = MIN(subsize,domain->prd[1]/comm->procgrid[1]);
      if (dimension == 3) 
	subsize = MIN(subsize,domain->prd[2]/comm->procgrid[2]);
    } else {
      subsize = 1.0/comm->procgrid[0]/length0;
      subsize = MIN(subsize,1.0/comm->procgrid[1]/length1);
      if (dimension == 3) 
	subsize = MIN(subsize,1.0/comm->procgrid[2]/length2);
    }
    dist_srd_reneigh = subsize - onemove;
  }

  // lo/hi = bbox on this proc which SRD particles must stay inside
  // lo/hi reneigh = bbox on this proc outside of which SRDs trigger a reneigh
  // for triclinic, these bbox are in lamda units

  if (triclinic == 0) {
    srdlo[0] = domain->sublo[0] - dist_srd;
    srdhi[0] = domain->subhi[0] + dist_srd;
    srdlo[1] = domain->sublo[1] - dist_srd;
    srdhi[1] = domain->subhi[1] + dist_srd;
    srdlo[2] = domain->sublo[2] - dist_srd;
    srdhi[2] = domain->subhi[2] + dist_srd;

    srdlo_reneigh[0] = domain->sublo[0] - dist_srd_reneigh;
    srdhi_reneigh[0] = domain->subhi[0] + dist_srd_reneigh;
    srdlo_reneigh[1] = domain->sublo[1] - dist_srd_reneigh;
    srdhi_reneigh[1] = domain->subhi[1] + dist_srd_reneigh;
    srdlo_reneigh[2] = domain->sublo[2] - dist_srd_reneigh;
    srdhi_reneigh[2] = domain->subhi[2] + dist_srd_reneigh;

  } else {
    srdlo[0] = domain->sublo_lamda[0] - dist_srd*length0;
    srdhi[0] = domain->subhi_lamda[0] + dist_srd*length0;
    srdlo[1] = domain->sublo_lamda[1] - dist_srd*length1;
    srdhi[1] = domain->subhi_lamda[1] + dist_srd*length1;
    srdlo[2] = domain->sublo_lamda[2] - dist_srd*length2;
    srdhi[2] = domain->subhi_lamda[2] + dist_srd*length2;

    srdlo_reneigh[0] = domain->sublo_lamda[0] - dist_srd_reneigh*length0;
    srdhi_reneigh[0] = domain->subhi_lamda[0] + dist_srd_reneigh*length0;
    srdlo_reneigh[1] = domain->sublo_lamda[1] - dist_srd_reneigh*length1;
    srdhi_reneigh[1] = domain->subhi_lamda[1] + dist_srd_reneigh*length1;
    srdlo_reneigh[2] = domain->sublo_lamda[2] - dist_srd_reneigh*length2;
    srdhi_reneigh[2] = domain->subhi_lamda[2] + dist_srd_reneigh*length2;
  }
}

/* ----------------------------------------------------------------------
   setup bins used for binning SRD particles for velocity reset
   gridsrd = desired bin size
   also setup bin shifting parameters
   also setup comm of bins that straddle processor boundaries
   called at beginning of each run
   called every reneighbor if box size changes, but not if box shape changes
------------------------------------------------------------------------- */

void FixSRD::setup_velocity_bins()
{
  // require integer # of bins across global domain

  nbin1x = static_cast<int> (domain->xprd/gridsrd + 0.5);
  nbin1y = static_cast<int> (domain->yprd/gridsrd + 0.5);
  nbin1z = static_cast<int> (domain->zprd/gridsrd + 0.5);
  if (dimension == 2) nbin1z = 1;

  if (nbin1x == 0) nbin1x = 1;
  if (nbin1y == 0) nbin1y = 1;
  if (nbin1z == 0) nbin1z = 1;

  if (triclinic == 0) {
    binsize1x = domain->xprd / nbin1x;
    binsize1y = domain->yprd / nbin1y;
    binsize1z = domain->zprd / nbin1z;
    bininv1x = 1.0/binsize1x;
    bininv1y = 1.0/binsize1y;
    bininv1z = 1.0/binsize1z;
  } else {
    binsize1x = 1.0 / nbin1x;
    binsize1y = 1.0 / nbin1y;
    binsize1z = 1.0 / nbin1z;
    bininv1x = nbin1x;
    bininv1y = nbin1y;
    bininv1z = nbin1z;
  }

  nbins1 = nbin1x*nbin1y*nbin1z;

  // setup two shifts, 0 = no shift, 1 = shift
  // initialize no shift case since static
  // shift case is dynamic, has to be initialized each time shift occurs
  // setup_velocity_shift allocates memory for vbin and sendlist/recvlist

  double *boxlo;
  if (triclinic == 0) boxlo = domain->boxlo;
  else boxlo = domain->boxlo_lamda;

  shifts[0].corner[0] = boxlo[0];
  shifts[0].corner[1] = boxlo[1];
  shifts[0].corner[2] = boxlo[2];
  setup_velocity_shift(0,0);

  shifts[1].corner[0] = boxlo[0];
  shifts[1].corner[1] = boxlo[1];
  shifts[1].corner[2] = boxlo[2];
  setup_velocity_shift(1,0);

  // allocate binhead based on max # of bins in either shift

  int max = shifts[0].nbins;
  max = MAX(max,shifts[1].nbins);

  if (max > maxbin1) {
    memory->destroy(binhead);
    maxbin1 = max;
    memory->create(binhead,max,"fix/srd:binhead");
  }

  // allocate sbuf,rbuf based on biggest bin message

  max = 0;
  for (int ishift = 0; ishift < 2; ishift++)
    for (int iswap = 0; iswap < 2*dimension; iswap++) {
      max = MAX(max,shifts[ishift].bcomm[iswap].nsend);
      max = MAX(max,shifts[ishift].bcomm[iswap].nrecv);
    }

  if (max > maxbuf) {
    memory->destroy(sbuf1);
    memory->destroy(sbuf2);
    memory->destroy(rbuf1);
    memory->destroy(rbuf2);
    maxbuf = max;
    memory->create(sbuf1,max*VBINSIZE,"fix/srd:sbuf");
    memory->create(sbuf2,max*VBINSIZE,"fix/srd:sbuf");
    memory->create(rbuf1,max*VBINSIZE,"fix/srd:rbuf");
    memory->create(rbuf2,max*VBINSIZE,"fix/srd:rbuf");
  }

  // commflag = 1 if any comm required due to bins overlapping proc boundaries

  shifts[0].commflag = 0;
  if (nbin1x % comm->procgrid[0]) shifts[0].commflag = 1;
  if (nbin1y % comm->procgrid[1]) shifts[0].commflag = 1;
  if (nbin1z % comm->procgrid[2]) shifts[0].commflag = 1;
  shifts[1].commflag = 1;
}

/* ----------------------------------------------------------------------
   setup velocity shift parameters
   set binlo[]/binhi[] and nbins,nbinx,nbiny,nbinz for this proc
   set bcomm[6] params based on bin overlaps with proc boundaries
   no comm of bins across non-periodic global boundaries
   set vbin owner flags for bins I am owner of
   ishift = 0, dynamic = 0:
     set all settings since are static
     allocate and set bcomm params and vbins
     do not comm bins that align with proc boundaries
   ishift = 1, dynamic = 0:
     set max bounds on bin counts and message sizes
     allocate and set bcomm params and vbins based on max bounds
     other settings will later change dynamically
   ishift = 1, dynamic = 1:
     set actual bin bounds and counts for specific shift
     set bcomm params and vbins (already allocated)
   called by setup_velocity_bins() and reset_velocities()
------------------------------------------------------------------------- */

void FixSRD::setup_velocity_shift(int ishift, int dynamic)
{
  int i,j,k,m,id,nsend;
  int *sendlist;
  BinComm *first,*second;
  BinAve *vbin;

  double *sublo,*subhi;
  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  int *binlo = shifts[ishift].binlo;
  int *binhi = shifts[ishift].binhi;
  double *corner = shifts[ishift].corner;
  int *procgrid = comm->procgrid;
  int *myloc = comm->myloc;
  
  binlo[0] = static_cast<int> ((sublo[0]-corner[0])*bininv1x);
  binlo[1] = static_cast<int> ((sublo[1]-corner[1])*bininv1y);
  binlo[2] = static_cast<int> ((sublo[2]-corner[2])*bininv1z);
  if (dimension == 2) shifts[ishift].binlo[2] = 0;
  
  binhi[0] = static_cast<int> ((subhi[0]-corner[0])*bininv1x);
  binhi[1] = static_cast<int> ((subhi[1]-corner[1])*bininv1y);
  binhi[2] = static_cast<int> ((subhi[2]-corner[2])*bininv1z);
  if (dimension == 2) shifts[ishift].binhi[2] = 0;

  if (ishift == 0) {
    if (myloc[0]*nbin1x % procgrid[0] == 0)
      binlo[0] = myloc[0]*nbin1x/procgrid[0];
    if (myloc[1]*nbin1y % procgrid[1] == 0)
      binlo[1] = myloc[1]*nbin1y/procgrid[1];
    if (myloc[2]*nbin1z % procgrid[2] == 0)
      binlo[2] = myloc[2]*nbin1z/procgrid[2];
    
    if ((myloc[0]+1)*nbin1x % procgrid[0] == 0)
      binhi[0] = (myloc[0]+1)*nbin1x/procgrid[0] - 1;
    if ((myloc[1]+1)*nbin1y % procgrid[1] == 0)
      binhi[1] = (myloc[1]+1)*nbin1y/procgrid[1] - 1;
    if ((myloc[2]+1)*nbin1z % procgrid[2] == 0)
      binhi[2] = (myloc[2]+1)*nbin1z/procgrid[2] - 1;
  }
  
  int nbinx = binhi[0] - binlo[0] + 1;
  int nbiny = binhi[1] - binlo[1] + 1;
  int nbinz = binhi[2] - binlo[2] + 1;

  // allow for one extra bin if shifting will occur
  
  if (ishift == 1 && dynamic == 0) {
    nbinx++;
    nbiny++;
    if (dimension == 3) nbinz++;
  }
  
  int nbins = nbinx*nbiny*nbinz;
  int nbinxy = nbinx*nbiny;
  int nbinsq = nbinx*nbiny;
  nbinsq = MAX(nbiny*nbinz,nbinsq);
  nbinsq = MAX(nbinx*nbinz,nbinsq);
  
  shifts[ishift].nbins = nbins;
  shifts[ishift].nbinx = nbinx;
  shifts[ishift].nbiny = nbiny;
  shifts[ishift].nbinz = nbinz;

  int reallocflag = 0;
  if (dynamic == 0 && nbinsq > shifts[ishift].maxbinsq) {
    shifts[ishift].maxbinsq = nbinsq;
    reallocflag = 1;
  }

  // bcomm neighbors
  // first = send in lo direction, recv from hi direction
  // second = send in hi direction, recv from lo direction
  
  if (dynamic == 0) {
    shifts[ishift].bcomm[0].sendproc = comm->procneigh[0][0];
    shifts[ishift].bcomm[0].recvproc = comm->procneigh[0][1];
    shifts[ishift].bcomm[1].sendproc = comm->procneigh[0][1];
    shifts[ishift].bcomm[1].recvproc = comm->procneigh[0][0];
    shifts[ishift].bcomm[2].sendproc = comm->procneigh[1][0];
    shifts[ishift].bcomm[2].recvproc = comm->procneigh[1][1];
    shifts[ishift].bcomm[3].sendproc = comm->procneigh[1][1];
    shifts[ishift].bcomm[3].recvproc = comm->procneigh[1][0];
    shifts[ishift].bcomm[4].sendproc = comm->procneigh[2][0];
    shifts[ishift].bcomm[4].recvproc = comm->procneigh[2][1];
    shifts[ishift].bcomm[5].sendproc = comm->procneigh[2][1];
    shifts[ishift].bcomm[5].recvproc = comm->procneigh[2][0];
  }
  
  // set nsend,nrecv and sendlist,recvlist for each swap in x,y,z
  // set nsend,nrecv = 0 if static bins align with proc boundary
  //   or to prevent dynamic bin swapping across non-periodic global boundary
  // allocate sendlist,recvlist only for dynamic = 0
  
  first = &shifts[ishift].bcomm[0];
  second = &shifts[ishift].bcomm[1];
  
  first->nsend = first->nrecv = second->nsend = second->nrecv = nbiny*nbinz;
  if (ishift == 0) {
    if (myloc[0]*nbin1x % procgrid[0] == 0)
      first->nsend = second->nrecv = 0;
    if ((myloc[0]+1)*nbin1x % procgrid[0] == 0)
      second->nsend = first->nrecv = 0;
  } else {
    if (domain->xperiodic == 0) {
      if (myloc[0] == 0) first->nsend = second->nrecv = 0;
      if (myloc[0] == procgrid[0]-1) second->nsend = first->nrecv = 0;
    }
  }

  if (reallocflag) {
    memory->destroy(first->sendlist);
    memory->destroy(first->recvlist);
    memory->destroy(second->sendlist);
    memory->destroy(second->recvlist);
    memory->create(first->sendlist,nbinsq,"fix/srd:sendlist");
    memory->create(first->recvlist,nbinsq,"fix/srd:sendlist");
    memory->create(second->sendlist,nbinsq,"fix/srd:sendlist");
    memory->create(second->recvlist,nbinsq,"fix/srd:sendlist");
  }

  m = 0;
  i = 0;
  for (j = 0; j < nbiny; j++) 
    for (k = 0; k < nbinz; k++) {
      id = k*nbinxy + j*nbinx + i;
      first->sendlist[m] = second->recvlist[m] = id;
      m++;
    }
  m = 0;
  i = nbinx-1;
  for (j = 0; j < nbiny; j++) 
    for (k = 0; k < nbinz; k++) {
      id = k*nbinxy + j*nbinx + i;
      second->sendlist[m] = first->recvlist[m] = id;
      m++;
    }
  
  first = &shifts[ishift].bcomm[2];
  second = &shifts[ishift].bcomm[3];
  
  first->nsend = first->nrecv = second->nsend = second->nrecv = nbinx*nbinz;
  if (ishift == 0) {
    if (myloc[1]*nbin1y % procgrid[1] == 0)
      first->nsend = second->nrecv = 0;
    if ((myloc[1]+1)*nbin1y % procgrid[1] == 0)
      second->nsend = first->nrecv = 0;
  } else {
    if (domain->yperiodic == 0) {
      if (myloc[1] == 0) first->nsend = second->nrecv = 0;
      if (myloc[1] == procgrid[1]-1) second->nsend = first->nrecv = 0;
    }
  }
  
  if (reallocflag) {
    memory->destroy(first->sendlist);
    memory->destroy(first->recvlist);
    memory->destroy(second->sendlist);
    memory->destroy(second->recvlist);
    memory->create(first->sendlist,nbinsq,"fix/srd:sendlist");
    memory->create(first->recvlist,nbinsq,"fix/srd:sendlist");
    memory->create(second->sendlist,nbinsq,"fix/srd:sendlist");
    memory->create(second->recvlist,nbinsq,"fix/srd:sendlist");
  }
  
  m = 0;
  j = 0;
  for (i = 0; i < nbinx; i++) 
    for (k = 0; k < nbinz; k++) {
      id = k*nbinxy + j*nbinx + i;
      first->sendlist[m] = second->recvlist[m] = id;
      m++;
    }
  m = 0;
  j = nbiny-1;
  for (i = 0; i < nbinx; i++) 
    for (k = 0; k < nbinz; k++) {
      id = k*nbinxy + j*nbinx + i;
      second->sendlist[m] = first->recvlist[m] = id;
      m++;
    }
  
  if (dimension == 3) {
    first = &shifts[ishift].bcomm[4];
    second = &shifts[ishift].bcomm[5];
    
    first->nsend = first->nrecv = second->nsend = second->nrecv = nbinx*nbiny;
    if (ishift == 0) {
      if (myloc[2]*nbin1z % procgrid[2] == 0)
	first->nsend = second->nrecv = 0;
      if ((myloc[2]+1)*nbin1z % procgrid[2] == 0)
	second->nsend = first->nrecv = 0;
    } else {
      if (domain->zperiodic == 0) {
	if (myloc[2] == 0) first->nsend = second->nrecv = 0;
	if (myloc[2] == procgrid[2]-1) second->nsend = first->nrecv = 0;
      }
    }
    
    if (reallocflag) {
      memory->destroy(first->sendlist);
      memory->destroy(first->recvlist);
      memory->destroy(second->sendlist);
      memory->destroy(second->recvlist);
      memory->create(first->sendlist,nbinx*nbiny,"fix/srd:sendlist");
      memory->create(first->recvlist,nbinx*nbiny,"fix/srd:sendlist");
      memory->create(second->sendlist,nbinx*nbiny,"fix/srd:sendlist");
      memory->create(second->recvlist,nbinx*nbiny,"fix/srd:sendlist");
    }
    
    m = 0;
    k = 0;
    for (i = 0; i < nbinx; i++) 
      for (j = 0; j < nbiny; j++) {
	id = k*nbinxy + j*nbinx + i;
	first->sendlist[m] = second->recvlist[m] = id;
	m++;
      }
    m = 0;
    k = nbinz-1;
    for (i = 0; i < nbinx; i++) 
      for (j = 0; j < nbiny; j++) {
	id = k*nbinxy + j*nbinx + i;
	second->sendlist[m] = first->recvlist[m] = id;
	m++;
      }
  }
  
  // allocate vbins, only for dynamic = 0
  
  if (dynamic == 0 && nbins > shifts[ishift].maxvbin) {
    memory->destroy(shifts[ishift].vbin);
    shifts[ishift].maxvbin = nbins;
    shifts[ishift].vbin = (BinAve *) 
      memory->smalloc(nbins*sizeof(BinAve),"fix/srd:vbin");
  }
  
  // for vbins I own, set owner = 1
  // if bin never sent to anyone, I own it
  // if bin sent to lower numbered proc, I do not own it
  // if bin sent to self, I do not own it on even swap (avoids double counting)
  
  vbin = shifts[ishift].vbin;
  for (i = 0; i < nbins; i++) vbin[i].owner = 1;
  for (int iswap = 0; iswap < 2*dimension; iswap++) {
    if (shifts[ishift].bcomm[iswap].sendproc > me) continue;
    if (shifts[ishift].bcomm[iswap].sendproc == me && iswap % 2 == 0) continue;
    nsend = shifts[ishift].bcomm[iswap].nsend;
    sendlist = shifts[ishift].bcomm[iswap].sendlist;
    for (m = 0; m < nsend; m++) vbin[sendlist[m]].owner = 0;
  }

  // if tstat and deformflag:
  // set xctr for all nbins in lamda units so can later compute vstream of bin

  if (tstat && deformflag) {
    m = 0;
    for (k = 0; k < nbinz; k++)
      for (j = 0; j < nbiny; j++)
	for (i = 0; i < nbinx; i++) {
	  vbin[m].xctr[0] = corner[0] + (i+binlo[0]+0.5)/nbin1x;
	  vbin[m].xctr[1] = corner[1] + (j+binlo[1]+0.5)/nbin1y;
	  vbin[m].xctr[2] = corner[2] + (k+binlo[2]+0.5)/nbin1z;
	  m++;
	}
  }
}

/* ----------------------------------------------------------------------
   setup bins used for big and SRD particle searching
   gridsearch = desired bin size
   bins are orthogonal whether simulation box is orthogonal or triclinic
   for orthogonal boxes, called at each setup since cutghost may change
   for triclinic boxes, called at every reneigh, since tilt may change
   sets the following:
     nbin2 xyz = # of bins in each dim
     binsize2 and bininv2 xyz = size of bins in each dim
     xyz blo2 = origin of bins
------------------------------------------------------------------------- */

void FixSRD::setup_search_bins()
{
  // subboxlo/hi = real space bbox which 
  //   owned/ghost big particles or walls can be in
  // start with bounding box for my sub-domain, add dist_ghost
  // for triclinic, need to:
  //   a) convert dist_ghost to lamda space via length012
  //   b) lo/hi = sub-domain big particle bbox in lamda space
  //   c) convert lo/hi to real space bounding box via domain->bbox()
  //   similar to neighbor::setup_bins() and comm::cutghost[] calculation

  double subboxlo[3],subboxhi[3];

  if (triclinic == 0) {
    subboxlo[0] = domain->sublo[0] - dist_ghost;
    subboxlo[1] = domain->sublo[1] - dist_ghost;
    subboxlo[2] = domain->sublo[2] - dist_ghost;
    subboxhi[0] = domain->subhi[0] + dist_ghost;
    subboxhi[1] = domain->subhi[1] + dist_ghost;
    subboxhi[2] = domain->subhi[2] + dist_ghost;
  } else {
    double *h_inv = domain->h_inv;
    double length0,length1,length2;
    length0 = sqrt(h_inv[0]*h_inv[0] + h_inv[5]*h_inv[5] + h_inv[4]*h_inv[4]);
    length1 = sqrt(h_inv[1]*h_inv[1] + h_inv[3]*h_inv[3]);
    length2 = h_inv[2];
    double lo[3],hi[3];
    lo[0] = domain->sublo_lamda[0] - dist_ghost*length0;
    lo[1] = domain->sublo_lamda[1] - dist_ghost*length1;
    lo[2] = domain->sublo_lamda[2] - dist_ghost*length2;
    hi[0] = domain->subhi_lamda[0] + dist_ghost*length0;
    hi[1] = domain->subhi_lamda[1] + dist_ghost*length1;
    hi[2] = domain->subhi_lamda[2] + dist_ghost*length2;
    domain->bbox(lo,hi,subboxlo,subboxhi);
  }

  // require integer # of bins for that volume

  nbin2x = static_cast<int> ((subboxhi[0] - subboxlo[0]) / gridsearch);
  nbin2y = static_cast<int> ((subboxhi[1] - subboxlo[1]) / gridsearch);
  nbin2z = static_cast<int> ((subboxhi[2] - subboxlo[2]) / gridsearch);
  if (dimension == 2) nbin2z = 1;

  if (nbin2x == 0) nbin2x = 1;
  if (nbin2y == 0) nbin2y = 1;
  if (nbin2z == 0) nbin2z = 1;

  binsize2x = (subboxhi[0] - subboxlo[0]) / nbin2x;
  binsize2y = (subboxhi[1] - subboxlo[1]) / nbin2y;
  binsize2z = (subboxhi[2] - subboxlo[2]) / nbin2z;
  bininv2x = 1.0/binsize2x;
  bininv2y = 1.0/binsize2y;
  bininv2z = 1.0/binsize2z;

  // add bins on either end due to extent of big particles
  // radmax = max distance from central bin that biggest particle overlaps
  // includes skin movement
  // nx,ny,nz = max # of bins to search away from central bin

  double radmax = 0.5*maxbigdiam + 0.5*neighbor->skin;

  int nx = static_cast<int> (radmax/binsize2x) + 1;
  int ny = static_cast<int> (radmax/binsize2y) + 1;
  int nz = static_cast<int> (radmax/binsize2z) + 1;
  if (dimension == 2) nz = 0;

  nbin2x += 2*nx;
  nbin2y += 2*ny;
  nbin2z += 2*nz;

  xblo2 = subboxlo[0] - nx*binsize2x;
  yblo2 = subboxlo[1] - ny*binsize2y;
  zblo2 = subboxlo[2] - nz*binsize2z;
  if (dimension == 2) zblo2 = domain->boxlo[2];

  // allocate bins
  // first deallocate previous bins if necessary

  nbins2 = nbin2x*nbin2y*nbin2z;
  if (nbins2 > maxbin2) {
    memory->destroy(nbinbig);
    memory->destroy(binbig);
    maxbin2 = nbins2;
    memory->create(nbinbig,nbins2,"fix/srd:nbinbig");
    memory->create(binbig,nbins2,ATOMPERBIN,"fix/srd:binbig");
  }
}

/* ----------------------------------------------------------------------
   compute stencil of bin offsets for a big particle overlapping search bins
------------------------------------------------------------------------- */

void FixSRD::setup_search_stencil()
{
  // radmax = max distance from central bin that any big particle overlaps
  // includes skin movement
  // nx,ny,nz = max # of bins to search away from central bin

  double radmax = 0.5*maxbigdiam + 0.5*neighbor->skin;
  double radsq = radmax*radmax;

  int nx = static_cast<int> (radmax/binsize2x) + 1;
  int ny = static_cast<int> (radmax/binsize2y) + 1;
  int nz = static_cast<int> (radmax/binsize2z) + 1;
  if (dimension == 2) nz = 0;

  // allocate stencil array
  // deallocate previous stencil if necessary

  int max = (2*nx+1) * (2*ny+1) * (2*nz+1);
  if (max > maxstencil) {
    memory->destroy(stencil);
    maxstencil = max;
    memory->create(stencil,max,4,"fix/srd:stencil");
  }

  // loop over all bins
  // add bin to stencil:
  // if distance of closest corners of bin and central bin is within cutoff

  nstencil = 0;
  for (int k = -nz; k <= nz; k++)
    for (int j = -ny; j <= ny; j++)
      for (int i = -nx; i <= nx; i++)
	if (bin_bin_distance(i,j,k) < radsq) {
	  stencil[nstencil][0] = i;
	  stencil[nstencil][1] = j;
	  stencil[nstencil][2] = k;
	  stencil[nstencil][3] = k*nbin2y*nbin2x + j*nbin2x + i;
	  nstencil++;
	}
}

/* ----------------------------------------------------------------------
   compute closest squared distance between point x and bin ibin
------------------------------------------------------------------------- */

double FixSRD::point_bin_distance(double *x, int i, int j, int k)
{
  double delx,dely,delz;

  double xlo = xblo2 + i*binsize2x;
  double xhi = xlo + binsize2x;
  double ylo = yblo2 + j*binsize2y;
  double yhi = ylo + binsize2y;
  double zlo = zblo2 + k*binsize2z;
  double zhi = zlo + binsize2z;

  if (x[0] < xlo) delx = xlo - x[0];
  else if (x[0] > xhi) delx = x[0] - xhi;
  else delx = 0.0;

  if (x[1] < ylo) dely = ylo - x[1];
  else if (x[1] > yhi) dely = x[1] - yhi;
  else dely = 0.0;

  if (x[2] < zlo) delz = zlo - x[2];
  else if (x[2] > zhi) delz = x[2] - zhi;
  else delz = 0.0;

  return (delx*delx + dely*dely + delz*delz);
}

/* ----------------------------------------------------------------------
   compute closest squared distance between 2 bins
   central bin (0,0,0) and bin (i,j,k)
------------------------------------------------------------------------- */

double FixSRD::bin_bin_distance(int i, int j, int k)
{
  double delx,dely,delz;

  if (i > 0) delx = (i-1)*binsize2x;
  else if (i == 0) delx = 0.0;
  else delx = (i+1)*binsize2x;

  if (j > 0) dely = (j-1)*binsize2y;
  else if (j == 0) dely = 0.0;
  else dely = (j+1)*binsize2y;

  if (k > 0) delz = (k-1)*binsize2z;
  else if (k == 0) delz = 0.0;
  else delz = (k+1)*binsize2z;

  return (delx*delx + dely*dely + delz*delz);
}

/* ---------------------------------------------------------------------- */

int FixSRD::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (torqueflag) {
    for (i = first; i < last; i++) {
      buf[m++] = flocal[i][0];
      buf[m++] = flocal[i][1];
      buf[m++] = flocal[i][2];
      buf[m++] = tlocal[i][0];
      buf[m++] = tlocal[i][1];
      buf[m++] = tlocal[i][2];
    }

  } else {
    for (i = first; i < last; i++) {
      buf[m++] = flocal[i][0];
      buf[m++] = flocal[i][1];
      buf[m++] = flocal[i][2];
    }
  }

  return comm_reverse;
}

/* ---------------------------------------------------------------------- */

void FixSRD::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;

  if (torqueflag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      flocal[j][0] += buf[m++];
      flocal[j][1] += buf[m++];
      flocal[j][2] += buf[m++];
      tlocal[j][0] += buf[m++];
      tlocal[j][1] += buf[m++];
      tlocal[j][2] += buf[m++];
    }

  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      flocal[j][0] += buf[m++];
      flocal[j][1] += buf[m++];
      flocal[j][2] += buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   SRD stats
------------------------------------------------------------------------- */

double FixSRD::compute_vector(int n)
{
  // only sum across procs one time

  if (stats_flag == 0) {
    stats[0] = ncheck;
    stats[1] = ncollide;
    stats[2] = nbounce;
    stats[3] = ninside;
    stats[4] = nrescale;
    stats[5] = nbins2;
    stats[6] = nbins1;
    stats[7] = srd_bin_count;
    stats[8] = srd_bin_temp;
    stats[9] = bouncemaxnum;
    stats[10] = bouncemax;
    stats[11] = reneighcount;

    MPI_Allreduce(stats,stats_all,10,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&stats[10],&stats_all[10],1,MPI_DOUBLE,MPI_MAX,world);
    if (stats_all[7] != 0.0) stats_all[8] /= stats_all[7];
    stats_all[6] /= nprocs;

    stats_flag = 1;
  }

  return stats_all[n];
}

/* ---------------------------------------------------------------------- */

void FixSRD::velocity_stats(int groupnum)
{
  int bitmask = group->bitmask[groupnum];

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double vone;
  double vave = 0.0;
  double vmax = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & bitmask) {
      vone = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      vave += vone;
      if (vone > vmax) vmax = vone;
    }

  double all;
  MPI_Allreduce(&vave,&all,1,MPI_DOUBLE,MPI_SUM,world);
  double count = group->count(groupnum);
  if (count != 0.0) vave = all/count;
  else vave = 0.0;

  MPI_Allreduce(&vmax,&all,1,MPI_DOUBLE,MPI_MAX,world);
  vmax = all;

  if (me == 0) {
    if (screen)
      fprintf(screen,"  ave/max %s velocity = %g %g\n",
	      group->names[groupnum],vave,vmax);
    if (logfile)
      fprintf(logfile,"  ave/max %s velocity = %g %g\n",
	      group->names[groupnum],vave,vmax);
  }
}

/* ---------------------------------------------------------------------- */

double FixSRD::newton_raphson(double t1, double t2)
{
  double f1,f2,df,tlo,thi;
  lineside(t1,f1,df);
  if (f1 < 0.0) {
    tlo = t1;
    thi = t2;
  } else {
    thi = t1;
    tlo = t2;
  }

  double f;
  double t = 0.5*(t1+t2);
  double dtold = fabs(t2-t1);
  double dt = dtold;
  lineside(t,f,df);

  double temp;
  for (int i = 0; i < MAXITER; i++) {
    if ((((t-thi)*df - f)*((t-tlo)*df - f) > 0.0) ||
	(fabs(2.0*f) > fabs(dtold*df))) {
      dtold = dt;
      dt = 0.5 * (thi-tlo);
      t = tlo + dt;
      if (tlo == t) return t;
    } else {
      dtold = dt;
      dt = f / df;
      temp = t;
      t -= dt;
      if (temp == t) return t;
    }
    if (fabs(dt) < TOLERANCE) return t;
    lineside(t,f,df);
    if (f < 0.0) tlo = t;
    else thi = t;
  }

  return t;
}

/* ---------------------------------------------------------------------- */

void FixSRD::lineside(double t, double &f, double &df)
{
  double p[2],c[2];

  p[0] = xs0[0] + (xs1[0]-xs0[0])*t;
  p[1] = xs0[1] + (xs1[1]-xs0[1])*t;
  c[0] = xb0[0] + (xb1[0]-xb0[0])*t;
  c[1] = xb0[1] + (xb1[1]-xb0[1])*t;
  double dtheta = theta1 - theta0;
  double theta = theta0 + dtheta*t;
  double cosT = cos(theta);
  double sinT = sin(theta);

  f = (p[1]-c[1]) * cosT - (p[0]-c[0]) * sinT;
  df = ((xs1[1]-xs0[1]) - (xb1[1]-xb0[1]))*cosT - (p[1]-c[1])*sinT*dtheta - 
    ((xs1[0]-xs0[0]) - (xb1[0]-xb0[0]))*sinT - (p[0]-c[0])*cosT*dtheta;
}

/* ---------------------------------------------------------------------- */

void FixSRD::triside(double t, double &f, double &df)
{
  double p[2],c[2];

  p[0] = xs0[0] + (xs1[0]-xs0[0])*t;
  p[1] = xs0[1] + (xs1[1]-xs0[1])*t;
  c[0] = xb0[0] + (xb1[0]-xb0[0])*t;
  c[1] = xb0[1] + (xb1[1]-xb0[1])*t;
  double dtheta = theta1 - theta0;
  double theta = theta0 + dtheta*t;
  double cosT = cos(theta);
  double sinT = sin(theta);

  f = (p[1]-c[1]) * cosT - (p[0]-c[0]) * sinT;
  df = ((xs1[1]-xs0[1]) - (xb1[1]-xb0[1]))*cosT - (p[1]-c[1])*sinT*dtheta - 
    ((xs1[0]-xs0[0]) - (xb1[0]-xb0[0]))*sinT - (p[0]-c[0])*cosT*dtheta;
}

/* ---------------------------------------------------------------------- */

double FixSRD::memory_usage()
{
  double bytes = 0.0;
  bytes += (shifts[0].nbins + shifts[1].nbins) * sizeof(BinAve);
  bytes += nmax * sizeof(int);
  if (bigexist) {
    bytes += nbins2 * sizeof(int);
    bytes += nbins2*ATOMPERBIN * sizeof(int);
  }
  bytes += nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   useful debugging functions
------------------------------------------------------------------------- */

double FixSRD::distance(int i, int j)
{
  double dx = atom->x[i][0] - atom->x[j][0];
  double dy = atom->x[i][1] - atom->x[j][1];
  double dz = atom->x[i][2] - atom->x[j][2];
  return sqrt(dx*dx + dy*dy + dz*dz);
}

/* ---------------------------------------------------------------------- */

void FixSRD::print_collision(int i, int j, int ibounce,
			      double t_remain, double dt,
			      double *xscoll, double *xbcoll, double *norm,
			      int type)
{
  double xsstart[3],xbstart[3];
  double **x = atom->x;
  double **v = atom->v;

  if (type != WALL) {
    printf("COLLISION between SRD %d and BIG %d\n",atom->tag[i],atom->tag[j]);
    printf("  bounce # = %d\n",ibounce+1);
    printf("  local indices: %d %d\n",i,j);
    printf("  timestep = %g\n",dt);
    printf("  time remaining post-collision = %g\n",t_remain);
    
    xsstart[0] = x[i][0] - dt*v[i][0];
    xsstart[1] = x[i][1] - dt*v[i][1];
    xsstart[2] = x[i][2] - dt*v[i][2];
    xbstart[0] = x[j][0] - dt*v[j][0];
    xbstart[1] = x[j][1] - dt*v[j][1];
    xbstart[2] = x[j][2] - dt*v[j][2];

    printf("  SRD start position = %g %g %g\n",
	   xsstart[0],xsstart[1],xsstart[2]);
    printf("  BIG start position = %g %g %g\n",
	   xbstart[0],xbstart[1],xbstart[2]);
    printf("  SRD coll  position = %g %g %g\n",
	   xscoll[0],xscoll[1],xscoll[2]);
    printf("  BIG coll  position = %g %g %g\n",
	   xbcoll[0],xbcoll[1],xbcoll[2]);
    printf("  SRD end   position = %g %g %g\n",x[i][0],x[i][1],x[i][2]);
    printf("  BIG end   position = %g %g %g\n",x[j][0],x[j][1],x[j][2]);

    printf("  SRD vel = %g %g %g\n",v[i][0],v[i][1],v[i][2]);
    printf("  BIG vel = %g %g %g\n",v[j][0],v[j][1],v[j][2]);
    printf("  surf norm = %g %g %g\n",norm[0],norm[1],norm[2]);
    
    double rstart = sqrt((xsstart[0]-xbstart[0])*(xsstart[0]-xbstart[0]) +
			 (xsstart[1]-xbstart[1])*(xsstart[1]-xbstart[1]) +
			 (xsstart[2]-xbstart[2])*(xsstart[2]-xbstart[2]));
    double rcoll = sqrt((xscoll[0]-xbcoll[0])*(xscoll[0]-xbcoll[0]) +
			(xscoll[1]-xbcoll[1])*(xscoll[1]-xbcoll[1]) +
			(xscoll[2]-xbcoll[2])*(xscoll[2]-xbcoll[2]));
    double rend = sqrt((x[i][0]-x[j][0])*(x[i][0]-x[j][0]) +
		       (x[i][1]-x[j][1])*(x[i][1]-x[j][1]) +
		       (x[i][2]-x[j][2])*(x[i][2]-x[j][2]));
    
    printf("  separation at start = %g\n",rstart);
    printf("  separation at coll  = %g\n",rcoll);
    printf("  separation at end   = %g\n",rend);

  } else {
    int dim = wallwhich[j] / 2;

    printf("COLLISION between SRD %d and WALL %d\n",atom->tag[i],j);
    printf("  bounce # = %d\n",ibounce+1);
    printf("  local indices: %d %d\n",i,j);
    printf("  timestep = %g\n",dt);
    printf("  time remaining post-collision = %g\n",t_remain);
    
    xsstart[0] = x[i][0] - dt*v[i][0];
    xsstart[1] = x[i][1] - dt*v[i][1];
    xsstart[2] = x[i][2] - dt*v[i][2];
    xbstart[0] = xbstart[1] = xbstart[2] = 0.0;
    xbstart[dim] = xwall[j] - dt*vwall[j];

    printf("  SRD start position = %g %g %g\n",
	   xsstart[0],xsstart[1],xsstart[2]);
    printf("  WALL start position = %g\n",xbstart[dim]);
    printf("  SRD coll  position = %g %g %g\n",
	   xscoll[0],xscoll[1],xscoll[2]);
    printf("  WALL coll position = %g\n",xbcoll[dim]);
    printf("  SRD end   position = %g %g %g\n",x[i][0],x[i][1],x[i][2]);
    printf("  WALL end  position = %g\n",xwall[j]);

    printf("  SRD vel = %g %g %g\n",v[i][0],v[i][1],v[i][2]);
    printf("  WALL vel = %g\n",vwall[j]);
    printf("  surf norm = %g %g %g\n",norm[0],norm[1],norm[2]);
    
    double rstart = xsstart[dim]-xbstart[dim];
    double rcoll = xscoll[dim]-xbcoll[dim];
    double rend = x[dim][0]-xwall[j];
    
    printf("  separation at start = %g\n",rstart);
    printf("  separation at coll  = %g\n",rcoll);
    printf("  separation at end   = %g\n",rend);
  }
}
