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
#include "fix_deposit.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDeposit::FixDeposit(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix deposit command");

  restart_global = 1;
  time_depend = 1;

  // required args

  ninsert = force->inumeric(FLERR,arg[3]);
  ntype = force->inumeric(FLERR,arg[4]);
  nfreq = force->inumeric(FLERR,arg[5]);
  seed = force->inumeric(FLERR,arg[6]);

  if (seed <= 0) error->all(FLERR,"Illegal fix deposit command");

  // set defaults

  iregion = -1;
  idregion = NULL;
  idnext = 0;
  globalflag = localflag = 0;
  lo = hi = deltasq = 0.0;
  nearsq = 0.0;
  maxattempt = 10;
  rateflag = 0;
  vxlo = vxhi = vylo = vyhi = vzlo = vzhi = 0.0;
  scaleflag = 1;
  targetflag = 0;

  // read options from end of input line

  options(narg-7,&arg[7]);

  // error checks on region and its extent being inside simulation box

  if (iregion == -1) error->all(FLERR,"Must specify a region in fix deposit");
  if (domain->regions[iregion]->bboxflag == 0)
    error->all(FLERR,"Fix deposit region does not support a bounding box");
  if (domain->regions[iregion]->dynamic_check())
    error->all(FLERR,"Fix deposit region cannot be dynamic");

  xlo = domain->regions[iregion]->extent_xlo;
  xhi = domain->regions[iregion]->extent_xhi;
  ylo = domain->regions[iregion]->extent_ylo;
  yhi = domain->regions[iregion]->extent_yhi;
  zlo = domain->regions[iregion]->extent_zlo;
  zhi = domain->regions[iregion]->extent_zhi;

  if (domain->triclinic == 0) {
    if (xlo < domain->boxlo[0] || xhi > domain->boxhi[0] ||
        ylo < domain->boxlo[1] || yhi > domain->boxhi[1] ||
        zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all(FLERR,"Deposition region extends outside simulation box");
  } else {
    if (xlo < domain->boxlo_bound[0] || xhi > domain->boxhi_bound[0] ||
        ylo < domain->boxlo_bound[1] || yhi > domain->boxhi_bound[1] ||
        zlo < domain->boxlo_bound[2] || zhi > domain->boxhi_bound[2])
      error->all(FLERR,"Deposition region extends outside simulation box");
  }

  // setup scaling

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // apply scaling to all input parameters with dist/vel units

  if (domain->dimension == 2) {
    lo *= yscale;
    hi *= yscale;
    rate *= yscale;
  } else {
    lo *= zscale;
    hi *= zscale;
    rate *= zscale;
  }
  deltasq *= xscale*xscale;
  nearsq *= xscale*xscale;
  vxlo *= xscale;
  vxhi *= xscale;
  vylo *= yscale;
  vyhi *= yscale;
  vzlo *= zscale;
  vzhi *= zscale;
  tx *= xscale;
  ty *= yscale;
  tz *= zscale;

  // maxtag_all = current max tag for all atoms

  if (idnext) {
    int *tag = atom->tag;
    int nlocal = atom->nlocal;

    int maxtag = 0;
    for (int i = 0; i < nlocal; i++) maxtag = MAX(maxtag,tag[i]);
    MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_INT,MPI_MAX,world);
  }

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  nfirst = next_reneighbor;
  ninserted = 0;
}

/* ---------------------------------------------------------------------- */

FixDeposit::~FixDeposit()
{
  delete random;
  delete [] idregion;
}

/* ---------------------------------------------------------------------- */

int FixDeposit::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDeposit::init()
{
  // set index and check validity of region

  iregion = domain->find_region(idregion);
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix deposit does not exist");
}

/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixDeposit::pre_exchange()
{
  int i,j;
  int flag,flagall;
  double coord[3],lamda[3],delx,dely,delz,rsq;
  double *newcoord;

  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  // compute current offset = bottom of insertion volume

  double offset = 0.0;
  if (rateflag) offset = (update->ntimestep - nfirst) * update->dt * rate;

  double *sublo,*subhi;
  if (domain->triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // attempt an insertion until successful

  int nfix = modify->nfix;
  Fix **fix = modify->fix;

  int success = 0;
  int attempt = 0;
  while (attempt < maxattempt) {
    attempt++;

    // choose random position for new atom within region

    coord[0] = xlo + random->uniform() * (xhi-xlo);
    coord[1] = ylo + random->uniform() * (yhi-ylo);
    coord[2] = zlo + random->uniform() * (zhi-zlo);
    while (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) == 0) {
      coord[0] = xlo + random->uniform() * (xhi-xlo);
      coord[1] = ylo + random->uniform() * (yhi-ylo);
      coord[2] = zlo + random->uniform() * (zhi-zlo);
    }

    // adjust vertical coord by offset

    if (domain->dimension == 2) coord[1] += offset;
    else coord[2] += offset;

    // if global, reset vertical coord to be lo-hi above highest atom
    // if local, reset vertical coord to be lo-hi above highest "nearby" atom
    // local computation computes lateral distance between 2 particles w/ PBC

    if (globalflag || localflag) {
      int dim;
      double max,maxall,delx,dely,delz,rsq;

      if (domain->dimension == 2) {
        dim = 1;
        max = domain->boxlo[1];
      } else {
        dim = 2;
        max = domain->boxlo[2];
      }

      double **x = atom->x;
      int nlocal = atom->nlocal;
      for (i = 0; i < nlocal; i++) {
        if (localflag) {
          delx = coord[0] - x[i][0];
          dely = coord[1] - x[i][1];
          delz = 0.0;
          domain->minimum_image(delx,dely,delz);
          if (domain->dimension == 2) rsq = delx*delx;
          else rsq = delx*delx + dely*dely;
          if (rsq > deltasq) continue;
        }
        if (x[i][dim] > max) max = x[i][dim];
      }

      MPI_Allreduce(&max,&maxall,1,MPI_DOUBLE,MPI_MAX,world);
      if (domain->dimension == 2)
        coord[1] = maxall + lo + random->uniform()*(hi-lo);
      else
        coord[2] = maxall + lo + random->uniform()*(hi-lo);
    }

    // now have final coord
    // if distance to any atom is less than near, try again

    double **x = atom->x;
    int nlocal = atom->nlocal;

    flag = 0;
    for (i = 0; i < nlocal; i++) {
      delx = coord[0] - x[i][0];
      dely = coord[1] - x[i][1];
      delz = coord[2] - x[i][2];
      domain->minimum_image(delx,dely,delz);
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < nearsq) flag = 1;
    }
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
    if (flagall) continue;

    // insertion will proceed
    // choose random velocity for new atom

    double vxtmp = vxlo + random->uniform() * (vxhi-vxlo);
    double vytmp = vylo + random->uniform() * (vyhi-vylo);
    double vztmp = vzlo + random->uniform() * (vzhi-vzlo);

    // if we have a sputter target change velocity vector accordingly
    if (targetflag) {
      double vel = sqrt(vxtmp*vxtmp + vytmp*vytmp + vztmp*vztmp);
      delx = tx - coord[0];
      dely = ty - coord[1];
      delz = tz - coord[2];
      double rsq = delx*delx + dely*dely + delz*delz;
      if (rsq > 0.0) {
        double rinv = sqrt(1.0/rsq);
        vxtmp = delx*rinv*vel;
        vytmp = dely*rinv*vel;
        vztmp = delz*rinv*vel;
      }
    }

    // check if new atom is in my sub-box or above it if I'm highest proc
    // if so, add to my list via create_atom()
    // initialize info about the atoms
    // set group mask to "all" plus fix group

    if (domain->triclinic) {
      domain->x2lamda(coord,lamda);
      newcoord = lamda;
    } else newcoord = coord;

    flag = 0;
    if (newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
        newcoord[1] >= sublo[1] && newcoord[1] < subhi[1] &&
        newcoord[2] >= sublo[2] && newcoord[2] < subhi[2]) flag = 1;
    else if (domain->dimension == 3 && newcoord[2] >= domain->boxhi[2] &&
             comm->myloc[2] == comm->procgrid[2]-1 &&
             newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
             newcoord[1] >= sublo[1] && newcoord[1] < subhi[1]) flag = 1;
    else if (domain->dimension == 2 && newcoord[1] >= domain->boxhi[1] &&
             comm->myloc[1] == comm->procgrid[1]-1 &&
             newcoord[0] >= sublo[0] && newcoord[0] < subhi[0]) flag = 1;

    if (flag) {
      atom->avec->create_atom(ntype,coord);
      int m = atom->nlocal - 1;
      atom->type[m] = ntype;
      atom->mask[m] = 1 | groupbit;
      atom->v[m][0] = vxtmp;
      atom->v[m][1] = vytmp;
      atom->v[m][2] = vztmp;
      for (j = 0; j < nfix; j++)
        if (fix[j]->create_attribute) fix[j]->set_arrays(m);
    }
    MPI_Allreduce(&flag,&success,1,MPI_INT,MPI_MAX,world);
    break;
  }

  // warn if not successful b/c too many attempts or no proc owned particle

  if (!success && comm->me == 0)
    error->warning(FLERR,"Particle deposition was unsuccessful",0);

  // reset global natoms
  // if idnext, set new atom ID to incremented maxtag_all
  // else set new atom ID to value beyond all current atoms
  // if global map exists, reset it now instead of waiting for comm
  // since adding an atom messes up ghosts

  if (success) {
    atom->natoms += 1;
    if (atom->tag_enable) {
      if (idnext) {
        maxtag_all++;
        if (atom->nlocal && atom->tag[atom->nlocal-1] == 0) 
          atom->tag[atom->nlocal-1] = maxtag_all;
      } else atom->tag_extend();
      if (atom->map_style) {
        atom->nghost = 0;
        atom->map_init();
        atom->map_set();
      }
    }
  }

  // next timestep to insert
  // next_reneighbor = 0 if done

  if (success) ninserted++;
  if (ninserted < ninsert) next_reneighbor += nfreq;
  else next_reneighbor = 0;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixDeposit::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix indent command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix deposit command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix deposit does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"id") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix deposit command");
      if (strcmp(arg[iarg+1],"max") == 0) idnext = 0;
      else if (strcmp(arg[iarg+1],"next") == 0) idnext = 1;
      else error->all(FLERR,"Illegal fix deposit command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"global") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix deposit command");
      globalflag = 1;
      localflag = 0;
      lo = force->numeric(FLERR,arg[iarg+1]);
      hi = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"local") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix deposit command");
      localflag = 1;
      globalflag = 0;
      lo = force->numeric(FLERR,arg[iarg+1]);
      hi = force->numeric(FLERR,arg[iarg+2]);
      deltasq = force->numeric(FLERR,arg[iarg+3])*force->numeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"near") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix deposit command");
      nearsq = force->numeric(FLERR,arg[iarg+1])*force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"attempt") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix deposit command");
      maxattempt = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"rate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix deposit command");
      rateflag = 1;
      rate = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vx") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix deposit command");
      vxlo = force->numeric(FLERR,arg[iarg+1]);
      vxhi = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix deposit command");
      vylo = force->numeric(FLERR,arg[iarg+1]);
      vyhi = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix deposit command");
      vzlo = force->numeric(FLERR,arg[iarg+1]);
      vzhi = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix deposit command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix deposit command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"target") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix deposit command");
      tx = force->numeric(FLERR,arg[iarg+1]);
      ty = force->numeric(FLERR,arg[iarg+2]);
      tz = force->numeric(FLERR,arg[iarg+3]);
      targetflag = 1;
      iarg += 4;
    } else error->all(FLERR,"Illegal fix deposit command");
  }
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixDeposit::write_restart(FILE *fp)
{
  int n = 0;
  double list[4];
  list[n++] = random->state();
  list[n++] = ninserted;
  list[n++] = nfirst;
  list[n++] = next_reneighbor;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixDeposit::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  ninserted = static_cast<int> (list[n++]);
  nfirst = static_cast<int> (list[n++]);
  next_reneighbor = static_cast<int> (list[n++]);

  random->reset(seed);
}
