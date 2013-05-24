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
#include "string.h"
#include "stdlib.h"
#include "fix_append_atoms.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "random_mars.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define BIG      1.0e30
#define EPSILON  1.0e-6

/* ---------------------------------------------------------------------- */

FixAppendAtoms::FixAppendAtoms(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  force_reneighbor = 1;
  next_reneighbor = -1;
  box_change = 1;
  time_depend = 1;

  if (narg < 4) error->all(FLERR,"Illegal fix append/atoms command");

  // default settings

  scaleflag = 1;
  spatflag=0;
  xloflag = xhiflag = yloflag = yhiflag = zloflag = zhiflag = 0;

  tempflag = 0;

  ranflag = 0;
  ranx = 0.0;
  rany = 0.0;
  ranz = 0.0;

  randomx = NULL;
  randomt = NULL;

  if (domain->lattice->nbasis == 0)
    error->all(FLERR,"Fix append/atoms requires a lattice be defined");

  nbasis = domain->lattice->nbasis;
  basistype = new int[nbasis];
  for (int i = 0; i < nbasis; i++) basistype[i] = 1;

  int iarg = 0;
  iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"xlo") == 0) {
      error->all(FLERR,"Only zhi currently implemented for fix append/atoms");
      xloflag = 1;
      iarg++;
      if (domain->boundary[0][0] != 3)
        error->all(FLERR,"Append boundary must be shrink/minimum");
    } else if (strcmp(arg[iarg],"xhi") == 0) {
      error->all(FLERR,"Only zhi currently implemented for fix append/atoms");
      xhiflag = 1;
      iarg++;
      if (domain->boundary[0][1] != 3)
        error->all(FLERR,"Append boundary must be shrink/minimum");
    } else if (strcmp(arg[iarg],"ylo") == 0) {
      error->all(FLERR,"Only zhi currently implemented for fix append/atoms");
      yloflag = 1;
      iarg++;
      if (domain->boundary[1][0] != 3)
        error->all(FLERR,"Append boundary must be shrink/minimum");
    } else if (strcmp(arg[iarg],"yhi") == 0) {
      error->all(FLERR,"Only zhi currently implemented for fix append/atoms");
      yhiflag = 1;
      iarg++;
      if (domain->boundary[1][1] != 3)
        error->all(FLERR,"Append boundary must be shrink/minimum");
    } else if (strcmp(arg[iarg],"zlo") == 0) {
      error->all(FLERR,"Only zhi currently implemented for fix append/atoms");
      zloflag = 1;
      iarg++;
      if (domain->boundary[2][0] != 3)
        error->all(FLERR,"Append boundary must be shrink/minimum");
    } else if (strcmp(arg[iarg],"zhi") == 0) {
      zhiflag = 1;
      iarg++;
      if (domain->boundary[2][1] != 3)
        error->all(FLERR,"Append boundary must be shrink/minimum");
    } else if (strcmp(arg[iarg],"freq") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix append/atoms command");
      freq = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"spatial") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix append/atoms command");
      if (strcmp(arg[iarg+1],"f_") == 0)
        error->all(FLERR,
                   "Bad fix ID in fix append/atoms command");
      spatflag = 1;
      int n = strlen(arg[iarg+1]);
      spatlead = atof(arg[iarg+2]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg+1][2]);
      n = strlen(suffix) + 1;
      spatialid = new char[n];
      strcpy(spatialid,suffix);
      delete [] suffix;
      iarg += 3;
    } else if (strcmp(arg[iarg],"basis") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix append/atoms command");
      int ibasis = atoi(arg[iarg+1]);
      int itype = atoi(arg[iarg+2]);
      if (ibasis <= 0 || ibasis > nbasis || itype <= 0 || itype > atom->ntypes)
        error->all(FLERR,"Invalid basis setting in fix append/atoms command");
      basistype[ibasis-1] = itype;
      iarg += 3;
    } else if (strcmp(arg[iarg],"size") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix append/atoms command");
      size = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix append/atoms command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix append/atoms command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"random") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix append/atoms command");
      ranflag = 1;
      ranx = atof(arg[iarg+1]);
      rany = atof(arg[iarg+2]);
      ranz = atof(arg[iarg+3]);
      xseed = atoi(arg[iarg+4]);
      if (xseed <= 0) error->all(FLERR,"Illegal fix append/atoms command");
      randomx = new RanMars(lmp,xseed + comm->me);
      iarg += 5;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix append/atoms command");
      tempflag = 1;
      t_target = atof(arg[iarg+1]);
      t_period = atof(arg[iarg+2]);
      tseed    = atoi(arg[iarg+3]);
      t_extent = atof(arg[iarg+4]);
      if (t_target <= 0) error->all(FLERR,"Illegal fix append/atoms command");
      if (t_period <= 0) error->all(FLERR,"Illegal fix append/atoms command");
      if (t_extent <= 0) error->all(FLERR,"Illegal fix append/atoms command");
      if (tseed <= 0) error->all(FLERR,"Illegal fix append/atoms command");
      randomt = new RanMars(lmp,tseed + comm->me);
      gfactor1 = new double[atom->ntypes+1];
      gfactor2 = new double[atom->ntypes+1];
      iarg += 5;
    } else error->all(FLERR,"Illegal fix append/atoms command");
  }

  if ((xloflag || xhiflag) && domain->xperiodic)
    error->all(FLERR,"Cannot use append/atoms in periodic dimension");
  if ((yloflag || yhiflag) && domain->yperiodic)
    error->all(FLERR,"Cannot use append/atoms in periodic dimension");
  if ((zloflag || zhiflag) && domain->zperiodic)
    error->all(FLERR,"Cannot use append/atoms in periodic dimension");

  if (domain->triclinic == 1)
    error->all(FLERR,"Cannot append atoms to a triclinic box");

  // setup scaling

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  if (xloflag || xhiflag) size *= xscale;
  if (yloflag || yhiflag) size *= yscale;
  if (zloflag || zhiflag) size *= zscale;

  if (ranflag) {
    ranx *= xscale;
    rany *= yscale;
    ranz *= zscale;
  }
}

/* ---------------------------------------------------------------------- */

FixAppendAtoms::~FixAppendAtoms()
{
  delete [] basistype;

  if (ranflag) delete randomx;
  if (tempflag) {
    delete randomt;
    delete [] gfactor1;
    delete [] gfactor2;
  }
}

/* ---------------------------------------------------------------------- */

int FixAppendAtoms::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAppendAtoms::initial_integrate(int vflag)
{
  if (update->ntimestep % freq == 0) next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixAppendAtoms::setup(int vflag)
{
  /*** CALL TO CREATE GROUP?  SEE POST_FORCE ***/
  post_force(vflag);
}


/* ---------------------------------------------------------------------- */

int FixAppendAtoms::get_spatial()
{
  if (update->ntimestep % freq == 0) {
    int ifix = modify->find_fix(spatialid);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for fix ave/spatial does not exist");
    Fix *fix = modify->fix[ifix];

    int failed = 0;
    int count = 0;
    while (failed < 2) {
      double tmp = fix->compute_vector(2*count);
      if (tmp == 0.0) failed++;
      else failed = 0;
      count++;
    }
    double *pos = new double[count-2];
    double *val = new double[count-2];
    for (int loop=0; loop < count-2; loop++) {
      pos[loop] = fix->compute_vector(2*loop);
      val[loop] = fix->compute_vector(2*loop+1);
    }

    // always ignore the first and last

    double binsize = 2.0;
    double min_energy=0.0;
    double max_energy=0.0;
    int header = static_cast<int> (size / binsize);
    advance = 0;

    for (int loop=1; loop <= header; loop++) {
        max_energy += val[loop];
    }
    for (int loop=count-2-2*header; loop <=count-3-header; loop++) {
      min_energy += val[loop];
    }
    max_energy /= header;
    min_energy /= header;

    double shockfront_min = 0.0;
    double shockfront_max = 0.0;
    double shockfront_loc = 0.0;
    int front_found1 = 0;
    for (int loop=count-3-header; loop > header; loop--) {
      if (front_found1 == 1) continue;
      if (val[loop] > min_energy + 0.1*(max_energy - min_energy)) {
        shockfront_max = pos[loop];
        front_found1=1;
      }
    }
    int front_found2 = 0;
    for (int loop=header+1; loop <=count-3-header; loop++) {
      if (val[loop] > min_energy + 0.6*(max_energy - min_energy)) {
        shockfront_min = pos[loop];
        front_found2=1;
      }
    }
    if      (front_found1 + front_found2 == 0) shockfront_loc = 0.0;
    else if (front_found1 + front_found2 == 1)
      shockfront_loc = shockfront_max + shockfront_min;
    else if (front_found1 == 1 && front_found2 == 1 &&
             shockfront_max-shockfront_min > spatlead/2.0)
      shockfront_loc = shockfront_max;
    else shockfront_loc = (shockfront_max + shockfront_min) / 2.0;
    if (comm->me == 0)
      printf("SHOCK: %g %g %g %g %g\n", shockfront_loc, shockfront_min,
             shockfront_max, domain->boxlo[2], domain->boxhi[2]);

    if (domain->boxhi[2] - shockfront_loc < spatlead) advance = 1;

    delete [] pos,val;
  }

  advance_sum = 0;
  MPI_Allreduce(&advance,&advance_sum,1,MPI_INT,MPI_SUM,world);

  if (advance_sum > 0) return 1;
  else return 0;
}

/* ---------------------------------------------------------------------- */

void FixAppendAtoms::post_force(int vflag)
{
  double **f = atom->f;
  double **v = atom->v;
  double **x = atom->x;
  int  *type = atom->type;
  int nlocal = atom->nlocal;

  double gamma1,gamma2;
  double tsqrt = sqrt(t_target);

  if (atom->mass) {
    if (tempflag) {
      for (int i = 1; i <= atom->ntypes; i++) {
        gfactor1[i] = -atom->mass[i] / t_period / force->ftm2v;
        gfactor2[i] = sqrt(atom->mass[i]) *
          sqrt(24.0*force->boltz/t_period/update->dt/force->mvv2e) /
          force->ftm2v;
      }
    }
    for (int i = 0; i < nlocal; i++) {
      // SET TEMP AHEAD OF SHOCK
      if (tempflag && x[i][2] >= domain->boxhi[2] - t_extent ) {
        gamma1 = gfactor1[type[i]];
        gamma2 = gfactor2[type[i]] * tsqrt;
        f[i][0] += gamma1*v[i][0] + gamma2*(randomt->uniform()-0.5);
        f[i][1] += gamma1*v[i][1] + gamma2*(randomt->uniform()-0.5);
        f[i][2] += gamma1*v[i][2] + gamma2*(randomt->uniform()-0.5);
      }
      // FREEZE ATOMS AT BOUNDARY
      if (x[i][2] >= domain->boxhi[2] - size) {
        f[i][0] = 0.0;
        f[i][1] = 0.0;
        f[i][2] = 0.0;
        v[i][0] = 0.0;
        v[i][1] = 0.0;
        v[i][2] = 0.0;
      }
    }
  } else {
    double *rmass = atom->rmass;
    double boltz = force->boltz;
    double dt = update->dt;
    double mvv2e = force->mvv2e;
    double ftm2v = force->ftm2v;

    for (int i = 0; i < nlocal; i++) {

      // set temp ahead of shock

      if (tempflag && x[i][2] >= domain->boxhi[2] - t_extent ) {
        gamma1 = -rmass[i] / t_period / ftm2v;
        gamma2 = sqrt(rmass[i]) * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
        gamma2 *= tsqrt;
        f[i][0] += gamma1*v[i][0] + gamma2*(randomt->uniform()-0.5);
        f[i][1] += gamma1*v[i][1] + gamma2*(randomt->uniform()-0.5);
        f[i][2] += gamma1*v[i][2] + gamma2*(randomt->uniform()-0.5);
      }

      // freeze atoms at boundary

      if (x[i][2] >= domain->boxhi[2] - size) {
        f[i][0] = 0.0;
        f[i][1] = 0.0;
        f[i][2] = 0.0;
        v[i][0] = 0.0;
        v[i][1] = 0.0;
        v[i][2] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAppendAtoms::pre_exchange()
{
  int ntimestep = update->ntimestep;
  int addnode = 0;

  if (ntimestep % freq == 0) {
    if (spatflag==1) if (get_spatial()==0) return;
    if (comm->myloc[2] == comm->procgrid[2]-1) {
      double bboxlo[3],bboxhi[3];

      bboxlo[0] = domain->sublo[0]; bboxhi[0] = domain->subhi[0];
      bboxlo[1] = domain->sublo[1]; bboxhi[1] = domain->subhi[1];
      bboxlo[2] = domain->subhi[2]; bboxhi[2] = domain->subhi[2]+size;

      double xmin,ymin,zmin,xmax,ymax,zmax;
      xmin = ymin = zmin = BIG;
      xmax = ymax = zmax = -BIG;

      domain->lattice->bbox(1,bboxlo[0],bboxlo[1],bboxlo[2],
                            xmin,ymin,zmin,xmax,ymax,zmax);
      domain->lattice->bbox(1,bboxhi[0],bboxlo[1],bboxlo[2],
                            xmin,ymin,zmin,xmax,ymax,zmax);
      domain->lattice->bbox(1,bboxlo[0],bboxhi[1],bboxlo[2],
                            xmin,ymin,zmin,xmax,ymax,zmax);
      domain->lattice->bbox(1,bboxhi[0],bboxhi[1],bboxlo[2],
                            xmin,ymin,zmin,xmax,ymax,zmax);
      domain->lattice->bbox(1,bboxlo[0],bboxlo[1],bboxhi[2],
                            xmin,ymin,zmin,xmax,ymax,zmax);
      domain->lattice->bbox(1,bboxhi[0],bboxlo[1],bboxhi[2],
                            xmin,ymin,zmin,xmax,ymax,zmax);
      domain->lattice->bbox(1,bboxlo[0],bboxhi[1],bboxhi[2],
                            xmin,ymin,zmin,xmax,ymax,zmax);
      domain->lattice->bbox(1,bboxhi[0],bboxhi[1],bboxhi[2],
                            xmin,ymin,zmin,xmax,ymax,zmax);

      int ilo,ihi,jlo,jhi,klo,khi;
      ilo = static_cast<int> (xmin);
      jlo = static_cast<int> (ymin);
      klo = static_cast<int> (zmin);
      ihi = static_cast<int> (xmax);
      jhi = static_cast<int> (ymax);
      khi = static_cast<int> (zmax);

      if (xmin < 0.0) ilo--;
      if (ymin < 0.0) jlo--;
      if (zmin < 0.0) klo--;

      double **basis = domain->lattice->basis;
      double x[3];
      double *sublo = domain->sublo;
      double *subhi = domain->subhi;
      double *mass = atom->mass;

      int i,j,k,m;
      for (k = klo; k <= khi; k++) {
        for (j = jlo; j <= jhi; j++) {
          for (i = ilo; i <= ihi; i++) {
            for (m = 0; m < nbasis; m++) {
              x[0] = i + basis[m][0];
              x[1] = j + basis[m][1];
              x[2] = k + basis[m][2];

              int flag = 0;
              // convert from lattice coords to box coords
              domain->lattice->lattice2box(x[0],x[1],x[2]);

              if (x[0] >= sublo[0] && x[0] < subhi[0] &&
                  x[1] >= sublo[1] && x[1] < subhi[1] &&
                  x[2] >= subhi[2] && x[2] < subhi[2]+size) flag = 1;
              else if (domain->dimension == 2 && x[1] >= domain->boxhi[1] &&
                       comm->myloc[1] == comm->procgrid[1]-1 &&
                       x[0] >= sublo[0] && x[0] < subhi[0]) flag = 1;

              if (flag) {
                if (ranflag) {
                  x[0] += ranx * 2.0*(randomx->uniform()-0.5);
                  x[1] += rany * 2.0*(randomx->uniform()-0.5);
                  x[2] += ranz * 2.0*(randomx->uniform()-0.5);
                }
                addnode++;
                atom->avec->create_atom(basistype[m],x);
              }
            }
          }
        }
      }
    }
    int addtotal = 0;
    MPI_Barrier(world);
    MPI_Allreduce(&addnode,&addtotal,1,MPI_INT,MPI_SUM,world);

    if (addtotal) {
      domain->reset_box();
      if (atom->tag_enable) {
        atom->tag_extend();
        atom->natoms += addtotal;
        if (atom->map_style) {
          atom->nghost = 0;
          atom->map_init();
          atom->map_set();
        }
      }
    }
  }
}
