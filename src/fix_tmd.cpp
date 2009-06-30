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
   Contributing authors: Paul Crozier (SNL)
                         Christian Burisch (Bochum Univeristy, Germany)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_tmd.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "group.h"
#include "respa.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define CHUNK 1000
#define MAXLINE 256

/* ---------------------------------------------------------------------- */

FixTMD::FixTMD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all("Illegal fix tmd command");

  rho_stop = atof(arg[3]);
  nfileevery = atoi(arg[5]);
  if (rho_stop < 0 || nfileevery < 0) error->all("Illegal fix tmd command");
  if (nfileevery && narg != 7) error->all("Illegal fix tmd command");

  MPI_Comm_rank(world,&me);

  // perform initial allocation of atom-based arrays
  // register with Atom class

  xf = NULL;
  xold = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // make sure an atom map exists before reading in target coordinates

  if (atom->map_style == 0) 
    error->all("Cannot use fix TMD unless atom map exists");

  // read from arg[4] and store coordinates of final target in xf

  readfile(arg[4]);

  // open arg[6] statistics file and write header

  if (nfileevery) {
    if (narg != 7) error->all("Illegal fix tmd command");
    if (me == 0) {
      fp = fopen(arg[6],"w");
      if (fp == NULL) {
	char str[128];
	sprintf(str,"Cannot open fix tmd file %s",arg[6]);
	error->one(str);
      }
      fprintf(fp,"%s %s\n","# Step rho_target rho_old gamma_back",
	      "gamma_forward lambda work_lambda work_analytical");
    }
  }

  masstotal = group->mass(igroup);

  // rho_start = initial rho
  // xold = initial x or 0.0 if not in group

  int *mask = atom->mask;
  int *type = atom->type;
  int *image = atom->image;
  double **x = atom->x;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double dx,dy,dz;
  int xbox,ybox,zbox;
  rho_start = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      dx = x[i][0] + xbox*xprd - xf[i][0];
      dy = x[i][1] + ybox*yprd - xf[i][1];
      dz = x[i][2] + zbox*zprd - xf[i][2];
      rho_start += mass[type[i]]*(dx*dx + dy*dy + dz*dz);
      xold[i][0] = x[i][0] + xbox*xprd;
      xold[i][1] = x[i][1] + ybox*yprd;
      xold[i][2] = x[i][2] + zbox*zprd;
    } else xold[i][0] = xold[i][1] = xold[i][2] = 0.0;
  }

  double rho_start_total;
  MPI_Allreduce(&rho_start,&rho_start_total,1,MPI_DOUBLE,MPI_SUM,world);
  rho_start = sqrt(rho_start_total/masstotal);
  rho_old = rho_start;

  work_lambda = 0.0;
  work_analytical = 0.0;
  previous_stat = 0;
}

/* ---------------------------------------------------------------------- */

FixTMD::~FixTMD()
{
  if (nfileevery && me == 0) fclose(fp);

  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored arrays

  memory->destroy_2d_double_array(xf);
  memory->destroy_2d_double_array(xold);
}

/* ---------------------------------------------------------------------- */

int FixTMD::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTMD::init()
{
  // check that no integrator fix comes after a TMD fix

  int flag = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"tmd") == 0) flag = 1;
    if (flag && strcmp(modify->fix[i]->style,"nve") == 0) flag = 2;
    if (flag && strcmp(modify->fix[i]->style,"nvt") == 0) flag = 2;
    if (flag && strcmp(modify->fix[i]->style,"npt") == 0) flag = 2;
    if (flag && strcmp(modify->fix[i]->style,"nph") == 0) flag = 2;
  }
  if (flag == 2) error->all("Fix tmd must come after integration fixes");

  // timesteps

  dtv = update->dt;
  dtf = update->dt * force->ftm2v;
  if (strcmp(update->integrate_style,"respa") == 0)
    step_respa = ((Respa *) update->integrate)->step;
}

/* ---------------------------------------------------------------------- */

void FixTMD::initial_integrate(int vflag)
{
  double a,b,c,d,e;
  double dx,dy,dz,dxkt,dykt,dzkt;
  double dxold,dyold,dzold,xback,yback,zback;
  double gamma_forward,gamma_back,gamma_max,lambda;
  double kt,fr,kttotal,frtotal,dtfm;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *image = atom->image;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int xbox,ybox,zbox;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  double rho_target = rho_start + delta * (rho_stop - rho_start);

  // compute the Lagrange multiplier

  a = b = e = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dxold = xold[i][0] - xf[i][0];
      dyold = xold[i][1] - xf[i][1];
      dzold = xold[i][2] - xf[i][2];
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      dx = x[i][0] + xbox*xprd - xf[i][0];
      dy = x[i][1] + ybox*yprd - xf[i][1];
      dz = x[i][2] + zbox*zprd - xf[i][2];
      a += mass[type[i]]*(dxold*dxold + dyold*dyold + dzold*dzold);
      b += mass[type[i]]*(dx   *dxold + dy   *dyold + dz   *dzold);
      e += mass[type[i]]*(dx   *dx    + dy   *dy    + dz   *dz);
    }
  }

  double abe[3],abetotal[3];
  abe[0] = a;  abe[1] = b;  abe[2] = e;
  MPI_Allreduce(abe,abetotal,3,MPI_DOUBLE,MPI_SUM,world);

  a = abetotal[0]/masstotal;
  b = 2.0*abetotal[1]/masstotal;
  e = abetotal[2]/masstotal;
  c = e - rho_old*rho_old;
  d = b*b - 4*a*c;  

  if (d < 0) d = 0;
  if (b >= 0) gamma_max = (-b - sqrt(d))/(2*a);
  else        gamma_max = (-b + sqrt(d))/(2*a);
  gamma_back = c/(a*gamma_max);
  if (a == 0.0) gamma_back = 0;

  c = e - rho_target*rho_target;
  d = b*b - 4*a*c;  
  if (d < 0) d = 0;
  if (b >= 0) gamma_max = (-b - sqrt(d))/(2*a);
  else        gamma_max = (-b + sqrt(d))/(2*a);
  gamma_forward = c/(a*gamma_max);
  if (a == 0.0) gamma_forward = 0;

  fr = kt = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dxold = xold[i][0] - xf[i][0];
      dyold = xold[i][1] - xf[i][1];
      dzold = xold[i][2] - xf[i][2];
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      xback = x[i][0] + xbox*xprd + gamma_back*dxold;
      yback = x[i][1] + ybox*yprd + gamma_back*dyold;
      zback = x[i][2] + zbox*zprd + gamma_back*dzold;
      dxkt = xback - xold[i][0];
      dykt = yback - xold[i][1];
      dzkt = zback - xold[i][2];
      kt += mass[type[i]]*(dxkt*dxkt + dykt*dykt + dzkt*dzkt);
      fr += f[i][0]*dxold + f[i][1]*dyold + f[i][2]*dzold;
    }
  }

  double r[2],rtotal[2];
  r[0] = fr;  r[1] = kt;
  MPI_Allreduce(r,rtotal,2,MPI_DOUBLE,MPI_SUM,world);
  frtotal = rtotal[0];
  kttotal = rtotal[1];

  // stat write of mean constraint force based on previous time step constraint

  if (nfileevery && me == 0) {   
    work_analytical += 
      (-frtotal - kttotal/dtv/dtf)*(rho_target - rho_old)/rho_old;
    lambda = gamma_back*rho_old*masstotal/dtv/dtf;
    work_lambda += lambda*(rho_target - rho_old);
    if (!(update->ntimestep % nfileevery) && 
	(previous_stat != update->ntimestep)) {
      fprintf(fp,"%d %g %g %g %g %g %g %g\n",
	      update->ntimestep,rho_target,rho_old,
              gamma_back,gamma_forward,lambda,work_lambda,work_analytical);
      fflush(fp);
      previous_stat = update->ntimestep;
    }
  }
  rho_old = rho_target;

  // apply the constraint and save constrained positions for next step

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      dxold = xold[i][0] - xf[i][0];
      x[i][0] += gamma_forward*dxold;
      v[i][0] += gamma_forward*dxold/dtv;
      f[i][0] += gamma_forward*dxold/dtv/dtfm;
      dyold = xold[i][1] - xf[i][1];
      x[i][1] += gamma_forward*dyold;
      v[i][1] += gamma_forward*dyold/dtv;
      f[i][1] += gamma_forward*dyold/dtv/dtfm;
      dzold = xold[i][2] - xf[i][2];
      x[i][2] += gamma_forward*dzold;
      v[i][2] += gamma_forward*dzold/dtv;
      f[i][2] += gamma_forward*dzold/dtv/dtfm;
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      xold[i][0] = x[i][0] + xbox*xprd;
      xold[i][1] = x[i][1] + ybox*yprd;
      xold[i][2] = x[i][2] + zbox*zprd;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTMD::initial_integrate_respa(int vflag, int ilevel, int flag)
{
  if (flag) return;             // only used by NPT,NPH

  dtv = step_respa[ilevel];
  dtf = step_respa[ilevel] * force->ftm2v;

  if (ilevel == 0) initial_integrate(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixTMD::memory_usage()
{
  double bytes = 2*atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based arrays
------------------------------------------------------------------------- */

void FixTMD::grow_arrays(int nmax)
{
  xf = memory->grow_2d_double_array(xf,nmax,3,"fix_tmd:xf");
  xold = memory->grow_2d_double_array(xold,nmax,3,"fix_tmd:xold");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixTMD::copy_arrays(int i, int j)
{
  xf[j][0] = xf[i][0];
  xf[j][1] = xf[i][1];
  xf[j][2] = xf[i][2];
  xold[j][0] = xold[i][0];
  xold[j][1] = xold[i][1];
  xold[j][2] = xold[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixTMD::pack_exchange(int i, double *buf)
{
  buf[0] = xf[i][0];
  buf[1] = xf[i][1];
  buf[2] = xf[i][2];
  buf[3] = xold[i][0];
  buf[4] = xold[i][1];
  buf[5] = xold[i][2];
  return 6;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixTMD::unpack_exchange(int nlocal, double *buf)
{
  xf[nlocal][0] = buf[0];
  xf[nlocal][1] = buf[1];
  xf[nlocal][2] = buf[2];
  xold[nlocal][0] = buf[3];
  xold[nlocal][1] = buf[4];
  xold[nlocal][2] = buf[5];
  return 6;
}

/* ----------------------------------------------------------------------
   read target coordinates from file, store with appropriate atom
------------------------------------------------------------------------- */

void FixTMD::readfile(char *file)
{
  if (me == 0) {
    if (screen) fprintf(screen,"Reading TMD target file %s ...\n",file);
    open(file);
  }

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  char *buffer = new char[CHUNK*MAXLINE];
  char *ptr,*next,*bufptr;
  int i,m,nlines,tag,imageflag,ix,iy,iz;
  double x,y,z,xprd,yprd,zprd;

  int firstline = 1;
  int ncount = 0;
  int eof = 0;
  xprd = yprd = zprd = -1.0;

  while (!eof) {
    if (me == 0) {
      m = 0;
      for (nlines = 0; nlines < CHUNK; nlines++) {
	ptr = fgets(&buffer[m],MAXLINE,fp);
	if (ptr == NULL) break;
	m += strlen(&buffer[m]);
      }
      if (ptr == NULL) eof = 1;
      buffer[m++] = '\n';
    }

    MPI_Bcast(&eof,1,MPI_INT,0,world);
    MPI_Bcast(&nlines,1,MPI_INT,0,world);
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);

    bufptr = buffer;
    for (i = 0; i < nlines; i++) {
      next = strchr(bufptr,'\n');
      *next = '\0';

      if (firstline) {
	if (strstr(bufptr,"xlo xhi")) {
	  double lo,hi;
	  sscanf(bufptr,"%lg %lg",&lo,&hi);
	  xprd = hi - lo;
	  bufptr = next + 1;
	  continue;
	} else if (strstr(bufptr,"ylo yhi")) {
	  double lo,hi;
	  sscanf(bufptr,"%lg %lg",&lo,&hi);
	  yprd = hi - lo;
	  bufptr = next + 1;
	  continue;
	} else if (strstr(bufptr,"zlo zhi")) {
	  double lo,hi;
	  sscanf(bufptr,"%lg %lg",&lo,&hi);
	  zprd = hi - lo;
	  bufptr = next + 1;
	  continue;
	} else if (atom->count_words(bufptr) == 4) {
	  if (xprd >= 0.0 || yprd >= 0.0 || zprd >= 0.0) 
	    error->all("Incorrect format in TMD target file");
	  imageflag = 0;
	  firstline = 0;
	} else if (atom->count_words(bufptr) == 7) {
	  if (xprd < 0.0 || yprd < 0.0 || zprd < 0.0) 
	    error->all("Incorrect format in TMD target file");
	  imageflag = 1;
	  firstline = 0;
	} else error->all("Incorrect format in TMD target file");
      }

      if (imageflag)
	sscanf(bufptr,"%d %lg %lg %lg %d %d %d",&tag,&x,&y,&z,&ix,&iy,&iz);
      else
	sscanf(bufptr,"%d %lg %lg %lg",&tag,&x,&y,&z);

      m = atom->map(tag);
      if (m >= 0 && m < nlocal && mask[m] & groupbit) {
	if (imageflag) {
	  xf[m][0] = x + ix*xprd;
	  xf[m][1] = y + iy*yprd;
	  xf[m][2] = z + iz*zprd;
	} else {
	  xf[m][0] = x;
	  xf[m][1] = y;
	  xf[m][2] = z;
	}
	ncount++;
      }

      bufptr = next + 1;
    }
  }

  // clean up

  delete [] buffer;

  if (me == 0) {
    if (compressed) pclose(fp);
    else fclose(fp);
  }

  // check that all atoms in group were listed in target file
  // set xf = 0.0 for atoms not in group

  int gcount = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) gcount++;
    else xf[i][0] = xf[i][1] = xf[i][2] = 0.0;

  int flag = 0;
  if (gcount != ncount) flag = 1;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all("TMD target file did not list all group atoms");
}

/* ----------------------------------------------------------------------
   proc 0 opens TMD data file
   test if gzipped
------------------------------------------------------------------------- */

void FixTMD::open(char *file)
{
  compressed = 0;
  char *suffix = file + strlen(file) - 3;
  if (suffix > file && strcmp(suffix,".gz") == 0) compressed = 1;
  if (!compressed) fp = fopen(file,"r");
  else {
#ifdef LAMMPS_GZIP
    char gunzip[128];
    sprintf(gunzip,"gunzip -c %s",file);
    fp = popen(gunzip,"r");
#else
    error->one("Cannot open gzipped file");
#endif
  }

  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(str);
  }
}

/* ---------------------------------------------------------------------- */

void FixTMD::reset_dt()
{
  dtv = update->dt;
  dtf = update->dt * force->ftm2v;
}
