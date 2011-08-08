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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "pppm_cg.h"

using namespace LAMMPS_NS;

#define OFFSET 16384
#define SMALLQ 0.00001
#if defined(FFT_SINGLE)
#define ZEROF 0.0f
#else
#define ZEROF 0.0
#endif

/* ---------------------------------------------------------------------- */

PPPMCG::PPPMCG(LAMMPS *lmp, int narg, char **arg) : PPPM(lmp, narg, arg)
{
  if ((narg < 1) || (narg > 2)) error->all("Illegal kspace_style pppm/cg command");

  if (narg == 2)
    smallq = atof(arg[1]);
  else
    smallq = SMALLQ;

  num_charged = -1;
  is_charged = NULL;
}

/* ----------------------------------------------------------------------
   free all memory 
------------------------------------------------------------------------- */

PPPMCG::~PPPMCG()
{
  memory->destroy(is_charged);
}

/* ----------------------------------------------------------------------
   compute the PPPM long-range force, energy, virial 
------------------------------------------------------------------------- */

void PPPMCG::compute(int eflag, int vflag)
{
  int i;

  // convert atoms from box to lamda coords
  
  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // extend size of per-atom arrays if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(part2grid);
    memory->destroy(is_charged);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm:part2grid");
    memory->create(is_charged,nmax,"pppm/cg:is_charged");
  }

  // one time setup message.

  if (num_charged < 0) {
    bigint charged_all, charged_num;
    double charged_frac, charged_fmax, charged_fmin;

    num_charged=0;
    for (i=0; i < atom->nlocal; ++i)
      if (fabs(atom->q[i]) > smallq)
        ++num_charged;

    // get fraction of charged particles per domain

    if (atom->nlocal > 0)
      charged_frac = static_cast<double>(num_charged) * 100.0 
                   / static_cast<double>(atom->nlocal);
    else
      charged_frac = 0.0;
    MPI_Reduce(&charged_frac,&charged_fmax,1,MPI_DOUBLE,MPI_MAX,0,world);
    MPI_Reduce(&charged_frac,&charged_fmin,1,MPI_DOUBLE,MPI_MIN,0,world);

    // get fraction of charged particles overall

    charged_num = num_charged;
    MPI_Reduce(&charged_num,&charged_all,1,MPI_LMP_BIGINT,MPI_SUM,0,world);
    charged_frac = static_cast<double>(charged_all) * 100.0 
                   / static_cast<double>(atom->natoms);

    if (me == 0) {
      if (screen) fprintf(screen,"Using pppm/cg optimization. Cutoff: %g. \n"
			"Total charged: %.1f%%.  Min./Max. charged / proc: %.1f%%/%.1f%%.\n",
			smallq, charged_frac, charged_fmin, charged_fmax);
      if (logfile) fprintf(logfile,"Using pppm/cg optimization. Cutoff: %g. \n"
			"Total charged: %.1f%%.  Max. charged / proc: %.1f%%.\n",
			smallq, charged_frac, charged_fmax);
    }
  }

  num_charged=0;
  for (i=0; i < atom->nlocal; ++i)
    if (fabs(atom->q[i]) > smallq) {
      is_charged[num_charged] = i;
      ++num_charged;
    }

  energy = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid

  particle_map();
  make_rho();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  
  poisson(eflag,vflag);

  // all procs communicate E-field values to fill ghost cells
  //   surrounding their 3d bricks

  fillbrick();

  // calculate the force on my particles

  fieldforce();

  // sum energy across procs and add in volume-dependent term

  if (eflag) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;
   
    energy *= 0.5*volume;
    energy -= g_ewald*qsqsum/1.772453851 +
      0.5*PI*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qqrd2e*scale;
  }

  // sum virial across procs

  if (vflag) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qqrd2e*scale*volume*virial_all[i];
  }

  // 2d slab correction

  if (slabflag) slabcorr(eflag);

  // convert atoms back from lamda to box coords
  
  if (triclinic) domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array 
------------------------------------------------------------------------- */

void PPPMCG::particle_map()
{
  int nx,ny,nz;

  double **x = atom->x;

  int flag = 0;
  for (int j = 0; j < num_charged; j++) {
    int i = is_charged[j];
    
    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    nx = static_cast<int> ((x[i][0]-boxlo[0])*delxinv+shift) - OFFSET;
    ny = static_cast<int> ((x[i][1]-boxlo[1])*delyinv+shift) - OFFSET;
    nz = static_cast<int> ((x[i][2]-boxlo[2])*delzinv+shift) - OFFSET;

    part2grid[i][0] = nx;
    part2grid[i][1] = ny;
    part2grid[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
	ny+nlower < nylo_out || ny+nupper > nyhi_out ||
	nz+nlower < nzlo_out || nz+nupper > nzhi_out) flag = 1;
  }

  if (flag) error->one("Out of range atoms - cannot compute PPPM");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid 
------------------------------------------------------------------------- */

void PPPMCG::make_rho()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density array

  FFT_SCALAR *vec = &density_brick[nzlo_out][nylo_out][nxlo_out];
  for (i = 0; i < ngrid; i++) vec[i] = ZEROF;

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;

  for (int j = 0; j < num_charged; j++) {
    int i = is_charged[j];

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    z0 = delvolinv * q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
	my = m+ny;
	x0 = y0*rho1d[1][m];
	for (l = nlower; l <= nupper; l++) {
	  mx = l+nx;
	  density_brick[mz][my][mx] += x0*rho1d[0][l];
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles 
------------------------------------------------------------------------- */

void PPPMCG::fieldforce()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  for (int j = 0; j < num_charged; j++) {
    i = is_charged[j];

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
	my = m+ny;
	y0 = z0*rho1d[1][m];
	for (l = nlower; l <= nupper; l++) {
	  mx = l+nx;
	  x0 = y0*rho1d[0][l];
	  ekx -= x0*vdx_brick[mz][my][mx];
	  eky -= x0*vdy_brick[mz][my][mx];
	  ekz -= x0*vdz_brick[mz][my][mx];
	}
      }
    }

    // convert E-field to force
    const double qfactor = qqrd2e*scale*q[i];
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    f[i][2] += qfactor*ekz;
  }
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if 
   adequate empty space is left between repeating slabs (J. Chem. Phys. 
   111, 3155).  Slabs defined here to be parallel to the xy plane. 
------------------------------------------------------------------------- */

void PPPMCG::slabcorr(int eflag)
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;

  double dipole = 0.0;
  for (int j = 0; j < num_charged; j++) {
    int i = is_charged[j];
    dipole += q[i]*x[i][2];
  }

  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // compute corrections
  
  double e_slabcorr = 2.0*PI*dipole_all*dipole_all/volume;
  
  if (eflag) energy += qqrd2e*scale * e_slabcorr;

  // add on force corrections

  double ffact = -4.0*PI*dipole_all/volume * qqrd2e * scale; 
  double **f = atom->f;

  for (int j = 0; j < num_charged; j++) {
    int i = is_charged[j];
    f[i][2] += q[i]*ffact;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local arrays 
------------------------------------------------------------------------- */

double PPPMCG::memory_usage()
{
  double bytes = PPPM::memory_usage();
  bytes += nmax * sizeof(int);
  return bytes;
}
