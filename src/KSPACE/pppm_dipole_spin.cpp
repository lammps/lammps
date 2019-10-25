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
   Contributing author: Julien Tranchida (SNL)
------------------------------------------------------------------------- */

#include "pppm_dipole_spin.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "gridcomm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "update.h"

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXORDER 7
#define OFFSET 16384
#define LARGE 10000.0
#define SMALL 0.00001
#define EPS_HOC 1.0e-7

enum{REVERSE_MU};
enum{FORWARD_MU,FORWARD_MU_PERATOM};

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

PPPMDipoleSpin::PPPMDipoleSpin(LAMMPS *lmp) : 
  PPPMDipole(lmp)
{
  dipoleflag = 0;
  spinflag = 1;
  
  hbar = force->hplanck/MY_2PI;         	// eV/(rad.THz)
  mub = 9.274e-4;                     		// in A.Ang^2
  mu_0 = 785.15;               			// in eV/Ang/A^2
  mub2mu0 = mub * mub * mu_0 / (4.0*MY_PI);	// in eV.Ang^3
  mub2mu0hbinv = mub2mu0 / hbar;        	// in rad.THz
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

PPPMDipoleSpin::~PPPMDipoleSpin()
{
  if (copymode) return;

  deallocate();
  if (peratom_allocate_flag) deallocate_peratom();
  fft1 = NULL;
  fft2 = NULL;
  remap = NULL;
  cg_dipole = NULL;
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void PPPMDipoleSpin::init()
{
  if (me == 0) {
    if (screen) fprintf(screen,"PPPMDipoleSpin initialization ...\n");
    if (logfile) fprintf(logfile,"PPPMDipoleSpin initialization ...\n");
  }

  // error check

  spinflag = atom->sp?1:0;
  
  triclinic_check();

  if (triclinic != domain->triclinic)
    error->all(FLERR,"Must redefine kspace_style after changing to triclinic box");

  if (domain->dimension == 2) error->all(FLERR,
                                         "Cannot use PPPMDipoleSpin with 2d simulation");
  if (comm->style != 0)
    error->universe_all(FLERR,"PPPMDipoleSpin can only currently be used with "
                        "comm_style brick");

  if (!atom->sp) error->all(FLERR,"Kspace style requires atom attribute sp");

  if (atom->sp && differentiation_flag == 1) error->all(FLERR,"Cannot (yet) use"
     " kspace_modify diff ad with spins");

  if (spinflag && strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"'metal' units have to be used with spins");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use nonperiodic boundaries with PPPMDipoleSpin");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab PPPMDipoleSpin");
  }

  if (order < 2 || order > MAXORDER) {
    char str[128];
    sprintf(str,"PPPMDipoleSpin order cannot be < 2 or > than %d",MAXORDER);
    error->all(FLERR,str);
  }

  // extract short-range Coulombic cutoff from pair style

  triclinic = domain->triclinic;
  if (triclinic)
    error->all(FLERR,"Cannot yet use triclinic cells with PPPMDipoleSpin");

  pair_check();

  int itmp = 0;
  double *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  // check the correct extract here
  if (p_cutoff == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;

  // kspace TIP4P not yet supported
  // qdist = offset only for TIP4P fictitious charge
  
  qdist = 0.0; 
  if (tip4pflag)
    error->all(FLERR,"Cannot yet use TIP4P with PPPMDipoleSpin");

  scale = 1.0;
  spsum_spsq();
  natoms_original = atom->natoms;

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  // is two_charge_force still relevant for spin systems?

  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;

  // free all arrays previously allocated

  deallocate();
  if (peratom_allocate_flag) deallocate_peratom();

  // setup FFT grid resolution and g_ewald
  // normally one iteration thru while loop is all that is required
  // if grid stencil does not extend beyond neighbor proc
  //   or overlap is allowed, then done
  // else reduce order and try again

  int (*procneigh)[2] = comm->procneigh;

  GridComm *cgtmp = NULL;
  int iteration = 0;
    
  while (order >= minorder) {
    if (iteration && me == 0)
      error->warning(FLERR,"Reducing PPPMDipoleSpin order b/c stencil extends "
                     "beyond nearest neighbor processor");

    compute_gf_denom();
    set_grid_global();
    set_grid_local();
    if (overlap_allowed) break;

    cgtmp = new GridComm(lmp,world,1,1,
                         nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                         nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                         procneigh[0][0],procneigh[0][1],procneigh[1][0],
                         procneigh[1][1],procneigh[2][0],procneigh[2][1]);
    cgtmp->ghost_notify();
    if (!cgtmp->ghost_overlap()) break;
    delete cgtmp;

    order--;
    iteration++;
  }

  if (order < minorder) error->all(FLERR,"PPPMDipoleSpin order < minimum allowed order");
  if (!overlap_allowed && cgtmp->ghost_overlap())
    error->all(FLERR,"PPPMDipoleSpin grid stencil extends "
               "beyond nearest neighbor processor");
  if (cgtmp) delete cgtmp;

  // adjust g_ewald

  if (!gewaldflag) adjust_gewald();

  // calculate the final accuracy

  double estimated_accuracy = final_accuracy_dipole();

  // print stats

  int ngrid_max,nfft_both_max;
  MPI_Allreduce(&ngrid,&ngrid_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nfft_both,&nfft_both_max,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {

#ifdef FFT_SINGLE
    const char fft_prec[] = "single";
#else
    const char fft_prec[] = "double";
#endif

    if (screen) {
      fprintf(screen,"  G vector (1/distance) = %g\n",g_ewald);
      fprintf(screen,"  grid = %d %d %d\n",nx_pppm,ny_pppm,nz_pppm);
      fprintf(screen,"  stencil order = %d\n",order);
      fprintf(screen,"  estimated absolute RMS force accuracy = %g\n",
              estimated_accuracy);
      fprintf(screen,"  estimated relative force accuracy = %g\n",
              estimated_accuracy/two_charge_force);
      fprintf(screen,"  using %s precision FFTs\n",fft_prec);
      fprintf(screen,"  3d grid and FFT values/proc = %d %d\n",
              ngrid_max,nfft_both_max);
    }
    if (logfile) {
      fprintf(logfile,"  G vector (1/distance) = %g\n",g_ewald);
      fprintf(logfile,"  grid = %d %d %d\n",nx_pppm,ny_pppm,nz_pppm);
      fprintf(logfile,"  stencil order = %d\n",order);
      fprintf(logfile,"  estimated absolute RMS force accuracy = %g\n",
              estimated_accuracy);
      fprintf(logfile,"  estimated relative force accuracy = %g\n",
              estimated_accuracy/two_charge_force);
      fprintf(logfile,"  using %s precision FFTs\n",fft_prec);
      fprintf(logfile,"  3d grid and FFT values/proc = %d %d\n",
              ngrid_max,nfft_both_max);
    }
  }

  // allocate K-space dependent memory
  // don't invoke allocate peratom(), will be allocated when needed

  allocate();
  cg_dipole->ghost_notify();
  cg_dipole->setup();

  // pre-compute Green's function denominator expansion
  // pre-compute 1d charge distribution coefficients

  compute_gf_denom();
  compute_rho_coeff();
}

/* ----------------------------------------------------------------------
   compute the PPPMDipoleSpin long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPMDipoleSpin::compute(int eflag, int vflag)
{
  int i,j;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  if (vflag_atom)
    error->all(FLERR,"Cannot (yet) compute per-atom virial "
                       "with kspace style pppm/dipole/spin");

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    cg_peratom_dipole->ghost_notify();
    cg_peratom_dipole->setup();
  }

  // if atom count has changed, update qsum and qsqsum

  if (atom->natoms != natoms_original) {
    spsum_spsq();
    natoms_original = atom->natoms;
  }

  // return if there are no spins

  if (musqsum == 0.0) return;

  // convert atoms from box to lamda coords

  boxlo = domain->boxlo;

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm_spin:part2grid");
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d on-grid density

  particle_map();
  make_rho_spin();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  cg_dipole->reverse_comm(this,REVERSE_MU);
  brick2fft_dipole();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  // also performs per-atom calculations via poisson_peratom()

  poisson_ik_dipole();

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  cg_dipole->forward_comm(this,FORWARD_MU);

  // extra per-atom energy/virial communication

  if (evflag_atom) {
    cg_peratom_dipole->forward_comm(this,FORWARD_MU_PERATOM);
  }

  // calculate the force on my particles

  fieldforce_ik_spin();

  // extra per-atom energy/virial communication

  if (evflag_atom) fieldforce_peratom_spin();

  // sum global energy across procs and add in volume-dependent term

  const double spscale = mub2mu0 * scale;
  const double g3 = g_ewald*g_ewald*g_ewald;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    energy *= 0.5*volume;
    energy -= musqsum*2.0*g3/3.0/MY_PIS;
    energy *= spscale;
  }

  // sum global virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*spscale*volume*virial_all[i];
  }

  // per-atom energy/virial
  // energy includes self-energy correction

  if (evflag_atom) {
    double **sp = atom->sp;
    double spx,spy,spz;
    int nlocal = atom->nlocal;
    int ntotal = nlocal;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
	spx = sp[i][0]*sp[i][3];
	spy = sp[i][1]*sp[i][3];
	spz = sp[i][2]*sp[i][3];
        eatom[i] *= 0.5;
        eatom[i] -= (spx*spx + spy*spy + spz*spz)*2.0*g3/3.0/MY_PIS;
        eatom[i] *= spscale;
      }
    }

    if (vflag_atom) {
      for (i = 0; i < ntotal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5*spscale;
    }
  }

  // 2d slab correction

  if (slabflag == 1) slabcorr();
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void PPPMDipoleSpin::make_rho_spin()
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz;
  FFT_SCALAR x0,y0,z0;
  FFT_SCALAR x1,y1,z1;
  FFT_SCALAR x2,y2,z2;

  // clear 3d density array

  memset(&(densityx_brick_dipole[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));
  memset(&(densityy_brick_dipole[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));
  memset(&(densityz_brick_dipole[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double **sp = atom->sp;
  double spx,spy,spz;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    spx = sp[i][0]*sp[i][3];
    spy = sp[i][1]*sp[i][3];
    spz = sp[i][2]*sp[i][3];
    z0 = delvolinv * spx;
    z1 = delvolinv * spy;
    z2 = delvolinv * spz;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      y1 = z1*rho1d[2][n];
      y2 = z2*rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        x1 = y1*rho1d[1][m];
        x2 = y2*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          densityx_brick_dipole[mz][my][mx] += x0*rho1d[0][l];
          densityy_brick_dipole[mz][my][mx] += x1*rho1d[0][l];
          densityz_brick_dipole[mz][my][mx] += x2*rho1d[0][l];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get magnetic field & force on my particles for ik
------------------------------------------------------------------------- */

void PPPMDipoleSpin::fieldforce_ik_spin()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz;
  FFT_SCALAR x0,y0,z0;
  FFT_SCALAR ex,ey,ez;
  FFT_SCALAR vxx,vyy,vzz,vxy,vxz,vyz;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double **sp = atom->sp;
  double spx,spy,spz;
  double **x = atom->x;
  double **f = atom->f;
  double **fm_long = atom->fm_long;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    ex = ey = ez = ZEROF;
    vxx = vyy = vzz = vxy = vxz = vyz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          ex -= x0*ux_brick_dipole[mz][my][mx];
          ey -= x0*uy_brick_dipole[mz][my][mx];
          ez -= x0*uz_brick_dipole[mz][my][mx];
          vxx -= x0*vdxx_brick_dipole[mz][my][mx];
          vyy -= x0*vdyy_brick_dipole[mz][my][mx];
          vzz -= x0*vdzz_brick_dipole[mz][my][mx];
          vxy -= x0*vdxy_brick_dipole[mz][my][mx];
          vxz -= x0*vdxz_brick_dipole[mz][my][mx];
          vyz -= x0*vdyz_brick_dipole[mz][my][mx];
        }
      }
    }

    // convert M-field and store mech. forces

    const double spfactor = mub2mu0 * scale;
    spx = sp[i][0]*sp[i][3];
    spy = sp[i][1]*sp[i][3];
    spz = sp[i][2]*sp[i][3];
    f[i][0] += spfactor*(vxx*spx + vxy*spy + vxz*spz);
    f[i][1] += spfactor*(vxy*spx + vyy*spy + vyz*spz);
    f[i][2] += spfactor*(vxz*spx + vyz*spy + vzz*spz);

    // store long-range mag. precessions
    
    const double spfactorh = mub2mu0hbinv * scale;
    fm_long[i][0] += spfactorh*ex;
    fm_long[i][1] += spfactorh*ey;
    fm_long[i][2] += spfactorh*ez;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get per-atom energy/virial
------------------------------------------------------------------------- */

void PPPMDipoleSpin::fieldforce_peratom_spin()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ux,uy,uz;
  FFT_SCALAR v0x,v1x,v2x,v3x,v4x,v5x;
  FFT_SCALAR v0y,v1y,v2y,v3y,v4y,v5y;
  FFT_SCALAR v0z,v1z,v2z,v3z,v4z,v5z;

  // loop over my charges, interpolate from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double **sp = atom->sp;
  double spx,spy,spz;
  double **x = atom->x;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    ux = uy = uz = ZEROF;
    v0x = v1x = v2x = v3x = v4x = v5x = ZEROF;
    v0y = v1y = v2y = v3y = v4y = v5y = ZEROF;
    v0z = v1z = v2z = v3z = v4z = v5z = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          if (eflag_atom) {
            ux += x0*ux_brick_dipole[mz][my][mx];
            uy += x0*uy_brick_dipole[mz][my][mx];
            uz += x0*uz_brick_dipole[mz][my][mx];
          }
          if (vflag_atom) {
            v0x += x0*v0x_brick_dipole[mz][my][mx];
            v1x += x0*v1x_brick_dipole[mz][my][mx];
            v2x += x0*v2x_brick_dipole[mz][my][mx];
            v3x += x0*v3x_brick_dipole[mz][my][mx];
            v4x += x0*v4x_brick_dipole[mz][my][mx];
            v5x += x0*v5x_brick_dipole[mz][my][mx];
            v0y += x0*v0y_brick_dipole[mz][my][mx];
            v1y += x0*v1y_brick_dipole[mz][my][mx];
            v2y += x0*v2y_brick_dipole[mz][my][mx];
            v3y += x0*v3y_brick_dipole[mz][my][mx];
            v4y += x0*v4y_brick_dipole[mz][my][mx];
            v5y += x0*v5y_brick_dipole[mz][my][mx];
            v0z += x0*v0z_brick_dipole[mz][my][mx];
            v1z += x0*v1z_brick_dipole[mz][my][mx];
            v2z += x0*v2z_brick_dipole[mz][my][mx];
            v3z += x0*v3z_brick_dipole[mz][my][mx];
            v4z += x0*v4z_brick_dipole[mz][my][mx];
            v5z += x0*v5z_brick_dipole[mz][my][mx];
          }
        }
      }
    }

    spx = sp[i][0]*sp[i][3];
    spy = sp[i][1]*sp[i][3];
    spz = sp[i][2]*sp[i][3];
    if (eflag_atom) eatom[i] += spx*ux + spy*uy + spz*uz;
    if (vflag_atom) {
      vatom[i][0] += spx*v0x + spy*v0y + spz*v0z;
      vatom[i][1] += spx*v1x + spy*v1y + spz*v1z;
      vatom[i][2] += spx*v2x + spy*v2y + spz*v2z;
      vatom[i][3] += spx*v3x + spy*v3y + spz*v3z;
      vatom[i][4] += spx*v4x + spy*v4y + spz*v4z;
      vatom[i][5] += spx*v5x + spy*v5y + spz*v5z;
    }
  }
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
------------------------------------------------------------------------- */

void PPPMDipoleSpin::slabcorr()
{
  // compute local contribution to global spin moment

  double spin = 0.0;
  double **sp = atom->sp;
  double spz;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) { 
    spz = sp[i][2]*sp[i][3];
    spin += spz;
  }

  // sum local contributions to get global spin moment

  double spin_all;
  MPI_Allreduce(&spin,&spin_all,1,MPI_DOUBLE,MPI_SUM,world);

  // compute corrections

  const double e_slabcorr = MY_2PI*(spin_all*spin_all/12.0)/volume;
  const double spscale = mub2mu0 * scale;

  if (eflag_global) energy += spscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = spscale * MY_2PI/volume/12.0;
    for (int i = 0; i < nlocal; i++) {
      spz = sp[i][2]*sp[i][3];
      eatom[i] += efact * spz * spin_all;
    }
  }

  // add on mag. force corrections

  double ffact = spscale * (-4.0*MY_PI/volume);
  double **fm_long = atom->fm_long;
  for (int i = 0; i < nlocal; i++) {
    fm_long[i][2] += ffact * spin_all;
  }
}

/* ----------------------------------------------------------------------
   compute spsum,spsqsum,sp2
   called initially, when particle count changes, when spins are changed
------------------------------------------------------------------------- */

void PPPMDipoleSpin::spsum_spsq()
{
  const int nlocal = atom->nlocal;

  musum = musqsum = mu2 = 0.0;
  if (atom->sp_flag) {
    double **sp = atom->sp;
    double spx, spy, spz;
    double spsum_local(0.0), spsqsum_local(0.0);

    // sum (direction x norm) of all spins

    for (int i = 0; i < nlocal; i++) {
      spx = sp[i][0]*sp[i][3];
      spy = sp[i][1]*sp[i][3];
      spz = sp[i][2]*sp[i][3];
      spsum_local += spx + spy + spz;
      spsqsum_local += spx*spx + spy*spy + spz*spz;
    }

    // store results into pppm_dipole quantities 

    MPI_Allreduce(&spsum_local,&musum,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&spsqsum_local,&musqsum,1,MPI_DOUBLE,MPI_SUM,world);

    //mu2 = musqsum * mub2mu0;
    mu2 = musqsum;
  }

  if (mu2 == 0 && comm->me == 0)
    error->all(FLERR,"Using kspace solver PPPMDipoleSpin on system with no spins");
}
