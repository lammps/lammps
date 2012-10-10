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
   Contributing authors: Paul Crozier, Stan Moore, Stephen Bond, (all SNL)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "msm.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "modify.h"
#include "fix.h"

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAX_LEVELS 10
#define OFFSET 16384
#define SMALL 0.00001
#define LARGE 10000.0

/* ---------------------------------------------------------------------- */

MSM::MSM(LAMMPS *lmp, int narg, char **arg) : KSpace(lmp, narg, arg)
{
  if (narg < 1) error->all(FLERR,"Illegal kspace_style msm command");

  msmflag = 1;

  accuracy_relative = atof(arg[0]);

  nfactors = 1;
  factors = new int[nfactors];
  factors[0] = 2;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  buf1 = buf2 = buf3 = buf4 = NULL;

  phi1d = dphi1d = NULL;

  nmax = 0;
  part2grid = NULL;

  g_direct = NULL;

  dgx_direct = NULL;
  dgy_direct = NULL;
  dgz_direct = NULL;

  v0_direct = v1_direct = v2_direct = NULL;
  v3_direct = v4_direct = v5_direct = NULL;

  levels = 0;
  
  order = 4;

  differentiation_flag = 1;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

MSM::~MSM()
{
  delete [] factors;
  deallocate();
  deallocate_peratom();
  memory->destroy(part2grid);
  memory->destroy(g_direct);
  memory->destroy(dgx_direct);
  memory->destroy(dgy_direct);
  memory->destroy(dgz_direct);
  memory->destroy(v0_direct);
  memory->destroy(v1_direct);
  memory->destroy(v2_direct);
  memory->destroy(v3_direct);
  memory->destroy(v4_direct);
  memory->destroy(v5_direct);
  deallocate_levels();
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void MSM::init()
{
  if (me == 0) {
    if (screen) fprintf(screen,"MSM initialization ...\n");
    if (logfile) fprintf(logfile,"MSM initialization ...\n");
  }

  // error check

  if (domain->triclinic)
    error->all(FLERR,"Cannot (yet) use MSM with triclinic box");

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot (yet) use MSM with 2d simulation");

  if (!atom->q_flag) error->all(FLERR,"Kspace style requires atom attribute q");

  if (slabflag == 1)
      error->all(FLERR,"Cannot use slab correction with MSM");

  if (domain->nonperiodic > 0)
    error->all(FLERR,"Cannot (yet) use nonperiodic boundaries with MSM");

  if (order < 4 || order > 10) {
    char str[128];
    sprintf(str,"MSM order must be 4, 6, 8, or 10");
    error->all(FLERR,str);
  }

  if (order%2 != 0) error->all(FLERR,"MSM order must be 4, 6, 8, or 10");

  // free all arrays previously allocated

  deallocate();
  deallocate_peratom();
  peratom_allocate_flag = 0;

  // extract short-range Coulombic cutoff from pair style

  qqrd2e = force->qqrd2e;
  scale = 1.0;

  pair_check();

  int itmp;
  double *p_cutoff = (double *) force->pair->extract("cut_msm",itmp);
  if (p_cutoff == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;

  // compute qsum & qsqsum and give error if not charge-neutral

  qsum = qsqsum = 0.0;
  for (int i = 0; i < atom->nlocal; i++) {
    qsum += atom->q[i];
    qsqsum += atom->q[i]*atom->q[i];
  }

  double tmp;
  MPI_Allreduce(&qsum,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsum = tmp;
  MPI_Allreduce(&qsqsum,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsqsum = tmp;
  q2 = qsqsum * force->qqrd2e / force->dielectric;

  if (qsqsum == 0.0)
    error->all(FLERR,"Cannot use kspace solver on system with no charge");

  // not yet sure of the correction needed for non-neutral systems

  if (fabs(qsum) > SMALL) {
    char str[128];
    sprintf(str,"System is not charge neutral, net charge = %g",qsum);
    error->all(FLERR,str);  
  }

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;

  // setup grid resolution

  set_grid();

  int flag_global = 0;

  if (nx_msm[0] >= OFFSET || ny_msm[0] >= OFFSET || nz_msm[0] >= OFFSET)
    error->all(FLERR,"MSM grid is too large");

  // loop over grid levels

  for (int n=0; n<levels; n++) {

    // global indices of MSM grid range from 0 to N-1
    // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
    //   global MSM grid that I own without ghost cells

    if (n == 0) {

      nxlo_in_d[n] = static_cast<int> (comm->xsplit[comm->myloc[0]] * nx_msm[n]);
      nxhi_in_d[n] = static_cast<int> (comm->xsplit[comm->myloc[0]+1] * nx_msm[n]) - 1;

      nylo_in_d[n] = static_cast<int> (comm->ysplit[comm->myloc[1]] * ny_msm[n]);
      nyhi_in_d[n] = static_cast<int> (comm->ysplit[comm->myloc[1]+1] * ny_msm[n]) - 1;

      nzlo_in_d[n] = static_cast<int> (comm->zsplit[comm->myloc[2]] * nz_msm[n]);
      nzhi_in_d[n] = static_cast<int> (comm->zsplit[comm->myloc[2]+1] * nz_msm[n]) - 1;

    } else {

      nxlo_in_d[n] = 0;
      nxhi_in_d[n] = nx_msm[n] - 1;

      nylo_in_d[n] = 0;
      nyhi_in_d[n] = ny_msm[n] - 1;

      nzlo_in_d[n] = 0;
      nzhi_in_d[n] = nz_msm[n] - 1;
    }

    // Use simple method of parallel communication for now

    nxlo_in[n] = 0;
    nxhi_in[n] = nx_msm[n] - 1;

    nylo_in[n] = 0;
    nyhi_in[n] = ny_msm[n] - 1;

    nzlo_in[n] = 0;
    nzhi_in[n] = nz_msm[n] - 1;

    // nlower,nupper = stencil size for mapping particles to MSM grid

    nlower = -(order-1)/2;
    nupper = order/2;

    // shift values for particle <-> grid mapping
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    // nlo_out,nhi_out = lower/upper limits of the 3d sub-brick of
    //   global MSM grid that my particles can contribute charge to
    // effectively nlo_in,nhi_in + ghost cells
    // nlo,nhi = global coords of grid pt to "lower left" of smallest/largest
    //           position a particle in my box can be at
    // dist[3] = particle position bound = subbox + skin/2.0
    // nlo_out,nhi_out = nlo,nhi + stencil size for particle mapping

    double *prd,*sublo,*subhi;

    //prd = domain->prd;
    //boxlo = domain->boxlo;
    //sublo = domain->sublo;
    //subhi = domain->subhi;

    // Use only one partition for now

    prd = domain->prd;
    boxlo = domain->boxlo;
    sublo = boxlo;
    subhi = domain->boxhi;

    double xprd = prd[0];
    double yprd = prd[1];
    double zprd = prd[2];

    double dist[3];
    double cuthalf = 0.0;
    if (n == 0) cuthalf = 0.5*neighbor->skin; // Only applies to finest grid
    dist[0] = dist[1] = dist[2] = cuthalf;

    int nlo,nhi;

    nlo = static_cast<int> ((sublo[0]-dist[0]-boxlo[0]) *
                            nx_msm[n]/xprd + OFFSET) - OFFSET;
    nhi = static_cast<int> ((subhi[0]+dist[0]-boxlo[0]) *
                            nx_msm[n]/xprd + OFFSET) - OFFSET;
    nxlo_out[n] = nlo + nlower;
    nxhi_out[n] = nhi + nupper;

    nlo = static_cast<int> ((sublo[1]-dist[1]-boxlo[1]) *
                            ny_msm[n]/yprd + OFFSET) - OFFSET;
    nhi = static_cast<int> ((subhi[1]+dist[1]-boxlo[1]) *
                            ny_msm[n]/yprd + OFFSET) - OFFSET;
    nylo_out[n] = nlo + nlower;
    nyhi_out[n] = nhi + nupper;

    nlo = static_cast<int> ((sublo[2]-dist[2]-boxlo[2]) *
                            nz_msm[n]/zprd + OFFSET) - OFFSET;
    nhi = static_cast<int> ((subhi[2]+dist[2]-boxlo[2]) *
                            nz_msm[n]/zprd + OFFSET) - OFFSET;
    nzlo_out[n] = nlo + nlower;
    nzhi_out[n] = nhi + nupper;

    // nlo_ghost,nhi_ghost = # of planes I will recv from 6 directions
    //   that overlay domain I own
    // proc in that direction tells me via sendrecv()
    // if no neighbor proc, value is from self since I have ghosts regardless

    int nplanes;
    MPI_Status status;

    nplanes = nxlo_in[n] - nxlo_out[n];
    if (comm->procneigh[0][0] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[0][0],0,
                   &nxhi_ghost[n],1,MPI_INT,comm->procneigh[0][1],0,
                   world,&status);
    else nxhi_ghost[n] = nplanes;

    nplanes = nxhi_out[n] - nxhi_in[n];
    if (comm->procneigh[0][1] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[0][1],0,
                   &nxlo_ghost[n],1,MPI_INT,comm->procneigh[0][0],
                   0,world,&status);
    else nxlo_ghost[n] = nplanes;

    nplanes = nylo_in[n] - nylo_out[n];
    if (comm->procneigh[1][0] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[1][0],0,
                   &nyhi_ghost[n],1,MPI_INT,comm->procneigh[1][1],0,
                   world,&status);
    else nyhi_ghost[n] = nplanes;

    nplanes = nyhi_out[n] - nyhi_in[n];
    if (comm->procneigh[1][1] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[1][1],0,
                   &nylo_ghost[n],1,MPI_INT,comm->procneigh[1][0],0,
                   world,&status);
    else nylo_ghost[n] = nplanes;

    nplanes = nzlo_in[n] - nzlo_out[n];
    if (comm->procneigh[2][0] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[2][0],0,
                   &nzhi_ghost[n],1,MPI_INT,comm->procneigh[2][1],0,
                   world,&status);
    else nzhi_ghost[n] = nplanes;

    nplanes = nzhi_out[n] - nzhi_in[n];
    if (comm->procneigh[2][1] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[2][1],0,
                   &nzlo_ghost[n],1,MPI_INT,comm->procneigh[2][0],0,
                   world,&status);
    else nzlo_ghost[n] = nplanes;

    int flag = 0;

    if (n == 0) {

      // test that ghost overlap is not bigger than my sub-domain

      if (nxlo_ghost[n] > nxhi_in[n]-nxlo_in[n]+1) flag = 1;
      if (nxhi_ghost[n] > nxhi_in[n]-nxlo_in[n]+1) flag = 1;
      if (nylo_ghost[n] > nyhi_in[n]-nylo_in[n]+1) flag = 1;
      if (nyhi_ghost[n] > nyhi_in[n]-nylo_in[n]+1) flag = 1;
      if (nzlo_ghost[n] > nzhi_in[n]-nzlo_in[n]+1) flag = 1;
      if (nzhi_ghost[n] > nzhi_in[n]-nzlo_in[n]+1) flag = 1;
    }

    int flag_all;
    MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);

    if (flag_all != 0) {
      char str[128];
      sprintf(str,"MSM parallel communication error, try reducing number of procs");
      error->all(FLERR,str);
    }

  // MSM grids for this proc, including ghosts

  ngrid[n] = (nxhi_out[n]-nxlo_out[n]+1) * (nyhi_out[n]-nylo_out[n]+1) *
    (nzhi_out[n]-nzlo_out[n]+1);
  }


  // buffer space for use in ghost_swap and fillbrick
  // idel = max # of ghost planes to send or recv in +/- dir of each dim
  // nx,ny,nz = owned planes (including ghosts) in each dim
  // nxx,nyy,nzz = max # of grid cells to send in each dim
  // nbuf = max in any dim, augment by 3x for components of vd_xyz in fillbrick

  int idelx,idely,idelz,nx,ny,nz,nxx,nyy,nzz;

  idelx = MAX(nxlo_ghost[0],nxhi_ghost[0]);
  idelx = MAX(idelx,nxhi_out[0]-nxhi_in[0]);
  idelx = MAX(idelx,nxlo_in[0]-nxlo_out[0]);

  idely = MAX(nylo_ghost[0],nyhi_ghost[0]);
  idely = MAX(idely,nyhi_out[0]-nyhi_in[0]);
  idely = MAX(idely,nylo_in[0]-nylo_out[0]);

  idelz = MAX(nzlo_ghost[0],nzhi_ghost[0]);
  idelz = MAX(idelz,nzhi_out[0]-nzhi_in[0]);
  idelz = MAX(idelz,nzlo_in[0]-nzlo_out[0]);

  nx = nxhi_out[0] - nxlo_out[0] + 1;
  ny = nyhi_out[0] - nylo_out[0] + 1;
  nz = nzhi_out[0] - nzlo_out[0] + 1;

  nxx = idelx * ny * nz;
  nyy = idely * nx * nz;
  nzz = idelz * nx * ny;

  nbuf = MAX(nxx,nyy);
  nbuf = MAX(nbuf,nzz);

  nbuf_peratom = 7*nbuf;
  if (!differentiation_flag) {
    nbuf *= 3;
  }

  double estimated_error = estimate_total_error();

  // print stats

  int ngrid_max,nbuf_max;

  // All processors have a copy of the complete grid at each level

  nbuf_max = nbuf;
  ngrid_max = ngrid[0];

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  brick buffer size/proc = %d %d\n",
                        ngrid_max,nbuf_max);
      fprintf(screen,"  estimated absolute RMS force accuracy = %g\n",
              estimated_error);
      fprintf(screen,"  estimated relative force accuracy = %g\n",
              estimated_error/two_charge_force);
    }
    if (logfile) {
      fprintf(logfile,"  brick buffer size/proc = %d %d\n",
                         ngrid_max,nbuf_max);
      fprintf(logfile,"  estimated absolute RMS force accuracy = %g\n",
              estimated_error);
      fprintf(logfile,"  estimated relative force accuracy = %g\n",
              estimated_error/two_charge_force);
    }
  }

  // allocate K-space dependent memory

  allocate();
}

/* ----------------------------------------------------------------------
   estimate 1d grid RMS force error for MSM 
------------------------------------------------------------------------- */

double MSM::estimate_1d_error(double h, double prd)
{
  double a = cutoff;
  int p = order - 1;
  
  double Mp,cprime,error_scaling; 
  Mp = cprime = error_scaling = 1;
  // Mp values from Table 5.1 of Hardy's thesis
  // cprime values from equation 4.17 of Hardy's thesis
  // error scaling from empirical fitting to convert to rms force errors
  if (p == 3) {
    Mp = 9;
    cprime = 1.0/6.0;
    if (differentiation_flag) error_scaling = 0.39189561;
    else error_scaling = 0.215328372;
  } else if (p == 5) {
    Mp = 825;
    cprime = 1.0/30.0;
    if (differentiation_flag) error_scaling = 0.150829428;
    else error_scaling = 0.10751471;
  } else if (p == 7) {
    Mp = 130095;
    cprime = 1.0/140.0;
    if (differentiation_flag) error_scaling = 0.049632967;
    else error_scaling = 0.047579461;
  } else if (p == 9) {
    Mp = 34096545;
    cprime = 1.0/630.0;
    if (differentiation_flag) error_scaling = 0.013520855;
    else error_scaling = 0.010403771;
  } else {
    error->all(FLERR,"MSM order must be 4, 6, 8, or 10");
  }
  
  // equation 4.1 from Hardy's thesis
  double C_p = 4.0*cprime*Mp/3.0;
  
  // use empirical parameters to convert to rms force errors
  C_p *= error_scaling;
  
  // equation 3.197 from Hardy's thesis
  double error_1d = C_p*pow(h,(p-1))/pow(a,(p+1));
  
  // include dependency of error on other terms
  error_1d *= q2*a/(prd*sqrt(atom->natoms));
  
  return error_1d;
}

/* ----------------------------------------------------------------------
   estimate 3d grid RMS force error
------------------------------------------------------------------------- */

double MSM::estimate_3d_error()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double error_x = estimate_1d_error(xprd/nx_msm[0],xprd);
  double error_y = estimate_1d_error(yprd/ny_msm[0],yprd);
  double error_z = estimate_1d_error(zprd/nz_msm[0],zprd);
  double error_3d =
   sqrt(error_x*error_x + error_y*error_y + error_z*error_z) / sqrt(3.0);
  return error_3d;
}

/* ----------------------------------------------------------------------
   estimate total RMS force error
------------------------------------------------------------------------- */

double MSM::estimate_total_error()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  bigint natoms = atom->natoms;

  double grid_error = estimate_3d_error();
  double q2_over_sqrt = q2 / sqrt(natoms*cutoff*xprd*yprd*zprd);
  double short_range_error = 0.0;
  double table_error =
   estimate_table_accuracy(q2_over_sqrt,short_range_error);
  double estimated_total_error = sqrt(grid_error*grid_error +
   short_range_error*short_range_error + table_error*table_error);

  return estimated_total_error;
}

/* ----------------------------------------------------------------------
   adjust MSM coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void MSM::setup()
{
  int i,j,k,l,m,n;
  double *prd;

  double a = cutoff;

  // volume-dependent factors

  prd = domain->prd;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  volume = xprd * yprd * zprd;

  // loop over grid levels

  for (int n=0; n<levels; n++) {

    delxinv[n] = nx_msm[n]/xprd;
    delyinv[n] = ny_msm[n]/yprd;
    delzinv[n] = nz_msm[n]/zprd;

    delvolinv[n] = delxinv[n]*delyinv[n]*delzinv[n];
  }

  nxhi_direct = static_cast<int> (2.0*a*delxinv[0]);
  nxlo_direct = -nxhi_direct;
  nyhi_direct = static_cast<int> (2.0*a*delyinv[0]);
  nylo_direct = -nyhi_direct;
  nzhi_direct = static_cast<int> (2.0*a*delzinv[0]);
  nzlo_direct = -nzhi_direct;

  nmax_direct = 8*(nxhi_direct+1)*(nyhi_direct+1)*(nzhi_direct+1);

  if (!peratom_allocate_flag) { // Timestep 0
    get_g_direct();
    get_dg_direct();
    get_virial_direct();
  } else {
    if (differentiation_flag) get_g_direct();
    if (!differentiation_flag) get_dg_direct();

    if (!differentiation_flag && eflag_either) get_g_direct();
    if (differentiation_flag && vflag_either) get_dg_direct();
    if (vflag_either) get_virial_direct();
  }
  boxlo = domain->boxlo;
}

/* ----------------------------------------------------------------------
   compute the MSM long-range force, energy, virial
------------------------------------------------------------------------- */

void MSM::compute(int eflag, int vflag)
{
  int i,j;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
    eflag_atom = vflag_atom = eflag_either = vflag_either = 0;

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    peratom_allocate_flag = 1;
  }

  // extend size of per-atom arrays if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"msm:part2grid");
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid (aninterpolation)

  particle_map();
  make_rho();


  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks

  ghost_swap(0);

  grid_swap(0,qgrid[0]);

  // Direct sum on finest grid level is parallel

  if (differentiation_flag) direct_ad(0);
  else direct(0);

  if (evflag_atom) direct_peratom(0);

  if (differentiation_flag || eflag_atom)
    grid_swap(0,egrid[0]);

  if (!differentiation_flag) {
    grid_swap(0,fxgrid[0]);
    grid_swap(0,fygrid[0]);
    grid_swap(0,fzgrid[0]);
  }

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;
  }

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = virial_all[i];
  }

  if (vflag_atom) {
    grid_swap(0,v0grid[0]);
    grid_swap(0,v1grid[0]);
    grid_swap(0,v2grid[0]);
    grid_swap(0,v3grid[0]);
    grid_swap(0,v4grid[0]);
    grid_swap(0,v5grid[0]);
  }

  restriction(0);

  // compute potential gradient on my MSM grid and
  //   portion of e_long on this proc's MSM grid
  // return gradients (electric fields) in 3d brick decomposition

  for (int n=1; n<levels; n++) {
    if (differentiation_flag) direct_ad(n);
    else direct(n);

    if (evflag_atom) direct_peratom(n);

    if (n < levels-1) restriction(n);
  }

  for (int n=levels-2; n>=0; n--)
    prolongation(n);

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  if (differentiation_flag || evflag_atom) fillbrick_ad_peratom(0);

  if (!differentiation_flag) fillbrick(0);

  // calculate the force on my particles (interpolation)

  if (differentiation_flag) fieldforce_ad();
  else fieldforce();

  // calculate the per-atom energy for my particles

  if (evflag_atom) fieldforce_peratom();

  const double qscale = force->qqrd2e * scale;

  // Total long-range energy

  if (eflag_global) {
    double e_self = qsqsum*gamma(0.0)/cutoff;  // Self-energy term
    energy -= e_self;
    energy *= 0.5*qscale;
  }

  // Total long-range virial

  if (vflag_global) {
    for (i = 0; i < 6; i++) virial[i] *= 0.5*qscale;
  }

  // per-atom energy/virial
  // energy includes self-energy correction

  if (evflag_atom) {
    double *q = atom->q;
    int nlocal = atom->nlocal;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
        eatom[i] -= q[i]*q[i]*gamma(0.0)/cutoff;
        eatom[i] *= 0.5*qscale;
      }
    }

    if (vflag_atom) {
      for (i = 0; i < nlocal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5*qscale;
    }
  }

}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of grid points
------------------------------------------------------------------------- */

void MSM::allocate()
{
  memory->create(buf1,nbuf,"msm:buf1");
  memory->create(buf2,nbuf,"msm:buf2");

  if (differentiation_flag) {
    memory->create(buf3,nbuf_peratom,"msm:buf3");
    memory->create(buf4,nbuf_peratom,"msm:buf4");
  }

  // summation coeffs

  memory->create2d_offset(phi1d,3,-order+1,order-1,"msm:phi1d");
  memory->create2d_offset(dphi1d,3,-order+1,order-1,"msm:dphi1d");

  // allocate grid levels

  for (int n=0; n<levels; n++) {
    memory->create3d_offset(qgrid[n],nzlo_out[n],nzhi_out[n],
            nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:qgrid");

    if (differentiation_flag) {
      memory->create3d_offset(egrid[n],nzlo_out[n],nzhi_out[n],
              nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:egrid");
    } else {
      memory->create3d_offset(fxgrid[n],nzlo_out[n],nzhi_out[n],
              nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:fxgrid");
      memory->create3d_offset(fygrid[n],nzlo_out[n],nzhi_out[n],
              nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:fygrid");
      memory->create3d_offset(fzgrid[n],nzlo_out[n],nzhi_out[n],
              nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:fzgrid");
    }
  }
}

/* ----------------------------------------------------------------------
   allocate per-atom memory that depends on # of grid points
------------------------------------------------------------------------- */

void MSM::allocate_peratom()
{
  // allocate grid levels

  for (int n=0; n<levels; n++) {
    if (!differentiation_flag)
      memory->create3d_offset(egrid[n],nzlo_out[n],nzhi_out[n],
              nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:egrid");

    memory->create3d_offset(v0grid[n],nzlo_out[n],nzhi_out[n],
            nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:v0grid");
    memory->create3d_offset(v1grid[n],nzlo_out[n],nzhi_out[n],
            nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:v1grid");
    memory->create3d_offset(v2grid[n],nzlo_out[n],nzhi_out[n],
            nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:v2grid");
    memory->create3d_offset(v3grid[n],nzlo_out[n],nzhi_out[n],
            nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:v3grid");
    memory->create3d_offset(v4grid[n],nzlo_out[n],nzhi_out[n],
            nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:v4grid");
    memory->create3d_offset(v5grid[n],nzlo_out[n],nzhi_out[n],
            nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:v5grid");
  }

  if (!differentiation_flag) {
    memory->create(buf3,nbuf_peratom,"msm:buf3");
    memory->create(buf4,nbuf_peratom,"msm:buf4");
  }
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of grid points
------------------------------------------------------------------------- */

void MSM::deallocate()
{
  memory->destroy(buf1);
  memory->destroy(buf2);

  if (differentiation_flag) {
    memory->destroy(buf3);
    memory->destroy(buf4);
  }

  memory->destroy2d_offset(phi1d,-order+1);
  memory->destroy2d_offset(dphi1d,-order+1);

  // deallocate grid levels

  for (int n=0; n<levels; n++) {
    if (qgrid[n])
      memory->destroy3d_offset(qgrid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);

    if (differentiation_flag) {
      if (egrid[n])
        memory->destroy3d_offset(egrid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
    } else {
      if (fxgrid[n])
        memory->destroy3d_offset(fxgrid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
      if (fygrid[n])
        memory->destroy3d_offset(fygrid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
      if (fzgrid[n])
        memory->destroy3d_offset(fzgrid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
    }
  }
}

/* ----------------------------------------------------------------------
   deallocate per-atom memory that depends on # of grid points
------------------------------------------------------------------------- */

void MSM::deallocate_peratom()
{
  if (!differentiation_flag) {
    memory->destroy(buf3);
    memory->destroy(buf4);
  }

  // deallocate grid levels

  for (int n=0; n<levels; n++) {
    if (!differentiation_flag)
      if (egrid[n])
        memory->destroy3d_offset(egrid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);

    if (v0grid[n])
      memory->destroy3d_offset(v0grid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
    if (v1grid[n])
      memory->destroy3d_offset(v1grid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
    if (v2grid[n])
      memory->destroy3d_offset(v2grid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
    if (v3grid[n])
      memory->destroy3d_offset(v3grid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
    if (v4grid[n])
      memory->destroy3d_offset(v4grid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
    if (v5grid[n])
      memory->destroy3d_offset(v5grid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
  }
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of grid levels
------------------------------------------------------------------------- */

void MSM::allocate_levels()
{
  ngrid = new int[levels];

  nx_msm = new int[levels];
  ny_msm = new int[levels];
  nz_msm = new int[levels];

  nxlo_in = new int[levels];
  nylo_in = new int[levels];
  nzlo_in = new int[levels];

  nxhi_in = new int[levels];
  nyhi_in = new int[levels];
  nzhi_in = new int[levels];

  nxlo_in_d = new int[levels];
  nylo_in_d = new int[levels];
  nzlo_in_d = new int[levels];

  nxhi_in_d = new int[levels];
  nyhi_in_d = new int[levels];
  nzhi_in_d = new int[levels];

  nxlo_out = new int[levels];
  nylo_out = new int[levels];
  nzlo_out = new int[levels];

  nxhi_out = new int[levels];
  nyhi_out = new int[levels];
  nzhi_out = new int[levels];

  nxlo_ghost = new int[levels];
  nylo_ghost = new int[levels];
  nzlo_ghost = new int[levels];

  nxhi_ghost = new int[levels];
  nyhi_ghost = new int[levels];
  nzhi_ghost = new int[levels];

  delxinv = new double[levels];
  delyinv = new double[levels];
  delzinv = new double[levels];
  delvolinv = new double[levels];

  qgrid = new double***[levels];
  egrid = new double***[levels];

  fxgrid = new double***[levels];
  fygrid = new double***[levels];
  fzgrid = new double***[levels];

  v0grid = new double***[levels];
  v1grid = new double***[levels];
  v2grid = new double***[levels];
  v3grid = new double***[levels];
  v4grid = new double***[levels];
  v5grid = new double***[levels];

  for (int n=0; n<levels; n++) {
    qgrid[n] = NULL;
    egrid[n] = NULL;

    fxgrid[n] = NULL;
    fygrid[n] = NULL;
    fzgrid[n] = NULL;

    v0grid[n] = NULL;
    v1grid[n] = NULL;
    v2grid[n] = NULL;
    v3grid[n] = NULL;
    v4grid[n] = NULL;
    v5grid[n] = NULL;
  }

}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of grid levels
------------------------------------------------------------------------- */

void MSM::deallocate_levels()
{
  delete [] ngrid;

  delete [] nx_msm;
  delete [] ny_msm;
  delete [] nz_msm;

  delete [] nxlo_in;
  delete [] nylo_in;
  delete [] nzlo_in;

  delete [] nxhi_in;
  delete [] nyhi_in;
  delete [] nzhi_in;

  delete [] nxlo_in_d;
  delete [] nylo_in_d;
  delete [] nzlo_in_d;

  delete [] nxhi_in_d;
  delete [] nyhi_in_d;
  delete [] nzhi_in_d;

  delete [] nxlo_out;
  delete [] nylo_out;
  delete [] nzlo_out;

  delete [] nxhi_out;
  delete [] nyhi_out;
  delete [] nzhi_out;

  delete [] nxlo_ghost;
  delete [] nylo_ghost;
  delete [] nzlo_ghost;

  delete [] nxhi_ghost;
  delete [] nyhi_ghost;
  delete [] nzhi_ghost;

  delete [] delxinv;
  delete [] delyinv;
  delete [] delzinv;
  delete [] delvolinv;

  delete [] qgrid;
  delete [] egrid;

  delete [] fxgrid;
  delete [] fygrid;
  delete [] fzgrid;

  delete [] v0grid;
  delete [] v1grid;
  delete [] v2grid;
  delete [] v3grid;
  delete [] v4grid;
  delete [] v5grid;
}

/* ----------------------------------------------------------------------
   set size of MSM grids
------------------------------------------------------------------------- */

void MSM::set_grid()
{
  if (accuracy_relative <= 0.0)
    error->all(FLERR,"KSpace accuracy must be > 0");

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int nx_max,ny_max,nz_max;

  if (!gridflag) {
    nx_max = ny_max = nz_max = 2;
    double hx = xprd/nx_max;
    double hy = yprd/ny_max;
    double hz = zprd/nz_max;

    double x_error = 2.0*accuracy;
    double y_error = 2.0*accuracy;
    double z_error = 2.0*accuracy;

    while (x_error > accuracy) {
      nx_max *= 2;
      hx = xprd/nx_max;
      x_error = estimate_1d_error(hx,xprd);
    }

    while (y_error > accuracy) {
      ny_max *= 2;
      hy = yprd/ny_max;
      y_error = estimate_1d_error(hy,yprd);
    }

    while (z_error > accuracy) {
      nz_max *= 2;
      hz = zprd/nz_max;
      z_error = estimate_1d_error(hz,zprd);
    }
  } else {
    nx_max = nx_msm_max;
    ny_max = ny_msm_max;
    nz_max = nz_msm_max;
  }

  // boost grid size until it is factorable

  int flag = 0;
  int xlevels,ylevels,zlevels;

  while (!factorable(nx_max,flag,xlevels)) nx_max++;
  while (!factorable(ny_max,flag,ylevels)) ny_max++;
  while (!factorable(nz_max,flag,zlevels)) nz_max++;

  if (flag)
    error->warning(FLERR,"Number of MSM mesh points increased to be a multiple of 2");

  // Find maximum number of levels

  levels = MAX(xlevels,ylevels);
  levels = MAX(levels,zlevels);

  if (levels > MAX_LEVELS)
    error->all(FLERR,"Too many MSM grid levels");

  allocate_levels();

  for (int n = 0; n < levels; n++) {

    if (xlevels-n-1 > 0)
      nx_msm[n] = static_cast<int> (pow(2.0,xlevels-n-1));
    else
      nx_msm[n] = 1;

    if (ylevels-n-1 > 0)
      ny_msm[n] = static_cast<int> (pow(2.0,ylevels-n-1));
    else
      ny_msm[n] = 1;

    if (zlevels-n-1 > 0)
      nz_msm[n] = static_cast<int> (pow(2.0,zlevels-n-1));
    else
      nz_msm[n] = 1;
  }

  // Output grid stats

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  grid = %d %d %d\n",nx_msm[0],ny_msm[0],nz_msm[0]);
      fprintf(screen,"  stencil order = %d\n",order);
    }
    if (logfile) {
      fprintf(logfile,"  grid = %d %d %d\n",nx_msm[0],ny_msm[0],nz_msm[0]);
      fprintf(logfile,"  stencil order = %d\n",order);
    }
  }

}

/* ----------------------------------------------------------------------
   check if all factors of n are in list of factors
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int MSM::factorable(int n, int &flag, int &levels)
{
  int i,norig;
  norig = n;
  levels = 1;

  while (n > 1) {
    for (i = 0; i < nfactors; i++) {
      if (n % factors[i] == 0) {
        n /= factors[i];
        levels++;
        break;
      }
    }
    if (i == nfactors) {
      flag = 1;
      return 0;
    }
  }

  return 1;
}

/* ----------------------------------------------------------------------
   MPI-Reduce so each processor has the full grid
------------------------------------------------------------------------- */
void MSM::grid_swap(int n, double*** &gridn)
{
  double ***gridn_all;
  memory->create3d_offset(gridn_all,nzlo_out[n],nzhi_out[n],nylo_out[n],nyhi_out[n],
                          nxlo_out[n],nxhi_out[n],"msm:grid_all");

  memset(&(gridn_all[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));

  MPI_Allreduce(&(gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),
                &(gridn_all[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),
                ngrid[n],MPI_DOUBLE,MPI_SUM,world);

  // Swap pointers between gridn and gridn_all to avoid need of copy operation

  double ***tmp;
  tmp = gridn;
  gridn = gridn_all;
  gridn_all = tmp;

  memory->destroy3d_offset(gridn_all,nzlo_out[n],nylo_out[n],nxlo_out[n]);
}

/* ----------------------------------------------------------------------
   ghost-swap to accumulate full density in brick decomposition
------------------------------------------------------------------------- */

void MSM::ghost_swap(int n)
{
  double ***qgridn = qgrid[n];

  int i,k,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my ghosts for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my real cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxhi_in[n]+1; ix <= nxhi_out[n]; ix++)
        buf1[k++] = qgridn[iz][iy][ix];

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxlo_in[n]; ix < nxlo_in[n]+nxlo_ghost[n]; ix++)
        qgridn[iz][iy][ix] += buf2[k++];

  // pack my ghosts for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my real cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxlo_out[n]; ix < nxlo_in[n]; ix++)
        buf1[k++] = qgridn[iz][iy][ix];

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxhi_in[n]-nxhi_ghost[n]+1; ix <= nxhi_in[n]; ix++)
        qgridn[iz][iy][ix] += buf2[k++];

  // pack my ghosts for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my real cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nyhi_in[n]+1; iy <= nyhi_out[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++)
        buf1[k++] = qgridn[iz][iy][ix];

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_in[n]; iy < nylo_in[n]+nylo_ghost[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++)
        qgridn[iz][iy][ix] += buf2[k++];

  // pack my ghosts for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my real cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy < nylo_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++)
        buf1[k++] = qgridn[iz][iy][ix];

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nyhi_in[n]-nyhi_ghost[n]+1; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++)
        qgridn[iz][iy][ix] += buf2[k++];

  // pack my ghosts for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my real cells

  k = 0;
  for (iz = nzhi_in[n]+1; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++)
        buf1[k++] = qgridn[iz][iy][ix];

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_in[n]; iz < nzlo_in[n]+nzlo_ghost[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++)
        qgridn[iz][iy][ix] += buf2[k++];

  // pack my ghosts for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my real cells

  k = 0;
  for (iz = nzlo_out[n]; iz < nzlo_in[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++)
        buf1[k++] = qgridn[iz][iy][ix];

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzhi_in[n]-nzhi_ghost[n]+1; iz <= nzhi_in[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++)
        qgridn[iz][iy][ix] += buf2[k++];

}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */

void MSM::particle_map()
{
  int nx,ny,nz;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int flag = 0;

  for (int i = 0; i < nlocal; i++) {

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    nx = static_cast<int> ((x[i][0]-boxlo[0])*delxinv[0]+OFFSET) - OFFSET;
    ny = static_cast<int> ((x[i][1]-boxlo[1])*delyinv[0]+OFFSET) - OFFSET;
    nz = static_cast<int> ((x[i][2]-boxlo[2])*delzinv[0]+OFFSET) - OFFSET;

    part2grid[i][0] = nx;
    part2grid[i][1] = ny;
    part2grid[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlower < nxlo_out[0] || nx+nupper > nxhi_out[0] ||
        ny+nlower < nylo_out[0] || ny+nupper > nyhi_out[0] ||
        nz+nlower < nzlo_out[0] || nz+nupper > nzhi_out[0]) flag = 1;
  }

  if (flag) error->one(FLERR,"Out of range atoms - cannot compute MSM");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void MSM::make_rho()
{
  //fprintf(screen,"MSM aninterpolation\n\n");

  int i,l,m,n,nn,nx,ny,nz,mx,my,mz;
  double dx,dy,dz,x0,y0,z0;

  // clear 3d density array

  double ***qgridn = qgrid[0];

  memset(&(qgridn[nzlo_out[0]][nylo_out[0]][nxlo_out[0]]),0,ngrid[0]*sizeof(double));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx - (x[i][0]-boxlo[0])*delxinv[0];
    dy = ny - (x[i][1]-boxlo[1])*delyinv[0];
    dz = nz - (x[i][2]-boxlo[2])*delzinv[0];

    compute_phis_and_dphis(dx,dy,dz);

    z0 = q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*phi1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*phi1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          qgridn[mz][my][mx] += x0*phi1d[0][l];
        }
      }
    }
  }

}

/* ----------------------------------------------------------------------
   MSM direct part procedure for intermediate grid levels
------------------------------------------------------------------------- */

void MSM::direct_ad(int n)
{
  //fprintf(screen,"Direct contribution on level %i\n\n",n);

  double ***egridn = egrid[n];
  double ***qgridn = qgrid[n];

  // bitmask for PBCs (only works for power of 2 numbers)

  int PBCx,PBCy,PBCz;

  PBCx = nx_msm[n]-1;
  PBCy = ny_msm[n]-1;
  PBCz = nz_msm[n]-1;

  // zero out electric potential

  memset(&(egridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));

  // Simple parallelization of direct sum

  int icx,icy,icz,ix,iy,iz,k;
  int jj,kk;
  double qtmp;
  double esum,v0sum,v1sum,v2sum,v3sum,v4sum,v5sum;

  for (icz = nzlo_in_d[n]; icz <= nzhi_in_d[n]; icz++) {
    for (icy = nylo_in_d[n]; icy <= nyhi_in_d[n]; icy++) {
      for (icx = nxlo_in_d[n]; icx <= nxhi_in_d[n]; icx++) {
        if (evflag) {
          esum = 0.0;
          v0sum = v1sum = v2sum = 0.0;
          v3sum = v4sum = v5sum = 0.0;
        }

  // do double loop over points on the intermediate grid level
  // for now, assume I own all points on the intermediate grid level

        k = 0;
        for (iz = nzlo_direct; iz <= nzhi_direct; iz++) {
          kk = (icz+iz)&PBCz;
          for (iy = nylo_direct; iy <= nyhi_direct; iy++) {
            jj = (icy+iy)&PBCy;
            for (ix = nxlo_direct; ix <= nxhi_direct; ix++) {
              qtmp = qgridn[kk][jj][(icx+ix)&PBCx];
              egridn[icz][icy][icx] += g_direct[n][k] * qtmp;

              if (evflag) {
                if (eflag_global) esum += g_direct[n][k] * qtmp;
                if (vflag_global) {
                  v0sum += v0_direct[n][k] * qtmp;
                  v1sum += v1_direct[n][k] * qtmp;
                  v2sum += v2_direct[n][k] * qtmp;
                  v3sum += v3_direct[n][k] * qtmp;
                  v4sum += v4_direct[n][k] * qtmp;
                  v5sum += v5_direct[n][k] * qtmp;
                }
              }
              k++;
            }
          }
        }

        if (evflag) {
          qtmp = qgridn[icz][icy][icx];
          if (eflag_global) energy += esum * qtmp;
          if (vflag_global) {
            virial[0] += v0sum * qtmp;
            virial[1] += v1sum * qtmp;
            virial[2] += v2sum * qtmp;
            virial[3] += v3sum * qtmp;
            virial[4] += v4sum * qtmp;
            virial[5] += v5sum * qtmp;
          }
        }

      }
    }
  }
}

/* ----------------------------------------------------------------------
   MSM direct part procedure for intermediate grid levels
------------------------------------------------------------------------- */

void MSM::direct(int n)
{
  double ***qgridn = qgrid[n];

  double ***fxgridn = fxgrid[n];
  double ***fygridn = fygrid[n];
  double ***fzgridn = fzgrid[n];

  // bitmask for PBCs (only works for power of 2 numbers)

  int PBCx,PBCy,PBCz;

  PBCx = nx_msm[n]-1;
  PBCy = ny_msm[n]-1;
  PBCz = nz_msm[n]-1;

  int icx,icy,icz,ix,iy,iz,k;

  // zero out forces

  memset(&(fxgridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
  memset(&(fygridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
  memset(&(fzgridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));

  // Simple parallelization of direct sum

  int jj,kk;
  double qtmp;
  double esum,v0sum,v1sum,v2sum,v3sum,v4sum,v5sum;

  for (icz = nzlo_in_d[n]; icz <= nzhi_in_d[n]; icz++) {
    for (icy = nylo_in_d[n]; icy <= nyhi_in_d[n]; icy++) {
      for (icx = nxlo_in_d[n]; icx <= nxhi_in_d[n]; icx++) {
        if (evflag) {
          esum = 0.0;
          v0sum = v1sum = v2sum = 0.0;
          v3sum = v4sum = v5sum = 0.0;
        }

  // do double loop over points on the intermediate grid level
  // for now, assume I own all points on the intermediate grid level

        k = 0;
        for (iz = nzlo_direct; iz <= nzhi_direct; iz++) {
          kk = (icz+iz)&PBCz;
          for (iy = nylo_direct; iy <= nyhi_direct; iy++) {
            jj = (icy+iy)&PBCy;
            for (ix = nxlo_direct; ix <= nxhi_direct; ix++) {
              qtmp = qgridn[kk][jj][(icx+ix)&PBCx];
              fxgridn[icz][icy][icx] += dgx_direct[n][k] * qtmp;
              fygridn[icz][icy][icx] += dgy_direct[n][k] * qtmp;
              fzgridn[icz][icy][icx] += dgz_direct[n][k] * qtmp;

              if (evflag) {
                if (eflag_global) esum += g_direct[n][k] * qtmp;
                if (vflag_global) {
                  v0sum += v0_direct[n][k] * qtmp;
                  v1sum += v1_direct[n][k] * qtmp;
                  v2sum += v2_direct[n][k] * qtmp;
                  v3sum += v3_direct[n][k] * qtmp;
                  v4sum += v4_direct[n][k] * qtmp;
                  v5sum += v5_direct[n][k] * qtmp;
                }
              }

              k++;
            }
          }
        }

        if (evflag) {
          qtmp = qgridn[icz][icy][icx];
          if (eflag_global) energy += esum * qtmp;
          if (vflag_global) {
            virial[0] += v0sum * qtmp;
            virial[1] += v1sum * qtmp;
            virial[2] += v2sum * qtmp;
            virial[3] += v3sum * qtmp;
            virial[4] += v4sum * qtmp;
            virial[5] += v5sum * qtmp;
          }
        }

      }
    }
  }
}

/* ----------------------------------------------------------------------
   MSM direct part procedure for intermediate grid levels
------------------------------------------------------------------------- */

void MSM::direct_peratom(int n)
{
  double ***egridn = egrid[n];
  double ***qgridn = qgrid[n];

  double ***v0gridn = v0grid[n];
  double ***v1gridn = v1grid[n];
  double ***v2gridn = v2grid[n];
  double ***v3gridn = v3grid[n];
  double ***v4gridn = v4grid[n];
  double ***v5gridn = v5grid[n];

  // bitmask for PBCs (only works for power of 2 numbers)

  int PBCx,PBCy,PBCz;

  PBCx = nx_msm[n]-1;
  PBCy = ny_msm[n]-1;
  PBCz = nz_msm[n]-1;

  int icx,icy,icz,ix,iy,iz,k;

  // zero out electric potential

  if (!differentiation_flag && eflag_atom)
    memset(&(egridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));

  // zero out virial

  if (vflag_atom) {
    memset(&(v0gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
    memset(&(v1gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
    memset(&(v2gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
    memset(&(v3gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
    memset(&(v4gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
    memset(&(v5gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
  }

  // Simple parallelization of direct sum

  int jj,kk;
  double qtmp;

  int ndiff_eflag_atom = !differentiation_flag && eflag_atom;

  for (icz = nzlo_in_d[n]; icz <= nzhi_in_d[n]; icz++) {
    for (icy = nylo_in_d[n]; icy <= nyhi_in_d[n]; icy++) {
      for (icx = nxlo_in_d[n]; icx <= nxhi_in_d[n]; icx++) {

  // do double loop over points on the intermediate grid level
  // for now, assume I own all points on the intermediate grid level

        k = 0;
        for (iz = nzlo_direct; iz <= nzhi_direct; iz++) {
          kk = (icz+iz)&PBCz;
          for (iy = nylo_direct; iy <= nyhi_direct; iy++) {
            jj = (icy+iy)&PBCy;
            for (ix = nxlo_direct; ix <= nxhi_direct; ix++) {
              qtmp = qgridn[kk][jj][(icx+ix)&PBCx];

              if (ndiff_eflag_atom)
                egridn[icz][icy][icx] += g_direct[n][k] * qtmp;

              if (vflag_atom) {
                v0gridn[icz][icy][icx] += v0_direct[n][k] * qtmp;
                v1gridn[icz][icy][icx] += v1_direct[n][k] * qtmp;
                v2gridn[icz][icy][icx] += v2_direct[n][k] * qtmp;
                v3gridn[icz][icy][icx] += v3_direct[n][k] * qtmp;
                v4gridn[icz][icy][icx] += v4_direct[n][k] * qtmp;
                v5gridn[icz][icy][icx] += v5_direct[n][k] * qtmp;
              }

              k++;
            }
          }
        }

      }
    }
  }
}

/* ----------------------------------------------------------------------
   MSM restriction  procedure for intermediate grid levels
------------------------------------------------------------------------- */

void MSM::restriction(int n)
{
  //fprintf(screen,"Restricting from level %i to %i\n\n",n,n+1);

  int p = order-1;

  double ***qgrid1 = qgrid[n];
  double ***qgrid2 = qgrid[n+1];

  // bitmask for PBCs (only works for power of 2 numbers)

  int PBCx,PBCy,PBCz;

  PBCx = nx_msm[n]-1;
  PBCy = ny_msm[n]-1;
  PBCz = nz_msm[n]-1;

  //restrict grid (going from grid n to grid n+1, i.e. to a coarser grid)

  for (int nu=-p; nu<=p; nu++) {
    phi1d[0][nu] = compute_phi(nu*delxinv[n+1]/delxinv[n]);
    phi1d[1][nu] = compute_phi(nu*delyinv[n+1]/delyinv[n]);
    phi1d[2][nu] = compute_phi(nu*delzinv[n+1]/delzinv[n]);
  }

  int ip,jp,kp,ic,jc,kc,i,j,k;
  int jj,kk;
  double phiz,phizy;

  // zero out charge on coarser grid

  memset(&(qgrid2[nzlo_out[n+1]][nylo_out[n+1]][nxlo_out[n+1]]),0,ngrid[n+1]*sizeof(double));

  for (kp = nzlo_in[n+1]; kp <= nzhi_in[n+1]; kp++)
    for (jp = nylo_in[n+1]; jp <= nyhi_in[n+1]; jp++)
      for (ip = nxlo_in[n+1]; ip <= nxhi_in[n+1]; ip++) {

        ic = static_cast<int> (ip*delxinv[n]/delxinv[n+1]);
        jc = static_cast<int> (jp*delyinv[n]/delyinv[n+1]);
        kc = static_cast<int> (kp*delzinv[n]/delzinv[n+1]);
        for (k=-p; k<=p; k++) { // Could make this faster by eliminating zeros
          kk = (kc+k)&PBCz;
          phiz = phi1d[2][k];
          for (j=-p; j<=p; j++) {
            jj = (jc+j)&PBCy;
            phizy = phi1d[1][j]*phiz;
            for (i=-p; i<=p; i++) {
              qgrid2[kp][jp][ip] += qgrid1[kk][jj][(ic+i)&PBCx] *
                phi1d[0][i]*phizy;
            }
          }
        }
      }

}

/* ----------------------------------------------------------------------
   MSM prolongation procedure for intermediate grid levels
------------------------------------------------------------------------- */

void MSM::prolongation(int n)
{
  //fprintf(screen,"Prolongating from level %i to %i\n\n",n+1,n);

  int p = order-1;

  double ***egrid1 = egrid[n];
  double ***egrid2 = egrid[n+1];

  double ***fxgrid1 = fxgrid[n];
  double ***fxgrid2 = fxgrid[n+1];
  double ***fygrid1 = fygrid[n];
  double ***fygrid2 = fygrid[n+1];
  double ***fzgrid1 = fzgrid[n];
  double ***fzgrid2 = fzgrid[n+1];

  double ***v0grid1 = v0grid[n];
  double ***v0grid2 = v0grid[n+1];
  double ***v1grid1 = v1grid[n];
  double ***v1grid2 = v1grid[n+1];
  double ***v2grid1 = v2grid[n];
  double ***v2grid2 = v2grid[n+1];
  double ***v3grid1 = v3grid[n];
  double ***v3grid2 = v3grid[n+1];
  double ***v4grid1 = v4grid[n];
  double ***v4grid2 = v4grid[n+1];
  double ***v5grid1 = v5grid[n];
  double ***v5grid2 = v5grid[n+1];

  // bitmask for PBCs (only works for power of 2 numbers)

  int PBCx,PBCy,PBCz;

  PBCx = nx_msm[n]-1;
  PBCy = ny_msm[n]-1;
  PBCz = nz_msm[n]-1;

  //prolongate grid (going from grid n to grid n-1, i.e. to a finer grid)

  for (int nu=-p; nu<=p; nu++) {
    phi1d[0][nu] = compute_phi(nu*delxinv[n+1]/delxinv[n]);
    phi1d[1][nu] = compute_phi(nu*delyinv[n+1]/delyinv[n]);
    phi1d[2][nu] = compute_phi(nu*delzinv[n+1]/delzinv[n]);
  }

  int ip,jp,kp,ic,jc,kc,i,j,k;
  int jj,kk,ii;
  double phiz,phizy,phi3d;

  for (kp = nzlo_in[n+1]; kp <= nzhi_in[n+1]; kp++)
    for (jp = nylo_in[n+1]; jp <= nyhi_in[n+1]; jp++)
      for (ip = nxlo_in[n+1]; ip <= nxhi_in[n+1]; ip++) {

        ic = static_cast<int> (ip*delxinv[n]/delxinv[n+1]);
        jc = static_cast<int> (jp*delyinv[n]/delyinv[n+1]);
        kc = static_cast<int> (kp*delzinv[n]/delzinv[n+1]);
        for (k=-p; k<=p; k++) { // Could make this faster by eliminating zeros or creating separate functions
          kk = (kc+k)&PBCz;
          phiz = phi1d[2][k];
          for (j=-p; j<=p; j++) {
            jj = (jc+j)&PBCy;
            phizy = phi1d[1][j]*phiz;
            for (i=-p; i<=p; i++) {
              ii = (ic+i)&PBCx;
              phi3d = phi1d[0][i]*phizy;

              if (differentiation_flag || eflag_atom) {
                egrid1[kk][jj][ii] += egrid2[kp][jp][ip] * phi3d;
              }

              if (!differentiation_flag) {
                fxgrid1[kk][jj][ii] += fxgrid2[kp][jp][ip] * phi3d;
                fygrid1[kk][jj][ii] += fygrid2[kp][jp][ip] * phi3d;
                fzgrid1[kk][jj][ii] += fzgrid2[kp][jp][ip] * phi3d;
              }

              if (vflag_atom) {
                v0grid1[kk][jj][ii] += v0grid2[kp][jp][ip] * phi3d;
                v1grid1[kk][jj][ii] += v1grid2[kp][jp][ip] * phi3d;
                v2grid1[kk][jj][ii] += v2grid2[kp][jp][ip] * phi3d;
                v3grid1[kk][jj][ii] += v3grid2[kp][jp][ip] * phi3d;
                v4grid1[kk][jj][ii] += v4grid2[kp][jp][ip] * phi3d;
                v5grid1[kk][jj][ii] += v5grid2[kp][jp][ip] * phi3d;
              }

            }
          }
        }
      }

}

/* ----------------------------------------------------------------------
   ghost-swap to fill ghost cells of my brick with field values
------------------------------------------------------------------------- */

void MSM::fillbrick_ad_peratom(int n)
{
  double ***egridn = egrid[n];

  double ***v0gridn = v0grid[n];
  double ***v1gridn = v1grid[n];
  double ***v2gridn = v2grid[n];
  double ***v3gridn = v3grid[n];
  double ***v4gridn = v4grid[n];
  double ***v5gridn = v5grid[n];

  int i,k,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my real cells for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzhi_in[n]-nzhi_ghost[n]+1; iz <= nzhi_in[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          buf3[k++] = egridn[iz][iy][ix];
        if (vflag_atom) {
          buf3[k++] = v0gridn[iz][iy][ix];
          buf3[k++] = v1gridn[iz][iy][ix];
          buf3[k++] = v2gridn[iz][iy][ix];
          buf3[k++] = v3gridn[iz][iy][ix];
          buf3[k++] = v4gridn[iz][iy][ix];
          buf3[k++] = v5gridn[iz][iy][ix];
        }
      }

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < k; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_DOUBLE,comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf3,k,MPI_DOUBLE,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz < nzlo_in[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          egridn[iz][iy][ix] = buf4[k++];
        if (vflag_atom) {
          v0gridn[iz][iy][ix] = buf4[k++];
          v1gridn[iz][iy][ix] = buf4[k++];
          v2gridn[iz][iy][ix] = buf4[k++];
          v3gridn[iz][iy][ix] = buf4[k++];
          v4gridn[iz][iy][ix] = buf4[k++];
          v5gridn[iz][iy][ix] = buf4[k++];
        }
      }

  // pack my real cells for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_in[n]; iz < nzlo_in[n]+nzlo_ghost[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          buf3[k++] = egridn[iz][iy][ix];
        if (vflag_atom) {
          buf3[k++] = v0gridn[iz][iy][ix];
          buf3[k++] = v1gridn[iz][iy][ix];
          buf3[k++] = v2gridn[iz][iy][ix];
          buf3[k++] = v3gridn[iz][iy][ix];
          buf3[k++] = v4gridn[iz][iy][ix];
          buf3[k++] = v5gridn[iz][iy][ix];
        }
      }

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < k; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_DOUBLE,comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf3,k,MPI_DOUBLE,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzhi_in[n]+1; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          egridn[iz][iy][ix] = buf4[k++];
        if (vflag_atom) {
          v0gridn[iz][iy][ix] = buf4[k++];
          v1gridn[iz][iy][ix] = buf4[k++];
          v2gridn[iz][iy][ix] = buf4[k++];
          v3gridn[iz][iy][ix] = buf4[k++];
          v4gridn[iz][iy][ix] = buf4[k++];
          v5gridn[iz][iy][ix] = buf4[k++];
        }
      }

  // pack my real cells for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nyhi_in[n]-nyhi_ghost[n]+1; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          buf3[k++] = egridn[iz][iy][ix];
        if (vflag_atom) {
          buf3[k++] = v0gridn[iz][iy][ix];
          buf3[k++] = v1gridn[iz][iy][ix];
          buf3[k++] = v2gridn[iz][iy][ix];
          buf3[k++] = v3gridn[iz][iy][ix];
          buf3[k++] = v4gridn[iz][iy][ix];
          buf3[k++] = v5gridn[iz][iy][ix];
        }
      }

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < k; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_DOUBLE,comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf3,k,MPI_DOUBLE,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy < nylo_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          egridn[iz][iy][ix] = buf4[k++];
        if (vflag_atom) {
          v0gridn[iz][iy][ix] = buf4[k++];
          v1gridn[iz][iy][ix] = buf4[k++];
          v2gridn[iz][iy][ix] = buf4[k++];
          v3gridn[iz][iy][ix] = buf4[k++];
          v4gridn[iz][iy][ix] = buf4[k++];
          v5gridn[iz][iy][ix] = buf4[k++];
        }
      }

  // pack my real cells for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_in[n]; iy < nylo_in[n]+nylo_ghost[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          buf3[k++] = egridn[iz][iy][ix];
        if (vflag_atom) {
          buf3[k++] = v0gridn[iz][iy][ix];
          buf3[k++] = v1gridn[iz][iy][ix];
          buf3[k++] = v2gridn[iz][iy][ix];
          buf3[k++] = v3gridn[iz][iy][ix];
          buf3[k++] = v4gridn[iz][iy][ix];
          buf3[k++] = v5gridn[iz][iy][ix];
        }
      }

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < k; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_DOUBLE,comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf3,k,MPI_DOUBLE,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nyhi_in[n]+1; iy <= nyhi_out[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          egridn[iz][iy][ix] = buf4[k++];
        if (vflag_atom) {
          v0gridn[iz][iy][ix] = buf4[k++];
          v1gridn[iz][iy][ix] = buf4[k++];
          v2gridn[iz][iy][ix] = buf4[k++];
          v3gridn[iz][iy][ix] = buf4[k++];
          v4gridn[iz][iy][ix] = buf4[k++];
          v5gridn[iz][iy][ix] = buf4[k++];
        }
      }

  // pack my real cells for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxhi_in[n]-nxhi_ghost[n]+1; ix <= nxhi_in[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          buf3[k++] = egridn[iz][iy][ix];
        if (vflag_atom) {
          buf3[k++] = v0gridn[iz][iy][ix];
          buf3[k++] = v1gridn[iz][iy][ix];
          buf3[k++] = v2gridn[iz][iy][ix];
          buf3[k++] = v3gridn[iz][iy][ix];
          buf3[k++] = v4gridn[iz][iy][ix];
          buf3[k++] = v5gridn[iz][iy][ix];
        }
      }

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < k; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_DOUBLE,comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf3,k,MPI_DOUBLE,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxlo_out[n]; ix < nxlo_in[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          egridn[iz][iy][ix] = buf4[k++];
        if (vflag_atom) {
          v0gridn[iz][iy][ix] = buf4[k++];
          v1gridn[iz][iy][ix] = buf4[k++];
          v2gridn[iz][iy][ix] = buf4[k++];
          v3gridn[iz][iy][ix] = buf4[k++];
          v4gridn[iz][iy][ix] = buf4[k++];
          v5gridn[iz][iy][ix] = buf4[k++];
        }
      }

  // pack my real cells for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxlo_in[n]; ix < nxlo_in[n]+nxlo_ghost[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          buf3[k++] = egridn[iz][iy][ix];
        if (vflag_atom) {
          buf3[k++] = v0gridn[iz][iy][ix];
          buf3[k++] = v1gridn[iz][iy][ix];
          buf3[k++] = v2gridn[iz][iy][ix];
          buf3[k++] = v3gridn[iz][iy][ix];
          buf3[k++] = v4gridn[iz][iy][ix];
          buf3[k++] = v5gridn[iz][iy][ix];
        }
      }

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < k; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_DOUBLE,comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf3,k,MPI_DOUBLE,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxhi_in[n]+1; ix <= nxhi_out[n]; ix++) {
        if (differentiation_flag || eflag_atom)
          egridn[iz][iy][ix] = buf4[k++];
        if (vflag_atom) {
          v0gridn[iz][iy][ix] = buf4[k++];
          v1gridn[iz][iy][ix] = buf4[k++];
          v2gridn[iz][iy][ix] = buf4[k++];
          v3gridn[iz][iy][ix] = buf4[k++];
          v4gridn[iz][iy][ix] = buf4[k++];
          v5gridn[iz][iy][ix] = buf4[k++];
        }
      }
}

/* ----------------------------------------------------------------------
   ghost-swap to fill ghost cells of my brick with field values
------------------------------------------------------------------------- */

void MSM::fillbrick(int n)
{
  double ***fxgridn = fxgrid[n];
  double ***fygridn = fygrid[n];
  double ***fzgridn = fzgrid[n];

  int i,k,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my real cells for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzhi_in[n]-nzhi_ghost[n]+1; iz <= nzhi_in[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        buf1[k++] = fxgridn[iz][iy][ix];
        buf1[k++] = fygridn[iz][iy][ix];
        buf1[k++] = fzgridn[iz][iy][ix];
      }

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz < nzlo_in[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        fxgridn[iz][iy][ix] = buf2[k++];
        fygridn[iz][iy][ix] = buf2[k++];
        fzgridn[iz][iy][ix] = buf2[k++];
      }

  // pack my real cells for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_in[n]; iz < nzlo_in[n]+nzlo_ghost[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        buf1[k++] = fxgridn[iz][iy][ix];
        buf1[k++] = fygridn[iz][iy][ix];
        buf1[k++] = fzgridn[iz][iy][ix];
      }

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzhi_in[n]+1; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        fxgridn[iz][iy][ix] = buf2[k++];
        fygridn[iz][iy][ix] = buf2[k++];
        fzgridn[iz][iy][ix] = buf2[k++];
      }

  // pack my real cells for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nyhi_in[n]-nyhi_ghost[n]+1; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        buf1[k++] = fxgridn[iz][iy][ix];
        buf1[k++] = fygridn[iz][iy][ix];
        buf1[k++] = fzgridn[iz][iy][ix];
      }

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy < nylo_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        fxgridn[iz][iy][ix] = buf2[k++];
        fygridn[iz][iy][ix] = buf2[k++];
        fzgridn[iz][iy][ix] = buf2[k++];
      }

  // pack my real cells for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_in[n]; iy < nylo_in[n]+nylo_ghost[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        buf1[k++] = fxgridn[iz][iy][ix];
        buf1[k++] = fygridn[iz][iy][ix];
        buf1[k++] = fzgridn[iz][iy][ix];
      }

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nyhi_in[n]+1; iy <= nyhi_out[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        fxgridn[iz][iy][ix] = buf2[k++];
        fygridn[iz][iy][ix] = buf2[k++];
        fzgridn[iz][iy][ix] = buf2[k++];
      }

  // pack my real cells for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxhi_in[n]-nxhi_ghost[n]+1; ix <= nxhi_in[n]; ix++) {
        buf1[k++] = fxgridn[iz][iy][ix];
        buf1[k++] = fygridn[iz][iy][ix];
        buf1[k++] = fzgridn[iz][iy][ix];
      }

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxlo_out[n]; ix < nxlo_in[n]; ix++) {
        fxgridn[iz][iy][ix] = buf2[k++];
        fygridn[iz][iy][ix] = buf2[k++];
        fzgridn[iz][iy][ix] = buf2[k++];
      }

  // pack my real cells for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxlo_in[n]; ix < nxlo_in[n]+nxlo_ghost[n]; ix++) {
        buf1[k++] = fxgridn[iz][iy][ix];
        buf1[k++] = fygridn[iz][iy][ix];
        buf1[k++] = fzgridn[iz][iy][ix];
      }

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < k; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf1,k,MPI_DOUBLE,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxhi_in[n]+1; ix <= nxhi_out[n]; ix++) {
        fxgridn[iz][iy][ix] = buf2[k++];
        fygridn[iz][iy][ix] = buf2[k++];
        fzgridn[iz][iy][ix] = buf2[k++];
      }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get force on my particles
------------------------------------------------------------------------- */

void MSM::fieldforce_ad()
{
  //fprintf(screen,"MSM interpolation\n\n");

  double ***egridn = egrid[0];
  double ***qgridn = qgrid[0];

  int i,l,m,n,nx,ny,nz,mx,my,mz;
  double dx,dy,dz;
  double phi_x,phi_y,phi_z;
  double dphi_x,dphi_y,dphi_z;
  double ekx,eky,ekz;


  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx - (x[i][0]-boxlo[0])*delxinv[0];
    dy = ny - (x[i][1]-boxlo[1])*delyinv[0];
    dz = nz - (x[i][2]-boxlo[2])*delzinv[0];

    compute_phis_and_dphis(dx,dy,dz);

    ekx = eky = ekz = 0.0;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      phi_z = phi1d[2][n];
      dphi_z = dphi1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        phi_y = phi1d[1][m];
        dphi_y = dphi1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          phi_x = phi1d[0][l];
          dphi_x = dphi1d[0][l];
          ekx += dphi_x*phi_y*phi_z*egridn[mz][my][mx];
          eky += phi_x*dphi_y*phi_z*egridn[mz][my][mx];
          ekz += phi_x*phi_y*dphi_z*egridn[mz][my][mx];
        }
      }
    }

    ekx *= delxinv[0];
    eky *= delyinv[0];
    ekz *= delzinv[0];

    // convert E-field to force

    const double qfactor = force->qqrd2e*scale*q[i];
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    f[i][2] += qfactor*ekz;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get my particles
------------------------------------------------------------------------- */

void MSM::fieldforce()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  double dx,dy,dz,x0,y0,z0;
  double ekx,eky,ekz;

  double ***fxgridn = fxgrid[0];
  double ***fygridn = fygrid[0];
  double ***fzgridn = fzgrid[0];

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx - (x[i][0]-boxlo[0])*delxinv[0];
    dy = ny - (x[i][1]-boxlo[1])*delyinv[0];
    dz = nz - (x[i][2]-boxlo[2])*delzinv[0];

    compute_phis_and_dphis(dx,dy,dz);

    ekx = eky = ekz = 0.0;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = phi1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*phi1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*phi1d[0][l];
          ekx -= x0*fxgridn[mz][my][mx];
          eky -= x0*fygridn[mz][my][mx];
          ekz -= x0*fzgridn[mz][my][mx];
        }
      }
    }

    // convert E-field to force

    const double qfactor = force->qqrd2e * scale * q[i];
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    f[i][2] += qfactor*ekz;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get per-atom energy/virial
------------------------------------------------------------------------- */

void MSM::fieldforce_peratom()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  double dx,dy,dz,x0,y0,z0;
  double u,v0,v1,v2,v3,v4,v5;

  double ***egridn = egrid[0];

  double ***v0gridn = v0grid[0];
  double ***v1gridn = v1grid[0];
  double ***v2gridn = v2grid[0];
  double ***v3gridn = v3grid[0];
  double ***v4gridn = v4grid[0];
  double ***v5gridn = v5grid[0];

  // loop over my charges, interpolate from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx - (x[i][0]-boxlo[0])*delxinv[0];
    dy = ny - (x[i][1]-boxlo[1])*delyinv[0];
    dz = nz - (x[i][2]-boxlo[2])*delzinv[0];

    compute_phis_and_dphis(dx,dy,dz);

    u = v0 = v1 = v2 = v3 = v4 = v5 = 0.0;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = phi1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*phi1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*phi1d[0][l];
          if (eflag_atom) u += x0*egridn[mz][my][mx];
          if (vflag_atom) {
            v0 += x0*v0gridn[mz][my][mx];
            v1 += x0*v1gridn[mz][my][mx];
            v2 += x0*v2gridn[mz][my][mx];
            v3 += x0*v3gridn[mz][my][mx];
            v4 += x0*v4gridn[mz][my][mx];
            v5 += x0*v5gridn[mz][my][mx];
          }
        }
      }
    }

    if (eflag_atom) eatom[i] += q[i]*u;
    if (vflag_atom) {
      vatom[i][0] += q[i]*v0;
      vatom[i][1] += q[i]*v1;
      vatom[i][2] += q[i]*v2;
      vatom[i][3] += q[i]*v3;
      vatom[i][4] += q[i]*v4;
      vatom[i][5] += q[i]*v5;
    }
  }
}

/* ----------------------------------------------------------------------
   charge assignment into phi1d
------------------------------------------------------------------------- */

void MSM::compute_phis_and_dphis(const double &dx, const double &dy, 
                                 const double &dz)
{
  double delx,dely,delz;

  for (int nu = nlower; nu <= nupper; nu++) {
    delx = dx + double(nu);
    dely = dy + double(nu);
    delz = dz + double(nu);

    phi1d[0][nu] = compute_phi(delx);
    phi1d[1][nu] = compute_phi(dely);
    phi1d[2][nu] = compute_phi(delz);
    dphi1d[0][nu] = compute_dphi(delx);
    dphi1d[1][nu] = compute_dphi(dely);
    dphi1d[2][nu] = compute_dphi(delz);
  }
}

/* ----------------------------------------------------------------------
   compute phi using interpolating polynomial
   see Eq 7 from Parallel Computing 35 (2009) 164–177
   and Hardy's thesis
------------------------------------------------------------------------- */

double MSM::compute_phi(const double &xi)
{
  double phi;
  double abs_xi = fabs(xi);
  double xi2 = xi*xi;

  if (order == 4) {
    if (abs_xi <= 1) {
      phi = (1.0 - abs_xi)*(1.0 + abs_xi - 1.5*xi2);
    } else if (abs_xi <= 2) {
      phi = -0.5*(abs_xi - 1.0)*(2.0 - abs_xi)*(2.0 - abs_xi);
    } else {
      phi = 0.0;
    }

  } else if (order == 6) {
    if (abs_xi <= 1) {
      phi = (1.0 - xi2)*(2.0 - abs_xi)*(6.0 + 3.0*abs_xi -
        5.0*xi2)/12.0;
    } else if (abs_xi <= 2) {
      phi = -(abs_xi - 1.0)*(2.0 - abs_xi)*(3.0 - abs_xi)*
        (4.0 + 9.0*abs_xi - 5.0*xi2)/24.0;
    } else if (abs_xi <= 3) {
      phi = (abs_xi - 1.0)*(abs_xi - 2.0)*(3.0 - abs_xi)*
        (3.0 - abs_xi)*(4.0 - abs_xi)/24.0;
    } else {
      phi = 0.0;
    }

  } else if (order == 8) {
    if (abs_xi <= 1) {
      phi = (1.0 - xi2)*(4.0 - xi2)*(3.0 - abs_xi)*
        (12.0 + 4.0*abs_xi - 7.0*xi2)/144.0;
    } else if (abs_xi <= 2) {
      phi = -(xi2 - 1.0)*(2.0 - abs_xi)*(3.0 - abs_xi)*
        (4.0 - abs_xi)*(10.0 + 12.0*abs_xi - 7.0*xi2)/240.0;
    } else if (abs_xi <= 3) {
      phi = (abs_xi - 1.0)*(abs_xi - 2.0)*(3.0 - abs_xi)*(4.0 - abs_xi)*
        (5.0 - abs_xi)*(6.0 + 20.0*abs_xi - 7.0*xi2)/720.0;
    } else if (abs_xi <= 4) {
      phi = -(abs_xi - 1.0)*(abs_xi - 2.0)*(abs_xi - 3.0)*(4.0 - abs_xi)*
        (4.0 - abs_xi)*(5.0 - abs_xi)*(6.0 - abs_xi)/720.0;
    } else {
      phi = 0.0;
    }

  } else if (order == 10) {
    if (abs_xi <= 1) {
      phi = (1.0 - xi2)*(4.0 - xi2)*(9.0 - xi2)*
        (4.0 - abs_xi)*(20.0 + 5.0*abs_xi - 9.0*xi2)/2880.0;
    } else if (abs_xi <= 2) {
      phi = -(xi2 - 1.0)*(4.0 - xi2)*(3.0 - abs_xi)*(4.0 - abs_xi)*
        (5.0 - abs_xi)*(6.0 + 5.0*abs_xi - 3.0*xi2)/1440.0;
    } else if (abs_xi <= 3) {
      phi = (xi2 - 1.0)*(abs_xi - 2.0)*(3.0 - abs_xi)*(4.0 - abs_xi)*
        (5.0 - abs_xi)*(6.0 - abs_xi)*(14.0 + 25.0*abs_xi - 9.0*xi2)/10080.0;
    } else if (abs_xi <= 4) {
      phi = -(abs_xi - 1.0)*(abs_xi - 2.0)*(abs_xi - 3.0)*(4.0 - abs_xi)*
        (5.0 - abs_xi)*(6.0 - abs_xi)*(7.0 - abs_xi)*
        (8.0 + 35.0*abs_xi - 9.0*xi2)/40320.0;
    } else if (abs_xi <= 5) {
      phi = (abs_xi - 1.0)*(abs_xi - 2.0)*(abs_xi - 3.0)*
        (abs_xi - 4.0)*(5.0 - abs_xi)*(5.0 - abs_xi)*(6.0 - abs_xi)*
        (7.0 - abs_xi)*(8.0 - abs_xi)/40320.0;
    } else {
      phi = 0.0;
    }
  }

  return phi;
}

/* ----------------------------------------------------------------------
   compute the derivative of phi
   phi is an interpolating polynomial
   see Eq 7 from Parallel Computing 35 (2009) 164–177
   and Hardy's thesis
------------------------------------------------------------------------- */

double MSM::compute_dphi(const double &xi)
{
  double dphi;
  double abs_xi = fabs(xi);

  if (order == 4) {
    double xi2 = xi*xi;
    double abs_xi2 = abs_xi*abs_xi;
    if (abs_xi == 0.0) {
      dphi = 0.0;
    } else if (abs_xi <= 1) {
      dphi = xi*(3*xi2 + 6*abs_xi2 - 10*abs_xi)/2.0/abs_xi;
    } else if (abs_xi <= 2) {
      dphi = xi*(2 - abs_xi)*(3*abs_xi - 4)/2.0/abs_xi;
    } else {
      dphi = 0.0;
    }

  } else if (order == 6) {
    double xi2 = xi*xi;
    double xi4 = xi2*xi2;
    double abs_xi2 = abs_xi*abs_xi;
    double abs_xi3 = abs_xi2*abs_xi;
    double abs_xi4 = abs_xi2*abs_xi2;
    if (abs_xi == 0.0) {
      dphi = 0.0;
    } else if (abs_xi <= 1) {
      dphi = xi*(46*xi2*abs_xi - 20*xi2*abs_xi2 - 5*xi4 + 5*xi2 +
        6*abs_xi3 + 10*abs_xi2 - 50*abs_xi)/12.0/abs_xi;
    } else if (abs_xi <= 2) {
      dphi = xi*(15*xi2*abs_xi2 - 60*xi2*abs_xi + 55*xi2 +
        10*abs_xi4 - 96*abs_xi3 + 260*abs_xi2 - 210*abs_xi + 10)/
        24.0/abs_xi;
    } else if (abs_xi <= 3) {
      dphi = -xi*(abs_xi - 3)*(5*abs_xi3 - 37*abs_xi2 +
      84*abs_xi - 58)/24.0/abs_xi;
    } else {
      dphi = 0.0;
    }

  } else if (order == 8) {
    double xi2 = xi*xi;
    double xi4 = xi2*xi2;
    double xi6 = xi4*xi2;
    double abs_xi2 = abs_xi*abs_xi;
    double abs_xi3 = abs_xi2*abs_xi;
    double abs_xi4 = abs_xi2*abs_xi2;
    double abs_xi5 = abs_xi4*abs_xi;
    double abs_xi6 = abs_xi5*abs_xi;
    if (abs_xi == 0.0) {
      dphi = 0.0;
    } else if (abs_xi <= 1) {
      dphi = xi*(7*xi6 + 42*xi4*abs_xi2 - 134*xi4*abs_xi - 35*xi4 -
        16*xi2*abs_xi3 - 140*xi2*abs_xi2 + 604*xi2*abs_xi + 28*xi2 +
        40*abs_xi3 + 56*abs_xi2 - 560*abs_xi)/144.0/abs_xi;
    } else if (abs_xi <= 2) {
      dphi = xi*(126*xi4*abs_xi - 21*xi4*abs_xi2 - 182*xi4 -
        28*xi2*abs_xi4 + 300*xi2*abs_xi3 - 1001*xi2*abs_xi2 +
        990*xi2*abs_xi + 154*xi2 + 24*abs_xi5 - 182*abs_xi4 +
        270*abs_xi3 + 602*abs_xi2 - 1260*abs_xi + 28)/240.0/abs_xi;
    } else if (abs_xi <= 3) {
      dphi = xi*(35*xi2*abs_xi4 - 420*xi2*abs_xi3 +
        1785*xi2*abs_xi2 - 3150*xi2*abs_xi + 1918*xi2 +
        14*abs_xi6 - 330*abs_xi5 + 2660*abs_xi4 -
        9590*abs_xi3 + 15806*abs_xi2 - 9940*abs_xi + 756)/720.0/abs_xi;
    } else if (abs_xi <= 4) {
      dphi = -xi*(abs_xi - 4)*(7*abs_xi5 - 122*abs_xi4 +
        807*abs_xi3 - 2512*abs_xi2 + 3644*abs_xi - 1944)/720.0/abs_xi;
    } else {
      dphi = 0.0;
    }

  } else if (order == 10) {
    double xi2 = xi*xi;
    double xi4 = xi2*xi2;
    double xi6 = xi4*xi2;
    double xi8 = xi6*xi2;
    double abs_xi2 = abs_xi*abs_xi;
    double abs_xi3 = abs_xi2*abs_xi;
    double abs_xi4 = abs_xi2*abs_xi2;
    double abs_xi5 = abs_xi4*abs_xi;
    double abs_xi6 = abs_xi5*abs_xi;
    double abs_xi7 = abs_xi6*abs_xi;
    double abs_xi8 = abs_xi7*abs_xi;
    if (abs_xi == 0.0) {
      dphi = 0.0;
    } else if (abs_xi <= 1) {
      dphi = xi*(298*xi6*abs_xi - 72*xi6*abs_xi2 - 9*xi8 +
        126*xi6 + 30*xi4*abs_xi3 + 756*xi4*abs_xi2 - 3644*xi4*abs_xi -
        441*xi4 - 280*xi2*abs_xi3 - 1764*xi2*abs_xi2 + 12026*xi2*abs_xi +
        324*xi2 + 490*abs_xi3 + 648*abs_xi2 - 10792*abs_xi)/2880.0/abs_xi;
    } else if (abs_xi <= 2) {
      dphi = xi*(9*xi6*abs_xi2 - 72*xi6*abs_xi + 141*xi6 +
        18*xi4*abs_xi4 - 236*xi4*abs_xi3 + 963*xi4*abs_xi2 -
        1046*xi4*abs_xi - 687*xi4 - 20*xi2*abs_xi5 + 156*xi2*abs_xi4 +
        168*xi2*abs_xi3 - 3522*xi2*abs_xi2 + 6382*xi2*abs_xi + 474*xi2 +
        50*abs_xi5 - 516*abs_xi4 + 1262*abs_xi3 + 1596*abs_xi2 -
        6344*abs_xi + 72)/1440.0/abs_xi;
    } else if (abs_xi <= 3) {
      dphi = xi*(720*xi4*abs_xi3 - 45*xi4*abs_xi4 - 4185*xi4*abs_xi2 +
        10440*xi4*abs_xi - 9396*xi4 - 36*xi2*abs_xi6 + 870*xi2*abs_xi5 -
        7965*xi2*abs_xi4 + 34540*xi2*abs_xi3 - 70389*xi2*abs_xi2 +
        51440*xi2*abs_xi + 6012*xi2 + 50*abs_xi7 - 954*abs_xi6 +
        6680*abs_xi5 - 19440*abs_xi4 + 11140*abs_xi3 + 49014*abs_xi2 -
        69080*abs_xi + 3384)/10080.0/abs_xi;
    } else if (abs_xi <= 4) {
      dphi = xi*(63*xi2*abs_xi6 - 1512*xi2*abs_xi5 + 14490*xi2*abs_xi4 -
        70560*xi2*abs_xi3 + 182763*xi2*abs_xi2 - 236376*xi2*abs_xi +
        117612*xi2 + 18*abs_xi8 - 784*abs_xi7 + 12600*abs_xi6 -
        101556*abs_xi5 + 451962*abs_xi4 - 1121316*abs_xi3 +
        1451628*abs_xi2 - 795368*abs_xi + 71856)/40320.0/abs_xi;
    } else if (abs_xi <= 5) {
      dphi = -xi*(abs_xi - 5)*(9*abs_xi7 - 283*abs_xi6 +
        3667*abs_xi5 - 25261*abs_xi4 + 99340*abs_xi3 -
        221416*abs_xi2 + 256552*abs_xi - 117648)/40320.0/abs_xi;
    } else {
      dphi = 0.0;
    }
  }

  return dphi;
}

/* ----------------------------------------------------------------------
   Compute direct interaction for each grid level
------------------------------------------------------------------------- */
void MSM::get_g_direct()
{
  if (g_direct) memory->destroy(g_direct);
  memory->create(g_direct,levels,nmax_direct,"msm:g_direct");

  double a = cutoff;

  int n,k,ix,iy,iz;
  double xdiff,ydiff,zdiff;
  double rsq,rho,two_n;

  two_n = 1.0;

  for (n=0; n<levels; n++) {

    k = 0;
    for (iz = nzlo_direct; iz <= nzhi_direct; iz++) {
      zdiff = iz/delzinv[n];
      for (iy = nylo_direct; iy <= nyhi_direct; iy++) {
        ydiff = iy/delyinv[n];
        for (ix = nxlo_direct; ix <= nxhi_direct; ix++) {
          xdiff = ix/delxinv[n];
          rsq = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff;
          rho = sqrt(rsq)/(two_n*a);
          g_direct[n][k] = gamma(rho)/(two_n*a) - gamma(rho/2.0)/(2.0*two_n*a);
          k++;
        }
      }
    }
    two_n *= 2.0;
  }
}

/* ----------------------------------------------------------------------
   Compute direct interaction for each grid level
------------------------------------------------------------------------- */
void MSM::get_dg_direct()
{
  if (dgx_direct) memory->destroy(dgx_direct);
  memory->create(dgx_direct,levels,nmax_direct,"msm:dgx_direct");
  if (dgy_direct) memory->destroy(dgy_direct);
  memory->create(dgy_direct,levels,nmax_direct,"msm:dgy_direct");
  if (dgz_direct) memory->destroy(dgz_direct);
  memory->create(dgz_direct,levels,nmax_direct,"msm:dgz_direct");

  double a = cutoff;
  double a_sq = cutoff*cutoff;

  int n,k,ix,iy,iz;
  double xdiff,ydiff,zdiff;
  double rsq,r,rho,two_n,two_nsq,dg;

  two_n = 1.0;

  for (n=0; n<levels; n++) {
    two_nsq = two_n * two_n;

    k = 0;
    for (iz = nzlo_direct; iz <= nzhi_direct; iz++) {
      zdiff = iz/delzinv[n];
      for (iy = nylo_direct; iy <= nyhi_direct; iy++) {
        ydiff = iy/delyinv[n];
        for (ix = nxlo_direct; ix <= nxhi_direct; ix++) {
          xdiff = ix/delxinv[n];
          rsq = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff;
          r = sqrt(rsq);
          if (r == 0) {
            dgx_direct[n][k] = 0.0;
            dgy_direct[n][k] = 0.0;
            dgz_direct[n][k] = 0.0;
          } else {
            rho = r/(two_n*a);
            dg = -(dgamma(rho)/(two_nsq*a_sq) -
              dgamma(rho/2.0)/(4.0*two_nsq*a_sq))/r;
            dgx_direct[n][k] = dg * xdiff;
            dgy_direct[n][k] = dg * ydiff;
            dgz_direct[n][k] = dg * zdiff;
          }
          k++;
        }
      }
    }
    two_n *= 2.0;
  }
}

/* ----------------------------------------------------------------------
   Compute direct interaction for each grid level
------------------------------------------------------------------------- */
void MSM::get_virial_direct()
{
  if (v0_direct) memory->destroy(v0_direct);
  memory->create(v0_direct,levels,nmax_direct,"msm:v0_direct");
  if (v1_direct) memory->destroy(v1_direct);
  memory->create(v1_direct,levels,nmax_direct,"msm:v1_direct");
  if (v2_direct) memory->destroy(v2_direct);
  memory->create(v2_direct,levels,nmax_direct,"msm:v2_direct");
  if (v3_direct) memory->destroy(v3_direct);
  memory->create(v3_direct,levels,nmax_direct,"msm:v3_direct");
  if (v4_direct) memory->destroy(v4_direct);
  memory->create(v4_direct,levels,nmax_direct,"msm:v4_direct");
  if (v5_direct) memory->destroy(v5_direct);
  memory->create(v5_direct,levels,nmax_direct,"msm:v5_direct");

  int n,k,ix,iy,iz;
  double xdiff,ydiff,zdiff;

  for (n=0; n<levels; n++) {

    k = 0;
    for (iz = nzlo_direct; iz <= nzhi_direct; iz++) {
      zdiff = iz/delzinv[n];
      for (iy = nylo_direct; iy <= nyhi_direct; iy++) {
        ydiff = iy/delyinv[n];
        for (ix = nxlo_direct; ix <= nxhi_direct; ix++) {
          xdiff = ix/delxinv[n];
          v0_direct[n][k] = dgx_direct[n][k] * xdiff;
          v1_direct[n][k] = dgy_direct[n][k] * ydiff;
          v2_direct[n][k] = dgz_direct[n][k] * zdiff;
          v3_direct[n][k] = dgx_direct[n][k] * ydiff;
          v4_direct[n][k] = dgx_direct[n][k] * zdiff;
          v5_direct[n][k] = dgy_direct[n][k] * zdiff;
          k++;
        }
      }
    }

  }
}
