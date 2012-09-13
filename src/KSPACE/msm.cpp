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

#define MAXORDER 7
#define MAX_LEVELS 10
#define OFFSET 16384
#define SMALL 0.00001
#define LARGE 10000.0
#define EPS_HOC 1.0e-7

/* ---------------------------------------------------------------------- */

MSM::MSM(LAMMPS *lmp, int narg, char **arg) : KSpace(lmp, narg, arg)
{
  if (narg < 1) error->all(FLERR,"Illegal kspace_style msm command");

  accuracy_relative = atof(arg[0]);

  nfactors = 1;
  factors = new int[nfactors];
  factors[0] = 2;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  buf1 = buf2 = NULL;

  phi1d = dphi1d = NULL;

  nmax = 0;
  part2grid = NULL;

  g_direct = NULL;

  levels = 0;

}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

MSM::~MSM()
{
  delete [] factors;
  deallocate();
  memory->destroy(part2grid);
  memory->destroy(g_direct);
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

  for (int i = 0; i < modify->nfix; i++) {
    if ((strcmp(modify->fix[i]->style,"npt") == 0) ||
        (strcmp(modify->fix[i]->style,"nph") == 0)) {
      error->all(FLERR,"Cannot (yet) use MSM with npt/nph simulation");
    }
  }

  order = 4;  // 4th order interpolation scheme has been implemented for MSM

  if (order > MAXORDER) {
    char str[128];
    sprintf(str,"MSM order cannot be greater than %d",MAXORDER);
    error->all(FLERR,str);
  }

  // free all arrays previously allocated

  deallocate();

  // extract short-range Coulombic cutoff from pair style

  qqrd2e = force->qqrd2e;
  scale = 1.0;

  if (force->pair == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  int itmp;
  double *p_cutoff = (double *) force->pair->extract("cut_msm",itmp);
  if (p_cutoff == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;

  if ((strcmp(force->kspace_style,"pppm/tip4p") == 0) ||
       (strcmp(force->kspace_style,"pppm/tip4p/proxy") == 0)) {
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  }

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
  if (fabs(qsum) > SMALL) {
    char str[128];
    sprintf(str,"System is not charge neutral, net charge = %g",qsum);
    error->all(FLERR,str);  // Not yet sure of the correction needed for non-neutral systems
  }

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;

  // setup grid resolution

  set_grid();

  int flag_global = 0;

  // loop over grid levels

  for (int n=0; n<levels; n++) {

    if (nx_msm[n] >= OFFSET || ny_msm[n] >= OFFSET || nz_msm[n] >= OFFSET)
      error->all(FLERR,"MSM grid is too large");

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
  }

  // Largest MSM grid for this proc, including ghosts

  ngrid = (nxhi_out[0]-nxlo_out[0]+1) * (nyhi_out[0]-nylo_out[0]+1) *
    (nzhi_out[0]-nzlo_out[0]+1);


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
  nbuf *= 3;

  double estimated_error = estimate_total_error();

  // print stats

  int ngrid_max,nbuf_max;

  // All processors have a copy of the complete grid at each level

  nbuf_max = nbuf;
  ngrid_max = ngrid;

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
   estimate 1d grid RMS force error for MSM from Sec 3.1 of Hardy's thesis
------------------------------------------------------------------------- */

double MSM::estimate_1d_error(double h, double prd)
{
  double a = cutoff;
  int p = order - 1;
  double error_1d = pow(h,(p-1))/pow(a,(p+1));
  error_1d *= 24.0*q2/sqrt(atom->natoms*a*prd*prd);
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

  get_g_direct();

  boxlo = domain->boxlo;
}

/* ----------------------------------------------------------------------
   compute the MSM long-range force, energy, virial
------------------------------------------------------------------------- */

void MSM::compute(int eflag, int vflag)
{
  int i;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  // extend size of per-atom arrays if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"msm:part2grid");
  }

  energy = 0.0;

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid (aninterpolation)

  particle_map();
  make_rho();


  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks

  ghost_swap(0);

  charge_swap(0);

  // Direct sum on finest grid level is parallel

  direct(eflag_global,vflag_global,0);

  potential_swap(0);

  restrict(eflag_global,vflag_global,0);

  // compute potential gradient on my MSM grid and
  //   portion of e_long on this proc's MSM grid
  // return gradients (electric fields) in 3d brick decomposition

  for (int n=1; n<levels; n++) {
    direct(eflag_global,vflag_global,n);

    if (n < levels-1) restrict(eflag_global,vflag_global,n);
  }

  for (int n=levels-2; n>=0; n--)
    prolongate(eflag_global,vflag_global,n);

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  fillbrick(0);

  // calculate the force on my particles (interpolation)

  fieldforce();

  // sum energy across procs and add in volume-dependent term

  if (eflag_global) {
    double e_self = qsqsum*gamma(0.0)/cutoff;  // Self-energy term
    energy -= e_self;
    double energy_all;
    energy *= 0.5*qqrd2e*scale;
  }

}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void MSM::allocate()
{
  memory->create(buf1,nbuf,"msm:buf1");
  memory->create(buf2,nbuf,"msm:buf2");

  // summation coeffs

  memory->create2d_offset(phi1d,3,nlower-2,nupper+2,"msm:phi1d");
  memory->create2d_offset(dphi1d,3,nlower-2,nupper+2,"msm:dphi1d");

  // allocate grid levels

  for (int n=0; n<levels; n++) {
    memory->create3d_offset(qgrid[n],nzlo_out[n],nzhi_out[n],
            nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:qgrid");
    memory->create3d_offset(egrid[n],nzlo_out[n],nzhi_out[n],
            nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:egrid");
  }

}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void MSM::deallocate()
{
  memory->destroy(buf1);
  memory->destroy(buf2);

  memory->destroy2d_offset(phi1d,nlower-2);
  memory->destroy2d_offset(dphi1d,nlower-2);

  // deallocate grid levels

  for (int n=0; n<levels; n++) {
    if (qgrid[n]) 
      memory->destroy3d_offset(qgrid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
    if (egrid[n]) 
      memory->destroy3d_offset(egrid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);
  }
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of grid levels
------------------------------------------------------------------------- */

void MSM::allocate_levels()
{
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

}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of grid levels
------------------------------------------------------------------------- */

void MSM::deallocate_levels()
{

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
   MPI-Reduce so each processor has all the info it needs
------------------------------------------------------------------------- */
void MSM::charge_swap(int n)
{
  double ***qgridn = qgrid[n];
  double ***qgridn_all;
  memory->create3d_offset(qgridn_all,nzlo_out[n],nzhi_out[n],nylo_out[n],nyhi_out[n],
                          nxlo_out[n],nxhi_out[n],"msm:qgrid_all");

  memset(&(qgridn_all[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid*sizeof(double));

  MPI_Allreduce(&(qgridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),
                &(qgridn_all[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),
                ngrid,MPI_DOUBLE,MPI_SUM,world);

  // Swap pointers between qgridn and qgridn_all to avoid need of copy operation

  double ***tmp;
  tmp = qgridn;
  qgrid[n] = qgridn_all;
  qgridn_all = tmp;

  memory->destroy3d_offset(qgridn_all,nzlo_out[n],nylo_out[n],nxlo_out[n]);

}

/* ----------------------------------------------------------------------
   MPI-Reduce so each processor has all the info it needs
------------------------------------------------------------------------- */
void MSM::potential_swap(int n)
{
  double ***egridn = egrid[n];
  double ***egridn_all;
  memory->create3d_offset(egridn_all,nzlo_out[n],nzhi_out[n],nylo_out[n],nyhi_out[n],
                          nxlo_out[n],nxhi_out[n],"msm:qgrid_all");

  memset(&(egridn_all[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid*sizeof(double));

  MPI_Allreduce(&(egridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),
                &(egridn_all[nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),
                ngrid,MPI_DOUBLE,MPI_SUM,world);

  // Swap pointers between egridn and egridn_all to avoid need of copy operation

  double ***tmp;
  tmp = egridn;
  egrid[n] = egridn_all;
  egridn_all = tmp;

  memory->destroy3d_offset(egridn_all,nzlo_out[n],nylo_out[n],nxlo_out[n]);

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

  memset(&(qgridn[nzlo_out[0]][nylo_out[0]][nxlo_out[0]]),0,ngrid*sizeof(double));

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

void MSM::direct(int eflag, int vflag, int n)
{
  //fprintf(screen,"Direct contribution on level %i\n\n",n);

  double ***egridn = egrid[n];
  double ***qgridn = qgrid[n];

  // bitmask for PBCs (only works for power of 2 numbers)

  int PBCx,PBCy,PBCz;

  PBCx = nx_msm[n]-1;
  PBCy = ny_msm[n]-1;
  PBCz = nz_msm[n]-1;

  // zero out electric field brick

  for (int icz = nzlo_in[n]; icz <= nzhi_in[n]; icz++)
    for (int icy = nylo_in[n]; icy <= nyhi_in[n]; icy++)
      for (int icx = nxlo_in[n]; icx <= nxhi_in[n]; icx++)
        egridn[icz][icy][icx] = 0.0;

  // Simple parallelization of direct sum

  for (int icz = nzlo_in_d[n]; icz <= nzhi_in_d[n]; icz++) {
    for (int icy = nylo_in_d[n]; icy <= nyhi_in_d[n]; icy++) {
      for (int icx = nxlo_in_d[n]; icx <= nxhi_in_d[n]; icx++) {


  // do double loop over points on the intermediate grid level
  // for now, assume I own all points on the intermediate grid level

        int k = 0;
        for (int iz = nzlo_direct; iz <= nzhi_direct; iz++) {
          for (int iy = nylo_direct; iy <= nyhi_direct; iy++) {
            for (int ix = nxlo_direct; ix <= nxhi_direct; ix++) {
              egridn[icz][icy][icx] += g_direct[n][k++]
                * qgridn[(icz+iz) & PBCz][(icy+iy) & PBCy][(icx+ix) & PBCx];
            }
          }
        }

      }
    }
  }

}

/* ----------------------------------------------------------------------
   MSM restrict  procedure for intermediate grid levels
------------------------------------------------------------------------- */

void MSM::restrict(int eflag, int vflag,int n)
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

  for (int kp = nzlo_in[n+1]; kp <= nzhi_in[n+1]; kp++)
    for (int jp = nylo_in[n+1]; jp <= nyhi_in[n+1]; jp++)
      for (int ip = nxlo_in[n+1]; ip <= nxhi_in[n+1]; ip++) {
        qgrid2[kp][jp][ip] = 0.0;

        int ic = static_cast<int> (ip*delxinv[n]/delxinv[n+1]);
        int jc = static_cast<int> (jp*delyinv[n]/delyinv[n+1]);
        int kc = static_cast<int> (kp*delzinv[n]/delzinv[n+1]);
        for (int k=-p; k<=p; k++) // Could make this faster by eliminating zeros
          for (int j=-p; j<=p; j++)
            for (int i=-p; i<=p; i++)
              qgrid2[kp][jp][ip] += 
               qgrid1[(kc+k)&PBCz][(jc+j)&PBCy][(ic+i)&PBCx] * 
               phi1d[0][i]*phi1d[1][j]*phi1d[2][k];
      }


}

/* ----------------------------------------------------------------------
   MSM prolongate procedure for intermediate grid levels
------------------------------------------------------------------------- */

void MSM::prolongate(int eflag, int vflag,int n)
{
  //fprintf(screen,"Prolongating from level %i to %i\n\n",n+1,n);

  int p = order-1;

  double ***egrid1 = egrid[n];
  double ***egrid2 = egrid[n+1];

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

  for (int kp = nzlo_in[n+1]; kp <= nzhi_in[n+1]; kp++)
    for (int jp = nylo_in[n+1]; jp <= nyhi_in[n+1]; jp++)
      for (int ip = nxlo_in[n+1]; ip <= nxhi_in[n+1]; ip++) {

        int ic = static_cast<int> (ip*delxinv[n]/delxinv[n+1]);
        int jc = static_cast<int> (jp*delyinv[n]/delyinv[n+1]);
        int kc = static_cast<int> (kp*delzinv[n]/delzinv[n+1]);
        for (int k=-p; k<=p; k++) // Could make this faster by eliminating zeros
          for (int j=-p; j<=p; j++)
            for (int i=-p; i<=p; i++)
              egrid1[(kc+k)&PBCz][(jc+j)&PBCy][(ic+i)&PBCx] += 
                egrid2[kp][jp][ip] * phi1d[0][i]*phi1d[1][j]*phi1d[2][k];
        }

}

/* ----------------------------------------------------------------------
   ghost-swap to fill ghost cells of my brick with field values
------------------------------------------------------------------------- */

void MSM::fillbrick(int n)
{
  double ***egridn = egrid[n];

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
        buf1[k++] = egridn[iz][iy][ix];
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
        egridn[iz][iy][ix] = buf2[k++];
      }

  // pack my real cells for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_in[n]; iz < nzlo_in[n]+nzlo_ghost[n]; iz++)
    for (iy = nylo_in[n]; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        buf1[k++] = egridn[iz][iy][ix];
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
        egridn[iz][iy][ix] = buf2[k++];
      }

  // pack my real cells for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nyhi_in[n]-nyhi_ghost[n]+1; iy <= nyhi_in[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        buf1[k++] = egridn[iz][iy][ix];
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
        egridn[iz][iy][ix] = buf2[k++];
      }

  // pack my real cells for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_in[n]; iy < nylo_in[n]+nylo_ghost[n]; iy++)
      for (ix = nxlo_in[n]; ix <= nxhi_in[n]; ix++) {
        buf1[k++] = egridn[iz][iy][ix];
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
        egridn[iz][iy][ix] = buf2[k++];
      }

  // pack my real cells for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxhi_in[n]-nxhi_ghost[n]+1; ix <= nxhi_in[n]; ix++) {
        buf1[k++] = egridn[iz][iy][ix];
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
      for (ix = nzlo_out[n]; ix < nxlo_in[n]; ix++) {
        egridn[iz][iy][ix] = buf2[k++];
      }

  // pack my real cells for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my ghost cells

  k = 0;
  for (iz = nzlo_out[n]; iz <= nzhi_out[n]; iz++)
    for (iy = nylo_out[n]; iy <= nyhi_out[n]; iy++)
      for (ix = nxlo_in[n]; ix < nxlo_in[n]+nxlo_ghost[n]; ix++) {
        buf1[k++] = egridn[iz][iy][ix];
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
        egridn[iz][iy][ix] = buf2[k++];
      }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get force on my particles
------------------------------------------------------------------------- */

void MSM::fieldforce()
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

  // Sum total long-range energy

  for (int kp=0; kp<nz_msm[0]; kp++)
    for (int jp=0; jp<ny_msm[0]; jp++)
      for (int ip=0; ip<nx_msm[0]; ip++) {
        energy += egridn[kp][jp][ip]*qgridn[kp][jp][ip];
      }

}

/* ----------------------------------------------------------------------
   charge assignment into phi1d
------------------------------------------------------------------------- */
void MSM::compute_phis_and_dphis(const double &dx, const double &dy, const double &dz)
{
  for (int nu = nlower; nu <= nupper; nu++) {
    double delx = dx + double(nu);
    double dely = dy + double(nu);
    double delz = dz + double(nu);

    //fprintf(screen,"delx = %f, phi = %f\n",delx,compute_phi(delx));

    phi1d[0][nu] = compute_phi(delx);
    phi1d[1][nu] = compute_phi(dely);
    phi1d[2][nu] = compute_phi(delz);
    dphi1d[0][nu] = compute_dphi(delx);
    dphi1d[1][nu] = compute_dphi(dely);
    dphi1d[2][nu] = compute_dphi(delz);
  }
}

/* ----------------------------------------------------------------------
   compute phi using quadratic interpolating polynomial
   see Eq 7 from Parallel Computing 35 (2009) 164–177
------------------------------------------------------------------------- */
double MSM::compute_phi(const double &xi)
{
   double phi;
   double abs_xi = fabs(xi);
   if (abs_xi <= 1) {
     phi = (1 - abs_xi)*(1 + abs_xi - 1.5*abs_xi*abs_xi);
   } else if (abs_xi <= 2) {
     phi = -0.5*(abs_xi - 1)*(2 - abs_xi)*(2 - abs_xi);
   } else {
     phi = 0.0;
   }
   return phi;
}

/* ----------------------------------------------------------------------
   compute the derivative of phi
   phi is a quadratic interpolating polynomial
   see Eq 7 from Parallel Computing 35 (2009) 164–177
------------------------------------------------------------------------- */

double MSM::compute_dphi(const double &xi)
{
   double dphi;
   double abs_xi = fabs(xi);
   if (abs_xi == 0.0) {
     dphi = 0.0;
   } else if (abs_xi <= 1) {
     dphi = 3*xi*xi*xi/(2*abs_xi) - 5*xi + 3*xi*abs_xi;
   } else if (abs_xi <= 2) {
     dphi = xi*(2 - abs_xi)*(3*abs_xi - 4)/(2*abs_xi);
   } else {
     dphi = 0.0;
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

  for (int n=0; n<levels; n++) {
    double two_to_n = pow(2.0,n);

    int k = 0;
    for (int iz = nzlo_direct; iz <= nzhi_direct; iz++) {
      double zdiff = iz/delzinv[n];
      for (int iy = nylo_direct; iy <= nyhi_direct; iy++) {
        double ydiff = iy/delyinv[n];
        for (int ix = nxlo_direct; ix <= nxhi_direct; ix++) {
          double xdiff = ix/delxinv[n];
          double rsq = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff;
          double rho = sqrt(rsq)/(two_to_n*a);
          g_direct[n][k++] = gamma(rho)/(two_to_n*a) - gamma(rho/2.0)/(2.0*two_to_n*a);
        }
      }
    }

  }
}
