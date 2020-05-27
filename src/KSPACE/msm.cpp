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

#include "msm.h"
#include <mpi.h>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "gridcomm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAX_LEVELS 10
#define OFFSET 16384
#define SMALL 0.00001

enum{REVERSE_RHO,REVERSE_AD,REVERSE_AD_PERATOM};
enum{FORWARD_RHO,FORWARD_AD,FORWARD_AD_PERATOM};
/* ---------------------------------------------------------------------- */

MSM::MSM(LAMMPS *lmp) : KSpace(lmp),
  factors(NULL), delxinv(NULL), delyinv(NULL), delzinv(NULL), nx_msm(NULL),
  ny_msm(NULL), nz_msm(NULL), nxlo_in(NULL), nylo_in(NULL), nzlo_in(NULL),
  nxhi_in(NULL), nyhi_in(NULL), nzhi_in(NULL), nxlo_out(NULL), nylo_out(NULL),
  nzlo_out(NULL), nxhi_out(NULL), nyhi_out(NULL), nzhi_out(NULL), ngrid(NULL),
  active_flag(NULL), alpha(NULL), betax(NULL), betay(NULL), betaz(NULL), peratom_allocate_flag(0),
  levels(0), world_levels(NULL), qgrid(NULL), egrid(NULL), v0grid(NULL), v1grid(NULL),
  v2grid(NULL), v3grid(NULL), v4grid(NULL), v5grid(NULL), g_direct(NULL),
  v0_direct(NULL), v1_direct(NULL), v2_direct(NULL), v3_direct(NULL), v4_direct(NULL),
  v5_direct(NULL), g_direct_top(NULL), v0_direct_top(NULL), v1_direct_top(NULL),
  v2_direct_top(NULL), v3_direct_top(NULL), v4_direct_top(NULL), v5_direct_top(NULL),
  phi1d(NULL), dphi1d(NULL), procneigh_levels(NULL), cg(NULL), cg_peratom(NULL),
  cg_all(NULL), cg_peratom_all(NULL), part2grid(NULL), boxlo(NULL)
{
  msmflag = 1;

  nfactors = 1;
  factors = new int[nfactors];
  factors[0] = 2;

  MPI_Comm_rank(world,&me);
  procneigh_levels = NULL;
  world_levels = NULL;
  active_flag = NULL;

  phi1d = dphi1d = NULL;

  nmax = 0;
  part2grid = NULL;

  g_direct = NULL;
  g_direct_top = NULL;

  v0_direct = v1_direct = v2_direct = NULL;
  v3_direct = v4_direct = v5_direct = NULL;

  v0_direct_top = v1_direct_top = v2_direct_top = NULL;
  v3_direct_top = v4_direct_top = v5_direct_top = NULL;

  cg_all = cg_peratom_all = NULL;
  cg = cg_peratom = NULL;

  ngrid = NULL;
  cg = NULL;
  cg_peratom = NULL;
  procneigh_levels = NULL;
  world_levels = NULL;
  active_flag = NULL;

  alpha = betax = betay = betaz = NULL;
  nx_msm = ny_msm = nz_msm = NULL;
  nxlo_in = nylo_in = nzlo_in = NULL;
  nxhi_in = nyhi_in = nzhi_in = NULL;
  nxlo_out = nylo_out = nzlo_out = NULL;
  nxhi_out = nyhi_out = nzhi_out = NULL;
  delxinv = delyinv = delzinv = NULL;
  qgrid = NULL;
  egrid = NULL;
  v0grid = v1grid = v2grid = v3grid = v4grid = v5grid = NULL;

  peratom_allocate_flag = 0;
  scalar_pressure_flag = 1;
  warn_nonneutral = 0;

  order = 10;
}

/* ---------------------------------------------------------------------- */

void MSM::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal kspace_style msm command");
  accuracy_relative = fabs(force->numeric(FLERR,arg[0]));
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

MSM::~MSM()
{
  delete [] factors;
  deallocate();
  if (peratom_allocate_flag) deallocate_peratom();
  memory->destroy(part2grid);
  memory->destroy(g_direct);
  memory->destroy(g_direct_top);
  memory->destroy(v0_direct);
  memory->destroy(v1_direct);
  memory->destroy(v2_direct);
  memory->destroy(v3_direct);
  memory->destroy(v4_direct);
  memory->destroy(v5_direct);
  memory->destroy(v0_direct_top);
  memory->destroy(v1_direct_top);
  memory->destroy(v2_direct_top);
  memory->destroy(v3_direct_top);
  memory->destroy(v4_direct_top);
  memory->destroy(v5_direct_top);
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

  triclinic_check();
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot (yet) use MSM with 2d simulation");
  if (comm->style != 0)
    error->universe_all(FLERR,"MSM can only currently be used with "
                        "comm_style brick");

  if (!atom->q_flag) error->all(FLERR,"Kspace style requires atom attribute q");

  if ((slabflag == 1) && (me == 0))
    error->warning(FLERR,"Slab correction not needed for MSM");

  if (order < 4 || order > 10) {
    char str[128];
    sprintf(str,"MSM order must be 4, 6, 8, or 10");
    error->all(FLERR,str);
  }

  if (order%2 != 0) error->all(FLERR,"MSM order must be 4, 6, 8, or 10");

  if (sizeof(FFT_SCALAR) != 8)
    error->all(FLERR,"Cannot (yet) use single precision with MSM "
               "(remove -DFFT_SINGLE from Makefile and re-compile)");

  // compute two charge force

  two_charge();

  // extract short-range Coulombic cutoff from pair style

  triclinic = domain->triclinic;
  pair_check();

  int itmp;
  double *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  if (p_cutoff == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;

  // compute qsum & qsqsum and error if not charge-neutral

  scale = 1.0;
  qqrd2e = force->qqrd2e;
  qsum_qsq();
  natoms_original = atom->natoms;

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;

  // setup MSM grid resolution

  set_grid_global();
  setup();

  double estimated_error = estimate_total_error();

  // output grid stats

  int ngrid_max;
  MPI_Allreduce(&ngrid[0],&ngrid_max,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  3d grid size/proc = %d\n",
                        ngrid_max);
      fprintf(screen,"  estimated absolute RMS force accuracy = %g\n",
              estimated_error);
      fprintf(screen,"  estimated relative force accuracy = %g\n",
              estimated_error/two_charge_force);
    }
    if (logfile) {
      fprintf(logfile,"  3d grid size/proc = %d\n",
                         ngrid_max);
      fprintf(logfile,"  estimated absolute RMS force accuracy = %g\n",
              estimated_error);
      fprintf(logfile,"  estimated relative force accuracy = %g\n",
              estimated_error/two_charge_force);
    }
  }

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  grid = %d %d %d\n",nx_msm[0],ny_msm[0],nz_msm[0]);
      fprintf(screen,"  order = %d\n",order);
    }
    if (logfile) {
      fprintf(logfile,"  grid = %d %d %d\n",nx_msm[0],ny_msm[0],nz_msm[0]);
      fprintf(logfile,"  order = %d\n",order);
    }
  }
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
    error_scaling = 0.39189561;
  } else if (p == 5) {
    Mp = 825;
    cprime = 1.0/30.0;
    error_scaling = 0.150829428;
  } else if (p == 7) {
    Mp = 130095;
    cprime = 1.0/140.0;
    error_scaling = 0.049632967;
  } else if (p == 9) {
    Mp = 34096545;
    cprime = 1.0/630.0;
    error_scaling = 0.013520855;
  } else {
    error->all(FLERR,"MSM order must be 4, 6, 8, or 10");
  }

  // equation 4.1 from Hardy's thesis
  C_p = 4.0*cprime*Mp/3.0;

  // use empirical parameters to convert to rms force errors
  C_p *= error_scaling;

  // equation 3.197 from Hardy's thesis
  double error_1d = C_p*pow(h,(p-1))/pow(a,(p+1));

  // include dependency of error on other terms
  error_1d *= q2*a/(prd*sqrt(double(atom->natoms)));

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
  double error_x = estimate_1d_error(h_x,xprd);
  double error_y = estimate_1d_error(h_y,yprd);
  double error_z = estimate_1d_error(h_z,zprd);
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
  double *prd;
  double a = cutoff;

  // volume-dependent factors

  prd = domain->prd;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  volume = xprd * yprd * zprd;

  // loop over grid levels and compute grid spacing

  for (int n=0; n<levels; n++) {
    if (triclinic == 0) {
      delxinv[n] = nx_msm[n]/xprd;
      delyinv[n] = ny_msm[n]/yprd;
      delzinv[n] = nz_msm[n]/zprd;
    } else { // use lamda (0-1) coordinates
      delxinv[n] = nx_msm[n];
      delyinv[n] = ny_msm[n];
      delzinv[n] = nz_msm[n];
    }
  }

  double ax = a;
  double ay = a;
  double az = a;

  // transform the interaction sphere in box coords to an
  // ellipsoid in lamda (0-1) coords to
  // get the direct sum interaction limits for a triclinic system

  if (triclinic) {
    double tmp[3];
    kspacebbox(a,&tmp[0]);
    ax = tmp[0];
    ay = tmp[1];
    az = tmp[2];
  }

  // direct sum interaction limits

  nxhi_direct = static_cast<int> (2.0*ax*delxinv[0]);
  nxlo_direct = -nxhi_direct;
  nyhi_direct = static_cast<int> (2.0*ay*delyinv[0]);
  nylo_direct = -nyhi_direct;
  nzhi_direct = static_cast<int> (2.0*az*delzinv[0]);
  nzlo_direct = -nzhi_direct;

  nmax_direct = 8*(nxhi_direct+1)*(nyhi_direct+1)*(nzhi_direct+1);

  deallocate();
  if (peratom_allocate_flag) deallocate_peratom();

  // compute direct sum interaction weights

  if (!peratom_allocate_flag) { // Timestep 0
    get_g_direct();
    get_virial_direct();
    if (domain->nonperiodic) {
      get_g_direct_top(levels-1);
      get_virial_direct_top(levels-1);
    }
  } else {
    get_g_direct();
    if (domain->nonperiodic) get_g_direct_top(levels-1);
    if (vflag_either && !scalar_pressure_flag) {
      get_virial_direct();
      if (domain->nonperiodic) get_virial_direct_top(levels-1);
    }
  }

  if (!triclinic)
    boxlo = domain->boxlo;
  else
    boxlo = domain->boxlo_lamda;

  // ghost grid points depend on direct sum interaction limits,
  // so need to re-compute local grid

  set_grid_local();

  // allocate K-space dependent memory
  // don't invoke allocate_peratom(), compute() will allocate when needed

  allocate();

  // setup commgrid

  cg_all->ghost_notify();
  cg_all->setup();
  for (int n=0; n<levels; n++) {
    if (!active_flag[n]) continue;
    cg[n]->ghost_notify();
    cg[n]->setup();
  }

}

/* ----------------------------------------------------------------------
   compute the MSM long-range force, energy, virial
------------------------------------------------------------------------- */

void MSM::compute(int eflag, int vflag)
{
  int i,j;

  // set energy/virial flags

  ev_init(eflag,vflag);

  if (scalar_pressure_flag && vflag_either) {
    if (vflag_atom)
      error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' to obtain "
        "per-atom virial with kspace_style MSM");

    // must switch on global energy computation if not already on

    if (eflag == 0 || eflag == 2) {
      eflag++;
      ev_setup(eflag,vflag);
    }
  }

  // if atom count has changed, update qsum and qsqsum

  if (atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }

  // return if there are no charges

  if (qsqsum == 0.0) return;

  // invoke allocate_peratom() if needed for first time

  if (vflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    cg_peratom_all->ghost_notify();
    cg_peratom_all->setup();
    for (int n=0; n<levels; n++) {
      if (!active_flag[n]) continue;
      cg_peratom[n]->ghost_notify();
      cg_peratom[n]->setup();
    }
  }

  // convert atoms from box to lamda coords

  if (triclinic)
    domain->x2lamda(atom->nlocal);

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"msm:part2grid");
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid (aninterpolation)

  particle_map();
  make_rho();

  // all procs reverse communicate charge density values from
  // their ghost grid points
  // to fully sum contribution in their 3d grid

  current_level = 0;
  cg_all->reverse_comm(this,REVERSE_RHO);

  // forward communicate charge density values to fill ghost grid points
  // compute direct sum interaction and then restrict to coarser grid

  for (int n=0; n<=levels-2; n++) {
    if (!active_flag[n]) continue;
    current_level = n;
    cg[n]->forward_comm(this,FORWARD_RHO);

    direct(n);
    restriction(n);
  }

  // compute direct interation for top grid level for non-periodic
  //   and for second from top grid level for periodic

  if (active_flag[levels-1]) {
    if (domain->nonperiodic) {
      current_level = levels-1;
      cg[levels-1]->forward_comm(this,FORWARD_RHO);
      direct_top(levels-1);
      cg[levels-1]->reverse_comm(this,REVERSE_AD);
      if (vflag_atom)
        cg_peratom[levels-1]->reverse_comm(this,REVERSE_AD_PERATOM);
    } else {
      // Here using MPI_Allreduce is cheaper than using commgrid
      grid_swap_forward(levels-1,qgrid[levels-1]);
      direct(levels-1);
      grid_swap_reverse(levels-1,egrid[levels-1]);
      current_level = levels-1;
      if (vflag_atom)
        cg_peratom[levels-1]->reverse_comm(this,REVERSE_AD_PERATOM);
    }
  }

  // prolongate energy/virial from coarser grid to finer grid
  // reverse communicate from ghost grid points to get full sum

  for (int n=levels-2; n>=0; n--) {
    if (!active_flag[n]) continue;
    prolongation(n);

    current_level = n;
    cg[n]->reverse_comm(this,REVERSE_AD);

    // extra per-atom virial communication

    if (vflag_atom)
      cg_peratom[n]->reverse_comm(this,REVERSE_AD_PERATOM);
  }

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  current_level = 0;
  cg_all->forward_comm(this,FORWARD_AD);

  // extra per-atom energy/virial communication

  if (vflag_atom)
    cg_peratom_all->forward_comm(this,FORWARD_AD_PERATOM);

  // calculate the force on my particles (interpolation)

  fieldforce();

  // calculate the per-atom energy/virial for my particles

  if (evflag_atom) fieldforce_peratom();

  // sum global energy across procs and add in self-energy term

  const double qscale = qqrd2e * scale;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    double e_self = qsqsum*gamma(0.0)/cutoff;
    energy -= e_self;
    energy *= 0.5*qscale;
  }

  // total long-range virial

  if (vflag_global && !scalar_pressure_flag) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*virial_all[i];
  }

  // fast compute of scalar pressure (if requested)

  if (scalar_pressure_flag && vflag_global)
    for (i = 0; i < 3; i++) virial[i] = energy/3.0;

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

  // convert atoms back from lamda to box coords

  if (triclinic)
    domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of grid points
------------------------------------------------------------------------- */

void MSM::allocate()
{
  // interpolation coeffs

  order_allocated = order;
  memory->create2d_offset(phi1d,3,-order,order,"msm:phi1d");
  memory->create2d_offset(dphi1d,3,-order,order,"msm:dphi1d");

  // commgrid using all processors for finest grid level

  int (*procneigh_all)[2] = comm->procneigh;


  cg_all = new GridComm(lmp,world,1,1,
                    nxlo_in[0],nxhi_in[0],nylo_in[0],nyhi_in[0],nzlo_in[0],nzhi_in[0],
                    nxlo_out_all,nxhi_out_all,nylo_out_all,nyhi_out_all,nzlo_out_all,nzhi_out_all,
                    nxlo_out[0],nxhi_out[0],nylo_out[0],nyhi_out[0],nzlo_out[0],nzhi_out[0],
                    procneigh_all[0][0],procneigh_all[0][1],procneigh_all[1][0],
                    procneigh_all[1][1],procneigh_all[2][0],procneigh_all[2][1]);

  // allocate memory for each grid level

  for (int n=0; n<levels; n++) {
    memory->create3d_offset(qgrid[n],nzlo_out[n],nzhi_out[n],
            nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:qgrid");

    memory->create3d_offset(egrid[n],nzlo_out[n],nzhi_out[n],
            nylo_out[n],nyhi_out[n],nxlo_out[n],nxhi_out[n],"msm:egrid");

    // create commgrid object for rho and electric field communication

    if (active_flag[n]) {
      int **procneigh = procneigh_levels[n];
      cg[n] = new GridComm(lmp,world_levels[n],1,1,
                        nxlo_in[n],nxhi_in[n],nylo_in[n],nyhi_in[n],nzlo_in[n],nzhi_in[n],
                        nxlo_out[n],nxhi_out[n],nylo_out[n],nyhi_out[n],nzlo_out[n],nzhi_out[n],
                        procneigh[0][0],procneigh[0][1],procneigh[1][0],
                        procneigh[1][1],procneigh[2][0],procneigh[2][1]);
    }
  }
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of grid points
------------------------------------------------------------------------- */

void MSM::deallocate()
{
  memory->destroy2d_offset(phi1d,-order_allocated);
  memory->destroy2d_offset(dphi1d,-order_allocated);

  if (cg_all) delete cg_all;

  for (int n=0; n<levels; n++) {
    if (qgrid[n])
      memory->destroy3d_offset(qgrid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);

    if (egrid[n])
      memory->destroy3d_offset(egrid[n],nzlo_out[n],nylo_out[n],nxlo_out[n]);

    if (world_levels)
      if (world_levels[n] != MPI_COMM_NULL)
          MPI_Comm_free(&world_levels[n]);

    if (cg)
      if (cg[n]) delete cg[n];
  }
}

/* ----------------------------------------------------------------------
   allocate per-atom virial memory that depends on # of grid points
------------------------------------------------------------------------- */

void MSM::allocate_peratom()
{
  peratom_allocate_flag = 1;

  // create commgrid object for per-atom virial using all processors

  int (*procneigh_all)[2] = comm->procneigh;

  cg_peratom_all =
    new GridComm(lmp,world,6,6,
                 nxlo_in[0],nxhi_in[0],nylo_in[0],nyhi_in[0],nzlo_in[0],nzhi_in[0],
                 nxlo_out_all,nxhi_out_all,nylo_out_all,nyhi_out_all,nzlo_out_all,nzhi_out_all,
                 nxlo_out[0],nxhi_out[0],nylo_out[0],nyhi_out[0],nzlo_out[0],nzhi_out[0],
                 procneigh_all[0][0],procneigh_all[0][1],procneigh_all[1][0],
                 procneigh_all[1][1],procneigh_all[2][0],procneigh_all[2][1]);

  // allocate memory for each grid level

  for (int n=0; n<levels; n++) {
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

    // create commgrid object for per-atom virial

    if (active_flag[n]) {
      int **procneigh = procneigh_levels[n];
      cg_peratom[n] =
        new GridComm(lmp,world_levels[n],6,6,
                     nxlo_in[n],nxhi_in[n],nylo_in[n],nyhi_in[n],nzlo_in[n],nzhi_in[n],
                     nxlo_out[n],nxhi_out[n],nylo_out[n],nyhi_out[n],nzlo_out[n],nzhi_out[n],
                     procneigh[0][0],procneigh[0][1],procneigh[1][0],
                     procneigh[1][1],procneigh[2][0],procneigh[2][1]);
    }
  }
}

/* ----------------------------------------------------------------------
   deallocate per-atom virial memory that depends on # of grid points
------------------------------------------------------------------------- */

void MSM::deallocate_peratom()
{
  peratom_allocate_flag = 0;

  if (cg_peratom_all) delete cg_peratom_all;

  for (int n=0; n<levels; n++) {
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

    if (cg_peratom)
      if (cg_peratom[n]) delete cg_peratom[n];
  }
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of grid levels
------------------------------------------------------------------------- */

void MSM::allocate_levels()
{
  ngrid = new int[levels];

  cg = new GridComm*[levels];
  cg_peratom = new GridComm*[levels];

  memory->create(procneigh_levels,levels,3,2,"msm:procneigh_levels");
  world_levels = new MPI_Comm[levels];
  active_flag = new int[levels];

  alpha = new int[levels];
  betax = new int[levels];
  betay = new int[levels];
  betaz = new int[levels];

  nx_msm = new int[levels];
  ny_msm = new int[levels];
  nz_msm = new int[levels];

  nxlo_in = new int[levels];
  nylo_in = new int[levels];
  nzlo_in = new int[levels];

  nxhi_in = new int[levels];
  nyhi_in = new int[levels];
  nzhi_in = new int[levels];

  nxlo_out = new int[levels];
  nylo_out = new int[levels];
  nzlo_out = new int[levels];

  nxhi_out = new int[levels];
  nyhi_out = new int[levels];
  nzhi_out = new int[levels];

  delxinv = new double[levels];
  delyinv = new double[levels];
  delzinv = new double[levels];

  qgrid = new double***[levels];
  egrid = new double***[levels];

  v0grid = new double***[levels];
  v1grid = new double***[levels];
  v2grid = new double***[levels];
  v3grid = new double***[levels];
  v4grid = new double***[levels];
  v5grid = new double***[levels];

  for (int n=0; n<levels; n++) {
    cg[n] = NULL;
    world_levels[n] = MPI_COMM_NULL;
    cg_peratom[n] = NULL;

    qgrid[n] = NULL;
    egrid[n] = NULL;

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

  memory->destroy(procneigh_levels);
  delete [] world_levels;
  delete [] active_flag;
  delete [] cg;
  delete [] cg_peratom;

  delete [] alpha;
  delete [] betax;
  delete [] betay;
  delete [] betaz;

  delete [] nx_msm;
  delete [] ny_msm;
  delete [] nz_msm;

  delete [] nxlo_in;
  delete [] nylo_in;
  delete [] nzlo_in;

  delete [] nxhi_in;
  delete [] nyhi_in;
  delete [] nzhi_in;

  delete [] nxlo_out;
  delete [] nylo_out;
  delete [] nzlo_out;

  delete [] nxhi_out;
  delete [] nyhi_out;
  delete [] nzhi_out;

  delete [] delxinv;
  delete [] delyinv;
  delete [] delzinv;

  delete [] qgrid;
  delete [] egrid;

  delete [] v0grid;
  delete [] v1grid;
  delete [] v2grid;
  delete [] v3grid;
  delete [] v4grid;
  delete [] v5grid;
}

/* ----------------------------------------------------------------------
   set total size of MSM grids
------------------------------------------------------------------------- */

void MSM::set_grid_global()
{
  if (accuracy_relative <= 0.0)
    error->all(FLERR,"KSpace accuracy must be > 0");

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int nx_max,ny_max,nz_max;
  double hx,hy,hz;

  if (adjust_cutoff_flag && !gridflag) {
    // seek to choose optimal Coulombic cutoff and number of grid levels
    //  (based on a cost estimate in Hardy's thesis)
    int p = order - 1;
    double hmin = 3072.0*(p+1)/(p-1)/
      (448.0*MY_PI + 56.0*MY_PI*order/2 + 1701.0);
    hmin = pow(hmin,1.0/6.0)*pow(xprd*yprd*zprd/atom->natoms,1.0/3.0);

    nx_max = static_cast<int>(xprd/hmin);
    ny_max = static_cast<int>(yprd/hmin);
    nz_max = static_cast<int>(zprd/hmin);

    nx_max = MAX(nx_max,2);
    ny_max = MAX(ny_max,2);
    nz_max = MAX(nz_max,2);

  } else if (!gridflag) {
    // Coulombic cutoff is set by user, choose grid to give requested error
    nx_max = ny_max = nz_max = 2;
    hx = xprd/nx_max;
    hy = yprd/ny_max;
    hz = zprd/nz_max;

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
    // cutoff and grid are set by user
    nx_max = nx_msm_max;
    ny_max = ny_msm_max;
    nz_max = nz_msm_max;
  }

  // scale grid for triclinic skew

  if (triclinic && !gridflag) {
    double tmp[3];
    tmp[0] = nx_max/xprd;
    tmp[1] = ny_max/yprd;
    tmp[2] = nz_max/zprd;
    lamda2xT(&tmp[0],&tmp[0]);
    nx_max = static_cast<int>(tmp[0]);
    ny_max = static_cast<int>(tmp[1]);
    nz_max = static_cast<int>(tmp[2]);
  }

  // boost grid size until it is factorable by 2
  // round up or down, depending on which is closer

  int flag = 0;
  int xlevels,ylevels,zlevels;

  while (!factorable(nx_max,flag,xlevels)) {
    double k = log(nx_max)/log(2.0);
    double r = k - floor(k);
    if (r > 0.5) nx_max++;
    else nx_max--;
  }
  while (!factorable(ny_max,flag,ylevels)) {
    double k = log(ny_max)/log(2.0);
    double r = k - floor(k);
    if (r > 0.5) ny_max++;
    else ny_max--;
  }
  while (!factorable(nz_max,flag,zlevels)) {
    double k = log(nz_max)/log(2.0);
    double r = k - floor(k);
    if (r > 0.5) nz_max++;
    else nz_max--;
  }

  if (flag && gridflag && me == 0)
    error->warning(FLERR,
                   "Number of MSM mesh points changed to be a multiple of 2");

  // adjust Coulombic cutoff to give desired error (if requested)

  if (adjust_cutoff_flag) {
    hx = xprd/nx_max;
    hy = yprd/ny_max;
    hz = zprd/nz_max;

    int p = order - 1;
    double Lx2 = xprd*xprd;
    double Ly2 = yprd*yprd;
    double Lz2 = zprd*zprd;
    double hx2pm2 = pow(hx,2.0*p-2.0);
    double hy2pm2 = pow(hy,2.0*p-2.0);
    double hz2pm2 = pow(hz,2.0*p-2.0);
    estimate_1d_error(1.0,1.0); // make sure that C_p is defined
    double k = q2*C_p/accuracy/sqrt(double(atom->natoms));
    double sum = hx2pm2/Lx2 + hy2pm2/Ly2 + hz2pm2/Lz2;

    cutoff = pow(k*k*sum/3.0,1.0/(2.0*p));
    int itmp;
    double *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
    *p_cutoff = cutoff;

    char str[128];
    sprintf(str,"Adjusting Coulombic cutoff for MSM, new cutoff = %g",cutoff);
    if (me == 0) error->warning(FLERR,str);
  }

  if (triclinic == 0) {
    h_x = xprd/nx_max;
    h_y = yprd/ny_max;
    h_z = zprd/nz_max;
  } else {
    double tmp[3];
    tmp[0] = nx_max;
    tmp[1] = ny_max;
    tmp[2] = nz_max;
    x2lamdaT(&tmp[0],&tmp[0]);
    h_x = 1.0/tmp[0];
    h_y = 1.0/tmp[1];
    h_z = 1.0/tmp[2];
  }

  // find maximum number of levels

  levels = MAX(xlevels,ylevels);
  levels = MAX(levels,zlevels);

  if (levels > MAX_LEVELS) error->all(FLERR,"Too many MSM grid levels");

  // need at least 2 MSM levels for periodic systems

  if (levels <= 1) {
    levels = xlevels = ylevels = zlevels = 2;
    nx_max = ny_max = nz_max = 2;
    if (gridflag)
      error->warning(FLERR,
             "MSM mesh too small, increasing to 2 points in each direction");
  }

  // omit top grid level for periodic systems

  if (!domain->nonperiodic) levels -= 1;

  deallocate_levels();
  allocate_levels();

  // find number of grid levels in each direction

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

  if (nx_msm[0] >= OFFSET || ny_msm[0] >= OFFSET || nz_msm[0] >= OFFSET)
    error->all(FLERR,"MSM grid is too large");

  // compute number of extra grid points needed for non-periodic boundary conditions

  if (domain->nonperiodic) {
    alpha[0] = -(order/2 - 1);
    betax[0] = nx_msm[0] + (order/2 - 1);
    betay[0] = ny_msm[0] + (order/2 - 1);
    betaz[0] = nz_msm[0] + (order/2 - 1);
    for (int n = 1; n < levels; n++) {
      alpha[n] = -((-alpha[n-1]+1)/2) - (order/2 - 1);
      betax[n] = ((betax[n-1]+1)/2) + (order/2 - 1);
      betay[n] = ((betay[n-1]+1)/2) + (order/2 - 1);
      betaz[n] = ((betaz[n-1]+1)/2) + (order/2 - 1);
    }
  }

  if (domain->nonperiodic) {
    alpha[0] = -(order/2 - 1);
    betax[0] = nx_msm[0] + (order/2 - 1);
    betay[0] = ny_msm[0] + (order/2 - 1);
    betaz[0] = nz_msm[0] + (order/2 - 1);
    for (int n = 1; n < levels; n++) {
      alpha[n] = -((-alpha[n-1]+1)/2) - (order/2 - 1);
      betax[n] = ((betax[n-1]+1)/2) + (order/2 - 1);
      betay[n] = ((betay[n-1]+1)/2) + (order/2 - 1);
      betaz[n] = ((betaz[n-1]+1)/2) + (order/2 - 1);
    }
  }

}

/* ----------------------------------------------------------------------
   set local subset of MSM grid that I own
   n xyz lo/hi in = 3d grid that I own (inclusive)
   n xyz lo/hi out = 3d grid + ghost cells in 6 directions (inclusive)
------------------------------------------------------------------------- */

void MSM::set_grid_local()
{
  // loop over grid levels

  for (int n=0; n<levels; n++) {

    // global indices of MSM grid range from 0 to N-1
    // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
    //   global MSM grid that I own without ghost cells

    nxlo_in[n] = static_cast<int> (comm->xsplit[comm->myloc[0]] * nx_msm[n]);
    nxhi_in[n] = static_cast<int> (comm->xsplit[comm->myloc[0]+1] * nx_msm[n]) - 1;

    nylo_in[n] = static_cast<int> (comm->ysplit[comm->myloc[1]] * ny_msm[n]);
    nyhi_in[n] = static_cast<int> (comm->ysplit[comm->myloc[1]+1] * ny_msm[n]) - 1;

    nzlo_in[n] = static_cast<int> (comm->zsplit[comm->myloc[2]] * nz_msm[n]);
    nzhi_in[n] = static_cast<int> (comm->zsplit[comm->myloc[2]+1] * nz_msm[n]) - 1;

    // nlower,nupper = stencil size for mapping (interpolating) particles to MSM grid

    nlower = -(order-1)/2;
    nupper = order/2;

    // lengths of box and processor sub-domains

    double *prd,*sublo,*subhi;

    if (!triclinic) {
      prd = domain->prd;
      sublo = domain->sublo;
      subhi = domain->subhi;
    } else {
      prd = domain->prd_lamda;
      sublo = domain->sublo_lamda;
      subhi = domain->subhi_lamda;
    }

    double xprd = prd[0];
    double yprd = prd[1];
    double zprd = prd[2];

    // shift values for particle <-> grid mapping
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    // nlo_out,nhi_out = lower/upper limits of the 3d sub-brick of
    //   global MSM grid that my particles can contribute charge to
    // effectively nlo_in,nhi_in + ghost cells
    // nlo,nhi = global coords of grid pt to "lower left" of smallest/largest
    //           position a particle in my box can be at
    // dist[3] = particle position bound = subbox + skin/2.0
    // nlo_out,nhi_out = nlo,nhi + stencil size for particle mapping

    double dist[3];
    double cuthalf = 0.0;
    if (n == 0) cuthalf = 0.5*neighbor->skin; // only applies to finest grid
    dist[0] = dist[1] = dist[2] = cuthalf;
    if (triclinic) kspacebbox(cuthalf,&dist[0]);

    int nlo,nhi;

    nlo = static_cast<int> ((sublo[0]-dist[0]-boxlo[0]) *
                            nx_msm[n]/xprd + OFFSET) - OFFSET;
    nhi = static_cast<int> ((subhi[0]+dist[0]-boxlo[0]) *
                            nx_msm[n]/xprd + OFFSET) - OFFSET;
    if (n == 0) {
      // use a smaller ghost region for interpolation
      nxlo_out_all = nlo + nlower;
      nxhi_out_all = nhi + nupper;
    }
    // a larger ghost region is needed for the direct sum and for restriction/prolongation
    nxlo_out[n] = nlo + MIN(-order,nxlo_direct);
    nxhi_out[n] = nhi + MAX(order,nxhi_direct);

    nlo = static_cast<int> ((sublo[1]-dist[1]-boxlo[1]) *
                            ny_msm[n]/yprd + OFFSET) - OFFSET;
    nhi = static_cast<int> ((subhi[1]+dist[1]-boxlo[1]) *
                            ny_msm[n]/yprd + OFFSET) - OFFSET;
    if (n == 0) {
      nylo_out_all = nlo + nlower;
      nyhi_out_all = nhi + nupper;
    }
    nylo_out[n] = nlo + MIN(-order,nylo_direct);
    nyhi_out[n] = nhi + MAX(order,nyhi_direct);

    nlo = static_cast<int> ((sublo[2]-dist[2]-boxlo[2]) *
                            nz_msm[n]/zprd + OFFSET) - OFFSET;
    nhi = static_cast<int> ((subhi[2]+dist[2]-boxlo[2]) *
                            nz_msm[n]/zprd + OFFSET) - OFFSET;
    if (n == 0) {
      nzlo_out_all = nlo + nlower;
      nzhi_out_all = nhi + nupper;
    }
    // a hemisphere is used for direct sum interactions,
    //   so no ghosting is needed for direct sum in the -z direction
    nzlo_out[n] = nlo - order;
    nzhi_out[n] = nhi + MAX(order,nzhi_direct);

    // add extra grid points for non-periodic boundary conditions

    if (domain->nonperiodic) {

      if (!domain->xperiodic) {
        if (nxlo_in[n] == 0)
          nxlo_in[n] = alpha[n];
        nxlo_out[n] = MAX(nxlo_out[n],alpha[n]);
        if (n == 0) nxlo_out_all = MAX(nxlo_out_all,alpha[0]);

        if (nxhi_in[n] == nx_msm[n] - 1)
          nxhi_in[n] = betax[n];
        nxhi_out[n] = MIN(nxhi_out[n],betax[n]);
        if (n == 0) nxhi_out_all = MIN(nxhi_out_all,betax[0]);
        if (nxhi_in[n] < 0)
          nxhi_in[n] = alpha[n] - 1;
      }

      if (!domain->yperiodic) {
        if (nylo_in[n] == 0)
          nylo_in[n] = alpha[n];
        nylo_out[n] = MAX(nylo_out[n],alpha[n]);
        if (n == 0) nylo_out_all = MAX(nylo_out_all,alpha[0]);

        if (nyhi_in[n] == ny_msm[n] - 1)
          nyhi_in[n] = betay[n];
        nyhi_out[n] = MIN(nyhi_out[n],betay[n]);
        if (n == 0) nyhi_out_all = MIN(nyhi_out_all,betay[0]);
        if (nyhi_in[n] < 0)
          nyhi_in[n] = alpha[n] - 1;
      }

      if (!domain->zperiodic) {
        if (nzlo_in[n] == 0)
          nzlo_in[n] = alpha[n];
        nzlo_out[n] = MAX(nzlo_out[n],alpha[n]);
        if (n == 0) nzlo_out_all = MAX(nzlo_out_all,alpha[0]);

        if (nzhi_in[n] == nz_msm[n] - 1)
          nzhi_in[n] = betaz[n];
        nzhi_out[n] = MIN(nzhi_out[n],betaz[n]);
        if (n == 0) nzhi_out_all = MIN(nzhi_out_all,betaz[0]);
        if (nzhi_in[n] < 0)
          nzhi_in[n] = alpha[n] - 1;
      }
    }

    // prevent inactive processors from participating in MPI communication routines

    set_proc_grid(n);

    // MSM grids for this proc, including ghosts

    ngrid[n] = (nxhi_out[n]-nxlo_out[n]+1) * (nyhi_out[n]-nylo_out[n]+1) *
      (nzhi_out[n]-nzlo_out[n]+1);
  }
}

/* ----------------------------------------------------------------------
   find active procs and neighbor procs for each grid level
------------------------------------------------------------------------- */

void MSM::set_proc_grid(int n)
{
    for (int i=0; i<3; i++)
      myloc[i] = comm->myloc[i];

  // size of inner MSM grid owned by this proc

  int nxgrid_in = nxhi_in[n]-nxlo_in[n]+1;
  int nygrid_in = nyhi_in[n]-nylo_in[n]+1;
  int nzgrid_in = nzhi_in[n]-nzlo_in[n]+1;
  int ngrid_in = nxgrid_in * nygrid_in * nzgrid_in;

  // check to see if this proc owns any inner grid points on this grid level
  //  if not, deactivate by setting active_flag = 0

  int cnt[3];

  cnt[0] = 0;
  if (myloc[1] == 0 && myloc[2] == 0)
    if (nxgrid_in > 0)
      cnt[0] = 1;

  cnt[1] = 0;
  if (myloc[0] == 0 && myloc[2] == 0)
    if (nygrid_in > 0)
      cnt[1] = 1;

  cnt[2] = 0;
  if (myloc[0] == 0 && myloc[1] == 0)
    if (nzgrid_in > 0)
      cnt[2] = 1;

  MPI_Allreduce(&cnt[0],&procgrid[0],3,MPI_INT,MPI_SUM,world);

  int color;

  if (ngrid_in > 0) {
    active_flag[n] = 1;
    color = 0;
  } else {
    active_flag[n] = 0;
    color = MPI_UNDEFINED;
  }

  // define a new MPI communicator for this grid level that only includes active procs

  MPI_Comm_split(world,color,me,&world_levels[n]);

  if (!active_flag[n]) return;

  int procneigh[3][2]; // my 6 neighboring procs, 0/1 = left/right

  // map processor IDs to new 3d processor grid
  // produces myloc, procneigh

  int periods[3];
  periods[0] = periods[1] = periods[2] = 1;
  MPI_Comm cartesian;

  MPI_Cart_create(world_levels[n],3,procgrid,periods,0,&cartesian);
  MPI_Cart_get(cartesian,3,procgrid,periods,myloc);
  MPI_Cart_shift(cartesian,0,1,&procneigh[0][0],&procneigh[0][1]);
  MPI_Cart_shift(cartesian,1,1,&procneigh[1][0],&procneigh[1][1]);
  MPI_Cart_shift(cartesian,2,1,&procneigh[2][0],&procneigh[2][1]);

  MPI_Comm_free(&cartesian);

  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      procneigh_levels[n][i][j] = procneigh[i][j];
}

/* ----------------------------------------------------------------------
   reset local grid arrays and communication stencils
   called by fix balance b/c it changed sizes of processor sub-domains
------------------------------------------------------------------------- */

void MSM::setup_grid()
{
  // free all arrays previously allocated
  // pre-compute volume-dependent coeffs
  // reset portion of global grid that each proc owns
  // reallocate MSM long-range dependent memory
  // don't invoke allocate_peratom(), compute() will allocate when needed

  setup();
}

/* ----------------------------------------------------------------------
   check if all factors of n are in list of factors
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int MSM::factorable(int n, int &flag, int &levels)
{
  int i;
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

  if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

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
   aninterpolation: interpolate charges from particles to grid
------------------------------------------------------------------------- */

void MSM::make_rho()
{
  //fprintf(screen,"MSM aninterpolation\n\n");

  int i,l,m,n,nx,ny,nz,mx,my,mz;
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

    compute_phis(dx,dy,dz);

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
   MSM direct sum procedure for intermediate grid levels, solve Poisson's
   equation to get energy, virial, etc.
------------------------------------------------------------------------- */

void MSM::direct(int n)
{
  //fprintf(screen,"Direct contribution on level %i\n\n",n);

  double ***qgridn = qgrid[n];
  double ***egridn = egrid[n];
  double ***v0gridn = v0grid[n];
  double ***v1gridn = v1grid[n];
  double ***v2gridn = v2grid[n];
  double ***v3gridn = v3grid[n];
  double ***v4gridn = v4grid[n];
  double ***v5gridn = v5grid[n];
  double *g_directn = g_direct[n];
  double *v0_directn = v0_direct[n];
  double *v1_directn = v1_direct[n];
  double *v2_directn = v2_direct[n];
  double *v3_directn = v3_direct[n];
  double *v4_directn = v4_direct[n];
  double *v5_directn = v5_direct[n];

  // zero out electric potential

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

  int icx,icy,icz,ix,iy,iz,zk,zyk,k;
  int ii,jj,kk;
  int imin,imax,jmin,jmax,kmax;
  double qtmp,qtmp2,gtmp;
  double esum,v0sum,v1sum,v2sum,v3sum,v4sum,v5sum;
  double **qk,**ek;
  double *qkj,*ekj;

  int nx = nxhi_direct - nxlo_direct + 1;
  int ny = nyhi_direct - nylo_direct + 1;

  // loop over inner grid points

  for (icz = nzlo_in[n]; icz <= nzhi_in[n]; icz++) {

    if (domain->zperiodic) {
      kmax = nzhi_direct;
    } else {
      kmax = MIN(nzhi_direct,betaz[n] - icz);
    }

    for (icy = nylo_in[n]; icy <= nyhi_in[n]; icy++) {

      if (domain->yperiodic) {
        jmin = nylo_direct;
        jmax = nyhi_direct;
      } else {
        jmin = MAX(nylo_direct,alpha[n] - icy);
        jmax = MIN(nyhi_direct,betay[n] - icy);
      }

      for (icx = nxlo_in[n]; icx <= nxhi_in[n]; icx++) {

        if (domain->xperiodic) {
          imin = nxlo_direct;
          imax = nxhi_direct;
        } else {
          imin = MAX(nxlo_direct,alpha[n] - icx);
          imax = MIN(nxhi_direct,betax[n] - icx);
        }

        qtmp = qgridn[icz][icy][icx]; // charge on center grid point

        esum = 0.0;
        if (vflag_either && !scalar_pressure_flag)
          v0sum = v1sum = v2sum = v3sum = v4sum = v5sum = 0.0;

        // use hemisphere to avoid double computation of pair-wise
        //   interactions in direct sum (no computations in -z direction)

        for (iz = 1; iz <= kmax; iz++) {
          kk = icz+iz;
          qk = qgridn[kk];
          ek = egridn[kk];
          zk = (iz + nzhi_direct)*ny;
          for (iy = jmin; iy <= jmax; iy++) {
            jj = icy+iy;
            qkj = qk[jj];
            ekj = ek[jj];
            zyk = (zk + iy + nyhi_direct)*nx;
            for (ix = imin; ix <= imax; ix++) {
              ii = icx+ix;
              qtmp2 = qkj[ii]; // charge on outer grid point
              k = zyk + ix + nxhi_direct;
              gtmp = g_directn[k];
              esum += gtmp * qtmp2;
              ekj[ii] += gtmp * qtmp;

              if (vflag_either && !scalar_pressure_flag) {
                v0sum += v0_directn[k] * qtmp2;
                v1sum += v1_directn[k] * qtmp2;
                v2sum += v2_directn[k] * qtmp2;
                v3sum += v3_directn[k] * qtmp2;
                v4sum += v4_directn[k] * qtmp2;
                v5sum += v5_directn[k] * qtmp2;
              }
            }
          }
        }

        // iz=0

        iz = 0;
        kk = icz+iz;
        qk = qgridn[kk];
        ek = egridn[kk];
        zk = (iz + nzhi_direct)*ny;
        for (iy = 1; iy <= jmax; iy++) {
          jj = icy+iy;
          qkj = qk[jj];
          ekj = ek[jj];
          zyk = (zk + iy + nyhi_direct)*nx;
          for (ix = imin; ix <= imax; ix++) {
            ii = icx+ix;
            qtmp2 = qkj[ii];
            k = zyk + ix + nxhi_direct;
            gtmp = g_directn[k];
            esum += gtmp * qtmp2;
            ekj[ii] += gtmp * qtmp;

            if (vflag_either && !scalar_pressure_flag) {
              v0sum += v0_directn[k] * qtmp2;
              v1sum += v1_directn[k] * qtmp2;
              v2sum += v2_directn[k] * qtmp2;
              v3sum += v3_directn[k] * qtmp2;
              v4sum += v4_directn[k] * qtmp2;
              v5sum += v5_directn[k] * qtmp2;
            }
          }
        }

        // iz=0, iy=0

        iz = 0;
        kk = icz+iz;
        qk = qgridn[kk];
        ek = egridn[kk];
        zk = (iz + nzhi_direct)*ny;
        iy = 0;
        jj = icy+iy;
        qkj = qk[jj];
        ekj = ek[jj];
        zyk = (zk + iy + nyhi_direct)*nx;
        for (ix = 1; ix <= imax; ix++) {
          ii = icx+ix;
          qtmp2 = qkj[ii];
          k = zyk + ix + nxhi_direct;
          gtmp = g_directn[k];
          esum += gtmp * qtmp2;
          ekj[ii] += gtmp * qtmp;

          if (vflag_either && !scalar_pressure_flag) {
            v0sum += v0_directn[k] * qtmp2;
            v1sum += v1_directn[k] * qtmp2;
            v2sum += v2_directn[k] * qtmp2;
            v3sum += v3_directn[k] * qtmp2;
            v4sum += v4_directn[k] * qtmp2;
            v5sum += v5_directn[k] * qtmp2;
          }
        }

        // iz=0, iy=0, ix=0

        iz = 0;
        zk = (iz + nzhi_direct)*ny;
        iy = 0;
        zyk = (zk + iy + nyhi_direct)*nx;
        ix = 0;
        k = zyk + ix + nxhi_direct;
        gtmp = g_directn[k];
        esum += 0.5 * gtmp * qtmp;
        egridn[icz][icy][icx] += 0.5 * gtmp * qtmp;

        // virial is zero for iz=0, iy=0, ix=0

        // accumulate per-atom energy/virial

        egridn[icz][icy][icx] += esum;

        if (vflag_atom && !scalar_pressure_flag) {
          v0gridn[icz][icy][icx] += v0sum;
          v1gridn[icz][icy][icx] += v1sum;
          v2gridn[icz][icy][icx] += v2sum;
          v3gridn[icz][icy][icx] += v3sum;
          v4gridn[icz][icy][icx] += v4sum;
          v5gridn[icz][icy][icx] += v5sum;
        }

        // accumulate total energy/virial

        if (evflag) {
          qtmp = qgridn[icz][icy][icx];
          if (eflag_global) energy += 2.0 * esum * qtmp;
          if (vflag_global && !scalar_pressure_flag) {
            virial[0] += 2.0 * v0sum * qtmp;
            virial[1] += 2.0 * v1sum * qtmp;
            virial[2] += 2.0 * v2sum * qtmp;
            virial[3] += 2.0 * v3sum * qtmp;
            virial[4] += 2.0 * v4sum * qtmp;
            virial[5] += 2.0 * v5sum * qtmp;
          }
        }

      }
    }
  }

  // compute per-atom virial (if requested)

  if (vflag_atom)
    direct_peratom(n);
}

/* ----------------------------------------------------------------------
   MSM direct sum procedure for intermediate grid levels, solve Poisson's
   equation to get per-atom virial, separate method used for performance
   reasons
------------------------------------------------------------------------- */

void MSM::direct_peratom(int n)
{
  //fprintf(screen,"Direct contribution on level %i\n\n",n);

  double ***qgridn = qgrid[n];
  double ***v0gridn = v0grid[n];
  double ***v1gridn = v1grid[n];
  double ***v2gridn = v2grid[n];
  double ***v3gridn = v3grid[n];
  double ***v4gridn = v4grid[n];
  double ***v5gridn = v5grid[n];

  int icx,icy,icz,ix,iy,iz,zk,zyk,k;
  int ii,jj,kk;
  int imin,imax,jmin,jmax,kmax;
  double qtmp;

  int nx = nxhi_direct - nxlo_direct + 1;
  int ny = nyhi_direct - nylo_direct + 1;

  // loop over inner grid points

  for (icz = nzlo_in[n]; icz <= nzhi_in[n]; icz++) {

    if (domain->zperiodic) {
      kmax = nzhi_direct;
    } else {
      kmax = MIN(nzhi_direct,betaz[n] - icz);
    }

    for (icy = nylo_in[n]; icy <= nyhi_in[n]; icy++) {

      if (domain->yperiodic) {
        jmin = nylo_direct;
        jmax = nyhi_direct;
      } else {
        jmin = MAX(nylo_direct,alpha[n] - icy);
        jmax = MIN(nyhi_direct,betay[n] - icy);
      }

      for (icx = nxlo_in[n]; icx <= nxhi_in[n]; icx++) {

        if (domain->xperiodic) {
          imin = nxlo_direct;
          imax = nxhi_direct;
        } else {
          imin = MAX(nxlo_direct,alpha[n] - icx);
          imax = MIN(nxhi_direct,betax[n] - icx);
        }

        qtmp = qgridn[icz][icy][icx]; // center grid point

        // use hemisphere to avoid double computation of pair-wise
        //   interactions in direct sum (no computations in -z direction)

        for (iz = 1; iz <= kmax; iz++) {
          kk = icz+iz;
          zk = (iz + nzhi_direct)*ny;
          for (iy = jmin; iy <= jmax; iy++) {
            jj = icy+iy;
            zyk = (zk + iy + nyhi_direct)*nx;
            for (ix = imin; ix <= imax; ix++) {
              ii = icx+ix;
              k = zyk + ix + nxhi_direct;
              v0gridn[kk][jj][ii] += v0_direct[n][k] * qtmp;
              v1gridn[kk][jj][ii] += v1_direct[n][k] * qtmp;
              v2gridn[kk][jj][ii] += v2_direct[n][k] * qtmp;
              v3gridn[kk][jj][ii] += v3_direct[n][k] * qtmp;
              v4gridn[kk][jj][ii] += v4_direct[n][k] * qtmp;
              v5gridn[kk][jj][ii] += v5_direct[n][k] * qtmp;
            }
          }
        }

        // iz=0

        iz = 0;
        kk = icz+iz;
        zk = (iz + nzhi_direct)*ny;
        for (iy = 1; iy <= jmax; iy++) {
          jj = icy+iy;
          zyk = (zk + iy + nyhi_direct)*nx;
          for (ix = imin; ix <= imax; ix++) {
            ii = icx+ix;
            k = zyk + ix + nxhi_direct;
            v0gridn[kk][jj][ii] += v0_direct[n][k] * qtmp;
            v1gridn[kk][jj][ii] += v1_direct[n][k] * qtmp;
            v2gridn[kk][jj][ii] += v2_direct[n][k] * qtmp;
            v3gridn[kk][jj][ii] += v3_direct[n][k] * qtmp;
            v4gridn[kk][jj][ii] += v4_direct[n][k] * qtmp;
            v5gridn[kk][jj][ii] += v5_direct[n][k] * qtmp;
          }
        }

        // iz=0, iy=0

        iz = 0;
        kk = icz+iz;
        zk = (iz + nzhi_direct)*ny;
        iy = 0;
        jj = icy+iy;
        zyk = (zk + iy + nyhi_direct)*nx;
        for (ix = 1; ix <= imax; ix++) {
          ii = icx+ix;
          k = zyk + ix + nxhi_direct;
          v0gridn[kk][jj][ii] += v0_direct[n][k] * qtmp;
          v1gridn[kk][jj][ii] += v1_direct[n][k] * qtmp;
          v2gridn[kk][jj][ii] += v2_direct[n][k] * qtmp;
          v3gridn[kk][jj][ii] += v3_direct[n][k] * qtmp;
          v4gridn[kk][jj][ii] += v4_direct[n][k] * qtmp;
          v5gridn[kk][jj][ii] += v5_direct[n][k] * qtmp;
       }

        // virial is zero for iz=0, iy=0, ix=0

      }
    }
  }
}

/* ----------------------------------------------------------------------
   MSM direct sum procedure for top grid level (nonperiodic systems only)
------------------------------------------------------------------------- */

void MSM::direct_top(int n)
{
  //fprintf(screen,"Direct contribution on level %i\n\n",n);

  double ***qgridn = qgrid[n];
  double ***egridn = egrid[n];
  double ***v0gridn = v0grid[n];
  double ***v1gridn = v1grid[n];
  double ***v2gridn = v2grid[n];
  double ***v3gridn = v3grid[n];
  double ***v4gridn = v4grid[n];
  double ***v5gridn = v5grid[n];

  // zero out electric potential

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

  int icx,icy,icz,ix,iy,iz,zk,zyk,k;
  int ii,jj,kk;
  int imin,imax,jmin,jmax,kmax;
  double qtmp,qtmp2,gtmp;
  double esum,v0sum,v1sum,v2sum,v3sum,v4sum,v5sum;
  double **qk,**ek;
  double *qkj,*ekj;

  int nx_top = betax[n] - alpha[n];
  int ny_top = betay[n] - alpha[n];
  int nz_top = betaz[n] - alpha[n];

  int nx = 2*nx_top + 1;
  int ny = 2*ny_top + 1;

  // loop over inner grid points

  for (icz = nzlo_in[n]; icz <= nzhi_in[n]; icz++) {

    if (domain->zperiodic) {
      kmax = nz_msm[n]-1;
    } else {
      kmax = betaz[n] - icz;
    }

    for (icy = nylo_in[n]; icy <= nyhi_in[n]; icy++) {

      if (domain->yperiodic) {
        jmin = 0;
        jmax = ny_msm[n]-1;
      } else {
        jmin = alpha[n] - icy;
        jmax = betay[n] - icy;
      }

      for (icx = nxlo_in[n]; icx <= nxhi_in[n]; icx++) {

        if (domain->xperiodic) {
          imin = 0;
          imax = nx_msm[n]-1;
        } else {
          imin = alpha[n] - icx;
          imax = betax[n] - icx;
        }

        qtmp = qgridn[icz][icy][icx];

        esum = 0.0;
        if (vflag_either && !scalar_pressure_flag)
          v0sum = v1sum = v2sum = v3sum = v4sum = v5sum = 0.0;

        // use hemisphere to avoid double computation of pair-wise
        //   interactions in direct sum (no computations in -z direction)

        for (iz = 1; iz <= kmax; iz++) {
          kk = icz+iz;
          qk = qgridn[kk];
          ek = egridn[kk];
          zk = (iz + nz_top)*ny;
          for (iy = jmin; iy <= jmax; iy++) {
            jj = icy+iy;
            qkj = qk[jj];
            ekj = ek[jj];
            zyk = (zk + iy + ny_top)*nx;
            for (ix = imin; ix <= imax; ix++) {
              ii = icx+ix;
              qtmp2 = qkj[ii];
              k = zyk + ix + nx_top;
              gtmp = g_direct_top[k];
              esum += gtmp * qtmp2;
              ekj[ii] += gtmp * qtmp;

              if (vflag_either && !scalar_pressure_flag) {
                v0sum += v0_direct_top[k] * qtmp2;
                v1sum += v1_direct_top[k] * qtmp2;
                v2sum += v2_direct_top[k] * qtmp2;
                v3sum += v3_direct_top[k] * qtmp2;
                v4sum += v4_direct_top[k] * qtmp2;
                v5sum += v5_direct_top[k] * qtmp2;
              }
            }
          }
        }

        // iz=0

        iz = 0;
        kk = icz+iz;
        qk = qgridn[kk];
        ek = egridn[kk];
        zk = (iz + nz_top)*ny;
        for (iy = 1; iy <= jmax; iy++) {
          jj = icy+iy;
          qkj = qk[jj];
          ekj = ek[jj];
          zyk = (zk + iy + ny_top)*nx;
          for (ix = imin; ix <= imax; ix++) {
            ii = icx+ix;
            qtmp2 = qkj[ii];
            k = zyk + ix + nx_top;
            gtmp = g_direct_top[k];
            esum += gtmp * qtmp2;
            ekj[ii] += gtmp * qtmp;

            if (vflag_either && !scalar_pressure_flag) {
              v0sum += v0_direct_top[k] * qtmp2;
              v1sum += v1_direct_top[k] * qtmp2;
              v2sum += v2_direct_top[k] * qtmp2;
              v3sum += v3_direct_top[k] * qtmp2;
              v4sum += v4_direct_top[k] * qtmp2;
              v5sum += v5_direct_top[k] * qtmp2;
            }
          }
        }

        // iz=0, iy=0

        iz = 0;
        kk = icz+iz;
        qk = qgridn[kk];
        ek = egridn[kk];
        zk = (iz + nz_top)*ny;
        iy = 0;
        jj = icy+iy;
        qkj = qk[jj];
        ekj = ek[jj];
        zyk = (zk + iy + ny_top)*nx;
        for (ix = 1; ix <= imax; ix++) {
          ii = icx+ix;
          qtmp2 = qkj[ii];
          k = zyk + ix + nx_top;
          gtmp = g_direct_top[k];
          esum += gtmp * qtmp2;
          ekj[ii] += gtmp * qtmp;

          if (vflag_either && !scalar_pressure_flag) {
            v0sum += v0_direct_top[k] * qtmp2;
            v1sum += v1_direct_top[k] * qtmp2;
            v2sum += v2_direct_top[k] * qtmp2;
            v3sum += v3_direct_top[k] * qtmp2;
            v4sum += v4_direct_top[k] * qtmp2;
            v5sum += v5_direct_top[k] * qtmp2;
          }
        }

        // iz=0, iy=0, ix=0

        iz = 0;
        zk = (iz + nz_top)*ny;
        iy = 0;
        zyk = (zk + iy + ny_top)*nx;
        ix = 0;
        ii = icx+ix;
        k = zyk + ix + nx_top;
        gtmp = g_direct_top[k];
        esum += 0.5 * gtmp * qtmp;
        egridn[icz][icy][icx] += 0.5 * gtmp * qtmp;

        if (vflag_either && !scalar_pressure_flag) {
          v0sum += v0_direct_top[k] * qtmp;
          v1sum += v1_direct_top[k] * qtmp;
          v2sum += v2_direct_top[k] * qtmp;
          v3sum += v3_direct_top[k] * qtmp;
          v4sum += v4_direct_top[k] * qtmp;
          v5sum += v5_direct_top[k] * qtmp;
        }

        // accumulate per-atom energy/virial

        egridn[icz][icy][icx] += esum;

        if (vflag_atom && !scalar_pressure_flag) {
          v0gridn[icz][icy][icx] += v0sum;
          v1gridn[icz][icy][icx] += v1sum;
          v2gridn[icz][icy][icx] += v2sum;
          v3gridn[icz][icy][icx] += v3sum;
          v4gridn[icz][icy][icx] += v4sum;
          v5gridn[icz][icy][icx] += v5sum;
        }

        // accumulate total energy/virial

        if (evflag) {
          qtmp = qgridn[icz][icy][icx];
          if (eflag_global) energy += 2.0 * esum * qtmp;
          if (vflag_global && !scalar_pressure_flag) {
            virial[0] += 2.0 * v0sum * qtmp;
            virial[1] += 2.0 * v1sum * qtmp;
            virial[2] += 2.0 * v2sum * qtmp;
            virial[3] += 2.0 * v3sum * qtmp;
            virial[4] += 2.0 * v4sum * qtmp;
            virial[5] += 2.0 * v5sum * qtmp;
          }
        }

      }
    }
  }

  // compute per-atom virial (if requested)

  if (vflag_atom)
    direct_peratom_top(n);
}
/* ----------------------------------------------------------------------
   MSM direct sum procedure for top grid level, solve Poisson's
   equation to get per-atom virial, separate method used for performance
   reasons
------------------------------------------------------------------------- */

void MSM::direct_peratom_top(int n)
{
  double ***qgridn = qgrid[n];
  double ***v0gridn = v0grid[n];
  double ***v1gridn = v1grid[n];
  double ***v2gridn = v2grid[n];
  double ***v3gridn = v3grid[n];
  double ***v4gridn = v4grid[n];
  double ***v5gridn = v5grid[n];

  int icx,icy,icz,ix,iy,iz,zk,zyk,k;
  int ii,jj,kk;
  int imin,imax,jmin,jmax,kmax;
  double qtmp;

  int nx_top = betax[n] - alpha[n];
  int ny_top = betay[n] - alpha[n];
  int nz_top = betaz[n] - alpha[n];

  int nx = 2*nx_top + 1;
  int ny = 2*ny_top + 1;

  // loop over inner grid points

  for (icz = nzlo_in[n]; icz <= nzhi_in[n]; icz++) {

    if (domain->zperiodic) {
      kmax = nz_msm[n]-1;
    } else {
      kmax = betaz[n] - icz;
    }

    for (icy = nylo_in[n]; icy <= nyhi_in[n]; icy++) {

      if (domain->yperiodic) {
        jmin = 0;
        jmax = ny_msm[n]-1;
      } else {
        jmin = alpha[n] - icy;
        jmax = betay[n] - icy;
      }

      for (icx = nxlo_in[n]; icx <= nxhi_in[n]; icx++) {

        if (domain->xperiodic) {
          imin = 0;
          imax = nx_msm[n]-1;
        } else {
          imin = alpha[n] - icx;
          imax = betax[n] - icx;
        }

        qtmp = qgridn[icz][icy][icx]; // center grid point

        // use hemisphere to avoid double computation of pair-wise
        //   interactions in direct sum (no computations in -z direction)

        for (iz = 1; iz <= kmax; iz++) {
          kk = icz+iz;
          zk = (iz + nz_top)*ny;
          for (iy = jmin; iy <= jmax; iy++) {
            jj = icy+iy;
            zyk = (zk + iy + ny_top)*nx;
            for (ix = imin; ix <= imax; ix++) {
              ii = icx+ix;
              k = zyk + ix + nx_top;
              v0gridn[kk][jj][ii] += v0_direct_top[k] * qtmp;
              v1gridn[kk][jj][ii] += v1_direct_top[k] * qtmp;
              v2gridn[kk][jj][ii] += v2_direct_top[k] * qtmp;
              v3gridn[kk][jj][ii] += v3_direct_top[k] * qtmp;
              v4gridn[kk][jj][ii] += v4_direct_top[k] * qtmp;
              v5gridn[kk][jj][ii] += v5_direct_top[k] * qtmp;
            }
          }
        }

        // iz=0

        iz = 0;
        kk = icz+iz;
        zk = (iz + nz_top)*ny;
        for (iy = 1; iy <= jmax; iy++) {
          jj = icy+iy;
          zyk = (zk + iy + ny_top)*nx;
          for (ix = imin; ix <= imax; ix++) {
            ii = icx+ix;
            k = zyk + ix + nx_top;
            v0gridn[kk][jj][ii] += v0_direct_top[k] * qtmp;
            v1gridn[kk][jj][ii] += v1_direct_top[k] * qtmp;
            v2gridn[kk][jj][ii] += v2_direct_top[k] * qtmp;
            v3gridn[kk][jj][ii] += v3_direct_top[k] * qtmp;
            v4gridn[kk][jj][ii] += v4_direct_top[k] * qtmp;
            v5gridn[kk][jj][ii] += v5_direct_top[k] * qtmp;
          }
        }

        // iz=0, iy=0

        iz = 0;
        kk = icz+iz;
        zk = (iz + nz_top)*ny;
        iy = 0;
        jj = icy+iy;
        zyk = (zk + iy + ny_top)*nx;
        for (ix = 1; ix <= imax; ix++) {
          ii = icx+ix;
          k = zyk + ix + nx_top;
          v0gridn[kk][jj][ii] += v0_direct_top[k] * qtmp;
          v1gridn[kk][jj][ii] += v1_direct_top[k] * qtmp;
          v2gridn[kk][jj][ii] += v2_direct_top[k] * qtmp;
          v3gridn[kk][jj][ii] += v3_direct_top[k] * qtmp;
          v4gridn[kk][jj][ii] += v4_direct_top[k] * qtmp;
          v5gridn[kk][jj][ii] += v5_direct_top[k] * qtmp;
       }

        // virial is zero for iz=0, iy=0, ix=0

      }
    }
  }
}

/* ----------------------------------------------------------------------
   MSM restriction procedure for intermediate grid levels, interpolate
   charges from finer grid to coarser grid
------------------------------------------------------------------------- */

void MSM::restriction(int n)
{
  //fprintf(screen,"Restricting from level %i to %i\n\n",n,n+1);

  const int p = order-1;

  double ***qgrid1 = qgrid[n];
  double ***qgrid2 = qgrid[n+1];

  int k = 0;
  int *index = new int[p+2];
  for (int nu=-p; nu<=p; nu++) {
    if (nu%2 == 0 && nu != 0) continue;
    phi1d[0][k] = compute_phi(nu*delxinv[n+1]/delxinv[n]);
    phi1d[1][k] = compute_phi(nu*delyinv[n+1]/delyinv[n]);
    phi1d[2][k] = compute_phi(nu*delzinv[n+1]/delzinv[n]);
    index[k] = nu;
    k++;
  }

  int ip,jp,kp,ic,jc,kc,i,j;
  int ii,jj,kk;
  double phiz,phizy,q2sum;

  // zero out charge on coarser grid

  memset(&(qgrid2[nzlo_out[n+1]][nylo_out[n+1]][nxlo_out[n+1]]),0,
         ngrid[n+1]*sizeof(double));

  for (kp = nzlo_in[n+1]; kp <= nzhi_in[n+1]; kp++)
    for (jp = nylo_in[n+1]; jp <= nyhi_in[n+1]; jp++)
      for (ip = nxlo_in[n+1]; ip <= nxhi_in[n+1]; ip++) {

        ic = ip * static_cast<int> (delxinv[n]/delxinv[n+1]);
        jc = jp * static_cast<int> (delyinv[n]/delyinv[n+1]);
        kc = kp * static_cast<int> (delzinv[n]/delzinv[n+1]);

        q2sum = 0.0;

        for (k=0; k<=p+1; k++) {
          kk = kc+index[k];
          if (!domain->zperiodic) {
            if (kk < alpha[n]) continue;
            if (kk > betaz[n]) break;
          }
          phiz = phi1d[2][k];
          for (j=0; j<=p+1; j++) {
            jj = jc+index[j];
            if (!domain->yperiodic) {
              if (jj < alpha[n]) continue;
              if (jj > betay[n]) break;
            }
            phizy = phi1d[1][j]*phiz;
            for (i=0; i<=p+1; i++) {
              ii = ic+index[i];
              if (!domain->xperiodic) {
                if (ii < alpha[n]) continue;
                if (ii > betax[n]) break;
              }
              q2sum += qgrid1[kk][jj][ii] *
                phi1d[0][i]*phizy;
            }
          }
        }
        qgrid2[kp][jp][ip] += q2sum;
      }
  delete[] index;
}

/* ----------------------------------------------------------------------
   MSM prolongation procedure for intermediate grid levels, interpolate
   per-atom energy/virial from coarser grid to finer grid
------------------------------------------------------------------------- */

void MSM::prolongation(int n)
{
  //fprintf(screen,"Prolongating from level %i to %i\n\n",n+1,n);

  const int p = order-1;

  double ***egrid1 = egrid[n];
  double ***egrid2 = egrid[n+1];

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

  int k = 0;
  int *index = new int[p+2];
  for (int nu=-p; nu<=p; nu++) {
    if (nu%2 == 0 && nu != 0) continue;
    phi1d[0][k] = compute_phi(nu*delxinv[n+1]/delxinv[n]);
    phi1d[1][k] = compute_phi(nu*delyinv[n+1]/delyinv[n]);
    phi1d[2][k] = compute_phi(nu*delzinv[n+1]/delzinv[n]);
    index[k] = nu;
    k++;
  }

  int ip,jp,kp,ic,jc,kc,i,j;
  int ii,jj,kk;
  double phiz,phizy,phi3d;
  double etmp2,v0tmp2,v1tmp2,v2tmp2,v3tmp2,v4tmp2,v5tmp2;

  for (kp = nzlo_in[n+1]; kp <= nzhi_in[n+1]; kp++)
    for (jp = nylo_in[n+1]; jp <= nyhi_in[n+1]; jp++)
      for (ip = nxlo_in[n+1]; ip <= nxhi_in[n+1]; ip++) {

        ic = ip * static_cast<int> (delxinv[n]/delxinv[n+1]);
        jc = jp * static_cast<int> (delyinv[n]/delyinv[n+1]);
        kc = kp * static_cast<int> (delzinv[n]/delzinv[n+1]);

        etmp2 = egrid2[kp][jp][ip];

        if (vflag_atom) {
          v0tmp2 = v0grid2[kp][jp][ip];
          v1tmp2 = v1grid2[kp][jp][ip];
          v2tmp2 = v2grid2[kp][jp][ip];
          v3tmp2 = v3grid2[kp][jp][ip];
          v4tmp2 = v4grid2[kp][jp][ip];
          v5tmp2 = v5grid2[kp][jp][ip];
        }

        for (k=0; k<=p+1; k++) {
          kk = kc+index[k];
          if (!domain->zperiodic) {
            if (kk < alpha[n]) continue;
            if (kk > betaz[n]) break;
          }
          phiz = phi1d[2][k];
          for (j=0; j<=p+1; j++) {
            jj = jc+index[j];
            if (!domain->yperiodic) {
              if (jj < alpha[n]) continue;
              if (jj > betay[n]) break;
            }
            phizy = phi1d[1][j]*phiz;
            for (i=0; i<=p+1; i++) {
              ii = ic+index[i];
              if (!domain->xperiodic) {
                if (ii < alpha[n]) continue;
                if (ii > betax[n]) break;
              }
              phi3d = phi1d[0][i]*phizy;

              egrid1[kk][jj][ii] += etmp2 * phi3d;

              if (vflag_atom) {
                v0grid1[kk][jj][ii] += v0tmp2 * phi3d;
                v1grid1[kk][jj][ii] += v1tmp2 * phi3d;
                v2grid1[kk][jj][ii] += v2tmp2 * phi3d;
                v3grid1[kk][jj][ii] += v3tmp2 * phi3d;
                v4grid1[kk][jj][ii] += v4tmp2 * phi3d;
                v5grid1[kk][jj][ii] += v5tmp2 * phi3d;
              }

            }
          }
        }

      }
  delete[] index;
}

/* ----------------------------------------------------------------------
   Use MPI_Allreduce to fill ghost grid values, for coarse grids this may
   be cheaper than using nearest-neighbor communication (commgrid), right
   now only works for periodic boundary conditions
------------------------------------------------------------------------- */
void MSM::grid_swap_forward(int n, double*** &gridn)
{
  double ***gridn_tmp;
  memory->create(gridn_tmp,nz_msm[n],ny_msm[n],nx_msm[n],"msm:grid_tmp");

  double ***gridn_all;
  memory->create(gridn_all,nz_msm[n],ny_msm[n],nx_msm[n],"msm:grid_all");

  int ngrid_in = nx_msm[n] * ny_msm[n] * nz_msm[n];

  memset(&(gridn_tmp[0][0][0]),0,ngrid_in*sizeof(double));
  memset(&(gridn_all[0][0][0]),0,ngrid_in*sizeof(double));

  // copy inner grid cell values from gridn into gridn_tmp

  int icx,icy,icz;

  for (icz = nzlo_in[n]; icz <= nzhi_in[n]; icz++)
    for (icy = nylo_in[n]; icy <= nyhi_in[n]; icy++)
      for (icx = nxlo_in[n]; icx <= nxhi_in[n]; icx++)
        gridn_tmp[icz][icy][icx] = gridn[icz][icy][icx];

  MPI_Allreduce(&(gridn_tmp[0][0][0]),
                &(gridn_all[0][0][0]),
                ngrid_in,MPI_DOUBLE,MPI_SUM,world_levels[n]);

  // bitmask for PBCs (only works for power of 2 numbers)

  int PBCx,PBCy,PBCz;

  PBCx = nx_msm[n]-1;
  PBCy = ny_msm[n]-1;
  PBCz = nz_msm[n]-1;

  // copy from gridn_all into gridn to fill ghost grid cell values

  for (icz = nzlo_out[n]; icz <= nzhi_out[n]; icz++)
    for (icy = nylo_out[n]; icy <= nyhi_out[n]; icy++)
      for (icx = nxlo_out[n]; icx <= nxhi_out[n]; icx++)
        gridn[icz][icy][icx] = gridn_all[icz&PBCz][icy&PBCy][icx&PBCx];

  memory->destroy(gridn_tmp);
  memory->destroy(gridn_all);
}

/* ----------------------------------------------------------------------
   Use MPI_Allreduce to get contribution from ghost grid cells, for coarse
   grids this may be cheaper than using nearest-neighbor communication
   (commgrid), right now only works for periodic boundary conditions
------------------------------------------------------------------------- */

void MSM::grid_swap_reverse(int n, double*** &gridn)
{
  double ***gridn_tmp;
  memory->create(gridn_tmp,nz_msm[n],ny_msm[n],nx_msm[n],"msm:grid_tmp");

  double ***gridn_all;
  memory->create(gridn_all,nz_msm[n],ny_msm[n],nx_msm[n],"msm:grid_all");

  int ngrid_in = nx_msm[n] * ny_msm[n] * nz_msm[n];

  memset(&(gridn_tmp[0][0][0]),0,ngrid_in*sizeof(double));
  memset(&(gridn_all[0][0][0]),0,ngrid_in*sizeof(double));

  // bitmask for PBCs (only works for power of 2 numbers)

  int icx,icy,icz;
  int PBCx,PBCy,PBCz;

  PBCx = nx_msm[n]-1;
  PBCy = ny_msm[n]-1;
  PBCz = nz_msm[n]-1;

  // copy ghost grid cell values from gridn into inner portion of gridn_tmp

  for (icz = nzlo_out[n]; icz <= nzhi_out[n]; icz++)
    for (icy = nylo_out[n]; icy <= nyhi_out[n]; icy++)
      for (icx = nxlo_out[n]; icx <= nxhi_out[n]; icx++)
        gridn_tmp[icz&PBCz][icy&PBCy][icx&PBCx] += gridn[icz][icy][icx];

  MPI_Allreduce(&(gridn_tmp[0][0][0]),
                &(gridn_all[0][0][0]),
                ngrid_in,MPI_DOUBLE,MPI_SUM,world_levels[n]);

  // copy inner grid cell values from gridn_all into gridn

  for (icz = nzlo_in[n]; icz <= nzhi_in[n]; icz++)
    for (icy = nylo_in[n]; icy <= nyhi_in[n]; icy++)
      for (icx = nxlo_in[n]; icx <= nxhi_in[n]; icx++)
        gridn[icz][icy][icx] = gridn_all[icz][icy][icx];

  memory->destroy(gridn_tmp);
  memory->destroy(gridn_all);
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc (used by commgrid)
------------------------------------------------------------------------- */

void MSM::pack_forward(int flag, double *buf, int nlist, int *list)
{
  int n = current_level;

  double ***qgridn = qgrid[n];
  double ***egridn = egrid[n];

  double ***v0gridn = v0grid[n];
  double ***v1gridn = v1grid[n];
  double ***v2gridn = v2grid[n];
  double ***v3gridn = v3grid[n];
  double ***v4gridn = v4grid[n];
  double ***v5gridn = v5grid[n];

  int k = 0;

  if (flag == FORWARD_RHO) {
    double *qsrc = &qgridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++) {
      buf[k++] = qsrc[list[i]];
    }
  } else if (flag == FORWARD_AD) {
    double *src = &egridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == FORWARD_AD_PERATOM) {
    double *v0src = &v0gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v1src = &v1gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v2src = &v2gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v3src = &v3gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v4src = &v4gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v5src = &v5gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++) {
      buf[k++] = v0src[list[i]];
      buf[k++] = v1src[list[i]];
      buf[k++] = v2src[list[i]];
      buf[k++] = v3src[list[i]];
      buf[k++] = v4src[list[i]];
      buf[k++] = v5src[list[i]];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void MSM::unpack_forward(int flag, double *buf, int nlist, int *list)
{
  int n = current_level;

  double ***qgridn = qgrid[n];
  double ***egridn = egrid[n];

  double ***v0gridn = v0grid[n];
  double ***v1gridn = v1grid[n];
  double ***v2gridn = v2grid[n];
  double ***v3gridn = v3grid[n];
  double ***v4gridn = v4grid[n];
  double ***v5gridn = v5grid[n];

  int k = 0;

  if (flag == FORWARD_RHO) {
    double *dest = &qgridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++) {
      dest[list[i]] = buf[k++];
    }
  } else if (flag == FORWARD_AD) {
    double *dest = &egridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[k++];
  } else if (flag == FORWARD_AD_PERATOM) {
    double *v0src = &v0gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v1src = &v1gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v2src = &v2gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v3src = &v3gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v4src = &v4gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v5src = &v5gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++) {
      v0src[list[i]] = buf[k++];
      v1src[list[i]] = buf[k++];
      v2src[list[i]] = buf[k++];
      v3src[list[i]] = buf[k++];
      v4src[list[i]] = buf[k++];
      v5src[list[i]] = buf[k++];
    }
  }
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void MSM::pack_reverse(int flag, double *buf, int nlist, int *list)
{
  int n = current_level;

  double ***qgridn = qgrid[n];
  double ***egridn = egrid[n];

  double ***v0gridn = v0grid[n];
  double ***v1gridn = v1grid[n];
  double ***v2gridn = v2grid[n];
  double ***v3gridn = v3grid[n];
  double ***v4gridn = v4grid[n];
  double ***v5gridn = v5grid[n];

  int k = 0;

  if (flag == REVERSE_RHO) {
    double *qsrc = &qgridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++) {
      buf[k++] = qsrc[list[i]];
    }
  } else if (flag == REVERSE_AD) {
    double *src = &egridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == REVERSE_AD_PERATOM) {
    double *v0src = &v0gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v1src = &v1gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v2src = &v2gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v3src = &v3gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v4src = &v4gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v5src = &v5gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++) {
      buf[k++] = v0src[list[i]];
      buf[k++] = v1src[list[i]];
      buf[k++] = v2src[list[i]];
      buf[k++] = v3src[list[i]];
      buf[k++] = v4src[list[i]];
      buf[k++] = v5src[list[i]];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void MSM::unpack_reverse(int flag, double *buf, int nlist, int *list)
{
  int n = current_level;

  double ***qgridn = qgrid[n];
  double ***egridn = egrid[n];

  double ***v0gridn = v0grid[n];
  double ***v1gridn = v1grid[n];
  double ***v2gridn = v2grid[n];
  double ***v3gridn = v3grid[n];
  double ***v4gridn = v4grid[n];
  double ***v5gridn = v5grid[n];

  int k = 0;

  if (flag == REVERSE_RHO) {
    double *dest = &qgridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++) {
      dest[list[i]] += buf[k++];
    }
  } else if (flag == REVERSE_AD) {
    double *dest = &egridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[k++];
  } else if (flag == REVERSE_AD_PERATOM) {
    double *v0src = &v0gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v1src = &v1gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v2src = &v2gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v3src = &v3gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v4src = &v4gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    double *v5src = &v5gridn[nzlo_out[n]][nylo_out[n]][nxlo_out[n]];
    for (int i = 0; i < nlist; i++) {
      v0src[list[i]] += buf[k++];
      v1src[list[i]] += buf[k++];
      v2src[list[i]] += buf[k++];
      v3src[list[i]] += buf[k++];
      v4src[list[i]] += buf[k++];
      v5src[list[i]] += buf[k++];
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get force on my particles
------------------------------------------------------------------------- */

void MSM::fieldforce()
{
  //fprintf(screen,"MSM interpolation\n\n");

  double ***egridn = egrid[0];

  int i,l,m,n,nx,ny,nz,mx,my,mz;
  double dx,dy,dz;
  double phi_x,phi_y,phi_z;
  double dphi_x,dphi_y,dphi_z;
  double ekx,eky,ekz,etmp;


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
          etmp = egridn[mz][my][mx];
          ekx += dphi_x*phi_y*phi_z*etmp;
          eky += phi_x*dphi_y*phi_z*etmp;
          ekz += phi_x*phi_y*dphi_z*etmp;
        }
      }
    }

    ekx *= delxinv[0];
    eky *= delyinv[0];
    ekz *= delzinv[0];

    // effectively divide by length for a triclinic system

    if (triclinic) {
      double tmp[3];
      tmp[0] = ekx;
      tmp[1] = eky;
      tmp[2] = ekz;
      x2lamdaT(&tmp[0],&tmp[0]);
      ekx = tmp[0];
      eky = tmp[1];
      ekz = tmp[2];
    }

    // convert E-field to force

    const double qfactor = qqrd2e*scale*q[i];
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
   charge assignment into phi1d (interpolation coefficients)
------------------------------------------------------------------------- */

void MSM::compute_phis(const double &dx, const double &dy,
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
  }
}

/* ----------------------------------------------------------------------
   charge assignment into phi1d and dphi1d (interpolation coefficients)
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
   see Eq 7 from Parallel Computing 35 (2009) 164-177
   and Hardy's thesis
------------------------------------------------------------------------- */

inline double MSM::compute_phi(const double &xi)
{
  double phi = 0.0;
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
   see Eq 7 from Parallel Computing 35 (2009) 164-177
   and Hardy's thesis
------------------------------------------------------------------------- */

inline double MSM::compute_dphi(const double &xi)
{
  double dphi = 0.0;
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
    double abs_xi3 = xi2*abs_xi;
    double abs_xi5 = xi4*abs_xi;
    if (abs_xi == 0.0) {
      dphi = 0.0;
    } else if (abs_xi <= 1) {
      dphi = xi*(49*xi6 - 175*xi4 + 84*xi2 - 150*abs_xi5 +
        644*abs_xi3 - 560*abs_xi)/144.0/abs_xi;
    } else if (abs_xi <= 2) {
      dphi = xi*(-49*xi6 - 1365*xi4 + 756*xi2 +
        450*abs_xi5 + 1260*abs_xi3 - 1260*abs_xi + 28)/240.0/abs_xi;
    } else if (abs_xi <= 3) {
      dphi = xi*(49*xi6 + 4445*xi4 + 17724*xi2 -
        750*abs_xi5 - 12740*abs_xi3 - 9940*abs_xi + 756)/720.0/abs_xi;
    } else if (abs_xi <= 4) {
      dphi = -xi*(abs_xi - 4)*(7*abs_xi5 - 122*xi4 +
        807*abs_xi3 - 2512*xi2 + 3644*abs_xi - 1944)/720.0/abs_xi;
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
   Compute direct interaction (energy) weights for intermediate grid levels
------------------------------------------------------------------------- */
void MSM::get_g_direct()
{
  if (g_direct) memory->destroy(g_direct);
  memory->create(g_direct,levels,nmax_direct,"msm:g_direct");

  double a = cutoff;

  int n,zk,zyk,k,ix,iy,iz;
  double xdiff,ydiff,zdiff;
  double dx,dy,dz;
  double tmp[3];
  double rsq,rho,two_n;

  two_n = 1.0;

  int nx = nxhi_direct - nxlo_direct + 1;
  int ny = nyhi_direct - nylo_direct + 1;

  for (n=0; n<levels; n++) {

    for (iz = nzlo_direct; iz <= nzhi_direct; iz++) {
      zdiff = iz/delzinv[n];
      zk = (iz + nzhi_direct)*ny;
      for (iy = nylo_direct; iy <= nyhi_direct; iy++) {
        ydiff = iy/delyinv[n];
        zyk = (zk + iy + nyhi_direct)*nx;
        for (ix = nxlo_direct; ix <= nxhi_direct; ix++) {
          xdiff = ix/delxinv[n];

          // transform grid point pair-wise distance from lamda (0-1) coords to box coords

          if (triclinic) {
            tmp[0] = xdiff;
            tmp[1] = ydiff;
            tmp[2] = zdiff;
            lamda2xvector(&tmp[0],&tmp[0]);
            dx = tmp[0];
            dy = tmp[1];
            dz = tmp[2];
          } else {
            dx = xdiff;
            dy = ydiff;
            dz = zdiff;
          }

          rsq = dx*dx + dy*dy + dz*dz;

          rho = sqrt(rsq)/(two_n*a);
          k = zyk + ix + nxhi_direct;
          g_direct[n][k] = gamma(rho)/(two_n*a) - gamma(rho/2.0)/(2.0*two_n*a);
        }
      }
    }
    two_n *= 2.0;
  }
}

/* ----------------------------------------------------------------------
   Compute direct interaction (virial) weights for intermediate grid levels
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

  double a = cutoff;
  double a_sq = cutoff*cutoff;

  int n,zk,zyk,k,ix,iy,iz;
  double xdiff,ydiff,zdiff;
  double dx,dy,dz;
  double tmp[3];
  double rsq,r,rho,two_n,two_nsq,dg;

  two_n = 1.0;

  int nx = nxhi_direct - nxlo_direct + 1;
  int ny = nyhi_direct - nylo_direct + 1;

  for (n=0; n<levels; n++) {
    two_nsq = two_n * two_n;

    for (iz = nzlo_direct; iz <= nzhi_direct; iz++) {
      zdiff = iz/delzinv[n];
      zk = (iz + nzhi_direct)*ny;
      for (iy = nylo_direct; iy <= nyhi_direct; iy++) {
        ydiff = iy/delyinv[n];
        zyk = (zk + iy + nyhi_direct)*nx;
        for (ix = nxlo_direct; ix <= nxhi_direct; ix++) {
          xdiff = ix/delxinv[n];

          if (triclinic) {
            tmp[0] = xdiff;
            tmp[1] = ydiff;
            tmp[2] = zdiff;
            lamda2xvector(&tmp[0],&tmp[0]);
            dx = tmp[0];
            dy = tmp[1];
            dz = tmp[2];
          } else {
            dx = xdiff;
            dy = ydiff;
            dz = zdiff;
          }

          rsq = dx*dx + dy*dy + dz*dz;
          k = zyk + ix + nxhi_direct;
          r = sqrt(rsq);
          if (r == 0) {
            v0_direct[n][k] = 0.0;
            v1_direct[n][k] = 0.0;
            v2_direct[n][k] = 0.0;
            v3_direct[n][k] = 0.0;
            v4_direct[n][k] = 0.0;
            v5_direct[n][k] = 0.0;
          } else {
            rho = r/(two_n*a);
            dg = -(dgamma(rho)/(two_nsq*a_sq) -
              dgamma(rho/2.0)/(4.0*two_nsq*a_sq))/r;
            v0_direct[n][k] = dg * dx * dx;
            v1_direct[n][k] = dg * dy * dy;
            v2_direct[n][k] = dg * dz * dz;
            v3_direct[n][k] = dg * dx * dy;
            v4_direct[n][k] = dg * dx * dz;
            v5_direct[n][k] = dg * dy * dz;
          }

        }
      }
    }
    two_n *= 2.0;
  }
}

/* ----------------------------------------------------------------------
   Compute direct interaction (energy) weights for top grid level
   (nonperiodic systems only)
------------------------------------------------------------------------- */

void MSM::get_g_direct_top(int n)
{
  int nx_top = betax[n] - alpha[n];
  int ny_top = betay[n] - alpha[n];
  int nz_top = betaz[n] - alpha[n];

  int nx = 2*nx_top + 1;
  int ny = 2*ny_top + 1;
  int nz = 2*nz_top + 1;

  int nmax_top = 8*(nx+1)*(ny*1)*(nz+1);

  if (g_direct_top) memory->destroy(g_direct_top);
  memory->create(g_direct_top,nmax_top,"msm:g_direct_top");

  double a = cutoff;

  int zk,zyk,k,ix,iy,iz;
  double xdiff,ydiff,zdiff;
  double dx,dy,dz;
  double tmp[3];
  double rsq,rho,two_n;

  two_n = pow(2.0,n);

  for (iz = -nz_top; iz <= nz_top; iz++) {
    zdiff = iz/delzinv[n];
    zk = (iz + nz_top)*ny;
    for (iy = -ny_top; iy <= ny_top; iy++) {
      ydiff = iy/delyinv[n];
      zyk = (zk + iy + ny_top)*nx;
      for (ix = -nx_top; ix <= nx_top; ix++) {
        xdiff = ix/delxinv[n];

        if (triclinic) {
          tmp[0] = xdiff;
          tmp[1] = ydiff;
          tmp[2] = zdiff;
          lamda2xvector(&tmp[0],&tmp[0]);
          dx = tmp[0];
          dy = tmp[1];
          dz = tmp[2];
        } else {
          dx = xdiff;
          dy = ydiff;
          dz = zdiff;
        }

        rsq = dx*dx + dy*dy + dz*dz;
        rho = sqrt(rsq)/(two_n*a);
        k = zyk + ix + nx_top;
        g_direct_top[k] = gamma(rho)/(two_n*a);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Compute direct interaction (virial) weights for top grid level
   (nonperiodic systems only)
------------------------------------------------------------------------- */

void MSM::get_virial_direct_top(int n)
{
  int nx_top = betax[n] - alpha[n];
  int ny_top = betay[n] - alpha[n];
  int nz_top = betaz[n] - alpha[n];

  int nx = 2*nx_top + 1;
  int ny = 2*ny_top + 1;
  int nz = 2*nz_top + 1;

  int nmax_top = 8*(nx+1)*(ny*1)*(nz+1);

  if (v0_direct_top) memory->destroy(v0_direct_top);
  memory->create(v0_direct_top,nmax_top,"msm:v0_direct_top");
  if (v1_direct_top) memory->destroy(v1_direct_top);
  memory->create(v1_direct_top,nmax_top,"msm:v1_direct_top");
  if (v2_direct_top) memory->destroy(v2_direct_top);
  memory->create(v2_direct_top,nmax_top,"msm:v2_direct_top");
  if (v3_direct_top) memory->destroy(v3_direct_top);
  memory->create(v3_direct_top,nmax_top,"msm:v3_direct_top");
  if (v4_direct_top) memory->destroy(v4_direct_top);
  memory->create(v4_direct_top,nmax_top,"msm:v4_direct_top");
  if (v5_direct_top) memory->destroy(v5_direct_top);
  memory->create(v5_direct_top,nmax_top,"msm:v5_direct_top");

  double a = cutoff;
  double a_sq = cutoff*cutoff;

  int zk,zyk,k,ix,iy,iz;
  double xdiff,ydiff,zdiff;
  double dx,dy,dz;
  double tmp[3];
  double rsq,r,rho,two_n,two_nsq,dg;

  two_n = pow(2.0,n);
  two_nsq = two_n * two_n;

  for (iz = -nz_top; iz <= nz_top; iz++) {
    zdiff = iz/delzinv[n];
    zk = (iz + nz_top)*ny;
    for (iy = -ny_top; iy <= ny_top; iy++) {
      ydiff = iy/delyinv[n];
      zyk = (zk + iy + ny_top)*nx;
      for (ix = -nx_top; ix <= nx_top; ix++) {
        xdiff = ix/delxinv[n];
        if (triclinic) {
          tmp[0] = xdiff;
          tmp[1] = ydiff;
          tmp[2] = zdiff;
          lamda2xvector(&tmp[0],&tmp[0]);
          dx = tmp[0];
          dy = tmp[1];
          dz = tmp[2];
        } else {
          dx = xdiff;
          dy = ydiff;
          dz = zdiff;
        }

        rsq = dx*dx + dy*dy + dz*dz;
        k = zyk + ix + nx_top;
        r = sqrt(rsq);
        if (r == 0) {
          v0_direct_top[k] = 0.0;
          v1_direct_top[k] = 0.0;
          v2_direct_top[k] = 0.0;
          v3_direct_top[k] = 0.0;
          v4_direct_top[k] = 0.0;
          v5_direct_top[k] = 0.0;
        } else {
          rho = r/(two_n*a);
          dg = -(dgamma(rho)/(two_nsq*a_sq))/r;
          v0_direct_top[k] = dg * dx * dx;
          v1_direct_top[k] = dg * dy * dy;
          v2_direct_top[k] = dg * dz * dz;
          v3_direct_top[k] = dg * dx * dy;
          v4_direct_top[k] = dg * dx * dz;
          v5_direct_top[k] = dg * dy * dz;
        }
      }
    }
  }
}
