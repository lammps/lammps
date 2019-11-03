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
   Contributing authors: Stan Moore (SNL), Julien Tranchida (SNL)
------------------------------------------------------------------------- */

#include "pppm_dipole.h"
#include <mpi.h>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "gridcomm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "fft3d_wrap.h"
#include "remap_wrap.h"
#include "memory.h"
#include "error.h"
#include "update.h"

#include "math_const.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

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

PPPMDipole::PPPMDipole(LAMMPS *lmp) : PPPM(lmp),
  densityx_brick_dipole(NULL), densityy_brick_dipole(NULL), 
  densityz_brick_dipole(NULL),
  vdxx_brick_dipole(NULL), vdyy_brick_dipole(NULL), vdzz_brick_dipole(NULL),
  vdxy_brick_dipole(NULL), vdxz_brick_dipole(NULL), vdyz_brick_dipole(NULL),
  ux_brick_dipole(NULL), uy_brick_dipole(NULL), uz_brick_dipole(NULL),
  v0x_brick_dipole(NULL), v1x_brick_dipole(NULL),
  v2x_brick_dipole(NULL), v3x_brick_dipole(NULL), v4x_brick_dipole(NULL), 
  v5x_brick_dipole(NULL), v0y_brick_dipole(NULL), v1y_brick_dipole(NULL), 
  v2y_brick_dipole(NULL), v3y_brick_dipole(NULL), v4y_brick_dipole(NULL), 
  v5y_brick_dipole(NULL), v0z_brick_dipole(NULL), v1z_brick_dipole(NULL), 
  v2z_brick_dipole(NULL), v3z_brick_dipole(NULL), v4z_brick_dipole(NULL), 
  v5z_brick_dipole(NULL), work3(NULL), work4(NULL), 
  densityx_fft_dipole(NULL), densityy_fft_dipole(NULL), 
  densityz_fft_dipole(NULL)
{
  dipoleflag = 1;
  group_group_enable = 0;

  cg_dipole = NULL;
  cg_peratom_dipole = NULL;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

PPPMDipole::~PPPMDipole()
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

void PPPMDipole::init()
{
  if (me == 0) {
    if (screen) fprintf(screen,"PPPMDipole initialization ...\n");
    if (logfile) fprintf(logfile,"PPPMDipole initialization ...\n");
  }

  // error check

  dipoleflag = atom->mu?1:0;
  qsum_qsq(0); // q[i] might not be declared ?

  if (dipoleflag && q2)
    error->all(FLERR,"Cannot (yet) use charges with Kspace style PPPMDipole");

  triclinic_check();

  if (triclinic != domain->triclinic)
    error->all(FLERR,"Must redefine kspace_style after changing to triclinic box");

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use PPPMDipole with 2d simulation");

  if (comm->style != 0)
    error->universe_all(FLERR,"PPPMDipole can only currently be used with "
                        "comm_style brick");

  if (!atom->mu) error->all(FLERR,"Kspace style requires atom attribute mu");

  if (atom->mu && differentiation_flag == 1) error->all(FLERR,"Cannot (yet) use kspace_modify diff"
       " ad with dipoles");

  if (dipoleflag && strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use nonperiodic boundaries with PPPMDipole");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab PPPMDipole");
  }

  if (order < 2 || order > MAXORDER) {
    char str[128];
    sprintf(str,"PPPMDipole order cannot be < 2 or > than %d",MAXORDER);
    error->all(FLERR,str);
  }

  // extract short-range Coulombic cutoff from pair style

  triclinic = domain->triclinic;
  if (triclinic)
    error->all(FLERR,"Cannot yet use triclinic cells with PPPMDipole");

  pair_check();

  int itmp = 0;
  double *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  if (p_cutoff == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;

  // kspace TIP4P not yet supported
  // qdist = offset only for TIP4P fictitious charge

  qdist = 0.0;
  if (tip4pflag)
    error->all(FLERR,"Cannot yet use TIP4P with PPPMDipole");

  // compute musum & musqsum and warn if no dipoles

  scale = 1.0;
  qqrd2e = force->qqrd2e;
  musum_musq();
  natoms_original = atom->natoms;

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

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
      error->warning(FLERR,"Reducing PPPMDipole order b/c stencil extends "
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

  if (order < minorder) error->all(FLERR,"PPPMDipole order < minimum allowed order");
  if (!overlap_allowed && cgtmp->ghost_overlap())
    error->all(FLERR,"PPPMDipole grid stencil extends "
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

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients

  compute_gf_denom();
  compute_rho_coeff();
}

/* ----------------------------------------------------------------------
   adjust PPPMDipole coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void PPPMDipole::setup()
{
  // perform some checks to avoid illegal boundaries with read_data

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use nonperiodic boundaries with PPPMDipole");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab PPPMDipole");
  }

  int i,j,k,n;
  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPMDipole
  // z dimension for 3d PPPMDipole is zprd since slab_volfactor = 1.0

  prd = domain->prd;
  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  delxinv = nx_pppm/xprd;
  delyinv = ny_pppm/yprd;
  delzinv = nz_pppm/zprd_slab;

  delvolinv = delxinv*delyinv*delzinv;

  double unitkx = (MY_2PI/xprd);
  double unitky = (MY_2PI/yprd);
  double unitkz = (MY_2PI/zprd_slab);

  // fkx,fky,fkz for my FFT grid pts

  double per;

  for (i = nxlo_fft; i <= nxhi_fft; i++) {
    per = i - nx_pppm*(2*i/nx_pppm);
    fkx[i] = unitkx*per;
  }

  for (i = nylo_fft; i <= nyhi_fft; i++) {
    per = i - ny_pppm*(2*i/ny_pppm);
    fky[i] = unitky*per;
  }

  for (i = nzlo_fft; i <= nzhi_fft; i++) {
    per = i - nz_pppm*(2*i/nz_pppm);
    fkz[i] = unitkz*per;
  }

  // virial coefficients

  double sqk,vterm;

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++) {
    for (j = nylo_fft; j <= nyhi_fft; j++) {
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        sqk = fkx[i]*fkx[i] + fky[j]*fky[j] + fkz[k]*fkz[k];
        if (sqk == 0.0) {
          vg[n][0] = 0.0;
          vg[n][1] = 0.0;
          vg[n][2] = 0.0;
          vg[n][3] = 0.0;
          vg[n][4] = 0.0;
          vg[n][5] = 0.0;
        } else {
          vterm = -2.0 * (1.0/sqk + 0.25/(g_ewald*g_ewald));
          vg[n][0] = 1.0 + vterm*fkx[i]*fkx[i];
          vg[n][1] = 1.0 + vterm*fky[j]*fky[j];
          vg[n][2] = 1.0 + vterm*fkz[k]*fkz[k];
          vg[n][3] = vterm*fkx[i]*fky[j];
          vg[n][4] = vterm*fkx[i]*fkz[k];
          vg[n][5] = vterm*fky[j]*fkz[k];
        }
        n++;
      }
    }
  }

  compute_gf_dipole();
}

/* ----------------------------------------------------------------------
   reset local grid arrays and communication stencils
   called by fix balance b/c it changed sizes of processor sub-domains
------------------------------------------------------------------------- */

void PPPMDipole::setup_grid()
{
  // free all arrays previously allocated

  deallocate();
  if (peratom_allocate_flag) deallocate_peratom();

  // reset portion of global grid that each proc owns

  set_grid_local();

  // reallocate K-space dependent memory
  // check if grid communication is now overlapping if not allowed
  // don't invoke allocate peratom(), will be allocated when needed

  allocate();

  cg_dipole->ghost_notify();
  if (overlap_allowed == 0 && cg_dipole->ghost_overlap())
    error->all(FLERR,"PPPMDipole grid stencil extends "
               "beyond nearest neighbor processor");
  cg_dipole->setup();

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients

  compute_gf_denom();
  compute_rho_coeff();

  // pre-compute volume-dependent coeffs

  setup();
}

/* ----------------------------------------------------------------------
   compute the PPPMDipole long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPMDipole::compute(int eflag, int vflag)
{
  int i,j;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  if (vflag_atom)
    error->all(FLERR,"Cannot (yet) compute per-atom virial "
                       "with kspace style pppm/dipole");

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    cg_peratom_dipole->ghost_notify();
    cg_peratom_dipole->setup();
  }

  // if atom count has changed, update qsum and qsqsum

  if (atom->natoms != natoms_original) {
    musum_musq();
    natoms_original = atom->natoms;
  }

  // return if there are no dipoles

  if (musqsum == 0.0) return;

  // convert atoms from box to lamda coords

  boxlo = domain->boxlo;

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm_dipole:part2grid");
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid

  particle_map();
  make_rho_dipole();

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

  fieldforce_ik_dipole();

  // extra per-atom energy/virial communication

  if (evflag_atom) fieldforce_peratom_dipole();

  // sum global energy across procs and add in volume-dependent term

  const double qscale = qqrd2e * scale;
  const double g3 = g_ewald*g_ewald*g_ewald;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    energy *= 0.5*volume;
    energy -= musqsum*2.0*g3/3.0/MY_PIS;
    energy *= qscale;
  }

  // sum global virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*volume*virial_all[i];
  }

  // per-atom energy/virial
  // energy includes self-energy correction

  if (evflag_atom) {
    double **mu = atom->mu;
    int nlocal = atom->nlocal;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
        eatom[i] *= 0.5;
        eatom[i] -= (mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2])*2.0*g3/3.0/MY_PIS;
        eatom[i] *= qscale;
      }
    }

    if (vflag_atom) {
      for (i = 0; i < nlocal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5*qscale;
    }
  }

  // 2d slab correction

  if (slabflag == 1) slabcorr();
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMDipole::allocate()
{
  memory->create3d_offset(densityx_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:densityx_brick_dipole");
  memory->create3d_offset(densityy_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:densityy_brick_dipole");
  memory->create3d_offset(densityz_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:densityz_brick_dipole");

  memory->create(densityx_fft_dipole,nfft_both,"pppm_dipole:densityy_fft_dipole");
  memory->create(densityy_fft_dipole,nfft_both,"pppm_dipole:densityy_fft_dipole");
  memory->create(densityz_fft_dipole,nfft_both,"pppm_dipole:densityz_fft_dipole");

  memory->create(greensfn,nfft_both,"pppm_dipole:greensfn");
  memory->create(work1,2*nfft_both,"pppm_dipole:work1");
  memory->create(work2,2*nfft_both,"pppm_dipole:work2");
  memory->create(work3,2*nfft_both,"pppm_dipole:work3");
  memory->create(work4,2*nfft_both,"pppm_dipole:work4");
  memory->create(vg,nfft_both,6,"pppm_dipole:vg");

  memory->create1d_offset(fkx,nxlo_fft,nxhi_fft,"pppm_dipole:fkx");
  memory->create1d_offset(fky,nylo_fft,nyhi_fft,"pppm_dipole:fky");
  memory->create1d_offset(fkz,nzlo_fft,nzhi_fft,"pppm_dipole:fkz");

  memory->create3d_offset(ux_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:ux_brick_dipole");
  memory->create3d_offset(uy_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:uy_brick_dipole");
  memory->create3d_offset(uz_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:uz_brick_dipole");

  memory->create3d_offset(vdxx_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:vdxx_brick_dipole");
  memory->create3d_offset(vdxy_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:vdxy_brick_dipole");
  memory->create3d_offset(vdyy_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:vdyy_brick_dipole");
  memory->create3d_offset(vdxz_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:vdxz_brick_dipole");
  memory->create3d_offset(vdyz_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:vdyz_brick_dipole");
  memory->create3d_offset(vdzz_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:vdzz_brick_dipole");

  // summation coeffs

  order_allocated = order;
  if (!gf_b) memory->create(gf_b,order,"pppm_dipole:gf_b");
  memory->create2d_offset(rho1d,3,-order/2,order/2,"pppm_dipole:rho1d");
  memory->create2d_offset(drho1d,3,-order/2,order/2,"pppm_dipole:drho1d");
  memory->create2d_offset(rho_coeff,order,(1-order)/2,order/2,"pppm_dipole:rho_coeff");
  memory->create2d_offset(drho_coeff,order,(1-order)/2,order/2,
                          "pppm_dipole:drho_coeff");

  // create 2 FFTs and a Remap
  // 1st FFT keeps data in FFT decompostion
  // 2nd FFT returns data in 3d brick decomposition
  // remap takes data from 3d brick to FFT decomposition

  int tmp;

  fft1 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   0,0,&tmp,collective_flag);

  fft2 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                   0,0,&tmp,collective_flag);

  remap = new Remap(lmp,world,
                    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                    nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                    1,0,0,FFT_PRECISION,collective_flag);

  // create ghost grid object for rho and electric field communication

  int (*procneigh)[2] = comm->procneigh;

  cg_dipole = new GridComm(lmp,world,9,3,
                           nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                           nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                           procneigh[0][0],procneigh[0][1],procneigh[1][0],
                           procneigh[1][1],procneigh[2][0],procneigh[2][1]);
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMDipole::deallocate()
{
  memory->destroy3d_offset(densityx_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(densityy_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(densityz_brick_dipole,nzlo_out,nylo_out,nxlo_out);

  memory->destroy3d_offset(ux_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(uy_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(uz_brick_dipole,nzlo_out,nylo_out,nxlo_out);

  memory->destroy3d_offset(vdxx_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(vdxy_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(vdyy_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(vdxz_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(vdyz_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(vdzz_brick_dipole,nzlo_out,nylo_out,nxlo_out);

  memory->destroy(densityx_fft_dipole);
  memory->destroy(densityy_fft_dipole);
  memory->destroy(densityz_fft_dipole);

  memory->destroy(greensfn);
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(work3);
  memory->destroy(work4);
  memory->destroy(vg);

  memory->destroy1d_offset(fkx,nxlo_fft);
  memory->destroy1d_offset(fky,nylo_fft);
  memory->destroy1d_offset(fkz,nzlo_fft);

  memory->destroy(gf_b);
  memory->destroy2d_offset(rho1d,-order_allocated/2);
  memory->destroy2d_offset(drho1d,-order_allocated/2);
  memory->destroy2d_offset(rho_coeff,(1-order_allocated)/2);
  memory->destroy2d_offset(drho_coeff,(1-order_allocated)/2);

  delete fft1;
  delete fft2;
  delete remap;
  delete cg_dipole;
}

/* ----------------------------------------------------------------------
   allocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMDipole::allocate_peratom()
{
  peratom_allocate_flag = 1;

  memory->create3d_offset(v0x_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v0x_brick_dipole");
  memory->create3d_offset(v1x_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v1x_brick_dipole");
  memory->create3d_offset(v2x_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v2x_brick_dipole");
  memory->create3d_offset(v3x_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v3x_brick_dipole");
  memory->create3d_offset(v4x_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v4x_brick_dipole");
  memory->create3d_offset(v5x_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v5x_brick_dipole");

  memory->create3d_offset(v0y_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v0y_brick_dipole");
  memory->create3d_offset(v1y_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v1y_brick_dipole");
  memory->create3d_offset(v2y_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v2y_brick_dipole");
  memory->create3d_offset(v3y_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v3y_brick_dipole");
  memory->create3d_offset(v4y_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v4y_brick_dipole");
  memory->create3d_offset(v5y_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v5y_brick_dipole");

  memory->create3d_offset(v0z_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v0z_brick_dipole");
  memory->create3d_offset(v1z_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v1z_brick_dipole");
  memory->create3d_offset(v2z_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v2z_brick_dipole");
  memory->create3d_offset(v3z_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v3z_brick_dipole");
  memory->create3d_offset(v4z_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v4z_brick_dipole");
  memory->create3d_offset(v5z_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm_dipole:v5z_brick_dipole");

  // create ghost grid object for rho and electric field communication

  int (*procneigh)[2] = comm->procneigh;

  cg_peratom_dipole =
    new GridComm(lmp,world,18,1,
                 nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                 nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                 procneigh[0][0],procneigh[0][1],procneigh[1][0],
                 procneigh[1][1],procneigh[2][0],procneigh[2][1]);
}

/* ----------------------------------------------------------------------
   deallocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMDipole::deallocate_peratom()
{
  peratom_allocate_flag = 0;

  memory->destroy3d_offset(v0x_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v1x_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v2x_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v3x_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v4x_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v5x_brick_dipole,nzlo_out,nylo_out,nxlo_out);

  memory->destroy3d_offset(v0y_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v1y_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v2y_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v3y_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v4y_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v5y_brick_dipole,nzlo_out,nylo_out,nxlo_out);

  memory->destroy3d_offset(v0z_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v1z_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v2z_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v3z_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v4z_brick_dipole,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v5z_brick_dipole,nzlo_out,nylo_out,nxlo_out);

  delete cg_peratom_dipole;
}

/* ----------------------------------------------------------------------
   set global size of PPPMDipole grid = nx,ny,nz_pppm
   used for charge accumulation, FFTs, and electric field interpolation
------------------------------------------------------------------------- */

void PPPMDipole::set_grid_global()
{
  // use xprd,yprd,zprd
  // adjust z dimension for 2d slab PPPMDipole
  // 3d PPPMDipole just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;

  // make initial g_ewald estimate
  // based on desired accuracy and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  double h;
  bigint natoms = atom->natoms;

  if (!gewaldflag) {
    if (accuracy <= 0.0)
      error->all(FLERR,"KSpace accuracy must be > 0");
    if (mu2 == 0.0)
     error->all(FLERR,"Must use kspace_modify gewald for systems with no dipoles");
    g_ewald = (1.35 - 0.15*log(accuracy))/cutoff;
    double g_ewald_new =
      find_gewald_dipole(g_ewald,cutoff,natoms,xprd*yprd*zprd,mu2);
    if (g_ewald_new > 0.0) g_ewald = g_ewald_new;
    else error->warning(FLERR,"PPPMDipole dipole Newton solver failed, "
                        "using old method to estimate g_ewald");
  }

  // set optimal nx_pppm,ny_pppm,nz_pppm based on order and accuracy
  // nz_pppm uses extended zprd_slab instead of zprd
  // reduce it until accuracy target is met

  if (!gridflag) {

    h = h_x = h_y = h_z = 4.0/g_ewald;
    int count = 0;
    while (1) {

      // set grid dimension

      nx_pppm = static_cast<int> (xprd/h_x);
      ny_pppm = static_cast<int> (yprd/h_y);
      nz_pppm = static_cast<int> (zprd_slab/h_z);

      if (nx_pppm <= 1) nx_pppm = 2;
      if (ny_pppm <= 1) ny_pppm = 2;
      if (nz_pppm <= 1) nz_pppm = 2;

      // set local grid dimension

      int npey_fft,npez_fft;
      if (nz_pppm >= nprocs) {
        npey_fft = 1;
        npez_fft = nprocs;
      } else procs2grid2d(nprocs,ny_pppm,nz_pppm,&npey_fft,&npez_fft);

      int me_y = me % npey_fft;
      int me_z = me / npey_fft;

      nxlo_fft = 0;
      nxhi_fft = nx_pppm - 1;
      nylo_fft = me_y*ny_pppm/npey_fft;
      nyhi_fft = (me_y+1)*ny_pppm/npey_fft - 1;
      nzlo_fft = me_z*nz_pppm/npez_fft;
      nzhi_fft = (me_z+1)*nz_pppm/npez_fft - 1;

      double df_kspace = compute_df_kspace_dipole();

      count++;

      // break loop if the accuracy has been reached or
      // too many loops have been performed

      if (df_kspace <= accuracy) break;
      if (count > 500) error->all(FLERR, "Could not compute grid size");
      h *= 0.95;
      h_x = h_y = h_z = h;
    }
  }

  // boost grid size until it is factorable

  while (!factorable(nx_pppm)) nx_pppm++;
  while (!factorable(ny_pppm)) ny_pppm++;
  while (!factorable(nz_pppm)) nz_pppm++;

  h_x = xprd/nx_pppm;
  h_y = yprd/ny_pppm;
  h_z = zprd_slab/nz_pppm;

  if (nx_pppm >= OFFSET || ny_pppm >= OFFSET || nz_pppm >= OFFSET)
    error->all(FLERR,"PPPMDipole grid is too large");
}


/* ----------------------------------------------------------------------
   compute estimated kspace force error for dipoles
------------------------------------------------------------------------- */

double PPPMDipole::compute_df_kspace_dipole()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;
  bigint natoms = atom->natoms;
  double qopt = compute_qopt_dipole();
  double df_kspace = sqrt(qopt/natoms)*mu2/(3.0*xprd*yprd*zprd_slab);
  return df_kspace;
}

/* ----------------------------------------------------------------------
   compute qopt for dipoles with ik differentiation
------------------------------------------------------------------------- */

double PPPMDipole::compute_qopt_dipole()
{
  double qopt = 0.0;
  const double * const prd = domain->prd;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double zprd_slab = zprd*slab_volfactor;
  const double unitkx = (MY_2PI/xprd);
  const double unitky = (MY_2PI/yprd);
  const double unitkz = (MY_2PI/zprd_slab);

  double snx,sny,snz;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2,dot1,dot2;
  double denominator;
  double u1,sqk;

  int k,l,m,nx,ny,nz,kper,lper,mper;

  const int nbx = 2;
  const int nby = 2;
  const int nbz = 2;

  const int twoorder = 2*order;

  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);
    snz = square(sin(0.5*unitkz*mper*zprd_slab/nz_pppm));

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);
      sny = square(sin(0.5*unitky*lper*yprd/ny_pppm));

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);
        snx = square(sin(0.5*unitkx*kper*xprd/nx_pppm));

        sqk = square(unitkx*kper) + square(unitky*lper) + square(unitkz*mper);

        if (sqk != 0.0) {
          denominator = gf_denom(snx,sny,snz);
          sum1 = 0.0;
          sum2 = 0.0;

          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);
            sx = exp(-0.25*square(qx/g_ewald));
            argx = 0.5*qx*xprd/nx_pppm;
            wx = powsinxx(argx,twoorder);

            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);
              sy = exp(-0.25*square(qy/g_ewald));
              argy = 0.5*qy*yprd/ny_pppm;
              wy = powsinxx(argy,twoorder);

              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);
                sz = exp(-0.25*square(qz/g_ewald));
                argz = 0.5*qz*zprd_slab/nz_pppm;
                wz = powsinxx(argz,twoorder);

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                dot2 = qx*qx + qy*qy + qz*qz;
                //dot1 = dot1*dot1*dot1; // power of 3 for dipole forces
                //dot2 = dot2*dot2*dot2;
                u1 = sx*sy*sz;
                const double w2 = wx*wy*wz;
                const double phi = u1*MY_4PI/dot2;
                const double top = dot1*dot1*dot1*w2*phi;
                sum1 += phi*phi*dot2*dot2*dot2;
                sum2 += top*top/sqk/sqk/sqk;
              }
            }
          }
          qopt += sum1 - sum2/denominator;
        }
      }
    }
  }
  double qopt_all;
  MPI_Allreduce(&qopt,&qopt_all,1,MPI_DOUBLE,MPI_SUM,world);
  return qopt_all;
}

/* ----------------------------------------------------------------------
   pre-compute modified (Hockney-Eastwood) Coulomb Green's function
------------------------------------------------------------------------- */

void PPPMDipole::compute_gf_dipole()
{
  const double * const prd = domain->prd;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double zprd_slab = zprd*slab_volfactor;
  const double unitkx = (MY_2PI/xprd);
  const double unitky = (MY_2PI/yprd);
  const double unitkz = (MY_2PI/zprd_slab);

  double snx,sny,snz;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,dot1,dot2;
  double denominator;
  double sqk;

  int k,l,m,n,nx,ny,nz,kper,lper,mper;

  int nbx = static_cast<int> ((g_ewald*xprd/(MY_PI*nx_pppm)) *
                              pow(-log(EPS_HOC),0.25));
  int nby = static_cast<int> ((g_ewald*yprd/(MY_PI*ny_pppm)) *
                              pow(-log(EPS_HOC),0.25));
  int nbz = static_cast<int> ((g_ewald*zprd_slab/(MY_PI*nz_pppm)) *
                              pow(-log(EPS_HOC),0.25));
  nbx = MAX(nbx,2);
  nby = MAX(nby,2);
  nbz = MAX(nbz,2);
  const int twoorder = 2*order;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);
    snz = square(sin(0.5*unitkz*mper*zprd_slab/nz_pppm));

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);
      sny = square(sin(0.5*unitky*lper*yprd/ny_pppm));

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);
        snx = square(sin(0.5*unitkx*kper*xprd/nx_pppm));

        sqk = square(unitkx*kper) + square(unitky*lper) + square(unitkz*mper);

        if (sqk != 0.0) {
          denominator = gf_denom(snx,sny,snz);
          sum1 = 0.0;

          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);
            sx = exp(-0.25*square(qx/g_ewald));
            argx = 0.5*qx*xprd/nx_pppm;
            wx = powsinxx(argx,twoorder);

            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);
              sy = exp(-0.25*square(qy/g_ewald));
              argy = 0.5*qy*yprd/ny_pppm;
              wy = powsinxx(argy,twoorder);

              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);
                sz = exp(-0.25*square(qz/g_ewald));
                argz = 0.5*qz*zprd_slab/nz_pppm;
                wz = powsinxx(argz,twoorder);

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                dot2 = qx*qx + qy*qy + qz*qz;
                const double u1 = sx*sy*sz;
                const double w2 = wx*wy*wz;
                const double phi = u1*MY_4PI/dot2;
                sum1 += dot1*dot1*dot1*w2*phi/sqk/sqk/sqk;
              }
            }
          }
          greensfn[n++] = sum1/denominator;
        } else greensfn[n++] = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   calculate f(x) for use in Newton-Raphson solver
------------------------------------------------------------------------- */

double PPPMDipole::newton_raphson_f()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  bigint natoms = atom->natoms;

  double df_rspace,df_kspace;
  double vol = xprd*yprd*zprd;
  double a = cutoff*g_ewald;
  double rg2 = a*a;
  double rg4 = rg2*rg2;
  double rg6 = rg4*rg2;
  double Cc = 4.0*rg4 + 6.0*rg2 + 3.0;
  double Dc = 8.0*rg6 + 20.0*rg4 + 30.0*rg2 + 15.0;
  df_rspace = (mu2/(sqrt(vol*powint(g_ewald,4)*powint(cutoff,9)*natoms)) *
      sqrt(13.0/6.0*Cc*Cc + 2.0/15.0*Dc*Dc - 13.0/15.0*Cc*Dc) * exp(-rg2));
  df_kspace = compute_df_kspace_dipole();

  return df_rspace - df_kspace;
}

/* ----------------------------------------------------------------------
   find g_ewald parameter for dipoles based on desired accuracy
   using a Newton-Raphson solver
------------------------------------------------------------------------- */

double PPPMDipole::find_gewald_dipole(double x, double Rc,
                              bigint natoms, double vol, double b2)
{
  double dx,tol;
  int maxit;

  maxit = 10000; //Maximum number of iterations
  tol = 0.00001; //Convergence tolerance

  //Begin algorithm

  for (int i = 0; i < maxit; i++) {
    dx = newton_raphson_f_dipole(x,Rc,natoms,vol,b2) / derivf_dipole(x,Rc,natoms,vol,b2);
    x = x - dx; //Update x
    if (fabs(dx) < tol) return x;
    if (x < 0 || x != x) // solver failed
      return -1;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   calculate f(x) objective function for dipoles
 ------------------------------------------------------------------------- */

double PPPMDipole::newton_raphson_f_dipole(double x, double Rc, bigint
natoms, double vol, double b2)
{
  double a = Rc*x;
  double rg2 = a*a;
  double rg4 = rg2*rg2;
  double rg6 = rg4*rg2;
  double Cc = 4.0*rg4 + 6.0*rg2 + 3.0;
  double Dc = 8.0*rg6 + 20.0*rg4 + 30.0*rg2 + 15.0;
  double f = (b2/(sqrt(vol*powint(x,4)*powint(Rc,9)*natoms)) *
    sqrt(13.0/6.0*Cc*Cc + 2.0/15.0*Dc*Dc - 13.0/15.0*Cc*Dc) *
    exp(-rg2)) - accuracy;

  return f;
}

/* ----------------------------------------------------------------------
   calculate numerical derivative f'(x) of objective function for dipoles
 ------------------------------------------------------------------------- */

double PPPMDipole::derivf_dipole(double x, double Rc,
                         bigint natoms, double vol, double b2)
{
  double h = 0.000001;  //Derivative step-size
  return (newton_raphson_f_dipole(x + h,Rc,natoms,vol,b2) - newton_raphson_f_dipole(x,Rc,natoms,vol,b2)) / h;
}

/* ----------------------------------------------------------------------
   calculate the final estimate of the accuracy
------------------------------------------------------------------------- */

double PPPMDipole::final_accuracy_dipole()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double vol = xprd*yprd*zprd;
  bigint natoms = atom->natoms;
  if (natoms == 0) natoms = 1; // avoid division by zero

  double df_kspace = compute_df_kspace_dipole();

  double a = cutoff*g_ewald;
  double rg2 = a*a;
  double rg4 = rg2*rg2;
  double rg6 = rg4*rg2;
  double Cc = 4.0*rg4 + 6.0*rg2 + 3.0;
  double Dc = 8.0*rg6 + 20.0*rg4 + 30.0*rg2 + 15.0;
  double df_rspace = (mu2/(sqrt(vol*powint(g_ewald,4)*powint(cutoff,9)*natoms)) *
    sqrt(13.0/6.0*Cc*Cc + 2.0/15.0*Dc*Dc - 13.0/15.0*Cc*Dc) *
    exp(-rg2));

  double estimated_accuracy = sqrt(df_kspace*df_kspace + df_rspace*df_rspace);

  return estimated_accuracy;
}

/* ----------------------------------------------------------------------
   pre-compute Green's function denominator expansion coeffs, Gamma(2n)
------------------------------------------------------------------------- */

void PPPMDipole::compute_gf_denom()
{
  if (gf_b) memory->destroy(gf_b);
  memory->create(gf_b,order,"pppm_dipole:gf_b");

  int k,l,m;

  for (l = 1; l < order; l++) gf_b[l] = 0.0;
  gf_b[0] = 1.0;

  for (m = 1; m < order; m++) {
    for (l = m; l > 0; l--)
      gf_b[l] = 4.0 * (gf_b[l]*(l-m)*(l-m-0.5)-gf_b[l-1]*(l-m-1)*(l-m-1));
    gf_b[0] = 4.0 * (gf_b[0]*(l-m)*(l-m-0.5));
  }

  bigint ifact = 1;
  for (k = 1; k < 2*order; k++) ifact *= k;
  double gaminv = 1.0/ifact;
  for (l = 0; l < order; l++) gf_b[l] *= gaminv;
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void PPPMDipole::make_rho_dipole()
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

  double **mu = atom->mu;
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

    z0 = delvolinv * mu[i][0];
    z1 = delvolinv * mu[i][1];
    z2 = delvolinv * mu[i][2];
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
   remap density from 3d brick decomposition to FFT decomposition
------------------------------------------------------------------------- */

void PPPMDipole::brick2fft_dipole()
{
  int n,ix,iy,iz;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        densityx_fft_dipole[n] = densityx_brick_dipole[iz][iy][ix];
        densityy_fft_dipole[n] = densityy_brick_dipole[iz][iy][ix];
        densityz_fft_dipole[n] = densityz_brick_dipole[iz][iy][ix];
        n++;
      }

  remap->perform(densityx_fft_dipole,densityx_fft_dipole,work1);
  remap->perform(densityy_fft_dipole,densityy_fft_dipole,work1);
  remap->perform(densityz_fft_dipole,densityz_fft_dipole,work1);
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for ik
------------------------------------------------------------------------- */

void PPPMDipole::poisson_ik_dipole()
{
  int i,j,k,n,ii;
  double eng;
  double wreal,wimg;

  // transform dipole density (r -> k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n] = densityx_fft_dipole[i];
    work1[n+1] = ZEROF;
    work2[n] = densityy_fft_dipole[i];
    work2[n+1] = ZEROF;
    work3[n] = densityz_fft_dipole[i];
    work3[n+1] = ZEROF;
    n += 2;
  }

  fft1->compute(work1,work1,1);
  fft1->compute(work2,work2,1);
  fft1->compute(work3,work3,1);

  // global energy and virial contribution

  double scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm);
  double s2 = scaleinv*scaleinv;

  if (eflag_global || vflag_global) {
    if (vflag_global) {
      n = 0;
      ii = 0;
      for (k = nzlo_fft; k <= nzhi_fft; k++)
        for (j = nylo_fft; j <= nyhi_fft; j++)
          for (i = nxlo_fft; i <= nxhi_fft; i++) {
            wreal = (work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
            wimg = (work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
            eng = s2 * greensfn[ii] * (wreal*wreal + wimg*wimg);
            for (int jj = 0; jj < 6; jj++) virial[jj] += eng*vg[ii][jj];
            virial[0] += 2.0*s2*greensfn[ii]*fkx[i]*(work1[n]*wreal + work1[n+1]*wimg);
            virial[1] += 2.0*s2*greensfn[ii]*fky[j]*(work2[n]*wreal + work2[n+1]*wimg);
            virial[2] += 2.0*s2*greensfn[ii]*fkz[k]*(work3[n]*wreal + work3[n+1]*wimg);
            virial[3] += 2.0*s2*greensfn[ii]*fky[j]*(work1[n]*wreal + work1[n+1]*wimg);
            virial[4] += 2.0*s2*greensfn[ii]*fkz[k]*(work1[n]*wreal + work1[n+1]*wimg);
            virial[5] += 2.0*s2*greensfn[ii]*fkz[k]*(work2[n]*wreal + work2[n+1]*wimg);
            if (eflag_global) energy += eng;
            ii++;
            n += 2;
          }
    } else {
      n = 0;
      ii = 0;
      for (k = nzlo_fft; k <= nzhi_fft; k++)
        for (j = nylo_fft; j <= nyhi_fft; j++)
          for (i = nxlo_fft; i <= nxhi_fft; i++) {
            wreal = (work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
            wimg = (work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
            energy +=
            s2 * greensfn[ii] * (wreal*wreal + wimg*wimg);
   ii++;
            n += 2;
          }
    }
  }

  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n]   *= scaleinv * greensfn[i];
    work1[n+1] *= scaleinv * greensfn[i];
    work2[n]   *= scaleinv * greensfn[i];
    work2[n+1] *= scaleinv * greensfn[i];
    work3[n]   *= scaleinv * greensfn[i];
    work3[n+1] *= scaleinv * greensfn[i];
    n += 2;
  }

  // extra FFTs for per-atom energy/virial

  if (vflag_atom) poisson_peratom_dipole();

  // compute electric potential
  // FFT leaves data in 3d brick decomposition

  // Ex

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        work4[n+1] = fkx[i]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        ux_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Ey

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        work4[n+1] = fky[j]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        uy_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Ez

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkz[k]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        work4[n+1] = fkz[k]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        uz_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vxx

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*fkx[i]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fkx[i]*fkx[i]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdxx_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vyy

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*fky[j]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fky[j]*fky[j]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdyy_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vzz

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkz[k]*fkz[k]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fkz[k]*fkz[k]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdzz_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vxy

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*fky[j]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fkx[i]*fky[j]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdxy_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vxz

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*fkz[k]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fkx[i]*fkz[k]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdxz_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vyz

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*fkz[k]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fky[j]*fkz[k]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdyz_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for per-atom energy/virial
------------------------------------------------------------------------- */

void PPPMDipole::poisson_peratom_dipole()
{
  int i,ii,j,k,n;

  // 18 components of virial in v0 thru v5

  if (!vflag_atom) return;

  // V0x

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*(vg[ii][0]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkx[i]*work1[n]);
        work4[n+1] = fkx[i]*(vg[ii][0]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkx[i]*work1[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v0x_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V0y

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*(vg[ii][0]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkx[i]*work1[n]);
        work4[n+1] = fky[j]*(vg[ii][0]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkx[i]*work1[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v0y_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V0z

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkz[k]*(vg[ii][0]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkx[i]*work1[n]);
        work4[n+1] = fkz[k]*(vg[ii][0]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkx[i]*work1[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v0z_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V1x

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*(vg[ii][1]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fky[j]*work2[n]);
        work4[n+1] = fkx[i]*(vg[ii][1]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fky[j]*work2[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v1x_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V1y

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*(vg[ii][1]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fky[j]*work2[n]);
        work4[n+1] = fky[j]*(vg[ii][1]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fky[j]*work2[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v1y_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V1z

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkz[k]*(vg[ii][1]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fky[j]*work2[n]);
        work4[n+1] = fkz[k]*(vg[ii][1]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fky[j]*work2[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v1z_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V2x

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*(vg[ii][2]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkz[k]*work3[n]);
        work4[n+1] = fkx[i]*(vg[ii][2]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkz[k]*work3[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v2x_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V2y

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*(vg[ii][2]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkz[k]*work3[n]);
        work4[n+1] = fky[j]*(vg[ii][2]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkz[k]*work3[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v2y_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V2z

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkz[k]*(vg[ii][2]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkz[k]*work3[n]);
        work4[n+1] = fkz[k]*(vg[ii][2]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkz[k]*work3[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v2z_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V3x

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*(vg[ii][3]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fky[j]*work1[n]);
        work4[n+1] = fkx[i]*(vg[ii][3]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fky[j]*work1[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v3x_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V3y

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*(vg[ii][3]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fky[j]*work1[n]);
        work4[n+1] = fky[j]*(vg[ii][3]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fky[j]*work1[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v3y_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V3z

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkz[k]*(vg[ii][3]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fky[j]*work1[n]);
        work4[n+1] = fkz[k]*(vg[ii][3]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fky[j]*work1[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v3z_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V4x

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*(vg[ii][4]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkz[k]*work1[n]);
        work4[n+1] = fkx[i]*(vg[ii][4]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkz[k]*work1[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v4x_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V4y

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*(vg[ii][4]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkz[k]*work1[n]);
        work4[n+1] = fky[j]*(vg[ii][4]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkz[k]*work1[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v4y_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V4z

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkz[k]*(vg[ii][4]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkz[k]*work1[n]);
        work4[n+1] = fkz[k]*(vg[ii][4]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkz[k]*work1[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v4z_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V5x

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*(vg[ii][5]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkz[k]*work2[n]);
        work4[n+1] = fkx[i]*(vg[ii][5]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkz[k]*work2[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v5x_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V5y

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*(vg[ii][5]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkz[k]*work2[n]);
        work4[n+1] = fky[j]*(vg[ii][5]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkz[k]*work2[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v5y_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // V5z

  n = 0;
  ii = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkz[k]*(vg[ii][5]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]) + 2.0*fkz[k]*work2[n]);
        work4[n+1] = fkz[k]*(vg[ii][5]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]) + 2.0*fkz[k]*work2[n+1]);
        n += 2;
        ii++;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v5z_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

void PPPMDipole::fieldforce_ik_dipole()
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


  double **mu = atom->mu;
  double **x = atom->x;
  double **f = atom->f;
  double **t = atom->torque;

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

    // convert E-field to torque

    const double mufactor = qqrd2e * scale;
    f[i][0] += mufactor*(vxx*mu[i][0] + vxy*mu[i][1] + vxz*mu[i][2]);
    f[i][1] += mufactor*(vxy*mu[i][0] + vyy*mu[i][1] + vyz*mu[i][2]);
    f[i][2] += mufactor*(vxz*mu[i][0] + vyz*mu[i][1] + vzz*mu[i][2]);

    t[i][0] += mufactor*(mu[i][1]*ez - mu[i][2]*ey);
    t[i][1] += mufactor*(mu[i][2]*ex - mu[i][0]*ez);
    t[i][2] += mufactor*(mu[i][0]*ey - mu[i][1]*ex);
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get per-atom energy/virial
------------------------------------------------------------------------- */

void PPPMDipole::fieldforce_peratom_dipole()
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

  double **mu = atom->mu;
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

    if (eflag_atom) eatom[i] += mu[i][0]*ux + mu[i][1]*uy + mu[i][2]*uz;
    if (vflag_atom) {
      vatom[i][0] += mu[i][0]*v0x + mu[i][1]*v0y + mu[i][2]*v0z;
      vatom[i][1] += mu[i][0]*v1x + mu[i][1]*v1y + mu[i][2]*v1z;
      vatom[i][2] += mu[i][0]*v2x + mu[i][1]*v2y + mu[i][2]*v2z;
      vatom[i][3] += mu[i][0]*v3x + mu[i][1]*v3y + mu[i][2]*v3z;
      vatom[i][4] += mu[i][0]*v4x + mu[i][1]*v4y + mu[i][2]*v4z;
      vatom[i][5] += mu[i][0]*v5x + mu[i][1]*v5y + mu[i][2]*v5z;
    }
  }
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void PPPMDipole::pack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;

  if (flag == FORWARD_MU) {
    FFT_SCALAR *src_ux = &ux_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_uy = &uy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_uz = &uz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vxx = &vdxx_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vyy = &vdyy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vzz = &vdzz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vxy = &vdxy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vxz = &vdxz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vyz = &vdyz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = src_ux[list[i]];
      buf[n++] = src_uy[list[i]];
      buf[n++] = src_uz[list[i]];
      buf[n++] = src_vxx[list[i]];
      buf[n++] = src_vyy[list[i]];
      buf[n++] = src_vzz[list[i]];
      buf[n++] = src_vxy[list[i]];
      buf[n++] = src_vxz[list[i]];
      buf[n++] = src_vyz[list[i]];
    }
  } else if (flag == FORWARD_MU_PERATOM) {
    FFT_SCALAR *v0xsrc = &v0x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1xsrc = &v1x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2xsrc = &v2x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3xsrc = &v3x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4xsrc = &v4x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5xsrc = &v5x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0ysrc = &v0y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1ysrc = &v1y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2ysrc = &v2y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3ysrc = &v3y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4ysrc = &v4y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5ysrc = &v5y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0zsrc = &v0z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1zsrc = &v1z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2zsrc = &v2z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3zsrc = &v3z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4zsrc = &v4z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5zsrc = &v5z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = v0xsrc[list[i]];
      buf[n++] = v1xsrc[list[i]];
      buf[n++] = v2xsrc[list[i]];
      buf[n++] = v3xsrc[list[i]];
      buf[n++] = v4xsrc[list[i]];
      buf[n++] = v5xsrc[list[i]];
      buf[n++] = v0ysrc[list[i]];
      buf[n++] = v1ysrc[list[i]];
      buf[n++] = v2ysrc[list[i]];
      buf[n++] = v3ysrc[list[i]];
      buf[n++] = v4ysrc[list[i]];
      buf[n++] = v5ysrc[list[i]];
      buf[n++] = v0zsrc[list[i]];
      buf[n++] = v1zsrc[list[i]];
      buf[n++] = v2zsrc[list[i]];
      buf[n++] = v3zsrc[list[i]];
      buf[n++] = v4zsrc[list[i]];
      buf[n++] = v5zsrc[list[i]];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void PPPMDipole::unpack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;

  if (flag == FORWARD_MU) {
    FFT_SCALAR *dest_ux = &ux_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_uy = &uy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_uz = &uz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vxx = &vdxx_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vyy = &vdyy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vzz = &vdzz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vxy = &vdxy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vxz = &vdxz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vyz = &vdyz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      dest_ux[list[i]] = buf[n++];
      dest_uy[list[i]] = buf[n++];
      dest_uz[list[i]] = buf[n++];
      dest_vxx[list[i]] = buf[n++];
      dest_vyy[list[i]] = buf[n++];
      dest_vzz[list[i]] = buf[n++];
      dest_vxy[list[i]] = buf[n++];
      dest_vxz[list[i]] = buf[n++];
      dest_vyz[list[i]] = buf[n++];
    }
  } else if (flag == FORWARD_MU_PERATOM) {
    FFT_SCALAR *v0xsrc = &v0x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1xsrc = &v1x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2xsrc = &v2x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3xsrc = &v3x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4xsrc = &v4x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5xsrc = &v5x_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0ysrc = &v0y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1ysrc = &v1y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2ysrc = &v2y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3ysrc = &v3y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4ysrc = &v4y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5ysrc = &v5y_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0zsrc = &v0z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1zsrc = &v1z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2zsrc = &v2z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3zsrc = &v3z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4zsrc = &v4z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5zsrc = &v5z_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      v0xsrc[list[i]] = buf[n++];
      v1xsrc[list[i]] = buf[n++];
      v2xsrc[list[i]] = buf[n++];
      v3xsrc[list[i]] = buf[n++];
      v4xsrc[list[i]] = buf[n++];
      v5xsrc[list[i]] = buf[n++];
      v0ysrc[list[i]] = buf[n++];
      v1ysrc[list[i]] = buf[n++];
      v2ysrc[list[i]] = buf[n++];
      v3ysrc[list[i]] = buf[n++];
      v4ysrc[list[i]] = buf[n++];
      v5ysrc[list[i]] = buf[n++];
      v0zsrc[list[i]] = buf[n++];
      v1zsrc[list[i]] = buf[n++];
      v2zsrc[list[i]] = buf[n++];
      v3zsrc[list[i]] = buf[n++];
      v4zsrc[list[i]] = buf[n++];
      v5zsrc[list[i]] = buf[n++];
    }
  }
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void PPPMDipole::pack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;
  if (flag == REVERSE_MU) {
    FFT_SCALAR *src_dipole0 = &densityx_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_dipole1 = &densityy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_dipole2 = &densityz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = src_dipole0[list[i]];
      buf[n++] = src_dipole1[list[i]];
      buf[n++] = src_dipole2[list[i]];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void PPPMDipole::unpack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;
  if (flag == REVERSE_MU) {
    FFT_SCALAR *dest_dipole0 = &densityx_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_dipole1 = &densityy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_dipole2 = &densityz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      dest_dipole0[list[i]] += buf[n++];
      dest_dipole1[list[i]] += buf[n++];
      dest_dipole2[list[i]] += buf[n++];
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

void PPPMDipole::slabcorr()
{
  // compute local contribution to global dipole moment

  double dipole = 0.0;
  double **mu = atom->mu;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) dipole += mu[i][2];

  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant

  if (eflag_atom || fabs(qsum) > SMALL) {

    error->all(FLERR,"Cannot (yet) use kspace slab correction with "
      "long-range dipoles and non-neutral systems or per-atom energy");
  }

  // compute corrections

  const double e_slabcorr = MY_2PI*(dipole_all*dipole_all/12.0)/volume;
  const double qscale = qqrd2e * scale;

  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = qscale * MY_2PI/volume/12.0;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * mu[i][2]*dipole_all;
  }

  // add on torque corrections

  if (atom->torque) {
    double ffact = qscale * (-4.0*MY_PI/volume);
    double **mu = atom->mu;
    double **torque = atom->torque;
    for (int i = 0; i < nlocal; i++) {
      torque[i][0] += ffact * dipole_all * mu[i][1];
      torque[i][1] += -ffact * dipole_all * mu[i][0];
    }
  }
}

/* ----------------------------------------------------------------------
   perform and time the 1d FFTs required for N timesteps
------------------------------------------------------------------------- */

int PPPMDipole::timing_1d(int n, double &time1d)
{
  double time1,time2;

  for (int i = 0; i < 2*nfft_both; i++) work1[i] = ZEROF;

  MPI_Barrier(world);
  time1 = MPI_Wtime();

  for (int i = 0; i < n; i++) {
    fft1->timing1d(work1,nfft_both,1);
    fft1->timing1d(work1,nfft_both,1);
    fft1->timing1d(work1,nfft_both,1);
    fft2->timing1d(work1,nfft_both,-1);
    fft2->timing1d(work1,nfft_both,-1);
    fft2->timing1d(work1,nfft_both,-1);
    fft2->timing1d(work1,nfft_both,-1);
    fft2->timing1d(work1,nfft_both,-1);
    fft2->timing1d(work1,nfft_both,-1);
    fft2->timing1d(work1,nfft_both,-1);
    fft2->timing1d(work1,nfft_both,-1);
    fft2->timing1d(work1,nfft_both,-1);
  }

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  time1d = time2 - time1;

  return 12;
}

/* ----------------------------------------------------------------------
   perform and time the 3d FFTs required for N timesteps
------------------------------------------------------------------------- */

int PPPMDipole::timing_3d(int n, double &time3d)
{
  double time1,time2;

  for (int i = 0; i < 2*nfft_both; i++) work1[i] = ZEROF;

  MPI_Barrier(world);
  time1 = MPI_Wtime();

  for (int i = 0; i < n; i++) {
    fft1->compute(work1,work1,1);
    fft1->compute(work1,work1,1);
    fft1->compute(work1,work1,1);
    fft2->compute(work1,work1,-1);
    fft2->compute(work1,work1,-1);
    fft2->compute(work1,work1,-1);
    fft2->compute(work1,work1,-1);
    fft2->compute(work1,work1,-1);
    fft2->compute(work1,work1,-1);
    fft2->compute(work1,work1,-1);
    fft2->compute(work1,work1,-1);
    fft2->compute(work1,work1,-1);
  }

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  time3d = time2 - time1;

  return 12;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double PPPMDipole::memory_usage()
{
  double bytes = nmax*3 * sizeof(double);
  int nbrick = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
    (nzhi_out-nzlo_out+1);
  bytes += 6 * nfft_both * sizeof(double);   // vg
  bytes += nfft_both * sizeof(double);       // greensfn
  bytes += nfft_both*5 * sizeof(FFT_SCALAR); // work*2*2
  bytes += 9 * nbrick * sizeof(FFT_SCALAR);  // ubrick*3 + vdbrick*6
  bytes += nfft_both*7 * sizeof(FFT_SCALAR); // density_ffx*3 + work*2*2

  if (peratom_allocate_flag)
    bytes += 21 * nbrick * sizeof(FFT_SCALAR);

  if (cg_dipole) bytes += cg_dipole->memory_usage();
  if (cg_peratom_dipole) bytes += cg_peratom_dipole->memory_usage();

  return bytes;
}

/* ----------------------------------------------------------------------
   compute musum,musqsum,mu2
   called initially, when particle count changes, when dipoles are changed
------------------------------------------------------------------------- */

void PPPMDipole::musum_musq()
{
  const int nlocal = atom->nlocal;

  musum = musqsum = mu2 = 0.0;
  if (atom->mu_flag) {
    double** mu = atom->mu;
    double musum_local(0.0), musqsum_local(0.0);

    for (int i = 0; i < nlocal; i++) {
      musum_local += mu[i][0] + mu[i][1] + mu[i][2];
      musqsum_local += mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2];
    }

    MPI_Allreduce(&musum_local,&musum,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&musqsum_local,&musqsum,1,MPI_DOUBLE,MPI_SUM,world);

    mu2 = musqsum * force->qqrd2e;
  }

  if (mu2 == 0 && comm->me == 0)
    error->all(FLERR,"Using kspace solver PPPMDipole on system with no dipoles");
}
