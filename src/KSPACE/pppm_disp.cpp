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
   Contributing authors: Rolf Isele-Holder (Aachen University)
                         Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "pppm_disp.h"
#include "math_const.h"
#include "atom.h"
#include "comm.h"
#include "gridcomm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "domain.h"
#include "fft3d_wrap.h"
#include "remap_wrap.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXORDER   7
#define OFFSET 16384
#define SMALL 0.00001
#define LARGE 10000.0
#define EPS_HOC 1.0e-7

enum{GEOMETRIC,ARITHMETIC,SIXTHPOWER};
enum{REVERSE_RHO, REVERSE_RHO_G, REVERSE_RHO_A, REVERSE_RHO_NONE};
enum{FORWARD_IK, FORWARD_AD, FORWARD_IK_PERATOM, FORWARD_AD_PERATOM,
     FORWARD_IK_G, FORWARD_AD_G, FORWARD_IK_PERATOM_G, FORWARD_AD_PERATOM_G,
     FORWARD_IK_A, FORWARD_AD_A, FORWARD_IK_PERATOM_A, FORWARD_AD_PERATOM_A,
     FORWARD_IK_NONE, FORWARD_AD_NONE, FORWARD_IK_PERATOM_NONE, FORWARD_AD_PERATOM_NONE};


#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

PPPMDisp::PPPMDisp(LAMMPS *lmp, int narg, char **arg) : KSpace(lmp, narg, arg)
{
  if (narg < 1) error->all(FLERR,"Illegal kspace_style pppm/disp command");

  triclinic_support = 0;
  pppmflag = dispersionflag = 1;
  accuracy_relative = fabs(force->numeric(FLERR,arg[0]));
  
  nfactors = 3;
  factors = new int[nfactors];
  factors[0] = 2;
  factors[1] = 3;
  factors[2] = 5;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  csumflag = 0;
  B = NULL;
  cii = NULL;
  csumi = NULL;
  peratom_allocate_flag = 0;

  density_brick = vdx_brick = vdy_brick = vdz_brick = NULL;
  density_fft = NULL;
  u_brick = v0_brick = v1_brick = v2_brick = v3_brick = 
    v4_brick = v5_brick = NULL;

  density_brick_g = vdx_brick_g = vdy_brick_g = vdz_brick_g = NULL;
  density_fft_g = NULL;
  u_brick_g = v0_brick_g = v1_brick_g = v2_brick_g = v3_brick_g = 
    v4_brick_g = v5_brick_g = NULL;

  density_brick_a0 = vdx_brick_a0 = vdy_brick_a0 = vdz_brick_a0 = NULL;
  density_fft_a0 = NULL;
  u_brick_a0 = v0_brick_a0 = v1_brick_a0 = v2_brick_a0 = v3_brick_a0 = 
    v4_brick_a0 = v5_brick_a0 = NULL;

  density_brick_a1 = vdx_brick_a1 = vdy_brick_a1 = vdz_brick_a1 = NULL;
  density_fft_a1 = NULL;
  u_brick_a1 = v0_brick_a1 = v1_brick_a1 = v2_brick_a1 = v3_brick_a1 = 
    v4_brick_a1 = v5_brick_a1 = NULL;

  density_brick_a2 = vdx_brick_a2 = vdy_brick_a2 = vdz_brick_a2 = NULL;
  density_fft_a2 = NULL;
  u_brick_a2 = v0_brick_a2 = v1_brick_a2 = v2_brick_a2 = v3_brick_a2 = 
    v4_brick_a2 = v5_brick_a2 = NULL;

  density_brick_a3 = vdx_brick_a3 = vdy_brick_a3 = vdz_brick_a3 = NULL;
  density_fft_a3 = NULL;
  u_brick_a3 = v0_brick_a3 = v1_brick_a3 = v2_brick_a3 = v3_brick_a3 = 
    v4_brick_a3 = v5_brick_a3 = NULL;

  density_brick_a4 = vdx_brick_a4 = vdy_brick_a4 = vdz_brick_a4 = NULL;
  density_fft_a4 = NULL;
  u_brick_a4 = v0_brick_a4 = v1_brick_a4 = v2_brick_a4 = v3_brick_a4 = 
    v4_brick_a4 = v5_brick_a4 = NULL;

  density_brick_a5 = vdx_brick_a5 = vdy_brick_a5 = vdz_brick_a5 = NULL;
  density_fft_a5 = NULL;
  u_brick_a5 = v0_brick_a5 = v1_brick_a5 = v2_brick_a5 = v3_brick_a5 = 
    v4_brick_a5 = v5_brick_a5 = NULL;

  density_brick_a6 = vdx_brick_a6 = vdy_brick_a6 = vdz_brick_a6 = NULL;
  density_fft_a6 = NULL;
  u_brick_a6 = v0_brick_a6 = v1_brick_a6 = v2_brick_a6 = v3_brick_a6 = 
    v4_brick_a6 = v5_brick_a6 = NULL;

  density_brick_none = vdx_brick_none = vdy_brick_none = vdz_brick_none = NULL;
  density_fft_none = NULL;
  u_brick_none = v0_brick_none = v1_brick_none = v2_brick_none = v3_brick_none = 
    v4_brick_none = v5_brick_none = NULL;

  greensfn = NULL;
  greensfn_6 = NULL;
  work1 = work2 = NULL;
  work1_6 = work2_6 = NULL;
  vg = NULL;
  vg2 = NULL;
  vg_6 = NULL;
  vg2_6 = NULL;
  fkx = fky = fkz = NULL;
  fkx2 = fky2 = fkz2 = NULL;
  fkx_6 = fky_6 = fkz_6 = NULL;
  fkx2_6 = fky2_6 = fkz2_6 = NULL;

  sf_precoeff1 = sf_precoeff2 = sf_precoeff3 = sf_precoeff4 = 
    sf_precoeff5 = sf_precoeff6 = NULL;
  sf_precoeff1_6 = sf_precoeff2_6 = sf_precoeff3_6 = sf_precoeff4_6 = 
    sf_precoeff5_6 = sf_precoeff6_6 = NULL;

  gf_b = NULL;
  gf_b_6 = NULL;
  rho1d = rho_coeff = NULL;
  drho1d = drho_coeff = NULL;
  rho1d_6 = rho_coeff_6 = NULL;
  drho1d_6 = drho_coeff_6 = NULL;
  fft1 = fft2 = NULL;
  fft1_6 = fft2_6 = NULL;
  remap = NULL;
  remap_6 = NULL;

  nmax = 0;
  part2grid = NULL;
  part2grid_6 = NULL;

  cg = NULL;
  cg_peratom = NULL;
  cg_6 = NULL;
  cg_peratom_6 = NULL;

  memset(function, 0, EWALD_FUNCS*sizeof(int));
}

/* ----------------------------------------------------------------------
   free all memory 
------------------------------------------------------------------------- */

PPPMDisp::~PPPMDisp()
{
  delete [] factors;
  delete [] B;
  B = NULL;
  delete [] cii;
  cii = NULL;
  delete [] csumi;
  csumi = NULL;
  deallocate();
  deallocate_peratom();
  memory->destroy(part2grid);
  memory->destroy(part2grid_6);
  part2grid = part2grid_6 = NULL;
}

/* ----------------------------------------------------------------------
   called once before run 
------------------------------------------------------------------------- */

void PPPMDisp::init()
{
  if (me == 0) {
    if (screen) fprintf(screen,"PPPMDisp initialization ...\n");
    if (logfile) fprintf(logfile,"PPPMDisp initialization ...\n");
  }

  triclinic_check();
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use PPPMDisp with 2d simulation");
  if (comm->style != 0) 
    error->universe_all(FLERR,"PPPMDisp can only currently be used with "
                        "comm_style brick");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use nonperiodic boundaries with PPPMDisp");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 || 
	domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab PPPMDisp");
  }
 
  if (order > MAXORDER || order_6 > MAXORDER) {
    char str[128];
    sprintf(str,"PPPMDisp coulomb order cannot be greater than %d",MAXORDER);
    error->all(FLERR,str);
  }

  // free all arrays previously allocated

  deallocate();
  deallocate_peratom(); 

  // check whether cutoff and pair style are set

  triclinic = domain->triclinic;
  pair_check();

  int tmp;
  Pair *pair = force->pair;
  int *ptr = pair ? (int *) pair->extract("ewald_order",tmp) : NULL;
  double *p_cutoff = pair ? (double *) pair->extract("cut_coul",tmp) : NULL;
  double *p_cutoff_lj = pair ? (double *) pair->extract("cut_LJ",tmp) : NULL;
  if (!(ptr||*p_cutoff||*p_cutoff_lj)) 
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;
  cutoff_lj = *p_cutoff_lj;

  double tmp2;
  MPI_Allreduce(&cutoff, &tmp2,1,MPI_DOUBLE,MPI_SUM,world); 

  // check out which types of potentials will have to be calculated

  int ewald_order = ptr ? *((int *) ptr) : 1<<1;
  int ewald_mix = ptr ? *((int *) pair->extract("ewald_mix",tmp)) : GEOMETRIC;
  memset(function, 0, EWALD_FUNCS*sizeof(int));
  for (int i=0; i<=EWALD_MAXORDER; ++i)			// transcribe order
    if (ewald_order&(1<<i)) {				// from pair_style
      int  k=0;
      char str[128];
      switch (i) {
	case 1:
	  k = 0; break;
	case 6:
	  if ((ewald_mix==GEOMETRIC || ewald_mix==SIXTHPOWER || 
               mixflag == 1) && mixflag!= 2) { k = 1; break; }
	  else if (ewald_mix==ARITHMETIC && mixflag!=2) { k = 2; break; }
	  else if (mixflag == 2) { k = 3; break; }
	default:
	  sprintf(str, "Unsupported order in kspace_style "
                  "pppm/disp, pair_style %s", force->pair_style);
	  error->all(FLERR,str);
      }
      function[k] = 1;
    }
 

  // warn, if function[0] is not set but charge attribute is set!

  if (!function[0] && atom->q_flag && me == 0) {
    char str[128];
    sprintf(str, "Charges are set, but coulombic solver is not used");
    error->warning(FLERR, str);
  }

  // compute qsum & qsqsum, if function[0] is set, warn if not charge-neutral

  scale = 1.0;
  qqrd2e = force->qqrd2e;
  natoms_original = atom->natoms;
 
  if (function[0]) qsum_qsq(0);

  // if kspace is TIP4P, extract TIP4P params from pair style
  // bond/angle are not yet init(), so insure equilibrium request is valid

  qdist = 0.0;
 
  if (tip4pflag) {
    int itmp;
    double *p_qdist = (double *) force->pair->extract("qdist",itmp);
    int *p_typeO = (int *) force->pair->extract("typeO",itmp);
    int *p_typeH = (int *) force->pair->extract("typeH",itmp);
    int *p_typeA = (int *) force->pair->extract("typeA",itmp);
    int *p_typeB = (int *) force->pair->extract("typeB",itmp);
    if (!p_qdist || !p_typeO || !p_typeH || !p_typeA || !p_typeB)
      error->all(FLERR,"KSpace style is incompatible with Pair style");
    qdist = *p_qdist;
    typeO = *p_typeO;
    typeH = *p_typeH;
    int typeA = *p_typeA;
    int typeB = *p_typeB;

    if (force->angle == NULL || force->bond == NULL)
      error->all(FLERR,"Bond and angle potentials must be defined for TIP4P");
    if (typeA < 1 || typeA > atom->nangletypes || 
	force->angle->setflag[typeA] == 0)
      error->all(FLERR,"Bad TIP4P angle type for PPPMDisp/TIP4P");
    if (typeB < 1 || typeB > atom->nbondtypes || 
	force->bond->setflag[typeB] == 0)
      error->all(FLERR,"Bad TIP4P bond type for PPPMDisp/TIP4P");
    double theta = force->angle->equilibrium_angle(typeA);
    double blen = force->bond->equilibrium_distance(typeB);
    alpha = qdist / (cos(0.5*theta) * blen);
  }

  // initialize the pair style to get the coefficients

  neighrequest_flag = 0;
  pair->init();
  neighrequest_flag = 1;
  init_coeffs();

  //if g_ewald and g_ewald_6 have not been specified, set some initial value
  //  to avoid problems when calculating the energies!

  if (!gewaldflag) g_ewald = 1;
  if (!gewaldflag_6) g_ewald_6 = 1;

  // set accuracy (force units) from accuracy_relative or accuracy_absolute
  
  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;

  int (*procneigh)[2] = comm->procneigh;

  int iteration = 0;
  if (function[0]) {
    GridComm *cgtmp = NULL;
    while (order >= minorder) {

      if (iteration && me == 0)
          error->warning(FLERR,"Reducing PPPMDisp Coulomb order "
                         "b/c stencil extends beyond neighbor processor");
      iteration++;

      // set grid for dispersion interaction and coulomb interactions
 
      set_grid();

      if (nx_pppm >= OFFSET || ny_pppm >= OFFSET || nz_pppm >= OFFSET)
      error->all(FLERR,"PPPMDisp Coulomb grid is too large");

      set_fft_parameters(nx_pppm, ny_pppm, nz_pppm,
                         nxlo_fft, nylo_fft, nzlo_fft,
                         nxhi_fft, nyhi_fft, nzhi_fft,
                         nxlo_in, nylo_in, nzlo_in,
                         nxhi_in, nyhi_in, nzhi_in,
                         nxlo_out, nylo_out, nzlo_out,
                         nxhi_out, nyhi_out, nzhi_out,
                         nlower, nupper,
                         ngrid, nfft, nfft_both,
                         shift, shiftone, order);

      if (overlap_allowed) break;

      cgtmp = new GridComm(lmp, world,1,1,
                           nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                           nxlo_out,nxhi_out,nylo_out,nyhi_out,
                           nzlo_out,nzhi_out,
                           procneigh[0][0],procneigh[0][1],procneigh[1][0],
                           procneigh[1][1],procneigh[2][0],procneigh[2][1]);
      cgtmp->ghost_notify();
      if (!cgtmp->ghost_overlap()) break;
      delete cgtmp;

      order--;
    }

    if (order < minorder)
      error->all(FLERR,
                 "Coulomb PPPMDisp order has been reduced below minorder");
    if (cgtmp) delete cgtmp;

    // adjust g_ewald
  
    if (!gewaldflag) adjust_gewald();

    // calculate the final accuracy
  
    double acc = final_accuracy();
  
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
        fprintf(screen,"  Coulomb G vector (1/distance)= %g\n",g_ewald);
        fprintf(screen,"  Coulomb grid = %d %d %d\n",nx_pppm,ny_pppm,nz_pppm);
        fprintf(screen,"  Coulomb stencil order = %d\n",order);
        fprintf(screen,"  Coulomb estimated absolute RMS force accuracy = %g\n",
                acc);
        fprintf(screen,"  Coulomb estimated relative force accuracy = %g\n",
                acc/two_charge_force);
        fprintf(screen,"  using %s precision FFTs\n",fft_prec);
        fprintf(screen,"  3d grid and FFT values/proc = %d %d\n",
		ngrid_max, nfft_both_max);
      }
      if (logfile) {
        fprintf(logfile,"  Coulomb G vector (1/distance) = %g\n",g_ewald);
        fprintf(logfile,"  Coulomb grid = %d %d %d\n",nx_pppm,ny_pppm,nz_pppm);
        fprintf(logfile,"  Coulomb stencil order = %d\n",order);
        fprintf(logfile,
                "  Coulomb estimated absolute RMS force accuracy = %g\n",
                acc);
        fprintf(logfile,"  Coulomb estimated relative force accuracy = %g\n",
                acc/two_charge_force);
        fprintf(logfile,"  using %s precision FFTs\n",fft_prec);
        fprintf(logfile,"  3d grid and FFT values/proc = %d %d\n",
		ngrid_max, nfft_both_max);
      }
    }
  }

  iteration = 0;
  if (function[1] + function[2] + function[3]) {
    GridComm *cgtmp = NULL;
    while (order_6 >= minorder) {

      if (iteration && me == 0)
          error->warning(FLERR,"Reducing PPPMDisp dispersion order "
                         "b/c stencil extends beyond neighbor processor");
      iteration++;

      set_grid_6();
   
      if (nx_pppm_6 >= OFFSET || ny_pppm_6 >= OFFSET || nz_pppm_6 >= OFFSET)
      error->all(FLERR,"PPPMDisp Dispersion grid is too large");

      set_fft_parameters(nx_pppm_6, ny_pppm_6, nz_pppm_6,
                         nxlo_fft_6, nylo_fft_6, nzlo_fft_6,
                         nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                         nxlo_in_6, nylo_in_6, nzlo_in_6,
                         nxhi_in_6, nyhi_in_6, nzhi_in_6,
                         nxlo_out_6, nylo_out_6, nzlo_out_6,
                         nxhi_out_6, nyhi_out_6, nzhi_out_6,
                         nlower_6, nupper_6,
                         ngrid_6, nfft_6, nfft_both_6,
                         shift_6, shiftone_6, order_6);

      if (overlap_allowed) break;

      cgtmp = new GridComm(lmp,world,1,1,
                           nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,
                           nzlo_in_6,nzhi_in_6,
                           nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,
                           nzlo_out_6,nzhi_out_6,
                           procneigh[0][0],procneigh[0][1],procneigh[1][0],
                           procneigh[1][1],procneigh[2][0],procneigh[2][1]);
      cgtmp->ghost_notify();
      if (!cgtmp->ghost_overlap()) break;
      delete cgtmp;
      order_6--;
    }

    if (order_6 < minorder) 
      error->all(FLERR,"Dispersion PPPMDisp order has been "
                 "reduced below minorder");
    if (cgtmp) delete cgtmp;

    // adjust g_ewald_6

    if (!gewaldflag_6 && accuracy_kspace_6 == accuracy_real_6) 
      adjust_gewald_6();

    // calculate the final accuracy

    double acc, acc_real, acc_kspace;
    final_accuracy_6(acc, acc_real, acc_kspace);


    // print stats

    int ngrid_max,nfft_both_max;
    MPI_Allreduce(&ngrid_6,&ngrid_max,1,MPI_INT,MPI_MAX,world);
    MPI_Allreduce(&nfft_both_6,&nfft_both_max,1,MPI_INT,MPI_MAX,world);

    if (me == 0) {
    #ifdef FFT_SINGLE
      const char fft_prec[] = "single";
    #else
      const char fft_prec[] = "double";
    #endif
  
      if (screen) {
        fprintf(screen,"  Dispersion G vector (1/distance)= %g\n",g_ewald_6);
        fprintf(screen,"  Dispersion grid = %d %d %d\n",
                nx_pppm_6,ny_pppm_6,nz_pppm_6);
        fprintf(screen,"  Dispersion stencil order = %d\n",order_6);
        fprintf(screen,"  Dispersion estimated absolute "
                "RMS force accuracy = %g\n",acc);
        fprintf(screen,"  Dispersion estimated absolute "
                "real space RMS force accuracy = %g\n",acc_real);
        fprintf(screen,"  Dispersion estimated absolute "
                "kspace RMS force accuracy = %g\n",acc_kspace);
        fprintf(screen,"  Dispersion estimated relative force accuracy = %g\n",
                acc/two_charge_force);
        fprintf(screen,"  using %s precision FFTs\n",fft_prec);
        fprintf(screen,"  3d grid and FFT values/proc dispersion = %d %d\n",
                          ngrid_max,nfft_both_max);
      }
      if (logfile) {
        fprintf(logfile,"  Dispersion G vector (1/distance) = %g\n",g_ewald_6);
        fprintf(logfile,"  Dispersion grid = %d %d %d\n",
                nx_pppm_6,ny_pppm_6,nz_pppm_6);
        fprintf(logfile,"  Dispersion stencil order = %d\n",order_6);
        fprintf(logfile,"  Dispersion estimated absolute "
                "RMS force accuracy = %g\n",acc);
        fprintf(logfile,"  Dispersion estimated absolute "
                "real space RMS force accuracy = %g\n",acc_real);
        fprintf(logfile,"  Dispersion estimated absolute "
                "kspace RMS force accuracy = %g\n",acc_kspace);
        fprintf(logfile,"  Disperion estimated relative force accuracy = %g\n",
                acc/two_charge_force);
        fprintf(logfile,"  using %s precision FFTs\n",fft_prec);
        fprintf(logfile,"  3d grid and FFT values/proc dispersion = %d %d\n",
                           ngrid_max,nfft_both_max);
      }
    }
  }

  // allocate K-space dependent memory

  allocate();

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients

  if (function[0]) {
    compute_gf_denom(gf_b, order);
    compute_rho_coeff(rho_coeff, drho_coeff, order);
    cg->ghost_notify();
    cg->setup();
    if (differentiation_flag == 1)
      compute_sf_precoeff(nx_pppm, ny_pppm, nz_pppm, order,
                          nxlo_fft, nylo_fft, nzlo_fft, 
                          nxhi_fft, nyhi_fft, nzhi_fft,
                          sf_precoeff1, sf_precoeff2, sf_precoeff3,
                          sf_precoeff4, sf_precoeff5, sf_precoeff6);
  }
  if (function[1] + function[2] + function[3]) {
    compute_gf_denom(gf_b_6, order_6);
    compute_rho_coeff(rho_coeff_6, drho_coeff_6, order_6);
    cg_6->ghost_notify();
    cg_6->setup();
    if (differentiation_flag == 1)
      compute_sf_precoeff(nx_pppm_6, ny_pppm_6, nz_pppm_6, order_6,
                          nxlo_fft_6, nylo_fft_6, nzlo_fft_6, 
                          nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                          sf_precoeff1_6, sf_precoeff2_6, sf_precoeff3_6,
                          sf_precoeff4_6, sf_precoeff5_6, sf_precoeff6_6);
  }

}

/* ----------------------------------------------------------------------
   adjust PPPM coeffs, called initially and whenever volume has changed 
------------------------------------------------------------------------- */

void PPPMDisp::setup()
{
  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPM
  // z dimension for 3d PPPM is zprd since slab_volfactor = 1.0

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

 // compute fkx,fky,fkz for my FFT grid pts

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  //compute the virial coefficients and green functions
  if (function[0]){

    delxinv = nx_pppm/xprd;
    delyinv = ny_pppm/yprd;
    delzinv = nz_pppm/zprd_slab;

    delvolinv = delxinv*delyinv*delzinv;

    double per;
    int i, j, k, n;

    for (i = nxlo_fft; i <= nxhi_fft; i++) {
      per = i - nx_pppm*(2*i/nx_pppm);
      fkx[i] = unitkx*per;
      j = (nx_pppm - i) % nx_pppm;
      per = j - nx_pppm*(2*j/nx_pppm);
      fkx2[i] = unitkx*per;
    }

    for (i = nylo_fft; i <= nyhi_fft; i++) {
      per = i - ny_pppm*(2*i/ny_pppm);
      fky[i] = unitky*per;
      j = (ny_pppm - i) % ny_pppm;
      per = j - ny_pppm*(2*j/ny_pppm);
      fky2[i] = unitky*per;
    }

    for (i = nzlo_fft; i <= nzhi_fft; i++) {
      per = i - nz_pppm*(2*i/nz_pppm);
      fkz[i] = unitkz*per;
      j = (nz_pppm - i) % nz_pppm;
      per = j - nz_pppm*(2*j/nz_pppm);
      fkz2[i] = unitkz*per;
    }

    double sqk,vterm;
    double gew2inv = 1/(g_ewald*g_ewald);
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
	    vterm = -2.0 * (1.0/sqk + 0.25*gew2inv);
	    vg[n][0] = 1.0 + vterm*fkx[i]*fkx[i];
	    vg[n][1] = 1.0 + vterm*fky[j]*fky[j];
	    vg[n][2] = 1.0 + vterm*fkz[k]*fkz[k];
	    vg[n][3] = vterm*fkx[i]*fky[j];
	    vg[n][4] = vterm*fkx[i]*fkz[k];
	    vg[n][5] = vterm*fky[j]*fkz[k];
            vg2[n][0] = vterm*0.5*(fkx[i]*fky[j] + fkx2[i]*fky2[j]);
            vg2[n][1] = vterm*0.5*(fkx[i]*fkz[k] + fkx2[i]*fkz2[k]);
            vg2[n][2] = vterm*0.5*(fky[j]*fkz[k] + fky2[j]*fkz2[k]);
  	  }
	  n++;
        }
      }
    }
    compute_gf();
    if (differentiation_flag == 1) compute_sf_coeff();
  }

  if (function[1] + function[2] + function[3]) {
    delxinv_6 = nx_pppm_6/xprd;
    delyinv_6 = ny_pppm_6/yprd;
    delzinv_6 = nz_pppm_6/zprd_slab;
    delvolinv_6 = delxinv_6*delyinv_6*delzinv_6;

    double per;
    int i, j, k, n;
    for (i = nxlo_fft_6; i <= nxhi_fft_6; i++) {
      per = i - nx_pppm_6*(2*i/nx_pppm_6);
      fkx_6[i] = unitkx*per;
      j = (nx_pppm_6 - i) % nx_pppm_6;
      per = j - nx_pppm_6*(2*j/nx_pppm_6);
      fkx2_6[i] = unitkx*per;
    }
    for (i = nylo_fft_6; i <= nyhi_fft_6; i++) {
      per = i - ny_pppm_6*(2*i/ny_pppm_6);
      fky_6[i] = unitky*per;
      j = (ny_pppm_6 - i) % ny_pppm_6;
      per = j - ny_pppm_6*(2*j/ny_pppm_6);
      fky2_6[i] = unitky*per;
    }
    for (i = nzlo_fft_6; i <= nzhi_fft_6; i++) {
      per = i - nz_pppm_6*(2*i/nz_pppm_6);
      fkz_6[i] = unitkz*per;
      j = (nz_pppm_6 - i) % nz_pppm_6;
      per = j - nz_pppm_6*(2*j/nz_pppm_6);
      fkz2_6[i] = unitkz*per;
    }
    double sqk,vterm;
    long double erft, expt,nom, denom;
    long double b, bs, bt;
    double rtpi = sqrt(MY_PI);
    double gewinv = 1/g_ewald_6;
    n = 0;
    for (k = nzlo_fft_6; k <= nzhi_fft_6; k++) {
      for (j = nylo_fft_6; j <= nyhi_fft_6; j++) {
        for (i = nxlo_fft_6; i <= nxhi_fft_6; i++) {
	  sqk = fkx_6[i]*fkx_6[i] + fky_6[j]*fky_6[j] + fkz_6[k]*fkz_6[k];
	  if (sqk == 0.0) {
	    vg_6[n][0] = 0.0;
	    vg_6[n][1] = 0.0;
	    vg_6[n][2] = 0.0;
	    vg_6[n][3] = 0.0;
	    vg_6[n][4] = 0.0;
	    vg_6[n][5] = 0.0;
	  } else {
            b = 0.5*sqrt(sqk)*gewinv;
            bs = b*b;
            bt = bs*b;
            erft = 2*bt*rtpi*erfc((double) b);
            expt = exp(-bs);
            nom = erft - 2*bs*expt;
            denom = nom + expt;
            if (denom == 0) vterm = 3.0/sqk;
            else vterm = 3.0*nom/(sqk*denom);
	    vg_6[n][0] = 1.0 + vterm*fkx_6[i]*fkx_6[i];
	    vg_6[n][1] = 1.0 + vterm*fky_6[j]*fky_6[j];
	    vg_6[n][2] = 1.0 + vterm*fkz_6[k]*fkz_6[k];
	    vg_6[n][3] = vterm*fkx_6[i]*fky_6[j];
	    vg_6[n][4] = vterm*fkx_6[i]*fkz_6[k];
	    vg_6[n][5] = vterm*fky_6[j]*fkz_6[k];
            vg2_6[n][0] = vterm*0.5*(fkx_6[i]*fky_6[j] + fkx2_6[i]*fky2_6[j]);
            vg2_6[n][1] = vterm*0.5*(fkx_6[i]*fkz_6[k] + fkx2_6[i]*fkz2_6[k]);
            vg2_6[n][2] = vterm*0.5*(fky_6[j]*fkz_6[k] + fky2_6[j]*fkz2_6[k]);
	  }
	  n++;
        }
      }
    }
    compute_gf_6();
    if (differentiation_flag == 1) compute_sf_coeff_6();
  }
}

/* ----------------------------------------------------------------------
   reset local grid arrays and communication stencils
   called by fix balance b/c it changed sizes of processor sub-domains
------------------------------------------------------------------------- */

void PPPMDisp::setup_grid()
{
  // free all arrays previously allocated

  deallocate();
  deallocate_peratom();

  // reset portion of global grid that each proc owns

  if (function[0])
    set_fft_parameters(nx_pppm, ny_pppm, nz_pppm,
                       nxlo_fft, nylo_fft, nzlo_fft,
                       nxhi_fft, nyhi_fft, nzhi_fft,
                       nxlo_in, nylo_in, nzlo_in,
                       nxhi_in, nyhi_in, nzhi_in,
                       nxlo_out, nylo_out, nzlo_out,
                       nxhi_out, nyhi_out, nzhi_out,
                       nlower, nupper,
                       ngrid, nfft, nfft_both,
                       shift, shiftone, order);

  if (function[1] + function[2] + function[3])
    set_fft_parameters(nx_pppm_6, ny_pppm_6, nz_pppm_6,
                       nxlo_fft_6, nylo_fft_6, nzlo_fft_6,
                       nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                       nxlo_in_6, nylo_in_6, nzlo_in_6,
                       nxhi_in_6, nyhi_in_6, nzhi_in_6,
                       nxlo_out_6, nylo_out_6, nzlo_out_6,
                       nxhi_out_6, nyhi_out_6, nzhi_out_6,
                       nlower_6, nupper_6,
                       ngrid_6, nfft_6, nfft_both_6,
                       shift_6, shiftone_6, order_6);

  // reallocate K-space dependent memory
  // check if grid communication is now overlapping if not allowed
  // don't invoke allocate_peratom(), compute() will allocate when needed

  allocate();

  if (function[0]) {
    cg->ghost_notify();
    if (overlap_allowed == 0 && cg->ghost_overlap())
      error->all(FLERR,"PPPM grid stencil extends "
                 "beyond nearest neighbor processor");
    cg->setup();
  }
  if (function[1] + function[2] + function[3]) {
    cg_6->ghost_notify();
    if (overlap_allowed == 0 && cg_6->ghost_overlap())
      error->all(FLERR,"PPPM grid stencil extends "
                 "beyond nearest neighbor processor");
    cg_6->setup();
  }

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients

  if (function[0]) {
    compute_gf_denom(gf_b, order);
    compute_rho_coeff(rho_coeff, drho_coeff, order);
    if (differentiation_flag == 1) 
      compute_sf_precoeff(nx_pppm, ny_pppm, nz_pppm, order,
                          nxlo_fft, nylo_fft, nzlo_fft, 
                          nxhi_fft, nyhi_fft, nzhi_fft,
                          sf_precoeff1, sf_precoeff2, sf_precoeff3,
                          sf_precoeff4, sf_precoeff5, sf_precoeff6);
  }
  if (function[1] + function[2] + function[3]) {
    compute_gf_denom(gf_b_6, order_6);
    compute_rho_coeff(rho_coeff_6, drho_coeff_6, order_6);
    if (differentiation_flag == 1)
      compute_sf_precoeff(nx_pppm_6, ny_pppm_6, nz_pppm_6, order_6,
                          nxlo_fft_6, nylo_fft_6, nzlo_fft_6, 
                          nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                          sf_precoeff1_6, sf_precoeff2_6, sf_precoeff3_6,
                          sf_precoeff4_6, sf_precoeff5_6, sf_precoeff6_6);
  }

  // pre-compute volume-dependent coeffs

  setup();
}

/* ----------------------------------------------------------------------
   compute the PPPM long-range force, energy, virial 
------------------------------------------------------------------------- */

void PPPMDisp::compute(int eflag, int vflag)
{

  int i;
  // convert atoms from box to lamda coords

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global = 
	 eflag_atom = vflag_atom = 0;

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    if (function[0]) {
      cg_peratom->ghost_notify();
      cg_peratom->setup();
    }
    if (function[1] + function[2] + function[3]) {
      cg_peratom_6->ghost_notify();
      cg_peratom_6->setup();
    }
    peratom_allocate_flag = 1;
  }
  
  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }
  // extend size of per-atom arrays if necessary

  if (atom->nlocal > nmax) {

    if (function[0]) memory->destroy(part2grid);
    if (function[1] + function[2] + function[3]) memory->destroy(part2grid_6);
    nmax = atom->nmax;
    if (function[0]) memory->create(part2grid,nmax,3,"pppm/disp:part2grid");
    if (function[1] + function[2] + function[3]) 
      memory->create(part2grid_6,nmax,3,"pppm/disp:part2grid_6");
  }


  energy = 0.0;
  energy_1 = 0.0;
  energy_6 = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial_6[i] = virial_1[i] = 0.0;

  // find grid points for all my particles
  // distribute partcles' charges/dispersion coefficients on the grid
  // communication between processors and remapping two fft
  // Solution of poissons equation in k-space and backtransformation
  // communication between processors
  // calculation of forces

  if (function[0]) {

    //perfrom calculations for coulomb interactions only

    particle_map_c(delxinv, delyinv, delzinv, shift, part2grid, nupper, nlower,
                 nxlo_out, nylo_out, nzlo_out, nxhi_out, nyhi_out, nzhi_out);

    make_rho_c();

    cg->reverse_comm(this,REVERSE_RHO);
 
    brick2fft(nxlo_in, nylo_in, nzlo_in, nxhi_in, nyhi_in, nzhi_in,
	      density_brick, density_fft, work1,remap); 
 
    if (differentiation_flag == 1) {

      poisson_ad(work1, work2, density_fft, fft1, fft2,
                 nx_pppm, ny_pppm, nz_pppm, nfft,
                 nxlo_fft, nylo_fft, nzlo_fft, nxhi_fft, nyhi_fft, nzhi_fft,
                 nxlo_in, nylo_in, nzlo_in, nxhi_in, nyhi_in, nzhi_in,
                 energy_1, greensfn, 
                 virial_1, vg,vg2,
                 u_brick, v0_brick, v1_brick, v2_brick, v3_brick, v4_brick, v5_brick);

      cg->forward_comm(this,FORWARD_AD);

      fieldforce_c_ad(); 

      if (vflag_atom) cg_peratom->forward_comm(this, FORWARD_AD_PERATOM);

    } else {
      poisson_ik(work1, work2, density_fft, fft1, fft2,
                 nx_pppm, ny_pppm, nz_pppm, nfft,
                 nxlo_fft, nylo_fft, nzlo_fft, nxhi_fft, nyhi_fft, nzhi_fft,
                 nxlo_in, nylo_in, nzlo_in, nxhi_in, nyhi_in, nzhi_in,
                 energy_1, greensfn, 
	         fkx, fky, fkz,fkx2, fky2, fkz2,
                 vdx_brick, vdy_brick, vdz_brick, virial_1, vg,vg2,
                 u_brick, v0_brick, v1_brick, v2_brick, v3_brick, v4_brick, v5_brick);

      cg->forward_comm(this, FORWARD_IK);

      fieldforce_c_ik(); 

      if (evflag_atom) cg_peratom->forward_comm(this, FORWARD_IK_PERATOM);
    }
    if (evflag_atom) fieldforce_c_peratom();
  }

  if (function[1]) {
    //perfrom calculations for geometric mixing
    particle_map(delxinv_6, delyinv_6, delzinv_6, shift_6, part2grid_6, nupper_6, nlower_6,
                 nxlo_out_6, nylo_out_6, nzlo_out_6, nxhi_out_6, nyhi_out_6, nzhi_out_6);
    make_rho_g();


    cg_6->reverse_comm(this, REVERSE_RHO_G);

    brick2fft(nxlo_in_6, nylo_in_6, nzlo_in_6, nxhi_in_6, nyhi_in_6, nzhi_in_6,
	      density_brick_g, density_fft_g, work1_6,remap_6);
 
    if (differentiation_flag == 1) {

      poisson_ad(work1_6, work2_6, density_fft_g, fft1_6, fft2_6,
                 nx_pppm_6, ny_pppm_6, nz_pppm_6, nfft_6,
                 nxlo_fft_6, nylo_fft_6, nzlo_fft_6, nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                 nxlo_in_6, nylo_in_6, nzlo_in_6, nxhi_in_6, nyhi_in_6, nzhi_in_6,
                 energy_6, greensfn_6, 
                 virial_6, vg_6, vg2_6,
                 u_brick_g, v0_brick_g, v1_brick_g, v2_brick_g, v3_brick_g, v4_brick_g, v5_brick_g);

      cg_6->forward_comm(this,FORWARD_AD_G);

      fieldforce_g_ad();

      if (vflag_atom) cg_peratom_6->forward_comm(this,FORWARD_AD_PERATOM_G);

    } else {
      poisson_ik(work1_6, work2_6, density_fft_g, fft1_6, fft2_6,
                 nx_pppm_6, ny_pppm_6, nz_pppm_6, nfft_6,
                 nxlo_fft_6, nylo_fft_6, nzlo_fft_6, nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                 nxlo_in_6, nylo_in_6, nzlo_in_6, nxhi_in_6, nyhi_in_6, nzhi_in_6,
                 energy_6, greensfn_6, 
	         fkx_6, fky_6, fkz_6,fkx2_6, fky2_6, fkz2_6,
                 vdx_brick_g, vdy_brick_g, vdz_brick_g, virial_6, vg_6, vg2_6,
                 u_brick_g, v0_brick_g, v1_brick_g, v2_brick_g, v3_brick_g, v4_brick_g, v5_brick_g);
 
      cg_6->forward_comm(this,FORWARD_IK_G);
 
      fieldforce_g_ik();


      if (evflag_atom) cg_peratom_6->forward_comm(this, FORWARD_IK_PERATOM_G);
    }
    if (evflag_atom) fieldforce_g_peratom();
  }

  if (function[2]) {
    //perform calculations for arithmetic mixing
    particle_map(delxinv_6, delyinv_6, delzinv_6, shift_6, part2grid_6, nupper_6, nlower_6,
                 nxlo_out_6, nylo_out_6, nzlo_out_6, nxhi_out_6, nyhi_out_6, nzhi_out_6);
    make_rho_a();

    cg_6->reverse_comm(this, REVERSE_RHO_A);

    brick2fft_a();

    if ( differentiation_flag == 1) {

      poisson_ad(work1_6, work2_6, density_fft_a3, fft1_6, fft2_6,
                 nx_pppm_6, ny_pppm_6, nz_pppm_6, nfft_6,
                 nxlo_fft_6, nylo_fft_6, nzlo_fft_6, nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                 nxlo_in_6, nylo_in_6, nzlo_in_6, nxhi_in_6, nyhi_in_6, nzhi_in_6,
                 energy_6, greensfn_6, 
                 virial_6, vg_6, vg2_6,
                 u_brick_a3, v0_brick_a3, v1_brick_a3, v2_brick_a3, v3_brick_a3, v4_brick_a3, v5_brick_a3);
      poisson_2s_ad(density_fft_a0, density_fft_a6,
                    u_brick_a0, v0_brick_a0, v1_brick_a0, v2_brick_a0, v3_brick_a0, v4_brick_a0, v5_brick_a0,
                    u_brick_a6, v0_brick_a6, v1_brick_a6, v2_brick_a6, v3_brick_a6, v4_brick_a6, v5_brick_a6);
      poisson_2s_ad(density_fft_a1, density_fft_a5,
                    u_brick_a1, v0_brick_a1, v1_brick_a1, v2_brick_a1, v3_brick_a1, v4_brick_a1, v5_brick_a1,
                    u_brick_a5, v0_brick_a5, v1_brick_a5, v2_brick_a5, v3_brick_a5, v4_brick_a5, v5_brick_a5);
      poisson_2s_ad(density_fft_a2, density_fft_a4,
                    u_brick_a2, v0_brick_a2, v1_brick_a2, v2_brick_a2, v3_brick_a2, v4_brick_a2, v5_brick_a2,
                    u_brick_a4, v0_brick_a4, v1_brick_a4, v2_brick_a4, v3_brick_a4, v4_brick_a4, v5_brick_a4);

      cg_6->forward_comm(this, FORWARD_AD_A);

      fieldforce_a_ad();

      if (evflag_atom) cg_peratom_6->forward_comm(this, FORWARD_AD_PERATOM_A);

    }  else {
    
      poisson_ik(work1_6, work2_6, density_fft_a3, fft1_6, fft2_6,
                 nx_pppm_6, ny_pppm_6, nz_pppm_6, nfft_6,
                 nxlo_fft_6, nylo_fft_6, nzlo_fft_6, nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                 nxlo_in_6, nylo_in_6, nzlo_in_6, nxhi_in_6, nyhi_in_6, nzhi_in_6,
                 energy_6, greensfn_6, 
	         fkx_6, fky_6, fkz_6,fkx2_6, fky2_6, fkz2_6,
                 vdx_brick_a3, vdy_brick_a3, vdz_brick_a3, virial_6, vg_6, vg2_6,
                 u_brick_a3, v0_brick_a3, v1_brick_a3, v2_brick_a3, v3_brick_a3, v4_brick_a3, v5_brick_a3);
      poisson_2s_ik(density_fft_a0, density_fft_a6,
                    vdx_brick_a0, vdy_brick_a0, vdz_brick_a0,
                    vdx_brick_a6, vdy_brick_a6, vdz_brick_a6,
                    u_brick_a0, v0_brick_a0, v1_brick_a0, v2_brick_a0, v3_brick_a0, v4_brick_a0, v5_brick_a0,
                    u_brick_a6, v0_brick_a6, v1_brick_a6, v2_brick_a6, v3_brick_a6, v4_brick_a6, v5_brick_a6);
      poisson_2s_ik(density_fft_a1, density_fft_a5,
                    vdx_brick_a1, vdy_brick_a1, vdz_brick_a1,
                    vdx_brick_a5, vdy_brick_a5, vdz_brick_a5,
                    u_brick_a1, v0_brick_a1, v1_brick_a1, v2_brick_a1, v3_brick_a1, v4_brick_a1, v5_brick_a1,
                    u_brick_a5, v0_brick_a5, v1_brick_a5, v2_brick_a5, v3_brick_a5, v4_brick_a5, v5_brick_a5);
      poisson_2s_ik(density_fft_a2, density_fft_a4,
                    vdx_brick_a2, vdy_brick_a2, vdz_brick_a2,
                    vdx_brick_a4, vdy_brick_a4, vdz_brick_a4,
                    u_brick_a2, v0_brick_a2, v1_brick_a2, v2_brick_a2, v3_brick_a2, v4_brick_a2, v5_brick_a2,
                    u_brick_a4, v0_brick_a4, v1_brick_a4, v2_brick_a4, v3_brick_a4, v4_brick_a4, v5_brick_a4);

      cg_6->forward_comm(this, FORWARD_IK_A);

      fieldforce_a_ik();

      if (evflag_atom) cg_peratom_6->forward_comm(this, FORWARD_IK_PERATOM_A);
    }
    if (evflag_atom) fieldforce_a_peratom();
  }

  if (function[3]) {
    //perfrom calculations if no mixing rule applies
    particle_map(delxinv_6, delyinv_6, delzinv_6, shift_6, part2grid_6, nupper_6, nlower_6,
                 nxlo_out_6, nylo_out_6, nzlo_out_6, nxhi_out_6, nyhi_out_6, nzhi_out_6);

    make_rho_none();

    cg_6->reverse_comm(this, REVERSE_RHO_NONE);

    brick2fft_none();

    if (differentiation_flag == 1) {

      int n = 0;
      for (int k = 0; k<nsplit_alloc/2; k++) {
        poisson_none_ad(n,n+1,density_fft_none[n],density_fft_none[n+1],
                        u_brick_none[n],u_brick_none[n+1],
                        v0_brick_none, v1_brick_none, v2_brick_none,
                        v3_brick_none, v4_brick_none, v5_brick_none);
        n += 2;
      }

      cg_6->forward_comm(this,FORWARD_AD_NONE);

      fieldforce_none_ad();

      if (vflag_atom) cg_peratom_6->forward_comm(this,FORWARD_AD_PERATOM_NONE);

    } else {
      int n = 0;
      for (int k = 0; k<nsplit_alloc/2; k++) {

        poisson_none_ik(n,n+1,density_fft_none[n], density_fft_none[n+1],
                        vdx_brick_none[n], vdy_brick_none[n], vdz_brick_none[n],
                        vdx_brick_none[n+1], vdy_brick_none[n+1], vdz_brick_none[n+1],
                        u_brick_none, v0_brick_none, v1_brick_none, v2_brick_none,
                        v3_brick_none, v4_brick_none, v5_brick_none);
        n += 2;
      }

      cg_6->forward_comm(this,FORWARD_IK_NONE);

      fieldforce_none_ik();


      if (evflag_atom) cg_peratom_6->forward_comm(this, FORWARD_IK_PERATOM_NONE);
    }
    if (evflag_atom) fieldforce_none_peratom();
  }

  // update qsum and qsqsum, if needed

  if (eflag_global || eflag_atom) {
    if (qsum_update_flag || (atom->natoms != natoms_original)) {
      qsum_qsq(0);
      natoms_original = atom->natoms;
    }
  }

  // sum energy across procs and add in volume-dependent term

  const double qscale = force->qqrd2e * scale;
  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy_1,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy_1 = energy_all;
    MPI_Allreduce(&energy_6,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy_6 = energy_all;
   
    energy_1 *= 0.5*volume;
    energy_6 *= 0.5*volume;
    
    energy_1 -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy_6 += - MY_PI*MY_PIS/(6*volume)*pow(g_ewald_6,3)*csumij +
      1.0/12.0*pow(g_ewald_6,6)*csum;
    energy_1 *= qscale;
  }

  // sum virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial_1,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*volume*virial_all[i];
    MPI_Allreduce(virial_6,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] += 0.5*volume*virial_all[i];
    if (function[1]+function[2]+function[3]){
      double a =  MY_PI*MY_PIS/(6*volume)*pow(g_ewald_6,3)*csumij;
      virial[0] -= a;
      virial[1] -= a;
      virial[2] -= a;
    }
  }

  if (eflag_atom) {
    if (function[0]) {
      double *q = atom->q;
      for (i = 0; i < atom->nlocal; i++) {
        eatom[i] -= qscale*g_ewald*q[i]*q[i]/MY_PIS + qscale*MY_PI2*q[i]*qsum / (g_ewald*g_ewald*volume); //coulomb self energy correction
      }
    }
    if (function[1] + function[2] + function[3]) {
      int tmp;
      for (i = 0; i < atom->nlocal; i++) {
        tmp = atom->type[i];
        eatom[i] += - MY_PI*MY_PIS/(6*volume)*pow(g_ewald_6,3)*csumi[tmp] +
                      1.0/12.0*pow(g_ewald_6,6)*cii[tmp];
      }
    }
  }
            
  if (vflag_atom) {
    if (function[1] + function[2] + function[3]) {
      int tmp;
      for (i = 0; i < atom->nlocal; i++) {
        tmp = atom->type[i];
        for (int n = 0; n < 3; n++) vatom[i][n] -= MY_PI*MY_PIS/(6*volume)*pow(g_ewald_6,3)*csumi[tmp]; //dispersion self virial correction
      }
    }
  }


  // 2d slab correction

  if (slabflag) slabcorr(eflag);
  if (function[0]) energy += energy_1;
  if (function[1] + function[2] + function[3]) energy += energy_6;

  // convert atoms back from lamda to box coords
  
  if (triclinic) domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   initialize coefficients needed for the dispersion density on the grids
------------------------------------------------------------------------- */

void PPPMDisp::init_coeffs()				// local pair coeffs
{
  int tmp;
  int n = atom->ntypes;
  int converged;
  delete [] B;
  B = NULL;
  if (function[3] + function[2]) {                     // no mixing rule or arithmetic
    if (function[2] && me == 0) {
      if (screen) fprintf(screen,"  Optimizing splitting of Dispersion coefficients\n");
      if (logfile) fprintf(logfile,"  Optimizing splitting of Dispersion coefficients\n");
    }

    // allocate data for eigenvalue decomposition
    double **A=NULL;
    double **Q=NULL;
    if ( n > 1 ) {
      // get dispersion coefficients
      double **b = (double **) force->pair->extract("B",tmp);
      memory->create(A,n,n,"pppm/disp:A");
      memory->create(Q,n,n,"pppm/disp:Q");
      // fill coefficients to matrix a
      for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
          A[i-1][j-1] = b[i][j];
      // transform q to a unity matrix
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          Q[i][j] = 0.0;
      for (int i = 0; i < n; i++)
        Q[i][i] = 1.0;
      // perfrom eigenvalue decomposition with QR algorithm
      converged = qr_alg(A,Q,n);
      if (function[3] && !converged) {
        error->all(FLERR,"Matrix factorization to split dispersion coefficients failed");
      }
      // determine number of used eigenvalues 
      //   based on maximum allowed number or cutoff criterion
      //   sort eigenvalues according to their size with bubble sort
      double t;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n-1-i; j++) {
          if (fabs(A[j][j]) < fabs(A[j+1][j+1])) {
            t = A[j][j];
            A[j][j] = A[j+1][j+1];
            A[j+1][j+1] = t;
            for (int k = 0; k < n; k++) {
              t = Q[k][j];
              Q[k][j] = Q[k][j+1];
              Q[k][j+1] = t;
            }
          }
        }
      }

      //   check which eigenvalue is the first that is smaller
      //   than a specified tolerance
      //   check how many are maximum allowed by the user
      double amax = fabs(A[0][0]);
      double acrit = amax*splittol;
      double bmax = 0;
      double err = 0;
      nsplit = 0;
      for (int i = 0; i < n; i++) {
        if (fabs(A[i][i]) > acrit) nsplit++;
        else {
          bmax = fabs(A[i][i]);
          break;
        }
      }

      err =  bmax/amax;
      if (err > 1.0e-4) {
        char str[128];
        sprintf(str,"Estimated error in splitting of dispersion coeffs is %g",err);
        error->warning(FLERR, str);
      }
      // set B
      B = new double[nsplit*n+nsplit];
      for (int i = 0; i< nsplit; i++) {
        B[i] = A[i][i];
        for (int j = 0; j < n; j++) {
          B[nsplit*(j+1) + i] = Q[j][i];
        }
      }

      nsplit_alloc = nsplit;
      if (nsplit%2 == 1) nsplit_alloc = nsplit + 1;
    } else
        nsplit = 1; // use geometric mixing

    // check if the function should preferably be [1] or [2] or [3]
    if (nsplit == 1) {
      if ( B ) delete [] B;
      function[3] = 0;
      function[2] = 0;
      function[1] = 1;
      if (me == 0) {
        if (screen) fprintf(screen,"  Using geometric mixing for reciprocal space\n");
        if (logfile) fprintf(logfile,"  Using geometric mixing for reciprocal space\n");
      }
    }
    if (function[2] && nsplit <= 6) {
      if (me == 0) {
        if (screen) fprintf(screen,"  Using %d instead of 7 structure factors\n",nsplit);
        if (logfile) fprintf(logfile,"  Using %d instead of 7 structure factors\n",nsplit);
      }
      function[3] = 1;
      function[2] = 0;
    }
    if (function[2] && (nsplit > 6)) {
      if (me == 0) {
        if (screen) fprintf(screen,"  Using 7 structure factors\n");
        if (logfile) fprintf(logfile,"  Using 7 structure factors\n");
      }
      if ( B ) delete [] B;
    }
    if (function[3]) {
      if (me == 0) {
        if (screen) fprintf(screen,"  Using %d structure factors\n",nsplit);
        if (logfile) fprintf(logfile,"  Using %d structure factors\n",nsplit);
      }
      if (nsplit > 9) error->warning(FLERR, "Simulations might be very slow because of large number of structure factors");
    }

    memory->destroy(A);
    memory->destroy(Q);
  }
  if (function[1]) {					// geometric 1/r^6
    double **b = (double **) force->pair->extract("B",tmp);
    B = new double[n+1];
    for (int i=0; i<=n; ++i) B[i] = sqrt(fabs(b[i][i]));
  }
  if (function[2]) {					// arithmetic 1/r^6
    //cannot use epsilon, because this has not been set yet
    double **epsilon = (double **) force->pair->extract("epsilon",tmp);  
    //cannot use sigma, because this has not been set yet
    double **sigma = (double **) force->pair->extract("sigma",tmp);  
    if (!(epsilon&&sigma))
      error->all(FLERR,"Epsilon or sigma reference not set by pair style in PPPMDisp");
    double eps_i, sigma_i, sigma_n, *bi = B = new double[7*n+7];
    double c[7] = {
      1.0, sqrt(6.0), sqrt(15.0), sqrt(20.0), sqrt(15.0), sqrt(6.0), 1.0};
    for (int i=0; i<=n; ++i) {
      eps_i = sqrt(epsilon[i][i]);
      sigma_i = sigma[i][i];
      sigma_n = 1.0;
      for (int j=0; j<7; ++j) {
        *(bi++) = sigma_n*eps_i*c[j]*0.25;
        sigma_n *= sigma_i;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Eigenvalue decomposition of a real, symmetric matrix with the QR
   method (includes transpformation to Tridiagonal Matrix + Wilkinson
   shift)
------------------------------------------------------------------------- */

int PPPMDisp::qr_alg(double **A, double **Q, int n)
{
  int converged = 0;
  double an1, an, bn1, d, mue;
  // allocate some memory for the required operations
  double **A0,**Qi,**C,**D,**E;
  // make a copy of A for convergence check
  memory->create(A0,n,n,"pppm/disp:A0");
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      A0[i][j] = A[i][j];

  // allocate an auxiliary matrix Qi
  memory->create(Qi,n,n,"pppm/disp:Qi");

  // alllocate an auxillary matrices for the matrix multiplication
  memory->create(C,n,n,"pppm/disp:C");
  memory->create(D,n,n,"pppm/disp:D");
  memory->create(E,n,n,"pppm/disp:E");

  // transform Matrix A to Tridiagonal form
  hessenberg(A,Q,n);

  // start loop for the matrix factorization
  int count = 0;
  int countmax = 100000;
  while (1) {
    // make a Wilkinson shift
    an1 = A[n-2][n-2];
    an = A[n-1][n-1];
    bn1 = A[n-2][n-1];
    d = (an1-an)/2;
    mue = an + d - copysign(1.,d)*sqrt(d*d + bn1*bn1);
    for (int i = 0; i < n; i++) 
      A[i][i] -= mue;

    // perform a QR factorization for a tridiagonal matrix A
    qr_tri(Qi,A,n);

    // update the matrices
    mmult(A,Qi,C,n);
    mmult(Q,Qi,C,n);

    // backward Wilkinson shift
    for (int i = 0; i < n; i++)
      A[i][i] += mue;

    // check the convergence
    converged = check_convergence(A,Q,A0,C,D,E,n);
    if (converged) break;
    count = count + 1;
    if (count == countmax) break;
  }
  
  // free allocated memory
  memory->destroy(Qi);
  memory->destroy(A0);
  memory->destroy(C);
  memory->destroy(D);
  memory->destroy(E);
  
  return converged;
}

/* ----------------------------------------------------------------------
   Transform a Matrix to Hessenberg form (for symmetric Matrices, the 
   result will be a tridiagonal matrix)
------------------------------------------------------------------------- */

void PPPMDisp::hessenberg(double **A, double **Q, int n)
{
  double r,a,b,c,s,x1,x2;
  for (int i = 0; i < n-1; i++) {
    for (int j = i+2; j < n; j++) {
      // compute coeffs for the rotation matrix
      a = A[i+1][i];
      b = A[j][i];
      r = sqrt(a*a + b*b);
      c = a/r;
      s = b/r;
      // update the entries of A with multiplication from the left
      for (int k = 0; k < n; k++) {
        x1 = A[i+1][k];
        x2 = A[j][k];
        A[i+1][k] = c*x1 + s*x2;
        A[j][k] = -s*x1 + c*x2;
      }
      // update the entries of A and Q with a multiplication from the right
      for (int k = 0; k < n; k++) {
        x1 = A[k][i+1];
        x2 = A[k][j];
        A[k][i+1] = c*x1 + s*x2;
        A[k][j] = -s*x1 + c*x2;
        x1 = Q[k][i+1];
        x2 = Q[k][j];
        Q[k][i+1] = c*x1 + s*x2;
        Q[k][j] = -s*x1 + c*x2;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   QR factorization for a tridiagonal matrix; Result of the factorization
   is stored in A and Qi
------------------------------------------------------------------------- */

void PPPMDisp::qr_tri(double** Qi,double** A,int n)
{
  double r,a,b,c,s,x1,x2;
  int j,k,k0,kmax;
  // make Qi a unity matrix
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      Qi[i][j] = 0.0;
  for (int i = 0; i < n; i++)
    Qi[i][i] = 1.0;
  // loop over main diagonal and first of diagonal of A
  for (int i = 0; i < n-1; i++) {
    j = i+1;
    // coefficients of the rotation matrix
    a = A[i][i];
    b = A[j][i];
    r = sqrt(a*a + b*b);
    c = a/r;
    s = b/r;
    // update the entries of A and Q
    k0 = (i-1>0)?i-1:0;   //min(i-1,0);
    kmax = (i+3<n)?i+3:n;  //min(i+3,n);
    for (k = k0; k < kmax; k++) {
      x1 = A[i][k];
      x2 = A[j][k];
      A[i][k] = c*x1 + s*x2;
      A[j][k] = -s*x1 + c*x2;
    }
    for (k = 0; k < n; k++) {
      x1 = Qi[k][i];
      x2 = Qi[k][j];
      Qi[k][i] = c*x1 + s*x2;
      Qi[k][j] = -s*x1 + c*x2;
    }
  }
}

/* ----------------------------------------------------------------------
   Multiply two matrices A and B, store the result in A; C provides
   some memory to store intermediate results
------------------------------------------------------------------------- */

void PPPMDisp::mmult(double** A, double** B, double** C, int n)
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      C[i][j] = 0.0;

  // perform matrix multiplication 
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
        C[i][j] += A[i][k] * B[k][j];
  // copy the result back to matrix A
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      A[i][j] = C[i][j];
}

/* ----------------------------------------------------------------------
   Check if the factorization has converged by comparing all elements of the
   original matrix and the new matrix
------------------------------------------------------------------------- */

int PPPMDisp::check_convergence(double** A,double** Q,double** A0,
                                double** C,double** D,double** E,int n)
{
  double eps = 1.0e-8;
  int converged = 1;
  double epsmax = -1;
  double Bmax = 0.0;
  double diff;
  // get the largest eigenvalue of the original matrix
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      Bmax = (Bmax>A0[i][j])?Bmax:A0[i][j];  //max(Bmax,A0[i][j]);
  double epsabs = eps*Bmax;
  
  // reconstruct the original matrix
  // store the diagonal elements in D
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      D[i][j] = 0.0;
  for (int i = 0; i < n; i++)
    D[i][i] = A[i][i];
  // store matrix Q in E
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      E[i][j] = Q[i][j];
  // E = Q*A
  mmult(E,D,C,n);
  // store transpose of Q in D
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      D[i][j] = Q[j][i];
  // E = Q*A*Q.t
  mmult(E,D,C,n);

  //compare the original matrix and the final matrix
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      diff = A0[i][j] - E[i][j];
      epsmax = (epsmax>fabs(diff))?epsmax:fabs(diff);//max(epsmax,fabs(diff));
    }
  }
  if (epsmax > epsabs) converged = 0;
  return converged;
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order 
------------------------------------------------------------------------- */

void PPPMDisp::allocate()
{

  int (*procneigh)[2] = comm->procneigh;

  if (function[0]) {
    memory->create(work1,2*nfft_both,"pppm/disp:work1");
    memory->create(work2,2*nfft_both,"pppm/disp:work2");

    memory->create1d_offset(fkx,nxlo_fft,nxhi_fft,"pppm/disp:fkx");
    memory->create1d_offset(fky,nylo_fft,nyhi_fft,"pppm/disp:fky");
    memory->create1d_offset(fkz,nzlo_fft,nzhi_fft,"pppm/disp:fkz");

    memory->create1d_offset(fkx2,nxlo_fft,nxhi_fft,"pppm/disp:fkx2");
    memory->create1d_offset(fky2,nylo_fft,nyhi_fft,"pppm/disp:fky2");
    memory->create1d_offset(fkz2,nzlo_fft,nzhi_fft,"pppm/disp:fkz2");


    memory->create(gf_b,order,"pppm/disp:gf_b");
    memory->create2d_offset(rho1d,3,-order/2,order/2,"pppm/disp:rho1d");
    memory->create2d_offset(rho_coeff,order,(1-order)/2,order/2,"pppm/disp:rho_coeff");
    memory->create2d_offset(drho1d,3,-order/2,order/2,"pppm/disp:rho1d");
    memory->create2d_offset(drho_coeff,order,(1-order)/2,order/2,"pppm/disp:drho_coeff");

    memory->create(greensfn,nfft_both,"pppm/disp:greensfn");
    memory->create(vg,nfft_both,6,"pppm/disp:vg");
    memory->create(vg2,nfft_both,3,"pppm/disp:vg2");

    memory->create3d_offset(density_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
  			    nxlo_out,nxhi_out,"pppm/disp:density_brick");
    if ( differentiation_flag == 1) {
      memory->create3d_offset(u_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
  		  	      nxlo_out,nxhi_out,"pppm/disp:u_brick");
      memory->create(sf_precoeff1,nfft_both,"pppm/disp:sf_precoeff1");
      memory->create(sf_precoeff2,nfft_both,"pppm/disp:sf_precoeff2");
      memory->create(sf_precoeff3,nfft_both,"pppm/disp:sf_precoeff3");
      memory->create(sf_precoeff4,nfft_both,"pppm/disp:sf_precoeff4");
      memory->create(sf_precoeff5,nfft_both,"pppm/disp:sf_precoeff5");
      memory->create(sf_precoeff6,nfft_both,"pppm/disp:sf_precoeff6");

    } else {
      memory->create3d_offset(vdx_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
  			      nxlo_out,nxhi_out,"pppm/disp:vdx_brick");
      memory->create3d_offset(vdy_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
			      nxlo_out,nxhi_out,"pppm/disp:vdy_brick");
      memory->create3d_offset(vdz_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
			      nxlo_out,nxhi_out,"pppm/disp:vdz_brick");
    }
    memory->create(density_fft,nfft_both,"pppm/disp:density_fft");

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

  if (differentiation_flag == 1)
    cg = new GridComm(lmp,world,1,1,
                      nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                      nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                      procneigh[0][0],procneigh[0][1],procneigh[1][0],
                      procneigh[1][1],procneigh[2][0],procneigh[2][1]);
  else
    cg = new GridComm(lmp,world,3,1,
                      nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                      nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                      procneigh[0][0],procneigh[0][1],procneigh[1][0],
                      procneigh[1][1],procneigh[2][0],procneigh[2][1]);
  }

  if (function[1]) {
    memory->create(work1_6,2*nfft_both_6,"pppm/disp:work1_6");
    memory->create(work2_6,2*nfft_both_6,"pppm/disp:work2_6");

    memory->create1d_offset(fkx_6,nxlo_fft_6,nxhi_fft_6,"pppm/disp:fkx_6");
    memory->create1d_offset(fky_6,nylo_fft_6,nyhi_fft_6,"pppm/disp:fky_6");
    memory->create1d_offset(fkz_6,nzlo_fft_6,nzhi_fft_6,"pppm/disp:fkz_6");

    memory->create1d_offset(fkx2_6,nxlo_fft_6,nxhi_fft_6,"pppm/disp:fkx2_6");
    memory->create1d_offset(fky2_6,nylo_fft_6,nyhi_fft_6,"pppm/disp:fky2_6");
    memory->create1d_offset(fkz2_6,nzlo_fft_6,nzhi_fft_6,"pppm/disp:fkz2_6");

    memory->create(gf_b_6,order_6,"pppm/disp:gf_b_6");
    memory->create2d_offset(rho1d_6,3,-order_6/2,order_6/2,"pppm/disp:rho1d_6");
    memory->create2d_offset(rho_coeff_6,order_6,(1-order_6)/2,order_6/2,"pppm/disp:rho_coeff_6");
    memory->create2d_offset(drho1d_6,3,-order_6/2,order_6/2,"pppm/disp:drho1d_6");
    memory->create2d_offset(drho_coeff_6,order_6,(1-order_6)/2,order_6/2,"pppm/disp:drho_coeff_6");

    memory->create(greensfn_6,nfft_both_6,"pppm/disp:greensfn_6");
    memory->create(vg_6,nfft_both_6,6,"pppm/disp:vg_6");
    memory->create(vg2_6,nfft_both_6,3,"pppm/disp:vg2_6");

    memory->create3d_offset(density_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  			    nxlo_out_6,nxhi_out_6,"pppm/disp:density_brick_g");
    if ( differentiation_flag == 1) {
      memory->create3d_offset(u_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_g");

      memory->create(sf_precoeff1_6,nfft_both_6,"pppm/disp:sf_precoeff1_6");
      memory->create(sf_precoeff2_6,nfft_both_6,"pppm/disp:sf_precoeff2_6");
      memory->create(sf_precoeff3_6,nfft_both_6,"pppm/disp:sf_precoeff3_6");
      memory->create(sf_precoeff4_6,nfft_both_6,"pppm/disp:sf_precoeff4_6");
      memory->create(sf_precoeff5_6,nfft_both_6,"pppm/disp:sf_precoeff5_6");
      memory->create(sf_precoeff6_6,nfft_both_6,"pppm/disp:sf_precoeff6_6");

    }  else {
      memory->create3d_offset(vdx_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdx_brick_g");
      memory->create3d_offset(vdy_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdy_brick_g");
      memory->create3d_offset(vdz_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdz_brick_g");
    }
    memory->create(density_fft_g,nfft_both_6,"pppm/disp:density_fft_g");


    int tmp;

    fft1_6 = new FFT3d(lmp,world,nx_pppm_6,ny_pppm_6,nz_pppm_6,
		     nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		     nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		     0,0,&tmp,collective_flag);

    fft2_6 = new FFT3d(lmp,world,nx_pppm_6,ny_pppm_6,nz_pppm_6,
		     nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		     nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
		     0,0,&tmp,collective_flag);

    remap_6 = new Remap(lmp,world,
		      nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
		      nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		      1,0,0,FFT_PRECISION,collective_flag);

    // create ghost grid object for rho and electric field communication

    if (differentiation_flag == 1)
      cg_6 = new GridComm(lmp,world,1,1,
                        nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                        nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                        procneigh[0][0],procneigh[0][1],procneigh[1][0],
                        procneigh[1][1],procneigh[2][0],procneigh[2][1]);
    else
      cg_6 = new GridComm(lmp,world,3,1,
                        nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                        nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                        procneigh[0][0],procneigh[0][1],procneigh[1][0],
                        procneigh[1][1],procneigh[2][0],procneigh[2][1]);
  }

  if (function[2]) {
    memory->create(work1_6,2*nfft_both_6,"pppm/disp:work1_6");
    memory->create(work2_6,2*nfft_both_6,"pppm/disp:work2_6");

    memory->create1d_offset(fkx_6,nxlo_fft_6,nxhi_fft_6,"pppm/disp:fkx_6");
    memory->create1d_offset(fky_6,nylo_fft_6,nyhi_fft_6,"pppm/disp:fky_6");
    memory->create1d_offset(fkz_6,nzlo_fft_6,nzhi_fft_6,"pppm/disp:fkz_6");

    memory->create1d_offset(fkx2_6,nxlo_fft_6,nxhi_fft_6,"pppm/disp:fkx2_6");
    memory->create1d_offset(fky2_6,nylo_fft_6,nyhi_fft_6,"pppm/disp:fky2_6");
    memory->create1d_offset(fkz2_6,nzlo_fft_6,nzhi_fft_6,"pppm/disp:fkz2_6");

    memory->create(gf_b_6,order_6,"pppm/disp:gf_b_6");
    memory->create2d_offset(rho1d_6,3,-order_6/2,order_6/2,"pppm/disp:rho1d_6");
    memory->create2d_offset(rho_coeff_6,order_6,(1-order_6)/2,order_6/2,"pppm/disp:rho_coeff_6");
    memory->create2d_offset(drho1d_6,3,-order_6/2,order_6/2,"pppm/disp:drho1d_6");
    memory->create2d_offset(drho_coeff_6,order_6,(1-order_6)/2,order_6/2,"pppm/disp:drho_coeff_6");

    memory->create(greensfn_6,nfft_both_6,"pppm/disp:greensfn_6");
    memory->create(vg_6,nfft_both_6,6,"pppm/disp:vg_6");
    memory->create(vg2_6,nfft_both_6,3,"pppm/disp:vg2_6");

    memory->create3d_offset(density_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  			    nxlo_out_6,nxhi_out_6,"pppm/disp:density_brick_a0");
    memory->create3d_offset(density_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  			    nxlo_out_6,nxhi_out_6,"pppm/disp:density_brick_a1");
    memory->create3d_offset(density_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  			    nxlo_out_6,nxhi_out_6,"pppm/disp:density_brick_a2");
    memory->create3d_offset(density_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  			    nxlo_out_6,nxhi_out_6,"pppm/disp:density_brick_a3");
    memory->create3d_offset(density_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  			    nxlo_out_6,nxhi_out_6,"pppm/disp:density_brick_a4");
    memory->create3d_offset(density_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  			    nxlo_out_6,nxhi_out_6,"pppm/disp:density_brick_a5");
    memory->create3d_offset(density_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  			    nxlo_out_6,nxhi_out_6,"pppm/disp:density_brick_a6");

    memory->create(density_fft_a0,nfft_both_6,"pppm/disp:density_fft_a0");
    memory->create(density_fft_a1,nfft_both_6,"pppm/disp:density_fft_a1");
    memory->create(density_fft_a2,nfft_both_6,"pppm/disp:density_fft_a2");
    memory->create(density_fft_a3,nfft_both_6,"pppm/disp:density_fft_a3");
    memory->create(density_fft_a4,nfft_both_6,"pppm/disp:density_fft_a4");
    memory->create(density_fft_a5,nfft_both_6,"pppm/disp:density_fft_a5");
    memory->create(density_fft_a6,nfft_both_6,"pppm/disp:density_fft_a6");


    if ( differentiation_flag == 1 ) {
      memory->create3d_offset(u_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a0");
      memory->create3d_offset(u_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a1");
      memory->create3d_offset(u_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a2");
      memory->create3d_offset(u_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a3");
      memory->create3d_offset(u_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a4");
      memory->create3d_offset(u_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a5");
      memory->create3d_offset(u_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a6");

      memory->create(sf_precoeff1_6,nfft_both_6,"pppm/disp:sf_precoeff1_6");
      memory->create(sf_precoeff2_6,nfft_both_6,"pppm/disp:sf_precoeff2_6");
      memory->create(sf_precoeff3_6,nfft_both_6,"pppm/disp:sf_precoeff3_6");
      memory->create(sf_precoeff4_6,nfft_both_6,"pppm/disp:sf_precoeff4_6");
      memory->create(sf_precoeff5_6,nfft_both_6,"pppm/disp:sf_precoeff5_6");
      memory->create(sf_precoeff6_6,nfft_both_6,"pppm/disp:sf_precoeff6_6");

    } else {

      memory->create3d_offset(vdx_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdx_brick_a0");
      memory->create3d_offset(vdy_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdy_brick_a0");
      memory->create3d_offset(vdz_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdz_brick_a0");

      memory->create3d_offset(vdx_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdx_brick_a1");
      memory->create3d_offset(vdy_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdy_brick_a1");
      memory->create3d_offset(vdz_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdz_brick_a1");

      memory->create3d_offset(vdx_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdx_brick_a2");
      memory->create3d_offset(vdy_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdy_brick_a2");
      memory->create3d_offset(vdz_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdz_brick_a2");

      memory->create3d_offset(vdx_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdx_brick_a3");
      memory->create3d_offset(vdy_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdy_brick_a3");
      memory->create3d_offset(vdz_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdz_brick_a3");

      memory->create3d_offset(vdx_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdx_brick_a4");
      memory->create3d_offset(vdy_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdy_brick_a4");
      memory->create3d_offset(vdz_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdz_brick_a4");

      memory->create3d_offset(vdx_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdx_brick_a5");
      memory->create3d_offset(vdy_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdy_brick_a5");
      memory->create3d_offset(vdz_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdz_brick_a5");

      memory->create3d_offset(vdx_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdx_brick_a6");
      memory->create3d_offset(vdy_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdy_brick_a6");
      memory->create3d_offset(vdz_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdz_brick_a6");
    }



    int tmp;

    fft1_6 = new FFT3d(lmp,world,nx_pppm_6,ny_pppm_6,nz_pppm_6,
		     nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		     nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		     0,0,&tmp,collective_flag);

    fft2_6 = new FFT3d(lmp,world,nx_pppm_6,ny_pppm_6,nz_pppm_6,
		     nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		     nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
		     0,0,&tmp,collective_flag);

    remap_6 = new Remap(lmp,world,
		      nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
		      nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		      1,0,0,FFT_PRECISION,collective_flag);

    // create ghost grid object for rho and electric field communication


    if (differentiation_flag == 1)
      cg_6 = new GridComm(lmp,world,7,7,
                        nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                        nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                        procneigh[0][0],procneigh[0][1],procneigh[1][0],
                        procneigh[1][1],procneigh[2][0],procneigh[2][1]);
    else
      cg_6 = new GridComm(lmp,world,21,7,
                        nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                        nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                        procneigh[0][0],procneigh[0][1],procneigh[1][0],
                        procneigh[1][1],procneigh[2][0],procneigh[2][1]);
  }  

  if (function[3]) {
    memory->create(work1_6,2*nfft_both_6,"pppm/disp:work1_6");
    memory->create(work2_6,2*nfft_both_6,"pppm/disp:work2_6");

    memory->create1d_offset(fkx_6,nxlo_fft_6,nxhi_fft_6,"pppm/disp:fkx_6");
    memory->create1d_offset(fky_6,nylo_fft_6,nyhi_fft_6,"pppm/disp:fky_6");
    memory->create1d_offset(fkz_6,nzlo_fft_6,nzhi_fft_6,"pppm/disp:fkz_6");

    memory->create1d_offset(fkx2_6,nxlo_fft_6,nxhi_fft_6,"pppm/disp:fkx2_6");
    memory->create1d_offset(fky2_6,nylo_fft_6,nyhi_fft_6,"pppm/disp:fky2_6");
    memory->create1d_offset(fkz2_6,nzlo_fft_6,nzhi_fft_6,"pppm/disp:fkz2_6");

    memory->create(gf_b_6,order_6,"pppm/disp:gf_b_6");
    memory->create2d_offset(rho1d_6,3,-order_6/2,order_6/2,"pppm/disp:rho1d_6");
    memory->create2d_offset(rho_coeff_6,order_6,(1-order_6)/2,order_6/2,"pppm/disp:rho_coeff_6");
    memory->create2d_offset(drho1d_6,3,-order_6/2,order_6/2,"pppm/disp:drho1d_6");
    memory->create2d_offset(drho_coeff_6,order_6,(1-order_6)/2,order_6/2,"pppm/disp:drho_coeff_6");

    memory->create(greensfn_6,nfft_both_6,"pppm/disp:greensfn_6");
    memory->create(vg_6,nfft_both_6,6,"pppm/disp:vg_6");
    memory->create(vg2_6,nfft_both_6,3,"pppm/disp:vg2_6");

    memory->create4d_offset(density_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  			    nxlo_out_6,nxhi_out_6,"pppm/disp:density_brick_none");
    if ( differentiation_flag == 1) {
      memory->create4d_offset(u_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_none");

      memory->create(sf_precoeff1_6,nfft_both_6,"pppm/disp:sf_precoeff1_6");
      memory->create(sf_precoeff2_6,nfft_both_6,"pppm/disp:sf_precoeff2_6");
      memory->create(sf_precoeff3_6,nfft_both_6,"pppm/disp:sf_precoeff3_6");
      memory->create(sf_precoeff4_6,nfft_both_6,"pppm/disp:sf_precoeff4_6");
      memory->create(sf_precoeff5_6,nfft_both_6,"pppm/disp:sf_precoeff5_6");
      memory->create(sf_precoeff6_6,nfft_both_6,"pppm/disp:sf_precoeff6_6");

    }  else {
      memory->create4d_offset(vdx_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdx_brick_none");
      memory->create4d_offset(vdy_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdy_brick_none");
      memory->create4d_offset(vdz_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
			      nxlo_out_6,nxhi_out_6,"pppm/disp:vdz_brick_none");
    }
    memory->create(density_fft_none,nsplit_alloc,nfft_both_6,"pppm/disp:density_fft_none");


    int tmp;

    fft1_6 = new FFT3d(lmp,world,nx_pppm_6,ny_pppm_6,nz_pppm_6,
		     nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		     nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		     0,0,&tmp,collective_flag);

    fft2_6 = new FFT3d(lmp,world,nx_pppm_6,ny_pppm_6,nz_pppm_6,
		     nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		     nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
		     0,0,&tmp,collective_flag);

    remap_6 = new Remap(lmp,world,
		      nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
		      nxlo_fft_6,nxhi_fft_6,nylo_fft_6,nyhi_fft_6,nzlo_fft_6,nzhi_fft_6,
		      1,0,0,FFT_PRECISION,collective_flag);

    // create ghost grid object for rho and electric field communication

    if (differentiation_flag == 1)
      cg_6 = new GridComm(lmp,world,nsplit_alloc,nsplit_alloc,
                        nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                        nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                        procneigh[0][0],procneigh[0][1],procneigh[1][0],
                        procneigh[1][1],procneigh[2][0],procneigh[2][1]);
    else
      cg_6 = new GridComm(lmp,world,3*nsplit_alloc,nsplit_alloc,
                        nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                        nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                        procneigh[0][0],procneigh[0][1],procneigh[1][0],
                        procneigh[1][1],procneigh[2][0],procneigh[2][1]);
  }

}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
   for per atom calculations 
------------------------------------------------------------------------- */

void PPPMDisp::allocate_peratom()
{

  int (*procneigh)[2] = comm->procneigh;

  if (function[0]) {

    if (differentiation_flag != 1)
      memory->create3d_offset(u_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
    	                      nxlo_out,nxhi_out,"pppm/disp:u_brick");

    memory->create3d_offset(v0_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
			    nxlo_out,nxhi_out,"pppm/disp:v0_brick");
    memory->create3d_offset(v1_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
  			    nxlo_out,nxhi_out,"pppm/disp:v1_brick");
    memory->create3d_offset(v2_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
  			    nxlo_out,nxhi_out,"pppm/disp:v2_brick");
    memory->create3d_offset(v3_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
  			    nxlo_out,nxhi_out,"pppm/disp:v3_brick");
    memory->create3d_offset(v4_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
  			    nxlo_out,nxhi_out,"pppm/disp:v4_brick");
    memory->create3d_offset(v5_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
  			    nxlo_out,nxhi_out,"pppm/disp:v5_brick");

    // create ghost grid object for rho and electric field communication

    if (differentiation_flag == 1)
      cg_peratom =
        new GridComm(lmp,world,6,1,
                     nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                     nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                     procneigh[0][0],procneigh[0][1],procneigh[1][0],
                     procneigh[1][1],procneigh[2][0],procneigh[2][1]);
    else
      cg_peratom =
        new GridComm(lmp,world,7,1,
                     nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                     nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                     procneigh[0][0],procneigh[0][1],procneigh[1][0],
                     procneigh[1][1],procneigh[2][0],procneigh[2][1]);

  }


  if (function[1]) {

    if ( differentiation_flag != 1 )
      memory->create3d_offset(u_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_g");

    memory->create3d_offset(v0_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v0_brick_g");
    memory->create3d_offset(v1_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v1_brick_g");
    memory->create3d_offset(v2_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v2_brick_g");
    memory->create3d_offset(v3_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v3_brick_g");
    memory->create3d_offset(v4_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v4_brick_g");
    memory->create3d_offset(v5_brick_g,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v5_brick_g");

    // create ghost grid object for rho and electric field communication

    if (differentiation_flag == 1)
      cg_peratom_6 =
        new GridComm(lmp,world,6,1,
                     nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                     nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                     procneigh[0][0],procneigh[0][1],procneigh[1][0],
                     procneigh[1][1],procneigh[2][0],procneigh[2][1]);
    else
      cg_peratom_6 =
        new GridComm(lmp,world,7,1,
                     nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                     nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                     procneigh[0][0],procneigh[0][1],procneigh[1][0],
                     procneigh[1][1],procneigh[2][0],procneigh[2][1]);

  }

  if (function[2]) {
   
    if ( differentiation_flag != 1 ) {
      memory->create3d_offset(u_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a0");
      memory->create3d_offset(u_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a1");
      memory->create3d_offset(u_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a2");
      memory->create3d_offset(u_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a3");
      memory->create3d_offset(u_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a4");
      memory->create3d_offset(u_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a5");
      memory->create3d_offset(u_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_a6");
    }

    memory->create3d_offset(v0_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v0_brick_a0");
    memory->create3d_offset(v1_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
    	                        nxlo_out_6,nxhi_out_6,"pppm/disp:v1_brick_a0");
    memory->create3d_offset(v2_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v2_brick_a0");
    memory->create3d_offset(v3_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v3_brick_a0");
    memory->create3d_offset(v4_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v4_brick_a0");
    memory->create3d_offset(v5_brick_a0,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v5_brick_a0");

    memory->create3d_offset(v0_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v0_brick_a1");
    memory->create3d_offset(v1_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
   	                        nxlo_out_6,nxhi_out_6,"pppm/disp:v1_brick_a1");
    memory->create3d_offset(v2_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v2_brick_a1");
    memory->create3d_offset(v3_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v3_brick_a1");
    memory->create3d_offset(v4_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  	  	                nxlo_out_6,nxhi_out_6,"pppm/disp:v4_brick_a1");
    memory->create3d_offset(v5_brick_a1,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v5_brick_a1");

    memory->create3d_offset(v0_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v0_brick_a2");
    memory->create3d_offset(v1_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v1_brick_a2");
    memory->create3d_offset(v2_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v2_brick_a2");
    memory->create3d_offset(v3_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v3_brick_a2");
    memory->create3d_offset(v4_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v4_brick_a2");
    memory->create3d_offset(v5_brick_a2,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v5_brick_a2");

    memory->create3d_offset(v0_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v0_brick_a3");
    memory->create3d_offset(v1_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v1_brick_a3");
    memory->create3d_offset(v2_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v2_brick_a3");
    memory->create3d_offset(v3_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  	  	                nxlo_out_6,nxhi_out_6,"pppm/disp:v3_brick_a3");
    memory->create3d_offset(v4_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v4_brick_a3");
    memory->create3d_offset(v5_brick_a3,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v5_brick_a3");

    memory->create3d_offset(v0_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v0_brick_a4");
    memory->create3d_offset(v1_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v1_brick_a4");
    memory->create3d_offset(v2_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v2_brick_a4");
    memory->create3d_offset(v3_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v3_brick_a4");
    memory->create3d_offset(v4_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v4_brick_a4");
    memory->create3d_offset(v5_brick_a4,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v5_brick_a4");

    memory->create3d_offset(v0_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v0_brick_a5");
    memory->create3d_offset(v1_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v1_brick_a5");
    memory->create3d_offset(v2_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v2_brick_a5");
    memory->create3d_offset(v3_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v3_brick_a5");
    memory->create3d_offset(v4_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v4_brick_a5");
    memory->create3d_offset(v5_brick_a5,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v5_brick_a5");

    memory->create3d_offset(v0_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  	  	                nxlo_out_6,nxhi_out_6,"pppm/disp:v0_brick_a6");
    memory->create3d_offset(v1_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v1_brick_a6");
    memory->create3d_offset(v2_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v2_brick_a6");
    memory->create3d_offset(v3_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v3_brick_a6");
    memory->create3d_offset(v4_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v4_brick_a6");
    memory->create3d_offset(v5_brick_a6,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	        nxlo_out_6,nxhi_out_6,"pppm/disp:v5_brick_a6");

    // create ghost grid object for rho and electric field communication

    if (differentiation_flag == 1)
      cg_peratom_6 =
        new GridComm(lmp,world,42,1,
                     nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                     nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                     procneigh[0][0],procneigh[0][1],procneigh[1][0],
                     procneigh[1][1],procneigh[2][0],procneigh[2][1]);
    else
      cg_peratom_6 =
        new GridComm(lmp,world,49,1,
                     nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                     nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                     procneigh[0][0],procneigh[0][1],procneigh[1][0],
                     procneigh[1][1],procneigh[2][0],procneigh[2][1]);

  }  

  if (function[3]) {

    if ( differentiation_flag != 1 )
      memory->create4d_offset(u_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	      nxlo_out_6,nxhi_out_6,"pppm/disp:u_brick_none");

    memory->create4d_offset(v0_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v0_brick_none");
    memory->create4d_offset(v1_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v1_brick_none");
    memory->create4d_offset(v2_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v2_brick_none");
    memory->create4d_offset(v3_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v3_brick_none");
    memory->create4d_offset(v4_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v4_brick_none");
    memory->create4d_offset(v5_brick_none,nsplit_alloc,nzlo_out_6,nzhi_out_6,nylo_out_6,nyhi_out_6,
  		  	    nxlo_out_6,nxhi_out_6,"pppm/disp:v5_brick_none");

    // create ghost grid object for rho and electric field communication

    if (differentiation_flag == 1)
      cg_peratom_6 =
        new GridComm(lmp,world,6*nsplit_alloc,1,
                     nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                     nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                     procneigh[0][0],procneigh[0][1],procneigh[1][0],
                     procneigh[1][1],procneigh[2][0],procneigh[2][1]);
    else
      cg_peratom_6 =
        new GridComm(lmp,world,7*nsplit_alloc,1,
                     nxlo_in_6,nxhi_in_6,nylo_in_6,nyhi_in_6,nzlo_in_6,nzhi_in_6,
                     nxlo_out_6,nxhi_out_6,nylo_out_6,nyhi_out_6,nzlo_out_6,nzhi_out_6,
                     procneigh[0][0],procneigh[0][1],procneigh[1][0],
                     procneigh[1][1],procneigh[2][0],procneigh[2][1]);

  }
}


/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order 
------------------------------------------------------------------------- */

void PPPMDisp::deallocate()
{
  memory->destroy3d_offset(density_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(vdx_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(vdy_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(vdz_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy(density_fft);
  density_brick = vdx_brick = vdy_brick = vdz_brick = NULL;
  density_fft = NULL;

  memory->destroy3d_offset(density_brick_g,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdx_brick_g,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdy_brick_g,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdz_brick_g,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy(density_fft_g);
  density_brick_g = vdx_brick_g = vdy_brick_g = vdz_brick_g = NULL;
  density_fft_g = NULL;

  memory->destroy3d_offset(density_brick_a0,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdx_brick_a0,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdy_brick_a0,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdz_brick_a0,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy(density_fft_a0);
  density_brick_a0 = vdx_brick_a0 = vdy_brick_a0 = vdz_brick_a0 = NULL;
  density_fft_a0 = NULL;

  memory->destroy3d_offset(density_brick_a1,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdx_brick_a1,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdy_brick_a1,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdz_brick_a1,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy(density_fft_a1);
  density_brick_a1 = vdx_brick_a1 = vdy_brick_a1 = vdz_brick_a1 = NULL;
  density_fft_a1 = NULL;

  memory->destroy3d_offset(density_brick_a2,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdx_brick_a2,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdy_brick_a2,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdz_brick_a2,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy(density_fft_a2);
  density_brick_a2 = vdx_brick_a2 = vdy_brick_a2 = vdz_brick_a2 = NULL;
  density_fft_a2 = NULL;

  memory->destroy3d_offset(density_brick_a3,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdx_brick_a3,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdy_brick_a3,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdz_brick_a3,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy(density_fft_a3);
  density_brick_a3 = vdx_brick_a3 = vdy_brick_a3 = vdz_brick_a3 = NULL;
  density_fft_a3 = NULL;
 
  memory->destroy3d_offset(density_brick_a4,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdx_brick_a4,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdy_brick_a4,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdz_brick_a4,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy(density_fft_a4);
  density_brick_a4 = vdx_brick_a4 = vdy_brick_a4 = vdz_brick_a4 = NULL;
  density_fft_a4 = NULL;

  memory->destroy3d_offset(density_brick_a5,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdx_brick_a5,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdy_brick_a5,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdz_brick_a5,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy(density_fft_a5);
  density_brick_a5 = vdx_brick_a5 = vdy_brick_a5 = vdz_brick_a5 = NULL;
  density_fft_a5 = NULL;

  memory->destroy3d_offset(density_brick_a6,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdx_brick_a6,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdy_brick_a6,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy3d_offset(vdz_brick_a6,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy(density_fft_a6);
  density_brick_a6 = vdx_brick_a6 = vdy_brick_a6 = vdz_brick_a6 = NULL;
  density_fft_a6 = NULL;

  memory->destroy4d_offset(density_brick_none,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy4d_offset(vdx_brick_none,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy4d_offset(vdy_brick_none,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy4d_offset(vdz_brick_none,nzlo_out_6,nylo_out_6,nxlo_out_6);
  memory->destroy(density_fft_none);
  density_brick_none = vdx_brick_none = vdy_brick_none = vdz_brick_none = NULL;
  density_fft_none = NULL;

  memory->destroy(sf_precoeff1);
  memory->destroy(sf_precoeff2);
  memory->destroy(sf_precoeff3);
  memory->destroy(sf_precoeff4);
  memory->destroy(sf_precoeff5);
  memory->destroy(sf_precoeff6);
  sf_precoeff1 = sf_precoeff2 = sf_precoeff3 = sf_precoeff4 = sf_precoeff5 = sf_precoeff6 = NULL;

  memory->destroy(sf_precoeff1_6);
  memory->destroy(sf_precoeff2_6);
  memory->destroy(sf_precoeff3_6);
  memory->destroy(sf_precoeff4_6);
  memory->destroy(sf_precoeff5_6);
  memory->destroy(sf_precoeff6_6);
  sf_precoeff1_6 = sf_precoeff2_6 = sf_precoeff3_6 = sf_precoeff4_6 = sf_precoeff5_6 = sf_precoeff6_6 = NULL;

  memory->destroy(greensfn);
  memory->destroy(greensfn_6);
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(work1_6);
  memory->destroy(work2_6);
  memory->destroy(vg);
  memory->destroy(vg2);
  memory->destroy(vg_6);
  memory->destroy(vg2_6);
  greensfn = greensfn_6 = NULL;
  work1 = work2 = work1_6 = work2_6 = NULL;
  vg = vg2 = vg_6 = vg2_6 = NULL;

  memory->destroy1d_offset(fkx,nxlo_fft);
  memory->destroy1d_offset(fky,nylo_fft);
  memory->destroy1d_offset(fkz,nzlo_fft);
  fkx = fky = fkz = NULL;

  memory->destroy1d_offset(fkx2,nxlo_fft);
  memory->destroy1d_offset(fky2,nylo_fft);
  memory->destroy1d_offset(fkz2,nzlo_fft);
  fkx2 = fky2 = fkz2 = NULL;

  memory->destroy1d_offset(fkx_6,nxlo_fft_6);
  memory->destroy1d_offset(fky_6,nylo_fft_6);
  memory->destroy1d_offset(fkz_6,nzlo_fft_6);
  fkx_6 = fky_6 = fkz_6 = NULL;

  memory->destroy1d_offset(fkx2_6,nxlo_fft_6);
  memory->destroy1d_offset(fky2_6,nylo_fft_6);
  memory->destroy1d_offset(fkz2_6,nzlo_fft_6);
  fkx2_6 = fky2_6 = fkz2_6 = NULL;


  memory->destroy(gf_b);
  memory->destroy2d_offset(rho1d,-order/2);
  memory->destroy2d_offset(rho_coeff,(1-order)/2);
  memory->destroy2d_offset(drho1d,-order/2);
  memory->destroy2d_offset(drho_coeff, (1-order)/2);
  gf_b = NULL;
  rho1d = rho_coeff = drho1d = drho_coeff = NULL;

  memory->destroy(gf_b_6);
  memory->destroy2d_offset(rho1d_6,-order_6/2);
  memory->destroy2d_offset(rho_coeff_6,(1-order_6)/2);
  memory->destroy2d_offset(drho1d_6,-order_6/2); 
  memory->destroy2d_offset(drho_coeff_6,(1-order_6)/2);
  gf_b_6 = NULL;
  rho1d_6 = rho_coeff_6 = drho1d_6 = drho_coeff_6 = NULL;

  delete fft1;
  delete fft2;
  delete remap;
  delete cg;
  fft1 = fft2 = NULL;
  remap = NULL;
  cg = NULL;

  delete fft1_6;
  delete fft2_6;
  delete remap_6;
  delete cg_6;
  fft1_6 = fft2_6 = NULL;
  remap_6 = NULL;
  cg_6 = NULL;
}


/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
   for per atom calculations 
------------------------------------------------------------------------- */

void PPPMDisp::deallocate_peratom()
{
  peratom_allocate_flag = 0;

  memory->destroy3d_offset(u_brick, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(v0_brick, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(v1_brick, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(v2_brick, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(v3_brick, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(v4_brick, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(v5_brick, nzlo_out, nylo_out, nxlo_out);
  u_brick = v0_brick = v1_brick = v2_brick = v3_brick = v4_brick = v5_brick = NULL;

  memory->destroy3d_offset(u_brick_g, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v0_brick_g, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v1_brick_g, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v2_brick_g, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v3_brick_g, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v4_brick_g, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v5_brick_g, nzlo_out_6, nylo_out_6, nxlo_out_6);
  u_brick_g = v0_brick_g = v1_brick_g = v2_brick_g = v3_brick_g = v4_brick_g = v5_brick_g = NULL;

  memory->destroy3d_offset(u_brick_a0, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v0_brick_a0, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v1_brick_a0, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v2_brick_a0, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v3_brick_a0, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v4_brick_a0, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v5_brick_a0, nzlo_out_6, nylo_out_6, nxlo_out_6);
  u_brick_a0 = v0_brick_a0 = v1_brick_a0 = v2_brick_a0 = v3_brick_a0 = v4_brick_a0 = v5_brick_a0 = NULL;

  memory->destroy3d_offset(u_brick_a1, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v0_brick_a1, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v1_brick_a1, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v2_brick_a1, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v3_brick_a1, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v4_brick_a1, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v5_brick_a1, nzlo_out_6, nylo_out_6, nxlo_out_6);
  u_brick_a1 = v0_brick_a1 = v1_brick_a1 = v2_brick_a1 = v3_brick_a1 = v4_brick_a1 = v5_brick_a1 = NULL;

  memory->destroy3d_offset(u_brick_a2, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v0_brick_a2, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v1_brick_a2, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v2_brick_a2, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v3_brick_a2, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v4_brick_a2, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v5_brick_a2, nzlo_out_6, nylo_out_6, nxlo_out_6);
  u_brick_a2 = v0_brick_a2 = v1_brick_a2 = v2_brick_a2 = v3_brick_a2 = v4_brick_a2 = v5_brick_a2 = NULL;

  memory->destroy3d_offset(u_brick_a3, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v0_brick_a3, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v1_brick_a3, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v2_brick_a3, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v3_brick_a3, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v4_brick_a3, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v5_brick_a3, nzlo_out_6, nylo_out_6, nxlo_out_6);
  u_brick_a3 = v0_brick_a3 = v1_brick_a3 = v2_brick_a3 = v3_brick_a3 = v4_brick_a3 = v5_brick_a3 = NULL;
 
  memory->destroy3d_offset(u_brick_a4, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v0_brick_a4, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v1_brick_a4, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v2_brick_a4, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v3_brick_a4, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v4_brick_a4, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v5_brick_a4, nzlo_out_6, nylo_out_6, nxlo_out_6);
  u_brick_a4 = v0_brick_a4 = v1_brick_a4 = v2_brick_a4 = v3_brick_a4 = v4_brick_a4 = v5_brick_a4 = NULL;
 
  memory->destroy3d_offset(u_brick_a5, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v0_brick_a5, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v1_brick_a5, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v2_brick_a5, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v3_brick_a5, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v4_brick_a5, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v5_brick_a5, nzlo_out_6, nylo_out_6, nxlo_out_6);
  u_brick_a5 = v0_brick_a5 = v1_brick_a5 = v2_brick_a5 = v3_brick_a5 = v4_brick_a5 = v5_brick_a5 = NULL;

  memory->destroy3d_offset(u_brick_a6, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v0_brick_a6, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v1_brick_a6, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v2_brick_a6, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v3_brick_a6, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v4_brick_a6, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy3d_offset(v5_brick_a6, nzlo_out_6, nylo_out_6, nxlo_out_6);
  u_brick_a6 = v0_brick_a6 = v1_brick_a6 = v2_brick_a6 = v3_brick_a6 = v4_brick_a6 = v5_brick_a6 = NULL;

  memory->destroy4d_offset(u_brick_none, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy4d_offset(v0_brick_none, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy4d_offset(v1_brick_none, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy4d_offset(v2_brick_none, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy4d_offset(v3_brick_none, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy4d_offset(v4_brick_none, nzlo_out_6, nylo_out_6, nxlo_out_6);
  memory->destroy4d_offset(v5_brick_none, nzlo_out_6, nylo_out_6, nxlo_out_6);
  u_brick_none = v0_brick_none = v1_brick_none = v2_brick_none = v3_brick_none = v4_brick_none = v5_brick_none = NULL;

  delete cg_peratom;
  delete cg_peratom_6;
  cg_peratom = cg_peratom_6 = NULL;
}

/* ----------------------------------------------------------------------
   set size of FFT grid (nx,ny,nz_pppm) and g_ewald
   for Coulomb interactions
------------------------------------------------------------------------- */

void PPPMDisp::set_grid()
{
  double q2 = qsqsum * force->qqrd2e;

  // use xprd,yprd,zprd even if triclinic so grid size is the same
  // adjust z dimension for 2d slab PPPM
  // 3d PPPM just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;
  
  // make initial g_ewald estimate
  // based on desired accuracy and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  double h, h_x,h_y,h_z;
  bigint natoms = atom->natoms;

  if (!gewaldflag) {
    g_ewald = accuracy*sqrt(natoms*cutoff*xprd*yprd*zprd) / (2.0*q2);
    if (g_ewald >= 1.0)  
      error->all(FLERR,"KSpace accuracy too large to estimate G vector");
    g_ewald = sqrt(-log(g_ewald)) / cutoff;
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

      //set local grid dimension
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

      double qopt = compute_qopt();
   
      double dfkspace = sqrt(qopt/natoms)*q2/(xprd*yprd*zprd_slab);

      count++;

      // break loop if the accuracy has been reached or too many loops have been performed
      if (dfkspace <= accuracy) break;
      if (count > 500) error->all(FLERR, "Could not compute grid size for Coulomb interaction");
      h *= 0.95;
      h_x = h_y = h_z = h;
    }
  }
  
  // boost grid size until it is factorable

  while (!factorable(nx_pppm)) nx_pppm++;
  while (!factorable(ny_pppm)) ny_pppm++;
  while (!factorable(nz_pppm)) nz_pppm++;
}

/* ----------------------------------------------------------------------
   set the FFT parameters 
------------------------------------------------------------------------- */

void PPPMDisp::set_fft_parameters(int& nx_p,int& ny_p,int& nz_p,
                                   int& nxlo_f,int& nylo_f,int& nzlo_f,
                                   int& nxhi_f,int& nyhi_f,int& nzhi_f,
                                   int& nxlo_i,int& nylo_i,int& nzlo_i,
                                   int& nxhi_i,int& nyhi_i,int& nzhi_i,
                                   int& nxlo_o,int& nylo_o,int& nzlo_o,
                                   int& nxhi_o,int& nyhi_o,int& nzhi_o,
		                   int& nlow, int& nupp,
                                   int& ng, int& nf, int& nfb,
		                   double& sft,double& sftone, int& ord)
{
  // global indices of PPPM grid range from 0 to N-1
  // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
  //   global PPPM grid that I own without ghost cells
  // for slab PPPM, assign z grid as if it were not extended

  nxlo_i = static_cast<int> (comm->xsplit[comm->myloc[0]] * nx_p);
  nxhi_i = static_cast<int> (comm->xsplit[comm->myloc[0]+1] * nx_p) - 1;

  nylo_i = static_cast<int> (comm->ysplit[comm->myloc[1]] * ny_p);
  nyhi_i = static_cast<int> (comm->ysplit[comm->myloc[1]+1] * ny_p) - 1;

  nzlo_i = static_cast<int> 
      (comm->zsplit[comm->myloc[2]] * nz_p/slab_volfactor);
  nzhi_i = static_cast<int> 
      (comm->zsplit[comm->myloc[2]+1] * nz_p/slab_volfactor) - 1;


  // nlow,nupp = stencil size for mapping particles to PPPM grid

  nlow = -(ord-1)/2;
  nupp = ord/2;

  // sft values for particle <-> grid mapping
  // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

  if (ord % 2) sft = OFFSET + 0.5;
  else sft = OFFSET;
  if (ord % 2) sftone = 0.0;
  else sftone = 0.5;

  // nlo_out,nhi_out = lower/upper limits of the 3d sub-brick of
  //   global PPPM grid that my particles can contribute charge to
  // effectively nlo_in,nhi_in + ghost cells
  // nlo,nhi = global coords of grid pt to "lower left" of smallest/largest
  //           position a particle in my box can be at
  // dist[3] = particle position bound = subbox + skin/2.0 + qdist
  //   qdist = offset due to TIP4P fictitious charge
  //   convert to triclinic if necessary
  // nlo_out,nhi_out = nlo,nhi + stencil size for particle mapping
  // for slab PPPM, assign z grid as if it were not extended

  double *prd,*sublo,*subhi;

  if (triclinic == 0) {
    prd = domain->prd;
    boxlo = domain->boxlo;
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    prd = domain->prd_lamda;
    boxlo = domain->boxlo_lamda;
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double dist[3];
  double cuthalf = 0.5*neighbor->skin + qdist;
  if (triclinic == 0) dist[0] = dist[1] = dist[2] = cuthalf;
  else {
    dist[0] = cuthalf/domain->prd[0];
    dist[1] = cuthalf/domain->prd[1];
    dist[2] = cuthalf/domain->prd[2];
  }
    
  int nlo,nhi;
    
  nlo = static_cast<int> ((sublo[0]-dist[0]-boxlo[0]) * 
                            nx_p/xprd + sft) - OFFSET;
  nhi = static_cast<int> ((subhi[0]+dist[0]-boxlo[0]) * 
                            nx_p/xprd + sft) - OFFSET;
  nxlo_o = nlo + nlow;
  nxhi_o = nhi + nupp;

  nlo = static_cast<int> ((sublo[1]-dist[1]-boxlo[1]) * 
                            ny_p/yprd + sft) - OFFSET;
  nhi = static_cast<int> ((subhi[1]+dist[1]-boxlo[1]) * 
                            ny_p/yprd + sft) - OFFSET;
  nylo_o = nlo + nlow;
  nyhi_o = nhi + nupp;

  nlo = static_cast<int> ((sublo[2]-dist[2]-boxlo[2]) * 
                            nz_p/zprd_slab + sft) - OFFSET;
  nhi = static_cast<int> ((subhi[2]+dist[2]-boxlo[2]) * 
                            nz_p/zprd_slab + sft) - OFFSET;
  nzlo_o = nlo + nlow;
  nzhi_o = nhi + nupp;

  // for slab PPPM, change the grid boundary for processors at +z end
  //   to include the empty volume between periodically repeating slabs
  // for slab PPPM, want charge data communicated from -z proc to +z proc,
  //   but not vice versa, also want field data communicated from +z proc to
  //   -z proc, but not vice versa
  // this is accomplished by nzhi_i = nzhi_o on +z end (no ghost cells)

  if (slabflag && (comm->myloc[2] == comm->procgrid[2]-1)) {
    nzhi_i = nz_p - 1;
    nzhi_o = nz_p - 1;
  }
  
  // decomposition of FFT mesh
  // global indices range from 0 to N-1
  // proc owns entire x-dimension, clump of columns in y,z dimensions
  // npey_fft,npez_fft = # of procs in y,z dims
  // if nprocs is small enough, proc can own 1 or more entire xy planes,
  //   else proc owns 2d sub-blocks of yz plane
  // me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
  // nlo_fft,nhi_fft = lower/upper limit of the section
  //   of the global FFT mesh that I own

  int npey_fft,npez_fft;
  if (nz_p >= nprocs) {
    npey_fft = 1;
    npez_fft = nprocs;
  } else procs2grid2d(nprocs,ny_p,nz_p,&npey_fft,&npez_fft);

  int me_y = me % npey_fft;
  int me_z = me / npey_fft;

  nxlo_f = 0;
  nxhi_f = nx_p - 1;
  nylo_f = me_y*ny_p/npey_fft;
  nyhi_f = (me_y+1)*ny_p/npey_fft - 1;
  nzlo_f = me_z*nz_p/npez_fft;
  nzhi_f = (me_z+1)*nz_p/npez_fft - 1;

  // PPPM grid for this proc, including ghosts

  ng = (nxhi_o-nxlo_o+1) * (nyhi_o-nylo_o+1) *
    (nzhi_o-nzlo_o+1);

  // FFT arrays on this proc, without ghosts
  // nfft = FFT points in FFT decomposition on this proc
  // nfft_brick = FFT points in 3d brick-decomposition on this proc
  // nfft_both = greater of 2 values

  nf = (nxhi_f-nxlo_f+1) * (nyhi_f-nylo_f+1) *
    (nzhi_f-nzlo_f+1);
  int nfft_brick = (nxhi_i-nxlo_i+1) * (nyhi_i-nylo_i+1) *
    (nzhi_i-nzlo_i+1);
  nfb = MAX(nf,nfft_brick);

}

/* ----------------------------------------------------------------------
   check if all factors of n are in list of factors
   return 1 if yes, 0 if no 
------------------------------------------------------------------------- */

int PPPMDisp::factorable(int n)
{
  int i;

  while (n > 1) {
    for (i = 0; i < nfactors; i++) {
      if (n % factors[i] == 0) {
	n /= factors[i];
	break;
      }
    }
    if (i == nfactors) return 0;
  }

  return 1;
}

/* ----------------------------------------------------------------------
   pre-compute Green's function denominator expansion coeffs, Gamma(2n) 
------------------------------------------------------------------------- */
void PPPMDisp::adjust_gewald()
{
  
  // Use Newton solver to find g_ewald

  double dx;
        
  // Begin algorithm
  
  for (int i = 0; i < LARGE; i++) {
    dx = f() / derivf(); 
    g_ewald -= dx; //Update g_ewald
    if (fabs(f()) < SMALL) return;
  }
   
  // Failed to converge
  
  char str[128];
  sprintf(str, "Could not compute g_ewald");
  error->all(FLERR, str);

}

/* ----------------------------------------------------------------------
 Calculate f(x)
 ------------------------------------------------------------------------- */

double PPPMDisp::f()
{
  double df_rspace, df_kspace;
  double q2 = qsqsum * force->qqrd2e;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;
  bigint natoms = atom->natoms;

  df_rspace = 2.0*q2*exp(-g_ewald*g_ewald*cutoff*cutoff) / 
       sqrt(natoms*cutoff*xprd*yprd*zprd);
   
  double qopt = compute_qopt();
  df_kspace = sqrt(qopt/natoms)*q2/(xprd*yprd*zprd_slab);
   
  return df_rspace - df_kspace;
}

/* ----------------------------------------------------------------------
 Calculate numerical derivative f'(x) using forward difference
 [f(x + h) - f(x)] / h
 ------------------------------------------------------------------------- */
            
double PPPMDisp::derivf()
{  
  double h = 0.000001;  //Derivative step-size
  double df,f1,f2,g_ewald_old;
  
  f1 = f();
  g_ewald_old = g_ewald;
  g_ewald += h;
  f2 = f();
  g_ewald = g_ewald_old;
  df = (f2 - f1)/h;
  
  return df;
} 

/* ----------------------------------------------------------------------
   Calculate the final estimator for the accuracy
------------------------------------------------------------------------- */

double PPPMDisp::final_accuracy()
{
  double df_rspace, df_kspace;
  double q2 = qsqsum * force->qqrd2e;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;
  bigint natoms = atom->natoms;
  df_rspace = 2.0*q2 * exp(-g_ewald*g_ewald*cutoff*cutoff) / 
             sqrt(natoms*cutoff*xprd*yprd*zprd);

  double qopt = compute_qopt();

  df_kspace = sqrt(qopt/natoms)*q2/(xprd*yprd*zprd_slab);

  double acc = sqrt(df_rspace*df_rspace + df_kspace*df_kspace);
  return acc;
}

/* ----------------------------------------------------------------------
   Calculate the final estimator for the Dispersion accuracy
------------------------------------------------------------------------- */

void PPPMDisp::final_accuracy_6(double& acc, double& acc_real, double& acc_kspace)
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;
  bigint natoms = atom->natoms;
  acc_real = lj_rspace_error();

  double qopt = compute_qopt_6();

  acc_kspace = sqrt(qopt/natoms)*csum/(xprd*yprd*zprd_slab);

  acc = sqrt(acc_real*acc_real + acc_kspace*acc_kspace);
  return;
}

/* ----------------------------------------------------------------------
   Compute qopt for Coulomb interactions
------------------------------------------------------------------------- */

double PPPMDisp::compute_qopt()
{
  double qopt;
  if (differentiation_flag == 1) {
    qopt = compute_qopt_ad();
  } else {
    qopt = compute_qopt_ik();
  }
  double qopt_all;
  MPI_Allreduce(&qopt,&qopt_all,1,MPI_DOUBLE,MPI_SUM,world);
  return qopt_all;
}

/* ----------------------------------------------------------------------
   Compute qopt for Dispersion interactions
------------------------------------------------------------------------- */

double PPPMDisp::compute_qopt_6()
{
  double qopt;
  if (differentiation_flag == 1) {
    qopt = compute_qopt_6_ad();
  } else {
    qopt = compute_qopt_6_ik();
  }
  double qopt_all;
  MPI_Allreduce(&qopt,&qopt_all,1,MPI_DOUBLE,MPI_SUM,world);
  return qopt_all;
}

/* ----------------------------------------------------------------------
   Compute qopt for the ik differentiation scheme and Coulomb interaction
------------------------------------------------------------------------- */

double PPPMDisp::compute_qopt_ik()
{
  double qopt = 0.0;
  int k,l,m;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int nx,ny,nz,kper,lper,mper;
  double sqk, u2;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2, sum3,dot1,dot2;

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);
      
        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) + 
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);
            sx = exp(-0.25*pow(qx/g_ewald,2.0));
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order);
            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);
              sy = exp(-0.25*pow(qy/g_ewald,2.0));
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order);
              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);
                sz = exp(-0.25*pow(qz/g_ewald,2.0));
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order);

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                dot2 = qx*qx+qy*qy+qz*qz;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += sx*sy*sz*sx*sy*sz/dot2*4.0*4.0*MY_PI*MY_PI;
                sum2 += u2*sx*sy*sz*4.0*MY_PI/dot2*dot1;
		sum3 += u2;
              }
            }
          }
	  sum2 *= sum2;
	  sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
        }
      }
    }
  }
  return qopt;
}

/* ----------------------------------------------------------------------
   Compute qopt for the ad differentiation scheme and Coulomb interaction
------------------------------------------------------------------------- */

double PPPMDisp::compute_qopt_ad()
{
  double qopt = 0.0;
  int k,l,m;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;


  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int nx,ny,nz,kper,lper,mper;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double u2, sqk;
  double sum1,sum2,sum3,sum4,dot2;

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);
      
        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) + 
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {
    
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          sum4 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);
            sx = exp(-0.25*pow(qx/g_ewald,2.0));
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order);
            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);
              sy = exp(-0.25*pow(qy/g_ewald,2.0));
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order);
              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);
                sz = exp(-0.25*pow(qz/g_ewald,2.0));
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order);

                dot2 = qx*qx+qy*qy+qz*qz;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += sx*sy*sz*sx*sy*sz/dot2*4.0*4.0*MY_PI*MY_PI;
                sum2 += sx*sy*sz * u2*4.0*MY_PI;
                sum3 += u2;
                sum4 += dot2*u2;
              }
            }
          }
          sum2 *= sum2;
          qopt += sum1 - sum2/(sum3*sum4);
        }
      }
    }
  }
  return qopt;
}

/* ----------------------------------------------------------------------
   Compute qopt for the ik differentiation scheme and Dispersion interaction
------------------------------------------------------------------------- */

double PPPMDisp::compute_qopt_6_ik()
{
  double qopt = 0.0;
  int k,l,m;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int nx,ny,nz,kper,lper,mper;
  double sqk, u2;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2, sum3;
  double dot1,dot2, rtdot2, term;
  double inv2ew = 2*g_ewald_6;
  inv2ew = 1.0/inv2ew;
  double rtpi = sqrt(MY_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) {
    mper = m - nz_pppm_6*(2*m/nz_pppm_6);

    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      lper = l - ny_pppm_6*(2*l/ny_pppm_6);

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - nx_pppm_6*(2*k/nx_pppm_6);
      
        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) + 
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm_6*nx);
            sx = exp(-qx*qx*inv2ew*inv2ew);
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm_6;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order_6);
            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm_6*ny);
              sy = exp(-qy*qy*inv2ew*inv2ew);
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm_6;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order_6);
              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm_6*nz);
                sz = exp(-qz*qz*inv2ew*inv2ew);
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm_6;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order_6);

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                dot2 = qx*qx+qy*qy+qz*qz;
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew*inv2ew)*sx*sy*sz +
		       2*dot2*rtdot2*inv2ew*inv2ew*inv2ew*rtpi*erfc(rtdot2*inv2ew);
                term *= g_ewald_6*g_ewald_6*g_ewald_6;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += term*term*MY_PI*MY_PI*MY_PI/9.0 * dot2;
                sum2 += -u2*term*MY_PI*rtpi/3.0*dot1;
		sum3 += u2;
              }
            }
          }
	  sum2 *= sum2;
	  sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
        }
      }
    }
  }
  return qopt;
}

/* ----------------------------------------------------------------------
   Compute qopt for the ad differentiation scheme and Dispersion interaction
------------------------------------------------------------------------- */

double PPPMDisp::compute_qopt_6_ad()
{
  double qopt = 0.0;
  int k,l,m;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int nx,ny,nz,kper,lper,mper;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double u2, sqk;
  double sum1,sum2,sum3,sum4;
  double dot2, rtdot2, term;
  double inv2ew = 2*g_ewald_6;
  inv2ew = 1/inv2ew;
  double rtpi = sqrt(MY_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) {
    mper = m - nz_pppm_6*(2*m/nz_pppm_6);

    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      lper = l - ny_pppm_6*(2*l/ny_pppm_6);

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - nx_pppm_6*(2*k/nx_pppm_6);
      
        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) + 
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {
    
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          sum4 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm_6*nx);
            sx = exp(-qx*qx*inv2ew*inv2ew);
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm_6;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order_6);
            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm_6*ny);
              sy = exp(-qy*qy*inv2ew*inv2ew);
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm_6;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order_6);
              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm_6*nz);
                sz = exp(-qz*qz*inv2ew*inv2ew);
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm_6;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order_6);

                dot2 = qx*qx+qy*qy+qz*qz;
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew*inv2ew)*sx*sy*sz +
		       2*dot2*rtdot2*inv2ew*inv2ew*inv2ew*rtpi*erfc(rtdot2*inv2ew);
                term *= g_ewald_6*g_ewald_6*g_ewald_6;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += term*term*MY_PI*MY_PI*MY_PI/9.0 * dot2;
                sum2 += -term*MY_PI*rtpi/3.0 * u2 * dot2;
                sum3 += u2;
                sum4 += dot2*u2;
              }
            }
          }
          sum2 *= sum2;
          qopt += sum1 - sum2/(sum3*sum4);
        }
      }
    }
  }
  return qopt;
}

/* ----------------------------------------------------------------------
   set size of FFT grid  and g_ewald_6
   for Dispersion interactions
------------------------------------------------------------------------- */

void PPPMDisp::set_grid_6()
{
  // Calculate csum
  if (!csumflag) calc_csum();
  if (!gewaldflag_6) set_init_g6();
  if (!gridflag_6) set_n_pppm_6();
  while (!factorable(nx_pppm_6)) nx_pppm_6++;
  while (!factorable(ny_pppm_6)) ny_pppm_6++;
  while (!factorable(nz_pppm_6)) nz_pppm_6++;
  
}

/* ----------------------------------------------------------------------
   Calculate the sum of the squared dispersion coefficients and other 
   related quantities required for the calculations
------------------------------------------------------------------------- */

void PPPMDisp::calc_csum()
{
  csumij = 0.0;
  csum = 0.0;

  int ntypes = atom->ntypes;   
  int i,j,k;

  delete [] cii;
  cii = new double[ntypes +1];
  for (i = 0; i<=ntypes; i++) cii[i] = 0.0;
  delete [] csumi; 
  csumi = new double[ntypes +1];
  for (i = 0; i<=ntypes; i++) csumi[i] = 0.0; 
  int *neach = new int[ntypes+1];
  for (i = 0; i<=ntypes; i++) neach[i] = 0; 

  //the following variables are needed to distinguish between arithmetic
  //  and geometric mixing

  if (function[1]) {
    for (i = 1; i <= ntypes; i++)
      cii[i] = B[i]*B[i];
    int tmp;
    for (i = 0; i < atom->nlocal; i++) {
      tmp = atom->type[i];
      neach[tmp]++;
      csum += B[tmp]*B[tmp];
    }
  }
  if (function[2]) {
    for (i = 1; i <= ntypes; i++)
      cii[i] = 64.0/20.0*B[7*i+3]*B[7*i+3];
    int tmp;
    for (i = 0; i < atom->nlocal; i++) {
      tmp = atom->type[i];
      neach[tmp]++;
      csum += 64.0/20.0*B[7*tmp+3]*B[7*tmp+3];
    }
  }
  if (function[3]) {
    for (i = 1; i <= ntypes; i++)
      for (j = 0; j < nsplit; j++)
        cii[i] += B[j]*B[nsplit*i + j]*B[nsplit*i + j];
    int tmp;
    for (i = 0; i < atom->nlocal; i++) {
      tmp = atom->type[i];
      neach[tmp]++;
      for (j = 0; j < nsplit; j++)
        csum += B[j]*B[nsplit*tmp + j]*B[nsplit*tmp + j];
    }
  }


  double tmp2;
  MPI_Allreduce(&csum,&tmp2,1,MPI_DOUBLE,MPI_SUM,world);
  csum = tmp2;
  csumflag = 1;

  int *neach_all = new int[ntypes+1];
  MPI_Allreduce(neach,neach_all,ntypes+1,MPI_INT,MPI_SUM,world);

  // copmute csumij and csumi
  double d1, d2;
  if (function[1]){
    for (i=1; i<=ntypes; i++) {
      for (j=1; j<=ntypes; j++) {
        csumi[i] += neach_all[j]*B[i]*B[j];
        d1 = neach_all[i]*B[i];
        d2 = neach_all[j]*B[j];
        csumij += d1*d2;
        //csumij += neach_all[i]*neach_all[j]*B[i]*B[j]; 
      }
    }
  }
  if (function[2]) {
    for (i=1; i<=ntypes; i++) {
      for (j=1; j<=ntypes; j++) {
        for (k=0; k<=6; k++) {
          csumi[i] += neach_all[j]*B[7*i + k]*B[7*(j+1)-k-1];
          d1 = neach_all[i]*B[7*i + k];
          d2 = neach_all[j]*B[7*(j+1)-k-1];
          csumij += d1*d2;
          //csumij += neach_all[i]*neach_all[j]*B[7*i + k]*B[7*(j+1)-k-1];
        }
      }
    }
  }
  if (function[3]) {
    for (i=1; i<=ntypes; i++) {
      for (j=1; j<=ntypes; j++) {
        for (k=0; k<nsplit; k++) {
	  csumi[i] += neach_all[j]*B[k]*B[nsplit*i+k]*B[nsplit*j+k];
	  d1 = neach_all[i]*B[nsplit*i+k];
	  d2 = neach_all[j]*B[nsplit*j+k];
          csumij += B[k]*d1*d2;
	}
      }
    }
  }    

  delete [] neach;
  delete [] neach_all;
}

/* ----------------------------------------------------------------------
   adjust g_ewald_6 to the new grid size
------------------------------------------------------------------------- */

void PPPMDisp::adjust_gewald_6()
{
  // Use Newton solver to find g_ewald_6
  double dx;

  // Start loop

  for (int i = 0; i <  LARGE; i++) {
    dx = f_6() / derivf_6();
    g_ewald_6 -= dx; //update g_ewald_6
    if (fabs(f_6()) < SMALL) return;
  }

  // Failed to converge

  char str[128];
  sprintf(str, "Could not adjust g_ewald_6");
  error->all(FLERR, str);

}

/* ----------------------------------------------------------------------
 Calculate f(x) for Dispersion interaction
 ------------------------------------------------------------------------- */

double PPPMDisp::f_6()
{
  double df_rspace, df_kspace;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  bigint natoms = atom->natoms;

  df_rspace = lj_rspace_error();
   
  double qopt = compute_qopt_6();
  df_kspace = sqrt(qopt/natoms)*csum/(xprd*yprd*zprd_slab);
   
  return df_rspace - df_kspace;
}

/* ----------------------------------------------------------------------
 Calculate numerical derivative f'(x) using forward difference
 [f(x + h) - f(x)] / h
 ------------------------------------------------------------------------- */
            
double PPPMDisp::derivf_6()
{  
  double h = 0.000001;  //Derivative step-size
  double df,f1,f2,g_ewald_old;
  
  f1 = f_6();
  g_ewald_old = g_ewald_6;
  g_ewald_6 += h;
  f2 = f_6();
  g_ewald_6 = g_ewald_old;
  df = (f2 - f1)/h;
  
  return df;
} 


/* ----------------------------------------------------------------------
   calculate an initial value for g_ewald_6
   ---------------------------------------------------------------------- */

void PPPMDisp::set_init_g6()
{
  // use xprd,yprd,zprd even if triclinic so grid size is the same
  // adjust z dimension for 2d slab PPPM
  // 3d PPPM just uses zprd since slab_volfactor = 1.0

  // make initial g_ewald estimate
  // based on desired error and real space cutoff
 
  // compute initial value for df_real with g_ewald_6 = 1/cutoff_lj
  // if df_real > 0, repeat divide g_ewald_6 by 2 until df_real < 0
  // else, repeat multiply g_ewald_6 by 2 until df_real > 0
  // perform bisection for the last two values of
  double df_real;
  double g_ewald_old; 
  double gmin, gmax;

  // check if there is a user defined accuracy
  double acc_rspace = accuracy;
  if (accuracy_real_6 > 0) acc_rspace = accuracy_real_6;

  g_ewald_6 = 1.0/cutoff_lj;
  df_real = lj_rspace_error() - acc_rspace;
  int counter = 0;
  if (df_real > 0) {
    while (df_real > 0 && counter < LARGE) {
      counter++;
      g_ewald_old = g_ewald_6;
      g_ewald_6 *= 2;
      df_real = lj_rspace_error() - acc_rspace;
    }
  }

  if (df_real < 0) {
    while (df_real < 0 && counter < LARGE) {
      counter++;
      g_ewald_old = g_ewald_6;
      g_ewald_6 *= 0.5;
      df_real = lj_rspace_error() - acc_rspace;
    }
  }

  if (counter >= LARGE-1) error->all(FLERR,"Cannot compute initial g_ewald_disp");

  gmin = MIN(g_ewald_6, g_ewald_old);
  gmax = MAX(g_ewald_6, g_ewald_old);
  g_ewald_6 = gmin + 0.5*(gmax-gmin);
  counter = 0;
  while (gmax-gmin > SMALL && counter < LARGE) {
    counter++;
    df_real = lj_rspace_error() -acc_rspace;
    if (df_real < 0) gmax = g_ewald_6;
    else gmin = g_ewald_6;
    g_ewald_6 = gmin + 0.5*(gmax-gmin);
  }
  if (counter >= LARGE-1) error->all(FLERR,"Cannot compute initial g_ewald_disp");

}

/* ----------------------------------------------------------------------
   calculate nx_pppm, ny_pppm, nz_pppm for dispersion interaction
   ---------------------------------------------------------------------- */

void PPPMDisp::set_n_pppm_6()
{
  bigint natoms = atom->natoms;

  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  double h, h_x,h_y,h_z;

  double acc_kspace = accuracy;
  if (accuracy_kspace_6 > 0.0) acc_kspace = accuracy_kspace_6;

  // initial value for the grid spacing
  h = h_x = h_y = h_z = 4.0/g_ewald_6;
  // decrease grid spacing untill required precision is obtained
  int count = 0;
  while(1) {
  
    // set grid dimension
    nx_pppm_6 = static_cast<int> (xprd/h_x);
    ny_pppm_6 = static_cast<int> (yprd/h_y);
    nz_pppm_6 = static_cast<int> (zprd_slab/h_z);

    if (nx_pppm_6 <= 1) nx_pppm_6 = 2;
    if (ny_pppm_6 <= 1) ny_pppm_6 = 2;
    if (nz_pppm_6 <= 1) nz_pppm_6 = 2;

    //set local grid dimension
    int npey_fft,npez_fft;
    if (nz_pppm_6 >= nprocs) {
      npey_fft = 1;
      npez_fft = nprocs;
    } else procs2grid2d(nprocs,ny_pppm_6,nz_pppm_6,&npey_fft,&npez_fft);

    int me_y = me % npey_fft;
    int me_z = me / npey_fft;

    nxlo_fft_6 = 0;
    nxhi_fft_6 = nx_pppm_6 - 1;
    nylo_fft_6 = me_y*ny_pppm_6/npey_fft;
    nyhi_fft_6 = (me_y+1)*ny_pppm_6/npey_fft - 1;
    nzlo_fft_6 = me_z*nz_pppm_6/npez_fft;
    nzhi_fft_6 = (me_z+1)*nz_pppm_6/npez_fft - 1;

    double qopt = compute_qopt_6();
 
    double df_kspace = sqrt(qopt/natoms)*csum/(xprd*yprd*zprd_slab);

    count++;

    // break loop if the accuracy has been reached or too many loops have been performed
    if (df_kspace <= acc_kspace) break;
    if (count > 500) error->all(FLERR, "Could not compute grid size for Dispersion");
    h *= 0.95;
    h_x = h_y = h_z = h;
  }
}

/* ----------------------------------------------------------------------
   calculate the real space error for dispersion interactions
   ---------------------------------------------------------------------- */

double PPPMDisp::lj_rspace_error()
{
  bigint natoms = atom->natoms;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;

  double deltaf;
  double rgs = (cutoff_lj*g_ewald_6);
  rgs *= rgs;
  double rgs_inv = 1.0/rgs;
  deltaf = csum/sqrt(natoms*xprd*yprd*zprd_slab*cutoff_lj)*sqrt(MY_PI)*pow(g_ewald_6, 5)*
    exp(-rgs)*(1+rgs_inv*(3+rgs_inv*(6+rgs_inv*6)));
  return deltaf;
}


/* ----------------------------------------------------------------------
   Compyute the modified (hockney-eastwood) coulomb green function
   ---------------------------------------------------------------------- */ 

void PPPMDisp::compute_gf()
{
  int k,l,m,n;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int kper,lper,mper;
  double snx,sny,snz,snx2,sny2,snz2;
  double sqk;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double numerator,denominator;


  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);
    qz = unitkz*mper;
    snz = sin(0.5*qz*zprd_slab/nz_pppm);
    snz2 = snz*snz;
    sz = exp(-0.25*pow(qz/g_ewald,2.0));
    wz = 1.0;
    argz = 0.5*qz*zprd_slab/nz_pppm;
    if (argz != 0.0) wz = pow(sin(argz)/argz,order);
    wz *= wz;

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);
      qy = unitky*lper;
      sny = sin(0.5*qy*yprd/ny_pppm);
      sny2 = sny*sny;
      sy = exp(-0.25*pow(qy/g_ewald,2.0));
      wy = 1.0;
      argy = 0.5*qy*yprd/ny_pppm;
      if (argy != 0.0) wy = pow(sin(argy)/argy,order);
      wy *= wy;

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);
        qx = unitkx*kper;
        snx = sin(0.5*qx*xprd/nx_pppm);
        snx2 = snx*snx;
        sx = exp(-0.25*pow(qx/g_ewald,2.0));
        wx = 1.0;
        argx = 0.5*qx*xprd/nx_pppm;
        if (argx != 0.0) wx = pow(sin(argx)/argx,order);
        wx *= wx;

        sqk = pow(qx,2.0) + pow(qy,2.0) + pow(qz,2.0);

        if (sqk != 0.0) {
          numerator = 4.0*MY_PI/sqk;
          denominator = gf_denom(snx2,sny2,snz2, gf_b, order);  
          greensfn[n++] = numerator*sx*sy*sz*wx*wy*wz/denominator;
        } else greensfn[n++] = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute self force coefficients for ad-differentiation scheme
   and Coulomb interaction 
------------------------------------------------------------------------- */

void PPPMDisp::compute_sf_precoeff(int nxp, int nyp, int nzp, int ord, 
                                    int nxlo_ft, int nylo_ft, int nzlo_ft,
                                    int nxhi_ft, int nyhi_ft, int nzhi_ft,
                                    double *sf_pre1, double *sf_pre2, double *sf_pre3,
                                    double *sf_pre4, double *sf_pre5, double *sf_pre6)
{

  int i,k,l,m,n;
  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPM
  // z dimension for 3d PPPM is zprd since slab_volfactor = 1.0

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int nx,ny,nz,kper,lper,mper;
  double argx,argy,argz;
  double wx0[5],wy0[5],wz0[5],wx1[5],wy1[5],wz1[5],wx2[5],wy2[5],wz2[5];
  double qx0,qy0,qz0,qx1,qy1,qz1,qx2,qy2,qz2;
  double u0,u1,u2,u3,u4,u5,u6;
  double sum1,sum2,sum3,sum4,sum5,sum6;

  int nb = 2;

  n = 0;
  for (m = nzlo_ft; m <= nzhi_ft; m++) {
    mper = m - nzp*(2*m/nzp);

    for (l = nylo_ft; l <= nyhi_ft; l++) {
      lper = l - nyp*(2*l/nyp);

      for (k = nxlo_ft; k <= nxhi_ft; k++) {
        kper = k - nxp*(2*k/nxp);
      
        sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = 0.0;
        for (i = -nb; i <= nb; i++) {

          qx0 = unitkx*(kper+nxp*i);
          qx1 = unitkx*(kper+nxp*(i+1));
          qx2 = unitkx*(kper+nxp*(i+2));
          wx0[i+2] = 1.0;
          wx1[i+2] = 1.0;
          wx2[i+2] = 1.0;
          argx = 0.5*qx0*xprd/nxp;
          if (argx != 0.0) wx0[i+2] = pow(sin(argx)/argx,ord);
          argx = 0.5*qx1*xprd/nxp;
          if (argx != 0.0) wx1[i+2] = pow(sin(argx)/argx,ord);
          argx = 0.5*qx2*xprd/nxp;
          if (argx != 0.0) wx2[i+2] = pow(sin(argx)/argx,ord);

          qy0 = unitky*(lper+nyp*i);
          qy1 = unitky*(lper+nyp*(i+1));
          qy2 = unitky*(lper+nyp*(i+2));
          wy0[i+2] = 1.0;
          wy1[i+2] = 1.0;
          wy2[i+2] = 1.0;
          argy = 0.5*qy0*yprd/nyp;
          if (argy != 0.0) wy0[i+2] = pow(sin(argy)/argy,ord);
          argy = 0.5*qy1*yprd/nyp;
          if (argy != 0.0) wy1[i+2] = pow(sin(argy)/argy,ord);
          argy = 0.5*qy2*yprd/nyp;
          if (argy != 0.0) wy2[i+2] = pow(sin(argy)/argy,ord);
   
          qz0 = unitkz*(mper+nzp*i);
          qz1 = unitkz*(mper+nzp*(i+1));
          qz2 = unitkz*(mper+nzp*(i+2));
          wz0[i+2] = 1.0;
          wz1[i+2] = 1.0;
          wz2[i+2] = 1.0;
          argz = 0.5*qz0*zprd_slab/nzp;
          if (argz != 0.0) wz0[i+2] = pow(sin(argz)/argz,ord);
          argz = 0.5*qz1*zprd_slab/nzp;
          if (argz != 0.0) wz1[i+2] = pow(sin(argz)/argz,ord);
           argz = 0.5*qz2*zprd_slab/nzp;
          if (argz != 0.0) wz2[i+2] = pow(sin(argz)/argz,ord);
        }
    
        for (nx = 0; nx <= 4; nx++) {
          for (ny = 0; ny <= 4; ny++) {
            for (nz = 0; nz <= 4; nz++) {
              u0 = wx0[nx]*wy0[ny]*wz0[nz];
              u1 = wx1[nx]*wy0[ny]*wz0[nz];
              u2 = wx2[nx]*wy0[ny]*wz0[nz];
              u3 = wx0[nx]*wy1[ny]*wz0[nz];
              u4 = wx0[nx]*wy2[ny]*wz0[nz];
              u5 = wx0[nx]*wy0[ny]*wz1[nz];
              u6 = wx0[nx]*wy0[ny]*wz2[nz];

              sum1 += u0*u1;
              sum2 += u0*u2;
              sum3 += u0*u3;
              sum4 += u0*u4;
              sum5 += u0*u5;
              sum6 += u0*u6;
            }
          }
        }
        
        // store values

        sf_pre1[n] = sum1;
        sf_pre2[n] = sum2;
        sf_pre3[n] = sum3;
        sf_pre4[n] = sum4;
        sf_pre5[n] = sum5;
        sf_pre6[n++] = sum6;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Compute the modified (hockney-eastwood) dispersion green function
   ---------------------------------------------------------------------- */

void PPPMDisp::compute_gf_6()
{
  double *prd;
  int k,l,m,n;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPM
  // z dimension for 3d PPPM is zprd since slab_volfactor = 1.0

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int kper,lper,mper;
  double sqk;
  double snx,sny,snz,snx2,sny2,snz2;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz;
  double qx,qy,qz;
  double rtsqk, term;
  double numerator,denominator;
  double inv2ew = 2*g_ewald_6;
  inv2ew = 1/inv2ew;
  double rtpi = sqrt(MY_PI);

  numerator = -MY_PI*rtpi*g_ewald_6*g_ewald_6*g_ewald_6/(3.0);

  n = 0;
  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) {
    mper = m - nz_pppm_6*(2*m/nz_pppm_6);
    qz = unitkz*mper;
    snz = sin(0.5*unitkz*mper*zprd_slab/nz_pppm_6);
    snz2 = snz*snz;
    sz = exp(-qz*qz*inv2ew*inv2ew);
    wz = 1.0;
    argz = 0.5*qz*zprd_slab/nz_pppm_6;
    if (argz != 0.0) wz = pow(sin(argz)/argz,order_6);
    wz *= wz;
              
    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      lper = l - ny_pppm_6*(2*l/ny_pppm_6);
      qy = unitky*lper;
      sny = sin(0.5*unitky*lper*yprd/ny_pppm_6);
      sny2 = sny*sny;
      sy = exp(-qy*qy*inv2ew*inv2ew);
      wy = 1.0;
      argy = 0.5*qy*yprd/ny_pppm_6;
      if (argy != 0.0) wy = pow(sin(argy)/argy,order_6);
      wy *= wy;

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
	kper = k - nx_pppm_6*(2*k/nx_pppm_6);
        qx = unitkx*kper;
	snx = sin(0.5*unitkx*kper*xprd/nx_pppm_6);
	snx2 = snx*snx;
        sx = exp(-qx*qx*inv2ew*inv2ew);
	wx = 1.0;
	argx = 0.5*qx*xprd/nx_pppm_6;
	if (argx != 0.0) wx = pow(sin(argx)/argx,order_6);
        wx *= wx;
      
	sqk = pow(qx,2.0) + pow(qy,2.0) + pow(qz,2.0);

        if (sqk != 0.0) {
	  denominator = gf_denom(snx2,sny2,snz2, gf_b_6, order_6); 
	  rtsqk = sqrt(sqk);
          term = (1-2*sqk*inv2ew*inv2ew)*sx*sy*sz +
                  2*sqk*rtsqk*inv2ew*inv2ew*inv2ew*rtpi*erfc(rtsqk*inv2ew);
	  greensfn_6[n++] = numerator*term*wx*wy*wz/denominator;
        } else greensfn_6[n++] = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute self force coefficients for ad-differentiation scheme
   and Coulomb interaction 
------------------------------------------------------------------------- */
void PPPMDisp::compute_sf_coeff()
{
  int i,k,l,m,n;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  for (i = 0; i <= 5; i++) sf_coeff[i] = 0.0;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    for (l = nylo_fft; l <= nyhi_fft; l++) {
      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        sf_coeff[0] += sf_precoeff1[n]*greensfn[n];
        sf_coeff[1] += sf_precoeff2[n]*greensfn[n];
        sf_coeff[2] += sf_precoeff3[n]*greensfn[n];
        sf_coeff[3] += sf_precoeff4[n]*greensfn[n];
        sf_coeff[4] += sf_precoeff5[n]*greensfn[n];
        sf_coeff[5] += sf_precoeff6[n]*greensfn[n];
        ++n;
      }
    }
  }

  // Compute the coefficients for the self-force correction

  double prex, prey, prez;
  prex = prey = prez = MY_PI/volume;
  prex *= nx_pppm/xprd;
  prey *= ny_pppm/yprd;
  prez *= nz_pppm/zprd_slab;
  sf_coeff[0] *= prex;
  sf_coeff[1] *= prex*2;
  sf_coeff[2] *= prey;
  sf_coeff[3] *= prey*2;
  sf_coeff[4] *= prez;
  sf_coeff[5] *= prez*2;

  // communicate values with other procs

  double tmp[6];
  MPI_Allreduce(sf_coeff,tmp,6,MPI_DOUBLE,MPI_SUM,world);
  for (n = 0; n < 6; n++) sf_coeff[n] = tmp[n];
}

/* ----------------------------------------------------------------------
   compute self force coefficients for ad-differentiation scheme
   and Dispersion interaction 
------------------------------------------------------------------------- */

void PPPMDisp::compute_sf_coeff_6()
{
  int i,k,l,m,n;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  for (i = 0; i <= 5; i++) sf_coeff_6[i] = 0.0;

  n = 0;
  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) {
    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        sf_coeff_6[0] += sf_precoeff1_6[n]*greensfn_6[n];
        sf_coeff_6[1] += sf_precoeff2_6[n]*greensfn_6[n];
        sf_coeff_6[2] += sf_precoeff3_6[n]*greensfn_6[n];
        sf_coeff_6[3] += sf_precoeff4_6[n]*greensfn_6[n];
        sf_coeff_6[4] += sf_precoeff5_6[n]*greensfn_6[n];
        sf_coeff_6[5] += sf_precoeff6_6[n]*greensfn_6[n];
        ++n;
      }
    }
  }

  
  // perform multiplication with prefactors
  
  double prex, prey, prez;
  prex = prey = prez = MY_PI/volume;
  prex *= nx_pppm_6/xprd;
  prey *= ny_pppm_6/yprd;
  prez *= nz_pppm_6/zprd_slab;
  sf_coeff_6[0] *= prex;
  sf_coeff_6[1] *= prex*2;
  sf_coeff_6[2] *= prey;
  sf_coeff_6[3] *= prey*2;
  sf_coeff_6[4] *= prez;
  sf_coeff_6[5] *= prez*2;
  
  // communicate values with other procs
  
  double tmp[6];
  MPI_Allreduce(sf_coeff_6,tmp,6,MPI_DOUBLE,MPI_SUM,world);
  for (n = 0; n < 6; n++) sf_coeff_6[n] = tmp[n];

}

/* ----------------------------------------------------------------------
   denominator for Hockney-Eastwood Green's function
     of x,y,z = sin(kx*deltax/2), etc

            inf                 n-1
   S(n,k) = Sum  W(k+pi*j)**2 = Sum b(l)*(z*z)**l
           j=-inf               l=0

          = -(z*z)**n /(2n-1)! * (d/dx)**(2n-1) cot(x)  at z = sin(x)
   gf_b = denominator expansion coeffs 
------------------------------------------------------------------------- */

double PPPMDisp::gf_denom(double x, double y, double z, double *g_b, int ord)
{
  double sx,sy,sz;
  sz = sy = sx = 0.0;
  for (int l = ord-1; l >= 0; l--) {
    sx = g_b[l] + sx*x;
    sy = g_b[l] + sy*y;
    sz = g_b[l] + sz*z;
  }
  double s = sx*sy*sz;
  return s*s;
}

/* ----------------------------------------------------------------------
   pre-compute Green's function denominator expansion coeffs, Gamma(2n) 
------------------------------------------------------------------------- */

void PPPMDisp::compute_gf_denom(double* gf, int ord)
{
  int k,l,m;
  
  for (l = 1; l < ord; l++) gf[l] = 0.0;
  gf[0] = 1.0;
  
  for (m = 1; m < ord; m++) {
    for (l = m; l > 0; l--) 
      gf[l] = 4.0 * (gf[l]*(l-m)*(l-m-0.5)-gf[l-1]*(l-m-1)*(l-m-1));
    gf[0] = 4.0 * (gf[0]*(l-m)*(l-m-0.5));
  }

  bigint ifact = 1;
  for (k = 1; k < 2*ord; k++) ifact *= k;
  double gaminv = 1.0/ifact;
  for (l = 0; l < ord; l++) gf[l] *= gaminv;
}

/* ----------------------------------------------------------------------
   ghost-swap to accumulate full density in brick decomposition 
   remap density from 3d brick decomposition to FFTdecomposition
   for coulomb interaction or dispersion interaction with geometric
   mixing
------------------------------------------------------------------------- */

void PPPMDisp::brick2fft(int nxlo_i, int nylo_i, int nzlo_i,
                         int nxhi_i, int nyhi_i, int nzhi_i,
                         FFT_SCALAR*** dbrick, FFT_SCALAR* dfft, FFT_SCALAR* work,
                         LAMMPS_NS::Remap* rmp)
{
  int n,ix,iy,iz;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  n = 0;
  for (iz = nzlo_i; iz <= nzhi_i; iz++)
    for (iy = nylo_i; iy <= nyhi_i; iy++)
      for (ix = nxlo_i; ix <= nxhi_i; ix++)
	dfft[n++] = dbrick[iz][iy][ix];

  rmp->perform(dfft,dfft,work);
}


/* ----------------------------------------------------------------------
   ghost-swap to accumulate full density in brick decomposition 
   remap density from 3d brick decomposition to FFTdecomposition
   for dispersion with arithmetic mixing rule
------------------------------------------------------------------------- */

void PPPMDisp::brick2fft_a()
{
  int n,ix,iy,iz;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  n = 0;
  for (iz = nzlo_in_6; iz <= nzhi_in_6; iz++)
    for (iy = nylo_in_6; iy <= nyhi_in_6; iy++)
      for (ix = nxlo_in_6; ix <= nxhi_in_6; ix++) {
        density_fft_a0[n] = density_brick_a0[iz][iy][ix];
        density_fft_a1[n] = density_brick_a1[iz][iy][ix];
        density_fft_a2[n] = density_brick_a2[iz][iy][ix];
        density_fft_a3[n] = density_brick_a3[iz][iy][ix];
        density_fft_a4[n] = density_brick_a4[iz][iy][ix];
        density_fft_a5[n] = density_brick_a5[iz][iy][ix];
        density_fft_a6[n++] = density_brick_a6[iz][iy][ix];
      }

  remap_6->perform(density_fft_a0,density_fft_a0,work1_6);
  remap_6->perform(density_fft_a1,density_fft_a1,work1_6);
  remap_6->perform(density_fft_a2,density_fft_a2,work1_6);
  remap_6->perform(density_fft_a3,density_fft_a3,work1_6);
  remap_6->perform(density_fft_a4,density_fft_a4,work1_6);
  remap_6->perform(density_fft_a5,density_fft_a5,work1_6);
  remap_6->perform(density_fft_a6,density_fft_a6,work1_6);

}

/* ----------------------------------------------------------------------
   ghost-swap to accumulate full density in brick decomposition 
   remap density from 3d brick decomposition to FFTdecomposition
   for dispersion with special case
------------------------------------------------------------------------- */

void PPPMDisp::brick2fft_none()
{
  int k,n,ix,iy,iz;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  for (k = 0; k<nsplit_alloc; k++) {
    n = 0;
    for (iz = nzlo_in_6; iz <= nzhi_in_6; iz++)
      for (iy = nylo_in_6; iy <= nyhi_in_6; iy++)
        for (ix = nxlo_in_6; ix <= nxhi_in_6; ix++) 
          density_fft_none[k][n++] = density_brick_none[k][iz][iy][ix];
  }

  for (k=0; k<nsplit_alloc; k++)
    remap_6->perform(density_fft_none[k],density_fft_none[k],work1_6);
}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array 
------------------------------------------------------------------------- */

void PPPMDisp::particle_map(double delx, double dely, double delz,
                             double sft, int** p2g, int nup, int nlow,
                             int nxlo, int nylo, int nzlo,
                             int nxhi, int nyhi, int nzhi)
{
  int nx,ny,nz;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  if (!isfinite(boxlo[0]) || !isfinite(boxlo[1]) || !isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    
    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    nx = static_cast<int> ((x[i][0]-boxlo[0])*delx+sft) - OFFSET;
    ny = static_cast<int> ((x[i][1]-boxlo[1])*dely+sft) - OFFSET;
    nz = static_cast<int> ((x[i][2]-boxlo[2])*delz+sft) - OFFSET;

    p2g[i][0] = nx;
    p2g[i][1] = ny;
    p2g[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlow < nxlo || nx+nup > nxhi ||
	ny+nlow < nylo || ny+nup > nyhi ||
	nz+nlow < nzlo || nz+nup > nzhi)
      flag = 1;
  }

  if (flag) error->one(FLERR,"Out of range atoms - cannot compute PPPMDisp");
}


void PPPMDisp::particle_map_c(double delx, double dely, double delz,
                               double sft, int** p2g, int nup, int nlow,
                               int nxlo, int nylo, int nzlo,
                               int nxhi, int nyhi, int nzhi)
{
  particle_map(delx, dely, delz, sft, p2g, nup, nlow,
               nxlo, nylo, nzlo, nxhi, nyhi, nzhi);
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid 
------------------------------------------------------------------------- */

void PPPMDisp::make_rho_c()
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density array

  memset(&(density_brick[nzlo_out][nylo_out][nxlo_out]),0,
	 ngrid*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);

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
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = dispersion "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid --- geometric mixing
------------------------------------------------------------------------- */

void PPPMDisp::make_rho_g()
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density array

  memset(&(density_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
	 ngrid_6*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  int type;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {

    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;

    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);
    type = atom->type[i];
    z0 = delvolinv_6 * B[type];
    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      y0 = z0*rho1d_6[2][n];
      for (m = nlower_6; m <= nupper_6; m++) {
	my = m+ny;
	x0 = y0*rho1d_6[1][m];
	for (l = nlower_6; l <= nupper_6; l++) {
	  mx = l+nx;
	  density_brick_g[mz][my][mx] += x0*rho1d_6[0][l];
	}
      }
    }
  }
}


/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = dispersion "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid --- arithmetic mixing
------------------------------------------------------------------------- */

void PPPMDisp::make_rho_a()
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0,w;

  // clear 3d density array

  memset(&(density_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
	 ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
	 ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
	 ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
	 ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
	 ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
	 ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
	 ngrid_6*sizeof(FFT_SCALAR));

  // loop over my particles, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  int type;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {

    //do the following for all 4 grids
    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;
    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);
    type = atom->type[i];
    z0 = delvolinv_6;
    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      y0 = z0*rho1d_6[2][n];
      for (m = nlower_6; m <= nupper_6; m++) {
	my = m+ny;
	x0 = y0*rho1d_6[1][m];
	for (l = nlower_6; l <= nupper_6; l++) {
	  mx = l+nx;
          w = x0*rho1d_6[0][l];
	  density_brick_a0[mz][my][mx] += w*B[7*type];
	  density_brick_a1[mz][my][mx] += w*B[7*type+1];
	  density_brick_a2[mz][my][mx] += w*B[7*type+2];
	  density_brick_a3[mz][my][mx] += w*B[7*type+3];
	  density_brick_a4[mz][my][mx] += w*B[7*type+4];
	  density_brick_a5[mz][my][mx] += w*B[7*type+5];
	  density_brick_a6[mz][my][mx] += w*B[7*type+6];
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = dispersion "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid --- case when mixing rules don't apply
------------------------------------------------------------------------- */

void PPPMDisp::make_rho_none()
{
  int k,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0,w;

  // clear 3d density array
  for (k = 0; k < nsplit_alloc; k++)
    memset(&(density_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
	   ngrid_6*sizeof(FFT_SCALAR));


  // loop over my particles, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  int type;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {

    //do the following for all 4 grids
    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;
    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);
    type = atom->type[i];
    z0 = delvolinv_6;
    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      y0 = z0*rho1d_6[2][n];
      for (m = nlower_6; m <= nupper_6; m++) {
	my = m+ny;
	x0 = y0*rho1d_6[1][m];
	for (l = nlower_6; l <= nupper_6; l++) {
	  mx = l+nx;
          w = x0*rho1d_6[0][l];
          for (k = 0; k < nsplit; k++)
	    density_brick_none[k][mz][my][mx] += w*B[nsplit*type + k];
	}
      }
    }
  }
}


/* ----------------------------------------------------------------------
   FFT-based Poisson solver for ik differentiation
------------------------------------------------------------------------- */

void PPPMDisp::poisson_ik(FFT_SCALAR* wk1, FFT_SCALAR* wk2,
                           FFT_SCALAR* dfft, LAMMPS_NS::FFT3d* ft1,LAMMPS_NS::FFT3d* ft2, 
                           int nx_p, int ny_p, int nz_p, int nft,
                           int nxlo_ft, int nylo_ft, int nzlo_ft,
                           int nxhi_ft, int nyhi_ft, int nzhi_ft,
                           int nxlo_i, int nylo_i, int nzlo_i,
                           int nxhi_i, int nyhi_i, int nzhi_i,
                           double& egy, double* gfn,
                           double* kx, double* ky, double* kz,
                           double* kx2, double* ky2, double* kz2,
                           FFT_SCALAR*** vx_brick, FFT_SCALAR*** vy_brick, FFT_SCALAR*** vz_brick,
                           double* vir, double** vcoeff, double** vcoeff2,
                           FFT_SCALAR*** u_pa, FFT_SCALAR*** v0_pa, FFT_SCALAR*** v1_pa, FFT_SCALAR*** v2_pa,
                           FFT_SCALAR*** v3_pa, FFT_SCALAR*** v4_pa, FFT_SCALAR*** v5_pa)


{
  int i,j,k,n;
  double eng;

  // transform charge/dispersion density (r -> k) 
  n = 0;
  for (i = 0; i < nft; i++) {
    wk1[n++] = dfft[i];
    wk1[n++] = ZEROF;
  }

  ft1->compute(wk1,wk1,1);

  // if requested, compute energy and virial contribution

  double scaleinv = 1.0/(nx_p*ny_p*nz_p);
  double s2 = scaleinv*scaleinv;

  if (eflag_global || vflag_global) {
    if (vflag_global) {
      n = 0;
      for (i = 0; i < nft; i++) {
	eng = s2 * gfn[i] * (wk1[n]*wk1[n] + wk1[n+1]*wk1[n+1]);
	for (j = 0; j < 6; j++) vir[j] += eng*vcoeff[i][j];
	if (eflag_global) egy += eng;
	n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nft; i++) {
	egy += 
	  s2 * gfn[i] * (wk1[n]*wk1[n] + wk1[n+1]*wk1[n+1]);
	n += 2;
      }
    }
  }

  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  n = 0;
  for (i = 0; i < nft; i++) {
    wk1[n++] *= scaleinv * gfn[i];
    wk1[n++] *= scaleinv * gfn[i];
  }

  // compute gradients of V(r) in each of 3 dims by transformimg -ik*V(k)
  // FFT leaves data in 3d brick decomposition
  // copy it into inner portion of vdx,vdy,vdz arrays

  // x & y direction gradient

  n = 0;
  for (k = nzlo_ft; k <= nzhi_ft; k++)
    for (j = nylo_ft; j <= nyhi_ft; j++)
      for (i = nxlo_ft; i <= nxhi_ft; i++) {
	wk2[n] = 0.5*(kx[i]-kx2[i])*wk1[n+1] + 0.5*(ky[j]-ky2[j])*wk1[n];
	wk2[n+1] = -0.5*(kx[i]-kx2[i])*wk1[n] + 0.5*(ky[j]-ky2[j])*wk1[n+1];
	n += 2;
      }

  ft2->compute(wk2,wk2,-1);

  n = 0;
  for (k = nzlo_i; k <= nzhi_i; k++)
    for (j = nylo_i; j <= nyhi_i; j++)
      for (i = nxlo_i; i <= nxhi_i; i++) {
	vx_brick[k][j][i] = wk2[n++];
	vy_brick[k][j][i] = wk2[n++];
      }

  if (!eflag_atom) {
    // z direction gradient only

    n = 0;
    for (k = nzlo_ft; k <= nzhi_ft; k++)
      for (j = nylo_ft; j <= nyhi_ft; j++)
        for (i = nxlo_ft; i <= nxhi_ft; i++) {
	  wk2[n] = kz[k]*wk1[n+1];
	  wk2[n+1] = -kz[k]*wk1[n];
	  n += 2;
        }

    ft2->compute(wk2,wk2,-1);


    n = 0;
    for (k = nzlo_i; k <= nzhi_i; k++)
      for (j = nylo_i; j <= nyhi_i; j++)
        for (i = nxlo_i; i <= nxhi_i; i++) {
	  vz_brick[k][j][i] = wk2[n];
	  n += 2;
        }

  }

  else {
    // z direction gradient & per-atom energy

    n = 0;
    for (k = nzlo_ft; k <= nzhi_ft; k++)
      for (j = nylo_ft; j <= nyhi_ft; j++)
        for (i = nxlo_ft; i <= nxhi_ft; i++) {
	  wk2[n] = 0.5*(kz[k]-kz2[k])*wk1[n+1] - wk1[n+1];
	  wk2[n+1] = -0.5*(kz[k]-kz2[k])*wk1[n] + wk1[n];
	  n += 2;
        }

    ft2->compute(wk2,wk2,-1);

    n = 0;
    for (k = nzlo_i; k <= nzhi_i; k++)
      for (j = nylo_i; j <= nyhi_i; j++)
        for (i = nxlo_i; i <= nxhi_i; i++) {
	  vz_brick[k][j][i] = wk2[n++];
	  u_pa[k][j][i] = wk2[n++];;
        }
  }

  if (vflag_atom) poisson_peratom(wk1, wk2, ft2, vcoeff, vcoeff2, nft,
                                  nxlo_i, nylo_i, nzlo_i, nxhi_i, nyhi_i, nzhi_i,
                                  v0_pa, v1_pa, v2_pa, v3_pa, v4_pa, v5_pa);
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for ad differentiation
------------------------------------------------------------------------- */

void PPPMDisp::poisson_ad(FFT_SCALAR* wk1, FFT_SCALAR* wk2,
                           FFT_SCALAR* dfft, LAMMPS_NS::FFT3d* ft1,LAMMPS_NS::FFT3d* ft2, 
                           int nx_p, int ny_p, int nz_p, int nft,
                           int nxlo_ft, int nylo_ft, int nzlo_ft,
                           int nxhi_ft, int nyhi_ft, int nzhi_ft,
                           int nxlo_i, int nylo_i, int nzlo_i,
                           int nxhi_i, int nyhi_i, int nzhi_i,
                           double& egy, double* gfn,
                           double* vir, double** vcoeff, double** vcoeff2,
                           FFT_SCALAR*** u_pa, FFT_SCALAR*** v0_pa, FFT_SCALAR*** v1_pa, FFT_SCALAR*** v2_pa,
                           FFT_SCALAR*** v3_pa, FFT_SCALAR*** v4_pa, FFT_SCALAR*** v5_pa)


{
  int i,j,k,n;
  double eng;

  // transform charge/dispersion density (r -> k) 
  n = 0;
  for (i = 0; i < nft; i++) {
    wk1[n++] = dfft[i];
    wk1[n++] = ZEROF;
  }

  ft1->compute(wk1,wk1,1);

  // if requested, compute energy and virial contribution

  double scaleinv = 1.0/(nx_p*ny_p*nz_p);
  double s2 = scaleinv*scaleinv;

  if (eflag_global || vflag_global) {
    if (vflag_global) {
      n = 0;
      for (i = 0; i < nft; i++) {
	eng = s2 * gfn[i] * (wk1[n]*wk1[n] + wk1[n+1]*wk1[n+1]);
	for (j = 0; j < 6; j++) vir[j] += eng*vcoeff[i][j];
	if (eflag_global) egy += eng;
	n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nft; i++) {
	egy += 
	  s2 * gfn[i] * (wk1[n]*wk1[n] + wk1[n+1]*wk1[n+1]);
	n += 2;
      }
    }
  }

  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  n = 0;
  for (i = 0; i < nft; i++) {
    wk1[n++] *= scaleinv * gfn[i];
    wk1[n++] *= scaleinv * gfn[i];
  }


  n = 0;
  for (k = nzlo_ft; k <= nzhi_ft; k++)
    for (j = nylo_ft; j <= nyhi_ft; j++)
      for (i = nxlo_ft; i <= nxhi_ft; i++) {
        wk2[n] = wk1[n];
	wk2[n+1] = wk1[n+1];
	n += 2;
      }

  ft2->compute(wk2,wk2,-1);


  n = 0;
  for (k = nzlo_i; k <= nzhi_i; k++)
    for (j = nylo_i; j <= nyhi_i; j++)
      for (i = nxlo_i; i <= nxhi_i; i++) {
	u_pa[k][j][i] = wk2[n++];
        n++;
      }


  if (vflag_atom) poisson_peratom(wk1, wk2, ft2, vcoeff, vcoeff2, nft,
                                  nxlo_i, nylo_i, nzlo_i, nxhi_i, nyhi_i, nzhi_i,
                                  v0_pa, v1_pa, v2_pa, v3_pa, v4_pa, v5_pa);

}

/* ----------------------------------------------------------------------
   Fourier Transform for per atom virial calculations
------------------------------------------------------------------------- */

void PPPMDisp:: poisson_peratom(FFT_SCALAR* wk1, FFT_SCALAR* wk2, LAMMPS_NS::FFT3d* ft2, 
                                 double** vcoeff, double** vcoeff2, int nft,
                                 int nxlo_i, int nylo_i, int nzlo_i,
                                 int nxhi_i, int nyhi_i, int nzhi_i,
                                 FFT_SCALAR*** v0_pa, FFT_SCALAR*** v1_pa, FFT_SCALAR*** v2_pa,
                                 FFT_SCALAR*** v3_pa, FFT_SCALAR*** v4_pa, FFT_SCALAR*** v5_pa)
{
 //v0 & v1 term
  int n, i, j, k;
  n = 0;
  for (i = 0; i < nft; i++) {
    wk2[n] = wk1[n]*vcoeff[i][0] - wk1[n+1]*vcoeff[i][1];
    wk2[n+1] = wk1[n+1]*vcoeff[i][0] +  wk1[n]*vcoeff[i][1];
    n += 2;
  }

  ft2->compute(wk2,wk2,-1); 

  n = 0;
  for (k = nzlo_i; k <= nzhi_i; k++)
    for (j = nylo_i; j <= nyhi_i; j++)
      for (i = nxlo_i; i <= nxhi_i; i++) {
        v0_pa[k][j][i] = wk2[n++];
        v1_pa[k][j][i] = wk2[n++];
      }

  //v2 & v3 term
   
  n = 0;
  for (i = 0; i < nft; i++) {
    wk2[n] = wk1[n]*vcoeff[i][2] - wk1[n+1]*vcoeff2[i][0];
    wk2[n+1] = wk1[n+1]*vcoeff[i][2] + wk1[n]*vcoeff2[i][0];
    n += 2;
  }

  ft2->compute(wk2,wk2,-1); 

  n = 0;
  for (k = nzlo_i; k <= nzhi_i; k++)
    for (j = nylo_i; j <= nyhi_i; j++)
      for (i = nxlo_i; i <= nxhi_i; i++) {
        v2_pa[k][j][i] = wk2[n++];
        v3_pa[k][j][i] = wk2[n++];
      }

  //v4 & v5 term
   
  n = 0;
  for (i = 0; i < nft; i++) {
    wk2[n] = wk1[n]*vcoeff2[i][1] - wk1[n+1]*vcoeff2[i][2];
    wk2[n+1] = wk1[n+1]*vcoeff2[i][1] + wk1[n]*vcoeff2[i][2];
    n += 2;
  }

  ft2->compute(wk2,wk2,-1); 

  n = 0;
  for (k = nzlo_i; k <= nzhi_i; k++)
    for (j = nylo_i; j <= nyhi_i; j++)
      for (i = nxlo_i; i <= nxhi_i; i++) {
        v4_pa[k][j][i] = wk2[n++];
        v5_pa[k][j][i] = wk2[n++];
      }	 
 
}

/* ----------------------------------------------------------------------
   Poisson solver for one mesh with 2 different dispersion densities 
   for ik scheme
------------------------------------------------------------------------- */

void PPPMDisp::poisson_2s_ik(FFT_SCALAR* dfft_1, FFT_SCALAR* dfft_2,
                              FFT_SCALAR*** vxbrick_1, FFT_SCALAR*** vybrick_1, FFT_SCALAR*** vzbrick_1,
                              FFT_SCALAR*** vxbrick_2, FFT_SCALAR*** vybrick_2, FFT_SCALAR*** vzbrick_2,
                              FFT_SCALAR*** u_pa_1, FFT_SCALAR*** v0_pa_1, FFT_SCALAR*** v1_pa_1, FFT_SCALAR*** v2_pa_1,
                              FFT_SCALAR*** v3_pa_1, FFT_SCALAR*** v4_pa_1, FFT_SCALAR*** v5_pa_1,
                              FFT_SCALAR*** u_pa_2, FFT_SCALAR*** v0_pa_2, FFT_SCALAR*** v1_pa_2, FFT_SCALAR*** v2_pa_2,
                              FFT_SCALAR*** v3_pa_2, FFT_SCALAR*** v4_pa_2, FFT_SCALAR*** v5_pa_2)

{
  int i,j,k,n;
  double eng;

  double scaleinv = 1.0/(nx_pppm_6*ny_pppm_6*nz_pppm_6);

  // transform charge/dispersion density (r -> k)
  // only one tansform required when energies and pressures do not
  //  need to be calculated 
  if (eflag_global + vflag_global == 0) {
    n = 0;
    for (i = 0; i < nfft_6; i++) {
      work1_6[n++] = dfft_1[i];
      work1_6[n++] = dfft_2[i];
    }
  
    fft1_6->compute(work1_6,work1_6,1);
  }
  // two transforms are required when energies and pressures are
  //   calculated
  else {
    n = 0;
    for (i = 0; i < nfft_6; i++) {
      work1_6[n] = dfft_1[i];
      work2_6[n++] = ZEROF;
      work1_6[n] = ZEROF;
      work2_6[n++] = dfft_2[i];
    }

    fft1_6->compute(work1_6,work1_6,1);
    fft1_6->compute(work2_6,work2_6,1);

    double s2 = scaleinv*scaleinv;

    if (vflag_global) {
      n = 0;
      for (i = 0; i < nfft_6; i++) {
	eng = 2 * s2 * greensfn_6[i] * (work1_6[n]*work2_6[n+1] - work1_6[n+1]*work2_6[n]);
	for (j = 0; j < 6; j++) virial_6[j] += eng*vg_6[i][j];
	if (eflag_global)energy_6 += eng;
	n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nfft_6; i++) {
	energy_6 += 
	  2 * s2 * greensfn_6[i] * (work1_6[n]*work2_6[n+1] - work1_6[n+1]*work2_6[n]);
	n += 2;
      }
    }
    // unify the two transformed vectors for efficient calculations later
    for ( i = 0; i < 2*nfft_6; i++) {
      work1_6[i] += work2_6[i];
    }
  }

  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work1_6[n++] *= scaleinv * greensfn_6[i];
    work1_6[n++] *= scaleinv * greensfn_6[i];
  }

  // compute gradients of V(r) in each of 3 dims by transformimg -ik*V(k)
  // FFT leaves data in 3d brick decomposition
  // copy it into inner portion of vdx,vdy,vdz arrays

  // x direction gradient

  n = 0;
  for (k = nzlo_fft_6; k <= nzhi_fft_6; k++)
    for (j = nylo_fft_6; j <= nyhi_fft_6; j++)
      for (i = nxlo_fft_6; i <= nxhi_fft_6; i++) {
	work2_6[n] = 0.5*(fkx_6[i]-fkx2_6[i])*work1_6[n+1];
	work2_6[n+1] = -0.5*(fkx_6[i]-fkx2_6[i])*work1_6[n];
	n += 2;
      }

  fft2_6->compute(work2_6,work2_6,-1);
  
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
	vxbrick_1[k][j][i] = work2_6[n++];
        vxbrick_2[k][j][i] = work2_6[n++];
      }

  // y direction gradient

  n = 0;
  for (k = nzlo_fft_6; k <= nzhi_fft_6; k++)
    for (j = nylo_fft_6; j <= nyhi_fft_6; j++)
      for (i = nxlo_fft_6; i <= nxhi_fft_6; i++) {
	work2_6[n] = 0.5*(fky_6[j]-fky2_6[j])*work1_6[n+1];
	work2_6[n+1] = -0.5*(fky_6[j]-fky2_6[j])*work1_6[n];
	n += 2;
      }

  fft2_6->compute(work2_6,work2_6,-1);

  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
	vybrick_1[k][j][i] = work2_6[n++];
        vybrick_2[k][j][i] = work2_6[n++];
      }

  // z direction gradient

  n = 0;
  for (k = nzlo_fft_6; k <= nzhi_fft_6; k++)
    for (j = nylo_fft_6; j <= nyhi_fft_6; j++)
      for (i = nxlo_fft_6; i <= nxhi_fft_6; i++) {
	work2_6[n] = 0.5*(fkz_6[k]-fkz2_6[k])*work1_6[n+1];
	work2_6[n+1] = -0.5*(fkz_6[k]-fkz2_6[k])*work1_6[n];
	n += 2;
      }

  fft2_6->compute(work2_6,work2_6,-1);

  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
	vzbrick_1[k][j][i] = work2_6[n++];
	vzbrick_2[k][j][i] = work2_6[n++];
      }

  //Per-atom energy
    
  if (eflag_atom) {
    n = 0;
    for (i = 0; i < nfft_6; i++) {
      work2_6[n] = work1_6[n];
      work2_6[n+1] = work1_6[n+1];
      n += 2;
    }
    
    fft2_6->compute(work2_6,work2_6,-1); 
    
    n = 0;
    for (k = nzlo_in_6; k <= nzhi_in_6; k++)
      for (j = nylo_in_6; j <= nyhi_in_6; j++)
        for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
          u_pa_1[k][j][i] = work2_6[n++];
          u_pa_2[k][j][i] = work2_6[n++];
        }
  } 

  if (vflag_atom) poisson_2s_peratom(v0_pa_1, v1_pa_1, v2_pa_1, v3_pa_1, v4_pa_1, v5_pa_1,
                                     v0_pa_2, v1_pa_2, v2_pa_2, v3_pa_2, v4_pa_2, v5_pa_2);
}


/* ----------------------------------------------------------------------
   Poisson solver for one mesh with 2 different dispersion densities 
   for ik scheme
------------------------------------------------------------------------- */

void PPPMDisp::poisson_none_ik(int n1, int n2,FFT_SCALAR* dfft_1, FFT_SCALAR* dfft_2,
                              FFT_SCALAR*** vxbrick_1, FFT_SCALAR*** vybrick_1, FFT_SCALAR*** vzbrick_1,
                              FFT_SCALAR*** vxbrick_2, FFT_SCALAR*** vybrick_2, FFT_SCALAR*** vzbrick_2,
                              FFT_SCALAR**** u_pa, FFT_SCALAR**** v0_pa, FFT_SCALAR**** v1_pa, FFT_SCALAR**** v2_pa,
                              FFT_SCALAR**** v3_pa, FFT_SCALAR**** v4_pa, FFT_SCALAR**** v5_pa)
{
  int i,j,k,n;
  double eng;

  double scaleinv = 1.0/(nx_pppm_6*ny_pppm_6*nz_pppm_6);

  // transform charge/dispersion density (r -> k)
  // only one tansform required when energies and pressures do not
  //  need to be calculated 
  if (eflag_global + vflag_global == 0) {
    n = 0;
    for (i = 0; i < nfft_6; i++) {
      work1_6[n++] = dfft_1[i];
      work1_6[n++] = dfft_2[i];
    }
  
    fft1_6->compute(work1_6,work1_6,1);
  }


  // two transforms are required when energies and pressures are
  //   calculated
  else {
    n = 0;
    for (i = 0; i < nfft_6; i++) {
      work1_6[n] = dfft_1[i];
      work2_6[n++] = ZEROF;
      work1_6[n] = ZEROF;
      work2_6[n++] = dfft_2[i];
    }
   

    fft1_6->compute(work1_6,work1_6,1);
    fft1_6->compute(work2_6,work2_6,1);

    double s2 = scaleinv*scaleinv;

    if (vflag_global) {
      n = 0;
      for (i = 0; i < nfft_6; i++) {
	eng = s2 * greensfn_6[i] * (B[n1]*(work1_6[n]*work1_6[n] + work1_6[n+1]*work1_6[n+1]) + B[n2]*(work2_6[n]*work2_6[n] + work2_6[n+1]*work2_6[n+1]));
	for (j = 0; j < 6; j++) virial_6[j] += eng*vg_6[i][j];
	if (eflag_global)energy_6 += eng;
	n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nfft_6; i++) {
	energy_6 += 
	  s2 * greensfn_6[i] * (B[n1]*(work1_6[n]*work1_6[n] + work1_6[n+1]*work1_6[n+1]) + B[n2]*(work2_6[n]*work2_6[n] + work2_6[n+1]*work2_6[n+1]));
	n += 2;
      }
    }
    // unify the two transformed vectors for efficient calculations later
    for ( i = 0; i < 2*nfft_6; i++) {
      work1_6[i] += work2_6[i];
    }
  }

  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work1_6[n++] *= scaleinv * greensfn_6[i];
    work1_6[n++] *= scaleinv * greensfn_6[i];
  }

  // compute gradients of V(r) in each of 3 dims by transformimg -ik*V(k)
  // FFT leaves data in 3d brick decomposition
  // copy it into inner portion of vdx,vdy,vdz arrays

  // x direction gradient

  n = 0;
  for (k = nzlo_fft_6; k <= nzhi_fft_6; k++)
    for (j = nylo_fft_6; j <= nyhi_fft_6; j++)
      for (i = nxlo_fft_6; i <= nxhi_fft_6; i++) {
	work2_6[n] = 0.5*(fkx_6[i]-fkx2_6[i])*work1_6[n+1];
	work2_6[n+1] = -0.5*(fkx_6[i]-fkx2_6[i])*work1_6[n];
	n += 2;
      }

  fft2_6->compute(work2_6,work2_6,-1);
  
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
	vxbrick_1[k][j][i] = B[n1]*work2_6[n++];
        vxbrick_2[k][j][i] = B[n2]*work2_6[n++];
      }

  // y direction gradient

  n = 0;
  for (k = nzlo_fft_6; k <= nzhi_fft_6; k++)
    for (j = nylo_fft_6; j <= nyhi_fft_6; j++)
      for (i = nxlo_fft_6; i <= nxhi_fft_6; i++) {
	work2_6[n] = 0.5*(fky_6[j]-fky2_6[j])*work1_6[n+1];
	work2_6[n+1] = -0.5*(fky_6[j]-fky2_6[j])*work1_6[n];
	n += 2;
      }

  fft2_6->compute(work2_6,work2_6,-1);

  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
	vybrick_1[k][j][i] = B[n1]*work2_6[n++];
        vybrick_2[k][j][i] = B[n2]*work2_6[n++];
      }

  // z direction gradient

  n = 0;
  for (k = nzlo_fft_6; k <= nzhi_fft_6; k++)
    for (j = nylo_fft_6; j <= nyhi_fft_6; j++)
      for (i = nxlo_fft_6; i <= nxhi_fft_6; i++) {
	work2_6[n] = 0.5*(fkz_6[k]-fkz2_6[k])*work1_6[n+1];
	work2_6[n+1] = -0.5*(fkz_6[k]-fkz2_6[k])*work1_6[n];
	n += 2;
      }

  fft2_6->compute(work2_6,work2_6,-1);

  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
	vzbrick_1[k][j][i] = B[n1]*work2_6[n++];
	vzbrick_2[k][j][i] = B[n2]*work2_6[n++];
      }

  //Per-atom energy
    
  if (eflag_atom) {
    n = 0;
    for (i = 0; i < nfft_6; i++) {
      work2_6[n] = work1_6[n];
      work2_6[n+1] = work1_6[n+1];
      n += 2;
    }
    
    fft2_6->compute(work2_6,work2_6,-1); 
    
    n = 0;
    for (k = nzlo_in_6; k <= nzhi_in_6; k++)
      for (j = nylo_in_6; j <= nyhi_in_6; j++)
        for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
          u_pa[n1][k][j][i] = B[n1]*work2_6[n++];
          u_pa[n2][k][j][i] = B[n2]*work2_6[n++];
        }
  } 

  if (vflag_atom) poisson_none_peratom(n1,n2,
                                       v0_pa[n1], v1_pa[n1], v2_pa[n1], v3_pa[n1], v4_pa[n1], v5_pa[n1],
                                       v0_pa[n2], v1_pa[n2], v2_pa[n2], v3_pa[n2], v4_pa[n2], v5_pa[n2]);
}

/* ----------------------------------------------------------------------
   Poisson solver for one mesh with 2 different dispersion densities 
   for ad scheme
------------------------------------------------------------------------- */

void PPPMDisp::poisson_2s_ad(FFT_SCALAR* dfft_1, FFT_SCALAR* dfft_2,
                              FFT_SCALAR*** u_pa_1, FFT_SCALAR*** v0_pa_1, FFT_SCALAR*** v1_pa_1, FFT_SCALAR*** v2_pa_1,
                              FFT_SCALAR*** v3_pa_1, FFT_SCALAR*** v4_pa_1, FFT_SCALAR*** v5_pa_1,
                              FFT_SCALAR*** u_pa_2, FFT_SCALAR*** v0_pa_2, FFT_SCALAR*** v1_pa_2, FFT_SCALAR*** v2_pa_2,
                              FFT_SCALAR*** v3_pa_2, FFT_SCALAR*** v4_pa_2, FFT_SCALAR*** v5_pa_2)

{
  int i,j,k,n;
  double eng;

  double scaleinv = 1.0/(nx_pppm_6*ny_pppm_6*nz_pppm_6);

  // transform charge/dispersion density (r -> k)
  // only one tansform required when energies and pressures do not
  //  need to be calculated 
  if (eflag_global + vflag_global == 0) {
    n = 0;
    for (i = 0; i < nfft_6; i++) {
      work1_6[n++] = dfft_1[i];
      work1_6[n++] = dfft_2[i];
    }
  
    fft1_6->compute(work1_6,work1_6,1);
  }
  // two transforms are required when energies and pressures are
  //   calculated
  else {
    n = 0;
    for (i = 0; i < nfft_6; i++) {
      work1_6[n] = dfft_1[i];
      work2_6[n++] = ZEROF;
      work1_6[n] = ZEROF;
      work2_6[n++] = dfft_2[i];
    }

    fft1_6->compute(work1_6,work1_6,1);
    fft1_6->compute(work2_6,work2_6,1);

    double s2 = scaleinv*scaleinv;

    if (vflag_global) {
      n = 0;
      for (i = 0; i < nfft_6; i++) {
	eng = 2 * s2 * greensfn_6[i] * (work1_6[n]*work2_6[n+1] - work1_6[n+1]*work2_6[n]);
	for (j = 0; j < 6; j++) virial_6[j] += eng*vg_6[i][j];
	if (eflag_global)energy_6 += eng;
	n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nfft_6; i++) {
	energy_6 += 
	  2 * s2 * greensfn_6[i] * (work1_6[n]*work2_6[n+1] - work1_6[n+1]*work2_6[n]);
	n += 2;
      }
    }
    // unify the two transformed vectors for efficient calculations later
    for ( i = 0; i < 2*nfft_6; i++) {
      work1_6[i] += work2_6[i];
    }
  }


  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work1_6[n++] *= scaleinv * greensfn_6[i];
    work1_6[n++] *= scaleinv * greensfn_6[i];
  }


  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n];
    work2_6[n+1] = work1_6[n+1];
    n += 2;
  }
    
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        u_pa_1[k][j][i] = work2_6[n++];
        u_pa_2[k][j][i] = work2_6[n++];
      } 

  if (vflag_atom) poisson_2s_peratom(v0_pa_1, v1_pa_1, v2_pa_1, v3_pa_1, v4_pa_1, v5_pa_1,
                                     v0_pa_2, v1_pa_2, v2_pa_2, v3_pa_2, v4_pa_2, v5_pa_2);
}

/* ----------------------------------------------------------------------
   Poisson solver for one mesh with 2 different dispersion densities 
   for ad scheme
------------------------------------------------------------------------- */

void PPPMDisp::poisson_none_ad(int n1, int n2, FFT_SCALAR* dfft_1, FFT_SCALAR* dfft_2,
                               FFT_SCALAR*** u_pa_1, FFT_SCALAR*** u_pa_2,
                               FFT_SCALAR**** v0_pa, FFT_SCALAR**** v1_pa, FFT_SCALAR**** v2_pa,
                               FFT_SCALAR**** v3_pa, FFT_SCALAR**** v4_pa, FFT_SCALAR**** v5_pa)
{
  int i,j,k,n;
  double eng;

  double scaleinv = 1.0/(nx_pppm_6*ny_pppm_6*nz_pppm_6);

  // transform charge/dispersion density (r -> k)
  // only one tansform required when energies and pressures do not
  //  need to be calculated 
  if (eflag_global + vflag_global == 0) {
    n = 0;
    for (i = 0; i < nfft_6; i++) {
      work1_6[n++] = dfft_1[i];
      work1_6[n++] = dfft_2[i];
    }
  
    fft1_6->compute(work1_6,work1_6,1);
  }
  // two transforms are required when energies and pressures are
  //   calculated
  else {
    n = 0;
    for (i = 0; i < nfft_6; i++) {
      work1_6[n] = dfft_1[i];
      work2_6[n++] = ZEROF;
      work1_6[n] = ZEROF;
      work2_6[n++] = dfft_2[i];
    }

    fft1_6->compute(work1_6,work1_6,1);
    fft1_6->compute(work2_6,work2_6,1);

    double s2 = scaleinv*scaleinv;

    if (vflag_global) {
      n = 0;
      for (i = 0; i < nfft_6; i++) {
	eng = s2 * greensfn_6[i] * (B[n1]*(work1_6[n]*work1_6[n] + work1_6[n+1]*work1_6[n+1]) + B[n2]*(work2_6[n]*work2_6[n] + work2_6[n+1]*work2_6[n+1]));
	for (j = 0; j < 6; j++) virial_6[j] += eng*vg_6[i][j];
	if (eflag_global)energy_6 += eng;
	n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nfft_6; i++) {
	energy_6 += 
	  s2 * greensfn_6[i] * (B[n1]*(work1_6[n]*work1_6[n] + work1_6[n+1]*work1_6[n+1]) + B[n2]*(work2_6[n]*work2_6[n] + work2_6[n+1]*work2_6[n+1]));
	n += 2;
      }
    }
    // unify the two transformed vectors for efficient calculations later
    for ( i = 0; i < 2*nfft_6; i++) {
      work1_6[i] += work2_6[i];
    }
  }


  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work1_6[n++] *= scaleinv * greensfn_6[i];
    work1_6[n++] *= scaleinv * greensfn_6[i];
  }


  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n];
    work2_6[n+1] = work1_6[n+1];
    n += 2;
  }
  
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        u_pa_1[k][j][i] = B[n1]*work2_6[n++];
        u_pa_2[k][j][i] = B[n2]*work2_6[n++];
      } 

  if (vflag_atom) poisson_none_peratom(n1,n2,
                                       v0_pa[n1], v1_pa[n1], v2_pa[n1], v3_pa[n1], v4_pa[n1], v5_pa[n1],
                                       v0_pa[n2], v1_pa[n2], v2_pa[n2], v3_pa[n2], v4_pa[n2], v5_pa[n2]);
}

/* ----------------------------------------------------------------------
   Fourier Transform for per atom virial calculations
------------------------------------------------------------------------- */

void PPPMDisp::poisson_2s_peratom(FFT_SCALAR*** v0_pa_1, FFT_SCALAR*** v1_pa_1, FFT_SCALAR*** v2_pa_1,
                                   FFT_SCALAR*** v3_pa_1, FFT_SCALAR*** v4_pa_1, FFT_SCALAR*** v5_pa_1,
                                   FFT_SCALAR*** v0_pa_2, FFT_SCALAR*** v1_pa_2, FFT_SCALAR*** v2_pa_2,
                                   FFT_SCALAR*** v3_pa_2, FFT_SCALAR*** v4_pa_2, FFT_SCALAR*** v5_pa_2)
{
  //Compute first virial term v0
  int n, i, j, k;

  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg_6[i][0];
    work2_6[n+1] = work1_6[n+1]*vg_6[i][0];
    n += 2;
  }
   
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v0_pa_1[k][j][i] = work2_6[n++];
        v0_pa_2[k][j][i] = work2_6[n++];
      }
	 
  //Compute second virial term v1  
  
  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg_6[i][1];
    work2_6[n+1] = work1_6[n+1]*vg_6[i][1];
    n += 2;
  }
    
  fft2_6->compute(work2_6,work2_6,-1); 
  
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v1_pa_1[k][j][i] = work2_6[n++];
        v1_pa_2[k][j][i] = work2_6[n++];
      }
	  
  //Compute third virial term v2
   
  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg_6[i][2];
    work2_6[n+1] = work1_6[n+1]*vg_6[i][2];
    n += 2;
  }
    
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v2_pa_1[k][j][i] = work2_6[n++];
        v2_pa_2[k][j][i] = work2_6[n++];
      }

  //Compute fourth virial term v3
   
  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg2_6[i][0];
    work2_6[n+1] = work1_6[n+1]*vg2_6[i][0];
    n += 2;
  }
    
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v3_pa_1[k][j][i] = work2_6[n++];
        v3_pa_2[k][j][i] = work2_6[n++];
      }

  //Compute fifth virial term v4
   
  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg2_6[i][1];
    work2_6[n+1] = work1_6[n+1]*vg2_6[i][1];
    n += 2;
  }
    
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v4_pa_1[k][j][i] = work2_6[n++];
        v4_pa_2[k][j][i] = work2_6[n++];
      }
   
  //Compute last virial term v5
   
  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg2_6[i][2];
    work2_6[n+1] = work1_6[n+1]*vg2_6[i][2];
    n += 2;
  }
    
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v5_pa_1[k][j][i] = work2_6[n++];
        v5_pa_2[k][j][i] = work2_6[n++];
      }
}

/* ----------------------------------------------------------------------
   Fourier Transform for per atom virial calculations
------------------------------------------------------------------------- */

void PPPMDisp::poisson_none_peratom(int n1, int n2,                              
                                 FFT_SCALAR*** v0_pa_1, FFT_SCALAR*** v1_pa_1, FFT_SCALAR*** v2_pa_1,
                                 FFT_SCALAR*** v3_pa_1, FFT_SCALAR*** v4_pa_1, FFT_SCALAR*** v5_pa_1,
                                 FFT_SCALAR*** v0_pa_2, FFT_SCALAR*** v1_pa_2, FFT_SCALAR*** v2_pa_2,
                                 FFT_SCALAR*** v3_pa_2, FFT_SCALAR*** v4_pa_2, FFT_SCALAR*** v5_pa_2)
{
  //Compute first virial term v0
  int n, i, j, k;

  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg_6[i][0];
    work2_6[n+1] = work1_6[n+1]*vg_6[i][0];
    n += 2;
  }
   
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v0_pa_1[k][j][i] = B[n1]*work2_6[n++];
        v0_pa_2[k][j][i] = B[n2]*work2_6[n++];
      }
	 
  //Compute second virial term v1  
  
  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg_6[i][1];
    work2_6[n+1] = work1_6[n+1]*vg_6[i][1];
    n += 2;
  }
    
  fft2_6->compute(work2_6,work2_6,-1); 
  
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v1_pa_1[k][j][i] = B[n1]*work2_6[n++];
        v1_pa_2[k][j][i] = B[n2]*work2_6[n++];
      }
	  
  //Compute third virial term v2
   
  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg_6[i][2];
    work2_6[n+1] = work1_6[n+1]*vg_6[i][2];
    n += 2;
  }
    
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v2_pa_1[k][j][i] = B[n1]*work2_6[n++];
        v2_pa_2[k][j][i] = B[n2]*work2_6[n++];
      }

  //Compute fourth virial term v3
   
  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg2_6[i][0];
    work2_6[n+1] = work1_6[n+1]*vg2_6[i][0];
    n += 2;
  }
    
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v3_pa_1[k][j][i] = B[n1]*work2_6[n++];
        v3_pa_2[k][j][i] = B[n2]*work2_6[n++];
      }

  //Compute fifth virial term v4
   
  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg2_6[i][1];
    work2_6[n+1] = work1_6[n+1]*vg2_6[i][1];
    n += 2;
  }
    
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v4_pa_1[k][j][i] = B[n1]*work2_6[n++];
        v4_pa_2[k][j][i] = B[n2]*work2_6[n++];
      }
   
  //Compute last virial term v5
   
  n = 0;
  for (i = 0; i < nfft_6; i++) {
    work2_6[n] = work1_6[n]*vg2_6[i][2];
    work2_6[n+1] = work1_6[n+1]*vg2_6[i][2];
    n += 2;
  }
    
  fft2_6->compute(work2_6,work2_6,-1); 
    
  n = 0;
  for (k = nzlo_in_6; k <= nzhi_in_6; k++)
    for (j = nylo_in_6; j <= nyhi_in_6; j++)
      for (i = nxlo_in_6; i <= nxhi_in_6; i++) {
        v5_pa_1[k][j][i] = B[n1]*work2_6[n++];
        v5_pa_2[k][j][i] = B[n2]*work2_6[n++];
      }
}
 
/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles 
   for ik scheme
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_c_ik()
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

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);

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

    const double qfactor = force->qqrd2e * scale * q[i];
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    if (slabflag != 2) f[i][2] += qfactor*ekz;
  }
}
/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
   for ad scheme 
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_c_ad()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz;
  FFT_SCALAR ekx,eky,ekz;
  double s1,s2,s3;
  double sf = 0.0;

  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double hx_inv = nx_pppm/xprd;
  double hy_inv = ny_pppm/yprd;
  double hz_inv = nz_pppm/zprd_slab;

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
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);
    compute_drho1d(dx,dy,dz, order, drho_coeff, drho1d);

    ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          ekx += drho1d[0][l]*rho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          eky += rho1d[0][l]*drho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          ekz += rho1d[0][l]*rho1d[1][m]*drho1d[2][n]*u_brick[mz][my][mx];
        }
      }
    }
    ekx *= hx_inv;
    eky *= hy_inv;
    ekz *= hz_inv;
    // convert E-field to force and substract self forces
    const double qfactor = force->qqrd2e * scale;

    s1 = x[i][0]*hx_inv;
    s2 = x[i][1]*hy_inv;
    s3 = x[i][2]*hz_inv;
    sf = sf_coeff[0]*sin(2*MY_PI*s1);
    sf += sf_coeff[1]*sin(4*MY_PI*s1);
    sf *= 2*q[i]*q[i];
    f[i][0] += qfactor*(ekx*q[i] - sf);

    sf = sf_coeff[2]*sin(2*MY_PI*s2);
    sf += sf_coeff[3]*sin(4*MY_PI*s2);
    sf *= 2*q[i]*q[i];
    f[i][1] += qfactor*(eky*q[i] - sf);


    sf = sf_coeff[4]*sin(2*MY_PI*s3);
    sf += sf_coeff[5]*sin(4*MY_PI*s3);
    sf *= 2*q[i]*q[i];
    if (slabflag != 2) f[i][2] += qfactor*(ekz*q[i] - sf);
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles 
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_c_peratom()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR u_pa,v0,v1,v2,v3,v4,v5;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);

    u_pa = v0 = v1 = v2 = v3 = v4 = v5 = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
	my = m+ny;
	y0 = z0*rho1d[1][m];
	for (l = nlower; l <= nupper; l++) {
	  mx = l+nx;
	  x0 = y0*rho1d[0][l];
	  if (eflag_atom) u_pa += x0*u_brick[mz][my][mx];	
	  if (vflag_atom) {
            v0 += x0*v0_brick[mz][my][mx];
            v1 += x0*v1_brick[mz][my][mx];
            v2 += x0*v2_brick[mz][my][mx];
            v3 += x0*v3_brick[mz][my][mx];
            v4 += x0*v4_brick[mz][my][mx];
            v5 += x0*v5_brick[mz][my][mx];
          }
	}
      }
    }

    // convert E-field to force

    const double qfactor = 0.5*force->qqrd2e * scale * q[i];

    if (eflag_atom) eatom[i] += u_pa*qfactor;
    if (vflag_atom) {
      vatom[i][0] += v0*qfactor;
      vatom[i][1] += v1*qfactor;
      vatom[i][2] += v2*qfactor;
      vatom[i][3] += v3*qfactor;
      vatom[i][4] += v4*qfactor;
      vatom[i][5] += v5*qfactor;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for geometric mixing rule 
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_g_ik()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  double **x = atom->x;
  double **f = atom->f;
  int type;
  double lj;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;

    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);

    ekx = eky = ekz = ZEROF;
    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      z0 = rho1d_6[2][n];
      for (m = nlower_6; m <= nupper_6; m++) {
	my = m+ny;
	y0 = z0*rho1d_6[1][m];
	for (l = nlower_6; l <= nupper_6; l++) {
	  mx = l+nx;
	  x0 = y0*rho1d_6[0][l];
	  ekx -= x0*vdx_brick_g[mz][my][mx];
	  eky -= x0*vdy_brick_g[mz][my][mx];
	  ekz -= x0*vdz_brick_g[mz][my][mx];
	}
      }
    }

    // convert E-field to force
    type = atom->type[i];
    lj = B[type];
    f[i][0] += lj*ekx;
    f[i][1] += lj*eky;
    if (slabflag != 2) f[i][2] += lj*ekz;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for geometric mixing rule for ad scheme
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_g_ad()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz;
  FFT_SCALAR ekx,eky,ekz;
  double s1,s2,s3;
  double sf = 0.0;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double hx_inv = nx_pppm_6/xprd;
  double hy_inv = ny_pppm_6/yprd;
  double hz_inv = nz_pppm_6/zprd_slab;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  double **x = atom->x;
  double **f = atom->f;
  int type;
  double lj;

  int nlocal = atom->nlocal;

 
  for (i = 0; i < nlocal; i++) {
    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;

    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);
    compute_drho1d(dx,dy,dz, order_6, drho_coeff_6, drho1d_6);


    ekx = eky = ekz = ZEROF;
    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      for (m = nlower_6; m <= nupper_6; m++) {
        my = m+ny;
        for (l = nlower_6; l <= nupper_6; l++) {
          mx = l+nx;
          ekx += drho1d_6[0][l]*rho1d_6[1][m]*rho1d_6[2][n]*u_brick_g[mz][my][mx];
          eky += rho1d_6[0][l]*drho1d_6[1][m]*rho1d_6[2][n]*u_brick_g[mz][my][mx];
          ekz += rho1d_6[0][l]*rho1d_6[1][m]*drho1d_6[2][n]*u_brick_g[mz][my][mx];
        }
      }
    }
    ekx *= hx_inv;
    eky *= hy_inv;
    ekz *= hz_inv;

    // convert E-field to force
    type = atom->type[i];
    lj = B[type];

    s1 = x[i][0]*hx_inv;
    s2 = x[i][1]*hy_inv;
    s3 = x[i][2]*hz_inv;

    sf = sf_coeff_6[0]*sin(2*MY_PI*s1);
    sf += sf_coeff_6[1]*sin(4*MY_PI*s1);
    sf *= 2*lj*lj;
    f[i][0] += ekx*lj - sf;

    sf = sf_coeff_6[2]*sin(2*MY_PI*s2);
    sf += sf_coeff_6[3]*sin(4*MY_PI*s2);
    sf *= 2*lj*lj;
    f[i][1] += eky*lj - sf;


    sf = sf_coeff_6[4]*sin(2*MY_PI*s3);
    sf += sf_coeff_6[5]*sin(4*MY_PI*s3);
    sf *= 2*lj*lj;
    if (slabflag != 2) f[i][2] += ekz*lj - sf;

  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for geometric mixing rule for per atom quantities
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_g_peratom()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR u_pa,v0,v1,v2,v3,v4,v5;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  double **x = atom->x;
  int type;
  double lj;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;

    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);

    u_pa = v0 = v1 = v2 = v3 = v4 = v5 = ZEROF;
    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      z0 = rho1d_6[2][n];
      for (m = nlower_6; m <= nupper_6; m++) {
	my = m+ny;
	y0 = z0*rho1d_6[1][m];
	for (l = nlower_6; l <= nupper_6; l++) {
	  mx = l+nx;
	  x0 = y0*rho1d_6[0][l];
	  if (eflag_atom) u_pa += x0*u_brick_g[mz][my][mx];	
	  if (vflag_atom) {
            v0 += x0*v0_brick_g[mz][my][mx];
            v1 += x0*v1_brick_g[mz][my][mx];
            v2 += x0*v2_brick_g[mz][my][mx];
            v3 += x0*v3_brick_g[mz][my][mx];
            v4 += x0*v4_brick_g[mz][my][mx];
            v5 += x0*v5_brick_g[mz][my][mx];
          }
	}
      }
    }

    // convert E-field to force
    type = atom->type[i];
    lj = B[type]*0.5;

    if (eflag_atom) eatom[i] += u_pa*lj;
    if (vflag_atom) {
      vatom[i][0] += v0*lj;
      vatom[i][1] += v1*lj;
      vatom[i][2] += v2*lj;
      vatom[i][3] += v3*lj;
      vatom[i][4] += v4*lj;
      vatom[i][5] += v5*lj;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for arithmetic mixing rule and ik scheme
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_a_ik()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx0, eky0, ekz0, ekx1, eky1, ekz1, ekx2, eky2, ekz2;
  FFT_SCALAR ekx3, eky3, ekz3, ekx4, eky4, ekz4, ekx5, eky5, ekz5;
  FFT_SCALAR ekx6, eky6, ekz6;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  double **x = atom->x;
  double **f = atom->f;
  int type;
  double lj0, lj1, lj2, lj3, lj4, lj5, lj6;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {

    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;
    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);
    ekx0 = eky0 = ekz0 = ZEROF;
    ekx1 = eky1 = ekz1 = ZEROF;
    ekx2 = eky2 = ekz2 = ZEROF;
    ekx3 = eky3 = ekz3 = ZEROF;
    ekx4 = eky4 = ekz4 = ZEROF;
    ekx5 = eky5 = ekz5 = ZEROF;
    ekx6 = eky6 = ekz6 = ZEROF;
    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      z0 = rho1d_6[2][n];
      for (m = nlower_6; m <= nupper_6; m++) {
	my = m+ny;
	y0 = z0*rho1d_6[1][m];
	for (l = nlower_6; l <= nupper_6; l++) {
	  mx = l+nx;
	  x0 = y0*rho1d_6[0][l];
	  ekx0 -= x0*vdx_brick_a0[mz][my][mx];
	  eky0 -= x0*vdy_brick_a0[mz][my][mx];
	  ekz0 -= x0*vdz_brick_a0[mz][my][mx];
	  ekx1 -= x0*vdx_brick_a1[mz][my][mx];
	  eky1 -= x0*vdy_brick_a1[mz][my][mx];
	  ekz1 -= x0*vdz_brick_a1[mz][my][mx];
          ekx2 -= x0*vdx_brick_a2[mz][my][mx];
	  eky2 -= x0*vdy_brick_a2[mz][my][mx];
	  ekz2 -= x0*vdz_brick_a2[mz][my][mx];
	  ekx3 -= x0*vdx_brick_a3[mz][my][mx];
	  eky3 -= x0*vdy_brick_a3[mz][my][mx];
	  ekz3 -= x0*vdz_brick_a3[mz][my][mx];
	  ekx4 -= x0*vdx_brick_a4[mz][my][mx];
	  eky4 -= x0*vdy_brick_a4[mz][my][mx];
	  ekz4 -= x0*vdz_brick_a4[mz][my][mx];
          ekx5 -= x0*vdx_brick_a5[mz][my][mx];
	  eky5 -= x0*vdy_brick_a5[mz][my][mx];
	  ekz5 -= x0*vdz_brick_a5[mz][my][mx];
          ekx6 -= x0*vdx_brick_a6[mz][my][mx];
	  eky6 -= x0*vdy_brick_a6[mz][my][mx];
	  ekz6 -= x0*vdz_brick_a6[mz][my][mx];
	}
      }
    }
    // convert D-field to force
    type = atom->type[i];
    lj0 = B[7*type+6];
    lj1 = B[7*type+5];
    lj2 = B[7*type+4];
    lj3 = B[7*type+3];
    lj4 = B[7*type+2];
    lj5 = B[7*type+1];
    lj6 = B[7*type];
    f[i][0] += lj0*ekx0 + lj1*ekx1 + lj2*ekx2 + lj3*ekx3 + lj4*ekx4 + lj5*ekx5 + lj6*ekx6;
    f[i][1] += lj0*eky0 + lj1*eky1 + lj2*eky2 + lj3*eky3 + lj4*eky4 + lj5*eky5 + lj6*eky6;
    if (slabflag != 2) f[i][2] += lj0*ekz0 + lj1*ekz1 + lj2*ekz2 + lj3*ekz3 + lj4*ekz4 + lj5*ekz5 + lj6*ekz6;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for arithmetic mixing rule for the ad scheme
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_a_ad()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx0, eky0, ekz0, ekx1, eky1, ekz1, ekx2, eky2, ekz2;
  FFT_SCALAR ekx3, eky3, ekz3, ekx4, eky4, ekz4, ekx5, eky5, ekz5;
  FFT_SCALAR ekx6, eky6, ekz6;

  double s1,s2,s3;
  double sf = 0.0;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double hx_inv = nx_pppm_6/xprd;
  double hy_inv = ny_pppm_6/yprd;
  double hz_inv = nz_pppm_6/zprd_slab;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  double **x = atom->x;
  double **f = atom->f;
  int type;
  double lj0, lj1, lj2, lj3, lj4, lj5, lj6;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {

    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;

    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);
    compute_drho1d(dx,dy,dz, order_6, drho_coeff_6, drho1d_6);

    ekx0 = eky0 = ekz0 = ZEROF;
    ekx1 = eky1 = ekz1 = ZEROF;
    ekx2 = eky2 = ekz2 = ZEROF;
    ekx3 = eky3 = ekz3 = ZEROF;
    ekx4 = eky4 = ekz4 = ZEROF;
    ekx5 = eky5 = ekz5 = ZEROF;
    ekx6 = eky6 = ekz6 = ZEROF;
    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      for (m = nlower_6; m <= nupper_6; m++) {
	my = m+ny;
	for (l = nlower_6; l <= nupper_6; l++) {
	  mx = l+nx;
          x0 = drho1d_6[0][l]*rho1d_6[1][m]*rho1d_6[2][n];
          y0 = rho1d_6[0][l]*drho1d_6[1][m]*rho1d_6[2][n];
          z0 = rho1d_6[0][l]*rho1d_6[1][m]*drho1d_6[2][n];

          ekx0 += x0*u_brick_a0[mz][my][mx];
          eky0 += y0*u_brick_a0[mz][my][mx];
          ekz0 += z0*u_brick_a0[mz][my][mx];

          ekx1 += x0*u_brick_a1[mz][my][mx];
          eky1 += y0*u_brick_a1[mz][my][mx];
          ekz1 += z0*u_brick_a1[mz][my][mx];

          ekx2 += x0*u_brick_a2[mz][my][mx];
          eky2 += y0*u_brick_a2[mz][my][mx];
          ekz2 += z0*u_brick_a2[mz][my][mx];

          ekx3 += x0*u_brick_a3[mz][my][mx];
          eky3 += y0*u_brick_a3[mz][my][mx];
          ekz3 += z0*u_brick_a3[mz][my][mx];

          ekx4 += x0*u_brick_a4[mz][my][mx];
          eky4 += y0*u_brick_a4[mz][my][mx];
          ekz4 += z0*u_brick_a4[mz][my][mx];

          ekx5 += x0*u_brick_a5[mz][my][mx];
          eky5 += y0*u_brick_a5[mz][my][mx];
          ekz5 += z0*u_brick_a5[mz][my][mx];

          ekx6 += x0*u_brick_a6[mz][my][mx];
          eky6 += y0*u_brick_a6[mz][my][mx];
          ekz6 += z0*u_brick_a6[mz][my][mx];
	}
      }
    }

    ekx0 *= hx_inv;
    eky0 *= hy_inv;
    ekz0 *= hz_inv;

    ekx1 *= hx_inv;
    eky1 *= hy_inv;
    ekz1 *= hz_inv;

    ekx2 *= hx_inv;
    eky2 *= hy_inv;
    ekz2 *= hz_inv;

    ekx3 *= hx_inv;
    eky3 *= hy_inv;
    ekz3 *= hz_inv;

    ekx4 *= hx_inv;
    eky4 *= hy_inv;
    ekz4 *= hz_inv;

    ekx5 *= hx_inv;
    eky5 *= hy_inv;
    ekz5 *= hz_inv;

    ekx6 *= hx_inv;
    eky6 *= hy_inv;
    ekz6 *= hz_inv;

    // convert D-field to force
    type = atom->type[i];
    lj0 = B[7*type+6];
    lj1 = B[7*type+5];
    lj2 = B[7*type+4];
    lj3 = B[7*type+3];
    lj4 = B[7*type+2];
    lj5 = B[7*type+1];
    lj6 = B[7*type];

    s1 = x[i][0]*hx_inv;
    s2 = x[i][1]*hy_inv;
    s3 = x[i][2]*hz_inv;

    sf = sf_coeff_6[0]*sin(2*MY_PI*s1);
    sf += sf_coeff_6[1]*sin(4*MY_PI*s1);
    sf *= 4*lj0*lj6 + 4*lj1*lj5 + 4*lj2*lj4 + 2*lj3*lj3;
    f[i][0] += lj0*ekx0 + lj1*ekx1 + lj2*ekx2 + lj3*ekx3 + lj4*ekx4 + lj5*ekx5 + lj6*ekx6 - sf;

    sf = sf_coeff_6[2]*sin(2*MY_PI*s2);
    sf += sf_coeff_6[3]*sin(4*MY_PI*s2);
    sf *= 4*lj0*lj6 + 4*lj1*lj5 + 4*lj2*lj4 + 2*lj3*lj3;
    f[i][1] += lj0*eky0 + lj1*eky1 + lj2*eky2 + lj3*eky3 + lj4*eky4 + lj5*eky5 + lj6*eky6 - sf;

    sf = sf_coeff_6[4]*sin(2*MY_PI*s3);
    sf += sf_coeff_6[5]*sin(4*MY_PI*s3);
    sf *= 4*lj0*lj6 + 4*lj1*lj5 + 4*lj2*lj4 + 2*lj3*lj3;
    if (slabflag != 2) f[i][2] += lj0*ekz0 + lj1*ekz1 + lj2*ekz2 + lj3*ekz3 + lj4*ekz4 + lj5*ekz5 + lj6*ekz6 - sf;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for arithmetic mixing rule for per atom quantities
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_a_peratom()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR u_pa0,v00,v10,v20,v30,v40,v50;
  FFT_SCALAR u_pa1,v01,v11,v21,v31,v41,v51;
  FFT_SCALAR u_pa2,v02,v12,v22,v32,v42,v52;
  FFT_SCALAR u_pa3,v03,v13,v23,v33,v43,v53;
  FFT_SCALAR u_pa4,v04,v14,v24,v34,v44,v54;
  FFT_SCALAR u_pa5,v05,v15,v25,v35,v45,v55;
  FFT_SCALAR u_pa6,v06,v16,v26,v36,v46,v56;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  double **x = atom->x;
  int type;
  double lj0, lj1, lj2, lj3, lj4, lj5, lj6;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {

    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;
    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);

    u_pa0 = v00 = v10 = v20 = v30 = v40 = v50 = ZEROF;
    u_pa1 = v01 = v11 = v21 = v31 = v41 = v51 = ZEROF;
    u_pa2 = v02 = v12 = v22 = v32 = v42 = v52 = ZEROF;
    u_pa3 = v03 = v13 = v23 = v33 = v43 = v53 = ZEROF;
    u_pa4 = v04 = v14 = v24 = v34 = v44 = v54 = ZEROF;
    u_pa5 = v05 = v15 = v25 = v35 = v45 = v55 = ZEROF;
    u_pa6 = v06 = v16 = v26 = v36 = v46 = v56 = ZEROF;
    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      z0 = rho1d_6[2][n];
      for (m = nlower_6; m <= nupper_6; m++) {
	my = m+ny;
	y0 = z0*rho1d_6[1][m];
	for (l = nlower_6; l <= nupper_6; l++) {
	  mx = l+nx;
	  x0 = y0*rho1d_6[0][l];
          if (eflag_atom) {
            u_pa0 += x0*u_brick_a0[mz][my][mx];
            u_pa1 += x0*u_brick_a1[mz][my][mx];
            u_pa2 += x0*u_brick_a2[mz][my][mx];
            u_pa3 += x0*u_brick_a3[mz][my][mx];
            u_pa4 += x0*u_brick_a4[mz][my][mx];
            u_pa5 += x0*u_brick_a5[mz][my][mx];
            u_pa6 += x0*u_brick_a6[mz][my][mx];
	  }
          if (vflag_atom) {
            v00 += x0*v0_brick_a0[mz][my][mx];
            v10 += x0*v1_brick_a0[mz][my][mx];
            v20 += x0*v2_brick_a0[mz][my][mx];
            v30 += x0*v3_brick_a0[mz][my][mx];
            v40 += x0*v4_brick_a0[mz][my][mx];
            v50 += x0*v5_brick_a0[mz][my][mx];
            v01 += x0*v0_brick_a1[mz][my][mx];
            v11 += x0*v1_brick_a1[mz][my][mx];
            v21 += x0*v2_brick_a1[mz][my][mx];
            v31 += x0*v3_brick_a1[mz][my][mx];
            v41 += x0*v4_brick_a1[mz][my][mx];
            v51 += x0*v5_brick_a1[mz][my][mx];
            v02 += x0*v0_brick_a2[mz][my][mx];
            v12 += x0*v1_brick_a2[mz][my][mx];
            v22 += x0*v2_brick_a2[mz][my][mx];
            v32 += x0*v3_brick_a2[mz][my][mx];
            v42 += x0*v4_brick_a2[mz][my][mx];
            v52 += x0*v5_brick_a2[mz][my][mx];
            v03 += x0*v0_brick_a3[mz][my][mx];
            v13 += x0*v1_brick_a3[mz][my][mx];
            v23 += x0*v2_brick_a3[mz][my][mx];
            v33 += x0*v3_brick_a3[mz][my][mx];
            v43 += x0*v4_brick_a3[mz][my][mx];
            v53 += x0*v5_brick_a3[mz][my][mx];
            v04 += x0*v0_brick_a4[mz][my][mx];
            v14 += x0*v1_brick_a4[mz][my][mx];
            v24 += x0*v2_brick_a4[mz][my][mx];
            v34 += x0*v3_brick_a4[mz][my][mx];
            v44 += x0*v4_brick_a4[mz][my][mx];
            v54 += x0*v5_brick_a4[mz][my][mx];
            v05 += x0*v0_brick_a5[mz][my][mx];
            v15 += x0*v1_brick_a5[mz][my][mx];
            v25 += x0*v2_brick_a5[mz][my][mx];
            v35 += x0*v3_brick_a5[mz][my][mx];
            v45 += x0*v4_brick_a5[mz][my][mx];
            v55 += x0*v5_brick_a5[mz][my][mx];
            v06 += x0*v0_brick_a6[mz][my][mx];
            v16 += x0*v1_brick_a6[mz][my][mx];
            v26 += x0*v2_brick_a6[mz][my][mx];
            v36 += x0*v3_brick_a6[mz][my][mx];
            v46 += x0*v4_brick_a6[mz][my][mx];
            v56 += x0*v5_brick_a6[mz][my][mx];
          }
	}
      }
    }
    // convert D-field to force
    type = atom->type[i];
    lj0 = B[7*type+6]*0.5;
    lj1 = B[7*type+5]*0.5;
    lj2 = B[7*type+4]*0.5;
    lj3 = B[7*type+3]*0.5;
    lj4 = B[7*type+2]*0.5;
    lj5 = B[7*type+1]*0.5;
    lj6 = B[7*type]*0.5;

 
    if (eflag_atom) 
      eatom[i] += u_pa0*lj0 + u_pa1*lj1 + u_pa2*lj2 + 
        u_pa3*lj3 + u_pa4*lj4 + u_pa5*lj5 + u_pa6*lj6;
    if (vflag_atom) {
      vatom[i][0] += v00*lj0 + v01*lj1 + v02*lj2 + v03*lj3 + 
        v04*lj4 + v05*lj5 + v06*lj6;
      vatom[i][1] += v10*lj0 + v11*lj1 + v12*lj2 + v13*lj3 + 
        v14*lj4 + v15*lj5 + v16*lj6;
      vatom[i][2] += v20*lj0 + v21*lj1 + v22*lj2 + v23*lj3 + 
        v24*lj4 + v25*lj5 + v26*lj6;
      vatom[i][3] += v30*lj0 + v31*lj1 + v32*lj2 + v33*lj3 + 
        v34*lj4 + v35*lj5 + v36*lj6;
      vatom[i][4] += v40*lj0 + v41*lj1 + v42*lj2 + v43*lj3 + 
        v44*lj4 + v45*lj5 + v46*lj6;
      vatom[i][5] += v50*lj0 + v51*lj1 + v52*lj2 + v53*lj3 + 
        v54*lj4 + v55*lj5 + v56*lj6;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for arithmetic mixing rule and ik scheme
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_none_ik()
{
  int i,k,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR *ekx, *eky, *ekz;

  ekx = new FFT_SCALAR[nsplit];
  eky = new FFT_SCALAR[nsplit];
  ekz = new FFT_SCALAR[nsplit];
  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  double **x = atom->x;
  double **f = atom->f;
  int type;
  double lj;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {

    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;
    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);
    for (k = 0; k < nsplit; k++)
      ekx[k] = eky[k] = ekz[k] = ZEROF;
    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      z0 = rho1d_6[2][n];
      for (m = nlower_6; m <= nupper_6; m++) {
	my = m+ny;
	y0 = z0*rho1d_6[1][m];
	for (l = nlower_6; l <= nupper_6; l++) {
	  mx = l+nx;
	  x0 = y0*rho1d_6[0][l];
          for (k = 0; k < nsplit; k++) {
	    ekx[k] -= x0*vdx_brick_none[k][mz][my][mx];
	    eky[k] -= x0*vdy_brick_none[k][mz][my][mx];
	    ekz[k] -= x0*vdz_brick_none[k][mz][my][mx];
          }
	}
      }
    }
    // convert D-field to force
    type = atom->type[i];
    for (k = 0; k < nsplit; k++) {
      lj = B[nsplit*type + k];
      f[i][0] += lj*ekx[k];
      f[i][1] +=lj*eky[k];
      if (slabflag != 2) f[i][2] +=lj*ekz[k];
    }
  }

  delete [] ekx;
  delete [] eky;
  delete [] ekz;
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for arithmetic mixing rule for the ad scheme
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_none_ad()
{
  int i,k,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR *ekx, *eky, *ekz;

  ekx = new FFT_SCALAR[nsplit];
  eky = new FFT_SCALAR[nsplit];
  ekz = new FFT_SCALAR[nsplit];


  double s1,s2,s3;
  double sf1,sf2,sf3;
  double sf = 0.0;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double hx_inv = nx_pppm_6/xprd;
  double hy_inv = ny_pppm_6/yprd;
  double hz_inv = nz_pppm_6/zprd_slab;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  double **x = atom->x;
  double **f = atom->f;
  int type;
  double lj;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {

    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;

    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);
    compute_drho1d(dx,dy,dz, order_6, drho_coeff_6, drho1d_6);

    for (k = 0; k < nsplit; k++)
      ekx[k] = eky[k] = ekz[k] = ZEROF;

    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      for (m = nlower_6; m <= nupper_6; m++) {
	my = m+ny;
	for (l = nlower_6; l <= nupper_6; l++) {
	  mx = l+nx;
          x0 = drho1d_6[0][l]*rho1d_6[1][m]*rho1d_6[2][n];
          y0 = rho1d_6[0][l]*drho1d_6[1][m]*rho1d_6[2][n];
          z0 = rho1d_6[0][l]*rho1d_6[1][m]*drho1d_6[2][n];
          
          for (k = 0; k < nsplit; k++) {
            ekx[k] += x0*u_brick_none[k][mz][my][mx];
            eky[k] += y0*u_brick_none[k][mz][my][mx];
            ekz[k] += z0*u_brick_none[k][mz][my][mx];
          }
	}
      }
    }

    for (k = 0; k < nsplit; k++) {
      ekx[k] *= hx_inv;
      eky[k] *= hy_inv;
      ekz[k] *= hz_inv;
    }

    // convert D-field to force
    type = atom->type[i];

    s1 = x[i][0]*hx_inv;
    s2 = x[i][1]*hy_inv;
    s3 = x[i][2]*hz_inv;

    sf1 = sf_coeff_6[0]*sin(2*MY_PI*s1);
    sf1 += sf_coeff_6[1]*sin(4*MY_PI*s1);

    sf2 = sf_coeff_6[2]*sin(2*MY_PI*s2);
    sf2 += sf_coeff_6[3]*sin(4*MY_PI*s2);

    sf3 = sf_coeff_6[4]*sin(2*MY_PI*s3);
    sf3 += sf_coeff_6[5]*sin(4*MY_PI*s3);

    for (k = 0; k < nsplit; k++) {
      lj = B[nsplit*type + k];

      sf = sf1*B[k]*2*lj*lj;
      f[i][0] += lj*ekx[k] - sf;


      sf = sf2*B[k]*2*lj*lj;
      f[i][1] += lj*eky[k] - sf;

      sf = sf3*B[k]*2*lj*lj;
      if (slabflag != 2) f[i][2] += lj*ekz[k] - sf;
    }
  }

  delete [] ekx;
  delete [] eky;
  delete [] ekz;
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for arithmetic mixing rule for per atom quantities
------------------------------------------------------------------------- */

void PPPMDisp::fieldforce_none_peratom()
{
  int i,k,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR *u_pa,*v0,*v1,*v2,*v3,*v4,*v5;
  
  u_pa = new FFT_SCALAR[nsplit];
  v0 = new FFT_SCALAR[nsplit];
  v1 = new FFT_SCALAR[nsplit];
  v2 = new FFT_SCALAR[nsplit];
  v3 = new FFT_SCALAR[nsplit];
  v4 = new FFT_SCALAR[nsplit];
  v5 = new FFT_SCALAR[nsplit];

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  double **x = atom->x;
  int type;
  double lj;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {

    nx = part2grid_6[i][0];
    ny = part2grid_6[i][1];
    nz = part2grid_6[i][2];
    dx = nx+shiftone_6 - (x[i][0]-boxlo[0])*delxinv_6;
    dy = ny+shiftone_6 - (x[i][1]-boxlo[1])*delyinv_6;
    dz = nz+shiftone_6 - (x[i][2]-boxlo[2])*delzinv_6;
    compute_rho1d(dx,dy,dz, order_6, rho_coeff_6, rho1d_6);

    for (k = 0; k < nsplit; k++) 
      u_pa[k] = v0[k] = v1[k] = v2[k] = v3[k] = v4[k] = v5[k] = ZEROF;
 
    for (n = nlower_6; n <= nupper_6; n++) {
      mz = n+nz;
      z0 = rho1d_6[2][n];
      for (m = nlower_6; m <= nupper_6; m++) {
	my = m+ny;
	y0 = z0*rho1d_6[1][m];
	for (l = nlower_6; l <= nupper_6; l++) {
	  mx = l+nx;
	  x0 = y0*rho1d_6[0][l];
          if (eflag_atom) {
            for (k = 0; k < nsplit; k++)
              u_pa[k] += x0*u_brick_none[k][mz][my][mx];
	  }
          if (vflag_atom) {
            for (k = 0; k < nsplit; k++) {
              v0[k] += x0*v0_brick_none[k][mz][my][mx];
              v1[k] += x0*v1_brick_none[k][mz][my][mx];
              v2[k] += x0*v2_brick_none[k][mz][my][mx];
              v3[k] += x0*v3_brick_none[k][mz][my][mx];
              v4[k] += x0*v4_brick_none[k][mz][my][mx];
              v5[k] += x0*v5_brick_none[k][mz][my][mx];
            }
          }
	}
      }
    }
    // convert D-field to force
    type = atom->type[i];
    for (k = 0; k < nsplit; k++) {
      lj = B[nsplit*type + k]*0.5;
 
      if (eflag_atom) {
        eatom[i] += u_pa[k]*lj;
      }
      if (vflag_atom) {
        vatom[i][0] += v0[k]*lj;
        vatom[i][1] += v1[k]*lj;
        vatom[i][2] += v2[k]*lj;
        vatom[i][3] += v3[k]*lj;
        vatom[i][4] += v4[k]*lj;
        vatom[i][5] += v5[k]*lj;
      }
    }
  }

  delete [] u_pa;
  delete [] v0;
  delete [] v1;
  delete [] v2;
  delete [] v3;
  delete [] v4;
  delete [] v5;
}

/* ----------------------------------------------------------------------
   pack values to buf to send to another proc
------------------------------------------------------------------------- */

void PPPMDisp::pack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;

  switch (flag) {

  // Coulomb interactions

  case FORWARD_IK: {
    FFT_SCALAR *xsrc = &vdx_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *ysrc = &vdy_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *zsrc = &vdz_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = xsrc[list[i]];
      buf[n++] = ysrc[list[i]];
      buf[n++] = zsrc[list[i]];
    }
    break;
  }

  case FORWARD_AD: {
    FFT_SCALAR *src = &u_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
    break;
  }

  case FORWARD_IK_PERATOM: {
    FFT_SCALAR *esrc = &u_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) buf[n++] = esrc[list[i]];
      if (vflag_atom) {
        buf[n++] = v0src[list[i]];
        buf[n++] = v1src[list[i]];
        buf[n++] = v2src[list[i]];
        buf[n++] = v3src[list[i]];
        buf[n++] = v4src[list[i]];
        buf[n++] = v5src[list[i]];
      }
    }
    break;
  }

  case FORWARD_AD_PERATOM: {
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = v0src[list[i]];
      buf[n++] = v1src[list[i]];
      buf[n++] = v2src[list[i]];
      buf[n++] = v3src[list[i]];
      buf[n++] = v4src[list[i]];
      buf[n++] = v5src[list[i]];
    }
    break;
  }

  // Dispersion interactions, geometric mixing

  case FORWARD_IK_G: {
    FFT_SCALAR *xsrc = &vdx_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ysrc = &vdy_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zsrc = &vdz_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = xsrc[list[i]];
      buf[n++] = ysrc[list[i]];
      buf[n++] = zsrc[list[i]];
    }
    break;
  }

  case FORWARD_AD_G: {
    FFT_SCALAR *src = &u_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
    break;
  }

  case FORWARD_IK_PERATOM_G: {
    FFT_SCALAR *esrc = &u_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src = &v0_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src = &v1_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src = &v2_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src = &v3_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src = &v4_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src = &v5_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) buf[n++] = esrc[list[i]];
      if (vflag_atom) {
        buf[n++] = v0src[list[i]];
        buf[n++] = v1src[list[i]];
        buf[n++] = v2src[list[i]];
        buf[n++] = v3src[list[i]];
        buf[n++] = v4src[list[i]];
        buf[n++] = v5src[list[i]];
      }
    }
    break;
  }

  case FORWARD_AD_PERATOM_G: {
    FFT_SCALAR *v0src = &v0_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src = &v1_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src = &v2_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src = &v3_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src = &v4_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src = &v5_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = v0src[list[i]];
      buf[n++] = v1src[list[i]];
      buf[n++] = v2src[list[i]];
      buf[n++] = v3src[list[i]];
      buf[n++] = v4src[list[i]];
      buf[n++] = v5src[list[i]];
    }
    break;
  }

  // Dispersion interactions, arithmetic mixing

  case FORWARD_IK_A: {
    FFT_SCALAR *xsrc0 = &vdx_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ysrc0 = &vdy_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zsrc0 = &vdz_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xsrc1 = &vdx_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ysrc1 = &vdy_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zsrc1 = &vdz_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xsrc2 = &vdx_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ysrc2 = &vdy_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zsrc2 = &vdz_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xsrc3 = &vdx_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ysrc3 = &vdy_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zsrc3 = &vdz_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xsrc4 = &vdx_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ysrc4 = &vdy_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zsrc4 = &vdz_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xsrc5 = &vdx_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ysrc5 = &vdy_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zsrc5 = &vdz_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xsrc6 = &vdx_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ysrc6 = &vdy_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zsrc6 = &vdz_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];

    for (int i = 0; i < nlist; i++) {
      buf[n++] = xsrc0[list[i]];
      buf[n++] = ysrc0[list[i]];
      buf[n++] = zsrc0[list[i]];

      buf[n++] = xsrc1[list[i]];
      buf[n++] = ysrc1[list[i]];
      buf[n++] = zsrc1[list[i]];

      buf[n++] = xsrc2[list[i]];
      buf[n++] = ysrc2[list[i]];
      buf[n++] = zsrc2[list[i]];

      buf[n++] = xsrc3[list[i]];
      buf[n++] = ysrc3[list[i]];
      buf[n++] = zsrc3[list[i]];

      buf[n++] = xsrc4[list[i]];
      buf[n++] = ysrc4[list[i]];
      buf[n++] = zsrc4[list[i]];

      buf[n++] = xsrc5[list[i]];
      buf[n++] = ysrc5[list[i]];
      buf[n++] = zsrc5[list[i]];

      buf[n++] = xsrc6[list[i]];
      buf[n++] = ysrc6[list[i]];
      buf[n++] = zsrc6[list[i]];
    }
    break;
  }

  case FORWARD_AD_A: {
    FFT_SCALAR *src0 = &u_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src1 = &u_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src2 = &u_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src3 = &u_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src4 = &u_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src5 = &u_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src6 = &u_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];

    for (int i = 0; i < nlist; i++) {
      buf[n++] = src0[list[i]];
      buf[n++] = src1[list[i]];
      buf[n++] = src2[list[i]];
      buf[n++] = src3[list[i]];
      buf[n++] = src4[list[i]];
      buf[n++] = src5[list[i]];
      buf[n++] = src6[list[i]];
    }
    break;
  }

  case FORWARD_IK_PERATOM_A: {
    FFT_SCALAR *esrc0 = &u_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src0 = &v0_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src0 = &v1_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src0 = &v2_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src0 = &v3_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src0 = &v4_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src0 = &v5_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc1 = &u_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src1 = &v0_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src1 = &v1_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src1 = &v2_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src1 = &v3_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src1 = &v4_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src1 = &v5_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc2 = &u_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src2 = &v0_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src2 = &v1_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src2 = &v2_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src2 = &v3_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src2 = &v4_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src2 = &v5_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc3 = &u_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src3 = &v0_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src3 = &v1_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src3 = &v2_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src3 = &v3_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src3 = &v4_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src3 = &v5_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc4 = &u_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src4 = &v0_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src4 = &v1_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src4 = &v2_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src4 = &v3_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src4 = &v4_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src4 = &v5_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc5 = &u_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src5 = &v0_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src5 = &v1_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src5 = &v2_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src5 = &v3_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src5 = &v4_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src5 = &v5_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc6 = &u_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src6 = &v0_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src6 = &v1_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src6 = &v2_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src6 = &v3_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src6 = &v4_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src6 = &v5_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];

    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) {
        buf[n++] = esrc0[list[i]];
        buf[n++] = esrc1[list[i]];
        buf[n++] = esrc2[list[i]];
        buf[n++] = esrc3[list[i]];
        buf[n++] = esrc4[list[i]];
        buf[n++] = esrc5[list[i]];
        buf[n++] = esrc6[list[i]];
      }
      if (vflag_atom) {
        buf[n++] = v0src0[list[i]];
        buf[n++] = v1src0[list[i]];
        buf[n++] = v2src0[list[i]];
        buf[n++] = v3src0[list[i]];
        buf[n++] = v4src0[list[i]];
        buf[n++] = v5src0[list[i]];

        buf[n++] = v0src1[list[i]];
        buf[n++] = v1src1[list[i]];
        buf[n++] = v2src1[list[i]];
        buf[n++] = v3src1[list[i]];
        buf[n++] = v4src1[list[i]];
        buf[n++] = v5src1[list[i]];

        buf[n++] = v0src2[list[i]];
        buf[n++] = v1src2[list[i]];
        buf[n++] = v2src2[list[i]];
        buf[n++] = v3src2[list[i]];
        buf[n++] = v4src2[list[i]];
        buf[n++] = v5src2[list[i]];

        buf[n++] = v0src3[list[i]];
        buf[n++] = v1src3[list[i]];
        buf[n++] = v2src3[list[i]];
        buf[n++] = v3src3[list[i]];
        buf[n++] = v4src3[list[i]];
        buf[n++] = v5src3[list[i]];

        buf[n++] = v0src4[list[i]];
        buf[n++] = v1src4[list[i]];
        buf[n++] = v2src4[list[i]];
        buf[n++] = v3src4[list[i]];
        buf[n++] = v4src4[list[i]];
        buf[n++] = v5src4[list[i]];

        buf[n++] = v0src5[list[i]];
        buf[n++] = v1src5[list[i]];
        buf[n++] = v2src5[list[i]];
        buf[n++] = v3src5[list[i]];
        buf[n++] = v4src5[list[i]];
        buf[n++] = v5src5[list[i]];

        buf[n++] = v0src6[list[i]];
        buf[n++] = v1src6[list[i]];
        buf[n++] = v2src6[list[i]];
        buf[n++] = v3src6[list[i]];
        buf[n++] = v4src6[list[i]];
        buf[n++] = v5src6[list[i]];
      }
    }
    break;
  }

  case FORWARD_AD_PERATOM_A: {
    FFT_SCALAR *v0src0 = &v0_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src0 = &v1_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src0 = &v2_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src0 = &v3_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src0 = &v4_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src0 = &v5_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src1 = &v0_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src1 = &v1_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src1 = &v2_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src1 = &v3_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src1 = &v4_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src1 = &v5_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src2 = &v0_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src2 = &v1_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src2 = &v2_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src2 = &v3_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src2 = &v4_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src2 = &v5_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src3 = &v0_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src3 = &v1_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src3 = &v2_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src3 = &v3_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src3 = &v4_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src3 = &v5_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src4 = &v0_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src4 = &v1_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src4 = &v2_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src4 = &v3_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src4 = &v4_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src4 = &v5_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src5 = &v0_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src5 = &v1_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src5 = &v2_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src5 = &v3_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src5 = &v4_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src5 = &v5_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src6 = &v0_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src6 = &v1_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src6 = &v2_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src6 = &v3_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src6 = &v4_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src6 = &v5_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];

    for (int i = 0; i < nlist; i++) {
      buf[n++] = v0src0[list[i]];
      buf[n++] = v1src0[list[i]];
      buf[n++] = v2src0[list[i]];
      buf[n++] = v3src0[list[i]];
      buf[n++] = v4src0[list[i]];
      buf[n++] = v5src0[list[i]];

      buf[n++] = v0src1[list[i]];
      buf[n++] = v1src1[list[i]];
      buf[n++] = v2src1[list[i]];
      buf[n++] = v3src1[list[i]];
      buf[n++] = v4src1[list[i]];
      buf[n++] = v5src1[list[i]];

      buf[n++] = v0src2[list[i]];
      buf[n++] = v1src2[list[i]];
      buf[n++] = v2src2[list[i]];
      buf[n++] = v3src2[list[i]];
      buf[n++] = v4src2[list[i]];
      buf[n++] = v5src2[list[i]];

      buf[n++] = v0src3[list[i]];
      buf[n++] = v1src3[list[i]];
      buf[n++] = v2src3[list[i]];
      buf[n++] = v3src3[list[i]];
      buf[n++] = v4src3[list[i]];
      buf[n++] = v5src3[list[i]];

      buf[n++] = v0src4[list[i]];
      buf[n++] = v1src4[list[i]];
      buf[n++] = v2src4[list[i]];
      buf[n++] = v3src4[list[i]];
      buf[n++] = v4src4[list[i]];
      buf[n++] = v5src4[list[i]];

      buf[n++] = v0src5[list[i]];
      buf[n++] = v1src5[list[i]];
      buf[n++] = v2src5[list[i]];
      buf[n++] = v3src5[list[i]];
      buf[n++] = v4src5[list[i]];
      buf[n++] = v5src5[list[i]];

      buf[n++] = v0src6[list[i]];
      buf[n++] = v1src6[list[i]];
      buf[n++] = v2src6[list[i]];
      buf[n++] = v3src6[list[i]];
      buf[n++] = v4src6[list[i]];
      buf[n++] = v5src6[list[i]];
    }
    break;
  }

  // Dispersion interactions, no mixing

  case FORWARD_IK_NONE: {
    for (int k = 0; k < nsplit_alloc; k++) {
      FFT_SCALAR *xsrc = &vdx_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *ysrc = &vdy_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *zsrc = &vdz_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      for (int i = 0; i < nlist; i++) {
        buf[n++] = xsrc[list[i]];
        buf[n++] = ysrc[list[i]];
        buf[n++] = zsrc[list[i]];
      }
    }
    break;
  }

  case FORWARD_AD_NONE: {
    for (int k = 0; k < nsplit_alloc; k++) {
      FFT_SCALAR *src = &u_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      for (int i = 0; i < nlist; i++)
        buf[n++] = src[list[i]];
    }
    break;
  }

  case FORWARD_IK_PERATOM_NONE: {
    for (int k = 0; k < nsplit_alloc; k++) {
      FFT_SCALAR *esrc = &u_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v0src = &v0_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v1src = &v1_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v2src = &v2_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v3src = &v3_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v4src = &v4_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v5src = &v5_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      for (int i = 0; i < nlist; i++) {
        if (eflag_atom) buf[n++] = esrc[list[i]];
        if (vflag_atom) {
          buf[n++] = v0src[list[i]];
          buf[n++] = v1src[list[i]];
          buf[n++] = v2src[list[i]];
          buf[n++] = v3src[list[i]];
          buf[n++] = v4src[list[i]];
          buf[n++] = v5src[list[i]];
        }
      }
    }
    break;
  }

  case FORWARD_AD_PERATOM_NONE: {
    for (int k = 0; k < nsplit_alloc; k++) {
      FFT_SCALAR *v0src = &v0_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v1src = &v1_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v2src = &v2_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v3src = &v3_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v4src = &v4_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v5src = &v5_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      for (int i = 0; i < nlist; i++) {
        buf[n++] = v0src[list[i]];
        buf[n++] = v1src[list[i]];
        buf[n++] = v2src[list[i]];
        buf[n++] = v3src[list[i]];
        buf[n++] = v4src[list[i]];
        buf[n++] = v5src[list[i]];
      }
    }
    break;
  }

  }
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void PPPMDisp::unpack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;

  switch (flag) {

  // Coulomb interactions

  case FORWARD_IK: {
    FFT_SCALAR *xdest = &vdx_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *ydest = &vdy_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *zdest = &vdz_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      xdest[list[i]] = buf[n++];
      ydest[list[i]] = buf[n++];
      zdest[list[i]] = buf[n++];
    }
    break;
  }

  case FORWARD_AD: {
    FFT_SCALAR *dest = &u_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[n++];
    break;
  }

  case FORWARD_IK_PERATOM: {
    FFT_SCALAR *esrc = &u_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) esrc[list[i]] = buf[n++];
      if (vflag_atom) {
        v0src[list[i]] = buf[n++];
        v1src[list[i]] = buf[n++];
        v2src[list[i]] = buf[n++];
        v3src[list[i]] = buf[n++];
        v4src[list[i]] = buf[n++];
        v5src[list[i]] = buf[n++];
      }
    }
    break;
  }

  case FORWARD_AD_PERATOM: {
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      v0src[list[i]] = buf[n++];
      v1src[list[i]] = buf[n++];
      v2src[list[i]] = buf[n++];
      v3src[list[i]] = buf[n++];
      v4src[list[i]] = buf[n++];
      v5src[list[i]] = buf[n++];
    }
    break;
  }

  // Disperion interactions, geometric mixing

  case FORWARD_IK_G: {
    FFT_SCALAR *xdest = &vdx_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ydest = &vdy_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zdest = &vdz_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++) {
      xdest[list[i]] = buf[n++];
      ydest[list[i]] = buf[n++];
      zdest[list[i]] = buf[n++];
    }
    break;
  }

  case FORWARD_AD_G: {
    FFT_SCALAR *dest = &u_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[n++];
    break;
  }

  case FORWARD_IK_PERATOM_G: {
    FFT_SCALAR *esrc = &u_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src = &v0_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src = &v1_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src = &v2_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src = &v3_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src = &v4_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src = &v5_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) esrc[list[i]] = buf[n++];
      if (vflag_atom) {
        v0src[list[i]] = buf[n++];
        v1src[list[i]] = buf[n++];
        v2src[list[i]] = buf[n++];
        v3src[list[i]] = buf[n++];
        v4src[list[i]] = buf[n++];
        v5src[list[i]] = buf[n++];
      }
    }
    break;
  }

  case FORWARD_AD_PERATOM_G: {
    FFT_SCALAR *v0src = &v0_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src = &v1_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src = &v2_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src = &v3_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src = &v4_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src = &v5_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++) {
      v0src[list[i]] = buf[n++];
      v1src[list[i]] = buf[n++];
      v2src[list[i]] = buf[n++];
      v3src[list[i]] = buf[n++];
      v4src[list[i]] = buf[n++];
      v5src[list[i]] = buf[n++];
    }
    break;
  }

  // Disperion interactions, arithmetic mixing

  case FORWARD_IK_A: {
    FFT_SCALAR *xdest0 = &vdx_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ydest0 = &vdy_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zdest0 = &vdz_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xdest1 = &vdx_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ydest1 = &vdy_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zdest1 = &vdz_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xdest2 = &vdx_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ydest2 = &vdy_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zdest2 = &vdz_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xdest3 = &vdx_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ydest3 = &vdy_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zdest3 = &vdz_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xdest4 = &vdx_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ydest4 = &vdy_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zdest4 = &vdz_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xdest5 = &vdx_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ydest5 = &vdy_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zdest5 = &vdz_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *xdest6 = &vdx_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *ydest6 = &vdy_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *zdest6 = &vdz_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];

    for (int i = 0; i < nlist; i++) {
      xdest0[list[i]] = buf[n++];
      ydest0[list[i]] = buf[n++];
      zdest0[list[i]] = buf[n++];

      xdest1[list[i]] = buf[n++];
      ydest1[list[i]] = buf[n++];
      zdest1[list[i]] = buf[n++];

      xdest2[list[i]] = buf[n++];
      ydest2[list[i]] = buf[n++];
      zdest2[list[i]] = buf[n++];

      xdest3[list[i]] = buf[n++];
      ydest3[list[i]] = buf[n++];
      zdest3[list[i]] = buf[n++];

      xdest4[list[i]] = buf[n++];
      ydest4[list[i]] = buf[n++];
      zdest4[list[i]] = buf[n++];

      xdest5[list[i]] = buf[n++];
      ydest5[list[i]] = buf[n++];
      zdest5[list[i]] = buf[n++];

      xdest6[list[i]] = buf[n++];
      ydest6[list[i]] = buf[n++];
      zdest6[list[i]] = buf[n++];
    }
    break;
  }

  case FORWARD_AD_A: {
    FFT_SCALAR *dest0 = &u_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest1 = &u_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest2 = &u_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest3 = &u_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest4 = &u_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest5 = &u_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest6 = &u_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];

    for (int i = 0; i < nlist; i++) {
      dest0[list[i]] = buf[n++];
      dest1[list[i]] = buf[n++];
      dest2[list[i]] = buf[n++];
      dest3[list[i]] = buf[n++];
      dest4[list[i]] = buf[n++];
      dest5[list[i]] = buf[n++];
      dest6[list[i]] = buf[n++];
    }
    break;
  }

  case FORWARD_IK_PERATOM_A: {
    FFT_SCALAR *esrc0 = &u_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src0 = &v0_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src0 = &v1_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src0 = &v2_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src0 = &v3_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src0 = &v4_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src0 = &v5_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc1 = &u_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src1 = &v0_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src1 = &v1_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src1 = &v2_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src1 = &v3_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src1 = &v4_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src1 = &v5_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc2 = &u_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src2 = &v0_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src2 = &v1_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src2 = &v2_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src2 = &v3_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src2 = &v4_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src2 = &v5_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc3 = &u_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src3 = &v0_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src3 = &v1_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src3 = &v2_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src3 = &v3_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src3 = &v4_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src3 = &v5_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc4 = &u_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src4 = &v0_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src4 = &v1_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src4 = &v2_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src4 = &v3_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src4 = &v4_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src4 = &v5_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc5 = &u_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src5 = &v0_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src5 = &v1_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src5 = &v2_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src5 = &v3_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src5 = &v4_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src5 = &v5_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *esrc6 = &u_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v0src6 = &v0_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src6 = &v1_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src6 = &v2_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src6 = &v3_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src6 = &v4_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src6 = &v5_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];

    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) {
        esrc0[list[i]] = buf[n++];
        esrc1[list[i]] = buf[n++];
        esrc2[list[i]] = buf[n++];
        esrc3[list[i]] = buf[n++];
        esrc4[list[i]] = buf[n++];
        esrc5[list[i]] = buf[n++];
        esrc6[list[i]] = buf[n++];
      }
      if (vflag_atom) {
        v0src0[list[i]] = buf[n++];
        v1src0[list[i]] = buf[n++];
        v2src0[list[i]] = buf[n++];
        v3src0[list[i]] = buf[n++];
        v4src0[list[i]] = buf[n++];
        v5src0[list[i]] = buf[n++];

        v0src1[list[i]] = buf[n++];
        v1src1[list[i]] = buf[n++];
        v2src1[list[i]] = buf[n++];
        v3src1[list[i]] = buf[n++];
        v4src1[list[i]] = buf[n++];
        v5src1[list[i]] = buf[n++];

        v0src2[list[i]] = buf[n++];
        v1src2[list[i]] = buf[n++];
        v2src2[list[i]] = buf[n++];
        v3src2[list[i]] = buf[n++];
        v4src2[list[i]] = buf[n++];
        v5src2[list[i]] = buf[n++];

        v0src3[list[i]] = buf[n++];
        v1src3[list[i]] = buf[n++];
        v2src3[list[i]] = buf[n++];
        v3src3[list[i]] = buf[n++];
        v4src3[list[i]] = buf[n++];
        v5src3[list[i]] = buf[n++];

        v0src4[list[i]] = buf[n++];
        v1src4[list[i]] = buf[n++];
        v2src4[list[i]] = buf[n++];
        v3src4[list[i]] = buf[n++];
        v4src4[list[i]] = buf[n++];
        v5src4[list[i]] = buf[n++];

        v0src5[list[i]] = buf[n++];
        v1src5[list[i]] = buf[n++];
        v2src5[list[i]] = buf[n++];
        v3src5[list[i]] = buf[n++];
        v4src5[list[i]] = buf[n++];
        v5src5[list[i]] = buf[n++];

        v0src6[list[i]] = buf[n++];
        v1src6[list[i]] = buf[n++];
        v2src6[list[i]] = buf[n++];
        v3src6[list[i]] = buf[n++];
        v4src6[list[i]] = buf[n++];
        v5src6[list[i]] = buf[n++];
      }
    }
    break;
  }

  case FORWARD_AD_PERATOM_A: {
    FFT_SCALAR *v0src0 = &v0_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src0 = &v1_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src0 = &v2_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src0 = &v3_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src0 = &v4_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src0 = &v5_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src1 = &v0_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src1 = &v1_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src1 = &v2_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src1 = &v3_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src1 = &v4_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src1 = &v5_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src2 = &v0_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src2 = &v1_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src2 = &v2_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src2 = &v3_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src2 = &v4_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src2 = &v5_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src3 = &v0_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src3 = &v1_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src3 = &v2_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src3 = &v3_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src3 = &v4_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src3 = &v5_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src4 = &v0_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src4 = &v1_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src4 = &v2_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src4 = &v3_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src4 = &v4_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src4 = &v5_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src5 = &v0_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src5 = &v1_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src5 = &v2_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src5 = &v3_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src5 = &v4_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src5 = &v5_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];

    FFT_SCALAR *v0src6 = &v0_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v1src6 = &v1_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v2src6 = &v2_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v3src6 = &v3_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v4src6 = &v4_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *v5src6 = &v5_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];

    for (int i = 0; i < nlist; i++) {
      v0src0[list[i]] = buf[n++];
      v1src0[list[i]] = buf[n++];
      v2src0[list[i]] = buf[n++];
      v3src0[list[i]] = buf[n++];
      v4src0[list[i]] = buf[n++];
      v5src0[list[i]] = buf[n++];

      v0src1[list[i]] = buf[n++];
      v1src1[list[i]] = buf[n++];
      v2src1[list[i]] = buf[n++];
      v3src1[list[i]] = buf[n++];
      v4src1[list[i]] = buf[n++];
      v5src1[list[i]] = buf[n++];

      v0src2[list[i]] = buf[n++];
      v1src2[list[i]] = buf[n++];
      v2src2[list[i]] = buf[n++];
      v3src2[list[i]] = buf[n++];
      v4src2[list[i]] = buf[n++];
      v5src2[list[i]] = buf[n++];

      v0src3[list[i]] = buf[n++];
      v1src3[list[i]] = buf[n++];
      v2src3[list[i]] = buf[n++];
      v3src3[list[i]] = buf[n++];
      v4src3[list[i]] = buf[n++];
      v5src3[list[i]] = buf[n++];

      v0src4[list[i]] = buf[n++];
      v1src4[list[i]] = buf[n++];
      v2src4[list[i]] = buf[n++];
      v3src4[list[i]] = buf[n++];
      v4src4[list[i]] = buf[n++];
      v5src4[list[i]] = buf[n++];

      v0src5[list[i]] = buf[n++];
      v1src5[list[i]] = buf[n++];
      v2src5[list[i]] = buf[n++];
      v3src5[list[i]] = buf[n++];
      v4src5[list[i]] = buf[n++];
      v5src5[list[i]] = buf[n++];

      v0src6[list[i]] = buf[n++];
      v1src6[list[i]] = buf[n++];
      v2src6[list[i]] = buf[n++];
      v3src6[list[i]] = buf[n++];
      v4src6[list[i]] = buf[n++];
      v5src6[list[i]] = buf[n++];
    }
    break;
  }

  // Disperion interactions, geometric mixing

  case FORWARD_IK_NONE: {
    for (int k = 0; k < nsplit_alloc; k++) {
      FFT_SCALAR *xdest = &vdx_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *ydest = &vdy_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *zdest = &vdz_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      for (int i = 0; i < nlist; i++) {
        xdest[list[i]] = buf[n++];
        ydest[list[i]] = buf[n++];
        zdest[list[i]] = buf[n++];
      }
    }
    break;
  }

  case FORWARD_AD_NONE: {
    for (int k = 0; k < nsplit_alloc; k++) {
      FFT_SCALAR *dest = &u_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      for (int i = 0; i < nlist; i++)
        dest[list[i]] = buf[n++];
    }
    break;
  }

  case FORWARD_IK_PERATOM_NONE: {
    for (int k = 0; k < nsplit_alloc; k++) {
      FFT_SCALAR *esrc = &u_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v0src = &v0_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v1src = &v1_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v2src = &v2_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v3src = &v3_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v4src = &v4_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v5src = &v5_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      for (int i = 0; i < nlist; i++) {
        if (eflag_atom) esrc[list[i]] = buf[n++];
        if (vflag_atom) {
          v0src[list[i]] = buf[n++];
          v1src[list[i]] = buf[n++];
          v2src[list[i]] = buf[n++];
          v3src[list[i]] = buf[n++];
          v4src[list[i]] = buf[n++];
          v5src[list[i]] = buf[n++];
        }
      }
    }
    break;
  }

  case FORWARD_AD_PERATOM_NONE: {
    for (int k = 0; k < nsplit_alloc; k++) {
      FFT_SCALAR *v0src = &v0_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v1src = &v1_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v2src = &v2_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v3src = &v3_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v4src = &v4_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      FFT_SCALAR *v5src = &v5_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      for (int i = 0; i < nlist; i++) {
        v0src[list[i]] = buf[n++];
        v1src[list[i]] = buf[n++];
        v2src[list[i]] = buf[n++];
        v3src[list[i]] = buf[n++];
        v4src[list[i]] = buf[n++];
        v5src[list[i]] = buf[n++];
      }
    }
    break;
  }

  }
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void PPPMDisp::pack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;

  //Coulomb interactions

  if (flag == REVERSE_RHO) {
    FFT_SCALAR *src = &density_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];

  //Dispersion interactions, geometric mixing

  } else if (flag == REVERSE_RHO_G) {
    FFT_SCALAR *src = &density_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];

  //Dispersion interactions, arithmetic mixing

  } else if (flag == REVERSE_RHO_A) {
    FFT_SCALAR *src0 = &density_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src1 = &density_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src2 = &density_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src3 = &density_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src4 = &density_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src5 = &density_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *src6 = &density_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = src0[list[i]];
      buf[n++] = src1[list[i]];
      buf[n++] = src2[list[i]];
      buf[n++] = src3[list[i]];
      buf[n++] = src4[list[i]];
      buf[n++] = src5[list[i]];
      buf[n++] = src6[list[i]];
    }

  //Dispersion interactions, no mixing

  } else if (flag == REVERSE_RHO_NONE) {
    for (int k = 0; k < nsplit_alloc; k++) {
      FFT_SCALAR *src = &density_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      for (int i = 0; i < nlist; i++) {
        buf[n++] = src[list[i]];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void PPPMDisp::unpack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;

  //Coulomb interactions

  if (flag == REVERSE_RHO) {
    FFT_SCALAR *dest = &density_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];

  //Dispersion interactions, geometric mixing

  } else if (flag == REVERSE_RHO_G) {
    FFT_SCALAR *dest = &density_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];

  //Dispersion interactions, arithmetic mixing

  } else if (flag == REVERSE_RHO_A) {
    FFT_SCALAR *dest0 = &density_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest1 = &density_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest2 = &density_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest3 = &density_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest4 = &density_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest5 = &density_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6];
    FFT_SCALAR *dest6 = &density_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6];
    for (int i = 0; i < nlist; i++) {
      dest0[list[i]] += buf[n++];
      dest1[list[i]] += buf[n++];
      dest2[list[i]] += buf[n++];
      dest3[list[i]] += buf[n++];
      dest4[list[i]] += buf[n++];
      dest5[list[i]] += buf[n++];
      dest6[list[i]] += buf[n++];
    }

  //Dispersion interactions, no mixing

  } else if (flag == REVERSE_RHO_NONE) {
    for (int k = 0; k < nsplit_alloc; k++) {
      FFT_SCALAR *dest = &density_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6];
      for (int i = 0; i < nlist; i++)
        dest[list[i]] += buf[n++];
    }
  }
}

/* ----------------------------------------------------------------------
   map nprocs to NX by NY grid as PX by PY procs - return optimal px,py 
------------------------------------------------------------------------- */

void PPPMDisp::procs2grid2d(int nprocs, int nx, int ny, int *px, int *py)
{
  // loop thru all possible factorizations of nprocs
  // surf = surface area of largest proc sub-domain
  // innermost if test minimizes surface area and surface/volume ratio

  int bestsurf = 2 * (nx + ny);
  int bestboxx = 0;
  int bestboxy = 0;

  int boxx,boxy,surf,ipx,ipy;

  ipx = 1;
  while (ipx <= nprocs) {
    if (nprocs % ipx == 0) {
      ipy = nprocs/ipx;
      boxx = nx/ipx;
      if (nx % ipx) boxx++;
      boxy = ny/ipy;
      if (ny % ipy) boxy++;
      surf = boxx + boxy;
      if (surf < bestsurf || 
	  (surf == bestsurf && boxx*boxy > bestboxx*bestboxy)) {
	bestsurf = surf;
	bestboxx = boxx;
	bestboxy = boxy;
	*px = ipx;
	*py = ipy;
      }
    }
    ipx++;
  }
}

/* ----------------------------------------------------------------------
   charge assignment into rho1d
   dx,dy,dz = distance of particle from "lower left" grid point 
------------------------------------------------------------------------- */

void PPPMDisp::compute_rho1d(const FFT_SCALAR &dx, const FFT_SCALAR &dy,
			      const FFT_SCALAR &dz, int ord, 
                             FFT_SCALAR **rho_c, FFT_SCALAR **r1d)
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  for (k = (1-ord)/2; k <= ord/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = ord-1; l >= 0; l--) {
      r1 = rho_c[l][k] + r1*dx;
      r2 = rho_c[l][k] + r2*dy;
      r3 = rho_c[l][k] + r3*dz;
    }
    r1d[0][k] = r1;
    r1d[1][k] = r2;
    r1d[2][k] = r3;
  }
}

/* ----------------------------------------------------------------------
   charge assignment into drho1d
   dx,dy,dz = distance of particle from "lower left" grid point
------------------------------------------------------------------------- */

void PPPMDisp::compute_drho1d(const FFT_SCALAR &dx, const FFT_SCALAR &dy,
                          const FFT_SCALAR &dz, int ord, 
                              FFT_SCALAR **drho_c, FFT_SCALAR **dr1d)
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  for (k = (1-ord)/2; k <= ord/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = ord-2; l >= 0; l--) {
      r1 = drho_c[l][k] + r1*dx;
      r2 = drho_c[l][k] + r2*dy;
      r3 = drho_c[l][k] + r3*dz;
    }
    dr1d[0][k] = r1;
    dr1d[1][k] = r2;
    dr1d[2][k] = r3;
  }
}

/* ----------------------------------------------------------------------
   generate coeffients for the weight function of order n

              (n-1)
  Wn(x) =     Sum    wn(k,x) , Sum is over every other integer
           k=-(n-1)
  For k=-(n-1),-(n-1)+2, ....., (n-1)-2,n-1
      k is odd integers if n is even and even integers if n is odd
              ---
             | n-1
             | Sum a(l,j)*(x-k/2)**l   if abs(x-k/2) < 1/2
  wn(k,x) = <  l=0
             |
             |  0                       otherwise
              ---
  a coeffients are packed into the array rho_coeff to eliminate zeros
  rho_coeff(l,((k+mod(n+1,2))/2) = a(l,k) 
------------------------------------------------------------------------- */

void PPPMDisp::compute_rho_coeff(FFT_SCALAR **coeff , FFT_SCALAR **dcoeff, 
                                 int ord)
{
  int j,k,l,m;
  FFT_SCALAR s;

  FFT_SCALAR **a;
  memory->create2d_offset(a,ord,-ord,ord,"pppm/disp:a");

  for (k = -ord; k <= ord; k++) 
    for (l = 0; l < ord; l++)
      a[l][k] = 0.0;
        
  a[0][0] = 1.0;
  for (j = 1; j < ord; j++) {
    for (k = -j; k <= j; k += 2) {
      s = 0.0;
      for (l = 0; l < j; l++) {
	a[l+1][k] = (a[l][k+1]-a[l][k-1]) / (l+1);
#ifdef FFT_SINGLE
	s += powf(0.5,(float) l+1) *
	  (a[l][k-1] + powf(-1.0,(float) l) * a[l][k+1]) / (l+1);
#else
	s += pow(0.5,(double) l+1) * 
	  (a[l][k-1] + pow(-1.0,(double) l) * a[l][k+1]) / (l+1);
#endif
      }
      a[0][k] = s;
    }
  }

  m = (1-ord)/2;
  for (k = -(ord-1); k < ord; k += 2) {
    for (l = 0; l < ord; l++)
      coeff[l][m] = a[l][k];
    for (l = 1; l < ord; l++)
      dcoeff[l-1][m] = l*a[l][k];
    m++;
  }

  memory->destroy2d_offset(a,-ord);
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
------------------------------------------------------------------------- */

void PPPMDisp::slabcorr(int eflag)
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  double zprd = domain->zprd;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i]*x[i][2];

  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant

  double dipole_r2 = 0.0;
  if (eflag_atom || fabs(qsum) > SMALL) {
    for (int i = 0; i < nlocal; i++)
      dipole_r2 += q[i]*x[i][2]*x[i][2];

    // sum local contributions

    double tmp;
    MPI_Allreduce(&dipole_r2,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    dipole_r2 = tmp;
  }

  // compute corrections

  const double e_slabcorr = MY_2PI*(dipole_all*dipole_all -
    qsum*dipole_r2 - qsum*qsum*zprd*zprd/12.0)/volume;
  const double qscale = force->qqrd2e * scale;

  if (eflag_global) energy_1 += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = qscale * MY_2PI/volume;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * q[i]*(x[i][2]*dipole_all - 0.5*(dipole_r2 +
        qsum*x[i][2]*x[i][2]) - qsum*zprd*zprd/12.0);
  }

  // add on force corrections

  double ffact = qscale * (-4.0*MY_PI/volume);
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) f[i][2] += ffact * q[i]*(dipole_all - qsum*x[i][2]);
}

/* ----------------------------------------------------------------------
   perform and time the 1d FFTs required for N timesteps
------------------------------------------------------------------------- */

int PPPMDisp::timing_1d(int n, double &time1d)
{
  double time1,time2;
  int mixing = 1;
  if (function[2]) mixing = 4;
  if (function[3]) mixing = nsplit_alloc/2;

  if (function[0]) for (int i = 0; i < 2*nfft_both; i++) work1[i] = ZEROF;
  if (function[1] + function[2] + function[3])
    for (int i = 0; i < 2*nfft_both_6; i++) work1_6[i] = ZEROF;

  MPI_Barrier(world);
  time1 = MPI_Wtime();

  if (function[0]) {
    for (int i = 0; i < n; i++) {
      fft1->timing1d(work1,nfft_both,1);
      fft2->timing1d(work1,nfft_both,-1);
      if (differentiation_flag != 1){
        fft2->timing1d(work1,nfft_both,-1);
        fft2->timing1d(work1,nfft_both,-1);
      }
    }
  }

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  time1d = time2 - time1;

  MPI_Barrier(world);
  time1 = MPI_Wtime();

  if (function[1] + function[2] + function[3]) {
    for (int i = 0; i < n; i++) {
      fft1_6->timing1d(work1_6,nfft_both_6,1);
      fft2_6->timing1d(work1_6,nfft_both_6,-1);
      if (differentiation_flag != 1){
        fft2_6->timing1d(work1_6,nfft_both_6,-1);
        fft2_6->timing1d(work1_6,nfft_both_6,-1);
      }
    }
  }

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  time1d += (time2 - time1)*mixing;

  if (differentiation_flag) return 2;
  return 4;
}

/* ----------------------------------------------------------------------
   perform and time the 3d FFTs required for N timesteps
------------------------------------------------------------------------- */

int PPPMDisp::timing_3d(int n, double &time3d)
{
  double time1,time2;
  int mixing = 1;
  if (function[2]) mixing = 4;
  if (function[3]) mixing = nsplit_alloc/2;

  if (function[0]) for (int i = 0; i < 2*nfft_both; i++) work1[i] = ZEROF;
  if (function[1] + function[2] + function[3]) 
    for (int i = 0; i < 2*nfft_both_6; i++) work1_6[i] = ZEROF;



  MPI_Barrier(world);
  time1 = MPI_Wtime();

  if (function[0]) {
    for (int i = 0; i < n; i++) {
      fft1->compute(work1,work1,1);
      fft2->compute(work1,work1,-1);
      if (differentiation_flag != 1) {
        fft2->compute(work1,work1,-1);
        fft2->compute(work1,work1,-1);
      }
    }
  }

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  time3d = time2 - time1;

  MPI_Barrier(world);
  time1 = MPI_Wtime();
  
  if (function[1] + function[2] + function[3]) {
    for (int i = 0; i < n; i++) {
      fft1_6->compute(work1_6,work1_6,1);
      fft2_6->compute(work1_6,work1_6,-1);
      if (differentiation_flag != 1) {
        fft2_6->compute(work1_6,work1_6,-1);
        fft2_6->compute(work1_6,work1_6,-1);
      }
    }
  }

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  time3d += (time2 - time1) * mixing;

  if (differentiation_flag) return 2;
  return 4;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays 
------------------------------------------------------------------------- */

double PPPMDisp::memory_usage()
{
  double bytes = nmax*3 * sizeof(double);
  int mixing = 1;
  int diff = 3;     //depends on differentiation
  int per = 7;      //depends on per atom calculations
  if (differentiation_flag) {
    diff = 1;
    per = 6;
  }
  if (!evflag_atom) per = 0;
  if (function[2]) mixing = 7;
  if (function[3]) mixing = nsplit_alloc;

  if (function[0]) {
    int nbrick = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) * 
      (nzhi_out-nzlo_out+1);
    bytes += (1 + diff +  per) * nbrick * sizeof(FFT_SCALAR);     //brick memory
    bytes += 6 * nfft_both * sizeof(double);      // vg
    bytes += nfft_both * sizeof(double);          // greensfn
    bytes += nfft_both * 3 * sizeof(FFT_SCALAR);    // density_FFT, work1, work2 
    bytes += cg->memory_usage();
  }

  if (function[1] + function[2] + function[3]) {
    int nbrick = (nxhi_out_6-nxlo_out_6+1) * (nyhi_out_6-nylo_out_6+1) * 
      (nzhi_out_6-nzlo_out_6+1);
    bytes += (1 + diff + per ) * nbrick * sizeof(FFT_SCALAR) * mixing;     // density_brick + vd_brick + per atom bricks
    bytes += 6 * nfft_both_6 * sizeof(double);      // vg
    bytes += nfft_both_6 * sizeof(double);          // greensfn
    bytes += nfft_both_6 * (mixing + 2) * sizeof(FFT_SCALAR);    // density_FFT, work1, work2 
    bytes += cg_6->memory_usage();
  }
  return bytes;
}
