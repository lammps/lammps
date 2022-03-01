/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Christian Tabedzki (U Penn) cat159@scarletmail.rutgers.edu
                         Andrew Santos (SNL) a.p.santos11@gmail.com
                         Zachariah Vicars (U Penn).
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------
    TILD - Theoretically Informed Langevan Dynamics
    This replicates the TILD code done by the Riggleman group, 
    previously known as Dynamical Mean Field Theory. 
-------------------------------------------------------------------- */

#include "tild.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "force.h"
#include "gridcomm.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include "remap_wrap.h"
#include "thermo.h"
#include "update.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <tuple>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;

#define SMALL 0.00001
#define OFFSET 16384
#define MAXORDER   7

enum{REVERSE_RHO_NONE};
enum{FORWARD_NONE, FORWARD_GRID_DEN, FORWARD_AVG_GRID_DEN};

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif


/* ------------------------------------------------------------ */

TILD::TILD(LAMMPS *lmp) : KSpace(lmp),
  factors(nullptr), fkx(nullptr), fky(nullptr), fkz(nullptr), vg(nullptr),
  work1(nullptr), work2(nullptr), rho1d(nullptr), rho_coeff(nullptr), drho_coeff(nullptr), 
  gradWtypex(nullptr), gradWtypey(nullptr), gradWtypez(nullptr), density_brick_types(nullptr),
  fft1(nullptr), fft2(nullptr), remap(nullptr), gc(nullptr), 
  gc_buf1(nullptr), gc_buf2(nullptr),
  part2grid(nullptr), boxlo(nullptr)
 {
  peratom_allocate_flag = 0;
  group_allocate_flag = 0;
  
  nstyles = 3;

  pppmflag = 0;
  group_group_enable = 0;
  tildflag = 1;
  triclinic = domain->triclinic;
  triclinic_support = 0;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  nfft_both = 0;
  nxhi_in = nxlo_in = nxhi_out = nxlo_out = 0;
  nyhi_in = nylo_in = nyhi_out = nylo_out = 0;
  nzhi_in = nzlo_in = nzhi_out = nzlo_out = 0;

  nfactors = 3;
  factors = new int[nfactors];
  factors[0] = 2;
  factors[1] = 3;
  factors[2] = 5;

  fkx = fky = fkz = nullptr;
  fkx2 = fky2 = fkz2 = nullptr;

  vg = nullptr;
  vg_hat = nullptr;

  work1 = work2 = nullptr;
  rho1d = rho_coeff = drho_coeff = nullptr;

  fft1 = fft2 = nullptr;
  remap = nullptr;
  gc = nullptr;
  gc_buf1 = gc_buf2 = nullptr;

  part2grid = nullptr;

  density_brick_types = nullptr;
  avg_density_brick_types = nullptr;
  density_fft_types = nullptr;
  density_hat_fft_types = nullptr;
  gradWtypex = nullptr;
  gradWtypey = nullptr;
  gradWtypez = nullptr;
  potent = nullptr;
  potent_hat = nullptr;
  grad_potent = nullptr;
  grad_potent_hat = nullptr;
  setflag = nullptr;
  potent_type_map = nullptr;
  chi = nullptr;
  a2 = nullptr;
  rp = nullptr;
  xi = nullptr;
  ktmp = nullptr;
  ktmp2 = nullptr;
  ktmpi = ktmpj = nullptr;
  ktmp2i = ktmp2j = nullptr;
  tmp = nullptr;

  rho0 = 0.0;
  nmax = 0;
  sub_flag  = 1;
  mix_flag  = 1;
  ave_grid_flag  = 0;
  nvalid_last = -1;
  nvalid = 0;
  nevery = 0;
  irepeat = 0;
  nrepeat = 0;
  peratom_freq = 0;
  write_grid_flag  = 0;
  set_rho0 = 0.0;
  subtract_rho0 = 0;
  norm_flag = 1;
  npergrid = 0; 

  int ntypes = atom->ntypes;
  memory->create(potent_type_map,nstyles+1,ntypes+1,ntypes+1,"tild:potent_type_map");  
  memory->create(density_flags,ntypes+1,"tild:density_flags");  // determines if to pass data types or not

  memory->create(chi,ntypes+1,ntypes+1,"tild:chi"); // interaction parameters
  memory->create(a2,ntypes+1,ntypes+1,"tild:a2"); // gaussian parameter
  memory->create(xi,ntypes+1,ntypes+1,"tild:xi"); // erfc parameter
  memory->create(rp,ntypes+1,ntypes+1,"tild:rp"); // erfc parameter
  memory->create(np_rho,ntypes+1,ntypes+1,"tild:np_rho"); // erfc parameter

  for (int i = 0; i <= ntypes; i++) {
    for (int j = 0; j <= ntypes; j++) {
      potent_type_map[0][i][j] = 1; // style 0 is 1 if no tild potential is used
      for (int k = 1; k <= nstyles; k++) {
        potent_type_map[k][i][j] = 0; // style is set to 1 if it is used by type-type interaction
      }
      chi[i][j] = 0.0;
      a2[i][j] = 0;
      xi[i][j] = 0;
      rp[i][j] = 0;
      np_rho[i][j] = 0;
    }
    density_flags[i] = 0;
  }
};

/* ---------------------------------------------------------------------- */

void TILD::settings(int narg, char **arg)
{
    if (narg < 1) error->all(FLERR,"Illegal kspace_style tild command");
    grid_res = fabs(utils::numeric(FLERR,arg[0],false,lmp));
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

TILD::~TILD(){

  delete [] factors;
  TILD::deallocate();
  memory->destroy(part2grid);
  part2grid = nullptr;

  // check whether cutoff and pair style are set

  triclinic = domain->triclinic;
  pair_check();
  memory->destroy(potent_type_map);
  memory->destroy(density_flags);
  memory->destroy(chi);
  memory->destroy(a2);
  memory->destroy(rp);
  memory->destroy(xi);

  potent_type_map = nullptr;
  density_flags = nullptr;
  chi = nullptr;
  a2 = nullptr;
  rp = nullptr;
  xi = nullptr;
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void TILD::init()
{
  if (me == 0) utils::logmesg(lmp,"TILD initialization...\n");

  triclinic_check();
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use TILD with 2d simulation");
  if (differentiation_flag != 0)
    error->all(FLERR, "Cannot use analytic differentiation with TILD");
  if (comm->style != 0)
    error->universe_all(FLERR,"TILD can only currently be used with "
                        "comm_style brick");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with TILD");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab TILD");
  }

  if (order < 2 || order > MAXORDER)
    error->all(FLERR,"TILD order cannot be < 2 or > {}",MAXORDER);


  triclinic = domain->triclinic;
  pair_check();

  // free all arrays previously allocated
  deallocate();

  // setup FFT grid resolution and g_ewald
  // normally one iteration thru while loop is all that is required
  // if grid stencil does not extend beyond neighbor proc
  //   or overlap is allowed, then done
  // else reduce order and try again

  GridComm *gctmp = nullptr;
  int iteration = 0;

  while (order >= minorder) {
    if (iteration && me == 0)
      error->warning(FLERR,"Reducing TILD order b/c stencil extends "
                     "beyond nearest neighbor processor");

    // if (stagger_flag && !differentiation_flag) compute_gf_denom();
    set_grid_global();
    set_grid_local();
    if (overlap_allowed) break;

    gctmp = new GridComm(lmp,world,nx_tild,ny_tild,nz_tild,
                         nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                         nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out);

    int tmp1,tmp2;
    gctmp->setup(tmp1,tmp2);
    if (gctmp->ghost_adjacent()) break;
    delete gctmp;

    order--;
    iteration++;
  }

  if (order < minorder) error->all(FLERR,"TILD order < minimum allowed order");
  if (!overlap_allowed && !gctmp->ghost_adjacent())
    error->all(FLERR,"TILD grid stencil extends "
               "beyond nearest neighbor processor");
  if (gctmp) delete gctmp;

  set_grid_global();
  set_grid_local();

  // print stats

  int ngrid_max,nfft_both_max;
  MPI_Allreduce(&ngrid,&ngrid_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nfft_both,&nfft_both_max,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {
    std::string mesg = fmt::format("  grid = {} {} {}\n",nx_tild,ny_tild,nz_tild);
    mesg += fmt::format("  stencil order = {}\n",order);
    mesg += fmt::format("  using " LMP_FFT_PREC " precision " LMP_FFT_LIB "\n");
    mesg += fmt::format("  3d grid and FFT values/proc = {} {}\n", ngrid_max, nfft_both_max);
    utils::logmesg(lmp,mesg);
  }

  // allocate k-space dependent memory

  allocate();
  rho0 = calculate_rho0();

  double volume = domain->xprd * domain->yprd * domain->zprd;

  compute_rho_coeff();
}

/* ----------------------------------------------------------------------
   adjust TILD coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void TILD::setup(){

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with TILD");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab TILD");
  }


  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab TILD
  // z dimension for 3d TILD is zprd since slab_volfactor = 1.0

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

  delxinv = nx_tild/xprd;
  delyinv = ny_tild/yprd;
  delzinv = nz_tild/zprd_slab;

  double per;
    int i, j;

    for (i = nxlo_fft; i <= nxhi_fft; i++) {
      per = i - nx_tild*(2*i/nx_tild);
      fkx[i] = unitkx*per;
      j = (nx_tild - i) % nx_tild;
      per = j - nx_tild*(2*j/nx_tild);
      fkx2[i] = unitkx*per;
    }

    for (i = nylo_fft; i <= nyhi_fft; i++) {
      per = i - ny_tild*(2*i/ny_tild);
      fky[i] = unitky*per;
      j = (ny_tild - i) % ny_tild;
      per = j - ny_tild*(2*j/ny_tild);
      fky2[i] = unitky*per;
    }

    for (i = nzlo_fft; i <= nzhi_fft; i++) {
      per = i - nz_tild*(2*i/nz_tild);
      fkz[i] = unitkz*per;
      j = (nz_tild - i) % nz_tild;
      per = j - nz_tild*(2*j/nz_tild);
      fkz2[i] = unitkz*per;
    }

  delvolinv = delxinv*delyinv*delzinv;

  if (sub_flag == 1) {
    subtract_rho0 = 1;
  } else {
    subtract_rho0 = 0;
  }  
  if (norm_flag == 1) {
    normalize_by_rho0 = 1;
  } else {
    normalize_by_rho0 = 0;
  }

  if (mix_flag == 1) {
    int ntypes = atom->ntypes;
    for (int itype = 1; itype <= ntypes; itype++) {
      for (int jtype = itype+1; jtype <= ntypes; jtype++) {
        if (potent_type_map[0][itype][itype] ==1 || potent_type_map[0][jtype][jtype] == 1 ) {
          potent_type_map[0][itype][jtype] = 1;
          for (int istyle = 1; istyle <= nstyles; istyle++) 
            potent_type_map[istyle][itype][jtype] = potent_type_map[istyle][itype][itype];
        } else {
          potent_type_map[0][itype][jtype] = 0;
          for (int istyle = 1; istyle <= nstyles; istyle++) 
            // assume it is of type istyle, but it does the convolution
            potent_type_map[istyle][itype][jtype] = -1; 
        }
      }
    }
  }

  manually_flip_density_flags();
  rho0 = calculate_rho0();
  init_cross_potentials();
  vir_func_init();

  return;
}

/* ----------------------------------------------------------------------
   Initialization of the cross-potential/virial functions
------------------------------------------------------------------------- */

void TILD::vir_func_init() {
  int z, y, x, n;
  int Dim = domain->dimension;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  //double zprd_slab = zprd * slab_volfactor;
  double k[Dim];
  double scale_inv = 1.0 / (nx_tild * ny_tild * nz_tild);
  double delx = xprd/nx_tild;
  double dely = yprd/ny_tild;
  double delz = zprd/nz_tild;
  int ntypes = atom->ntypes;

  // Initializating for shape defined interactions
  int loc = -1;
  for (int itype = 1; itype <= ntypes; itype++) {
    for (int jtype = itype; jtype <= ntypes; jtype++) {
      // Skip if type cross-interaction does not use density/tild
      loc++;
      if ( potent_type_map[0][itype][jtype] == 1) continue;

      n = 0;
      for (z = nzlo_fft; z <= nzhi_fft; z++) {
        // PBC
        if (double(z) < double(nz_tild) / 2.0)
          k[2] = double(z) * delz;
        else
          k[2] = -double(nz_tild - z) * delz;

        for (y = nylo_fft; y <= nyhi_fft; y++) {
          if (double(y) < double(ny_tild) / 2.0)
            k[1] = double(y) * dely;
          else
            k[1] = -double(ny_tild - y) * dely;

          for (x = nxlo_fft; x <= nxhi_fft; x++) {
            if (double(x) < double(nx_tild) / 2.0)
              k[0] = double(x) * delx;
            else
              k[0] = -double(nx_tild - x) * delx;

            vg[loc][0][n] = k[0] * -grad_potent[loc][0][n];
            vg[loc][1][n] = k[1] * -grad_potent[loc][1][n];
            vg[loc][2][n] = k[2] * -grad_potent[loc][2][n];
            vg[loc][3][n] = k[1] * -grad_potent[loc][0][n];
            vg[loc][4][n] = k[2] * -grad_potent[loc][0][n];
            vg[loc][5][n] = k[2] * -grad_potent[loc][1][n];

            n++;

            }
          }
      }

      for (int i = 0; i < 6; i++){
        n = 0;
        for (int j = 0; j < nfft; j++) {
          ktmp[n++] = vg[loc][i][j];
          ktmp[n++] = ZEROF;
        }
  
        fft1->compute(ktmp, ktmp2, FFT3d::FORWARD);

        for (int j = 0; j < 2 * nfft; j++) {
          ktmp2[j] *= scale_inv;
          vg_hat[loc][i][j] = ktmp2[j];
        }
      }
    }
  }
  
  // Initializating for user-defined cross-interactions
  for (auto& potential: cross_iter){
    int i_index = min(potential.i,potential.j);
    int j_index = max(potential.i,potential.j);
    loc = (i_index*(2*ntypes - i_index - 1))/2 + j_index;

    if (potential.type != NONE){
      n = 0;
      for (z = nzlo_fft; z <= nzhi_fft; z++) {
        // PBC
        if (double(z) < double(nz_tild) / 2.0)
          k[2] = double(z) * delz;
        else
          k[2] = -double(nz_tild - z) * delz;

        for (y = nylo_fft; y <= nyhi_fft; y++) {
          if (double(y) < double(ny_tild) / 2.0)
            k[1] = double(y) * dely;
          else
            k[1] = -double(ny_tild - y) * dely;

          for (x = nxlo_fft; x <= nxhi_fft; x++) {
            if (double(x) < double(nx_tild) / 2.0)
              k[0] = double(x) * delx;
            else
              k[0] = -double(nx_tild - x) * delx;

            vg[loc][0][n] = k[0] * -grad_potent[loc][0][n];
            vg[loc][1][n] = k[1] * -grad_potent[loc][1][n];
            vg[loc][2][n] = k[2] * -grad_potent[loc][2][n];
            vg[loc][3][n] = k[1] * -grad_potent[loc][0][n];
            vg[loc][4][n] = k[2] * -grad_potent[loc][0][n];
            vg[loc][5][n] = k[2] * -grad_potent[loc][1][n];

            n++;

    }
  }
}

      for (int i = 0; i < 6; i++){
        n = 0;
        for (int j = 0; j < nfft; j++) {
          ktmp[n++] = vg[loc][i][j];
          ktmp[n++] = ZEROF;
        }
  
        fft1->compute(ktmp, ktmp2, FFT3d::FORWARD);

        for (int j = 0; j < 2 * nfft; j++) {
          ktmp2[j] *= scale_inv;
          vg_hat[loc][i][j] = ktmp2[j];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reset local grid arrays and communication stencils
   called by fix balance b/c it changed sizes of processor sub-domains
------------------------------------------------------------------------- */

void TILD::setup_grid()
{
  // free all arrays previously allocated

  deallocate();

  // reset portion of global grid that each proc owns

  set_grid_local();

  // reallocate K-space dependent memory
  // check if grid communication is now overlapping if not allowed

  allocate();

  if (!overlap_allowed && !gc->ghost_adjacent())
    error->all(FLERR,"TILD grid stencil extends "
               "beyond nearest neighbor processor");

  // pre-compute 1d density distribution coefficients
  compute_rho_coeff();

  // pre-compute volume-dependent coeffs for portion of grid I now own
  setup();
}

/* ----------------------------------------------------------------------
   precompute the fft of the gridded densities
------------------------------------------------------------------------- */

void TILD::precompute_density_hat_fft() {
  int n = 0;
  double scale_inv = 1.0 / (nx_tild * ny_tild * nz_tild);
  int ntypes = atom->ntypes;

  for ( int ktype = 1; ktype <= ntypes; ktype++) {
    if (density_flags[ktype] != 1) continue;
    n = 0;
    for (int k = 0; k < nfft; k++) {
      work1[n++] = density_fft_types[ktype][k];
      work1[n++] = ZEROF;
    }

    // FFT the density to k-space
    fft1->compute(work1, work1, FFT3d::FORWARD);

    for (int k = 0; k < 2*nfft; k++) {
      work1[k] *= scale_inv;
      density_hat_fft_types[ktype][k] = work1[k];
    }
  }
}

/* ----------------------------------------------------------------------
   compute the forces and energies of the system
------------------------------------------------------------------------- */

void TILD::compute(int eflag, int vflag){
  
  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  ev_init(eflag,vflag);

  if (atom->nmax > nmax) {

    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm/disp:part2grid");
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid

  particle_map();
  make_rho();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  int ntypes = atom->ntypes;
  gc->reverse_comm(GridComm::KSPACE,this,(ntypes+1)*3,sizeof(FFT_SCALAR),REVERSE_RHO_NONE,
                          gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  brick2fft();

  // compute potential gradient on my FFT grid 
  // returns gradients in 3d brick decomposition 
  // Question about whether it should perform per_atom calculations 

  accumulate_gradient();

  // all procs communicate gradient potential values
  // to fill ghost cells surrounding their 3d bricks
  
  gc->forward_comm(GridComm::KSPACE,this,(ntypes+1)*3,sizeof(FFT_SCALAR),FORWARD_NONE,
                          gc_buf1,gc_buf2,MPI_FFT_SCALAR);


  // PPPM has  extra per-atom energy/virial communication here. 
  // Again, I don't know if we should do per-atom calculations.
  // if so, the head LAMMPS developers will be able to tell us

  // calculate the force on my particles

  fieldforce_param();

  // check whether grid densities should be written out at this time 
  // Serious question about whether this should be a fix instead of this 
  // Answer: This should be a fix. That being said, currently easier to implement here.

  if ( write_grid_flag == 1 ) {
    if ( (update->ntimestep % grid_data_output_freq) ==  0) {
      write_grid_data(grid_data_filename, 0);
    }
  }
  if ( ave_grid_flag == 1 ) ave_grid();

  // sum global energy across procs and add in volume-dependent term

  if (eflag_global){
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all; 
  }

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    if (normalize_by_rho0 == 1)
    for (int i = 0; i < 6; i++) {
      // change number density to tild density for this contribution to virial
      virial[i] = virial_all[i] / rho0 * atom->natoms; 
    }
  }

}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double TILD::memory_usage()
{
  int ntypes = atom->ntypes;
  int ntypecross = ((ntypes + 1)*ntypes);
  int Dim = domain->dimension;
  double bytes = (ntypes + 1) * (ntypes + 1) * (nstyles + 1) * sizeof(int); // potent_type_map
  bytes += (ntypes + 1) * (nstyles + 1) * sizeof(int); // density_flags

  // parameters and chis of the sim, may need to be changed if new potential are added
  bytes += (double) (ntypes + 1) * (ntypes + 1) * 4 * sizeof(double); 
  bytes += (double) (ntypes + 1) * (ntypes + 1) * sizeof(int); // setflag

  int nbrick = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) *
    (nzhi_out - nzlo_out + 1); // brick size

  bytes += (double) (ntypes + 1) * 5 * nbrick * sizeof(FFT_SCALAR); // density brick types

  bytes += (double) (ntypecross + 1) * 6 * (2+1) * nfft_both * sizeof(FFT_SCALAR); // vg vg_hat
  bytes += (double) (ntypes + 1) * (2+1) * nfft_both * sizeof(FFT_SCALAR); // density_fft_types
  bytes += (double) (ntypecross+1) * (2+1) * sizeof(FFT_SCALAR); // potent and hat
  bytes += (double) Dim * (ntypecross+1) * (2+1) * sizeof(FFT_SCALAR); // grad potent and hat
  
  // fkx, fky, fkz, fkx2, fky2, fkz2 
  bytes += (double) 6 * npergrid * sizeof(FFT_SCALAR);
  // rho1d, rho_coeff, drho_coeff
  bytes += (double) 3 * (order+1) * sizeof(FFT_SCALAR);

  bytes += (double) nfft_both * 2 * 8 * sizeof(FFT_SCALAR); // all works, ktmps
  bytes += (double) nfft * sizeof(FFT_SCALAR); // just the tmp

  // two GridComm bufs 
  bytes += (double)(ngc_buf1 + ngc_buf2) * (ntypes + 1) * npergrid * sizeof(FFT_SCALAR);

  return bytes;
}


/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void TILD::allocate()
{

  int (*procneigh)[2] = comm->procneigh;
  
  int ntypes = atom->ntypes;
  int ntypecross = ((ntypes+1)*ntypes);

  // style coeffs
  memory->create(setflag,ntypes+1,ntypes+1,"pair:setflag");
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++)
      setflag[i][j] = 0;


  memory->create4d_offset(density_brick_types,ntypes+1,
                          nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"tild:density_brick_types");
  memory->create4d_offset(avg_density_brick_types,ntypes+1,
                          nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"tild:avg_density_brick_types");
  memory->create4d_offset(gradWtypex,ntypes+1, 
                          nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"tild:gradWtypex");
  memory->create4d_offset(gradWtypey,ntypes+1, 
                          nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"tild:gradWtypey");
  memory->create4d_offset(gradWtypez,ntypes+1, 
                          nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"tild:gradWtypez");

  memory->create(work1,2*nfft_both,"tild:work1");
  memory->create(work2,2*nfft_both,"tild:work2");
  memory->create(ktmp,2*nfft_both,"tild:ktmp");
  memory->create(ktmpi,2*nfft_both,"tild:ktmpi");
  memory->create(ktmpj,2*nfft_both,"tild:ktmpj");
  memory->create(ktmp2,2*nfft_both,"tild:ktmp2");
  memory->create(ktmp2i,2*nfft_both,"tild:ktmp2i");
  memory->create(ktmp2j,2*nfft_both,"tild:ktmp2j");
  memory->create(tmp,nfft,"tild:tmp");
  memory->create(vg,ntypecross+1,6,nfft_both,"tild:vg");
  memory->create(vg_hat,ntypecross+1,6,2*nfft_both,"tild:vg_hat");
  memory->create(density_fft_types,ntypes+1,nfft_both,"tild:density_fft_types");
  memory->create(density_hat_fft_types,ntypes+1,2*nfft_both,"tild:density_hat_fft_types");
  memory->create(potent,ntypecross+1,nfft_both,"tild:potent"); // Voignot 
  memory->create(potent_hat,ntypecross+1,2*nfft_both,"tild:potent_hat");
  memory->create(grad_potent,ntypecross+1,domain->dimension,nfft_both,"tild:grad_potent");
  memory->create(grad_potent_hat,ntypecross+1,domain->dimension,2*nfft_both,"tild:grad_potent_hat");

  if (triclinic == 0) {
    memory->create1d_offset(fkx,nxlo_fft,nxhi_fft,"tild:fkx");
    memory->create1d_offset(fky,nylo_fft,nyhi_fft,"tild:fky");
    memory->create1d_offset(fkz,nzlo_fft,nzhi_fft,"tild:fkz");
    memory->create1d_offset(fkx2,nxlo_fft,nxhi_fft,"tild:fkx2");
    memory->create1d_offset(fky2,nylo_fft,nyhi_fft,"tild:fky2");
    memory->create1d_offset(fkz2,nzlo_fft,nzhi_fft,"tild:fkz2");
  } else {
    memory->create(fkx,nfft_both,"tild:fkx");
    memory->create(fky,nfft_both,"tild:fky");
    memory->create(fkz,nfft_both,"tild:fkz");
    memory->create(fkx2,nfft_both,"tild:fkx2");
    memory->create(fky2,nfft_both,"tild:fky2");
    memory->create(fkz2,nfft_both,"tild:fkz2");
  }

  // summation coeffs

  memory->create2d_offset(rho1d,3,-order/2,order/2,"tild:rho1d");
  memory->create2d_offset(rho_coeff,order,(1-order)/2,order/2,"tild:rho_coeff");
  memory->create2d_offset(drho_coeff,order,(1-order)/2,order/2,"tild:drho_coeff");

  // create 2 FFTs and a Remap
  // 1st FFT keeps data in FFT decompostion
  // 2nd FFT returns data in 3d brick decomposition
  // remap takes data from 3d brick to FFT decomposition

  int tmp;

  fft1 = new FFT3d(lmp,world,nx_tild,ny_tild,nz_tild,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   0,0,&tmp,collective_flag);

  fft2 = new FFT3d(lmp,world,nx_tild,ny_tild,nz_tild,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                   0,0,&tmp,collective_flag);

  remap = new Remap(lmp,world,
                    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                    nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                    1,0,0,FFT_PRECISION,collective_flag);

  // create ghost grid object for rho and electric field communication
  // also create 2 bufs for ghost grid cell comm, passed to GridComm methods

  gc = new GridComm(lmp,world,nx_tild,ny_tild,nz_tild,
                    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                    nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out);

  gc->setup(ngc_buf1,ngc_buf2);

  if (differentiation_flag) npergrid = 1;
  else npergrid = 3;

  memory->create(gc_buf1,(ntypes+1) * npergrid*ngc_buf1,"tild:gc_buf1");
  memory->create(gc_buf2,(ntypes+1) * npergrid*ngc_buf2,"tild:gc_buf2");
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void TILD::deallocate()
{
  memory->destroy(setflag);
  setflag = nullptr;

  memory->destroy4d_offset(density_brick_types,nzlo_out,nylo_out,nxlo_out);
  memory->destroy4d_offset(avg_density_brick_types,nzlo_out,nylo_out,nxlo_out);
  memory->destroy(density_fft_types);
  memory->destroy(density_hat_fft_types);
  density_brick_types = nullptr;
  avg_density_brick_types = nullptr;
  density_fft_types = density_hat_fft_types = nullptr;
  memory->destroy(ktmp);
  memory->destroy(ktmpi);
  memory->destroy(ktmpj);
  memory->destroy(ktmp2);
  memory->destroy(ktmp2i);
  memory->destroy(ktmp2j);
  memory->destroy(tmp);
  tmp = ktmp = ktmpi = ktmpj = ktmp2 = ktmp2i = ktmp2j = tmp = nullptr;

  memory->destroy(vg);
  memory->destroy(vg_hat);
  memory->destroy(potent);
  memory->destroy(potent_hat);
  memory->destroy(grad_potent);
  memory->destroy(grad_potent_hat);
  vg_hat = vg = nullptr;
  potent_hat = potent = nullptr;
  grad_potent_hat = grad_potent = nullptr;

  memory->destroy4d_offset(gradWtypex, nzlo_out,nylo_out, nxlo_out);
  memory->destroy4d_offset(gradWtypey, nzlo_out,nylo_out, nxlo_out);
  memory->destroy4d_offset(gradWtypez, nzlo_out,nylo_out, nxlo_out);

  gradWtypex = nullptr;
  gradWtypey = nullptr;
  gradWtypez = nullptr;
  memory->destroy(work1);
  memory->destroy(work2);
  work1 = work2 = nullptr;

  if (triclinic == 0) {
    memory->destroy1d_offset(fkx,nxlo_fft);
    memory->destroy1d_offset(fky,nylo_fft);
    memory->destroy1d_offset(fkz,nzlo_fft);
    memory->destroy1d_offset(fkx2,nxlo_fft);
    memory->destroy1d_offset(fky2,nylo_fft);
    memory->destroy1d_offset(fkz2,nzlo_fft);
  } else {
    memory->destroy(fkx);
    memory->destroy(fky);
    memory->destroy(fkz);
    memory->destroy(fkx2);
    memory->destroy(fky2);
    memory->destroy(fkz2);
  }

  int order_allocated = order;
  memory->destroy2d_offset(rho1d,-order_allocated/2);
  memory->destroy2d_offset(rho_coeff,(1-order_allocated)/2);
  memory->destroy2d_offset(drho_coeff,(1-order_allocated)/2);
  rho1d = rho_coeff = drho_coeff = nullptr;

  delete fft1;
  delete fft2;
  delete remap;
  delete gc;
  memory->destroy(gc_buf1);
  memory->destroy(gc_buf2);
  fft1 = fft2 = nullptr;
  remap = nullptr;
}

/* ----------------------------------------------------------------------
   set size of FFT grid (nx,ny,nz_tild) and g_ewald
   for TILD interactions
------------------------------------------------------------------------- */

void TILD::set_grid_global()
{
  // use xprd,yprd,zprd even if triclinic so grid size is the same
  // adjust z dimension for 2d slab TILD
  // 3d PPPM just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;


  double h, h_x,h_y,h_z;

  // set optimal nx_tild,ny_tild,nz_tild based on order and accuracy
  // nz_tild uses extended zprd_slab instead of zprd
  // reduce it until accuracy target is met

  if (!gridflag) {
    h = h_x = h_y = h_z = grid_res;

    // set grid dimension
    nx_tild = static_cast<int> (xprd/h_x);
    ny_tild = static_cast<int> (yprd/h_y);
    nz_tild = static_cast<int> (zprd_slab/h_z);

    if (nx_tild <= 1) nx_tild = 2;
    if (ny_tild <= 1) ny_tild = 2;
    if (nz_tild <= 1) nz_tild = 2;
  } else {
    nx_tild = nx_pppm;
    ny_tild = ny_pppm;
    nz_tild = nz_pppm;
  }

  if (triclinic) {
    double tmp[3];
    tmp[0] = nx_tild / xprd;
    tmp[1] = ny_tild / yprd;
    tmp[2] = nz_tild / zprd;
    lamda2xT(&tmp[0], &tmp[0]);
    nx_tild = static_cast<int>(tmp[0]) + 1;
    ny_tild = static_cast<int>(tmp[1]) + 1;
    nz_tild = static_cast<int>(tmp[2]) + 1;
  }

  // boost grid size until it is factorable

  while (!factorable(nx_tild)) nx_tild++;
  while (!factorable(ny_tild)) ny_tild++;
  while (!factorable(nz_tild)) nz_tild++;

  if (nx_tild >= OFFSET || ny_tild >= OFFSET || nz_tild >= OFFSET)
    error->all(FLERR, "TILD grid is too large");
}

/* ----------------------------------------------------------------------
   Initialize in real space the cross-potentials
------------------------------------------------------------------------- */

void TILD::init_cross_potentials(){
  

  int Dim = domain->dimension;
  int ntypes = atom->ntypes;
  double scale_inv = 1.0/ nx_tild/ ny_tild/ nz_tild;
  int n = 0;

  // Loop over potental styles
  int loc = -1;
  for (int itype = 1; itype <= ntypes; itype++) {
    // Skip if type cross-interaction does not use density/tild

    for (int jtype = itype; jtype <= ntypes; jtype++) {
      loc++;
      // Skip if type cross-interaction does not use density/tild
      if ( potent_type_map[0][itype][jtype] == 1) continue;


      // If both parameters are Gaussian, just do analytical convolution
      if (potent_type_map[1][itype][jtype] == 1 or (mix_flag == 1 && potent_type_map[1][itype][itype] == 1 && potent_type_map[1][jtype][jtype] == 1) ){
        // mixing
        double a2_mix;
        if (mix_flag == 1) {
          a2_mix = a2[itype][itype] + a2[jtype][jtype];
        } else {
          a2_mix = a2[itype][jtype] + a2[itype][jtype];
        }
        const double* p = &a2_mix;
        init_potential(potent[loc], 1, p);

        init_potential_ft(potent_hat[loc], 1, p);

      } 
      // Computational Convolution
      else {
        if (mix_flag == 1) {
          calc_work(work1, itype, itype);
          if ( itype == jtype) {
            for (int i = 0; i < 2*nfft; i++) work2[i] = work1[i];
          } else {
            calc_work(work2, jtype, jtype);
          }
        } else {
          calc_work(work1, itype, jtype);
          for (int i = 0; i < 2*nfft; i++) work2[i] = work1[i];
        }

        // Convolution of potentials
        n = 0;
        for (int i = 0; i < nfft; i++) {
          complex_multiply(work1, work2, ktmp2, n);
          potent_hat[loc][n] = ktmp2[n];
          potent_hat[loc][n+1] = ktmp2[n+1];
          n += 2;
        }

        // Inverse fourier transform
        fft1->compute(ktmp2, ktmp, FFT3d::BACKWARD);

        // Storing the real space result
        n = 0;
        for (int j = 0; j < nfft; j++){
          potent[loc][j] = ktmp[n];
          n += 2;
        }
      }

      get_k_alias(potent_hat[loc], grad_potent_hat[loc]);
      for (int i=0; i < Dim; i++){
        for (int j = 0; j < 2*nfft; j++) {
          work1[j] = grad_potent_hat[loc][i][j];
        }
        //fft2->compute(work1, work2, -1);
        fft1->compute(work1, work2, FFT3d::BACKWARD);
        n = 0;
        for (int j = 0; j < nfft; j++) {
          grad_potent[loc][i][j] = -work2[n]; // still not sure about minus sign ..
          n+=2;
        }
      } 

    }
  }

  // Implementing the manually specified interactions

  for (auto& potential: cross_iter){
    int i_index = min(potential.i,potential.j);
    int j_index = max(potential.i,potential.j);
    loc = (i_index*(2*ntypes - i_index - 1))/2 + j_index;

    calc_cross_work(potential);

    if (potential.type != NONE){
      get_k_alias(potent_hat[loc], grad_potent_hat[loc]);
      for (int i = 0; i < Dim; i++) {
        for (int j = 0; j < 2 * nfft; j++) { 
          work1[j] = grad_potent_hat[loc][i][j]; 
        }

        fft1->compute(work1, work2, FFT3d::BACKWARD);
        n = 0;
        for (int j = 0; j < nfft; j++) {
          grad_potent[loc][i][j] = -work2[n];
          n += 2;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Determine which cross-potential should be used
------------------------------------------------------------------------- */

int TILD::get_style( const int i, const int j) {
  for (int istyle = 1; istyle <= nstyles; istyle++) { 
    if ( potent_type_map[istyle][i][j] == 1 ) return istyle;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   Larger cross potential startup functoin 
------------------------------------------------------------------------- */

void TILD::calc_work(FFT_SCALAR *wk, const int itype, const int jtype){
  double scale_inv = 1.0/ nx_tild/ ny_tild/ nz_tild;

  // needs work is this right for cross terms of the same potential?
  
  double params[4];
  for (int i = 0; i < 4; i++) params[i] = 0.0;

  int style = get_style(itype, jtype);
  if (style == 1) {
    params[0] = a2[itype][jtype];
    init_potential_ft(wk, style, params);
  } else if (style == 2) {
    params[0] = rp[itype][jtype];
    params[1] = xi[itype][jtype];
    params[2] = np_rho[itype][jtype];
    if (set_rho0 > 0.0f) params[2] = modified_rho0;
    init_potential(tmp, style, params);

    int j = 0;
    for (int i = 0; i < nfft; i++) {
      ktmp[j++] = tmp[i];
      ktmp[j++] = ZEROF;
    }

    fft1->compute(ktmp, wk, FFT3d::FORWARD);

    for (int i = 0; i < 2 * nfft; i++) {
      wk[i] *= scale_inv;
    }
  }
}

/* ----------------------------------------------------------------------
   Calculate the cross_work for the manually specified interactions
------------------------------------------------------------------------- */

void TILD::calc_cross_work(const Interaction& intrxn){
  int ntypes = atom->ntypes;
  int Dim = domain->dimension;
  int n = 0;

  int i_index = min(intrxn.i,intrxn.j);
  int j_index = max(intrxn.i,intrxn.j);

  const int loc = (i_index*(2*ntypes - i_index - 1))/2 + j_index;
  if (intrxn.type == DELETE) error->all(FLERR, "DELETE for calc_cross_work is invalid and shouldn't exist.");


  if (intrxn.type == NONE){
    
    
    for (int i = 0; i < Dim; i++) {
      n = 0;
      for (int j = 0; j < nfft; j++) {
        grad_potent[loc][i][j] = 0.0f;
        grad_potent_hat[loc][i][n] = 0.0f;
        grad_potent_hat[loc][i][n+1] = 0.0f;
        n += 2;
      }
    }

    n = 0;

    for (int j = 0; j < nfft; j++) {
      potent[loc][j] = 0.0f;
      potent_hat[loc][n] = 0.0f;
      potent_hat[loc][n+1] = 0.0f;
      n += 2;
    }
    return;
  }


  int m,l,k;
  double xprd=domain->xprd;
  double yprd=domain->yprd;
  double zprd=domain->zprd;
  double zper,yper,xper;
  double mdr2;
  double scale_inv = 1.0/(nx_tild * ny_tild * nz_tild);

  double vole = domain->xprd * domain->yprd * domain->zprd; // Note: the factor of V comes from the FFT

  double pref;

  if (intrxn.type == GAUSSIAN){

    init_potential(potent[loc],1, intrxn.parameters.data());
    init_potential_ft(potent_hat[loc], 1, intrxn.parameters.data());

  } else if (intrxn.type == ERFC){
    vector<double> params;
    params.push_back(intrxn.parameters[0]);
    params.push_back(intrxn.parameters[1]);
    params.push_back(1);

    init_potential(tmp, 2, params.data());

    int j = 0;
    for (int i = 0; i < nfft; i++) {
      ktmp[j++] = tmp[i];
      ktmp[j++] = ZEROF;
    }

    fft1->compute(ktmp, work1, FFT3d::FORWARD);

    FFT_SCALAR temp;
    for (int nn = 0; nn < 2*nfft; nn += 2) {
      temp = -work1[nn+1] * work1[nn+1] + work1[nn] * work1[nn];
      work1[nn+1] = (work1[nn+1] * work1[nn] + work1[nn] * work1[nn+1]) * scale_inv * scale_inv;
      work1[nn] = temp * scale_inv * scale_inv;
      potent_hat[loc][nn] = work1[nn];
      potent_hat[loc][nn+1] = work1[nn+1];
    }

    fft1->compute(work1, work1, FFT3d::BACKWARD);

    n = 0;
    for (int j = 0; j < nfft; j++){
      potent[loc][j] = work1[n];
      n += 2;
    }

  } else if (intrxn.type == GAUSSIAN_ERFC){

    init_potential_ft(ktmp2,1, intrxn.parameters.data());
    
    vector<double> params;
    params.push_back(intrxn.parameters[1]);
    params.push_back(intrxn.parameters[2]);
    params.push_back(1);

    init_potential(tmp, 2, params.data());

    int j = 0;
    for (int i = 0; i < nfft; i++) {
      ktmp[j++] = tmp[i];
      ktmp[j++] = ZEROF;
    }

    fft1->compute(ktmp, work2, FFT3d::FORWARD);
    for (int i = 0; i < 2 * nfft; i++) {
      ktmp[i] = work2[i] * scale_inv;
    }

    // Convolution of potentials
    n = 0;
    for (int i = 0; i < nfft; i++) {
      complex_multiply(ktmp, ktmp2, work1, n);
      potent_hat[loc][n] = work1[n];
      potent_hat[loc][n + 1] = work1[n + 1];
      n += 2;
    }

    fft1->compute(work1, work1, FFT3d::BACKWARD);

    n = 0;
    for (int j = 0; j < nfft; j++){
      potent[loc][j] = work1[n];
      n += 2;
    }
  }

}

/* ----------------------------------------------------------------------
   Initialize potentials in fourier space when possible
------------------------------------------------------------------------- */

void TILD::init_potential_ft(FFT_SCALAR *wk1, const int type, const double* parameters){

  int n = 0;
  double xprd=domain->xprd;
  double yprd=domain->yprd;
  double zprd=domain->zprd;
  double zper,yper,xper;
  double k2;
  double factor = 4.0 * MY_PI * MY_PI;

  if (type == 1) {
                                                                // should be 3/2 right?
    for (int z = nzlo_fft; z <= nzhi_fft; z++) {
      zper = static_cast<double>(z) / zprd;
      if (z >= nz_tild / 2.0) {
        zper -= static_cast<double>(nz_tild) / zprd;
      }
      double zper2 = factor * zper * zper; 

      for (int y = nylo_fft; y <= nyhi_fft; y++) {
        yper = static_cast<double>(y) / yprd;
        if (y >= ny_tild / 2.0) {
          yper -= static_cast<double>(ny_tild) / yprd;
        }
        double yper2 = factor * yper * yper; 

        for (int x = nxlo_fft; x <= nxhi_fft; x++) {
          xper = static_cast<double>(x) / xprd;
          if (x >= nx_tild / 2.0) {
            xper -= static_cast<double>(nx_tild) / xprd;
          }

          k2 = (factor * xper * xper) + yper2 + zper2;
          wk1[n++] = exp(-k2 * 0.5 * parameters[0]);
          wk1[n++] = ZEROF;
        
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Initialize potentials in real-space when fourier space is not possible
------------------------------------------------------------------------- */

void TILD::init_potential(FFT_SCALAR *wk1, const int type, const double* parameters){

  int m,l,k;
  int n = 0;
  double xprd=domain->xprd;
  double yprd=domain->yprd;
  double zprd=domain->zprd;
  double zper,yper,xper;
  double mdr2;

  double vole = domain->xprd * domain->yprd * domain->zprd; // Note: the factor of V comes from the FFT
  int Dim = domain->dimension;

  if (type == 1) {
                                                                // should be 3/2 right?
    double pref = vole / (pow( sqrt(2.0 * MY_PI * (parameters[0]) ), Dim));
    for (m = nzlo_fft; m <= nzhi_fft; m++) {
      zper = zprd * (static_cast<double>(m) / nz_tild);
      if (zper >= zprd / 2.0) {
        zper = zprd - zper;
      }

      for (l = nylo_fft; l <= nyhi_fft; l++) {
        yper = yprd * (static_cast<double>(l) / ny_tild);
        if (yper >= yprd / 2.0) {
          yper = yprd - yper;
        }

        for (k = nxlo_fft; k <= nxhi_fft; k++) {
          xper = xprd * (static_cast<double>(k) / nx_tild);
          if (xper >= xprd / 2.0) {
            xper = xprd - xper;
          }

          mdr2 = xper * xper + yper * yper + zper * zper;
          wk1[n++] = exp(-mdr2 * 0.5 / parameters[0]) * pref;
        }
      }
    }
  }
  else if (type == 2) {
    for (m = nzlo_fft; m <= nzhi_fft; m++) {
      zper = zprd * (static_cast<double>(m) / nz_tild);
      if (zper >= zprd / 2.0) {
        zper = zprd - zper;
      }

      for (l = nylo_fft; l <= nyhi_fft; l++) {
        yper = yprd * (static_cast<double>(l) / ny_tild);
        if (yper >= yprd / 2.0) {
          yper = yprd - yper;
        }

        for (k = nxlo_fft; k <= nxhi_fft; k++) {
          xper = xprd * (static_cast<double>(k) / nx_tild);
          if (xper >= xprd / 2.0) {
            xper = xprd - xper;
          }
          mdr2 = xper * xper + yper * yper + zper * zper;
          wk1[n++] = parameters[2]* 0.5 * (1.0 - erf((sqrt(mdr2) - parameters[0])/parameters[1])) * vole;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   check if all factors of n are in list of factors
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int TILD::factorable(int n)
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
   Get the kspace version of the inputted function/potential
   return the grid points in k-space
------------------------------------------------------------------------- */

void TILD::get_k_alias(FFT_SCALAR* wk1, FFT_SCALAR **out){
  int Dim = domain->dimension; 
  double k[Dim];
  int x, y, z;
  int n=0;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];

  int has_nyquist = 0;
  if ( nz_tild % 2 == 0 || ny_tild % 2 == 0 || nx_tild % 2 == 0 ) {
     has_nyquist = 1;
  }

  for (z = nzlo_fft; z <= nzhi_fft; z++) {
    if (nz_tild % 2 == 0 && z == nz_tild / 2)
      k[2] = 0.0;
    else if (double(z) < double(nz_tild) / 2.0)
      k[2] = 2 * MY_PI * double(z) / zprd;
    else
      k[2] = 2 * MY_PI * double(z - nz_tild) / zprd;

    for (y = nylo_fft; y <= nyhi_fft; y++) {
      if (ny_tild % 2 == 0 && y == ny_tild / 2) {
        k[1] = 0.0;
      } else if (double(y) < double(ny_tild) / 2.0)
        k[1] = 2 * MY_PI * double(y) / yprd;
      else
        k[1] = 2 * MY_PI * double(y - ny_tild) / yprd;
      for (x = nxlo_fft; x <= nxhi_fft; x++) {
        if (nx_tild % 2 == 0 && x == nx_tild / 2) {
          k[0] = 0.0;
        } else if (double(x) < double(nx_tild) / 2.0)
          k[0] = 2 * MY_PI * double(x) / xprd;
        else
          k[0] = 2 * MY_PI * double(x - nx_tild) / xprd;
        out[0][n] = wk1[n + 1] * k[0];
        out[0][n + 1] = -wk1[n] * k[0];
        out[1][n] = wk1[n + 1] * k[1];
        out[1][n + 1] = -wk1[n] * k[1];
        out[2][n] = wk1[n + 1] * k[2];
        out[2][n + 1] = -wk1[n] * k[2];
        n += 2;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */

void TILD::particle_map()
{
  int nx,ny,nz;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int flag = 0;

  if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

  for (int i = 0; i < nlocal; i++) {

    // order = even:
    //   (nx,ny,nz) = global index of grid pt to "lower left" of charge
    // order = odd:
    //   (nx,ny,nz) = global index of grid pt closest to charge due to shift
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
        nz+nlower < nzlo_out || nz+nupper > nzhi_out)
      flag = 1;
  }

  if (flag) error->one(FLERR,"Out of range atoms - cannot compute TILD");
}

/* ----------------------------------------------------------------------
   function for parsing TILD specific settings 
------------------------------------------------------------------------- */

int TILD::modify_param(int narg, char** arg)
{
  int ntypes = atom->ntypes;
  int iarg = 0;
  if (strcmp(arg[iarg], "tild") == 0) {
    iarg += 1;

    while (iarg < narg) {

      if (strcmp(arg[iarg], "prefactor") == 0) {
        if (domain->box_exist == 0)
          error->all(FLERR, "TILD command before simulation box is defined");

        if (iarg + 4 > narg) error->all(FLERR, "Illegal kspace_modify tild command");

        int ilo, ihi, jlo, jhi;
        utils::bounds(FLERR, arg[iarg + 1], 1, ntypes, ilo, ihi, error);
        utils::bounds(FLERR, arg[iarg + 2], 1, ntypes, jlo, jhi, error);
        double lprefactor = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        for (int i = ilo; i <= ihi; i++) {
          for (int j = MAX(jlo, i); j <= jhi; j++) { chi[i][j] = lprefactor; }
        }
        iarg += 4;
      } else if (strcmp(arg[iarg], "shape") == 0) {
        if (iarg + 4 > narg) error->all(FLERR, "Illegal kspace_modify tild command");

        int ilo, ihi, jlo, jhi;

        //Read in types 1 and 2
        utils::bounds(FLERR, arg[iarg + 1], 1, ntypes, ilo, ihi, error);
        utils::bounds(FLERR, arg[iarg + 2], 1, ntypes, jlo, jhi, error);

        //Obtain Gaussian width parameter
        if (strcmp(arg[iarg + 3], "gaussian") == 0) {
          if (iarg + 5 > narg) error->all(FLERR, "Illegal kspace_modify tild command");

          double a2_one = utils::numeric(FLERR, arg[iarg + 4], false, lmp) *
              utils::numeric(FLERR, arg[iarg + 4], false, lmp);
          for (int i = ilo; i <= ihi; i++) {
            for (int j = MAX(jlo, i); j <= jhi; j++) {
              potent_type_map[1][i][j] = 1;
              potent_type_map[0][i][j] = 0;
              a2[i][j] = a2_one;
            }
          }
          iarg += 5;
        } else if (strcmp(arg[iarg + 3], "erfc") == 0) {

          //Obtain Erfc Radius, width, and density for Erf
          if (iarg + 7 > narg) error->all(FLERR, "Illegal kspace_modify tild command");
          double rp_one = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
          double xi_one = utils::numeric(FLERR, arg[iarg + 5], false, lmp);
          double np_rho_one = utils::numeric(FLERR, arg[iarg + 6], false, lmp);
          for (int i = ilo; i <= ihi; i++) {
            for (int j = MAX(jlo, i); j <= jhi; j++) {
              potent_type_map[2][i][j] = 1;
              potent_type_map[0][i][j] = 0;
              rp[i][j] = rp_one;
              xi[i][j] = xi_one;
              np_rho[i][j] = np_rho_one;
            }
          }
          iarg += 7;
        } else if (strcmp(arg[iarg + 3], "none") == 0) {
          for (int i = ilo; i <= ihi; i++) {
            for (int j = MAX(jlo, i); j <= jhi; j++) {
              potent_type_map[0][i][j] = 1;
              for (int istyle = 1; istyle <= nstyles; istyle++) potent_type_map[istyle][i][j] = 0;
            }
          }
          iarg += 4;
        } else
          error->all(FLERR, "Illegal kspace_modify tild shape density function argument");

      } else if (strcmp(arg[iarg], "cross-interaction") == 0) {
        if (iarg + 3 > narg) error->all(FLERR, "Illegal kspace_modify tild command");
        int ilo, ihi, jlo, jhi;
        utils::bounds(FLERR, arg[iarg + 1], 1, ntypes, ilo, ihi, error);
        utils::bounds(FLERR, arg[iarg + 2], 1, ntypes, jlo, jhi, error);

        cross_type tmp_type;
        int nparams;
        int temp_arg = 0;
        if (strcmp(arg[iarg + 3], "gaussian") == 0) {
          tmp_type = GAUSSIAN;
          nparams = 1;
          temp_arg = 5;
        } else if (strcmp(arg[iarg + 3], "erfc") == 0) {
          tmp_type = ERFC;
          nparams = 2;
          temp_arg = 6;
        } else if (strcmp(arg[iarg + 3], "gaussian_erfc") == 0) {
          tmp_type = GAUSSIAN_ERFC;
          nparams = 3;
          temp_arg = 7;
        } else if (strcmp(arg[iarg + 3], "none") == 0) {
          tmp_type = NONE;
          nparams = 0;
          temp_arg = 4;
        } else if (strcmp(arg[iarg + 3], "delete") == 0) {
          tmp_type = DELETE;
          temp_arg = 4;
        } else {
          error->all(FLERR, "Illegal tild cross-interaction specified.");
        }

        if (tmp_type != DELETE) {
          for (int i = ilo; i <= ihi; i++) {
            for (int j = jlo; j <= jhi; j++) {
              Interaction temp;
              temp.i = i - 1;
              temp.j = j - 1;
              temp.type = tmp_type;

              for (int p = 0; p < nparams; p++) {
                temp.parameters.push_back(utils::numeric(FLERR, arg[iarg + 4 + p], false, lmp));
              }

              auto it = cross_iter.begin();
              int idx = 0;
              while (it != cross_iter.end()) {
                if ((it->i == i - 1) && (it->j == j - 1)) cross_iter[idx] = temp;
                it++;
              }

              if (it == cross_iter.end()) cross_iter.emplace_back(temp);
            }
          }
        } else {
          for (int ii = ilo; ii <= ihi; ii++) {
            for (int jj = jlo, ii; jj <= jhi; jj++) {
              auto it = cross_iter.begin();
              while (it != cross_iter.end()) {
                if ((it->i == ii - 1) && (it->j == jj - 1)) { cross_iter.erase(it--); }
                it++;
              }
            }
          }
        }
        iarg += temp_arg;

      } else if (strcmp(arg[iarg], "mix") == 0) {
        if (iarg + 2 > narg) error->all(FLERR, "Illegal kspace_modify tild command");
        mix_flag = 1;
        if (strcmp(arg[iarg + 1], "convolution") == 0)
          mix_flag = 1;
        else if (strcmp(arg[iarg + 1], "define") == 0)
          mix_flag = 0;
        else
          error->all(FLERR, "Illegal kspace_modify tild mix argument");
        iarg += 2;
      } else if (strcmp(arg[iarg], "set_rho0") == 0) {
        if (iarg + 2 > narg) error->all(FLERR, "Illegal kspace_modify tild command");
        set_rho0 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        iarg += 2;
      } else if (strcmp(arg[iarg], "subtract_rho0") == 0) {
        if (iarg + 2 > narg) error->all(FLERR, "Illegal kspace_modify tild command");
        if (strcmp(arg[iarg + 1], "yes") == 0)
          sub_flag = 1;
        else if (strcmp(arg[iarg + 1], "no") == 0)
          sub_flag = 0;
        else
          error->all(FLERR, "Illegal kspace_modify tild subtract_rho0 argument");
        iarg += 2;
      } else if (strcmp(arg[iarg], "normalize_by_rho0") == 0) {
        if (iarg + 2 > narg) error->all(FLERR, "Illegal kspace_modify tild command");
        if (strcmp(arg[iarg + 1], "yes") == 0)
          norm_flag = 1;
        else if (strcmp(arg[iarg + 1], "no") == 0)
          norm_flag = 0;
        else
          error->all(FLERR, "Illegal kspace_modify tild normalize_by_rho0 argument");
        iarg += 2;
      } else if (strcmp(arg[iarg], "write_grid_data") == 0) {
        if (iarg + 3 > narg) error->all(FLERR, "Illegal kspace_modify tild command");
        grid_data_output_freq = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
        if (grid_data_output_freq <= 0)
          write_grid_flag = 0;
        else
          write_grid_flag = 1;
        strcpy(grid_data_filename, arg[iarg + 2]);
        iarg += 3;
      } else if (strcmp(arg[iarg], "ave/grid") == 0) {
        if (iarg + 5 > narg) error->all(FLERR, "Illegal kspace_modify tild command");

        ave_grid_flag = 1;
        nevery = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
        nrepeat = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
        peratom_freq = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
        strcpy(ave_grid_filename, arg[iarg + 4]);
        nvalid = nextvalid();
        if (nevery <= 0 || nrepeat <= 0 || peratom_freq <= 0)
          error->all(FLERR, "Illegal fix tild/ave/grid command");
        if (peratom_freq % nevery || nrepeat * nevery > peratom_freq)
          error->all(FLERR, "Illegal kspace_modify tild/ave/grid command");
        iarg += 5;
      } else
        error->all(FLERR, "Illegal kspace_modify tild command");
    }
  }
  return iarg;
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void TILD::pack_forward_grid(int flag, void *vbuf, int nlist, int *list)
{
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;
  int n = 0;

  if (flag == FORWARD_NONE){
    for (int ktype = 1; ktype <= atom->ntypes; ktype++) {
      if ( density_flags[ktype] != 1 ) continue;
      //for (int j = 0; j < Dim; j++) {
      FFT_SCALAR *xsrc = &gradWtypex[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *ysrc = &gradWtypey[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *zsrc = &gradWtypez[ktype][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++) {
        buf[n++] = xsrc[list[i]];
        buf[n++] = ysrc[list[i]];
        buf[n++] = zsrc[list[i]];
      }
    }
  } else if (flag == FORWARD_GRID_DEN){
    for (int ktype = 1; ktype <= atom->ntypes; ktype++) {
      if ( density_flags[ktype] != 1 ) continue;
      FFT_SCALAR *srcx = &density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *srcy = &density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *srcz = &density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++) {
        buf[n++] = srcx[list[i]];
        buf[n++] = srcy[list[i]];
        buf[n++] = srcz[list[i]];
      }
    }
  } else if (flag == FORWARD_AVG_GRID_DEN){
    for (int ktype = 1; ktype <= atom->ntypes; ktype++) {
      if ( density_flags[ktype] != 1 ) continue;
      FFT_SCALAR *srcx = &avg_density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *srcy = &avg_density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *srcz = &avg_density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++) {
        buf[n++] = srcx[list[i]];
        buf[n++] = srcy[list[i]];
        buf[n++] = srcz[list[i]];
      }
    }
  }
}


/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void TILD::unpack_forward_grid(int flag, void *vbuf, int nlist, int *list)
{
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;

  int n = 0;
  //int Dim = domain->dimension;

  if (flag == FORWARD_NONE){
    for (int ktype = 1; ktype <= atom->ntypes; ktype++) {
      if ( density_flags[ktype] != 1 ) continue;
      FFT_SCALAR *destx = &gradWtypex[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *desty = &gradWtypey[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *destz = &gradWtypez[ktype][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++) {
        destx[list[i]] = buf[n++];
        desty[list[i]] = buf[n++];
        destz[list[i]] = buf[n++];
      }
    }
  } else if (flag == FORWARD_GRID_DEN){
    for (int ktype = 1; ktype <= atom->ntypes; ktype++) {
      if ( density_flags[ktype] != 1 ) continue;
      FFT_SCALAR *destx = &density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *desty = &density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *destz = &density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++) {
        destx[list[i]] = buf[n++];
        desty[list[i]] = buf[n++];
        destz[list[i]] = buf[n++];
      }
    }
  } else if (flag == FORWARD_AVG_GRID_DEN){
    for (int ktype = 1; ktype <= atom->ntypes; ktype++) {
      if ( density_flags[ktype] != 1 ) continue;
      FFT_SCALAR *destx = &avg_density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *desty = &avg_density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *destz = &avg_density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++) {
        destx[list[i]] = buf[n++];
        desty[list[i]] = buf[n++];
        destz[list[i]] = buf[n++];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void TILD::pack_reverse_grid(int flag, void *vbuf, int nlist, int *list) {
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;

  int n = 0;
  if (flag == REVERSE_RHO_NONE) {
    for (int ktype = 1; ktype <= atom->ntypes; ktype++) {
      if ( density_flags[ktype] != 1 ) continue;
      FFT_SCALAR *src = &density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++)
        buf[n++] = src[list[i]];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void TILD::unpack_reverse_grid(int flag, void *vbuf, int nlist, int *list)
{
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;  
  int n = 0;
  if (flag == REVERSE_RHO_NONE) {
    for (int ktype = 1; ktype <= atom->ntypes; ktype++) {
      if ( density_flags[ktype] != 1 ) continue;
      FFT_SCALAR *dest = &density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++)
        dest[list[i]] += buf[n++];
    }
  }
}

/* ----------------------------------------------------------------------
   map nprocs to NX by NY grid as PX by PY procs - return optimal px,py
------------------------------------------------------------------------- */

void TILD::procs2grid2d(int nprocs, int nx, int ny, int *px, int *py)
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
   set local subset of TILD/FFT grid that I own
   n xyz lo/hi in = 3d brick that I own (inclusive)
   n xyz lo/hi out = 3d brick + ghost cells in 6 directions (inclusive)
   n xyz lo/hi fft = FFT columns that I own (all of x dim, 2d decomp in yz)
------------------------------------------------------------------------- */

void TILD::set_grid_local()
{
  // global indices of TILD grid range from 0 to N-1
  // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
  //   global TILD grid that I own without ghost cells
  // for slab TILD, assign z grid as if it were not extended
  // both non-tiled and tiled proc layouts use 0-1 fractional sumdomain info

  if (comm->layout != Comm::LAYOUT_TILED) {
    nxlo_in = static_cast<int> (comm->xsplit[comm->myloc[0]] * nx_tild);
    nxhi_in = static_cast<int> (comm->xsplit[comm->myloc[0]+1] * nx_tild) - 1;

    nylo_in = static_cast<int> (comm->ysplit[comm->myloc[1]] * ny_tild);
    nyhi_in = static_cast<int> (comm->ysplit[comm->myloc[1]+1] * ny_tild) - 1;

    nzlo_in = static_cast<int>
      (comm->zsplit[comm->myloc[2]] * nz_tild/slab_volfactor);
    nzhi_in = static_cast<int>
      (comm->zsplit[comm->myloc[2]+1] * nz_tild/slab_volfactor) - 1;

  } else {
    nxlo_in = static_cast<int> (comm->mysplit[0][0] * nx_tild);
    nxhi_in = static_cast<int> (comm->mysplit[0][1] * nx_tild) - 1;

    nylo_in = static_cast<int> (comm->mysplit[1][0] * ny_tild);
    nyhi_in = static_cast<int> (comm->mysplit[1][1] * ny_tild) - 1;

    nzlo_in = static_cast<int> (comm->mysplit[2][0] * nz_tild/slab_volfactor);
    nzhi_in = static_cast<int> (comm->mysplit[2][1] * nz_tild/slab_volfactor) - 1;
  }

  // nlower,nupper = stencil size for mapping particles to TILD grid

  nlower = -(order-1)/2;
  nupper = order/2;

  // shift values for particle <-> grid mapping
  // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

  if (order % 2) shift = OFFSET + 0.5;
  else shift = OFFSET;
  if (order % 2) shiftone = 0.0;
  else shiftone = 0.5;

  // nlo_out,nhi_out = lower/upper limits of the 3d sub-brick of
  //   global TILD grid that my particles can contribute charge to
  // effectively nlo_in,nhi_in + ghost cells
  // nlo,nhi = index of global grid pt to "lower left" of smallest/largest
  //           position a particle in my box can be at
  // dist[3] = particle position bound = subbox + skin/2.0 
  //   convert to triclinic if necessary
  // nlo_out,nhi_out = nlo,nhi + stencil size for particle mapping
  // for slab TILD, assign z grid as if it were not extended

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

  double dist[3] = {0.0,0.0,0.0};
  double cuthalf = 0.5*neighbor->skin;
  if (triclinic == 0) dist[0] = dist[1] = dist[2] = cuthalf;
  else kspacebbox(cuthalf,&dist[0]);

  int nlo,nhi;
  nlo = nhi = 0;

  nlo = static_cast<int> ((sublo[0]-dist[0]-boxlo[0]) *
                            nx_tild/xprd + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[0] + dist[0]-boxlo[0]) *
                            nx_tild/xprd + shift) - OFFSET;
  nxlo_out = nlo + nlower;
  nxhi_out = nhi + nupper;

  nlo = static_cast<int> ((sublo[1]-dist[1]-boxlo[1]) *
                            ny_tild/yprd + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[1]+dist[1]-boxlo[1]) *
                            ny_tild/yprd + shift) - OFFSET;
  nylo_out = nlo + nlower;
  nyhi_out = nhi + nupper;

  nlo = static_cast<int> ((sublo[2]-dist[2]-boxlo[2]) *
                            nz_tild/zprd_slab + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[2]+dist[2]-boxlo[2]) *
                            nz_tild/zprd_slab + shift) - OFFSET;
  nzlo_out = nlo + nlower;
  nzhi_out = nhi + nupper;

  if (stagger_flag) {
    nxhi_out++;
    nyhi_out++;
    nzhi_out++;
  }

  // for slab TILD, change the grid boundary for processors at +z end
  //   to include the empty volume between periodically repeating slabs
  // for slab TILD, want charge data communicated from -z proc to +z proc,
  //   but not vice versa, also want field data communicated from +z proc to
  //   -z proc, but not vice versa
  // this is accomplished by nzhi_in = nzhi_out on +z end (no ghost cells)
  // also insure no other procs use ghost cells beyond +z limit
  // differnet logic for non-tiled vs tiled decomposition

  if (slabflag == 1) {
    if (comm->layout != Comm::LAYOUT_TILED) {
      if (comm->myloc[2] == comm->procgrid[2]-1) nzhi_in = nzhi_out = nz_tild - 1;
    } else {
      if (comm->mysplit[2][1] == 1.0) nzhi_in = nzhi_out = nz_tild - 1;
    }
    nzhi_out = MIN(nzhi_out,nz_tild-1);
  }

  // x-pencil decomposition of FFT mesh
  // global indices range from 0 to N-1
  // each proc owns entire x-dimension, clumps of columns in y,z dimensions
  // npey_fft,npez_fft = # of procs in y,z dims
  // if nprocs is small enough, proc can own 1 or more entire xy planes,
  //   else proc owns 2d sub-blocks of yz plane
  // me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
  // nlo_fft,nhi_fft = lower/upper limit of the section
  //   of the global FFT mesh that I own in x-pencil decomposition

  int npey_fft,npez_fft;
  if (nz_tild >= nprocs) {
    npey_fft = 1;
    npez_fft = nprocs;
  } else procs2grid2d(nprocs,ny_tild,nz_tild,&npey_fft,&npez_fft);

  int me_y = me % npey_fft;
  int me_z = me / npey_fft;

  nxlo_fft = 0;
  nxhi_fft = nx_tild - 1;
  nylo_fft = me_y*ny_tild/npey_fft;
  nyhi_fft = (me_y+1)*ny_tild/npey_fft - 1;
  nzlo_fft = me_z*nz_tild/npez_fft;
  nzhi_fft = (me_z+1)*nz_tild/npez_fft - 1;

  // ngrid = count of TILD grid pts owned by this proc, including ghosts

  ngrid = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
    (nzhi_out-nzlo_out+1);

  // count of FFT grids pts owned by this proc, without ghosts
  // nfft = FFT points in x-pencil FFT decomposition on this proc
  // nfft_brick = FFT points in 3d brick-decomposition on this proc
  // nfft_both = greater of 2 values

  nfft = (nxhi_fft-nxlo_fft+1) * (nyhi_fft-nylo_fft+1) *
    (nzhi_fft-nzlo_fft+1);
  int nfft_brick = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1) *
    (nzhi_in-nzlo_in+1);
  nfft_both = MAX(nfft,nfft_brick);
}


/* ----------------------------------------------------------------------
   create discretized density on section of global grid due to my particles
   density(x,y,z) = density at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void TILD::make_rho()
{
  int nx,ny,nz,mx,my,mz;
  int ntypes = atom->ntypes;
  FFT_SCALAR dx,dy,dz,x0,y0,z0,w;

  for (int k = 0; k <= ntypes; k++) {
    memset(&(density_brick_types[k][nzlo_out][nylo_out][nxlo_out]),0,
           ngrid*sizeof(FFT_SCALAR));
  }

  // loop over my particles, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  for (int i = 0; i < nlocal; i++) {
    
    // Skip if the particle type has no TILD interaction
    if (density_flags[type[i]] != 1) continue;

    // do the following for all grids
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];

    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);
    z0 = delvolinv;
    // loop over the stencil
    for (int n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      for (int m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        for (int l = nlower; l <= nupper; l++) {
          mx = l+nx;
          w = x0*rho1d[0][l];

          density_brick_types[type[i]][mz][my][mx] += w;
          density_brick_types[0][mz][my][mx] += w;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   remap density from 3d brick decomposition to FFT decomposition
------------------------------------------------------------------------- */

void TILD::brick2fft()
{
  int ix,iy,iz;
  int ntypes = atom->ntypes;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values


  for (int k = 1; k <= ntypes; k++) {
    int n = 0;
    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++){
          density_fft_types[k][n++] = density_brick_types[k][iz][iy][ix];
        }
  }

  for (int k = 1; k <= ntypes; k++) {
    remap->perform(density_fft_types[k],density_fft_types[k],work1); 
  }
}


/* ----------------------------------------------------------------------
   density assignment into rho1d
   dx,dy,dz = distance of particle from "lower left" grid point
------------------------------------------------------------------------- */

void TILD::compute_rho1d(const FFT_SCALAR &dx, const FFT_SCALAR &dy,
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
   calculates the gradient of the density fields and convolves it with 
   the interaction potential
------------------------------------------------------------------------- */

void TILD::accumulate_gradient() {
  int Dim = domain->dimension;
  int n = 0;
  int ntypes = atom->ntypes;
  
  precompute_density_hat_fft(); 

  for ( int ktype = 1; ktype <= ntypes; ktype++) {
    memset(&(gradWtypex[ktype][nzlo_out][nylo_out][nxlo_out]), 0,
           ngrid * sizeof(FFT_SCALAR));
    memset(&(gradWtypey[ktype][nzlo_out][nylo_out][nxlo_out]), 0,
           ngrid * sizeof(FFT_SCALAR));
    memset(&(gradWtypez[ktype][nzlo_out][nylo_out][nxlo_out]), 0,
           ngrid * sizeof(FFT_SCALAR));
  }

  // This part is for the specified chi interactions.

  FFT_SCALAR tmp_sub = 0.0;
  if (subtract_rho0 == 1)
    tmp_sub = rho0;

  int loc = -1;
  for (int itype = 1; itype <= ntypes; itype++) {
    for (int jtype = itype; jtype <= ntypes; jtype++) {
      loc++;
      if ( potent_type_map[0][itype][jtype] == 1) continue;

      double tmp_chi = chi[itype][jtype];
      if (normalize_by_rho0 == 1) tmp_chi /= rho0;

      if ( tmp_chi == 0 ) continue;

      int diff_type = 1;
      if (itype == jtype) {
        diff_type = 0;
      }

      if (eflag_global || vflag_global) {
        ev_calculation(loc, itype, jtype);
      }

      for (int i = 0; i < Dim; i++) {

        n = 0;
        for (int k = 0; k < nfft; k++) {
          complex_multiply(grad_potent_hat[loc][i], density_hat_fft_types[itype], ktmp2i, n);
          if (diff_type) complex_multiply(grad_potent_hat[loc][i], density_hat_fft_types[jtype], ktmp2j, n);
          n += 2;
        }

        fft2->compute(ktmp2i, ktmpi, FFT3d::BACKWARD);
        if (diff_type) fft2->compute(ktmp2j, ktmpj, FFT3d::BACKWARD);

        n = 0;
        int j = 0;
        for (int k = nzlo_in; k <= nzhi_in; k++)  {
          for (int m = nylo_in; m <= nyhi_in; m++) {
            for (int o = nxlo_in; o <= nxhi_in; o++) {
              if ( i == 0 ) {
                gradWtypex[jtype][k][m][o] += ktmpi[n] * tmp_chi;
                if (diff_type) gradWtypex[itype][k][m][o] += ktmpj[n] * tmp_chi;
              } else if ( i == 1 ) {
                gradWtypey[jtype][k][m][o] += ktmpi[n] * tmp_chi;
                if (diff_type) gradWtypey[itype][k][m][o] += ktmpj[n] * tmp_chi;
              } else {
                gradWtypez[jtype][k][m][o] += ktmpi[n] * tmp_chi;
                if (diff_type) gradWtypez[itype][k][m][o] += ktmpj[n] * tmp_chi;
              }
              n += 2;
              j++;
            }
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Applies the from the grid onto each particle
------------------------------------------------------------------------- */

void TILD::fieldforce_param(){
  int nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz;

  // loop over my charges, interpolate second desnity field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of density
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double **x = atom->x;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  // Convert field to force per particle

  int *type = atom->type;

  for (int i = 0; i < nlocal; i++) {

    int temp_type = type[i];
    if (temp_type == -1) continue;
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);

    ekx = eky = ekz = ZEROF;
    for (int n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (int m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (int l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          ekx += x0 * gradWtypex[temp_type][mz][my][mx];
          eky += x0 * gradWtypey[temp_type][mz][my][mx];
          ekz += x0 * gradWtypez[temp_type][mz][my][mx];
        }
      }
    }

    // convert field to force
    
    f[i][0] += ekx;
    f[i][1] += eky;
    f[i][2] += ekz;
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

void TILD::compute_rho_coeff()
{
  int j,k,l,m;
  FFT_SCALAR s;

  FFT_SCALAR **a;
  memory->create2d_offset(a,order,-order,order,"tild:a");

  for (k = -order; k <= order; k++)
    for (l = 0; l < order; l++)
      a[l][k] = 0.0;

  a[0][0] = 1.0;
  for (j = 1; j < order; j++) {
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

  m = (1-order)/2;
  for (k = -(order-1); k < order; k += 2) {
    for (l = 0; l < order; l++)
      rho_coeff[l][m] = a[l][k];
    for (l = 1; l < order; l++)
      drho_coeff[l-1][m] = l*a[l][k];
    m++;
  }

  memory->destroy2d_offset(a,-order);
}

/* ----------------------------------------------------------------------
   Abstraction to hide the annoying multiplication.
------------------------------------------------------------------------- */

void TILD::complex_multiply(FFT_SCALAR *in1, FFT_SCALAR  *in2, FFT_SCALAR  *out, int n){
  out[n] = ((in1[n] * in2[n]) - (in1[n + 1] * in2[n + 1]));
  out[n + 1] = ((in1[n + 1] * in2[n]) + (in1[n] * in2[n + 1]));
}

/* ----------------------------------------------------------------------
   Calculation of the energy and the virials
------------------------------------------------------------------------- */

void TILD::ev_calculation(const int loc, const int itype, const int jtype) {
  int n = 0;
  double eng, engk;
  double vtmp, vtmpk;
  double scale_inv = 1.0/(nx_tild *ny_tild * nz_tild);
  double V = domain->xprd * domain->yprd * domain->zprd;

  double tmp_rho_div = 1.0;
  if (normalize_by_rho0 == 1) tmp_rho_div = rho0;

  double type_factor = 1.0;
  if (itype == jtype) type_factor = 0.5;

  double factor = scale_inv / tmp_rho_div * type_factor * chi[itype][jtype];

  if (eflag_global) {
    // convolve itype-jtype interaction potential and itype density
    n = 0;
    for (int k = 0; k < nfft; k++) {
      complex_multiply(density_hat_fft_types[itype], potent_hat[loc], ktmpi, n);
      n += 2;
    }
    // IFFT the convolution
    fft1->compute(ktmpi, ktmp2i, FFT3d::BACKWARD);

    // chi and kappa energy contribution
    n = 0;
    eng = 0.0;
    for (int k = 0; k < nfft; k++) {
      // rho_i * u_ij * rho_j * chi * prefactor
      eng += ktmp2i[n] * density_fft_types[jtype][k];
      n += 2;
    }
    energy += eng * factor * V; 
  }

  // pressure tensor calculation
  if (vflag_global) {
    // loop over stress tensor
    for (int i = 0; i < 6; i++) {
      n = 0;
      for (int k = 0; k < nfft; k++) {
        complex_multiply(density_hat_fft_types[itype], vg_hat[loc][i], ktmpi, n);
        n += 2;
      }
      fft1->compute(ktmpi, ktmp2i, FFT3d::BACKWARD);
      n=0;
      vtmp = 0.0;
      vtmpk = 0.0;
      for (int k = 0; k < nfft; k++) {
        vtmp += ktmp2i[n] * density_fft_types[jtype][k];
        n += 2;
      }
      virial[i] += vtmp * factor;
      // diagonal IGEOS/chain-ruled on-diagonal part is called by including the kinetic energy term
    }
  }
}

/* ----------------------------------------------------------------------
  Calculates the rho0 needed for this system instead of the rho0 that 
  comes from the particle densities. 
  Each Gaussian particle particle contributes a 1/V and erfc particles 
    contributes 4/3*pi*r^3/V.
------------------------------------------------------------------------- */

double TILD::calculate_rho0(){
  MPI_Bcast(&set_rho0, 1, MPI_DOUBLE, 0, world);
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int ntypes = atom->ntypes;
  int particles_not_tild = 0;
  vector<int> count_per_type (ntypes+1,0);
  vector<int> count_per_type_local (ntypes+1,0);
  double lmass = 0, lmass_all = 0;
  int tild_particles_non_gaussian = 0;
  int particles_gaussian = 0;
  double rho = 0;

  for (int i = 0; i < nlocal; i++) {
    if ( potent_type_map[0][type[i]][type[i]] == 1) {
      particles_not_tild++;
    } else {
      count_per_type_local[type[i]]++;
      if (potent_type_map[1][type[i]][type[i]] != 1 ){
        tild_particles_non_gaussian++;
      }
    }
  }

  int temp_int;
  MPI_Allreduce(&tild_particles_non_gaussian, &temp_int, 1, MPI_INT, MPI_SUM, world);
  tild_particles_non_gaussian = temp_int;
  MPI_Allreduce(&particles_not_tild, &temp_int, 1, MPI_INT, MPI_SUM, world);
  particles_not_tild = temp_int;

  MPI_Allreduce(count_per_type_local.data(), count_per_type.data(), ntypes + 1, MPI_INT, MPI_SUM, world);

  double vole = domain->xprd * domain->yprd * domain->zprd;

  for (int itype = 1; itype <= ntypes; itype++) {
    // gaussian potential, has no radius
    if (potent_type_map[1][itype][itype] == 1) {
      lmass += count_per_type[itype];
    } 
  }

  double temp_double;
  MPI_Allreduce(&lmass, &temp_double, 1, MPI_DOUBLE, MPI_SUM, world);
  lmass = temp_double;

  particles_gaussian = lmass;

  double total_sphere_vol = 0;
    for ( int itype = 1; itype <= ntypes; itype++) {
      if (potent_type_map[2][itype][itype] == 1) {
      total_sphere_vol += count_per_type[itype] * (4.0 * MY_PI / 3.0) * rp[itype][itype] * rp[itype][itype] * rp[itype][itype];
      }
    }

  if (set_rho0 <= 0.0f){
    for ( int itype = 1; itype <= ntypes; itype++) {
      if (potent_type_map[2][itype][itype] == 1) {
        lmass += count_per_type[itype] * (4.0 * MY_PI / 3.0) * rp[itype][itype] * rp[itype][itype] * rp[itype][itype] * np_rho[itype][itype];
      }
    }
  } else {
    MPI_Allreduce(&total_sphere_vol, &temp_double, 1, MPI_DOUBLE, MPI_SUM, world);
    total_sphere_vol = temp_double;
    modified_rho0 = (set_rho0 * vole - lmass) /total_sphere_vol;
    for ( int itype = 1; itype <= ntypes; itype++) {
      if (potent_type_map[2][itype][itype] == 1) {
        lmass += count_per_type[itype] * (4.0 * MY_PI / 3.0) * rp[itype][itype] * rp[itype][itype] * rp[itype][itype] * modified_rho0;
      }
    }
  }

  MPI_Allreduce(&lmass, &lmass_all, 1, MPI_DOUBLE, MPI_SUM, world);

  rho = lmass_all / vole;
  if (comm->me == 0) {
    std::string mesg = 
      fmt::format("  Found {} particles without a TILD potential\n", particles_not_tild);
    mesg += fmt::format("  Calculated rho = {:.6f}; user set rho0 = {:.6f} for TILD potential.\n",
                        rho, set_rho0);

    utils::logmesg(lmp, mesg);
    mesg.clear();

    if (set_rho0 > 0) {

      // If TILD particles are only Gaussian and a user shape is defined.
      if ((particles_gaussian == lmass && particles_gaussian > 0)) {
        mesg +=
            fmt::format("  For pure Gaussian shapes and a user_defined rho0, adopting user_rho0 "
                        "= {:.6f} for TILD potential.\n",
                        set_rho0);
        rho = set_rho0;
      } else if (total_sphere_vol > 0 && set_rho0 > 0) {
        if (modified_rho0 != 0.f) {
          std::string mesg_mod = fmt::format(
              "Modified NP rho set to {:.6f} for all particles to achieve user specified rho.\n", modified_rho0);
          error->warning(FLERR, mesg_mod);
          if (modified_rho0 < 0) {
            std::string mesg_mod = fmt::format("Modified NP rho is below zero and is unphysical.\n");
            error->all(FLERR, mesg_mod);
          }
        }
      }
      utils::logmesg(lmp, mesg);
      mesg.clear();
    }


    if (rho <= 0.f) {
      if (set_rho0 > 0.f) {
        error->warning(FLERR, "Calculated rho <= 0; using user specificed rho0");
        mesg += fmt::format("  Using user set rho0 = {:.6f} for TILD potential.\n", set_rho0);
        utils::logmesg(lmp, mesg);
        rho = set_rho0;
      }
    }
  }

  MPI_Bcast(&rho, 1, MPI_DOUBLE, 0, world);

  if (rho <= 0.f) {
    if (set_rho0 <= 0.f) {
      std::string mesg = fmt::format("  Both densities are not positive. Aborting run.");
      error->all(FLERR, mesg);
    }
  }
  return rho;
}

/* ----------------------------------------------------------------------
   Write out the gridded densities out to the system
------------------------------------------------------------------------- */

void TILD::write_grid_data( char *filename, const int avg ) {
  int ntypes = atom->ntypes;
  if (comm->me == 0) {
    otp = fopen( filename, "w" ) ;

    // header
    fprintf( otp, "# x y z") ;
    for (int itype = 1; itype <= ntypes; itype++)
      fprintf( otp, " rho_%d", itype) ;
    fprintf( otp, "\n") ;
  }

  // communication buffer for all my Atom info
  // max_size = largest buffer needed by any proc

  int ncol = ntypes + 3;

  int sendrow = nfft_both;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"TILD:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"TILD:buf");

  // pack my atom data into buf

  if (avg == 1) pack_avg_grid_data(buf);
  else pack_grid_data(buf);

  // write one chunk of atoms per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      for (int n = 0; n < recvrow; n++) {
        fprintf( otp, "%lf %lf %lf", buf[n][0], buf[n][1], buf[n][2]);
        for (int itype = 1; itype <= ntypes; itype++) {
          fprintf( otp, " %1.16e", buf[n][2+itype]);
        }
        fprintf( otp, "\n" ) ;
      }
    } 
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
  if (me == 0) fclose( otp ) ;
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void TILD::pack_grid_data(double **buf)
{
  int ntypes = atom->ntypes;
  double fx = domain->xprd/nx_tild;
  double fy = domain->yprd/ny_tild;
  double fz = domain->zprd/nz_tild;
  int n = 0;
  //int n = nzlo_in + nylo_in + nxlo_in;
  for (int iz = nzlo_in; iz <= nzhi_in; iz++) {
    for (int iy = nylo_in; iy <= nyhi_in; iy++) {
      for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf[n][0] = ix * fx;
        buf[n][1] = iy * fy;
        buf[n][2] = iz * fz;
        for (int itype = 1; itype <= ntypes; itype++) {
          buf[n][2+itype] = density_brick_types[itype][iz][iy][ix];
        }
        n++;
      }
    }
  }  
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void TILD::pack_avg_grid_data(double **buf)
{
  int ntypes = atom->ntypes;
  double fx = domain->xprd/nx_tild;
  double fy = domain->yprd/ny_tild;
  double fz = domain->zprd/nz_tild;
  int n = 0;
  //int n = nzlo_in + nylo_in + nxlo_in;
  for (int iz = nzlo_in; iz <= nzhi_in; iz++) {
    for (int iy = nylo_in; iy <= nyhi_in; iy++) {
      for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf[n][0] = ix * fx;
        buf[n][1] = iy * fy;
        buf[n][2] = iz * fz;
        for (int itype = 1; itype <= ntypes; itype++) {
          buf[n][2+itype] = avg_density_brick_types[itype][iz][iy][ix];
        }
        n++;
      }
    }
  }  
}

/* ----------------------------------------------------------------------
   Sums up density types to later normalize
------------------------------------------------------------------------- */

void TILD::sum_grid_data()
{
  int ntypes = atom->ntypes;
  for (int iz = nzlo_in; iz <= nzhi_in; iz++) {
    for (int iy = nylo_in; iy <= nyhi_in; iy++) {
      for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
        for (int itype = 1; itype <= ntypes; itype++) {
          avg_density_brick_types[itype][iz][iy][ix] += density_brick_types[itype][iz][iy][ix];
        }
      }
    }
  }  
}

/* ----------------------------------------------------------------------
   Used to normalize the gridded density by the amount of time steps
------------------------------------------------------------------------- */

void TILD::multiply_ave_grid_data(const double factor)
{
  int ntypes = atom->ntypes;
  for (int iz = nzlo_in; iz <= nzhi_in; iz++) {
    for (int iy = nylo_in; iy <= nyhi_in; iy++) {
      for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
        for (int itype = 1; itype <= ntypes; itype++) {
          avg_density_brick_types[itype][iz][iy][ix] *= factor;
        }
      }
    }
  }  
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint TILD::nextvalid()
{
  bigint nvalid = (update->ntimestep/peratom_freq)*peratom_freq + peratom_freq;
  if (nvalid-peratom_freq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += peratom_freq;
  return nvalid;
}

/* ----------------------------------------------------------------------
   write out the average gridded density if it is time to do so
------------------------------------------------------------------------- */

void TILD::ave_grid()
{
  // skip if not step which requires doing something
  // error check if timestep was reset in an invalid manner

  bigint ntimestep = update->ntimestep;
  if (ntimestep < nvalid_last || ntimestep > nvalid)
    error->all(FLERR,"Invalid timestep reset for fix ave/atom");
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // zero if first step

  if (irepeat == 0)
    multiply_ave_grid_data( 0.0 ); 

  // accumulate results of attributes,computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  sum_grid_data(); 

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;

  if (irepeat < nrepeat) {
    nvalid += nevery;
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+peratom_freq - (nrepeat-1)*nevery;

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  
  multiply_ave_grid_data( 1.0 / repeat );
  write_grid_data(ave_grid_filename, 1);
}

/* ----------------------------------------------------------------------
  Turn on density_flag for each particle type that has a TILD interaction
  Updates potent_type_map
------------------------------------------------------------------------- */
void TILD::manually_flip_density_flags(){

  int ntypes = atom->ntypes;
  for (int i = 0; i <= ntypes; i++){
    density_flags[i] = 0;
  }

  for (int i = 0; i <= ntypes; i++){
    for (int j = i; j <= ntypes; j++){
      for (int sty = 1; sty < nstyles; sty++){
        if (potent_type_map[sty][i][j] != 0){
          density_flags[i] = density_flags [j] = 1;
        }
      }
    }
  }

  for (auto& intrxn : cross_iter) {
    if (intrxn.type != NONE ){
      density_flags[intrxn.i+1] = density_flags[intrxn.j+1] = 1;
      potent_type_map[0][intrxn.i+1][intrxn.j+1] = potent_type_map[0][intrxn.j+1][intrxn.i+1] = 0;
      potent_type_map[nstyles][intrxn.i+1][intrxn.j+1] = potent_type_map[nstyles][intrxn.j+1][intrxn.i+1] = 1;
    }
  }
}
