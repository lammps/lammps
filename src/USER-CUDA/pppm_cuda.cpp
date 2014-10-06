/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

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
   Contributing authors: Roy Pollock (LLNL), Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "pppm_cuda.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "domain.h"
#include "fft3d_wrap_cuda.h"
#include "remap_wrap.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include <ctime> //crmadd
#include "cuda_wrapper_cu.h"
#include "pppm_cuda_cu.h"
#include "user_cuda.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXORDER 7
#define OFFSET 4096
#define SMALL 0.00001
#define LARGE 10000.0
#define EPS_HOC 1.0e-7


void printArray(double* data,int nx, int ny, int nz)
{
  for(int i=0;i<nz;i++)
  for(int j=0;j<ny;j++)
  {
          printf("%i %i\n",i,j);
          for(int k=0;k<nx;k++)
          printf("%e ",data[2*(i*ny*nx+j*nx+k)]);
          printf("\n\n");
  }
}
void printArray(double*** data,int nx, int ny, int nz)
{
  for(int i=0;i<nx;i++)
  for(int j=0;j<ny;j++)
  {
          printf("%i %i\n",i,j);
          for(int k=0;k<nz;k++)
          printf("%e ",data[i][j][k]);
          printf("\n\n");
  }
}
/* ---------------------------------------------------------------------- */

PPPMCuda::PPPMCuda(LAMMPS *lmp, int narg, char **arg) : 
  PPPMOld(lmp, (narg==2?1:narg), arg)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  if ((narg > 3)||(narg<1)) error->all(FLERR,"Illegal kspace_style pppm/cuda command");
  #ifndef FFT_CUFFT
  error->all(FLERR,"Using kspace_style pppm/cuda without cufft is not possible. Compile with cufft=1 to include cufft. Aborting.");
  #endif

  triclinic_support = 0;
  accuracy_relative = fabs(force->numeric(FLERR,arg[0]));

  nfactors = 3;
  factors = new int[nfactors];
  factors[0] = 2;
  factors[1] = 3;
  factors[2] = 5;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  density_brick = vdx_brick = vdy_brick = vdz_brick = vdx_brick_tmp = NULL;
  density_fft = NULL;
  greensfn = NULL;
  work1 = work2 = NULL;
  vg = NULL;
  fkx = fky = fkz = NULL;
  buf1 = buf2 = NULL;

  gf_b = NULL;
  rho1d = rho_coeff = NULL;

  fft1c = fft2c = NULL;
  remap = NULL;

  density_brick_int=NULL;
  density_intScale=1000000;
  cu_vdx_brick = cu_vdy_brick = cu_vdz_brick = NULL;
  cu_density_brick = NULL;
  cu_density_brick_int = NULL;
  cu_density_fft = NULL;
  cu_energy=NULL;
  cu_greensfn = NULL;
  cu_work1 = cu_work2 = cu_work3 = NULL;
  cu_vg = NULL;
  cu_fkx = cu_fky = cu_fkz = NULL;

  cu_flag = NULL;
  cu_debugdata = NULL;
  cu_rho_coeff = NULL;
  cu_virial = NULL;

  cu_gf_b = NULL;

  cu_slabbuf = NULL;
  slabbuf = NULL;

  nmax = 0;
  part2grid = NULL;
  cu_part2grid = NULL;
  adev_data_array=NULL;
  poissontime=0;
  old_nmax=0;
  cu_pppm_grid_n=NULL;
  cu_pppm_grid_ids=NULL;

  pppm_grid_nmax=0;
  pppm2partgrid=new int[3];
  pppm_grid=new int[3];
  firstpass=true;
  scale = 1.0;
}


/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

PPPMCuda::~PPPMCuda()
{
  delete [] slabbuf;
  delete cu_slabbuf;

  delete [] factors;
  factors=NULL;
  deallocate();
  delete cu_part2grid;
  cu_part2grid=NULL;
  memory->destroy(part2grid);
  part2grid = NULL;
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void PPPMCuda::init()
{

  cuda->shared_data.pppm.cudable_force=1;

    //if(cuda->finished_run) {PPPM::init(); return;}

  if (me == 0) {
    if (screen) fprintf(screen,"PPPMCuda initialization ...\n");
    if (logfile) fprintf(logfile,"PPPMCuda initialization ...\n");
  }

  // error check

  if (domain->dimension == 2) error->all(FLERR,"Cannot use PPPMCuda with 2d simulation");
  if (comm->style != 0) 
    error->universe_all(FLERR,"PPPMCuda can only currently be used with "
                        "comm_style brick");

  if (!atom->q_flag) error->all(FLERR,"Kspace style requires atom attribute q");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use nonperiodic boundaries with PPPMCuda");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab PPPMCuda");
  }

  if (order < 2 || order > MAXORDER) {
    char str[128];
    sprintf(str,"PPPMCuda order cannot be smaller than 2 or greater than %d",MAXORDER);
    error->all(FLERR,str);
  }
  // free all arrays previously allocated

  deallocate();

  // extract short-range Coulombic cutoff from pair style

  triclinic_check();

  if (force->pair == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  int itmp=0;
  double *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  if (p_cutoff == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;

  // if kspace is TIP4P, extract TIP4P params from pair style

  qdist = 0.0;

  if (strcmp(force->kspace_style,"pppm/tip4p") == 0) {
    if (force->pair == NULL)
      error->all(FLERR,"KSpace style is incompatible with Pair style");
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
    double theta = force->angle->equilibrium_angle(typeA);
    double blen = force->bond->equilibrium_distance(typeB);
    alpha = qdist / (2.0 * cos(0.5*theta) * blen);
  }

  // compute qsum & qsqsum and warn if not charge-neutral

  scale = 1.0;
  qqrd2e = force->qqrd2e;
  qsum_qsq(0);
  natoms_original = atom->natoms;

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;

  // setup FFT grid resolution and g_ewald
  // normally one iteration thru while loop is all that is required
  // if grid stencil extends beyond neighbor proc, reduce order and try again

  int iteration = 0;

  while (order > 1) {
    if (iteration && me == 0)
      error->warning(FLERR,"Reducing PPPMCuda order b/c stencil extends "
                     "beyond neighbor processor");
    iteration++;

    set_grid();

    if (nx_pppm >= OFFSET || ny_pppm >= OFFSET || nz_pppm >= OFFSET)
      error->all(FLERR,"PPPMCuda grid is too large");

    // global indices of PPPMCuda grid range from 0 to N-1
    // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
    //   global PPPMCuda grid that I own without ghost cells
    // for slab PPPMCuda, assign z grid as if it were not extended

    nxlo_in = comm->myloc[0]*nx_pppm / comm->procgrid[0];
    nxhi_in = (comm->myloc[0]+1)*nx_pppm / comm->procgrid[0] - 1;
    nylo_in = comm->myloc[1]*ny_pppm / comm->procgrid[1];
    nyhi_in = (comm->myloc[1]+1)*ny_pppm / comm->procgrid[1] - 1;
    nzlo_in = comm->myloc[2] *
      (static_cast<int> (nz_pppm/slab_volfactor)) / comm->procgrid[2];
    nzhi_in = (comm->myloc[2]+1) *
      (static_cast<int> (nz_pppm/slab_volfactor)) / comm->procgrid[2] - 1;

    // nlower,nupper = stencil size for mapping particles to PPPMCuda grid

    nlower = -(order-1)/2;
    nupper = order/2;

    // shift values for particle <-> grid mapping
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    if (order % 2) shift = OFFSET + 0.5;
    else shift = OFFSET;
    if (order % 2) shiftone = 0.0;
    else shiftone = 0.5;

    // nlo_out,nhi_out = lower/upper limits of the 3d sub-brick of
    //   global PPPMCuda grid that my particles can contribute charge to
    // effectively nlo_in,nhi_in + ghost cells
    // nlo,nhi = global coords of grid pt to "lower left" of smallest/largest
    //           position a particle in my box can be at
    // dist[3] = particle position bound = subbox + skin/2.0 + qdist
    //   qdist = offset due to TIP4P fictitious charge
    //   convert to triclinic if necessary
    // nlo_out,nhi_out = nlo,nhi + stencil size for particle mapping
    // for slab PPPMCuda, assign z grid as if it were not extended


    triclinic = domain->triclinic;
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
                            nx_pppm/xprd + shift) - OFFSET;
    nhi = static_cast<int> ((subhi[0]+dist[0]-boxlo[0]) *
                            nx_pppm/xprd + shift) - OFFSET;
    nxlo_out = nlo + nlower;
    nxhi_out = nhi + nupper;

    nlo = static_cast<int> ((sublo[1]-dist[1]-boxlo[1]) *
                            ny_pppm/yprd + shift) - OFFSET;
    nhi = static_cast<int> ((subhi[1]+dist[1]-boxlo[1]) *
                            ny_pppm/yprd + shift) - OFFSET;
    nylo_out = nlo + nlower;
    nyhi_out = nhi + nupper;

    nlo = static_cast<int> ((sublo[2]-dist[2]-boxlo[2]) *
                            nz_pppm/zprd_slab + shift) - OFFSET;
    nhi = static_cast<int> ((subhi[2]+dist[2]-boxlo[2]) *
                            nz_pppm/zprd_slab + shift) - OFFSET;
    nzlo_out = nlo + nlower;
    nzhi_out = nhi + nupper;

    // for slab PPPMCuda, change the grid boundary for processors at +z end
    //   to include the empty volume between periodically repeating slabs
    // for slab PPPMCuda, want charge data communicated from -z proc to +z proc,
    //   but not vice versa, also want field data communicated from +z proc to
    //   -z proc, but not vice versa
    // this is accomplished by nzhi_in = nzhi_out on +z end (no ghost cells)

    if (slabflag && ((comm->myloc[2]+1) == (comm->procgrid[2]))) {
      nzhi_in =  nz_pppm - 1;
      nzhi_out = nz_pppm - 1;
    }

    // nlo_ghost,nhi_ghost = # of planes I will recv from 6 directions
    //   that overlay domain I own
    // proc in that direction tells me via sendrecv()
    // if no neighbor proc, value is from self since I have ghosts regardless

    int nplanes;
    MPI_Status status;

    nplanes = nxlo_in - nxlo_out;
    if (comm->procneigh[0][0] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[0][0],0,
                   &nxhi_ghost,1,MPI_INT,comm->procneigh[0][1],0,
                   world,&status);
    else nxhi_ghost = nplanes;

    nplanes = nxhi_out - nxhi_in;
    if (comm->procneigh[0][1] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[0][1],0,
                   &nxlo_ghost,1,MPI_INT,comm->procneigh[0][0],
                   0,world,&status);
    else nxlo_ghost = nplanes;

    nplanes = nylo_in - nylo_out;
    if (comm->procneigh[1][0] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[1][0],0,
                   &nyhi_ghost,1,MPI_INT,comm->procneigh[1][1],0,
                   world,&status);
    else nyhi_ghost = nplanes;

    nplanes = nyhi_out - nyhi_in;
    if (comm->procneigh[1][1] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[1][1],0,
                   &nylo_ghost,1,MPI_INT,comm->procneigh[1][0],0,
                   world,&status);
    else nylo_ghost = nplanes;

    nplanes = nzlo_in - nzlo_out;
    if (comm->procneigh[2][0] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[2][0],0,
                   &nzhi_ghost,1,MPI_INT,comm->procneigh[2][1],0,
                   world,&status);
    else nzhi_ghost = nplanes;

    nplanes = nzhi_out - nzhi_in;
    if (comm->procneigh[2][1] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[2][1],0,
                   &nzlo_ghost,1,MPI_INT,comm->procneigh[2][0],0,
                   world,&status);
    else nzlo_ghost = nplanes;

    // test that ghost overlap is not bigger than my sub-domain

    int flag = 0;
    if (nxlo_ghost > nxhi_in-nxlo_in+1) flag = 1;
    if (nxhi_ghost > nxhi_in-nxlo_in+1) flag = 1;
    if (nylo_ghost > nyhi_in-nylo_in+1) flag = 1;
    if (nyhi_ghost > nyhi_in-nylo_in+1) flag = 1;
    if (nzlo_ghost > nzhi_in-nzlo_in+1) flag = 1;
    if (nzhi_ghost > nzhi_in-nzlo_in+1) flag = 1;

    int flag_all;
    MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);

    if (flag_all == 0) break;
    order--;
  }

  if (order == 0) error->all(FLERR,"PPPMCuda order has been reduced to 0");



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

  // PPPMCuda grid for this proc, including ghosts

  ngrid = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
    (nzhi_out-nzlo_out+1);

  // FFT arrays on this proc, without ghosts
  // nfft = FFT points in FFT decomposition on this proc
  // nfft_brick = FFT points in 3d brick-decomposition on this proc
  // nfft_both = greater of 2 values

  nfft = (nxhi_fft-nxlo_fft+1) * (nyhi_fft-nylo_fft+1) *
    (nzhi_fft-nzlo_fft+1);
  int nfft_brick = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1) *
    (nzhi_in-nzlo_in+1);
  nfft_both = MAX(nfft,nfft_brick);

  // buffer space for use in brick2fft and fillbrick
  // idel = max # of ghost planes to send or recv in +/- dir of each dim
  // nx,ny,nz = owned planes (including ghosts) in each dim
  // nxx,nyy,nzz = max # of grid cells to send in each dim
  // nbuf = max in any dim, augment by 3x for components of vd_xyz in fillbrick

  int idelx,idely,idelz,nx,ny,nz,nxx,nyy,nzz;

  idelx = MAX(nxlo_ghost,nxhi_ghost);
  idelx = MAX(idelx,nxhi_out-nxhi_in);
  idelx = MAX(idelx,nxlo_in-nxlo_out);

  idely = MAX(nylo_ghost,nyhi_ghost);
  idely = MAX(idely,nyhi_out-nyhi_in);
  idely = MAX(idely,nylo_in-nylo_out);

  idelz = MAX(nzlo_ghost,nzhi_ghost);
  idelz = MAX(idelz,nzhi_out-nzhi_in);
  idelz = MAX(idelz,nzlo_in-nzlo_out);

  nx = nxhi_out - nxlo_out + 1;
  ny = nyhi_out - nylo_out + 1;
  nz = nzhi_out - nzlo_out + 1;

  nxx = idelx * ny * nz;
  nyy = idely * nx * nz;
  nzz = idelz * nx * ny;

  nbuf = MAX(nxx,nyy);
  nbuf = MAX(nbuf,nzz);
  nbuf *= 3;

  // print stats

  int ngrid_max,nfft_both_max,nbuf_max;
  MPI_Allreduce(&ngrid,&ngrid_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nfft_both,&nfft_both_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nbuf,&nbuf_max,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  brick FFT buffer size/proc = %d %d %d\n",
                        ngrid_max,nfft_both_max,nbuf_max);
    if (logfile) fprintf(logfile,"  brick FFT buffer size/proc = %d %d %d\n",
                         ngrid_max,nfft_both_max,nbuf_max);
  }
  cuda_shared_pppm* ap=&(cuda->shared_data.pppm);

   ap->density_intScale=density_intScale;
   ap->nxlo_in=nxlo_in;
   ap->nxhi_in=nxhi_in;
   ap->nxlo_out=nxlo_out;
   ap->nxhi_out=nxhi_out;
   ap->nylo_in=nylo_in;
   ap->nyhi_in=nyhi_in;
   ap->nylo_out=nylo_out;
   ap->nyhi_out=nyhi_out;
   ap->nzlo_in=nzlo_in;
   ap->nzhi_in=nzhi_in;
   ap->nzlo_out=nzlo_out;
   ap->nzhi_out=nzhi_out;
   ap->nxlo_in=nxlo_fft;
   ap->nxhi_in=nxhi_fft;
   ap->nylo_in=nylo_fft;
   ap->nyhi_in=nyhi_fft;
   ap->nzlo_in=nzlo_fft;
   ap->nzhi_in=nzhi_fft;
   ap->nx_pppm=nx_pppm;
   ap->ny_pppm=ny_pppm;
   ap->nz_pppm=nz_pppm;
   ap->qqrd2e=qqrd2e;
   ap->order=order;
   ap->nmax=nmax;
   ap->nlocal=atom->nlocal;
   ap->delxinv=delxinv;
   ap->delyinv=delyinv;
   ap->delzinv=delzinv;
   ap->nlower=nlower;
   ap->nupper=nupper;
   ap->shiftone=shiftone;

  // allocate K-space dependent memory


  allocate();

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients

  compute_gf_denom();
  compute_rho_coeff();
}

/* ----------------------------------------------------------------------
   adjust PPPMCuda coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void PPPMCuda::setup()
{
  double *prd;
  cu_gf_b->upload();
  // volume-dependent factors
  // adjust z dimension for 2d slab PPPMCuda
  // z dimension for 3d PPPMCuda is zprd since slab_volfactor = 1.0

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  delxinv = nx_pppm/xprd;
  delyinv = ny_pppm/yprd;
  delzinv = nz_pppm/zprd_slab;

  delvolinv = delxinv*delyinv*delzinv;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  // fkx,fky,fkz for my FFT grid pts
  Cuda_PPPM_Setup_fkxyz_vg(nx_pppm, ny_pppm,nz_pppm,unitkx,unitky,unitkz,g_ewald);



  // modified (Hockney-Eastwood) Coulomb Green's function

  int nbx = static_cast<int> ((g_ewald*xprd/(MY_PI*nx_pppm)) *
                              pow(-log(EPS_HOC),0.25));
  int nby = static_cast<int> ((g_ewald*yprd/(MY_PI*ny_pppm)) *
                              pow(-log(EPS_HOC),0.25));
  int nbz = static_cast<int> ((g_ewald*zprd_slab/(MY_PI*nz_pppm)) *
                              pow(-log(EPS_HOC),0.25));
  Cuda_PPPM_setup_greensfn(nx_pppm,ny_pppm,nz_pppm,unitkx,unitky,unitkz,g_ewald,
nbx,nby,nbz,xprd,yprd,zprd_slab);


#ifdef FFT_CUFFT
  cu_vdx_brick->upload();
  cu_vdy_brick->upload();
  cu_vdz_brick->upload();
#endif
  cu_rho_coeff->upload();
  cu_density_brick->memset_device(0);
  pppm_device_init_setup(&cuda->shared_data,shiftone,delxinv,delyinv,delzinv,nlower,nupper);
}

/* ----------------------------------------------------------------------
   compute the PPPMCuda long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPMCuda::compute(int eflag, int vflag)
{
  cuda_shared_atom*   cu_atom   = & cuda->shared_data.atom;

  int i;
  my_times starttime;
  my_times endtime;
  my_times starttotal;
  my_times endtotal;
  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // extend size of per-atom arrays if necessary

  if ((cu_atom->update_nmax)||(old_nmax==0)) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm:part2grid");
         delete cu_part2grid;
         delete [] adev_data_array;
         adev_data_array=new dev_array[1];
         cu_part2grid = new cCudaData<int  , int   , yx > ((int*)part2grid,adev_data_array, nmax,3);

          pppm_device_update(&cuda->shared_data,cu_part2grid->dev_data(),atom->nlocal,atom->nmax);
    old_nmax=nmax;
  }
  if(cu_atom->update_nlocal) {pppm_update_nlocal(cu_atom->nlocal);}

  energy = 0.0;
  if (vflag)
  {
          for (i = 0; i < 6; i++) virial[i] = 0.0;
          cu_virial->memset_device(0);
  }
  if(eflag) cu_energy->memset_device(0);
  my_gettime(CLOCK_REALTIME,&starttotal);

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid


  my_gettime(CLOCK_REALTIME,&starttime);

  particle_map();

  my_gettime(CLOCK_REALTIME,&endtime);
  cuda->shared_data.cuda_timings.pppm_particle_map+=(endtime.tv_sec-starttime.tv_sec+1.0*(endtime.tv_nsec-starttime.tv_nsec)/1000000000);

  //cu_part2grid->download();
  my_gettime(CLOCK_REALTIME,&starttime);
  make_rho();
  my_gettime(CLOCK_REALTIME,&endtime);
  cuda->shared_data.cuda_timings.pppm_make_rho+=(endtime.tv_sec-starttime.tv_sec+1.0*(endtime.tv_nsec-starttime.tv_nsec)/1000000000);

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  int nprocs=comm->nprocs;

  my_gettime(CLOCK_REALTIME,&starttime);

  if(nprocs>1)
  {
    cu_density_brick->download();
    brick2fft();
  }
  else
  {
     #ifdef FFT_CUFFT
     pppm_initfftdata(&cuda->shared_data,(PPPM_CFLOAT*)cu_density_brick->dev_data(),(FFT_CFLOAT*)cu_work2->dev_data());
     #endif
  }

  my_gettime(CLOCK_REALTIME,&endtime);
  cuda->shared_data.cuda_timings.pppm_brick2fft+=(endtime.tv_sec-starttime.tv_sec+1.0*(endtime.tv_nsec-starttime.tv_nsec)/1000000000);

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition

  my_gettime(CLOCK_REALTIME,&starttime);
  poisson(eflag,vflag);
  my_gettime(CLOCK_REALTIME,&endtime);
  cuda->shared_data.cuda_timings.pppm_poisson+=(endtime.tv_sec-starttime.tv_sec+1.0*(endtime.tv_nsec-starttime.tv_nsec)/1000000000);

  // all procs communicate E-field values to fill ghost cells
  //   surrounding their 3d bricks

  // not necessary since all the calculations are done on one proc

  // calculate the force on my particles


  my_gettime(CLOCK_REALTIME,&starttime);
  fieldforce();
  my_gettime(CLOCK_REALTIME,&endtime);
  cuda->shared_data.cuda_timings.pppm_fieldforce+=(endtime.tv_sec-starttime.tv_sec+1.0*(endtime.tv_nsec-starttime.tv_nsec)/1000000000);

  // sum energy across procs and add in volume-dependent term
  // reset qsum and qsqsum if atom count has changed

  my_gettime(CLOCK_REALTIME,&endtotal);
  cuda->shared_data.cuda_timings.pppm_compute+=(endtotal.tv_sec-starttotal.tv_sec+1.0*(endtotal.tv_nsec-starttotal.tv_nsec)/1000000000);

  if (eflag) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    if (atom->natoms != natoms_original) {
      qsum_qsq(0);
      natoms_original = atom->natoms;
    }

    energy *= 0.5*volume;
    energy -= g_ewald*qsqsum/1.772453851 +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qqrd2e;
  }

  // sum virial across procs

  if (vflag) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qqrd2e*volume*virial_all[i];
  }

  // 2d slab correction

  if (slabflag) slabcorr(eflag);

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);

  if(firstpass) firstpass=false;
}


/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */


void PPPMCuda::allocate()
{

  struct dev_array* dev_tmp=new struct dev_array[20];
  int n_cudata=0;


  memory->create3d_offset(density_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:density_brick");
  memory->create3d_offset(density_brick_int,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:density_brick_int");


  cu_density_brick = new cCudaData<double, PPPM_CFLOAT, x> ((double*) &(density_brick[nzlo_out][nylo_out][nxlo_out]), & (dev_tmp[n_cudata++]),
                                     (nzhi_out-nzlo_out+1)*(nyhi_out-nylo_out+1)*(nxhi_out-nxlo_out+1));

  cu_density_brick_int = new cCudaData<int, int, x> ((int*) &(density_brick_int[nzlo_out][nylo_out][nxlo_out]), & (dev_tmp[n_cudata++]),
                                     (nzhi_out-nzlo_out+1)*(nyhi_out-nylo_out+1)*(nxhi_out-nxlo_out+1));

  memory->create3d_offset(vdx_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:vdx_brick");
  memory->create3d_offset(vdx_brick_tmp,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:vdx_brick_tmp");

  cu_vdx_brick = new cCudaData<double, PPPM_CFLOAT, x> ((double*) &(vdx_brick[nzlo_out][nylo_out][nxlo_out]), & (dev_tmp[n_cudata++]),
                                     (nzhi_out-nzlo_out+1)*(nyhi_out-nylo_out+1)*(nxhi_out-nxlo_out+1));

  memory->create3d_offset(vdy_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:vdy_brick");
  cu_vdy_brick = new cCudaData<double, PPPM_CFLOAT, x> ((double*) &(vdy_brick[nzlo_out][nylo_out][nxlo_out]), & (dev_tmp[n_cudata++]),
                                     (nzhi_out-nzlo_out+1)*(nyhi_out-nylo_out+1)*(nxhi_out-nxlo_out+1));

  memory->create3d_offset(vdz_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:vdz_brick");
  cu_vdz_brick = new cCudaData<double, PPPM_CFLOAT, x> ((double*) &(vdz_brick[nzlo_out][nylo_out][nxlo_out]), & (dev_tmp[n_cudata++]),
                                     (nzhi_out-nzlo_out+1)*(nyhi_out-nylo_out+1)*(nxhi_out-nxlo_out+1));

  memory->create(density_fft,nfft_both,"pppm:density_fft");

  cu_density_fft = new cCudaData<double, PPPM_CFLOAT, x> (density_fft, & (dev_tmp[n_cudata++]),nfft_both);

  cu_energy = new cCudaData<double, ENERGY_CFLOAT, x> (NULL, &(dev_tmp[n_cudata++]),ny_pppm*nz_pppm);
  cu_virial = new cCudaData<double, ENERGY_CFLOAT, x> (NULL, &(dev_tmp[n_cudata++]),ny_pppm*nz_pppm*6);

  memory->create(greensfn,nfft_both,"pppm:greensfn");
  cu_greensfn = new cCudaData<double, PPPM_CFLOAT, x> (greensfn, & (dev_tmp[n_cudata++]) , nx_pppm*ny_pppm*nz_pppm);

  memory->create(work1,2*nx_pppm*ny_pppm*nz_pppm,"pppm:work1");
  memory->create(work2,2*nx_pppm*ny_pppm*nz_pppm,"pppm:work2");
  memory->create(work3,2*nx_pppm*ny_pppm*nz_pppm,"pppm:work3");

  cu_work1 = new cCudaData<double, FFT_CFLOAT, x> (work1, & (dev_tmp[n_cudata++]) , 2*nx_pppm*ny_pppm*nz_pppm);
  cu_work2 = new cCudaData<double, FFT_CFLOAT, x> (work2, & (dev_tmp[n_cudata++]) , 2*nx_pppm*ny_pppm*nz_pppm);
  cu_work3 = new cCudaData<double, FFT_CFLOAT, x> (work3, & (dev_tmp[n_cudata++]) , 2*nx_pppm*ny_pppm*nz_pppm);


  memory->create(fkx,nx_pppm,"pppmcuda:fkx");
  cu_fkx = new cCudaData<double, PPPM_CFLOAT, x> (fkx, & (dev_tmp[n_cudata++]) , nx_pppm);
  memory->create(fky,ny_pppm,"pppmcuda:fky");
  cu_fky = new cCudaData<double, PPPM_CFLOAT, x> (fky, & (dev_tmp[n_cudata++]) , ny_pppm);
  memory->create(fkz,nz_pppm,"pppmcuda:fkz");
  cu_fkz = new cCudaData<double, PPPM_CFLOAT, x> (fkz, & (dev_tmp[n_cudata++]) , nz_pppm);

  memory->create(vg,nfft_both,6,"pppm:vg");

  cu_vg = new cCudaData<double, PPPM_CFLOAT, xy> ((double*)vg, & (dev_tmp[n_cudata++]) , nfft_both,6);

  memory->create(buf1,nbuf,"pppm:buf1");
  memory->create(buf2,nbuf,"pppm:buf2");


  // summation coeffs


  gf_b = new double[order];
  cu_gf_b = new cCudaData<double,PPPM_CFLOAT,x> (gf_b, &(dev_tmp[n_cudata++]) , order);
  memory->create2d_offset(rho1d,3,-order/2,order/2,"pppm:rho1d");
  memory->create2d_offset(rho_coeff,order,(1-order)/2,order/2,"pppm:rho_coeff");

  cu_rho_coeff = new cCudaData<double, PPPM_CFLOAT, x> ((double*) &(rho_coeff[0][(1-order)/2]), & (dev_tmp[n_cudata++]) , order*(order/2-(1-order)/2+1));

  debugdata=new PPPM_CFLOAT[100];
  cu_debugdata = new cCudaData<PPPM_CFLOAT, PPPM_CFLOAT, x> (debugdata,& (dev_tmp[n_cudata++]),100);
  cu_flag = new cCudaData<int, int, x> (&global_flag,& (dev_tmp[n_cudata++]),3);

  // create 2 FFTs and a Remap
  // 1st FFT keeps data in FFT decompostion
  // 2nd FFT returns data in 3d brick decomposition
  // remap takes data from 3d brick to FFT decomposition

  int tmp;




  fft1c = new FFT3dCuda(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   0,0,&tmp,true);

  fft2c = new FFT3dCuda(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                   0,0,&tmp,false);


#ifdef FFT_CUFFT
  fft1c->set_cudata(cu_work2->dev_data(),cu_work1->dev_data());
  fft2c->set_cudata(cu_work2->dev_data(),cu_work3->dev_data());
#endif

  remap = new Remap(lmp,world,
                    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                    nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                    1,0,0,2,0);


pppm_device_init(cu_density_brick->dev_data(), cu_vdx_brick->dev_data(), cu_vdy_brick->dev_data(), cu_vdz_brick->dev_data(), cu_density_fft->dev_data(),cu_energy->dev_data(),cu_virial->dev_data()
            , cu_work1->dev_data(), cu_work2->dev_data(), cu_work3->dev_data(), cu_greensfn->dev_data(), cu_fkx->dev_data(), cu_fky->dev_data(), cu_fkz->dev_data(), cu_vg->dev_data()
            ,nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,nx_pppm,ny_pppm,nz_pppm
            ,nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,cu_gf_b->dev_data()
            ,qqrd2e,order,cu_rho_coeff->dev_data(),cu_debugdata->dev_data(),cu_density_brick_int->dev_data(),slabflag
         );
}



/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
 ---------------------------------------------------------------------- */

void PPPMCuda::deallocate()
{
  memory->destroy3d_offset(density_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(vdx_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(vdy_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(vdz_brick,nzlo_out,nylo_out,nxlo_out);

  density_brick = vdx_brick = vdy_brick = vdz_brick = NULL;

  memory->destroy(density_fft);
  memory->destroy(greensfn);
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(vg);

  density_fft = NULL;
  greensfn = NULL;
  work1 = NULL;
  work2 = NULL;
  vg = NULL;

  memory->destroy(fkx);
  memory->destroy(fky);
  memory->destroy(fkz);

  fkx = NULL;
  fky = NULL;
  fkz = NULL;

  delete cu_density_brick;
  delete cu_density_brick_int;
  delete cu_vdx_brick;
  delete cu_vdy_brick;
  delete cu_vdz_brick;
  delete cu_density_fft;
  delete cu_energy;
  delete cu_virial;
#ifdef FFT_CUFFT
  delete cu_greensfn;
  delete cu_gf_b;
  delete cu_vg;
  delete cu_work1;
  delete cu_work2;
  delete cu_work3;
  delete cu_fkx;
  delete cu_fky;
  delete cu_fkz;
#endif

  delete cu_flag;
  delete cu_debugdata;
  delete cu_rho_coeff;


  cu_vdx_brick = cu_vdy_brick = cu_vdz_brick = NULL;
  cu_density_brick = NULL;
  cu_density_brick_int = NULL;
  cu_density_fft = NULL;
  cu_energy=NULL;
  cu_virial=NULL;
#ifdef FFT_CUFFT
  cu_greensfn = NULL;
  cu_gf_b = NULL;
  cu_work1 = cu_work2 = cu_work3 = NULL;
  cu_vg = NULL;
  cu_fkx = cu_fky = cu_fkz = NULL;
#endif

  cu_flag = NULL;
  cu_debugdata = NULL;
  cu_rho_coeff = NULL;
  cu_part2grid = NULL;

  memory->destroy(buf1);
  memory->destroy(buf2);

  delete [] gf_b;
  gf_b = NULL;
  memory->destroy2d_offset(rho1d,-order/2); rho1d = NULL;
  memory->destroy2d_offset(rho_coeff,(1-order)/2); rho_coeff = NULL;

  delete fft1c;
  fft1c = NULL;

  delete fft2c;
  fft2c = NULL;
  delete remap;
  remap = NULL;
  buf1 = NULL;
  buf2 = NULL;
}

/* ----------------------------------------------------------------------
   set size of FFT grid (nx,ny,nz_pppm) and g_ewald
-------------------------------------------------------------------------*/

void PPPMCuda::set_grid()
{
  // see JCP 109, pg 7698 for derivation of coefficients
  // higher order coefficients may be computed if needed

  double **acons;
  memory->create(acons,8,7,"pppm:acons");

  acons[1][0] = 2.0 / 3.0;
  acons[2][0] = 1.0 / 50.0;
  acons[2][1] = 5.0 / 294.0;
  acons[3][0] = 1.0 / 588.0;
  acons[3][1] = 7.0 / 1440.0;
  acons[3][2] = 21.0 / 3872.0;
  acons[4][0] = 1.0 / 4320.0;
  acons[4][1] = 3.0 / 1936.0;
  acons[4][2] = 7601.0 / 2271360.0;
  acons[4][3] = 143.0 / 28800.0;
  acons[5][0] = 1.0 / 23232.0;
  acons[5][1] = 7601.0 / 13628160.0;
  acons[5][2] = 143.0 / 69120.0;
  acons[5][3] = 517231.0 / 106536960.0;
  acons[5][4] = 106640677.0 / 11737571328.0;
  acons[6][0] = 691.0 / 68140800.0;
  acons[6][1] = 13.0 / 57600.0;
  acons[6][2] = 47021.0 / 35512320.0;
  acons[6][3] = 9694607.0 / 2095994880.0;
  acons[6][4] = 733191589.0 / 59609088000.0;
  acons[6][5] = 326190917.0 / 11700633600.0;
  acons[7][0] = 1.0 / 345600.0;
  acons[7][1] = 3617.0 / 35512320.0;
  acons[7][2] = 745739.0 / 838397952.0;
  acons[7][3] = 56399353.0 / 12773376000.0;
  acons[7][4] = 25091609.0 / 1560084480.0;
  acons[7][5] = 1755948832039.0 / 36229939200000.0;
  acons[7][6] = 4887769399.0 / 37838389248.0;

  bigint natoms = atom->natoms;

  // use xprd,yprd,zprd even if triclinic so grid size is the same
  // adjust z dimension for 2d slab PPPMCuda
  // 3d PPPMCuda just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;

  // make initial g_ewald estimate
  // based on desired error and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  double h_x,h_y,h_z;

  if (!gewaldflag)
    g_ewald = sqrt(-log(accuracy*sqrt(natoms*cutoff*xprd*yprd*zprd) /
                        (2.0*q2))) / cutoff;

  // set optimal nx_pppm,ny_pppm,nz_pppm based on order and precision
  // nz_pppm uses extended zprd_slab instead of zprd
  // h = 1/g_ewald is upper bound on h such that h*g_ewald <= 1
  // reduce it until precision target is met

  if (!gridflag) {
    double err;
    h_x = h_y = h_z = 1/g_ewald;

    nx_pppm = static_cast<int> (xprd/h_x + 1);
    ny_pppm = static_cast<int> (yprd/h_y + 1);
    nz_pppm = static_cast<int> (zprd_slab/h_z + 1);

    err = rms(h_x,xprd,natoms,q2,acons);
    while (err > accuracy) {
      err = rms(h_x,xprd,natoms,q2,acons);
      nx_pppm++;
      h_x = xprd/nx_pppm;
    }

    err = rms(h_y,yprd,natoms,q2,acons);
    while (err > accuracy) {
      err = rms(h_y,yprd,natoms,q2,acons);
      ny_pppm++;
      h_y = yprd/ny_pppm;
    }

    err = rms(h_z,zprd_slab,natoms,q2,acons);
    while (err > accuracy) {
      err = rms(h_z,zprd_slab,natoms,q2,acons);
      nz_pppm++;
      h_z = zprd_slab/nz_pppm;
    }
  }

  // boost grid size until it is factorable

  while (!factorable(nx_pppm)) nx_pppm++;
  while (!factorable(ny_pppm)) ny_pppm++;
  while (!factorable(nz_pppm)) nz_pppm++;


  // adjust g_ewald for new grid size

  h_x = xprd/nx_pppm;
  h_y = yprd/ny_pppm;
  h_z = zprd_slab/nz_pppm;

  if (!gewaldflag) {
    double gew1,gew2,dgew,f,fmid,hmin,rtb;
    int ncount;

    gew1 = 0.0;
    g_ewald = gew1;
    f = diffpr(h_x,h_y,h_z,q2,acons);

    hmin = MIN(h_x,MIN(h_y,h_z));
    gew2 = 10/hmin;
    g_ewald = gew2;
    fmid = diffpr(h_x,h_y,h_z,q2,acons);

    if (f*fmid >= 0.0) error->all(FLERR,"Cannot compute PPPMCuda G");
    rtb = f < 0.0 ? (dgew=gew2-gew1,gew1) : (dgew=gew1-gew2,gew2);
    ncount = 0;
    while (fabs(dgew) > SMALL && fmid != 0.0) {
      dgew *= 0.5;
      g_ewald = rtb + dgew;
      fmid = diffpr(h_x,h_y,h_z,q2,acons);
      if (fmid <= 0.0) rtb = g_ewald;
      ncount++;
      if (ncount > LARGE) error->all(FLERR,"Cannot compute PPPMCuda G");
    }
  }

  // final RMS precision

  double lprx = rms(h_x,xprd,natoms,q2,acons);
  double lpry = rms(h_y,yprd,natoms,q2,acons);
  double lprz = rms(h_z,zprd_slab,natoms,q2,acons);
  double lpr = sqrt(lprx*lprx + lpry*lpry + lprz*lprz) / sqrt(3.0);
  double spr = 2.0*q2 * exp(-g_ewald*g_ewald*cutoff*cutoff) /
    sqrt(natoms*cutoff*xprd*yprd*zprd_slab);

  // free local memory

  memory->destroy(acons);

  // print info

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  G vector = %g\n",g_ewald);
      fprintf(screen,"  grid = %d %d %d\n",nx_pppm,ny_pppm,nz_pppm);
      fprintf(screen,"  stencil order = %d\n",order);
      fprintf(screen,"  absolute RMS force accuracy = %g\n",MAX(lpr,spr));
      fprintf(screen,"  relative force accuracy = %g\n",
              MAX(lpr,spr)/two_charge_force);
    }
    if (logfile) {
      fprintf(logfile,"  G vector = %g\n",g_ewald);
      fprintf(logfile,"  grid = %d %d %d\n",nx_pppm,ny_pppm,nz_pppm);
      fprintf(logfile,"  stencil order = %d\n",order);
      fprintf(logfile,"  absolute RMS force accuracy = %g\n",MAX(lpr,spr));
      fprintf(logfile,"  relative force accuracy = %g\n",
              MAX(lpr,spr)/two_charge_force);
    }
  }
}


/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */


void PPPMCuda::particle_map()
{
  MYDBG(printf("# CUDA PPPMCuda::particle_map() ... start\n");)
  int flag = 0;

    cu_flag->memset_device(0);
    flag=cuda_particle_map(&cuda->shared_data,cu_flag->dev_data());
    if(flag)
    {
      cu_debugdata->download();
      printf("Out of range atom: ");
       printf("ID: %i ",atom->tag[int(debugdata[0])]);
       printf("x: %e ",debugdata[7]);
       printf("y: %e ",debugdata[8]);
       printf("z: %e ",debugdata[9]);
       printf("nx: %e ",debugdata[4]);
       printf("ny: %e ",debugdata[5]);

      printf("\n");
      //printf("debugdata: cpu: %e %e %e %i\n",boxlo[0],boxlo[1],boxlo[2],atom->nlocal);
      cuda->cu_x->download();
            int nx,ny,nz;

            double **x = atom->x;
      int nlocal = atom->nlocal;
            for (int i = 0; i < nlocal; i++) {
        nx = static_cast<int> ((x[i][0]-boxlo[0])*delxinv+shift) - OFFSET;
        ny = static_cast<int> ((x[i][1]-boxlo[1])*delyinv+shift) - OFFSET;
              nz = static_cast<int> ((x[i][2]-boxlo[2])*delzinv+shift) - OFFSET;

            if(i==1203)printf("Outside Atom: %i %e %e %e (%i %i %i)\n",i,x[i][0],x[i][1],x[i][2],nx,ny,nz);
            if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
                ny+nlower < nylo_out || ny+nupper > nyhi_out ||
                nz+nlower < nzlo_out || nz+nupper > nzhi_out || i==1203) {printf("Outside Atom: %i %e %e %e (%i %i %i)\n",i,x[i][0],x[i][1],x[i][2],nx,ny,nz); }
            }

    }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Out of range atoms - cannot compute PPPMCuda!");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */


void PPPMCuda::make_rho()
{
    cuda_make_rho(&cuda->shared_data,cu_flag->dev_data(),&density_intScale,nxhi_out,nxlo_out,nyhi_out,nylo_out,nzhi_out,nzlo_out,cu_density_brick->dev_data(),cu_density_brick_int->dev_data());
}


/* ----------------------------------------------------------------------
   FFT-based Poisson solver
------------------------------------------------------------------------- */
void PPPMCuda::poisson(int eflag, int vflag)
{

#ifndef FFT_CUFFT
    PPPMOld::poisson();
    return;
#endif
#ifdef FFT_CUFFT
  my_times starttime;
  my_times endtime;


  my_gettime(CLOCK_REALTIME,&starttime);
  fft1c->compute(density_fft,work1,1);

  my_gettime(CLOCK_REALTIME,&endtime);
  poissontime+=(endtime.tv_sec-starttime.tv_sec+1.0*(endtime.tv_nsec-starttime.tv_nsec)/1000000000);



  if (eflag || vflag) {
    poisson_energy(nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,vflag);
    ENERGY_CFLOAT gpuvirial[6];
    energy+=sum_energy(cu_virial->dev_data(),cu_energy->dev_data(),nx_pppm,ny_pppm,nz_pppm,vflag,gpuvirial);
    if(vflag)
    {
      for(int j=0;j<6;j++) virial[j]+=gpuvirial[j];
    }
  }


  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  poisson_scale(nx_pppm,ny_pppm,nz_pppm);

   // compute gradients of V(r) in each of 3 dims by transformimg -ik*V(k)
  // FFT leaves data in 3d brick decomposition
  // copy it into inner portion of vdx,vdy,vdz arrays

  // x direction gradient


  poisson_xgrad(nx_pppm,ny_pppm,nz_pppm);


  my_gettime(CLOCK_REALTIME,&starttime);
  fft2c->compute(work2,work2,-1);
  my_gettime(CLOCK_REALTIME,&endtime);
  poissontime+=(endtime.tv_sec-starttime.tv_sec+1.0*(endtime.tv_nsec-starttime.tv_nsec)/1000000000);

  poisson_vdx_brick(nxhi_out,nxlo_out,nyhi_out,nylo_out,nzhi_out,nzlo_out,nx_pppm,ny_pppm,nz_pppm);


  // y direction gradient

  poisson_ygrad(nx_pppm,ny_pppm,nz_pppm);

  my_gettime(CLOCK_REALTIME,&starttime);
  fft2c->compute(work2,work2,-1);
  my_gettime(CLOCK_REALTIME,&endtime);
  poissontime+=(endtime.tv_sec-starttime.tv_sec+1.0*(endtime.tv_nsec-starttime.tv_nsec)/1000000000);

  poisson_vdy_brick(nxhi_out,nxlo_out,nyhi_out,nylo_out,nzhi_out,nzlo_out,nx_pppm,ny_pppm,nz_pppm);

  // z direction gradient

  poisson_zgrad(nx_pppm,ny_pppm,nz_pppm);

  my_gettime(CLOCK_REALTIME,&starttime);
  fft2c->compute(work2,work2,-1);
  my_gettime(CLOCK_REALTIME,&endtime);
  poissontime+=(endtime.tv_sec-starttime.tv_sec+1.0*(endtime.tv_nsec-starttime.tv_nsec)/1000000000);

  poisson_vdz_brick(nxhi_out,nxlo_out,nyhi_out,nylo_out,nzhi_out,nzlo_out,nx_pppm,ny_pppm,nz_pppm);
 #endif
}

/*----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
-------------------------------------------------------------------------*/

void PPPMCuda::fieldforce()
{
  cuda_fieldforce(& cuda->shared_data,cu_flag);
  return;
}

/* ----------------------------------------------------------------------
   perform and time the 4 FFTs required for N timesteps
------------------------------------------------------------------------- */

int PPPMCuda::timing_1d(int n, double &time1d)
{
  time1d = cuda->shared_data.cuda_timings.pppm_poisson/update->nsteps/4*n;
  return 4;
}

int PPPMCuda::timing_3d(int n, double &time3d)
{
  time3d = cuda->shared_data.cuda_timings.pppm_poisson/update->nsteps*n;
  return 4;
}

void PPPMCuda::slabcorr(int eflag)
{
  // compute local contribution to global dipole moment
  if(slabbuf==NULL)
  {
          slabbuf=new ENERGY_CFLOAT[(atom->nmax+31)/32];
          cu_slabbuf = new cCudaData<ENERGY_CFLOAT,ENERGY_CFLOAT, x> (slabbuf, (atom->nmax+31)/32);
  }
  if(unsigned((atom->nlocal+31)/32)*sizeof(ENERGY_CFLOAT)>=unsigned(cu_slabbuf->dev_size()))
  {
          delete [] slabbuf;
          delete cu_slabbuf;
          slabbuf=new ENERGY_CFLOAT[(atom->nmax+31)/32];
          cu_slabbuf = new cCudaData<ENERGY_CFLOAT,ENERGY_CFLOAT, x> (slabbuf, (atom->nmax+31)/32);
  }


  double dipole = cuda_slabcorr_energy(&cuda->shared_data,slabbuf,(ENERGY_CFLOAT*) cu_slabbuf->dev_data());

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // compute corrections

  double e_slabcorr = 2.0*MY_PI*dipole_all*dipole_all/volume;

  //if (eflag) energy += qqrd2e*scale * e_slabcorr;
  // need to add a correction to make non-neutral systems and per-atom energy translationally invariant
  if (eflag || fabs(qsum) > SMALL)
    error->all(FLERR,"Cannot (yet) use slab correction with kspace_style pppm/cuda for non-neutral systems or to get per-atom energy. Aborting.");

  double ffact = -4.0*MY_PI*dipole_all/volume;

  cuda_slabcorr_force(&cuda->shared_data,ffact);
}
