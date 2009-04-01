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
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "pppm.h"
#include "atom.h"
#include "comm.h"
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

#define MAXORDER 7
#define OFFSET 4096
#define SMALL 0.00001
#define LARGE 10000.0
#define EPS_HOC 1.0e-7

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PPPM::PPPM(LAMMPS *lmp, int narg, char **arg) : KSpace(lmp, narg, arg)
{
  if (narg != 1) error->all("Illegal kspace_style pppm command");

  precision = atof(arg[0]);
  PI = 4.0*atan(1.0);
  
  nfactors = 3;
  factors = new int[nfactors];
  factors[0] = 2;
  factors[1] = 3;
  factors[2] = 5;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  density_brick = vdx_brick = vdy_brick = vdz_brick = NULL;
  density_fft = NULL;
  greensfn = NULL;
  work1 = work2 = NULL;
  vg = NULL;
  fkx = fky = fkz = NULL;
  buf1 = buf2 = NULL;

  gf_b = NULL;
  rho1d = rho_coeff = NULL;

  fft1 = fft2 = NULL;
  remap = NULL;

  nmax = 0;
  part2grid = NULL;
}

/* ----------------------------------------------------------------------
   free all memory 
------------------------------------------------------------------------- */

PPPM::~PPPM()
{
  delete [] factors;
  deallocate();
  memory->destroy_2d_int_array(part2grid);
}

/* ----------------------------------------------------------------------
   called once before run 
------------------------------------------------------------------------- */

void PPPM::init()
{
  if (me == 0) {
    if (screen) fprintf(screen,"PPPM initialization ...\n");
    if (logfile) fprintf(logfile,"PPPM initialization ...\n");
  }

  // error check

  if (domain->triclinic)
    error->all("Cannot (yet) use PPPM with triclinic box");
  if (domain->dimension == 2) error->all("Cannot use PPPM with 2d simulation");

  if (!atom->q_flag) error->all("Kspace style requires atom attribute q");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all("Cannot use nonperiodic boundaries with PPPM");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 || 
	domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all("Incorrect boundaries with slab PPPM");
  }

  if (order > MAXORDER) {
    char str[128];
    sprintf(str,"PPPM order cannot be greater than %d",MAXORDER);
    error->all(str);
  }

  // free all arrays previously allocated

  deallocate();

  // extract short-range Coulombic cutoff from pair style

  qqrd2e = force->qqrd2e;

  if (force->pair == NULL)
    error->all("KSpace style is incompatible with Pair style");
  double *p_cutoff = (double *) force->pair->extract("cut_coul");
  if (p_cutoff == NULL)
    error->all("KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;

  // if kspace is TIP4P, extract TIP4P params from pair style

  qdist = 0.0;

  if (strcmp(force->kspace_style,"pppm/tip4p") == 0) {
    if (force->pair == NULL)
      error->all("KSpace style is incompatible with Pair style");
    double *p_qdist = (double *) force->pair->extract("qdist");
    int *p_typeO = (int *) force->pair->extract("typeO");
    int *p_typeH = (int *) force->pair->extract("typeH");
    int *p_typeA = (int *) force->pair->extract("typeA");
    int *p_typeB = (int *) force->pair->extract("typeB");
    if (!p_qdist || !p_typeO || !p_typeH || !p_typeA || !p_typeB)
      error->all("KSpace style is incompatible with Pair style");
    qdist = *p_qdist;
    typeO = *p_typeO;
    typeH = *p_typeH;
    int typeA = *p_typeA;
    int typeB = *p_typeB;

    if (force->angle == NULL || force->bond == NULL)
      error->all("Bond and angle potentials must be defined for TIP4P");
    double theta = force->angle->equilibrium_angle(typeA);
    double blen = force->bond->equilibrium_distance(typeB);
    alpha = qdist / (2.0 * cos(0.5*theta) * blen);
  }

  // compute qsum & qsqsum and warn if not charge-neutral

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

  if (qsqsum == 0.0)
    error->all("Cannot use kspace solver on system with no charge");
  if (fabs(qsum) > SMALL && me == 0) {
    char str[128];
    sprintf(str,"System is not charge neutral, net charge = %g",qsum);
    error->warning(str);
  }

  // setup FFT grid resolution and g_ewald
  // normally one iteration thru while loop is all that is required
  // if grid stencil extends beyond neighbor proc, reduce order and try again

  int iteration = 0;

  while (order > 0) {

    if (iteration && me == 0)
      error->warning("Reducing PPPM order b/c stencil extends "
		     "beyond neighbor processor");
    iteration++;

    set_grid();

    if (nx_pppm >= OFFSET || ny_pppm >= OFFSET || nz_pppm >= OFFSET)
      error->all("PPPM grid is too large");

    // global indices of PPPM grid range from 0 to N-1
    // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
    //   global PPPM grid that I own without ghost cells
    // for slab PPPM, assign z grid as if it were not extended

    nxlo_in = comm->myloc[0]*nx_pppm / comm->procgrid[0];
    nxhi_in = (comm->myloc[0]+1)*nx_pppm / comm->procgrid[0] - 1;
    nylo_in = comm->myloc[1]*ny_pppm / comm->procgrid[1];
    nyhi_in = (comm->myloc[1]+1)*ny_pppm / comm->procgrid[1] - 1;
    nzlo_in = comm->myloc[2] * 
      (static_cast<int> (nz_pppm/slab_volfactor)) / comm->procgrid[2];
    nzhi_in = (comm->myloc[2]+1) * 
      (static_cast<int> (nz_pppm/slab_volfactor)) / comm->procgrid[2] - 1;

    // nlower,nupper = stencil size for mapping particles to PPPM grid

    nlower = -(order-1)/2;
    nupper = order/2;

    // shift values for particle <-> grid mapping
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    if (order % 2) shift = OFFSET + 0.5;
    else shift = OFFSET;
    if (order % 2) shiftone = 0.0;
    else shiftone = 0.5;

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

    // for slab PPPM, change the grid boundary for processors at +z end
    //   to include the empty volume between periodically repeating slabs
    // for slab PPPM, want charge data communicated from -z proc to +z proc,
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

  if (order == 0) error->all("PPPM order has been reduced to 0");

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

  // PPPM grid for this proc, including ghosts

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

  // allocate K-space dependent memory

  allocate();

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients

  compute_gf_denom();
  compute_rho_coeff();
}

/* ----------------------------------------------------------------------
   adjust PPPM coeffs, called initially and whenever volume has changed 
------------------------------------------------------------------------- */

void PPPM::setup()
{
  int i,j,k,l,m,n;
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
    
  delxinv = nx_pppm/xprd;
  delyinv = ny_pppm/yprd;
  delzinv = nz_pppm/zprd_slab;

  delvolinv = delxinv*delyinv*delzinv;

  double unitkx = (2.0*PI/xprd);
  double unitky = (2.0*PI/yprd);
  double unitkz = (2.0*PI/zprd_slab);

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

  // modified (Hockney-Eastwood) Coulomb Green's function

  int nx,ny,nz,kper,lper,mper;
  double snx,sny,snz,snx2,sny2,snz2;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,dot1,dot2;
  double numerator,denominator;

  int nbx = static_cast<int> ((g_ewald*xprd/(PI*nx_pppm)) * 
			      pow(-log(EPS_HOC),0.25));
  int nby = static_cast<int> ((g_ewald*yprd/(PI*ny_pppm)) * 
			      pow(-log(EPS_HOC),0.25));
  int nbz = static_cast<int> ((g_ewald*zprd_slab/(PI*nz_pppm)) * 
			      pow(-log(EPS_HOC),0.25));

  double form = 1.0;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);
    snz = sin(0.5*unitkz*mper*zprd_slab/nz_pppm);
    snz2 = snz*snz;

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);
      sny = sin(0.5*unitky*lper*yprd/ny_pppm);
      sny2 = sny*sny;

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
	kper = k - nx_pppm*(2*k/nx_pppm);
	snx = sin(0.5*unitkx*kper*xprd/nx_pppm);
	snx2 = snx*snx;
      
	sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) + 
	  pow(unitkz*mper,2.0);

	if (sqk != 0.0) {
	  numerator = form*12.5663706/sqk;
	  denominator = gf_denom(snx2,sny2,snz2);  
	  sum1 = 0.0;
	  for (nx = -nbx; nx <= nbx; nx++) {
	    qx = unitkx*(kper+nx_pppm*nx);
	    sx = exp(-.25*pow(qx/g_ewald,2.0));
	    wx = 1.0;
	    argx = 0.5*qx*xprd/nx_pppm;
	    if (argx != 0.0) wx = pow(sin(argx)/argx,order);
	    for (ny = -nby; ny <= nby; ny++) {
	      qy = unitky*(lper+ny_pppm*ny);
	      sy = exp(-.25*pow(qy/g_ewald,2.0));
	      wy = 1.0;
	      argy = 0.5*qy*yprd/ny_pppm;
	      if (argy != 0.0) wy = pow(sin(argy)/argy,order);
	      for (nz = -nbz; nz <= nbz; nz++) {
		qz = unitkz*(mper+nz_pppm*nz);
		sz = exp(-.25*pow(qz/g_ewald,2.0));
		wz = 1.0;
		argz = 0.5*qz*zprd_slab/nz_pppm;
		if (argz != 0.0) wz = pow(sin(argz)/argz,order);

		dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
		dot2 = qx*qx+qy*qy+qz*qz;
		sum1 += (dot1/dot2) * sx*sy*sz * pow(wx*wy*wz,2.0);
	      }
	    }
	  }
	  greensfn[n++] = numerator*sum1/denominator;
	} else greensfn[n++] = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute the PPPM long-range force, energy, virial 
------------------------------------------------------------------------- */

void PPPM::compute(int eflag, int vflag)
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
    memory->destroy_2d_int_array(part2grid);
    nmax = atom->nmax;
    part2grid = memory->create_2d_int_array(nmax,3,"pppm:part2grid");
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
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order 
------------------------------------------------------------------------- */

void PPPM::allocate()
{
  density_brick = 
    memory->create_3d_double_array(nzlo_out,nzhi_out,nylo_out,nyhi_out,
				   nxlo_out,nxhi_out,"pppm:density_brick");
  vdx_brick =
    memory->create_3d_double_array(nzlo_out,nzhi_out,nylo_out,nyhi_out,
				   nxlo_out,nxhi_out,"pppm:vdx_brick");
  vdy_brick = 
    memory->create_3d_double_array(nzlo_out,nzhi_out,nylo_out,nyhi_out,
				   nxlo_out,nxhi_out,"pppm:vdy_brick");
  vdz_brick = 
    memory->create_3d_double_array(nzlo_out,nzhi_out,nylo_out,nyhi_out,
				   nxlo_out,nxhi_out,"pppm:vdz_brick");

  density_fft = 
    (double *) memory->smalloc(nfft_both*sizeof(double),"pppm:density_fft");
  greensfn = 
    (double *) memory->smalloc(nfft_both*sizeof(double),"pppm:greensfn");
  work1 = (double *) memory->smalloc(2*nfft_both*sizeof(double),"pppm:work1");
  work2 = (double *) memory->smalloc(2*nfft_both*sizeof(double),"pppm:work2");
  vg = memory->create_2d_double_array(nfft_both,6,"pppm:vg");

  fkx = memory->create_1d_double_array(nxlo_fft,nxhi_fft,"pppm:fkx");
  fky = memory->create_1d_double_array(nylo_fft,nyhi_fft,"pppm:fky");
  fkz = memory->create_1d_double_array(nzlo_fft,nzhi_fft,"pppm:fkz");

  buf1 = (double *) memory->smalloc(nbuf*sizeof(double),"pppm:buf1");
  buf2 = (double *) memory->smalloc(nbuf*sizeof(double),"pppm:buf2");

  // summation coeffs

  gf_b = new double[order];
  rho1d = memory->create_2d_double_array(3,-order/2,order/2,"pppm:rho1d");
  rho_coeff = memory->create_2d_double_array(order,(1-order)/2,order/2,
					     "pppm:rho_coeff");

  // create 2 FFTs and a Remap
  // 1st FFT keeps data in FFT decompostion
  // 2nd FFT returns data in 3d brick decomposition
  // remap takes data from 3d brick to FFT decomposition

  int tmp;

  fft1 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
		   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
		   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
		   0,0,&tmp);

  fft2 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
		   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
		   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
		   0,0,&tmp);

  remap = new Remap(lmp,world,
		    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
		    nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
		    1,0,0,2);
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order 
------------------------------------------------------------------------- */

void PPPM::deallocate()
{
  memory->destroy_3d_double_array(density_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy_3d_double_array(vdx_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy_3d_double_array(vdy_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy_3d_double_array(vdz_brick,nzlo_out,nylo_out,nxlo_out);

  memory->sfree(density_fft);
  memory->sfree(greensfn);
  memory->sfree(work1);
  memory->sfree(work2);
  memory->destroy_2d_double_array(vg);

  memory->destroy_1d_double_array(fkx,nxlo_fft);
  memory->destroy_1d_double_array(fky,nylo_fft);
  memory->destroy_1d_double_array(fkz,nzlo_fft);

  memory->sfree(buf1);
  memory->sfree(buf2);

  delete [] gf_b;
  memory->destroy_2d_double_array(rho1d,-order/2);
  memory->destroy_2d_double_array(rho_coeff,(1-order)/2);

  delete fft1;
  delete fft2;
  delete remap;
}

/* ----------------------------------------------------------------------
   set size of FFT grid (nx,ny,nz_pppm) and g_ewald 
------------------------------------------------------------------------- */

void PPPM::set_grid()
{
  // see JCP 109, pg 7698 for derivation of coefficients
  // higher order coefficients may be computed if needed

  double **acons = memory->create_2d_double_array(8,7,"pppm:acons");

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

  double q2 = qsqsum / force->dielectric;
  double natoms = atom->natoms;

  // use xprd,yprd,zprd even if triclinic so grid size is the same
  // adjust z dimension for 2d slab PPPM
  // 3d PPPM just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;
  
  // make initial g_ewald estimate
  // based on desired error and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  double hx,hy,hz;

  if (!gewaldflag)
    g_ewald = sqrt(-log(precision*sqrt(natoms*cutoff*xprd*yprd*zprd) / 
			(2.0*q2))) / cutoff;

  // set optimal nx_pppm,ny_pppm,nz_pppm based on order and precision
  // nz_pppm uses extended zprd_slab instead of zprd
  // h = 1/g_ewald is upper bound on h such that h*g_ewald <= 1
  // reduce it until precision target is met

  if (!gridflag) {
    double err;
    hx = hy = hz = 1/g_ewald;  

    nx_pppm = static_cast<int> (xprd/hx + 1);
    ny_pppm = static_cast<int> (yprd/hy + 1);
    nz_pppm = static_cast<int> (zprd_slab/hz + 1);

    err = rms(hx,xprd,natoms,q2,acons);
    while (err > precision) {
      err = rms(hx,xprd,natoms,q2,acons);
      nx_pppm++;
      hx = xprd/nx_pppm;
    }

    err = rms(hy,yprd,natoms,q2,acons);
    while (err > precision) {
      err = rms(hy,yprd,natoms,q2,acons);
      ny_pppm++;
      hy = yprd/ny_pppm;
    }

    err = rms(hz,zprd_slab,natoms,q2,acons);
    while (err > precision) {
      err = rms(hz,zprd_slab,natoms,q2,acons);
      nz_pppm++;
      hz = zprd_slab/nz_pppm;
    }
  }

  // boost grid size until it is factorable

  while (!factorable(nx_pppm)) nx_pppm++;
  while (!factorable(ny_pppm)) ny_pppm++;
  while (!factorable(nz_pppm)) nz_pppm++;

  // adjust g_ewald for new grid size

  hx = xprd/nx_pppm;
  hy = yprd/ny_pppm;
  hz = zprd_slab/nz_pppm;

  if (!gewaldflag) {
    double gew1,gew2,dgew,f,fmid,hmin,rtb;
    int ncount;

    gew1 = 0.0;
    g_ewald = gew1;
    f = diffpr(hx,hy,hz,q2,acons);

    hmin = MIN(hx,MIN(hy,hz));
    gew2 = 10/hmin;
    g_ewald = gew2;
    fmid = diffpr(hx,hy,hz,q2,acons);

    if (f*fmid >= 0.0) error->all("Cannot compute PPPM G");
    rtb = f < 0.0 ? (dgew=gew2-gew1,gew1) : (dgew=gew1-gew2,gew2);
    ncount = 0;
    while (fabs(dgew) > SMALL && fmid != 0.0) {
      dgew *= 0.5;
      g_ewald = rtb + dgew;
      fmid = diffpr(hx,hy,hz,q2,acons);      
      if (fmid <= 0.0) rtb = g_ewald;
      ncount++;
      if (ncount > LARGE) error->all("Cannot compute PPPM G");
    }
  }

  // final RMS precision

  double lprx = rms(hx,xprd,natoms,q2,acons);
  double lpry = rms(hy,yprd,natoms,q2,acons);
  double lprz = rms(hz,zprd_slab,natoms,q2,acons);
  double lpr = sqrt(lprx*lprx + lpry*lpry + lprz*lprz) / sqrt(3.0);
  double spr = 2.0*q2 * exp(-g_ewald*g_ewald*cutoff*cutoff) / 
    sqrt(natoms*cutoff*xprd*yprd*zprd_slab);

  // free local memory

  memory->destroy_2d_double_array(acons);

  // print info

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  G vector = %g\n",g_ewald);
      fprintf(screen,"  grid = %d %d %d\n",nx_pppm,ny_pppm,nz_pppm);
      fprintf(screen,"  stencil order = %d\n",order);
      fprintf(screen,"  RMS precision = %g\n",MAX(lpr,spr));
    }
    if (logfile) {
      fprintf(logfile,"  G vector = %g\n",g_ewald);
      fprintf(logfile,"  grid = %d %d %d\n",nx_pppm,ny_pppm,nz_pppm);
      fprintf(logfile,"  stencil order = %d\n",order);
      fprintf(logfile,"  RMS precision = %g\n",MAX(lpr,spr));
    }
  }
}

/* ----------------------------------------------------------------------
   check if all factors of n are in list of factors
   return 1 if yes, 0 if no 
------------------------------------------------------------------------- */

int PPPM::factorable(int n)
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
   compute RMS precision for a dimension
------------------------------------------------------------------------- */

double PPPM::rms(double h, double prd, double natoms,
		 double q2, double **acons)
{
  double sum = 0.0;
  for (int m = 0; m < order; m++) 
    sum += acons[order][m] * pow(h*g_ewald,2.0*m);
  double value = q2 * pow(h*g_ewald,order) *
    sqrt(g_ewald*prd*sqrt(2.0*PI)*sum/natoms) / (prd*prd);
  return value;
}

/* ----------------------------------------------------------------------
   compute difference in real-space and kspace RMS precision
------------------------------------------------------------------------- */

double PPPM::diffpr(double hx, double hy, double hz, double q2, double **acons)
{
  double lprx,lpry,lprz,kspace_prec,real_prec;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double natoms = atom->natoms;

  lprx = rms(hx,xprd,natoms,q2,acons);
  lpry = rms(hy,yprd,natoms,q2,acons);
  lprz = rms(hz,zprd*slab_volfactor,natoms,q2,acons);
  kspace_prec = sqrt(lprx*lprx + lpry*lpry + lprz*lprz) / sqrt(3.0);
  real_prec = 2.0*q2 * exp(-g_ewald*g_ewald*cutoff*cutoff) / 
   sqrt(natoms*cutoff*xprd*yprd*zprd);
  double value = kspace_prec - real_prec;
  return value;
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

double PPPM::gf_denom(double x, double y, double z)
{
  double sx,sy,sz;
  sz = sy = sx = 0.0;
  for (int l = order-1; l >= 0; l--) {
    sx = gf_b[l] + sx*x;
    sy = gf_b[l] + sy*y;
    sz = gf_b[l] + sz*z;
  }
  double s = sx*sy*sz;
  return s*s;
}

/* ----------------------------------------------------------------------
   pre-compute Green's function denominator expansion coeffs, Gamma(2n) 
------------------------------------------------------------------------- */

void PPPM::compute_gf_denom()
{
  int k,l,m;
  
  for (l = 1; l < order; l++) gf_b[l] = 0.0;
  gf_b[0] = 1.0;
  
  for (m = 1; m < order; m++) {
    for (l = m; l > 0; l--) 
      gf_b[l] = 4.0 * (gf_b[l]*(l-m)*(l-m-0.5)-gf_b[l-1]*(l-m-1)*(l-m-1));
    gf_b[0] = 4.0 * (gf_b[0]*(l-m)*(l-m-0.5));
  }

  int ifact = 1;
  for (k = 1; k < 2*order; k++) ifact *= k;
  double gaminv = 1.0/ifact;
  for (l = 0; l < order; l++) gf_b[l] *= gaminv;
}

/* ----------------------------------------------------------------------
   ghost-swap to accumulate full density in brick decomposition 
   remap density from 3d brick decomposition to FFT decomposition
------------------------------------------------------------------------- */

void PPPM::brick2fft()
{
  int i,n,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my ghosts for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in+1; ix <= nxhi_out; ix++)
	buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix < nxlo_in+nxlo_ghost; ix++)
	density_brick[iz][iy][ix] += buf2[n++];

  // pack my ghosts for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_out; ix < nxlo_in; ix++)
	buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in-nxhi_ghost+1; ix <= nxhi_in; ix++)
	density_brick[iz][iy][ix] += buf2[n++];

  // pack my ghosts for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in+1; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy < nylo_in+nylo_ghost; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	density_brick[iz][iy][ix] += buf2[n++];

  // pack my ghosts for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy < nylo_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in-nyhi_ghost+1; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	density_brick[iz][iy][ix] += buf2[n++];

  // pack my ghosts for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzhi_in+1; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_in; iz < nzlo_in+nzlo_ghost; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	density_brick[iz][iy][ix] += buf2[n++];

  // pack my ghosts for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz < nzlo_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzhi_in-nzhi_ghost+1; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	density_brick[iz][iy][ix] += buf2[n++];

  // remap from 3d brick decomposition to FFT decomposition
  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	density_fft[n++] = density_brick[iz][iy][ix];

  remap->perform(density_fft,density_fft,work1);
}

/* ----------------------------------------------------------------------
   ghost-swap to fill ghost cells of my brick with field values
------------------------------------------------------------------------- */

void PPPM::fillbrick()
{
  int i,n,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my real cells for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzhi_in-nzhi_ghost+1; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
	buf1[n++] = vdx_brick[iz][iy][ix];
	buf1[n++] = vdy_brick[iz][iy][ix];
	buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz < nzlo_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
	vdx_brick[iz][iy][ix] = buf2[n++];
	vdy_brick[iz][iy][ix] = buf2[n++];
	vdz_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_in; iz < nzlo_in+nzlo_ghost; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
	buf1[n++] = vdx_brick[iz][iy][ix];
	buf1[n++] = vdy_brick[iz][iy][ix];
	buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzhi_in+1; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
	vdx_brick[iz][iy][ix] = buf2[n++];
	vdy_brick[iz][iy][ix] = buf2[n++];
	vdz_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in-nyhi_ghost+1; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
	buf1[n++] = vdx_brick[iz][iy][ix];
	buf1[n++] = vdy_brick[iz][iy][ix];
	buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy < nylo_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
	vdx_brick[iz][iy][ix] = buf2[n++];
	vdy_brick[iz][iy][ix] = buf2[n++];
	vdz_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy < nylo_in+nylo_ghost; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
	buf1[n++] = vdx_brick[iz][iy][ix];
	buf1[n++] = vdy_brick[iz][iy][ix];
	buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in+1; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
	vdx_brick[iz][iy][ix] = buf2[n++];
	vdy_brick[iz][iy][ix] = buf2[n++];
	vdz_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in-nxhi_ghost+1; ix <= nxhi_in; ix++) {
	buf1[n++] = vdx_brick[iz][iy][ix];
	buf1[n++] = vdy_brick[iz][iy][ix];
	buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_out; ix < nxlo_in; ix++) {
	vdx_brick[iz][iy][ix] = buf2[n++];
	vdy_brick[iz][iy][ix] = buf2[n++];
	vdz_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix < nxlo_in+nxlo_ghost; ix++) {
	buf1[n++] = vdx_brick[iz][iy][ix];
	buf1[n++] = vdy_brick[iz][iy][ix];
	buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_DOUBLE,comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf1,n,MPI_DOUBLE,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in+1; ix <= nxhi_out; ix++) {
	vdx_brick[iz][iy][ix] = buf2[n++];
	vdy_brick[iz][iy][ix] = buf2[n++];
	vdz_brick[iz][iy][ix] = buf2[n++];
      }
}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array 
------------------------------------------------------------------------- */

void PPPM::particle_map()
{
  int nx,ny,nz;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    
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
	nz+nlower < nzlo_out || nz+nupper > nzhi_out) flag++;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all("Out of range atoms - cannot compute PPPM");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid 
------------------------------------------------------------------------- */

void PPPM::make_rho()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  double dx,dy,dz,x0,y0,z0;

  // clear 3d density array

  double *vec = &density_brick[nzlo_out][nylo_out][nxlo_out];
  for (i = 0; i < ngrid; i++) vec[i] = 0.0;

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
   FFT-based Poisson solver 
------------------------------------------------------------------------- */

void PPPM::poisson(int eflag, int vflag)
{
  int i,j,k,n;
  double eng;

  // transform charge density (r -> k) 

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] = density_fft[i];
    work1[n++] = 0.0;
  }
 
  fft1->compute(work1,work1,1);

  // if requested, compute energy and virial contribution

  double scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm);
  double s2 = scaleinv*scaleinv;

  if (eflag || vflag) {
    if (vflag) {
      n = 0;
      for (i = 0; i < nfft; i++) {
	eng = s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
	for (j = 0; j < 6; j++) virial[j] += eng*vg[i][j];
	energy += eng;
	n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nfft; i++) {
	energy += 
	  s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
	n += 2;
      }
    }
  }

  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] *= scaleinv * greensfn[i];
    work1[n++] *= scaleinv * greensfn[i];
  }

  // compute gradients of V(r) in each of 3 dims by transformimg -ik*V(k)
  // FFT leaves data in 3d brick decomposition
  // copy it into inner portion of vdx,vdy,vdz arrays

  // x direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
	work2[n] = fkx[i]*work1[n+1];
	work2[n+1] = -fkx[i]*work1[n];
	n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
	vdx_brick[k][j][i] = work2[n];
	n += 2;
      }

  // y direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
	work2[n] = fky[j]*work1[n+1];
	work2[n+1] = -fky[j]*work1[n];
	n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
	vdy_brick[k][j][i] = work2[n];
	n += 2;
      }

  // z direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
	work2[n] = fkz[k]*work1[n+1];
	work2[n+1] = -fkz[k]*work1[n];
	n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
	vdz_brick[k][j][i] = work2[n];
	n += 2;
      }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles 
------------------------------------------------------------------------- */

void PPPM::fieldforce()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  double dx,dy,dz,x0,y0,z0;
  double ek[3];

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

    compute_rho1d(dx,dy,dz);

    ek[0] = ek[1] = ek[2] = 0.0;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
	my = m+ny;
	y0 = z0*rho1d[1][m];
	for (l = nlower; l <= nupper; l++) {
	  mx = l+nx;
	  x0 = y0*rho1d[0][l];
	  ek[0] -= x0*vdx_brick[mz][my][mx];;
	  ek[1] -= x0*vdy_brick[mz][my][mx];;
	  ek[2] -= x0*vdz_brick[mz][my][mx];;
	}
      }
    }

    // convert E-field to force

    f[i][0] += qqrd2e*q[i]*ek[0];
    f[i][1] += qqrd2e*q[i]*ek[1];
    f[i][2] += qqrd2e*q[i]*ek[2];
  }
}

/* ----------------------------------------------------------------------
   map nprocs to NX by NY grid as PX by PY procs - return optimal px,py 
------------------------------------------------------------------------- */

void PPPM::procs2grid2d(int nprocs, int nx, int ny, int *px, int *py)
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

void PPPM::compute_rho1d(double dx, double dy, double dz)
{
  int k,l;

  for (k = (1-order)/2; k <= order/2; k++) {
    rho1d[0][k] = 0.0;
    rho1d[1][k] = 0.0;
    rho1d[2][k] = 0.0;
    for (l = order-1; l >= 0; l--) {
      rho1d[0][k] = rho_coeff[l][k] + rho1d[0][k]*dx;
      rho1d[1][k] = rho_coeff[l][k] + rho1d[1][k]*dy;
      rho1d[2][k] = rho_coeff[l][k] + rho1d[2][k]*dz;
    }
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

void PPPM::compute_rho_coeff()
{
  int j,k,l,m;
  double s;

  double **a = memory->create_2d_double_array(order,-order,order,"pppm:a");

  for (k = -order; k <= order; k++) 
    for (l = 0; l < order; l++)
      a[l][k] = 0.0;
        
  a[0][0] = 1.0;
  for (j = 1; j < order; j++) {
    for (k = -j; k <= j; k += 2) {
      s = 0.0;
      for (l = 0; l < j; l++) {
	a[l+1][k] = (a[l][k+1]-a[l][k-1]) / (l+1);
	s += pow(0.5,(double) l+1) * 
	  (a[l][k-1] + pow(-1.0,(double) l) * a[l][k+1]) / (l+1);
      }
      a[0][k] = s;
    }
  }

  m = (1-order)/2;
  for (k = -(order-1); k < order; k += 2) {
    for (l = 0; l < order; l++)
      rho_coeff[l][m] = a[l][k];
    m++;
  }

  memory->destroy_2d_double_array(a,-order);
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if 
   adequate empty space is left between repeating slabs (J. Chem. Phys. 
   111, 3155).  Slabs defined here to be parallel to the xy plane. 
------------------------------------------------------------------------- */

void PPPM::slabcorr(int eflag)
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i]*x[i][2];
  
  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // compute corrections
  
  double e_slabcorr = 2.0*PI*dipole_all*dipole_all/volume;
  
  if (eflag) energy += qqrd2e*e_slabcorr;

  // add on force corrections

  double ffact = -4.0*PI*dipole_all/volume; 
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) f[i][2] += qqrd2e*q[i]*ffact;
}

/* ----------------------------------------------------------------------
   perform and time the 4 FFTs required for N timesteps
------------------------------------------------------------------------- */

void PPPM::timing(int n, double &time3d, double &time1d)
{
  double time1,time2;

  for (int i = 0; i < 2*nfft_both; i++) work1[i] = 0.0;

  MPI_Barrier(world);
  time1 = MPI_Wtime();

  for (int i = 0; i < n; i++) {
    fft1->compute(work1,work1,1);
    fft2->compute(work1,work1,-1);
    fft2->compute(work1,work1,-1);
    fft2->compute(work1,work1,-1);
  }

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  time3d = time2 - time1;

  MPI_Barrier(world);
  time1 = MPI_Wtime();

  for (int i = 0; i < n; i++) {
    fft1->timing1d(work1,nfft_both,1);
    fft2->timing1d(work1,nfft_both,-1);
    fft2->timing1d(work1,nfft_both,-1);
    fft2->timing1d(work1,nfft_both,-1);
  }

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  time1d = time2 - time1;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays 
------------------------------------------------------------------------- */

double PPPM::memory_usage()
{
  double bytes = nmax*3 * sizeof(double);
  int nbrick = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) * 
    (nzhi_out-nzlo_out+1);
  bytes += 4 * nbrick * sizeof(double);
  bytes += 6 * nfft_both * sizeof(double);
  bytes += nfft_both*6 * sizeof(double);
  bytes += 2 * nbuf * sizeof(double);
  return bytes;
}
