// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "amoeba_convolution.h"

#include "comm.h"
#include "domain.h"
#include "fft3d_wrap.h"
#include "gridcomm.h"
#include "memory.h"
#include "neighbor.h"
#include "remap_wrap.h"
#include "update.h"

using namespace LAMMPS_NS;

// DEBUG

#define DEBUG_AMOEBA 0
#if DEBUG_AMOEBA
char *labels[7] =
  {(char *) "MPOLE_GRID", (char *) "POLAR_GRID",
   (char *) "POLAR_GRIDC", (char *) "DISP_GRID",
   (char *) "INDUCE_GRID", (char *) "INDUCE_GRIDC"};

enum{GRIDBRICK_OUT,GRIDBRICK_IN,FFT,CFFT1,CFFT2};
#endif
// END DEBUG

enum{MPOLE_GRID,POLAR_GRID,POLAR_GRIDC,DISP_GRID,INDUCE_GRID,INDUCE_GRIDC};

//#define SCALE 1
#define SCALE 0

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ----------------------------------------------------------------------
   partition an FFT grid across processors
   both for a brick and FFT x pencil decomposition
   nx,nz,nz = global FFT grid size
   order = size of stencil in each dimension that maps atoms to grid
   adapted from PPPM::set_grid_local()
------------------------------------------------------------------------- */

AmoebaConvolution::AmoebaConvolution(LAMMPS *lmp, Pair *pair,
                                     int nx_caller, int ny_caller, int nz_caller,
                                     int order_caller, int which_caller) :
  Pointers(lmp)
{
  amoeba = pair;
  nx = nx_caller;
  ny = ny_caller;
  nz = nz_caller;
  order = order_caller;
  which = which_caller;

  flag3d = 1;
  if (which == POLAR_GRIDC || which == INDUCE_GRIDC) flag3d = 0;

  nfft_global = (bigint) nx * ny * nz;

  // global indices of grid range from 0 to N-1
  // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
  //   global grid that I own without ghost cells
  // both non-tiled and tiled proc layouts use 0-1 fractional subdomain info

  if (comm->layout != Comm::LAYOUT_TILED) {
    nxlo_in = static_cast<int> (comm->xsplit[comm->myloc[0]] * nx);
    nxhi_in = static_cast<int> (comm->xsplit[comm->myloc[0]+1] * nx) - 1;
    nylo_in = static_cast<int> (comm->ysplit[comm->myloc[1]] * ny);
    nyhi_in = static_cast<int> (comm->ysplit[comm->myloc[1]+1] * ny) - 1;
    nzlo_in = static_cast<int> (comm->zsplit[comm->myloc[2]] * nz);
    nzhi_in = static_cast<int> (comm->zsplit[comm->myloc[2]+1] * nz) - 1;

  } else {
    nxlo_in = static_cast<int> (comm->mysplit[0][0] * nx);
    nxhi_in = static_cast<int> (comm->mysplit[0][1] * nx) - 1;
    nylo_in = static_cast<int> (comm->mysplit[1][0] * ny);
    nyhi_in = static_cast<int> (comm->mysplit[1][1] * ny) - 1;
    nzlo_in = static_cast<int> (comm->mysplit[2][0] * nz);
    nzhi_in = static_cast<int> (comm->mysplit[2][1] * nz) - 1;
  }

  // nlower,nupper = stencil size for mapping particles to FFT grid

  int nlower = -(order-1)/2;
  int nupper = order/2;

  // nlo_out,nhi_out = lower/upper limits of the 3d sub-brick of
  //   global grid that my particles can contribute charge to
  // effectively nlo_in,nhi_in + ghost cells
  // nlo,nhi = global coords of grid pt to "lower left" of smallest/largest
  //           position a particle in my box can be at
  // dist[3] = particle position bound = subbox + skin/2.0
  //   convert to triclinic if necessary
  // nlo_out,nhi_out = nlo,nhi + stencil size for particle mapping

  double *prd,*boxlo,*sublo,*subhi;
  int triclinic = domain->triclinic;

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

  double dist[3] = {0.0,0.0,0.0};
  double cuthalf = 0.5*neighbor->skin;
  if (triclinic == 0) dist[0] = dist[1] = dist[2] = cuthalf;
  else kspacebbox(cuthalf,&dist[0]);

  int nlo,nhi;

  nlo = static_cast<int> ((sublo[0]-dist[0]-boxlo[0]) * nx/xprd);
  nhi = static_cast<int> ((subhi[0]+dist[0]-boxlo[0]) * nx/xprd);
  nxlo_out = nlo + nlower;
  nxhi_out = nhi + nupper;

  nlo = static_cast<int> ((sublo[1]-dist[1]-boxlo[1]) * ny/yprd);
  nhi = static_cast<int> ((subhi[1]+dist[1]-boxlo[1]) * ny/yprd);
  nylo_out = nlo + nlower;
  nyhi_out = nhi + nupper;

  nlo = static_cast<int> ((sublo[2]-dist[2]-boxlo[2]) * nz/zprd);
  nhi = static_cast<int> ((subhi[2]+dist[2]-boxlo[2]) * nz/zprd);
  nzlo_out = nlo + nlower;
  nzhi_out = nhi + nupper;

  // x-pencil decomposition of FFT mesh
  // global indices range from 0 to N-1
  // each proc owns entire x-dimension, clumps of columns in y,z dimensions
  // npey_fft,npez_fft = # of procs in y,z dims
  // if nprocs is small enough, proc can own 1 or more entire xy planes,
  //   else proc owns 2d sub-blocks of yz plane
  // me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
  // nlo_fft,nhi_fft = lower/upper limit of the section
  //   of the global FFT mesh that I own in x-pencil decomposition

  int me = comm->me;
  int nprocs = comm->nprocs;

  int npey_fft,npez_fft;
  if (nz >= nprocs) {
    npey_fft = 1;
    npez_fft = nprocs;
  } else procs2grid2d(nprocs,ny,nz,npey_fft,npez_fft);

  int me_y = me % npey_fft;
  int me_z = me / npey_fft;

  nxlo_fft = 0;
  nxhi_fft = nx - 1;
  nylo_fft = me_y*ny/npey_fft;
  nyhi_fft = (me_y+1)*ny/npey_fft - 1;
  nzlo_fft = me_z*nz/npez_fft;
  nzhi_fft = (me_z+1)*nz/npez_fft - 1;

  // grid sizes
  // nbrick_owned = owned grid points in brick decomp
  // nbrick_ghosts = owned + ghost grid points in grid decomp
  // nfft_owned = owned grid points in FFT decomp
  // ngrid_either = max of nbrick_onwed and nfft_owned
  // nfft = total FFT grid points

  nbrick_owned = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1) *
    (nzhi_in-nzlo_in+1);
  nbrick_ghosts = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
    (nzhi_out-nzlo_out+1);
  nfft_owned = (nxhi_fft-nxlo_fft+1) * (nyhi_fft-nylo_fft+1) *
    (nzhi_fft-nzlo_fft+1);

  ngrid_either = MAX(nbrick_owned,nfft_owned);

  // instantiate FFT, GridComm, and Remap

  int tmp;

  fft1 = new FFT3d(lmp,world,nx,ny,nz,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
       1,0,&tmp,0);
  //       0,0,&tmp,0);

  fft2 = new FFT3d(lmp,world,nx,ny,nz,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                   //1,0,&tmp,0);
                   0,0,&tmp,0);

  gc = new GridComm(lmp,world,nx,ny,nz,
                    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                    nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out);

  int nqty = flag3d ? 1 : 2;
  remap = new Remap(lmp,world,
                    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                    nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                    nqty,0,0,FFT_PRECISION,0);

  // memory allocations

  if (flag3d) {
    memory->create3d_offset(grid_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"amoeba:grid_brick");
    grid_brick_start = &grid_brick[nzlo_out][nylo_out][nxlo_out];
    cgrid_brick = nullptr;
  } else {
    memory->create4d_offset_last(cgrid_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                                 nxlo_out,nxhi_out,2,"amoeba:cgrid_brick");
    grid_brick_start = &cgrid_brick[nzlo_out][nylo_out][nxlo_out][0];
    grid_brick = nullptr;
  }

  memory->create(grid_fft,ngrid_either,"amoeba:grid_fft");
  memory->create(cfft,2*ngrid_either,"amoeba:cfft");

  int ngc_buf1,ngc_buf2;
  gc->setup(ngc_buf1,ngc_buf2);
  memory->create(gc_buf1,nqty*ngc_buf1,"amoeba:gc_buf1");
  memory->create(gc_buf2,nqty*ngc_buf2,"amoeba:gc_buf2");

  memory->create(remap_buf,nqty*nfft_owned,"amoeba:remap_buf");
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

AmoebaConvolution::~AmoebaConvolution()
{
  memory->destroy3d_offset(grid_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy4d_offset_last(cgrid_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy(grid_fft);
  memory->destroy(cfft);
  memory->destroy(gc_buf1);
  memory->destroy(gc_buf2);
  memory->destroy(remap_buf);

  delete fft1;
  delete fft2;
  delete gc;
  delete remap;
}

/* ----------------------------------------------------------------------
   zero brick grid, including ghosts
   can be 3d real or 4d complex array
   return pointer to data in brick grid, caller casts to 3d or 4d
------------------------------------------------------------------------- */

void *AmoebaConvolution::zero()
{
  if (flag3d) return zero_3d();
  return zero_4d();
}

/* ---------------------------------------------------------------------- */

void *AmoebaConvolution::zero_3d()
{
  if (!grid_brick) return nullptr;
  memset(&(grid_brick[nzlo_out][nylo_out][nxlo_out]),0,
         nbrick_ghosts*sizeof(FFT_SCALAR));
  return (void *) grid_brick;
}

/* ---------------------------------------------------------------------- */

void *AmoebaConvolution::zero_4d()
{
  if (!cgrid_brick) return nullptr;
  memset(&(cgrid_brick[nzlo_out][nylo_out][nxlo_out][0]),0,
         2*nbrick_ghosts*sizeof(FFT_SCALAR));
  return (void *) cgrid_brick;
}

/* ----------------------------------------------------------------------
   perform pre-convolution grid operations
   can be 3d real or 4d complex array
   return pointer to complex cfft vector
------------------------------------------------------------------------- */

FFT_SCALAR *AmoebaConvolution::pre_convolution()
{
  if (flag3d) return pre_convolution_3d();
  return pre_convolution_4d();
}

/* ----------------------------------------------------------------------
   perform pre-convolution grid operations for 3d grid_brick array
------------------------------------------------------------------------- */

FFT_SCALAR *AmoebaConvolution::pre_convolution_3d()
{
  int ix,iy,iz,n;

  // reverse comm for 3d brick grid + ghosts

#if DEBUG_AMOEBA
  debug_scalar(GRIDBRICK_OUT,"PRE Convo / PRE GridComm");
#endif

  gc->reverse_comm(GridComm::PAIR,amoeba,1,sizeof(FFT_SCALAR),which,
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);

#if DEBUG_AMOEBA
  debug_scalar(GRIDBRICK_IN,"PRE Convo / POST GridComm");
  debug_file(GRIDBRICK_IN,"pre.convo.post.gridcomm");
#endif

  // copy owned 3d brick grid values to FFT grid

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        grid_fft[n++] = grid_brick[iz][iy][ix];

  // remap FFT grid from brick to x pencil partitioning

  remap->perform(grid_fft,grid_fft,remap_buf);

#if DEBUG_AMOEBA
  debug_scalar(FFT,"PRE Convo / POST Remap");
  debug_file(FFT,"pre.convo.post.remap");
#endif

  // copy real values into complex grid

  n = 0;
  for (int i = 0; i < nfft_owned; i++) {
    cfft[n++] = grid_fft[i];
    cfft[n++] = ZEROF;
  }

  // perform forward FFT

  fft1->compute(cfft,cfft,FFT3d::FORWARD);

  if (SCALE) {
    double scale = 1.0/nfft_global;
    for (int i = 0; i < 2*nfft_owned; i++) cfft[i] *= scale;
  }

#if DEBUG_AMOEBA
  debug_scalar(CFFT1,"PRE Convo / POST FFT");
  debug_file(CFFT1,"pre.convo.post.fft");
#endif
  return cfft;
}

/* ----------------------------------------------------------------------
   perform pre-convolution grid operations for 4d cgrid_brick array
------------------------------------------------------------------------- */

FFT_SCALAR *AmoebaConvolution::pre_convolution_4d()
{
  int ix,iy,iz,n;

  // reverse comm for 4d brick grid + ghosts

#if DEBUG_AMOEBA
  debug_scalar(GRIDBRICK_OUT,"PRE Convo / PRE GridComm");
#endif

  gc->reverse_comm(GridComm::PAIR,amoeba,2,sizeof(FFT_SCALAR),which,
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);

#if DEBUG_AMOEBA
  debug_scalar(GRIDBRICK_IN,"PRE Convo / POST GridComm");
  debug_file(GRIDBRICK_IN,"pre.convo.post.gridcomm");
#endif
  // copy owned 4d brick grid values to FFT grid

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        cfft[n++] = cgrid_brick[iz][iy][ix][0];
        cfft[n++] = cgrid_brick[iz][iy][ix][1];
      }

  // remap FFT grid from brick to x pencil partitioning
  // NOTE: could just setup FFT to start from brick decomp and skip remap

  remap->perform(cfft,cfft,remap_buf);

#if DEBUG_AMOEBA
  debug_scalar(FFT,"PRE Convo / POST Remap");
  debug_file(FFT,"pre.convo.post.remap");
#endif
  // perform forward FFT

  fft1->compute(cfft,cfft,FFT3d::FORWARD);

  if (SCALE) {
    double scale = 1.0/nfft_global;
    for (int i = 0; i < 2*nfft_owned; i++) cfft[i] *= scale;
  }

#if DEBUG_AMOEBA
  debug_scalar(CFFT1,"PRE Convo / POST FFT");
  debug_file(CFFT1,"pre.convo.post.fft");
#endif
  return cfft;
}

/* ----------------------------------------------------------------------
   perform post-convolution grid operations
   can be 3d real or 4d complex array
   return pointer to data in brick grid, caller casts to 3d or 4d
------------------------------------------------------------------------- */

void *AmoebaConvolution::post_convolution()
{
  if (flag3d) return post_convolution_3d();
  return post_convolution_4d();
}

/* ----------------------------------------------------------------------
   perform post-convolution grid operations for 3d grid_brick array
------------------------------------------------------------------------- */

void *AmoebaConvolution::post_convolution_3d()
{
  int ix,iy,iz,n;

  // perform backward FFT
#if DEBUG_AMOEBA
  debug_scalar(CFFT1,"POST Convo / PRE FFT");
  debug_file(CFFT1,"post.convo.pre.fft");
#endif
  fft2->compute(cfft,cfft,FFT3d::BACKWARD);

#if DEBUG_AMOEBA
  debug_scalar(CFFT2,"POST Convo / POST FFT");
  debug_file(CFFT2,"post.convo.post.fft");
#endif
  // copy real portion of 1d complex values into 3d real grid

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        grid_brick[iz][iy][ix] = cfft[n];
        n += 2;
      }

  // forward comm to populate ghost grid values

#if DEBUG_AMOEBA
  debug_scalar(GRIDBRICK_IN,"POST Convo / PRE gridcomm");
  debug_file(GRIDBRICK_IN,"post.convo.pre.gridcomm");
#endif
  gc->forward_comm(GridComm::PAIR,amoeba,1,sizeof(FFT_SCALAR),which,
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);

  return (void *) grid_brick;
}

/* ----------------------------------------------------------------------
   perform post-convolution grid operations for 4d cgrid_brick array
------------------------------------------------------------------------- */

void *AmoebaConvolution::post_convolution_4d()
{
  int ix,iy,iz,n;

  // perform backward FFT

#if DEBUG_AMOEBA
  debug_scalar(CFFT1,"POST Convo / PRE FFT");
  debug_file(CFFT1,"post.convo.pre.fft");
#endif
  fft2->compute(cfft,cfft,FFT3d::BACKWARD);

#if DEBUG_AMOEBA
  debug_scalar(CFFT2,"POST Convo / POST FFT");
  debug_file(CFFT2,"post.convo.post.fft");
#endif
  // copy 1d complex values into 4d complex grid

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        cgrid_brick[iz][iy][ix][0] = cfft[n++];
        cgrid_brick[iz][iy][ix][1] = cfft[n++];
      }

  // forward comm to populate ghost grid values

#if DEBUG_AMOEBA
  debug_scalar(GRIDBRICK_IN,"POST Convo / PRE gridcomm");
  debug_file(GRIDBRICK_IN,"post.convo.pre.gridcomm");
#endif
  gc->forward_comm(GridComm::PAIR,amoeba,2,sizeof(FFT_SCALAR),which,
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);

  return (void *) cgrid_brick;
}

/* ----------------------------------------------------------------------
   convert a sphere in box coords to an ellipsoid in lamda (0-1)
   coords and return the tight (axis-aligned) bounding box, does not
   preserve vector magnitude
   see http://www.loria.fr/~shornus/ellipsoid-bbox.html and
   http://yiningkarlli.blogspot.com/2013/02/
     bounding-boxes-for-ellipsoidsfigure.html
------------------------------------------------------------------------- */

void AmoebaConvolution::kspacebbox(double r, double *b)
{
  double *h = domain->h;
  double lx,ly,lz,xy,xz,yz;

  lx = h[0]; ly = h[1]; lz = h[2];
  yz = h[3]; xz = h[4]; xy = h[5];

  b[0] = r*sqrt(ly*ly*lz*lz + ly*ly*xz*xz - 2.0*ly*xy*xz*yz + lz*lz*xy*xy +
         xy*xy*yz*yz)/(lx*ly*lz);
  b[1] = r*sqrt(lz*lz + yz*yz)/(ly*lz);
  b[2] = r/lz;
}

/* ----------------------------------------------------------------------
   map nprocs to NX by NY grid as PX by PY procs - return optimal px,py
   copy of PPPM::procs2grid2d()
------------------------------------------------------------------------- */

void AmoebaConvolution::procs2grid2d(int nprocs, int nx, int ny, int &px, int &py)
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
        px = ipx;
        py = ipy;
      }
    }
    ipx++;
  }
}

#if DEBUG_AMOEBA
/* ----------------------------------------------------------------------
   output a scalar value to screen
   array = which array is being summed over
---------------------------------------------------------------------- */

void AmoebaConvolution::debug_scalar(int array, const char *label)
{
  double sum = 0.0;

  if (array == GRIDBRICK_OUT) {
    if (flag3d) {
      for (int iz = nzlo_out; iz <= nzhi_out; iz++)
        for (int iy = nylo_out; iy <= nyhi_out; iy++)
          for (int ix = nxlo_out; ix <= nxhi_out; ix++)
            sum += grid_brick[iz][iy][ix];
    } else {
      for (int iz = nzlo_out; iz <= nzhi_out; iz++)
        for (int iy = nylo_out; iy <= nyhi_out; iy++)
          for (int ix = nxlo_out; ix <= nxhi_out; ix++) {
            sum += cgrid_brick[iz][iy][ix][0];
            sum += cgrid_brick[iz][iy][ix][1];
          }
    }
  }

  if (array == GRIDBRICK_IN) {
    if (flag3d) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            sum += grid_brick[iz][iy][ix];
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
            sum += cgrid_brick[iz][iy][ix][0];
            sum += cgrid_brick[iz][iy][ix][1];
          }
    }
  }

  if (array == FFT) {
    if (flag3d) {
      for (int i = 0; i < nfft_owned; i++)
        sum += grid_fft[i];
    } else {
      for (int i = 0; i < 2*nfft_owned; i++)
        sum += cfft[i];
    }
  }

  if (array == CFFT1) {
    for (int i = 0; i < 2*nfft_owned; i++)
      sum += cfft[i];
  }

  if (array == CFFT2) {
    for (int i = 0; i < 2*nbrick_owned; i++)
      sum += cfft[i];
  }

  /*
  double sumall;
  MPI_Allreduce(&sum,&sumall,1,MPI_DOUBLE,MPI_SUM,world);
  if (comm->me == 0) printf("%s: %s: %12.8g\n",labels[which],label,sumall);
  */
}

/* ----------------------------------------------------------------------
   dump grid values to a file
   array = which array is being output
---------------------------------------------------------------------- */

void AmoebaConvolution::debug_file(int array, const char *label)
{
  FILE *fp;

  int me = comm->me;
  int nprocs = comm->nprocs;

  // open file

  char fname[128];
  sprintf(fname,"tmp.%s.%s",labels[which],label);
  if (me == 0) fp = fopen(fname,"w");

  // file header
  // ncol = # of columns, including grid cell ID

  bigint ntot = nx * ny * nz;

  int ncol;
  char *columns;

  if (array == CFFT1 || array == CFFT2 || !flag3d) {
    ncol = 3;
    columns = (char *) "id real imag";
  } else {
    ncol = 2;
    columns = (char *) "id value";
  }

  char boundstr[9];          // encoding of boundary flags
  domain->boundary_string(boundstr);

  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,BIGINT_FORMAT "\n",ntot);
    fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[0],domain->boxhi[0]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[1],domain->boxhi[1]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[2],domain->boxhi[2]);
    fprintf(fp,"ITEM: ATOMS %s\n",columns);
  }

  // pack my values
  // ngrid = # of grid cells I own

  int ngrid;
  if (array == GRIDBRICK_IN) ngrid = nbrick_owned;
  else if (array == FFT) ngrid = nfft_owned;
  else if (array == CFFT1) ngrid = nfft_owned;
  else if (array == CFFT2) ngrid = nbrick_owned;

  int ngridmax;
  MPI_Allreduce(&ngrid,&ngridmax,1,MPI_INT,MPI_MAX,world);

  double *buf,*buf2;
  memory->create(buf,ncol*ngridmax,"amoeba:buf");
  memory->create(buf2,ncol*ngridmax,"amoeba:buf2");

  ngrid = 0;

  if (array == GRIDBRICK_IN) {
    if (flag3d) {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
            int id = iz*ny*nx + iy*nx + ix + 1;
            buf[ncol*ngrid] = id;
            buf[ncol*ngrid+1] = grid_brick[iz][iy][ix];
            ngrid++;
          }
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
            int id = iz*ny*nx + iy*nx + ix + 1;
            buf[ncol*ngrid] = id;
            buf[ncol*ngrid+1] = cgrid_brick[iz][iy][ix][0];
            buf[ncol*ngrid+2] = cgrid_brick[iz][iy][ix][1];
            ngrid++;
          }
    }
  }

  if (array == FFT) {
    if (flag3d) {
      int m = 0;
      for (int iz = nzlo_fft; iz <= nzhi_fft; iz++)
        for (int iy = nylo_fft; iy <= nyhi_fft; iy++)
          for (int ix = nxlo_fft; ix <= nxhi_fft; ix++) {
            int id = iz*ny*nx + iy*nx + ix + 1;
            buf[ncol*ngrid] = id;
            buf[ncol*ngrid+1] = grid_fft[m++];
            ngrid++;
          }
    } else {
      int m = 0;
      for (int iz = nzlo_fft; iz <= nzhi_fft; iz++)
        for (int iy = nylo_fft; iy <= nyhi_fft; iy++)
          for (int ix = nxlo_fft; ix <= nxhi_fft; ix++) {
            int id = iz*ny*nx + iy*nx + ix + 1;
            buf[ncol*ngrid] = id;
            buf[ncol*ngrid+1] = cfft[m++];
            buf[ncol*ngrid+2] = cfft[m++];
            ngrid++;
          }
    }
  }

  if (array == CFFT1) {
    int m = 0;
    for (int iz = nzlo_fft; iz <= nzhi_fft; iz++)
      for (int iy = nylo_fft; iy <= nyhi_fft; iy++)
        for (int ix = nxlo_fft; ix <= nxhi_fft; ix++) {
          int id = iz*ny*nx + iy*nx + ix + 1;
          buf[ncol*ngrid] = id;
          buf[ncol*ngrid+1] = cfft[m++];
          buf[ncol*ngrid+2] = cfft[m++];
          ngrid++;
        }
  }

  if (array == CFFT2) {
    int m = 0;
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
          int id = iz*ny*nx + iy*nx + ix + 1;
          buf[ncol*ngrid] = id;
          buf[ncol*ngrid+1] = cfft[m++];
          buf[ncol*ngrid+2] = cfft[m++];
          ngrid++;
        }
  }

  // proc 0 outputs values
  // pings other procs, send/recv of their values

  int tmp,nlines;
  MPI_Request request;
  MPI_Status status;

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,ngridmax*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&nlines);
        nlines /= ncol;
      } else nlines = ngrid;

      int n = 0;
      for (int m = 0; m < nlines; m++) {
        if (ncol == 2)
          fprintf(fp,"%d %12.8g\n",(int) buf[n],buf[n+1]);
        else if (ncol == 3)
          fprintf(fp,"%d %12.8g %12.8g\n",(int ) buf[n],buf[n+1],buf[n+2]);
        n += ncol;
      }
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(buf,ngrid*ncol,MPI_DOUBLE,0,0,world);
  }

  // close file

  if (me == 0) fclose(fp);

  // clean up

  memory->destroy(buf);
  memory->destroy(buf2);
}
#endif
