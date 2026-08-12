// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   FFT based version of compute xrd.

   compute xrd evaluates F(K) = sum_j f_j(|K|/2) exp(2 pi i K.r_j) directly for
   every reciprocal lattice node, which costs O(nRows*N).  The same sum is a
   non-uniform Fourier transform, so it can instead be evaluated by spreading
   the atoms onto a uniform mesh with a Kaiser-Bessel window, taking one FFT,
   and dividing out the analytic Fourier transform of that window.

   compute xrd samples K = (i*dK[0], j*dK[1], k*dK[2]).  exp(2 pi i m x dK) has
   period 1/dK in x, so folding the atom coordinates into a cell of edge length
   1/dK[d] (here called the diffraction cell) is exact rather than an
   approximation, whether that cell is larger than the simulation box (c < 1) or
   smaller than it (c > 1).  A single uniform mesh spanning the diffraction cell
   therefore reproduces the requested nodes exactly.

   The atomic scattering factors depend on the atom type and on |K|, so one
   transform is done per distinct element and the results are combined with the
   scattering factor of that element at each node.

   Contributing author: derived from compute xrd by Shawn P. Coleman
------------------------------------------------------------------------- */

#include "compute_xrd_fft.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "update.h"

#include "compute_xrd_consts.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_PI;

static constexpr double MIN_OVERSAMPLE = 1.25;

/* ---------------------------------------------------------------------- */

ComputeXRDFFT::ComputeXRDFFT(LAMMPS *lmp, int narg, char **arg) :
  ComputeXRD(lmp, narg, arg),
  density_all(nullptr), density_slab(nullptr), work1(nullptr), recvcounts(nullptr),
  slot_of_type(nullptr), ztype_of_slot(nullptr), own_row(nullptr), own_idx(nullptr),
  own_deconv(nullptr), own_lp(nullptr), own_asf(nullptr), Fre(nullptr), Fim(nullptr),
  Iloc(nullptr), Iall(nullptr), fft1(nullptr)
{
  nprocs = comm->nprocs;
  setup_done = 0;
  fft_comm = MPI_COMM_NULL;

  order = nufft_order;
  oversample = nufft_oversample;

  if ((order < 3) || (order % 2 == 0))
    error->all(FLERR,"Compute XRD/FFT: order must be an odd number of 3 or larger");
  if (oversample < MIN_OVERSAMPLE)
    error->all(FLERR,"Compute XRD/FFT: oversample must be {} or larger",MIN_OVERSAMPLE);

  nlower = -(order-1)/2;
  nupper = (order-1)/2;

  nown = 0;
  nslot = 0;
}

/* ---------------------------------------------------------------------- */

ComputeXRDFFT::~ComputeXRDFFT()
{
  deallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeXRDFFT::deallocate()
{
  delete fft1;
  fft1 = nullptr;

  if (fft_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&fft_comm);
    fft_comm = MPI_COMM_NULL;
  }

  memory->destroy(density_all);
  memory->destroy(density_slab);
  memory->destroy(work1);
  memory->destroy(recvcounts);
  memory->destroy(slot_of_type);
  memory->destroy(ztype_of_slot);
  memory->destroy(own_row);
  memory->destroy(own_idx);
  memory->destroy(own_deconv);
  memory->destroy(own_lp);
  memory->destroy(own_asf);
  memory->destroy(Fre);
  memory->destroy(Fim);
  memory->destroy(Iloc);
  memory->destroy(Iall);

  density_all = density_slab = work1 = nullptr;
  recvcounts = slot_of_type = ztype_of_slot = own_row = own_idx = nullptr;
  own_deconv = own_lp = own_asf = Fre = Fim = Iloc = Iall = nullptr;

  setup_done = 0;
}

/* ----------------------------------------------------------------------
   enumerate the reciprocal lattice nodes exactly as compute xrd does, then
   build the FFT mesh that covers them
------------------------------------------------------------------------- */

void ComputeXRDFFT::init()
{
  // fills store_tmp with the (k,j,i) index of every row and column 0 of the
  // output array with the corresponding 2*theta.  Reusing it verbatim is what
  // guarantees that this style produces the same rows in the same order.

  ComputeXRD::init();

  deallocate();
  set_grid();
  setup_mesh();

#if defined(FFT_SINGLE)
  if (me == 0)
    error->warning(FLERR,"Compute XRD/FFT with single precision FFTs limits the accuracy "
                   "of weak diffuse intensities");
#endif

  if ((me == 0) && echo) {
    double mb = 1.0/1024.0/1024.0;
    double mem = (double)nmesh[0]*nmesh[1]*nmesh[2]*sizeof(FFT_SCALAR)*mb;
    utils::logmesg(lmp,"-----\nCompute XRD/FFT id:{}, # of relp:{}\n"
                   "Kaiser-Bessel order {}, FFT mesh {}x{}x{}, {} element transform(s)\n"
                   "Mesh memory per MPI rank = {:.4} Mbytes\n-----\n",
                   id,size_array_rows,order,nmesh[0],nmesh[1],nmesh[2],nslot,mem);
  }
}

/* ----------------------------------------------------------------------
   size the FFT mesh and pick the Kaiser-Bessel shape parameter
------------------------------------------------------------------------- */

void ComputeXRDFFT::set_grid()
{
  for (int d = 0; d < 3; d++) {

    // the mesh must resolve modes out to Knmax[d] with the requested
    // oversampling.  the 2*Knmax+1 floor keeps modes +Knmax and -Knmax from
    // aliasing onto the same mesh point.

    int nmin = (int) ceil(2.0*oversample*Knmax[d]);
    if (nmin < 2*Knmax[d]+1) nmin = 2*Knmax[d]+1;
    if (nmin < order+1) nmin = order+1;

    int n = nmin;
    while (!factorable(n)) n++;
    nmesh[d] = n;
  }

  bigint ntotal = (bigint)nmesh[0]*(bigint)nmesh[1]*(bigint)nmesh[2];
  if (ntotal > MAXSMALLINT)
    error->all(FLERR,"Compute XRD/FFT: FFT mesh of {}x{}x{} is too large; reduce the 2Theta "
               "range, increase the c values, or reduce oversample",
               nmesh[0],nmesh[1],nmesh[2]);

  // Beatty et al, IEEE Trans Med Imaging 24, 799 (2005), optimal Kaiser-Bessel
  // shape parameter for a stencil of width order at oversampling sigma

  for (int d = 0; d < 3; d++) {
    double sigma = (0.5*nmesh[d])/Knmax[d];
    double arg = (order*order/(sigma*sigma))*(sigma-0.5)*(sigma-0.5) - 0.8;
    if (arg < 0.0) arg = 0.0;
    kb_beta[d] = MY_PI*sqrt(arg);
  }

  // one mesh unit in each dimension corresponds to 1/(nmesh*dK) in distance
  // units, so a coordinate maps to mesh units by multiplying with dK*nmesh

  for (int d = 0; d < 3; d++) mesh_scale[d] = dK[d]*nmesh[d];
}

/* ----------------------------------------------------------------------
   allocate the mesh, set up the FFT, and cache everything that only depends
   on the mesh geometry
------------------------------------------------------------------------- */

void ComputeXRDFFT::setup_mesh()
{
  int nx = nmesh[0];
  int ny = nmesh[1];
  int nz = nmesh[2];
  int ntotal = nx*ny*nz;

  // the mesh is split into z slabs.  with more ranks than mesh planes the
  // trailing ranks own nothing: they still spread their own atoms, but they
  // take no part in the transform.  keeping the mesh independent of the number
  // of ranks is what makes the result reproducible across rank counts.

  nzlo_fft = me*nz/nprocs;
  nzhi_fft = (me+1)*nz/nprocs - 1;
  nslab = nzhi_fft - nzlo_fft + 1;
  nfft = nx*ny*nslab;

  memory->create(density_all,ntotal,"xrd/fft:density_all");
  memory->create(density_slab,MAX(nfft,1),"xrd/fft:density_slab");
  memory->create(work1,MAX(2*nfft,1),"xrd/fft:work1");
  memory->create(recvcounts,nprocs,"xrd/fft:recvcounts");

  for (int p = 0; p < nprocs; p++) {
    int zlo = p*nz/nprocs;
    int zhi = (p+1)*nz/nprocs - 1;
    recvcounts[p] = nx*ny*(zhi-zlo+1);
  }

  // input and output layouts are identical, which lets FFT3d skip the initial
  // remap: each rank already holds complete x and y pencils.  the transform
  // runs on a communicator holding only the ranks that own a slab, since not
  // every FFT backend accepts an empty block.

  MPI_Comm_split(world,(nslab > 0) ? 0 : MPI_UNDEFINED,me,&fft_comm);

  if (nslab > 0) {
    int tmp;
    fft1 = new FFT3d(lmp,fft_comm,nx,ny,nz,
                     0,nx-1,0,ny-1,nzlo_fft,nzhi_fft,
                     0,nx-1,0,ny-1,nzlo_fft,nzhi_fft,
                     0,0,&tmp,0,0);
  }

  // group the atom types by element: types sharing a row of ASFXRD have the
  // same scattering factor and can share a transform

  memory->create(slot_of_type,ntypes,"xrd/fft:slot_of_type");
  memory->create(ztype_of_slot,ntypes,"xrd/fft:ztype_of_slot");

  nslot = 0;
  for (int i = 0; i < ntypes; i++) {
    int slot = -1;
    for (int s = 0; s < nslot; s++)
      if (ztype_of_slot[s] == ztype[i]) { slot = s; break; }
    if (slot < 0) {
      slot = nslot;
      ztype_of_slot[nslot++] = ztype[i];
    }
    slot_of_type[i] = slot;
  }

  // cache the rows whose FFT mode this rank owns

  nown = 0;
  for (int n = 0; n < size_array_rows; n++) {
    int k = store_tmp[3*n];
    int iz = ((k % nz) + nz) % nz;
    if ((iz >= nzlo_fft) && (iz <= nzhi_fft)) nown++;
  }

  memory->create(own_row,MAX(nown,1),"xrd/fft:own_row");
  memory->create(own_idx,MAX(nown,1),"xrd/fft:own_idx");
  memory->create(own_deconv,MAX(nown,1),"xrd/fft:own_deconv");
  memory->create(own_lp,MAX(nown,1),"xrd/fft:own_lp");
  memory->create(own_asf,MAX(nslot*nown,1),"xrd/fft:own_asf");
  memory->create(Fre,MAX(nown,1),"xrd/fft:Fre");
  memory->create(Fim,MAX(nown,1),"xrd/fft:Fim");
  memory->create(Iloc,size_array_rows,"xrd/fft:Iloc");
  memory->create(Iall,size_array_rows,"xrd/fft:Iall");

  int a = 0;
  for (int n = 0; n < size_array_rows; n++) {
    int k = store_tmp[3*n];
    int j = store_tmp[3*n+1];
    int i = store_tmp[3*n+2];

    int iz = ((k % nz) + nz) % nz;
    if ((iz < nzlo_fft) || (iz > nzhi_fft)) continue;
    int iy = ((j % ny) + ny) % ny;
    int ix = ((i % nx) + nx) % nx;

    own_row[a] = n;
    own_idx[a] = ((iz-nzlo_fft)*ny + iy)*nx + ix;

    // the window transform must be evaluated at the signed mode index: it is
    // even in m but not periodic in m, so the folded index would be wrong

    own_deconv[a] = 1.0/(kb_window(i,nx,kb_beta[0]) * kb_window(j,ny,kb_beta[1]) *
                         kb_window(k,nz,kb_beta[2]));

    double K0 = i*dK[0];
    double K1 = j*dK[1];
    double K2 = k*dK[2];
    double dinv2 = K0*K0 + K1*K1 + K2*K2;
    double SinTheta_lambda = 0.5*sqrt(dinv2);

    if (LP == 1) {
      double SinTheta = SinTheta_lambda*lambda;
      double ang = asin(SinTheta);
      double Cos2Theta = cos(2.0*ang);
      double CosTheta = cos(ang);
      own_lp[a] = (1 + Cos2Theta*Cos2Theta) / (CosTheta * SinTheta * SinTheta);
    } else own_lp[a] = 1.0;

    for (int s = 0; s < nslot; s++) {
      double f = 0.0;
      for (int C = 0; C < 8; C += 2)
        f += ASFXRD[ztype_of_slot[s]][C] *
          exp(-1*ASFXRD[ztype_of_slot[s]][C+1] * SinTheta_lambda * SinTheta_lambda);
      f += ASFXRD[ztype_of_slot[s]][8];
      own_asf[s*nown + a] = f;
    }

    a++;
  }

  setup_done = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeXRDFFT::compute_array()
{
  invoked_array = update->ntimestep;

  double t0 = platform::walltime();

  int natoms = group->count(igroup);
  if (natoms == 0) natoms = 1;

  for (int a = 0; a < nown; a++) Fre[a] = Fim[a] = 0.0;

  int ntotal = nmesh[0]*nmesh[1]*nmesh[2];

  for (int s = 0; s < nslot; s++) {

    memset(density_all,0,ntotal*sizeof(FFT_SCALAR));
    spread(s);

    MPI_Reduce_scatter(density_all,density_slab,recvcounts,MPI_FFT_SCALAR,MPI_SUM,world);

    for (int i = 0; i < nfft; i++) {
      work1[2*i] = density_slab[i];
      work1[2*i+1] = (FFT_SCALAR) 0.0;
    }


    if (nslab > 0) fft1->compute(work1,work1,FFT3d::FORWARD);

    // the FFT uses the exp(-i...) convention while compute xrd is defined with
    // exp(+i...), so this yields the complex conjugate of F.  only |F|^2 is
    // reported, so the conjugation is not observable.

    double *asf = &own_asf[s*nown];
    for (int a = 0; a < nown; a++) {
      int idx = own_idx[a];
      double w = own_deconv[a]*asf[a];
      Fre[a] += w*(double)work1[2*idx];
      Fim[a] += w*(double)work1[2*idx+1];
    }
  }

  // every mode is owned by exactly one rank, so each row is contributed once

  for (int n = 0; n < size_array_rows; n++) Iloc[n] = 0.0;
  for (int a = 0; a < nown; a++)
    Iloc[own_row[a]] = own_lp[a]*(Fre[a]*Fre[a] + Fim[a]*Fim[a])/natoms;

  MPI_Allreduce(Iloc,Iall,size_array_rows,MPI_DOUBLE,MPI_SUM,world);

  for (int n = 0; n < size_array_rows; n++) array[n][1] = Iall[n];

  if ((me == 0) && echo)
    utils::logmesg(lmp,"-----\nCompute XRD/FFT id:{} Elapsed time: {:0.2f} s\n-----\n",
                   id,platform::walltime()-t0);
}

/* ----------------------------------------------------------------------
   spread all group atoms of element slot s onto the mesh with unit weight

   every rank spreads onto a full copy of the mesh and the copies are summed by
   the MPI_Reduce_scatter in compute_array().  the mesh is not aligned with the
   simulation box -- when the diffraction cell is smaller than the box, two
   points that are far apart in space fold onto the same mesh cell -- so a
   spatial decomposition of the mesh is not possible.
------------------------------------------------------------------------- */

void ComputeXRDFFT::spread(int s)
{
  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int nx = nmesh[0];
  int ny = nmesh[1];
  int nz = nmesh[2];

  auto *rho0 = new double[order];
  auto *rho1 = new double[order];
  auto *rho2 = new double[order];
  auto *mx = new int[order];
  auto *my = new int[order];
  auto *mz = new int[order];

  double invhalf = 2.0/order;

  for (int ii = 0; ii < nlocal; ii++) {
    if (!(mask[ii] & groupbit)) continue;
    if (slot_of_type[type[ii]-1] != s) continue;

    for (int d = 0; d < 3; d++) {
      int n = nmesh[d];

      // fold into the diffraction cell.  atoms may sit far outside the box in a
      // non-periodic direction, and the cell may be much smaller than the box,
      // so the reduction has to handle arbitrarily large offsets.

      double u = x[ii][d]*mesh_scale[d];
      u -= n*floor(u/n);
      if ((u >= n) || (u < 0.0)) u = 0.0;

      int ngrid = (int) floor(u + 0.5);
      double delta = ngrid - u;

      double *rho = (d == 0) ? rho0 : ((d == 1) ? rho1 : rho2);
      int *m = (d == 0) ? mx : ((d == 1) ? my : mz);

      for (int k = nlower; k <= nupper; k++) {
        double a = (k + delta)*invhalf;
        double sa = 1.0 - a*a;
        if (sa < 0.0) sa = 0.0;
        rho[k-nlower] = bessel_i0(kb_beta[d]*sqrt(sa));
        m[k-nlower] = (((ngrid + k) % n) + n) % n;
      }
    }

    for (int kk = 0; kk < order; kk++) {
      int zoff = mz[kk]*ny;
      double z0 = rho2[kk];
      for (int jj = 0; jj < order; jj++) {
        int yoff = (zoff + my[jj])*nx;
        double y0 = z0*rho1[jj];
        for (int i = 0; i < order; i++)
          density_all[yoff + mx[i]] += (FFT_SCALAR) (y0*rho0[i]);
      }
    }
  }

  delete [] rho0;
  delete [] rho1;
  delete [] rho2;
  delete [] mx;
  delete [] my;
  delete [] mz;
}

/* ----------------------------------------------------------------------
   Fourier transform of the Kaiser-Bessel window at integer mode m of an
   n point mesh:  W*sinh(s)/s with s = sqrt(beta^2 - (pi*W*m/n)^2), W = order,
   continued as W*sin(s)/s where the argument turns imaginary
------------------------------------------------------------------------- */

double ComputeXRDFFT::kb_window(int m, int n, double beta)
{
  double arg = MY_PI*order*((double) m)/n;
  double t = beta*beta - arg*arg;
  if (t > 0.0) {
    double s = sqrt(t);
    return order*sinh(s)/s;
  } else if (t < 0.0) {
    double s = sqrt(-t);
    return order*sin(s)/s;
  }
  return (double) order;
}

/* ----------------------------------------------------------------------
   modified Bessel function of the first kind, order zero

   I0(x) = sum_k (x^2/4)^k / (k!)^2, summed until the terms stop contributing.
   every term is positive, so there is no cancellation and the result is
   accurate to full double precision over the range of shape parameters used
   here.  the commonly used Abramowitz & Stegun 9.8.1/9.8.2 polynomials are
   only accurate to about 2e-7, which would put a floor of that order on the
   spreading weights and prevent the larger stencils from delivering the
   accuracy they should.
------------------------------------------------------------------------- */

double ComputeXRDFFT::bessel_i0(double x)
{
  double t = 0.25*x*x;
  double term = 1.0;
  double sum = 1.0;

  for (int k = 1; k < 200; k++) {
    term *= t/((double)k*(double)k);
    sum += term;
    if (term < 1.0e-17*sum) break;
  }
  return sum;
}

/* ----------------------------------------------------------------------
   only meshes whose dimensions have factors of 2, 3, and 5 are allowed
------------------------------------------------------------------------- */

int ComputeXRDFFT::factorable(int n)
{
  static const int factors[3] = {2,3,5};

  while (n > 1) {
    int found = 0;
    for (int i = 0; i < 3; i++) {
      if (n % factors[i] == 0) {
        n /= factors[i];
        found = 1;
        break;
      }
    }
    if (!found) return 0;
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

double ComputeXRDFFT::memory_usage()
{
  if (!setup_done) return (double)size_array_rows*2 * sizeof(double);

  double bytes = (double)size_array_rows*2 * sizeof(double);        // array
  bytes += (double)size_array_rows*3 * sizeof(int);                 // store_tmp
  bytes += (double)size_array_rows*2 * sizeof(double);              // Iloc,Iall
  bytes += (double)nmesh[0]*nmesh[1]*nmesh[2] * sizeof(FFT_SCALAR); // density_all
  bytes += (double)nfft*3 * sizeof(FFT_SCALAR);                     // density_slab,work1
  bytes += (double)nown*(3+nslot) * sizeof(double);                 // deconv,lp,asf,F
  bytes += (double)nown*2 * sizeof(int);                            // own_row,own_idx
  return bytes;
}
