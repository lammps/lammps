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

   compute xrd samples K = i*g0 + j*g1 + k*g2, with g the reciprocal lattice
   vectors of the cell scaled by the c parameters.  The phase of node index d
   is then 2 pi times the projection g_d.r, which is periodic in that
   projection, so folding it into one period (here called the diffraction cell)
   is exact rather than an approximation, whether that cell is larger than the
   simulation box (c < 1) or smaller than it (c > 1).  A single uniform mesh
   spanning the diffraction cell therefore reproduces the requested nodes
   exactly.  For a triclinic cell the projections are oblique, which is what
   mesh_vec accounts for; the mesh itself is uniform in the node indices, just
   as the kspace styles mesh a triclinic cell uniformly in lamda coordinates.

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
  foot_all(nullptr), recvcounts(nullptr), density_own(nullptr), density_slab(nullptr),
  work1(nullptr), sendbuf(nullptr), recvbuf(nullptr),
  kb_cheb(nullptr), kb_tcheb(nullptr), rho0(nullptr), rho1(nullptr), rho2(nullptr),
  mx(nullptr), my(nullptr), mz(nullptr), slot_atoms(nullptr), slot_start(nullptr),
  sstart(nullptr), scount(nullptr),
  destlist(nullptr), requests(nullptr),
  slot_of_type(nullptr), ztype_of_slot(nullptr), own_row(nullptr), own_idx(nullptr),
  own_deconv(nullptr), own_lp(nullptr), own_asf(nullptr), Fre(nullptr), Fim(nullptr),
  Iloc(nullptr), Iall(nullptr), fft1(nullptr)
{
  nprocs = comm->nprocs;
  setup_done = 0;
  fft_comm = MPI_COMM_NULL;
  nfoot = 0;
  all_full = 0;
  maxbucket = 0;
  maxsend = maxrecv = 0;
  foot_lo[0] = foot_lo[1] = foot_lo[2] = 0;
  foot_n[0] = foot_n[1] = foot_n[2] = 0;

  // the stencil width and the oversampling belong to this style alone, so they
  // are claimed from the words ComputeXRD did not recognize rather than taught
  // to it.  what nothing claims is rejected here.

  order = 7;
  oversample = 2.0;

  int kept = 0;
  for (int i = 0; i < nunclaimed; i++) {
    int iarg = unclaimed[i];

    if (strcmp(arg[iarg],"order") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"compute xrd/fft order",error);
      order = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if ((order < 3) || (order % 2 == 0))
        error->all(FLERR,"Compute XRD/FFT: order must be an odd number of 3 or larger");
      i++;    // the value was recorded as unclaimed as well
      continue;

    } else if (strcmp(arg[iarg],"oversample") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"compute xrd/fft oversample",error);
      oversample = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (oversample < MIN_OVERSAMPLE)
        error->all(FLERR,"Compute XRD/FFT: oversample must be {} or larger",MIN_OVERSAMPLE);
      i++;
      continue;
    }

    unclaimed[kept++] = iarg;
  }

  nunclaimed = kept;
  reject_unclaimed(arg);

  nlower = -(order-1)/2;
  nupper = (order-1)/2;

  nown = 0;
  nslot = 0;

  inv_ksq[0] = 0.0;
  for (int k = 1; k < BESSEL_TERMS; k++) inv_ksq[k] = 1.0/((double)k*(double)k);
}

/* ---------------------------------------------------------------------- */

ComputeXRDFFT::~ComputeXRDFFT()
{
  if (copymode) return;

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

  memory->destroy(density_own);
  memory->destroy(density_slab);
  memory->destroy(work1);
  memory->destroy(sendbuf);
  memory->destroy(recvbuf);
  memory->destroy(foot_all);
  memory->destroy(recvcounts);
  memory->destroy(slot_atoms);
  memory->destroy(slot_start);
  memory->destroy(kb_cheb);
  memory->destroy(kb_tcheb);
  memory->destroy(rho0);
  memory->destroy(rho1);
  memory->destroy(rho2);
  memory->destroy(mx);
  memory->destroy(my);
  memory->destroy(mz);
  memory->destroy(sstart);
  memory->destroy(scount);
  memory->destroy(destlist);
  delete [] requests;
  requests = nullptr;
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

  // Memory::destroy() takes its pointer by reference and nulls it, so only the
  // two allocated another way are cleared here

  nfoot = 0;
  maxsend = maxrecv = 0;
  maxbucket = 0;
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

  // the mesh, the transform and the spreading window follow from the node set,
  // the stencil width and the oversampling, none of which can change once the
  // compute is defined.  only their scaling with the box is refreshed, so that
  // a script with many run commands does not pay to free and rebuild a
  // communicator and an FFT plan before each one.

  if (setup_done) {
    refresh_scaling();
    return;
  }

  deallocate();
  set_grid();
  setup_mesh();

#if defined(FFT_SINGLE)
  if (me == 0)
    error->warning(FLERR,"Compute XRD/FFT with single precision FFTs limits the accuracy "
                   "of weak diffuse intensities");
#endif

  if ((me == 0) && echo) {

    // how much of the mesh a rank actually holds depends on where its atoms
    // are, which is not known until the first invocation, so what is reported
    // here is the whole mesh: the most a rank can be asked for

    double mb = 1.0/1024.0/1024.0;
    double mem = (double)nmesh[0]*nmesh[1]*nmesh[2]*sizeof(FFT_SCALAR)*mb;
    utils::logmesg(lmp,"-----\nCompute XRD/FFT id:{}, # of relp:{}\n"
                   "Kaiser-Bessel order {}, FFT mesh {}x{}x{}, {} element transform(s)\n"
                   "Mesh divided {}x{}x{} over the MPI ranks\n"
                   "Whole mesh {:.4} Mbytes, at most that much on one MPI rank\n-----\n",
                   id,size_array_rows,order,nmesh[0],nmesh[1],nmesh[2],nslot,
                   pgrid[0],pgrid[1],pgrid[2],mem);
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

    // in double, since a large oversample would overflow the count of mesh
    // points and silently give back a mesh smaller than the one asked for

    double nreq = ceil(2.0*oversample*Knmax[d]);
    if (nreq > MAXSMALLINT)
      error->all(FLERR,"Compute XRD/FFT: oversample {} asks for more than {} mesh points in "
                 "one dimension; reduce the 2Theta range, increase the c values, or reduce "
                 "oversample",oversample,MAXSMALLINT);

    int nmin = (int) nreq;
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

  // a rank stores the part of the mesh its own atoms reach, see set_footprint().
  // that is the whole mesh when those atoms wrap it, so this is the most a rank
  // can need.  warn before the mesh, rather than the number of atoms, becomes
  // the limit on the calculation

  double mesh_mb = (double)ntotal*sizeof(FFT_SCALAR)/1024.0/1024.0;
  if ((mesh_mb > 512.0) && (me == 0))
    error->warning(FLERR,"Compute XRD/FFT mesh of {}x{}x{} needs up to {:.4} Mbytes on an MPI "
                   "rank.  A rank stores only the part of the mesh its own atoms reach, which "
                   "shrinks as ranks are added, but reaches the whole mesh when the "
                   "diffraction cell is small.  Increase the c values to explore reciprocal "
                   "space more coarsely",nmesh[0],nmesh[1],nmesh[2],mesh_mb);

  // Beatty et al, IEEE Trans Med Imaging 24, 799 (2005), optimal Kaiser-Bessel
  // shape parameter for a stencil of width order at oversampling sigma

  for (int d = 0; d < 3; d++) {
    double sigma = (0.5*nmesh[d])/Knmax[d];
    double arg = (order*order/(sigma*sigma))*(sigma-0.5)*(sigma-0.5) - 0.8;
    if (arg < 0.0) arg = 0.0;
    kb_beta[d] = MY_PI*sqrt(arg);
    kb_norm[d] = 1.0/bessel_i0(kb_beta[d]);
  }

  build_window();

}

/* ----------------------------------------------------------------------
   expand the spreading window in Chebyshev polynomials of the offset of an
   atom from its nearest mesh point, once per dimension per stencil point
------------------------------------------------------------------------- */

void ComputeXRDFFT::build_window()
{
  memory->grow(kb_cheb,3*order*KB_NCHEB,"xrd/fft:kb_cheb");
  memory->grow(kb_tcheb,KB_NCHEB,"xrd/fft:kb_tcheb");

  double invhalf = 2.0/order;
  auto *f = new double[KB_NCHEB];

  for (int d = 0; d < 3; d++) {
    for (int j = 0; j < order; j++) {
      int k = j + nlower;

      // the window at the Chebyshev points of the offset, which runs over
      // half a mesh spacing either side of the nearest point

      for (int m = 0; m < KB_NCHEB; m++) {
        double u = cos(MY_PI*(m+0.5)/KB_NCHEB);
        double a = (k + 0.5*u)*invhalf;
        double sa = 1.0 - a*a;
        if (sa < 0.0) sa = 0.0;
        f[m] = bessel_i0(kb_beta[d]*sqrt(sa))*kb_norm[d];
      }

      double *c = &kb_cheb[(d*order + j)*KB_NCHEB];
      for (int i = 0; i < KB_NCHEB; i++) {
        double sum = 0.0;
        for (int m = 0; m < KB_NCHEB; m++)
          sum += f[m]*cos(MY_PI*i*(m+0.5)/KB_NCHEB);
        c[i] = 2.0*sum/KB_NCHEB;
      }
      c[0] *= 0.5;
    }
  }

  delete [] f;
}

/* ----------------------------------------------------------------------
   split the ranks over the three dimensions of the mesh

   the surface of a brick is what has to be received, so the grid is chosen to
   make the bricks as cubic as the number of ranks allows.  with no more ranks
   than the mesh has planes the slab decomposition is both cubic enough and
   contiguous, and it lets the reduction fall to a single MPI call, so it is
   preferred there.
------------------------------------------------------------------------- */

void ComputeXRDFFT::set_pgrid()
{
  if (nprocs <= nmesh[2]) {
    pgrid[0] = pgrid[1] = 1;
    pgrid[2] = nprocs;
    is_slab = 1;

  } else {
    double best = -1.0;
    pgrid[0] = pgrid[1] = 1;
    pgrid[2] = nprocs;

    for (int px = 1; px <= nprocs; px++) {
      if (nprocs % px) continue;
      int rest = nprocs/px;
      for (int py = 1; py <= rest; py++) {
        if (rest % py) continue;
        int pz = rest/py;

        // a rank with no mesh of its own still has to send, so grids that
        // leave ranks empty are allowed, but they waste one

        double sx = (double)nmesh[0]/px;
        double sy = (double)nmesh[1]/py;
        double sz = (double)nmesh[2]/pz;
        if ((sx < 1.0) || (sy < 1.0) || (sz < 1.0)) continue;

        // surface of one brick, which is what its owner receives

        double surf = sx*sy + sy*sz + sx*sz;
        if ((best < 0.0) || (surf < best)) {
          best = surf;
          pgrid[0] = px;
          pgrid[1] = py;
          pgrid[2] = pz;
        }
      }
    }
    // with a prime number of ranks larger than every dimension of the mesh
    // there is no grid that gives all of them a piece of it, and the seed
    // above stands: the ranks past the end own nothing and only contribute
    // their own atoms

    is_slab = ((pgrid[0] == 1) && (pgrid[1] == 1)) ? 1 : 0;
  }

  pme[0] = me % pgrid[0];
  pme[1] = (me/pgrid[0]) % pgrid[1];
  pme[2] = me/(pgrid[0]*pgrid[1]);
}

/* ----------------------------------------------------------------------
   the buffers the transform works out of, and the transform itself

   these are the only parts of the setup that depend on where the mesh lives,
   so a derived style that keeps the mesh somewhere else replaces this and
   inherits the rest of setup_mesh()
------------------------------------------------------------------------- */

void ComputeXRDFFT::allocate_mesh()
{
  memory->create(density_slab,MAX(nfft,1),"xrd/fft:density_slab");
  memory->create(work1,MAX(2*nfft,1),"xrd/fft:work1");
  memory->create(rho0,order,"xrd/fft:rho0");
  memory->create(rho1,order,"xrd/fft:rho1");
  memory->create(rho2,order,"xrd/fft:rho2");
  memory->create(mx,order,"xrd/fft:mx");
  memory->create(my,order,"xrd/fft:my");
  memory->create(mz,order,"xrd/fft:mz");

  if (nfft > 0) {
    int tmp;
    fft1 = new FFT3d(lmp,fft_comm,nmesh[0],nmesh[1],nmesh[2],
                     fftlo[0],ffthi[0],fftlo[1],ffthi[1],fftlo[2],ffthi[2],
                     fftlo[0],ffthi[0],fftlo[1],ffthi[1],fftlo[2],ffthi[2],
                     0,0,&tmp,0,0);
  }
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

  set_pgrid();

  for (int d = 0; d < 3; d++) {
    fftlo[d] = pme[d]*nmesh[d]/pgrid[d];
    ffthi[d] = (pme[d]+1)*nmesh[d]/pgrid[d] - 1;
    fftn[d] = ffthi[d] - fftlo[d] + 1;
    if (fftn[d] < 0) fftn[d] = 0;
  }

  nfft = fftn[0]*fftn[1]*fftn[2];

  memory->create(foot_all,6*nprocs,"xrd/fft:foot_all");
  memory->create(recvcounts,nprocs,"xrd/fft:recvcounts");
  memory->create(sstart,nprocs,"xrd/fft:sstart");
  memory->create(scount,nprocs,"xrd/fft:scount");
  memory->create(destlist,nprocs,"xrd/fft:destlist");
  requests = new MPI_Request[nprocs];

  for (int p = 0; p < nprocs; p++) scount[p] = 0;

  // a slab decomposition leaves the blocks of the mesh contiguous and in rank
  // order, which is what MPI_Reduce_scatter needs for the shortcut in
  // fold_reduce()

  if (is_slab)
    for (int p = 0; p < nprocs; p++) {
      int zlo = p*nz/nprocs;
      int zhi = (p+1)*nz/nprocs - 1;
      recvcounts[p] = nx*ny*(zhi-zlo+1);
    }

  // input and output layouts are the same brick, so the modes come back in the
  // layout the intensities are read from.  the transform runs on a
  // communicator holding only the ranks that own part of the mesh, since not
  // every FFT backend accepts an empty block.

  MPI_Comm_split(world,(nfft > 0) ? 0 : MPI_UNDEFINED,me,&fft_comm);

  allocate_mesh();

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
    int iz = ((store_tmp[3*n] % nz) + nz) % nz;
    int iy = ((store_tmp[3*n+1] % ny) + ny) % ny;
    int ix = ((store_tmp[3*n+2] % nx) + nx) % nx;
    if ((ix < fftlo[0]) || (ix > ffthi[0])) continue;
    if ((iy < fftlo[1]) || (iy > ffthi[1])) continue;
    if ((iz < fftlo[2]) || (iz > ffthi[2])) continue;
    nown++;
  }

  memory->create(own_row,MAX(nown,1),"xrd/fft:own_row");
  memory->create(own_idx,MAX(nown,1),"xrd/fft:own_idx");
  memory->create(own_deconv,MAX(nown,1),"xrd/fft:own_deconv");
  memory->create(own_lp,MAX(nown,1),"xrd/fft:own_lp");
  bigint nasf = (bigint)nslot*nown;
  if (nasf > MAXSMALLINT)
    error->one(FLERR,"Compute XRD/FFT: too many reciprocal lattice nodes for {} elements; "
               "reduce the 2Theta range or increase the c values",nslot);
  memory->create(own_asf,(int)MAX(nasf,(bigint)1),"xrd/fft:own_asf");
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
    int iy = ((j % ny) + ny) % ny;
    int ix = ((i % nx) + nx) % nx;
    if ((ix < fftlo[0]) || (ix > ffthi[0])) continue;
    if ((iy < fftlo[1]) || (iy > ffthi[1])) continue;
    if ((iz < fftlo[2]) || (iz > ffthi[2])) continue;

    own_row[a] = n;
    own_idx[a] = ((iz-fftlo[2])*fftn[1] + (iy-fftlo[1]))*fftn[0] + (ix-fftlo[0]);

    // the window transform must be evaluated at the signed mode index: it is
    // even in m but not periodic in m, so the folded index would be wrong.
    // it depends only on the mode and the mesh, so a box change cannot alter it

    own_deconv[a] = 1.0/(kb_window(i,nx,0) * kb_window(j,ny,1) * kb_window(k,nz,2));
    a++;
  }

  refresh_scaling();

  setup_done = 1;
}

/* ----------------------------------------------------------------------
   recompute everything that depends on the length of the reciprocal lattice
   vectors.  the mesh dimensions and the spreading kernel are fixed, as the
   kspace styles fix their grid, but the mesh is rescaled with the box.
------------------------------------------------------------------------- */

void ComputeXRDFFT::refresh_scaling()
{
  for (int d = 0; d < 3; d++)
    for (int e = 0; e < 3; e++) mesh_vec[d][e] = rlv[d][e]*nmesh[d];

  for (int a = 0; a < nown; a++) {
    int n = own_row[a];
    double K[3];
    kvector(rlv,store_tmp[3*n+2],store_tmp[3*n+1],store_tmp[3*n],K);
    double dinv2 = K[0]*K[0] + K[1]*K[1] + K[2]*K[2];
    double SinTheta_lambda = 0.5*sqrt(dinv2);

    // a node driven beyond the limit of the Ewald sphere no longer diffracts

    if (4 < dinv2*lambda*lambda) {
      own_lp[a] = 0.0;
      for (int s = 0; s < nslot; s++) own_asf[(bigint)s*nown + a] = 0.0;
      continue;
    }

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
      own_asf[(bigint)s*nown + a] = f;
    }
  }
}

/* ----------------------------------------------------------------------
   find the part of the mesh this rank spreads into

   the mesh index of an atom is an affine function of its coordinate, folded
   into the mesh.  before the fold the image of a set of atoms that lie in a
   box is one interval per dimension, so this rank only ever touches a
   contiguous run of mesh points, which wraps at most once after the fold.
   The run is found from the atoms themselves rather than from the subdomain,
   so it holds whatever the domain decomposition looks like, whether the box
   has been load balanced, and for atoms that have left the box in a
   non-periodic direction.

   Its length grows with the ratio of the box to the diffraction cell, so for
   coarse reciprocal spacing on few ranks it reaches the whole mesh and every
   rank holds a full copy, as it did before.
------------------------------------------------------------------------- */

int ComputeXRDFFT::minmax_u(double *umin, double *umax)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nany = 0;

  for (int d = 0; d < 3; d++) {
    umin[d] = 0.0;
    umax[d] = 0.0;
  }

  for (int ii = 0; ii < nlocal; ii++) {
    if (!(mask[ii] & groupbit)) continue;
    for (int d = 0; d < 3; d++) {
      double u = x[ii][0]*mesh_vec[d][0] + x[ii][1]*mesh_vec[d][1] +
                 x[ii][2]*mesh_vec[d][2];
      if (!nany) umin[d] = umax[d] = u;
      else if (u < umin[d]) umin[d] = u;
      else if (u > umax[d]) umax[d] = u;
    }
    nany = 1;
  }

  return nany;
}

/* ---------------------------------------------------------------------- */

void ComputeXRDFFT::grow_density_own()
{
  memory->grow(density_own,MAX(nfoot,(bigint)1),"xrd/fft:density_own");
}

/* ---------------------------------------------------------------------- */

void ComputeXRDFFT::set_footprint()
{
  double umin[3], umax[3];
  int nany = minmax_u(umin,umax);

  nfoot = 1;

  for (int d = 0; d < 3; d++) {
    int n = nmesh[d];

    if (!nany) {
      foot_lo[d] = 0;
      foot_n[d] = 0;
      nfoot = 0;
      continue;
    }

    // the interval of mesh points the stencil of these atoms can reach.  the
    // bounds are taken in double, so they are only good to a few ulp of the
    // unfolded index; one extra point on each side absorbs that.

    double glo = floor(umin[d] + 0.5) + nlower - 1.0;
    double ghi = floor(umax[d] + 0.5) + nupper + 1.0;
    double len = ghi - glo + 1.0;

    if (len >= (double) n) {
      foot_lo[d] = 0;
      foot_n[d] = n;
    } else {
      foot_n[d] = (int) len;
      double f = glo - n*floor(glo/n);
      int lo = (int) f;
      if (lo < 0) lo = 0;
      if (lo >= n) lo = n-1;
      foot_lo[d] = lo;
    }
    nfoot *= foot_n[d];
  }

  // the footprint of a rank owning atoms is at least the stencil width, so it
  // cannot overflow an int unless the whole mesh does, which set_grid() rejects

  grow_density_own();

  // a rank that already reaches most of the mesh saves little by sending the
  // pieces of it separately, and MPI reduces a whole mesh better than this can.
  // when that is true of every rank, round the footprints up to the whole mesh
  // and hand the sum to MPI_Reduce_scatter instead.  the threshold only picks
  // between two ways of computing the same thing.

  bigint ntotal = (bigint)nmesh[0]*nmesh[1]*nmesh[2];
  int nearly = (nany && (nfoot*4 > ntotal*3)) ? 1 : 0;
  MPI_Allreduce(&nearly,&all_full,1,MPI_INT,MPI_MIN,world);

  if (all_full) {
    nfoot = ntotal;
    for (int d = 0; d < 3; d++) {
      foot_lo[d] = 0;
      foot_n[d] = nmesh[d];
    }
    grow_density_own();
  }

  bucket_atoms();

  // every rank needs the others' footprints to know what it will be sent

  int mine[6];
  for (int d = 0; d < 3; d++) {
    mine[d] = foot_lo[d];
    mine[3+d] = foot_n[d];
  }
  MPI_Allgather(mine,6,MPI_INT,foot_all,6,MPI_INT,world);
}

/* ----------------------------------------------------------------------
   list the group atoms of each element together

   spread() is called once per element, and walked every atom each time to pick
   out the ones of that element, so a system of several elements tested every
   atom several times over to spread each one once.
------------------------------------------------------------------------- */

void ComputeXRDFFT::bucket_atoms()
{
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  memory->grow(slot_start,nslot+1,"xrd/fft:slot_start");
  for (int s = 0; s <= nslot; s++) slot_start[s] = 0;

  // count, then place; the counts become the starts by a running sum

  for (int ii = 0; ii < nlocal; ii++) {
    if (!(mask[ii] & groupbit)) continue;
    slot_start[slot_of_type[type[ii]-1] + 1]++;
  }

  for (int s = 0; s < nslot; s++) slot_start[s+1] += slot_start[s];

  if (slot_start[nslot] > maxbucket) {
    maxbucket = slot_start[nslot];
    memory->grow(slot_atoms,maxbucket,"xrd/fft:slot_atoms");
  }

  auto *fill = new int[nslot];
  for (int s = 0; s < nslot; s++) fill[s] = slot_start[s];

  for (int ii = 0; ii < nlocal; ii++) {
    if (!(mask[ii] & groupbit)) continue;
    slot_atoms[fill[slot_of_type[type[ii]-1]]++] = ii;
  }

  delete [] fill;
}

/* ----------------------------------------------------------------------
   sum every rank's footprint into the z slab decomposition of the mesh

   the element index is the message tag.  a rank whose brick is empty takes no
   part in the transform, which is the only thing that would hold it with the
   others, so it can reach the next element while they are still receiving this
   one; without the tag those messages would be taken for this element's.

   this replaces a reduction over a full copy of the mesh held by every rank.
   A rank sends each destination only the mesh planes that destination owns,
   so what it moves is the size of its footprint rather than the size of the
   mesh, and that shrinks as ranks are added.
------------------------------------------------------------------------- */

void ComputeXRDFFT::fold_reduce(int tag)
{
  // a slab decomposition leaves every rank holding the same contiguous blocks
  // the reduction wants, so when they also hold the whole mesh MPI can do all
  // of it in one call

  if (all_full && is_slab) {
    MPI_Reduce_scatter(density_own,density_slab,recvcounts,MPI_FFT_SCALAR,MPI_SUM,world);
    return;
  }

  for (int i = 0; i < nfft; i++) density_slab[i] = (FFT_SCALAR) 0.0;

  if (nfoot > maxsend) {
    maxsend = nfoot;
    memory->grow(sendbuf,(int)maxsend,"xrd/fft:sendbuf");
  }

  // the ranks whose bricks this footprint reaches.  along each dimension the
  // footprint is one run of mesh points, or two once it wraps, and each run
  // covers a consecutive range of the grid, so the destinations are found
  // without looking at every rank

  int ndest = 0;

  for (int qz = 0; qz < pgrid[2]; qz++) {
    for (int qy = 0; qy < pgrid[1]; qy++) {
      for (int qx = 0; qx < pgrid[0]; qx++) {
        int q = (qz*pgrid[1] + qy)*pgrid[0] + qx;
        bigint n = brick_count(q);
        if (n == 0) continue;
        destlist[ndest++] = q;
        scount[q] = (int) n;
      }
    }
  }

  int nout = 0;
  for (int i = 0; i < ndest; i++) {
    int q = destlist[i];
    sstart[q] = nout;
    nout += scount[q];
  }

  int nsend = 0;
  int selfoff = -1;

  for (int i = 0; i < ndest; i++) {
    int q = destlist[i];
    pack_brick(q,&sendbuf[sstart[q]]);
    if (q == me) selfoff = sstart[q];
    else MPI_Isend(&sendbuf[sstart[q]],scount[q],MPI_FFT_SCALAR,q,tag,world,
                   &requests[nsend++]);
  }

  for (int i = 0; i < ndest; i++) scount[destlist[i]] = 0;

  // this rank's own contribution needs no message

  if (selfoff >= 0) unpack_brick(&sendbuf[selfoff],me);

  // count the messages to expect, then take them as they arrive

  int nrecv = 0;
  bigint biggest = 0;

  if (nfft > 0) {
    for (int q = 0; q < nprocs; q++) {
      if (q == me) continue;
      const int *lo = &foot_all[6*q];
      const int *ln = &foot_all[6*q+3];
      if ((ln[0] == 0) || (ln[1] == 0) || (ln[2] == 0)) continue;
      bigint n = 1;
      for (int d = 0; d < 3; d++) {
        int seg[2][3];
        int ns = segments(lo[d],ln[d],nmesh[d],fftlo[d],ffthi[d],seg);
        int len = 0;
        for (int t = 0; t < ns; t++) len += seg[t][2];
        n *= len;
      }
      if (n == 0) continue;
      nrecv++;
      if (n > biggest) biggest = n;
    }
  }

  if (biggest > maxrecv) {
    maxrecv = biggest;
    memory->grow(recvbuf,(int)maxrecv,"xrd/fft:recvbuf");
  }

  for (int i = 0; i < nrecv; i++) {
    MPI_Status status;
    MPI_Recv(recvbuf,(int)maxrecv,MPI_FFT_SCALAR,MPI_ANY_SOURCE,tag,world,&status);
    unpack_brick(recvbuf,status.MPI_SOURCE);
  }

  if (nsend) MPI_Waitall(nsend,requests,MPI_STATUSES_IGNORE);
}

/* ----------------------------------------------------------------------
   where a footprint and a brick overlap along one dimension

   the footprint holds mesh points (lo + i) mod n for i below its length, which
   is one run of the mesh or two once it wraps.  each run is intersected with
   the brick, giving at most two pieces, returned as (index within the
   footprint, mesh point, length).
------------------------------------------------------------------------- */

int ComputeXRDFFT::segments(int flo, int flen, int n, int blo, int bhi, int (*seg)[3])
{
  int ns = 0;
  if ((flen <= 0) || (bhi < blo)) return 0;

  // run one: mesh points flo .. min(flo+flen,n)-1, starting the footprint
  // run two, when the footprint wraps: mesh points 0 .. flo+flen-n-1

  int glo[2], ghi[2], loc[2];
  int nrun = 0;

  glo[nrun] = flo;
  ghi[nrun] = ((flo + flen) < n) ? flo + flen - 1 : n - 1;
  loc[nrun] = 0;
  nrun++;

  if (flo + flen > n) {
    glo[nrun] = 0;
    ghi[nrun] = flo + flen - n - 1;
    loc[nrun] = n - flo;
    nrun++;
  }

  for (int r = 0; r < nrun; r++) {
    int gs = (glo[r] > blo) ? glo[r] : blo;
    int ge = (ghi[r] < bhi) ? ghi[r] : bhi;
    if (ge < gs) continue;
    seg[ns][0] = loc[r] + (gs - glo[r]);
    seg[ns][1] = gs;
    seg[ns][2] = ge - gs + 1;
    ns++;
  }

  return ns;
}

/* ----------------------------------------------------------------------
   how many mesh points this rank's footprint contributes to rank q's brick
------------------------------------------------------------------------- */

bigint ComputeXRDFFT::brick_count(int q)
{
  if ((foot_n[0] == 0) || (foot_n[1] == 0) || (foot_n[2] == 0)) return 0;

  bigint total = 1;
  for (int d = 0; d < 3; d++) {
    int blo = qgrid_lo(q,d);
    int bhi = qgrid_hi(q,d);
    int seg[2][3];
    int ns = segments(foot_lo[d],foot_n[d],nmesh[d],blo,bhi,seg);
    int len = 0;
    for (int t = 0; t < ns; t++) len += seg[t][2];
    if (len == 0) return 0;
    total *= len;
  }
  return total;
}

/* ----------------------------------------------------------------------
   copy the part of this rank's footprint that belongs to rank q into buf,
   and add a part received from rank q of the mesh into this rank's brick.
   both walk the overlap in the same order, so no description of the pieces
   has to travel with them
------------------------------------------------------------------------- */

void ComputeXRDFFT::pack_brick(int q, FFT_SCALAR *buf)
{
  int sz[2][3], sy[2][3], sx[2][3];
  int nz = segments(foot_lo[2],foot_n[2],nmesh[2],qgrid_lo(q,2),qgrid_hi(q,2),sz);
  int ny = segments(foot_lo[1],foot_n[1],nmesh[1],qgrid_lo(q,1),qgrid_hi(q,1),sy);
  int nx = segments(foot_lo[0],foot_n[0],nmesh[0],qgrid_lo(q,0),qgrid_hi(q,0),sx);

  const int fx = foot_n[0], fy = foot_n[1];
  bigint m = 0;

  for (int a = 0; a < nz; a++)
    for (int k = 0; k < sz[a][2]; k++) {
      bigint zoff = (bigint)(sz[a][0] + k)*fy;
      for (int b = 0; b < ny; b++)
        for (int j = 0; j < sy[b][2]; j++) {
          const FFT_SCALAR *row = &density_own[(zoff + sy[b][0] + j)*fx];
          for (int c = 0; c < nx; c++) {
            memcpy(&buf[m],&row[sx[c][0]],sx[c][2]*sizeof(FFT_SCALAR));
            m += sx[c][2];
          }
        }
    }
}

/* ---------------------------------------------------------------------- */

void ComputeXRDFFT::unpack_brick(const FFT_SCALAR *buf, int q)
{
  const int *lo = &foot_all[6*q];
  const int *ln = &foot_all[6*q+3];
  if ((ln[0] == 0) || (ln[1] == 0) || (ln[2] == 0)) return;

  int sz[2][3], sy[2][3], sx[2][3];
  int nz = segments(lo[2],ln[2],nmesh[2],fftlo[2],ffthi[2],sz);
  int ny = segments(lo[1],ln[1],nmesh[1],fftlo[1],ffthi[1],sy);
  int nx = segments(lo[0],ln[0],nmesh[0],fftlo[0],ffthi[0],sx);

  bigint m = 0;

  for (int a = 0; a < nz; a++)
    for (int k = 0; k < sz[a][2]; k++) {
      bigint zoff = (bigint)(sz[a][1] + k - fftlo[2])*fftn[1];
      for (int b = 0; b < ny; b++)
        for (int j = 0; j < sy[b][2]; j++) {
          FFT_SCALAR *row = &density_slab[(zoff + sy[b][1] + j - fftlo[1])*fftn[0]];
          for (int c = 0; c < nx; c++) {
            FFT_SCALAR *dst = &row[sx[c][1] - fftlo[0]];
            for (int i = 0; i < sx[c][2]; i++) dst[i] += buf[m++];
          }
        }
    }
}

/* ---------------------------------------------------------------------- */

void ComputeXRDFFT::compute_array()
{
  invoked_array = update->ntimestep;

  // rescale the mesh and everything derived from |K| if the box has changed

  if (update_reciprocal()) refresh_scaling();

  double t0 = platform::walltime();

  bigint natoms = group->count(igroup);
  if (natoms == 0) natoms = 1;

  for (int a = 0; a < nown; a++) Fre[a] = Fim[a] = 0.0;

  // which part of the mesh this rank spreads into depends on where its atoms
  // are, so it is found once per invocation and shared by all the elements

  set_footprint();

  for (int s = 0; s < nslot; s++) {

    if (nfoot) memset(density_own,0,(size_t)nfoot*sizeof(FFT_SCALAR));
    spread(s);

    fold_reduce(s);

    for (int i = 0; i < nfft; i++) {
      work1[2*i] = density_slab[i];
      work1[2*i+1] = (FFT_SCALAR) 0.0;
    }


    if (nfft > 0) fft1->compute(work1,work1,FFT3d::FORWARD);

    // the FFT uses the exp(-i...) convention while compute xrd is defined with
    // exp(+i...), so this yields the complex conjugate of F.  only |F|^2 is
    // reported, so the conjugation is not observable.

    double *asf = &own_asf[(bigint)s*nown];
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

   the mesh is not aligned with the simulation box: when the diffraction cell is
   smaller than the box, two points far apart in space fold onto the same mesh
   cell.  A rank therefore spreads into the run of mesh points its own atoms
   reach, found by set_footprint(), and fold_reduce() then sums those runs into
   the bricks of the mesh the ranks own.
------------------------------------------------------------------------- */

void ComputeXRDFFT::spread(int s)
{
  double **x = atom->x;

  if (!nfoot) return;

  // hoisted out of the loops below: the mesh dimensions, the start of this
  // rank's footprint and the stride of a row and a plane of it

  const int nx = nmesh[0], ny = nmesh[1], nz = nmesh[2];
  const int fx = foot_n[0], fy = foot_n[1];
  const int lox = foot_lo[0], loy = foot_lo[1], loz = foot_lo[2];

  for (int b = slot_start[s]; b < slot_start[s+1]; b++) {
    const int ii = slot_atoms[b];

    const double xi = x[ii][0], yi = x[ii][1], zi = x[ii][2];

    for (int d = 0; d < 3; d++) {
      const int n = (d == 0) ? nx : ((d == 1) ? ny : nz);
      const int flen = (d == 0) ? fx : ((d == 1) ? fy : foot_n[2]);
      const int flo = (d == 0) ? lox : ((d == 1) ? loy : loz);

      // fold into the diffraction cell.  atoms may sit far outside the box in a
      // non-periodic direction, and the cell may be much smaller than the box,
      // so the reduction has to handle arbitrarily large offsets.

      double u = xi*mesh_vec[d][0] + yi*mesh_vec[d][1] + zi*mesh_vec[d][2];
      u -= n*floor(u/n);
      if ((u >= n) || (u < 0.0)) u = 0.0;

      const int ngrid = (int) floor(u + 0.5);
      const double delta = ngrid - u;

      double *rho = (d == 0) ? rho0 : ((d == 1) ? rho1 : rho2);
      int *m = (d == 0) ? mx : ((d == 1) ? my : mz);

      // the window at this offset, from its Chebyshev expansion.  the
      // polynomials are built once for the offset and then shared by every
      // point of the stencil.  the expansion carries the normalization of the
      // window to a peak of one: unnormalized it peaks at I0(beta), which is
      // 1e12 at order 13, and the product over the three dimensions would put
      // 1e36 in every mesh cell and overflow a single precision mesh.  the same
      // factor divides out of the window transform in kb_window().

      const double uu = 2.0*delta;
      const double uu2 = 2.0*uu;
      kb_tcheb[0] = 1.0;
      kb_tcheb[1] = uu;
      for (int i = 2; i < KB_NCHEB; i++)
        kb_tcheb[i] = uu2*kb_tcheb[i-1] - kb_tcheb[i-2];

      const double *cd = &kb_cheb[d*order*KB_NCHEB];
      for (int j = 0; j < order; j++) {
        const double *cj = cd + j*KB_NCHEB;
        double w = 0.0;
        for (int i = 0; i < KB_NCHEB; i++) w += cj[i]*kb_tcheb[i];
        rho[j] = w;
      }

      // index within this rank's footprint rather than within the mesh.
      // set_footprint() sized the footprint from these same atoms, so the
      // stencil cannot reach past either end of it

      int g = ngrid + nlower - flo;
      g %= n;
      if (g < 0) g += n;

      for (int j = 0; j < order; j++) {
        if (g >= flen)
          error->one(FLERR,"Compute XRD/FFT: atom outside the mesh footprint of its "
                     "MPI rank; this is an internal error, please report it");
        m[j] = g;
        if (++g == n) g = 0;
      }
    }

    for (int kk = 0; kk < order; kk++) {
      const bigint zoff = (bigint)(mz[kk]*fy);
      const double z0 = rho2[kk];
      for (int jj = 0; jj < order; jj++) {
        FFT_SCALAR *row = &density_own[(zoff + my[jj])*fx];
        const double y0 = z0*rho1[jj];
        for (int i = 0; i < order; i++)
          row[mx[i]] += (FFT_SCALAR) (y0*rho0[i]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Fourier transform of the Kaiser-Bessel window at integer mode m of an
   n point mesh:  W*sinh(s)/s with s = sqrt(beta^2 - (pi*W*m/n)^2), W = order,
   continued as W*sin(s)/s where the argument turns imaginary

   scaled by the same 1/I0(beta) as the window itself, see spread().  that
   factor is kb_norm[d], so it is not recomputed here
------------------------------------------------------------------------- */

double ComputeXRDFFT::kb_window(int m, int n, int d)
{
  const double beta = kb_beta[d];
  double arg = MY_PI*order*((double) m)/n;
  double t = beta*beta - arg*arg;
  double w;

  if (t > 0.0) {
    double s = sqrt(t);
    w = order*sinh(s)/s;
  } else if (t < 0.0) {
    double s = sqrt(-t);
    w = order*sin(s)/s;
  } else w = (double) order;

  // the same 1/I0(beta) that normalizes the window itself, already computed

  return w*kb_norm[d];
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

  for (int k = 1; k < BESSEL_TERMS; k++) {
    term *= t*inv_ksq[k];
    sum += term;
    if (term < 1.0e-17*sum) return sum;
  }

  // the widest stencils have not converged in that many terms, so finish the
  // sum rather than truncate it

  for (int k = BESSEL_TERMS; k < 400; k++) {
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
  bytes += (double)nfoot * sizeof(FFT_SCALAR);                      // density_own
  bytes += (double)(maxsend+maxrecv) * sizeof(FFT_SCALAR);          // exchange buffers
  bytes += (double)nfft*3 * sizeof(FFT_SCALAR);                     // density_slab,work1
  bytes += (double)nprocs*7 * sizeof(int);                          // foot_all,recvcounts
  bytes += (double)nown*(4+nslot) * sizeof(double);                 // deconv,lp,asf,Fre,Fim
  bytes += (double)nown*2 * sizeof(int);                            // own_row,own_idx
  return bytes;
}
