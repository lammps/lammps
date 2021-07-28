// clang-format off
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
   Contributing author: Ray Shan (Sandia)
     Based on fix qeq/reax by H. Metin Aktulga
------------------------------------------------------------------------- */

#include "fix_qeq.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "pair.h"
#include "suffix.h"
#include "text_file_reader.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXLINE 1024

class parser_error : public std::exception {
  std::string message;
public:
  parser_error(const std::string &mesg) { message = mesg; }
  const char *what() const noexcept { return message.c_str(); }
};

/* ---------------------------------------------------------------------- */

FixQEq::FixQEq(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), list(nullptr), chi(nullptr), eta(nullptr),
  gamma(nullptr), zeta(nullptr), zcore(nullptr), chizj(nullptr), shld(nullptr),
  s(nullptr), t(nullptr), s_hist(nullptr), t_hist(nullptr), Hdia_inv(nullptr), b_s(nullptr),
  b_t(nullptr), p(nullptr), q(nullptr), r(nullptr), d(nullptr),
  qf(nullptr), q1(nullptr), q2(nullptr), qv(nullptr)
{
  if (narg < 8) error->all(FLERR,"Illegal fix qeq command");

  scalar_flag = 1;
  extscalar = 0;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  cutoff = utils::numeric(FLERR,arg[4],false,lmp);
  tolerance = utils::numeric(FLERR,arg[5],false,lmp);
  maxiter = utils::inumeric(FLERR,arg[6],false,lmp);
  maxwarn = 1;
  matvecs = 0;

  // check for sane arguments
  if ((nevery <= 0) || (cutoff <= 0.0) || (tolerance <= 0.0) || (maxiter <= 0))
    error->all(FLERR,"Illegal fix qeq command");

  alpha = 0.20;
  swa = 0.0;
  swb = cutoff;

  shld = nullptr;

  nlocal = n_cap = 0;
  nall = nmax = 0;
  m_fill = m_cap = 0;
  pack_flag = 0;
  s = nullptr;
  t = nullptr;
  nprev = 5;

  Hdia_inv = nullptr;
  b_s = nullptr;
  b_t = nullptr;

  // CG
  p = nullptr;
  q = nullptr;
  r = nullptr;
  d = nullptr;

  // H matrix
  H.firstnbr = nullptr;
  H.numnbrs = nullptr;
  H.jlist = nullptr;
  H.val = nullptr;

  // others
  cutoff_sq = cutoff*cutoff;
  chizj = nullptr;
  qf = nullptr;
  q1 = nullptr;
  q2 = nullptr;
  streitz_flag = 0;
  reax_flag = 0;
  qv = nullptr;

  comm_forward = comm_reverse = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  s_hist = t_hist = nullptr;
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  for (int i = 0; i < atom->nmax; i++)
    for (int j = 0; j < nprev; ++j)
      s_hist[i][j] = t_hist[i][j] = atom->q[i];

  if (strcmp(arg[7],"coul/streitz") == 0) {
    streitz_flag = 1;
  } else if (utils::strmatch(arg[7],"^reax..")) {
    reax_flag = 1;
  } else {
    read_file(arg[7]);
  }
}

/* ---------------------------------------------------------------------- */

FixQEq::~FixQEq()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,Atom::GROW);

  memory->destroy(s_hist);
  memory->destroy(t_hist);

  deallocate_storage();
  deallocate_matrix();

  memory->destroy(shld);

  if (!streitz_flag && !reax_flag) {
    memory->destroy(chi);
    memory->destroy(eta);
    memory->destroy(gamma);
    memory->destroy(zeta);
    memory->destroy(zcore);
  }
}

/* ---------------------------------------------------------------------- */

int FixQEq::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixQEq::allocate_storage()
{
  nmax = atom->nmax;

  memory->create(s,nmax,"qeq:s");
  memory->create(t,nmax,"qeq:t");

  memory->create(Hdia_inv,nmax,"qeq:Hdia_inv");
  memory->create(b_s,nmax,"qeq:b_s");
  memory->create(b_t,nmax,"qeq:b_t");

  memory->create(p,nmax,"qeq:p");
  memory->create(q,nmax,"qeq:q");
  memory->create(r,nmax,"qeq:r");
  memory->create(d,nmax,"qeq:d");

  memory->create(chizj,nmax,"qeq:chizj");
  memory->create(qf,nmax,"qeq:qf");
  memory->create(q1,nmax,"qeq:q1");
  memory->create(q2,nmax,"qeq:q2");

  memory->create(qv,nmax,"qeq:qv");
}

/* ---------------------------------------------------------------------- */

void FixQEq::deallocate_storage()
{
  memory->destroy(s);
  memory->destroy(t);

  memory->destroy(Hdia_inv);
  memory->destroy(b_s);
  memory->destroy(b_t);

  memory->destroy(p);
  memory->destroy(q);
  memory->destroy(r);
  memory->destroy(d);

  memory->destroy(chizj);
  memory->destroy(qf);
  memory->destroy(q1);
  memory->destroy(q2);

  memory->destroy(qv);
}

/* ---------------------------------------------------------------------- */

void FixQEq::reallocate_storage()
{
  deallocate_storage();
  allocate_storage();
  init_storage();
}

/* ---------------------------------------------------------------------- */

void FixQEq::allocate_matrix()
{
  int i,ii,inum,m;
  int *ilist, *numneigh;

  int mincap;
  double safezone;

  mincap = MIN_CAP;
  safezone = SAFE_ZONE;

  nlocal = atom->nlocal;
  n_cap = MAX((int)(nlocal * safezone), mincap);
  nall = atom->nlocal + atom->nghost;

  // determine the total space for the H matrix

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;

  m = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    m += numneigh[i];
  }
  m_cap = MAX((int)(m * safezone), mincap * MIN_NBRS);

  H.n = n_cap;
  H.m = m_cap;
  memory->create(H.firstnbr,n_cap,"qeq:H.firstnbr");
  memory->create(H.numnbrs,n_cap,"qeq:H.numnbrs");
  memory->create(H.jlist,m_cap,"qeq:H.jlist");
  memory->create(H.val,m_cap,"qeq:H.val");
}

/* ---------------------------------------------------------------------- */

void FixQEq::deallocate_matrix()
{
  memory->destroy(H.firstnbr);
  memory->destroy(H.numnbrs);
  memory->destroy(H.jlist);
  memory->destroy(H.val);
}

/* ---------------------------------------------------------------------- */

void FixQEq::reallocate_matrix()
{
  deallocate_matrix();
  allocate_matrix();
}

/* ---------------------------------------------------------------------- */

double FixQEq::compute_scalar()
{
  return matvecs;
}

/* ---------------------------------------------------------------------- */

void FixQEq::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixQEq::setup_pre_force(int vflag)
{
  if (force->newton_pair == 0)
    error->all(FLERR,"QEQ with 'newton pair off' not supported");

  if (force->pair) {
    if (force->pair->suffix_flag & (Suffix::INTEL|Suffix::GPU))
      error->all(FLERR,"QEQ is not compatiple with suffix version "
                 "of pair style");
  }

  deallocate_storage();
  allocate_storage();

  init_storage();

  deallocate_matrix();
  allocate_matrix();

  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEq::setup_pre_force_respa(int vflag, int ilevel)
{
  if (ilevel < nlevels_respa-1) return;
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEq::init_storage()
{
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;

  for (int i = 0; i < nall; i++) {
    Hdia_inv[i] = 1. / eta[atom->type[i]];
    b_s[i] = -chi[atom->type[i]];
    b_t[i] = -1.0;
    s[i] = t[i] = atom->q[i];

    chizj[i] = 0.0;
    qf[i] = 0.0;
    q1[i] = 0.0;
    q2[i] = 0.0;

    qv[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixQEq::pre_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEq::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

int FixQEq::CG(double *b, double *x)
{
  int  loop, i, ii, inum, *ilist;
  double tmp, alfa, beta, b_norm;
  double sig_old, sig_new;

  inum = list->inum;
  ilist = list->ilist;

  pack_flag = 1;
  sparse_matvec(&H, x, q);
  comm->reverse_comm_fix(this);

  vector_sum(r , 1.,  b, -1., q, inum);

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      d[i] = r[i] * Hdia_inv[i];
    else d[i] = 0.0;
  }

  b_norm = parallel_norm(b, inum);
  sig_new = parallel_dot(r, d, inum);

  for (loop = 1; loop < maxiter && sqrt(sig_new)/b_norm > tolerance; ++loop) {
    comm->forward_comm_fix(this);
    sparse_matvec(&H, d, q);
    comm->reverse_comm_fix(this);

    tmp = parallel_dot(d, q, inum);
    alfa = sig_new / tmp;

    vector_add(x, alfa, d, inum);
    vector_add(r, -alfa, q, inum);

    for (ii = 0; ii < inum; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit)
        p[i] = r[i] * Hdia_inv[i];
    }

    sig_old = sig_new;
    sig_new = parallel_dot(r, p, inum);

    beta = sig_new / sig_old;
    vector_sum(d, 1., p, beta, d, inum);
  }

  if ((comm->me == 0) && maxwarn && (loop >= maxiter))
    error->warning(FLERR,"Fix qeq CG convergence failed ({}) after {} "
                   "iterations at step {}",sqrt(sig_new)/b_norm,loop,
                   update->ntimestep);
  return loop;
}


/* ---------------------------------------------------------------------- */

void FixQEq::sparse_matvec(sparse_matrix *A, double *x, double *b)
{
  int i, j, itr_j;

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nlocal; ++i) {
    if (atom->mask[i] & groupbit)
      b[i] = eta[atom->type[i]] * x[i];
  }

  for (i = nlocal; i < nall; ++i) {
    if (atom->mask[i] & groupbit)
      b[i] = 0;
  }

  for (i = 0; i < nlocal; ++i) {
    if (atom->mask[i] & groupbit) {
      for (itr_j=A->firstnbr[i]; itr_j<A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
        j = A->jlist[itr_j];
        b[i] += A->val[itr_j] * x[j];
        b[j] += A->val[itr_j] * x[i];
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixQEq::calculate_Q()
{
  int i, k, inum, ii;
  int *ilist;
  double u, s_sum, t_sum;
  double *q = atom->q;

  inum = list->inum;
  ilist = list->ilist;

  s_sum = parallel_vector_acc(s, inum);
  t_sum = parallel_vector_acc(t, inum);
  u = s_sum / t_sum;

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      q[i] = s[i] - u * t[i];

      for (k = 4; k > 0; --k) {
        s_hist[i][k] = s_hist[i][k-1];
        t_hist[i][k] = t_hist[i][k-1];
      }
      s_hist[i][0] = s[i];
      t_hist[i][0] = t[i];
    }
  }

  pack_flag = 4;
  comm->forward_comm_fix(this); //Dist_vector(atom->q);
}

/* ---------------------------------------------------------------------- */

int FixQEq::pack_forward_comm(int n, int *list, double *buf,
                          int /*pbc_flag*/, int * /*pbc*/)
{
  int m;

  if (pack_flag == 1)
    for (m = 0; m < n; m++) buf[m] = d[list[m]];
  else if (pack_flag == 2)
    for (m = 0; m < n; m++) buf[m] = s[list[m]];
  else if (pack_flag == 3)
    for (m = 0; m < n; m++) buf[m] = t[list[m]];
  else if (pack_flag == 4)
    for (m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  else m = 0;

  return m;
}

/* ---------------------------------------------------------------------- */

void FixQEq::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if (pack_flag == 1)
    for (m = 0, i = first; m < n; m++, i++) d[i] = buf[m];
  else if (pack_flag == 2)
    for (m = 0, i = first; m < n; m++, i++) s[i] = buf[m];
  else if (pack_flag == 3)
    for (m = 0, i = first; m < n; m++, i++) t[i] = buf[m];
  else if (pack_flag == 4)
    for (m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
}

/* ---------------------------------------------------------------------- */

int FixQEq::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  for (m = 0, i = first; m < n; m++, i++) buf[m] = q[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void FixQEq::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m;

  for (m = 0; m < n; m++) q[list[m]] += buf[m];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixQEq::memory_usage()
{
  double bytes;

  bytes = (double)atom->nmax*nprev*2 * sizeof(double); // s_hist & t_hist
  bytes += (double)atom->nmax*11 * sizeof(double); // storage
  bytes += (double)n_cap*2 * sizeof(int); // matrix...
  bytes += (double)m_cap * sizeof(int);
  bytes += (double)m_cap * sizeof(double);

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

void FixQEq::grow_arrays(int nmax)
{
  memory->grow(s_hist,nmax,nprev,"qeq:s_hist");
  memory->grow(t_hist,nmax,nprev,"qeq:t_hist");
}

/* ----------------------------------------------------------------------
   copy values within fictitious charge arrays
------------------------------------------------------------------------- */

void FixQEq::copy_arrays(int i, int j, int /*delflag*/)
{
  for (int m = 0; m < nprev; m++) {
    s_hist[j][m] = s_hist[i][m];
    t_hist[j][m] = t_hist[i][m];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixQEq::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nprev; m++) buf[m] = s_hist[i][m];
  for (int m = 0; m < nprev; m++) buf[nprev+m] = t_hist[i][m];
  return nprev*2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixQEq::unpack_exchange(int n, double *buf)
{
  for (int m = 0; m < nprev; m++) s_hist[n][m] = buf[m];
  for (int m = 0; m < nprev; m++) t_hist[n][m] = buf[nprev+m];
  return nprev*2;
}

/* ---------------------------------------------------------------------- */

double FixQEq::parallel_norm(double *v, int n)
{
  int  i;
  double my_sum, norm_sqr;

  int ii;
  int *ilist;

  ilist = list->ilist;

  my_sum = 0.0;
  norm_sqr = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_sum += v[i]*v[i];
  }

  MPI_Allreduce(&my_sum, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, world);

  return sqrt(norm_sqr);
}

/* ---------------------------------------------------------------------- */

double FixQEq::parallel_dot(double *v1, double *v2, int n)
{
  int  i;
  double my_dot, res;

  int ii;
  int *ilist;

  ilist = list->ilist;

  my_dot = 0.0;
  res = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_dot += v1[i] * v2[i];
  }

  MPI_Allreduce(&my_dot, &res, 1, MPI_DOUBLE, MPI_SUM, world);

  return res;
}

/* ---------------------------------------------------------------------- */

double FixQEq::parallel_vector_acc(double *v, int n)
{
  int  i;
  double my_acc, res;

  int ii;
  int *ilist;

  ilist = list->ilist;

  my_acc = 0.0;
  res = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_acc += v[i];
  }

  MPI_Allreduce(&my_acc, &res, 1, MPI_DOUBLE, MPI_SUM, world);

  return res;
}

/* ---------------------------------------------------------------------- */

void FixQEq::vector_sum(double* dest, double c, double* v,
                                double d, double* y, int k)
{
  int kk;
  int *ilist;

  ilist = list->ilist;

  for (--k; k>=0; --k) {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] = c * v[kk] + d * y[kk];
  }
}

/* ---------------------------------------------------------------------- */

void FixQEq::vector_add(double* dest, double c, double* v, int k)
{
  int kk;
  int *ilist;

  ilist = list->ilist;

  for (--k; k>=0; --k) {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] += c * v[kk];
  }
}

/* ---------------------------------------------------------------------- */

void FixQEq::read_file(char *file)
{
  const int ntypes = atom->ntypes;

  memory->create(chi,ntypes+1,"qeq:chi");
  memory->create(eta,ntypes+1,"qeq:eta");
  memory->create(gamma,ntypes+1,"qeq:gamma");
  memory->create(zeta,ntypes+1,"qeq:zeta");
  memory->create(zcore,ntypes+1,"qeq:zcore");

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  if (comm->me == 0) {
    int *setflag = new int[ntypes+1];
    for (int n=0; n <= ntypes; ++n) {
      setflag[n] = 0;
      chi[n] = eta[n] = gamma[n] = zeta[n] = zcore[n] = 0.0;
    }

    try {
      int nlo,nhi;
      double val;

      FILE *fp = utils::open_potential(file,lmp,nullptr);
      if (fp == nullptr)
        throw parser_error(fmt::format("Cannot open fix qeq parameter file {}:"
                                       " {}", file,utils::getsyserror()));
      TextFileReader reader(fp, "qeq parameter");

      while (1) {
        auto values = reader.next_values(0);

        if (values.count() == 0) continue;
        if (values.count() < 6)
          throw parser_error("Invalid qeq parameter file");

        auto word = values.next_string();
        utils::bounds(FLERR,word,1,ntypes,nlo,nhi,nullptr);
        if ((nlo < 0) || (nhi < 0))
          throw parser_error("Invalid atom type range");

        val = values.next_double();
        for (int n=nlo; n <= nhi; ++n) chi[n] = val;
        val = values.next_double();
        for (int n=nlo; n <= nhi; ++n) eta[n] = val;
        val = values.next_double();
        for (int n=nlo; n <= nhi; ++n) gamma[n] = val;
        val = values.next_double();
        for (int n=nlo; n <= nhi; ++n) zeta[n] = val;
        val = values.next_double();
        for (int n=nlo; n <= nhi; ++n) zcore[n] = val;
        for (int n=nlo; n <= nhi; ++n) setflag[n] = 1;
      }
    } catch (EOFException &e) {
      ; // catch and ignore to exit loop
    } catch (std::exception &e) {
      error->one(FLERR,e.what());
    }

    for (int n=1; n <= ntypes; ++n)
      if (setflag[n] == 0)
        error->one(FLERR,fmt::format("Parameters for atom type {} missing in "
                                     "qeq parameter file", n));
    delete[] setflag;
  }

  MPI_Bcast(chi,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(eta,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(gamma,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(zeta,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(zcore,ntypes+1,MPI_DOUBLE,0,world);
}
