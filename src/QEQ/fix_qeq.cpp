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
   Contributing author: Ray Shan (Sandia)
     Based on fix qeq/reax by H. Metin Aktulga
------------------------------------------------------------------------- */

#include "fix_qeq.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "respa.h"
#include "text_file_reader.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double QSUMSMALL = 0.00001;
static constexpr int MIN_CAP = 50;
static constexpr double SAFE_ZONE = 1.2;
static constexpr bigint MIN_NBRS = 100;

namespace {
  class qeq_parser_error : public std::exception {
    std::string message;
  public:
    explicit qeq_parser_error(const std::string &mesg) { message = mesg; }
    [[nodiscard]] const char *what() const noexcept override { return message.c_str(); }
  };
}

/* ---------------------------------------------------------------------- */

FixQEq::FixQEq(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), list(nullptr), chi(nullptr), eta(nullptr),
  gamma(nullptr), zeta(nullptr), zcore(nullptr), qmin(nullptr), qmax(nullptr), omega(nullptr), chizj(nullptr), shld(nullptr),
  s(nullptr), t(nullptr), s_hist(nullptr), t_hist(nullptr), xls_hist(nullptr), xlt_hist(nullptr),
  Hdia_inv(nullptr), b_s(nullptr), b_t(nullptr), p(nullptr), q(nullptr), r(nullptr), d(nullptr),
  qf(nullptr), q1(nullptr), q2(nullptr), qv(nullptr)
{
  if (narg < 8) utils::missing_cmd_args(FLERR, "fix " + std::string(style), error);

  scalar_flag = 1;
  extscalar = 0;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  cutoff = utils::numeric(FLERR,arg[4],false,lmp);
  tolerance = utils::numeric(FLERR,arg[5],false,lmp);
  maxiter = utils::inumeric(FLERR,arg[6],false,lmp);
  maxwarn = 1;
  matvecs = 0;
  imax = maxiter;

  // extended-Lagrangian charge propagation is selected via the fix style name

  xl.flag = utils::strmatch(style,"/xlmd") ? 1 : 0;

  // check for sane arguments
  if ((nevery <= 0) || (cutoff <= 0.0) || (tolerance <= 0.0) || (maxiter <= 0))
    error->all(FLERR,"Illegal fix qeq command");

  // must have charges

  if (!atom->q_flag) error->all(FLERR, "Fix {} requires atom attribute q", style);

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
  maxexchange = 2*nprev;

  // initialize the per-atom history of atoms created during a run
  create_attribute = 1;

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
  ctip_flag = 0;
  qv = nullptr;

  comm_forward = comm_reverse = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  s_hist = t_hist = nullptr;
  FixQEq::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  for (int i = 0; i < atom->nmax; i++)
    for (int j = 0; j < nprev; ++j)
      s_hist[i][j] = t_hist[i][j] = atom->q[i];

  // QEq parameters from a pair style are only supported by the matching fix styles

  if (strcmp(arg[7],"coul/streitz") == 0) {
    if (!utils::strmatch(style,"^qeq/slater"))
      error->all(FLERR, 7, "Fix {} cannot use QEq parameters from pair style coul/streitz", style);
    streitz_flag = 1;
  } else if (utils::strmatch(arg[7],"^reax..")) {
    if (!utils::strmatch(style,"^qeq/shielded"))
      error->all(FLERR, 7, "Fix {} cannot use QEq parameters from pair style reaxff", style);
    reax_flag = 1;
  } else if (strcmp(arg[7],"coul/ctip") == 0) {
    if (!utils::strmatch(style,"^qeq/ctip"))
      error->all(FLERR, 7, "Fix {} cannot use QEq parameters from pair style coul/ctip", style);
    ctip_flag = 1;
  } else {
    read_file(arg[7]);
  }
}

/* ---------------------------------------------------------------------- */

FixQEq::~FixQEq()
{
  // unregister callbacks to this fix from Atom class
  if (modify->get_fix_by_id(id)) atom->delete_callback(id,Atom::GROW);

  memory->destroy(s_hist);
  memory->destroy(t_hist);
  memory->destroy(xls_hist);
  memory->destroy(xlt_hist);

  deallocate_storage();
  deallocate_matrix();

  memory->destroy(shld);

  if (!streitz_flag && !reax_flag && !ctip_flag) {
    memory->destroy(chi);
    memory->destroy(eta);
    memory->destroy(gamma);
    memory->destroy(zeta);
    memory->destroy(zcore);
    memory->destroy(qmin);
    memory->destroy(qmax);
    memory->destroy(omega);
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
  int i,ii,inum;
  int *ilist, *numneigh;
  bigint m;

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
  auto m_cap_big = (bigint)MAX(m * safezone, mincap * MIN_NBRS);
  if (m_cap_big > MAXSMALLINT)
    error->one(FLERR,"Too many neighbors in fix {}",style);
  m_cap = m_cap_big;

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

void FixQEq::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix {} requires atom attribute q", style);

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix {} group has no atoms", style);

  if ((comm->me == 0) && (!modify->get_fix_by_style("^efield").empty()))
    error->warning(FLERR,"Fix efield is ignored during charge equilibration");

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;

  // compute net charge and print warning if too large

  double qsum_local = 0.0, qsum = 0.0;
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit)
      qsum_local += atom->q[i];
  }
  MPI_Allreduce(&qsum_local,&qsum,1,MPI_DOUBLE,MPI_SUM,world);

  if ((comm->me == 0) && (fabs(qsum) > QSUMSMALL))
    error->warning(FLERR,"Fix {} group is not charge neutral, net charge = {:.8}{}",
                   style, qsum, utils::errorurl(29));
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

  deallocate_storage();
  allocate_storage();

  init_storage();

  deallocate_matrix();
  allocate_matrix();

  // restart the extended-Lagrangian warm-up at the beginning of every run, so
  // that stale or uninitialized auxiliary history (after read_restart, minimize,
  // replicate, delete_atoms, ...) is never used for the charge propagation

  xl.nseed = 0;
  xl.bypass = 1;
  pre_force(vflag);
  xl.bypass = 0;
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

  int *mask = atom->mask;
  for (int i = 0; i < nall; i++) {
    Hdia_inv[i] = 1. / eta[atom->type[i]];
    b_s[i] = -chi[atom->type[i]];
    b_t[i] = -1.0;
    // s is initialized to the current charge so that atoms outside the fix
    // group contribute their (fixed) charge to the electric field felt by the
    // group atoms through the sparse matrix-vector product. The t vector,
    // however, only encodes the charge-neutrality response of the group, so
    // atoms outside the group must contribute 0 to it -- otherwise their fixed
    // charge pollutes the neutralization and the equilibrated charges come out
    // with the wrong sign (see GitHub issue #3543).
    s[i] = atom->q[i];
    t[i] = (mask[i] & groupbit) ? atom->q[i] : 0.0;

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
  comm->reverse_comm(this);

  vector_sum(r , 1.,  b, -1., q, inum);

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      d[i] = r[i] * Hdia_inv[i];
    else d[i] = 0.0;
  }

  b_norm = parallel_norm(b, inum);
  sig_new = parallel_dot(r, d, inum);

  for (loop = 1; loop < imax && sqrt(sig_new)/b_norm > tolerance; ++loop) {
    comm->forward_comm(this);
    sparse_matvec(&H, d, q);
    comm->reverse_comm(this);

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

  // no warning for a deliberately truncated solve in extended-Lagrangian mode

  if ((comm->me == 0) && maxwarn && (imax == maxiter) && (loop >= imax))
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
    else
      b[i] = 0.0;
  }

  for (i = nlocal; i < nall; ++i) {
    b[i] = 0.0;
  }

  for (i = 0; i < nlocal; ++i) {
    // include the rows of atoms outside the group for their fixed charges
    for (itr_j=A->firstnbr[i]; itr_j<A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
      j = A->jlist[itr_j];
      b[i] += A->val[itr_j] * x[j];
      b[j] += A->val[itr_j] * x[i];
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
  comm->forward_comm(this); //Dist_vector(atom->q);
}

/* ----------------------------------------------------------------------
   parse keywords shared between qeq fix styles starting at position iarg.
   returns the number of arguments consumed, 0 if the keyword is unknown
------------------------------------------------------------------------- */

int FixQEq::parse_common_keyword(int narg, char **arg, int iarg)
{
  if (strcmp(arg[iarg],"warn") == 0) {
    if (iarg+2 > narg)
      utils::missing_cmd_args(FLERR, std::string("fix ") + style + " warn", error);
    maxwarn = utils::logical(FLERR,arg[iarg+1],false,lmp);
    return 2;
  }

  // keywords below apply to the extended-Lagrangian (xlmd) fix styles only

  if (!xl.flag) return 0;
  return xl.parse_keyword(lmp, style, narg, arg, iarg);
}

/* ----------------------------------------------------------------------
   validate extended-Lagrangian settings and set up the auxiliary storage.
   must be called at the end of the constructor of styles supporting xlmd
------------------------------------------------------------------------- */

void FixQEq::finalize_xl()
{
  if (!xl.flag) return;

  xl.finalize(lmp, style, nevery);
  maxexchange += 2*xl.nhist;

  // allocate and initialize the auxiliary history arrays. they stay in sync
  // with the atoms through the grow/copy/exchange callbacks that were already
  // registered in the base class constructor

  FixQEq::grow_arrays(atom->nmax);
  for (int i = 0; i < atom->nmax; i++) {
    xl.fill(xls_hist[i], atom->q[i]);
    xl.fill(xlt_hist[i], 0.0);
  }
}

/* ----------------------------------------------------------------------
   shared solver driver for the CG based qeq styles: solve H s = b_s and
   H t = b_t either to convergence or, in extended-Lagrangian mode, with
   a fixed small number of CG iterations seeded from the auxiliary variables
------------------------------------------------------------------------- */

int FixQEq::solve_st()
{
  const bool active = xl.ready(lmp);

  if (active) xl_predict();
  imax = active ? xl.ncg + 1 : maxiter;
  int nsolve = CG(b_s,s);
  nsolve += CG(b_t,t);
  imax = maxiter;
  if (xl.flag) xl_update(active);
  return nsolve/2;
}

/* ----------------------------------------------------------------------
   use the auxiliary variables as the initial guess for the truncated solve
------------------------------------------------------------------------- */

void FixQEq::xl_predict()
{
  int inum = list->inum;
  int *ilist = list->ilist;

  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      s[i] = xls_hist[i][0];
      t[i] = xlt_hist[i][0];
    }
  }

  pack_flag = 2;
  comm->forward_comm(this);
  pack_flag = 3;
  comm->forward_comm(this);
}

/* ----------------------------------------------------------------------
   advance the auxiliary variables theta_s and theta_t. when active, use
   the time-reversible (dissipative) Verlet propagation following
   Nomura et al., CPC 192, 91 (2015) for the plain Verlet case and
   Niklasson et al., JCP 130, 214109 (2009) for the dissipation term.
   otherwise (re-)seed the history from the fully converged solution
------------------------------------------------------------------------- */

void FixQEq::xl_update(bool active)
{
  int inum = list->inum;
  int *ilist = list->ilist;
  const bool restart = !active && xl.restart(lmp);

  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      if (active) {
        const double ths = xl.propagate(s[i], xls_hist[i]);
        const double tht = xl.propagate(t[i], xlt_hist[i]);
        xl.push(xls_hist[i], ths);
        xl.push(xlt_hist[i], tht);
      } else if (restart) {
        xl.fill(xls_hist[i], s[i]);
        xl.fill(xlt_hist[i], t[i]);
      } else {
        xl.push(xls_hist[i], s[i]);
        xl.push(xlt_hist[i], t[i]);
      }
    }
  }
  xl.advance(lmp, active, restart);
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
  bytes += (double)atom->nmax*xl.nhist*2 * sizeof(double); // xls_hist & xlt_hist
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
  if (xl.nhist) {
    memory->grow(xls_hist,nmax,xl.nhist,"qeq:xls_hist");
    memory->grow(xlt_hist,nmax,xl.nhist,"qeq:xlt_hist");
  }
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
  for (int m = 0; m < xl.nhist; m++) {
    xls_hist[j][m] = xls_hist[i][m];
    xlt_hist[j][m] = xlt_hist[i][m];
  }
}

/* ----------------------------------------------------------------------
   initialize history for an atom created during the run. grow_arrays()
   leaves the new storage uninitialized and would feed garbage into the
   initial guess extrapolation and the auxiliary variable propagation
------------------------------------------------------------------------- */

void FixQEq::set_arrays(int i)
{
  for (int m = 0; m < nprev; m++) s_hist[i][m] = t_hist[i][m] = 0.0;
  for (int m = 0; m < xl.nhist; m++) xls_hist[i][m] = xlt_hist[i][m] = 0.0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixQEq::pack_exchange(int i, double *buf)
{
  int n = 0;
  for (int m = 0; m < nprev; m++) buf[n++] = s_hist[i][m];
  for (int m = 0; m < nprev; m++) buf[n++] = t_hist[i][m];
  for (int m = 0; m < xl.nhist; m++) buf[n++] = xls_hist[i][m];
  for (int m = 0; m < xl.nhist; m++) buf[n++] = xlt_hist[i][m];
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixQEq::unpack_exchange(int i, double *buf)
{
  int n = 0;
  for (int m = 0; m < nprev; m++) s_hist[i][m] = buf[n++];
  for (int m = 0; m < nprev; m++) t_hist[i][m] = buf[n++];
  for (int m = 0; m < xl.nhist; m++) xls_hist[i][m] = buf[n++];
  for (int m = 0; m < xl.nhist; m++) xlt_hist[i][m] = buf[n++];
  return n;
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
  // fix qeq/ctip reads three more columns per atom type
  const bool ctip_format = utils::strmatch(style, "^qeq/ctip");

  const int ntypes = atom->ntypes;

  memory->create(chi,ntypes+1,"qeq:chi");
  memory->create(eta,ntypes+1,"qeq:eta");
  memory->create(gamma,ntypes+1,"qeq:gamma");
  memory->create(zeta,ntypes+1,"qeq:zeta");
  memory->create(zcore,ntypes+1,"qeq:zcore");
  memory->create(qmin,ntypes+1,"qeq:qmin");
  memory->create(qmax,ntypes+1,"qeq:qmax");
  memory->create(omega,ntypes+1,"qeq:omega");

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  if (comm->me == 0) {
    int *setflag = new int[ntypes+1];
    for (int n=0; n <= ntypes; ++n) {
      setflag[n] = 0;
      chi[n] = eta[n] = gamma[n] = zeta[n] = zcore[n] = qmin[n] = qmax[n] = omega[n] = 0.0;
    }

    FILE *fp = nullptr;
    try {
      int nlo,nhi;
      double val;

      fp = utils::open_potential(file,lmp,nullptr);
      if (fp == nullptr)
        throw qeq_parser_error(fmt::format("Cannot open fix qeq parameter file {}: {}",
                                           file,utils::getsyserror()));
      TextFileReader reader(fp, "qeq parameter");

      while (true) {
        auto values = reader.next_values(0);

        if (values.count() == 0) continue;
        if (ctip_format) {
          if (values.count() < 9)
            throw qeq_parser_error(fmt::format("Invalid qeq parameter file for {}", style));
        } else {
          if (values.count() < 6)
            throw qeq_parser_error(fmt::format("Invalid qeq parameter file for {}", style));
        }

        auto word = values.next_string();
        utils::bounds(FLERR,word,1,ntypes,nlo,nhi,nullptr);
        if ((nlo < 0) || (nhi < 0))
          throw qeq_parser_error(fmt::format("Invalid atom type range: {}",word));

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
        if (ctip_format) {
          val = values.next_double();
          for (int n=nlo; n <= nhi; ++n) qmin[n] = val;
          val = values.next_double();
          for (int n=nlo; n <= nhi; ++n) qmax[n] = val;
          val = values.next_double();
          for (int n=nlo; n <= nhi; ++n) omega[n] = val;
        }
        for (int n=nlo; n <= nhi; ++n) setflag[n] = 1;
      }
    } catch (EOFException &) {
      fclose(fp);
    } catch (std::exception &e) {
      error->one(FLERR,e.what());
    }

    for (int n=1; n <= ntypes; ++n)
      if (setflag[n] == 0)
        error->one(FLERR, "Parameters for atom type {} missing in qeq parameter file", n);
    delete[] setflag;
  }

  MPI_Bcast(chi,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(eta,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(gamma,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(zeta,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(zcore,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(qmin,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(qmax,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(omega,ntypes+1,MPI_DOUBLE,0,world);
}
