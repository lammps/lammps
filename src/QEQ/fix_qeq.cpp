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
   Contributing author: Ray Shan (Sandia)
     Based on fix qeq/reax by H. Metin Aktulga
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "fix_qeq.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "kspace.h"
#include "group.h"
#include "pair.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

FixQEq::FixQEq(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), list(NULL), chi(NULL), eta(NULL),
  gamma(NULL), zeta(NULL), zcore(NULL), chizj(NULL), shld(NULL),
  s(NULL), t(NULL), s_hist(NULL), t_hist(NULL), Hdia_inv(NULL), b_s(NULL),
  b_t(NULL), p(NULL), q(NULL), r(NULL), d(NULL),
  qf(NULL), q1(NULL), q2(NULL), qv(NULL)
{
  if (narg < 8) error->all(FLERR,"Illegal fix qeq command");

  nevery = force->inumeric(FLERR,arg[3]);
  cutoff = force->numeric(FLERR,arg[4]);
  tolerance = force->numeric(FLERR,arg[5]);
  maxiter = force->inumeric(FLERR,arg[6]);

  // check for sane arguments
  if ((nevery <= 0) || (cutoff <= 0.0) || (tolerance <= 0.0) || (maxiter <= 0))
    error->all(FLERR,"Illegal fix qeq command");

  alpha = 0.20;
  swa = 0.0;
  swb = cutoff;

  shld = NULL;

  nlocal = n_cap = 0;
  nall = nmax = 0;
  m_fill = m_cap = 0;
  pack_flag = 0;
  s = NULL;
  t = NULL;
  nprev = 5;

  Hdia_inv = NULL;
  b_s = NULL;
  b_t = NULL;

  // CG
  p = NULL;
  q = NULL;
  r = NULL;
  d = NULL;

  // H matrix
  H.firstnbr = NULL;
  H.numnbrs = NULL;
  H.jlist = NULL;
  H.val = NULL;

  // others
  cutoff_sq = cutoff*cutoff;
  chizj = NULL;
  qf = NULL;
  q1 = NULL;
  q2 = NULL;
  streitz_flag = 0;
  reax_flag = 0;
  qv = NULL;

  comm_forward = comm_reverse = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  s_hist = t_hist = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  for( int i = 0; i < atom->nmax; i++ )
    for (int j = 0; j < nprev; ++j )
      s_hist[i][j] = t_hist[i][j] = atom->q[i];

  if (strcmp(arg[7],"coul/streitz") == 0) {
    streitz_flag = 1;
  } else if (strcmp(arg[7],"reax/c") == 0) {
    reax_flag = 1;
  } else {
    read_file(arg[7]);
  }

}

/* ---------------------------------------------------------------------- */

FixQEq::~FixQEq()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);

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

  memory->destroy( Hdia_inv );
  memory->destroy( b_s );
  memory->destroy( b_t );

  memory->destroy( p );
  memory->destroy( q );
  memory->destroy( r );
  memory->destroy( d );

  memory->destroy( chizj );
  memory->destroy( qf );
  memory->destroy( q1 );
  memory->destroy( q2 );

  memory->destroy( qv );
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
  n_cap = MAX( (int)(nlocal * safezone), mincap );
  nall = atom->nlocal + atom->nghost;

  // determine the total space for the H matrix

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;

  m = 0;
  for( ii = 0; ii < inum; ii++ ) {
    i = ilist[ii];
    m += numneigh[i];
  }
  m_cap = MAX( (int)(m * safezone), mincap * MIN_NBRS );

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
  memory->destroy( H.firstnbr );
  memory->destroy( H.numnbrs );
  memory->destroy( H.jlist );
  memory->destroy( H.val );
}

/* ---------------------------------------------------------------------- */

void FixQEq::reallocate_matrix()
{
  deallocate_matrix();
  allocate_matrix();
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

  for( int i = 0; i < nall; i++ ) {
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

int FixQEq::CG( double *b, double *x )
{
  int  loop, i, ii, inum, *ilist;
  double tmp, alfa, beta, b_norm;
  double sig_old, sig_new;

  inum = list->inum;
  ilist = list->ilist;

  pack_flag = 1;
  sparse_matvec( &H, x, q );
  comm->reverse_comm_fix( this );

  vector_sum( r , 1.,  b, -1., q, inum );

  for( ii = 0; ii < inum; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      d[i] = r[i] * Hdia_inv[i];
    else d[i] = 0.0;
  }

  b_norm = parallel_norm( b, inum );
  sig_new = parallel_dot( r, d, inum);

  for( loop = 1; loop < maxiter && sqrt(sig_new)/b_norm > tolerance; ++loop ) {
    comm->forward_comm_fix(this);
    sparse_matvec( &H, d, q );
    comm->reverse_comm_fix(this);

    tmp = parallel_dot( d, q, inum);
    alfa = sig_new / tmp;

    vector_add( x, alfa, d, inum );
    vector_add( r, -alfa, q, inum );

    for( ii = 0; ii < inum; ++ii ) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit)
        p[i] = r[i] * Hdia_inv[i];
    }

    sig_old = sig_new;
    sig_new = parallel_dot( r, p, inum);

    beta = sig_new / sig_old;
    vector_sum( d, 1., p, beta, d, inum );
  }

  if (loop >= maxiter && comm->me == 0) {
    char str[128];
    sprintf(str,"Fix qeq CG convergence failed (%g) after %d iterations "
            "at " BIGINT_FORMAT " step",sqrt(sig_new)/b_norm,loop,update->ntimestep);
    error->warning(FLERR,str);
  }

  return loop;
}


/* ---------------------------------------------------------------------- */

void FixQEq::sparse_matvec( sparse_matrix *A, double *x, double *b )
{
  int i, j, itr_j;

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;

  for( i = 0; i < nlocal; ++i ) {
    if (atom->mask[i] & groupbit)
      b[i] = eta[ atom->type[i] ] * x[i];
  }

  for( i = nlocal; i < nall; ++i ) {
    if (atom->mask[i] & groupbit)
      b[i] = 0;
  }

  for( i = 0; i < nlocal; ++i ) {
    if (atom->mask[i] & groupbit) {
      for( itr_j=A->firstnbr[i]; itr_j<A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
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

  s_sum = parallel_vector_acc( s, inum );
  t_sum = parallel_vector_acc( t, inum);
  u = s_sum / t_sum;

  for( ii = 0; ii < inum; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      q[i] = s[i] - u * t[i];

      for( k = 4; k > 0; --k ) {
        s_hist[i][k] = s_hist[i][k-1];
        t_hist[i][k] = t_hist[i][k-1];
      }
      s_hist[i][0] = s[i];
      t_hist[i][0] = t[i];
    }
  }

  pack_flag = 4;
  comm->forward_comm_fix( this ); //Dist_vector( atom->q );
}

/* ---------------------------------------------------------------------- */

int FixQEq::pack_forward_comm(int n, int *list, double *buf,
                          int /*pbc_flag*/, int * /*pbc*/)
{
  int m;

  if (pack_flag == 1)
    for(m = 0; m < n; m++) buf[m] = d[list[m]];
  else if( pack_flag == 2 )
    for(m = 0; m < n; m++) buf[m] = s[list[m]];
  else if( pack_flag == 3 )
    for(m = 0; m < n; m++) buf[m] = t[list[m]];
  else if( pack_flag == 4 )
    for(m = 0; m < n; m++) buf[m] = atom->q[list[m]];

  return m;
}

/* ---------------------------------------------------------------------- */

void FixQEq::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if (pack_flag == 1)
    for(m = 0, i = first; m < n; m++, i++) d[i] = buf[m];
  else if( pack_flag == 2)
    for(m = 0, i = first; m < n; m++, i++) s[i] = buf[m];
  else if( pack_flag == 3)
    for(m = 0, i = first; m < n; m++, i++) t[i] = buf[m];
  else if( pack_flag == 4)
    for(m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
}

/* ---------------------------------------------------------------------- */

int FixQEq::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  for(m = 0, i = first; m < n; m++, i++) buf[m] = q[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void FixQEq::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m;

  for(m = 0; m < n; m++) q[list[m]] += buf[m];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixQEq::memory_usage()
{
  double bytes;

  bytes = atom->nmax*nprev*2 * sizeof(double); // s_hist & t_hist
  bytes += atom->nmax*11 * sizeof(double); // storage
  bytes += n_cap*2 * sizeof(int); // matrix...
  bytes += m_cap * sizeof(int);
  bytes += m_cap * sizeof(double);

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

double FixQEq::parallel_norm( double *v, int n )
{
  int  i;
  double my_sum, norm_sqr;

  int ii;
  int *ilist;

  ilist = list->ilist;

  my_sum = 0.0;
  norm_sqr = 0.0;
  for( ii = 0; ii < n; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_sum += v[i]*v[i];
  }

  MPI_Allreduce( &my_sum, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, world );

  return sqrt( norm_sqr );
}

/* ---------------------------------------------------------------------- */

double FixQEq::parallel_dot( double *v1, double *v2, int n)
{
  int  i;
  double my_dot, res;

  int ii;
  int *ilist;

  ilist = list->ilist;

  my_dot = 0.0;
  res = 0.0;
  for( ii = 0; ii < n; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_dot += v1[i] * v2[i];
  }

  MPI_Allreduce( &my_dot, &res, 1, MPI_DOUBLE, MPI_SUM, world );

  return res;
}

/* ---------------------------------------------------------------------- */

double FixQEq::parallel_vector_acc( double *v, int n )
{
  int  i;
  double my_acc, res;

  int ii;
  int *ilist;

  ilist = list->ilist;

  my_acc = 0.0;
  res = 0.0;
  for( ii = 0; ii < n; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_acc += v[i];
  }

  MPI_Allreduce( &my_acc, &res, 1, MPI_DOUBLE, MPI_SUM, world );

  return res;
}

/* ---------------------------------------------------------------------- */

void FixQEq::vector_sum( double* dest, double c, double* v,
                                double d, double* y, int k )
{
  int kk;
  int *ilist;

  ilist = list->ilist;

  for( --k; k>=0; --k ) {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] = c * v[kk] + d * y[kk];
  }
}

/* ---------------------------------------------------------------------- */

void FixQEq::vector_add( double* dest, double c, double* v, int k )
{
  int kk;
  int *ilist;

  ilist = list->ilist;

  for( --k; k>=0; --k ) {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] += c * v[kk];
  }
}

/* ---------------------------------------------------------------------- */

void FixQEq::read_file(char *file)
{
  int i;
  int params_per_line = 6;
  char **words = new char*[params_per_line+1];

  int ntypes = atom->ntypes;
  int *setflag = new int[ntypes+1];
  for (i=0; i <= ntypes; ++i) setflag[i] = 0;

  memory->create(chi,ntypes+1,"qeq:chi");
  memory->create(eta,ntypes+1,"qeq:eta");
  memory->create(gamma,ntypes+1,"qeq:gamma");
  memory->create(zeta,ntypes+1,"qeq:zeta");
  memory->create(zcore,ntypes+1,"qeq:zcore");

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open fix qeq parameter file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int n,nwords,eof,nlo,nhi;
  char line[MAXLINE],*ptr;

  eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // must have 6 parameters per line.

    if (nwords < 6)
      error->all(FLERR,"Invalid fix qeq parameter file");

    // words = ptrs to first 6 words in line

    for (n=0, words[n] = strtok(line," \t\n\r\f");
         n < 6;
         words[++n] = strtok(NULL," \t\n\r\f"));

    force->bounds(FLERR,words[0],ntypes,nlo,nhi);
    for (n=nlo; n <=nhi; ++n) {
      chi[n]     = force->numeric(FLERR,words[1]);
      eta[n]     = force->numeric(FLERR,words[2]);
      gamma[n]   = force->numeric(FLERR,words[3]);
      zeta[n]    = force->numeric(FLERR,words[4]);
      zcore[n]   = force->numeric(FLERR,words[5]);
      setflag[n] = 1;
    }
  }

  // check if all types are set
  for (n=1; n <= ntypes; ++n)
    if (setflag[n] == 0)
      error->all(FLERR,"Invalid fix qeq parameter file");

  delete [] words;
  delete [] setflag;
}
