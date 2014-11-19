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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
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
  Fix(lmp, narg, arg)
{
  if (narg < 8) error->all(FLERR,"Illegal fix qeq command");
  
  nevery = force->inumeric(FLERR,arg[3]);
  cutoff = force->numeric(FLERR,arg[4]);
  tolerance = force->numeric(FLERR,arg[5]);
  maxiter = force->inumeric(FLERR,arg[6]);

  alpha = 0.20;
  swa = 0.0;
  swb = cutoff;

  shld = NULL;

  n = n_cap = 0;
  N = nmax = 0;
  m_fill = m_cap = 0;
  pack_flag = 0;
  s = NULL;
  t = NULL;
  nprev = 5;

  Hdia_inv = NULL;
  b_s = NULL;
  b_t = NULL;
  b_prc = NULL;
  b_prm = NULL;

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

  comm_forward = comm_reverse = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  s_hist = t_hist = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  for( int i = 0; i < atom->nmax; i++ )
    for (int j = 0; j < nprev; ++j )
      s_hist[i][j] = t_hist[i][j] = 0;

  read_file(arg[7]);

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

  memory->destroy(chi);
  memory->destroy(eta);
  memory->destroy(gamma);
  memory->destroy(zeta);
  memory->destroy(zcore);
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
    memory->create(b_prc,nmax,"qeq:b_prc");
    memory->create(b_prm,nmax,"qeq:b_prm");

    memory->create(p,nmax,"qeq:p");
    memory->create(q,nmax,"qeq:q");
    memory->create(r,nmax,"qeq:r");
    memory->create(d,nmax,"qeq:d");

    memory->create(chizj,nmax,"qeq:chizj");
    memory->create(qf,nmax,"qeq:qf");
    memory->create(q1,nmax,"qeq:q1");
    memory->create(q2,nmax,"qeq:q2");
}

/* ---------------------------------------------------------------------- */

void FixQEq::deallocate_storage()
{
    memory->destroy(s);
    memory->destroy(t);

    memory->destroy( Hdia_inv );
    memory->destroy( b_s );
    memory->destroy( b_t );
    memory->destroy( b_prc );
    memory->destroy( b_prm );

    memory->destroy( p );
    memory->destroy( q );
    memory->destroy( r );
    memory->destroy( d );

    memory->destroy( chizj );
    memory->destroy( qf );
    memory->destroy( q1 );
    memory->destroy( q2 );
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

  n = atom->nlocal;
  n_cap = MAX( (int)(n * safezone), mincap );
  N = atom->nlocal + atom->nghost;

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

void FixQEq::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixQEq::setup_pre_force(int vflag)
{
  neighbor->build_one(list);

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

void FixQEq::min_setup_pre_force(int vflag)
{
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEq::init_storage()
{
  int N = atom->nlocal + atom->nghost;

  for( int i = 0; i < N; i++ ) {
    Hdia_inv[i] = 1. / eta[atom->type[i]];
    b_s[i] = -chi[atom->type[i]];
    b_t[i] = -1.0;
    b_prc[i] = 0;
    b_prm[i] = 0;
    s[i] = t[i] = 0;

    chizj[i] = 0.0;

    qf[i] = 0.0;
    q1[i] = 0.0;
    q2[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixQEq::pre_force_respa(int vflag, int ilevel, int iloop)
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
  int  i, j;
  double tmp, alfa, beta, b_norm;
  double sig_old, sig_new;

  int nn, jj;
  int *ilist;

  nn = list->inum;
  ilist = list->ilist;

  pack_flag = 1;
  sparse_matvec( &H, x, q );
  comm->reverse_comm_fix( this ); //Coll_Vector( q );

  vector_sum( r , 1.,  b, -1., q, nn );

  for( jj = 0; jj < nn; ++jj ) {
    j = ilist[jj];
    if (atom->mask[j] & groupbit)
      d[j] = r[j] * Hdia_inv[j]; //pre-condition
  }

  b_norm = parallel_norm( b, nn );
  sig_new = parallel_dot( r, d, nn);

  for( i = 1; i < maxiter && sqrt(sig_new) / b_norm > tolerance; ++i ) {
    comm->forward_comm_fix(this); //Dist_vector( d );
    sparse_matvec( &H, d, q );
    comm->reverse_comm_fix(this); //Coll_vector( q );

    tmp = parallel_dot( d, q, nn);
    alfa = sig_new / tmp;

    vector_add( x, alfa, d, nn );
    vector_add( r, -alfa, q, nn );

    // pre-conditioning
    for( jj = 0; jj < nn; ++jj ) {
      j = ilist[jj];
      if (atom->mask[j] & groupbit)
        p[j] = r[j] * Hdia_inv[j];
    }

    sig_old = sig_new;
    sig_new = parallel_dot( r, p, nn);

    beta = sig_new / sig_old;
    vector_sum( d, 1., p, beta, d, nn );

  }

  if (i >= maxiter && comm->me == 0) {
    char str[128];
    sprintf(str,"Fix qeq CG convergence failed (%g) after %d iterations "
            "at " BIGINT_FORMAT " step",sqrt(sig_new) / b_norm,i,update->ntimestep);
    error->warning(FLERR,str);
  }

  return i;
}


/* ---------------------------------------------------------------------- */

void FixQEq::sparse_matvec( sparse_matrix *A, double *x, double *b )
{
  int i, j, itr_j;
  int nn, NN;
  int *ilist;

  nn = atom->nlocal;
  NN = atom->nlocal + atom->nghost;
  ilist = list->ilist;

  for( i = 0; i < nn; ++i ) {
    if (atom->mask[i] & groupbit)
      b[i] = eta[ atom->type[i] ] * x[i];
  }

  for( i = nn; i < NN; ++i ) {
    if (atom->mask[i] & groupbit)
      b[i] = 0;
  }

  for( i = 0; i < nn; ++i ) {
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
  int i, k;
  double u, s_sum, t_sum;
  double *q = atom->q;

  int nn, ii;
  int *ilist;

  nn = list->inum;
  ilist = list->ilist;

  s_sum = parallel_vector_acc( s, nn );
  t_sum = parallel_vector_acc( t, nn);
  u = s_sum / t_sum;

  for( ii = 0; ii < nn; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      q[i] = s[i] - u * t[i];

      /* backup s & t */
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
                          int pbc_flag, int *pbc)
{
  int m;

  if( pack_flag == 1)
    for(m = 0; m < n; m++) buf[m] = d[list[m]];
  else if( pack_flag == 2 )
    for(m = 0; m < n; m++) buf[m] = s[list[m]];
  else if( pack_flag == 3 )
    for(m = 0; m < n; m++) buf[m] = t[list[m]];
  else if( pack_flag == 4 )
    for(m = 0; m < n; m++) buf[m] = atom->q[list[m]];

  return n;
}

/* ---------------------------------------------------------------------- */

void FixQEq::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if( pack_flag == 1)
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
  return n;
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

void FixQEq::copy_arrays(int i, int j, int delflag)
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

int FixQEq::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nprev; m++) s_hist[nlocal][m] = buf[m];
  for (int m = 0; m < nprev; m++) t_hist[nlocal][m] = buf[nprev+m];
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
/* ---------------------------------------------------------------------- */

void FixQEq::read_file(char *file)
{
  int i,itype,ntypes;
  int params_per_line = 6;
  char **words = new char*[params_per_line+1];

  ntypes = atom->ntypes;

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
      sprintf(str,"Cannot open Tersoff potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int n,nwords,ielement,eof;
  char line[MAXLINE],*ptr;

  eof = ielement = 0;

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

    ielement ++;
    if (ielement > ntypes)
      error->all(FLERR,"Invalid fix qeq parameter file");

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    itype = atoi(words[0]);
    chi[itype]   = atof(words[1]);
    eta[itype]   = atof(words[2]);
    gamma[itype] = atof(words[3]);
    zeta[itype]  = atof(words[4]);
    zcore[itype] = atof(words[5]);
  }
  delete [] words;
}
