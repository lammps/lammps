/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "pair_reax_c.h"
#if defined(PURE_REAX)
#include "basic_comm.h"
#include "vector.h"
#elif defined(LAMMPS_REAX)
#include "reaxc_basic_comm.h"
#include "reaxc_vector.h"
#endif

#if defined(PURE_REAX)
void real_packer( void *dummy, mpi_out_data *out_buf )
{
  int i;
  real *buf = (real*) dummy;
  real *out = (real*) out_buf->out_atoms;

  for( i = 0; i < out_buf->cnt; ++i )
    out[i] = buf[ out_buf->index[i] ];
}


void rvec_packer( void *dummy, mpi_out_data *out_buf )
{
  int i;
  rvec *buf = (rvec*) dummy;
  rvec *out = (rvec*)out_buf->out_atoms;

  for( i = 0; i < out_buf->cnt; ++i )
    memcpy( out[i], buf[ out_buf->index[i] ], sizeof(rvec) );
}


void rvec2_packer( void *dummy, mpi_out_data *out_buf )
{
  int i;
  rvec2 *buf = (rvec2*) dummy;
  rvec2 *out = (rvec2*) out_buf->out_atoms;

  for( i = 0; i < out_buf->cnt; ++i )
    memcpy( out[i], buf[ out_buf->index[i] ], sizeof(rvec2) );
}


void Dist( reax_system* system, mpi_datatypes *mpi_data,
           void *buf, MPI_Datatype type, int scale, dist_packer pack )
{
  int d;
  mpi_out_data *out_bufs;
  MPI_Comm comm;
  MPI_Request req1, req2;
  MPI_Status stat1, stat2;
  neighbor_proc *nbr1, *nbr2;

#if defined(DEBUG)
  fprintf( stderr, "p%d dist: entered\n", system->my_rank );
#endif
  comm = mpi_data->comm_mesh3D;
  out_bufs = mpi_data->out_buffers;

  for( d = 0; d < 3; ++d ) {
    /* initiate recvs */
    nbr1 = &(system->my_nbrs[2*d]);
    if( nbr1->atoms_cnt )
      MPI_Irecv( buf + nbr1->atoms_str*scale, nbr1->atoms_cnt, type,
                 nbr1->rank, 2*d+1,comm, &req1 );

    nbr2 = &(system->my_nbrs[2*d+1]);
    if( nbr2->atoms_cnt )
      MPI_Irecv( buf + nbr2->atoms_str*scale, nbr2->atoms_cnt, type,
                 nbr2->rank, 2*d, comm, &req2 );

    /* send both messages in dimension d */
    if( out_bufs[2*d].cnt ) {
      pack( buf, out_bufs + (2*d) );
      MPI_Send( out_bufs[2*d].out_atoms, out_bufs[2*d].cnt, type,
                nbr1->rank, 2*d, comm );
    }

    if( out_bufs[2*d+1].cnt ) {
      pack( buf, out_bufs + (2*d+1) );
      MPI_Send( out_bufs[2*d+1].out_atoms, out_bufs[2*d+1].cnt, type,
                nbr2->rank, 2*d+1, comm );
    }

    if( nbr1->atoms_cnt ) MPI_Wait( &req1, &stat1 );
    if( nbr2->atoms_cnt ) MPI_Wait( &req2, &stat2 );
  }

#if defined(DEBUG)
  fprintf( stderr, "p%d dist: done\n", system->my_rank );
#endif
}


void real_unpacker( void *dummy_in, void *dummy_buf, mpi_out_data *out_buf )
{
  int i;
  real *in = (real*) dummy_in;
  real *buf = (real*) dummy_buf;

  for( i = 0; i < out_buf->cnt; ++i )
    buf[ out_buf->index[i] ] += in[i];
}


void rvec_unpacker( void *dummy_in, void *dummy_buf, mpi_out_data *out_buf )
{
  int i;
  rvec *in = (rvec*) dummy_in;
  rvec *buf = (rvec*) dummy_buf;

  for( i = 0; i < out_buf->cnt; ++i ) {
    rvec_Add( buf[ out_buf->index[i] ], in[i] );
#if defined(DEBUG)
    fprintf( stderr, "rvec_unpacker: cnt=%d  i =%d  index[i]=%d\n",
             out_buf->cnt, i, out_buf->index[i] );
#endif
  }
}


void rvec2_unpacker( void *dummy_in, void *dummy_buf, mpi_out_data *out_buf )
{
  int i;
  rvec2 *in = (rvec2*) dummy_in;
  rvec2 *buf = (rvec2*) dummy_buf;

  for( i = 0; i < out_buf->cnt; ++i ) {
    buf[ out_buf->index[i] ][0] += in[i][0];
    buf[ out_buf->index[i] ][1] += in[i][1];
  }
}


void Coll( reax_system* system, mpi_datatypes *mpi_data,
           void *buf, MPI_Datatype type, int scale, coll_unpacker unpack )
{
  int d;
  void *in1, *in2;
  mpi_out_data *out_bufs;
  MPI_Comm comm;
  MPI_Request req1, req2;
  MPI_Status stat1, stat2;
  neighbor_proc *nbr1, *nbr2;

#if defined(DEBUG)
  fprintf( stderr, "p%d coll: entered\n", system->my_rank );
#endif
  comm = mpi_data->comm_mesh3D;
  in1 = mpi_data->in1_buffer;
  in2 = mpi_data->in2_buffer;
  out_bufs = mpi_data->out_buffers;

  for( d = 2; d >= 0; --d ) {
    /* initiate recvs */
    nbr1 = &(system->my_nbrs[2*d]);
    if( out_bufs[2*d].cnt )
      MPI_Irecv(in1, out_bufs[2*d].cnt, type, nbr1->rank, 2*d+1, comm, &req1);

    nbr2 = &(system->my_nbrs[2*d+1]);
    if( out_bufs[2*d+1].cnt )
      MPI_Irecv(in2, out_bufs[2*d+1].cnt, type, nbr2->rank, 2*d, comm, &req2);

    /* send both messages in dimension d */
    if( nbr1->atoms_cnt )
      MPI_Send( buf + nbr1->atoms_str*scale, nbr1->atoms_cnt, type,
                nbr1->rank, 2*d, comm );

    if( nbr2->atoms_cnt )
      MPI_Send( buf + nbr2->atoms_str*scale, nbr2->atoms_cnt, type,
                nbr2->rank, 2*d+1, comm );

#if defined(DEBUG)
    fprintf( stderr, "p%d coll[%d] nbr1: str=%d cnt=%d recv=%d\n",
             system->my_rank, d, nbr1->atoms_str, nbr1->atoms_cnt,
             out_bufs[2*d].cnt );
    fprintf( stderr, "p%d coll[%d] nbr2: str=%d cnt=%d recv=%d\n",
             system->my_rank, d, nbr2->atoms_str, nbr2->atoms_cnt,
             out_bufs[2*d+1].cnt );
#endif

    if( out_bufs[2*d].cnt ) {
      MPI_Wait( &req1, &stat1 );
      unpack( in1, buf, out_bufs + (2*d) );
    }

    if( out_bufs[2*d+1].cnt ) {
      MPI_Wait( &req2, &stat2 );
      unpack( in2, buf, out_bufs + (2*d+1) );
    }
  }

#if defined(DEBUG)
  fprintf( stderr, "p%d coll: done\n", system->my_rank );
#endif
}
#endif /*PURE_REAX*/

/*****************************************************************************/
real Parallel_Norm( real *v, int n, MPI_Comm comm )
{
  int  i;
  real my_sum, norm_sqr;

  my_sum = 0;
  for( i = 0; i < n; ++i )
    my_sum += SQR( v[i] );

  MPI_Allreduce( &my_sum, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, comm );

  return sqrt( norm_sqr );
}



real Parallel_Dot( real *v1, real *v2, int n, MPI_Comm comm )
{
  int  i;
  real my_dot, res;

  my_dot = 0;
  for( i = 0; i < n; ++i )
    my_dot += v1[i] * v2[i];

  MPI_Allreduce( &my_dot, &res, 1, MPI_DOUBLE, MPI_SUM, comm );

  return res;
}



real Parallel_Vector_Acc( real *v, int n, MPI_Comm comm )
{
  int  i;
  real my_acc, res;

  my_acc = 0;
  for( i = 0; i < n; ++i )
    my_acc += v[i];

  MPI_Allreduce( &my_acc, &res, 1, MPI_DOUBLE, MPI_SUM, comm );

  return res;
}


/*****************************************************************************/
#if defined(TEST_FORCES)
void Coll_ids_at_Master( reax_system *system, storage *workspace,
                         mpi_datatypes *mpi_data )
{
  int i;
  rc_tagint *id_list;

  MPI_Gather( &system->n, 1, MPI_INT, workspace->rcounts, 1, MPI_INT,
              MASTER_NODE, mpi_data->world );

  if( system->my_rank == MASTER_NODE ){
    workspace->displs[0] = 0;
    for( i = 1; i < system->wsize; ++i )
      workspace->displs[i] = workspace->displs[i-1] + workspace->rcounts[i-1];
  }

  id_list = (int*) malloc( system->n * sizeof(int) );
  for( i = 0; i < system->n; ++i )
    id_list[i] = system->my_atoms[i].orig_id;

  MPI_Gatherv( id_list, system->n, MPI_INT,
               workspace->id_all, workspace->rcounts, workspace->displs,
               MPI_INT, MASTER_NODE, mpi_data->world );

  free( id_list );

#if defined(DEBUG)
  if( system->my_rank == MASTER_NODE ) {
    for( i = 0 ; i < system->bigN; ++i )
      fprintf( stderr, "id_all[%d]: %d\n", i, workspace->id_all[i] );
  }
#endif
}


void Coll_rvecs_at_Master( reax_system *system, storage *workspace,
                           mpi_datatypes *mpi_data, rvec* v )
{
  MPI_Gatherv( v, system->n, mpi_data->mpi_rvec,
               workspace->f_all, workspace->rcounts, workspace->displs,
               mpi_data->mpi_rvec, MASTER_NODE, mpi_data->world );
}

#endif
