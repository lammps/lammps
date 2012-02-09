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

#ifndef LMP_IRREGULAR_H
#define LMP_IRREGULAR_H

#include "pointers.h"

namespace LAMMPS_NS {

class Irregular : protected Pointers {
 public:
  Irregular(class LAMMPS *);
  ~Irregular();
  void migrate_atoms();
  int create_data(int, int *);
  void exchange_data(char *, int, char *);
  void destroy_data();
  bigint memory_usage();

 private:
  int me,nprocs;
  int triclinic;
  int map_style;
  int uniform;
  double *xsplit,*ysplit,*zsplit;   // ptrs to comm
  int *procgrid;                    // ptr to comm
  int ***grid2proc;                 // ptr to comm
  double *boxlo;                    // ptr to domain 
  double *prd;                      // ptr to domain

  int maxsend,maxrecv;              // size of buffers in # of doubles
  double *buf_send,*buf_recv;

  // plan for irregular communication of atoms
  // no params refer to atoms copied to self

  struct PlanAtom {
    int nsend;                 // # of messages to send
    int nrecv;                 // # of messages to recv
    int sendmax;               // # of doubles in largest send message
    int *proc_send;            // procs to send to
    int *length_send;          // # of doubles to send to each proc
    int *num_send;             // # of atoms to send to each proc
    int *index_send;           // list of which atoms to send to each proc
    int *offset_send;          // where each atom starts in send buffer
    int *proc_recv;            // procs to recv from
    int *length_recv;          // # of doubles to recv from each proc
    MPI_Request *request;      // MPI requests for posted recvs
    MPI_Status *status;        // MPI statuses for WaitAll
  };

  // plan for irregular communication of datums
  // only 2 self params refer to atoms copied to self

  struct PlanData {            // plan for irregular communication of data
    int nsend;                 // # of messages to send
    int nrecv;                 // # of messages to recv
    int sendmax;               // # of datums in largest send message
    int *proc_send;            // procs to send to
    int *num_send;             // # of datums to send to each proc
    int *index_send;           // list of which datums to send to each proc
    int *proc_recv;            // procs to recv from
    int *num_recv;             // # of datums to recv from each proc
    int num_self;              // # of datums to copy to self
    int *index_self;           // list of which datums to copy to self
    MPI_Request *request;      // MPI requests for posted recvs
    MPI_Status *status;        // MPI statuses for WaitAll
  };

  PlanAtom *aplan;
  PlanData *dplan;

  int create_atom(int, int *, int *);
  void exchange_atom(double *, int *, double *);
  void destroy_atom();
  int coord2proc(double *);
  int binary(double, int, double *);

  void grow_send(int,int);          // reallocate send buffer
  void grow_recv(int);              // free/allocate recv buffer
};

}

#endif
/* ERROR/WARNING messages:

*/
