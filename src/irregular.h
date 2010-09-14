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
  struct PlanAtom {            // plan for irregular communication of atoms
    int nsend;                 // # of messages to send
    int nrecv;                 // # of messages to recv
    int sendmax;               // # of doubles in largest send message
    int *proc_send;            // procs to send to
    int *length_send;          // # of doubles to send to each proc
    int *num_send;             // # of datums to send to each proc
    int *index_send;           // list of which datums to send to each proc
    int *offset_send;          // where each datum starts in send buffer
    int *proc_recv;            // procs to recv from
    int *length_recv;          // # of doubles to recv from each proc
    MPI_Request *request;      // MPI requests for posted recvs
    MPI_Status *status;        // MPI statuses for WaitAll
  };

  struct PlanData {            // plan for irregular communication of data
    int nsend;                 // # of messages to send
    int nrecv;                 // # of messages to recv
    int sendmax;               // # of datums in largest send message
    int *proc_send;            // procs to send to
    int *num_send;             // # of datums to send to each proc
    int *index_send;           // list of which datums to send to each proc
    int *proc_recv;            // procs to recv from
    int *num_recv;             // # of doubles to recv from each proc
    int num_self;
    int *index_self;
    MPI_Request *request;      // MPI requests for posted recvs
    MPI_Status *status;        // MPI statuses for WaitAll
  };

  Irregular(class LAMMPS *);
  ~Irregular();
  void migrate_atoms();
  struct PlanData *create_data(int, int *, int *);
  void exchange_data(PlanData *, char *, int, char *);
  void destroy_data(PlanData *);

 private:
  int me,nprocs;
  int triclinic;
  int map_style;
  int *procgrid;
  int ***grid2proc;
  int maxsend,maxrecv;
  double *buf_send,*buf_recv;

  struct PlanAtom *create_atom(int, int *, int *, int *);
  void exchange_atom(PlanAtom *, double *, int *, double *);
  void destroy_atom(PlanAtom *);
  int coord2proc(double *);

  void grow_send(int,int);          // reallocate send buffer
  void grow_recv(int);              // free/allocate recv buffer
};

}

#endif
