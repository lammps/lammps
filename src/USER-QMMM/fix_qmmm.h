/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


#ifdef FIX_CLASS

FixStyle(qmmm,FixQMMM)

#else

#ifndef LMP_FIX_QMMM_H
#define LMP_FIX_QMMM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixQMMM : public Fix {
 public:
  FixQMMM(class LAMMPS *, int, char **);
  ~FixQMMM();
  int setmask();
  void init();
//  void setup_pre_force(int);
//  void setup(int);
//  void pre_force(int);          // to send out updated positions
//  void post_force(int);         // to receive and update forces
  double memory_usage();

 protected:
  MPI_Comm qm_comm;   // intra communicator with QM subsystem
  MPI_Comm mm_comm;   // intra communicator with MM subsystem
  void  *comm_buf;    // message buffer for internal communication
  void  *qm_idmap;    // hash for mapping QM atom indices to consistent order.
  int   *qm_remap;    // list of the hash keys for reverse mapping.
  void  *mm_idmap;    // hash for mapping MM atom indices to consistent order.
  int   *mm_remap;    // list of the hash keys for reverse mapping.
  double *qm_coord;   // QM system coordinates
  double *mm_coord;   // MM system coordinates used for electrostatic coupling
  int    *mm_type;    // system atom types used for electrostatic coupling
  double *qm_force;   // QM force data buffer
  double *mm_force;   // MM slave force data buffer
  double qmmm_fscale; // scale factor for forces. in case VMD's units are off.

  int    num_qm;      // total number of QM atoms controlled by this fix
  int    num_mm;      // total number of MM atoms for electrostatic coupling
  int    mm_group;    // group of MM atoms for electrostatic coupling
  int    comm_mode;   // QM/MM communication method (MPI or shmemq)
  int    qmmm_mode;   // QM/MM coupling mode (mechanical or electrostatic)
  int    qmmm_role;   // role in QM/MM coupling (MM master or MM slave)
  int    size_one;    // size of one element in communication buffer
  int    do_init;     // flag for one time initialization
};

}

#endif
#endif

// Local Variables:
// mode: c++
// c-basic-offset: 2
// fill-column: 76
// indent-tabs-mode: nil
// End:
