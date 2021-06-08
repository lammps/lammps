/* -*- c++ -*- ----------------------------------------------------------
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
   The FixIMD class contains code from VMD and NAMD which is copyrighted
   by the Board of Trustees of the University of Illinois and is free to
   use with LAMMPS according to point 2 of the UIUC license (10% clause):

" Licensee may, at its own expense, create and freely distribute
complimentary works that interoperate with the Software, directing others to
the TCBG server to license and obtain the Software itself. Licensee may, at
its own expense, modify the Software to make derivative works.  Except as
explicitly provided below, this License shall apply to any derivative work
as it does to the original Software distributed by Illinois.  Any derivative
work should be clearly marked and renamed to notify users that it is a
modified version and not the original Software distributed by Illinois.
Licensee agrees to reproduce the copyright notice and other proprietary
markings on any derivative work and to include in the documentation of such
work the acknowledgement:

 "This software includes code developed by the Theoretical and Computational
  Biophysics Group in the Beckman Institute for Advanced Science and
  Technology at the University of Illinois at Urbana-Champaign."

Licensee may redistribute without restriction works with up to 1/2 of their
non-comment source code derived from at most 1/10 of the non-comment source
code developed by Illinois and contained in the Software, provided that the
above directions for notice and acknowledgement are observed.  Any other
distribution of the Software or any derivative work requires a separate
license with Illinois.  Licensee may contact Illinois (vmd@ks.uiuc.edu) to
negotiate an appropriate license for such distribution."
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(imd,FixIMD);
// clang-format on
#else

#ifndef LMP_FIX_IMD_H
#define LMP_FIX_IMD_H

#include "fix.h"

#if defined(LAMMPS_ASYNC_IMD)
#include <pthread.h>
#endif

/* prototype for c wrapper that calls the real worker */
extern "C" void *fix_imd_ioworker(void *);

namespace LAMMPS_NS {

class FixIMD : public Fix {
 public:
  FixIMD(class LAMMPS *, int, char **);
  ~FixIMD();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  double memory_usage();

 protected:
  int imd_port;
  void *localsock;
  void *clientsock;

  int num_coords;       // total number of atoms controlled by this fix
  int size_one;         // bytes per atom in communication buffer.
  int maxbuf;           // size of atom communication buffer.
  void *comm_buf;       // communication buffer
  void *idmap;          // hash for mapping atom indices to consistent order.
  tagint *rev_idmap;    // list of the hash keys for reverse mapping.

  int imd_forces;       // number of forces communicated via IMD.
  void *force_buf;      // force data buffer
  double imd_fscale;    // scale factor for forces. in case VMD's units are off.

  int imd_inactive;     // true if IMD connection stopped.
  int imd_terminate;    // true if IMD requests termination of run.
  int imd_trate;        // IMD transmission rate.

  int unwrap_flag;    // true if coordinates need to be unwrapped before sending
  int nowait_flag;    // true if LAMMPS should not wait with the execution for VMD.
  int connect_msg;    // flag to indicate whether a "listen for connection message" is needed.

  int me;               // my MPI rank in this "world".
  int nlevels_respa;    // flag to determine respa levels.

  int msglen;
  char *msgdata;

#if defined(LAMMPS_ASYNC_IMD)
  int buf_has_data;               // flag to indicate to the i/o thread what to do.
  pthread_mutex_t write_mutex;    // mutex for sending coordinates to i/o thread
  pthread_cond_t write_cond;      // conditional variable for coordinate i/o
  pthread_mutex_t read_mutex;     // mutex for accessing data receieved by i/o thread
  pthread_t iothread;             // thread id for i/o thread.
  pthread_attr_t iot_attr;        // i/o thread attributes.
 public:
  void ioworker(void);
#endif

 protected:
  int reconnect();
};

}    // namespace LAMMPS_NS

#endif
#endif

// Local Variables:
// mode: c++
// compile-command: "make -j4 openmpi"
// c-basic-offset: 2
// fill-column: 76
// indent-tabs-mode: nil
// End:
