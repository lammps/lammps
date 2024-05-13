/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_FIX_EVENT_H
#define LMP_FIX_EVENT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEvent : public Fix {
 public:
  FixEvent(class LAMMPS *, int, char **);
  ~FixEvent() override = 0;    // use destructor to make base class virtual
  int setmask() override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  void write_restart(FILE *) override {}
  void restart(char *) override {}

  // methods specific to FixEvent

  void store_event();              // store quenched atoms
  void restore_event();            // restore quenched atoms
  void store_state_quench();       // store hot atoms prior to quench
  void restore_state_quench();     // restore hot atoms after quench
  void store_state_dephase();      // store atoms before dephase iteration
  void restore_state_dephase();    // restore atoms if dephase had event

 private:
  double **xevent;        // atom coords at last event
  double **xold;          // atom coords for reset/restore
  double **vold;          // atom vels for reset/restore
  imageint *imageold;     // image flags for reset/restore
  double **xorig;         // original atom coords for reset/restore
  double **vorig;         // original atom vels for reset/restore
  imageint *imageorig;    // original image flags for reset/restore
};

}    // namespace LAMMPS_NS

#endif
