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

#ifdef FIX_CLASS
// clang-format off
FixStyle(atc,FixATC);
// clang-format on
#else

#ifndef FIX_ATC_H
#define FIX_ATC_H

#include "fix.h"

namespace ATC {
class ATC_Method;
}
namespace LAMMPS_NS {
class NeighList;

/**
   *  @class FixATC
   *  @brief Class for an atom-to-continuum (ATC) LAMMPS fix.
   */

class FixATC : public Fix {
 public:
  /** constructor & destructor */
  FixATC(class LAMMPS *, int, char **);
  ~FixATC() override;

  /** initialization functions */
  void init() override;
  void init_list(int id, NeighList *ptr) override;
  void setup(int vflag) override;
  void min_setup(int vflag) override;

  /** setmask: tell LAMMPS which fix methods to call */
  int setmask() override;

  /** initial_integrate */
  void initial_integrate(int vflag) override;

  /** after first integrate phase */
  void post_integrate() override;

  /** final_integrate */
  void final_integrate() override;

  /** end of step for run or minimize */
  void end_of_step() override;

  /** pre_exchange is used to modify fix-specific data
       and is called before domain->pbc() and comm->exchange().  */
  void setup_pre_exchange() override;
  void pre_exchange() override;
  void min_pre_exchange() override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;

  /** pack_exchange called from atom_vec->pack_exchange()
       and packs fix-specific data for a given real (local)
       atom being sent to another processor.  */
  int pack_exchange(int, double *) override;

  /** unpack_exchange called from atom_vec->unpack_exchange()
       and unpacks fix-specific data for a given real (local)
       atom received from another processor. */
  int unpack_exchange(int, double *) override;

  /** pack_comm called from comm->forward_comm and
       packs fix-specific data for a given ghost atom
       from exchange with another proc */
  int pack_forward_comm(int, int *, double *, int, int *) override;

  /** unpack_comm called from comm->forward_comm and
       unpacks fix-specific data for a given ghost atom
       from exchange with another proc */
  void unpack_forward_comm(int, int, double *) override;

  /** pre_neighbor is used to modify fix-specific data
       and is called before neighbor list is built in
       neighbor->build().  */
  void pre_neighbor() override;
  void setup_pre_neighbor() override;

  /** pre/post_force is used to modify fix-specific data
        and is before/after the various force computations. */
  void pre_force(int vflag) override;
  void post_force(int vflag) override;

  /** post_run is called after a run completes */
  void post_run() override;

  /** min_pre_force is called before forces are calculated in minimize */
  void min_pre_force(int vflag) override;

  /** min_post_force is called after forces are calculated in minimize */
  void min_post_force(int vflag) override;

  /** modify atc parameters (parser) */
  int modify_param(int narg, char **arg) override;

  /** calls ATC_Method to handle restarting/checkpointing */
  /** these four methods are for writing per-atom quantities */
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;
  /** these two methods are for writing all other quantities */
  void write_restart(FILE *) override;
  void restart(char *) override;

  /** accessor function for ATC_Method class pointer */
  const ATC::ATC_Method *atc() { return atc_; }

 protected:
  LAMMPS *lammps_;

  /** functions for "thermo" output */
  double compute_scalar() override;
  double compute_vector(int n) override;
  double compute_array(int irow, int icol) override;
  double dtv, dtf;
  ATC::ATC_Method *atc_;
};
}    // namespace LAMMPS_NS

#endif
#endif
