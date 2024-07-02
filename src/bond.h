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

#ifndef LMP_BOND_H
#define LMP_BOND_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class Bond : protected Pointers {
  friend class ThrOMP;
  friend class FixOMP;

 public:
  static int instance_total;    // # of Bond classes ever instantiated

  int allocated;
  int *setflag;
  int partial_flag;          // 1 if bond type can be set to 0 and deleted
  int writedata;             // 1 if writes coeffs to data file
  double energy;             // accumulated energies
  double virial[6];          // accumulated virial: xx,yy,zz,xy,xz,yz
  double *eatom, **vatom;    // accumulated per-atom energy/virial

  int born_matrix_enable;

  int comm_forward;        // size of forward communication (0 if none)
  int comm_reverse;        // size of reverse communication (0 if none)
  int comm_reverse_off;    // size of reverse comm even if newton off

  int reinitflag;    // 0 if not compatible with fix adapt
                     // extract() method may still need to be added

  int single_extra;    // number of extra single values calculated
  double *svector;     // vector of extra single quantities

  // KOKKOS host/device flag and data masks

  ExecutionSpace execution_space;
  unsigned int datamask_read, datamask_modify;
  int copymode, kokkosable;

  Bond(class LAMMPS *);
  ~Bond() override;
  virtual void init();
  virtual void init_style() {}
  virtual void compute(int, int) = 0;
  virtual void settings(int, char **);
  virtual void coeff(int, char **) = 0;
  virtual double equilibrium_distance(int) = 0;
  virtual void write_restart(FILE *) = 0;
  virtual void read_restart(FILE *) = 0;
  virtual void write_restart_settings(FILE *){};
  virtual void read_restart_settings(FILE *){};
  virtual void write_data(FILE *) {}
  virtual double single(int, double, int, int, double &) = 0;
  virtual double memory_usage();
  virtual void *extract(const char *, int &) { return nullptr; }
  virtual void reinit();

  virtual int pack_forward_comm(int, int *, double *, int, int *) { return 0; }
  virtual void unpack_forward_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) { return 0; }
  virtual void unpack_reverse_comm(int, int *, double *) {}

  virtual void born_matrix(int /*btype*/, double /*rsq*/, int /*at1*/, int /*at2*/, double &du,
                           double &du2)
  {
    du = 0.0;
    du2 = 0.0;
  }

  void write_file(int, char **);

 protected:
  int instance_me;    // which Bond class instantiation I am
  int suffix_flag;    // suffix compatibility flag

  int evflag;
  int eflag_either, eflag_global, eflag_atom;
  int vflag_either, vflag_global, vflag_atom;
  int maxeatom, maxvatom;

  void ev_init(int eflag, int vflag, int alloc = 1)
  {
    if (eflag || vflag)
      ev_setup(eflag, vflag, alloc);
    else
      evflag = eflag_either = eflag_global = eflag_atom = vflag_either = vflag_global = vflag_atom =
          0;
  }
  void ev_setup(int, int, int alloc = 1);
  void ev_tally(int, int, int, int, double, double, double, double, double);
  void ev_tally_xyz(int, int, int, int, double, double, double, double, double, double, double);
};

}    // namespace LAMMPS_NS

#endif
