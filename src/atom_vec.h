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

#ifndef LMP_ATOM_VEC_H
#define LMP_ATOM_VEC_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class AtomVec : protected Pointers {
 public:
  int molecular;                       // 0 = atomic, 1 = molecular system
  int bonds_allow,angles_allow;        // 1 if bonds, angles are used
  int dihedrals_allow,impropers_allow; // 1 if dihedrals, impropers used
  int mass_type;                       // 1 if per-type masses
  int dipole_type;                     // 1 if per-type dipole moments

  int comm_x_only;                     // 1 if only exchange x in forward comm
  int comm_f_only;                     // 1 if only exchange f in reverse comm

  int size_forward;                    // # of values per atom in comm
  int size_reverse;                    // # in reverse comm
  int size_border;                     // # in border comm
  int size_velocity;                   // # of velocity based quantities
  int size_data_atom;                  // number of values in Atom line
  int size_data_vel;                   // number of values in Velocity line
  int size_data_bonus;                 // number of values in Bonus line
  int xcol_data;                       // column (1-N) where x is in Atom line

  int cudable;                         // 1 if atom style is CUDA-enabled
  int *maxsend;                        // CUDA-specific variable

  AtomVec(class LAMMPS *);
  virtual ~AtomVec() {}
  virtual void settings(int, char **);
  virtual void init();

  virtual void grow(int) = 0;
  virtual void grow_reset() = 0;
  virtual void copy(int, int, int) = 0;
  virtual void clear_bonus() {}

  virtual int pack_comm(int, int *, double *, int, int *) = 0;
  virtual int pack_comm_vel(int, int *, double *, int, int *) = 0;
  virtual int pack_comm_hybrid(int, int *, double *) {return 0;}
  virtual void unpack_comm(int, int, double *) = 0;
  virtual void unpack_comm_vel(int, int, double *) = 0;
  virtual int unpack_comm_hybrid(int, int, double *) {return 0;}

  virtual int pack_reverse(int, int, double *) = 0;
  virtual int pack_reverse_hybrid(int, int, double *) {return 0;}
  virtual void unpack_reverse(int, int *, double *) = 0;
  virtual int unpack_reverse_hybrid(int, int *, double *) {return 0;}

  virtual int pack_border(int, int *, double *, int, int *) = 0;
  virtual int pack_border_vel(int, int *, double *, int, int *) = 0;
  virtual int pack_border_hybrid(int, int *, double *) {return 0;}
  virtual void unpack_border(int, int, double *) = 0;
  virtual void unpack_border_vel(int, int, double *) = 0;
  virtual int unpack_border_hybrid(int, int, double *) {return 0;}

  virtual int pack_exchange(int, double *) = 0;
  virtual int unpack_exchange(double *) = 0;

  virtual int size_restart() = 0;
  virtual int pack_restart(int, double *) = 0;
  virtual int unpack_restart(double *) = 0;

  virtual void write_restart_settings(FILE *) {}
  virtual void read_restart_settings(FILE *) {}

  virtual void create_atom(int, double *) = 0;

  virtual void data_atom(double *, tagint, char **) = 0;
  virtual void data_atom_bonus(int, char **) {}
  virtual int data_atom_hybrid(int, char **) {return 0;}
  virtual void data_vel(int, char **);
  virtual int data_vel_hybrid(int, char **) {return 0;}

  virtual void pack_data(double **) = 0;
  virtual int pack_data_hybrid(int, double *) {return 0;}
  virtual void write_data(FILE *, int, double **) = 0;
  virtual int write_data_hybrid(FILE *, double *) {return 0;}
  virtual void pack_vel(double **);
  virtual int pack_vel_hybrid(int, double *) {return 0;}
  virtual void write_vel(FILE *, int, double **);
  virtual int write_vel_hybrid(FILE *, double *) {return 0;}

  int pack_bond(int **);
  void write_bond(FILE *, int, int **, int);
  int pack_angle(int **);
  void write_angle(FILE *, int, int **, int);
  void pack_dihedral(int **);
  void write_dihedral(FILE *, int, int **, int);
  void pack_improper(int **);
  void write_improper(FILE *, int, int **, int);

  virtual bigint memory_usage() = 0;

 protected:
  int nmax;                             // local copy of atom->nmax
  int deform_vremap;                    // local copy of domain properties
  int deform_groupbit;
  double *h_rate;
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid atom_style command

Self-explanatory.

E: USER-CUDA package requires a cuda enabled atom_style

Self-explanatory.

*/
