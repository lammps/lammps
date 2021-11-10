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

#ifndef LMP_ATOM_VEC_H
#define LMP_ATOM_VEC_H

#include "pointers.h"    // IWYU pragma: export
#include <vector>

namespace LAMMPS_NS {

class AtomVec : protected Pointers {
 public:
  enum { PER_ATOM = 0, PER_TYPE = 1 };
  int molecular;    // 0 = atomic, 1 = molecular system, 2 = molecular template system
  int bonds_allow, angles_allow;           // 1 if bonds, angles are used
  int dihedrals_allow, impropers_allow;    // 1 if dihedrals, impropers used
  int mass_type;                           // 1 if per-type masses
  int dipole_type;                         // 1 if per-type dipole moments
  int forceclearflag;                      // 1 if has forceclear() method

  int comm_x_only;    // 1 if only exchange x in forward comm
  int comm_f_only;    // 1 if only exchange f in reverse comm

  int size_forward;      // # of values per atom in comm
  int size_reverse;      // # in reverse comm
  int size_border;       // # in border comm
  int size_velocity;     // # of velocity based quantities
  int size_data_atom;    // number of values in Atom line
  int size_data_vel;     // number of values in Velocity line
  int xcol_data;         // column (1-N) where x is in Atom line
  int maxexchange;       // max size of exchanged atom
                         // only needs to be set if size > BUFEXTRA

  int bonus_flag;                // 1 if stores bonus data
  int size_forward_bonus;        // # in forward bonus comm
  int size_border_bonus;         // # in border bonus comm
  int size_restart_bonus_one;    // # in restart bonus comm
  int size_data_bonus;           // number of values in Bonus line

  class Molecule **onemols;    // list of molecules for style template
  int nset;                    // # of molecules in list

  int kokkosable;    // 1 if atom style is KOKKOS-enabled

  int nargcopy;      // copy of command-line args for atom_style command
  char **argcopy;    // used when AtomVec is realloced (restart,replicate)

  // additional list of peratom fields operated on by different methods
  // set or created by child styles

  char *fields_grow, *fields_copy;
  char *fields_comm, *fields_comm_vel, *fields_reverse;
  char *fields_border, *fields_border_vel;
  char *fields_exchange, *fields_restart;
  char *fields_create, *fields_data_atom, *fields_data_vel;

  // methods

  AtomVec(class LAMMPS *);
  virtual ~AtomVec();

  void store_args(int, char **);
  virtual void process_args(int, char **);
  virtual void init();

  virtual void force_clear(int, size_t) {}
  virtual bigint roundup(bigint);

  virtual void grow(int);
  virtual void grow_pointers() {}
  virtual void copy(int, int, int);

  virtual void copy_bonus(int, int, int) {}
  virtual void clear_bonus() {}

  virtual int pack_comm(int, int *, double *, int, int *);
  virtual int pack_comm_vel(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  virtual void unpack_comm_vel(int, int, double *);

  virtual int pack_comm_bonus(int, int *, double *) { return 0; }
  virtual void unpack_comm_bonus(int, int, double *) {}

  virtual int pack_reverse(int, int, double *);
  virtual void unpack_reverse(int, int *, double *);

  virtual int pack_border(int, int *, double *, int, int *);
  virtual int pack_border_vel(int, int *, double *, int, int *);
  virtual void unpack_border(int, int, double *);
  virtual void unpack_border_vel(int, int, double *);

  virtual int pack_border_bonus(int, int *, double *) { return 0; }
  virtual int unpack_border_bonus(int, int, double *) { return 0; }

  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(double *);

  virtual int pack_exchange_bonus(int, double *) { return 0; }
  virtual int unpack_exchange_bonus(int, double *) { return 0; }

  virtual int size_restart();
  virtual int pack_restart(int, double *);
  virtual int unpack_restart(double *);

  virtual void pack_restart_pre(int) {}
  virtual void pack_restart_post(int) {}
  virtual void unpack_restart_init(int) {}

  virtual int size_restart_bonus() { return 0; }
  virtual int pack_restart_bonus(int, double *) { return 0; }
  virtual int unpack_restart_bonus(int, double *) { return 0; }

  virtual void create_atom(int, double *);
  virtual void create_atom_post(int) {}

  virtual void data_atom(double *, imageint, char **);
  virtual void data_atom_post(int) {}
  virtual void data_atom_bonus(int, char **) {}
  virtual void data_body(int, int, int, int *, double *) {}

  virtual void data_bonds_post(int, int, tagint, tagint, tagint) {}

  virtual void pack_data(double **);
  virtual void write_data(FILE *, int, double **);
  virtual void pack_data_pre(int) {}
  virtual void pack_data_post(int) {}

  virtual void data_vel(int, char **);
  virtual void pack_vel(double **);
  virtual void write_vel(FILE *, int, double **);

  virtual int pack_bond(tagint **);
  virtual void write_bond(FILE *, int, tagint **, int);
  virtual int pack_angle(tagint **);
  virtual void write_angle(FILE *, int, tagint **, int);
  virtual int pack_dihedral(tagint **);
  virtual void write_dihedral(FILE *, int, tagint **, int);
  virtual int pack_improper(tagint **);
  virtual void write_improper(FILE *, int, tagint **, int);

  virtual int pack_data_bonus(double *, int) { return 0; }
  virtual void write_data_bonus(FILE *, int, double *, int) {}

  virtual int property_atom(char *) { return -1; }
  virtual void pack_property_atom(int, double *, int, int) {}

  virtual double memory_usage();
  virtual double memory_usage_bonus() { return 0; }

  // old hybrid functions, needed by Kokkos package

  virtual int pack_comm_hybrid(int, int *, double *) { return 0; }
  virtual int unpack_comm_hybrid(int, int, double *) { return 0; }
  virtual int pack_reverse_hybrid(int, int, double *) { return 0; }
  virtual int unpack_reverse_hybrid(int, int *, double *) { return 0; }
  virtual int pack_border_hybrid(int, int *, double *) { return 0; }
  virtual int unpack_border_hybrid(int, int, double *) { return 0; }
  virtual int data_atom_hybrid(int, char **) { return 0; }
  virtual int data_vel_hybrid(int, char **) { return 0; }
  virtual int pack_data_hybrid(int, double *) { return 0; }
  virtual int write_data_hybrid(FILE *, double *) { return 0; }
  virtual int pack_vel_hybrid(int, double *) { return 0; }
  virtual int write_vel_hybrid(FILE *, double *) { return 0; }

 protected:
  int nmax;             // local copy of atom->nmax
  int deform_vremap;    // local copy of domain properties
  int deform_groupbit;
  double *h_rate;

  tagint *tag;    // peratom fields common to all styles
  int *type, *mask;
  imageint *image;
  double **x, **v, **f;

  // standard list of peratom fields always operated on by different methods
  // common to all styles, so not listed in field strings

  const char *default_grow, *default_copy;
  const char *default_comm, *default_comm_vel, *default_reverse;
  const char *default_border, *default_border_vel;
  const char *default_exchange, *default_restart;
  const char *default_create, *default_data_atom, *default_data_vel;

  struct Method {
    std::vector<void *> pdata;
    std::vector<int> datatype;
    std::vector<int> cols;
    std::vector<int *> maxcols;
    std::vector<int> collength;
    std::vector<void *> plength;
    std::vector<int> index;

    void resize(int nfield);
  };

  Method mgrow, mcopy;
  Method mcomm, mcomm_vel, mreverse, mborder, mborder_vel, mexchange, mrestart;
  Method mcreate, mdata_atom, mdata_vel;

  int ngrow, ncopy;
  int ncomm, ncomm_vel, nreverse, nborder, nborder_vel, nexchange, nrestart;
  int ncreate, ndata_atom, ndata_vel;

  // thread info for fields that are duplicated over threads
  // used by fields in grow() and memory_usage()

  bool *threads;

  // local methods

  void grow_nmax();
  int grow_nmax_bonus(int);
  void setup_fields();
  int process_fields(char *, const char *, Method *);
  void init_method(int, Method *);
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Invalid atom_style command

Self-explanatory.

E: KOKKOS package requires a kokkos enabled atom_style

Self-explanatory.

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

*/
