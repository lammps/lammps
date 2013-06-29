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

#ifndef LMP_PAIR_H
#define LMP_PAIR_H

#include "pointers.h"

namespace LAMMPS_NS {

class Pair : protected Pointers {
  friend class AngleSDK;
  friend class AngleSDKOMP;
  friend class BondQuartic;
  friend class BondQuarticOMP;
  friend class DihedralCharmm;
  friend class DihedralCharmmOMP;
  friend class FixGPU;
  friend class FixOMP;
  friend class ThrOMP;

 public:
  double eng_vdwl,eng_coul;      // accumulated energies
  double virial[6];              // accumulated virial
  double *eatom,**vatom;         // accumulated per-atom energy/virial

  double cutforce;               // max cutoff for all atom pairs
  double **cutsq;                // cutoff sq for each atom pair
  int **setflag;                 // 0/1 = whether each i,j has been set

  int comm_forward;              // size of forward communication (0 if none)
  int comm_reverse;              // size of reverse communication (0 if none)
  int comm_reverse_off;          // size of reverse comm even if newton off

  int single_enable;             // 1 if single() routine exists
  int restartinfo;               // 1 if pair style writes restart info
  int respa_enable;              // 1 if inner/middle/outer rRESPA routines
  int one_coeff;                 // 1 if allows only one coeff * * call
  int manybody_flag;             // 1 if a manybody potential
  int no_virial_fdotr_compute;   // 1 if does not invoke virial_fdotr_compute()
  int writedata;                 // 1 if writes coeffs to data file
  int ghostneigh;                // 1 if pair style needs neighbors of ghosts
  double **cutghost;             // cutoff for each ghost pair

  int ewaldflag;                 // 1 if compatible with Ewald solver
  int pppmflag;                  // 1 if compatible with PPPM solver
  int msmflag;                   // 1 if compatible with MSM solver
  int dispersionflag;            // 1 if compatible with LJ/dispersion solver
  int tip4pflag;                 // 1 if compatible with TIP4P solver
  int dipoleflag;                // 1 if compatible with dipole solver

  int tail_flag;                 // pair_modify flag for LJ tail correction
  double etail,ptail;            // energy/pressure tail corrections
  double etail_ij,ptail_ij;

  int evflag;                    // energy,virial settings
  int eflag_either,eflag_global,eflag_atom;
  int vflag_either,vflag_global,vflag_atom;

  int ncoultablebits;            // size of Coulomb table, accessed by KSpace
  int ndisptablebits;            // size of dispersion table
  double tabinnersq;
  double tabinnerdispsq;
  double *rtable,*drtable,*ftable,*dftable,*ctable,*dctable;
  double *etable,*detable,*ptable,*dptable,*vtable,*dvtable;
  double *rdisptable, *drdisptable, *fdisptable, *dfdisptable;
  double *edisptable, *dedisptable;
  int ncoulshiftbits,ncoulmask;
  int ndispshiftbits, ndispmask;

  int nextra;                    // # of extra quantities pair style calculates
  double *pvector;               // vector of extra pair quantities

  int single_extra;              // number of extra single values calculated
  double *svector;               // vector of extra single quantities

  class NeighList *list;         // standard neighbor list used by most pairs
  class NeighList *listhalf;     // half list used by some pairs
  class NeighList *listfull;     // full list used by some pairs
  class NeighList *listgranhistory;  // granular history list used by some pairs
  class NeighList *listinner;    // rRESPA lists used by some pairs
  class NeighList *listmiddle;
  class NeighList *listouter;

  unsigned int datamask;
  unsigned int datamask_ext;

  int compute_flag;              // 0 if skip compute()

  Pair(class LAMMPS *);
  virtual ~Pair();

  // top-level Pair methods

  void init();
  void reinit();
  double mix_energy(double, double, double, double);
  double mix_distance(double, double);
  void write_file(int, char **);
  void init_bitmap(double, double, int, int &, int &, int &, int &);
  virtual void modify_params(int, char **);
  void compute_dummy(int, int);

  // need to be public, so can be called by pair_style reaxc

  void v_tally(int, double *, double *);
  void ev_tally(int, int, int, int, double, double, double,
                double, double, double);
  void ev_tally3(int, int, int, double, double,
                 double *, double *, double *, double *);
  void v_tally3(int, int, int, double *, double *, double *, double *);
  void v_tally4(int, int, int, int, double *, double *, double *,
                double *, double *, double *);
  void ev_tally_xyz(int, int, int, int, double, double,
                    double, double, double, double, double, double);

  // general child-class methods

  virtual void compute(int, int) = 0;
  virtual void compute_inner() {}
  virtual void compute_middle() {}
  virtual void compute_outer(int, int) {}

  virtual double single(int, int, int, int,
                        double, double, double, 
			double& fforce) {
    fforce = 0.0;
    return 0.0;
  }

  virtual void settings(int, char **) = 0;
  virtual void coeff(int, char **) = 0;

  virtual void init_style();
  virtual void init_list(int, class NeighList *);
  virtual double init_one(int, int) {return 0.0;}

  virtual void init_tables(double, double *);
  virtual void init_tables_disp(double);
  virtual void free_tables();
  virtual void free_disp_tables();

  virtual void write_restart(FILE *) {}
  virtual void read_restart(FILE *) {}
  virtual void write_restart_settings(FILE *) {}
  virtual void read_restart_settings(FILE *) {}
  virtual void write_data(FILE *) {}
  virtual void write_data_all(FILE *) {}

  virtual int pack_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_reverse_comm(int, int *, double *) {}
  virtual double memory_usage();

  // specific child-class methods for certain Pair styles

  virtual void *extract(const char *, int &) {return NULL;}
  virtual void swap_eam(double *, double **) {}
  virtual void reset_dt() {}
  virtual void min_xf_pointers(int, double **, double **) {}
  virtual void min_xf_get(int) {}
  virtual void min_x_set(int) {}

  virtual unsigned int data_mask() {return datamask;}
  virtual unsigned int data_mask_ext() {return datamask_ext;}

 protected:
  enum{GEOMETRIC,ARITHMETIC,SIXTHPOWER};   // mixing options

  int allocated;               // 0/1 = whether arrays are allocated
  int suffix_flag;             // suffix compatibility flag

                                       // pair_modify settings
  int offset_flag,mix_flag;            // flags for offset and mixing
  double tabinner;                     // inner cutoff for Coulomb table
  double tabinner_disp;                 // inner cutoff for dispersion table

  // custom data type for accessing Coulomb tables

  typedef union {int i; float f;} union_int_float_t;

  double THIRD;

  int vflag_fdotr;
  int maxeatom,maxvatom;

  virtual void ev_setup(int, int);
  void ev_unset();
  void ev_tally_full(int, double, double, double, double, double, double);
  void ev_tally_xyz_full(int, double, double,
                         double, double, double, double, double, double);
  void ev_tally4(int, int, int, int, double,
                 double *, double *, double *, double *, double *, double *);
  void ev_tally_tip4p(int, int *, double *, double, double);
  void v_tally2(int, int, double, double *);
  void v_tally_tensor(int, int, int, int,
                      double, double, double, double, double, double);
  void virial_fdotr_compute();

  FILE *open_potential(const char *);
  const char *potname(const char *);

  inline int sbmask(int j) {
    return j >> SBBITS & 3;
  }
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Too many total bits for bitmapped lookup table

Table size specified via pair_modify command is too large.  Note that
a value of N generates a 2^N size table.

E: Cannot have both pair_modify shift and tail set to yes

These 2 options are contradictory.

E: Cannot use pair tail corrections with 2d simulations

The correction factors are only currently defined for 3d systems.

W: Using pair tail corrections with nonperiodic system

This is probably a bogus thing to do, since tail corrections are
computed by integrating the density of a periodic system out to
infinity.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Pair style does not support pair_write

The pair style does not have a single() function, so it can
not be invoked by pair write.

E: Invalid atom types in pair_write command

Atom types must range from 1 to Ntypes inclusive.

E: Invalid style in pair_write command

Self-explanatory.  Check the input script.

E: Invalid cutoffs in pair_write command

Inner cutoff must be larger than 0.0 and less than outer cutoff.

E: Cannot open pair_write file

The specified output file for pair energies and forces cannot be
opened.  Check that the path and name are correct.

E: Bitmapped lookup tables require int/float be same size

Cannot use pair tables on this machine, because of word sizes.  Use
the pair_modify command with table 0 instead.

W: Table inner cutoff >= outer cutoff

You specified an inner cutoff for a Coulombic table that is longer
than the global cutoff.  Probably not what you wanted.

E: Too many exponent bits for lookup table

Table size specified via pair_modify command does not work with your
machine's floating point representation.

E: Too many mantissa bits for lookup table

Table size specified via pair_modify command does not work with your
machine's floating point representation.

E: Too few bits for lookup table

Table size specified via pair_modify command does not work with your
machine's floating point representation.

*/
