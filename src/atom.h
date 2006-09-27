/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef ATOM_H
#define ATOM_H

#include "lammps.h"

class Atom : public LAMMPS {
 public:
  char *style;
  double PI;

  // atom counts

  double natoms;                // total # of atoms in system, could be 0
  int nlocal,nghost;            // # of owned and ghost atoms on this proc
  int nmax;                     // max # of owned+ghost in arrays on this proc

  // # of atom and topology types

  int ntypes,nbondtypes,nangletypes,ndihedraltypes,nimpropertypes;

  // molecular counts

  int nbonds,nangles,ndihedrals,nimpropers;
  int bond_per_atom,angle_per_atom,dihedral_per_atom,improper_per_atom;

  // atom styles

  int style_angle,style_atomic,style_bond,style_charge,style_dipole,style_dpd;
  int style_eam,style_full,style_granular,style_molecular,style_peri;
  int style_hybrid;

  // flags set by atom styles

  int molecular;                       // 0 = atomic, 1 = molecular system
  int bonds_allow,angles_allow;        // 0/1 if bonds, angles are used
  int dihedrals_allow,impropers_allow; // 0/1 if dihedrals, impropers used
  int charge_allow;                    // 0/1 if charges used
  int mass_allow;                      // 0/1 if per-atom rmass array
  int mass_require;                    // 0/1 if per-type masses
  int dipole_require;                  // 0/1 if per-type dipole moments
  int tag_enable;                      // 0/1 if atom ID tags are defined

  // communication sizes - set by each atom style

  int size_comm,size_reverse,size_border;
  int size_atom_valid,size_atom_actual;

  // arrays of length ntypes and flags if set

  double *mass;
  int *mass_setflag;
  double *dipole;
  int *dipole_setflag;

  // arrays common to all atom classes

  int *tag,*type,*mask,*image;
  double **x,**v,**f;

  // arrays specific to some atom classes

  double *q,**mu;
  double **omega,**torque;
  double *radius,*density,*rmass,*vfrac;
  double **phix,**phiv,**phia;

  int *molecule;

  int maxspecial;
  int **nspecial,**special;

  int *num_bond;
  int **bond_type;
  int **bond_atom;

  int *num_angle;
  int **angle_type;
  int **angle_atom1,**angle_atom2,**angle_atom3;

  int *num_dihedral;
  int **dihedral_type;
  int **dihedral_atom1,**dihedral_atom2,**dihedral_atom3,**dihedral_atom4;

  int *num_improper;
  int **improper_type;
  int **improper_atom1,**improper_atom2,**improper_atom3,**improper_atom4;
 
  // array for extra peratom info in restart file destined for fix & diag 

  double **extra;

  // callback ptrs for atom arrays managed by fix classes

  int nextra_grow,nextra_restart;             // # of callbacks of each type
  int *extra_grow,*extra_restart;             // index of fix to callback to
  int nextra_grow_max,nextra_restart_max;     // size of callback lists
  int nextra_store;

  // data for global to local ID mapping

  int map_style;                  // default or user-specified style of map
                                  // 0 = none, 1 = array, 2 = hash
  int map_tag_max;
  int *map_array;

  struct HashElem {
    int global;                   // key to search on = global ID
    int local;                    // value associated with key = local index
    int next;                     // next entry in this bucket, -1 if last
  };
  int map_nhash;                  // # of entries hash table can hold
  int map_nused;                  // # of actual entries in hash table
  int map_free;                   // ptr to 1st unused entry in hash table
  int map_nbucket;                // # of hash buckets
  int *map_bucket;                // ptr to 1st entry in each bucket
  HashElem *map_hash;             // hash table
  int *primes;                    // table of prime #s for hashing
  int nprimes;                    // # of primes

  // functions common to all atom styles

  Atom(int, char **);
  virtual ~Atom();
  void set_style(char *);
  int check_style(char *);
  int style2arg(char **&);
  char *style2word(char *);
  void settings(Atom *);

  void init();
  void grow(int);

  void modify_params(int, char **);
  void tag_extend();
  int tag_consecutive();

  int parse_data(char *);
  int count_words(char *);

  void unpack_data(int, char *);
  void create_one(int, double, double, double);
  virtual void unpack_vels(int, char *);     // can be overwritten by child
  void unpack_bonds(int, char *);
  void unpack_angles(int, char *);
  void unpack_dihedrals(int, char *);
  void unpack_impropers(int, char *);

  void allocate_type_arrays();
  void set_mass(char *);
  void set_mass(int, double);
  void set_mass(int, char **);
  void set_mass(double *);
  void check_mass();
  void set_dipole(char *);
  void set_dipole(int, char **);
  void set_dipole(double *);
  void check_dipole();

  void add_callback(int);
  void delete_callback(char *, int);
  void update_callback(int);

  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);

  int memory_usage();

  // functions for global ID to local index atom mapping
  // map lookup function is inlined right here for efficiency
  
  inline int map(int global) {
    if (map_style == 1) return map_array[global];
    else return map_find_hash(global);
  };

  void map_init();
  void map_clear();
  void map_set();
  void map_one(int, int);
  void map_delete();
  int map_find_hash(int);

  // pure virtual functions, must be defined in child class

  virtual void copy(int, int) = 0;

  virtual void pack_comm(int, int *, double *, int *) = 0;
  virtual void unpack_comm(int, int, double *) = 0;
  virtual void pack_reverse(int, int, double *) = 0;
  virtual void unpack_reverse(int, int *, double *) = 0;

  virtual void pack_border(int, int *, double *, int *) = 0;
  virtual void unpack_border(int, int, double *) = 0;
  virtual int pack_exchange(int, double *) = 0;
  virtual int unpack_exchange(double *) = 0;
};

#endif
