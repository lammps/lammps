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

#ifndef LMP_ATOM_H
#define LMP_ATOM_H

#include "pointers.h"
#include <map>
#include <string>

namespace LAMMPS_NS {

class Atom : protected Pointers {
 public:
  char *atom_style;
  class AtomVec *avec;

  // atom counts

  bigint natoms;                // total # of atoms in system, could be 0
                                // natoms may not be current if atoms lost
  int nlocal,nghost;            // # of owned and ghost atoms on this proc
  int nmax;                     // max # of owned+ghost in arrays on this proc
  int tag_enable;               // 0/1 if atom ID tags are defined
  int molecular;                // 0 = atomic, 1 = standard molecular system,
                                // 2 = molecule template system
  bigint nellipsoids;           // number of ellipsoids
  bigint nlines;                // number of lines
  bigint ntris;                 // number of triangles
  bigint nbodies;               // number of bodies

  bigint nbonds,nangles,ndihedrals,nimpropers;
  int ntypes,nbondtypes,nangletypes,ndihedraltypes,nimpropertypes;
  int bond_per_atom,angle_per_atom,dihedral_per_atom,improper_per_atom;
  int extra_bond_per_atom,extra_angle_per_atom;
  int extra_dihedral_per_atom,extra_improper_per_atom;

  int firstgroup;               // store atoms in this group first, -1 if unset
  int nfirst;                   // # of atoms in first group on this proc
  char *firstgroupname;         // group-ID to store first, NULL if unset

  // per-atom arrays
  // customize by adding new array

  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;

  tagint *molecule;
  int *molindex,*molatom;

  double *q,**mu;
  double **omega,**angmom,**torque;
  double *radius,*rmass;
  int *ellipsoid,*line,*tri,*body;

  // SPIN package

  double *emag;
  double **sp;
  double **fm;
  double **fm_long;

  // PERI package

  double *vfrac,*s0;
  double **x0;

  // USER-EFF and USER-AWPMD packages

  int *spin;
  double *eradius,*ervel,*erforce,*ervelforce;
  double *cs,*csforce,*vforce;
  int *etag;

  // USER-SPH package

  double *rho,*drho,*e,*de,*cv;
  double **vest;

  // USER-SMD package

  double *contact_radius;
  double **smd_data_9;
  double **smd_stress;
  double *eff_plastic_strain;
  double *eff_plastic_strain_rate;
  double *damage;

  // USER-DPD package

  double *uCond,*uMech,*uChem,*uCGnew,*uCG;
  double *duChem;
  double *dpdTheta;
  int nspecies_dpd;

  // USER-MESO package

  double **cc, **cc_flux;        // cc = chemical concentration
  double *edpd_temp,*edpd_flux;  // temperature and heat flux
  double *edpd_cv;               // heat capacity
  int cc_species;

  // molecular info

  int **nspecial;               // 0,1,2 = cumulative # of 1-2,1-3,1-4 neighs
  tagint **special;             // IDs of 1-2,1-3,1-4 neighs of each atom
  int maxspecial;               // special[nlocal][maxspecial]

  int *num_bond;
  int **bond_type;
  tagint **bond_atom;

  int *num_angle;
  int **angle_type;
  tagint **angle_atom1,**angle_atom2,**angle_atom3;

  int *num_dihedral;
  int **dihedral_type;
  tagint **dihedral_atom1,**dihedral_atom2,**dihedral_atom3,**dihedral_atom4;

  int *num_improper;
  int **improper_type;
  tagint **improper_atom1,**improper_atom2,**improper_atom3,**improper_atom4;

  // custom arrays used by fix property/atom

  int **ivector;
  double **dvector;
  char **iname,**dname;
  int nivector,ndvector;

  // atom style and per-atom array existence flags
  // customize by adding new flag

  int sphere_flag,ellipsoid_flag,line_flag,tri_flag,body_flag;
  int peri_flag,electron_flag;
  int ecp_flag;
  int wavepacket_flag,sph_flag;

  int molecule_flag,molindex_flag,molatom_flag;
  int q_flag,mu_flag;
  int rmass_flag,radius_flag,omega_flag,torque_flag,angmom_flag;
  int vfrac_flag,spin_flag,eradius_flag,ervel_flag,erforce_flag;
  int cs_flag,csforce_flag,vforce_flag,ervelforce_flag,etag_flag;
  int rho_flag,e_flag,cv_flag,vest_flag;
  int dpd_flag,edpd_flag,tdpd_flag;

  //USER-SPIN package

  int sp_flag;

  // USER-SMD package

  int smd_flag;
  int contact_radius_flag;
  int smd_data_9_flag;
  int smd_stress_flag;
  int x0_flag;
  int eff_plastic_strain_flag;
  int eff_plastic_strain_rate_flag;
  int damage_flag;

  // Peridynamics scale factor, used by dump cfg

  double pdscale;

  // molecule templates
  // each template can be a set of consecutive molecules
  // each with same ID (stored in molecules)
  // 1st molecule in template stores nset = # in set

  int nmolecule;
  class Molecule **molecules;

  // extra peratom info in restart file destined for fix & diag

  double **extra;

  // per-type arrays

  double *mass;
  int *mass_setflag;

  // callback ptrs for atom arrays managed by fix classes

  int nextra_grow,nextra_restart,nextra_border;  // # of callbacks of each type
  int *extra_grow,*extra_restart,*extra_border;  // index of fix to callback to
  int nextra_grow_max,nextra_restart_max;        // size of callback lists
  int nextra_border_max;
  int nextra_store;

  int map_style;                  // style of atom map: 0=none, 1=array, 2=hash
  int map_user;                   // user requested map style:
                                  // 0 = no request, 1=array, 2=hash, 3=yes
  tagint map_tag_max;             // max atom ID that map() is setup for

  // spatial sorting of atoms

  int sortfreq;             // sort atoms every this many steps, 0 = off
  bigint nextsort;          // next timestep to sort on
  double userbinsize;       // requested sort bin size

  // indices of atoms with same ID

  int *sametag;      // sametag[I] = next atom with same ID, -1 if no more

  // AtomVec factory types and map

  typedef AtomVec *(*AtomVecCreator)(LAMMPS *);
  typedef std::map<std::string,AtomVecCreator> AtomVecCreatorMap;
  AtomVecCreatorMap *avec_map;

  // functions

  Atom(class LAMMPS *);
  ~Atom();

  void settings(class Atom *);
  void create_avec(const char *, int, char **, int);
  virtual class AtomVec *new_avec(const char *, int, int &);
  void init();
  void setup();

  class AtomVec *style_match(const char *);
  void modify_params(int, char **);
  void tag_check();
  void tag_extend();
  int tag_consecutive();

  void bonus_check();

  int parse_data(const char *);
  int count_words(const char *);
  int count_words(const char *, char *);

  void deallocate_topology();

  void data_atoms(int, char *, tagint, tagint, int, int, double *);
  void data_vels(int, char *, tagint);
  void data_bonds(int, char *, int *, tagint, int);
  void data_angles(int, char *, int *, tagint, int);
  void data_dihedrals(int, char *, int *, tagint, int);
  void data_impropers(int, char *, int *, tagint, int);
  void data_bonus(int, char *, class AtomVec *, tagint);
  void data_bodies(int, char *, class AtomVecBody *, tagint);
  void data_fix_compute_variable(int, int);

  virtual void allocate_type_arrays();
  void set_mass(const char *, int, const char *, int);
  void set_mass(const char *, int, int, double);
  void set_mass(const char *, int, int, char **);
  void set_mass(double *);
  void check_mass(const char *, int);

  int radius_consistency(int, double &);
  int shape_consistency(int, double &, double &, double &);

  void add_molecule(int, char **);
  int find_molecule(char *);
  void add_molecule_atom(class Molecule *, int, int, tagint);

  void first_reorder();
  virtual void sort();

  void add_callback(int);
  void delete_callback(const char *, int);
  void update_callback(int);

  int find_custom(const char *, int &);
  virtual int add_custom(const char *, int);
  virtual void remove_custom(int, int);

  virtual void sync_modify(ExecutionSpace, unsigned int, unsigned int) {}

  void *extract(char *);

  inline int* get_map_array() {return map_array;};
  inline int get_map_size() {return map_tag_max+1;};
  inline int get_map_maxarray() {return map_maxarray+1;};

  bigint memory_usage();
  int memcheck(const char *);

  // functions for global to local ID mapping
  // map lookup function inlined for efficiency
  // return -1 if no map defined

  inline int map(tagint global) {
    if (map_style == 1) return map_array[global];
    else if (map_style == 2) return map_find_hash(global);
    else return -1;
  };

  void map_init(int check = 1);
  void map_clear();
  void map_set();
  void map_one(tagint, int);
  int map_style_set();
  void map_delete();
  int map_find_hash(tagint);

 protected:

  // global to local ID mapping

  int *map_array;       // direct map via array that holds map_tag_max
  int map_maxarray;     // allocated size of map_array (1 larger than this)

  struct HashElem {     // hashed map
    tagint global;      // key to search on = global ID
    int local;          // value associated with key = local index
    int next;           // next entry in this bucket, -1 if last
  };
  int map_nhash;        // # of entries hash table can hold
  int map_nused;        // # of actual entries in hash table
  int map_free;         // ptr to 1st unused entry in hash table
  int map_nbucket;      // # of hash buckets
  int *map_bucket;      // ptr to 1st entry in each bucket
  HashElem *map_hash;   // hash table

  int max_same;         // allocated size of sametag

  // spatial sorting of atoms

  int nbins;                      // # of sorting bins
  int nbinx,nbiny,nbinz;          // bins in each dimension
  int maxbin;                     // max # of bins
  int maxnext;                    // max size of next,permute
  int *binhead;                   // 1st atom in each bin
  int *next;                      // next atom in bin
  int *permute;                   // permutation vector
  double bininvx,bininvy,bininvz; // inverse actual bin sizes
  double bboxlo[3],bboxhi[3];     // bounding box of my sub-domain

  int memlength;                  // allocated size of memstr
  char *memstr;                   // string of array names already counted

  void setup_sort_bins();
  int next_prime(int);

 private:
  template <typename T> static AtomVec *avec_creator(LAMMPS *);
};

}

#endif

/* ERROR/WARNING messages:

E: Atom IDs must be used for molecular systems

Atom IDs are used to identify and find partner atoms in bonds.

E: Unrecognized atom style

The choice of atom style is unknown.

E: Could not find atom_modify first group ID

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Atom_modify id command after simulation box is defined

The atom_modify id command cannot be used after a read_data,
read_restart, or create_box command.

E: Atom_modify map command after simulation box is defined

The atom_modify map command cannot be used after a read_data,
read_restart, or create_box command.

E: Atom_modify sort and first options cannot be used together

Self-explanatory.

E: One or more Atom IDs is negative

Atom IDs must be positive integers.

E: One or more atom IDs is too big

The limit on atom IDs is set by the SMALLBIG, BIGBIG, SMALLSMALL
setting in your Makefile.  See Section_start 2.2 of the manual for
more details.

E: One or more atom IDs is zero

Either all atoms IDs must be zero or none of them.

E: Non-zero atom IDs with atom_modify id = no

Self-explanatory.

E: All atom IDs = 0 but atom_modify id = yes

Self-explanatory.

E: Duplicate atom IDs exist

Self-explanatory.

E: New atom IDs exceed maximum allowed ID

See the setting for tagint in the src/lmptype.h file.

E: Incorrect atom format in data file

Number of values per atom line in the data file is not consistent with
the atom style.

E: Incorrect format of ... section in data file

Number or type of values per line in the given section of the data file
is not consistent with the requirements for this section.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

E: Incorrect velocity format in data file

Each atom style defines a format for the Velocity section
of the data file.  The read-in lines do not match.

E: Invalid atom ID in Velocities section of data file

Atom IDs must be positive integers and within range of defined
atoms.

E: Invalid atom ID in Bonds section of data file

Atom IDs must be positive integers and within range of defined
atoms.

E: Invalid bond type in Bonds section of data file

Bond type must be positive integer and within range of specified bond
types.

E: Invalid atom ID in Angles section of data file

Atom IDs must be positive integers and within range of defined
atoms.

E: Invalid angle type in Angles section of data file

Angle type must be positive integer and within range of specified angle
types.

E: Invalid atom ID in Dihedrals section of data file

Atom IDs must be positive integers and within range of defined
atoms.

E: Invalid dihedral type in Dihedrals section of data file

Dihedral type must be positive integer and within range of specified
dihedral types.

E: Invalid atom ID in Impropers section of data file

Atom IDs must be positive integers and within range of defined
atoms.

E: Invalid improper type in Impropers section of data file

Improper type must be positive integer and within range of specified
improper types.

E: Incorrect bonus data format in data file

See the read_data doc page for a description of how various kinds of
bonus data must be formatted for certain atom styles.

E: Invalid atom ID in Bonus section of data file

Atom IDs must be positive integers and within range of defined
atoms.

E: Invalid atom ID in Bodies section of data file

Atom IDs must be positive integers and within range of defined
atoms.

E: Reuse of molecule template ID

The template IDs must be unique.

E: Atom sort did not operate correctly

This is an internal LAMMPS error.  Please report it to the
developers.

E: Too many atom sorting bins

This is likely due to an immense simulation box that has blown up
to a large size.

U: Cannot set mass for this atom style

This atom style does not support mass settings for each atom type.
Instead they are defined on a per-atom basis in the data file.

U: Invalid mass line in data file

Self-explanatory.

U: Invalid type for mass set

Mass command must set a type from 1-N where N is the number of atom
types.

U: Invalid mass value

Self-explanatory.

U: All masses are not set

For atom styles that define masses for each atom type, all masses must
be set in the data file or by the mass command before running a
simulation.  They must also be set before using the velocity
command.

*/
