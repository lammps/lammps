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

#ifndef LMP_ATOM_H
#define LMP_ATOM_H

#include "pointers.h"

#include <map>
#include <set>

namespace LAMMPS_NS {

// forward declarations

class AtomVec;
class Molecule;

class Atom : protected Pointers {
 public:
  char *atom_style;
  AtomVec *avec;
  enum { DOUBLE, INT, BIGINT };
  enum { GROW = 0, RESTART = 1, BORDER = 2 };
  enum { ATOMIC = 0, MOLECULAR = 1, TEMPLATE = 2 };
  enum { ATOM = 0, BOND = 1, ANGLE = 2, DIHEDRAL = 3, IMPROPER = 4 };
  enum { NUMERIC = 0, LABELS = 1 };
  enum { MAP_NONE = 0, MAP_ARRAY = 1, MAP_HASH = 2, MAP_YES = 3 };

  // atom counts

  bigint natoms;         // total # of atoms in system, could be 0
                         // natoms may not be current if atoms lost
  int nlocal, nghost;    // # of owned and ghost atoms on this proc
  int nmax;              // max # of owned+ghost in arrays on this proc
  int tag_enable;        // 0/1 if atom ID tags are defined
  int molecular;         // 0 = atomic, 1 = standard molecular system,
                         // 2 = molecule template system
  bigint nellipsoids;    // number of ellipsoids
  bigint nlines;         // number of lines
  bigint ntris;          // number of triangles
  bigint nbodies;        // number of bodies

  // system properties

  bigint nbonds, nangles, ndihedrals, nimpropers;
  int ntypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes;
  int bond_per_atom, angle_per_atom, dihedral_per_atom, improper_per_atom;
  int extra_bond_per_atom, extra_angle_per_atom;
  int extra_dihedral_per_atom, extra_improper_per_atom;

  int firstgroup;          // store atoms in this group first, -1 if unset
  int nfirst;              // # of atoms in first group on this proc
  char *firstgroupname;    // group-ID to store first, null pointer if unset

  // --------------------------------------------------------------------
  // 1st customization section: customize by adding new per-atom variable
  // per-atom vectors and arrays

  tagint *tag;
  int *type, *mask;
  imageint *image;
  double **x, **v, **f;

  // charged and dipolar particles

  double *rmass;
  double *q, **mu;

  // finite-size particles

  double *radius;
  double **omega, **angmom, **torque;
  int *ellipsoid, *line, *tri, *body;
  double **quat;
  double *temperature, *heatflow;

  // molecular systems

  tagint *molecule;
  int *molindex, *molatom;

  int **nspecial;      // 0,1,2 = cumulative # of 1-2,1-3,1-4 neighs
  tagint **special;    // IDs of 1-2,1-3,1-4 neighs of each atom
  int maxspecial;      // special[nlocal][maxspecial]

  int *num_bond;
  int **bond_type;
  tagint **bond_atom;

  int *num_angle;
  int **angle_type;
  tagint **angle_atom1, **angle_atom2, **angle_atom3;

  int *num_dihedral;
  int **dihedral_type;
  tagint **dihedral_atom1, **dihedral_atom2, **dihedral_atom3, **dihedral_atom4;

  int *num_improper;
  int **improper_type;
  tagint **improper_atom1, **improper_atom2, **improper_atom3, **improper_atom4;

  // PERI package

  double *vfrac, *s0;
  double **x0;

  // SPIN package

  double **sp, **fm, **fm_long;

  // EFF and AWPMD packages

  int *spin;
  double *eradius, *ervel, *erforce;
  double *ervelforce;
  double **cs, **csforce, **vforce;
  int *etag;

  // CG-DNA package

  tagint *id5p;

  // DPD-REACT package

  double *uCond, *uMech, *uChem, *uCGnew, *uCG;
  double *duChem;
  double *dpdTheta;
  int nspecies_dpd;

  // MESO package

  double **cc, **cc_flux;           // cc = chemical concentration
  double *edpd_temp, *edpd_flux;    // temperature and heat flux
  double *vest_temp;
  double *edpd_cv;    // heat capacity
  int cc_species;

  // MACHDYN package

  double *contact_radius;
  double **smd_data_9;
  double **smd_stress;
  double *eff_plastic_strain;
  double *eff_plastic_strain_rate;
  double *damage;

  // SPH package

  double *rho, *drho, *esph, *desph, *cv;
  double **vest;

  // AMOEBA package

  int *nspecial15;       // # of 1-5 neighs
  tagint **special15;    // IDs of 1-5 neighs of each atom
  int maxspecial15;      // special15[nlocal][maxspecial15]

  // DIELECTRIC package

  double *area, *ed, *em, *epsilon, *curvature, *q_scaled;

  // end of customization section
  // --------------------------------------------------------------------

  // --------------------------------------------------------------------
  // 2nd customization section: customize by adding new flags
  // identical list as Atom::set_atomflag_defaults()
  // most are existence flags for per-atom vectors and arrays
  // 1 if variable is used, 0 if not

  int labelmapflag, types_style;
  int sphere_flag, ellipsoid_flag, line_flag, tri_flag, body_flag;
  int peri_flag, electron_flag;
  int wavepacket_flag, sph_flag;

  int molecule_flag, molindex_flag, molatom_flag;
  int q_flag, mu_flag;
  int rmass_flag, radius_flag, omega_flag, torque_flag, angmom_flag, quat_flag;
  int temperature_flag, heatflow_flag;
  int vfrac_flag, spin_flag, eradius_flag, ervel_flag, erforce_flag;
  int cs_flag, csforce_flag, vforce_flag, ervelforce_flag, etag_flag;
  int rho_flag, esph_flag, cv_flag, vest_flag;
  int dpd_flag, edpd_flag, tdpd_flag;
  int mesont_flag;

  // SPIN package

  int sp_flag;

  // MACHDYN package

  int x0_flag;
  int smd_flag, damage_flag;
  int contact_radius_flag, smd_data_9_flag, smd_stress_flag;
  int eff_plastic_strain_flag, eff_plastic_strain_rate_flag;

  // AMOEBA package

  int nspecial15_flag;

  // Peridynamics scale factor, used by dump cfg

  double pdscale;

  // DIELECTRIC package

  int dielectric_flag;

  // end of customization section
  // --------------------------------------------------------------------

  // per-atom data struct describing all per-atom vectors/arrays

  struct PerAtom {
    std::string name;
    void *address;
    void *address_length;
    int *address_maxcols;
    int datatype;
    int cols;
    int collength;
    int threadflag;
  };

  std::vector<PerAtom> peratom;

  // custom vectors and arrays used by fix property/atom

  int **ivector, ***iarray;
  double **dvector, ***darray;
  int *icols, *dcols;
  char **ivname, **dvname, **ianame, **daname;
  int nivector, ndvector, niarray, ndarray;

  // molecule templates
  // each template can be a set of consecutive molecules
  // each with same ID (stored in molecules)
  // 1st molecule in template stores nset = # in set

  int nmolecule;
  Molecule **molecules;

  // type label maps

  class LabelMap *lmap;

  // extra peratom info in restart file destined for fix & diag

  double **extra;

  // per-type arrays

  double *mass;
  int *mass_setflag;

  // callback ptrs for atom arrays managed by fix classes

  int nextra_grow, nextra_restart, nextra_border;    // # of callbacks of each type
  int *extra_grow, *extra_restart, *extra_border;    // index of fix to callback to
  int nextra_grow_max, nextra_restart_max;           // size of callback lists
  int nextra_border_max;
  int nextra_store;

  int map_style;                    // style of atom map: 0=none, 1=array, 2=hash
  int map_user;                     // user requested map style:
                                    // 0 = no request, 1=array, 2=hash, 3=yes
  tagint map_tag_max;               // max atom ID that map() is setup for
  std::set<tagint> *unique_tags;    // set to ensure that bodies have unique tags

  // spatial sorting of atoms

  int sortfreq;          // sort atoms every this many steps, 0 = off
  bigint nextsort;       // next timestep to sort on
  double userbinsize;    // requested sort bin size

  // indices of atoms with same ID

  int *sametag;    // sametag[I] = next atom with same ID, -1 if no more

  // true if image flags were reset to 0 during data_atoms()

  bool reset_image_flag[3];

  // AtomVec factory types and map

  typedef AtomVec *(*AtomVecCreator)(LAMMPS *);
  typedef std::map<std::string, AtomVecCreator> AtomVecCreatorMap;
  AtomVecCreatorMap *avec_map;

  // --------------------------------------------------------------------
  // functions

  Atom(class LAMMPS *);
  ~Atom() override;

  void settings(class Atom *);
  void peratom_create();
  void add_peratom(const std::string &, void *, int, int, int threadflag = 0);
  void add_peratom_change_columns(const std::string &, int);
  void add_peratom_vary(const std::string &, void *, int, int *, void *, int collength = 0);
  void create_avec(const std::string &, int, char **, int);
  virtual AtomVec *new_avec(const std::string &, int, int &);

  virtual void init();
  void setup();

  std::string get_style();
  AtomVec *style_match(const char *);
  void modify_params(int, char **);
  void tag_check();
  void tag_extend();
  int tag_consecutive();

  void bonus_check();

  int parse_data(const char *);

  void deallocate_topology();

  void data_atoms(int, char *, tagint, tagint, int, int, double *, int, int *);
  void data_vels(int, char *, tagint);
  void data_bonds(int, char *, int *, tagint, int, int, int *);
  void data_angles(int, char *, int *, tagint, int, int, int *);
  void data_dihedrals(int, char *, int *, tagint, int, int, int *);
  void data_impropers(int, char *, int *, tagint, int, int, int *);
  void data_bonus(int, char *, AtomVec *, tagint);
  void data_bodies(int, char *, AtomVec *, tagint);
  void data_fix_compute_variable(int, int);

  virtual void allocate_type_arrays();
  void set_mass(const char *, int, const char *, int, int, int *);
  void set_mass(const char *, int, int, double);
  void set_mass(const char *, int, int, char **);
  void set_mass(double *);
  void check_mass(const char *, int);

  int radius_consistency(int, double &);
  int shape_consistency(int, double &, double &, double &);

  void add_molecule(int, char **);
  int find_molecule(const char *);
  std::vector<Molecule *> get_molecule_by_id(const std::string &);
  void add_molecule_atom(Molecule *, int, int, tagint);

  void add_label_map();

  void first_reorder();
  virtual void sort();

  void add_callback(int);
  void delete_callback(const char *, int);
  void update_callback(int);

  int find_custom(const char *, int &, int &);
  virtual int add_custom(const char *, int, int);
  virtual void remove_custom(int, int, int);

  virtual void sync_modify(ExecutionSpace, unsigned int, unsigned int) {}

  void *extract(const char *);
  int extract_datatype(const char *);

  inline int *get_map_array() { return map_array; };
  inline int get_map_size() { return map_tag_max + 1; };
  inline int get_max_same() { return max_same; };
  inline int get_map_maxarray() { return map_maxarray + 1; };

  // NOTE: placeholder method until KOKKOS/AtomVec is refactored
  int memcheck(const char *) { return 1; }

  double memory_usage();

  // functions for global to local ID mapping
  // map lookup function inlined for efficiency
  // return -1 if no map defined

  inline int map(tagint global)
  {
    if (map_style == 1)
      return map_array[global];
    else if (map_style == 2)
      return map_find_hash(global);
    else
      return -1;
  };

  virtual void map_init(int check = 1);
  virtual void map_clear();
  virtual void map_set();
  void map_one(tagint, int);
  int map_style_set();
  virtual void map_delete();
  int map_find_hash(tagint);

 protected:
  // global to local ID mapping

  int *map_array;      // direct map via array that holds map_tag_max
  int map_maxarray;    // allocated size of map_array (1 larger than this)

  struct HashElem {    // hashed map
    tagint global;     // key to search on = global ID
    int local;         // value associated with key = local index
    int next;          // next entry in this bucket, -1 if last
  };
  int map_nhash;         // # of entries hash table can hold
  int map_nused;         // # of actual entries in hash table
  int map_free;          // ptr to 1st unused entry in hash table
  int map_nbucket;       // # of hash buckets
  int *map_bucket;       // ptr to 1st entry in each bucket
  HashElem *map_hash;    // hash table

  int max_same;    // allocated size of sametag

  // spatial sorting of atoms

  int nbins;                           // # of sorting bins
  int nbinx, nbiny, nbinz;             // bins in each dimension
  int maxbin;                          // max # of bins
  int maxnext;                         // max size of next,permute
  int *binhead;                        // 1st atom in each bin
  int *next;                           // next atom in bin
  int *permute;                        // permutation vector
  double bininvx, bininvy, bininvz;    // inverse actual bin sizes
  double bboxlo[3], bboxhi[3];         // bounding box of my sub-domain

  void set_atomflag_defaults();
  void setup_sort_bins();
  int next_prime(int);
};

}    // namespace LAMMPS_NS

#endif
