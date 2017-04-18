/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2006 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: molfile_plugin.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.108 $       $Date: 2016/02/26 03:17:01 $
 *
 ***************************************************************************/

/** @file
 * API for C extensions to define a way to load structure, coordinate,
 * trajectory, and volumetric data files
 */

#ifndef MOL_FILE_PLUGIN_H
#define MOL_FILE_PLUGIN_H

#include "vmdplugin.h"

#if defined(DESRES_READ_TIMESTEP2)
/* includes needed for large integer types used for frame counts */
#include <sys/types.h>
typedef ssize_t molfile_ssize_t;      /**< for frame counts */
#endif

/**
 * Define a common plugin type to be used when registering the plugin.
 */
#define MOLFILE_PLUGIN_TYPE "mol file reader"

/**
 * File converter plugins use the same API  but register under a different
 * type so that regular file readers can have priority.
 */
#define MOLFILE_CONVERTER_PLUGIN_TYPE "mol file converter"

/* File plugin symbolic constants for better code readability */
#define MOLFILE_SUCCESS           0   /**< succeeded in reading file      */
#define MOLFILE_EOF              -1   /**< end of file                    */
#define MOLFILE_ERROR            -1   /**< error reading/opening a file   */
#define MOLFILE_NOSTRUCTUREDATA  -2   /**< no structure data in this file */

#define MOLFILE_NUMATOMS_UNKNOWN -1   /**< unknown number of atoms       */
#define MOLFILE_NUMATOMS_NONE     0   /**< no atoms in this file type    */

/**
 * Maximum string size macro
 */
#define MOLFILE_BUFSIZ           81   /**< maximum chars in string data  */
#define MOLFILE_BIGBUFSIZ        4096 /**< maximum chars in long strings */

#define MOLFILE_MAXWAVEPERTS     25   /**< maximum number of wavefunctions
                                       *   per timestep */

/**
 * Hard-coded direct-I/O page size constants for use by both VMD
 * and the plugins that want to use direct, unbuffered I/O for high
 * performance with SSDs etc.  We use two constants to define the
 * range of hardware page sizes that we can support, so that we can
 * add support for larger 8KB or 16KB page sizes in the future
 * as they become more prevalent in high-end storage systems.
 *
 * At present, VMD uses a hard-coded 4KB page size to reduce memory
 * fragmentation, but these constants will make it easier to enable the
 * use of larger page sizes in the future if it becomes necessary.
 */
#define MOLFILE_DIRECTIO_MIN_BLOCK_SIZE 4096
#define MOLFILE_DIRECTIO_MAX_BLOCK_SIZE 4096


/**
 * File level comments, origin information, and annotations.
 */
typedef struct {
  char database[81];   /**< database of origin, if any        */
  char accession[81];  /**< database accession code, if any   */
  char date[81];       /**< date/time stamp for this data     */
  char title[81];      /**< brief title for this data         */
  int remarklen;       /**< length of remarks string          */
  char *remarks;       /**< free-form remarks about data      */
} molfile_metadata_t;


/*
 * Struct for specifying atoms in a molecular structure.  The first
 * six components are required, the rest are optional and their presence is
 * indicating by setting the corresponding bit in optsflag.  When omitted,
 * the application (for read_structure) or plugin (for write_structure)
 * must be able to supply default values if the missing parameters are
 * part of its internal data structure.
 * Note that it is not possible to specify coordinates with this structure.
 * This is intentional; all coordinate I/O is done with the read_timestep and
 * write_timestep functions.
 */

/**
 * Per-atom attributes and information.
 */
typedef struct {
  /* these fields absolutely must be set or initialized to empty */
  char name[16];      /**< required atom name string             */
  char type[16];      /**< required atom type string             */
  char resname[8];    /**< required residue name string          */
  int resid;          /**< required integer residue ID           */
  char segid[8];      /**< required segment name string, or ""   */
#if 0 && vmdplugin_ABIVERSION > 17
  /* The new PDB file formats allows for much larger structures, */
  /* which can therefore require longer chain ID strings.  The   */
  /* new PDBx/mmCIF file formats do not have length limits on    */
  /* fields, so PDB chains could be arbitrarily long strings     */
  /* in such files.  At present, we know we need at least 3-char */
  /* chains for existing PDBx/mmCIF files.                       */
  char chain[4];      /**< required chain name, or ""            */
#else
  char chain[2];      /**< required chain name, or ""            */
#endif
  /* rest are optional; use optflags to specify what's present   */
  char altloc[2];     /**< optional PDB alternate location code  */
  char insertion[2];  /**< optional PDB insertion code           */
  float occupancy;    /**< optional occupancy value              */
  float bfactor;      /**< optional B-factor value               */
  float mass;         /**< optional mass value                   */
  float charge;       /**< optional charge value                 */
  float radius;       /**< optional radius value                 */
  int atomicnumber;   /**< optional element atomic number        */

#if 0
  char complex[16];
  char assembly[16];
  int qmregion;
  int qmregionlink;
  int qmlayer;
  int qmlayerlink;
  int qmfrag;
  int qmfraglink;
  string qmecp;
  int qmadapt;
  int qmect;          /**< boolean */
  int qmparam;
  int autoparam;
#endif

#if defined(DESRES_CTNUMBER)
  int ctnumber;       /**< mae ct block, 0-based, including meta */
#endif
} molfile_atom_t;

/*@{*/
/** Plugin optional data field availability flag */
#define MOLFILE_NOOPTIONS     0x0000 /**< no optional data                 */
#define MOLFILE_INSERTION     0x0001 /**< insertion codes provided         */
#define MOLFILE_OCCUPANCY     0x0002 /**< occupancy data provided          */
#define MOLFILE_BFACTOR       0x0004 /**< B-factor data provided           */
#define MOLFILE_MASS          0x0008 /**< Atomic mass provided             */
#define MOLFILE_CHARGE        0x0010 /**< Atomic charge provided           */
#define MOLFILE_RADIUS        0x0020 /**< Atomic VDW radius provided       */
#define MOLFILE_ALTLOC        0x0040 /**< Multiple conformations present   */
#define MOLFILE_ATOMICNUMBER  0x0080 /**< Atomic element number provided   */
#define MOLFILE_BONDSSPECIAL  0x0100 /**< Only non-standard bonds provided */
#if defined(DESRES_CTNUMBER)
#define MOLFILE_CTNUMBER      0x0200 /**< ctnumber provided */
#endif
#define MOLFILE_BADOPTIONS    0xFFFFFFFF /**< Detect badly behaved plugins */

/*@}*/

/*@{*/
/** Flags indicating availability of optional data fields
 *  for QM timesteps
 */
#define MOLFILE_QMTS_NOOPTIONS     0x0000 /**< no optional data               */
#define MOLFILE_QMTS_GRADIENT      0x0001 /**< energy gradients provided      */
#define MOLFILE_QMTS_SCFITER       0x0002
/*@}*/

typedef struct molfile_timestep_metadata {
  unsigned int count;                  /**< total # timesteps; -1 if unknown */
  unsigned int avg_bytes_per_timestep; /** bytes per timestep                */
  int has_velocities;                  /**< if timesteps have velocities     */
} molfile_timestep_metadata_t;

/*
 * Per-timestep atom coordinates and periodic cell information
 */
typedef struct {
  float *coords;        /**< coordinates of all atoms, arranged xyzxyzxyz   */
  float *velocities;    /**< space for velocities of all atoms; same layout */
                        /**< NULL unless has_velocities is set              */

  /*@{*/
  /**
   * Unit cell specification of the form A, B, C, alpha, beta, gamma.
   * notes: A, B, C are side lengths of the unit cell
   * alpha = angle between b and c
   *  beta = angle between a and c
   * gamma = angle between a and b
   */
  float A, B, C, alpha, beta, gamma;
  /*@}*/

  double physical_time; /**< physical time point associated with this frame */

#if defined(DESRES_READ_TIMESTEP2)
  /* HACK to support generic trajectory information */
  double total_energy;
  double potential_energy;
  double kinetic_energy;
  double extended_energy;
  double force_energy;
  double total_pressure;
#endif

} molfile_timestep_t;


/**
 * Metadata for volumetric datasets, read initially and used for subsequent
 * memory allocations and file loading.
 */
typedef struct {
  char dataname[256];   /**< name of volumetric data set                    */
  float origin[3];      /**< origin: origin of volume (x=0, y=0, z=0 corner */

  /*
   * x/y/z axis:
   * These the three cell sides, providing both direction and length
   * (not unit vectors) for the x, y, and z axes.  In the simplest
   * case, these would be <size,0,0> <0,size,0> and <0,0,size) for
   * an orthogonal cubic volume set.  For other cell shapes these
   * axes can be oriented non-orthogonally, and the parallelpiped
   * may have different side lengths, not just a cube/rhombus.
   */
  float xaxis[3];       /**< direction (and length) for X axis              */
  float yaxis[3];       /**< direction (and length) for Y axis              */
  float zaxis[3];       /**< direction (and length) for Z axis              */

  /*
   * x/y/z size:
   * Number of grid cells along each axis.  This is _not_ the
   * physical size of the box, this is the number of voxels in each
   * direction, independent of the shape of the volume set.
   */
  int xsize;            /**< number of grid cells along the X axis           */
  int ysize;            /**< number of grid cells along the Y axis           */
  int zsize;            /**< number of grid cells along the Z axis           */

#if vmdplugin_ABIVERSION > 16
  int has_scalar;       /**< flag indicating presence of scalar volume       */
  int has_gradient;     /**< flag indicating presence of vector volume       */
  int has_variance;     /**< flag indicating presence of variance map        */
#endif
  int has_color;        /**< flag indicating presence of voxel color data    */
} molfile_volumetric_t;


#if vmdplugin_ABIVERSION > 16
/**
 * Volumetric dataset read/write structure with both flag/parameter sets
 * and VMD-allocated pointers for fields to be used by the plugin.
 */
typedef struct {
  int setidx;           /**< volumetric dataset index to load/save */
  float *scalar;        /**< scalar density/potential field data   */
  float *gradient;      /**< gradient vector field                 */
  float *variance;      /**< variance map indicating signal/noise  */
  float *rgb3f;         /**< RGB floating point color texture map  */
  unsigned char *rgb3u; /**< RGB unsigned byte color texture map   */
} molfile_volumetric_readwrite_t;
#endif


/**************************************************************
 **************************************************************
 ****                                                      ****
 ****          Data structures for QM files                ****
 ****                                                      ****
 **************************************************************
 **************************************************************/

/* macros for the convergence status of a QM calculation. */
#define MOLFILE_QMSTATUS_UNKNOWN       -1 /* don't know yet */
#define MOLFILE_QMSTATUS_OPT_CONV       0 /* optimization converged */
#define MOLFILE_QMSTATUS_SCF_NOT_CONV   1 /* SCF convergence failed */
#define MOLFILE_QMSTATUS_OPT_NOT_CONV   2 /* optimization not converged */
#define MOLFILE_QMSTATUS_FILE_TRUNCATED 3 /* file was truncated */

/* macros describing the SCF method (SCFTYP in GAMESS) */
#define MOLFILE_SCFTYPE_UNKNOWN -1 /* no info about the method  */
#define MOLFILE_SCFTYPE_NONE     0 /* calculation didn't make use of SCF */
#define MOLFILE_SCFTYPE_RHF      1 /* restricted Hartree-Fock   */
#define MOLFILE_SCFTYPE_UHF      2 /* unrestricted Hartree-Fock */
#define MOLFILE_SCFTYPE_ROHF     3 /* restricted open-shell Hartree-Fock */
#define MOLFILE_SCFTYPE_GVB      4 /* generalized valence bond orbitals  */
#define MOLFILE_SCFTYPE_MCSCF    5 /* multi-configuration SCF   */
#define MOLFILE_SCFTYPE_FF       6 /* classical force-field based sim.   */

/* macros describing the type of calculation (RUNTYP in GAMESS) */
#define MOLFILE_RUNTYPE_UNKNOWN    0  /* single point run */
#define MOLFILE_RUNTYPE_ENERGY     1  /* single point run */
#define MOLFILE_RUNTYPE_OPTIMIZE   2  /* geometry optimization */
#define MOLFILE_RUNTYPE_SADPOINT   3  /* saddle point search */
#define MOLFILE_RUNTYPE_HESSIAN    4  /* Hessian/frequency calculation */
#define MOLFILE_RUNTYPE_SURFACE    5  /* potential surface scan */
#define MOLFILE_RUNTYPE_GRADIENT   6  /* energy gradient calculation */
#define MOLFILE_RUNTYPE_MEX        7  /* minimum energy crossing */
#define MOLFILE_RUNTYPE_DYNAMICS   8  /* Any type of molecular dynamics
                                       * e.g. Born-Oppenheimer, Car-Parinello,
                                       * or classical MD */
#define MOLFILE_RUNTYPE_PROPERTIES 9  /* Properties were calculated from a
                                       * wavefunction that was read from file */


/**
 * Sizes of various QM-related, timestep independent data arrays
 * which must be allocated by the caller (VMD) so that the plugin
 * can fill in the arrays with data.
 */
typedef struct {
  /* hessian data */
  int nimag;                    /**< number of imaginary modes */
  int nintcoords;               /**< number internal coordinates */
  int ncart;                    /**< number cartesian coordinates */

  /* orbital/basisset data */
  int num_basis_funcs;          /**< number of uncontracted basis functions in basis array */
  int num_basis_atoms;          /**< number of atoms in basis set */
  int num_shells;               /**< total number of atomic shells */
  int wavef_size;               /**< size of the wavefunction
                                 *   i.e. size of secular eq. or
                                 *   # of cartesian contracted
                                 *   gaussian basis functions */

  /* everything else */
  int have_sysinfo;
  int have_carthessian;         /**< hessian in cartesian coords available  */
  int have_inthessian;          /**< hessian in internal coords available  */
  int have_normalmodes;         /**< normal modes available  */
} molfile_qm_metadata_t;


/**
 * QM run info. Parameters that stay unchanged during a single file.
 */
typedef struct {
  int nproc;             /**< number of processors used. */
  int memory;            /**< amount of memory used in Mbyte. */
  int runtype;           /**< flag indicating the calculation method. */
  int scftype;           /**< SCF type: RHF, UHF, ROHF, GVB or MCSCF wfn. */
  int status;            /**< indicates wether SCF and geometry optimization
                          *   have converged properly. */
  int num_electrons;     /**< number of electrons.    XXX: can be fractional in some DFT codes */
  int totalcharge;       /**< total charge of system. XXX: can be fractional in some DFT codes */
  int num_occupied_A;    /**< number of occupied alpha orbitals */
  int num_occupied_B;    /**< number of occupied beta orbitals */

  double *nuc_charge;    /**< array(natom) containing the nuclear charge of atom i */

  char basis_string[MOLFILE_BUFSIZ];    /**< basis name as "nice" string. */
  char runtitle[MOLFILE_BIGBUFSIZ];     /**< title of run.                */
  char geometry[MOLFILE_BUFSIZ];        /**< type of provided geometry,   XXX: remove?
                                         * e.g. UNIQUE, ZMT, CART, ...    */
  char version_string[MOLFILE_BUFSIZ];  /**< QM code version information. */
} molfile_qm_sysinfo_t;


/**
 * Data for QM basis set
 */
typedef struct {
  int *num_shells_per_atom; /**< number of shells per atom */
  int *num_prim_per_shell;  /**< number of shell primitives shell */

  float *basis;             /**< contraction coeffients and exponents for
                             *   the basis functions in the form
                             *   {exp(1), c-coeff(1), exp(2), c-coeff(2), ...};
                             *   array size = 2*num_basis_funcs
                             *   The basis must NOT be normalized. */
  int *atomic_number;       /**< atomic numbers (chem. element) of atoms in basis set */
  int *angular_momentum;    /**< 3 ints per wave function coefficient do describe the
                             *   cartesian components of the angular momentum.
                             *   E.g. S={0 0 0}, Px={1 0 0}, Dxy={1 1 0}, or Fyyz={0 2 1}.
                             */
  int *shell_types;         /**< type for each shell in basis */
} molfile_qm_basis_t;


/**
 * Data from QM Hessian/normal mode runs
 *
 * A noteworthy comment from one of Axel's emails:
 * The molfile_qm_hessian_t, I'd rename to molfile_hessian_t (one
 * can do vibrational analysis without QM) and would make this a
 * completely separate entity. This could then be also used to
 * read in data from, say, principal component analysis or normal
 * mode analysis and VMD could contain code to either project a
 * trajectory on the contained eigenvectors or animate them and
 * so on. There is a bunch of possible applications...
 */
typedef struct {
  double *carthessian;  /**< hessian matrix in cartesian coordinates (ncart)*(ncart)
                         *   as a single array of doubles (row(1), ...,row(natoms)) */
  int    *imag_modes;   /**< list(nimag) of imaginary modes */
  double *inthessian;   /**< hessian matrix in internal coordinates
                         *   (nintcoords*nintcoords) as a single array of
                         *   doubles (row(1), ...,row(nintcoords)) */
  float *wavenumbers;   /**< array(ncart) of wavenumbers of normal modes */
  float *intensities;   /**< array(ncart) of intensities of normal modes */
  float *normalmodes;   /**< matrix(ncart*ncart) of normal modes  */
} molfile_qm_hessian_t;


/**
 * QM related information that is timestep independent
 */
typedef struct {
  molfile_qm_sysinfo_t run;             /* system info  */
  molfile_qm_basis_t   basis;           /* basis set info */
  molfile_qm_hessian_t hess;            /* hessian info */
} molfile_qm_t;



/**
 *  Enumeration of all of the wavefunction types that can be read
 *  from QM file reader plugins.
 *
 *  CANON    = canonical (i.e diagonalized) wavefunction
 *  GEMINAL  = GVB-ROHF geminal pairs
 *  MCSCFNAT = Multi-Configuration SCF natural orbitals
 *  MCSCFOPT = Multi-Configuration SCF optimized orbitals
 *  CINATUR  = Configuration-Interaction natural orbitals
 *  BOYS     = Boys localization
 *  RUEDEN   = Ruedenberg localization
 *  PIPEK    = Pipek-Mezey population localization
 *
 *  NBO related localizations:
 *  --------------------------
 *  NAO      = Natural Atomic Orbitals
 *  PNAO     = pre-orthogonal NAOs
 *  NBO      = Natural Bond Orbitals
 *  PNBO     = pre-orthogonal NBOs
 *  NHO      = Natural Hybrid Orbitals
 *  PNHO     = pre-orthogonal NHOs
 *  NLMO     = Natural Localized Molecular Orbitals
 *  PNLMO    = pre-orthogonal NLMOs
 *
 *  UNKNOWN  = Use this for any type not listed here
 *             You can use the string field for description
 */
enum molfile_qm_wavefunc_type {
  MOLFILE_WAVE_CANON,    MOLFILE_WAVE_GEMINAL,
  MOLFILE_WAVE_MCSCFNAT, MOLFILE_WAVE_MCSCFOPT,
  MOLFILE_WAVE_CINATUR,
  MOLFILE_WAVE_PIPEK,  MOLFILE_WAVE_BOYS, MOLFILE_WAVE_RUEDEN,
  MOLFILE_WAVE_NAO,    MOLFILE_WAVE_PNAO, MOLFILE_WAVE_NHO,
  MOLFILE_WAVE_PNHO,   MOLFILE_WAVE_NBO,  MOLFILE_WAVE_PNBO,
  MOLFILE_WAVE_PNLMO,  MOLFILE_WAVE_NLMO, MOLFILE_WAVE_MOAO,
  MOLFILE_WAVE_NATO,   MOLFILE_WAVE_UNKNOWN
};


/**
 *  Enumeration of all of the supported QM related charge
 *  types
 */
enum molfile_qm_charge_type {
  MOLFILE_QMCHARGE_UNKNOWN,
  MOLFILE_QMCHARGE_MULLIKEN, MOLFILE_QMCHARGE_LOWDIN,
  MOLFILE_QMCHARGE_ESP, MOLFILE_QMCHARGE_NPA
};



/**
 * Sizes of various QM-related, per-timestep data arrays
 * which must be allocated by the caller (VMD) so that the plugin
 * can fill in the arrays with data.
 */
typedef struct molfile_qm_timestep_metadata {
  unsigned int count;                  /**< total # timesteps; -1 if unknown */
  unsigned int avg_bytes_per_timestep; /**< bytes per timestep                */
  int has_gradient;                    /**< if timestep contains gradient    */
  int num_scfiter;                     /**< # scf iterations for this ts     */
  int num_orbitals_per_wavef[MOLFILE_MAXWAVEPERTS]; /**< # orbitals for each wavefunction */
  int has_orben_per_wavef[MOLFILE_MAXWAVEPERTS]; /**< orbital energy flags */
  int has_occup_per_wavef[MOLFILE_MAXWAVEPERTS]; /**< orbital occupancy flags */
  int num_wavef ;                      /**< # wavefunctions in this ts     */
  int wavef_size;                      /**< size of one wavefunction
                                        *   (# of gaussian basis fctns)    */
  int num_charge_sets;                 /**< # of charge values per atom */
} molfile_qm_timestep_metadata_t;


/**
 * QM wavefunction
 */
typedef struct {
  int   type;               /**< MOLFILE_WAVE_CANON, MOLFILE_WAVE_MCSCFNAT, ... */
  int   spin;               /**< 1 for alpha, -1 for beta */
  int   excitation;         /**< 0 for ground state, 1,2,3,... for excited states */
  int   multiplicity;       /**< spin multiplicity of the state, zero if unknown */
  char info[MOLFILE_BUFSIZ]; /**< string for additional type info */

  double energy;            /**< energy of the electronic state.
                             *   i.e. HF-SCF energy, CI state energy,
                             *   MCSCF energy, etc. */

  float *wave_coeffs;       /**< expansion coefficients for wavefunction in the
                             *   form {orbital1(c1),orbital1(c2),.....,orbitalM(cN)} */
  float *orbital_energies;  /**< list of orbital energies for wavefunction */
  float *occupancies;       /**< orbital occupancies */
  int   *orbital_ids;       /**< orbital ID numbers; If NULL then VMD will
                             *   assume 1,2,3,...num_orbs.     */
} molfile_qm_wavefunction_t;


/**
 * QM per trajectory timestep info
 * Note that each timestep can contain multiple wavefunctions.
 */
typedef struct {
  molfile_qm_wavefunction_t *wave; /**< array of wavefunction objects */
  float  *gradient;         /**< force on each atom (=gradient of energy) */

  double *scfenergies;      /**< energies from the SCF cycles */
  double *charges;          /**< per-atom charges */
  int    *charge_types;     /**< type of each charge set */
} molfile_qm_timestep_t;


/**************************************************************
 **************************************************************/




/**
 *  Enumeration of all of the supported graphics objects that can be read
 *  from graphics file reader plugins.
 */
enum molfile_graphics_type {
  MOLFILE_POINT,  MOLFILE_TRIANGLE, MOLFILE_TRINORM, MOLFILE_NORMS,
  MOLFILE_LINE,   MOLFILE_CYLINDER, MOLFILE_CAPCYL,  MOLFILE_CONE,
  MOLFILE_SPHERE, MOLFILE_TEXT,     MOLFILE_COLOR,   MOLFILE_TRICOLOR
};

/**
 *  Individual graphics object/element data
 */
typedef struct {
  int type;             /* One of molfile_graphics_type */
  int style;            /* A general style parameter    */
  float size;           /* A general size parameter     */
  float data[9];        /* All data for the element     */
} molfile_graphics_t;


/*
 * Types for raw graphics elements stored in files.  Data for each type
 * should be stored by the plugin as follows:

type        data                                     style       size
----        ----                                     -----       ----
point       x, y, z                                              pixel size
triangle    x1,y1,z1,x2,y2,z2,x3,y3,z3
trinorm     x1,y1,z1,x2,y2,z2,x3,y3,z3
            the next array element must be NORMS
tricolor    x1,y1,z1,x2,y2,z2,x3,y3,z3
            the next array elements must be NORMS
            the following element must be COLOR, with three RGB triples
norms       x1,y1,z1,x2,y2,z2,x3,y3,z3
line        x1,y1,z1,x2,y2,z2                        0=solid     pixel width
                                                     1=stippled
cylinder    x1,y1,z1,x2,y2,z2                        resolution  radius
capcyl      x1,y1,z1,x2,y2,z2                        resolution  radius
sphere      x1,y1,z1                                 resolution  radius
text        x, y, z, up to 24 bytes of text                      pixel size
color       r, g, b
*/


/**
 * Main file reader API.  Any function in this struct may be NULL
 * if not implemented by the plugin; the application checks this to determine
 * what functionality is present in the plugin.
 */
typedef struct {
  /**
   * Required header
   */
  vmdplugin_HEAD

  /**
   * Filename extension for this file type.  May be NULL if no filename
   * extension exists and/or is known.  For file types that match several
   * common extensions, list them in a comma separated list such as:
   *  "pdb,ent,foo,bar,baz,ban"
   * The comma separated list will be expanded when filename extension matching
   * is performed.  If multiple plugins solicit the same filename extensions,
   * the one that lists the extension earliest in its list is selected. In the
   * case of a "tie", the first one tried/checked "wins".
   */
  const char *filename_extension;

  /**
   * Try to open the file for reading.  Return an opaque handle, or NULL on
   * failure. Set the number of atoms; if the number of atoms cannot be
   * determined, set natoms to MOLFILE_NUMATOMS_UNKNOWN.
   * Filetype should be the name under which this plugin was registered;
   * this is provided so that plugins can provide the same function pointer
   * to handle multiple file types.
   */
  void *(* open_file_read)(const char *filepath, const char *filetype,
      int *natoms);

  /**
   * Read molecular structure from the given file handle.  atoms is allocated
   * by the caller and points to space for natoms.
   * On success, place atom information in the passed-in pointer.
   * optflags specifies which optional fields in the atoms will be set by
   * the plugin.
   */
  int (*read_structure)(void *, int *optflags, molfile_atom_t *atoms);

  /**
   * Read bond information for the molecule.  On success the arrays from
   * and to should point to the (one-based) indices of bonded atoms.
   * Each unique bond should be specified only once, so file formats that list
   * bonds twice will need post-processing before the results are returned to
   * the caller.
   * If the plugin provides bond information, but the file loaded doesn't
   * actually contain any bond info, the nbonds parameter should be
   * set to 0 and from/to should be set to NULL to indicate that no bond
   * information was actually present, and automatic bond search should be
   * performed.
   *
   * If the plugin provides bond order information, the bondorder array
   * will contain the bond order for each from/to pair.  If not, the bondorder
   * pointer should be set to NULL, in which case the caller will provide a
   * default bond order value of 1.0.
   *
   * If the plugin provides bond type information, the bondtype array
   * will contain the bond type index for each from/to pair. These numbers
   * are consecutive integers starting from 0.
   * the bondtypenames list, contains the corresponding names, if available,
   * as a NULL string terminated list. nbondtypes is provided for convenience
   * and consistency checking.
   *
   * These arrays must be freed by the plugin in the close_file_read function.
   * This function can be called only after read_structure().
   * Return MOLFILE_SUCCESS if no errors occur.
   */
  int (*read_bonds)(void *, int *nbonds, int **from, int **to, float **bondorder,
                    int **bondtype, int *nbondtypes, char ***bondtypename);

  /**
   * XXX this function will be augmented and possibly superceded by a
   *     new QM-capable version named read_timestep(), when finished.
   *
   * Read the next timestep from the file.  Return MOLFILE_SUCCESS, or
   * MOLFILE_EOF on EOF.  If the molfile_timestep_t argument is NULL, then
   * the frame should be skipped.  Otherwise, the application must prepare
   * molfile_timestep_t by allocating space in coords for the corresponding
   * number of coordinates.
   * The natoms parameter exists because some coordinate file formats
   * (like CRD) cannot determine for themselves how many atoms are in a
   * timestep; the app must therefore obtain this information elsewhere
   * and provide it to the plugin.
   */
  int (* read_next_timestep)(void *, int natoms, molfile_timestep_t *);

  /**
   * Close the file and release all data.  The handle cannot be reused.
   */
  void (* close_file_read)(void *);

  /**
   * Open a coordinate file for writing using the given header information.
   * Return an opaque handle, or NULL on failure.  The application must
   * specify the number of atoms to be written.
   * filetype should be the name under which this plugin was registered.
   */
  void *(* open_file_write)(const char *filepath, const char *filetype,
      int natoms);

  /**
   * Write structure information.  Return success.
   */
  int (* write_structure)(void *, int optflags, const molfile_atom_t *atoms);

  /**
   * Write a timestep to the coordinate file.  Return MOLFILE_SUCCESS if no
   * errors occur.  If the file contains structure information in each
   * timestep (like a multi-entry PDB), it will have to cache the information
   * from the initial calls from write_structure.
   */
  int (* write_timestep)(void *, const molfile_timestep_t *);

  /**
   * Close the file and release all data.  The handle cannot be reused.
   */
  void (* close_file_write)(void *);

  /**
   * Retrieve metadata pertaining to volumetric datasets in this file.
   * Set nsets to the number of volumetric data sets, and set *metadata
   * to point to an array of molfile_volumetric_t.  The array is owned by
   * the plugin and should be freed by close_file_read().  The application
   * may call this function any number of times.
   */
  int (* read_volumetric_metadata)(void *, int *nsets,
        molfile_volumetric_t **metadata);

  /**
   * Read the specified volumetric data set into the space pointed to by
   * datablock.  The set is specified with a zero-based index.  The space
   * allocated for the datablock must be equal to
   * xsize * ysize * zsize.  No space will be allocated for colorblock
   * unless has_color is nonzero; in that case, colorblock should be
   * filled in with three RGB floats per datapoint.
   */
  int (* read_volumetric_data)(void *, int set, float *datablock,
        float *colorblock);
#if vmdplugin_ABIVERSION > 16
  int (* read_volumetric_data_ex)(void *, molfile_volumetric_readwrite_t *v);
#endif

  /**
   * Read raw graphics data stored in this file.   Return the number of data
   * elements and the data itself as an array of molfile_graphics_t in the
   * pointer provided by the application.  The plugin is responsible for
   * freeing the data when the file is closed.
   */
  int (* read_rawgraphics)(void *, int *nelem, const molfile_graphics_t **data);

  /**
   * Read molecule metadata such as what database (if any) this file/data
   * came from, what the accession code for the database is, textual remarks
   * and other notes pertaining to the contained structure/trajectory/volume
   * and anything else that's informative at the whole file level.
   */
  int (* read_molecule_metadata)(void *, molfile_metadata_t **metadata);

  /**
   * Write bond information for the molecule.  The arrays from
   * and to point to the (one-based) indices of bonded atoms.
   * Each unique bond will be specified only once by the caller.
   * File formats that list bonds twice will need to emit both the
   * from/to and to/from versions of each.
   * This function must be called before write_structure().
   *
   * Like the read_bonds() routine, the bondorder pointer is set to NULL
   * if the caller doesn't have such information, in which case the
   * plugin should assume a bond order of 1.0 if the file format requires
   * bond order information.
   *
   * Support for bond types follows the bondorder rules. bondtype is
   * an integer array of the size nbonds that contains the bond type
   * index (consecutive integers starting from 0) and bondtypenames
   * contain the corresponding strings, in case the naming/numbering
   * scheme is different from the index numbers.
   * if the pointers are set to NULL, then this information is not available.
   * bondtypenames can only be used of bondtypes is also given.
   * Return MOLFILE_SUCCESS if no errors occur.
   */
  int (* write_bonds)(void *, int nbonds, int *from, int *to, float *bondorder,
                     int *bondtype, int nbondtypes, char **bondtypename);

  /**
   * Write the specified volumetric data set into the space pointed to by
   * datablock.  The * allocated for the datablock must be equal to
   * xsize * ysize * zsize.  No space will be allocated for colorblock
   * unless has_color is nonzero; in that case, colorblock should be
   * filled in with three RGB floats per datapoint.
   */
  int (* write_volumetric_data)(void *, molfile_volumetric_t *metadata,
                                float *datablock, float *colorblock);
#if vmdplugin_ABIVERSION > 16
  int (* write_volumetric_data_ex)(void *, molfile_volumetric_t *metadata,
                                   molfile_volumetric_readwrite_t *v);
#endif

  /**
   * Read in Angles, Dihedrals, Impropers, and Cross Terms and optionally types.
   * (Cross terms pertain to the CHARMM/NAMD CMAP feature)
   */
  int (* read_angles)(void *handle, int *numangles, int **angles, int **angletypes,
                      int *numangletypes, char ***angletypenames, int *numdihedrals,
                      int **dihedrals, int **dihedraltypes, int *numdihedraltypes,
                      char ***dihedraltypenames, int *numimpropers, int **impropers,
                      int **impropertypes, int *numimpropertypes, char ***impropertypenames,
                      int *numcterms, int **cterms, int *ctermcols, int *ctermrows);

  /**
   * Write out Angles, Dihedrals, Impropers, and Cross Terms
   * (Cross terms pertain to the CHARMM/NAMD CMAP feature)
   */
  int (* write_angles)(void *handle, int numangles, const int *angles, const int *angletypes,
                       int numangletypes, const char **angletypenames, int numdihedrals,
                       const int *dihedrals, const int *dihedraltypes, int numdihedraltypes,
                       const char **dihedraltypenames, int numimpropers,
                       const int *impropers, const int *impropertypes, int numimpropertypes,
                       const char **impropertypenames, int numcterms,  const int *cterms,
                       int ctermcols, int ctermrows);


  /**
   * Retrieve metadata pertaining to timestep independent
   * QM datasets in this file.
   *
   * The metadata are the sizes of the QM related data structure
   * arrays that will be populated by the plugin when
   * read_qm_rundata() is called. Since the allocation of these
   * arrays is done by VMD rather than the plugin, VMD needs to
   * know the sizes beforehand. Consequently read_qm_metadata()
   * has to be called before read_qm_rundata().
   */
  int (* read_qm_metadata)(void *, molfile_qm_metadata_t *metadata);


  /**
   * Read timestep independent QM data.
   *
   * Typical data that are defined only once per trajectory are
   * general info about the calculation (such as the used method),
   * the basis set and normal modes.
   * The data structures to be populated must have been allocated
   * before by VMD according to sizes obtained through
   * read_qm_metadata().
   */
  int (* read_qm_rundata)(void *, molfile_qm_t *qmdata);


  /**
   * Read the next timestep from the file.  Return MOLFILE_SUCCESS, or
   * MOLFILE_EOF on EOF.  If the molfile_timestep_t or molfile_qm_metadata_t
   * arguments are NULL, then the coordinate or qm data should be skipped.
   * Otherwise, the application must prepare molfile_timestep_t and
   * molfile_qm_timestep_t by allocating space for the corresponding
   * number of coordinates, orbital wavefunction coefficients, etc.
   * Since it is common for users to want to load only the final timestep
   * data from a QM run, the application may provide any combination of
   * valid, or NULL pointers for the molfile_timestep_t and
   * molfile_qm_timestep_t parameters, depending on what information the
   * user is interested in.
   * The natoms and qm metadata parameters exist because some file formats
   * cannot determine for themselves how many atoms etc are in a
   * timestep; the app must therefore obtain this information elsewhere
   * and provide it to the plugin.
   */
  int (* read_timestep)(void *, int natoms, molfile_timestep_t *,
                        molfile_qm_metadata_t *, molfile_qm_timestep_t *);

  int (* read_timestep_metadata)(void *, molfile_timestep_metadata_t *);
  int (* read_qm_timestep_metadata)(void *, molfile_qm_timestep_metadata_t *);

#if defined(DESRES_READ_TIMESTEP2)
  /**
    * Read a specified timestep!
    */
  int (* read_timestep2)(void *, molfile_ssize_t index, molfile_timestep_t *);

  /**
    * write up to count times beginning at index start into the given
    * space.  Return the number read, or -1 on error.
    */
  molfile_ssize_t (* read_times)( void *,
                                  molfile_ssize_t start,
                                  molfile_ssize_t count,
                                  double * times );
#endif

  /**
   *  Console output, READ-ONLY function pointer.
   *  Function pointer that plugins can use for printing to the host
   *  application's text console.  This provides a clean way for plugins
   *  to send message strings back to the calling application, giving the
   *  caller the ability to prioritize, buffer, and redirect console messages
   *  to an appropriate output channel, window, etc.  This enables the use of
   *  graphical consoles like TkCon without losing console output from plugins.
   *  If the function pointer is NULL, no console output service is provided
   *  by the calling application, and the output should default to stdout
   *  stream.  If the function pointer is non-NULL, all output will be
   *  subsequently dealt with by the calling application.
   *
   *  XXX this should really be put into a separate block of
   *      application-provided read-only function pointers for any
   *      application-provided services
   */
  int (* cons_fputs)(const int, const char*);

} molfile_plugin_t;

#endif

