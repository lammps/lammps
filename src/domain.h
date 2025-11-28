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

#ifndef LMP_DOMAIN_H
#define LMP_DOMAIN_H

#include "pointers.h"

#include <cmath>
#include <map>
#include <unordered_set>

namespace LAMMPS_NS {
class Region;

/** \class LAMMPS_NS::Domain
 * \brief Simulation box and domain decomposition management for LAMMPS
 *
 * The Domain class manages the simulation box geometry and properties in LAMMPS.
 * It handles both orthogonal and triclinic (tilted) simulation boxes, and provides
 * the geometric framework for domain decomposition across MPI processes.
 *
 * Key responsibilities include:
 *
 * - Defining the global simulation box dimensions and boundaries
 * - Managing periodic, fixed, and shrink-wrap boundary conditions
 * - Handling triclinic boxes (both restricted and general forms)
 * - Coordinate transformations between Cartesian and lamda (fractional) coordinates
 * - Minimum image convention for periodic systems
 * - Image flag manipulation for atoms crossing periodic boundaries
 * - Managing simulation regions (geometric shapes for selections)
 * - Managing the lattice for atom creation
 *
 * The class supports both 2D and 3D simulations and provides methods for
 * coordinate remapping, unwrapping, and periodic image calculations.
 *
 * \sa LAMMPS_NS::Region, LAMMPS_NS::Lattice, LAMMPS_NS::Pointers */
class Domain : protected Pointers {
 public:
  int box_exist;                          /**< 0 = box not yet created, 1 = box exists */
  int dimension;                          /**< Simulation dimensionality: 2 = 2D, 3 = 3D */
  int nonperiodic;                        /**< 0 = periodic all dims, 1 = periodic/fixed, 2 = shrink-wrap */
  int xperiodic, yperiodic, zperiodic;    /**< Periodicity flags: 0 = non-periodic, 1 = periodic */
  int periodicity[3];                     /**< Periodicity flags as array [x,y,z] */

  int boundary[3][2];    /**< Boundary settings for 6 faces: 0=periodic, 1=fixed, 2=shrink-wrap, 3=shrink-wrap/min */

  int triclinic;            /**< 0 = orthogonal box, 1 = triclinic (tilted) box */
  int triclinic_general;    /**< 1 if general triclinic mapping stored, 0 if not */

  // orthogonal box

  double xprd, yprd, zprd;                   /**< Global box dimensions in x, y, z */
  double xprd_half, yprd_half, zprd_half;    /**< Half of global box dimensions */
  double prd[3];                             /**< Global box dimensions as array */
  double prd_half[3];                        /**< Half dimensions as array */

  // restricted triclinic box
  // xyz prd,xyz prd_half and prd,prd_half = same as if not tilted

  double prd_lamda[3];         /**< Box dimensions in lamda coords = (1,1,1) */
  double prd_half_lamda[3];    /**< Half box in lamda coords = (0.5,0.5,0.5) */

  // orthogonal box global bounds

  double boxlo[3], boxhi[3];    /**< Global box lower and upper bounds */

  // restricted triclinic box
  // boxlo/hi = same as if not tilted

  double boxlo_lamda[3], boxhi_lamda[3];    /**< Lamda box bounds = (0,1) */
  double boxlo_bound[3], boxhi_bound[3];    /**< Bounding box of tilted domain */
  double corners[8][3];                     /**< 8 corner points of the box */

  // orthogonal box & restricted triclinic box

  double minxlo, minxhi;    /**< Minimum x bounds when shrink-wrapping */
  double minylo, minyhi;    /**< Minimum y bounds when shrink-wrapping */
  double minzlo, minzhi;    /**< Minimum z bounds when shrink-wrapping */

  // orthogonal box

  double sublo[3], subhi[3];    /**< This processor's subdomain bounds */

  // restricted triclinic box
  // sublo/hi = undefined

  double sublo_lamda[3], subhi_lamda[3];    /**< Subdomain bounds in lamda coords */

  // restricted triclinic box

  double xy, xz, yz;                /**< 3 tilt factors for triclinic box */
  double h[6], h_inv[6];            /**< Shape matrix and inverse in Voigt ordering (xx,yy,zz,yz,xz,xy) */
  double h_rate[6], h_ratelo[3];    /**< Rate of box size/shape change */

  // general triclinic box
  // boxlo = lower left corner

  double avec[3], bvec[3], cvec[3];    /**< ABC edge vectors of general triclinic box */
  double rotate_g2r[3][3];             /**< Rotation matrix: general to restricted triclinic */
  double rotate_r2g[3][3];             /**< Rotation matrix: restricted to general triclinic */

  // box flags

  int box_change;           /**< 1 if any box change flag is set, 0 otherwise */
  int box_change_size;      /**< 1 if box size changes during run, 0 if not */
  int box_change_shape;     /**< 1 if box shape changes during run, 0 if not */
  int box_change_domain;    /**< 1 if processor subdomains change, 0 if not */

  int deform_flag;        /**< 1 if fix deform exists, 0 otherwise */
  int deform_vremap;      /**< 1 if fix deform remaps velocities, 0 otherwise */
  int deform_groupbit;    /**< Atom group bitmask for velocity remapping */

  class Lattice *lattice;    /**< Pointer to user-defined lattice */

  int copymode;    /**< 1 if this is a copy, don't deallocate in destructor */

  /* Velocity remapping modes */
  enum { NO_REMAP, X_REMAP, V_REMAP };

  using RegionCreator = Region *(*) (LAMMPS *, int, char **);    /**< Function pointer type for region creation */
  using RegionCreatorMap = std::map<std::string, RegionCreator>; /**< Map type for region style factory */
  RegionCreatorMap *region_map;    /**< Factory map for region styles */

  /** Constructor for Domain class
   * \param lmp Pointer to the main LAMMPS instance */
  Domain(class LAMMPS *);

  /** Destructor - cleans up regions, lattice, and region map */
  ~Domain() override;

  /** Initialize domain before a run */
  virtual void init();

  /** Set initial box dimensions
   * \param expandflag 1 to expand box for shrink-wrap, 0 otherwise (default 1) */
  void set_initial_box(int expandflag = 1);

  /** Set global box parameters */
  virtual void set_global_box();

  /** Set lamda box for triclinic */
  virtual void set_lamda_box();

  /** Set local (processor) box boundaries */
  virtual void set_local_box();

  /** Reset box after atoms have moved */
  virtual void reset_box();

  /** Apply periodic boundary conditions to atoms */
  virtual void pbc();

  /** Check if image flags are consistent */
  void image_check();

  /** Check if box is too small for cutoffs */
  void box_too_small_check();

  /** Check if subboxes are too small for communication
   * \param thresh Threshold distance */
  void subbox_too_small_check(double);

  /** Apply minimum image convention to displacement
   * \param file Source file for error message
   * \param line Line number for error message
   * \param dx X component of displacement
   * \param dy Y component of displacement
   * \param dz Z component of displacement */
  void minimum_image(const std::string &, int, double &, double &, double &) const;

  /** Apply minimum image convention (array version)
   * \param file Source file for error message
   * \param line Line number for error message
   * \param delta Displacement vector */
  void minimum_image(const std::string &file, int line, double *delta) const
  {
    minimum_image(file, line, delta[0], delta[1], delta[2]);
  }

  /** Apply minimum image for large displacements
   * \param file Source file for error message
   * \param line Line number for error message
   * \param dx X component of displacement
   * \param dy Y component of displacement
   * \param dz Z component of displacement */
  void minimum_image_big(const std::string &, int, double &, double &, double &) const;

  /** Apply minimum image for large displacements (array version)
   * \param file Source file for error message
   * \param line Line number for error message
   * \param delta Displacement vector */
  void minimum_image_big(const std::string &file, int line, double *delta) const
  {
    minimum_image_big(file, line, delta[0], delta[1], delta[2]);
  }

  /** Find closest periodic image of atom j to atom i
   * \param i Index of reference atom
   * \param j Index of atom to find image of
   * \return Index of closest image */
  int closest_image(int, int);

  /** Find closest periodic image of atom j to a point
   * \param pos Reference position
   * \param j Index of atom to find image of
   * \return Index of closest image */
  int closest_image(const double *const, int);

  /** Calculate position of closest image
   * \param xi Reference position
   * \param xj Atom position
   * \param xjimage Output: position of closest image */
  void closest_image(const double *const, const double *const, double *const);

  /** Remap atom into periodic box, update image flags
   * \param x Atom position (modified in place)
   * \param image Image flags (modified in place) */
  void remap(double *, imageint &);

  /** Remap atom into periodic box (no image update)
   * \param x Atom position (modified in place) */
  void remap(double *);

  /** Remap all atoms into periodic box */
  void remap_all();

  /** Remap atom near a reference position
   * \param x Atom position (modified in place)
   * \param xref Reference position */
  void remap_near(double *, double *);

  /** Unwrap atom position using image flags (inverse operation)
   * \param x Atom position (modified in place)
   * \param image Image flags */
  void unmap_inv(double *x, imageint);

  /** Unwrap atom position using image flags
   * \param x Atom position (modified in place)
   * \param image Image flags */
  void unmap(double *, imageint);

  /** Unwrap atom position into separate array
   * \param x Atom position
   * \param image Image flags
   * \param xnew Output: unwrapped position */
  void unmap(const double *, imageint, double *);

  /** Flip image flags when box is flipped
   * \param flipx Flip in x direction
   * \param flipy Flip in y direction
   * \param flipz Flip in z direction */
  void image_flip(int, int, int);

  /** Check if this processor owns an atom
   * \param id Atom ID
   * \param x Atom position
   * \param image Image flags
   * \param flag Ownership check mode
   * \return 1 if owned, 0 if not */
  int ownatom(int, double *, imageint *, int);

  /** Define a general triclinic box
   * \param avec A edge vector
   * \param bvec B edge vector
   * \param cvec C edge vector
   * \param origin Box origin */
  void define_general_triclinic(double *, double *, double *, double *);

  /** Calculate rotation matrix from general to restricted triclinic
   * \param avec A edge vector
   * \param bvec B edge vector
   * \param cvec C edge vector
   * \param rotate Output: rotation matrix
   * \param avec_new Output: rotated A vector
   * \param bvec_new Output: rotated B vector
   * \param cvec_new Output: rotated C vector */
  void general_to_restricted_rotation(double *, double *, double *, double[3][3], double *,
                                      double *, double *);

  /** Convert coordinates from general to restricted triclinic
   * \param x Coordinates (modified in place) */
  void general_to_restricted_coords(double *);

  /** Convert coordinates from restricted to general triclinic
   * \param x Coordinates (modified in place) */
  void restricted_to_general_coords(double *);

  /** Convert coordinates from restricted to general triclinic (with output)
   * \param x Input coordinates
   * \param xnew Output coordinates */
  void restricted_to_general_coords(double *, double *);

  /** Convert vector from general to restricted triclinic
   * \param v Vector (modified in place) */
  void general_to_restricted_vector(double *);

  /** Convert vector from restricted to general triclinic
   * \param v Vector (modified in place) */
  void restricted_to_general_vector(double *);

  /** Convert vector from restricted to general triclinic (with output)
   * \param v Input vector
   * \param vnew Output vector */
  void restricted_to_general_vector(double *, double *x);

  /** Set the lattice
   * \param narg Number of arguments
   * \param arg Argument array */
  void set_lattice(int, char **);

  /** Add a region
   * \param narg Number of arguments
   * \param arg Argument array */
  void add_region(int, char **);

  /** Delete a region by pointer
   * \param region Pointer to Region instance */
  void delete_region(Region *);

  /** Delete a region by ID
   * \param id Region ID */
  void delete_region(const std::string &);

  /** Get a region by its ID
   * \param id Region ID
   * \return Pointer to Region, or nullptr if not found */
  Region *get_region_by_id(const std::string &) const;

  /** Get all regions matching a style
   * \param style Style name
   * \return Vector of matching Region pointers */
  const std::vector<Region *> get_region_by_style(const std::string &) const;

  /** Get list of all regions
   * \return Vector of all Region pointers */
  const std::vector<Region *> get_region_list();

  /** Set boundary conditions
   * \param narg Number of arguments
   * \param arg Argument array
   * \param flag Flag for error checking */
  void set_boundary(int, char **, int);

  /** Print box info to screen/log
   * \param prefix String prefix for output */
  void print_box(const std::string &);

  /** Convert boundary settings to string
   * \param str Output string buffer */
  void boundary_string(char *);

  /** Convert lamda to Cartesian coords for all atoms
   * \param flag Processing flag */
  virtual void lamda2x(int);

  /** Convert lamda to Cartesian for range of atoms
   * \param lo First atom index
   * \param hi Last atom index */
  virtual void lamda2x(int, int);

  /** Convert Cartesian to lamda coords for all atoms
   * \param flag Processing flag */
  virtual void x2lamda(int);

  /** Convert Cartesian to lamda for range of atoms
   * \param lo First atom index
   * \param hi Last atom index */
  virtual void x2lamda(int, int);

  /** Convert lamda to Cartesian (single point)
   * \param lamda Input lamda coords
   * \param x Output Cartesian coords */
  virtual void lamda2x(double *, double *);

  /** Convert Cartesian to lamda (single point)
   * \param x Input Cartesian coords
   * \param lamda Output lamda coords */
  virtual void x2lamda(double *, double *);

  /** Check if point is inside global box
   * \param x Point coordinates
   * \return 1 if inside, 0 if outside */
  int inside(double *);

  /** Check if point is inside non-periodic dimensions
   * \param x Point coordinates
   * \return 1 if inside, 0 if outside */
  int inside_nonperiodic(double *);

  /** Convert Cartesian to lamda with bounds
   * \param x Input coords
   * \param lamda Output lamda coords
   * \param boxlo Box lower bounds
   * \param h Shape matrix */
  void x2lamda(double *, double *, double *, double *);

  /** Calculate bounding box of a region
   * \param xlo Output: lower bounds
   * \param xhi Output: upper bounds
   * \param lo Input: minimum bounds
   * \param hi Input: maximum bounds */
  void bbox(double *, double *, double *, double *);

  /** Calculate corners of the simulation box */
  void box_corners();

  /** Calculate corners of this processor's subbox */
  void subbox_corners();

  /** Calculate corners in lamda coordinates
   * \param lo Lamda lower bounds
   * \param hi Lamda upper bounds */
  void lamda_box_corners(double *, double *);

  // minimum image convention check
  // return 1 if any distance > 1/2 of box size
  // indicates a special neighbor is actually not in a bond,
  //   but is a far-away image that should be treated as an unbonded neighbor
  // inline since called from neighbor build inner loop

  /** Check if minimum image applies to displacement
   * \param dx X displacement
   * \param dy Y displacement
   * \param dz Z displacement
   * \return 1 if any component exceeds half box, 0 otherwise */
  inline int minimum_image_check(double dx, double dy, double dz) const
  {
    if (xperiodic && fabs(dx) > xprd_half) return 1;
    if (yperiodic && fabs(dy) > yprd_half) return 1;
    if (zperiodic && fabs(dz) > zprd_half) return 1;
    return 0;
  }

 protected:
  double small[3];    /**< Fractions of box lengths for tolerance */
  std::unordered_set<Region *> regions;    /**< Set of active regions */
};

}    // namespace LAMMPS_NS

#endif
