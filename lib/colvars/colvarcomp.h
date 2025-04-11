// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARCOMP_H
#define COLVARCOMP_H

// Declaration of colvar::cvc base class and derived ones.
//
// Future cvc's could be declared on additional header files.
// After the declaration of a new derived class, its metric
// functions must be reimplemented as well.
// If the new cvc has no symmetry or periodicity,
// this can be done straightforwardly by using the macro:
// simple_scalar_dist_functions (derived_class)

#include <memory>

#include "colvarmodule.h"
#include "colvaratoms.h"
#include "colvar.h"
#include "colvar_geometricpath.h"


/// \brief Colvar component (base class for collective variables)
///
/// A \link colvar::cvc \endlink object (or an object of a
/// cvc-derived class) implements the calculation of a collective
/// variable, its gradients and any other related physical quantities
/// that depend on microscopic degrees of freedom.
///
/// No restriction is set to what kind of calculation a \link colvar::cvc \endlink
/// object performs (usually an analytical function of atomic coordinates).
/// The only constraints are that: \par
///
/// - The value is calculated by the \link calc_value() \endlink
///   method, and is an object of \link colvarvalue \endlink class.  This
///   provides a transparent way to treat scalar and non-scalar variables
///   alike, and allows an automatic selection of the applicable algorithms.
///
/// - The object provides an implementation \link apply_force() \endlink to
///   apply forces to atoms.  Typically, one or more \link colvarmodule::atom_group
///   \endlink objects are used, but this is not a requirement for as long as
///   the \link colvar::cvc \endlink object communicates with the simulation program.
///
/// <b> If you wish to implement a new collective variable component, you
/// should write your own class by inheriting directly from \link
/// colvar::cvc \endlink, or one of its derived classes (for instance,
/// \link colvar::distance \endlink is frequently used, because it provides
/// useful data and function members for any colvar based on two
/// atom groups).</b>
///
/// The steps are: \par
/// 1. Declare the new class as a derivative of \link colvar::cvc \endlink
///    in the file \link colvarcomp.h \endlink
/// 2. Implement the new class in a file named colvarcomp_<something>.cpp
/// 3. Declare the name of the new class inside the \link colvar \endlink class
///    in \link colvar.h \endlink (see "list of available components")
/// 4. Add a call for the new class in colvar::init_components()
////   (file: colvar.cpp)
///

class colvar::cvc
  : public colvarparse, public colvardeps
{
public:

  /// \brief The name of the object (helps to identify this
  /// cvc instance when debugging)
  std::string name;

  /// String identifier for the type of collective variable
  std::string function_type() const;

  /// Keyword used in the input to denote this CVC
  std::string config_key;

  /// \brief Coefficient in the polynomial combination (default: 1.0)
  cvm::real sup_coeff = 1.0;

  /// \brief Exponent in the polynomial combination (default: 1)
  int       sup_np = 1;

  /// \brief Period of the values of this CVC (default: 0.0, non periodic)
  cvm::real period = 0.0;

  /// \brief If the component is periodic, wrap around this value (default: 0.0)
  cvm::real wrap_center = 0.0;

  /// Constructor
  cvc();

  /// Destructor
  virtual ~cvc();

  /// An init function should be defined for every class inheriting from cvc
  /// \param conf Contents of the configuration file pertaining to this \link
  /// cvc \endlink
  virtual int init(std::string const &conf);

  /// \brief Initialize dependency tree
  virtual int init_dependencies();

  /// \brief After construction, set data related to dependency handling
  int setup();

  /// \brief Implementation of the feature list accessor for colvar
  virtual const std::vector<feature *> &features() const
  {
    return cvc_features;
  }

  virtual std::vector<feature *> &modify_features()
  {
    return cvc_features;
  }

  static void delete_features() {
    for (size_t i=0; i < cvc_features.size(); i++) {
      delete cvc_features[i];
    }
    cvc_features.clear();
  }

  /// \brief Get vector of vectors of atom IDs for all atom groups
  virtual std::vector<std::vector<int> > get_atom_lists();

  /// \brief Obtain data needed for the calculation for the backend
  virtual void read_data();

  /// \brief Calculate the variable
  virtual void calc_value() = 0;

  /// \brief Calculate the atomic gradients, to be reused later in
  /// order to apply forces
  virtual void calc_gradients() {}

  /// \brief Calculate the atomic fit gradients
  void calc_fit_gradients();

  /// \brief Calculate finite-difference gradients alongside the analytical ones, for each Cartesian component
  virtual void debug_gradients();

  /// \brief Calculate atomic gradients and add them to the corresponding item in gradient vector
  /// May be overridden by CVCs that do not store their gradients in the classic way, see dihedPC
  virtual void collect_gradients(std::vector<int> const &atom_ids, std::vector<cvm::rvector> &atomic_gradients);

  /// \brief Calculate the total force from the system using the
  /// inverse atomic gradients
  virtual void calc_force_invgrads();

  /// \brief Calculate the divergence of the inverse atomic gradients
  virtual void calc_Jacobian_derivative();


  /// \brief Return the previously calculated value
  colvarvalue const & value() const;

  /// \brief Return the previously calculated total force
  colvarvalue const & total_force() const;

  /// \brief Return the previously calculated divergence of the
  /// inverse atomic gradients
  colvarvalue const & Jacobian_derivative() const;

  /// \brief Apply the collective variable force, by communicating the
  /// atomic forces to the simulation program (\b Note: the \link ft
  /// \endlink member is not altered by this function)
  ///
  /// Note: multiple calls to this function within the same simulation
  /// step will add the forces altogether \param cvforce The
  /// collective variable force, usually coming from the biases and
  /// eventually manipulated by the parent \link colvar \endlink
  /// object
  virtual void apply_force(colvarvalue const &cvforce);

  /// Square distance between x1 and x2 (can be redefined to transparently
  /// implement metrics in multi-dimensional spaces with or without
  /// constraints, symmetries and periodicities).  The default implementation
  /// assumes scalar numbers and no symmetries or periodicities.
  ///
  /// If symmetries or periodicities are present, the
  /// colvar::cvc::dist2() should be redefined to return the
  /// "closest distance" value and colvar::cvc::dist2_lgrad(),
  /// colvar::cvc::dist2_rgrad() to return its gradients.
  ///
  /// If constraints are present (and not already implemented by any
  /// of the \link colvarvalue \endlink types), the
  /// colvar::cvc::dist2_lgrad() and
  /// colvar::cvc::dist2_rgrad() functions should be redefined
  /// to provide a gradient which is compatible with the constraint,
  /// i.e. already deprived of its component normal to the constraint
  /// hypersurface.
  ///
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;

  /// \brief Gradient(with respect to x1) of the square distance (can
  /// be redefined to transparently implement constraints, symmetries
  /// and periodicities)
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;

  /// \brief Gradient(with respect to x2) of the square distance (can
  /// be redefined to transparently implement constraints, symmetries
  /// and periodicities)
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;

  /// \brief Wrap value (for periodic/symmetric cvcs)
  virtual void wrap(colvarvalue &x_unwrapped) const;

  /// \brief Pointers to all atom groups, to let colvars collect info
  /// e.g. atomic gradients
  std::vector<cvm::atom_group *> atom_groups;

  /// \brief Store a pointer to new atom group, and list as child for dependencies
  void register_atom_group(cvm::atom_group *ag);

  /// Pointer to the gradient of parameter param_name
  virtual colvarvalue const *get_param_grad(std::string const &param_name);

  /// Set the named parameter to the given value
  virtual int set_param(std::string const &param_name, void const *new_value);

  /// \brief Whether or not this CVC will be computed in parallel whenever possible
  bool b_try_scalable = true;

  /// Forcibly set value of CVC - useful for driving an external coordinate,
  /// eg. lambda dynamics
  inline void set_value(colvarvalue const &new_value) {
    x = new_value;
  }

protected:

  /// Set the value of \link function_type \endlink and its dependencies
  int set_function_type(std::string const &type);

  /// Update the description string based on name and type
  int update_description();

  /// Parse a group definition
  cvm::atom_group *parse_group(std::string const &conf, char const *group_key,
                               bool optional = false);

  /// \brief Parse options pertaining to total force calculation
  virtual int init_total_force_params(std::string const &conf);

  /// \brief Implementation of the feature list for colvar
  static std::vector<feature *> cvc_features;

  /// Record the type of this class as well as those it is derived from
  std::vector<std::string> function_types;

  /// \brief Cached value
  colvarvalue x;

  /// \brief Value at the previous step
  colvarvalue x_old;

  /// \brief Calculated total force (\b Note: this is calculated from
  /// the total atomic forces read from the program, subtracting fromt
  /// the "internal" forces of the system the "external" forces from
  /// the colvar biases)
  colvarvalue ft;

  /// \brief Calculated Jacobian derivative (divergence of the inverse
  /// gradients): serves to calculate the phase space correction
  colvarvalue jd;

  /// \brief Set data types for a scalar distance (convenience function)
  void init_as_distance();

  /// \brief Set data types for a bounded angle (0° to 180°)
  void init_as_angle();

  /// \brief Set data types for a periodic angle (-180° to 180°)
  void init_as_periodic_angle();

  /// \brief Set two scalar boundaries (convenience function)
  void init_scalar_boundaries(cvm::real lb, cvm::real ub);

  /// \brief Location of the lower boundary (not defined by user choice)
  colvarvalue lower_boundary;

  /// \brief Location of the upper boundary (not defined by user choice)
  colvarvalue upper_boundary;

  /// \brief CVC-specific default colvar width (default: not provided)
  cvm::real width = 0.0;
};


inline colvarvalue const & colvar::cvc::value() const
{
  return x;
}


inline colvarvalue const & colvar::cvc::total_force() const
{
  return ft;
}


inline colvarvalue const & colvar::cvc::Jacobian_derivative() const
{
  return jd;
}



/// \brief Colvar component: distance between the centers of mass of
/// two groups (colvarvalue::type_scalar type, range [0:*))

class colvar::distance
  : public colvar::cvc
{
protected:
  /// First atom group
  cvm::atom_group  *group1 = nullptr;
  /// Second atom group
  cvm::atom_group  *group2 = nullptr;
  /// Vector distance, cached to be recycled
  cvm::rvector     dist_v;
public:
  distance();
  virtual ~distance() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
};



// \brief Colvar component: distance vector between centers of mass
// of two groups (\link colvarvalue::type_3vector \endlink type,
// range (-*:*)x(-*:*)x(-*:*))
class colvar::distance_vec
  : public colvar::distance
{
public:
  distance_vec();
  virtual ~distance_vec() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  /// Redefined to deal with multiple dimensions
  virtual cvm::real dist2(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual void wrap(colvarvalue &x_unwrapped) const;
};


/// \brief Colvar component: distance unit vector (direction) between
/// centers of mass of two groups (colvarvalue::type_unit3vector type,
/// range [-1:1]x[-1:1]x[-1:1])
class colvar::distance_dir
  : public colvar::distance
{
public:
  distance_dir();
  virtual ~distance_dir() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  /// Redefined to deal with multiple dimensions
  virtual cvm::real dist2(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual void wrap(colvarvalue &x_unwrapped) const;
};



/// \brief Colvar component: projection of the distance vector along
/// an axis(colvarvalue::type_scalar type, range (-*:*))
class colvar::distance_z
  : public colvar::cvc
{
protected:
  /// Main atom group
  cvm::atom_group  *main = nullptr;
  /// Reference atom group
  cvm::atom_group  *ref1 = nullptr;
  /// Optional, second ref atom group
  cvm::atom_group  *ref2 = nullptr;
  /// Vector on which the distance vector is projected
  cvm::rvector axis;
  /// Norm of the axis
  cvm::real axis_norm;
  /// Vector distance, cached to be recycled
  cvm::rvector     dist_v;
  /// Flag: using a fixed axis vector?
  bool fixed_axis = true;
public:
  distance_z();
  virtual ~distance_z() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
};



/// \brief Colvar component: projection of the distance vector on a
/// plane (colvarvalue::type_scalar type, range [0:*))
class colvar::distance_xy
  : public colvar::distance_z
{
protected:
  /// Components of the distance vector orthogonal to the axis
  cvm::rvector dist_v_ortho;
  /// Vector distances
  cvm::rvector v12, v13;
public:
  distance_xy();
  virtual ~distance_xy() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
};


/// \brief Colvar component: polar coordinate phi of a group
/// (colvarvalue::type_scalar type, range [-180:180])
class colvar::polar_phi
  : public colvar::cvc
{
protected:
  cvm::atom_group *atoms = nullptr;
  cvm::real r, theta, phi;

public:
  polar_phi();
  virtual ~polar_phi() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
};


/// \brief Colvar component: polar coordinate theta of a group
/// (colvarvalue::type_scalar type, range [0:180])
class colvar::polar_theta
  : public colvar::cvc
{
public:
  polar_theta();
  virtual ~polar_theta() {}
  virtual int init(std::string const &conf);
protected:
  cvm::atom_group  *atoms = nullptr;
  cvm::real r, theta, phi;
public:
  virtual void calc_value();
  virtual void calc_gradients();
};


/// \brief Colvar component: average distance between two groups of atoms, weighted as the sixth power,
/// as in NMR refinements(colvarvalue::type_scalar type, range (0:*))
class colvar::distance_inv
  : public colvar::cvc
{
protected:
  /// First atom group
  cvm::atom_group  *group1 = nullptr;
  /// Second atom group
  cvm::atom_group  *group2 = nullptr;
  /// Components of the distance vector orthogonal to the axis
  int exponent = 6;
public:
  distance_inv();
  virtual ~distance_inv() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
};



/// \brief Colvar component: N1xN2 vector of pairwise distances
/// (colvarvalue::type_vector type, range (0:*) for each component)
class colvar::distance_pairs
  : public colvar::cvc
{
protected:
  /// First atom group
  cvm::atom_group  *group1 = nullptr;
  /// Second atom group
  cvm::atom_group  *group2 = nullptr;
public:
  distance_pairs();
  virtual ~distance_pairs() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  /// Redefined to deal with multiple dimensions
  virtual cvm::real dist2(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual void wrap(colvarvalue &x_unwrapped) const;
};


/// \brief Colvar component:  dipole magnitude of a molecule
class colvar::dipole_magnitude
  : public colvar::cvc
{
protected:
  /// Dipole atom group
  cvm::atom_group  *atoms = nullptr;
  cvm::atom_pos dipoleV;
public:
  dipole_magnitude();
  virtual ~dipole_magnitude() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
};



/// \brief Colvar component: Radius of gyration of an atom group
/// (colvarvalue::type_scalar type, range [0:*))
class colvar::gyration
  : public colvar::cvc
{
protected:
  /// Atoms involved
  cvm::atom_group  *atoms = nullptr;
public:
  gyration();
  virtual ~gyration() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
};



/// \brief Colvar component: moment of inertia of an atom group
/// (colvarvalue::type_scalar type, range [0:*))
class colvar::inertia
  : public colvar::gyration
{
public:
  inertia();
  virtual ~inertia() {}
  virtual void calc_value();
  virtual void calc_gradients();
};



/// \brief Colvar component: moment of inertia of an atom group
/// around a user-defined axis (colvarvalue::type_scalar type, range [0:*))
class colvar::inertia_z
  : public colvar::inertia
{
protected:
  /// Vector on which the inertia tensor is projected
  cvm::rvector axis;
public:
  inertia_z();
  virtual ~inertia_z() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
};



/// \brief Colvar component: projection of 3N coordinates onto an
/// eigenvector(colvarvalue::type_scalar type, range (-*:*))
class colvar::eigenvector
  : public colvar::cvc
{
protected:

  /// Atom group
  cvm::atom_group  *           atoms = nullptr;

  /// Reference coordinates
  std::vector<cvm::atom_pos>  ref_pos;

  /// Eigenvector (of a normal or essential mode): will always have zero center
  std::vector<cvm::rvector>   eigenvec;

  /// Inverse square norm of the eigenvector
  cvm::real                   eigenvec_invnorm2;

public:

  eigenvector();
  virtual ~eigenvector() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
};



/// \brief Colvar component: angle between the centers of mass of
/// three groups (colvarvalue::type_scalar type, range [0:PI])
class colvar::angle
  : public colvar::cvc
{
protected:

  /// Atom group
  cvm::atom_group  *group1 = nullptr;
  /// Atom group
  cvm::atom_group  *group2 = nullptr;
  /// Atom group
  cvm::atom_group  *group3 = nullptr;

  /// Inter site vectors
  cvm::rvector r21, r23;
  /// Inter site vector norms
  cvm::real r21l, r23l;
  /// Derivatives wrt group centers of mass
  cvm::rvector dxdr1, dxdr3;

  /// Compute total force on first site only to avoid unwanted
  /// coupling to other colvars (see e.g. Ciccotti et al., 2005)
  /// (or to allow dummy atoms)
  bool b_1site_force = false;
public:

  angle();
  /// \brief Initialize the three groups after three atoms
  angle(cvm::atom const &a1, cvm::atom const &a2, cvm::atom const &a3);
  virtual ~angle() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
};



/// \brief Colvar component: angle between the dipole of a molecule and an axis
/// formed by two groups of atoms(colvarvalue::type_scalar type, range [0:PI])
class colvar::dipole_angle
  : public colvar::cvc
{
protected:

  /// Dipole atom group
  cvm::atom_group  *group1 = nullptr;
  /// Atom group
  cvm::atom_group  *group2 = nullptr;
  /// Atom group
  cvm::atom_group  *group3 = nullptr;

  /// Inter site vectors
  cvm::rvector r21, r23;
  /// Inter site vector norms
  cvm::real r21l, r23l;
  /// Derivatives wrt group centers of mass
  cvm::rvector dxdr1, dxdr3;

  /// Compute total force on first site only to avoid unwanted
  /// coupling to other colvars (see e.g. Ciccotti et al., 2005)
  /// (or to allow dummy atoms)
  bool b_1site_force = false;
public:

  dipole_angle();
  virtual ~dipole_angle() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
};



/// \brief Colvar component: dihedral between the centers of mass of
/// four groups (colvarvalue::type_scalar type, range [-PI:PI])
class colvar::dihedral
  : public colvar::cvc
{
protected:

  /// Atom group
  cvm::atom_group  *group1 = nullptr;
  /// Atom group
  cvm::atom_group  *group2 = nullptr;
  /// Atom group
  cvm::atom_group  *group3 = nullptr;
  /// Atom group
  cvm::atom_group  *group4 = nullptr;
  /// Inter site vectors
  cvm::rvector r12, r23, r34;

  /// \brief Compute total force on first site only to avoid unwanted
  /// coupling to other colvars (see e.g. Ciccotti et al., 2005)
  bool b_1site_force = false;

public:

  /// \brief Initialize the four groups after four atoms
  dihedral(cvm::atom const &a1, cvm::atom const &a2, cvm::atom const &a3, cvm::atom const &a4);
  dihedral();
  virtual ~dihedral() {}
  virtual int init(std::string  const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
};



/// \brief Colvar component: coordination number between two groups
/// (colvarvalue::type_scalar type, range [0:N1*N2])
class colvar::coordnum
  : public colvar::cvc
{
protected:
  /// First atom group
  cvm::atom_group  *group1 = nullptr;
  /// Second atom group
  cvm::atom_group  *group2 = nullptr;
  /// \brief "Cutoff" for isotropic calculation (default)
  cvm::real     r0;
  /// \brief "Cutoff vector" for anisotropic calculation
  cvm::rvector  r0_vec;
  /// \brief Whether r/r0 or \vec{r}*\vec{1/r0_vec} should be used
  bool b_anisotropic = false;
  /// Integer exponent of the function numerator
  int en = 6;
  /// Integer exponent of the function denominator
  int ed = 12;

  /// If true, group2 will be treated as a single atom
  bool b_group2_center_only = false;

  /// Tolerance for the pair list
  cvm::real tolerance = 0.0;

  /// Frequency of update of the pair list
  int pairlist_freq = 100;

  /// Pair list
  bool *pairlist = nullptr;

public:

  coordnum();
  virtual ~coordnum();
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();

  enum {
    ef_null = 0,
    ef_gradients = 1,
    ef_anisotropic = (1<<8),
    ef_use_pairlist = (1<<9),
    ef_rebuild_pairlist = (1<<10)
  };

  /// \brief Calculate a coordination number through the function
  /// (1-x**n)/(1-x**m), where x = |A1-A2|/r0 \param r0, r0_vec "cutoff" for
  /// the coordination number (scalar or vector depending on user choice)
  /// \param en Numerator exponent \param ed Denominator exponent \param First
  /// atom \param Second atom \param pairlist_elem pointer to pair flag for
  /// this pair \param tolerance A pair is defined as having a larger
  /// coordination than this number
  template<int flags>
  static cvm::real switching_function(cvm::real const &r0,
                                      cvm::rvector const &r0_vec,
                                      int en,
                                      int ed,
                                      cvm::atom &A1,
                                      cvm::atom &A2,
                                      bool **pairlist_elem,
                                      cvm::real tolerance);

  /// Workhorse function
  template<int flags> int compute_coordnum();

  /// Workhorse function
  template<int flags> void main_loop(bool **pairlist_elem);

};



/// \brief Colvar component: self-coordination number within a group
/// (colvarvalue::type_scalar type, range [0:N*(N-1)/2])
class colvar::selfcoordnum
  : public colvar::cvc
{
protected:

  /// Selected atoms
  cvm::atom_group  *group1 = nullptr;
  /// \brief "Cutoff" for isotropic calculation (default)
  cvm::real     r0;
  /// Integer exponent of the function numerator
  int en = 6;
  /// Integer exponent of the function denominator
  int ed = 12;
  cvm::real tolerance = 0.0;
  int pairlist_freq = 100;

  bool *pairlist = nullptr;

public:

  selfcoordnum();
  ~selfcoordnum();
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();

  /// Main workhorse function
  template<int flags> int compute_selfcoordnum();
};



/// \brief Colvar component: coordination number between two groups
/// (colvarvalue::type_scalar type, range [0:N1*N2])
class colvar::groupcoordnum
  : public colvar::distance
{
protected:
  /// \brief "Cutoff" for isotropic calculation (default)
  cvm::real     r0;
  /// \brief "Cutoff vector" for anisotropic calculation
  cvm::rvector  r0_vec;
  /// \brief Wheter dist/r0 or \vec{dist}*\vec{1/r0_vec} should ne be
  /// used
  bool b_anisotropic = false;
  /// Integer exponent of the function numerator
  int en = 6;
  /// Integer exponent of the function denominator
  int ed = 12;
public:
  groupcoordnum();
  virtual ~groupcoordnum() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
};



/// \brief Colvar component: hydrogen bond, defined as the product of
/// a colvar::coordnum and 1/2*(1-cos((180-ang)/ang_tol))
/// (colvarvalue::type_scalar type, range [0:1])
class colvar::h_bond
  : public colvar::cvc
{
protected:
  /// \brief "Cutoff" distance between acceptor and donor
  cvm::real     r0;
  /// Integer exponent of the function numerator
  int en = 6;
  /// Integer exponent of the function denominator
  int ed = 8;
public:
  /// Constructor for atoms already allocated
  h_bond(cvm::atom const &acceptor,
         cvm::atom const &donor,
         cvm::real r0, int en, int ed);
  h_bond();
  virtual ~h_bond() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
};


/// \brief Colvar component: alpha helix content of a contiguous
/// segment of 5 or more residues, implemented as a sum of Ca-Ca-Ca
/// angles and hydrogen bonds (colvarvalue::type_scalar type, range
/// [0:1])
class colvar::alpha_angles
  : public colvar::cvc
{
protected:

  /// Reference Calpha-Calpha angle (default: 88 degrees)
  cvm::real theta_ref = 88.0;

  /// Tolerance on the Calpha-Calpha angle
  cvm::real theta_tol = 15.0;

  /// List of Calpha-Calpha angles
  std::vector<angle *> theta;

  /// List of hydrogen bonds
  std::vector<h_bond *>   hb;

  /// Contribution of the HB terms
  cvm::real hb_coeff = 0.5;

  /// Cutoff for HB
  cvm::real r0;

  /// Integer exponent of the HB numerator
  int en = 6;

  /// Integer exponent of the HB denominator
  int ed = 8;

public:

  alpha_angles();
  virtual ~alpha_angles();
  virtual int init(std::string const &conf);
  void calc_value();
  void calc_gradients();
  /// Re-implementation of cvc::collect_gradients() to carry over atomic gradients of sub-cvcs
  void collect_gradients(std::vector<int> const &atom_ids, std::vector<cvm::rvector> &atomic_gradients);
  void apply_force(colvarvalue const &force);
};



/// \brief Colvar component: dihedPC
/// Projection of the config onto a dihedral principal component
/// See e.g. Altis et al., J. Chem. Phys 126, 244111 (2007)
/// Based on a set of 'dihedral' cvcs
class colvar::dihedPC
  : public colvar::cvc
{
protected:

  std::vector<dihedral *> theta;
  std::vector<cvm::real> coeffs;

public:

  dihedPC();
  virtual ~dihedPC();
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  /// Re-implementation of cvc::collect_gradients() to carry over atomic gradients of sub-cvcs
  virtual void collect_gradients(std::vector<int> const &atom_ids, std::vector<cvm::rvector> &atomic_gradients);
  virtual void apply_force(colvarvalue const &force);
};



/// \brief Colvar component: orientation in space of an atom group,
/// with respect to a set of reference coordinates
/// (colvarvalue::type_quaternion type, range
/// [-1:1]x[-1:1]x[-1:1]x[-1:1])
class colvar::orientation
  : public colvar::cvc
{
protected:

  /// Atom group
  cvm::atom_group  *          atoms = nullptr;
  /// Center of geometry of the group
  cvm::atom_pos              atoms_cog;

  /// Reference coordinates
  std::vector<cvm::atom_pos> ref_pos;

  /// Shifted atomic positions
  std::vector<cvm::atom_pos> shifted_pos;

  /// Rotation object
  cvm::rotation              rot;

  /// \brief This is used to remove jumps in the sign of the
  /// quaternion, which may be annoying in the colvars trajectory
  cvm::quaternion            ref_quat;

  /// Rotation derivative
  struct rotation_derivative_impl_;
  std::unique_ptr<rotation_derivative_impl_> rot_deriv_impl;

public:

  orientation();
  virtual ~orientation();
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  /// Redefined to use quaternion metrics
  virtual cvm::real dist2(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to use quaternion metrics
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to use quaternion metrics
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to use quaternion metrics
  virtual void wrap(colvarvalue &x_unwrapped) const;
};


/// Colvar component: angle of rotation with respect to a set of reference coordinates
/// (colvarvalue::type_scalar type, range [0:PI))
/// This is also used to derive all other sub-rotation variables that return a scalar value
class colvar::orientation_angle
  : public colvar::orientation
{
public:

  orientation_angle();
  virtual ~orientation_angle() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  /// Redefined to use scalar metrics
  virtual cvm::real dist2(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to use scalar metrics
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to use scalar metrics
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to use scalar metrics
  virtual void wrap(colvarvalue &x_unwrapped) const;
};



/// \brief Colvar component: cosine of the angle of rotation with respect to a set
/// of reference coordinates (colvarvalue::type_scalar type, range
/// [-1:1])
class colvar::orientation_proj
  : public colvar::orientation_angle
{
public:

  orientation_proj();
  virtual ~orientation_proj() {}
  virtual void calc_value();
  virtual void calc_gradients();
};



/// \brief Colvar component: projection of the orientation vector onto
/// a predefined axis (colvarvalue::type_scalar type, range [-1:1])
class colvar::tilt
  : public colvar::orientation_proj
{
protected:

  cvm::rvector axis;

public:

  tilt();
  virtual ~tilt() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
};



/// \brief Colvar component: angle of rotation around a predefined
/// axis (colvarvalue::type_scalar type, range [-PI:PI])
class colvar::spin_angle
  : public colvar::tilt
{
public:

  spin_angle();
  virtual ~spin_angle() {}
  virtual void calc_value();
  virtual void calc_gradients();
};


class colvar::euler_phi
  : public colvar::orientation_angle
{
public:
  euler_phi();
  virtual ~euler_phi() {}
  virtual void calc_value();
  virtual void calc_gradients();
};


class colvar::euler_psi
  : public colvar::orientation_angle
{
public:
  euler_psi();
  virtual ~euler_psi() {}
  virtual void calc_value();
  virtual void calc_gradients();
};


class colvar::euler_theta
  : public colvar::orientation_angle
{
public:
  euler_theta();
  virtual ~euler_theta() {}
  virtual void calc_value();
  virtual void calc_gradients();
};


/// \brief Colvar component: root mean square deviation (RMSD) of a
/// group with respect to a set of reference coordinates; uses \link
/// colvar::orientation \endlink to calculate the rotation matrix
/// (colvarvalue::type_scalar type, range [0:*))
class colvar::rmsd
  : public colvar::cvc
{
protected:

  /// Atom group
  cvm::atom_group  *atoms = nullptr;

  /// Reference coordinates (for RMSD calculation only)
  /// Includes sets with symmetry permutations (n_permutations * n_atoms)
  std::vector<cvm::atom_pos>  ref_pos;

  /// Number of permutations of symmetry-related atoms
  size_t n_permutations = 1;

  /// Index of the permutation yielding the smallest RMSD (0 for identity)
  size_t best_perm_index = 0;

  /// Permutation RMSD input parsing
  int init_permutation(std::string const &conf);

public:
  rmsd();
  virtual ~rmsd() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
};



// \brief Colvar component: flat vector of Cartesian coordinates
// Mostly useful to compute scripted colvar values
class colvar::cartesian
  : public colvar::cvc
{
protected:
  /// Atom group
  cvm::atom_group  *atoms = nullptr;
  /// Which Cartesian coordinates to include
  std::vector<size_t> axes;
public:
  cartesian();
  virtual ~cartesian() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  /// Redefined to deal with multiple dimensions
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to deal with multiple dimensions
  virtual void wrap(colvarvalue &x_unwrapped) const;
};


// \brief Colvar component: alch_lambda
// To communicate value with back-end in lambda-dynamics
class colvar::alch_lambda
  : public colvar::cvc
{
protected:
  // No atom groups needed
public:
  alch_lambda();
  virtual ~alch_lambda() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
};


// \brief Colvar component: alch_Flambda
// To communicate force on lambda with back-end in lambda-dynamics
class colvar::alch_Flambda
  : public colvar::cvc
{
protected:
  // No atom groups needed
public:
  alch_Flambda();
  virtual ~alch_Flambda() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
};


class colvar::CartesianBasedPath
  : public colvar::cvc
{
protected:
    virtual void computeDistanceBetweenReferenceFrames(std::vector<cvm::real>& result);
    virtual void computeDistanceToReferenceFrames(std::vector<cvm::real>& result);
    /// Selected atoms
    cvm::atom_group *atoms = nullptr;
    /// Fitting options
    bool has_user_defined_fitting = false;
    /// Reference frames
    std::vector<std::vector<cvm::atom_pos>> reference_frames;
    std::vector<std::vector<cvm::atom_pos>> reference_fitting_frames;
    /// Atom groups for RMSD calculation together with reference frames
    std::vector<cvm::atom_group*> comp_atoms;
    /// Total number of reference frames
    size_t total_reference_frames = 0;
public:
    CartesianBasedPath();
    virtual ~CartesianBasedPath();
    virtual int init(std::string const &conf);
    virtual void calc_value() = 0;
    /// Redefined to raise error because this is an abstract type
    virtual void apply_force(colvarvalue const &force);
};

/// \brief Colvar component: alternative path collective variable using geometry, variable s
/// For more information see https://plumed.github.io/doc-v2.5/user-doc/html/_p_a_t_h.html
/// Diaz Leines, G.; Ensing, B. Path Finding on High-Dimensional Free Energy Landscapes. Phys. Rev. Lett. 2012, 109 (2), 020601. https://doi.org/10.1103/PhysRevLett.109.020601.
class colvar::gspath
  : public colvar::CartesianBasedPath, public GeometricPathCV::GeometricPathBase<cvm::atom_pos, cvm::real, GeometricPathCV::path_sz::S>
{
private:
    // Optimal rotation for compute v3
    cvm::rotation rot_v3;
protected:
    virtual void prepareVectors();
    virtual void updateDistanceToReferenceFrames();
public:
    gspath();
    virtual ~gspath() {}
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};



/// \brief Colvar component: alternative path collective variable using geometry, variable z
/// This should be merged with gspath in the same class by class inheritance or something else
class colvar::gzpath
  : public colvar::CartesianBasedPath, public GeometricPathCV::GeometricPathBase<cvm::atom_pos, cvm::real, GeometricPathCV::path_sz::Z>
{
private:
    // Optimal rotation for compute v3, v4
    cvm::rotation rot_v3;
    cvm::rotation rot_v4;
protected:
    virtual void prepareVectors();
    virtual void updateDistanceToReferenceFrames();
public:
    gzpath();
    virtual ~gzpath() {}
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};

/// Current only linear combination of sub-CVCs is available
class colvar::linearCombination
  : public colvar::cvc
{
protected:
    /// Sub-colvar components
    std::vector<colvar::cvc*> cv;
    /// If all sub-cvs use explicit gradients then we also use it
    bool use_explicit_gradients;
protected:
    cvm::real getPolynomialFactorOfCVGradient(size_t i_cv) const;
public:
    linearCombination();
    virtual ~linearCombination();
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
  /// Redefined to allow arbitrary dimensions
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  /// Redefined to allow arbitrary dimensions
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to allow arbitrary dimensions
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to allow arbitrary dimensions
  virtual void wrap(colvarvalue &x_unwrapped) const;
};


/// custom expression of colvars
class colvar::customColvar
  : public colvar::linearCombination
{
protected:
    bool use_custom_function = false;
#ifdef LEPTON
    /// Vector of evaluators for custom functions using Lepton
    std::vector<Lepton::CompiledExpression *> value_evaluators;
    /// Vector of evaluators for gradients of custom functions
    std::vector<Lepton::CompiledExpression *> gradient_evaluators;
    /// Vector of references to cvc values to be passed to Lepton evaluators
    std::vector<double *> value_eval_var_refs;
    std::vector<double *> grad_eval_var_refs;
    /// Unused value that is written to when a variable simplifies out of a Lepton expression
    double dev_null = 0.0;
#endif
public:
    customColvar();
    virtual ~customColvar();
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};


class colvar::CVBasedPath
  : public colvar::cvc
{
protected:
    /// Sub-colvar components
    std::vector<colvar::cvc*> cv;
    /// Reference colvar values from path
    std::vector<std::vector<colvarvalue>> ref_cv;
    /// If all sub-cvs use explicit gradients then we also use it
    bool use_explicit_gradients;
    /// Total number of reference frames
    size_t total_reference_frames = 0;
protected:
    virtual void computeDistanceToReferenceFrames(std::vector<cvm::real>& result);
    /// Helper function to determine the distance between reference frames
    virtual void computeDistanceBetweenReferenceFrames(std::vector<cvm::real>& result) const;
    cvm::real getPolynomialFactorOfCVGradient(size_t i_cv) const;
public:
    CVBasedPath();
    virtual ~CVBasedPath();
    virtual int init(std::string const &conf);
    virtual void calc_value() = 0;
    /// Redefined to raise error because this is an abstract type
    virtual void apply_force(colvarvalue const &force);
  /// Redefined to use the metric of the returned colvarvalue (defined at runtime)
  virtual cvm::real dist2(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to use the metric of the returned colvarvalue (defined at runtime)
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to use the metric of the returned colvarvalue (defined at runtime)
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1, colvarvalue const &x2) const;
  /// Redefined to use the metric of the returned colvarvalue (defined at runtime)
  virtual void wrap(colvarvalue &x_unwrapped) const;
};


/// \brief Colvar component: alternative path collective variable using geometry, variable s
/// Allow any combination of existing (scalar) CVs
/// For more information see https://plumed.github.io/doc-v2.5/user-doc/html/_p_a_t_h.html
/// Diaz Leines, G.; Ensing, B. Path Finding on High-Dimensional Free Energy Landscapes. Phys. Rev. Lett. 2012, 109 (2), 020601. https://doi.org/10.1103/PhysRevLett.109.020601.
class colvar::gspathCV
  : public colvar::CVBasedPath, public GeometricPathCV::GeometricPathBase<colvarvalue, cvm::real, GeometricPathCV::path_sz::S>
{
protected:
    virtual void updateDistanceToReferenceFrames();
    virtual void prepareVectors();
public:
    gspathCV();
    virtual ~gspathCV();
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};



class colvar::gzpathCV
  : public colvar::CVBasedPath, public GeometricPathCV::GeometricPathBase<colvarvalue, cvm::real, GeometricPathCV::path_sz::Z>
{
protected:
    virtual void updateDistanceToReferenceFrames();
    virtual void prepareVectors();
public:
    gzpathCV();
    virtual ~gzpathCV();
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};

struct ArithmeticPathImpl;

class colvar::aspath
  : public colvar::CartesianBasedPath
{
private:
    std::unique_ptr<ArithmeticPathImpl> impl_;
    friend struct ArithmeticPathImpl;
public:
    aspath();
    virtual ~aspath();
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};

class colvar::azpath
  : public colvar::CartesianBasedPath
{
private:
    std::unique_ptr<ArithmeticPathImpl> impl_;
    friend struct ArithmeticPathImpl;
public:
    azpath();
    virtual ~azpath();
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};

class colvar::aspathCV
  : public colvar::CVBasedPath
{
private:
    std::unique_ptr<ArithmeticPathImpl> impl_;
    friend struct ArithmeticPathImpl;
public:
    aspathCV();
    virtual ~aspathCV();
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};


class colvar::azpathCV
  : public colvar::CVBasedPath
{
private:
    std::unique_ptr<ArithmeticPathImpl> impl_;
    friend struct ArithmeticPathImpl;
public:
    azpathCV();
    virtual ~azpathCV();
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};

// forward declaration
namespace neuralnetworkCV {
    class neuralNetworkCompute;
}


class colvar::neuralNetwork
  : public linearCombination
{
protected:
    /// actual computation happens in neuralnetworkCV::neuralNetworkCompute
    std::unique_ptr<neuralnetworkCV::neuralNetworkCompute> nn;
    /// the index of nn output components
    size_t m_output_index = 0;
public:
    neuralNetwork();
    virtual ~neuralNetwork();
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
  /// Redefined to allow arbitrary dimensions
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  /// Redefined to allow arbitrary dimensions
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to allow arbitrary dimensions
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to allow arbitrary dimensions
  virtual void wrap(colvarvalue &x_unwrapped) const;
};


// \brief Colvar component: total value of a scalar map
// (usually implemented as a grid by the simulation engine)
class colvar::map_total
  : public colvar::cvc
{
public:

  map_total();
  virtual ~map_total() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);

protected:

  /// String identifier of the map object (as used by the simulation engine)
  std::string volmap_name;

  /// Numeric identifier of the map object (as used by the simulation engine)
  int volmap_id = -1;

  /// Index of the map objet in the proxy arrays
  int volmap_index = -1;

  /// Group of atoms selected internally (optional)
  cvm::atom_group *atoms = nullptr;

  /// Weights assigned to each atom (default: uniform weights)
  std::vector<cvm::real> atom_weights;
};



#endif
