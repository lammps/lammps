// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
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


#include "colvarmodule.h"
#include "colvar.h"
#include "colvaratoms.h"


/// \brief Colvar component (base class for collective variables)
///
/// A \link cvc \endlink object (or an object of a
/// cvc-derived class) implements the calculation of a collective
/// variable, its gradients and any other related physical quantities
/// that depend on microscopic degrees of freedom.
///
/// No restriction is set to what kind of calculation a \link cvc \endlink
/// object performs (usually an analytical function of atomic coordinates).
/// The only constraints are that: \par
///
/// - The value is calculated by the \link calc_value() \endlink
///   method, and is an object of \link colvarvalue \endlink class.  This
///   provides a transparent way to treat scalar and non-scalar variables
///   alike, and allows an automatic selection of the applicable algorithms.
///
/// - The object provides an implementation \link apply_force() \endlink to
///   apply forces to atoms.  Typically, one or more \link cvm::atom_group
///   \endlink objects are used, but this is not a requirement for as long as
///   the \link cvc \endlink object communicates with the simulation program.
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

  /// \brief Description of the type of collective variable
  ///
  /// Normally this string is set by the parent \link colvar \endlink
  /// object within its constructor, when all \link cvc \endlink
  /// objects are initialized; therefore the main "config string"
  /// constructor does not need to define it.  If a \link cvc
  /// \endlink is initialized and/or a different constructor is used,
  /// this variable definition should be set within the constructor.
  std::string function_type;

  /// Keyword used in the input to denote this CVC
  std::string config_key;

  /// \brief Coefficient in the polynomial combination (default: 1.0)
  cvm::real sup_coeff;
  /// \brief Exponent in the polynomial combination (default: 1)
  int       sup_np;

  /// \brief Is this a periodic component?
  bool b_periodic;

  /// \brief Period of this cvc value, (default: 0.0, non periodic)
  cvm::real period;

  /// \brief If the component is periodic, wrap around this value (default: 0.0)
  cvm::real wrap_center;

  /// \brief Constructor
  ///
  /// Calls the init() function of the class
  cvc(std::string const &conf);

  /// An init function should be defined for every class inheriting from cvc
  /// \param conf Contents of the configuration file pertaining to this \link
  /// cvc \endlink
  virtual int init(std::string const &conf);

  /// \brief Within the constructor, make a group parse its own
  /// options from the provided configuration string
  /// Returns reference to new group
  cvm::atom_group *parse_group(std::string const &conf,
                               char const *group_key,
                               bool optional = false);

  /// \brief Parse options pertaining to total force calculation
  virtual int init_total_force_params(std::string const &conf);

  /// \brief After construction, set data related to dependency handling
  int setup();

  /// \brief Default constructor (used when \link cvc \endlink
  /// objects are declared within other ones)
  cvc();

  /// Destructor
  virtual ~cvc();

  /// \brief Implementation of the feature list for colvar
  static std::vector<feature *> cvc_features;

  /// \brief Implementation of the feature list accessor for colvar
  virtual const std::vector<feature *> &features()
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
  virtual void apply_force(colvarvalue const &cvforce) = 0;

  /// \brief Square distance between x1 and x2 (can be redefined to
  /// transparently implement constraints, symmetries and
  /// periodicities)
  ///
  /// colvar::cvc::dist2() and the related functions are
  /// declared as "const" functions, but not "static", because
  /// additional parameters defining the metrics (e.g. the
  /// periodicity) may be specific to each colvar::cvc object.
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
  /// Finally, another useful application, if you are performing very
  /// many operations with these functions, could be to override the
  /// \link colvarvalue \endlink member functions and access directly
  /// its member data.  For instance: to define dist2(x1,x2) as
  /// (x2.real_value-x1.real_value)*(x2.real_value-x1.real_value) in
  /// case of a scalar \link colvarvalue \endlink type.
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
  virtual void wrap(colvarvalue &x) const;

  /// \brief Pointers to all atom groups, to let colvars collect info
  /// e.g. atomic gradients
  std::vector<cvm::atom_group *> atom_groups;

  /// \brief Store a pointer to new atom group, and list as child for dependencies
  inline void register_atom_group(cvm::atom_group *ag) {
    atom_groups.push_back(ag);
    add_child((colvardeps *) ag);
  }

  /// \brief Whether or not this CVC will be computed in parallel whenever possible
  bool b_try_scalable;

protected:

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
  cvm::atom_group  *group1;
  /// Second atom group
  cvm::atom_group  *group2;
  /// Vector distance, cached to be recycled
  cvm::rvector     dist_v;
  /// Use absolute positions, ignoring PBCs when present
  bool b_no_PBC;
public:
  distance(std::string const &conf);
  distance();
  virtual ~distance() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



// \brief Colvar component: distance vector between centers of mass
// of two groups (\link colvarvalue::type_3vector \endlink type,
// range (-*:*)x(-*:*)x(-*:*))
class colvar::distance_vec
  : public colvar::distance
{
public:
  distance_vec(std::string const &conf);
  distance_vec();
  virtual ~distance_vec() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  /// Redefined to handle the box periodicity
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  /// Redefined to handle the box periodicity
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to handle the box periodicity
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: distance unit vector (direction) between
/// centers of mass of two groups (colvarvalue::type_unit3vector type,
/// range [-1:1]x[-1:1]x[-1:1])
class colvar::distance_dir
  : public colvar::distance
{
public:
  distance_dir(std::string const &conf);
  distance_dir();
  virtual ~distance_dir() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  /// Redefined to override the distance ones
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  /// Redefined to override the distance ones
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to override the distance ones
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: projection of the distance vector along
/// an axis(colvarvalue::type_scalar type, range (-*:*))
class colvar::distance_z
  : public colvar::cvc
{
protected:
  /// Main atom group
  cvm::atom_group  *main;
  /// Reference atom group
  cvm::atom_group  *ref1;
  /// Optional, second ref atom group
  cvm::atom_group  *ref2;
  /// Use absolute positions, ignoring PBCs when present
  bool b_no_PBC;
  /// Vector on which the distance vector is projected
  cvm::rvector axis;
  /// Norm of the axis
  cvm::real axis_norm;
  /// Vector distance, cached to be recycled
  cvm::rvector     dist_v;
  /// Flag: using a fixed axis vector?
  bool fixed_axis;
public:
  distance_z(std::string const &conf);
  distance_z();
  virtual ~distance_z() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// \brief Redefined to make use of the user-provided period
  virtual void wrap(colvarvalue &x) const;
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
  distance_xy(std::string const &conf);
  distance_xy();
  virtual ~distance_xy() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};


/// \brief Colvar component: polar coordinate phi of a group
/// (colvarvalue::type_scalar type, range [-180:180])
class colvar::polar_phi
  : public colvar::cvc
{
public:
  polar_phi(std::string const &conf);
  polar_phi();
  virtual ~polar_phi() {}
protected:
  cvm::atom_group  *atoms;
  cvm::real r, theta, phi;
public:
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  /// Redefined to handle the 2*PI periodicity
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual void wrap(colvarvalue &x) const;
};


/// \brief Colvar component: polar coordinate theta of a group
/// (colvarvalue::type_scalar type, range [0:180])
class colvar::polar_theta
  : public colvar::cvc
{
public:
  polar_theta(std::string const &conf);
  polar_theta();
  virtual ~polar_theta() {}
protected:
  cvm::atom_group  *atoms;
  cvm::real r, theta, phi;
public:
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  /// Redefined to override the distance ones
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  /// Redefined to override the distance ones
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to override the distance ones
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};

/// \brief Colvar component: average distance between two groups of atoms, weighted as the sixth power,
/// as in NMR refinements(colvarvalue::type_scalar type, range (0:*))
class colvar::distance_inv
  : public colvar::cvc
{
protected:
  /// First atom group
  cvm::atom_group  *group1;
  /// Second atom group
  cvm::atom_group  *group2;
  /// Components of the distance vector orthogonal to the axis
  int exponent;
  /// Use absolute positions, ignoring PBCs when present
  bool b_no_PBC;
public:
  distance_inv(std::string const &conf);
  distance_inv();
  virtual ~distance_inv() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: N1xN2 vector of pairwise distances
/// (colvarvalue::type_vector type, range (0:*) for each component)
class colvar::distance_pairs
  : public colvar::cvc
{
protected:
  /// First atom group
  cvm::atom_group  *group1;
  /// Second atom group
  cvm::atom_group  *group2;
  /// Use absolute positions, ignoring PBCs when present
  bool b_no_PBC;
public:
  distance_pairs(std::string const &conf);
  distance_pairs();
  virtual ~distance_pairs() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
};



/// \brief Colvar component: Radius of gyration of an atom group
/// (colvarvalue::type_scalar type, range [0:*))
class colvar::gyration
  : public colvar::cvc
{
protected:
  /// Atoms involved
  cvm::atom_group  *atoms;
public:
  /// Constructor
  gyration(std::string const &conf);
  gyration();
  virtual ~gyration() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: moment of inertia of an atom group
/// (colvarvalue::type_scalar type, range [0:*))
class colvar::inertia
  : public colvar::gyration
{
public:
  /// Constructor
  inertia(std::string const &conf);
  inertia();
  virtual ~inertia() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
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
  /// Constructor
  inertia_z(std::string const &conf);
  inertia_z();
  virtual ~inertia_z() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: projection of 3N coordinates onto an
/// eigenvector(colvarvalue::type_scalar type, range (-*:*))
class colvar::eigenvector
  : public colvar::cvc
{
protected:

  /// Atom group
  cvm::atom_group  *           atoms;

  /// Reference coordinates
  std::vector<cvm::atom_pos>  ref_pos;

  /// Geometric center of the reference coordinates
  cvm::atom_pos                ref_pos_center;

  /// Eigenvector (of a normal or essential mode): will always have zero center
  std::vector<cvm::rvector>   eigenvec;

  /// Inverse square norm of the eigenvector
  cvm::real                   eigenvec_invnorm2;

public:

  /// Constructor
  eigenvector(std::string const &conf);
  virtual ~eigenvector() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: angle between the centers of mass of
/// three groups (colvarvalue::type_scalar type, range [0:PI])
class colvar::angle
  : public colvar::cvc
{
protected:

  /// Atom group
  cvm::atom_group  *group1;
  /// Atom group
  cvm::atom_group  *group2;
  /// Atom group
  cvm::atom_group  *group3;

  /// Inter site vectors
  cvm::rvector r21, r23;
  /// Inter site vector norms
  cvm::real r21l, r23l;
  /// Derivatives wrt group centers of mass
  cvm::rvector dxdr1, dxdr3;

  /// Compute total force on first site only to avoid unwanted
  /// coupling to other colvars (see e.g. Ciccotti et al., 2005)
  /// (or to allow dummy atoms)
  bool b_1site_force;
public:

  /// Initialize by parsing the configuration
  angle(std::string const &conf);
  /// \brief Initialize the three groups after three atoms
  angle(cvm::atom const &a1, cvm::atom const &a2, cvm::atom const &a3);
  angle();
  virtual ~angle() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: angle between the dipole of a molecule and an axis
/// formed by two groups of atoms(colvarvalue::type_scalar type, range [0:PI])
class colvar::dipole_angle
  : public colvar::cvc
{
protected:

  /// Dipole atom group
  cvm::atom_group  *group1;
  /// Atom group
  cvm::atom_group  *group2;
  /// Atom group
  cvm::atom_group  *group3;

  /// Inter site vectors
  cvm::rvector r21, r23;
  /// Inter site vector norms
  cvm::real r21l, r23l;
  /// Derivatives wrt group centers of mass
  cvm::rvector dxdr1, dxdr3;

  /// Compute total force on first site only to avoid unwanted
  /// coupling to other colvars (see e.g. Ciccotti et al., 2005)
  /// (or to allow dummy atoms)
  bool b_1site_force;
public:

  /// Initialize by parsing the configuration
  dipole_angle (std::string const &conf);
  /// \brief Initialize the three groups after three atoms
  dipole_angle (cvm::atom const &a1, cvm::atom const &a2, cvm::atom const &a3);
  dipole_angle();
  virtual ~dipole_angle() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force (colvarvalue const &force);
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
};



/// \brief Colvar component: dihedral between the centers of mass of
/// four groups (colvarvalue::type_scalar type, range [-PI:PI])
class colvar::dihedral
  : public colvar::cvc
{
protected:

  /// Atom group
  cvm::atom_group  *group1;
  /// Atom group
  cvm::atom_group  *group2;
  /// Atom group
  cvm::atom_group  *group3;
  /// Atom group
  cvm::atom_group  *group4;
  /// Inter site vectors
  cvm::rvector r12, r23, r34;

  /// \brief Compute total force on first site only to avoid unwanted
  /// coupling to other colvars (see e.g. Ciccotti et al., 2005)
  bool b_1site_force;

public:

  /// Initialize by parsing the configuration
  dihedral(std::string  const &conf);
  /// \brief Initialize the four groups after four atoms
  dihedral(cvm::atom const &a1, cvm::atom const &a2, cvm::atom const &a3, cvm::atom const &a4);
  dihedral();
  virtual ~dihedral() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force(colvarvalue const &force);

  /// Redefined to handle the 2*PI periodicity
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual void wrap(colvarvalue &x) const;
};



/// \brief Colvar component: coordination number between two groups
/// (colvarvalue::type_scalar type, range [0:N1*N2])
class colvar::coordnum
  : public colvar::cvc
{
protected:
  /// First atom group
  cvm::atom_group  *group1;
  /// Second atom group
  cvm::atom_group  *group2;
  /// \brief "Cutoff" for isotropic calculation (default)
  cvm::real     r0;
  /// \brief "Cutoff vector" for anisotropic calculation
  cvm::rvector  r0_vec;
  /// \brief Whether r/r0 or \vec{r}*\vec{1/r0_vec} should be used
  bool b_anisotropic;
  /// Integer exponent of the function numerator
  int en;
  /// Integer exponent of the function denominator
  int ed;

  /// \brief If true, group2 will be treated as a single atom, stored in this
  /// accessory group
  cvm::atom_group *group2_center;

  /// Tolerance for the pair list
  cvm::real tolerance;

  /// Frequency of update of the pair list
  int pairlist_freq;

  /// Pair list
  bool *pairlist;

public:

  coordnum(std::string const &conf);
  coordnum();
  ~coordnum();

  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;

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

  /// Main workhorse function
  template<int flags> int compute_coordnum();

};



/// \brief Colvar component: self-coordination number within a group
/// (colvarvalue::type_scalar type, range [0:N*(N-1)/2])
class colvar::selfcoordnum
  : public colvar::cvc
{
protected:

  /// Selected atoms
  cvm::atom_group  *group1;
  /// \brief "Cutoff" for isotropic calculation (default)
  cvm::real     r0;
  /// Integer exponent of the function numerator
  int en;
  /// Integer exponent of the function denominator
  int ed;
  cvm::real tolerance;
  int pairlist_freq;
  bool *pairlist;

public:

  selfcoordnum(std::string const &conf);
  selfcoordnum();
  ~selfcoordnum();
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);

  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;

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
  bool b_anisotropic;
  /// Integer exponent of the function numerator
  int en;
  /// Integer exponent of the function denominator
  int ed;
public:
  /// Constructor
  groupcoordnum(std::string const &conf);
  groupcoordnum();
  virtual ~groupcoordnum() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);

  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
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
  int en;
  /// Integer exponent of the function denominator
  int ed;
public:
  h_bond(std::string const &conf);
  /// Constructor for atoms already allocated
  h_bond(cvm::atom const &acceptor,
         cvm::atom const &donor,
         cvm::real r0, int en, int ed);
  h_bond();
  virtual ~h_bond();
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);

  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: alpha helix content of a contiguous
/// segment of 5 or more residues, implemented as a sum of phi/psi
/// dihedral angles and hydrogen bonds (colvarvalue::type_scalar type,
/// range [0:1])
// class colvar::alpha_dihedrals
//   : public colvar::cvc
// {
// protected:

//   /// Alpha-helical reference phi value
//   cvm::real phi_ref;

//   /// Alpha-helical reference psi value
//   cvm::real psi_ref;

//   /// List of phi dihedral angles
//   std::vector<dihedral *> phi;

//   /// List of psi dihedral angles
//   std::vector<dihedral *> psi;

//   /// List of hydrogen bonds
//   std::vector<h_bond *>   hb;

// public:

//   alpha_dihedrals (std::string const &conf);
//   alpha_dihedrals();
//   virtual ~alpha_dihedrals() {}
//   virtual void calc_value();
//   virtual void calc_gradients();
//   virtual void apply_force (colvarvalue const &force);
//   virtual cvm::real dist2 (colvarvalue const &x1,
//                            colvarvalue const &x2) const;
//   virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
//                                    colvarvalue const &x2) const;
//   virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
//                                    colvarvalue const &x2) const;
// };



/// \brief Colvar component: alpha helix content of a contiguous
/// segment of 5 or more residues, implemented as a sum of Ca-Ca-Ca
/// angles and hydrogen bonds (colvarvalue::type_scalar type, range
/// [0:1])
class colvar::alpha_angles
  : public colvar::cvc
{
protected:

  /// Reference Calpha-Calpha angle (default: 88 degrees)
  cvm::real theta_ref;

  /// Tolerance on the Calpha-Calpha angle
  cvm::real theta_tol;

  /// List of Calpha-Calpha angles
  std::vector<angle *> theta;

  /// List of hydrogen bonds
  std::vector<h_bond *>   hb;

  /// Contribution of the hb terms
  cvm::real hb_coeff;

public:

  alpha_angles(std::string const &conf);
  alpha_angles();
  virtual ~alpha_angles();
  void calc_value();
  void calc_gradients();
  void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
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

  dihedPC(std::string const &conf);
  dihedPC();
  virtual  ~dihedPC();
  void calc_value();
  void calc_gradients();
  void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
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
  cvm::atom_group  *          atoms;
  /// Center of geometry of the group
  cvm::atom_pos              atoms_cog;

  /// Reference coordinates
  std::vector<cvm::atom_pos> ref_pos;

  /// Rotation object
  cvm::rotation              rot;

  /// \brief This is used to remove jumps in the sign of the
  /// quaternion, which may be annoying in the colvars trajectory
  cvm::quaternion            ref_quat;

public:

  orientation(std::string const &conf);
  orientation();
  virtual ~orientation() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: angle of rotation with respect to a set
/// of reference coordinates (colvarvalue::type_scalar type, range
/// [0:PI))
class colvar::orientation_angle
  : public colvar::orientation
{
public:

  orientation_angle(std::string const &conf);
  orientation_angle();
  virtual ~orientation_angle() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: cosine of the angle of rotation with respect to a set
/// of reference coordinates (colvarvalue::type_scalar type, range
/// [-1:1])
class colvar::orientation_proj
  : public colvar::orientation
{
public:

  orientation_proj(std::string const &conf);
  orientation_proj();
  virtual ~orientation_proj() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: projection of the orientation vector onto
/// a predefined axis (colvarvalue::type_scalar type, range [-1:1])
class colvar::tilt
  : public colvar::orientation
{
protected:

  cvm::rvector axis;

public:

  tilt(std::string const &conf);
  tilt();
  virtual ~tilt() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



/// \brief Colvar component: angle of rotation around a predefined
/// axis (colvarvalue::type_scalar type, range [-PI:PI])
class colvar::spin_angle
  : public colvar::orientation
{
protected:

  cvm::rvector axis;

public:

  spin_angle(std::string const &conf);
  spin_angle();
  virtual ~spin_angle() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  /// Redefined to handle the 2*PI periodicity
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual void wrap(colvarvalue &x) const;
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
  cvm::atom_group  *atoms;

  /// Reference coordinates (for RMSD calculation only)
  std::vector<cvm::atom_pos>  ref_pos;

public:

  /// Constructor
  rmsd(std::string const &conf);
  virtual ~rmsd() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};



// \brief Colvar component: flat vector of Cartesian coordinates
// Mostly useful to compute scripted colvar values
class colvar::cartesian
  : public colvar::cvc
{
protected:
  /// Atom group
  cvm::atom_group  *atoms;
  /// Which Cartesian coordinates to include
  std::vector<size_t> axes;
public:
  cartesian(std::string const &conf);
  cartesian();
  virtual ~cartesian() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
};


// metrics functions for cvc implementations

// simple definitions of the distance functions; these are useful only
// for optimization (the type check performed in the default
// colvarcomp functions is skipped)

// definitions assuming the scalar type

#define simple_scalar_dist_functions(TYPE)                              \
                                                                        \
                                                                        \
  cvm::real colvar::TYPE::dist2(colvarvalue const &x1,                  \
                                colvarvalue const &x2) const            \
  {                                                                     \
    return (x1.real_value - x2.real_value)*(x1.real_value - x2.real_value); \
  }                                                                     \
                                                                        \
                                                                        \
  colvarvalue colvar::TYPE::dist2_lgrad(colvarvalue const &x1,          \
                                        colvarvalue const &x2) const    \
  {                                                                     \
    return 2.0 * (x1.real_value - x2.real_value);                       \
  }                                                                     \
                                                                        \
                                                                        \
  colvarvalue colvar::TYPE::dist2_rgrad(colvarvalue const &x1,          \
                                        colvarvalue const &x2) const    \
  {                                                                     \
    return this->dist2_lgrad(x2, x1);                                   \
  }                                                                     \
                                                                        \

#endif
