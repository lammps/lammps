#ifndef ATC_TYPEDEFS_H
#define ATC_TYPEDEFS_H

#include <set>
using std::pair;
using std::set;

#ifdef NEW_LAMMPS
#include "lmptype.h"
#endif

#include "Array.h"
#include "Array2D.h"

using namespace ATC_matrix;

#include "MatrixLibrary.h"
#include "DependencyManager.h"

namespace ATC
{
  /** physical constants */
  static const double kBeV_ = 8.617343e-5;// [eV/K]

  /** unsigned ints, when needed */
  typedef int INDEX; 

  /** elementset integral */
  enum ElementsetOperationType {
    ELEMENTSET_TOTAL=0,
    ELEMENTSET_AVERAGE
  };
  /** faceset integral */
  enum FacesetIntegralType {
    BOUNDARY_INTEGRAL=0,
    CONTOUR_INTEGRAL
  };

  /** nodeset operation */
  enum NodesetOperationType {
    NODESET_SUM=0,
    NODESET_AVERAGE
  };

  /** boundary integration */
  enum BoundaryIntegrationType {
    NO_QUADRATURE=0,
    FE_QUADRATURE,
    FE_INTERPOLATION
  };
  /** domain integration */
  enum IntegrationDomainType {
    FULL_DOMAIN=0,
    ATOM_DOMAIN,
    FE_DOMAIN,
    FULL_DOMAIN_ATOMIC_QUADRATURE_SOURCE 
  };
  /** domain decomposition */
  enum DomainDecompositionType {
    REPLICATED_MEMORY=0,
    DISTRIBUTED_MEMORY 
  };
  /** atomic weight specification */
  enum AtomicWeightType {
    USER=0,
    LATTICE,
    ELEMENT,
    REGION,
    GROUP,
    MULTISCALE,
    NODE,
    NODE_ELEMENT,
    READ_IN
  };
  /** geometry location with respect to MD domain */
  enum GeometryType {
    FE_ONLY = 0,
    MD_ONLY,
    BOUNDARY
  };
  /** enumerated type for atomic reference frame */
  enum AtomToElementMapType {
    LAGRANGIAN=0,
    EULERIAN
  };
  /* enumerated type for coupling matrix structure */
  enum MatrixStructure {
    FULL=0,     // contributions from all nodes
    LOCALIZED,  // contributions only from nodes with sources
    LUMPED      // row-sum lumped version of full matrix
  };
  /* enumerated type for distinguishing ghost from internal atoms */
  enum AtomType {
    INTERNAL=0,
    GHOST,
    ALL,
    PROC_GHOST,
    NO_ATOMS,
    NUM_ATOM_TYPES
  };
  /** field types */
  enum FieldName { 
      TIME=-2,
      POSITION=-1,
      TEMPERATURE=0, // Intrinsic Fields
      DISPLACEMENT,
      VELOCITY,
      MASS_DENSITY,
      CHARGE_DENSITY,
      SPECIES_CONCENTRATION,
      ELECTRON_DENSITY, // Extrinsic Fields
      ELECTRON_VELOCITY,
      ELECTRON_TEMPERATURE,
      ELECTRIC_POTENTIAL,
      ELECTRON_WAVEFUNCTION,
      ELECTRON_WAVEFUNCTIONS, 
      ELECTRON_WAVEFUNCTION_ENERGIES, 
      FERMI_ENERGY, 
      MOMENTUM,
      PROJECTED_VELOCITY,
      KINETIC_TEMPERATURE,
      THERMAL_ENERGY,
      KINETIC_ENERGY,
      STRESS,
      HEAT_FLUX,
      CHARGE_FLUX,
      SPECIES_FLUX,
      INTERNAL_ENERGY,
      REFERENCE_POTENTIAL_ENERGY,
      POTENTIAL_ENERGY,
      ENERGY,
      NUMBER_DENSITY,
      ESHELBY_STRESS,
      CAUCHY_BORN_STRESS,
      CAUCHY_BORN_ENERGY,
      CAUCHY_BORN_ESHELBY_STRESS,
      TRANSFORMED_STRESS,
      VACANCY_CONCENTRATION,
      ROTATION,
      STRETCH,
      DIPOLE_MOMENT,
      QUADRUPOLE_MOMENT,
      CAUCHY_BORN_ELASTIC_DEFORMATION_GRADIENT,
      DISLOCATION_DENSITY,
      NUM_TOTAL_FIELDS 
  };
  const int NUM_FIELDS = ELECTRON_WAVEFUNCTION+1; 

#define NDIM 3
  static const int FieldSizes[NUM_TOTAL_FIELDS] = {
    1, // TEMPERATURE
    NDIM, // DISPLACEMENT
    NDIM, // VELOCITY
    1, // MASS_DENSITY
    1, // CHARGE_DENSITY
    0, // SPECIES_CONCENTRATION - VARIABLE
    1, // ELECTRON_DENSITY
    NDIM, // ELECTRON_VELOCITY
    1, // ELECTRON_TEMPERATURE
    1, // ELECTRIC_POTENTIAL
    1, // ELECTRON_WAVEFUNCTION ?
    0, // ELECTRON_WAVEFUNCTIONS - VARIABLE
    0, // ELECTRON_WAVEFUNCTION_ENERGIES - VARIABLE
    1, // FERMI_ENERGY
    NDIM, // MOMENTUM
    NDIM, // PROJECTED_VELOCITY
    1, // KINETIC_TEMPERATURE
    1, // THERMAL_ENERGY
    1, // KINETIC_ENERGY
    NDIM*NDIM, // STRESS
    NDIM, // HEAT_FLUX
    NDIM, // CHARGE_FLUX
    0, // SPECIES_FLUX - VARIABLE
    1, // INTERNAL_ENERGY
    1, // REFERENCE_POTENTIAL_ENERGY
    1, // POTENTIAL_ENERGY
    1, // ENERGY
    1, // NUMBER_DENSITY
    NDIM*NDIM, // ESHELBY_STRESS
    NDIM*NDIM, // CAUCHY_BORN_STRESS,
    1, //  CAUCHY_BORN_ENERGY,
    NDIM*NDIM, //  CAUCHY_BORN_ESHELBY_STRESS,
    NDIM*NDIM, //  TRANSFORMED_STRESS,
    1, //  VACANCY_CONCENTRATION,
    NDIM*NDIM, //  ROTATION,
    NDIM*NDIM, //  STRETCH,
    NDIM, // DIPOLE_MOMENT,
    NDIM, // QUADRUPOLE_MOMENT,
    NDIM*NDIM, //  CAUCHY_BORN_ELASTIC_DEFORMATION_GRADIENT,
    NDIM*NDIM //  DISLOCATION_DENSITY
  };

  enum hardyNormalization { 
    NO_NORMALIZATION=0,
    VOLUME_NORMALIZATION, NUMBER_NORMALIZATION, MASS_NORMALIZATION
  };

  /** enums for FE Element and Interpolate classes */
  enum FeEltGeometry   {HEXA, TETRA};
  enum FeIntQuadrature {NODAL, GAUSS1, GAUSS2, GAUSS3, FACE};

  /** field name enum to string */
  inline FeIntQuadrature string_to_FIQ(const string &str) 
  {
    if (str == "nodal")
      return NODAL;
    else if (str == "gauss1")
      return GAUSS1;
    else if (str == "gauss2")
      return GAUSS2;
    else if (str == "gauss3")
      return GAUSS3;
    else if (str == "face")
      return FACE;
    else
      throw ATC_Error("Bad quadrature input" + str + ".");
  }

  /** field name enum to string */
  inline string field_to_string(const FieldName index) 
  {
    switch (index) {
    case TEMPERATURE:
      return "temperature";
    case DISPLACEMENT:
      return "displacement";
    case VELOCITY:
      return "velocity";
    case MASS_DENSITY:
      return "mass_density";
    case CHARGE_DENSITY:
      return "charge_density";
    case ELECTRON_DENSITY:
      return "electron_density";
    case ELECTRON_VELOCITY:
      return "electron_velocity";
    case ELECTRON_TEMPERATURE:
      return "electron_temperature";
    case ELECTRIC_POTENTIAL:
      return "electric_potential";
    case ELECTRON_WAVEFUNCTION:
      return "electron_wavefunction";
    case ELECTRON_WAVEFUNCTIONS:
      return "electron_wavefunctions";
    case ELECTRON_WAVEFUNCTION_ENERGIES:
      return "electron_wavefunction_energies";
    case FERMI_ENERGY:
      return "fermi_energy";
    case MOMENTUM:
      return "momentum";
    case PROJECTED_VELOCITY:
      return "projected_velocity";
    case KINETIC_TEMPERATURE:
      return "kinetic_temperature";
    case THERMAL_ENERGY:
      return "thermal_energy";
    case KINETIC_ENERGY:
      return "kinetic_energy";
    case STRESS:
      return "stress";
    case ESHELBY_STRESS:
      return "eshelby_stress";
    case CAUCHY_BORN_STRESS:
      return "cauchy_born_stress";
    case CAUCHY_BORN_ENERGY:
      return "cauchy_born_energy";
    case CAUCHY_BORN_ESHELBY_STRESS:
      return "cauchy_born_eshelby_stress";
    case HEAT_FLUX:
      return "heat_flux";
    case CHARGE_FLUX:
      return "charge_flux";
    case SPECIES_FLUX:
      return "species_flux";
    case INTERNAL_ENERGY:
      return "internal_energy";
    case POTENTIAL_ENERGY:
      return "potential_energy";
    case REFERENCE_POTENTIAL_ENERGY:
      return "reference_potential_energy";
    case ENERGY:
      return "energy";
    case NUMBER_DENSITY:
      return "number_density";
    case TRANSFORMED_STRESS:
      return "transformed_stress";
    case VACANCY_CONCENTRATION:
      return "vacancy_concentration";
    case SPECIES_CONCENTRATION:
      return "species_concentration";
    case ROTATION:
      return "rotation";
    case STRETCH:
      return "stretch";
    case DIPOLE_MOMENT:
      return "dipole_moment";
    case QUADRUPOLE_MOMENT:
      return "quadrupole_moment";
    default:
      throw ATC_Error("field not found in field_to_string");
    }
  };

  /** string to field enum */
  inline FieldName string_to_field(const string & name) 
  {
    if      (name=="temperature")
      return TEMPERATURE;
    else if (name=="displacement")
      return DISPLACEMENT;
    else if (name=="velocity")
      return VELOCITY;
    else if (name=="mass_density")
      return MASS_DENSITY;
    else if (name=="charge_density")
      return CHARGE_DENSITY;
    else if (name=="electron_density")
      return ELECTRON_DENSITY;
    else if (name=="electron_velocity")
      return ELECTRON_VELOCITY;
    else if (name=="electron_temperature")
      return ELECTRON_TEMPERATURE;
    else if (name=="electric_potential")
      return ELECTRIC_POTENTIAL;
    else if (name=="electron_wavefunction")
      return ELECTRON_WAVEFUNCTION;
    else if (name=="electron_wavefunctions")
      return ELECTRON_WAVEFUNCTIONS;
    else if (name=="electron_wavefunction_energies")
      return ELECTRON_WAVEFUNCTION_ENERGIES;
    else if (name=="fermi_energy")
      return FERMI_ENERGY;
    else if (name=="momentum")
      return MOMENTUM;
    else if (name=="projected_velocity")
      return PROJECTED_VELOCITY;
    else if (name=="kinetic_temperature")
      return KINETIC_TEMPERATURE; // temperature from total KE
    else if (name=="thermal_energy")
      return THERMAL_ENERGY;
    else if (name=="kinetic_energy")
      return KINETIC_ENERGY;
    else if (name=="stress")
      return STRESS;
    else if (name=="eshelby_stress")
      return ESHELBY_STRESS;
    else if (name=="cauchy_born_stress")
      return CAUCHY_BORN_STRESS;
    else if (name=="cauchy_born_energy")
      return CAUCHY_BORN_ENERGY;
    else if (name=="cauchy_born_eshelby_stress")
      return CAUCHY_BORN_ESHELBY_STRESS;
    else if (name=="heat_flux")
      return HEAT_FLUX;
    else if (name=="charge_flux")
      return CHARGE_FLUX;
    else if (name=="species_flux")
      return SPECIES_FLUX;
    else if (name=="internal_energy")
      return INTERNAL_ENERGY;
    else if (name=="reference_potential_energy")
      return REFERENCE_POTENTIAL_ENERGY;
    else if (name=="potential_energy")
      return POTENTIAL_ENERGY;
    else if (name=="energy")
      return ENERGY;
    else if (name=="number_density")
      return NUMBER_DENSITY;
    else if (name=="transformed_stress")
      return TRANSFORMED_STRESS;
    else if (name=="vacancy_concentration")
      return VACANCY_CONCENTRATION;
    else if (name=="species_concentration")
      return SPECIES_CONCENTRATION;
    else if (name=="rotation")
      return ROTATION;
    else if (name=="stretch")
      return STRETCH;
    else if (name=="dipole_moment")
      return DIPOLE_MOMENT;
    else if (name=="quadrupole_moment")
      return QUADRUPOLE_MOMENT;
    else
      throw ATC_Error(name + " is not a valid field");
  };

  inline bool is_intrinsic(const FieldName & field_enum) 
  {
    if  (field_enum==TEMPERATURE
      || field_enum==DISPLACEMENT
      || field_enum==VELOCITY
      || field_enum==MASS_DENSITY
      || field_enum==CHARGE_DENSITY
      || field_enum==SPECIES_CONCENTRATION
      || field_enum==KINETIC_TEMPERATURE
      || field_enum==POTENTIAL_ENERGY
      || field_enum==REFERENCE_POTENTIAL_ENERGY
     )   return true;
    else return false;
  };

  inline string field_to_intrinsic_name(const FieldName index) 
  {
    if (is_intrinsic(index)) {
      return "NodalAtomic"+ATC_Utility::to_cap(field_to_string(index));
    }
    else {
      throw ATC_Error("field "+field_to_string(index)+" is not an intrinsic field");
    }
  }
  inline string field_to_restriction_name(const FieldName index) 
  {
    if (is_intrinsic(index)) {
      return "Restricted"+ATC_Utility::to_cap(field_to_string(index));
    }
    else {
      throw ATC_Error("field "+field_to_string(index)+" is not an intrinsic field");
    }
  }
  inline string field_to_prolongation_name(const FieldName index) 
  {
    return "Prolonged"+ATC_Utility::to_cap(field_to_string(index));
  }


  /** types of temperature definitions */
  enum TemperatureDefType {
    NONE = 0,
    KINETIC,
    TOTAL
  };

  /** types of ghost boundary conditions in momentum */
  enum BoundaryDynamicsType {
    NO_BOUNDARY_DYNAMICS=0,
    PRESCRIBED,
    DAMPED_HARMONIC,
    COUPLED
  };

  /** string to temperature definition enum */
  inline bool string_to_temperature_def(const string & name, TemperatureDefType & index) {
    if      (name=="none")
      index = NONE;
    else if (name=="kinetic")
      index = KINETIC;
    else if (name=="total")
      index = TOTAL;
    else {
      throw ATC_Error("temperature operator type "+name+" not valid");
      return false;
    }

    return true;
  };

  /** solver types */
  enum SolverType { DIRECT=0, ITERATIVE};
  enum DirichletType {DIRICHLET_PENALTY=0, DIRICHLET_CONDENSE};

  /** physics types */
  enum PhysicsType
  {
    NO_PHYSICS=0, // for post-processing only
    THERMAL,
    ELASTIC,
    SHEAR,
    THERMO_ELASTIC,
    SPECIES // aka Mass
  };
 
  /** rhs types */
  enum FluxType 
  {
    FLUX = 0,           // has a source weighted by gradient of shape function
    SOURCE,             // has a source term weighted by the shape function
    PRESCRIBED_SOURCE,  // has a prescribed source term
    ROBIN_SOURCE,       // has a Robin source term
    EXTRINSIC_SOURCE,   // has an extrinsic source term
    NUM_FLUX
  };

  /** stiffness/ derivative of rhs types */
  enum StiffnessType 
  {
    BB_STIFFNESS = 0,    
    NN_STIFFNESS,    
    BN_STIFFNESS,    
    NB_STIFFNESS,    
    NUM_STIFFNESS
  };

  /** LAMMPS atom type identifiers */
  enum IdType {
    ATOM_TYPE=0,
    ATOM_GROUP
  };

  /** molecule size identifiers */
  enum MolSize {
    MOL_SMALL=0,
    MOL_LARGE
  };

  /** basic */
  typedef pair<int, int> PAIR;

  /** typedefs for compact set of bc values */
  typedef set < pair < int, double> > BC_SET; // node, value
  typedef vector< BC_SET > BCS;  // dof (node, value)
  typedef set<int> NSET;  // nodeset
  typedef set<PAIR> FSET;  // faceset
  typedef set<int> ESET;  // elementset

  /** typedefs for N and B integrand functions */
  typedef set<FieldName> ARG_NAMES; 
  typedef map<FieldName, DenseMatrix<double> > ARGS; 
  typedef MatrixDependencyManager<DenseMatrix, double>  FIELD;
  typedef vector<MatrixDependencyManager<DenseMatrix, double> >  GRAD_FIELD;
  typedef map<FieldName, MatrixDependencyManager<DenseMatrix, double> > FIELDS;
  typedef map<FieldName, MatrixDependencyManager<DenseMatrix, double> * > FIELD_POINTERS;
  typedef map<FieldName, DenseMatrix<double> > FIELD_MATS;
  typedef map<string, MatrixDependencyManager<DenseMatrix, double> > TAG_FIELDS;
  typedef map<FieldName, vector<MatrixDependencyManager<DenseMatrix, double> > > GRAD_FIELDS;
  typedef map<FieldName, vector<DenseMatrix<double> > > GRAD_FIELD_MATS;
  typedef map<FieldName, MatrixDependencyManager<DiagonalMatrix, double> > MASS_MATS;
  typedef map<FieldName, MatrixDependencyManager<SparseMatrix, double> > CON_MASS_MATS;
  typedef MatrixDependencyManager<DenseMatrix, double> DENS_MAN;
  typedef MatrixDependencyManager<SparseMatrix, double> SPAR_MAN;
  typedef MatrixDependencyManager<ParSparseMatrix, double> PAR_SPAR_MAN;
  typedef MatrixDependencyManager<DiagonalMatrix, double> DIAG_MAN;
  typedef MatrixDependencyManager<ParDiagonalMatrix, double> PAR_DIAG_MAN;

  /** typedefs for WeakEquation evaluation */
  typedef Array2D<bool>  RHS_MASK;

  /** typedefs for input/output */
  typedef map<string, const Matrix<double>*> OUTPUT_LIST;
  typedef map<string, Matrix<double>*> RESTART_LIST;

  typedef  pair<int, int>  ID_PAIR;
  typedef vector< pair<int, int> > ID_LIST;

  /** misc typedefs */
  class XT_Function;
  class UXT_Function;
  typedef map<FieldName, map<PAIR, Array<XT_Function*> > > SURFACE_SOURCE;
  typedef map<FieldName, map<PAIR, Array<UXT_Function*> > > ROBIN_SURFACE_SOURCE;
  typedef map<FieldName, Array2D<XT_Function *> > VOLUME_SOURCE;
  typedef map<string, MatrixDependencyManager<DenseMatrix, double> > ATOMIC_DATA;
  
  /** typedefs for FE_Mesh */
  typedef map<string, set<int > > NODE_SET_MAP;
  typedef map<string, set<int > > ELEMENT_SET_MAP;
  typedef map<string, set<PAIR> > FACE_SET_MAP;

  /** string to index */
  // inline vs. static is used to avoid compiler warnings that the function isn't used
  // the compiler seems to just check if static functions are used in the file they're
  // declared in rather than all the files that include the header,
  // same for arrays (but not primitives, e.g. ints) hopefully this also speeds up the code
  inline bool string_to_index(const string & dim, int & index, int & sgn)
  {
    char dir;
    if (dim.empty()) return false;
    sgn = (dim[0] == '-') ? -1 : 1;
    dir = dim[dim.size()-1]; // dir is last character
    if      (dir == 'x') index = 0;
    else if (dir == 'y') index = 1;
    else if (dir == 'z') index = 2;
    else return false;
    return true;
  };

  /** string to index */
  inline string index_to_string(const int &index)
  {
    if      (index==0) return "x";
    else if (index==1) return "y";
    else if (index==2) return "z";
    return "unknown";
  };

  /** string to index */
  inline bool string_to_index(const string &dim, int &index)
  {
    if (dim=="x")
      index = 0;
    else if (dim=="y")
      index = 1;
    else if (dim=="z")
      index = 2;
    else
      return false;

    return true;
  };

  inline string print_mask(const Array2D<bool> & rhsMask) 
  {
    string msg;
    for (int i = 0; i < NUM_FIELDS; i++) {       
      FieldName field = (FieldName) i;
      string name = field_to_string(field);
      if (rhsMask(field,FLUX)
          || rhsMask(field,SOURCE)           
          || rhsMask(field,PRESCRIBED_SOURCE)
          || rhsMask(field,ROBIN_SOURCE)
          || rhsMask(field,EXTRINSIC_SOURCE))  {
        msg = "RHS_MASK: " + name;
        if (rhsMask(field,FLUX))              msg += " flux";
        if (rhsMask(field,SOURCE))            msg += " source";
        if (rhsMask(field,PRESCRIBED_SOURCE)) msg += " prescribed_src";
        if (rhsMask(field,ROBIN_SOURCE))      msg += " robin_src";
        if (rhsMask(field,EXTRINSIC_SOURCE))  msg += " extrinsic_src";
      }
    }
    return msg;
  }
}

#endif
