#ifndef ATC_TYPEDEFS_H
#define ATC_TYPEDEFS_H

#include <set>
using std::pair;
using std::set;

#include "Array.h"
#include "Array2D.h"
#include "XT_Function.h"
#include "MatrixLibrary.h"
#include "ATC_Error.h"

namespace ATC
{
  /** Material types */
  enum FieldName { // Intrinsic Fields
      TEMPERATURE=0,
      DISPLACEMENT,
      VELOCITY,
      MASS_DENSITY,
      CHARGE_DENSITY,
      ELECTRON_DENSITY, // Extrinsic Fields
      ELECTRON_VELOCITY,
      ELECTRON_TEMPERATURE,
      ELECTRIC_POTENTIAL,
      NUM_FIELDS
  };

  /** solver types */
  enum SolverType { DIRECT=0, ITERATIVE}; 

  /** physics types */
  enum PhysicsType 
  {
    NO_PHYSICS=0, // for post-processing only
    THERMAL,
    ELASTIC,
    THERMO_ELASTIC,
    SPECIES
  };
 
  /** rhs types */
  enum FluxType 
  {
    FLUX = 0,           // has a source weighted by gradient of shape function
    SOURCE,             // has a source term weighted by the shape function
    PRESCRIBED_SOURCE,  // has a prescribed source term
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

  /** typedefs for N and B integrand functions */
  typedef DenseMatrix<double>  FIELD;
  typedef vector<DenseMatrix<double> >  GRAD_FIELD;
  typedef map<FieldName, DenseMatrix<double> > FIELDS;
  typedef map<FieldName, vector<DenseMatrix<double> > > GRAD_FIELDS;

  /** typedefs for input/output */
  typedef map<string, Matrix<double>*> OUTPUT_LIST;

  /** misc typedefs */
  typedef pair<int, int> PAIR;
  typedef map<FieldName, map<PAIR, Array<XT_Function*> > > SURFACE_SOURCE;
  typedef map<FieldName, Array2D<XT_Function *> > VOLUME_SOURCE;
  typedef vector<SparseMatrix<double> > GRAD_SHPFCN;

  /** typedefs for FE_Mesh */
  typedef map<string, set<int > > NODE_SET_MAP;
  typedef map<string, set<int > > ELEMENT_SET_MAP;
  typedef map<string, set<PAIR> > FACE_SET_MAP;

  /** string to index */
  static bool string_to_index(const string & dim, int & index, int & sgn)
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
  static string index_to_string(const int &index)
  {
    if      (index==0) return "x";
    else if (index==1) return "y";
    else if (index==2) return "z";
    return "unknown";
  };

  /** string to index */
  static bool string_to_index(const string &dim, int &index)
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


  /** field name enum to string */
  static string field_to_string(const FieldName index) 
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
    default:
      throw ATC_Error(0,"field not found in field_to_string");
    }
  };

  /** string to field enum */
  static FieldName string_to_field(const string & name) 
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
    else
      throw ATC_Error(0,name + " is not a valid field");
  };

  static bool is_intrinsic(const FieldName & field_enum) 
  {
    if  (field_enum==TEMPERATURE
      || field_enum==DISPLACEMENT
      || field_enum==VELOCITY
      || field_enum==MASS_DENSITY
      || field_enum==CHARGE_DENSITY)  return true;
    else                        return false;
  };


  static void print_mask(const Array2D<bool> & rhsMask) 
  {
    for (int i = 0; i < NUM_FIELDS; i++) {       
      FieldName field = (FieldName) i;
      string name = field_to_string(field);
      if (rhsMask(field,FLUX)
          || rhsMask(field,SOURCE)           
          || rhsMask(field,PRESCRIBED_SOURCE)
          || rhsMask(field,EXTRINSIC_SOURCE))  {
        cout << "RHS_MASK: " << name;
        if (rhsMask(field,FLUX))  cout << " flux";
        if (rhsMask(field,SOURCE)) cout << " source";
        if (rhsMask(field,PRESCRIBED_SOURCE)) cout << " prescribed_src";
        if (rhsMask(field,EXTRINSIC_SOURCE)) cout << " extrinsic_src";
        cout << "\n";
      }
    }
  }


}

#endif
