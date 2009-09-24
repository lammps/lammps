/**
 *  @class  PhysicsModel
 *  @brief  An adaptor for the FE_Engine of the specific weak form of 
 *  the continuum PDE for the FE_Engine.
 *  It is assumed that the PDE fits this template:
 *  DENSITY(FIELDS) FIELD_RATE 
 *    = DIV FLUX(FIELDS, GRAD_FIELDS) + SOURCE(FIELDS,GRAD_FIELDS)
 *    + PRESCRIBED_SOURCE(X,t) + EXTRINSIC_SOURCE(FIELDS,GRAD_FIELDS)
 */


#ifndef PHYSICS_MODEL_H
#define PHYSICS_MODEL_H


// included headers
#include <map>
#include "Array2D.h"
#include "MatrixLibrary.h"
#include "ATC_Error.h"
#include "Material.h"
#include "ATC_TypeDefs.h"
#include "StringManip.h"

using namespace std;
using namespace ATC_STRING;

namespace ATC 
{

  // Forward declarations
  class ATC_Transfer;

  class PhysicsModel 
  {

  public:

    // constructor 
    PhysicsModel(string fileName,
                 ATC_Transfer * atcTransfer) 
    : atcTransfer_(atcTransfer) 
    { 
      parse_material_file(fileName); 
    }

    // destructor
    virtual ~PhysicsModel()
    {
      vector< Material* >::iterator iter;
      for (iter = materials_.begin(); iter != materials_.end(); iter++) {
        Material * mat = *iter;
        if (mat) delete mat;
      }
    }

    /** parse material file */
    void parse_material_file(string fileName)
    {
      vector< Material* >::iterator iter;
      for (iter = materials_.begin(); iter != materials_.end(); iter++) {
        Material * mat = *iter;
          if (mat) delete mat;
      }
      fstream  fileId(fileName.c_str(), std::ios::in);
      if (!fileId.is_open()) throw ATC_Error(0,"cannot open material file");
      vector<string> line;
      int index = 0;
      while(fileId.good()) {      
        get_command_line(fileId, line);
        if (line.size() == 0 || line[0] == "#") continue;
        if (line[0] == "material") {
          string tag = line[1];
          Material * mat = new Material(tag,fileId);
          materials_.push_back(mat);
          materialNameToIndexMap_[tag] = index++;
        }
      }
      if (int(materials_.size()) == 0) {
        throw ATC_Error(0,"No materials were defined"); }
      cout << " ATC:: " << int(materials_.size()) << " materials defined\n";
      fileId.close();
    }

    /** initialize */
    virtual void initialize(void) = 0; 

    // set timescale parameters based on a given lengthscale
    virtual void set_timescales(const double lengthscale) {};

    /** access number of materials */
    int get_nMaterials(void) const
    {
      return materials_.size();
    }

    /** access material index from name */
    int material_index(const string & name) const
    {
      string tag = name;
      to_lower(tag); // this is an artifact of StringManip parsing
      map<string,int>::const_iterator iter;
      iter = materialNameToIndexMap_.find(tag);
      if (iter ==  materialNameToIndexMap_.end()) {
        throw ATC_Error(0,"No material named "+name+" found");
      }
      int index = iter->second;
      return index;
    }

    /** access to parameter values */
    bool parameter_value(const string& name, double& value, 
                         const int imat = 0) const 
    {
      // search owned parameters
      value = 0.0;
      map<string,double>::const_iterator it = parameterValues_.find(name);
      if (it != parameterValues_.end()) {
          value = it->second;
          return true;
      }
      // interogate material models
      bool found = materials_[imat]->get_parameter(name,value);
      return found;
    }

    /** return fields ids and length */
    virtual void get_num_fields(map<FieldName,int> & fieldSizes, 
                                Array2D<bool> & fieldMask) const = 0;

    /** is the material model linear */
    virtual bool is_linear(FieldName name) { 
      vector< Material* >::iterator iter;
      for (iter = materials_.begin(); iter != materials_.end(); iter++) {
        Material * mat = *iter;
        bool linear = mat->linear_flux(name) 
                   && mat->linear_source(name) 
                   && mat->constant_density(name);
        if (! linear) return linear;
      }
      return true;
    }

    /** is rhs linear */
    virtual bool has_linear_rhs(FieldName name) { 
      vector< Material* >::iterator iter;
      for (iter = materials_.begin(); iter != materials_.end(); iter++) {
        Material * mat = *iter;
        bool constant = mat->linear_flux(name) && mat->linear_source(name);
        if (! constant) return constant;
      }
      return true;
    }

    /** is mass matrix constant */
    virtual bool has_constant_mass(FieldName name) { 
      vector< Material* >::iterator iter;
      for (iter = materials_.begin(); iter != materials_.end(); iter++) {
        Material * mat = *iter;
        bool constant = mat->constant_density(name);
        if (! constant) return constant;
      }
      return true;
    }

    /** energy or other preserved quantity */
    virtual void E_integrand(const Array<FieldName> &mask, 
                             const FIELDS &fields, 
                             const GRAD_FIELDS &grad_fields,
                             FIELDS &capacity,
                             const int matIndex = 0)  const
    {
      throw ATC_Error(0,"E_integrand not implemented for this PhysicsModel");
    }

    /** heat/momentum/energy/mass capacity used in the LHS mass matrix */
    virtual void M_integrand(const Array<FieldName> &mask, 
                             const FIELDS &fields, 
                             FIELDS &capacity,
                             const int matIndex = 0)  const
    {
      throw ATC_Error(0,"M_integrand not implemented for this PhysicsModel");
    }
    // flux that is integrated with N as its weight 
    virtual void N_integrand(const Array2D<bool> &mask, 
                             const FIELDS &fields, 
                             const GRAD_FIELDS &grad_fields,
                             FIELDS &flux,
                             const int matIndex = 0)  const
    {
      throw ATC_Error(0,"N_integrand not implemented for this PhysicsModel");
    }

    /** flux that is integrated with Grad N as its weight */
    virtual void B_integrand(const Array2D<bool> & mask, 
                             const FIELDS &fields,
                             const GRAD_FIELDS &grad_fields,
                             GRAD_FIELDS &flux,
                             const int matIndex = 0) const
    {
      throw ATC_Error(0,"B_integrand not implemented for this PhysicsModel");
    }

    /** has a integrand for the N weighted integral */
    virtual bool has_N_integrand() const { return false; }

    /** has a integrand for the B=grad_x N weighted integral */
    virtual bool has_B_integrand() const { return false; }

  protected:

    /** associated ATC Transfer object */
    ATC_Transfer * atcTransfer_;

    // parameter values
    map<string, double> parameterValues_;

    // material models
    vector<Material *> materials_;
    map<string,int> materialNameToIndexMap_;
  };

};
#endif
