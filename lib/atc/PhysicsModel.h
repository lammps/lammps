#ifndef PHYSICS_MODEL_H
#define PHYSICS_MODEL_H

#include <map>
#include <vector>
#include <string>
#include "Array2D.h"
#include "MatrixLibrary.h"
#include "ATC_Error.h"
#include "Material.h"
#include "WeakEquation.h"
#include "NonLinearSolver.h"
#include "ATC_TypeDefs.h"
#include "Utility.h"

namespace ATC 
{
  class ATC_Coupling;
  //-------------------------------------------------------------------
  // @class PhysicsModel
  //-------------------------------------------------------------------
/**
 *  @brief  An adaptor for the FE_Engine of the specific weak form of 
 *  the continuum PDE for the FE_Engine.
 *  It is assumed that the PDE fits this template:
 *  DENSITY(FIELDS) FIELD_RATE 
 *    = DIV FLUX(FIELDS, GRAD_FIELDS) + SOURCE(FIELDS,GRAD_FIELDS)
 *    + PRESCRIBED_SOURCE(X,t) + EXTRINSIC_SOURCE(FIELDS,GRAD_FIELDS)
 *  Also it is important to understand that the physics model only handles
 *  extrinsic fields or surrogates of intrinsic fields
 */
  class PhysicsModel 
  {

  public:

    // constructor 
    PhysicsModel(std::string fileName);

    // destructor
    virtual ~PhysicsModel();

    /** parse material file */
    void parse_material_file(std::string fileName);

    /** initialize */
    void initialize(void); 

    // set timescale parameters based on a given lengthscale
    virtual void set_timescales(const double lengthscale) {};

    /** access number of materials */
    int nMaterials(void) const { return materials_.size(); }

    /** access material index from name */
    int material_index(const std::string & name) const;

    /** access material from index */
    const Material * material(const int index) const {return materials_[index];}

    /** access to parameter values */
    bool parameter_value(const std::string& name, double& value, 
                         const int imat = 0) const ;

    /** return fields ids and length */
    void num_fields(std::map<FieldName,int> & fieldSizes, 
                        Array2D<bool> & rhsMask) const;

    /** is the material model linear */
    bool is_linear(FieldName name) const { 
      std::vector< Material* >::const_iterator iter;
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
    bool has_linear_rhs(FieldName name) const { 
      std::vector< Material* >::const_iterator iter;
      for (iter = materials_.begin(); iter != materials_.end(); iter++) {
        Material * mat = *iter;
        bool constant = mat->linear_flux(name) && mat->linear_source(name);
        if (! constant) return constant;
      }
      return true;
    }

    /** is mass matrix constant */
    bool has_constant_mass(FieldName name) const { 
      std::vector< Material* >::const_iterator iter;
      for (iter = materials_.begin(); iter != materials_.end(); iter++) {
        Material * mat = *iter;
        bool constant = mat->constant_density(name);
        if (! constant) return constant;
      }
      return true;
    }

    /** access to weak equations */
    const WeakEquation * weak_equation(FieldName field) const
    { 
      std::map<FieldName,WeakEquation *>::const_iterator itr = weakEqns_.find(field);
      if (itr == weakEqns_.end()) return NULL;
      return (weakEqns_.find(field))->second;
    }

    /** requires ics */
    bool is_dynamic(FieldName field) const 
    {
      return (weak_equation(field)->type() == WeakEquation::DYNAMIC_PDE);
    }

    /** query null weak equations per material */
    bool null(FieldName field, int matID) const
    { 
      return null_(field,matID);
    }

  protected:
    /** parameter values */
    std::map<std::string, double> parameterValues_;

    /** material models */
    std::vector<Material *> materials_; 
    std::map<std::string,int> materialNameToIndexMap_;// maps tag to index

    /** weak equations */
    std::map<FieldName,WeakEquation *> weakEqns_;

    /** null weak equations per material */
    Array2D<int> null_;

    /** type tag */
    std::string type_;
  };


  // note that these classes do not use inheritance other than from the
  // generic base class above. Inheritance is meant to come from the
  // weak equations that they contain

  //-------------------------------------------------------------------
  // @class   PhysicsModelThermal
  //-------------------------------------------------------------------
  class PhysicsModelThermal : public PhysicsModel
  {
  public: 
    PhysicsModelThermal(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelElastic
  //-------------------------------------------------------------------
  class PhysicsModelElastic : public PhysicsModel
  {
  public: 
    PhysicsModelElastic(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelThemoMechanical
  //-------------------------------------------------------------------
  class PhysicsModelThermoElastic : public PhysicsModel
  {
  public: 
    PhysicsModelThermoElastic(std::string filename);
  };
  
  //-------------------------------------------------------------------
  // @class   PhysicsModelShear
  //-------------------------------------------------------------------
  class PhysicsModelShear : public PhysicsModel
  {
  public: 
    PhysicsModelShear(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelThemoShear
  //-------------------------------------------------------------------
  class PhysicsModelThermoShear : public PhysicsModel
  {
  public: 
    PhysicsModelThermoShear(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelSpecies
  //-------------------------------------------------------------------
  class PhysicsModelSpecies : public PhysicsModel
  {
  public: 
    PhysicsModelSpecies(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelTwoTemperature
  //-------------------------------------------------------------------
  class PhysicsModelTwoTemperature : public PhysicsModel
  {
  public: 
    PhysicsModelTwoTemperature(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelDriftDiffusion
  //-------------------------------------------------------------------
  class PhysicsModelDriftDiffusion : public PhysicsModel
  {
  public: 
    PhysicsModelDriftDiffusion(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelDriftDiffusionEquilibrium
  //-------------------------------------------------------------------
  class PhysicsModelDriftDiffusionEquilibrium : public PhysicsModel
  {
  public: 
    PhysicsModelDriftDiffusionEquilibrium(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelDriftDiffusionSchrodinger
  //-------------------------------------------------------------------
  class PhysicsModelDriftDiffusionSchrodinger : public PhysicsModel
  {
  public: 
    PhysicsModelDriftDiffusionSchrodinger(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelDriftDiffusionConvection
  //-------------------------------------------------------------------
  class PhysicsModelDriftDiffusionConvection : public PhysicsModel
  {
  public: 
    PhysicsModelDriftDiffusionConvection(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelDriftDiffusionEquilibrium
  //-------------------------------------------------------------------
  class PhysicsModelDriftDiffusionConvectionEquilibrium : public PhysicsModel
  {
  public: 
    PhysicsModelDriftDiffusionConvectionEquilibrium(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelDriftDiffusionSchrodinger
  //-------------------------------------------------------------------
  class PhysicsModelDriftDiffusionConvectionSchrodinger : public PhysicsModel
  {
  public: 
    PhysicsModelDriftDiffusionConvectionSchrodinger(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelDriftDiffusionSchrodingerSlice
  //-------------------------------------------------------------------
  class PhysicsModelDriftDiffusionSchrodingerSlice : public PhysicsModel
  {
  public: 
    PhysicsModelDriftDiffusionSchrodingerSlice(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelElectrostatic
  //-------------------------------------------------------------------
  class PhysicsModelElectrostatic : public PhysicsModel
  {
  public: 
    PhysicsModelElectrostatic(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelElectrostaticEquilibrium
  //-------------------------------------------------------------------
  class PhysicsModelElectrostaticEquilibrium : public PhysicsModel
  {
  public: 
    PhysicsModelElectrostaticEquilibrium(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class   PhysicsModelElectrostatic
  //-------------------------------------------------------------------
  class PhysicsModelSpeciesElectrostatic : public PhysicsModel
  {
  public: 
    PhysicsModelSpeciesElectrostatic(std::string filename);
  };
  //-------------------------------------------------------------------
  // @class PhysicsModelTangentOperator
  // @brief adaptor to NonLinearSolver to solve rhs(x,y) = 0 for x
  //-------------------------------------------------------------------
  class PhysicsModelTangentOperator :  public TangentOperator
  {
    public:
      PhysicsModelTangentOperator(ATC_Coupling * atc,
                               const PhysicsModel * physicsModel,
                               Array2D<bool> & rhsMask,
                               IntegrationDomainType integrationType,
                               FIELDS & rhs,
                               FIELDS & fields,
                               FieldName fieldName,
                               const int dof=0);
      PhysicsModelTangentOperator(ATC_Coupling * atc,
                               const PhysicsModel * physicsModel,
                               Array2D<bool> & rhsMask,
                               IntegrationDomainType integrationType,
                               FieldName fieldName,
                               const int dof=0);
      ~PhysicsModelTangentOperator(){};
      void function(const VECTOR & x, DENS_VEC & r);
      void tangent(const VECTOR & x, DENS_VEC & r, MATRIX & K);
    private:
      ATC_Coupling * atc_;
      const PhysicsModel * physicsModel_;
      Array2D<bool> rhsMask_;
      IntegrationDomainType integrationType_;
      FIELDS & rhs_;
      FIELDS & fields_;
      FieldName fieldName_;
      int dof_;
      SPAR_MAT stiffness_; 
  };

};
#endif
