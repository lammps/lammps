#include "PhysicsModel.h"
#include "WeakEquation.h"
#include "WeakEquationDiffusion.h"
#include "WeakEquationChargeDiffusion.h"
#include "WeakEquationMassDiffusion.h"
#include "WeakEquationElectronContinuity.h"
#include "WeakEquationElectronMomentum.h"
#include "WeakEquationElectronTemperature.h"
#include "WeakEquationMomentum.h"
#include "WeakEquationPhononTemperature.h"
#include "WeakEquationPoisson.h"
#include "WeakEquationSchrodinger.h"
#include "ATC_Coupling.h" // for tangent operator
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>

using ATC_Utility::command_line;
using ATC_Utility::to_lower;
using std::fstream;
using std::stringstream;
using std::pair;
using std::string;
using std::map;
using std::set;
using std::vector;

namespace ATC {

//---------------------------------------------------------------------
// PhysicsModel
//---------------------------------------------------------------------
PhysicsModel::PhysicsModel(string fileName)
{
  parse_material_file(fileName);
}

PhysicsModel::~PhysicsModel()
{
  {
    vector< Material* >::iterator iter;
    for (iter = materials_.begin(); iter != materials_.end(); iter++) {
          Material * mat = *iter;
      if (mat) delete mat;
    }
  }
  {
    map<FieldName, WeakEquation* >::iterator iter;
    for (iter = weakEqns_.begin(); iter != weakEqns_.end(); iter++) {
      WeakEquation * weakEq = iter->second;
      if (weakEq) delete weakEq;
    }
  }
}

void PhysicsModel::parse_material_file(string fileName)
{
  vector< Material* >::iterator iter;
  for (iter = materials_.begin(); iter != materials_.end(); iter++) {
    Material * mat = *iter;
    if (mat) delete mat;
  }
  LammpsInterface::UnitsType lammpsUnits = ATC::LammpsInterface::instance()->units_style();
  fstream  fileId(fileName.c_str(), std::ios::in);
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  int index = 0;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0 || line[0] == "#") continue;
    if (line[0] == "material") {
      string tag = line[1];
      Material * mat = new Material(tag,fileId);
      materials_.push_back(mat);
      materialNameToIndexMap_[tag] = index++;
      if (line.size() > 2) {
        string units = line[2];
        stringstream ss;
        ss << "WARNING: material units " << units << " do not match lammps";
        if      (units == "SI") {
          if (lammpsUnits != LammpsInterface::SI) 
            ATC::LammpsInterface::instance()->print_msg_once(ss.str());
        }
        else if (units == "real") {
          if (lammpsUnits != LammpsInterface::REAL ) 
            ATC::LammpsInterface::instance()->print_msg_once(ss.str());
        }
        else if (units == "metal") {
          if (lammpsUnits != LammpsInterface::METAL ) 
            ATC::LammpsInterface::instance()->print_msg_once(ss.str());
        } 
        else {
          throw ATC_Error("unknown units in material file"); 
        }
      }
      else {
        throw ATC_Error("units need to be specified in material file");
      }
    }
  }
  if (int(materials_.size()) == 0) {
    throw ATC_Error("No materials were defined"); }
  stringstream ss;

  ss << int(materials_.size()) << " materials defined from " << fileName;
  ATC::LammpsInterface::instance()->print_msg_once(ss.str());
  fileId.close();
}

void PhysicsModel::initialize(void)
{
  // initialize materials
  vector< Material* >::const_iterator iter;
  for (iter = materials_.begin(); iter != materials_.end(); iter++) {
    Material * mat = *iter;
     mat->initialize();
  }

  // set up null (weakEq, material) registry
  null_.reset(ATC::NUM_FIELDS, materials_.size());
  null_ = false;

  // initialize weak equations
  map<FieldName, WeakEquation * >::const_iterator weak;
  for (weak = weakEqns_.begin(); weak!=weakEqns_.end(); weak++) {
    FieldName fieldName = weak->first;
    WeakEquation * weakEq = weak->second;
    set<string> needs= weakEq->needs_material_functions();    
    vector< Material* >::iterator iter;
    for (iter = materials_.begin(); iter != materials_.end(); iter++) {
      Material * mat = *iter;
      if (! (mat->check_registry(needs)) ) {
        string tag = mat->label();
        int matId = materialNameToIndexMap_[tag];
        null_(fieldName,matId) = true;
        stringstream ss;
        ss << "WARNING: physics model: [" << type_ << "], material: [" << tag
           << "] does not provide all interfaces for <" 
           << field_to_string(fieldName) 
           << "> physics and will be treated as null ";
        ATC::LammpsInterface::instance()->print_msg_once(ss.str());
        // if (noNull_) 
        //throw ATC_Error("material does not provide all interfaces for physics");
      }
    }
  }
}

int PhysicsModel::material_index(const string & name) const
{
  string tag = name;
  to_lower(tag); // this is an artifact of StringManip parsing
  map<string,int>::const_iterator iter;
  iter = materialNameToIndexMap_.find(tag);
  if (iter ==  materialNameToIndexMap_.end()) {
    throw ATC_Error("No material named "+name+" found");
  }
  int index = iter->second;
  return index;
}

bool PhysicsModel::parameter_value(const string& name, double& value,
                     const int imat) const
{
  // search owned parameters
  value = 0.0;
  map<string,double>::const_iterator it = parameterValues_.find(name);
  if (it != parameterValues_.end()) {
    value = it->second;
    return true;
  }
  // interogate material models
  bool found = materials_[imat]->parameter(name,value);
  return found;
}

void PhysicsModel::num_fields(map<FieldName,int> & fieldSizes,
                            Array2D<bool> & rhsMask) const
{
  map<FieldName, WeakEquation * >::const_iterator itr;
  for (itr = weakEqns_.begin(); itr!=weakEqns_.end(); itr++) {
    FieldName field = itr->first;
    WeakEquation * weakEq = itr->second;
    int size = weakEq->field_size();
    fieldSizes[field] = size;
    rhsMask(field,FLUX)   = weakEq->has_B_integrand();
    rhsMask(field,SOURCE) = weakEq->has_N_integrand();
  } 
}

//---------------------------------------------------------------------
// PhysicsModelThermal
//---------------------------------------------------------------------
PhysicsModelThermal::PhysicsModelThermal(string filename):
  PhysicsModel(filename)
{
  weakEqns_[TEMPERATURE] = new WeakEquationPhononTemperature();
}
//---------------------------------------------------------------------
// PhysicsModelElastic
//---------------------------------------------------------------------
PhysicsModelElastic::PhysicsModelElastic(string filename):
  PhysicsModel(filename)
{
  weakEqns_[VELOCITY] = new WeakEquationMomentum();
}
//---------------------------------------------------------------------
// PhysicsModelThermoElastic
//---------------------------------------------------------------------
PhysicsModelThermoElastic::PhysicsModelThermoElastic(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "thermo-elastic";
  weakEqns_[VELOCITY]    = new WeakEquationMomentum();
  weakEqns_[TEMPERATURE] = new WeakEquationPhononTemperature();
}
//---------------------------------------------------------------------
// PhysicsModelShear
//---------------------------------------------------------------------
PhysicsModelShear::PhysicsModelShear(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "shear";
  weakEqns_[VELOCITY] = new WeakEquationMomentumDiffusion();
}
//---------------------------------------------------------------------
// PhysicsModelThermoShear
//---------------------------------------------------------------------
PhysicsModelThermoShear::PhysicsModelThermoShear(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "thermo-shear";
  weakEqns_[VELOCITY]    = new WeakEquationMomentumDiffusion();
  weakEqns_[TEMPERATURE] = new WeakEquationPhononTemperature();
}
//---------------------------------------------------------------------
// PhysicsModelSpecies
//---------------------------------------------------------------------
PhysicsModelSpecies::PhysicsModelSpecies(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "species";

  weakEqns_[MASS_DENSITY]   = new WeakEquationMassDiffusion();
  weakEqns_[SPECIES_CONCENTRATION] = new WeakEquationDiffusion();
}
//---------------------------------------------------------------------
// PhysicsModelSpeciesElectrostatic
//---------------------------------------------------------------------
PhysicsModelSpeciesElectrostatic::PhysicsModelSpeciesElectrostatic(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "species electrostatic";

  weakEqns_[MASS_DENSITY]         = new WeakEquationMassDiffusion();
  weakEqns_[CHARGE_DENSITY] = new WeakEquationChargeDiffusion();
  weakEqns_[SPECIES_CONCENTRATION]= new WeakEquationDiffusion();
  weakEqns_[ELECTRIC_POTENTIAL]   = new WeakEquationPoissonConstantRHS();
}
//---------------------------------------------------------------------
// PhysicsModelTwoTemperature
//---------------------------------------------------------------------
PhysicsModelTwoTemperature::PhysicsModelTwoTemperature(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "two-temperature";
  weakEqns_[TEMPERATURE]          = new WeakEquationPhononTemperatureExchange();
  weakEqns_[ELECTRON_TEMPERATURE] = new WeakEquationElectronTemperature();
}
//---------------------------------------------------------------------
// PhysicsModelDriftDiffusion
//---------------------------------------------------------------------
PhysicsModelDriftDiffusion::PhysicsModelDriftDiffusion(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "drift-diffusion";
  weakEqns_[TEMPERATURE]          = new WeakEquationPhononTemperatureExchange();
  weakEqns_[ELECTRON_TEMPERATURE] = new WeakEquationElectronTemperatureJouleHeating();
  weakEqns_[ELECTRON_DENSITY]     = new WeakEquationElectronContinuity();
  weakEqns_[ELECTRIC_POTENTIAL]   = new WeakEquationPoissonConstantRHS();
}
//---------------------------------------------------------------------
// PhysicsModelDriftDiffusionEquilibrium
//---------------------------------------------------------------------
PhysicsModelDriftDiffusionEquilibrium::PhysicsModelDriftDiffusionEquilibrium(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "equilibrium drift-diffusion";
  weakEqns_[TEMPERATURE]          = new WeakEquationPhononTemperatureExchange();
  weakEqns_[ELECTRON_TEMPERATURE] = new WeakEquationElectronTemperatureJouleHeating();
  weakEqns_[ELECTRON_DENSITY]     = new WeakEquationElectronEquilibrium();
  weakEqns_[ELECTRIC_POTENTIAL]   = new WeakEquationPoisson();
}
//---------------------------------------------------------------------
// PhysicsModelDriftDiffusionSchrodinger
//---------------------------------------------------------------------
PhysicsModelDriftDiffusionSchrodinger::PhysicsModelDriftDiffusionSchrodinger(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "schrodinger drift-diffusion";
  weakEqns_[TEMPERATURE]          = new WeakEquationPhononTemperatureExchange();
  weakEqns_[ELECTRON_TEMPERATURE] = new WeakEquationElectronTemperatureJouleHeating();
  weakEqns_[ELECTRON_DENSITY]     = new WeakEquationElectronEquilibrium();
  weakEqns_[ELECTRIC_POTENTIAL]   = new WeakEquationPoisson();
  weakEqns_[ELECTRON_WAVEFUNCTION]= new WeakEquationSchrodinger();
}
//---------------------------------------------------------------------
// PhysicsModelDriftDiffusionSchrodingerSlice
//---------------------------------------------------------------------
PhysicsModelDriftDiffusionSchrodingerSlice::PhysicsModelDriftDiffusionSchrodingerSlice(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "schrodinger drift-diffusion slice";
  weakEqns_[TEMPERATURE]          = new WeakEquationPhononTemperatureExchange();
  weakEqns_[ELECTRON_TEMPERATURE] = new WeakEquationElectronTemperatureJouleHeating();
  weakEqns_[ELECTRON_DENSITY]     = new WeakEquationElectronEquilibrium();
  weakEqns_[ELECTRIC_POTENTIAL]   = new WeakEquationPoissonConstantRHS();
  weakEqns_[ELECTRON_WAVEFUNCTION]= new WeakEquationSchrodinger();
}
//---------------------------------------------------------------------
// PhysicsModelDriftDiffusionConvection
//---------------------------------------------------------------------
PhysicsModelDriftDiffusionConvection::PhysicsModelDriftDiffusionConvection(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "convection drift-diffusion";
  weakEqns_[TEMPERATURE]             = new WeakEquationPhononTemperatureExchange();
  weakEqns_[ELECTRON_TEMPERATURE]    = new WeakEquationElectronTemperatureConvection();
  weakEqns_[ELECTRON_DENSITY]        = new WeakEquationElectronContinuity();
  weakEqns_[ELECTRON_VELOCITY]       = new WeakEquationElectronMomentumDDM();
  weakEqns_[ELECTRIC_POTENTIAL]      = new WeakEquationPoissonConstantRHS();
}
//---------------------------------------------------------------------
// PhysicsModelDriftDiffusionConvectionEquilibrium
//---------------------------------------------------------------------
PhysicsModelDriftDiffusionConvectionEquilibrium::PhysicsModelDriftDiffusionConvectionEquilibrium(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "equilibrium convection drift-diffusion";
  weakEqns_[TEMPERATURE]             = new WeakEquationPhononTemperatureExchange();
  weakEqns_[ELECTRON_TEMPERATURE]    = new WeakEquationElectronTemperatureConvection();
  weakEqns_[ELECTRON_DENSITY]        = new WeakEquationElectronEquilibrium();
  weakEqns_[ELECTRON_VELOCITY]       = new WeakEquationElectronMomentumDDM();
  weakEqns_[ELECTRIC_POTENTIAL]      = new WeakEquationPoisson();
}
//---------------------------------------------------------------------
// PhysicsModelDriftDiffusionConvectionSchrodinger
//---------------------------------------------------------------------
PhysicsModelDriftDiffusionConvectionSchrodinger::PhysicsModelDriftDiffusionConvectionSchrodinger(string filename):
  PhysicsModel(filename)
{
  PhysicsModel::type_ = "schrodinger convection drift-diffusion";
  weakEqns_[TEMPERATURE]             = new WeakEquationPhononTemperatureExchange();
  weakEqns_[ELECTRON_TEMPERATURE]    = new WeakEquationElectronTemperatureConvection();
  weakEqns_[ELECTRON_DENSITY]        = new WeakEquationElectronEquilibrium();
  weakEqns_[ELECTRON_VELOCITY]       = new WeakEquationElectronMomentumDDM();
  weakEqns_[ELECTRIC_POTENTIAL]      = new WeakEquationPoissonConstantRHS();
  weakEqns_[ELECTRON_WAVEFUNCTION]   = new WeakEquationSchrodinger();
}
//---------------------------------------------------------------------
// PhysicsModelElectrostatic
//---------------------------------------------------------------------
PhysicsModelElectrostatic::PhysicsModelElectrostatic(string filename):
  PhysicsModel(filename) 
{
  PhysicsModel::type_ = "electrostatic";
  weakEqns_[VELOCITY]           = new WeakEquationMomentumElectrostatic();
  weakEqns_[ELECTRON_DENSITY]   = new WeakEquationElectronContinuity();
  weakEqns_[ELECTRIC_POTENTIAL] = new WeakEquationPoissonConstantRHS();
}
//---------------------------------------------------------------------
// PhysicsModelElectrostaticEquilibrium
//---------------------------------------------------------------------
PhysicsModelElectrostaticEquilibrium::PhysicsModelElectrostaticEquilibrium(string filename):
  PhysicsModel(filename) 
{
  PhysicsModel::type_ = "equilibrium electrostatic";
  weakEqns_[VELOCITY]           = new WeakEquationMomentumElectrostatic();
  weakEqns_[ELECTRON_DENSITY]   = new WeakEquationElectronEquilibrium();
  weakEqns_[ELECTRIC_POTENTIAL] = new WeakEquationPoisson();
}

//-------------------------------------------------------------------
// PhysicsModelTangentOperator
//-------------------------------------------------------------------
PhysicsModelTangentOperator::PhysicsModelTangentOperator(ATC_Coupling * atc,
                               const PhysicsModel * physicsModel,
                               Array2D<bool> & rhsMask,
                               IntegrationDomainType integrationType,
                               FIELDS & rhs,
                               FIELDS & fields,
                               FieldName fieldName,
                               const int dof)
         :TangentOperator(),
          atc_(atc),
          physicsModel_(physicsModel),
          rhsMask_(rhsMask),
          integrationType_(integrationType),
          rhs_(rhs),
          fields_(fields),
          fieldName_(fieldName),
          dof_(dof)
 {
 };

PhysicsModelTangentOperator::PhysicsModelTangentOperator(ATC_Coupling * atc,
                               const PhysicsModel * physicsModel,
                               Array2D<bool> & rhsMask,
                               IntegrationDomainType integrationType,
                               FieldName fieldName,
                               const int dof)
         :TangentOperator(),
          atc_(atc),
          physicsModel_(physicsModel),
          rhsMask_(rhsMask),
          integrationType_(integrationType),
          rhs_(atc_->rhs()),
          fields_(atc_->fields()),
          fieldName_(fieldName),
          dof_(dof)
 {
 };

void PhysicsModelTangentOperator::function(const VECTOR & x, DENS_VEC & r)
{
  CLON_VEC f = column(fields_[fieldName_].set_quantity(),dof_); 
  f = x;
  atc_->compute_rhs_vector(rhsMask_, fields_, rhs_, integrationType_, physicsModel_);
  CLON_VEC rhsv = column(rhs_[fieldName_].quantity(),dof_);
  r = rhsv;
}
void PhysicsModelTangentOperator::tangent(const VECTOR & x, DENS_VEC & r, 
                                          MATRIX & K)
{
  function(x,r);
  pair<FieldName,FieldName> row_col(fieldName_,fieldName_);
  const FIELDS & fields = fields_;
  atc_->compute_rhs_tangent(row_col, rhsMask_, fields , stiffness_, integrationType_, physicsModel_);
  K = stiffness_.dense_copy(); 
}

}; // end namespace
