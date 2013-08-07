// ATC Headers
#include "ExtrinsicModel.h"
#include "ExtrinsicModelTwoTemperature.h"
#include "ExtrinsicModelDriftDiffusion.h"
#include "ExtrinsicModelElectrostatic.h"
#include "ATC_Error.h"
#include "TimeIntegrator.h"
#include "ATC_Coupling.h"
#include "LammpsInterface.h"
#include "PrescribedDataManager.h"
#include "PhysicsModel.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelManager
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ExtrinsicModelManager::ExtrinsicModelManager(ATC_Coupling * atc) :
    atc_(atc)
  {
    // do nothing
  }
  
  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ExtrinsicModelManager::~ExtrinsicModelManager()
  {
    vector<ExtrinsicModel *>::iterator imodel;
    for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
      delete *(imodel);
  }

  //--------------------------------------------------------
  //  modify
  //--------------------------------------------------------
  bool ExtrinsicModelManager::modify(int narg, char **arg)
  {
    bool foundMatch = false; 
        
    // loop over models with command
    vector<ExtrinsicModel *>::iterator imodel;
    for(imodel=extrinsicModels_.begin(); 
        imodel!=extrinsicModels_.end(); imodel++) {
      foundMatch = (*imodel)->modify(narg,arg);
      if (foundMatch) break;
    }
    
    return foundMatch;
  }

  //--------------------------------------------------------
  //  create_model
  //--------------------------------------------------------
  void ExtrinsicModelManager::create_model(ExtrinsicModelType modelType,
                                           string matFileName)
  {
    string typeName;
    bool validModel = model_to_string(modelType,typeName);
    if (!validModel) {
      throw ATC_Error("Could not create extrinsic model");
      return;
    }
    ExtrinsicModel * myModel;
    if      (modelType==TWO_TEMPERATURE) {
      stringstream ss;
      ss << "creating two_temperature extrinsic model";
      ATC::LammpsInterface::instance()->print_msg_once(ss.str());
      myModel = new ExtrinsicModelTwoTemperature
        (this,modelType,matFileName);
    }
    else if (modelType==DRIFT_DIFFUSION 
         ||  modelType==DRIFT_DIFFUSION_EQUILIBRIUM
         ||  modelType==DRIFT_DIFFUSION_SCHRODINGER
         || modelType==DRIFT_DIFFUSION_SCHRODINGER_SLICE) 
    {
      stringstream ss;
      ss << "creating drift_diffusion extrinsic model";
      ATC::LammpsInterface::instance()->print_msg_once(ss.str());
      myModel = new ExtrinsicModelDriftDiffusion
        (this,modelType,matFileName);
    }
    else if (modelType==CONVECTIVE_DRIFT_DIFFUSION 
         ||  modelType==CONVECTIVE_DRIFT_DIFFUSION_EQUILIBRIUM
         ||  modelType==CONVECTIVE_DRIFT_DIFFUSION_SCHRODINGER) {
      stringstream ss;
      ss << "creating convective_drift_diffusion extrinsic model";
      ATC::LammpsInterface::instance()->print_msg_once(ss.str());
      myModel = new ExtrinsicModelDriftDiffusionConvection
        (this,modelType,matFileName);
    }
    else if (modelType==ELECTROSTATIC || modelType==ELECTROSTATIC_EQUILIBRIUM) {
      stringstream ss;
      ss << "creating electrostatic extrinsic model";
      ATC::LammpsInterface::instance()->print_msg_once(ss.str());
      myModel = new ExtrinsicModelElectrostaticMomentum
        (this,modelType,matFileName);
    }
    else if (modelType==FEM_EFIELD) {
      stringstream ss;
      ss << "creating fem_efield extrinsic model";
      ATC::LammpsInterface::instance()->print_msg_once(ss.str());
      myModel = new ExtrinsicModelElectrostatic
        (this,modelType,matFileName);
    }
    extrinsicModels_.push_back(myModel);

     // add new fields to fields data
     map<FieldName,int> fieldSizes;
     myModel->num_fields(fieldSizes);
     atc_->add_fields(fieldSizes);
  }

  //--------------------------------------------------------
  //  initialize
  //--------------------------------------------------------
  void ExtrinsicModelManager::construct_transfers()
  {
    vector<ExtrinsicModel *>::iterator imodel;
    for(imodel=extrinsicModels_.begin(); 
        imodel!=extrinsicModels_.end(); imodel++) {
      // initialize models
      (*imodel)->construct_transfers();
    }
  }

  //--------------------------------------------------------
  //  initialize
  //--------------------------------------------------------
  void ExtrinsicModelManager::initialize()
  {
    vector<ExtrinsicModel *>::iterator imodel;
    for(imodel=extrinsicModels_.begin(); 
        imodel!=extrinsicModels_.end(); imodel++) {
      // initialize models
      (*imodel)->initialize();
    }
  }

  //--------------------------------------------------------
  //  get_model : access to a particular type of model
  //--------------------------------------------------------
  const ExtrinsicModel * ExtrinsicModelManager::model(const ExtrinsicModelType type) const {
      vector<ExtrinsicModel *>::const_iterator imodel;
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++) {
        if ((*imodel)->model_type()==type) return *imodel;
      }
      return NULL;
    }


  //--------------------------------------------------------
  //  size_vector
  //--------------------------------------------------------
  int ExtrinsicModelManager::size_vector(int intrinsicSize)
  {
    int extrinsicSize = 0;
    vector<ExtrinsicModel *>::iterator imodel;
    for(imodel=extrinsicModels_.begin(); 
        imodel!=extrinsicModels_.end(); imodel++) {
      // query all models for LAMMPS display
      int currentSize = intrinsicSize + extrinsicSize;
      extrinsicSize += (*imodel)->size_vector(currentSize);
    }
    
    return extrinsicSize;
  }

  //--------------------------------------------------------
  //  compute_scalar
  //--------------------------------------------------------
  double ExtrinsicModelManager::compute_scalar(void)
  {
    double value = 0.;
    vector<ExtrinsicModel *>::iterator imodel;
    for(imodel=extrinsicModels_.begin(); 
        imodel!=extrinsicModels_.end(); imodel++) {
      value += (*imodel)->compute_scalar(); // sum
    }
    return value;
  }

  //--------------------------------------------------------
  //  compute_vector
  //--------------------------------------------------------
  double ExtrinsicModelManager::compute_vector(int n)
  {
    double value = 0.;
    vector<ExtrinsicModel *>::iterator imodel;
    for(imodel=extrinsicModels_.begin(); 
        imodel!=extrinsicModels_.end(); imodel++) {
      // query all models for LAMMPS display
      if ((*imodel)->compute_vector(n,value))
        break;
    }
    return value;
  }

  //--------------------------------------------------------
  //  finish
  //--------------------------------------------------------
  void ExtrinsicModelManager::finish()
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  pre_init_integrate
  //--------------------------------------------------------
  void ExtrinsicModelManager::pre_init_integrate(ExtrinsicModelType modelType)
  {
    vector<ExtrinsicModel *>::iterator imodel;
    if (modelType == NUM_MODELS) {// execute all the models
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        (*imodel)->pre_init_integrate();
    }
    else { // execute only requested type of model
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        if ((*imodel)->model_type() == modelType)
          (*imodel)->pre_init_integrate();
    }
  }

  //--------------------------------------------------------
  //  mid_init_integrate
  //--------------------------------------------------------
  void ExtrinsicModelManager::mid_init_integrate(ExtrinsicModelType modelType)
  {
    vector<ExtrinsicModel *>::iterator imodel;
    if (modelType == NUM_MODELS) {// execute all the models
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        (*imodel)->mid_init_integrate();
    }
    else { // execute only requested type of model
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        if ((*imodel)->model_type() == modelType)
          (*imodel)->mid_init_integrate();
    }
  }

  //--------------------------------------------------------
  //  post_init_integrate
  //--------------------------------------------------------
  void ExtrinsicModelManager::post_init_integrate(ExtrinsicModelType modelType)
  {
    vector<ExtrinsicModel *>::iterator imodel;
    if (modelType == NUM_MODELS) {// execute all the models
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        (*imodel)->post_init_integrate();
    }
    else { // execute only requested type of model
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        if ((*imodel)->model_type() == modelType)
          (*imodel)->post_init_integrate();
    }
  }

  //--------------------------------------------------------
  //  post_force
  //--------------------------------------------------------
  void ExtrinsicModelManager::post_force(ExtrinsicModelType modelType)
  {
    vector<ExtrinsicModel *>::iterator imodel;
    if (modelType == NUM_MODELS) {// execute all the models
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        (*imodel)->post_force();
    }
    else { // execute only requested type of model
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        if ((*imodel)->model_type() == modelType)
          (*imodel)->post_force();
    }
  }

  //--------------------------------------------------------
  //  pre_final_integrate
  //--------------------------------------------------------
  void ExtrinsicModelManager::pre_final_integrate(ExtrinsicModelType modelType)
  {
    vector<ExtrinsicModel *>::iterator imodel;
    if (modelType == NUM_MODELS) {// execute all the models
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        (*imodel)->pre_final_integrate();
    }
    else { // execute only requested type of model
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        if ((*imodel)->model_type() == modelType)
          (*imodel)->pre_final_integrate();
    }
  }

  //--------------------------------------------------------
  //  post_final_integrate
  //--------------------------------------------------------
  void ExtrinsicModelManager::post_final_integrate(ExtrinsicModelType modelType)
  {
    vector<ExtrinsicModel *>::iterator imodel;
    if (modelType == NUM_MODELS) {// execute all the models
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        (*imodel)->post_final_integrate();
    }
    else { // execute only requested type of model
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        if ((*imodel)->model_type() == modelType)
          (*imodel)->post_final_integrate();
    }
  }

  //--------------------------------------------------------
  //  set_sources
  //--------------------------------------------------------
  void ExtrinsicModelManager::set_sources(FIELDS & fields, FIELDS & sources, ExtrinsicModelType modelType)
  {
    vector<ExtrinsicModel *>::iterator imodel;
    if (modelType == NUM_MODELS) {// execute all the models
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        (*imodel)->set_sources(fields,sources);
    }
    else { // execute only requested type of model
      for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
        if ((*imodel)->model_type() == modelType)
          (*imodel)->set_sources(fields,sources);
    }
  }

  //--------------------------------------------------------
  //  output
  //--------------------------------------------------------
  void ExtrinsicModelManager::output(OUTPUT_LIST & outputData)
  {
    vector<ExtrinsicModel *>::iterator imodel;
    for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
      (*imodel)->output(outputData);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModel
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ExtrinsicModel::ExtrinsicModel(ExtrinsicModelManager * modelManager,
                                 ExtrinsicModelType modelType,
                                 string matFileName) :
    atc_(modelManager->atc()),
    modelManager_(modelManager),
    modelType_(modelType),
    physicsModel_(NULL)
  {
    rhsMaskIntrinsic_.reset(NUM_FIELDS,NUM_FLUX);
    rhsMaskIntrinsic_ = false;
  }
  
  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ExtrinsicModel::~ExtrinsicModel()
  {
    if (physicsModel_) delete physicsModel_;
  }

  //--------------------------------------------------------
  // initialize
  //--------------------------------------------------------
  void ExtrinsicModel::initialize(void)
  {
    physicsModel_->initialize();
  }

  //--------------------------------------------------------
  //  get_num_fields
  //  - sets dict of fields
  //--------------------------------------------------------
  void ExtrinsicModel::num_fields(map<FieldName,int> & fieldSizes)
  {
    physicsModel_->num_fields(fieldSizes,atc_->fieldMask_); 
  }

};
