// ATC_Transfer Headers
#include "ExtrinsicModel.h"
#include "ExtrinsicModelTwoTemperature.h"
#include "ATC_Error.h"
#include "TimeIntegrator.h"
#include "ATC_Transfer.h"
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
  ExtrinsicModelManager::ExtrinsicModelManager(ATC_Transfer * atcTransfer) :
    atcTransfer_(atcTransfer)
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
    FieldName thisField;
    int thisIndex;
    int argIndx = 0;

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
      throw ATC_Error(0,"Could not create extrinsic model");
      return;
    }
    ExtrinsicModel * myModel;
    if      (modelType==TWO_TEMPERATURE) {
      cout << "ATC: creating two_temperature extrinsic model \n";
      myModel = new ExtrinsicModelTwoTemperature
	(this,modelType,matFileName);
    }
    extrinsicModels_.push_back(myModel);

     // add new fields to fields data
     map<FieldName,int> fieldSizes;
     myModel->get_num_fields(fieldSizes);
     atcTransfer_->add_fields(fieldSizes);
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
        if ((*imodel)->get_model_type() == modelType)
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
        if ((*imodel)->get_model_type() == modelType)
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
        if ((*imodel)->get_model_type() == modelType)
          (*imodel)->post_init_integrate();
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
        if ((*imodel)->get_model_type() == modelType)
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
        if ((*imodel)->get_model_type() == modelType)
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
        if ((*imodel)->get_model_type() == modelType)
          (*imodel)->set_sources(fields,sources);
    }
  }

  //--------------------------------------------------------
  //  output
  //--------------------------------------------------------
  void ExtrinsicModelManager::output(double dt,OUTPUT_LIST & outputData)
  {
    vector<ExtrinsicModel *>::iterator imodel;
    for(imodel=extrinsicModels_.begin(); imodel!=extrinsicModels_.end(); imodel++)
      (*imodel)->output(dt,outputData);
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
    modelManager_(modelManager),
    modelType_(modelType),
    atc_(modelManager->get_atc_transfer()),
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
  //  get_num_fields
  //  - sets dict of fields
  //--------------------------------------------------------
  void ExtrinsicModel::get_num_fields(map<FieldName,int> & fieldSizes)
  {
    physicsModel_->get_num_fields(fieldSizes,atc_->fieldMask_); // NOTE clunky
  }

};
