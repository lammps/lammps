#include "PrescribedDataManager.h"
#include "FE_Engine.h"
#include "ATC_Error.h"

#include <set>

namespace ATC {

//-------------------------------------------------------------------------
//  PrescribedDataManager
//-------------------------------------------------------------------------
  PrescribedDataManager::PrescribedDataManager
  (FE_Engine * feEngine, 
   const map<FieldName,int> & fieldSize) :
    feEngine_(feEngine), fieldSizes_(fieldSize)
  {
    // construct & initialize internal data
    nNodes_ =  feEngine_->get_nNodes();
    nElems_ =  feEngine_->get_nElements();
    map<FieldName,int>::const_iterator field;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      int thisSize = field->second;
      // nodal ics & essential bcs
      ics_[thisField].reset(nNodes_,thisSize);
      bcs_[thisField].reset(nNodes_,thisSize);
      for (int inode = 0; inode < nNodes_ ; ++inode) {
        for (int idof = 0; idof < thisSize ; ++idof) {
          ics_[thisField](inode,idof) = NULL;          
          bcs_[thisField](inode,idof) = NULL;          
        }
      }
      // element based sources
      elementSources_[thisField].reset(nElems_,thisSize);
      for (int ielem = 0; ielem < nElems_ ; ++ielem) {
        for (int idof = 0; idof < thisSize ; ++idof) {
          elementSources_[thisField](ielem,idof) = NULL;          
        }
      }
    }
  }

//-------------------------------------------------------------------------
// ~PrescribedDataManager
//-------------------------------------------------------------------------
  PrescribedDataManager::~PrescribedDataManager()
  {
  }

//-------------------------------------------------------------------------
// add_field
//-------------------------------------------------------------------------
  void PrescribedDataManager::add_field(FieldName fieldName, int size)
  {
    // check to see if field exists
    if (fieldSizes_.find(fieldName) == fieldSizes_.end()) return;

    // construct & initialize internal data
    nNodes_ =  feEngine_->get_nNodes();
    nElems_ =  feEngine_->get_nElements();
  
    // nodal ics & essential bcs
    ics_[fieldName].reset(nNodes_,size);
    bcs_[fieldName].reset(nNodes_,size);
    for (int inode = 0; inode < nNodes_ ; ++inode) {
      for (int idof = 0; idof < size ; ++idof) {
        ics_[fieldName](inode,idof) = NULL;          
        bcs_[fieldName](inode,idof) = NULL;          
      }
    }

    // element based sources
    elementSources_[fieldName].reset(nElems_,size);
    for (int ielem = 0; ielem < nElems_ ; ++ielem) {
      for (int idof = 0; idof < size ; ++idof) {
        elementSources_[fieldName](ielem,idof) = NULL;          
      }
    }
  }

//-------------------------------------------------------------------------
// remove_field
//-------------------------------------------------------------------------
  void PrescribedDataManager::remove_field(FieldName fieldName)
  {
    // check to see if field exists
    if (fieldSizes_.find(fieldName) == fieldSizes_.end())
      return;
    
    // delete field in maps
    fieldSizes_.erase(fieldName);
    ics_.erase(fieldName);
    bcs_.erase(fieldName);
    elementSources_.erase(fieldName);
  }

//-------------------------------------------------------------------------
//  fix_initial_field
//-------------------------------------------------------------------------
  void PrescribedDataManager::fix_initial_field
  (const string nodesetName, 
   const FieldName thisField, 
   const int thisIndex, 
   const XT_Function * f)
  {
    using std::set;
    set<int> nodeSet = (feEngine_->get_feMesh())->get_nodeset(nodesetName);
    set<int>::const_iterator iset;
    for (iset = nodeSet.begin(); iset != nodeSet.end(); iset++) {
      int inode = *iset;
      ics_[thisField](inode,thisIndex) = (XT_Function*) f; // NOTE const cast?
    }
  }
//-------------------------------------------------------------------------
//  fix_field
//-------------------------------------------------------------------------
  void PrescribedDataManager::fix_field
  (const string nodesetName, 
   const FieldName thisField, 
   const int thisIndex, 
   const XT_Function * f)
  {
    using std::set;
    // fix fields 
    set<int> nodeSet = (feEngine_->get_feMesh())->get_nodeset(nodesetName);
    set<int>::const_iterator iset;
    for (iset = nodeSet.begin(); iset != nodeSet.end(); iset++) {
      int inode = *iset;
      bcs_[thisField](inode,thisIndex) = (XT_Function*) f;
    }
  }
//-------------------------------------------------------------------------
//  unfix_field
//-------------------------------------------------------------------------
  void PrescribedDataManager::unfix_field
  (const string nodesetName, 
   const FieldName thisField, 
   const int thisIndex)
  {
    using std::set;
    set<int> nodeSet = (feEngine_->get_feMesh())->get_nodeset(nodesetName);
    set<int>::const_iterator iset;
    for (iset = nodeSet.begin(); iset != nodeSet.end(); iset++) {
      int inode = *iset;
      bcs_[thisField](inode,thisIndex) = NULL;
    }
  }
//-------------------------------------------------------------------------
//  fix_flux
//-------------------------------------------------------------------------
  void PrescribedDataManager::fix_flux
  (const string facesetName, 
   const FieldName thisField, 
   const int thisIndex, 
   const XT_Function * f)
  {
    const set< pair <int,int> > * fset 
      = & ( (feEngine_->get_feMesh())->get_faceset(facesetName));
    set< pair<int,int> >::const_iterator iset;
    for (iset = fset->begin(); iset != fset->end(); iset++) {
      pair<int,int>  face = *iset;
      // allocate, if necessary
      Array < XT_Function * > & dof  = faceSources_[thisField][face];
      if (dof.get_length() == 0) {
        int ndof = (fieldSizes_.find(thisField))->second;
        dof.reset(ndof);
        for(int i = 0; i < ndof; i++)  dof(i) = NULL; 
      }
      dof(thisIndex) = (XT_Function*) f;
    }
  }
//-------------------------------------------------------------------------
//  unfix_flux
//-------------------------------------------------------------------------
  void PrescribedDataManager::unfix_flux
  (const string facesetName, 
   const FieldName thisField, 
   const int thisIndex) 
  {
    const set< pair <int,int> > * fset 
      = & ( (feEngine_->get_feMesh())->get_faceset(facesetName));
    set< pair<int,int> >::const_iterator iset;
    for (iset = fset->begin(); iset != fset->end(); iset++) {
      pair<int,int>  face = *iset;
      Array < XT_Function * > & dof  = faceSources_[thisField][face];
      dof(thisIndex) = NULL;
    }
  }
//-------------------------------------------------------------------------
//  fix_source
//-------------------------------------------------------------------------
  void PrescribedDataManager::fix_source
  (const string elemsetName, 
   const FieldName thisField, 
   const int thisIndex, 
   const XT_Function *f)
  {
    using std::set;
    set<int> elemSet = (feEngine_->get_feMesh())->get_elementset(elemsetName);
    set<int>::const_iterator iset;
    for (iset = elemSet.begin(); iset != elemSet.end(); iset++) {
      int ielem = *iset;
      // fix source
      elementSources_[thisField](ielem,thisIndex) = (XT_Function*) f;
    }
  }
//-------------------------------------------------------------------------
//  unfix_source
//-------------------------------------------------------------------------
  void PrescribedDataManager::unfix_source
  (const string elemsetName, 
   const FieldName thisField, 
   const int thisIndex)
  {
    using std::set;
    set<int> elemSet = (feEngine_->get_feMesh())->get_elementset(elemsetName);
    set<int>::const_iterator iset;
    for (iset = elemSet.begin(); iset != elemSet.end(); iset++) {
      int ielem = *iset;
      elementSources_[thisField](ielem,thisIndex) = NULL;
    }
  }
//-------------------------------------------------------------------------
//  set_initial_conditions
//-------------------------------------------------------------------------
  void PrescribedDataManager::set_initial_conditions(const double t,  
                                                     FIELDS &fields,
                                                     FIELDS &dot_fields,
                                                     FIELDS &ddot_fields,
                                                     FIELDS &dddot_fields)
  {
    map<FieldName,int>::const_iterator field;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      int thisSize = field->second;
      for (int inode = 0; inode < nNodes_ ; ++inode) {
        for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
          XT_Function *f = ics_[thisField](inode,thisIndex);
          if (!f)  f = bcs_[thisField](inode,thisIndex);
          if (f) 
          {
            DENS_VEC coords(3); 
            coords = (feEngine_->get_feMesh())->nodal_coordinates(inode);
            double *x = coords.get_ptr();
            fields      [thisField](inode,thisIndex) = f->f(x,t);
            dot_fields  [thisField](inode,thisIndex) = f->dfdt(x,t); 
            ddot_fields [thisField](inode,thisIndex) = f->ddfdt(x,t); 
            dddot_fields[thisField](inode,thisIndex) = f->dddfdt(x,t); 
          }
          else throw ATC_Error(0,"all initial conditions have not been defined");
        }
      }
    }
  }
//-------------------------------------------------------------------------
//  set_fixed_fields
//-------------------------------------------------------------------------
  void PrescribedDataManager::set_fixed_fields(const double t, 
                                               FIELDS &fields,
                                               FIELDS &dot_fields,
                                               FIELDS &ddot_fields,
                                               FIELDS &dddot_fields)
  {
    map<FieldName,int>::const_iterator field;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      int thisSize = field->second;
      for (int inode = 0; inode < nNodes_ ; ++inode) {
        for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
          XT_Function * f = bcs_[thisField](inode,thisIndex);
          if (f) {
            DENS_VEC coords(3); 
            coords = (feEngine_->get_feMesh())->nodal_coordinates(inode);
            double * x = coords.get_ptr();
            FIELD * myField = & fields[thisField];
            fields      [thisField](inode,thisIndex) = f->f(x,t);
            dot_fields  [thisField](inode,thisIndex) = f->dfdt(x,t); 
            ddot_fields [thisField](inode,thisIndex) = f->ddfdt(x,t); 
            dddot_fields[thisField](inode,thisIndex) = f->dddfdt(x,t); 
          }
        }
      }
    }
  }
//-------------------------------------------------------------------------
//  set_fixed_field
//-------------------------------------------------------------------------
  void PrescribedDataManager::set_fixed_field(
    const double t, 
    const FieldName & fieldName,
    DENS_MAT & fieldMatrix)
  {
    map<FieldName,int>::iterator fieldSizeIter = fieldSizes_.find(fieldName);
    if (fieldSizeIter == fieldSizes_.end()) {
      throw ATC_Error(0, "Unrecognized FieldName in PrescribedDataManager::set_fixed_field()");
    }
    int thisSize = fieldSizeIter->second;
    for (int inode = 0; inode < nNodes_ ; ++inode) {
      for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
        XT_Function * f = bcs_[fieldName](inode,thisIndex);
        if (f) {
          DENS_VEC coords(3); 
          coords = (feEngine_->get_feMesh())->nodal_coordinates(inode);
          fieldMatrix(inode,thisIndex) = f->f(coords.get_ptr(),t);
        }
      }
    }
  }
//-------------------------------------------------------------------------
//  set_fixed_dfield
//-------------------------------------------------------------------------
  void PrescribedDataManager::set_fixed_dfield(
    const double t, 
    const FieldName & fieldName,
    DENS_MAT & dfieldMatrix)
  {
    map<FieldName,int>::iterator fieldSizeIter = fieldSizes_.find(fieldName);
    if (fieldSizeIter == fieldSizes_.end()) {
      throw ATC_Error(0, "Unrecognized FieldName in PrescribedDataManager::set_fixed_dfield()");
    }
    int thisSize = fieldSizeIter->second;
    for (int inode = 0; inode < nNodes_ ; ++inode) {
      for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
        XT_Function * f = bcs_[fieldName](inode,thisIndex);
        if (f) {
          DENS_VEC coords(3); 
          coords = (feEngine_->get_feMesh())->nodal_coordinates(inode);
          dfieldMatrix(inode,thisIndex) = f->dfdt(coords.get_ptr(),t);
        }
      }
    }
  }
//-------------------------------------------------------------------------
//  set_sources
//-------------------------------------------------------------------------
  void PrescribedDataManager::set_sources
  (double t,  
   FIELDS & sources)
  {
    // zero
    Array<bool> fieldMask(NUM_FIELDS);
    fieldMask = false;
    map<FieldName,int>::const_iterator field;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      fieldMask(thisField) = true;
      int thisSize = field->second;
      sources[thisField].reset(nNodes_,thisSize);
    }
    // compute boundary fluxes
    feEngine_->add_fluxes(fieldMask,t,faceSources_,sources);
  
    // compute internal sources
    feEngine_->add_sources(fieldMask,t,elementSources_,sources);

    // mask out nodes with essential bcs
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      int thisSize = field->second;
      for (int inode = 0; inode < nNodes_ ; ++inode) {
        for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
          XT_Function * f = bcs_[thisField](inode,thisIndex);
          if (f) {
            sources[thisField](inode,thisIndex) = 0.0;
          }
        }
      }
    }
  
  }

//-------------------------------------------------------------------------
//  print
//-------------------------------------------------------------------------
  void PrescribedDataManager::print(void) 
  {
    // print and check consistency
    enum dataType {FREE=0,FIELD,SOURCE};
    Array2D < int > bcTypes;
    Array <int> conn;
    map<FieldName,int>::const_iterator field;
    XT_Function * f;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      int thisFieldSize = field->second;
      string fieldName = field_to_string(thisField);
      int thisSize = field->second;
      bcTypes.reset(nNodes_,thisSize);
      for (int inode = 0; inode < nNodes_ ; ++inode) {
        for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
          f = bcs_[thisField](inode,thisIndex);
          if (f) { bcTypes(inode,thisIndex) = FIELD; }
          else   { bcTypes(inode,thisIndex) = FREE; }
        }
      }
      // FIXED has higher precidence than SOURCE
      for (int ielem = 0; ielem < nElems_ ; ++ielem) {
        for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
          f = elementSources_[thisField](ielem,thisIndex);
          if (f) { 
            feEngine_->element_connectivity(ielem,conn);
            for (int i = 0; i < conn.get_length() ; ++i) {
              int inode = conn(i);
              if (bcTypes(inode,thisIndex) != FIELD) 
              { bcTypes(inode,thisIndex) = SOURCE; }
            }
          }
        }
      }
      map < pair<int,int>, Array < XT_Function * > > & fset 
        = faceSources_[thisField];
      map < pair<int,int>, Array < XT_Function * > > ::const_iterator iset;
      for (iset = fset.begin(); iset != fset.end(); iset++) {
        pair<int,int> face = iset->first;
        Array < XT_Function * > fs = iset->second;
        for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
          f = fs(thisIndex);
          if (f) { 
            feEngine_->face_connectivity(face,conn);
            for (int i = 0; i < conn.get_length() ; ++i) {
              int inode = conn(i);
              if (bcTypes(inode,thisIndex) != FIELD) 
              { bcTypes(inode,thisIndex) = SOURCE; }
            }
          }
        }
      }
      for (int inode = 0; inode < nNodes_ ; ++inode) {
        for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
          cout << "node: " << inode << " "  <<  fieldName;
          if (thisFieldSize > 1) { cout << " " << thisIndex; }
          f = ics_[thisField](inode,thisIndex);
          if (f) { cout << " IC"; }
          if      (bcTypes(inode,thisIndex) == FIELD ) {cout << " FIXED"; }
          else if (bcTypes(inode,thisIndex) == SOURCE) {cout << " SOURCE"; }
          cout << "\n";
        }
      }
    }
  }

} // end namespace
