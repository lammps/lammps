#ifndef PRESCRIBED_DATA_MANAGER_H
#define PRESCRIBED_DATA_MANAGER_H

// manager for initial conditions, essential/natural "boundary" conditions 
// and sources 

// to do: 
//        handle no sources, time-independent sources
//        prescribed surface sources

#include <vector>
#include <map>
#include <string>

#include "XT_Function.h"
#include "PhysicsModel.h"
#include "FE_Element.h"
#include "Array.h"
#include "Array2D.h"
#include "FE_Engine.h"

//clase FE_Engine


namespace ATC {
  using std::vector;
  using std::pair;
  using std::map;

  class PrescribedDataManager {

  public:

    /** exclusive conditions: free | fixed field | flux or domain source */
    //enum Bc_Type {FREE=0,FIELD,SOURCE};

    PrescribedDataManager(FE_Engine * feEngine, 
                          const map<FieldName,int> & fieldSize);
    ~PrescribedDataManager();

    /** add/remove a field */
    void add_field(FieldName fieldName, int size);
    void remove_field(FieldName fieldName);

    /** direct access to ics */
    map < FieldName, Array2D < XT_Function * > > *  
      get_ics(void) { return & ics_; }
    const Array2D < XT_Function * > *  
      get_ics(FieldName fieldName) { return & ics_[fieldName]; }
    /** direct access to bcs */
    map < FieldName, Array2D < XT_Function * > > *  
      get_bcs(void) { return & bcs_; }
    const Array2D < XT_Function * > *  
      get_bcs(FieldName fieldName) { return & bcs_[fieldName]; }
    
    /** query initial state */
    bool is_initially_fixed(const int node, 
                            const FieldName thisField,
                            const int thisIndex=0)  const
    { 
      return ((ics_.find(thisField)->second))(node,thisIndex) ? true : false ; 
    }
    /** query state */
    bool is_fixed(const int node, 
                  const FieldName thisField,
                  const int thisIndex=0)  const
    { 
      return ((bcs_.find(thisField)->second))(node,thisIndex) ? true : false ; 
    }
  
    /** set initial field values */
    void fix_initial_field (const string nodesetName, 
                            const FieldName thisField, 
                            const int thisIndex, 
                            const XT_Function * f);
    /** un/set field values at fixed nodes */
    void fix_field (const string nodesetName, 
                    const FieldName thisField, 
                    const int thisIndex, 
                    const XT_Function * f);
    void unfix_field (const string nodesetName, 
                      const FieldName thisField, 
                      const int thisIndex); 
    /** un/set fluxes */
    void fix_flux  (const string facesetName,
                    const FieldName thisField, 
                    const int thisIndex, 
                    const XT_Function * f);
    void unfix_flux(const string facesetName,
                    const FieldName thisField, 
                    const int thisIndex);
    /** un/set sources */
    void fix_source(const string nodesetName,
                    const FieldName thisField, 
                    const int thisIndex, 
                    const XT_Function * f);
    void unfix_source(const string nodesetName,
                      const FieldName thisField, 
                      const int thisIndex);
    /** get initial conditions  */
    void set_initial_conditions(const double time, 
                                FIELDS & fields,
                                FIELDS & dot_fields,
                                FIELDS & ddot_fields,
                                FIELDS & dddot_fields);
    /** get "boundary" conditions on fields */
    void set_fixed_fields(const double time, 
                          FIELDS & fields,
                          FIELDS & dot_fields,
                          FIELDS & ddot_fields,
                          FIELDS & dddot_fields);
    /** get "boundary" conditions on a single field */
    void set_fixed_field(const double time, 
                         const FieldName & fieldName,
                         DENS_MAT & fieldMatrix);
    /** get "boundary" conditions on a single time derivative field */
    void set_fixed_dfield(const double time, 
                          const FieldName & fieldName,
                          DENS_MAT & dfieldMatrix);
    /** get "sources" (flux and sources: divided by leading coef of ODE) */
    void set_sources(const double time, 
                     FIELDS & sources);

    /** debugging status output */
    void print(void);

  private:
    /** number of nodes */
    int nNodes_;

    /** number of elements */
    int nElems_;

    /** names and sizes of fields */
    map<FieldName,int> fieldSizes_;

    /** access to all the FE computations */
    FE_Engine * feEngine_;

    // node numbering & dof numbering : contiguous
    // fieldname & bc_type : types/enums
    /** ics : XT_Function * f = ics_[field](inode,idof) */
    map < FieldName, Array2D < XT_Function * > > ics_;

    /** bcs: essential bcs XT_Function * f = bcs_[field][face](idof) */
    map < FieldName, Array2D < XT_Function * > > bcs_;

    /** sources :  XT_Function * f = faceSources_[field][face](idof) */
    map < FieldName, map < pair <int, int>, Array < XT_Function * > > > 
      faceSources_;
    /** sources :  XT_Function * f = elementSources_[field](ielem,idof) */
    map < FieldName, Array2D < XT_Function * > > elementSources_;
  };

}

#endif
