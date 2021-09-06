#ifndef PRESCRIBED_DATA_MANAGER_H
#define PRESCRIBED_DATA_MANAGER_H

#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>

#include "ATC_TypeDefs.h"
#include "Function.h"
#include "PhysicsModel.h"
#include "FE_Element.h"
#include "Array.h"
#include "Array2D.h"
#include "FE_Engine.h"


namespace ATC {

  /**
   *  @class  PrescribedDataManager
   *  @brief  Base class for managing initial conditions, essential/natural "boundary" conditions and sources
   */
  class PrescribedDataManager {

  public:

    /** exclusive conditions: free | fixed field | flux or domain source */
    //enum Bc_Type {FREE=0,FIELD,SOURCE};

    PrescribedDataManager(FE_Engine * feEngine,
                          const std::map<FieldName,int> & fieldSize);
    ~PrescribedDataManager();

    /** add/remove a field */
    void add_field(FieldName fieldName, int size);
    void remove_field(FieldName fieldName);

    /** direct access to ics */
    std::map < FieldName, Array2D < XT_Function * > > *
      ics(void) { return & ics_; }
    const Array2D < XT_Function * > *
      ics(FieldName fieldName) { return & ics_[fieldName]; }
    /** direct access to bcs */

    const std::map < FieldName, BCS > &  bcs(void) const
    {
      return bcValues_;
    }
    /** */
    const BCS & bcs(const FieldName fieldName) const
    {
      return (bcValues_.find(fieldName))->second;
    }
    /** */
    void bcs
      (const FieldName fieldName, const std::set<int> nodes, BCS & bcs,
      bool local = false) const;
    /** */
    std::map < FieldName, Array2D < XT_Function * > > *
      bc_functions(void) { return & bcs_; }
    /** */
    const Array2D < XT_Function * > *
      bc_functions(FieldName fieldName) { return & bcs_[fieldName]; }
    /** */
    ROBIN_SURFACE_SOURCE * robin_functions(void) { return & faceSourcesRobin_; }
    bool has_robin_source(FieldName fieldName) const {
       return ((faceSourcesRobin_.find(fieldName)->second).size() > 0) ;
    }
    /** */
    const std::map<PAIR, Array<UXT_Function*> > *
      robin_functions(FieldName fieldName) { return & faceSourcesRobin_[fieldName]; }
    /** */
    OPEN_SURFACE * open_faces(void) { return & facesOpen_; }
    bool has_open_face(FieldName fieldName) const {
       return ((facesOpen_.find(fieldName)->second).size() > 0) ;
    }
    /** */
    const std::set<PAIR> *
      open_faces(FieldName fieldName) { return & facesOpen_[fieldName]; }

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

    /** */
    std::set<int> fixed_nodes(
                  const FieldName thisField,
                  const int thisIndex=0)  const
    {
      std::set<int> fixed;
      const Array2D < XT_Function *> & bcs = bcs_.find(thisField)->second;
      for (int node = 0; node < bcs.nRows() ; node++) {
        if (bcs(node,thisIndex)) fixed.insert(node);
      }
      return fixed;
    }
    /** */
    void fixed_nodes(
                     const FieldName thisField,
                     std::set<int> & fixed,
                     const int thisIndex=0)  const
    {
      const Array2D < XT_Function *> & bcs = bcs_.find(thisField)->second;
      for (int node = 0; node < bcs.nRows() ; node++) {
        if (bcs(node,thisIndex)) fixed.insert(node);
      }
    }
    /** to determine whether a solution is needed */
    bool all_fixed(
                  const FieldName thisField,
                  const int thisIndex=-1)  const
    {
       if (thisIndex < 0) {
         // static_casts are to iterface with std::vector without compiler warngings
         bool allFixed = (fixed_nodes(thisField,0).size() == static_cast<unsigned>(nNodes_) );
         int ndof = (fieldSizes_.find(thisField)->second);
         for (int j = 1; j < ndof; ++j) {
           allFixed = allFixed && (fixed_nodes(thisField,j).size() == static_cast<unsigned>(nNodes_));
         }
         return allFixed;
       }
       else {
         return (fixed_nodes(thisField,thisIndex).size() == static_cast<unsigned>(nNodes_) );
       }
    }
    /** a function to determine if the tangent is invertible */
    bool none_fixed(
                  const FieldName thisField,
                  const int thisIndex=0)  const
    {
       return (fixed_nodes(thisField,thisIndex).size() == 0 )
           && (faceSourcesRobin_.size() == 0)
           && (facesOpen_.size() == 0);
    }

    /** */
    std::set<int> flux_face_nodes(
                  const FieldName thisField,
                  const int thisIndex=0)  const
    {
      std::set<int> fluxes;
      //list of nodes to insert.
      //1 for nodes to insert, 0 for nodes not to insert.
      int *toInsert = new int[nNodes_];
      for (int i = 0; i < nNodes_; ++i) toInsert[i] = 0;

      const std::map < std::pair <int, int>, Array < XT_Function * > > & sources = faceSources_.find(thisField)->second;
      std::map < std::pair <int, int>, Array < XT_Function * > >::const_iterator fset_iter;
      for (fset_iter = sources.begin(); fset_iter != sources.end(); fset_iter++) {
        int ielem = fset_iter->first.first;
        // if this is not our element, do not do calculations
        if (!feEngine_->fe_mesh()->is_owned_elt(ielem)) continue;
        const Array <XT_Function*> &fs = fset_iter->second;
        if (fs(thisIndex)) {
          Array<int> nodes;
          feEngine_->face_connectivity(fset_iter->first,nodes);
          for (int node = 0; node < nodes.size(); node++) {
            toInsert[nodes(node)] = 1;
          }
        }
      }
      // gather partial results
      LammpsInterface::instance()->logical_or(MPI_IN_PLACE, toInsert, nNodes_);
      // insert selected elements into fluxes
      for (int node = 0; node < nNodes_; ++node) {
        if (toInsert[node]) fluxes.insert(node);
      }
      delete[] toInsert;
      return fluxes;
    }

    /** */
    void flux_face_nodes(
                  const FieldName thisField,
                  std::set<int> fluxes,
                  const int thisIndex=0)  const
    {
      //list of nodes to insert.
      //1 for nodes to insert, 0 for nodes not to insert.
      int *toInsert = new int[nNodes_];
      for (int i = 0; i < nNodes_; ++i) toInsert[i] = 0;

      const std::map < std::pair <int, int>, Array < XT_Function * > > & sources = faceSources_.find(thisField)->second;
      std::map < std::pair <int, int>, Array < XT_Function * > >::const_iterator fset_iter;
      for (fset_iter = sources.begin(); fset_iter != sources.end(); fset_iter++) {
        int ielem = fset_iter->first.first;
        // if this is not our element, do not do calculations
        if (!feEngine_->fe_mesh()->is_owned_elt(ielem)) continue;
        const Array <XT_Function*> &fs = fset_iter->second;
        if (fs(thisIndex)) {
          Array<int> nodes;
          feEngine_->face_connectivity(fset_iter->first,nodes);
          for (int node = 0; node < nodes.size(); node++) {
            toInsert[nodes(node)] = 1;
          }
        }
      }
      // gather partial results
      LammpsInterface::instance()->logical_or(MPI_IN_PLACE, toInsert, nNodes_);
      // insert selected elements into fluxes
      for (int node = 0; node < nNodes_; ++node) {
        if (toInsert[node]) fluxes.insert(node);
      }
      delete[] toInsert;
    }

    /** */
    std::set<int> flux_element_nodes(
                  const FieldName thisField,
                  const int thisIndex=0)  const
    {
      std::set<int> fluxes;
      //list of nodes to insert.
      //1 for nodes to insert, 0 for nodes not to insert.
      int *toInsert = new int[nNodes_];
      for (int i = 0; i < nNodes_; ++i) toInsert[i] = 0;

      const Array2D < XT_Function *> & sources = elementSources_.find(thisField)->second;
      for (int element = 0; element < sources.nRows() ; element++) {
        // if this is not our element, do not do calculations
        if (!feEngine_->fe_mesh()->is_owned_elt(element)) continue;
        if (sources(element,thisIndex)) {
          Array<int> nodes;
          feEngine_->element_connectivity(element,nodes);
          for (int node = 0; node < nodes.size(); node++) {
            toInsert[nodes(node)] = 1;
          }
        }
      }
      // gather partial results
      LammpsInterface::instance()->logical_or(MPI_IN_PLACE, toInsert, nNodes_);
      // insert selected elements into fluxes
      for (int node = 0; node < nNodes_; ++node) {
        if (toInsert[node]) fluxes.insert(node);
      }
      delete[] toInsert;
      return fluxes;
    }

    /** */
    void flux_element_nodes(
                  const FieldName thisField,
                  std::set<int> fluxes,
                  const int thisIndex=0)  const
    {
      //list of nodes to insert.
      //1 for nodes to insert, 0 for nodes not to insert.
      int *toInsert = new int[nNodes_];
      for (int i = 0; i < nNodes_; ++i) toInsert[i] = 0;

      const Array2D < XT_Function *> & sources = elementSources_.find(thisField)->second;
      for (int element = 0; element < sources.nRows() ; element++) {
        // if this is not our element, do not do calculations
        if (!feEngine_->fe_mesh()->is_owned_elt(element)) continue;
        if (sources(element,thisIndex)) {
          Array<int> nodes;
          feEngine_->element_connectivity(element,nodes);
          for (int node = 0; node < nodes.size(); node++) {
            toInsert[nodes(node)] = 1;
          }
        }
      }
      // gather partial results
      LammpsInterface::instance()->logical_or(MPI_IN_PLACE, toInsert, nNodes_);
      // insert selected elements into fluxes
      for (int node = 0; node < nNodes_; ++node) {
        if (toInsert[node]) fluxes.insert(node);
      }
      delete[] toInsert;
    }

    /** */
    bool no_fluxes(
                   const FieldName thisField,
                   const int thisIndex=0)  const
    {
      return ((flux_element_nodes(thisField,thisIndex).size() == 0) &&
              (flux_face_nodes(thisField,thisIndex).size() == 0));
    }

    /** set initial field values */
    void fix_initial_field (const std::string nodesetName,
                            const FieldName thisField,
                            const int thisIndex,
                            const XT_Function * f);
    /** un/set field values at fixed nodesets */
    void fix_field (const std::set<int> nodeset,
                    const FieldName thisField,
                    const int thisIndex,
                    const XT_Function * f);
    /** un/set field values at fixed nodesets */
    void fix_field (const std::string nodesetName,
                    const FieldName thisField,
                    const int thisIndex,
                    const XT_Function * f);
    void unfix_field (const std::string nodesetName,
                      const FieldName thisField,
                      const int thisIndex);
    /** un/set field values at fixed nodes */
    void fix_field (const int nodeId,
                    const FieldName thisField,
                    const int thisIndex,
                    const XT_Function * f);
    void unfix_field (const int nodeId,
                      const FieldName thisField,
                      const int thisIndex);
    /** un/set fluxes */
    void fix_flux  (const std::string facesetName,
                    const FieldName thisField,
                    const int thisIndex,
                    const XT_Function * f);
    void unfix_flux(const std::string facesetName,
                    const FieldName thisField,
                    const int thisIndex);
    void fix_robin (const std::string facesetName,
                    const FieldName thisField,
                    const int thisIndex,
                    const UXT_Function * f);
    void unfix_robin(const std::string facesetName,
                    const FieldName thisField,
                    const int thisIndex);
    void fix_open  (const std::string facesetName,
                    const FieldName thisField);
    void unfix_open(const std::string facesetName,
                    const FieldName thisField);
    /** un/set sources */
    void fix_source(const std::string elemsetName,
                    const FieldName thisField,
                    const int thisIndex,
                    const XT_Function * f);
    void fix_source( const FieldName thisField,
                    const int thisIndex,
                    const std::set<std::pair<int,double> >  & source);
    void unfix_source(const std::string elemsetName,
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
    /** number of unique nodes */
    int nNodes_;

    /** number of elements */
    int nElems_;

    /** names and sizes of fields */
    std::map<FieldName,int> fieldSizes_;

    /** access to all the FE computations */
    FE_Engine * feEngine_;

    // node numbering & dof numbering : contiguous
    // fieldname & bc_type : types/enums
    /** ics : XT_Function * f = ics_[field](inode,idof) */
    std::map < FieldName, Array2D < XT_Function * > > ics_;

    /** bcs: essential bcs XT_Function * f = bcs_[field][face](idof) */
    std::map < FieldName, Array2D < XT_Function * > > bcs_;

    /** sources :  XT_Function * f = faceSources_[field][face](idof) */
    std::map < FieldName, std::map < std::pair <int, int>, Array < XT_Function * > > >
      faceSources_;
    /** sources :  UXT_Function * f = faceSourcesRobin_[field][face](idof) */
    std::map < FieldName, std::map < std::pair <int, int>, Array < UXT_Function * > > >
      faceSourcesRobin_;
    /** sources :  facesOpen_[field][face] */
    std::map < FieldName, std::set < std::pair <int, int> > >
      facesOpen_;
    /** sources :  XT_Function * f = elementSources_[field](ielem,idof) */
    std::map < FieldName, Array2D < XT_Function * > > elementSources_;
    FIELDS nodalSources_;

    /** values of bcs in a compact set */
    std::map < FieldName, BCS > bcValues_;
  };

}

#endif
