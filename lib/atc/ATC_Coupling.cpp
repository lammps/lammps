// ATC headers
#include "ATC_Coupling.h"
#include "FE_Engine.h"
#include "Array.h"
#include "Array2D.h"
#include "ATC_Error.h"
#include "PrescribedDataManager.h"
#include "AtomicRegulator.h"
#include "TimeIntegrator.h"
#include "PhysicsModel.h"
#include "AtomToMoleculeTransfer.h"
#include "MoleculeSet.h"
#include "FieldManager.h"

using std::string;
using std::map;
using std::pair;
using std::set;

namespace ATC {
  //--------------------------------------------------
  ATC_Coupling::ATC_Coupling(string groupName, double ** & perAtomArray, LAMMPS_NS::Fix * thisFix) :
    ATC_Method(groupName, perAtomArray, thisFix),
    consistentInitialization_(false),
    equilibriumStart_(false),
    useFeMdMassMatrix_(false),
    trackCharge_(false),
    temperatureDef_(NONE),
    prescribedDataMgr_(NULL),
    physicsModel_(NULL),
    extrinsicModelManager_(this),
    atomicRegulator_(NULL),
    atomQuadForInternal_(true),
    elementMask_(NULL),
    elementMaskMass_(NULL),
    elementMaskMassMd_(NULL),
    internalToMask_(NULL),
    internalElement_(NULL),
    ghostElement_(NULL),
    nodalGeometryType_(NULL),
    bndyIntType_(NO_QUADRATURE),
    bndyFaceSet_(NULL),
    atomicWeightsMask_(NULL),
    shpFcnMask_(NULL),
    shpFcnDerivsMask_(NULL),
    sourceIntegration_(FULL_DOMAIN)
  {
    // size the field mask
    fieldMask_.reset(NUM_FIELDS,NUM_FLUX); 
    fieldMask_ = false;
    // default: no consistent mass matrices
    useConsistentMassMatrix_.reset(NUM_FIELDS);
    useConsistentMassMatrix_ = false;
    mdMassNormalization_ = true;
    // check to see if lammps has any charges
    if (lammpsInterface_->atom_charge()) trackCharge_ = true;
    // default: perform velocity verlet
    integrateInternalAtoms_ = true;
  }
  //--------------------------------------------------
  ATC_Coupling::~ATC_Coupling()
  {
    interscaleManager_.clear(); 
    if (feEngine_) { delete feEngine_; feEngine_ = NULL; } 
    if (physicsModel_) delete physicsModel_;
    if (atomicRegulator_) delete atomicRegulator_;
    if (prescribedDataMgr_) delete prescribedDataMgr_;
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      delete _tiIt_->second;
    }
  }
  //--------------------------------------------------
  // Interactions with LAMMPS fix commands
  // parse input command and pass on to finite element engine 
  //   or physics specific transfers if necessary
  //   revert to physics-specific transfer if no command matches input
  // first keyword is unique to particular class
  // base class keyword matching must apply to ALL physics
  // order:  derived, base, owned objects
  //--------------------------------------------------
  bool ATC_Coupling::modify(int narg, char **arg)
  {
    FieldName thisField;
    int thisIndex;
    int argIdx=0;

    bool match = false; 
    
    // gateways to other modules e.g. extrinsic, control, mesh
    // pass off to extrinsic
    if (strcmp(arg[argIdx],"extrinsic")==0) {
      argIdx++;
      match = extrinsicModelManager_.modify(narg-argIdx,&arg[argIdx]);
    }
    // catch special case
    if ((strcmp(arg[argIdx],"control")==0)
      &&(strcmp(arg[argIdx+1],"charge")==0)) {
      match = extrinsicModelManager_.modify(narg-argIdx,&arg[argIdx]);
    }
    // parsing handled here
    else {
      /*! \page man_initial fix_modify AtC initial 
        \section syntax
        fix_modify AtC initial <field> <nodeset> <constant | function>
        - <field> = field name valid for type of physics, temperature | electron_temperature
        - <nodeset> = name of set of nodes to apply initial condition
        - <constant | function> = value or name of function followed by its
          parameters
        \section examples
        <TT> fix_modify atc initial temperature groupNAME 10. </TT>
        \section description
        Sets the initial values for the specified field at the specified nodes.
        \section restrictions
        keyword 'all' reserved in nodeset name
        \section default 
        none
      */
      // set initial conditions
      if (strcmp(arg[argIdx],"initial")==0) {
        argIdx++;
        parse_field(arg,argIdx,thisField,thisIndex);
        string nsetName(arg[argIdx++]);
        XT_Function * f = NULL;
        // parse constant
        if (narg == argIdx+1) {
          f = XT_Function_Mgr::instance()->constant_function(atof(arg[argIdx]));
        }
        // parse function
        else {
          f = XT_Function_Mgr::instance()->function(&(arg[argIdx]),narg-argIdx);
        }
        prescribedDataMgr_->fix_initial_field(nsetName,thisField,thisIndex,f);
        match = true;
      }

      /*! \page man_fix_nodes fix_modify AtC fix 
        \section syntax
        fix_modify AtC fix <field> <nodeset> <constant | function>
        - <field> = field name valid for type of physics
        - <nodeset> = name of set of nodes to apply boundary condition
        - <constant | function> = value or name of function followed by its
          parameters
        \section examples
        <TT> fix_modify AtC fix temperature groupNAME 10. </TT> \n
        <TT> fix_modify AtC fix temperature groupNAME 0 0 0 10.0 0 0 1.0 </TT> \n
        \section description
        Creates a constraint on the values of the specified field at specified nodes.
        \section restrictions
        keyword 'all' reserved in nodeset name
        \section related 
        see \ref man_unfix_nodes
        \section default 
        none
      */
      // fix and unfix nodes
      else if (strcmp(arg[argIdx],"fix")==0) {
        argIdx++;
        parse_field(arg,argIdx,thisField,thisIndex);
        string nsetName(arg[argIdx++]);
        XT_Function * f = NULL;
        // fix current value 
        if (narg == argIdx) {
          set<int> nodeSet = (feEngine_->fe_mesh())->nodeset(nsetName);
          set<int>::const_iterator iset;
          const DENS_MAT & field =(fields_.find(thisField)->second).quantity(); 
          for (iset = nodeSet.begin(); iset != nodeSet.end(); iset++) {
            int inode = *iset;
            double v = field(inode,thisIndex);
            f = XT_Function_Mgr::instance()->constant_function(v);
            set<int> one; one.insert(inode);
            prescribedDataMgr_->fix_field(one,thisField,thisIndex,f);
          }
         }
        // parse constant
        else if (narg == argIdx+1) {
          f = XT_Function_Mgr::instance()->constant_function(atof(arg[argIdx]));
          prescribedDataMgr_->fix_field(nsetName,thisField,thisIndex,f);
        }
        // parse function
        else {
          f = XT_Function_Mgr::instance()->function(&(arg[argIdx]),narg-argIdx);
          prescribedDataMgr_->fix_field(nsetName,thisField,thisIndex,f);
        }
        match = true;
      }

      /*! \page man_unfix_nodes fix_modify AtC unfix 
        \section syntax
        fix_modify AtC unfix <field> <nodeset> 
        - <field> = field name valid for type of physics
        - <nodeset> = name of set of nodes 
        \section examples
        <TT> fix_modify AtC unfix temperature groupNAME </TT>
        \section description
        Removes constraint on field values for specified nodes.
        \section restrictions
        keyword 'all' reserved in nodeset name
        \section related 
        see \ref man_fix_nodes
        \section default 
        none
      */
      else if (strcmp(arg[argIdx],"unfix")==0) {
        argIdx++;
        parse_field(arg,argIdx,thisField,thisIndex);
        string nsetName(arg[argIdx++]);
        prescribedDataMgr_->unfix_field(nsetName,thisField,thisIndex);
        match = true;
      }

    /*! \page man_source fix_modify AtC source
      \section syntax
       fix_modify AtC source <field> <element_set> <value | function>
        - <field> = field name valid for type of physics
        - <element_set> = name of set of elements 
      \section examples
      <TT> fix_modify atc source temperature middle temporal_ramp 10. 0. </TT>
      \section description
      Add domain sources to the mesh. The units are consistent with LAMMPS's 
      units for mass, length and time and are defined by the PDE being solved,
      e.g. for thermal transfer the balance equation is for energy and source
      is energy per time.
      \section restrictions 
      keyword 'all' reserved in element_set name
      \section related
      see \ref man_remove_source
      \section default
      none
    */
      else if (strcmp(arg[argIdx],"source")==0) {
        argIdx++;
        parse_field(arg,argIdx,thisField,thisIndex);
        string esetName(arg[argIdx++]);
        XT_Function * f = NULL;
        // parse constant
        if (narg == argIdx+1) {
          f = XT_Function_Mgr::instance()->constant_function(atof(arg[argIdx]));
        }
        // parse function
        else {
          f = XT_Function_Mgr::instance()->function(&(arg[argIdx]),narg-argIdx);
        }
        prescribedDataMgr_->fix_source(esetName,thisField,thisIndex,f);
        fieldMask_(thisField,PRESCRIBED_SOURCE) = true;
        match = true;
      }

    /*! \page man_remove_source fix_modify AtC remove_source
      \section syntax
      fix_modify AtC remove_source <field> <element_set>
        - <field> = field name valid for type of physics
        - <element_set> = name of set of elements 
      \section examples
      <TT> fix_modify atc remove_source temperature groupNAME </TT>
      \section description
      Remove a domain source.
      \section restrictions
      keyword 'all' reserved in element_set name
      \section related
      see \ref man_source
      \section default
    */
      else if (strcmp(arg[argIdx],"remove_source")==0) {
        argIdx++;
        parse_field(arg,argIdx,thisField,thisIndex);
        string esetName(arg[argIdx++]);
        prescribedDataMgr_->unfix_source(esetName,thisField,thisIndex);
        fieldMask_(thisField,PRESCRIBED_SOURCE) = false;
        match = true;
      }
  
    /*! \page man_fix_flux fix_modify AtC fix_flux
      \section syntax
       fix_modify AtC fix_flux <field> <face_set> <value | function>
        - <field> = field name valid for type of physics, temperature | electron_temperature
        - <face_set> = name of set of element faces
      \section examples
       <TT> fix_modify atc fix_flux temperature faceSet 10.0 </TT> \n
       
      \section description
       Command for fixing normal fluxes e.g. heat_flux. 
       This command only prescribes the normal component of the physical flux, e.g. heat (energy) flux.
       The units are in AtC units, i.e. derived from the LAMMPS length, time, and mass scales.
      \section restrictions 
      Only normal fluxes (Neumann data) can be prescribed.
      \section related
      see \ref man_unfix_flux
      \section default
    */
      else if (strcmp(arg[argIdx],"fix_flux")==0) {
        argIdx++;
        parse_field(arg,argIdx,thisField,thisIndex);
        string fsetName(arg[argIdx++]);
        XT_Function * f = NULL;
        // parse constant
        if (narg == argIdx+1) {
          f = XT_Function_Mgr::instance()->constant_function(atof(arg[argIdx]));
        }
        // parse function
        else {
          f = XT_Function_Mgr::instance()->function(&(arg[argIdx]),narg-argIdx);
        }
        prescribedDataMgr_->fix_flux(fsetName,thisField,thisIndex,f);
        fieldMask_(thisField,PRESCRIBED_SOURCE) = true;
        match = true;
      }

    /*! \page man_unfix_flux fix_modify AtC unfix_flux
      \section syntax
      fix_modify AtC fix_flux <field> <face_set> <value | function>
        - <field> = field name valid for type of physics, temperature | electron_temperature
        - <face_set> = name of set of element faces
      \section examples
       <TT> fix_modify atc unfix_flux temperature faceSet  </TT> \n
       
      \section description
       Command for removing prescribed normal fluxes e.g. heat_flux, stress. 
      \section restrictions 
      \section related
      see \ref man_unfix_flux
      \section default
    */
      else if (strcmp(arg[argIdx],"unfix_flux")==0) {
        argIdx++;
        parse_field(arg,argIdx,thisField,thisIndex);
        string fsetName(arg[argIdx++]);
        prescribedDataMgr_->unfix_flux(fsetName,thisField,thisIndex);
        fieldMask_(thisField,PRESCRIBED_SOURCE) = false;
        match = true;
      }

      
    /*! \page man_fe_md_boundary fix_modify AtC fe_md_boundary 
      \section syntax
      fix_modify AtC fe_md_boundary <faceset | interpolate | no_boundary> [args]
      \section examples
       <TT> fix_modify atc fe_md_boundary interpolate </TT> \n
      \section description
      Specifies different methods for computing fluxes between between the MD and FE integration regions.  Faceset defines a faceset separating the MD and FE regions and uses finite element face quadrature to compute the flux.  Interpolate uses a reconstruction scheme to approximate the flux, which is more robust but less accurate if the MD/FE boundary does correspond to a faceset.  No boundary results in no fluxes between the systems being computed.
      \section restrictions 
      If faceset is used, all the AtC non-boundary atoms must lie within and completely fill the domain enclosed by the faceset.
      \section related
      see \man_boundary_faceset for how to specify the faceset name.
      \section default
      Interpolate.
    */
      else if (strcmp(arg[argIdx],"fe_md_boundary")==0) {
        bndyIntType_ = FE_INTERPOLATION;// default
        if(strcmp(arg[argIdx],"faceset")==0) {
          argIdx++;
          bndyIntType_ = FE_QUADRATURE;
          string name(arg[argIdx++]);
          bndyFaceSet_ = & ( (feEngine_->fe_mesh())->faceset(name));
        } 
        else if (strcmp(arg[argIdx],"interpolate")==0) {
          argIdx++;
          bndyIntType_ = FE_INTERPOLATION;
        }
        else if (strcmp(arg[argIdx],"no_boundary")==0) { 
          bndyIntType_ = NO_QUADRATURE; 
        }
        else {
          throw ATC_Error("Bad boundary integration type");
        }
      }



    /*! \page man_boundary_faceset fix_modify AtC boundary_faceset 
      \section syntax
      fix_modify AtC boundary_faceset <is | add> [args]
      \section examples
      fix_modify AtC boundary_faceset is obndy
      \section description
      This command species the faceset name when using a faceset to compute the MD/FE boundary fluxes.  The faceset must already exist.
      \section restrictions
      This is only valid when fe_md_boundary is set to faceset.
      \section related
      \man_fe_md_boundary
      \section default
    */
      else if (strcmp(arg[argIdx],"boundary_faceset")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"is")==0) { // replace existing faceset
          argIdx++;
          boundaryFaceNames_.clear();
          string name(arg[argIdx++]);
          boundaryFaceNames_.insert(name);
          match = true;
        }
        else if (strcmp(arg[argIdx],"add")==0) { // add this faceset to list
          argIdx++;
          string name(arg[argIdx]);
          boundaryFaceNames_.insert(name);
          match = true;
        }
      }

      /*! \page man_internal_quadrature fix_modify AtC internal_quadrature 
        \section syntax
        fix_modify atc internal_quadrature <on | off> [region]
        \section examples
        <TT> fix_modify atc internal_quadrature off </TT>
        \section description
        Command to use or not use atomic quadrature on internal elements
        fully filled with atoms. By turning the internal quadrature off
        these elements do not contribute to the governing PDE and the fields
        at the internal nodes follow the weighted averages of the atomic data.
        \section optional
        Optional region tag specifies which finite element nodes will be treated
        as being within the MD region.  This option is only valid with
        internal_quadrature off.
        \section restrictions
        \section related 
        \section default 
        on
      */
      else if (strcmp(arg[argIdx],"internal_quadrature")==0) {
        if (initialized_) {
          throw ATC_Error("Cannot change internal_quadrature method after first run");
        }
        argIdx++;
        if (strcmp(arg[argIdx],"on")==0) {
          argIdx++;
          atomQuadForInternal_ = true;
          match = true;
        }
        else if (strcmp(arg[argIdx],"off")==0) {
          argIdx++;
          if (argIdx == narg) {
            atomQuadForInternal_ = false;
            regionID_ = -1;
            match = true;
          }
          else {
            for (regionID_ = 0; regionID_ < lammpsInterface_->nregion(); regionID_++) 
              if (strcmp(arg[argIdx],lammpsInterface_->region_name(regionID_)) == 0) break;
            if (regionID_ < lammpsInterface_->nregion()) {
              atomQuadForInternal_ = false;
              match = true;
            }
            else {
              throw ATC_Error("Region " + string(arg[argIdx]) + " does not exist");
            }
          }
        }
        if (match) {
          needReset_ = true;
        }
      }

      else if (strcmp(arg[argIdx],"fix_robin")==0) {
        argIdx++;
        parse_field(arg,argIdx,thisField,thisIndex);
        string fsetName(arg[argIdx++]);
        UXT_Function * f = NULL;
        // parse linear
        if (narg == argIdx+2) {
          f = UXT_Function_Mgr::instance()->linear_function(atof(arg[argIdx]),atof(arg[argIdx+1]));
        }
        // parse function
        else {
          throw ATC_Error("unimplemented function");
        }
        prescribedDataMgr_->fix_robin(fsetName,thisField,thisIndex,f);
        fieldMask_(thisField,ROBIN_SOURCE) = true;
        match = true;
      }
      else if (strcmp(arg[argIdx],"unfix_robin")==0) {
        argIdx++;
        parse_field(arg,argIdx,thisField,thisIndex);
        string fsetName(arg[argIdx++]);
        prescribedDataMgr_->unfix_robin(fsetName,thisField,thisIndex);
        fieldMask_(thisField,ROBIN_SOURCE) = false;
        match = true;
      }

      /*! \page man_atomic_charge fix_modify AtC atomic_charge
      \section syntax
      fix_modify AtC <include | omit> atomic_charge
        - <include | omit> = switch to activiate/deactiviate inclusion of intrinsic atomic charge in ATC
      \section examples
       <TT> fix_modify atc compute include atomic_charge </TT>
      \section description
       Determines whether AtC tracks the total charge as a finite element field
      \section restrictions
      Required for:  electrostatics
      \section related
      \section default
      if the atom charge is defined, default is on, otherwise default is off
    */
      else if (strcmp(arg[argIdx],"include")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"atomic_charge")==0) {
          trackCharge_ = true;
          match = true;
          needReset_ = true;
        }
      }
      else if (strcmp(arg[argIdx],"omit")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"atomic_charge")==0) {
          trackCharge_ = false;
          match = true;
          needReset_ = true;
        }
      }

      /*! \page man_source_integration fix_modify AtC source_integration
      \section syntax
      fix_modify AtC source_integration < fe | atom>
      \section examples
       <TT> fix_modify atc source_integration atom </TT>
      \section description
      \section restrictions
      \section related
      \section default
      Default is fe
    */
      else if (strcmp(arg[argIdx],"source_integration")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"fe")==0) {
          sourceIntegration_ = FULL_DOMAIN;
        }
        else {
          sourceIntegration_ = FULL_DOMAIN_ATOMIC_QUADRATURE_SOURCE;
        }
        match = true;
      }

      /*! \page man_consistent_fe_initialization fix_modify AtC consistent_fe_initialization
      \section syntax
       fix_modify AtC consistent_fe_initialization <on | off>
        - <on|off> = switch to activiate/deactiviate the intial setting of FE intrinsic field to match the projected MD field
      \section examples
       <TT> fix_modify atc consistent_fe_initialization on </TT>
      \section description
      Determines whether AtC initializes FE intrinsic fields (e.g., temperature) to match the projected MD values.  This is particularly useful for fully overlapping simulations.
      \section restrictions
      Can be used with:  thermal, two_temperature.
      Cannot be used with time filtering on. Does not include boundary nodes.
      \section related
      \section default
      Default is off
    */
      else if (strcmp(arg[argIdx],"consistent_fe_initialization")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"on")==0) {
          if (timeFilterManager_.filter_dynamics())
            throw ATC_Error("Consistent FE initialization cannot be used with time filtering");
          consistentInitialization_ = true;
          match = true;
        }
        else if (strcmp(arg[argIdx],"off")==0) {
          consistentInitialization_ = false;
          match = true;
        }
      }

      // switch for equilibrium filtering start
      /*! \page man_equilibrium_start fix_modify AtC equilibrium_start
        \section syntax
        fix_modify AtC equilibrium_start <on|off>
        
        \section examples
        <TT> fix_modify atc equilibrium_start on </TT> \n
        
        \section description
        Starts filtered calculations assuming they start in equilibrium, i.e. perfect finite element force balance.
        
        \section restrictions
        only needed before filtering is begun
        
        \section related
        see \ref man_time_filter
        
        \section default
        on
      */
      else if (strcmp(arg[argIdx],"equilibrium_start")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"on")==0) {
          equilibriumStart_ = true;
          match = true;
        }
        else if (strcmp(arg[argIdx],"off")==0) {
          equilibriumStart_ = false;
          match = true;
        }
      }

      /*! \page man_mass_matrix fix_modify AtC mass_matrix
        \section syntax 
        fix_modify AtC mass_matrix <fe | md_fe>
        - <fe | md_fe> = activiate/deactiviate using the FE mass matrix in the MD region
        \section examples
        <TT> fix_modify atc mass_matrix fe </TT>
        \section description
        Determines whether AtC uses the FE mass matrix based on Gaussian quadrature or based on atomic quadrature in the MD region.  This is useful for fully overlapping simulations to improve efficiency.
        \section restrictions
        Should not be used unless the FE region is contained within the MD region, otherwise the method will be unstable and inaccurate
        \section related
        \section default
        Default is off
      */

      else if (strcmp(arg[argIdx],"mass_matrix")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"fe")==0) {
          useFeMdMassMatrix_ = true;
          match = true;
        }
        else {
          useFeMdMassMatrix_ = false;
          match = true;
        }
        if (match) {
          needReset_ = true;
        }
      }

    /*! \page man_material fix_modify AtC material
      \section syntax
      fix_modify AtC material [elementset_name] [material_id]  \n
      \section examples
      <TT> fix_modify AtC material gap_region 2</TT>
      \section description
      Sets the material model in elementset_name to be of type material_id.
      \section restrictions 
      The element set must already be created and the material must be specified in the material file given the the atc fix on construction
      \section related
      \section default
      All elements default to the first material in the material file.
    */
    else if (strcmp(arg[argIdx],"material")==0) {
      argIdx++;
      string elemsetName(arg[argIdx++]);
      int matId = physicsModel_->material_index(arg[argIdx++]);
      using std::set;
      set<int> elemSet = (feEngine_->fe_mesh())->elementset(elemsetName);
      if(elementToMaterialMap_.size() == 0) {
        throw ATC_Error("need mesh before material command");
      }
      // set elementToMaterialMap
      set<int>::const_iterator iset;
      for (iset = elemSet.begin(); iset != elemSet.end(); iset++) {
        int ielem = *iset;
        
        // and the tag a string
        elementToMaterialMap_(ielem) = matId;
      }
      match = true;
      needReset_ = true;
    }

    } // end else 
    // no match, call base class parser
    if (!match) {
      match = ATC_Method::modify(narg, arg);
    }
    return match; // return to FixATC
  }

  //--------------------------------------------------
   /** PDE type */
   WeakEquation::PDE_Type ATC_Coupling::pde_type(const FieldName fieldName) const
   {
     const WeakEquation * weakEq = physicsModel_->weak_equation(fieldName);
     if (weakEq == NULL) return WeakEquation::PROJECTION_PDE; 
     return weakEq->type();
   }
  //--------------------------------------------------
   /** is dynamic PDE */
   bool ATC_Coupling::is_dynamic(const FieldName fieldName) const
   {
     const WeakEquation * weakEq = physicsModel_->weak_equation(fieldName);
     if (weakEq == NULL) return false; 
     return (physicsModel_->weak_equation(fieldName)->type() == WeakEquation::DYNAMIC_PDE);
   }


  //--------------------------------------------------
  /** allow FE_Engine to construct data manager after mesh is constructed */
  void ATC_Coupling::construct_prescribed_data_manager (void) {
    prescribedDataMgr_ = new PrescribedDataManager(feEngine_,fieldSizes_);
  }

  //--------------------------------------------------
  // pack_fields
  //   bundle all allocated field matrices into a list
  //   for output needs
  //--------------------------------------------------
  void ATC_Coupling::pack_fields(RESTART_LIST & data)
  {
    ATC_Method::pack_fields(data);
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->pack_fields(data);
    }
  }

  //--------------------------------------------------------------
  // create_physics_model
  // - method to create physics model 
  //--------------------------------------------------------------
  void ATC_Coupling::create_physics_model(const PhysicsType & physicsType,
                                          string matFileName)
  {
    if (physicsModel_) {
      throw ATC_Error("Attempted to create PhysicsModel multiple times in ATC_Coupling");
    }
    // Create PhysicsModel based on physicsType
    switch (physicsType) {
    case NO_PHYSICS :
      break;
    case THERMAL :
      physicsModel_ = new PhysicsModelThermal(matFileName);
      break;
    case ELASTIC :
      physicsModel_ = new PhysicsModelElastic(matFileName);
      break;
    case SHEAR:
      physicsModel_ = new PhysicsModelShear(matFileName);
      break;
    case SPECIES :
      physicsModel_ = new PhysicsModelSpecies(matFileName);
      break;
    case THERMO_ELASTIC :
      physicsModel_ = new PhysicsModelThermoElastic(matFileName);
      break;
    default:
      throw ATC_Error("Unknown physics type in ATC_Coupling::create_physics_model");
    }
  }

  //--------------------------------------------------------
  void ATC_Coupling::construct_methods()
  {
    ATC_Method::construct_methods();

    // construct needed time filters for mass matrices
    if (timeFilterManager_.need_reset()) {
      init_filter();
      map<FieldName,int>::const_iterator field;
      for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
        FieldName thisField = field->first;
        // fill in mass matrix time filters if needed
        if (!massMatTimeFilters_[thisField])
          massMatTimeFilters_[thisField] = timeFilterManager_.construct(TimeFilterManager::INSTANTANEOUS);
      }
    }
  }
  //-------------------------------------------------------------------
  void ATC_Coupling::init_filter()
  {
    if (timeFilterManager_.need_reset()) {
      map<FieldName,int>::const_iterator field;
      for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
        FieldName thisField = field->first;
        int thisSize = field->second;
        (nodalAtomicFieldsRoc_[thisField].set_quantity()).reset(nNodes_,thisSize);
      }
    }
  }
  //--------------------------------------------------------
  void ATC_Coupling::set_fixed_nodes()
  {
    // set fields
    prescribedDataMgr_->set_fixed_fields(time(),
      fields_,dot_fields_,ddot_fields_,dddot_fields_);


    // set related data
    map<FieldName,int>::const_iterator field;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      int thisSize = field->second;
      DENS_MAT & rhs(rhs_[thisField].set_quantity());
      for (int inode = 0; inode < nNodes_ ; ++inode) {
        for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
          if (prescribedDataMgr_->is_fixed(inode,thisField,thisIndex)) {
            rhs(inode,thisIndex) = 0.;
          }
        }
      }
    }
  }
  //--------------------------------------------------------
  void ATC_Coupling::set_initial_conditions()
  {
    // set fields 
    prescribedDataMgr_->set_initial_conditions(time(),
      fields_,dot_fields_,ddot_fields_,dddot_fields_);

    // set (all) related data
    map<FieldName,int>::const_iterator field;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      int thisSize = field->second;
      DENS_MAT & rhs(rhs_[thisField].set_quantity());
      for (int inode = 0; inode < nNodes_ ; ++inode) {
        for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
          rhs(inode,thisIndex) = 0.;
        }
      }
    }
  }
  //--------------------------------------------------------
  void ATC_Coupling::set_sources()
  {
    // set fields
    prescribedDataMgr_->set_sources(time(),sources_);

  }
  //-----------------------------------------------------------------
  // this is w_a source_a
  void ATC_Coupling::compute_sources_at_atoms(const RHS_MASK & rhsMask,
                                              const FIELDS & fields,
                                              const PhysicsModel * physicsModel,
                                              FIELD_MATS & atomicSources)
  {
    if (shpFcnMask_) {
      feEngine_->compute_source(rhsMask,
                                fields,
                                physicsModel,
                                atomMaterialGroupsMask_,
                                atomicWeightsMask_->quantity(),
                                shpFcnMask_->quantity(),
                                shpFcnDerivsMask_->quantity(),
                                atomicSources);
    }
    else {
      for (FIELDS::const_iterator field = fields.begin(); 
           field != fields.end(); field++) {
        FieldName thisFieldName = field->first;
        FIELDS::const_iterator fieldItr = fields.find(thisFieldName);
        const DENS_MAT & field = (fieldItr->second).quantity();
        atomicSources[thisFieldName].reset(field.nRows(),field.nCols());
      }
    }
  }
  //-----------------------------------------------------------------
  
  void ATC_Coupling::compute_atomic_sources(const RHS_MASK & fieldMask,
                                            const FIELDS & fields,
                                            FIELDS & atomicSources)
  {

    for (FIELDS::const_iterator field = fields.begin();
         field != fields.end(); field++) {
      FieldName thisFieldName = field->first;
      if (is_intrinsic(thisFieldName)) {
        atomicSources[thisFieldName] = 0.;
        if (fieldMask(thisFieldName,FLUX)) {
          atomicSources[thisFieldName] = boundaryFlux_[thisFieldName];
        }
        if (fieldMask(thisFieldName,PRESCRIBED_SOURCE)) {
          atomicSources[thisFieldName] -= fluxMask_*(sources_[thisFieldName].quantity());
        } 


        // add in sources from extrinsic models
        if (fieldMask(thisFieldName,EXTRINSIC_SOURCE))
          atomicSources[thisFieldName] -= fluxMask_*(extrinsicSources_[thisFieldName].quantity());

      }
    }
  }
  //-----------------------------------------------------------------
  void ATC_Coupling::masked_atom_domain_rhs_tangent(
    const pair<FieldName,FieldName> row_col,
    const RHS_MASK & rhsMask,
    const FIELDS & fields, 
    SPAR_MAT & stiffness,
    const PhysicsModel * physicsModel)
  {
    if (shpFcnMask_) {
      feEngine_->compute_tangent_matrix(rhsMask, row_col,
                                        fields, physicsModel, atomMaterialGroupsMask_,
                                        atomicWeightsMask_->quantity(), shpFcnMask_->quantity(),
                                        shpFcnDerivsMask_->quantity(),stiffness);
    }
    else {
      stiffness.reset(nNodes_,nNodes_);
    }
  }
  //-----------------------------------------------------------------
  void ATC_Coupling::compute_rhs_tangent(
    const pair<FieldName,FieldName> row_col,
    const RHS_MASK & rhsMask,
    const FIELDS & fields, 
    SPAR_MAT & stiffness,
    const IntegrationDomainType integrationType,
    const PhysicsModel * physicsModel)
  {

    if (integrationType  == FULL_DOMAIN_ATOMIC_QUADRATURE_SOURCE) {
      RHS_MASK rhsMaskFE = rhsMask;
      RHS_MASK rhsMaskMD = rhsMask; rhsMaskMD = false;
      for (FIELDS::const_iterator field = fields.begin(); 
           field != fields.end(); field++) {
        FieldName thisFieldName = field->first;
        if ( rhsMaskFE(thisFieldName,SOURCE) ) {
          rhsMaskFE(thisFieldName,SOURCE) = false;
          rhsMaskMD(thisFieldName,SOURCE) = true;
        }
      }
      feEngine_->compute_tangent_matrix(rhsMaskFE, row_col,
        fields , physicsModel, elementToMaterialMap_, stiffness);
      masked_atom_domain_rhs_tangent(row_col,
                                     rhsMaskMD,
                                     fields,
                                     stiffnessAtomDomain_,
                                     physicsModel);
      stiffness += stiffnessAtomDomain_; 

    }
    else {
      feEngine_->compute_tangent_matrix(rhsMask, row_col,
        fields , physicsModel, elementToMaterialMap_, stiffness);
    }
    ROBIN_SURFACE_SOURCE & robinFcn = *(prescribedDataMgr_->robin_functions()); 
    feEngine_->add_robin_tangent(rhsMask, fields, time(), robinFcn, stiffness);
  }
  //-----------------------------------------------------------------
  void ATC_Coupling::compute_rhs_vector(const RHS_MASK & rhsMask,
                                        const FIELDS & fields, 
                                              FIELDS & rhs,
                                        const IntegrationDomainType domain,
                                        const PhysicsModel * physicsModel)
  {
    if (!physicsModel) physicsModel = physicsModel_;
    

    // compute FE contributions
    
    evaluate_rhs_integral(rhsMask,fields,rhs,domain,physicsModel);

    for (int n = 0; n < rhsMask.nRows(); n++) {
      FieldName thisFieldName = FieldName(n);
      if (rhsMask(thisFieldName,PRESCRIBED_SOURCE)) {
        if (is_intrinsic(thisFieldName)) {
          rhs[thisFieldName] += fluxMaskComplement_*(sources_[thisFieldName].quantity());
        } 
        else {
          rhs[thisFieldName] +=                     sources_[thisFieldName].quantity();
        }
      }

      // add in sources from extrinsic models
      if (rhsMask(thisFieldName,EXTRINSIC_SOURCE)) {
        if (is_intrinsic(thisFieldName)) {
          rhs[thisFieldName] += fluxMaskComplement_*(extrinsicSources_[thisFieldName].quantity());
        }
        else {
          rhs[thisFieldName] +=                     extrinsicSources_[thisFieldName].quantity();
        }
      }
      
    }
    ROBIN_SURFACE_SOURCE & robinFcn = *(prescribedDataMgr_->robin_functions()); 
    feEngine_->add_robin_fluxes(rhsMask, fields, time(), robinFcn, rhs);
  }
  //-----------------------------------------------------------------
  void ATC_Coupling::masked_atom_domain_rhs_integral(
    const Array2D<bool> & rhsMask,
    const FIELDS & fields, FIELDS & rhs,
    const PhysicsModel * physicsModel)
  {
    if (shpFcnMask_) {
      feEngine_->compute_rhs_vector(rhsMask,
                                    fields,
                                    physicsModel,
                                    atomMaterialGroupsMask_,
                                    atomicWeightsMask_->quantity(),
                                    shpFcnMask_->quantity(),
                                    shpFcnDerivsMask_->quantity(),
                                    rhs);
    }
    else {
      for (FIELDS::const_iterator field = fields.begin(); 
           field != fields.end(); field++) {
        FieldName thisFieldName = field->first;
        FIELDS::const_iterator fieldItr = fields.find(thisFieldName);
        const DENS_MAT & field = (fieldItr->second).quantity();
        (rhs[thisFieldName].set_quantity()).reset(field.nRows(),field.nCols());
      }
    }
  }
  //-----------------------------------------------------------------
  void ATC_Coupling::evaluate_rhs_integral(
    const Array2D<bool> & rhsMask,
    const FIELDS & fields, FIELDS & rhs,
    const IntegrationDomainType integrationType,
    const PhysicsModel * physicsModel)
  {
    
    if (!physicsModel) physicsModel = physicsModel_;

    
    if      (integrationType == FE_DOMAIN ) {
      feEngine_->compute_rhs_vector(rhsMask, 
                                    fields, 
                                    physicsModel,
                                    elementToMaterialMap_,
                                    rhs, 
                                    &(elementMask_->quantity()));
      masked_atom_domain_rhs_integral(rhsMask,
                                      fields,
                                      rhsAtomDomain_,
                                      physicsModel);
      for (FIELDS::const_iterator field = fields.begin(); 
           field != fields.end(); field++) {
        FieldName thisFieldName = field->first;
        rhs[thisFieldName] -= rhsAtomDomain_[thisFieldName].quantity();
      }
    }
    else if (integrationType == ATOM_DOMAIN) {
      
      masked_atom_domain_rhs_integral(rhsMask,
                                      fields,
                                      rhs,
                                      physicsModel);
    }
    else if (integrationType  == FULL_DOMAIN_ATOMIC_QUADRATURE_SOURCE) {
      RHS_MASK rhsMaskFE = rhsMask;
      RHS_MASK rhsMaskMD = rhsMask; rhsMaskMD = false;
      for (FIELDS::const_iterator field = fields.begin(); 
           field != fields.end(); field++) {
        FieldName thisFieldName = field->first;
        if ( rhsMaskFE(thisFieldName,SOURCE) ) {
          rhsMaskFE(thisFieldName,SOURCE) = false;
          rhsMaskMD(thisFieldName,SOURCE) = true;
        }
      }
      feEngine_->compute_rhs_vector(rhsMaskFE,
                                    fields,
                                    physicsModel,
                                    elementToMaterialMap_,
                                    rhs);
      masked_atom_domain_rhs_integral(rhsMaskMD,
                                      fields,
                                      rhsAtomDomain_,
                                      physicsModel);
      for (FIELDS::const_iterator field = fields.begin(); 
           field != fields.end(); field++) {
        FieldName thisFieldName = field->first;

        if ( ((rhs[thisFieldName].quantity()).size() > 0)
         && ((rhsAtomDomain_[thisFieldName].quantity()).size() > 0) )
          rhs[thisFieldName] += rhsAtomDomain_[thisFieldName].quantity();
      }
    }
    else { // domain == FULL_DOMAIN
      feEngine_->compute_rhs_vector(rhsMask,
                                    fields,
                                    physicsModel,
                                    elementToMaterialMap_,
                                    rhs);
    }
  }
  //--------------------------------------------------
  void ATC_Coupling::initialize()
  { 
    // initialize physics model
    if (physicsModel_) physicsModel_->initialize();

    ATC_Method::initialize();
    
    // initialized_ is set to true by derived class initialize()
    // STEP 6 - data initialization continued:  set initial conditions
    if (!initialized_) {
      // Apply integration masking and new ICs
      // initialize schedule derivatives
      try {
        set_initial_conditions();
      }
      catch (ATC::ATC_Error& atcError) {
        if (!useRestart_)
          throw;
      }
    }
   
    // initialize and fix computational geometry, this can be changed in the future for Eulerian calculations that fill and empty elements which is why it is outside a !initialized_ guard
    internalElement_->unfix_quantity();
    if (ghostElement_) ghostElement_->unfix_quantity();
    internalElement_->quantity();
    if (ghostElement_) ghostElement_->quantity();
    nodalGeometryType_->quantity();
    internalElement_->fix_quantity();
    if (ghostElement_) ghostElement_->fix_quantity();
    reset_flux_mask();

    // setup grouping of atoms by material
    reset_atom_materials();

    // reset time filters if needed
    if (timeFilterManager_.need_reset()) {
      if ((!initialized_) || (atomToElementMapType_ == EULERIAN)) {
        map<FieldName,int>::const_iterator field;
        for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
          FieldName thisField = field->first;
          if (is_intrinsic(thisField) && is_dynamic(thisField)) { 
            compute_mass_matrix(thisField);
            if (!useConsistentMassMatrix_(thisField) && !useFeMdMassMatrix_) {
              massMatsMd_[thisField] = massMatsMdInstantaneous_[thisField].quantity();
              massMatsAq_[thisField] = massMatsAqInstantaneous_[thisField].quantity();
              update_mass_matrix(thisField);
            }
          }
        }
      }
    }

    
    // prepare computes for first timestep
    lammpsInterface_->computes_addstep(lammpsInterface_->ntimestep()+1);
  }
  //-------------------------------------------------------------------
  void ATC_Coupling::construct_time_integration_data()
  {
    if (!initialized_) {

      map<FieldName,int>::const_iterator field;
      for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
        FieldName thisField = field->first;
        int thisSize = field->second;
        
        // Allocate fields, initialize to default values, set up initial schedule
        
        fields_[thisField].reset(nNodes_,thisSize);
        dot_fields_[thisField].reset(nNodes_,thisSize);
        ddot_fields_[thisField].reset(nNodes_,thisSize);
        dddot_fields_[thisField].reset(nNodes_,thisSize);
        
        // Allocate restricted fields
        if (is_intrinsic(thisField)) {
          nodalAtomicFields_[thisField].reset(nNodes_,thisSize);
          nodalAtomicFieldsRoc_[thisField].reset(nNodes_,thisSize);
        }

        // Dimension finite element rhs matrix
        rhs_[thisField].reset(nNodes_,thisSize);
        rhsAtomDomain_[thisField].reset(nNodes_,thisSize);
        
        sources_[thisField].reset(nNodes_,thisSize);
        extrinsicSources_[thisField].reset(nNodes_,thisSize);
        boundaryFlux_[thisField].reset(nNodes_,thisSize);
        
        if (is_intrinsic(thisField) && is_dynamic(thisField)) {
          massMats_[thisField].reset(nNodes_,nNodes_); // PARALLELIZE
          massMatsFE_[thisField].reset(nNodes_,nNodes_);
          massMatsAq_[thisField].reset(nNodes_,nNodes_);
          massMatsMd_[thisField].reset(nNodes_,nNodes_);
          massMatsMdInstantaneous_[thisField].reset(nNodes_,nNodes_);
          massMatsAqInstantaneous_[thisField].reset(nNodes_,nNodes_);
          massMatsInv_[thisField].reset(nNodes_,nNodes_);  // PARALLELIZE
          massMatsMdInv_[thisField].reset(nNodes_,nNodes_); // PARALLELIZE
        }
        else {
          // no MD mass matrices needed, regular matrices computed in extrinsic model
          if (useConsistentMassMatrix_(thisField)) {
            // compute FE mass matrix in full domain
            
            consistentMassMats_[thisField].reset(nNodes_,nNodes_); // PARALLELIZE
            consistentMassMatsInv_[thisField].reset(nNodes_,nNodes_); // PARALLELIZE
          }
          else {
            massMats_[thisField].reset(nNodes_,nNodes_); // PARALLELIZE
            massMatsInv_[thisField].reset(nNodes_,nNodes_); // PARALLELIZE
          }
        }
      }
    }
  }
  //--------------------------------------------------------
  //  create_full_element_mask
  //    constructs element mask which only masks out 
  //    null elements
  //--------------------------------------------------------
  MatrixDependencyManager<DenseMatrix, bool> * ATC_Coupling::create_full_element_mask()
  {
    MatrixDependencyManager<DenseMatrix, bool> * elementMaskMan = new MatrixDependencyManager<DenseMatrix, bool>(feEngine_->num_elements(),1);
    DenseMatrix<bool> & elementMask(elementMaskMan->set_quantity());
    elementMask = true;
      
    const set<int> & nullElements = feEngine_->null_elements();
    set<int>::const_iterator iset;
    for (iset = nullElements.begin(); iset != nullElements.end(); iset++) {
      int ielem = *iset;
      elementMask(ielem,0) = false;
    }

    return elementMaskMan;
  }
  //--------------------------------------------------------
  //  create_element_set_mask
  //    constructs element mask based on an element set,
  //    uses ints for MPI communication later
  //--------------------------------------------------------
  MatrixDependencyManager<DenseMatrix, int> * ATC_Coupling::create_element_set_mask(const string & elementSetName)
  {
    MatrixDependencyManager<DenseMatrix, int> * elementMaskMan = new MatrixDependencyManager<DenseMatrix, int>(feEngine_->num_elements(),1);
    DenseMatrix<int> & elementMask(elementMaskMan->set_quantity());
    elementMask = false;

    const set<int> & elementSet((feEngine_->fe_mesh())->elementset(elementSetName));
    set<int>::const_iterator iset;
    for (iset = elementSet.begin(); iset != elementSet.end(); ++iset) {
      int ielem = *iset;
      elementMask(ielem,0) = true;
    }
      
    const set<int> & nullElements = feEngine_->null_elements();
    for (iset = nullElements.begin(); iset != nullElements.end(); iset++) {
      int ielem = *iset;
      elementMask(ielem,0) = false;
    }

    return elementMaskMan;
  }
  //--------------------------------------------------------
  //  set_computational_geometry
  //    constructs needed transfer operators which define
  //    hybrid atom/FE computational geometry
  //--------------------------------------------------------
  void ATC_Coupling::set_computational_geometry()
  {
    ATC_Method::set_computational_geometry();

    // does element contain internal atoms
    if (internalElementSet_.size()) {
      // set up elements and maps based on prescribed element sets
      internalElement_ = create_element_set_mask(internalElementSet_);
    }
    else {
      internalElement_ = new AtomTypeElement(this,atomElement_);
    }
    interscaleManager_.add_dense_matrix_int(internalElement_,
                                            "ElementHasInternal");

    if (groupbitGhost_) {
      atomGhostElement_ = new AtomToElementMap(this,
                                               atomGhostCoarseGrainingPositions_,
                                               GHOST);
      interscaleManager_.add_per_atom_int_quantity(atomGhostElement_,
                                                   "AtomGhostElement");
      
      // does element contain ghost atoms
      ghostElement_ = new AtomTypeElement(this,atomGhostElement_);
      interscaleManager_.add_dense_matrix_int(ghostElement_,
                                              "ElementHasGhost");
    }
    
    // element masking for approximate right-hand side FE atomic quadrature
    if (atomQuadForInternal_) {
      elementMask_ = create_full_element_mask();
    }
    else {
      if (internalElementSet_.size()) {
        // when geometry is based on elements, there are no mixed elements
        elementMask_ = new MatrixDependencyManager<DenseMatrix, bool>;
        (elementMask_->set_quantity()).reset(feEngine_->num_elements(),1,false);
      }
      else {
        elementMask_ = new ElementMask(this);
      }
      internalToMask_ = new AtomToElementset(this,elementMask_);
      interscaleManager_.add_per_atom_int_quantity(internalToMask_,
                                                   "InternalToMaskMap");
    }
    interscaleManager_.add_dense_matrix_bool(elementMask_,
                                             "ElementMask");

    if (useFeMdMassMatrix_) {
      if (atomQuadForInternal_) {
        elementMaskMass_ = elementMask_;
      }
      else {
        elementMaskMass_ = create_full_element_mask();
        interscaleManager_.add_dense_matrix_bool(elementMaskMass_,
                                                 "NonNullElementMask");
      }

      elementMaskMassMd_ = new AtomElementMask(this);
      interscaleManager_.add_dense_matrix_bool(elementMaskMassMd_,
                                               "InternalElementMask");
    }

    // assign element and node types for computational geometry
    if (internalElementSet_.size()) {
      nodalGeometryType_ = new NodalGeometryTypeElementSet(this);
    }
    else {
      nodalGeometryType_ = new NodalGeometryType(this);
    }
    interscaleManager_.add_dense_matrix_int(nodalGeometryType_,
                                            "NodalGeometryType");
  }
  //--------------------------------------------------------
  //  construct_interpolant
  //    constructs: interpolatn, accumulant, weights, and spatial derivatives
  //--------------------------------------------------------
  void ATC_Coupling::construct_interpolant()
  {
    // finite element shape functions for interpolants
    PerAtomShapeFunction * atomShapeFunctions = new PerAtomShapeFunction(this);
    interscaleManager_.add_per_atom_sparse_matrix(atomShapeFunctions,"Interpolant");
    shpFcn_ = atomShapeFunctions;

    // use shape functions for accumulants if no kernel function is provided
    if (!kernelFunction_) {
      accumulant_ = shpFcn_;
    }
    else {
      if (kernelOnTheFly_) throw ATC_Error("ATC_Coupling::construct_transfers - on the fly kernel evaluations not supported");
      PerAtomKernelFunction * atomKernelFunctions = new PerAtomKernelFunction(this);
      interscaleManager_.add_per_atom_sparse_matrix(atomKernelFunctions,
                                                    "Accumulant");
      accumulant_ = atomKernelFunctions;
      accumulantWeights_ = new AccumulantWeights(accumulant_);
      mdMassNormalization_ = false;
    }
    
    this->create_atom_volume();

    // masked atom weights
    if (atomQuadForInternal_) {
      atomicWeightsMask_ = atomVolume_;
    }
    else {
      atomicWeightsMask_ = new MappedDiagonalMatrix(this,
                                                    atomVolume_,
                                                    internalToMask_);
      interscaleManager_.add_diagonal_matrix(atomicWeightsMask_,
                                             "AtomWeightsMask");
    }
    // nodal volumes for mass matrix, relies on atomVolumes constructed in base class construct_transfers
    nodalAtomicVolume_ = new AdmtfShapeFunctionRestriction(this,atomVolume_,shpFcn_);
    interscaleManager_.add_dense_matrix(nodalAtomicVolume_,"NodalAtomicVolume");

    // shape function derivatives, masked shape function and derivatives if needed for FE quadrature in atomic domain
    if (atomQuadForInternal_) {
      shpFcnDerivs_ = new PerAtomShapeFunctionGradient(this);
      interscaleManager_.add_vector_sparse_matrix(shpFcnDerivs_,
                                                  "InterpolantGradient");

      shpFcnMask_ = shpFcn_;
      shpFcnDerivsMask_ = shpFcnDerivs_;
    }
    else {
      bool hasMaskedElt = false;
      const DenseMatrix<bool> & elementMask(elementMask_->quantity());
      for (int i = 0; i < elementMask.size(); ++i) {
        if (elementMask(i,0)) {
          hasMaskedElt = true;
          break;
        }
      }
      if (hasMaskedElt) {
        shpFcnDerivs_ = new PerAtomShapeFunctionGradient(this);
        interscaleManager_.add_vector_sparse_matrix(shpFcnDerivs_,
                                                    "InterpolantGradient");

        shpFcnMask_ = new RowMappedSparseMatrix(this,
                                                shpFcn_,
                                                internalToMask_);
        interscaleManager_.add_sparse_matrix(shpFcnMask_,
                                             "ShapeFunctionMask");
        shpFcnDerivsMask_ = new RowMappedSparseMatrixVector(this,
                                                            shpFcnDerivs_,
                                                            internalToMask_);
        interscaleManager_.add_vector_sparse_matrix(shpFcnDerivsMask_,"ShapeFunctionGradientMask");
      }
    }
  }
  //--------------------------------------------------------
  //  construct_molecule_transfers
  //--------------------------------------------------------
  void ATC_Coupling::construct_molecule_transfers()
  {
    
    map<string,pair<MolSize,int> >::const_iterator molecule;
    PerAtomQuantity<double> * atomProcGhostCoarseGrainingPositions = interscaleManager_.per_atom_quantity("AtomicProcGhostCoarseGrainingPositions");
    FundamentalAtomQuantity * mass = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                                                  PROC_GHOST);
    for (molecule = moleculeIds_.begin(); molecule != moleculeIds_.end(); molecule++) {
      const string moleculeName = molecule->first;
      int groupbit = (molecule->second).second;
      SmallMoleculeSet * smallMoleculeSet = new SmallMoleculeSet(this,groupbit);
      smallMoleculeSet->initialize();
      interscaleManager_.add_small_molecule_set(smallMoleculeSet,moleculeName);
      SmallMoleculeCentroid * moleculeCentroid = 
        new SmallMoleculeCentroid(this,mass,smallMoleculeSet,atomProcGhostCoarseGrainingPositions);
      interscaleManager_.add_dense_matrix(moleculeCentroid,"MoleculeCentroid"+moleculeName);

      // shape function at molecular coordinates
      PointToElementMap * elementMapMol = 
        new PointToElementMap(this,moleculeCentroid);
      interscaleManager_.add_dense_matrix_int(elementMapMol,
                                              "ElementMap"+moleculeName);
      InterpolantSmallMolecule * shpFcnMol = new InterpolantSmallMolecule(this,
        elementMapMol, moleculeCentroid, smallMoleculeSet);
      interscaleManager_.add_sparse_matrix(shpFcnMol,
                                           "ShapeFunction"+moleculeName);
    }
  }
  //--------------------------------------------------------
  //  construct_transfers
  //    constructs needed transfer operators
  //--------------------------------------------------------
  void ATC_Coupling::construct_transfers()
  {
    ATC_Method::construct_transfers();
    


    extrinsicModelManager_.construct_transfers();
  }
  //--------------------------------------------------
  void ATC_Coupling::delete_mass_mat_time_filter(FieldName thisField)
  {
  }
  //--------------------------------------------------
  void ATC_Coupling::set_mass_mat_time_filter(FieldName thisField,TimeFilterManager::FilterIntegrationType filterIntegrationType)
  {
    massMatTimeFilters_[thisField] = timeFilterManager_.construct(filterIntegrationType);
  }
  //--------------------------------------------------------------
  /** method to trigger construction of mesh data after mesh construction */
  //--------------------------------------------------------------
  void ATC_Coupling::initialize_mesh_data(void)
  {
    int nelts = feEngine_->fe_mesh()->num_elements();
    elementToMaterialMap_.reset(nelts);
    elementToMaterialMap_ = 0;

    construct_prescribed_data_manager();
    meshDataInitialized_ = true;
  }
  //--------------------------------------------------------
  
  void ATC_Coupling::reset_flux_mask(void)
  {
    int i;
    // this is exact only for uniform meshes and certain types of atomic weights
    // \int_{\Omega_MD} N_I dV = \sum_\alpha N_I\alpha V_\alpha
    fluxMask_.reset((invNodeVolumes_.quantity()) 
              * (nodalAtomicVolume_->quantity()));

    DIAG_MAT id(fluxMask_.nRows(),fluxMask_.nCols()); 
    id = 1.0;
    fluxMaskComplement_ = id + -1.0*fluxMask_;

    // set flux masks for nodes we can tell by geometry
    const INT_ARRAY & nodeType(nodalGeometryType_->quantity());
    for (i = 0; i < nNodes_; ++i) {
      if (nodeType(i,0)==MD_ONLY) {
        fluxMask_(i,i) = 1.;
        fluxMaskComplement_(i,i) = 0.;
      }
      else if (nodeType(i,0)==FE_ONLY) {
        fluxMask_(i,i) = 0.;
        fluxMaskComplement_(i,i) = 1.;
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Coupling::compute_mass_matrix(FieldName thisField, PhysicsModel * physicsModel)
  {

    if (!physicsModel) physicsModel = physicsModel_;
    if (useConsistentMassMatrix_(thisField)) {
      // compute FE mass matrix in full domain
      
      Array<FieldName> massMask(1);
      massMask(0) = thisField; 
      
      feEngine_->compute_mass_matrix(massMask,fields_,physicsModel,
                                     elementToMaterialMap_,consistentMassMats_,
                                     &(elementMask_->quantity()));
      // brute force computation of inverse
      consistentMassMatsInv_[thisField] = inv((consistentMassMats_[thisField].quantity()).dense_copy());
    }
    else { // lumped mass matrix
      // compute FE mass matrix in full domain
      Array<FieldName> massMask(1);
      massMask(0) = thisField;

      if (useFeMdMassMatrix_) {
        feEngine_->compute_lumped_mass_matrix(massMask,fields_,physicsModel,
                                              elementToMaterialMap_,massMats_,
                                              &(elementMaskMass_->quantity()));
        const DIAG_MAT & myMassMat(massMats_[thisField].quantity());
        DIAG_MAT & myMassMatInv(massMatsInv_[thisField].set_quantity());
        DIAG_MAT & myMassMatMdInv(massMatsMdInv_[thisField].set_quantity());

        feEngine_->compute_lumped_mass_matrix(massMask,fields_,physicsModel,
                                              elementToMaterialMap_,massMatsMd_,
                                              &(elementMaskMassMd_->quantity()));
        const DIAG_MAT & myMassMatMd(massMatsMd_[thisField].quantity());
        // compute inverse mass matrices since we're using lumped masses
        for (int iNode = 0; iNode < nNodes_; iNode++) {
          
          if (fabs(myMassMat(iNode,iNode))>0)
            myMassMatInv(iNode,iNode) = 1./myMassMat(iNode,iNode);
          else
            myMassMatInv(iNode,iNode) = 0.;

          if (fabs(myMassMatMd(iNode,iNode))>0)
            myMassMatMdInv(iNode,iNode) = 1./myMassMatMd(iNode,iNode);
          else
            myMassMatMdInv(iNode,iNode) = 0.;
        }
      }
      else {
        feEngine_->compute_lumped_mass_matrix(massMask,fields_,physicsModel,
                                              elementToMaterialMap_,massMatsFE_,
                                              &(elementMask_->quantity()));
        // fully remove contributions from internal nodes
        
        DIAG_MAT & myMassMatFE(massMatsFE_[thisField].set_quantity());
        if (!atomQuadForInternal_) {
          const INT_ARRAY & nodeType(nodalGeometryType_->quantity());
          for (int iNode = 0; iNode < nNodes_; iNode++)
            if (nodeType(iNode,0)==MD_ONLY) {
              myMassMatFE(iNode,iNode) = 0.;
            }
        }
        
        // atomic quadrature for FE mass matrix in atomic domain
        if (shpFcnMask_) {
          feEngine_->compute_lumped_mass_matrix(massMask,fields_,physicsModel,atomMaterialGroupsMask_,
                                                atomicWeightsMask_->quantity(),shpFcnMask_->quantity(),
                                                massMatsAqInstantaneous_);
        }
        else {
          (massMatsAqInstantaneous_[thisField].set_quantity()).reset(nNodes_,nNodes_);
        }
        
        // set up mass MD matrices
        compute_md_mass_matrix(thisField,massMatsMdInstantaneous_[thisField].set_quantity());
      }
    }
  }
  //--------------------------------------------------------
  void ATC_Coupling::update_mass_matrix(FieldName thisField)
  {
    DIAG_MAT & myMassMat(massMats_[thisField].set_quantity());
    DIAG_MAT & myMassMatInv(massMatsInv_[thisField].set_quantity());
    DIAG_MAT & myMassMatMDInv(massMatsMdInv_[thisField].set_quantity());
    const DIAG_MAT & myMassMatMD(massMatsMd_[thisField].quantity());

    myMassMat = massMatsFE_[thisField].quantity();
    // remove contributions from overlap by approximate quadrature
    myMassMat -= massMatsAq_[thisField].quantity();
    // add contributions from atomic region
    myMassMat += myMassMatMD;

    // compute inverse mass matrices since we're using lumped masses
    for (int iNode = 0; iNode < nNodes_; iNode++) {
      
      if (fabs(myMassMatMD(iNode,iNode))>0) { 
        myMassMatMDInv(iNode,iNode) = 1./myMassMatMD(iNode,iNode);
      }
      else
        myMassMatMDInv(iNode,iNode) = 0.;
      
      if (fabs(myMassMat(iNode,iNode))>0) {
        myMassMatInv(iNode,iNode) = 1./myMassMat(iNode,iNode);
      }
      else
        myMassMatInv(iNode,iNode) = 0.;
    }
  }
  //--------------------------------------------------------
  void ATC_Coupling::reset_atom_materials()
  {
    int nMaterials = physicsModel_->nMaterials();
    atomMaterialGroups_.reset(nMaterials);
    atomMaterialGroupsMask_.reset(nMaterials);

    for (int i = 0; i < nMaterials; i++) {
      atomMaterialGroups_(i).clear();
      atomMaterialGroupsMask_(i).clear();
    }

    const INT_ARRAY & atomToElementMap(atomElement_->quantity());
    for (int i = 0; i < nLocal_; i++) {
      atomMaterialGroups_(elementToMaterialMap_(atomToElementMap(i,0))).insert(i);
    }
    if (atomQuadForInternal_) {
      for (int i = 0; i < nLocal_; i++) {
        atomMaterialGroupsMask_(elementToMaterialMap_(atomToElementMap(i,0))).insert(i);
      }
    }
    else {
      const INT_ARRAY & map(internalToMask_->quantity());
      for (int i = 0; i < nLocal_; i++) {
        int idx = map(i,0);
        if (idx > -1) {
          atomMaterialGroupsMask_(elementToMaterialMap_(atomToElementMap(i,0))).insert(idx);
        }
      }
    }
  }

  //--------------------------------------------------------
  //  pre_init_integrate
  //    time integration before the lammps atomic
  //    integration of the Verlet step 1
  //--------------------------------------------------------
  void ATC_Coupling::pre_init_integrate()
  {
    ATC_Method::pre_init_integrate();
    double dt = lammpsInterface_->dt();

    // Perform any initialization, no actual integration
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->pre_initial_integrate1(dt);
    }

    // Apply controllers to atom velocities, if needed
    atomicRegulator_->apply_pre_predictor(dt,lammpsInterface_->ntimestep());

    // predict nodal fields and time derivatives
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->pre_initial_integrate2(dt);
    }
    extrinsicModelManager_.pre_init_integrate();
  }

  //--------------------------------------------------------
  //  mid_init_integrate
  //    time integration between the velocity update and
  //    the position lammps update of Verlet step 1
  //--------------------------------------------------------
  void ATC_Coupling::mid_init_integrate()
  {
    double dt = lammpsInterface_->dt();

    // Compute nodal velocity at n+1/2, if needed
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->mid_initial_integrate1(dt);
    }

    atomicRegulator_->apply_mid_predictor(dt,lammpsInterface_->ntimestep());

    extrinsicModelManager_.mid_init_integrate();
  }

  ///--------------------------------------------------------
  //  post_init_integrate
  //    time integration after the lammps atomic updates of
  //    Verlet step 1
  //--------------------------------------------------------
  void ATC_Coupling::post_init_integrate()
  {
    double dt = lammpsInterface_->dt();
  
    // Compute nodal velocity at n+1
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_initial_integrate1(dt);
    }

    // Update kinetostat quantities if displacement is being regulated
    atomicRegulator_->apply_post_predictor(dt,lammpsInterface_->ntimestep());

    // Update extrisic model
    extrinsicModelManager_.post_init_integrate();

    // fixed values, non-group bcs handled through FE
    set_fixed_nodes();
      
    update_time(0.5);

    // ghost update, if needed
    ATC_Method::post_init_integrate();

    // Apply time filtering to mass matrices, if needed
    if (timeFilterManager_.filter_dynamics() && !useFeMdMassMatrix_) {
      map<FieldName,int>::const_iterator field;
      for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
        FieldName thisField = field->first;
        if (!useConsistentMassMatrix_(thisField) && is_intrinsic(thisField)) {
          massMatTimeFilters_[thisField]->apply_pre_step1(massMatsAq_[thisField].set_quantity(),
                                                          massMatsAqInstantaneous_[thisField].quantity(),dt);
          massMatTimeFilters_[thisField]->apply_pre_step1(massMatsMd_[thisField].set_quantity(),
                                                          massMatsMdInstantaneous_[thisField].quantity(),dt);
        }
      }
    }
  }

  
  //--------------------------------------------------------
  void ATC_Coupling::pre_neighbor()
  {
    ATC_Method::pre_neighbor();
    reset_atom_materials();
  }

  //--------------------------------------------------------
  void ATC_Coupling::pre_exchange()
  {
    ATC_Method::pre_exchange();
  }

  //--------------------------------------------------------
  void ATC_Coupling::post_force()
  {
    ATC_Method::post_force();

    if ( (atomToElementMapType_ == EULERIAN) && (step() % atomToElementMapFrequency_ == 0) ) {
      reset_atom_materials();

      if (!useFeMdMassMatrix_) {
        map<FieldName,int>::const_iterator field;
        for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
          FieldName thisField = field->first;
          if (is_intrinsic(thisField) && is_dynamic(thisField)) {
            compute_mass_matrix(thisField);  
          }
        }
      }
    }

    if (atomToElementMapType_ == EULERIAN && !useFeMdMassMatrix_) {
      if (timeFilterManager_.filter_dynamics() || (step() % atomToElementMapFrequency_ == 0)) {
        double dt = lammpsInterface_->dt();
        map<FieldName,int>::const_iterator field;
        for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
          FieldName thisField = field->first;
          if (is_intrinsic(thisField) && is_dynamic(thisField)) {
            massMatTimeFilters_[thisField]->apply_post_step1(massMatsAq_[thisField].set_quantity(),
                                                             massMatsAqInstantaneous_[thisField].quantity(),dt);
            massMatTimeFilters_[thisField]->apply_post_step1(massMatsMd_[thisField].set_quantity(),
                                                             massMatsMdInstantaneous_[thisField].quantity(),dt);
            update_mass_matrix(thisField); 
          }
        }
      }
    }

    // apply extrinsic model
    extrinsicModelManager_.post_force();
  }

  //=================================================================
  //
  //=================================================================
  void ATC_Coupling::finish()
  {
    ATC_Method::finish();
    // Time integrator
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->finish();
    }
  }



  //=================================================================
  //
  //=================================================================
  void ATC_Coupling::compute_boundary_flux(const Array2D<bool> & rhsMask,
                                           const FIELDS & fields, 
                                           FIELDS & rhs,
                                           const Array< set <int> > atomMaterialGroups,
                                           const VectorDependencyManager<SPAR_MAT * > * shpFcnDerivs,
                                           const SPAR_MAN * shpFcn,
                                           const DIAG_MAN * atomicWeights,
                                           
                                           const MatrixDependencyManager<DenseMatrix, bool> * elementMask,
                                           const RegulatedNodes * nodeSet)
  {
    if (bndyIntType_ == FE_QUADRATURE) {
      feEngine_->compute_boundary_flux(rhsMask,
                                       fields,
                                       physicsModel_,
                                       elementToMaterialMap_,
                                       (* bndyFaceSet_),
                                       rhs);
    }
    else if (bndyIntType_ == FE_INTERPOLATION) {
      if (elementMask) {
        feEngine_->compute_boundary_flux(rhsMask,
                                         fields,
                                         physicsModel_,
                                         elementToMaterialMap_,
                                         atomMaterialGroups,
                                         atomicWeights->quantity(),
                                         shpFcn->quantity(),
                                         shpFcnDerivs->quantity(),
                                         fluxMask_,
                                         rhs,
                                         &elementMask->quantity(),
                                         &nodeSet->quantity());
      }
      else {
        feEngine_->compute_boundary_flux(rhsMask,
                                         fields,
                                         physicsModel_,
                                         elementToMaterialMap_,
                                         atomMaterialGroups_,
                                         atomVolume_->quantity(),
                                         shpFcn_->quantity(),
                                         shpFcnDerivs_->quantity(),
                                         fluxMask_,
                                         rhs);
      }
    }
    else if (bndyIntType_ == NO_QUADRATURE) {
      FIELDS::const_iterator field;
      for (field = fields.begin(); field != fields.end(); field++) {
        FieldName thisFieldName = field->first;

if (thisFieldName >= rhsMask.nRows()) break;
        if (rhsMask(thisFieldName,FLUX)) {
          int nrows = (field->second).nRows();
          int ncols = (field->second).nCols();
          rhs[thisFieldName].reset(nrows,ncols);
        }
      }
    }
  }

  //-----------------------------------------------------------------
  void ATC_Coupling::compute_flux(const Array2D<bool> & rhsMask,
                                  const FIELDS & fields, 
                                  GRAD_FIELD_MATS & flux,
                                  const PhysicsModel * physicsModel) 
  {
    if (! physicsModel) { physicsModel = physicsModel_; }
    feEngine_->compute_flux(rhsMask,
                            fields,
                            physicsModel,
                            elementToMaterialMap_,
                            flux);
  }


  //--------------------------------------------------------

  void ATC_Coupling::nodal_projection(const FieldName & fieldName,
                                      const PhysicsModel * physicsModel,
                                      FIELD & field)
  {
    FIELDS rhs;
    rhs[fieldName].reset(nNodes_,field.nCols());
    Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX);
    rhsMask = false;
    rhsMask(fieldName,SOURCE) = true;
    compute_rhs_vector(rhsMask, fields_, rhs, sourceIntegration_, physicsModel);
    const DENS_MAT & B(rhs[fieldName].quantity());

    field = (invNodeVolumes_.quantity())*B;
  }

  // parse_boundary_integration
  //   parses the boundary integration to determine
  //   the type of boundary integration being used
  //--------------------------------------------------
  
  
  BoundaryIntegrationType ATC_Coupling::parse_boundary_integration(int narg,
                                                                   char **arg,
                                                                   const set< pair<int,int> > * boundaryFaceSet)
  {
    
    int argIndex = 0;
    BoundaryIntegrationType myBoundaryIntegrationType = FE_INTERPOLATION;// default
      if (narg > 0) {
        if(strcmp(arg[argIndex],"faceset")==0) {
          argIndex++;
          myBoundaryIntegrationType = FE_QUADRATURE;
          string name(arg[argIndex]);
          boundaryFaceSet = & ( (feEngine_->fe_mesh())->faceset(name));
          set_boundary_face_set(boundaryFaceSet);
        } 
        else if (strcmp(arg[argIndex],"interpolate")==0) {
          myBoundaryIntegrationType = FE_INTERPOLATION;
        }
        else if (strcmp(arg[argIndex],"no_boundary")==0) { 
          myBoundaryIntegrationType = NO_QUADRATURE; 
        }
        else {
          throw ATC_Error("Bad boundary integration type");
        }
      }
    set_boundary_integration_type(myBoundaryIntegrationType);
    return myBoundaryIntegrationType;
  }

}; // namespace ATC
