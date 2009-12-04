// ATC_Transfer headers
#include "ATC_Transfer.h"
#include "FE_Engine.h"
#include "Array.h"
#include "Array2D.h"
#include "ATC_Error.h"
#include "CG.h"
#include "XT_Function.h"
#include "PrescribedDataManager.h"
#include "TimeIntegrator.h"
#include "PhysicsModel.h"
#include "PhysicsModelThermal.h"
#include "PhysicsModelTwoTemperature.h"

namespace ATC {

  ATC_Transfer::ATC_Transfer(void) :
    lammpsInterface_(LammpsInterface::instance()),
    physicsModel_(NULL),
    feEngine_(NULL),
    prescribedDataMgr_(NULL),
    simTime_(0.0),
    stepCounter_(0),
    sampleCounter_(0),
    outputFrequency_(0),
    sampleFrequency_(1),
    outputFrequencyAtom_(0),
    nLocal_(0),
    nLocalTotal_(0),
    nLocalGhost_(0),
    nLocalMask_(0),
    nInternal_(0),
    nGhost_(0),
    nLocalLambda_(0),
    equilibriumStart_(false),
    mdOutputManager_(), 
    bndyIntType_(NO_QUADRATURE), 
    bndyFaceSet_(NULL),
    globalSwitch_(0),
    resetSwitch_(false),
    atomToElementMapType_(LAGRANGIAN),
    atomToElementMapFrequency_(0),
    neighborResetFrequency_(0),
    scalarFlag_(0),
    vectorFlag_(0),
    sizeVector_(0),
    globalFreq_(0),
    extScalar_(0),
    extVector_(0),
    extList_(NULL),
    timeFilterManager_(this),
    groupbit_(0),
    groupbitGhost_(0),
    atomSwitch_(false),
    initialized_(false),
    atomQuadForInternal_(true),
    useLocalizedLambda_(false),
    useLumpedLambda_(false),
    regionID_(-1),
    atomWeightType_(LATTICE),
    xref_(NULL),
    readXref_(false),
    timeIntegrator_(NULL),
    extrinsicModelManager_(this),
    useRestart_(false),
    trackCharge_(false)
  {
    // Defaults
    grow_arrays(lammpsInterface_->nmax());
    fieldMask_.reset(NUM_FIELDS,NUM_FLUX);
    fieldMask_ = false;

    // fe_engine
    feEngine_ = new FE_Engine(this);

    // check to see if lammps has any charges
    double * atomicCharge = lammpsInterface_->atom_charge();
    if (atomicCharge)
      trackCharge_ = true;
  }

  ATC_Transfer::~ATC_Transfer()
  {
    lammpsInterface_->destroy_2d_double_array(xref_);
    if (feEngine_) delete feEngine_;
    if (physicsModel_) delete physicsModel_;
    if (prescribedDataMgr_) delete prescribedDataMgr_;
  }

  //--------------------------------------------------
  // pack_fields
  //   bundle all allocated field matrices into a list
  //   for output needs
  //--------------------------------------------------
  void ATC_Transfer::pack_fields(OUTPUT_LIST & data)
  {
    map<FieldName,int>::const_iterator field;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
        FieldName thisField = field->first;
        string fieldName = field_to_string(thisField);
        string matrixName;
        
        // copy all fields from ATC_Transfer.h
        matrixName = "fields_" + fieldName;
        data[matrixName] = & fields_[thisField];
        matrixName = "dot_fields_" + fieldName;
        data[matrixName] = & dot_fields_[thisField];
        matrixName = "ddot_fields_" + fieldName;
        data[matrixName] = &ddot_fields_[thisField];
        matrixName = "dddot_fields_" + fieldName;
        data[matrixName] = & dddot_fields_[thisField];
        matrixName = "dot_fieldsMD_" + fieldName;
        data[matrixName] = & dot_fieldsMD_[thisField];
        matrixName = "ddot_fieldsMD_" + fieldName;
        data[matrixName] = & ddot_fieldsMD_[thisField];
        matrixName = "dot_dot_fieldsMD_";
        data[matrixName] = & dot_dot_fieldsMD_[thisField];
  
        matrixName = "fieldNdOld_" + fieldName;
        data[matrixName] = & fieldNdOld_[thisField];
        matrixName = "fieldNdFiltered_" + fieldName;
        data[matrixName] = & fieldNdFiltered_[thisField];
        matrixName = "fieldRateNdOld_";
        data[matrixName] = & fieldRateNdOld_[thisField];
        matrixName = "fieldRateNdFiltered_" + fieldName;
        data[matrixName] = & fieldRateNdFiltered_[thisField];
        matrixName = "dot_fieldRateNdOld_" + fieldName;
        data[matrixName] = & dot_fieldRateNdOld_[thisField];
        matrixName = "dot_fieldRateNdFiltered_" + fieldName;
        data[matrixName] = & dot_fieldRateNdFiltered_[thisField];

        matrixName = "auxStorage_" + fieldName;
        data[matrixName] = & auxStorage_[thisField];
    }
  }
  
  //--------------------------------------------------
  // write_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_Transfer::write_restart_data(string fileName, OUTPUT_LIST & data)
  {
    pack_fields(data);
    feEngine_->write_restart_file(fileName,&data);
  }

  //--------------------------------------------------
  // read_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_Transfer::read_restart_data(string fileName, OUTPUT_LIST & data)
  {
    pack_fields(data);
    feEngine_->read_restart_file(fileName,&data);
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
  bool ATC_Transfer::modify(int narg, char **arg)
  {
    FieldName thisField;
    int thisIndex;
    int argIdx;

    bool match = false; 
    
    if (strcmp(arg[0],"transfer")==0) {
      /*! \page man_transfer_output fix_modify AtC transfer output
        \section syntax
        fix_modify AtC transfer output <filename_prefix> <frequency>
        [text | vector_components | tensor_components ]
        - filename_prefix (string) = prefix for data files 
        - frequency (integer) = frequency of output in time-steps 
        - options (keyword/s): \n
        text = creates text output as well as binary Ensight output \n
        vector_components = outputs vectors as scalar components \n
        tensor_components = outputs tensor as scalar components 
        (use this for Paraview)\n
        \section examples
        <TT> fix_modify AtC transfer output heatFE 100 </TT> \n
        <TT> fix_modify AtC transfer output hardyFE 1 text tensor_components </TT> \n
        \section description
        Creates (binary, "gold" format) Ensight output of nodal/mesh data 
        which is transfer/physics specific.
        \section restrictions
        \section related 
        see \ref man_fix_atc
        \section default 
        none
      */
      if (strcmp(arg[1],"output")==0) {
        if (strcmp(arg[2],"nodeset")==0) {
          string nset = arg[3];
          if (strcmp(arg[4],"sum")==0) {
            string field = arg[5];
            pair < string, FieldName >  id(nset,string_to_field(field));
            nsetCompute_[id] = 0.0;
            match = true;
          }
        }
        else {
          outputPrefix_ = arg[2];
          outputFrequency_ = atoi(arg[3]);
          bool text_output = false, vect_comp = false, tensor_comp = false;
          int rank = lammpsInterface_->comm_rank();
          for (int i = 4; i<narg; ++i) {
            if      (strcmp(arg[i],"text")==0) text_output = true;
            else if (strcmp(arg[i],"vector_components")==0) vect_comp = true;
            else if (strcmp(arg[i],"tensor_components")==0) tensor_comp = true;
          }
          if (outputFrequency_>0) {
            if (text_output) {
              feEngine_->initialize_output(rank,outputPrefix_,GNUPLOT);
              if (rank == 0) 
                cout << " ATC:: Warning : text output can create _LARGE_ files\n";
            } 
            else feEngine_->initialize_output(rank,outputPrefix_);
            if (vect_comp) 
              feEngine_->output_manager()
                ->set_option(OUTPUT_VECTOR_COMPONENTS,true);
            if (tensor_comp) 
              feEngine_->output_manager()
                ->set_option(OUTPUT_TENSOR_COMPONENTS,true);
          }
          match = true;
        }
      }
          

      /*! \page man_transfer_atomic_output fix_modify AtC transfer atomic_output
        \section syntax
        fix_modify AtC transfer atomic_output <filename_prefix> <frequency> 
        [text]
        - filename_prefix (string) = prefix for data files
        - frequency (integer) = frequency of output in time-steps
        - option "text" = creates a text version of the output as well as
        the binary Ensight output
        \section examples
        <TT> fix_modify AtC transfer atomic_output heatMD 100 </TT> \n
        <TT> fix_modify AtC transfer atomic_output nanoMD 1 text </TT> \n
        \section description
        Creates (binary, "gold" format) Ensight output of atom (point) data 
        which is transfer/physics specific.
        \section restrictions
        \section related 
        \section default 
        none
      */
      else if (strcmp(arg[1],"atomic_output")==0) {
        outputPrefixAtom_ = arg[2];
        outputFrequencyAtom_ = atoi(arg[3]);
        if (outputFrequencyAtom_>0) {
          if (narg == 5 && strcmp(arg[4],"text")==0) {
            mdOutputManager_.initialize(outputPrefixAtom_,GNUPLOT);
            cout << " ATC:: Warning : text output can create _LARGE_ files\n";
          } 
          else mdOutputManager_.initialize(outputPrefixAtom_);
        }
        match = true;
      }

      /*! \page man_transfer_internal fix_modify AtC transfer internal 
        \section syntax
        fix_modify AtC transfer internal type <atom-type-id>
        - <atom-type-id> = type id for non-fixed atoms internal to the FE mesh
        \section examples
        <TT> fix_modify AtC transfer internal type internal_atoms </TT>
        \section description
        Command to define the atoms to couple to finite elements. 
        This definition is required for a AtC simulation.
        \section restrictions
        \section related 
        see \ref man_transfer_boundary
        \section default 
        none
      */
      else if ((strcmp(arg[1],"internal")==0) && (strcmp(arg[2],"type")==0)) {
        int igroup = lammpsInterface_->find_group(arg[3]);
        groupbit_ |= lammpsInterface_->group_bit(igroup);
        igroups_.insert(igroup);
        match = true;
      }

      /*! \page man_transfer_boundary fix_modify AtC transfer boundary 
        \section syntax
        fix_modify AtC transfer boundary type <atom-type-id>
        - <atom-type-id> = type id for atoms that represent a ficticious
        boundary internal to the FE mesh
        \section examples
        <TT> fix_modify AtC transfer boundary type ghost_atoms </TT>
        \section description
        Command to define the atoms that represent the ficticious 
        boundary internal to the FE mesh. For fully overlapped MD/FE 
        domains with periodic boundary conditions no boundary atoms should
        be defined.
        \section restrictions
        \section related 
        see \ref man_transfer_internal
        \section default 
        none
      */
      else if ((strcmp(arg[1],"boundary")==0) && (strcmp(arg[2],"type")==0)) {
        int igroup = lammpsInterface_->find_group(arg[3]);
        groupbitGhost_ |= lammpsInterface_->group_bit(igroup);
        igroupsGhost_.insert(igroup);
        match = true;
      }

      /*! \page man_internal_quadrature fix_modify AtC transfer internal_quadrature 
        \section syntax
        fix_modify atc transfer internal_quadrature < on | off >
        \section examples
        <TT> fix_modify atc transfer internal_quadrature off </TT>
        \section description
        Command use or not use atomic quadrature on internal elements
        fully filled with atoms. By turning the internal quadrature off
        these elements do not contribute to the governing PDE and the fields
        at the internal nodes follow the weighted averages of the atomic data.
        \section restrictions
        \section related 
        \section default 
        on
      */
      else if (strcmp(arg[1],"internal_quadrature")==0) {
        if (strcmp(arg[2],"on")==0) {
          atomQuadForInternal_ = true;
          match = true;
        }
        else if (strcmp(arg[2],"off")==0) {
          atomQuadForInternal_ = false;
          regionID_ = -1;
          match = true;
        }
        else if (strcmp(arg[3],"off")==0) {
          for (regionID_ = 0; regionID_ < lammpsInterface_->nregion(); regionID_++) 
            if (strcmp(arg[2],lammpsInterface_->region_name(regionID_)) == 0) break;
          if (regionID_ < lammpsInterface_->nregion()) {
            atomQuadForInternal_ = false;
            match = true;
          }
          else
            regionID_ = -1;
        }
      }


      /*! \page man_initial fix_modify AtC transfer initial 
        \section syntax
        fix_modify AtC transfer initial <field> <nodeset> <constant | function>
        - <field> = field name valid for type of physics, temperature | electron_temperature
        - <nodeset> = name of set of nodes to apply boundary condition
        - <constant | function> = value or name of function followed by its
          parameters
        \section examples
        <TT> fix_modify atc transfer initial temperature groupNAME 10. </TT>
        \section description
        Sets the initial values for the specified field at the specified nodes.
        \section restrictions
        keyword 'all' reserved in nodeset name
        \section related 
        see \ref man_transfer_internal
        \section default 
        none
      */
      // set initial conditions
      else if (strcmp(arg[1],"initial")==0) {
        argIdx = 2;
        get_field(arg,argIdx,thisField,thisIndex);
        string nsetName(arg[argIdx++]);
        XT_Function * f = NULL;
        // parse constant
        if (narg == argIdx+1) {
          f = XT_Function_Mgr::instance()->get_constant_function(atof(arg[argIdx]));
        }
        // parse function
        else {
          f = XT_Function_Mgr::instance()->get_function(&(arg[argIdx]),narg-argIdx);
        }
        prescribedDataMgr_->fix_initial_field(nsetName,thisField,thisIndex,f);
        match = true;
      }

      /*! \page man_fix_nodes fix_modify AtC transfer fix 
        \section syntax
        fix_modify AtC transfer fix <field> <nodeset> <constant | function>
        - <field> = field name valid for type of physics
        - <nodeset> = name of set of nodes to apply boundary condition
        - <constant | function> = value or name of function followed by its
          parameters
        \section examples
        <TT> fix_modify AtC transfer fix temperature groupNAME 10. </TT> \n
        <TT> fix_modify AtC transfer fix temperature groupNAME 0 0 0 10.0 0 0 1.0 </TT> \n
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
      else if (strcmp(arg[1],"fix")==0) {
        argIdx = 2;
        get_field(arg,argIdx,thisField,thisIndex);
        string nsetName(arg[argIdx++]);
        XT_Function * f = NULL;
        // parse constant
        if (narg == argIdx+1) {
          f = XT_Function_Mgr::instance()->get_constant_function(atof(arg[argIdx]));
        }
        // parse function
        else {
          f = XT_Function_Mgr::instance()->get_function(&(arg[argIdx]),narg-argIdx);
        }
        prescribedDataMgr_->fix_field(nsetName,thisField,thisIndex,f);
        match = true;
      }

      /*! \page man_unfix_nodes fix_modify AtC transfer unfix 
        \section syntax
        fix_modify AtC transfer unfix <field> <nodeset> 
        - <field> = field name valid for type of physics
        - <nodeset> = name of set of nodes 
        \section examples
        <TT> fix_modify AtC transfer unfix temperature groupNAME </TT>
        \section description
        Removes constraint on field values for specified nodes.
        \section restrictions
        keyword 'all' reserved in nodeset name
        \section related 
        see \ref man_fix_nodes
        \section default 
        none
      */
      else if (strcmp(arg[1],"unfix")==0) {
        argIdx = 2;
        get_field(arg,argIdx,thisField,thisIndex);
        string nsetName(arg[argIdx++]);
        prescribedDataMgr_->unfix_field(nsetName,thisField,thisIndex);
        match = true;
      }

    /*! \page man_source fix_modify AtC transfer source
      \section syntax
       fix_modify AtC transfer source <field> <element_set> <value | function>
        - <field> = field name valid for type of physics
        - <element_set> = name of set of elements 
      \section examples
      <TT> fix_modify atc transfer source temperature middle temporal_ramp 10. 0. </TT>
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
      else if (strcmp(arg[1],"source")==0) {
        argIdx = 2;
        get_field(arg,argIdx,thisField,thisIndex);
        string esetName(arg[argIdx++]);
        XT_Function * f = NULL;
        // parse constant
        if (narg == argIdx+1) {
          f = XT_Function_Mgr::instance()->get_constant_function(atof(arg[argIdx]));
        }
        // parse function
        else {
          f = XT_Function_Mgr::instance()->get_function(&(arg[argIdx]),narg-argIdx);
        }
        prescribedDataMgr_->fix_source(esetName,thisField,thisIndex,f);
        fieldMask_(thisField,PRESCRIBED_SOURCE) = true;
        match = true;
      }

    /*! \page man_remove_source fix_modify AtC transfer remove_source
      \section syntax
      fix_modify AtC transfer remove_source <field> <element_set>
        - <field> = field name valid for type of physics
        - <element_set> = name of set of elements 
      \section examples
      <TT> fix_modify atc transfer remove_source temperature groupNAME </TT>
      \section description
      Remove a domain source.
      \section restrictions
      keyword 'all' reserved in element_set name
      \section related
      see \ref man_source
      \section default
    */
      else if (strcmp(arg[1],"remove_source")==0) {
        argIdx = 2;
        get_field(arg,argIdx,thisField,thisIndex);
        string esetName(arg[argIdx++]);
        prescribedDataMgr_->unfix_source(esetName,thisField,thisIndex);
        fieldMask_(thisField,PRESCRIBED_SOURCE) = false;
        match = true;
      }
  
    /*! \page man_fix_flux fix_modify AtC transfer fix_flux
      \section syntax
       fix_modify AtC transfer fix_flux <field> <face_set> <value | function>
        - <field> = field name valid for type of physics, temperature | electron_temperature
        - <face_set> = name of set of element faces
      \section examples
       <TT> fix_modify atc transfer fix_flux temperature faceSet 10.0 </TT> \n
       
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
      else if (strcmp(arg[1],"fix_flux")==0) {
        argIdx = 2;
        get_field(arg,argIdx,thisField,thisIndex);
        string fsetName(arg[argIdx++]);
        XT_Function * f = NULL;
        // parse constant
        if (narg == argIdx+1) {
          f = XT_Function_Mgr::instance()->get_constant_function(atof(arg[argIdx]));
        }
        // parse function
        else {
          f = XT_Function_Mgr::instance()->get_function(&(arg[argIdx]),narg-argIdx);
        }
        prescribedDataMgr_->fix_flux(fsetName,thisField,thisIndex,f);
        fieldMask_(thisField,PRESCRIBED_SOURCE) = true;
        match = true;
      }

    /*! \page man_unfix_flux fix_modify AtC transfer unfix_flux
      \section syntax
      fix_modify AtC transfer fix_flux <field> <face_set> <value | function>
        - <field> = field name valid for type of physics, temperature | electron_temperature
        - <face_set> = name of set of element faces
      \section examples
       <TT> fix_modify atc transfer unfix_flux temperature faceSet  </TT> \n
       
      \section description
       Command for removing prescribed normal fluxes e.g. heat_flux, stress. 
      \section restrictions 
      \section related
      see \ref man_unfix_flux
      \section default
    */
      else if (strcmp(arg[1],"unfix_flux")==0) {
        argIdx = 2;
        get_field(arg,argIdx,thisField,thisIndex);
        string fsetName(arg[argIdx++]);
        prescribedDataMgr_->unfix_flux(fsetName,thisField,thisIndex);
        fieldMask_(thisField,PRESCRIBED_SOURCE) = false;
        match = true;
      }

      // pass off to time filter
      else if (strcmp(arg[1],"filter")==0) {
        match = timeFilterManager_.modify(narg-2,&arg[2]);
        if (match && timeIntegrator_) timeIntegrator_->force_reset();
      }

    }

    // pass off to extrinsic model
    else if (strcmp(arg[0],"extrinsic")==0) {
      match = extrinsicModelManager_.modify(narg-1,&arg[1]);
    }

    /*! \page man_atom_element_map fix_modify AtC transfer atom_element_map
      \section syntax
      fix_modify AtC transfer atom_element_map  <eulerian|lagrangian> <frequency> \n
      - frequency (int) : frequency of updating atom-to-continuum maps based on the
      current configuration - only for eulerian 
      \section examples
      <TT> fix_modify atc transfer atom_element_map eulerian 100 </TT>
      \section description
      Changes frame of reference from eulerian to lagrangian and sets the
      frequency for which the map from atoms to elements is reformed and 
      all the attendant data is recalculated.
      \section restrictions 
      Cannot change map type after initialization.
      \section related
      \section default
      lagrangian
    */
    else if (strcmp(arg[0],"atom_element_map")==0) {
      if (strcmp(arg[1],"eulerian")==0) {
        atomToElementMapType_ = EULERIAN;
        atomToElementMapFrequency_ = atoi(arg[2]);
      } 
      else {
        atomToElementMapType_ = LAGRANGIAN;
        atomToElementMapFrequency_ = 0;
      } 
      match = true;
    }

    /*! \page man_neighbor_reset_frequency fix_modify AtC transfer neighbor_reset_frequency
        \section syntax
        fix_modify AtC transfer neighbor_reset_frequency <frequency>
        - frequency (int) : reset neighbor map for AtC analysis every "frequency"  steps
        \section examples
        <TT> fix_modify AtC transfer neighbor_reset_frequency 1000 </TT>
        \section description
        This command allows the control of how often the lists 
        of which atoms pairs are in the support of a node
        are reformed
        \section restrictions
        Only currently relevant to Hardy post-processing see \ref man_fix_atc
        \section related
        This command is different from \ref man_atom_element_map 
        in that it only concerns neighbor lists internal to ATC not
        the map of atom locations to elements which define the shape function
        values at atoms.
        \section default
        The default is for the internal neighbor lists to be recalculated on
        the next output step (or sample step if time filtering is used).
      */
    else if (strcmp(arg[0],"neighbor_reset_frequency")==0) {
      neighborResetFrequency_ = atoi(arg[1]);
      match = true;
    }

    /*! \page man_read_restart fix_modify AtC transfer read_restart
      \section syntax
      fix_modify AtC transfer read_restart [file_name]  \n
      \section examples
      <TT> fix_modify AtC transfer read_restart ATC_state </TT> \n
      \section description
      Reads the current state of the fields from a named text-based restart 
      file.
      \section restrictions 
      The restart file only contains fields and their time derivatives.
      The reference positions of the atoms and the commands that initialize
      the fix are not saved e.g. an identical mesh containing the same atoms 
      will have to be recreated. 
      \section related
      see write_restart \ref man_write_restart
      \section default
      none
    */
    else if (strcmp(arg[0],"read_restart")==0) {
      restartFileName_  = arg[1];
      useRestart_ = true;
      match = true;
    }
    /*! \page man_write_restart fix_modify AtC transfer write_restart
      \section syntax
      fix_modify AtC transfer write_restart [file_name]  \n
      \section examples
      <TT> fix_modify AtC transfer write_restart restart.mydata </TT> \n
      \section description
      Dumps the current state of the fields to a named text-based restart file.
      This done when the command is invoked and not repeated, unlike the 
      similar lammps command.
      \section restrictions 
      The restart file only contains fields and their time derivatives.
      The reference positions of the atoms and the commands that initialize
      the fix are not saved e.g. an identical mesh containing the same atoms 
      will have to be recreated. 
      \section related
      see read_restart \ref man_read_restart
      \section default
      none
    */
    else if (strcmp(arg[0],"write_restart")==0) {
      string restartFileName(arg[1]);
      OUTPUT_LIST data;
      write_restart_data(restartFileName,data);
      match = true;
    }


    // pass off to fe engine
    if (!match) {
      match = feEngine_->modify(narg, arg);
    }
  
    return match;
  
  }

  //--------------------------------------------------
  // parse_boundary_integration
  //   parses the boundary integration to determine
  //   the type of boundary integration being used
  //--------------------------------------------------
  ATC::ATC_Transfer::BoundaryIntegrationType ATC_Transfer::parse_boundary_integration(int narg,
                                                                                      char **arg,
                                                                                      const set< pair<int,int> > * boundaryFaceSet)
  {
    int argIndex = 0;
    BoundaryIntegrationType myBoundaryIntegrationType;
    if (narg > 0) {
      if(strcmp(arg[argIndex],"faceset")==0) {
        argIndex++;
        myBoundaryIntegrationType = FE_QUADRATURE;
        string name(arg[argIndex]);
        boundaryFaceSet = & ( (feEngine_->get_feMesh())->get_faceset(name));
        set_boundary_face_set(boundaryFaceSet);
      } 
      else if (strcmp(arg[argIndex],"interpolate")==0) {
        myBoundaryIntegrationType = FE_INTERPOLATION;
      }
      else if (strcmp(arg[argIndex],"no_boundary")==0) {
        myBoundaryIntegrationType = NO_QUADRATURE;
      }
      else {
        throw ATC_Error(0,"Bad boundary integration type");
      }
    }
    else { // default is interpolation
      myBoundaryIntegrationType = FE_INTERPOLATION;
    }
    
    set_boundary_integration_type(myBoundaryIntegrationType);
    return myBoundaryIntegrationType;
  }

  //--------------------------------------------------
  // helper function for parser
  // handles : "displacement x"  and "temperature" by indexing argIdx
  // for fluxes : only normal fluxes can be prescribed
  //--------------------------------------------------
  void ATC_Transfer::get_field(/*const*/ char ** args, int & argIdx,
                               FieldName & thisField, int & thisIndex)
  {
    string thisName = args[argIdx++];
    thisField = string_to_field(thisName);
    map<FieldName,int>::const_iterator iter = fieldSizes_.find(thisField);
    if (iter == fieldSizes_.end()) {
      throw ATC_Error(0,"Bad field name: "+thisName);
    }
    string thisDim = args[argIdx];
    thisIndex = 0;
    if (string_to_index(thisDim,thisIndex)) {
      if ( !( thisIndex < fieldSizes_[thisField]) ) {
        throw ATC_Error(0,"Bad field dimension "+thisDim);
      }
      argIdx++;
    }
  }

  void ATC_Transfer::init_filter()
  {
    if (timeFilterManager_.need_reset()) {
      map<FieldName,int>::const_iterator field;
      for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
        FieldName thisField = field->first;
        int thisSize = field->second;
        
        // Allocate restricted fields
        fieldRateNdFiltered_[thisField].reset(nNodes_,thisSize);
        fieldRateNdOld_[thisField].reset(nNodes_,thisSize);
        dot_fieldRateNdFiltered_[thisField].reset(nNodes_,thisSize);
        dot_fieldRateNdOld_[thisField].reset(nNodes_,thisSize);
      }
    }
  }

  void ATC_Transfer::set_xref()
  {
    double **x = lammpsInterface_->xatom();
    for (int i = 0; i < lammpsInterface_->nmax(); i++) {
      for (int j = 0; j < nsd_; j++) {
        xref_[i][j] = x[i][j];
      }
    }
    Xprd_ = lammpsInterface_->domain_xprd();
    Yprd_ = lammpsInterface_->domain_yprd();
    Zprd_ = lammpsInterface_->domain_zprd();
    XY_ = lammpsInterface_->domain_xy();
    XZ_ = lammpsInterface_->domain_xz();
    YZ_ = lammpsInterface_->domain_yz();

  }

  void ATC_Transfer::initialize()
  { 
    // initialized_ is set to true by derived class initialize()
    int i, j;
    if (!initialized_) {
      if (!trackCharge_) {
        fieldSizes_.erase(CHARGE_DENSITY);
        prescribedDataMgr_->remove_field(CHARGE_DENSITY);
        for (i = 0; i < NUM_FLUX; i++)
          fieldMask_(CHARGE_DENSITY,i) = false;
      }

      // Initialize finite elements shape function, mapping variables
      feEngine_->initialize();
      
      nNodes_ = feEngine_->get_nNodes(); 
      nsd_ = feEngine_->get_nsd();
      if (nsd_!=lammpsInterface_->dimension()) 
        throw ATC_Error(0,"Spatial dimensions inconsistent between LAMMPS and ATC");
    
      // set reference position of atoms
      set_xref();    

    } 
 
    // determine which elements are strictly in the MD region and mask out for intrinsic fields
    // start by finding local atom counts
    reset_nlocal();
    int nelts = feEngine_->get_feMesh()->get_nElements();
    elementMask_.reset(nelts);
    if (atomQuadForInternal_) {
      elementMask_ = true;
    }
    else {
      // determine which elements contain internal atoms
      Array<int> hasInternal(nelts);
      Array<int> hasInternal_local(nelts);
      hasInternal_local = 0;
      for (int i = 0; i < nLocal_; ++i) {
        DENS_VEC coords(3);
        coords(0) = atomicCoords_(0,i);
        coords(1) = atomicCoords_(1,i);
        coords(2) = atomicCoords_(2,i);
        int eltID = feEngine_->get_feMesh()->map_to_element(coords);
        hasInternal_local(eltID) = 1;
      }
      // swap contributions
      lammpsInterface_->logical_or(hasInternal_local.get_ptr(),
                                   hasInternal.get_ptr(),hasInternal.size());
      // determine which elements contain ghost atoms
      Array<int> hasGhost(nelts);
      Array<int> hasGhost_local(nelts);
      hasGhost_local = 0;
      for (int i = 0; i < nLocalGhost_; ++i) {
        DENS_VEC coords(3);
        coords(0) = ghostAtomCoords_(0,i);
        coords(1) = ghostAtomCoords_(1,i);
        coords(2) = ghostAtomCoords_(2,i);
        int eltID = feEngine_->get_feMesh()->map_to_element(coords);
        hasGhost_local(eltID) = 1;
      }
      //swap contributions
      lammpsInterface_->logical_or(hasGhost_local.get_ptr(),
                                   hasGhost.get_ptr(),hasGhost.size());
      for (int i = 0; i < nelts; ++i)
        elementMask_(i) = !hasInternal(i) || hasGhost(i);
    }
    set<int> & nullElements = feEngine_->null_elements();
    set<int>::const_iterator iset;
    for (iset = nullElements.begin(); iset != nullElements.end(); iset++) {
      int ielem = *iset;
      elementMask_(ielem) = false;
    }
  
    // set atomic weights
    reset_shape_functions(); // creates atom->element map
    set_atomic_weights();

    // group atoms by material
    int nMatls = physicsModel_->get_nMaterials();
    atomMaterialGroups_.reset(nMatls);
    for (int i = 0; i < nLocal_; i++) {
      int matId = elementToMaterialMap_(atomToElementMap_(i));
      atomMaterialGroups_(matId).insert(i);// not internalToAtom_
    }
    if (atomToElementMapType_ == EULERIAN && nMatls >1 ) {
      throw ATC_Error(0," multiple materials & eulerian map not supported");
    }
    
    // compute necessary restructuring matrices
    maskMat_.reset(nNodes_,nNodes_);
    if (!useLocalizedLambda_) {
      for (i = 0; i < nNodes_; ++i)
        if (shpWeight_(i,i) > 0.)
          maskMat_(i,i) = 1.;
    }
    else { // select nodes for thermostat
      // 1) their shape functions have both mask and ghost atoms in their support
      DENS_VEC hasInternalLocal(nNodes_);
      DENS_VEC hasInternal(nNodes_);
      if (nLocal_>0) hasInternalLocal = shpFcn_.col_sum();    
      lammpsInterface_->allsum(hasInternalLocal.get_ptr(),hasInternal.get_ptr(),nNodes_);
     
      DENS_VEC hasGhostLocal(nNodes_);
      DENS_VEC hasGhost(nNodes_);
      if (nLocalGhost_>0) hasGhostLocal = shpFcnGhost_.col_sum();
      lammpsInterface_->allsum(hasGhostLocal.get_ptr(),hasGhost.get_ptr(),nNodes_);
      
      for (i = 0; i < nNodes_; ++i) {
        if (hasInternal(i) > 0. && hasGhost(i) > 0.) {
          maskMat_(i,i) = 1.;
        }
      }

      // 2) they are in a specified boundary faceset
      set<string>::const_iterator iter;
      for (iter = boundaryFaceNames_.begin(); iter != boundaryFaceNames_.end(); iter++) {
        set<int> nodeSet;
        feEngine_->get_feMesh()->faceset_to_nodeset(*iter,nodeSet);
        set<int>::const_iterator nodeIter;
        for (nodeIter = nodeSet.begin(); nodeIter != nodeSet.end(); nodeIter++) {
          maskMat_(*nodeIter,*nodeIter) = 1.;
        }
      }
    }
    
    // determine non-zero values and row-column indices
    // also used in creating overlap map and mapping from atoms to elements
    // note we assume the overlap map is time-independent
    nodeToOverlapMap_.reset(nNodes_);
    nNodeOverlap_ = 0;
    for (i = 0; i < nNodes_; i++) {
      if (maskMat_(i,i) > 0.) {
        nodeToOverlapMap_(i) = nNodeOverlap_;
        nNodeOverlap_++;
      }
      else {
        nodeToOverlapMap_(i) = -1;
      }
    }
  
    overlapToNodeMap_.reset(nNodeOverlap_);
    int counter = 0;
    for (i = 0; i < nNodes_; i++) {
      if (nodeToOverlapMap_(i)!=-1) {
        overlapToNodeMap_(counter) = i;
        counter++;
      }
    }
    
    // reset local overlap shape function storage
    reset_NhatOverlap();
    
    // Shape sparse matrix used in Hoover-type  thermostats
    // first get local pattern sized nNodeOverlap X nNodeOverlap
    //compute_consistent_md_mass_matrix(NhatOverlap_,tmp);
    DENS_MAT tmpLocal(nNodeOverlap_,nNodeOverlap_);
    DENS_MAT tmp(nNodeOverlap_,nNodeOverlap_);
    if (nLocalLambda_>0) {
      DENS_VEC T(nLocalLambda_);
      T = 1.0;
      for (i = 0; i < nNodeOverlap_; i++) {
        for (j = 0; j < nNodeOverlap_; j++) {
          for (int iatom = 0; iatom < nLocalLambda_; iatom++) {
            tmpLocal(i,j) += NhatOverlap_(iatom,i)*T(iatom)*NhatOverlap_(iatom,j);
          }
        }
      }
    }
    // second accumulate total pattern
    lammpsInterface_->allsum(tmpLocal.get_ptr(), tmp.get_ptr(), tmp.size());
    // third extract non-zero entries & construct
    M_T_Template.reset(nNodeOverlap_,nNodeOverlap_);
    for (i = 0; i < nNodeOverlap_; i++) {
      for (j = 0; j < nNodeOverlap_; j++) {
        if (abs(tmp(i,j))>0) {
          M_T_Template.add(i,j,0.);
        }  
      }
    }
    M_T_Template.compress();
    
    // \int_\Omega N_I dV
    NodeVolumes_.reset(nNodes_,nNodes_);
    invNodeVolumes_.reset(nNodes_,nNodes_);
    feEngine_->compute_lumped_mass_matrix(NodeVolumes_);
    invNodeVolumes_ = NodeVolumes_.inv();

    // this is exact only for uniform meshes and certain types of atomic weights
    // \int_{\Omega_MD} N_I dV = \sum_\alpha N_I\alpha V_\alpha
    //     DENS_VEC atom_weights_vec(atomicWeights_.size(), 
    //                                     atomicWeights_.get_ptr()); 
    DENS_VEC subsetNodeVolume(NodeVolumes_.size());
    DENS_VEC unitVec(atomicWeights_.size());
    unitVec = 1.0;
    restrict_unscaled(unitVec, subsetNodeVolume);
    fluxMask_.reset(invNodeVolumes_* subsetNodeVolume);
    DIAG_MAT id(fluxMask_.nRows(),fluxMask_.nCols()); 
    id = 1.0;
    fluxMaskComplement_ = id + -1.0*fluxMask_;

    // set flux masks for nodes we can tell by geometry
    DENS_VEC hasInternalLocal(nNodes_);
    DENS_VEC hasInternal(nNodes_);
    if (nLocal_>0) hasInternalLocal = shpFcn_.col_sum();    
    lammpsInterface_->allsum(hasInternalLocal.get_ptr(),hasInternal.get_ptr(),nNodes_);
    
    DENS_VEC hasGhostLocal(nNodes_);
    DENS_VEC hasGhost(nNodes_);
    if (nLocalGhost_>0) hasGhostLocal = shpFcnGhost_.col_sum();
    lammpsInterface_->allsum(hasGhostLocal.get_ptr(),hasGhost.get_ptr(),nNodes_);
    
    for (i = 0; i < nNodes_; ++i) {
      if (hasInternal(i) > 0.) {
        if (hasGhost(i) == 0.) {
          fluxMask_(i,i) = 1.;
          fluxMaskComplement_(i,i) = 0.;
        }
      }
      else {
        fluxMask_(i,i) = 0.;
        fluxMaskComplement_(i,i) = 1.;
      }
    }

    // Initialize fields and set initial conditions
    if (!initialized_) {
      nodeType_.reset(nNodes_);
      for (int iNode = 0; iNode < nNodes_; ++iNode)
        nodeType_(iNode) = check_shape_function_type(iNode);
      
      map<FieldName,int>::const_iterator field;
      for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
        FieldName thisField = field->first;
        int thisSize = field->second;
      
        // Allocate fields, initialize to default values, set up initial schedule
        fields_[thisField].reset(nNodes_,thisSize);
        dot_fields_[thisField].reset(nNodes_,thisSize);
        ddot_fields_[thisField].reset(nNodes_,thisSize);
        dddot_fields_[thisField].reset(nNodes_,thisSize);
          
        dot_fieldsMD_[thisField].reset(nNodes_,thisSize);
        ddot_fieldsMD_[thisField].reset(nNodes_,thisSize);
        dot_dot_fieldsMD_[thisField].reset(nNodes_,thisSize);

        // Allocate restricted fields
        fieldRateNdFiltered_[thisField].reset(nNodes_,thisSize);
        fieldRateNdOld_[thisField].reset(nNodes_,thisSize);
        dot_fieldRateNdFiltered_[thisField].reset(nNodes_,thisSize);
        dot_fieldRateNdOld_[thisField].reset(nNodes_,thisSize);
    
        fieldNdFiltered_[thisField].reset(nNodes_,thisSize);
        fieldNdOld_[thisField].reset(nNodes_,thisSize);

        // Allocate auxilliary storage
        auxStorage_[thisField].reset(nNodes_,thisSize);
      
        // Dimension finite element rhs matrix
        rhs_[thisField].reset(nNodes_,thisSize);
        rhsAtomDomain_[thisField].reset(nNodes_,thisSize);
        sources_[thisField].reset(nNodes_,thisSize);
        extrinsicSources_[thisField].reset(nNodes_,thisSize);
        boundaryFlux_[thisField].reset(nNodes_,thisSize);

        if (is_intrinsic(thisField)) {
          // compute FE mass matrix in full domain
          massMatInv_[thisField].reset(nNodes_,nNodes_);
          Array<FieldName> massMask(1);
          massMask(0) = thisField;
          feEngine_->compute_lumped_mass_matrix(massMask,fields_,physicsModel_,
                                                elementToMaterialMap_,massMats_,
                                                &elementMask_);

          // atomic quadrature for FE mass matrix in atomic domain
          map<FieldName,DIAG_MAT> massMatAq;
          feEngine_->compute_lumped_mass_matrix(massMask,fields_,physicsModel_,
                                                atomMaterialGroups_,atomicWeightsMask_,
                                                shpFcnMasked_,massMatAq);
          // remove contributions from overlap by approximate quadrature
          // and fully remove contributions from internal
          massMats_[thisField] -= massMatAq[thisField];
          for (int iNode = 0; iNode < nNodes_; iNode++)
            if (nodeType_(iNode)==MD_ONLY)
              massMats_[thisField](iNode,iNode) = 0.;
          
          // set up mass MD matrices
          massMatMDInv_[thisField].reset(nNodes_,nNodes_);
          compute_md_mass_matrix(thisField,massMatsMD_);
          massMats_[thisField] += massMatsMD_[thisField];
          
          // compute inverse mass matrices since we're using lumped masses
          for (int iNode = 0; iNode < nNodes_; iNode++) {
            if (fabs(massMatsMD_[thisField](iNode,iNode))>0) { 
              massMatMDInv_[thisField](iNode,iNode) = 1./massMatsMD_[thisField](iNode,iNode);
            }
            if (fabs(massMats_[thisField](iNode,iNode))>0) {
              massMatInv_[thisField](iNode,iNode) = 1./massMats_[thisField](iNode,iNode);
            }
          }
        }
        else {
          // no MD mass matrices needed, regular matrices computed in extrinsic model
          massMats_[thisField].reset(nNodes_,nNodes_);
          massMatInv_[thisField].reset(nNodes_,nNodes_);
        }
      }

      // compute inverse consistent mass for the continuity equation
      SPAR_MAT consistentMass;
      feEngine_->compute_mass_matrix(consistentMass);
      consistentMass.compress();
      consistentMassInverse_ = inv(consistentMass.dense_copy());

      
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

    // set lists of fixed nodes per field
    map<FieldName,int>::const_iterator field;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      int thisSize = field->second;
      isFixedNode_[thisField].reset(nNodes_,thisSize);
      for (int i = 0; i < nNodes_; i++)
        for (int j = 0; j < thisSize; j++)
          isFixedNode_[thisField](i,j) = prescribedDataMgr_->is_fixed(i,thisField,j);
    }

    return;
  }

  void ATC_Transfer::pre_final_integrate()
  {
    // atomSwtich is triggered by a lammps call to pre_exchange 
    // (see ATC_Transfer.h)
    int localSwitch = 0;
    if (atomSwitch_ || resetSwitch_) localSwitch = 1;
    lammpsInterface_->int_allmax(&localSwitch,&globalSwitch_);  
    if (globalSwitch_>0) {
      reset_nlocal();
      reset_shape_functions();
      set_atomic_weights();
      reset_NhatOverlap();
      atomSwitch_ = false;
      globalSwitch_ = 0;
    }
    else if (atomToElementMapType_ == EULERIAN
             && stepCounter_ % atomToElementMapFrequency_ == 0 ) {
      reset_coordinates();
      reset_shape_functions();
      set_atomic_weights();
      reset_NhatOverlap();
    }
    
  }

  void ATC_Transfer::init_integrate_velocity() 
  {
    int nlocal = lammpsInterface_->nlocal();
    int *type = lammpsInterface_->atom_type();
    int *mask = lammpsInterface_->atom_mask();
    double *mass = lammpsInterface_->atom_mass();
    double **v = lammpsInterface_->vatom();
    double **f = lammpsInterface_->fatom();
    double dtf = 0.5 * lammpsInterface_->dt() * lammpsInterface_->ftm2v();

    if (mass) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit_) {
          double dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
      }
  
    } else {
      double *rmass = lammpsInterface_->atom_rmass();
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit_) {
          double dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
      }
    }
  }

  void ATC_Transfer::init_integrate_position() 
  {
    int nlocal = lammpsInterface_->nlocal();
    double **x = lammpsInterface_->xatom();
    double **v = lammpsInterface_->vatom();
    double dtv = lammpsInterface_->dt();
    int *mask  = lammpsInterface_->atom_mask();

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit_) {
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    }
  }

  void ATC_Transfer::final_integrate() 
  {
    int nlocal = lammpsInterface_->nlocal();
    int *type = lammpsInterface_->atom_type();
    int *mask = lammpsInterface_->atom_mask();
    double *mass = lammpsInterface_->atom_mass();
    double **v = lammpsInterface_->vatom();
    double **f = lammpsInterface_->fatom();
    double dtf = 0.5 * lammpsInterface_->dt() * lammpsInterface_->ftm2v();

    if (mass) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit_) {
          double dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
      }
 
    } else {
      double *rmass = lammpsInterface_->atom_rmass();
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit_) {
          double dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
      }
    }
  }


  void ATC_Transfer::finish()
  {
    // possibly resetting flags and deallocating corresponding memory

    // FE Engine
    feEngine_->finish();

    return;
  }

  //--------------------------------------------------------------
  // create_physics_model
  // - method to create physics model 
  //--------------------------------------------------------------
  void ATC_Transfer::create_physics_model(const PhysicsType & physicsType,
                                          string matFileName)
  {
    if (physicsModel_) {
      throw ATC_Error(0,"Attempted to create PhysicsModel multiple times in ATC_Transfer");
    }
    // Create PhysicsModel based on physicsType
    switch (physicsType) {
    case NO_PHYSICS :
      break;
    case THERMAL :
      physicsModel_ = new PhysicsModelThermal(matFileName, this);
      break;
    default:
      throw ATC_Error(0,"Unknown physics type in ATC_Transfer::create_physics_model");
    }

  }

  //--------------------------------------------------------------
  /** method to trigger construction of mesh data after mesh construction */
  //--------------------------------------------------------------
  void ATC_Transfer::initialize_mesh_data(void)
  {
    int nelts = feEngine_->get_feMesh()->get_nElements();
    elementToMaterialMap_.reset(nelts);
    elementToMaterialMap_ = 0;

    construct_prescribed_data_manager();
  }

  //--------------------------------------------------------------
  /** method to add new fields to the included list */
  //--------------------------------------------------------------
  void ATC_Transfer::add_fields(map<FieldName,int> & newFieldSizes)
  {
    map<FieldName,int>::const_iterator field;
    for (field = newFieldSizes.begin(); field!=newFieldSizes.end(); field++) {
      FieldName thisField = field->first;
      int thisSize = field->second;
      if (fieldSizes_.find(thisField)==fieldSizes_.end()) {
          fieldSizes_[thisField] = thisSize;
      }
    }
  }
  
}; // namespace ATC
