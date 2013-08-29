// ATC headers
#include "ATC_Method.h"
#include "LammpsInterface.h"
#include "FE_Engine.h"
#include "Array.h"
#include "Array2D.h"
#include "ATC_Error.h"
#include "Function.h"
#include "PrescribedDataManager.h"
#include "TimeIntegrator.h"
#include "PhysicsModel.h"
#include "PerAtomQuantityLibrary.h"
#include "TransferLibrary.h"
#include "KernelFunction.h"
#include "Utility.h"
#include "FieldManager.h"

#include <fstream>
#include <sstream>
#include <iostream>

using ATC_Utility::sgn;
using ATC_Utility::to_string;
using ATC_Utility::is_dbl;

using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::string;
using std::map;
using std::set;
using std::vector;
using std::pair;

namespace ATC {

  ATC_Method::ATC_Method(string groupName, double ** & perAtomArray, LAMMPS_NS::Fix * thisFix) :
    nodalAtomicVolume_(NULL),
    needReset_(true),
    lammpsInterface_(LammpsInterface::instance()),
    interscaleManager_(this),
    timeFilterManager_(this),
    integrateInternalAtoms_(false),
    atomTimeIntegrator_(NULL),
    ghostManager_(this),
    feEngine_(NULL),
    initialized_(false),
    meshDataInitialized_(false),
    localStep_(0),
    sizeComm_(8), // 3 positions + 1 material id * 2 for output
    atomCoarseGrainingPositions_(NULL),
    atomGhostCoarseGrainingPositions_(NULL),
    atomProcGhostCoarseGrainingPositions_(NULL),
    atomReferencePositions_(NULL),
    nsd_(lammpsInterface_->dimension()),
    xref_(NULL),
    readXref_(false),
    needXrefProcessorGhosts_(false),
    trackDisplacement_(false),
    needsAtomToElementMap_(true),
    atomElement_(NULL),
    atomGhostElement_(NULL),
    internalElementSet_(""),
    atomMasses_(NULL),
    atomPositions_(NULL),
    atomVelocities_(NULL),
    atomForces_(NULL),
    parallelConsistency_(true),
    outputNow_(false),
    outputTime_(true),
    outputFrequency_(0),
    sampleFrequency_(0),
    sampleCounter_(0),
    peScale_(1./(lammpsInterface_->mvv2e())),
    keScale_(1.),
    scalarFlag_(0),
    vectorFlag_(0),
    sizeVector_(0),
    scalarVectorFreq_(0),
    sizePerAtomCols_(4), 
    perAtomOutput_(NULL),
    perAtomArray_(perAtomArray),
    extScalar_(0),
    extVector_(0),
    extList_(NULL),
    thermoEnergyFlag_(0),
    atomVolume_(NULL),
    atomicWeightsWriteFlag_(false),
    atomicWeightsWriteFrequency_(0),
    atomWeightType_(LATTICE),
    domainDecomposition_(REPLICATED_MEMORY),
    groupbit_(0),
    groupbitGhost_(0),
    needProcGhost_(false),
    groupTag_(groupName),
    nLocal_(0),
    nLocalTotal_(0),
    nLocalGhost_(0),
    atomToElementMapType_(LAGRANGIAN),
    atomToElementMapFrequency_(0),
    regionID_(-1),
    mdMassNormalization_(false),
    kernelBased_(false),
    kernelOnTheFly_(false),
    kernelFunction_(NULL),
    bondOnTheFly_(false),
    accumulant_(NULL),
    accumulantMol_(NULL),
    accumulantMolGrad_(NULL),
    accumulantWeights_(NULL),
    accumulantInverseVolumes_(&invNodeVolumes_),
    accumulantBandwidth_(0),
    useRestart_(false),
    hasRefPE_(false),
    setRefPE_(false),
    setRefPEvalue_(false),
    refPEvalue_(0.),
    readRefPE_(false),
    nodalRefPotentialEnergy_(NULL),
    simTime_(0.0),
    stepCounter_(0)
  {
    lammpsInterface_->print_msg_once("version "+version());    
    lammpsInterface_->set_fix_pointer(thisFix);
    interscaleManager_.set_lammps_data_prefix();
    grow_arrays(lammpsInterface_->nmax());
    feEngine_ = new FE_Engine(lammpsInterface_->world());

    
    lammpsInterface_->create_compute_pe_peratom();
  }

  ATC_Method::~ATC_Method()
  {
    lammpsInterface_->destroy_2d_double_array(xref_);
    lammpsInterface_->destroy_2d_double_array(perAtomOutput_);
    if (atomTimeIntegrator_) delete atomTimeIntegrator_;
    if (feEngine_) delete feEngine_;
  }

  //--------------------------------------------------
  // pack_fields
  //   bundle all allocated field matrices into a list
  //   for output needs
  //--------------------------------------------------
  void ATC_Method::pack_fields(RESTART_LIST & data)
  {
    map<FieldName,int>::const_iterator field;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      string fieldName = field_to_string(thisField);
      string matrixName;
      // copy all fields from ATC_Method.h
      matrixName = "fields_" + fieldName;
      data[matrixName] = & fields_[thisField].set_quantity();
      matrixName = "dot_fields_" + fieldName;
      data[matrixName] = & dot_fields_[thisField].set_quantity();
      matrixName = "ddot_fields_" + fieldName;
      data[matrixName] = & ddot_fields_[thisField].set_quantity();
      matrixName = "dddot_fields_" + fieldName;
      data[matrixName] = & dddot_fields_[thisField].set_quantity();
      matrixName = "NodalAtomicFieldsRoc_" + fieldName;
      data[matrixName] = & nodalAtomicFieldsRoc_[thisField].set_quantity();
    }
  }
  
  //--------------------------------------------------
  // write_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_Method::write_restart_data(string fileName, RESTART_LIST & data)
  {
    pack_fields(data); 
    feEngine_->write_restart_file(fileName,&data);
  }

  //--------------------------------------------------
  // read_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_Method::read_restart_data(string fileName, RESTART_LIST & data)
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
  bool ATC_Method::modify(int narg, char **arg)
  {
    int argIdx=0;
    bool match = false; 
    
    // gateways to other modules e.g. extrinsic, control, mesh
    // pass off to fe engine
    if (strcmp(arg[argIdx],"mesh")==0) {
      match = feEngine_->modify(narg, arg); 
      if (feEngine_->has_mesh()  && !meshDataInitialized_) 
        this->initialize_mesh_data();
    }
    // pass off to time filter
    else if (strcmp(arg[argIdx],"filter")==0) { 
      argIdx++;
      match = timeFilterManager_.modify(narg-argIdx,&arg[argIdx]);

        // consistentInitialization_ = false; 
    }
    // pass off to kernel function manager
    else if (strcmp(arg[argIdx],"kernel")==0) {
      argIdx++;

      if (kernelFunction_) {
        //delete kernelFunction_;
        //resetKernelFunction_ = true;
      }
      kernelFunction_ = KernelFunctionMgr::instance()->function(&arg[argIdx],narg-argIdx); 
      if (kernelFunction_) match = true;
      else ATC_Error("no matching kernel found");
      feEngine_->set_kernel(kernelFunction_);
      
      accumulantMol_=&kernelAccumulantMol_; // KKM add
      accumulantMolGrad_=&kernelAccumulantMolGrad_; // KKM add
    }
    // pass off to ghost manager
    else if (strcmp(arg[argIdx],"boundary_dynamics")==0) {
      argIdx++;
      match = ghostManager_.modify(narg-argIdx,&arg[argIdx]);
    }

    // parsing handled here
    else {
      if (strcmp(arg[argIdx],"parallel_consistency")==0) {
        argIdx++;
        //if (!kernelFunction_)          { throw ATC_Error("on_the_fly requires a kernel function"); }
        if (strcmp(arg[argIdx],"off")==0) parallelConsistency_ = false;
        else                              parallelConsistency_ = true;
        match = true;  
      }
     /*! \page man_hardy_on_the_fly fix_modify AtC on_the_fly 
        \section syntax
        fix_modify AtC on_the_fly <bond | kernel> <optional on | off> \n        - bond | kernel (keyword) = specifies on-the-fly calculation of bond or 
kernel         matrix elements \n
        - on | off (keyword) =  activate or discontinue on-the-fly mode \n 
        \section examples
        <TT> fix_modify AtC on_the_fly bond on </TT> \n        <TT> fix_modify AtC on_the_fly kernel </TT> \n
        <TT> fix_modify AtC on_the_fly kernel off </TT> \n
        \section description
        Overrides normal mode of pre-calculating and storing bond pair-to-node a
nd 
        kernel atom-to-node matrices. If activated, will calculate elements of t
hese
        matrices during repeated calls of field computations (i.e. "on-the-fly") and not store them for
        future use.   \n        on flag is optional - if omitted, on_the_fly will be activated for the s
pecified 
        matrix. Can be deactivated using off flag. \n 
        \section restrictions
        Must be used with the hardy/field type of AtC fix 
        ( see \ref man_fix_atc )
        \section related
        \section default
        By default, on-the-fly calculation is not active (i.e. off). However, code does a memory allocation check to determine if it can store all needed bond and kernel matrix ele ments. If this allocation fails, on-the-fly is activated. \n
      */

      else if (strcmp(arg[argIdx],"on_the_fly")==0) {
        argIdx++;
        //if (!kernelFunction_)          { throw ATC_Error("on_the_fly requires a kernel function"); }
        if (strcmp(arg[argIdx],"bond")==0) {
          argIdx++;
          bondOnTheFly_ = true;
          if (narg > argIdx && strcmp(arg[argIdx],"off")==0) bondOnTheFly_ = false;
        }
        else if (strcmp(arg[argIdx],"kernel")==0) {
          argIdx++;
          kernelOnTheFly_ = true;
          if (narg > argIdx && strcmp(arg[argIdx],"off")==0) kernelOnTheFly_ = false;
        }
        else { throw ATC_Error("unsupported on_the_fly type"); }
        match = true;
      }

      /*! \page man_output fix_modify AtC output
        \section syntax
        fix_modify AtC output <filename_prefix> <frequency>
        [text | full_text | binary | vector_components | tensor_components ]
        fix_modify AtC output index [step | time ]
        - filename_prefix (string) = prefix for data files 
        - frequency (integer) = frequency of output in time-steps 
        - options (keyword/s): \n
        text = creates text output of index, step and nodal variable values for unique nodes \n
        full_text = creates text output index, nodal id, step, nodal coordinates and nodal variable values for unique and image nodes \n
        binary = creates binary Ensight output \n
        vector_components = outputs vectors as scalar components \n
        tensor_components = outputs tensor as scalar components 
        (use this for Paraview)\n
         
        \section examples
        <TT> fix_modify AtC output heatFE 100 </TT> \n
        <TT> fix_modify AtC output hardyFE 1 text tensor_components </TT> \n
        <TT> fix_modify AtC output hardyFE 10 text binary tensor_components </TT> \n
        <TT> fix_modify AtC output index step </TT> \n
        \section description
        Creates text and/or binary (Ensight, "gold" format) output of nodal/mesh data 
        which is transfer/physics specific. Output indexed by step or time is possible.
        \section restrictions
        \section related 
        see \ref man_fix_atc
        \section default 
        no default format
        output indexed by time
      */
      else if (strcmp(arg[argIdx],"output")==0) {
        argIdx++;
      /*! \page man_output_nodeset fix_modify AtC output nodeset
        \section syntax
        fix_modify AtC output nodeset <nodeset_name> <operation>
        - nodeset_name (string) = name of nodeset to be operated on
        - operation (keyword/s): \n
        sum = creates nodal sum over nodes in specified nodeset \n
        \section examples
        <TT> fix_modify AtC output nodeset nset1 sum </TT> \n
        \section description
        Performs operation over the nodes belonging to specified nodeset
        and outputs resulting variable values to GLOBALS file.
        \section restrictions
        \section related 
        see \ref man_fix_atc
        \section default 
        none
      */
        if (strcmp(arg[argIdx],"nodeset")==0) {
          argIdx++;
          string nset = arg[argIdx++];
          if       (strcmp(arg[argIdx],"sum")==0) {
            argIdx++;
            string field = arg[argIdx];
            pair < string, FieldName >  id(nset,string_to_field(field));
            nsetData_[id] = NODESET_SUM;
            match = true;
          }
          else if (strcmp(arg[argIdx],"average")==0) {
            argIdx++;
            string field = arg[argIdx];
            pair < string, FieldName >  id(nset,string_to_field(field));
            nsetData_[id] = NODESET_AVERAGE;
            match = true;
          }
        }


      /*! \page man_boundary_integral fix_modify AtC output boundary_integral 
        \section syntax
        fix_modify AtC output boundary_integral [field] faceset [name]
        - field (string) : name of hardy field
        - name (string)  : name of faceset
        \section examples
        <TT> fix_modify AtC output boundary_integral stress faceset loop1 </TT> \n
        \section description
        Calculates a surface integral of the given field dotted with the
        outward normal of the faces and puts output in the "GLOBALS" file
        \section restrictions
        Must be used with the hardy/field type of AtC fix 
        ( see \ref man_fix_atc )
        \section related
        \section default
        none
      */

      /*! \page man_contour_integral fix_modify AtC output contour_integral
        \section syntax
        fix_modify AtC output contour_integral [field] faceset [name] <axis [x | y | z
]>
        - field (string) : name of hardy field
        - name (string)  : name of faceset
        - axis (string)  : x or y or z
        \section examples
        <TT> fix_modify AtC output contour_integral stress faceset loop1 </TT> \n
        \section description
        Calculates a surface integral of the given field dotted with the
        outward normal of the faces and puts output in the "GLOBALS" file
        \section restrictions
        Must be used with the hardy/field type of AtC fix 
        ( see \ref man_fix_atc )
        \section related
        \section default
        none
      */

        else if ( (strcmp(arg[argIdx],"boundary_integral")==0) 
               || (strcmp(arg[argIdx],"contour_integral")==0) ) {
          FacesetIntegralType type = BOUNDARY_INTEGRAL;
          if  (strcmp(arg[argIdx],"contour_integral")==0) 
                              type = CONTOUR_INTEGRAL;
          argIdx++;
          string field(arg[argIdx++]);
          if(strcmp(arg[argIdx],"faceset")==0) {
            argIdx++;
            string name(arg[argIdx++]);
            pair <string,string> pair_name(name,field);
            fsetData_[pair_name] = type;
            match = true;
          }
        } // end "boundary_integral" || "contour_integral"

      /*! \page man_output_elementset fix_modify AtC output elementset
        \section syntax
        fix_modify AtC output volume_integral <eset_name> <field> {`
        - set_name (string) = name of elementset to be integrated over
        - fieldname (string) = name of field to integrate
        csum = creates nodal sum over nodes in specified nodeset \n
        \section examples
        <TT> fix_modify AtC output eset1 mass_density </TT> \n
        \section description
        Performs volume integration of specified field over elementset
        and outputs resulting variable values to GLOBALS file.
        \section restrictions
        \section related 
        see \ref man_fix_atc
        \section default 
        none
      */

        else if ( (strcmp(arg[argIdx],"volume_integral")==0) ) {
          argIdx++;
          string name(arg[argIdx++]);
          string field(arg[argIdx++]);
          pair <string,FieldName> pair_name(name,string_to_field(field));
          if (++argIdx < narg) { // keyword average
            esetData_[pair_name] = ELEMENTSET_AVERAGE; 
          }
          else {
            esetData_[pair_name] = ELEMENTSET_TOTAL; 
          }
          match = true;
        }

        else if (strcmp(arg[argIdx],"now")==0) {
          argIdx++;
          double dt = 1.0;
          if (argIdx < narg) {
            dt = atof(arg[argIdx++]);
          }
          update_time(dt);
          update_step();
          outputNow_ = true;
          this->output();
          outputNow_ = false;
          match = true;
        }
        else 
          if (strcmp(arg[argIdx],"index")==0) {
            argIdx++;
            if (strcmp(arg[argIdx],"step")==0) { outputTime_ = false; }
            else                               { outputTime_ = true; }
          match = true;
        }
        else {
          outputPrefix_ = arg[argIdx++];
          outputFrequency_ = atoi(arg[argIdx++]);
          bool ensight_output = false, full_text_output = false;
          bool text_output = false, vect_comp = false, tensor_comp = false;
          int rank = lammpsInterface_->comm_rank();
          for (int i = argIdx; i<narg; ++i) {
            if      (strcmp(arg[i],"full_text")==0) full_text_output = true;
            else if (strcmp(arg[i],"text")==0)           text_output = true;
            else if (strcmp(arg[i],"binary")==0)      ensight_output = true;
            else if (strcmp(arg[i],"vector_components")==0) vect_comp = true;
            else if (strcmp(arg[i],"tensor_components")==0) tensor_comp = true;
            else { throw ATC_Error(" output: unknown keyword ");  }
          }
          if (outputFrequency_>0) {
            set<int> otypes;
            if (full_text_output || text_output) {
              lammpsInterface_->print_msg_once("Warning : text output can create _LARGE_ files");
            }
            if (full_text_output) otypes.insert(FULL_GNUPLOT);
            if (text_output)      otypes.insert(GNUPLOT);
            if (ensight_output)   otypes.insert(ENSIGHT);
            if (ntracked() > 0) {
               string fstem = field_to_string(SPECIES_CONCENTRATION);
               string istem = field_to_intrinsic_name(SPECIES_CONCENTRATION);
               vector<string> tnames = tracked_names();
               vector<string> fnames;
               vector<string> inames;
               for (unsigned int i = 0; i < tnames.size(); i++) {
                 fnames.push_back(fstem+tnames[i]);
                 inames.push_back(istem+tnames[i]);
               }
               feEngine_->add_field_names(fstem,fnames);
               feEngine_->add_field_names(istem,inames);
            }
            feEngine_->initialize_output(rank,outputPrefix_,otypes);
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
    // add a species for tracking
    /*! \page man_add_species fix_modify AtC add_species
      \section syntax
      fix_modify_AtC add_species <TAG> <group|type> <ID> \n
      - <TAG> = tag for tracking a species \n
      - group|type = LAMMPS defined group or type of atoms \n
      - <ID> = name of group or type number \n
      \section examples
      <TT> fix_modify AtC add_species gold type 1 </TT> \n
      <TT> group GOLDGROUP type 1 </TT> \n
      <TT> fix_modify AtC add_species gold group GOLDGROUP </TT>
      \section description
      Associates a tag with all atoms of a specified type or within a specified group. \n 
      \section restrictions
      \section related
      \section default
      No defaults for this command. 
    */
    else if (strcmp(arg[argIdx],"add_species")==0) {
      argIdx++;
      string speciesTag = arg[argIdx];
      string tag = arg[argIdx];
      argIdx++;
      if (strcmp(arg[argIdx],"group")==0) {
        if (narg-argIdx == 2) {
          string name = arg[++argIdx];
          int id = lammpsInterface_->group_bit(name);
          groupList_.push_back(id); 
          groupNames_.push_back(tag); 
        } 
        else {
          while (++argIdx < narg) {
            string name = arg[argIdx];
            int id = lammpsInterface_->group_bit(name);
            string tag = speciesTag+"-"+name;
            groupList_.push_back(id); 
            groupNames_.push_back(tag); 
          }
        }
      }
      else if (strcmp(arg[argIdx],"type")==0) {
        if (narg-argIdx == 2) {
          int id = atoi(arg[++argIdx]);
          typeList_.push_back(id); 
          typeNames_.push_back(tag); 
        } 
        else {
          while (++argIdx < narg) {
            int id = atoi(arg[argIdx]);
            string tag = speciesTag+"_"+to_string(id);
            typeList_.push_back(id); 
            typeNames_.push_back(tag); 
          }
        }
      }
      else {
        throw ATC_Error("ATC_Method: add_species only handles groups or types"); }
      match = true;
    }

    // remove species from tracking
    
    /*! \page man_remove_species fix_modify AtC remove_species
      \section syntax
      fix_modify_AtC delete_species <TAG> \n
      
      - <TAG> = tag for tracking a species \n
      \section examples
      <TT> fix_modify AtC remove_species gold </TT> \n
      \section description
      Removes tag designated for tracking a specified species. \n 
      \section restrictions
      \section related
      \section default
      No defaults for this command. 
    */
    else if (strcmp(arg[argIdx],"delete_species")==0) {
      argIdx++;
      string tag = arg[argIdx++];
      if (strcmp(arg[argIdx],"group")==0) {
        for (unsigned int j = 0; j < groupList_.size(); j++) {
          if (tag == groupNames_[j]) {
            groupList_.erase(groupList_.begin()+j); 
            groupNames_.erase(groupNames_.begin()+j); 
            break;
          }
        }
      }
      else if (strcmp(arg[argIdx],"type")==0) {
        for (unsigned int j = 0; j < typeList_.size(); j++) {
          if (tag == typeNames_[j]) {
            typeList_.erase(typeList_.begin()+j); 
            typeNames_.erase(typeNames_.begin()+j); 
            break;
          }
        }
      }
      else {
        throw ATC_Error("ATC_Method: delete_species only handles groups or types"); }
      match = true;
      
    }

    // add a molecule for tracking
    /*! \page man_add_molecule fix_modify AtC add_molecule
      \section syntax
      fix_modify_AtC add_molecule <small|large> <TAG> <GROUP_NAME> \n
      
      - small|large = can be small if molecule size < cutoff radius, must be large otherwise \n
      - <TAG> = tag for tracking a species \n
      - <GROUP_NAME> = name of group that tracking will be applied to \n   
      \section examples
      <TT> group WATERGROUP type 1 2 </TT> \n
      <TT> fix_modify AtC add_molecule small water WATERGROUP </TT> \n
      \section description
      Associates a tag with all molecules corresponding to a specified group. \n 
      \section restrictions
      \section related
      \section default
      No defaults for this command. 
    */
    else if (strcmp(arg[argIdx],"add_molecule")==0) {
      argIdx++;
      MolSize size;
      if (strcmp(arg[argIdx],"small")==0) {
        size = MOL_SMALL;
        //needXrefProcessorGhosts_ = true;
        needProcGhost_ = true;
      }
      else
        throw ATC_Error("ATC_CouplingMass:  Bad molecule size in add_molecule");
      argIdx++;
      string moleculeTag = arg[argIdx];
      
      argIdx++;
      int groupBit = lammpsInterface_->group_bit(arg[argIdx]);
      moleculeIds_[moleculeTag] = pair<MolSize,int>(size,groupBit);
      match = true;
    }

    // remove molecule from tracking
    /*! \page man_remove_molecule fix_modify AtC remove_molecule
      \section syntax
      fix_modify_AtC remove_molecule <TAG> \n
      
      - <TAG> = tag for tracking a molecule type \n
      \section examples
      <TT> fix_modify AtC remove_molecule water </TT> \n
      \section description
      Removes tag designated for tracking a specified set of molecules. \n 
      \section restrictions
      \section related
      \section default
      No defaults for this command. 
    */
    else if (strcmp(arg[argIdx],"remove_molecule")==0) {
      argIdx++;
      string moleculeTag = arg[argIdx];
      moleculeIds_.erase(moleculeTag);
      
      taggedDensMan_.erase(moleculeTag);
    }
          
      /*! \page man_boundary fix_modify AtC boundary 
        \section syntax
        fix_modify AtC boundary type <atom-type-id>
        - <atom-type-id> = type id for atoms that represent a ficticious
        boundary internal to the FE mesh
        \section examples
        <TT> fix_modify AtC boundary type ghost_atoms </TT>
        \section description
        Command to define the atoms that represent the ficticious 
        boundary internal to the FE mesh. For fully overlapped MD/FE 
        domains with periodic boundary conditions no boundary atoms should
        be defined.
        \section restrictions
        \section default 
        none
      */
      else if (strcmp(arg[argIdx],"boundary")==0) {
        argIdx++;
        groupTagGhost_ = arg[argIdx++];
        match = true;
      }

      /*! \page man_internal_atom_integrate fix_modify AtC internal_atom_integrate
        \section syntax
        fix_modify AtC internal_atom_integrate <on | off>
        <TT> fix_modify AtC internal_atom_integrate on </TT>
        \section description
        Has AtC perform time integration for the atoms in the group on which it operates.  This does not include boundary atoms.
        \section restrictions
        AtC must be created before any fixes doing time integration.  It must be on for coupling methods which impose constraints on velocities during the first verlet step, e.g. control momentum glc_velocity.
        \section default
        on for coupling methods, off for post-processors
        off
       */
      else if (strcmp(arg[argIdx],"internal_atom_integrate")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"off")==0) {
          integrateInternalAtoms_ = false;
          match = true;
        }
        else {
          integrateInternalAtoms_ = true;
          match = true;
        }
      }

      /*! \page man_internal_element_set fix_modify AtC internal_element_set
        \section syntax
        fix_modify AtC internal_element_set <element-set-name>
        - <element-set-name> = name of element set defining internal region, or off
        \section examples
        <TT> fix_modify AtC internal_element_set myElementSet </TT>
        <TT> fix_modify AtC internal_element_set off </TT>
        \section description
        Enables AtC to base the region for internal atoms to be an element set.
        If no ghost atoms are used, all the AtC atoms must be constrained to remain
        in this element set by the user, e.g., with walls.  If boundary atoms are 
        used in conjunction with Eulerian atom maps
        AtC will partition all atoms of a boundary or internal type to be of type internal
        if they are in the internal region or to be of type boundary otherwise.
        \section restrictions
        If boundary atoms are used in conjunction with Eulerian atom maps, the Eulerian
        reset frequency must be an integer multiple of the Lammps reneighbor frequency
        \section related
        see \ref atom_element_map_type and \ref boundary
        \section default
        off
       */
      else if (strcmp(arg[argIdx],"internal_element_set")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"off")==0) {
          internalElementSet_ = "";
          match = true;
        }
        else {
          internalElementSet_ = string(arg[argIdx]);
          const set<int> & elementSet((feEngine_->fe_mesh())->elementset(internalElementSet_)); // check it exists and is not trivial
          if (elementSet.size()==0) throw ATC_Error("internal_element_set - element set " + internalElementSet_ + " has no elements");
          match = true;
        }
      }

    /*! \page man_atom_weight fix_modify AtC atom_weight 
      \section syntax
      fix_modify AtC atom_weight <method> <arguments>
        - <method> = \n
          value: atoms in specified group assigned constant value given \n 
          lattice: volume per atom for specified lattice type (e.g. fcc) and parameter \n
          element: element volume divided among atoms within element \n
          region: volume per atom determined based on the atom count in the MD regions and their volumes. Note: meaningful only if atoms completely fill all the regions. \n
          group: volume per atom determined based on the atom count in a group and its volume\n
          read_in: list of values for atoms are read-in from specified file \n
      \section examples
       <TT> fix_modify atc atom_weight constant myatoms 11.8 </TT> \n
       <TT> fix_modify atc atom_weight lattice </TT> \n
       <TT> fix_modify atc atom_weight read-in atm_wt_file.txt </TT> \n
      \section description
       Command for assigning the value of atomic weights used for atomic integration in
       atom-continuum coupled simulations. 
      \section restrictions 
      Use of lattice option requires a lattice type and parameter is already specified.
      \section related
      \section default
      lattice
    */
      else if (strcmp(arg[argIdx],"atom_weight")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"constant")==0) { 
          argIdx++;
          atomWeightType_ = USER;
          int groupbit = -1;
          if (strcmp(arg[argIdx],"all")==0) {
          }
          else {
            groupbit = lammpsInterface_->group_bit(arg[argIdx]);
          }
          argIdx++;
          double value = atof(arg[argIdx]);
          Valpha_[groupbit] = value;
          match = true;
        }
        else if (strcmp(arg[argIdx],"lattice")==0) {
          atomWeightType_ = LATTICE;
          match = true;
        }
        else if (strcmp(arg[argIdx],"element")==0) {
          atomWeightType_ = ELEMENT;
          match = true;
        }
        else if (strcmp(arg[argIdx],"region")==0) {
          atomWeightType_ = REGION;
          match = true;
        }
        else if (strcmp(arg[argIdx],"group")==0) {
          atomWeightType_ = GROUP;
          match = true;
        }
        else if (strcmp(arg[argIdx],"multiscale")==0) {
          atomWeightType_ = MULTISCALE;
          match = true;
        }
        else if (strcmp(arg[argIdx],"node")==0) {
          atomWeightType_ = NODE;
          match = true;
        }
        else if (strcmp(arg[argIdx],"node_element")==0) {
          atomWeightType_ = NODE_ELEMENT;
          match = true;
        }
        else if (strcmp(arg[argIdx],"read_in")==0) {
          atomWeightType_ = READ_IN;
          argIdx++;
          atomicWeightsFile_ = arg[argIdx];
          match = true;
        }
        if (match) {
          needReset_ = true;
        }
      }

    /*! \page man_decomposition fix_modify AtC decomposition 
      \section syntax
      fix_modify AtC decomposition <type> 
        - <type> = \n
          replicated_memory: nodal information replicated on each processor \n 
          distributed_memory: only owned nodal information on processor  \n
      \section examples
       <TT> fix_modify atc decomposition distributed_memory </TT> \n
      \section description
       Command for assigning the distribution of work and memory for parallel runs.
      \section restrictions 
      replicated_memory is appropriate for simulations were the number of nodes << number of atoms
      \section related
      \section default
      replicated_memory
    */
      else if (strcmp(arg[argIdx],"decomposition")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"replicated_memory")==0) { 
          domainDecomposition_ = REPLICATED_MEMORY;
          match = true;
        }
        else if (strcmp(arg[argIdx],"distributed_memory")==0) {
          domainDecomposition_ = DISTRIBUTED_MEMORY;
          match = true;
        }
      }

    /*! \page man_write_atom_weights fix_modify AtC write_atom_weights
      \section syntax
      fix_modify AtC write_atom_weights <filename> <frequency>
        - <filename> = name of file that atomic weights are written to \n
        - <frequency> = how often writes will occur \n
      \section examples
       <TT> fix_modify atc write_atom_weights atm_wt_file.txt 10 </TT> \n
      \section description
       Command for writing the values of atomic weights to a specified file.
      \section restrictions 
      \section related
      \section default
    */
      else if (strcmp(arg[argIdx],"write_atom_weights")==0) {
        argIdx++;
        atomicWeightsFile_ = arg[argIdx];
        argIdx++;
        atomicWeightsWriteFrequency_ = atoi(arg[argIdx]);
        atomicWeightsWriteFlag_ = true;
        match = true;
      }


      /*! \page man_reset_time fix_modify AtC reset_time 
      \section syntax
      fix_modify AtC reset_time <value> 
      \section examples
       <TT> fix_modify atc reset_time 0.0 </TT> \n
      \section description
      Resets the simulation time counter. 
      \section restrictions
      \section related
      \section default
      */
      else if (strcmp(arg[argIdx],"reset_time")==0) {
        argIdx++;
        set_time();
        if (narg > argIdx) {
          double time = atof(arg[argIdx]);
          set_time(time);
        }
        match = true;
      }

      /*! \page man_reset_time fix_modify AtC kernel_bandwidth
      \section syntax
      fix_modify AtC kernel_bandwidth <value> 
      \section examples
       <TT> fix_modify atc reset_time 8 </TT> \n
      \section description
      Sets a maximum parallel bandwidth for the kernel functions during parallel communication.  If the command is not issued, the default will be to assume the bandwidth of the kernel matrix corresponds to the number of sampling locations.
      \section restrictions
      Only is used if kernel functions are being used.
      \section related
      \section default
      Number of sample locations.
      */
      else if (strcmp(arg[argIdx],"kernel_bandwidth")==0) {
        argIdx++;
        accumulantBandwidth_ = atoi(arg[argIdx]);
        match = true;
      }

      /*! \page man_reset_atomic_reference_positions fix_modify AtC reset_atomic_reference_positions 
      \section syntax
      fix_modify AtC reset_atomic_reference_positions
      \section examples
       <TT> fix_modify atc reset_atomic_reference_positions
      \section description
      Resets the atomic positions ATC uses to perform point to field operations.
      In can be used to use perfect lattice sites in ATC but a thermalized or
      deformed lattice in LAMMPS.
      \section restrictions
      \section related
      \section default
      Default is off
      */
      else if (strcmp(arg[argIdx],"reset_atomic_reference_positions")==0) {
        argIdx++;
        xRefFile_ = arg[argIdx];
        readXref_ = true;
        match = true;
      }

      /*! \page man_set fix_modify AtC set 
        \section syntax
        fix_modify AtC set reference_potential_energy <value_or_filename(optional)>
        - value (double) : optional user specified zero point for PE in native LAMMPS energy units \n
        - filename (string) : optional user specified string for file of nodal PE values to be read-in 
        \section examples
        <TT> fix_modify AtC set reference_potential_energy </TT> \n
        <TT> fix_modify AtC set reference_potential_energy -0.05 </TT> \n
        <TT> fix_modify AtC set reference_potential_energy myPEvalues </TT> \n
        \section description
        Used to set various quantities for the post-processing algorithms.
        It sets the zero point for the potential energy density using 
        the value provided for all nodes, or from the current 
        configuration of the lattice if no value is provided, or 
        values provided within the specified filename.
        \section restrictions
        Must be used with the hardy/field type of AtC fix 
        ( see \ref man_fix_atc )
        \section related
        \section default
        Defaults to lammps zero point i.e. isolated atoms
      */
      else if (strcmp(arg[argIdx],"set")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"reference_potential_energy")==0) {
          argIdx++;
          setRefPE_ = true;
          if (narg > argIdx) {
            string a(arg[argIdx]);
            if (is_dbl(a)) {
              double value = atof(arg[argIdx]);
              refPEvalue_ = value;
              setRefPEvalue_ = true;
            }
            else {
              nodalRefPEfile_ = arg[argIdx];
              readRefPE_ = true;
            }
          } 
          match = true;
        }
      } // end "set"



      /*! \page man_atom_element_map fix_modify AtC atom_element_map
      \section syntax
      fix_modify AtC atom_element_map  <eulerian|lagrangian> <frequency> \n
      - frequency (int) : frequency of updating atom-to-continuum maps based on the
      current configuration - only for eulerian 
      \section examples
      <TT> fix_modify atc atom_element_map eulerian 100 </TT>
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
      else if (strcmp(arg[argIdx],"atom_element_map")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"eulerian")==0) {
          atomToElementMapType_ = EULERIAN;
          argIdx++;
          atomToElementMapFrequency_ = atoi(arg[argIdx]);
        } 
        else {
          atomToElementMapType_ = LAGRANGIAN;
          atomToElementMapFrequency_ = 0;
        } 
        match = true;
        needReset_ = true;
      }

      /*! \page man_read_restart fix_modify AtC read_restart
      \section syntax
      fix_modify AtC read_restart [file_name]  \n
      \section examples
      <TT> fix_modify AtC read_restart ATC_state </TT> \n
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
      else if (strcmp(arg[argIdx],"read_restart")==0) {
        argIdx++;
        restartFileName_  = arg[argIdx];
        useRestart_ = true;
        match = true;
      }

      /*! \page man_write_restart fix_modify AtC write_restart
      \section syntax
      fix_modify AtC write_restart [file_name]  \n
      \section examples
      <TT> fix_modify AtC write_restart restart.mydata </TT> \n
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
      else if (strcmp(arg[argIdx],"write_restart")==0) {
        argIdx++;
        string restartFileName(arg[argIdx]);
        RESTART_LIST data;
        write_restart_data(restartFileName,data);
        match = true;
      }

    } // end else 

    return match; // return to FixATC
  
  }

  //--------------------------------------------------
  // helper function for parser
  // handles : "displacement x"  and "temperature" by indexing argIdx
  // for fluxes : only normal fluxes can be prescribed
  //--------------------------------------------------
  void ATC_Method::parse_field(/*const*/ char ** args, int & argIdx,
                           FieldName & thisField, int & thisIndex)
  {
    string thisName = args[argIdx++];
    if (args[argIdx] == NULL) {
      throw ATC_Error("Need to give field '"+thisName+"' more args");
    }
    thisField = string_to_field(thisName);
    map<FieldName,int>::const_iterator iter = fieldSizes_.find(thisField);
    if (iter == fieldSizes_.end()) {
      throw ATC_Error("Bad field name: "+thisName);
    }
    string thisDim = args[argIdx];
    thisIndex = 0;
    if (string_to_index(thisDim,thisIndex)) {
      if ( !( thisIndex < fieldSizes_[thisField]) ) {
        throw ATC_Error("Bad field dimension "+thisDim);
      }
      argIdx++;
    }
  }

  //-------------------------------------------------------------------
  // this sets the peratom output
  
  void ATC_Method::update_peratom_output()
  {
    perAtomArray_ = perAtomOutput_; 
    // copy values
    for (int i = 0; i < lammpsInterface_->nlocal(); i++) {
      for (int j = 0; j < nsd_; j++) {
        perAtomOutput_[i][j] = xref_[i][j];
      }
      for (int j = nsd_; j < sizePerAtomCols_; j++) {
        perAtomOutput_[i][j] = 0;
      }
    }
    int indx = nsd_;
    if (atomVolume_->nRows() > 0) { // kernel Hardy does not compute these
      const DIAG_MAT & myAtomicWeights(atomVolume_->quantity());
      for (int i = 0; i < nLocal_; i++) {
        double wg = myAtomicWeights(i,i);
        if (wg > 0) {
          int ii = internalToAtom_(i);
          perAtomOutput_[ii][indx] = 1./wg;
        }
      }
    }
  }

  void ATC_Method::adjust_xref_pbc()
  {
    
    int nlocal = lammpsInterface_->nlocal();
    int xperiodic = lammpsInterface_->xperiodic();
    int yperiodic = lammpsInterface_->yperiodic();
    int zperiodic = lammpsInterface_->zperiodic();
    double **x = lammpsInterface_->xatom();
    double boxxlo,boxxhi;
    double boxylo,boxyhi;
    double boxzlo,boxzhi;

    lammpsInterface_->box_bounds(boxxlo,boxxhi,
                                 boxylo,boxyhi,
                                 boxzlo,boxzhi);
//  bool changed = false;
    for (int i = 0; i < nlocal; i++) {
      if (xperiodic) {
        if (x[i][0] < boxxlo) {
          xref_[i][0] += Xprd_;
//        changed = true;
        }
        if (x[i][0] >= boxxhi) { 
          xref_[i][0] -= Xprd_;
//        changed = true;
        } 
      } 

      if (yperiodic) {
        if (x[i][1] < boxylo) {
          xref_[i][1] += Yprd_;
//        changed = true;
        }
        if (x[i][1] >= boxyhi) {
          xref_[i][1] -= Yprd_;
//        changed = true;
        } 
      } 

      if (zperiodic) {
        if (x[i][2] < boxzlo) {
          xref_[i][2] += Zprd_;
//        changed = true;
        }
        if (x[i][2] >= boxzhi) {
          xref_[i][2] -= Zprd_;
//        changed = true;
        } 
      } 
    }

    // propagate reset if needed
    if (atomToElementMapType_ == LAGRANGIAN) {
      if (atomCoarseGrainingPositions_) {
        atomCoarseGrainingPositions_->force_reset();
      }
    }
    else if (atomReferencePositions_) {
      atomReferencePositions_->force_reset();
    }

  }
  //-------------------------------------------------------------------
  void ATC_Method::initialize()
  { 
    feEngine_->partition_mesh();
    // initialized_ is set to true by derived class initialize()
    // localStep_ is a counter within a run or minimize
    localStep_ = 0;
    // STEP 1)  get basic information data from Lammps/fix
    // 1a)  group ids for all internal atoms
    groupbit_ = lammpsInterface_->group_bit(groupTag_);

    // 1b) group ids for ghost atoms
    groupbitGhost_ = 0;
    if (!groupTagGhost_.empty()) {
      groupbitGhost_ = lammpsInterface_->group_bit(groupTagGhost_);
    }

    // 1c) periodicity and box bounds/lengths
    if (!initialized_) {
      
      lammpsInterface_->box_periodicity(periodicity[0],
                                            periodicity[1],
                                            periodicity[2]);
      lammpsInterface_->box_bounds(box_bounds[0][0],box_bounds[1][0],
                                       box_bounds[0][1],box_bounds[1][1],
                                       box_bounds[0][2],box_bounds[1][2]);
      for (int k = 0; k < nsd_; k++) {
        box_length[k] = box_bounds[1][k] - box_bounds[0][k]; 
      }

      lammpsInterface_->set_reference_box();

      // get periodicity data from lammps for parallel exchange to adjust for periodicity
      Xprd_ = lammpsInterface_->domain_xprd();
      Yprd_ = lammpsInterface_->domain_yprd();
      Zprd_ = lammpsInterface_->domain_zprd();
//    box_length[0] = Xprd_;
//    box_length[1] = Yprd_;
//    box_length[2] = Zprd_;
      XY_ = lammpsInterface_->domain_xy();
      XZ_ = lammpsInterface_->domain_xz();
      YZ_ = lammpsInterface_->domain_yz();
    }

    // STEP 2 computational geometry
    // 2a) get basic information from continuum/FE
    this->set_continuum_data();

    // STEP 2b) set up data structures for computational geometry
    if (this->reset_methods()) {
      // clear memory manager
      interscaleManager_.clear_temporary_data();
      atomVolume_ = NULL;

      // reference positions and energy
      if (!initialized_) {
        double **x = lammpsInterface_->xatom();
        for (int i = 0; i < lammpsInterface_->nmax(); i++) {
          for (int j = 0; j < nsd_; j++) {
            xref_[i][j] = x[i][j];
          }
        }
        
        // re-write non-ghosts' xref with values from a file
        if (readXref_) {
          bool success = read_atomic_ref_positions(xRefFile_.c_str());
          if (!success)
            throw ATC_Error("Error reading atomic reference positions");
          readXref_ = false;
        }

        // set up maps from lammps to atc indexing
        reset_nlocal();
      }

      this->set_computational_geometry();
    }

    // 2c) basic data regarding atomic system, e.g. atom coordinates
    if (atomToElementMapType_ == EULERIAN) {
      reset_coordinates();
    }

    // STEP 3) set up variables which will be integrated in time
    this->construct_time_integration_data();

    // STEP 4) instantiate all the various specific algorithms and methods
    this->construct_methods();

    // STEP 5) construct dependency-managed data
    // 5b) all other transfer operators
    // needs to be done before every run in case options have changed or the atoms have been changed by the user
    
    if (this->reset_methods()) {
      // construct all the needed data structures
      this->construct_transfers();

      // allocate all space needed for lammps arrays
      interscaleManager_.grow_arrays(lammpsInterface_->nmax());
    }
    // reset all computes invoked flags and lammps data
    interscaleManager_.lammps_force_reset();
 
    // STEP 6) initialize data
    // 6b) size quantities which use pack_comm
    interscaleManager_.size_comm_quantities();

    // 6c) set coarse-graining functions and atomic weights
    if (!initialized_) {
      // FE_Engine allocates all required memory
      // assume initial atomic position is the reference position for now
      
      // \int_\Omega N_I dV : static if the mesh is
      NodeVolumes_.reset(nNodes_,nNodes_);
      invNodeVolumes_.reset(nNodes_,nNodes_);
      feEngine_->compute_lumped_mass_matrix(NodeVolumes_);
      invNodeVolumes_.set_quantity() = NodeVolumes_.inv(); 
    }
    atomVolume_->set_reset();

    // 6d) reference values
    this->set_reference_potential_energy();

    // 6e) atomic output for 0th step
    update_peratom_output();

    // clear need for resets
    needReset_ = false;

  }
  //-------------------------------------------------------------------
  void ATC_Method::set_continuum_data()
  {
    // initialize finite element engine and get basic properties
    if (!initialized_) {
      feEngine_->initialize();
      if (nsd_!=feEngine_->nsd()) {
        throw ATC_Error("Spatial dimensions inconsistent between LAMMPS and ATC");
      }
      nNodes_ = feEngine_->num_nodes();
    }
  }

  //-------------------------------------------------------------------
  void ATC_Method::set_computational_geometry()
  {
    // set positions used for coarse-graining operators
    
    
    
    
    if (!initialized_) {
      if (atomToElementMapType_ == EULERIAN) {
        FundamentalAtomQuantity * atomPositionsAll = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_POSITION,ALL);
        ClonedAtomQuantity<double> * myAtomPositions =
          new ClonedAtomQuantity<double>(this,atomPositionsAll,INTERNAL);
        atomCoarseGrainingPositions_ = myAtomPositions;
        interscaleManager_.add_per_atom_quantity(myAtomPositions,
                                                 "AtomicCoarseGrainingPositions");

        if (trackDisplacement_) {
          XrefWrapper * myAtomReferencePositions = new XrefWrapper(this);
          atomReferencePositions_ = myAtomReferencePositions;
          interscaleManager_.add_per_atom_quantity(myAtomReferencePositions,
                                                   "AtomicReferencePositions");
          atomReferencePositions_->set_memory_type(PERSISTENT);
        }

        if (groupbitGhost_) {
          myAtomPositions = new ClonedAtomQuantity<double>(this,atomPositionsAll,GHOST);
          atomGhostCoarseGrainingPositions_ = myAtomPositions;
          interscaleManager_.add_per_atom_quantity(myAtomPositions,
                                                   "AtomicGhostCoarseGrainingPositions");
        }
        if(needProcGhost_){
          FundamentalAtomQuantity * atomPositionsAll = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_POSITION,PROC_GHOST);
          ClonedAtomQuantity<double> * myAtomPositions =
            new ClonedAtomQuantity<double>(this,atomPositionsAll,PROC_GHOST);
          atomProcGhostCoarseGrainingPositions_ = myAtomPositions;
          interscaleManager_.add_per_atom_quantity(myAtomPositions,
                                                   "AtomicProcGhostCoarseGrainingPositions");
        }  
      }
      else {
        XrefWrapper * myAtomPositions = new XrefWrapper(this);
        atomCoarseGrainingPositions_ = myAtomPositions;
        interscaleManager_.add_per_atom_quantity(myAtomPositions,
                                                 "AtomicCoarseGrainingPositions");
        atomReferencePositions_ = atomCoarseGrainingPositions_;

        if (groupbitGhost_) {
          myAtomPositions = new XrefWrapper(this,GHOST);
          atomGhostCoarseGrainingPositions_ = myAtomPositions;
          interscaleManager_.add_per_atom_quantity(myAtomPositions,
                                                   "AtomicGhostCoarseGrainingPositions");
        }
        if (needProcGhost_) {
          XrefWrapper * myAtomPositions = new XrefWrapper(this);
          atomProcGhostCoarseGrainingPositions_ = myAtomPositions;
          interscaleManager_.add_per_atom_quantity(myAtomPositions,
                                                   "AtomicProcGhostCoarseGrainingPositions");
        }
      }
      atomCoarseGrainingPositions_->set_memory_type(PERSISTENT);
      if (atomGhostCoarseGrainingPositions_) atomGhostCoarseGrainingPositions_->set_memory_type(PERSISTENT);
      if (atomProcGhostCoarseGrainingPositions_) atomProcGhostCoarseGrainingPositions_->set_memory_type(PERSISTENT);
    }

    // Add in atom to element map if using shape functions
    if (needsAtomToElementMap_) {
      atomElement_ = new AtomToElementMap(this);
      interscaleManager_.add_per_atom_int_quantity(atomElement_,"AtomElement");
    }
  }

  //-------------------------------------------------------------------
  void ATC_Method::construct_methods()
  {
    
    if (this->reset_methods()) {
      if (atomTimeIntegrator_) delete atomTimeIntegrator_;
      if (integrateInternalAtoms_) {
        atomTimeIntegrator_ = new AtomTimeIntegratorType(this,INTERNAL);
      }
      else {
        atomTimeIntegrator_ = new AtomTimeIntegrator();
      }

      // set up integration schemes for ghosts
      ghostManager_.construct_methods();
    }
  }

  //-------------------------------------------------------------------
  void ATC_Method::construct_transfers()
  {
    this->construct_interpolant();

    this->construct_molecule_transfers();

    atomTimeIntegrator_->construct_transfers();
    ghostManager_.construct_transfers();



  }
  //-------------------------------------------------------------------
  PerAtomDiagonalMatrix<double> * ATC_Method::create_atom_volume()
  {
    if (atomVolume_) {
      return atomVolume_;
    }
    else {
      // set variables to compute atomic weights
      DENS_MAN * nodalVolume(NULL);
      switch (atomWeightType_) {
      case USER:
        atomVolume_ = new AtomVolumeUser(this,Valpha_);
        break;
      case LATTICE:
        atomVolume_ = new AtomVolumeLattice(this);
        break;
      case ELEMENT:
        atomVolume_ = new AtomVolumeElement(this);
        break;
      case REGION:
        atomVolume_ = new AtomVolumeRegion(this);
        break;
      case GROUP:
        atomVolume_ = new AtomVolumeGroup(this,Valpha_);
        break;
      case MULTISCALE:
        if (!shpFcn_) {
          throw ATC_Error("ATC_Method::create_atom_volume - Multiscale algorithm requires an interpolant");
        }
        nodalVolume = new NodalAtomVolume(this,shpFcn_);
        interscaleManager_.add_dense_matrix(nodalVolume,"NodalAtomVolume");
        atomVolume_ = new FtaShapeFunctionProlongationDiagonalMatrix(this,nodalVolume,shpFcn_);
        break;
      case NODE:
        if (!shpFcn_) {
          throw ATC_Error("ATC_Method::create_atom_volume - Node algorithm requires an interpolant");
        }
        nodalVolume = new NodalVolume(this,shpFcn_);
        interscaleManager_.add_dense_matrix(nodalVolume,"NodalVolume");
        atomVolume_ = new FtaShapeFunctionProlongationDiagonalMatrix(this,nodalVolume,shpFcn_);
        break;
      case NODE_ELEMENT:
        if (!shpFcn_) {
          throw ATC_Error("ATC_Method::create_atom_volume - Node-Element algorithm requires an interpolant");
        }
        nodalVolume = new NodalAtomVolumeElement(this,shpFcn_);
        interscaleManager_.add_dense_matrix(nodalVolume,"NodalAtomVolumeElement");
        atomVolume_ = new FtaShapeFunctionProlongationDiagonalMatrix(this,nodalVolume,shpFcn_);
        break;
      case READ_IN:
        atomVolume_ = new AtomVolumeFile(this,atomicWeightsFile_);
        break;
      }
      if (atomVolume_) {
        interscaleManager_.add_per_atom_diagonal_matrix(atomVolume_,"AtomVolume");
      }
      else {
        throw ATC_Error("ATC_Method::create_atom_volume - bad option for atom volume algorithm");
      }

      return atomVolume_;
    }
  }
  //--------------------------------------------------------
  void ATC_Method::init_integrate()
  {
    atomTimeIntegrator_->init_integrate_velocity(dt());
    ghostManager_.init_integrate_velocity(dt());
    // account for other fixes doing time integration
    interscaleManager_.fundamental_force_reset(LammpsInterface::ATOM_VELOCITY);

    atomTimeIntegrator_->init_integrate_position(dt());
    ghostManager_.init_integrate_position(dt());
    // account for other fixes doing time integration
    interscaleManager_.fundamental_force_reset(LammpsInterface::ATOM_POSITION);
  }
  //-------------------------------------------------------------------
  void ATC_Method::post_init_integrate()
  {
    ghostManager_.post_init_integrate();
  }
  //-------------------------------------------------------------------
  void ATC_Method::pre_exchange()
  {
    adjust_xref_pbc();
    // call interscale manager to sync atc per-atom data with lammps array ahead of parallel communication
    interscaleManager_.prepare_exchange();

    // change types based on moving from internal region to ghost region
    if ((atomToElementMapType_ == EULERIAN) && (step() % atomToElementMapFrequency_ == 0)) {
      ghostManager_.pre_exchange();
    }
  }
  //-------------------------------------------------------------------
  void ATC_Method::setup_pre_exchange()
  {
    adjust_xref_pbc();
    // call interscale manager to sync atc per-atom data with lammps array ahead of parallel communication
    interscaleManager_.prepare_exchange();
  }
  //-------------------------------------------------------------------
  void ATC_Method::pre_neighbor()
  {
    // reset quantities arising from atom exchange
    reset_nlocal();

    interscaleManager_.post_exchange();

    // forward_comm should go here
  }
  //-------------------------------------------------------------------
  void ATC_Method::min_post_force()
  {
    post_force();
  }
  //-------------------------------------------------------------------
  void ATC_Method::post_force()
  {
    // this resets allow for the possibility of other fixes modifying positions and velocities, e.g. walls, but reduces efficiency
    interscaleManager_.lammps_force_reset();
  }
  //--------------------------------------------------------
  void ATC_Method::final_integrate()
  {
    atomTimeIntegrator_->final_integrate(dt());
    ghostManager_.final_integrate(dt());
    // account for other fixes doing time integration
    interscaleManager_.fundamental_force_reset(LammpsInterface::ATOM_VELOCITY);
  }
  //-------------------------------------------------------------------
  void ATC_Method::post_final_integrate()
  {
    if (atomicWeightsWriteFlag_ && (step() % atomicWeightsWriteFrequency_ == 0)) {
      write_atomic_weights(atomicWeightsFile_,atomVolume_->quantity());
    }
  }
  //-------------------------------------------------------------------
  void ATC_Method::end_of_step()
  {
    localStep_ += 1;
  }
  //--------------------------------------------------------------
  void ATC_Method::finish()
  {
    // FE Engine
    if (feEngine_) feEngine_->finish();
    feEngine_->departition_mesh();
  }

  //--------------------------------------------------------------
  /** method to add new fields to the included list */
  //--------------------------------------------------------------
  void ATC_Method::add_fields(map<FieldName,int> & newFieldSizes)
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
 
//-------------------------------------------------------------------
  void ATC_Method::set_reference_potential_energy(void)
  {
    if (setRefPE_) {
      if (setRefPEvalue_) {
        nodalRefPotentialEnergy_->set_quantity() = refPEvalue_;
        setRefPEvalue_ = false;
      }
      else if (readRefPE_) {
        if (LammpsInterface::instance()->rank_zero()) {
          stringstream ss;
          ss << "reading reference potential energy from " << nodalRefPEfile_;
          LammpsInterface::instance()->print_msg(ss.str());
        }
        (nodalRefPotentialEnergy_->set_quantity()).from_file(nodalRefPEfile_);
        readRefPE_ = false;
      }
      else {
        hasRefPE_ = false;
        SPAR_MAN * referenceAccumulant = interscaleManager_.sparse_matrix("ReferenceAccumulant");
        if (referenceAccumulant) {
          referenceAccumulant->set_quantity() = accumulant_->quantity();
        }
        DIAG_MAN * referenceAccumulantInverseVolumes = interscaleManager_.diagonal_matrix("ReferenceAccumulantInverseVolumes");
        if (referenceAccumulantInverseVolumes) {
          referenceAccumulantInverseVolumes->set_quantity() = accumulantInverseVolumes_->quantity();
        }
        PAQ * atomicRefPe = interscaleManager_.per_atom_quantity("AtomicReferencePotential");
        if (!atomicRefPe) {
          throw ATC_Error("ATC_Method::set_reference_potential_energy - atomic reference PE object was not created during construct_transfers");
        }
        PAQ* pe = interscaleManager_.per_atom_quantity("AtomicPotentialEnergy");
        if (!pe) {
          throw ATC_Error("ATC_Method::set_reference_potential_energy - atomic PE object was not created during construct_transfers");
        }
        atomicRefPe->set_quantity() = pe->quantity();
        atomicRefPe->fix_quantity();
      }
      setRefPE_ = false;
      hasRefPE_ = true;
    }
  }
//-------------------------------------------------------------------
 

  //=================================================================
  // memory management and processor information exchange
  //=================================================================


  //-----------------------------------------------------------------
  // number of doubles 
  //-----------------------------------------------------------------
  int ATC_Method::doubles_per_atom() const
  {
    
    int doubles = 4;
    doubles += interscaleManager_.memory_usage();
    return doubles;
  }

  //-----------------------------------------------------------------
  // memory usage of local atom-based arrays 
  //-----------------------------------------------------------------
  int ATC_Method::memory_usage()
  {
    int bytes = doubles_per_atom();
    bytes *= lammpsInterface_->nmax() * sizeof(double);
    return bytes;
  }

  //-----------------------------------------------------------------
  // allocate local atom-based arrays 
  //-----------------------------------------------------------------
  void ATC_Method::grow_arrays(int nmax)
  {
    xref_ =
      lammpsInterface_->grow_2d_double_array(xref_,nmax,3,"fix_atc:xref");

    perAtomOutput_ = 
      lammpsInterface_->grow_2d_double_array(perAtomOutput_,nmax,sizePerAtomCols_,"fix_atc:perAtomOutput");
    interscaleManager_.grow_arrays(nmax);
  }

  //-----------------------------------------------------------------
  // copy values within local atom-based arrays 
  //-----------------------------------------------------------------
  void ATC_Method::copy_arrays(int i, int j)
  {
    xref_[j][0] = xref_[i][0];
    xref_[j][1] = xref_[i][1];
    xref_[j][2] = xref_[i][2];

    for (int ii = 0 ; ii < sizePerAtomCols_ ; ii++ ) {
      perAtomOutput_[j][ii] = perAtomOutput_[i][ii];
    }
    interscaleManager_.copy_arrays(i,j);
  }

  //-----------------------------------------------------------------
  // pack values in local atom-based arrays for exchange with another proc 
  //-----------------------------------------------------------------
  int ATC_Method::pack_exchange(int i, double *buf)
  {
    buf[0] = xref_[i][0]; 
    buf[1] = xref_[i][1]; 
    buf[2] = xref_[i][2]; 
    
    int j = 4;
    for (int ii = 0 ; ii < sizePerAtomCols_ ; ii++ ) {
      buf[j++] = perAtomOutput_[i][ii];
    }
    int interscaleSizeComm = interscaleManager_.pack_exchange(i,&buf[j]);
    return sizeComm_ + interscaleSizeComm;
  }

  //-----------------------------------------------------------------
  // unpack values in local atom-based arrays from exchange with another proc 
  //-----------------------------------------------------------------
  int ATC_Method::unpack_exchange(int nlocal, double *buf)
  {
    xref_[nlocal][0] = buf[0];
    xref_[nlocal][1] = buf[1];
    xref_[nlocal][2] = buf[2];

    int j = 4;
    for (int ii = 0 ; ii < sizePerAtomCols_ ; ii++ ) {
      perAtomOutput_[nlocal][ii] = buf[j++];
    }
    int interscaleSizeComm = interscaleManager_.unpack_exchange(nlocal,&buf[j]);
    return sizeComm_ + interscaleSizeComm;
  }

  //-----------------------------------------------------------------
  // pack values in local atom-based arrays from exchange with another proc 
  //-----------------------------------------------------------------
  int ATC_Method::pack_comm(int n, int *list, double *buf, 
                            int pbc_flag, int *pbc)
  {
    int i,j,m;
    double dx = 0,dy = 0,dz = 0;

    int * num_bond = lammpsInterface_->num_bond();
    int ** bond_atom = lammpsInterface_->bond_atom();
  
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = xref_[j][0];
        buf[m++] = xref_[j][1];
        buf[m++] = xref_[j][2];
        
        if (num_bond) {
          buf[m++] = num_bond[j];
          for (int ii = 0; ii < lammpsInterface_->bond_per_atom(); ii++) {
            buf[m++] = bond_atom[j][ii];
          }
        }
      }
    } 
    else {
      if (lammpsInterface_->domain_triclinic() == 0) {
        dx = pbc[0]*Xprd_;
        dy = pbc[1]*Yprd_;
        dz = pbc[2]*Zprd_;
      } 
      else {
        dx = pbc[0]*Xprd_ + pbc[5]*XY_ + pbc[4]*XZ_;
        dy = pbc[1]*Yprd_ + pbc[3]*YZ_;
        dz = pbc[2]*Zprd_;
      }
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = xref_[j][0] + dx;
        buf[m++] = xref_[j][1] + dy;
        buf[m++] = xref_[j][2] + dz;

        if (num_bond) {
          buf[m++] = num_bond[j];
          for (int ii = 0; ii < lammpsInterface_->bond_per_atom(); ii++) {
            buf[m++] = bond_atom[j][ii];
          }
        }

      }
    }

    int mySize = 3;
    if (num_bond)
      mySize += 1 + lammpsInterface_->bond_per_atom();
    return mySize;
  }

  //-----------------------------------------------------------------
  // unpack values in local atom-based arrays from exchange with another proc 
  //-----------------------------------------------------------------
  void ATC_Method::unpack_comm(int n, int first, double *buf)
  {
    int i,m,last;

    int * num_bond = lammpsInterface_->num_bond();
    int ** bond_atom = lammpsInterface_->bond_atom();

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      xref_[i][0] = buf[m++];
      xref_[i][1] = buf[m++];
      xref_[i][2] = buf[m++];
      
      if (num_bond) {
        num_bond[i] = static_cast<int>(buf[m++]);
        for (int ii = 0; ii < lammpsInterface_->bond_per_atom(); ii++) {
          bond_atom[i][ii] = static_cast<int>(buf[m++]);
        }
      }

    }
  }

  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  void ATC_Method::reset_nlocal()
  {
    nLocalTotal_ = lammpsInterface_->nlocal();
    const int * mask = lammpsInterface_->atom_mask();
    nLocal_ = 0;
    nLocalGhost_ = 0;   

    for (int i = 0; i < nLocalTotal_; ++i) {
      if (mask[i] & groupbit_) nLocal_++;
      if (mask[i] & groupbitGhost_) nLocalGhost_++;
    }

    // set up internal & ghost maps
    
    if (nLocal_>0) {
      // set map
      internalToAtom_.resize(nLocal_);
      int j = 0;
      // construct internalToAtom map 
      //  : internal index -> local lammps atom index
      for (int i = 0; i < nLocalTotal_; ++i) {
        if (mask[i] & groupbit_) internalToAtom_(j++) = i;
      }
#ifdef EXTENDED_ERROR_CHECKING
      stringstream ss;
      ss << "Nlocal = " << nLocal_ << " but only found " << j << "atoms";
      if (j!=nLocal_) throw ATC_Error(ss.str());
#endif
      // construct reverse map
      atomToInternal_.clear();
      for (int i = 0; i < nLocal_; ++i) {
        atomToInternal_[internalToAtom_(i)] = i;
      }
    }
    if (nLocalGhost_>0) {
      // set map
      ghostToAtom_.resize(nLocalGhost_);
      int j = 0;
      for (int i = 0; i < nLocalTotal_; ++i) {
        if (mask[i] & groupbitGhost_) ghostToAtom_(j++) = i;
      }
    }

    //WIP_JAT this should not be needed at all, but a memory problem with sparse matrices requires them to be reset (possibly related to note in SparseMatrix-inl.h::_delete())
    interscaleManager_.reset_nlocal();

  }

  //-------------------------------------------------------------------
  void ATC_Method::reset_coordinates()
  {
    // update coarse graining positions for internal and ghost atoms
    atomCoarseGrainingPositions_->unfix_quantity();
    atomCoarseGrainingPositions_->quantity();
    atomCoarseGrainingPositions_->fix_quantity();
    if (atomGhostCoarseGrainingPositions_) {
      atomGhostCoarseGrainingPositions_->unfix_quantity();
      atomGhostCoarseGrainingPositions_->quantity();
      atomGhostCoarseGrainingPositions_->fix_quantity();
    }
     if (atomProcGhostCoarseGrainingPositions_) {
      atomProcGhostCoarseGrainingPositions_->unfix_quantity();
      atomProcGhostCoarseGrainingPositions_->quantity();
      atomProcGhostCoarseGrainingPositions_->fix_quantity();
    } 
  }

  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  void ATC_Method::write_atomic_weights(const string filename, const DIAG_MAT & atomicVolumeMatrix)
  {
    int nlocal = lammpsInterface_->nlocal();
    int nlocalmax;
    LammpsInterface::instance()->int_allmax(&nlocal,&nlocalmax);           
    int natoms = int(lammpsInterface_->natoms());
    ofstream out;
    const char* fname = &filename[0];
 
    // create tag to local id map for this processor
    map <int,int> id2tag;
    map <int,int>::const_iterator itr;
    int * atom_tag = lammpsInterface_->atom_tag();
    for (int i = 0; i < nlocal; ++i) {
      id2tag[i] = atom_tag[i];
    }

    int comm_rank = LammpsInterface::instance()->comm_rank();
    int nprocs;
    LammpsInterface::instance()->int_allmax(&comm_rank,&nprocs);
    nprocs += 1;

    if (comm_rank == 0) {
      out.open(fname);
      // print header lines
      out << "Atomic Weights for LAMMPS/atc analysis\n";
      out << " \n";
      out << natoms << " Atoms in system\n";
      out << " \n";
      // print atomic weights from proc 0
      for(int i = 0; i < nlocal; i++) {
        out << id2tag[i] << "  " << atomicVolumeMatrix(i,i) << "\n";
      }
    }
                
    if (nprocs > 1) {
      int max_size,send_size;
      send_size = nlocal;
      LammpsInterface::instance()->int_allmax(&send_size,&max_size);

      if (comm_rank == 0) {
        int intbuf[max_size];
        double buf[max_size];
        for (int iproc = 1; iproc < nprocs; iproc++) {
          LammpsInterface::instance()->int_recv(intbuf,max_size,iproc);
          LammpsInterface::instance()->recv(buf,max_size,iproc);
          for (int i = 0; i < max_size; i++) {
            out << intbuf[i] << "  " << buf[i] << "\n";
          }  
        }
      } else {
        int intbuf[send_size];
        double buf[send_size];
        for (int i = 0; i < send_size; i++) {
          intbuf[i] = id2tag[i];
          buf[i] = atomicVolumeMatrix(i,i);
        }
        LammpsInterface::instance()->int_send(intbuf,send_size);
        LammpsInterface::instance()->send(buf,send_size);
      }
    }
                
    if (comm_rank == 0) { 
      out.close();
    }           
  }
        
  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  void ATC_Method::compute_consistent_md_mass_matrix(const SPAR_MAT & shapeFunctionMatrix,
                                                       SPAR_MAT & mdMassMatrix) const
  {
    
    int nCols = shapeFunctionMatrix.nCols();
    DENS_MAT massMatrixLocal(nCols,nCols);
    DENS_MAT denseMdMassMatrix(nCols,nCols);
    if (nLocal_>0)
      massMatrixLocal = shapeFunctionMatrix.transMat(shapeFunctionMatrix);
    
    lammpsInterface_->allsum(massMatrixLocal.ptr(),
                             denseMdMassMatrix.ptr(),
                             denseMdMassMatrix.size());
    mdMassMatrix.reset(denseMdMassMatrix,1.e-10);
  }

  //=================================================================
  // Interscale operators
  //=================================================================
  // in the spirit of the current design of ATC: atoms local, nodes global
  
  
  
  
  bool ATC_Method::nodal_influence(const int groupbit,
                              set<int> & nset, set<int> & aset, double tol)
  {
    int nghost = nodal_influence(groupbit,nset,aset,true,tol);
    int local_nghost = nghost;
    lammpsInterface_->int_allsum(&local_nghost,&nghost);
    if (nghost == 0) {
       nodal_influence(groupbit,nset,aset,false,tol);
    }
    return (nghost > 0);
  }
  int ATC_Method::nodal_influence(const int groupbit,
        set<int> & nset, set<int> & aset, bool ghost, double tol)
  {
    Array<int> & amap = (ghost) ? ghostToAtom_ : internalToAtom_;
    int natoms = (ghost) ? nLocalGhost_ : nLocal_;
    DENS_MAT influence(nNodes_,1);
    DENS_MAT atomInfluence(natoms,1);
    const int *mask = lammpsInterface_->atom_mask();
    for (int i = 0; i < natoms; i++) {  
      if (mask[amap(i)] & groupbit) {
         atomInfluence(i,0) = 1;
         aset.insert(i);
      }
    }
    // relies on shape functions
    if (ghost) {
      restrict_volumetric_quantity(atomInfluence,influence,(interscaleManager_.per_atom_sparse_matrix("InterpolantGhost"))->quantity());
    }
    else {
    restrict_volumetric_quantity(atomInfluence,influence); 
    }

    DENS_MAT localInfluence = influence;
    lammpsInterface_->allsum(localInfluence.ptr(),
                             influence.ptr(),
                             influence.size());

    for (int i = 0; i < nNodes_; i++) {  
      if (fabs(influence(i,0)) > tol)  { nset.insert(i);  }
    }
    return aset.size();
  }


  //--------------------------------------------------------
  void ATC_Method::restrict_volumetric_quantity(const MATRIX & atomData,
                                                MATRIX & nodeData,
                                                const SPAR_MAT & shpFcn)
  {
    // computes nodeData = N*DeltaVAtom*atomData where N are the shape functions
    DENS_MAT workNodeArray(nodeData.nRows(),nodeData.nCols()); 
    //DENS_MAT workNodeArray;


    if (atomData.nRows() > 0) { // or shpFcn_??? 
      workNodeArray = shpFcn.transMat(atomData);
    }
    int count = nodeData.nRows()*nodeData.nCols();
    lammpsInterface_->allsum(workNodeArray.ptr(),nodeData.ptr(),count);
    return;
  }


  //--------------------------------------------------------
  void ATC_Method::restrict_volumetric_quantity(const MATRIX & atomData,
                                                MATRIX & nodeData)
  {
    restrict_volumetric_quantity(atomData,nodeData,shpFcn_->quantity());
    return;
  }


  //--------------------------------------------------------
  void ATC_Method::prolong(const MATRIX & nodeData,
                             MATRIX & atomData)
  {
    // computes the finite element interpolation at atoms atomData = N*nodeData
    if (nLocal_>0) {
      atomData = (shpFcn_->quantity())*nodeData;
    }
    return;
  }

  //========================================================
  // FE functions
  //========================================================

  //--------------------------------------------------------
  void ATC_Method::output()
  {
//  if (lammpsInterface_->comm_rank() == 0) {
      compute_nodeset_output();
      compute_faceset_output();
      compute_elementset_output();
//  }
  }
  //--------------------------------------------------------
  void ATC_Method::compute_nodeset_output(void) 
  {
    map< pair <string, FieldName>, NodesetOperationType >::const_iterator iter;
    for (iter = nsetData_.begin(); iter != nsetData_.end();iter++){
      pair <string, FieldName> id = iter->first;
      string nsetName = id.first;
      FieldName field = id.second;
      double sum = 0.0;
      const set<int> nset = feEngine_->fe_mesh()->nodeset(nsetName);
      const DENS_MAT & thisField = (fields_.find(field)->second).quantity();
      set< int >::const_iterator itr;
      for (itr = nset.begin(); itr != nset.end();itr++){
        int node = *itr;
        sum += thisField(node,0);
      }
      string name = nsetName + "_" + field_to_string(field);
      if (iter->second == NODESET_AVERAGE) {
        sum /= nset.size();
        name = "average_"+name;
      }
      feEngine_->add_global(name, sum);
    }
  }
  //--------------------------------------------------------
  void ATC_Method::compute_faceset_output(void)
  {
    map < pair<string,string>, FacesetIntegralType >::const_iterator iter;
    DENS_MAT values;
    for (iter = fsetData_.begin(); iter !=  fsetData_.end(); iter++) {
      string bndyName  = (iter->first).first;
      string fieldName = (iter->first).second;
      const set< PAIR > & faceSet = (feEngine_->fe_mesh())->faceset(bndyName);
      ATOMIC_DATA::iterator data_iter = filteredData_.find(fieldName);
      if (data_iter == filteredData_.end()) {
        string msg = "Specified fieldName "+fieldName+
          " not found in filteredData_ while attempting surface integration";
        throw ATC_Error(msg);
      }
      const DENS_MAT & data =  ((data_iter->second).quantity());
      string stem = bndyName + "_" + fieldName + "_";
      bool tf = false;
      if (iter->second == CONTOUR_INTEGRAL) {
        stem = "contour_"+stem;
        tf = true;
      }
      feEngine_->field_surface_flux(data,faceSet,values,tf);
      for (int i = 0; i < values.nRows() ; ++i ) {
        string name = stem + to_string(i+1);
        feEngine_->add_global(name, values(i,0));
      }
    }
  }
  //--------------------------------------------------------
  void ATC_Method::compute_elementset_output(void) 
  {
    map< pair <string, FieldName>, ElementsetOperationType >::const_iterator iter;
    for (iter = esetData_.begin(); iter != esetData_.end();iter++){
      pair <string, FieldName> id = iter->first;
      string esetName = id.first;
      FieldName field = id.second;
      const ESET eset = feEngine_->fe_mesh()->elementset(esetName);
      const DENS_MAT & thisField = (fields_.find(field)->second).quantity();
      DENS_VEC total = feEngine_->integrate(thisField,eset);
      string name = esetName + "_" + field_to_string(field);
      if (iter->second == ELEMENTSET_AVERAGE) {
        DENS_MAT ones(nNodes_,0); ones = 1;
        DENS_VEC V = feEngine_->integrate(ones,eset);
        total /= V[0];
        name = "average_"+name;
      }
      if (total.size() == 1) {
        feEngine_->add_global(name, total[0]);
      }
      else {
        for (int i = 0; i < total.size(); i++) {
          feEngine_->add_global(name+to_string(i), total[i]);
        }
      }
    }
  }


  //=================================================================
  // 
  //=================================================================
  //--------------------------------------------------------
  bool ATC_Method::read_atomic_ref_positions(const char * filename)
  {
    int nlocal = lammpsInterface_->nlocal();
    ifstream in;
    const int lineSize = 256;
    char line[lineSize];

    // create tag to local id map for this processor
    map <int,int> tag2id;
    map <int,int>::const_iterator itr;
    int * atom_tag = lammpsInterface_->atom_tag();
    for (int i = 0; i < nlocal; ++i) {
      tag2id[atom_tag[i]] = i;
    }

    // get number of atoms
    int natoms = 0;
    if (LammpsInterface::instance()->rank_zero()) {
      in.open(filename);
      string msg;
      string name = filename;
      msg = "no "+name+" reference file found";
      if (! in.good()) throw ATC_Error(msg);
      in.getline(line,lineSize); // header
      in.getline(line,lineSize); // blank line
      in >> natoms;
      in.close();
      stringstream ss; 
      ss << "found " << natoms << " atoms in reference file";
      LammpsInterface::instance()->print_msg(ss.str());
    }
    LammpsInterface::instance()->int_broadcast(&natoms);

    // read atoms and assign
    if (LammpsInterface::instance()->rank_zero()) {
      in.open(filename); 
      while(in.good()) {
        in.getline(line,lineSize);
        string str(line);
        int pos = str.find("Atoms");           
        if (pos > -1) {
          in.getline(line,lineSize); // blank line
          break;
        }
      }
    }
    int nread = 0,type = -1, tag = -1, count = 0;
    double x[3]={0,0,0};
    while (nread < natoms) {
      if (LammpsInterface::instance()->rank_zero()) {
         in.getline(line,lineSize);
         stringstream ss (line,stringstream::in | stringstream::out);
         ss >> tag >> type >> x[0] >> x[1] >> x[2]; 
         nread++;
      }
      LammpsInterface::instance()->int_broadcast(&nread);
      LammpsInterface::instance()->int_broadcast(&tag);
      LammpsInterface::instance()->broadcast(x,3);
      itr = tag2id.find(tag);
      if (itr != tag2id.end()) {
        int iatom = itr->second;
        xref_[iatom][0] = x[0];
        xref_[iatom][1] = x[1];
        xref_[iatom][2] = x[2];
        count++;
      }
    }
    if (LammpsInterface::instance()->rank_zero()) {
      in.close();
      stringstream ss; 
      ss << "read  " << natoms << " reference positions";
      LammpsInterface::instance()->print_msg(ss.str());
    }
    if (count != nlocal) 
       throw ATC_Error("reset "+to_string(count)+" atoms vs "+to_string(nlocal));


    return true;
  }  

//--------------------------------------------------------
  void ATC_Method::remap_ghost_ref_positions(void)
  {
    
    int nlocal = lammpsInterface_->nlocal();
    int nghost = lammpsInterface_->nghost();

    double box_bounds[2][3];
    lammpsInterface_->box_bounds(box_bounds[0][0],box_bounds[1][0],
                                    box_bounds[0][1],box_bounds[1][1],
                                    box_bounds[0][2],box_bounds[1][2]);
    double xlo = box_bounds[0][0], xhi = box_bounds[1][0];
    double ylo = box_bounds[0][1], yhi = box_bounds[1][1];
    double zlo = box_bounds[0][2], zhi = box_bounds[1][2];

    double box_length[3];
    for (int k = 0; k < 3; k++) {
      box_length[k] = box_bounds[1][k] - box_bounds[0][k];
    }
    double Lx = box_length[0], Ly = box_length[1], Lz = box_length[2];

    // create tag to local id map for this processor
    map <int,int> tag2id;
    map <int,int>::const_iterator itr;
    int * atom_tag = lammpsInterface_->atom_tag();
    for (int i = 0; i < nlocal; ++i) {
      tag2id[atom_tag[i]] = i;
    }
    
    // loop over ghosts
    double ** x = lammpsInterface_->xatom();
    for (int j = nlocal; j < nlocal + nghost; j++) {
      int tag = atom_tag[j];
      int i = tag2id[tag];
      //itr = tag2id.find(tag);
      //if (itr != tag2id.end()) 
      double* xj = x[j];
      double* Xj = xref_[j];
      //double Xj[3];
      double* Xi = xref_[i];
      // the assumption is that xref_[j] has been shuffled 
      // so make an image of xref_[i] that is close to x[j]
      if (xj[0] <= xlo) Xj[0] = Xi[0] -Lx; 
      if (xj[0] >= xhi) Xj[0] = Xi[0] +Lx; 
      if (xj[1] <= ylo) Xj[1] = Xi[1] -Ly; 
      if (xj[1] >= yhi) Xj[1] = Xi[1] +Ly; 
      if (xj[2] <= zlo) Xj[2] = Xi[2] -Lz;
      if (xj[2] >= zhi) Xj[2] = Xi[2] +Lz; 
    }
  }
};
