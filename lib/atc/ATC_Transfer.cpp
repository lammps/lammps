// ATC_Transfer headers
#include "ATC_Transfer.h"
#include "ATC_Error.h"
#include "FE_Engine.h"
#include "LammpsInterface.h"
#include "Quadrature.h"
#include "VoigtOperations.h"
#include "TransferLibrary.h"
#include "Stress.h"
#include "KernelFunction.h"
#include "PerPairQuantity.h"
#include "FieldManager.h"
#define ESHELBY_VIRIAL
#include "LinearSolver.h"


// Other Headers
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <fstream>
#include <sstream>
#include <exception>



// PLAN:
//* energies
//* filters - make filterFields class
//* output directly
//* enum, tagged, computes, mat(field to field) functions
//* grads & rates
//* on-the-fly
// * remove derived classes

using namespace std;
using namespace ATC_Utility;
using namespace voigt3;
//using ATC_Utility::to_string;
//using voigt3::vector_to_matrix;
//using voigt3::vector_to_symm_matrix;
//using voigt3::matrix_to_vector;
//using voigt3::symm_matrix_to_vector;

namespace ATC {

  const int numFields_ = 16;
  FieldName indices_[numFields_] = {
    CHARGE_DENSITY,
    MASS_DENSITY,
    SPECIES_CONCENTRATION,
    NUMBER_DENSITY,
    MOMENTUM,
    VELOCITY,
    PROJECTED_VELOCITY,
    DISPLACEMENT,
    POTENTIAL_ENERGY,
    KINETIC_ENERGY,
    KINETIC_TEMPERATURE,
    TEMPERATURE,
    CHARGE_FLUX,
    SPECIES_FLUX,
    THERMAL_ENERGY};
    //ELECTRIC_POTENTIAL};

  ATC_Transfer::ATC_Transfer(string groupName,
                             double ** & perAtomArray,
                             LAMMPS_NS::Fix * thisFix,
                             string matParamFile)
    : ATC_Method(groupName,perAtomArray,thisFix),
      xPointer_(NULL), 
      outputStepZero_(true),
      neighborReset_(false),
      pairMap_(NULL),
      bondMatrix_(NULL),
      pairVirial_(NULL),
      pairHeatFlux_(NULL),
      nComputes_(0),
      hasPairs_(true),
      hasBonds_(false),
      resetKernelFunction_(false),
      dxaExactMode_(true),
      cauchyBornStress_(NULL)
  {
    nTypes_ = lammpsInterface_->ntypes();

    peScale_=1.;
    keScale_= lammpsInterface_->mvv2e();
    // if surrogate model of md (no physics model created)
    if (matParamFile != "none") {
      fstream  fileId(matParamFile.c_str(), std::ios::in);
      if (!fileId.is_open()) throw ATC_Error("cannot open material file");
      CbData cb;
      LammpsInterface *lmp = LammpsInterface::instance();
      lmp->lattice(cb.cell_vectors, cb.basis_vectors);
      cb.inv_atom_volume = 1.0 / lmp->volume_per_atom();
      cb.e2mvv           = 1.0 / lmp->mvv2e();
      cb.atom_mass       = lmp->atom_mass(1); 
      cb.boltzmann       = lmp->boltz();
      cb.hbar            = lmp->hbar();
      cauchyBornStress_ = new StressCauchyBorn(fileId, cb);
    }

    // Defaults
    set_time();
 
    outputFlags_.reset(NUM_TOTAL_FIELDS);
    outputFlags_ = false;
    fieldFlags_.reset(NUM_TOTAL_FIELDS);
    fieldFlags_ = false;
    gradFlags_.reset(NUM_TOTAL_FIELDS);
    gradFlags_ = false;
    rateFlags_.reset(NUM_TOTAL_FIELDS);
    rateFlags_ = false;

    outputFields_.resize(NUM_TOTAL_FIELDS);
    for (int i = 0; i < NUM_TOTAL_FIELDS; i++) { outputFields_[i] = NULL; }

    // Hardy requires ref positions for processor ghosts for bond list
    
    //needXrefProcessorGhosts_ = true;
  }

  //-------------------------------------------------------------------
  ATC_Transfer::~ATC_Transfer()
  {
    interscaleManager_.clear(); 
    if (cauchyBornStress_) delete cauchyBornStress_; 
  }

  //-------------------------------------------------------------------
  // called before the beginning of a "run"
  void ATC_Transfer::initialize() 
  {
    ATC_Method::initialize();

    if (!initialized_) { 
      if (cauchyBornStress_) cauchyBornStress_->initialize();
    }

    if (!initialized_ || ATC::LammpsInterface::instance()->atoms_sorted() || resetKernelFunction_) {
      // initialize kernel funciton matrix N_Ia
      if (! kernelOnTheFly_) {
        try{
          if (!moleculeIds_.empty()) compute_kernel_matrix_molecule(); //KKM add
        }
        catch(bad_alloc&) {
          ATC::LammpsInterface::instance()->print_msg("kernel will be computed on-the-fly");
          kernelOnTheFly_ = true;
        }
      }
      resetKernelFunction_ = false;
    }

    // initialize bond matrix B_Iab
    if ((! bondOnTheFly_) 
       && ( ( fieldFlags_(STRESS) 
           || fieldFlags_(ESHELBY_STRESS) 
           || fieldFlags_(HEAT_FLUX) ) ) ) {
      try {
        compute_bond_matrix(); 
      } 
      catch(bad_alloc&) { 
        ATC::LammpsInterface::instance()->print_msg("stress/heat_flux will be computed on-the-fly");
        
        bondOnTheFly_ = true;
      }
    }

    // set sample frequency to output if sample has not be specified
    if (sampleFrequency_ == 0) sampleFrequency_ = outputFrequency_;

    // output for step 0 
    if (!initialized_) {
      if (outputFrequency_ > 0) {
        // initialize filtered data
        compute_fields();
        for(int index=0; index < NUM_TOTAL_FIELDS; ++index) {
          if(fieldFlags_(index)) {
            string name = field_to_string((FieldName) index);
            filteredData_[name] = hardyData_[name];
            timeFilters_(index)->initialize(filteredData_[name].quantity());
          }
          if (rateFlags_(index)) {
            string name = field_to_string((FieldName) index);
            string rate_field = name + "_rate";
            filteredData_[rate_field] = hardyData_[rate_field];
          }
          if (gradFlags_(index)) {
            string name = field_to_string((FieldName) index);
            string grad_field = name + "_gradient";
            filteredData_[grad_field] = hardyData_[grad_field];
          }
        }
        int index = NUM_TOTAL_FIELDS;
        map <string,int>::const_iterator iter;
        for (iter = computes_.begin(); iter != computes_.end(); iter++) {
          string tag = iter->first;
          filteredData_[tag] = hardyData_[tag];
          timeFilters_(index)->initialize(filteredData_[tag].quantity());
#ifdef ESHELBY_VIRIAL
          if (tag == "virial" && fieldFlags_(ESHELBY_STRESS)) {
            filteredData_["eshelby_virial"] = hardyData_["eshelby_virial"];
          }
#endif
          index++;
        }
        output();
      }
    }

    initialized_ = true;

    lammpsInterface_->computes_addstep(lammpsInterface_->ntimestep()+sampleFrequency_);

    
    //remap_ghost_ref_positions();
    update_peratom_output(); 
  }

  //-------------------------------------------------------------------
  void ATC_Transfer::set_continuum_data()
  {
    ATC_Method::set_continuum_data();
    if (!initialized_) {
      nNodesGlobal_ = feEngine_->fe_mesh()->num_nodes();
    }
  }

  //-------------------------------------------------------------------
  void ATC_Transfer::construct_time_integration_data()
  {
    if (!initialized_) {
      // ground state for PE
      nodalRefPotentialEnergy_.reset(nNodes_,1);
      
      // size arrays for requested/required  fields
      for(int index=0; index < NUM_TOTAL_FIELDS; ++index) {
        if (fieldFlags_(index)) {

          int size = FieldSizes[index];
          if (atomToElementMapType_ == EULERIAN) {
            if (index == STRESS) size=6;
            if (index == CAUCHY_BORN_STRESS) size=6;
          }
          if (size == 0) {
            if (index == SPECIES_CONCENTRATION) size=typeList_.size()+groupList_.size();
          }
          string name = field_to_string((FieldName) index);
          hardyData_   [name].reset(nNodes_,size);
          filteredData_[name].reset(nNodes_,size);
        }
      }

      // size arrays for projected compute fields
      map <string,int>::const_iterator iter;
      for (iter = computes_.begin(); iter != computes_.end(); iter++) {
        string tag = iter->first;
        COMPUTE_POINTER cmpt = lammpsInterface_->compute_pointer(tag);
        int ncols = lammpsInterface_->compute_ncols_peratom(cmpt);
        hardyData_        [tag].reset(nNodes_,ncols);
        filteredData_[tag].reset(nNodes_,ncols);
#ifdef ESHELBY_VIRIAL
        if (tag == "virial" && fieldFlags_(ESHELBY_STRESS)) {
          string esh = "eshelby_virial";
          int size = FieldSizes[ESHELBY_STRESS];
          hardyData_   [esh].reset(nNodes_,size);
          filteredData_[esh].reset(nNodes_,size);
        }
#endif
      }
    }
  }
  //--------------------------------------------------------
  //  set_computational_geometry
  //    constructs needed transfer operators which define
  //    hybrid atom/FE computational geometry
  //--------------------------------------------------------
  void ATC_Transfer::set_computational_geometry()
  {
    ATC_Method::set_computational_geometry();
  }

  //-------------------------------------------------------------------
  // constructs quantities
  void ATC_Transfer::construct_transfers()
  {

    // set pointer to positions
    // REFACTOR use method's handling of xref/xpointer
    set_xPointer(); 

    ATC_Method::construct_transfers();

    if (!(kernelOnTheFly_)) {
      // finite element shape functions for interpolants 
      PerAtomShapeFunction * atomShapeFunctions = new PerAtomShapeFunction(this);
      interscaleManager_.add_per_atom_sparse_matrix(atomShapeFunctions,"Interpolant");
      shpFcn_ = atomShapeFunctions;
    }

    
    this->create_atom_volume();

    // accumulants
    if (kernelFunction_) {
      // kernel-based accumulants
      if (kernelOnTheFly_) {
        ConstantQuantity<double> * atomCount = new ConstantQuantity<double>(this,1.);
        interscaleManager_.add_per_atom_quantity(atomCount,"AtomCount");
        OnTheFlyKernelAccumulation * myWeights 
           = new OnTheFlyKernelAccumulation(this,
             atomCount, kernelFunction_, atomCoarseGrainingPositions_);
        interscaleManager_.add_dense_matrix(myWeights,
                                            "KernelInverseWeights");
        accumulantWeights_ = new OnTheFlyKernelWeights(myWeights);
      }
      else {
        PerAtomKernelFunction * atomKernelFunctions = new PerAtomKernelFunction(this);
        interscaleManager_.add_per_atom_sparse_matrix(atomKernelFunctions,
                                                      "Accumulant");
        accumulant_ = atomKernelFunctions;
        accumulantWeights_ = new AccumulantWeights(accumulant_);
      }
      accumulantInverseVolumes_ = new KernelInverseVolumes(this,kernelFunction_);
      interscaleManager_.add_diagonal_matrix(accumulantInverseVolumes_,
                                            "AccumulantInverseVolumes");
      interscaleManager_.add_diagonal_matrix(accumulantWeights_,
                                             "AccumulantWeights");
    }
    else {
      // mesh-based accumulants
      if (kernelOnTheFly_) {
        ConstantQuantity<double> * atomCount = new ConstantQuantity<double>(this,1.);
        interscaleManager_.add_per_atom_quantity(atomCount,"AtomCount");
        OnTheFlyMeshAccumulation * myWeights 
           = new OnTheFlyMeshAccumulation(this,
             atomCount, atomCoarseGrainingPositions_);
        interscaleManager_.add_dense_matrix(myWeights,
                                            "KernelInverseWeights");
        accumulantWeights_ = new OnTheFlyKernelWeights(myWeights);
      } else {
        accumulant_ = shpFcn_;
        accumulantWeights_ = new AccumulantWeights(accumulant_);
        interscaleManager_.add_diagonal_matrix(accumulantWeights_,
                                              "AccumulantWeights");
      }
    }

    
    bool needsBondMatrix =  (! bondOnTheFly_ ) &&
             (fieldFlags_(STRESS)
           || fieldFlags_(ESHELBY_STRESS)
           || fieldFlags_(HEAT_FLUX));
    if (needsBondMatrix) {
      if (hasPairs_ && hasBonds_) {
        pairMap_ = new PairMapBoth(lammpsInterface_,groupbit_);
      }
      else if (hasBonds_) {
        pairMap_ = new PairMapBond(lammpsInterface_,groupbit_);
      }
      else if (hasPairs_) {
        pairMap_ = new PairMapNeighbor(lammpsInterface_,groupbit_);
      }
    }
    if (pairMap_) interscaleManager_.add_pair_map(pairMap_,"PairMap");
    //if (pairMap_ && !initialized_) interscaleManager_.add_pair_map(pairMap_,"PairMap");


    //const  PerAtomQuantity<double> *  x0= interscaleManager_.per_atom_quantity("AtomicReferencePositions");
    //const  PerAtomQuantity<double> *  x0= interscaleManager_.per_atom_quantity("AtomicCoarseGrainingPositions");
    //const  PerAtomQuantity<double> *  x0= interscaleManager_.per_atom_quantity("AtomicReferencePositions");
      
    if ( fieldFlags_(STRESS) || fieldFlags_(ESHELBY_STRESS) || fieldFlags_(HEAT_FLUX) ) {

      const FE_Mesh * fe_mesh = feEngine_->fe_mesh();
      if (!kernelBased_) {
        bondMatrix_ = new BondMatrixPartitionOfUnity(lammpsInterface_,*pairMap_,xPointer_,fe_mesh,accumulantInverseVolumes_); 
      }
      else {
        bondMatrix_ = new BondMatrixKernel(lammpsInterface_,*pairMap_,xPointer_,fe_mesh,kernelFunction_);
      }
    }
    if (bondMatrix_) interscaleManager_.add_sparse_matrix(bondMatrix_,"BondMatrix");

    if ( fieldFlags_(STRESS) || fieldFlags_(ESHELBY_STRESS) ) {
      if (atomToElementMapType_ == LAGRANGIAN) {
//      pairVirial_ = new PairVirialLagrangian(lammpsInterface_,*pairMap_,x0);
        pairVirial_ = new PairVirialLagrangian(lammpsInterface_,*pairMap_,xref_);
      }
      else if (atomToElementMapType_ == EULERIAN) {
        pairVirial_ = new PairVirialEulerian(lammpsInterface_,*pairMap_);
      }
      else {
        throw ATC_Error("no atom to element map specified");
      }
    }
    if (pairVirial_) interscaleManager_.add_dense_matrix(pairVirial_,"PairVirial");

    if ( fieldFlags_(HEAT_FLUX) ) {
      if (atomToElementMapType_ == LAGRANGIAN) {
        pairHeatFlux_ = new PairPotentialHeatFluxLagrangian(lammpsInterface_,*pairMap_,xref_);
      }
      else if (atomToElementMapType_ == EULERIAN) {
        pairHeatFlux_ = new PairPotentialHeatFluxEulerian(lammpsInterface_,*pairMap_);
      }
      else {
        throw ATC_Error("no atom to element map specified");
      }
    }
    if (pairHeatFlux_) interscaleManager_.add_dense_matrix(pairHeatFlux_,"PairHeatFlux");

    // gradient matrix
    if (gradFlags_.has_member(true)) {
      NativeShapeFunctionGradient * gradientMatrix = new NativeShapeFunctionGradient(this);
      interscaleManager_.add_vector_sparse_matrix(gradientMatrix,"GradientMatrix");
      gradientMatrix_ = gradientMatrix;
    }

    // molecule centroid, molecule charge, dipole moment and quadrupole moment calculations KKM add
    if (!moleculeIds_.empty()) {
      map<string,pair<MolSize,int> >::const_iterator molecule;
      InterscaleManager & interscaleManager = this->interscale_manager(); // KKM add, may be we do not need this as interscaleManager_ already exists.  
      PerAtomQuantity<double> * atomProcGhostCoarseGrainingPositions_ = interscaleManager.per_atom_quantity("AtomicProcGhostCoarseGrainingPositions");
      FundamentalAtomQuantity * mass = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,PROC_GHOST);
      molecule = moleculeIds_.begin();
      int groupbit = (molecule->second).second;
      smallMoleculeSet_ = new SmallMoleculeSet(this,groupbit);
      smallMoleculeSet_->initialize(); // KKM add, why should we? 
      interscaleManager_.add_small_molecule_set(smallMoleculeSet_,"MoleculeSet");
      moleculeCentroid_ = new SmallMoleculeCentroid(this,mass,smallMoleculeSet_,atomProcGhostCoarseGrainingPositions_);
      interscaleManager_.add_dense_matrix(moleculeCentroid_,"MoleculeCentroid");
      AtomToSmallMoleculeTransfer<double> * moleculeMass =
        new AtomToSmallMoleculeTransfer<double>(this,mass,smallMoleculeSet_);
      interscaleManager_.add_dense_matrix(moleculeMass,"MoleculeMass");
      FundamentalAtomQuantity * atomicCharge = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_CHARGE,PROC_GHOST);
      AtomToSmallMoleculeTransfer<double> * moleculeCharge =
        new AtomToSmallMoleculeTransfer<double>(this,atomicCharge,smallMoleculeSet_);
      interscaleManager_.add_dense_matrix(moleculeCharge,"MoleculeCharge");
      dipoleMoment_ = new SmallMoleculeDipoleMoment(this,atomicCharge,smallMoleculeSet_,atomProcGhostCoarseGrainingPositions_,moleculeCentroid_);
      interscaleManager_.add_dense_matrix(dipoleMoment_,"DipoleMoment");
      quadrupoleMoment_ = new SmallMoleculeQuadrupoleMoment(this,atomicCharge,smallMoleculeSet_,atomProcGhostCoarseGrainingPositions_,moleculeCentroid_);
      interscaleManager_.add_dense_matrix(quadrupoleMoment_,"QuadrupoleMoment");
    }

    FieldManager fmgr(this);

//  for(int index=0; index < NUM_TOTAL_FIELDS; ++index) 
    for(int i=0; i < numFields_; ++i) {
      FieldName index = indices_[i];
      if (fieldFlags_(index)) {
        outputFields_[index] = fmgr.nodal_atomic_field(index);
      }
    }

// WIP REJ
    if (fieldFlags_(ELECTRIC_POTENTIAL)) {
      restrictedCharge_ = fmgr.restricted_atom_quantity(CHARGE_DENSITY);
    }

    // computes
    map <string,int>::const_iterator iter;
    for (iter = computes_.begin(); iter != computes_.end(); iter++) {
      string tag = iter->first;
      ComputedAtomQuantity * c = new ComputedAtomQuantity(this, tag);
      interscaleManager_.add_per_atom_quantity(c,tag);
      int projection = iter->second;
      DIAG_MAN * w = NULL;
      if      (projection == VOLUME_NORMALIZATION ) 
         { w = accumulantInverseVolumes_; }
      else if (projection == NUMBER_NORMALIZATION ) 
         { w = accumulantWeights_; }
      if (kernelFunction_ && kernelOnTheFly_) {
        OnTheFlyKernelAccumulationNormalized * C = new OnTheFlyKernelAccumulationNormalized(this, c, kernelFunction_, atomCoarseGrainingPositions_, w);
        interscaleManager_.add_dense_matrix(C,tag);
        outputFieldsTagged_[tag] = C;
      }
      else {
        AtfProjection * C =  new AtfProjection(this, c, accumulant_, w);
        interscaleManager_.add_dense_matrix(C,tag);
        outputFieldsTagged_[tag] = C;
      }
    }
  
  }

  //-------------------------------------------------------------------
  // sets initial values of filtered quantities
  void ATC_Transfer::construct_methods()
  {
    if ((!initialized_) || timeFilterManager_.need_reset()) {
      timeFilters_.reset(NUM_TOTAL_FIELDS+nComputes_);
      sampleCounter_ = 0;
      
      // for filtered  fields
      for(int index=0; index < NUM_TOTAL_FIELDS; ++index) {
        if (fieldFlags_(index)) {
          string name = field_to_string((FieldName) index);
          filteredData_[name]  = 0.0;
          timeFilters_(index) = timeFilterManager_.construct(); 
        }
      }
      
      // for filtered projected computes
      
      //       lists/accessing of fields ( & computes)
      map <string,int>::const_iterator iter;
      int index = NUM_TOTAL_FIELDS;
      for (iter = computes_.begin(); iter != computes_.end(); iter++) {
        string tag = iter->first;
        filteredData_[tag] = 0.0;
        timeFilters_(index) = timeFilterManager_.construct();
        index++;
      }
    }
  }


  //-------------------------------------------------------------------
  // called after the end of a "run"
  void ATC_Transfer::finish() 
  {
    // base class
    ATC_Method::finish();
  }

  //-------------------------------------------------------------------
  // this is the parser
  bool ATC_Transfer::modify(int narg, char **arg)
  {
    bool match = false;

    int argIdx = 0;
    // check to see if it is a transfer class command
      /*! \page man_hardy_fields fix_modify AtC fields 
        \section syntax
        fix_modify AtC fields <all | none> \n
        fix_modify AtC fields <add | delete> <list_of_fields> \n
        - all | none (keyword) = output all or no fields  \n
        - add | delete (keyword) = add or delete the listed output fields \n
        - fields (keyword) =  \n
        density : mass per unit volume \n
        displacement : displacement vector \n
        momentum : momentum per unit volume \n
        velocity : defined by momentum divided by density \n
        projected_velocity : simple kernel estimation of atomic velocities \n
        temperature :  temperature derived from the relative atomic kinetic energy (as done by ) \n
        kinetic_temperature : temperature derived from the full kinetic energy  \n
        number_density : simple kernel estimation of number of atoms per unit volume \n
        stress : 
        Cauchy stress tensor for eulerian analysis (atom_element_map), or
        1st Piola-Kirchhoff stress tensor for lagrangian analysis    \n
        transformed_stress : 
        1st Piola-Kirchhoff stress tensor for eulerian analysis (atom_element_map), or 
                             Cauchy stress tensor for lagrangian analysis    \n
        heat_flux : spatial heat flux vector for eulerian, 
        or referential heat flux vector for lagrangian \n
        potential_energy : potential energy per unit volume \n
        kinetic_energy : kinetic energy per unit volume \n
        thermal_energy : thermal energy (kinetic energy - continuum kinetic energy) per unit volume \n
        internal_energy : total internal energy (potential + thermal) per unit volume \n
        energy : total energy (potential + kinetic) per unit volume \n
        number_density : number of atoms per unit volume \n
        eshelby_stress: configurational stress (energy-momentum) tensor defined by Eshelby
          [References: Philos. Trans. Royal Soc. London A, Math. Phys. Sci., Vol. 244, 
          No. 877 (1951) pp. 87-112; J. Elasticity, Vol. 5, Nos. 3-4 (1975) pp. 321-335] \n
        vacancy_concentration: volume fraction of vacancy content \n
        type_concentration: volume fraction of a specific atom type \n
        \section examples
        <TT> fix_modify AtC fields add velocity temperature </TT>
        \section description
        Allows modification of the fields calculated and output by the 
        transfer class. The commands are cumulative, e.g.\n
        <TT> fix_modify AtC fields none </TT> \n 
        followed by \n 
        <TT> fix_modify AtC fields add velocity temperature </TT> \n
        will only output the velocity and temperature fields.
        \section restrictions
        Must be used with the hardy/field type of AtC fix, see \ref man_fix_atc.
        Currently, the stress and heat flux formulas are only correct for 
        central force potentials, e.g. Lennard-Jones and EAM 
        but not Stillinger-Weber.
        \section related
        See \ref man_hardy_gradients , \ref man_hardy_rates  and \ref man_hardy_computes
        \section default
        By default, no fields are output
      */
      if (strcmp(arg[argIdx],"fields")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"all")==0) { 
          outputFlags_ = true;
          match = true;
        } 
        else if (strcmp(arg[argIdx],"none")==0) { 
          outputFlags_ = false;
          match = true;
        } 
        else if (strcmp(arg[argIdx],"add")==0) { 
          argIdx++;
          for (int i = argIdx; i < narg; ++i) {
            FieldName field_name = string_to_field(arg[i]);
            outputFlags_(field_name) = true; 
          }
          match = true;
        } 
        else if (strcmp(arg[argIdx],"delete")==0) { 
          argIdx++;
          for (int i = argIdx; i < narg; ++i) {
            FieldName field_name = string_to_field(arg[i]);
            outputFlags_(field_name) = false; 
          }
          match = true;
        } 
        check_field_dependencies();
        if (fieldFlags_(DISPLACEMENT)) { trackDisplacement_ = true; }
      }

      /*! \page man_hardy_gradients fix_modify AtC gradients
        \section syntax
        fix_modify AtC gradients <add | delete> <list_of_fields> \n
        - add | delete (keyword) = add or delete the calculation of gradients for the listed output fields \n
        - fields (keyword) =  \n
        gradients can be calculated for all fields listed in \ref man_hardy_fields
 
        \section examples
        <TT> fix_modify AtC gradients add temperature velocity stress </TT> \n
        <TT> fix_modify AtC gradients delete velocity </TT> \n
        \section description
        Requests calculation and ouput of gradients of the fields from the 
        transfer class. These gradients will be with regard to spatial or material
        coordinate for eulerian or lagrangian analysis, respectively, as specified by 
        atom_element_map (see \ref man_atom_element_map )
        \section restrictions
        Must be used with the hardy/field type of AtC fix 
        ( see \ref man_fix_atc )
        \section related
        \section default
        No gradients are calculated by default
      */
      else if (strcmp(arg[argIdx],"gradients")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"add")==0) { 
          argIdx++;
          FieldName field_name;
          for (int i = argIdx; i < narg; ++i) {
            field_name = string_to_field(arg[i]);
            gradFlags_(field_name) = true; 
          }
          match = true;
        } 
        else if (strcmp(arg[argIdx],"delete")==0) { 
          argIdx++;
          FieldName field_name;
          for (int i = argIdx; i < narg; ++i) {
            field_name = string_to_field(arg[i]);
            gradFlags_(field_name) = false; 
          }
          match = true;
        } 
      }

      /*! \page man_hardy_rates fix_modify AtC rates 
        \section syntax
        fix_modify AtC rates <add | delete> <list_of_fields> \n
        - add | delete (keyword) = add or delete the calculation of rates (time derivatives) for the listed output fields \n
        - fields (keyword) =  \n
        rates can be calculated for all fields listed in \ref man_hardy_fields
 
        \section examples
        <TT> fix_modify AtC rates add temperature velocity stress </TT> \n
        <TT> fix_modify AtC rates delete stress </TT> \n
        \section description
        Requests calculation and ouput of rates (time derivatives) of the fields from the 
        transfer class. For eulerian analysis (see \ref man_atom_element_map ), these rates
        are the partial time derivatives of the nodal fields, not the full (material) time
        derivatives. \n
        \section restrictions
        Must be used with the hardy/field type of AtC fix 
        ( see \ref man_fix_atc )
        \section related
        \section default
        No rates are calculated by default
      */
      else if (strcmp(arg[argIdx],"rates")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"add")==0) { 
          argIdx++;
          FieldName field_name;
          for (int i = argIdx; i < narg; ++i) {
            field_name = string_to_field(arg[i]);
            rateFlags_(field_name) = true; 
          }
          match = true;
        } 
        else if (strcmp(arg[argIdx],"delete")==0) { 
          argIdx++;
          FieldName field_name;
          for (int i = argIdx; i < narg; ++i) {
            field_name = string_to_field(arg[i]);
            rateFlags_(field_name) = false;
          }
          match = true;
        } 
      }


      /*! \page man_pair_interactions fix_modify AtC pair_interactions on|off 
        \section syntax
        fix_modify AtC pair_interactions on|off \n
        fix_modify AtC bond_interactions on|off \n
      
        \section examples
        <TT> fix_modify AtC bond_interactions on </TT> \n
        \section description
        include bonds and/or pairs in the stress and heat flux computations
        \section restrictions
        \section related
        \section default
         pair interactions: on, bond interactions: off
      */
      if (strcmp(arg[argIdx],"pair_interactions")==0) {    // default true
        argIdx++;
        if (strcmp(arg[argIdx],"on")==0) { hasPairs_ = true; }      
        else                             { hasPairs_ = false;}
        match = true;
      }  
      if (strcmp(arg[argIdx],"bond_interactions")==0) {    // default false
        argIdx++;
        if (strcmp(arg[argIdx],"on")==0) { hasBonds_ = true; }      
        else                             { hasBonds_ = false;}
        match = true;
      }  
 
      /*! \page man_hardy_computes fix_modify AtC computes 
        \section syntax
        fix_modify AtC computes <add | delete> [per-atom compute id] <volume | number> \n
        - add | delete (keyword) = add or delete the calculation of an equivalent continuum field
        for the specified per-atom compute as volume or number density quantity \n
        - per-atom compute id =  name/id for per-atom compute, 
        fields can be calculated for all per-atom computes available from LAMMPS \n
        - volume | number (keyword) = field created is a per-unit-volume quantity
        or a per-atom quantity as weighted by kernel functions \n 
 
        \section examples
        <TT> compute virial all stress/atom </TT> \n
        <TT> fix_modify AtC computes add virial volume </TT> \n
        <TT> fix_modify AtC computes delete virial </TT> \n
        \n
        <TT> compute centrosymmetry all centro/atom </TT> \n
        <TT> fix_modify AtC computes add centrosymmetry number </TT> \n
        \section description
        Calculates continuum fields corresponding to specified per-atom computes created by LAMMPS \n
        \section restrictions
        Must be used with the hardy/field type of AtC fix ( see \ref man_fix_atc ) \n
        Per-atom compute must be specified before corresponding continuum field can be requested \n
        \section related
        See manual page for compute 
        \section default
        No defaults exist for this command
      */
      else if (strcmp(arg[argIdx],"computes")==0) {
        argIdx++;
        if (strcmp(arg[argIdx],"add")==0) { 
          argIdx++;
          string tag(arg[argIdx++]);
          int normalization = NO_NORMALIZATION;
          if (narg > argIdx) {
            if      (strcmp(arg[argIdx],"volume")==0) { 
              normalization = VOLUME_NORMALIZATION;
            }
            else if (strcmp(arg[argIdx],"number")==0) { 
              normalization = NUMBER_NORMALIZATION;
            }
            else if (strcmp(arg[argIdx],"mass")==0) { 
              normalization = MASS_NORMALIZATION;
              throw ATC_Error("mass normalized not implemented");
            }
          }
          computes_[tag] = normalization;
          nComputes_++;
          match = true;
        } 
        else if (strcmp(arg[argIdx],"delete")==0) { 
          argIdx++;
          string tag(arg[argIdx]);
          if (computes_.find(tag) != computes_.end()) {
            computes_.erase(tag);
            nComputes_--;
          }
          else {
            throw ATC_Error(tag+" compute is not in list"); 
          }
          match = true;
        } 
      }


      /*! \page man_hardy_dxa_exact_mode fix_modify AtC dxa_exact_mode 
        \section syntax
        fix_modify AtC dxa_exact_mode <optional on | off> \n
        - on | off (keyword) =  use "exact"/serial mode for DXA-based
             calculation of dislocation density, or not \n

        \section examples
        <TT> fix_modify AtC dxa_exact_mode </TT> \n
        <TT> fix_modify AtC dxa_exact_mode on</TT> \n
        <TT> fix_modify AtC dxa_exact_mode off</TT> \n
        \section description
        Overrides normal "exact"/serial mode for DXA code to extract dislocation segments,
        as opposed to an "inexact" mode that's more efficient for parallel computation of
        large systems. \n 
        \section restrictions
        Must be used with the hardy/field type of AtC fix
        ( see \ref man_fix_atc )
        \section related
        \section default
        By default, the DXA "exact"/serial mode is used (i.e. on). \n
      */
      else if (strcmp(arg[argIdx],"dxa_exact_mode")==0) {
        argIdx++;
        dxaExactMode_ = true;
        if (narg > argIdx && strcmp(arg[argIdx],"off")==0) dxaExactMode_ = false;
        match = true; 
      }

      /*! \page man_sample_frequency fix_modify AtC sample_frequency
        \section syntax
        fix_modify AtC sample_frequency [freq]
        - freq (int) : frequency to sample field in number of steps
        \section examples
        <TT> fix_modify AtC sample_frequency 10
        \section description
        Calculates a surface integral of the given field dotted with the
        outward normal of the faces and puts output in the "GLOBALS" file
        \section restrictions
        Must be used with the hardy/field AtC fix ( see \ref man_fix_atc ) 
        and is only relevant when time filters are being used.
        \section related
        \section default
        none
      */
      else if (strcmp(arg[argIdx],"sample_frequency")==0) {
        argIdx++;
        int value = outputFrequency_; // default to output frequency
        if (narg > 1) {
          if (atoi(arg[argIdx]) > 0) value = atoi(arg[argIdx]);
        }
        sampleFrequency_ = value;
        match = true;
      } // end "sample_frequency"

    // no match, call base class parser
    if (!match) {
      match = ATC_Method::modify(narg, arg);
    }

    return match;
  }

  //-------------------------------------------------------------------
  // called at the beginning of a timestep
  void ATC_Transfer::pre_init_integrate()
  {
    ATC_Method::pre_init_integrate();
  }

  //-------------------------------------------------------------------
  // called at the begining of second half timestep
  // REFACTOR move this to post_neighbor
  void ATC_Transfer::pre_final_integrate()
  {
    // update time 
    update_time(); // time uses step if dt = 0


    
    if ( neighborReset_ && sample_now() ) {
      if (! kernelOnTheFly_ ) {
        if (!moleculeIds_.empty()) compute_kernel_matrix_molecule(); //KKM add
      }
      neighborReset_ = false;
    }
  }

  //-------------------------------------------------------------------
  // called at the end of second half timestep
  void ATC_Transfer::post_final_integrate()
  {
    // compute spatially smoothed quantities
    double dt = lammpsInterface_->dt();
    if ( sample_now() ) {
      
      bool needsBond =  (! bondOnTheFly_ ) &&
             (fieldFlags_(STRESS)
           || fieldFlags_(ESHELBY_STRESS)
           || fieldFlags_(HEAT_FLUX));

      if ( needsBond ) {
        if (pairMap_->need_reset()) {
//        ATC::LammpsInterface::instance()->print_msg("Recomputing bond matrix due to atomReset_ value");
          compute_bond_matrix(); 
        }
      }
      time_filter_pre (dt);
      compute_fields();
      time_filter_post(dt);
      lammpsInterface_->computes_addstep(lammpsInterface_->ntimestep()+sampleFrequency_);
    }

    // output
    if ( output_now() && !outputStepZero_ ) output();
    outputStepZero_ = false;

  }

  //-------------------------------------------------------------------
  void ATC_Transfer::compute_bond_matrix(void)
  {
     bondMatrix_->reset();
  }
  //-------------------------------------------------------------------
  void ATC_Transfer::compute_fields(void)
  {
    
    // keep per-atom computes fresh. JAZ and REJ not sure why; 
    // need to confer with JAT. (JAZ, 4/5/12)
    interscaleManager_.lammps_force_reset();

    // (1) direct quantities
    for(int i=0; i < numFields_; ++i) {
      FieldName index = indices_[i];
      if (fieldFlags_(index)) {
        hardyData_[field_to_string(index)].set_quantity() 
          = (outputFields_[index])->quantity();
      }
    }
    if (fieldFlags_(INTERNAL_ENERGY)) 
      compute_internal_energy(hardyData_["internal_energy"].set_quantity());
    if (fieldFlags_(ENERGY))
      compute_energy(hardyData_["energy"].set_quantity());
    if (fieldFlags_(STRESS)) 
      compute_stress(hardyData_["stress"].set_quantity());
    if (fieldFlags_(HEAT_FLUX)) 
      compute_heatflux(hardyData_["heat_flux"].set_quantity());
// molecule data
    if (fieldFlags_(DIPOLE_MOMENT))
      compute_dipole_moment(hardyData_["dipole_moment"].set_quantity()); 
    if (fieldFlags_(QUADRUPOLE_MOMENT))
      compute_quadrupole_moment(hardyData_["quadrupole_moment"].set_quantity());
    if (fieldFlags_(DISLOCATION_DENSITY)) 
      compute_dislocation_density(hardyData_["dislocation_density"].set_quantity());

    // (2) derived quantities
    // compute: gradients
    if (gradFlags_.has_member(true)) {
      for(int index=0; index < NUM_TOTAL_FIELDS; ++index) {
        if (gradFlags_(index)) {
          string field= field_to_string((FieldName) index);
          string grad_field = field + "_gradient";
          if (hardyData_.find(field) == hardyData_.end() ) {
            throw ATC_Error("field " + field + " needs to be defined for gradient");
          }
          gradient_compute(hardyData_[field].quantity(), hardyData_[grad_field].set_quantity());
        }
      }
    }
    // compute: eshelby stress 
    if (fieldFlags_(ESHELBY_STRESS)) {
      {
      compute_eshelby_stress(hardyData_["eshelby_stress"].set_quantity(),
                             hardyData_["internal_energy"].quantity(),
                             hardyData_["stress"].quantity(),
                             hardyData_["displacement_gradient"].quantity());
      }
    }
    if (fieldFlags_(CAUCHY_BORN_ESHELBY_STRESS)) {
      DENS_MAT & H = hardyData_["displacement_gradient"].set_quantity();
      DENS_MAT E(H.nRows(),1);
      ATOMIC_DATA::const_iterator tfield = hardyData_.find("temperature");
      const DENS_MAT *temp = tfield==hardyData_.end() ? NULL : &((tfield->second).quantity());
      //DENS_MAT & T = hardyData_["temperature"];
      //cauchy_born_entropic_energy(H,E,T); E += hardyData_["internal_energy"];
      cauchy_born_energy(H, E, temp);
      E -= nodalRefPotentialEnergy_;

      compute_eshelby_stress(hardyData_["cauchy_born_eshelby_stress"].set_quantity(),
                             E,hardyData_["stress"].quantity(),
                             hardyData_["displacement_gradient"].quantity());
    }
    // compute: cauchy born stress 
    if (fieldFlags_(CAUCHY_BORN_STRESS)) {
      ATOMIC_DATA::const_iterator tfield = hardyData_.find("temperature");
      const DENS_MAT *temp = tfield==hardyData_.end() ? NULL : &((tfield->second).quantity());
      cauchy_born_stress(hardyData_["displacement_gradient"].quantity(),
                         hardyData_["cauchy_born_stress"].set_quantity(), temp);
    }
    // compute: cauchy born energy 
    if (fieldFlags_(CAUCHY_BORN_ENERGY)) {
      ATOMIC_DATA::const_iterator tfield = hardyData_.find("temperature");
      const DENS_MAT *temp = tfield==hardyData_.end() ? NULL : &((tfield->second).quantity());
      cauchy_born_energy(hardyData_["displacement_gradient"].quantity(), 
                         hardyData_["cauchy_born_energy"].set_quantity(), temp);
    }
    // 1st PK transformed to cauchy (lag) or cauchy transformed to 1st PK (eul)
    if (fieldFlags_(TRANSFORMED_STRESS)) {
      compute_transformed_stress(hardyData_["transformed_stress"].set_quantity(),
                                 hardyData_["stress"].quantity(),
                                 hardyData_["displacement_gradient"].quantity());
    }
    if (fieldFlags_(VACANCY_CONCENTRATION)) {
      compute_vacancy_concentration(hardyData_["vacancy_concentration"].set_quantity(),
                                    hardyData_["displacement_gradient"].quantity(),
                                    hardyData_["number_density"].quantity());
    }
    if (fieldFlags_(ELECTRIC_POTENTIAL)) {
      compute_electric_potential(
        hardyData_[field_to_string(ELECTRIC_POTENTIAL)].set_quantity());
    }
    // compute: rotation and/or stretch from deformation gradient 
    if (fieldFlags_(ROTATION) || fieldFlags_(STRETCH)) {
      compute_polar_decomposition(hardyData_["rotation"].set_quantity(),
                                  hardyData_["stretch"].set_quantity(),
                                  hardyData_["displacement_gradient"].quantity());
    }
    // compute: rotation and/or stretch from deformation gradient 
    if (fieldFlags_(CAUCHY_BORN_ELASTIC_DEFORMATION_GRADIENT)) {
      compute_elastic_deformation_gradient2(hardyData_["elastic_deformation_gradient"].set_quantity(),
                                           hardyData_["stress"].quantity(),
                                           hardyData_["displacement_gradient"].quantity());
    }

    // (3) computes
    lammpsInterface_->computes_clearstep();
    map <string,int>::const_iterator iter;
    for (iter = computes_.begin(); iter != computes_.end(); iter++) {
      string tag = iter->first;
      COMPUTE_POINTER cmpt = lammpsInterface_->compute_pointer(tag);
      int projection = iter->second;
      int ncols = lammpsInterface_->compute_ncols_peratom(cmpt);;
      DENS_MAT atomicData(nLocal_,ncols);
      if (ncols == 1) {
        double * atomData = lammpsInterface_->compute_vector_peratom(cmpt);
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          atomicData(i,0) = atomData[atomIdx];
        }
      }
      else {
        double ** atomData = lammpsInterface_->compute_array_peratom(cmpt);
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          for (int k = 0; k < ncols; k++) {
            atomicData(i,k) = atomData[atomIdx][k];
          }
        }
      }
      // REFACTOR -- make dep manage
      if      (projection == NO_NORMALIZATION) {
        project(atomicData,hardyData_[tag].set_quantity());
      }
      else if (projection == VOLUME_NORMALIZATION) {
        project_volume_normalized(atomicData,hardyData_[tag].set_quantity());
      }
      else if (projection == NUMBER_NORMALIZATION) {
        project_count_normalized(atomicData,hardyData_[tag].set_quantity());
      }
      else if (projection == MASS_NORMALIZATION) {
        throw ATC_Error("unimplemented normalization");
      }
      else {
        throw ATC_Error("unimplemented normalization");
      }
#ifdef ESHELBY_VIRIAL
      if (tag == "virial" && fieldFlags_(ESHELBY_STRESS)) {
        if (atomToElementMapType_ == LAGRANGIAN) {
          DENS_MAT tmp = hardyData_[tag].quantity();
          DENS_MAT & myData(hardyData_[tag].set_quantity());
          myData.reset(nNodes_,FieldSizes[STRESS]);
          DENS_MAT F(3,3),FT(3,3),FTINV(3,3),CAUCHY(3,3),PK1(3,3);
          const DENS_MAT& H(hardyData_["displacement_gradient"].quantity());
          for (int k = 0; k < nNodes_; k++ ) {
            vector_to_symm_matrix(k,tmp,CAUCHY);
            vector_to_matrix(k,H,F);
            F(0,0) += 1.0; F(1,1) += 1.0; F(2,2) += 1.0;
            FT = F.transpose();
            FTINV = inv(FT);
            
            //       volumes are already reference volumes.
            PK1 = CAUCHY*FTINV; 
            matrix_to_vector(k,PK1,myData);
          }
        }
        compute_eshelby_stress(hardyData_["eshelby_virial"].set_quantity(),
           hardyData_["internal_energy"].quantity(),hardyData_[tag].quantity(),
           hardyData_["displacement_gradient"].quantity());
      }
#endif
    }
    
  }// end of compute_fields routine

  //-------------------------------------------------------------------
  void ATC_Transfer::time_filter_pre(double dt)
  {
    sampleCounter_++;
    string name;
    double delta_t = dt*sampleFrequency_;
    for(int index=0; index < NUM_TOTAL_FIELDS; ++index) {
      if (fieldFlags_(index)) {
        name = field_to_string((FieldName) index);
        timeFilters_(index)->apply_pre_step1(filteredData_[name].set_quantity(),
                                             hardyData_[name].quantity(), delta_t);
      }
    }
    map <string,int>::const_iterator iter;
    int index = NUM_TOTAL_FIELDS;
    for (iter = computes_.begin(); iter != computes_.end(); iter++) {
      string tag = iter->first;
      timeFilters_(index)->apply_pre_step1(filteredData_[tag].set_quantity(),
                                           hardyData_[tag].quantity(), delta_t);
      index++;
    }
  }

  //-------------------------------------------------------------------
  void ATC_Transfer::time_filter_post(double dt)
  {
    sampleCounter_++;
    string name;
    double delta_t = dt*sampleFrequency_;
    for(int index=0; index < NUM_TOTAL_FIELDS; ++index) {
      if (fieldFlags_(index)) {
        name = field_to_string((FieldName) index);
        timeFilters_(index)->apply_post_step2(filteredData_[name].set_quantity(),
                                              hardyData_[name].quantity(), delta_t);
      }
    }
    if (rateFlags_.has_member(true)) {
      for(int index=0; index < NUM_TOTAL_FIELDS; ++index) {
        if (rateFlags_(index)) {
          string field= field_to_string((FieldName) index);
          string rate_field = field + "_rate";
          timeFilters_(index)->rate(hardyData_[rate_field].set_quantity(),
                                    filteredData_[field].quantity(),
                                    hardyData_[field].quantity(), delta_t);
        }
      }
    }
    for(int index=0; index < NUM_TOTAL_FIELDS; ++index) {
      if (rateFlags_(index)) {
        name = field_to_string((FieldName) index);
        string rate_field = name + "_rate";
        filteredData_[rate_field] = hardyData_[rate_field];
      }
      if (gradFlags_(index)) {
        name = field_to_string((FieldName) index);
        string grad_field = name + "_gradient";
        filteredData_[grad_field] = hardyData_[grad_field];
      }
    }
    
    //       lists/accessing of fields ( & computes)
    map <string,int>::const_iterator iter;
    int index = NUM_TOTAL_FIELDS;
    for (iter = computes_.begin(); iter != computes_.end(); iter++) {
      string tag = iter->first;
      timeFilters_(index)->apply_post_step2(filteredData_[tag].set_quantity(),
                                            hardyData_[tag].quantity(), delta_t);
#ifdef ESHELBY_VIRIAL
      if (tag == "virial" && fieldFlags_(ESHELBY_STRESS)) {
        filteredData_["eshelby_virial"] = hardyData_["eshelby_virial"];
      }
#endif
      index++;
    }
  }

  //-------------------------------------------------------------------
  void ATC_Transfer::output()
  {
    feEngine_->departition_mesh();
    
    for(int index=0; index < NUM_TOTAL_FIELDS; ++index) {
      if (outputFlags_(index)) {
        FieldName fName = (FieldName) index;
        string name= field_to_string(fName);
        fields_[fName] = filteredData_[name];
      }
    }
    
    ATC_Method::output();
    if (lammpsInterface_->comm_rank() == 0) {
      // data
      OUTPUT_LIST output_data;
#ifdef REFERENCE_PE_OUTPUT
      output_data["reference_potential_energy"] = & nodalRefPotentialEnergy_;
#endif
      for(int index=0; index < NUM_TOTAL_FIELDS; ++index) {
        if (outputFlags_(index)) {
          string name= field_to_string((FieldName) index);
          output_data[name]       = & ( filteredData_[name].set_quantity());
        }
        if (rateFlags_(index)) {
          string name= field_to_string((FieldName) index);
          string rate_name = name + "_rate";
          output_data[rate_name] = & ( filteredData_[rate_name].set_quantity());
        }
        if (gradFlags_(index)) {
          string name= field_to_string((FieldName) index);
          string grad_name = name + "_gradient";
          output_data[grad_name] = & ( filteredData_[grad_name].set_quantity());
        }
      }
      
      //       lists/accessing of fields ( & computes)
      map <string,int>::const_iterator iter;
      for (iter = computes_.begin(); iter != computes_.end(); iter++) {
        string tag = iter->first;
        output_data[tag]       = & ( filteredData_[tag].set_quantity());
#ifdef ESHELBY_VIRIAL
        if (tag == "virial" && fieldFlags_(ESHELBY_STRESS)) {
          output_data["eshelby_virial"] = & ( filteredData_["eshelby_virial"].set_quantity() );
        }
#endif
      }

      // output
      feEngine_->write_data(output_index(), & output_data); 
    }
    feEngine_->partition_mesh();
  }

/////// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/////// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //-------------------------------------------------------------------
  // computes nodeData = N*atomData 
  void ATC_Transfer::project(const DENS_MAT & atomData,
                                  DENS_MAT & nodeData)
  {
    if (! kernelOnTheFly_ ) {
      nodeData.reset(nNodes_,atomData.nCols(),true);
      DENS_MAT workNodeArray(nodeData.nRows(),nodeData.nCols());
      if (nLocal_>0) workNodeArray = (accumulant_->quantity()).transMat(atomData);
      int count = nodeData.nRows()*nodeData.nCols();
      lammpsInterface_->allsum(workNodeArray.ptr(),nodeData.ptr(),count);
    }
    else {
      compute_projection(atomData,nodeData);
    }
  }

  //-------------------------------------------------------------------
  // computes nodeData = N*molData specially for molecules
  void ATC_Transfer::project_molecule(const DENS_MAT & molData,
                                  DENS_MAT & nodeData)
  {
    if (! kernelOnTheFly_ ) {
      nodeData.reset(nNodes_,molData.nCols(),true);
      DENS_MAT workNodeArray(nodeData.nRows(),nodeData.nCols());
      if (nLocal_>0) workNodeArray = (accumulantMol_->quantity()).transMat(molData);
      int count = nodeData.nRows()*nodeData.nCols();
      lammpsInterface_->allsum(workNodeArray.ptr(),nodeData.ptr(),count);
    }
    else {
      compute_projection(molData,nodeData);
    }
  }

  //-------------------------------------------------------------------
  // computes nodeData = gradient of N*molData specially for molecules
  void ATC_Transfer::project_molecule_gradient(const DENS_MAT & molData,
                                  DENS_MAT & nodeData)
  {
    if (! kernelOnTheFly_ ) {
      nodeData.reset(nNodes_,molData.nCols(),true);
      DENS_MAT workNodeArray(nodeData.nRows(),nodeData.nCols());
      if (nLocal_>0) workNodeArray = (accumulantMolGrad_->quantity()).transMat(molData);
      int count = nodeData.nRows()*nodeData.nCols();
      lammpsInterface_->allsum(workNodeArray.ptr(),nodeData.ptr(),count);
    }
    else {
      compute_projection(molData,nodeData);
    }
  }

  //-------------------------------------------------------------------
  // count normalized
  void ATC_Transfer::project_count_normalized(const DENS_MAT & atomData,
                                  DENS_MAT & nodeData)
  {
    DENS_MAT tmp; 
    project(atomData,tmp); 
    nodeData = (accumulantWeights_->quantity())*tmp;
  }

  //-------------------------------------------------------------------
  // volume normalized
  void ATC_Transfer::project_volume_normalized(const DENS_MAT & atomData,
                                        DENS_MAT & nodeData)
  {
    DENS_MAT tmp;
    project(atomData,tmp); 
    nodeData = (accumulantInverseVolumes_->quantity())*tmp;
  }

  //-------------------------------------------------------------------
  // volume normalized molecule
  void ATC_Transfer::project_volume_normalized_molecule(const DENS_MAT & molData,
                                        DENS_MAT & nodeData)
  {
    DENS_MAT tmp; 
    project_molecule(molData,tmp); 
    nodeData = (accumulantInverseVolumes_->quantity())*tmp;
  }

  //-------------------------------------------------------------------
  // volume normalized molecule_gradient
  void ATC_Transfer::project_volume_normalized_molecule_gradient(const DENS_MAT & molData,
                                        DENS_MAT & nodeData)
  {
    DENS_MAT tmp; 
    project_molecule_gradient(molData,tmp); 
    nodeData = (accumulantInverseVolumes_->quantity())*tmp;
  }


  //-------------------------------------------------------------------
  void ATC_Transfer::gradient_compute(const DENS_MAT & inNodeData,
                                           DENS_MAT & outNodeData)
  {
    int nrows = inNodeData.nRows();
    int ncols = inNodeData.nCols();
    outNodeData.reset(nrows,ncols*nsd_);
    int index = 0;
    for (int n = 0; n < ncols; n++) { //output v1,1 v1,2 v1,3 ...
      for (int m = 0; m < nsd_; m++) {
        CLON_VEC inData(inNodeData,CLONE_COL,n);
        CLON_VEC outData(outNodeData,CLONE_COL,index);
        outData = (*((gradientMatrix_->quantity())[m]))*inData;
        ++index;
      }
    }
  }



  //-------------------------------------------------------------------
  void ATC_Transfer::compute_force_matrix()
  {
    atomicForceMatrix_ = pairVirial_->quantity();
  }

  //-------------------------------------------------------------------
  // computes "virial" part of heat flux
  // This is correct ONLY for pair potentials. 
  void ATC_Transfer::compute_heat_matrix()
  {
    atomicHeatMatrix_ = pairHeatFlux_->quantity();
  }

  //-------------------------------------------------------------------
  // set xPointer_ to xref or xatom depending on Lagrangian/Eulerian analysis 
  void ATC_Transfer::set_xPointer()
  {
    xPointer_ = xref_;
    if (atomToElementMapType_ == EULERIAN) {
      xPointer_ = lammpsInterface_->xatom();
    }
  }

  //-------------------------------------------------------------------
// SOON TO BE OBSOLETE
  // check consistency of fieldFlags_
  void ATC_Transfer::check_field_dependencies()
  {
    fieldFlags_ = outputFlags_;
    if (fieldFlags_(TRANSFORMED_STRESS))  {
      fieldFlags_(STRESS) = true;
      fieldFlags_(DISPLACEMENT) = true;
    }
    if (fieldFlags_(ESHELBY_STRESS))  {
      fieldFlags_(STRESS) = true;
      fieldFlags_(INTERNAL_ENERGY) = true;
      fieldFlags_(DISPLACEMENT) = true;
    }
    if (fieldFlags_(CAUCHY_BORN_STRESS)
        || fieldFlags_(CAUCHY_BORN_ENERGY) 
        || fieldFlags_(CAUCHY_BORN_ESHELBY_STRESS)
        || fieldFlags_(CAUCHY_BORN_ELASTIC_DEFORMATION_GRADIENT))  {
      if (! (cauchyBornStress_) ) {
        throw ATC_Error("can't compute cauchy-born stress w/o cauchy born model");
      }
    }
    if (fieldFlags_(CAUCHY_BORN_ELASTIC_DEFORMATION_GRADIENT))  {
      fieldFlags_(STRESS) = true;
    }
    if (fieldFlags_(CAUCHY_BORN_STRESS)
        || fieldFlags_(CAUCHY_BORN_ENERGY)) {
      fieldFlags_(TEMPERATURE)  = true;
      fieldFlags_(DISPLACEMENT) = true;
    }
    if (fieldFlags_(CAUCHY_BORN_ESHELBY_STRESS))  {
      fieldFlags_(TEMPERATURE)  = true;
      fieldFlags_(DISPLACEMENT) = true;
      fieldFlags_(STRESS) = true;
    }
    if (fieldFlags_(VACANCY_CONCENTRATION)) {
      fieldFlags_(DISPLACEMENT) = true;
      fieldFlags_(NUMBER_DENSITY) = true;
    }
    if (fieldFlags_(INTERNAL_ENERGY)) {
      fieldFlags_(POTENTIAL_ENERGY) = true;
      fieldFlags_(THERMAL_ENERGY) = true;
    }
    if (fieldFlags_(ENERGY)) {
      fieldFlags_(POTENTIAL_ENERGY) = true;
      fieldFlags_(KINETIC_ENERGY) = true;
    }
    if (fieldFlags_(TEMPERATURE) || fieldFlags_(HEAT_FLUX) ||
        fieldFlags_(KINETIC_ENERGY) || fieldFlags_(THERMAL_ENERGY) || 
        fieldFlags_(ENERGY) || fieldFlags_(INTERNAL_ENERGY) ||
        fieldFlags_(KINETIC_ENERGY) || (fieldFlags_(STRESS) &&
        atomToElementMapType_ == EULERIAN) ) {
      fieldFlags_(VELOCITY) = true;
      fieldFlags_(MASS_DENSITY) = true;
    }

    if (fieldFlags_(VELOCITY)) {
      fieldFlags_(MASS_DENSITY) = true;
      fieldFlags_(MOMENTUM) = true;
    }
    if (fieldFlags_(DISPLACEMENT)) {
      fieldFlags_(MASS_DENSITY) = true;
    }
    if (fieldFlags_(TEMPERATURE) ) {
      fieldFlags_(NUMBER_DENSITY) = true;
    }

    if (fieldFlags_(ROTATION) || 
        fieldFlags_(STRETCH)) {
      fieldFlags_(DISPLACEMENT) = true;
    }
    if (fieldFlags_(ESHELBY_STRESS)
       || fieldFlags_(CAUCHY_BORN_STRESS) 
       || fieldFlags_(CAUCHY_BORN_ENERGY) 
       || fieldFlags_(CAUCHY_BORN_ESHELBY_STRESS) 
       || fieldFlags_(CAUCHY_BORN_ELASTIC_DEFORMATION_GRADIENT) 
       || fieldFlags_(VACANCY_CONCENTRATION)
       || fieldFlags_(ROTATION)
       || fieldFlags_(STRETCH) ) {
        gradFlags_(DISPLACEMENT) = true;
    }

    // check whether single_enable==0 for stress/heat flux calculation
    if (fieldFlags_(STRESS) || fieldFlags_(HEAT_FLUX)) {
      if (lammpsInterface_->single_enable()==0) {
        throw ATC_Error("Calculation of  stress field not possible with selected pair type.");
      }
    }
    
  } 

//============== THIN WRAPPERS ====================================
// OBSOLETE
// HARDY COMPUTES
// ***************UNCONVERTED**************************

  //-------------------------------------------------------------------
  // total energy
  
  void ATC_Transfer::compute_energy(DENS_MAT & energy)
  {
    PerAtomQuantity<double> * atomicPotentialEnergy = interscaleManager_.per_atom_quantity("AtomicPotentialEnergy");
    atomicScalar_=atomicPotentialEnergy->quantity();
    double mvv2e = lammpsInterface_->mvv2e();
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ** vatom = lammpsInterface_->vatom();
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      double ma =  mass ? mass[type[atomIdx]]: rmass[atomIdx];
      ma *= mvv2e; // convert mass to appropriate units
      // compute kinetic energy per atom
      double* v = vatom[atomIdx];
      double atomKE = 0.0;
      for (int k = 0; k < nsd_; k++) { atomKE += v[k]*v[k]; }
      atomKE *= 0.5*ma;
      // add up total energy per atom
      atomicScalar_(i,0) += atomKE;
    }
    project_volume_normalized(atomicScalar_, energy);
  }
  // internal energy
  
  void ATC_Transfer::compute_internal_energy(DENS_MAT & energy)
  {
    PerAtomQuantity<double> * atomicPotentialEnergy = interscaleManager_.per_atom_quantity("AtomicPotentialEnergy");
    PerAtomQuantity<double> * atomicProlongedVelocity = interscaleManager_.per_atom_quantity("ProlongedVelocity");
    atomicScalar_=atomicPotentialEnergy->quantity();
    atomicVector_=atomicProlongedVelocity->quantity();
    double mvv2e = lammpsInterface_->mvv2e();
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ** vatom = lammpsInterface_->vatom();
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      double ma =  mass ? mass[type[atomIdx]]: rmass[atomIdx];
      ma *= mvv2e; // convert mass to appropriate units
      // compute kinetic energy per atom
      double* v = vatom[atomIdx];
      double atomKE = 0.0;
      for (int k = 0; k < nsd_; k++) { 
        atomKE += (v[k]-atomicVector_(i,k))*(v[k]-atomicVector_(i,k)); 
      }
      atomKE *= 0.5*ma;
      // add up total energy per atom
      atomicScalar_(i,0) += atomKE;
    }
    project_volume_normalized(atomicScalar_, energy);
  }
  //-------------------------------------------------------------------
  // MOLECULE
  //-------------------------------------------------------------------
  void ATC_Transfer::compute_dipole_moment(DENS_MAT & dipole_moment)
  {
    const DENS_MAT & molecularVector(dipoleMoment_->quantity());
    project_volume_normalized_molecule(molecularVector,dipole_moment); // KKM add
   //
  }
  //-------------------------------------------------------------------
  void ATC_Transfer::compute_quadrupole_moment(DENS_MAT & quadrupole_moment)
  {
    const DENS_MAT & molecularVector(quadrupoleMoment_->quantity());
    project_volume_normalized_molecule_gradient(molecularVector,quadrupole_moment); // KKM add
   //
  }
  //-------------------------------------------------------------------
  void ATC_Transfer::compute_stress(DENS_MAT & stress)
  {
    // table of bond functions already calculated in initialize function
    // get conversion factor for nktV to p units
    double nktv2p = lammpsInterface_->nktv2p();

    // calculate kinetic energy tensor part of stress for Eulerian analysis
    if (atomToElementMapType_ == EULERIAN && nLocal_>0) {
      compute_kinetic_stress(stress); 
    }
    else {
      // zero stress table for Lagrangian analysis or if nLocal_ = 0
      stress.zero();
    }
    // add-in potential part of stress tensor
    int nrows = stress.nRows();
    int ncols = stress.nCols();
    DENS_MAT local_potential_hardy_stress(nrows,ncols);
    if (nLocal_>0) {
      if (bondOnTheFly_) {
        compute_potential_stress(local_potential_hardy_stress);
      }
      else {
        // compute table of force & position dyad
        compute_force_matrix();
        // calculate force part of stress tensor
        local_potential_hardy_stress = atomicBondMatrix_*atomicForceMatrix_;
        local_potential_hardy_stress *= 0.5; 
      }
    }
    // global summation of potential part of stress tensor
    DENS_MAT potential_hardy_stress(nrows,ncols);
    int count = nrows*ncols;
    lammpsInterface_->allsum(local_potential_hardy_stress.ptr(),
                             potential_hardy_stress.ptr(), count);
    stress += potential_hardy_stress;
    stress = nktv2p*stress;
  }

  //-------------------------------------------------------------------
  // kinetic energy portion of stress
  void ATC_Transfer::compute_kinetic_stress(DENS_MAT& stress)
  {
    const DENS_MAT& density = hardyData_["mass_density"].quantity();
    const DENS_MAT& velocity = hardyData_["velocity"].quantity();

    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ** vatom    = lammpsInterface_->vatom();
    double mvv2e = lammpsInterface_->mvv2e(); // [MV^2]-->[Energy]

    atomicTensor_.reset(nLocal_,6);
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      double ma =  mass ? mass[type[atomIdx]]: rmass[atomIdx];
      ma *= mvv2e; // convert mass to appropriate units
      double* v = vatom[atomIdx];
      atomicTensor_(i,0) -= ma*v[0]*v[0];
      atomicTensor_(i,1) -= ma*v[1]*v[1];
      atomicTensor_(i,2) -= ma*v[2]*v[2];
      atomicTensor_(i,3) -= ma*v[0]*v[1];
      atomicTensor_(i,4) -= ma*v[0]*v[2];
      atomicTensor_(i,5) -= ma*v[1]*v[2];
    }
    project_volume_normalized(atomicTensor_, stress);

    for (int i = 0; i < nNodes_; i++) {
      double rho_i = mvv2e*density(i,0);
      stress(i,0) += rho_i*velocity(i,0)*velocity(i,0);
      stress(i,1) += rho_i*velocity(i,1)*velocity(i,1);
      stress(i,2) += rho_i*velocity(i,2)*velocity(i,2);
      stress(i,3) += rho_i*velocity(i,0)*velocity(i,1);
      stress(i,4) += rho_i*velocity(i,0)*velocity(i,2);
      stress(i,5) += rho_i*velocity(i,1)*velocity(i,2);
    }
  }

  //-------------------------------------------------------------------
  void ATC_Transfer::compute_heatflux(DENS_MAT & flux)
  {
    // calculate kinetic part of heat flux
    if (atomToElementMapType_ == EULERIAN && nLocal_>0) {
      compute_kinetic_heatflux(flux);
    }
    else {
      flux.zero(); // zero stress table for Lagrangian analysis 
    }
    // add potential part of heat flux vector
    int nrows = flux.nRows();
    int ncols = flux.nCols();
    DENS_MAT local_hardy_heat(nrows,ncols);
    if (nLocal_>0) {
      if (bondOnTheFly_) {
        compute_potential_heatflux(local_hardy_heat);
      }
      else {
        // calculate force/potential-derivative part of heat flux
        compute_heat_matrix();
        local_hardy_heat = atomicBondMatrix_*atomicHeatMatrix_;
      }
    }
    // global summation of potential part of heat flux vector
    DENS_MAT hardy_heat(nrows,ncols);
    int count = nrows*ncols;
    lammpsInterface_->allsum(local_hardy_heat.ptr(),
                             hardy_heat.ptr(), count);
    flux += hardy_heat;
  }

  //-------------------------------------------------------------------
  // compute kinetic part of heat flux
  void ATC_Transfer::compute_kinetic_heatflux(DENS_MAT& flux)
  {
    const DENS_MAT& velocity = hardyData_["velocity"].quantity();
    const DENS_MAT& energy = hardyData_["mass_density"].quantity();
    const DENS_MAT& stress = hardyData_["stress"].quantity();

    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ** vatom    = lammpsInterface_->vatom();
    double mvv2e = lammpsInterface_->mvv2e();
    double * atomPE = lammpsInterface_->compute_pe_peratom();
    double atomKE, atomEnergy;
    atomicVector_.reset(nLocal_,3);
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      double ma = mass ? mass[type[atomIdx]]: rmass[atomIdx];
      ma *= mvv2e; // convert mass to appropriate units
      double* v = vatom[atomIdx];
      atomKE = 0.0;
      for (int k = 0; k < nsd_; k++) { atomKE += v[k]*v[k]; }
      atomKE *= 0.5*ma;
      atomEnergy = atomKE + atomPE[atomIdx];
      for (int j = 0; j < nsd_; j++) {
        atomicVector_(i,j) += atomEnergy*v[j];
      }
    }
    project_volume_normalized(atomicVector_,flux);

    // - e^0_I v_I + \sigma^T_I v_I
    for (int i = 0; i < nNodes_; i++) {
      double e_i = energy(i,0);
      flux(i,0) += (e_i + stress(i,0))*velocity(i,0) 
                 + stress(i,3)*velocity(i,1)+ stress(i,4)*velocity(i,2);
      flux(i,1) += (e_i + stress(i,1))*velocity(i,1) 
                 + stress(i,3)*velocity(i,0)+ stress(i,5)*velocity(i,2);
      flux(i,2) += (e_i + stress(i,2))*velocity(i,2) 
                 + stress(i,4)*velocity(i,0)+ stress(i,5)*velocity(i,1);
    }
  }
  //--------------------------------------------------------------------
  void ATC_Transfer::compute_electric_potential(DENS_MAT & phi)
  {
    // Poisson solve with insulating bcs
    const DENS_MAT & rho = (restrictedCharge_->quantity());
    SPAR_MAT K;
    feEngine_->stiffness_matrix(K);
    double permittivity = lammpsInterface_->dielectric(); 
    permittivity *= LammpsInterface::instance()->epsilon0();
    K *= permittivity;
    BC_SET bcs;
    bcs.insert(pair<int,int>(0,0));
    LinearSolver solver(K,bcs);
    CLON_VEC x = column(phi,0);
    CLON_VEC b = column(rho,0);
    solver.solve(x,b);
//x.print("x:phi");
//b.print("b:rho");
    //LinearSolver solver(K,AUTO_SOLVE,true);
  }
  //--------------------------------------------------------------------
  void ATC_Transfer::compute_vacancy_concentration(DENS_MAT & Cv,
    const DENS_MAT & H, const DENS_MAT & rhoN)
  {
    int * type = lammpsInterface_->atom_type();
    DENS_MAT new_rho(nNodes_,1);
    DENS_MAT atomCnt(nLocal_,1);
    double atomic_weight_sum = 0.0;
    double site_weight_sum = 0.0;
    int number_atoms = 0;
    const DIAG_MAT & myAtomicWeights(atomVolume_->quantity());

    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      if (type[atomIdx] != 13) { 
        atomCnt(i,0) = myAtomicWeights(i,i);
        atomic_weight_sum += myAtomicWeights(i,i);
        number_atoms++;
      }
      site_weight_sum += myAtomicWeights(i,i);
    }
    project_volume_normalized(atomCnt, new_rho);
    DENS_MAT F(3,3);
    for (int i = 0; i < nNodes_; i++) {
      if (atomToElementMapType_ == EULERIAN) {
        // for Eulerian analysis: F = (1-H)^{-1}
        DENS_MAT G(3,3);
        vector_to_matrix(i,H,G);
        G *= -1.;
        G(0,0) += 1.0; G(1,1) += 1.0; G(2,2) += 1.0;
        F = inv(G);
      }
      else if (atomToElementMapType_ == LAGRANGIAN) {
        // for Lagrangian analysis: F = (1+H)
        vector_to_matrix(i,H,F);
        F(0,0) += 1.0; F(1,1) += 1.0; F(2,2) += 1.0;
      }
      double J = det(F);
      Cv(i,0) = 1.0 - J*new_rho(i,0);
    }
  }

  //--------------------------------------------------------------------
  void ATC_Transfer::compute_eshelby_stress(DENS_MAT & M,
    const DENS_MAT & E, const DENS_MAT & S, const DENS_MAT & H)
  {
    // eshelby stress:M, energy:E, stress:S, displacement gradient: H
    // eshelby stress = W I - F^T.P = W I - C.S  [energy]
    // symmetric if isotropic S = a_0 I + a_1 C + a_2 C^2
    M.reset(nNodes_,FieldSizes[ESHELBY_STRESS]);
    double nktv2p = lammpsInterface_->nktv2p();
    DENS_MAT P(3,3),F(3,3),FT(3,3),FTP(3,3),ESH(3,3);
    for (int i = 0; i < nNodes_; i++) {
      double W = E(i,0);
      ESH.identity();
      ESH *= W;

      // copy to local
      if (atomToElementMapType_ == LAGRANGIAN) {
        // Stress notation convention:: 0:11 1:12 2:13  3:21 4:22 5:23  6:31 7:32 8:33
        vector_to_matrix(i,S,P);


        vector_to_matrix(i,H,F);
#ifndef H_BASED
        F(0,0) += 1.0; F(1,1) += 1.0; F(2,2) += 1.0;
#endif
        FT = F.transpose(); 
      }
      else if (atomToElementMapType_ == EULERIAN) {
        vector_to_symm_matrix(i,S,P);
        vector_to_matrix(i,H,F);
        FT = F.transpose();
      }
      FTP = (1.0/nktv2p)*FT*P;
      ESH -= FTP;
      if (atomToElementMapType_ == EULERIAN) {
        // For Eulerian analysis, M = F^T*(w-H^T.CauchyStress)
        DENS_MAT Q(3,3);
        Q.identity();
        // Q stores (1-H)
        Q -= FT.transpose();
        DENS_MAT F(3,3);
        F = inv(Q); 
        FT = F.transpose();
        ESH = FT*ESH; 
      }
      // copy to global
      matrix_to_vector(i,ESH,M);
    }
  }
  //---------------------------------------------------------------------------
  // Computes the Cauchy Born stress tensor, T given displacement gradient, H
  // and optional temperature argument (passed by pointer), TEMP
  //---------------------------------------------------------------------------
  void ATC_Transfer::cauchy_born_stress(const DENS_MAT &H, DENS_MAT &T, const DENS_MAT *temp)
  {
    FIELD_MATS uField;      // uField should contain temperature.
    DENS_MAT_VEC tField;
    GRAD_FIELD_MATS hField;
    DENS_MAT_VEC &h = hField[DISPLACEMENT];
    h.assign(nsd_, DENS_MAT(nNodes_,nsd_));
    tField.assign(nsd_, DENS_MAT(nNodes_,nsd_));
    // each row is the CB stress at a node stored in voigt form 
    T.reset(nNodes_,FieldSizes[CAUCHY_BORN_STRESS]);
    const double nktv2p = lammpsInterface_->nktv2p();
    const double fact = -lammpsInterface_->mvv2e()*nktv2p;

    // reshape H (#nodes,9) into h [3](#nodes,3) displacement gradient
    vector_to_dens_mat_vec(H,h);

    // if temperature is provided, then set it
    if (temp) uField[TEMPERATURE] = *temp;

    // Computes the stress at each node.
    cauchyBornStress_->stress(uField, hField, tField);

    // reshapes the stress, T to a (#node,6) DenseMatrix.
    DENS_MAT S(nNodes_,6);
    symm_dens_mat_vec_to_vector(tField,S);
    S *= fact;
    
    // tField/S holds the 2nd P-K stress tensor. Transform to
    // Cauchy for EULERIAN analysis, transform to 1st P-K
    // for LAGRANGIAN analysis.
    DENS_MAT PK2(3,3),G(3,3),F(3,3),FT(3,3),STRESS(3,3);
    for (int i = 0; i < nNodes_; i++) {

      vector_to_symm_matrix(i,S,PK2);

      if (atomToElementMapType_ == EULERIAN) {

        // for Eulerian analysis: F = (1-H)^{-1}
        vector_to_matrix(i,H,G);
        G *= -1.;

        G(0,0) += 1.0; G(1,1) += 1.0; G(2,2) += 1.0;
        F = inv(G);
        FT  = transpose(F);
        double J = det(F);
        STRESS = F*PK2*FT;
        STRESS *= 1/J; 
        symm_matrix_to_vector(i,STRESS,T);
      }
      else{
        // for Lagrangian analysis: F = 1 + H
        vector_to_matrix(i,H,F);

        F(0,0) += 1.0; F(1,1) += 1.0; F(2,2) += 1.0;
        STRESS = F*PK2;
        matrix_to_vector(i,STRESS,T);
      }
 
    }
  }
  //---------------------------------------------------------------------------
  // Computes the Cauchy Born energy density, E given displacement gradient, H
  // and optional temperature argument (passed by pointer), TEMP
  //---------------------------------------------------------------------------
  void ATC_Transfer::cauchy_born_energy(const DENS_MAT &H, DENS_MAT &E, const DENS_MAT *temp)
  {
    FIELD_MATS uField;      // uField should contain temperature.
    GRAD_FIELD_MATS hField;
    DENS_MAT_VEC &h = hField[DISPLACEMENT];
    h.assign(nsd_, DENS_MAT(nNodes_,nsd_));

    // reshape H (#nodes,9) into h [3](#nodes,3) displacement gradient
    vector_to_dens_mat_vec(H,h);

    // if temperature is provided, then set it
    if (temp) uField[TEMPERATURE] = *temp;

    // Computes the free/potential energy at each node.
    cauchyBornStress_->elastic_energy(uField, hField, E);

    // convert back to energy units for  ( ATC coupling uses MLT units)
    double mvv2e = lammpsInterface_->mvv2e(); // [MV^2]-->[Energy]
    E *= mvv2e;

    // for Eulerian analysis, convert energy density to per-unit deformed volume
    if (atomToElementMapType_ == EULERIAN) {
      DENS_MAT G(3,3),F(3,3);
      for (int i = 0; i < nNodes_; i++) {
        // for Eulerian analysis: F = (1-H)^{-1}
        vector_to_matrix(i,H,G);
        G *= -1.;

        G(0,0) += 1.0; G(1,1) += 1.0; G(2,2) += 1.0;
        F = inv(G);
        double J = det(F);
        E(i,0) *= 1/J;
      }
    }

    // subtract zero point energy
    E -= nodalRefPotentialEnergy_;
  }
  //---------------------------------------------------------------------------
  // Computes the M/LH entropic energy density
  //---------------------------------------------------------------------------
  void ATC_Transfer::cauchy_born_entropic_energy(const DENS_MAT &H, DENS_MAT &E, const DENS_MAT &T)
  {
    FIELD_MATS uField;      // uField should contain temperature.
    uField[TEMPERATURE] = T; 
    GRAD_FIELD_MATS hField;
    DENS_MAT_VEC &h = hField[DISPLACEMENT];
    h.assign(nsd_, DENS_MAT(nNodes_,nsd_));

    // reshape H (#nodes,9) into h [3](#nodes,3) displacement gradient
    vector_to_dens_mat_vec(H,h);

    // Computes the free/potential energy at each node.
    cauchyBornStress_->entropic_energy(uField, hField, E);

    // convert back to energy units for  ( ATC coupling uses MLT units)
    double mvv2e = lammpsInterface_->mvv2e(); // [MV^2]-->[Energy]
    E *= mvv2e;

    // for Eulerian analysis, convert energy density to per-unit deformed volume
    if (atomToElementMapType_ == EULERIAN) {
      DENS_MAT G(3,3),F(3,3);
      for (int i = 0; i < nNodes_; i++) {
        // for Eulerian analysis: F = (1-H)^{-1}
        vector_to_matrix(i,H,G);
        G *= -1.;

        G(0,0) += 1.0; G(1,1) += 1.0; G(2,2) += 1.0;
        F = inv(G);
        double J = det(F);
        E(i,0) *= 1/J;
      }
    }

  }
  //--------------------------------------------------------------------
  void ATC_Transfer::compute_transformed_stress(DENS_MAT & stress,
    const DENS_MAT & T, const DENS_MAT & H)
  {
      stress.reset(nNodes_,FieldSizes[TRANSFORMED_STRESS]);
      DENS_MAT S(3,3),FT(3,3),P(3,3);
      for (int i = 0; i < nNodes_; i++) {
        if (atomToElementMapType_ == EULERIAN) {
          vector_to_symm_matrix(i,T,P);
          // for Eulerian analysis: F^T = (1-H)^{-T}
          DENS_MAT G(3,3);
          vector_to_matrix(i,H,G);
          G *= -1.;

          G(0,0) += 1.0; G(1,1) += 1.0; G(2,2) += 1.0;
          FT = inv(G.transpose());
        }
        else{
          vector_to_matrix(i,T,P);
          // for Lagrangian analysis: F^T = (1+H)^T
          DENS_MAT F(3,3);
          vector_to_matrix(i,H,F);

          F(0,0) += 1.0; F(1,1) += 1.0; F(2,2) += 1.0;
          FT = F.transpose(); 
        }
        //
        double J = det(FT);
        FT *= 1/J;
        if (atomToElementMapType_ == EULERIAN) {
          FT = inv(FT); 
        }
        S = P*FT;
        matrix_to_vector(i,S,stress);
      }
  }
  //--------------------------------------------------------------------
  void ATC_Transfer::compute_polar_decomposition(DENS_MAT & rotation,
    DENS_MAT & stretch, const DENS_MAT & H)
  {
    DENS_MAT F(3,3),R(3,3),U(3,3);
    for (int i = 0; i < nNodes_; i++) { 
      vector_to_matrix(i,H,F);
      F(0,0) += 1.0; F(1,1) += 1.0; F(2,2) += 1.0;
      if (atomToElementMapType_ == EULERIAN) { 
        polar_decomposition(F,R,U,false); } // F = V R
      else  {
        polar_decomposition(F,R,U); } // F = R U
      // copy to local
      if ( fieldFlags_(ROTATION) ) {

        matrix_to_vector(i,R,rotation);
      }
      if ( fieldFlags_(STRETCH) ) {
        matrix_to_vector(i,U,stretch);
      }
    }
  }
  //--------------------------------------------------------------------
  void ATC_Transfer::compute_elastic_deformation_gradient(DENS_MAT & Fe,
    const DENS_MAT & P, const DENS_MAT & H)
  
  {
    // calculate Fe for every node
    const double nktv2p = lammpsInterface_->nktv2p();
    const double fact = 1.0/ ( lammpsInterface_->mvv2e()*nktv2p );
    for (int i = 0; i < nNodes_; i++) { 
      DENS_VEC Pv = global_vector_to_vector(i,P);
      Pv *= fact;
      CBElasticTangentOperator tangent(cauchyBornStress_, Pv);
      NonLinearSolver solver(&tangent);
      DENS_VEC Fv = global_vector_to_vector(i,H);  // pass in initial guess
      add_identity_voigt_unsymmetric(Fv);
      solver.solve(Fv);
      vector_to_global_vector(i,Fv,Fe);
    }
  }
  //--------------------------------------------------------------------
  void ATC_Transfer::compute_elastic_deformation_gradient2(DENS_MAT & Fe,
    const DENS_MAT & P, const DENS_MAT & H)
  {
    // calculate Fe for every node
    const double nktv2p = lammpsInterface_->nktv2p();
    const double fact = 1.0/ ( lammpsInterface_->mvv2e()*nktv2p );
    DENS_MAT F(3,3),R(3,3),U(3,3),PP(3,3),S(3,3);
    for (int i = 0; i < nNodes_; i++) { 
      // get F = RU
      vector_to_matrix(i,H,F);
      F(0,0) += 1.0; F(1,1) += 1.0; F(2,2) += 1.0;
      if (atomToElementMapType_ == EULERIAN) { 
        polar_decomposition(F,R,U,false); } // F = V R
      else  {
        polar_decomposition(F,R,U); } // F = R U
      // get S
      vector_to_matrix(i,P,PP);
      //S = PP*transpose(F);
      S = inv(F)*PP;
      
      S += S.transpose(); S *= 0.5; // symmetrize
      DENS_VEC Sv = to_voigt(S);
      Sv *= fact;
      // solve min_U || S - S_CB(U) ||
      CB2ndElasticTangentOperator tangent(cauchyBornStress_, Sv);
      NonLinearSolver solver(&tangent);
      //DENS_VEC Uv = to_voigt_unsymmetric(U);  // pass in initial guess
      DENS_VEC Uv = to_voigt(U);  // pass in initial guess
      //DENS_VEC Uv(6); Uv(0)=1;Uv(1)=1;Uv(2)=1;Uv(3)=0;Uv(4)=0;Uv(5)=0;
      solver.solve(Uv);
      DENS_MAT Ue = from_voigt(Uv);
      DENS_MAT FFe = R*Ue;
      matrix_to_vector(i,FFe,Fe);
    }
  }
// ===========   Analytical solutions ==========================
} // end namespace ATC


