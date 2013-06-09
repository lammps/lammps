// ATC_Transfer headers ?
#include "ATC_TransferHardy.h"
#include "ATC_Error.h"
#include "FE_Engine.h"
#include "ATC_HardyKernel.h"
#include "LammpsInterface.h"
#include "Quadrature.h"

// Other Headers
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <fstream>
#include <exception>

static int sgn(double x) { return (int) ((x>0) - (x<0)); } 
static int rnd(double x) { return (int) (x+sgn(x)*0.5); } 

using namespace std;

namespace ATC {

  ATC_TransferHardy::ATC_TransferHardy(string groupName,
                                       string matParamFile)
    : ATC_Transfer(),
      setRefPE_(false),
      setRefPEvalue_(false),
      xPointer_(NULL),
      kernelFunction_(NULL),
      cauchyBornStress_(NULL),
      bondOnTheFly_(false),
      kernelOnTheFly_(false),
      useAtomicShapeFunctions_(true)
  {
    // assign default "internal" group
    int igroup = lammpsInterface_->find_group(groupName.c_str());
    groupbit_ |= lammpsInterface_->group_bit(igroup);
    igroups_.insert(igroup);

    nTypes_ = lammpsInterface_->ntypes();


    // Defaults
    simTime_ = 0.;
 
    fieldFlags_.reset(NUM_HARDY_FIELDS);
    fieldFlags_ = true;
    gradFlags_.reset(NUM_HARDY_FIELDS);
    gradFlags_ = false;
    rateFlags_.reset(NUM_HARDY_FIELDS);
    rateFlags_ = false;

    fieldSizes_.reset(NUM_HARDY_FIELDS);
    fieldSizes_ = 1;
    fieldSizes_(HARDY_DISPLACEMENT) =3;
    fieldSizes_(HARDY_MOMENTUM) =3;
    fieldSizes_(HARDY_VELOCITY) =3;
    fieldSizes_(HARDY_PROJECTED_VELOCITY) =3;
    fieldSizes_(HARDY_HEAT_FLUX) =3;
    fieldSizes_(HARDY_STRESS) =9; // LAGRANGIAN
    fieldSizes_(HARDY_ESHELBY_STRESS) =9; 
    fieldSizes_(HARDY_CAUCHY_BORN_STRESS) =6; 
    fieldSizes_(HARDY_TRANSFORMED_STRESS) =9; 
    fieldSizes_(HARDY_TYPE_CONCENTRATION) =nTypes_; 

    atomicOutputMask_.insert("displacement");
    atomicOutputMask_.insert("velocity");
    atomicOutputMask_.insert("temperature");
    atomicOutputMask_.insert("potential_energy");
    atomicOutputMask_.insert("force");
    atomicOutputMask_.insert("virial");
    atomicOutputMask_.insert("centrosymmetry");

    set_line_quadrature(line_ngauss,line_xg,line_wg);

    int computeID;
    // create compute for pe/atom
    computeID = lammpsInterface_->atomPE_create();
    if (lammpsInterface_->comm_rank() == 0 ) {
      cout << "atomPE compute created with ID: " << computeID << "\n" << flush;
    }

  }

  //-------------------------------------------------------------------
  ATC_TransferHardy::~ATC_TransferHardy()
  {
    // Deallocate 
    if (kernelFunction_) delete kernelFunction_; 
  }

  //-------------------------------------------------------------------
  // called before the beginning of a "run"
  void ATC_TransferHardy::initialize() 
  {
    // Base class initalizations
    if (!initialized_) {
      feEngine_->initialize();
      nNodes_ = feEngine_->get_nNodes();
      nsd_ = feEngine_->get_nsd();
      // set ref positions
      set_xref();

      // set nlocal & internal to atom map
      reset_nlocal();
      if (useAtomicShapeFunctions_) { 
        try {
          reset_shape_functions();
        }
        catch(bad_alloc&) { 
          if (lammpsInterface_->comm_rank() == 0) {
            cout << " ATC:: insufficient memory for reset_shape_functions() to execute\n";
            throw;
          }
        }
        set_atomic_weights();
      }

      NodeVolumes_.reset(nNodes_,nNodes_);
      invNodeVolumes_.reset(nNodes_,nNodes_);
      feEngine_->compute_lumped_mass_matrix(NodeVolumes_);
      invNodeVolumes_ = inv(NodeVolumes_);
    }

    // positions
    set_xPointer();

    // set up workspaces 
    if (!initialized_) { 

      nNodesGlobal_ = feEngine_->get_feMesh()->get_nNodes();

      // get periodicity and box bounds/lengths
      lammpsInterface_->get_box_periodicity(periodicity[0],
                                            periodicity[1],periodicity[2]);
      lammpsInterface_->get_box_bounds(box_bounds[0][0],box_bounds[1][0],
                                       box_bounds[0][1],box_bounds[1][1],
                                       box_bounds[0][2],box_bounds[1][2]);

      for (int k = 0; k < 3; k++) {
        box_length[k] = box_bounds[1][k] - box_bounds[0][k]; 
      } 

      if (atomToElementMapType_ == EULERIAN) {
        fieldSizes_(HARDY_STRESS) =6;
      }

      // PE compute
      lammpsInterface_->atomPE_init();

      // ground state for PE
      nodalRefPotentialEnergy_.reset(nNodes_,1);

      // copy requested fieldFlags_ settings to outputFlags_
      outputFlags_ = fieldFlags_;

      // perform consistency check on fieldFlags
      check_fieldFlags_consistency();  

      // check whether single_enable==0 for stress/heat flux calculation
      if (fieldFlags_(HARDY_STRESS) || fieldFlags_(HARDY_HEAT_FLUX)) {
        if (lammpsInterface_->single_enable()==0) {
          throw ATC_Error(0,"Calculation of Hardy stress field not possible with selected pair type.");
        }
      }

      // size arrays for requested/required Hardy fields
      for(int index=0; index < NUM_HARDY_FIELDS; ++index) {
        string name;
        int size;
        if (fieldFlags_(index)) {
          hardy_field_to_string(index,name);
          size = fieldSizes_(index);
          hardyData_        [name].reset(nNodes_,size);
          hardyDataOld_     [name].reset(nNodes_,size);
          filteredHardyData_[name].reset(nNodes_,size);
        }
      }

      // compute table of Hardy bond functions for processors that own atoms 
      if ((! bondOnTheFly_) 
         && ( ( fieldFlags_(HARDY_STRESS) 
             || fieldFlags_(HARDY_ESHELBY_STRESS) 
             || fieldFlags_(HARDY_HEAT_FLUX) ) ) ) {
        try {
          compute_bond_matrix(); 
        } 
        catch(bad_alloc&) { 
          if (lammpsInterface_->comm_rank() == 0) {
            cout << "\n ATC:: stress/heat_flux will be computed on-the-fly\n";
            bondOnTheFly_ = true;
          }
        }
      }

      // compute shape function matrix for kernel functions
      if (kernelFunction_ && (! kernelOnTheFly_)) {
        try{ 
          compute_kernel_matrix();
        }
        catch(bad_alloc&) { 
          if (lammpsInterface_->comm_rank() == 0) {
            cout << "\n ATC:: kernel will be computed on-the-fly\n";
            kernelOnTheFly_ = true;
          }
        }
      }

      // compute table of gradient of Hardy functions
      if (fieldFlags_(HARDY_ESHELBY_STRESS)
       || fieldFlags_(HARDY_CAUCHY_BORN_STRESS) 
       || fieldFlags_(HARDY_VACANCY_CONCENTRATION)) {
        gradFlags_(HARDY_DISPLACEMENT) = true;
      }
      gradientTable_ = vector<SPAR_MAT>(nsd_);
      if (gradFlags_.has_member(true)) {
        compute_gradient_matrix(); 
      }

      // construct & initialize filter
      init_filter();
    }
    initialized_ = true;

    // addstep needs to be done _every_ initialize not just first
    lammpsInterface_->atomPE_addstep(lammpsInterface_->ntimestep()+1);

    if (lammpsInterface_->comm_rank() == 0) {
      cout << " ATC:: conversion factor for energy/vol -> stress " 
           << lammpsInterface_->nktv2p() << "\n";
      cout << " ATC:: cutoff radius "  
           << lammpsInterface_->pair_cutoff() << "\n";
    }
  }

  //-------------------------------------------------------------------
  // sets initial values of filtered quantities
  void ATC_TransferHardy::init_filter()
  {
    timeFilters_.reset(NUM_HARDY_FIELDS);
    sampleCounter_ = 0;
    for(int index=0; index < NUM_HARDY_FIELDS; ++index) {
      string name;
      hardy_field_to_string(index,name);
      filteredHardyData_[name]  = 0.0;
      timeFilters_(index) = timeFilterManager_.construct(); 
    }
  }


  //-------------------------------------------------------------------
  // called after the end of a "run"
  void ATC_TransferHardy::finish() 
  {
    // base class
    ATC_Transfer::finish();
  }

  //-------------------------------------------------------------------
  // this is the parser
  bool ATC_TransferHardy::modify(int narg, char **arg)
  {
    bool match = false;
    double val, dval;

    // check to see if it is a transfer class command
    if (strcmp(arg[0],"transfer")==0) {
      /*! \page man_hardy_fields fix_modify AtC transfer fields 
        \section syntax
        fix_modify AtC transfer fields <all | none> \n
        fix_modify AtC transfer fields <add | delete> <list_of_fields> \n
        - all | none (keyword) = output all or no fields  \n
        - add | delete (keyword) = add or delete the listed output fields \n
        - fields (keyword) =  \n
        density : mass per unit volume \n
        displacement : displacement vector \n
        momentum : momentum per unit volume \n
        velocity : defined by momentum divided by density \n
        projected_velocity : simple kernel estimation of atomic velocities \n
        temperature :  temperature derived from the relative atomic kinetic energy (as done by Hardy) \n
        kinetic_temperature : temperature derived from the full kinetic energy  \n
        energy : total energy (potential + kinetic) per unit volume \n
        number_density : simple kernel estimation of number of atoms per unit volume \n
        stress : 
        Cauchy stress tensor for eulerian analysis (atom_element_map), or
        1st Piola-Kirchhoff stress tensor for lagrangian analysis    \n
        transformed_stress : 
        1st Piola-Kirchhoff stress tensor for eulerian analysis (atom_element_map), or 
                             Cauchy stress tensor for lagrangian analysis    \n
        heat_flux : spatial heat flux vector for eulerian, 
        or referential heat flux vector for lagrangian \n
        \section examples
        <TT> fix_modify AtC transfer fields add velocity temperature </TT>
        \section description
        Allows modification of the fields calculated and output by the Hardy
        transfer class. The commands are cumulative, e.g.\n
        <TT> fix_modify AtC transfer fields none </TT> \n 
        followed by \n 
        <TT> fix_modify AtC transfer fields add velocity temperature </TT> \n
        will only output the velocity and temperature fields.
        \section restrictions
        Must be used with the hardy AtC transfer, see \ref man_fix_atc.
        Currently, the stress and heat flux formulas are only correct for 
        central force potentials, e.g. Lennard-Jones and EAM 
        but not Stillinger-Weber.
        \section related
        See \ref man_hardy_gradients , \ref man_hardy_rates  and \ref man_hardy_computes
        \section default
        All fields are output by default
      */
      if (strcmp(arg[1],"fields")==0) {
        if (strcmp(arg[2],"all")==0) { 
          fieldFlags_ = true;
          match = true;
        } 
        else if (strcmp(arg[2],"none")==0) { 
          fieldFlags_ = false;
          match = true;
        } 
        else if (strcmp(arg[2],"add")==0) { 
          hardyFieldName field_name;
          for (int i = 3; i < narg; ++i) {
            if (string_to_hardy_field(arg[i],field_name)) {
              fieldFlags_(field_name) = true; }
            else { throw ATC_Error(0,"unsupported Hardy field"); }
          }
          match = true;
        } 
        else if (strcmp(arg[2],"delete")==0) { 
          hardyFieldName field_name;
          for (int i = 3; i < narg; ++i) {
            if (string_to_hardy_field(arg[i],field_name)) {
              fieldFlags_(field_name) = false; }
            else { throw ATC_Error(0,"unsupported Hardy field"); }
          }
          match = true;
        } 
      }

      /*! \page man_hardy_gradients fix_modify AtC transfer gradients
        \section syntax
        fix_modify AtC transfer gradients <add | delete> <list_of_fields> \n
        - add | delete (keyword) = add or delete the calculation of gradients for the listed output fields \n
        - fields (keyword) =  \n
        gradients can be calculated for all fields listed in \ref man_hardy_fields
 
        \section examples
        <TT> fix_modify AtC transfer gradients add temperature velocity stress </TT> \n
        <TT> fix_modify AtC transfer gradients delete velocity </TT> \n
        \section description
        Requests calculation and ouput of gradients of the fields from the Hardy
        transfer class. These gradients will be with regard to spatial or material
        coordinate for eulerian or lagrangian analysis, respectively, as specified by 
        atom_element_map (see \ref man_atom_element_map )
        \section restrictions
        Must be used with the hardy AtC transfer
        ( see \ref man_fix_atc )
        \section related
        \section default
        No gradients are calculated by default
      */
      else if (strcmp(arg[1],"gradients")==0) {
        if (strcmp(arg[2],"add")==0) { 
          hardyFieldName field_name;
          for (int i = 3; i < narg; ++i) {
            if (string_to_hardy_field(arg[i],field_name)) {
              gradFlags_(field_name) = true; }
            else { throw ATC_Error(0,"unsupported Hardy field"); }
          }
          match = true;
        } 
        else if (strcmp(arg[2],"delete")==0) { 
          hardyFieldName field_name;
          for (int i = 3; i < narg; ++i) {
            if (string_to_hardy_field(arg[i],field_name)) {
              gradFlags_(field_name) = false; }
            else { throw ATC_Error(0,"unsupported Hardy field"); }
          }
          match = true;
        } 
      }

      /*! \page man_hardy_rates fix_modify AtC transfer rates 
        \section syntax
        fix_modify AtC transfer rates <add | delete> <list_of_fields> \n
        - add | delete (keyword) = add or delete the calculation of rates (time derivatives) for the listed output fields \n
        - fields (keyword) =  \n
        rates can be calculated for all fields listed in \ref man_hardy_fields
 
        \section examples
        <TT> fix_modify AtC transfer rates add temperature velocity stress </TT> \n
        <TT> fix_modify AtC transfer rates delete stress </TT> \n
        \section description
        Requests calculation and ouput of rates (time derivatives) of the fields from the Hardy
        transfer class. For eulerian analysis (see \ref man_atom_element_map ), these rates
        are the partial time derivatives of the nodal fields, not the full (material) time
        derivatives. \n
        \section restrictions
        Must be used with the hardy AtC transfer
        ( see \ref man_fix_atc )
        \section related
        \section default
        No rates are calculated by default
      */
      else if (strcmp(arg[1],"rates")==0) {
        if (strcmp(arg[2],"add")==0) { 
          hardyFieldName field_name;
          for (int i = 3; i < narg; ++i) {
            if (string_to_hardy_field(arg[i],field_name)) {
              rateFlags_(field_name) = true; }
            else { throw ATC_Error(0,"unsupported Hardy field"); }
          }
          match = true;
        } 
        else if (strcmp(arg[2],"delete")==0) { 
          hardyFieldName field_name;
          for (int i = 3; i < narg; ++i) {
            if (string_to_hardy_field(arg[i],field_name)) {
              rateFlags_(field_name) = false; }
            else { throw ATC_Error(0,"unsupported Hardy field"); }
          }
          match = true;
        } 
      }

      /*! \page man_hardy_computes fix_modify AtC transfer computes 
        \section syntax
        fix_modify AtC transfer computes <add | delete> [per-atom compute id] <volume | number> \n
        - add | delete (keyword) = add or delete the calculation of an equivalent continuum field
        for the specified per-atom compute as volume or number density quantity \n
        - per-atom compute id =  name/id for per-atom compute, 
        fields can be calculated for all per-atom computes available from LAMMPS \n
        - volume | number (keyword) = field created is a per-unit-volume quantity
        or a per-atom quantity as weighted by kernel functions \n 
 
        \section examples
        <TT> compute virial all stress/atom </TT> \n
        <TT> fix_modify AtC transfer computes add virial volume </TT> \n
        <TT> fix_modify AtC transfer computes delete virial </TT> \n
        \n
        <TT> compute centrosymmetry all centro/atom </TT> \n
        <TT> fix_modify AtC transfer computes add centrosymmetry number </TT> \n
        \section description
        Calculates continuum fields corresponding to specified per-atom computes created by LAMMPS \n
        \section restrictions
        Must be used with the hardy AtC transfer ( see \ref man_fix_atc ) \n
        Per-atom compute must be specified before corresponding continuum field can be requested \n
        \section related
        See manual page for compute 
        \section default
        No defaults exist for this command
      */
      else if (strcmp(arg[1],"computes")==0) {
        if (strcmp(arg[2],"add")==0) { 
          string tag(arg[3]);
          int icompute = lammpsInterface_->find_compute(tag.c_str());
          if (icompute < 0) 
            throw ATC_Error(0,"Could not find compute "+tag);
          int normalization = NO_NORMALIZATION;
          if (narg > 4) {
            if      (strcmp(arg[4],"volume")==0) { 
              normalization = VOLUME_NORMALIZATION;
            }
            else if (strcmp(arg[4],"number")==0) { 
              normalization = NUMBER_NORMALIZATION;
            }
          }
          computes_[tag] = normalization;
          match = true;
        } 
        else if (strcmp(arg[2],"delete")==0) { 
          string tag(arg[3]);
          if (computes_.find(tag) != computes_.end()) {
            computes_.erase(tag);
          }
          else {
            throw ATC_Error(0,tag+" compute is not in list"); 
          }
          match = true;
        } 
      }

      /*! \page man_hardy_kernel fix_modify AtC kernel 
        \section syntax
        fix_modify AtC transfer kernel <type> <parameters>
        - type (keyword) = mesh, step, cell, cubic_cylinder, cubic_sphere, quartic_cylinder, quartic_sphere \n
        - parameters :\n
        mesh = none\n
        step = radius (double) \n
        cell = hx, hy, hz (double) or h (double) \n
        cubic_cylinder = radius (double) \n
        cubic_sphere = radius (double) \n
        quartic_cylinder = radius (double) \n
        quartic_sphere = radius (double) \n
        \section examples
        fix_modify AtC transfer kernel cell 1.0 1.0 1.0
        \section description
      
        \section restrictions
        Must be used with the hardy AtC transfer \n
        For cylinder kernel types, cylindrical axis is assumed to be in z-direction \n
        ( see \ref man_fix_atc )
        \section related
        \section default
        Default to the mesh based kernel
      */
      else if (strcmp(arg[1],"kernel")==0) {
        if (kernelFunction_) delete kernelFunction_;
        if (strcmp(arg[2],"mesh")==0) { 
          match = true;
          useAtomicShapeFunctions_ = true;
        }
        else if (strcmp(arg[2],"step")==0) { 
          double parameters[1] = {atof(arg[3])}; 
          kernelFunction_ = new ATC_HardyKernelStep(1,parameters);
          match = true;
        }
        else if (strcmp(arg[2],"cell")==0) { 
          double parameters[3];
          parameters[0] = parameters[1] = parameters[2] = atof(arg[3]);
          if (narg > 5) {
            for (int i = 1; i < 3; i++) { parameters[i] = atof(arg[3+i]); }
          } 
          kernelFunction_ = new ATC_HardyKernelCell(2,parameters);
          match = true;
        }
        else if (strcmp(arg[2],"cubic_cylinder")==0) { 
          double parameters[1] = {atof(arg[3])}; // cutoff radius
          kernelFunction_ = new ATC_HardyKernelCubicCyl(1,parameters);
          match = true;
        }
        else if (strcmp(arg[2],"cubic_sphere")==0) { 
          double parameters[1] = {atof(arg[3])}; // cutoff radius
          kernelFunction_ = new ATC_HardyKernelCubicSphere(1,parameters);
          match = true;
        }
        else if (strcmp(arg[2],"quartic_cylinder")==0) { 
          double parameters[1] = {atof(arg[3])}; // cutoff radius
          kernelFunction_ = new ATC_HardyKernelQuarticCyl(1,parameters);
          match = true;
        }
        else if (strcmp(arg[2],"quartic_sphere")==0) { 
          double parameters[1] = {atof(arg[3])}; // cutoff radius
          kernelFunction_ = new ATC_HardyKernelQuarticSphere(1,parameters);
          match = true;
        }
      }

      /*! \page man_hardy_on_the_fly fix_modify AtC transfer on_the_fly 
        \section syntax
        fix_modify AtC transfer on_the_fly <bond | kernel> <optional on | off> \n
        - bond | kernel (keyword) = specifies on-the-fly calculation of bond or kernel 
        matrix elements \n
        - on | off (keyword) =  activate or discontinue on-the-fly mode \n
 
        \section examples
        <TT> fix_modify AtC transfer on_the_fly bond on </TT> \n
        <TT> fix_modify AtC transfer on_the_fly kernel </TT> \n
        <TT> fix_modify AtC transfer on_the_fly kernel off </TT> \n
        \section description
        Overrides normal mode of pre-calculating and storing bond pair-to-node and 
        kernel atom-to-node matrices. If activated, will calculate elements of these
        matrices during repeated calls of field computations (i.e. "on-the-fly") and not store them for
        future use.   \n
        on flag is optional - if omitted, on_the_fly will be activated for the specified 
        matrix. Can be deactivated using off flag. \n 
        \section restrictions
        Must be used with the hardy AtC transfer
        ( see \ref man_fix_atc )
        \section related
        \section default
        By default, on-the-fly calculation is not active (i.e. off). However, code does a memory allocation
        check to determine if it can store all needed bond and kernel matrix elements. If this allocation
        fails, on-the-fly is activated. \n
      */
      else if (strcmp(arg[1],"on_the_fly")==0) {
        if (strcmp(arg[2],"bond")==0) {
          bondOnTheFly_ = true;
          if (narg > 3 && strcmp(arg[3],"off")==0) bondOnTheFly_ = false;
        }
        else if (strcmp(arg[2],"kernel")==0) {
          kernelOnTheFly_ = true;
          if (narg > 3 && strcmp(arg[3],"off")==0) kernelOnTheFly_ = false;
        }
        else { throw ATC_Error(0,"unsupported on_the_fly type"); }
        match = true;
      }

      /*! \page man_hardy_set fix_modify AtC set 
        \section syntax
        fix_modify AtC transfer set reference_potential_energy <value>
        - value (double) : optional user specified zero point for PE
        \section examples
        <TT> fix_modify AtC transfer set reference_potential_energy </TT> \n
        <TT> fix_modify AtC transfer set reference_potential_energy -0.05 </TT> \n
        \section description
        Used to set various quantities for the post-processing algorithms.
        Currently it only 
        sets the zero point for the potential energy density using 
        the value provided for all nodes, or from the current 
        configuration of the lattice if no value is provided
        \section restrictions
        Must be used with the hardy AtC transfer
        ( see \ref man_fix_atc )
        \section related
        \section default
        Defaults to lammps zero point i.e. isolated atoms
      */
      else if (strcmp(arg[1],"set")==0) {
        if (strcmp(arg[2],"reference_potential_energy")==0) {
          if (narg > 3) {
            double value = atof(arg[3]);
            nodalRefPEvalue_ = value;
            setRefPEvalue_ = true;
          } 
          else { // NOTE make this part of initialize
            setRefPE_ = true;
          }
          match = true;
        }
      } // end "set"

      /*! \page man_boundary_integral fix_modify AtC boundary_integral 
        \section syntax
        fix_modify AtC transfer boundary_integral [field] faceset [name]
        - field (string) : name of hardy field
        - name (string)  : name of faceset
        \section examples
        <TT> fix_modify AtC transfer boundary_integral stress faceset loop1 </TT> \n
        \section description
        Calculates a surface integral of the given field dotted with the
        outward normal of the faces and puts output in the "GLOBALS" file
        \section restrictions
        Must be used with the hardy AtC transfer
        ( see \ref man_fix_atc )
        \section related
        \section default
        none
      */
      else if (strcmp(arg[1],"boundary_integral")==0) {
        hardyFieldName field;
//      if (string_to_hardy_field(arg[2],field)) { }
//      else { throw ATC_Error(0,"unsupported Hardy field"); }
        if(strcmp(arg[3],"faceset")==0) {
          string name(arg[4]); 
          string field(arg[2]); 
          const set< pair<int,int> > * faceSet
            = & ( (feEngine_->get_feMesh())->get_faceset(name));
          pair <string,string> pair_name(name,field);
          bndyIntegralData_[pair_name] = faceSet;
          match = true;
        }
      } // end "boundary_integral"

      /*! \page man_contour_integral fix_modify AtC contour_integral 
        \section syntax
        fix_modify AtC transfer contour_integral [field] faceset [name] <axis [x | y | z]>
        - field (string) : name of hardy field
        - name (string)  : name of faceset
        - axis (string)  : x or y or z
        \section examples
        <TT> fix_modify AtC transfer contour_integral stress faceset loop1 </TT> \n
        \section description
        Calculates a surface integral of the given field dotted with the
        outward normal of the faces and puts output in the "GLOBALS" file
        \section restrictions
        Must be used with the hardy AtC transfer
        ( see \ref man_fix_atc )
        \section related
        \section default
        none
      */
      else if (strcmp(arg[1],"contour_integral")==0) {
        hardyFieldName field;
//      if (string_to_hardy_field(arg[2],field)) { }
//      else { throw ATC_Error(0,"unsupported Hardy field"); }
        if(strcmp(arg[3],"faceset")==0) {
          string name(arg[4]); 
          string field(arg[2]); 
          const set< pair<int,int> > * faceSet
            = & ( (feEngine_->get_feMesh())->get_faceset(name));
          pair <string,string> pair_name(name,field);
          contourIntegralData_[pair_name] = faceSet;
          match = true;
        }
      } // end "contour_integral"
    } 

    // no match, call base class parser
    if (!match) {
      match = ATC_Transfer::modify(narg, arg);
    }

    return match;
  }

  //-------------------------------------------------------------------
  // called at the beginning of a timestep
  void ATC_TransferHardy::pre_init_integrate()
  {
    // output initial configuration
    if (stepCounter_ == 0 && outputFrequency_ > 0) {
      double dt = lammpsInterface_->dt();
      time_filter_pre (dt);
      compute_fields();
      // initialize filtered data
      for(int index=0; index < NUM_HARDY_FIELDS; ++index) {
        string name;
        hardy_field_to_string(index,name);
        filteredHardyData_[name] = hardyData_[name];
      }
      time_filter_post(dt);
      compute_boundary_integrals();
      output();
    }
  }
  //-------------------------------------------------------------------
  // called at the end of first half of a timestep
  void ATC_TransferHardy::post_init_integrate()
  {
  }

  //-------------------------------------------------------------------
  // called at the begining of second half timestep
  void ATC_TransferHardy::pre_final_integrate()
  {
    // cases to recompute transfer matrices:
    // (1) if atoms have changed processors
    bool atomSwitch = atomSwitch_;
    int localSwitch = atomSwitch_, globalSwitch = 0;
    lammpsInterface_->int_allmax(&localSwitch,&globalSwitch);
    atomSwitch = globalSwitch;
    // (2) if eulerian and we are at the atom_to_element_map reset frequency
    //  which is used to track the convection of atoms  
    bool eulerianReset = (atomToElementMapType_ == EULERIAN)
                  && (stepCounter_ % atomToElementMapFrequency_ ==0);

    bool needsKernel = ( kernelFunction_ && (! kernelOnTheFly_ ) );

    bool needsBond =  (! bondOnTheFly_ ) &&
             (fieldFlags_(HARDY_STRESS)
           || fieldFlags_(HARDY_ESHELBY_STRESS)
           || fieldFlags_(HARDY_HEAT_FLUX));

    // (3) if we suspect neighbor lists changing and are 
    // at the neighbor reset frequency
    int ago = lammpsInterface_->neighbor_ago(); 
    int interval = outputFrequency_;
    if (timeFilterManager_.filter_dynamics() && sampleFrequency_ < interval) 
      { interval = sampleFrequency_; } 
    bool neighborReset = (ago <= interval);
    if (neighborResetFrequency_ > 0) {
      neighborReset = (neighborReset 
                       || (stepCounter_ % neighborResetFrequency_ == 0) ); 
    }

    if (atomSwitch || eulerianReset) {
      reset_nlocal();
      if (needsKernel) {
        compute_kernel_matrix();
      } 
      if (useAtomicShapeFunctions_) {
        reset_shape_functions();
        set_atomic_weights();
      }
      if (needsBond) {
        compute_bond_matrix(); 
      }
      if (gradFlags_.has_member(true)) { 
        compute_gradient_matrix(); 
      }
    }
    else if (neighborReset && needsBond) {
      check_pair_map(); // will recompute bond matrix if necessary
    }

  }

  //-------------------------------------------------------------------
  // called at the end of second half timestep
  void ATC_TransferHardy::post_final_integrate()
  {
    double dt = lammpsInterface_->dt();
    simTime_ += dt;
    if (dt == 0.0) simTime_ = stepCounter_;
    ++stepCounter_;

    bool output_now = ( (outputFrequency_ > 0)
      && ((stepCounter_ % outputFrequency_ == 0)) );
    bool sample_now = ( (output_now) 
      || (timeFilterManager_.filter_dynamics() && (stepCounter_ % sampleFrequency_ == 0)));

    // compute spatially smoothed quantities
    if (sample_now) {
      time_filter_pre (dt, output_now); // NOTE drop output_now
      compute_fields();
      time_filter_post(dt, output_now);
      compute_boundary_integrals();
    }

    // output
    if ( output_now ) output();
    if ( (outputFrequencyAtom_ > 0)
         && ((stepCounter_ % outputFrequencyAtom_ == 0)) )
      atomic_output();

    // "1" should be outputFrequency_ or sampleFrequency_
    lammpsInterface_->atomPE_addstep(lammpsInterface_->ntimestep()+1);
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_fields(void)
  {
    if (setRefPE_) { 
      compute_potential_energy(nodalRefPotentialEnergy_);
      setRefPE_ = false;
    }

    if (setRefPEvalue_) {
      nodalRefPotentialEnergy_ = nodalRefPEvalue_;
      setRefPEvalue_ = false;
    }


// NOTE remove this
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ** v    = lammpsInterface_->vatom();
    double ** x    = lammpsInterface_->xatom();
    int atomIdx;

    bool needs_velocity =  fieldFlags_(HARDY_VELOCITY);
    bool needs_projected_velocity =  fieldFlags_(HARDY_PROJECTED_VELOCITY);
    bool needs_momentum =  fieldFlags_(HARDY_MOMENTUM);
    bool needs_density  =  fieldFlags_(HARDY_DENSITY);
    bool needs_displacement  =  fieldFlags_(HARDY_DISPLACEMENT);
    bool needs_temperature  =  fieldFlags_(HARDY_TEMPERATURE);
    bool needs_kinetic_temperature  =  fieldFlags_(HARDY_KINETIC_TEMPERATURE);
    bool needs_stress  =  fieldFlags_(HARDY_STRESS);
    bool needs_eshelby_stress  =  fieldFlags_(HARDY_ESHELBY_STRESS);
    bool needs_cauchy_born_stress  =  fieldFlags_(HARDY_CAUCHY_BORN_STRESS);
    bool needs_heat  =  fieldFlags_(HARDY_HEAT_FLUX);
    bool needs_energy  =  fieldFlags_(HARDY_ENERGY);
    bool needs_number_density  =  fieldFlags_(HARDY_NUMBER_DENSITY);
    bool needs_transformed_stress = fieldFlags_(HARDY_TRANSFORMED_STRESS);
    bool needs_vacancy_concentration = fieldFlags_(HARDY_VACANCY_CONCENTRATION);
    bool needs_type_concentration = fieldFlags_(HARDY_TYPE_CONCENTRATION);

    // (1) direct quantities
    //compute : kernel norm is a _weighted_ number density
    if (needs_number_density) { 
      compute_number_density(hardyData_["number_density"]);
    }
    // compute: density
    if (needs_density) {
      compute_mass_density(hardyData_["density"]);
    }
    // compute: displacement
    if (needs_displacement) { 
      compute_displacement(hardyData_["displacement"],hardyData_["density"]);
    }
    // compute: momentum
    if (needs_momentum) {
      compute_momentum(hardyData_["momentum"]);
    }
    //compute: projected velocity
    if (needs_projected_velocity) { 
      compute_projected_velocity(hardyData_["projected_velocity"]);
    }
    // compute: velocity
    if (needs_velocity) { 
      compute_velocity(hardyData_["velocity"],
                       hardyData_["density"], hardyData_["momentum"]);
// NOTE make this only for mesh
      compute_variation_velocity(uVariationVelocity_, hardyData_["velocity"]);
    }

    // compute : temperature & (full) kinetic temperature
    if (needs_temperature) { 
      compute_temperature(hardyData_["temperature"]);
    }
    if (needs_kinetic_temperature) { 
      compute_kinetic_temperature(hardyData_["kinetic_temperature"]);
    }
    // compute: stress
    if (needs_stress) {
      compute_stress(hardyData_["stress"]);
    }
    // compute: heat flux 
    if (needs_heat) {
      compute_heatflux(hardyData_["heat_flux"]);
    }
    // compute: energy
    if (needs_energy) {
      compute_total_energy(hardyData_["energy"]);
    }

    // (2) derived quantities
    // compute: gradients
    if (gradFlags_.has_member(true)) {
      for(int index=0; index < NUM_HARDY_FIELDS; ++index) {
        if (gradFlags_(index)) {
          string field;
          hardy_field_to_string(index,field);
          string grad_field = field + "_gradient";
          if (hardyData_.find(field) == hardyData_.end() ) {
            throw ATC_Error(0,"field " + field + " needs to be defined for gradient");
          }
          gradient_compute(hardyData_[field], hardyData_[grad_field]);
        }
      }
    }

    // compute:eshelby_stress
    if (fieldFlags_(HARDY_ESHELBY_STRESS)) {
      if (! (gradFlags_(DISPLACEMENT)) ) {
        throw ATC_Error(0,"displacement_gradient needs to be defined for eshelby stress");
      }
      compute_eshelby_stress(hardyData_["eshelby_stress"],
                             hardyData_["energy"],hardyData_["stress"],
                             hardyData_["displacement_gradient"]);
    }

    // (3) computes
    map <string,int>::const_iterator iter;
    for (iter = computes_.begin(); iter != computes_.end(); iter++) {
      string tag = iter->first;
      int projection = iter->second;
      int ncols = lammpsInterface_->compute_ncols(tag.c_str());;
      DENS_MAT atomicData(nLocal_,ncols);
      if (ncols == 1) {
        double * atomData = lammpsInterface_->compute_scalar_data(tag.c_str());
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          atomicData(i,0) = atomData[atomIdx];
        }
      }
      else {
        double ** atomData = lammpsInterface_->compute_vector_data(tag.c_str());
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          for (int k = 0; k < ncols; k++) {
            atomicData(i,k) = atomData[atomIdx][k];
          }
        }
      }
      if (projection == VOLUME_NORMALIZATION) {
        project_volume_normalized(atomicData, hardyData_[tag]);
      }
      else if (projection == NUMBER_NORMALIZATION) {
        project_count_normalized(atomicData, hardyData_[tag]);
      }
      else {
        ATC_Transfer::restrict_unscaled(atomicData, hardyData_[tag]);
      }
    }
    

  }// end of compute_fields routine

  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_boundary_integrals(void)
  {
    if (! bndyIntegralData_.empty())  {
      map < pair<string,string>, const set< PAIR > * >::const_iterator iter;
      DENS_MAT values;
      for (iter = bndyIntegralData_.begin(); 
           iter !=  bndyIntegralData_.end(); iter++) {
        string bndyName  = (iter->first).first;
        string fieldName = (iter->first).second;
        const set< PAIR > & faceSet = *(iter->second);
        const DENS_MAT & data =  (hardyData_.find(fieldName)->second);
        feEngine_->field_surface_flux(data,faceSet,values);
        for (int i = 0; i < values.nRows() ; ++i ) {
          string name = bndyName + "_" + fieldName + "_" 
                                 + ATC_STRING::tostring(i+1);
          feEngine_->add_global(name, values(i,0));
        }
      }
    }
    if (! contourIntegralData_.empty())  {
      map < pair<string,string>, const set< PAIR > * >::const_iterator iter;
      DENS_MAT values;
      for (iter = contourIntegralData_.begin(); 
           iter !=  contourIntegralData_.end(); iter++) {
        string bndyName  = (iter->first).first;
        string fieldName = (iter->first).second;
        const set< PAIR > & faceSet = *(iter->second);
        const DENS_MAT & data =  (hardyData_.find(fieldName)->second);
        feEngine_->field_surface_flux(data,faceSet,values,true);
        for (int i = 0; i < values.nRows() ; ++i ) {
          string name = "contour_"+bndyName + "_" + fieldName + "_" 
                                 + ATC_STRING::tostring(i+1);
          feEngine_->add_global(name, values(i,0));
        }
      }
    }
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_potential_energy(DENS_MAT & nodalPE)
  {
    DENS_MAT atomicEnergy(nLocal_,1);

    // compute pair energy per atom
    double * atomPE = lammpsInterface_->atomPE_compute();

    // add up potential energy per atom
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      atomicEnergy(i,0) += atomPE[atomIdx];
    }
    project_volume_normalized(atomicEnergy, nodalPE );
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_number_density(DENS_MAT & density)
  {
     DENS_MAT atomCnt(nLocal_,1);
     atomCnt = 1;
     project_volume_normalized(atomCnt, density);
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_mass_density(DENS_MAT & density)
  {
    atomicDensity_.reset(nLocal_,1);
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ma;
    for (int i = 0; i < nLocal_; i++) {
       int atomIdx = internalToAtom_(i);
       if (mass) ma = mass[type[atomIdx]];
       else ma = rmass[atomIdx];
       // density
       atomicDensity_(i,0) = ma;
     }
     project_volume_normalized(atomicDensity_, density);
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_displacement(DENS_MAT & displacement,
                                        const  DENS_MAT & density)
  {
    atomicDisplacement_.reset(nLocal_,nsd_);
    double ** x    = lammpsInterface_->xatom();
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ma;
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      if (mass) ma = mass[type[atomIdx]];
      else ma = rmass[atomIdx];
      for (int j = 0; j < nsd_; j++) {
        // displacement
        atomicDisplacement_(i,j) = x[atomIdx][j] - xref_[atomIdx][j];
        // the cheap version of unwrapping atomic positions
        if ((bool) periodicity[j]) {
          double u = atomicDisplacement_(i,j);
          if (u >=  0.5*box_length[j]) { u -= box_length[j]; }
          if (u <= -0.5*box_length[j]) { u += box_length[j]; }
          atomicDisplacement_(i,j) = u;
        }
        atomicDisplacement_(i,j) *= ma;
      }
    }
    project_volume_normalized(atomicDisplacement_, displacement);
    for (int i = 0; i < nNodes_; i++) {
      double rho_i = density(i,0);
      if (rho_i > 0.0) {
        for (int j = 0; j < nsd_; j++) {
          displacement(i,j) = displacement(i,j)/rho_i;
        }
      }
    }
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_momentum(DENS_MAT & momentum)
  {
    atomicMomentum_.reset(nLocal_,nsd_);
    double ** v    = lammpsInterface_->vatom();
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ma;
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      if (mass) ma = mass[type[atomIdx]];
      else ma = rmass[atomIdx];
      // momentum
      for (int j = 0; j < nsd_; j++) {
        atomicMomentum_(i,j) = ma*(v[atomIdx][j]);
      }
    }
    project_volume_normalized(atomicMomentum_, momentum);
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_projected_velocity(DENS_MAT & velocity)
  {
    atomicVelocity_.reset(nLocal_,nsd_);
    double ** v    = lammpsInterface_->vatom();
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      for (int j = 0; j < nsd_; j++) {
        // velocity
        atomicVelocity_(i,j) += v[atomIdx][j];
      }
    }
    project_count_normalized(atomicVelocity_, velocity);
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_velocity(DENS_MAT & velocity,
                                    const  DENS_MAT & density,
                                    const  DENS_MAT & momentum)
  {
    velocity.reset(nNodes_,nsd_);
    for (int i = 0; i < nNodes_; i++) {
      double rho_i = density(i,0);
      if (rho_i > 0.0) {
        for (int j = 0; j < nsd_; j++) {
          velocity(i,j) = momentum(i,j)/rho_i;
        }
      }
    }
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_variation_velocity(DENS_MAT & velocity,
                                               const DENS_MAT & vI)
  {
    if (nLocal_>0) {
    // interpolate nodal velocities to the atoms
    vbar_.reset(nLocal_,nsd_);
    double ** v    = lammpsInterface_->vatom();
    // use of prolong assumes atom system contained within mesh 
    if (useAtomicShapeFunctions_)  {
      ATC_Transfer::prolong(vI,vbar_);
    }
    // compute and store variation velocities of atoms
    uVariationVelocity_.reset(nLocal_,nsd_);
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      for (int j = 0; j < nsd_; j++) {
        velocity(i,j) = v[atomIdx][j] - vbar_(i,j);
      }
    }
    }
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_temperature(DENS_MAT & temperature)
  {
    atomicTemperature_.reset(nLocal_,1);
    double kB = lammpsInterface_->kBoltzmann();
    double Tcoef = 1./(nsd_*kB);
    double ** v    = lammpsInterface_->vatom();
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ma;
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      if (mass) ma = mass[type[atomIdx]];
      else ma = rmass[atomIdx];
      for (int j = 0; j < nsd_; j++) {
        // variation temperature
        atomicTemperature_(i,0) += Tcoef*ma*(uVariationVelocity_(i,j)*uVariationVelocity_(i,j));
      }
    }
    project_count_normalized(atomicTemperature_, temperature);
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_kinetic_temperature(DENS_MAT & temperature)
  {
    atomicKineticTemperature_.reset(nLocal_,1);
    double kB = lammpsInterface_->kBoltzmann();
    double Tcoef = 1./(nsd_*kB);
    double ** v    = lammpsInterface_->vatom();
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ma;
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      if (mass) ma = mass[type[atomIdx]];
      else ma = rmass[atomIdx];
      for (int j = 0; j < nsd_; j++) {
        // full temperature
        atomicKineticTemperature_(i,0) += Tcoef*ma*(v[atomIdx][j]*v[atomIdx][j]);
      }
    }
    project_count_normalized(atomicKineticTemperature_, temperature);
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_stress(DENS_MAT & stress)
  {
    // table of bond functions already calculated in initialize function
    // get conversion factor for mvv to e units
    double mvv2e = lammpsInterface_->mvv2e();
    // get conversion factor for nktV to p units
    double nktv2p = lammpsInterface_->nktv2p();
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ma;

    // calculate kinetic energy tensor part of stress for Eulerian analysis
    DENS_MAT & vI = uVariationVelocity_;
    if (atomToElementMapType_ == EULERIAN && nLocal_>0) {
      atomicStress_.reset(nLocal_,6);
      for (int i = 0; i < nLocal_; i++) {
        int atomIdx = internalToAtom_(i);
        if (mass) ma = mass[type[atomIdx]];
        else ma = rmass[atomIdx];
        // convert mass to appropriate units
        ma = mvv2e*ma;
        atomicStress_(i,0) -= ma*vI(i,0)*vI(i,0);
        atomicStress_(i,1) -= ma*vI(i,1)*vI(i,1);
        atomicStress_(i,2) -= ma*vI(i,2)*vI(i,2);
        atomicStress_(i,3) -= ma*vI(i,0)*vI(i,1);
        atomicStress_(i,4) -= ma*vI(i,0)*vI(i,2);
        atomicStress_(i,5) -= ma*vI(i,1)*vI(i,2);
      }
      project_volume_normalized(atomicStress_, stress);
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
        local_potential_hardy_stress = atomicBondTable_*atomicForceTable_;
        local_potential_hardy_stress *= 0.5; // NOTE: work around for matrix multiplication
      }
    }
      // global summation of potential part of stress tensor
    DENS_MAT potential_hardy_stress(nrows,ncols);
    int count = nrows*ncols;
    lammpsInterface_->allsum(local_potential_hardy_stress.get_ptr(),
                             potential_hardy_stress.get_ptr(), count);
    stress += potential_hardy_stress;
    stress = nktv2p*stress;
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_heatflux(DENS_MAT & flux)
  {
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ma;
    // calculate kinetic part of heat flux
    if (atomToElementMapType_ == EULERIAN && nLocal_>0) {
      // get conversion factor for mvv to e units
      double mvv2e = lammpsInterface_->mvv2e();
      // compute pair energy per atom
      double * atomPE = lammpsInterface_->atomPE_compute();
      double atomKE, atomEnergy;
      atomicHeat_.reset(nLocal_,3);
      for (int i = 0; i < nLocal_; i++) {
        int atomIdx = internalToAtom_(i);
        if (mass) ma = mass[type[atomIdx]];
        else ma = rmass[atomIdx];
        // convert mass to appropriate units
        ma = mvv2e*ma;
            atomKE = 0.0;
        for (int j = 0; j < nsd_; j++) {
          atomKE += 0.5*ma*(uVariationVelocity_(i,j)*uVariationVelocity_(i,j));
        }
        atomEnergy = atomKE + atomPE[atomIdx];
        for (int j = 0; j < nsd_; j++) {
          atomicHeat_(i,j) += atomEnergy*uVariationVelocity_(i,j);
        }
      }
      project_volume_normalized(atomicHeat_,flux);
    }
    else {
      // zero stress table for Lagrangian analysis or if nLocal_ = 0
      flux.zero();
    }
    // add-in potential part of heat flux vector
    int nrows = flux.nRows();
    int ncols = flux.nCols();
    DENS_MAT local_hardy_heat(nrows,ncols);
    if (nLocal_>0) {
      if (bondOnTheFly_) {
        compute_potential_heatflux(local_hardy_heat);
      }
      else {
        // compute table of heat vectors
        compute_heat_matrix();
        // calculate force/potential-derivative part of heat flux
        local_hardy_heat = atomicBondTable_*atomicHeatTable_;
      }
    }
    // global summation of potential part of heat flux vector
    DENS_MAT hardy_heat(nrows,ncols);
    int count = nrows*ncols;
    lammpsInterface_->allsum(local_hardy_heat.get_ptr(),
                             hardy_heat.get_ptr(), count);
    flux += hardy_heat;
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_total_energy(DENS_MAT & energy)
  {
    atomicEnergy_.reset(nLocal_,1);
    // compute pair energy per atom
    double * atomPE = lammpsInterface_->atomPE_compute();
    // get conversion factor for mvv to e units
    double mvv2e = lammpsInterface_->mvv2e();
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double ma;
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      if (mass) ma = mass[type[atomIdx]];
      else ma = rmass[atomIdx];
      // convert mass to appropriate units
      ma = mvv2e*ma;
      // compute kinetic energy per atom
      double atomKE = 0.0;
      for (int k = 0; k < nsd_; k++) {
        atomKE += 0.5*ma*(uVariationVelocity_(i,k)*uVariationVelocity_(i,k));
      }
      // add up total energy per atom
      atomicEnergy_(i,0) += atomPE[atomIdx] + atomKE;
    }
    project_volume_normalized(atomicEnergy_, energy);
    // subtract zero point energy
    energy -= nodalRefPotentialEnergy_;
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::time_filter_pre(double dt, bool output_now)
  {
    sampleCounter_++;
    string name;
    double delta_t = dt*sampleFrequency_;
    for(int index=0; index < NUM_HARDY_FIELDS; ++index) {
      hardy_field_to_string(index,name);
      timeFilters_(index)->apply_pre_step1(filteredHardyData_[name],
                                    hardyData_[name], delta_t);
    }
    // NOTE add computes_ here
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::time_filter_post(double dt, bool output_now)
  {
    sampleCounter_++;
    string name;
    double delta_t = dt*sampleFrequency_;
    for(int index=0; index < NUM_HARDY_FIELDS; ++index) {
      hardy_field_to_string(index,name);
      timeFilters_(index)->apply_post_step2(filteredHardyData_[name],
                                    hardyData_[name], delta_t);
    }
    if (rateFlags_.has_member(true)) {
      for(int index=0; index < NUM_HARDY_FIELDS; ++index) {
        if (rateFlags_(index)) {
          string field;
          hardy_field_to_string(index,field);
          string rate_field = field + "_rate";
          timeFilters_(index)->rate(hardyData_[rate_field],
                            filteredHardyData_[field],
                            hardyData_[field], delta_t);
        }
      }
    }
    for(int index=0; index < NUM_HARDY_FIELDS; ++index) {
      hardy_field_to_string(index,name);
      if (rateFlags_(index)) {
        string rate_field = name + "_rate";
        filteredHardyData_[rate_field] = hardyData_[rate_field];
      }
      if (gradFlags_(index)) {
        string grad_field = name + "_gradient";
        filteredHardyData_[grad_field] = hardyData_[grad_field];
      }
    }
    // NOTE add computes_ here
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::output()
  {
    if (lammpsInterface_->comm_rank() == 0) {
      // data
      OUTPUT_LIST output_data;
      for(int index=0; index < NUM_HARDY_FIELDS; ++index) {
        string name;
        hardy_field_to_string(index,name);
        if (outputFlags_(index)) {
          output_data[name]       = & ( filteredHardyData_[name]);
        }
        if (rateFlags_(index)) {
          string rate_name = name + "_rate";
          output_data[rate_name] = & ( filteredHardyData_[rate_name]);
        }
        if (gradFlags_(index)) {
          string grad_name = name + "_gradient";
          output_data[grad_name] = & ( filteredHardyData_[grad_name]);
        }
      }
      map <string,int>::const_iterator iter;
      for (iter = computes_.begin(); iter != computes_.end(); iter++) {
        string tag = iter->first;
        output_data[tag]       = & ( hardyData_[tag]); // NOTE should be filtered
      }

      // output
      feEngine_->write_data(simTime_, & output_data); 
    }
  }
  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_kernel_matrix(void)
  {
    if (nLocal_>0) {
      set_xPointer();
      kernelShpFcn_.reset(nNodes_,nLocal_);
      DENS_VEC xI(nsd_),xa(nsd_),xaI(nsd_);
      if (lammpsInterface_->comm_rank() == 0 ) {
        cout << "ATC:: computing kernel matrix " << flush;
      }
      int heartbeatFreq = (nNodes_ <= 10 ? 1 : (int) nNodes_ / 10);
      for (int i = 0; i < nNodes_; i++) {
        if (lammpsInterface_->comm_rank() == 0 && i % heartbeatFreq == 0 ) {
          cout << "." << flush;
        }
        // Hardy point
        xI = (feEngine_->get_feMesh())->nodal_coordinates(i);
        //xI = (feEngine_->get_feMesh())->global_coordinates(i);
        int ii = i;
        //int ii = feEngine_->get_feMesh()->map_global_to_unique(i);
        for (int j = 0; j < nLocal_; j++) {
          // atom location
          int lammps_j = internalToAtom_(j); 
          xa.copy(xPointer_[lammps_j],3);
          // difference vector
          xaI = xa - xI;
          periodicity_correction(xaI);
          // compute kernel value & add it to the matrix
          double val = kernelFunction_->value(xaI);
          if (val > 0) kernelShpFcn_.add(ii,j,val);
        }
      }
      kernelShpFcn_.compress();
      if (lammpsInterface_->comm_rank() == 0) {
        cout << "done\n";
      }
    }
  }

  //-------------------------------------------------------------------
  // count normalized
  void ATC_TransferHardy::project_count_normalized(const DENS_MAT & atomData,
                                  DENS_MAT & nodeData)
  {
    int ncols_atomData = atomData.nCols();
    nodeData.reset(nNodes_,ncols_atomData);
    nodeData.zero();
    // kernel functions
    if (kernelFunction_ ) {
      project_volume_normalized(atomData,nodeData);
      for (int i = 0; i < nNodes_; i++) {
        double numDens = hardyData_["number_density"](i,0);
        if (numDens > 0) { 
          double invNumDens = 1.0/numDens;
          for (int j = 0; j <  ncols_atomData; j++) {
            nodeData(i,j) *= invNumDens;
          }
        } 
        else {
          for (int j = 0; j < ncols_atomData; j++) {
            nodeData(i,j) = 0;
          }
        }
      }
    }
    // mesh-based kernel functions
    else {
      ATC_Transfer::restrict(atomData,nodeData);
    }
  }

  //-------------------------------------------------------------------
  // volume normalized
  void ATC_TransferHardy::project_volume_normalized(const DENS_MAT & atomData,
                                        DENS_MAT & nodeData)
  {
    int ncols_atomData = atomData.nCols();
    DENS_MAT workNodeArray(nNodes_,ncols_atomData);
    workNodeArray.zero();
    nodeData.reset(workNodeArray.nRows(),workNodeArray.nCols());
    nodeData.zero();
    // kernel functions
    if (kernelFunction_ ) {
      if (nLocal_>0) {
        double invVol = kernelFunction_->inv_vol();
        // on the fly calculation
        if (kernelOnTheFly_) {
          set_xPointer();
          DENS_VEC xI(nsd_),xa(nsd_),xaI(nsd_);
          double val;
          for (int i = 0; i < nNodes_; i++) {
            xI = (feEngine_->get_feMesh())->nodal_coordinates(i);
            for (int j = 0; j < nLocal_; j++) {
              int lammps_j = internalToAtom_(j);
              xa.copy(xPointer_[lammps_j],3);
              xaI = xa - xI;
              periodicity_correction(xaI);
              val = kernelFunction_->value(xaI);
              if (val > 0) { 
                for (int k=0; k < ncols_atomData; k++) {
                  workNodeArray(i,k) += val*atomData(j,k);
                }
              } 
            } 
          } 
        }
        // matrix calculation
        else {
          workNodeArray = kernelShpFcn_*atomData;
        }
        workNodeArray *= invVol;
      }
    }
    // mesh-based kernel functions
    else {
      // computes nodeData = N*atomData w/ N: volume normalized shape functions 
      if (nLocal_>0)
        workNodeArray = invNodeVolumes_*shpFcn_.transMat(atomData);
    }  
    // accumulate across processors
    int count = workNodeArray.nRows()*workNodeArray.nCols();
    lammpsInterface_->allsum(workNodeArray.get_ptr(),nodeData.get_ptr(),count);
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::gradient_compute(const DENS_MAT & inNodeData,
                                           DENS_MAT & outNodeData)
  {
    int nrows = inNodeData.nRows();
    int ncols = inNodeData.nCols();
    outNodeData.reset(nrows,ncols*nsd_);
    for (int j = 0; j < nNodes_; j++) {
      for (int i = 0; i < nNodes_; i++) {
        int index = 0;
        for (int n = 0; n < ncols; n++) { //output v1,1 v1,2 v1,3 ...
          for (int m = 0; m < nsd_; m++) {
            outNodeData(j,index) += gradientTable_[m](j,i)*inNodeData(i,n);
            index++;
          }
        }
      }
    }
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_pair_map()
  {
    // neighbor lists for "internal" group 
    // a: a "real" atom in the "internal" group 
    // b: any neighbor to a (it should project into the mesh for mesh-based)
    int inum = lammpsInterface_->neighbor_list_inum();
    if (inum != nLocalTotal_)
      throw ATC_Error(0,"size mismatch in neighbor list");
    int *ilist = lammpsInterface_->neighbor_list_ilist();
    int *numneigh = lammpsInterface_->neighbor_list_numneigh();
    int **firstneigh = lammpsInterface_->neighbor_list_firstneigh();
    int * mask = lammpsInterface_->atom_mask();

    pairMap_.clear(); 
    int pair_index = 0;
    pair< int,int > pair_ij;
    for (int i = 0; i < nLocalTotal_; i++) { // NOTE why nLocal__Total__?
      int lammps_i = ilist[i];
      // filter out atoms not in the internal group
      if (mask[lammps_i] & groupbit_) {
        //int tag_i = (lammpsInterface_->atom_tag())[lammps_i];
        for (int j = 0; j < numneigh[lammps_i]; j++) {
          int lammps_j = firstneigh[lammps_i][j];
          lammps_j &= NEIGHMASK;
          //int tag_j = (lammpsInterface_->atom_tag())[lammps_j];
          //if (lammps_i > lammps_j) continue; // skip a > b pairs
          //if (tag_i > tag_j) continue; // skip a > b pairs
          pair_ij.first  = lammps_i; // alpha 
          pair_ij.second = lammps_j; // beta
          pairMap_[pair_ij] = pair_index;
          pair_index++;
        }
      }
    }
    nPairs_ = pair_index;
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_bond_matrix()
  {
    set_xPointer();
    if (lammpsInterface_->comm_rank() == 0) {
      cout << " ATC:: computing bond matrix " << flush;
    }
    // compute pair map
    compute_pair_map();

    // neighbor lists
    int *numneigh = lammpsInterface_->neighbor_list_numneigh();
    int **firstneigh = lammpsInterface_->neighbor_list_firstneigh();

    Array<bool> latticePeriodicity(3);
    latticePeriodicity(0) = (bool) periodicity[0];
    latticePeriodicity(1) = (bool) periodicity[1];
    latticePeriodicity(2) = (bool) periodicity[2];

 
    double lam1,lam2;
    double bond_value; 
    pair< int,int > pair_jk;
    map< pair< int,int >,int >::iterator pairMapIterator;
    int pair_index;


    atomicBondTable_.reset(nNodes_,nPairs_);

    // process differently for mesh vs translation-invariant kernels
    if (kernelFunction_ ) {
      int heartbeatFreq = (nNodes_ <= 10 ? 1 : (int) nNodes_ / 10);
      //int heartbeatFreq = (nNodesGlobal_ <= 10 ? 1 : (int) nNodesGlobal_ / 10);
      // "normal" kernel functions
      DENS_VEC xa(nsd_),xI(nsd_),xaI(nsd_),xb(nsd_),xbI(nsd_),xba(nsd_);
      double kernel_inv_vol = kernelFunction_->inv_vol();
      for (int i = 0; i < nNodes_; i++) {
      //for (int i = 0; i < nNodesGlobal_; i++) {
        if (lammpsInterface_->comm_rank() == 0 && i % heartbeatFreq == 0 ) {
          cout << "." << flush;
        }
        int ii = i;
        //int ii = feEngine_->get_feMesh()->map_global_to_unique(i);
        // Hardy point
        xI = (feEngine_->get_feMesh())->nodal_coordinates(i);
        //xI = (feEngine_->get_feMesh())->global_coordinates(i);
        for (pairMapIterator = pairMap_.begin(); 
             pairMapIterator != pairMap_.end(); pairMapIterator++){
          int lammps_j = (pairMapIterator->first).first ; 
          int lammps_k = (pairMapIterator->first).second; 
          xa.copy(xPointer_[lammps_j],3);
          xaI = xa - xI;
          periodicity_correction(xaI);
          xb.copy(xPointer_[lammps_k],3);
          xba = xb - xa;
          xbI = xba + xaI;
          // intersect segment with kernel support
          kernelFunction_->bond_intercepts(xaI,xbI,lam1,lam2);
          if (lam1 < lam2) {
            bond_value 
              = kernel_inv_vol*(kernelFunction_->bond(xaI,xbI,lam1,lam2)); 
            pair_index = pairMapIterator->second;
            atomicBondTable_.add(ii,pair_index,bond_value);
          } // if lam1 < lam2
        } // pair map
      }// end nodes loop
      atomicBondTable_.compress();
    }
    else {
      // mesh-based kernel functions
      int heartbeatFreq = (int) pairMap_.size() / 10;
      int i=0;
      int nodes_per_element = feEngine_->get_feMesh()->get_nNodesPerElement();
      Array<int> node_list(nodes_per_element);
      DENS_VEC shp(nodes_per_element);
      DENS_VEC xa(nsd_),xb(nsd_),xab(nsd_),xlambda(nsd_);
      double **xatom = lammpsInterface_->xatom();
      for (pairMapIterator = pairMap_.begin(); 
           pairMapIterator != pairMap_.end(); pairMapIterator++){
        if (lammpsInterface_->comm_rank() == 0 && i++ % heartbeatFreq == 0 ) {
          cout << "." << flush;
        }
        int lammps_j = (pairMapIterator->first).first ; 
        int lammps_k = (pairMapIterator->first).second; 
        pair_index = pairMapIterator->second;
        xa.copy(xPointer_[lammps_j],3);
        xb.copy(xPointer_[lammps_k],3);
        lam1 = 0.0; lam2 = 1.0;
        double del_lambda = 0.5*(lam2 - lam1);
        double avg_lambda = 0.5*(lam2 + lam1);
        xab = xa - xb;
        for (int i = 0; i < line_ngauss; i++) {
          double lambda = del_lambda*line_xg[i] +avg_lambda;
          xlambda = lambda*xab + xb;
          // NOTE doesn't work for the case that fe_region \subset md_region
          int dummyEltID;
          feEngine_->shape_functions(xlambda,shp,dummyEltID,
                                     node_list,latticePeriodicity);
          // accumulate to nodes whose support overlaps the integration point
          for (int I = 0; I < nodes_per_element; I++) {
            // Use factor of 0.5 to account for line integration 
            //   domain of 0 to 1 instead of -1 to 1 
            int inode = node_list(I);
            double val = invNodeVolumes_(inode,inode)*shp(I)*line_wg[i]*0.5;
            atomicBondTable_.add(inode,pair_index,val);
          }
        }
      }
      atomicBondTable_.compress();
    }// end if kernelFunction_
    if (lammpsInterface_->comm_rank() == 0) {
      cout << "done\n";
    }
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_gradient_matrix()
  {
    feEngine_->compute_gradient_matrix(gradientTable_);
  }

  //-------------------------------------------------------------------
  // checks if atoms have changed neighbors
  void ATC_TransferHardy::check_pair_map()
  { 
    // (1) compute fresh temporary pairMap
    // (2) check if tmpPairMap == pairMap_

    // neighbor lists
    int *numneigh = lammpsInterface_->neighbor_list_numneigh();
    int **firstneigh = lammpsInterface_->neighbor_list_firstneigh();

    pair< int,int > pair_ij;
    map< pair< int,int >,int >::iterator pairMapIterator;

    for (int i = 0; i < nLocal_; i++) {
      int lammps_i = internalToAtom_(i);
      for (int k = 0; k < numneigh[lammps_i]; k++) {
        int j = firstneigh[lammps_i][k];
        j &= NEIGHMASK;
        pair_ij.first = lammps_i;
        pair_ij.second = j;
        pairMapIterator = pairMap_.find(pair_ij);
        // if pair is not found in pairMap_, recompute pair map and bond table
        // NOTE need to check if pairMap_ has a stale pair <<<<
        // construct new pair map and use "=="
        if (pairMapIterator == pairMap_.end()) {
          cout << "ATC::check_pair_map has found that pair map and bond table need to be recomputed" << endl;
          compute_bond_matrix();
          return;
        }
      }
    }
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_force_matrix()
  {
    // reset force table
    if (atomToElementMapType_ == LAGRANGIAN) {
      atomicForceTable_.reset(nPairs_,9);
    }
    else if (atomToElementMapType_ == EULERIAN) {
      atomicForceTable_.reset(nPairs_,6);
    }
    // calculate force table
    double **xatom = lammpsInterface_->xatom();
    map< pair< int,int >,int >::iterator pairMapIterator;
    for (pairMapIterator = pairMap_.begin();
         pairMapIterator != pairMap_.end(); pairMapIterator++){
      int lammps_i = (pairMapIterator->first).first ;
      int lammps_j = (pairMapIterator->first).second;
      int pair_index = pairMapIterator->second;
      double * xi = xatom[lammps_i];
      double * xj = xatom[lammps_j];
      double delx = xi[0] - xj[0];
      double dely = xi[1] - xj[1];
      double delz = xi[2] - xj[2];
      double rsq = delx*delx + dely*dely + delz*delz;
      double fforce = 0;
      lammpsInterface_->pair_force(lammps_i,lammps_j,rsq,fforce);
      if (atomToElementMapType_ == LAGRANGIAN) {
        double delX = xref_[lammps_i][0] - xref_[lammps_j][0]; 
        double delY = xref_[lammps_i][1] - xref_[lammps_j][1]; 
        double delZ = xref_[lammps_i][2] - xref_[lammps_j][2]; 
        atomicForceTable_(pair_index,0)=-delx*fforce*delX;
        atomicForceTable_(pair_index,1)=-delx*fforce*delY;
        atomicForceTable_(pair_index,2)=-delx*fforce*delZ;
        atomicForceTable_(pair_index,3)=-dely*fforce*delX;
        atomicForceTable_(pair_index,4)=-dely*fforce*delY;
        atomicForceTable_(pair_index,5)=-dely*fforce*delZ;
        atomicForceTable_(pair_index,6)=-delz*fforce*delX;
        atomicForceTable_(pair_index,7)=-delz*fforce*delY;
        atomicForceTable_(pair_index,8)=-delz*fforce*delZ;
      }
      else if (atomToElementMapType_ == EULERIAN) {
        atomicForceTable_(pair_index,0)=-delx*delx*fforce;
        atomicForceTable_(pair_index,1)=-dely*dely*fforce;
        atomicForceTable_(pair_index,2)=-delz*delz*fforce;
        atomicForceTable_(pair_index,3)=-delx*dely*fforce;
        atomicForceTable_(pair_index,4)=-delx*delz*fforce;
        atomicForceTable_(pair_index,5)=-dely*delz*fforce;
      }
    }
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_heat_matrix()
  {
    // reset heat flux and energy tables
    atomicHeatTable_.reset(nPairs_,3);
    // calculate force table
    map< pair< int,int >,int >::iterator pairMapIterator;
    double ** xatom    = lammpsInterface_->xatom();
    double fforce;
    for (pairMapIterator = pairMap_.begin();
         pairMapIterator != pairMap_.end(); pairMapIterator++){
      int lammps_i = (pairMapIterator->first).first ;
      int lammps_j = (pairMapIterator->first).second;
      int pair_index = pairMapIterator->second;
      double * xi = xatom[lammps_i];
      double * xj = xatom[lammps_j];
      double energy_i = 0.0;
      double delx = xi[0] - xj[0];
      double dely = xi[1] - xj[1];
      double delz = xi[2] - xj[2];
      double rsq = delx*delx + dely*dely + delz*delz;
      double pair_pe=lammpsInterface_->pair_force(lammps_i,lammps_j,rsq,fforce);
      fforce *= 0.5;
      // This term is correct ONLY for pair potentials. 
      //  fforce = - (1/x_{ij})*(d(phi_i)/d(x_{ij}) + d(phi_j)/d(x_{ij}))
      //  What we need for the heat flux calculation is the portion of this that belongs only to 
      //  atom j: - (1/x_{ij})*d(phi_j)/d(x_{ij}).
      int i = lammps_i; // NOTE is this right?
      fforce *= (delx*uVariationVelocity_(i,0) +
                 dely*uVariationVelocity_(i,1) + 
                 delz*uVariationVelocity_(i,2)); 
      // NOTE REJ don't think this is right, where is energy_i from above?
      if (atomToElementMapType_ == LAGRANGIAN) {
        double delX = xref_[lammps_i][0] - xref_[lammps_j][0]; 
        double delY = xref_[lammps_i][1] - xref_[lammps_j][1]; 
        double delZ = xref_[lammps_i][2] - xref_[lammps_j][2]; 
        atomicHeatTable_(pair_index,0)=delX*fforce;
        atomicHeatTable_(pair_index,1)=delY*fforce;
        atomicHeatTable_(pair_index,2)=delZ*fforce;
      }
      else if (atomToElementMapType_ == EULERIAN) {
        atomicHeatTable_(pair_index,0)=delx*fforce;
        atomicHeatTable_(pair_index,1)=dely*fforce;
        atomicHeatTable_(pair_index,2)=delz*fforce;
      }
    }
  }

  //-------------------------------------------------------------------
  // on-the-fly calculation of stress
  void ATC_TransferHardy::compute_potential_stress(DENS_MAT& stress)
  {
    set_xPointer();
    stress.zero();
    // neighbor lists
    int *numneigh = lammpsInterface_->neighbor_list_numneigh();
    int **firstneigh = lammpsInterface_->neighbor_list_firstneigh();
    double ** xatom    = lammpsInterface_->xatom();
    Array<bool> latticePeriodicity(3);
    latticePeriodicity(0) = (bool) periodicity[0];
    latticePeriodicity(1) = (bool) periodicity[1];
    latticePeriodicity(2) = (bool) periodicity[2];
    double lam1,lam2;
    double bond_value; 
    // process differently for mesh vs translation-invariant kernels
    if (lammpsInterface_->comm_rank() == 0) {
      cout << "ATC:: computing potential stress: " << flush;
    }
    int heartbeatFreq = (nNodes_ <= 10 ? 1 : (int) nNodes_ / 10);
    if (kernelFunction_ ) {
      // "normal" kernel functions
      DENS_VEC xa(nsd_),xI(nsd_),xaI(nsd_),xb(nsd_),xbI(nsd_),xba(nsd_);
      double kernel_inv_vol = kernelFunction_->inv_vol();
      for (int i = 0; i < nNodes_; i++) {
      //for (int i = 0; i < nNodesGlobal_; i++) {
        if (lammpsInterface_->comm_rank() == 0 && i % heartbeatFreq == 0 ) {
          cout << "." << flush;
        }
        // Hardy point
        xI = (feEngine_->get_feMesh())->nodal_coordinates(i);
        //xI = (feEngine_->get_feMesh())->global_coordinates(i);
        int inode = i;
        //int inode = feEngine_->get_feMesh()->map_global_to_unique(i);
        for (int j = 0; j < nLocal_; j++) {
          // second (neighbor) atom location
          int lammps_j = internalToAtom_(j); 
          xa.copy(xPointer_[lammps_j],3);
          double * xj = xatom[lammps_j];
          // difference vector
          xaI = xa - xI;
          periodicity_correction(xaI);
          for (int k = 0; k < numneigh[lammps_j]; ++k) {
            int lammps_k = firstneigh[lammps_j][k];
            lammps_k &= NEIGHMASK;
            // first atom location
            xb.copy(xPointer_[lammps_k],3);
            double * xk = xatom[lammps_k];
            // difference vector
            xba = xb - xa;
            xbI = xba + xaI;
            kernelFunction_->bond_intercepts(xaI,xbI,lam1,lam2);
            // compute virial
            if (lam1 < lam2) {
              bond_value 
                = kernel_inv_vol*(kernelFunction_->bond(xaI,xbI,lam1,lam2)); 
              double delx = xatom[lammps_j][0] - xatom[lammps_k][0];
              double dely = xatom[lammps_j][1] - xatom[lammps_k][1];
              double delz = xatom[lammps_j][2] - xatom[lammps_k][2];
              double rsq = delx*delx + dely*dely + delz*delz;
              double fforce = 0;
              lammpsInterface_->pair_force(lammps_j,lammps_k,rsq,fforce);
              fforce *= 0.5; // dbl count
              if (atomToElementMapType_ == LAGRANGIAN) {
                double delX = xref_[lammps_j][0] - xref_[lammps_k][0]; 
                double delY = xref_[lammps_j][1] - xref_[lammps_k][1]; 
                double delZ = xref_[lammps_j][2] - xref_[lammps_k][2]; 
                stress(inode,0) +=-delx*fforce*delX*bond_value;
                stress(inode,1) +=-delx*fforce*delY*bond_value;
                stress(inode,2) +=-delx*fforce*delZ*bond_value;
                stress(inode,3) +=-dely*fforce*delX*bond_value;
                stress(inode,4) +=-dely*fforce*delY*bond_value;
                stress(inode,5) +=-dely*fforce*delZ*bond_value;
                stress(inode,6) +=-delz*fforce*delX*bond_value;
                stress(inode,7) +=-delz*fforce*delY*bond_value;
                stress(inode,8) +=-delz*fforce*delZ*bond_value;
              }
              else { //if (atomToElementMapType_ == EULERIAN) {
                stress(inode,0) +=-delx*delx*fforce*bond_value;
                stress(inode,1) +=-dely*dely*fforce*bond_value;
                stress(inode,2) +=-delz*delz*fforce*bond_value;
                stress(inode,3) +=-delx*dely*fforce*bond_value;
                stress(inode,4) +=-delx*delz*fforce*bond_value;
                stress(inode,5) +=-dely*delz*fforce*bond_value;
              }
            }
          }
        }
      }
    }
    else {
      // mesh-based kernel functions
      int nodes_per_element = feEngine_->get_feMesh()->get_nNodesPerElement();
      Array<int> node_list(nodes_per_element);
      DENS_VEC shp(nodes_per_element);
      DENS_VEC xa(nsd_),xb(nsd_),xab(nsd_),xlambda(nsd_);
      for (int j = 0; j < nLocal_; j++) {
        if (lammpsInterface_->comm_rank() == 0 && j % heartbeatFreq == 0 ) {
          cout << "." << flush;
        }
        // first atom location
        int lammps_j = internalToAtom_(j); 
        xa.copy(xPointer_[lammps_j],3);
        for (int k = 0; k < numneigh[lammps_j]; ++k) {
          int lammps_k = firstneigh[lammps_j][k];
          lammps_k &= NEIGHMASK;
          //if (lammps_k < lammps_j) continue; // full neighbor list
          // second (neighbor) atom location
          xb.copy(xPointer_[lammps_k],3);
          double delx = xatom[lammps_j][0] - xatom[lammps_k][0];
          double dely = xatom[lammps_j][1] - xatom[lammps_k][1];
          double delz = xatom[lammps_j][2] - xatom[lammps_k][2];
          double rsq = delx*delx + dely*dely + delz*delz;
          double fforce = 0;
          lammpsInterface_->pair_force(lammps_j,lammps_k,rsq,fforce);
          fforce *= 0.5; // 1/2 sum_ab = sum_(ab)
          double virial[9];
          if (atomToElementMapType_ == LAGRANGIAN) {
            double delX = xref_[lammps_j][0] - xref_[lammps_k][0]; 
            double delY = xref_[lammps_j][1] - xref_[lammps_k][1]; 
            double delZ = xref_[lammps_j][2] - xref_[lammps_k][2]; 
            virial[0] =-delx*fforce*delX;
            virial[1] =-delx*fforce*delY;
            virial[2] =-delx*fforce*delZ;
            virial[3] =-dely*fforce*delX;
            virial[4] =-dely*fforce*delY;
            virial[5] =-dely*fforce*delZ;
            virial[6] =-delz*fforce*delX;
            virial[7] =-delz*fforce*delY;
            virial[8] =-delz*fforce*delZ;
          }
          else {//if (atomToElementMapType_ == EULERIAN) {
            virial[0] =-delx*delx*fforce;
            virial[1] =-dely*dely*fforce;
            virial[2] =-delz*delz*fforce;
            virial[3] =-delx*dely*fforce;
            virial[4] =-delx*delz*fforce;
            virial[5] =-dely*delz*fforce;
          }
          lam1 = 0.0; lam2 = 1.0;
          double del_lambda = 0.5*(lam2 - lam1);
          double avg_lambda = 0.5*(lam2 + lam1);
          xab = xa - xb;
          for (int i = 0; i < line_ngauss; i++) {
            double lambda = del_lambda*line_xg[i] +avg_lambda;
            xlambda = lambda*xab + xb;
            // NOTE doesn't work for the case that fe_region \subset md_region
            int dummyEltID;
            feEngine_->shape_functions(xlambda,shp,dummyEltID,
                                       node_list,latticePeriodicity);
            // accumulate to nodes whose support overlaps the integration point
            for (int I = 0; I < nodes_per_element; I++) {
              // Use factor of 0.5 to account for line integration 
              //   domain of 0 to 1 instead of -1 to 1 
              int inode = node_list(I);
              double bond_value 
                = invNodeVolumes_(inode,inode)*shp(I)*line_wg[i]*0.5;
              stress(inode,0) += virial[0]*bond_value;
              stress(inode,1) += virial[1]*bond_value;
              stress(inode,2) += virial[2]*bond_value;
              stress(inode,3) += virial[3]*bond_value;
              stress(inode,4) += virial[4]*bond_value;
              stress(inode,5) += virial[5]*bond_value;
              if (atomToElementMapType_ == LAGRANGIAN) {
                stress(inode,6) += virial[6]*bond_value;
                stress(inode,7) += virial[7]*bond_value;
                stress(inode,8) += virial[8]*bond_value;
              }
            }
          }
        }
      }
    }// end if kernelFunction_
    if (lammpsInterface_->comm_rank() == 0) {
      cout << "done\n" << flush;
    }
  }

  //-------------------------------------------------------------------
  // on-the-fly calculation of the heat flux
  void ATC_TransferHardy::compute_potential_heatflux(DENS_MAT& flux)
  {
    set_xPointer();
    flux.zero();
    // neighbor lists
    int *numneigh = lammpsInterface_->neighbor_list_numneigh();
    int **firstneigh = lammpsInterface_->neighbor_list_firstneigh();
    double ** xatom    = lammpsInterface_->xatom();
    Array<bool> latticePeriodicity(3);
    latticePeriodicity(0) = (bool) periodicity[0];
    latticePeriodicity(1) = (bool) periodicity[1];
    latticePeriodicity(2) = (bool) periodicity[2];
    double lam1,lam2;
    double bond_value; 
    // process differently for mesh vs translation-invariant kernels
    if (kernelFunction_ ) {
      // "normal" kernel functions
      DENS_VEC xa(nsd_),xI(nsd_),xaI(nsd_),xb(nsd_),xbI(nsd_),xba(nsd_);
      double kernel_inv_vol = kernelFunction_->inv_vol();
      for (int i = 0; i < nNodes_; i++) {
      //for (int i = 0; i < nNodesGlobal_; i++) {
        int inode = i;
        //int inode = feEngine_->get_feMesh()->map_global_to_unique(i);
        // Hardy point
        xI = (feEngine_->get_feMesh())->nodal_coordinates(i);
        //xI = (feEngine_->get_feMesh())->global_coordinates(i);
        for (int j = 0; j < nLocal_; j++) {
          // second (neighbor) atom location
          int lammps_j = internalToAtom_(j); 
          xa.copy(xPointer_[lammps_j],3);
          double * xj = xatom[lammps_j];
          // difference vector
          xaI = xa - xI;
          periodicity_correction(xaI);
          for (int k = 0; k < numneigh[lammps_j]; ++k) {
            int lammps_k = firstneigh[lammps_j][k];
            lammps_k &= NEIGHMASK;
            // first atom location
            xb.copy(xPointer_[lammps_k],3);
            double * xk = xatom[lammps_k];
            // difference vector
            xba = xb - xa;
            xbI = xba + xaI;
            kernelFunction_->bond_intercepts(xaI,xbI,lam1,lam2);
            // compute virial
            if (lam1 < lam2) {
              bond_value 
                = kernel_inv_vol*(kernelFunction_->bond(xaI,xbI,lam1,lam2)); 
              double delx = xatom[lammps_j][0] - xatom[lammps_k][0];
              double dely = xatom[lammps_j][1] - xatom[lammps_k][1];
              double delz = xatom[lammps_j][2] - xatom[lammps_k][2];
              double rsq = delx*delx + dely*dely + delz*delz;
              double fforce = 0;
              double pair_pe = lammpsInterface_->pair_force(lammps_j,lammps_k,rsq,fforce);
              fforce *= 0.5; // dbl count
              fforce *= (delx*uVariationVelocity_(j,0) +
                         dely*uVariationVelocity_(j,1) +
                         delz*uVariationVelocity_(j,2));
              if (atomToElementMapType_ == LAGRANGIAN) {
                double delX = xref_[lammps_j][0] - xref_[lammps_k][0]; 
                double delY = xref_[lammps_j][1] - xref_[lammps_k][1]; 
                double delZ = xref_[lammps_j][2] - xref_[lammps_k][2]; 
                flux(inode,0) +=fforce*delX*bond_value;
                flux(inode,1) +=fforce*delY*bond_value;
                flux(inode,2) +=fforce*delZ*bond_value;
              }
              else { //if (atomToElementMapType_ == EULERIAN) {
                flux(inode,0) +=fforce*delx*bond_value;
                flux(inode,1) +=fforce*dely*bond_value;
                flux(inode,2) +=fforce*delz*bond_value;
              }
            }
          }
        }
      }
    }
    else {
      // mesh-based kernel functions
      int nodes_per_element = feEngine_->get_feMesh()->get_nNodesPerElement();
      Array<int> node_list(nodes_per_element);
      DENS_VEC shp(nodes_per_element);
      DENS_VEC xa(nsd_),xb(nsd_),xab(nsd_),xlambda(nsd_);
      for (int j = 0; j < nLocal_; j++) {
        // first atom location
        int lammps_j = internalToAtom_(j); 
        xa.copy(xPointer_[lammps_j],3);
        for (int k = 0; k < numneigh[lammps_j]; ++k) {
          int lammps_k = firstneigh[lammps_j][k];
          lammps_k &= NEIGHMASK;
          //if (lammps_k < lammps_j) continue; // full neighbor list
          // second (neighbor) atom location
          xb.copy(xPointer_[lammps_k],3);
          double delx = xatom[lammps_j][0] - xatom[lammps_k][0];
          double dely = xatom[lammps_j][1] - xatom[lammps_k][1];
          double delz = xatom[lammps_j][2] - xatom[lammps_k][2];
          double rsq = delx*delx + dely*dely + delz*delz;
          double fforce = 0;
          double pair_pe = lammpsInterface_->pair_force(lammps_j,lammps_k,rsq,fforce);
          fforce *= 0.5; // 1/2 sum_ab = sum_(ab)
          fforce *= (delx*uVariationVelocity_(j,0) +
                     dely*uVariationVelocity_(j,1) +
                     delz*uVariationVelocity_(j,2));
          double flux_vec[3];
          if (atomToElementMapType_ == LAGRANGIAN) {
            double delX = xref_[lammps_j][0] - xref_[lammps_k][0]; 
            double delY = xref_[lammps_j][1] - xref_[lammps_k][1]; 
            double delZ = xref_[lammps_j][2] - xref_[lammps_k][2]; 
            flux_vec[0] =fforce*delX;
            flux_vec[1] =fforce*delY;
            flux_vec[2] =fforce*delZ;
          }
          else {//if (atomToElementMapType_ == EULERIAN) {
            flux_vec[0] =fforce*delx;
            flux_vec[1] =fforce*dely;
            flux_vec[2] =fforce*delz;
          }
          lam1 = 0.0; lam2 = 1.0;
          double del_lambda = 0.5*(lam2 - lam1);
          double avg_lambda = 0.5*(lam2 + lam1);
          xab = xa - xb;
          for (int i = 0; i < line_ngauss; i++) {
            double lambda = del_lambda*line_xg[i] +avg_lambda;
            xlambda = lambda*xab + xb;
            // NOTE doesn't work for the case that fe_region \subset md_region
            int dummyEltID;
            feEngine_->shape_functions(xlambda,shp,dummyEltID,
                                       node_list,latticePeriodicity);
            // accumulate to nodes whose support overlaps the integration point
            for (int I = 0; I < nodes_per_element; I++) {
              // Use factor of 0.5 to account for line integration 
              //   domain of 0 to 1 instead of -1 to 1 
              int inode = node_list(I);
              double bond_value 
                = invNodeVolumes_(inode,inode)*shp(I)*line_wg[i]*0.5;
              flux(inode,0) += flux_vec[0]*bond_value;
              flux(inode,1) += flux_vec[1]*bond_value;
              flux(inode,2) += flux_vec[2]*bond_value;
            }
          }
        }
      }
    }// end if kernelFunction_
  }

  //-------------------------------------------------------------------
  void ATC_TransferHardy::compute_eshelby_stress(DENS_MAT & M,
    const DENS_MAT & E, const DENS_MAT & S, const DENS_MAT & H)
  {
    // eshelby stess:M, energy:E, stress:S, displacement gradient: H
    // eshelby stress = W I - F^T.P = W I - C.S  [energy]
    // symmetric if isotropic S = a_0 I + a_1 C + a_2 C^2
    M.reset(nNodes_,fieldSizes_(HARDY_ESHELBY_STRESS));
    double nktv2p = lammpsInterface_->nktv2p();
    DENS_MAT P(3,3),FT(3,3),FTP(3,3),ESH(3,3);
    for (int i = 0; i < nNodes_; i++) {
      double W = E(i,0);
      ESH = 0.0;
      ESH(0,0) = W; ESH(1,1) = W; ESH(2,2) = W;
      // NOTE make matrix function to make this easier? see VoigtOp HACK
      // copy to local
      if (atomToElementMapType_ == LAGRANGIAN) {
        // Stress notation convention:: 0:11 1:12 2:13  3:21 4:22 5:23  6:31 7:32 8:33
        P(0,0) = S(i,0); P(0,1) = S(i,1); P(0,2) = S(i,2);
        P(1,0) = S(i,3); P(1,1) = S(i,4); P(1,2) = S(i,5);
        P(2,0) = S(i,6); P(2,1) = S(i,7); P(2,2) = S(i,8);
// NOTE can just remove \int P N dA in post-processing
//#define H_BASED
#ifdef H_BASED
        FT(0,0) = H(i,0); FT(1,0) = H(i,1); FT(2,0) = H(i,2);
        FT(0,1) = H(i,3); FT(1,1) = H(i,4); FT(2,1) = H(i,5);
        FT(0,2) = H(i,6); FT(1,2) = H(i,7); FT(2,2) = H(i,8);
#else
        FT(0,0) = 1+H(i,0); FT(1,0) = H(i,1);   FT(2,0) = H(i,2);
        FT(0,1) = H(i,3);   FT(1,1) = 1+H(i,4); FT(2,1) = H(i,5);
        FT(0,2) = H(i,6);   FT(1,2) = H(i,7);   FT(2,2) = 1+H(i,8);
#endif
      }
      else if (atomToElementMapType_ == EULERIAN) {
        // Stress notation convention:: 0:11 1:22 2:33 3:12=21 4:13=31 5:23=32
        P(0,0) = S(i,0); P(0,1) = S(i,3); P(0,2) = S(i,4);
        P(1,0) = S(i,3); P(1,1) = S(i,1); P(1,2) = S(i,5);
        P(2,0) = S(i,4); P(2,1) = S(i,5); P(2,2) = S(i,2);
        FT(0,0) = H(i,0); FT(1,0) = H(i,1); FT(2,0) = H(i,2);
        FT(0,1) = H(i,3); FT(1,1) = H(i,4); FT(2,1) = H(i,5);
        FT(0,2) = H(i,6); FT(1,2) = H(i,7); FT(2,2) = H(i,8);
      }
      FTP = (1.0/nktv2p)*FT*P;
      ESH -= FTP;
      if (atomToElementMapType_ == EULERIAN) {
        // For Eulerian analysis, M = F^T*(w-H^T.CauchyStress)
        DENS_MAT Q(3,3);
        // Q stores (1-H)
        Q(0,0) = 1; Q(1,1) = 1; Q(2,2) = 1;
        Q(0,1) = 0; Q(1,0) = 0;
        Q(0,2) = 0; Q(2,0) = 0;
        Q(1,2) = 0; Q(2,1) = 0;
        Q -= FT.transpose();
        DENS_MAT F; F.reset(3,3);
        F = inv(Q); 
        FT = F.transpose();
        ESH = FT*ESH; 
      }
      // copy to global
      M(i,0) = ESH(0,0); M(i,1) = ESH(0,1); M(i,2) = ESH(0,2);
      M(i,3) = ESH(1,0); M(i,4) = ESH(1,1); M(i,5) = ESH(1,2);
      M(i,6) = ESH(2,0); M(i,7) = ESH(2,1); M(i,8) = ESH(2,2);
#ifdef EXTENDED_ERROR_CHECKING
      if (atomToElementMapType_ == LAGRANGIAN) {
        cout << i << " W: " << W << "\n";
#ifdef H_BASED
        FTP.print("H^T.P");
#else
        FTP.print("F^T.P");
#endif
      }
      else if (atomToElementMapType_ == EULERIAN) {
        cout << i << " w: " << W << "\n";
        FTP.print("H^T.CauchyStress");
      }
      ESH.print("Eshelby stress");
#endif
    }
  }

  //-------------------------------------------------------------------
  // correction of node-pair distances due to periodicity 
  void ATC_TransferHardy::periodicity_correction(DENS_VEC& xaI)
  {
    for (int m = 0; m < nsd_; m++) {
      xaI(m) -= periodicity[m]*box_length[m]*rnd(xaI(m)/box_length[m]);
    } 
  }

  //-------------------------------------------------------------------
  // set xPointer_ to xref or xatom depending on Lagrangian/Eulerian analysis 
  void ATC_TransferHardy::set_xPointer()
  {
    xPointer_ = xref_;
    if (atomToElementMapType_ == EULERIAN) {
      xPointer_ = lammpsInterface_->xatom();
    }
  }

  //-------------------------------------------------------------------
  // check consistency of fieldFlags_
  void ATC_TransferHardy::check_fieldFlags_consistency()
  {
    if (fieldFlags_(HARDY_TRANSFORMED_STRESS))  {
      fieldFlags_(HARDY_STRESS) = true;
      fieldFlags_(HARDY_DISPLACEMENT) = true;
    }
    if (fieldFlags_(HARDY_ESHELBY_STRESS))  {
      fieldFlags_(HARDY_STRESS) = true;
      fieldFlags_(HARDY_ENERGY) = true;
      fieldFlags_(HARDY_DISPLACEMENT) = true;
    }
    if (fieldFlags_(HARDY_CAUCHY_BORN_STRESS))  {
      if (! (cauchyBornStress_) ) {
        throw ATC_Error(0,"can't compute cauchy-born stress w/o cauchy born model");
      }
      fieldFlags_(HARDY_DISPLACEMENT) = true;
    }
    if (fieldFlags_(HARDY_VACANCY_CONCENTRATION)) {
      fieldFlags_(HARDY_DISPLACEMENT) = true;
      fieldFlags_(HARDY_NUMBER_DENSITY) = true;
    }
    if (fieldFlags_(HARDY_TEMPERATURE) || fieldFlags_(HARDY_HEAT_FLUX) ||
        fieldFlags_(HARDY_ENERGY) || (fieldFlags_(HARDY_STRESS) &&
        atomToElementMapType_ == EULERIAN) ) {
      fieldFlags_(HARDY_VELOCITY) = true;
    }
    if (fieldFlags_(HARDY_VELOCITY)) {
      fieldFlags_(HARDY_DENSITY) = true;
      fieldFlags_(HARDY_MOMENTUM) = true;
    }
    if (fieldFlags_(HARDY_DISPLACEMENT)) {
      fieldFlags_(HARDY_DENSITY) = true;
    }
    if (fieldFlags_(HARDY_TEMPERATURE) || 
        fieldFlags_(HARDY_TYPE_CONCENTRATION)) {
      fieldFlags_(HARDY_NUMBER_DENSITY) = true;
    }
    
  } 

} // end namespace ATC
