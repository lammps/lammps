/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

#include "fix_atc.h"
#include <cstdio>
#include <cstring>
#include <sstream>
#include "fix_nve.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "pointers.h"
#include "comm.h"
#include "group.h"

#include "ATC_Method.h"
#include "ATC_Transfer.h"
#include "ATC_TransferKernel.h"
#include "ATC_TransferPartitionOfUnity.h"
#include "ATC_CouplingEnergy.h"
#include "ATC_CouplingMomentum.h"
#include "ATC_CouplingMass.h"
#include "ATC_CouplingMomentumEnergy.h"
#include "LammpsInterface.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using std::string;

#ifdef LAMMPS_BIGBIG
#error "The USER-ATC package is not compatible with -DLAMMPS_BIGBIG"
#endif

// main page of doxygen documentation
/*! \mainpage AtC : Atom-to-Continuum methods
    fix commands:
    - \ref man_fix_atc (links to all related commands)
*/

/* ------------------------------------------------------------------------- */

FixATC::FixATC(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
  lammps_(lmp), atc_(NULL)
{
  // ID GROUP atc PHYSICSTYPE [PARAMETERFILE]
  if (narg < 4 || narg > 5) lmp->error->all(FLERR,"Illegal fix atc command");

  // Set LAMMPS pointer on LammpsInterface
  ATC::LammpsInterface::instance()->set_lammps(lmp);

  /*! \page man_fix_atc fix atc command
    \section syntax
    fix <fixID> <group> atc <type> <parameter_file>
    - fixID = name of fix
    - group = name of group fix is to be applied
    - type\n
    = thermal : thermal coupling with fields: temperature  \n
    = two_temperature : electron-phonon coupling with field: temperature and electron_temperature  \n
    = hardy :  on-the-fly post-processing using kernel localization functions (see "related" section for possible fields) \n
    = field :  on-the-fly post-processing using mesh-based localization functions (see "related" section for possible fields) \n
    - parameter_file = name of the file with material parameters. \n
    note: Neither hardy nor field requires a parameter file
    \section examples
    <TT> fix AtC internal atc thermal Ar_thermal.dat </TT> \n
    <TT> fix AtC internal atc two_temperature Ar_ttm.mat </TT> \n
    <TT> fix AtC internal atc hardy  </TT> \n
    <TT> fix AtC internal atc field  </TT> \n
    \section description
    This fix is the beginning to creating a coupled FE/MD simulation and/or
    an on-the-fly estimation of continuum fields. The coupled versions of this
    fix do Verlet integration and the /post-processing does not.
    After instantiating this fix, several other fix_modify commands will be
    needed to set up the problem, e.g. define the finite element mesh and
    prescribe initial and boundary conditions.

    The following coupling example is typical, but non-exhaustive:\n

<TT>
     # ... commands to create and initialize the MD system \n

     # initial fix to designate coupling type and group to apply it to \n
     #            tag group          physics     material_file \n
     fix          AtC internal   atc thermal     Ar_thermal.mat\n \n
     # create a uniform 12 x 2 x 2 mesh that covers region contain the group \n
     #                                 nx ny nz region   periodicity \n
     fix_modify   AtC mesh create 12 2  2  mdRegion f p p\n \n
     # specify the control method for the type of coupling \n
     #                         physics         control_type \n
     fix_modify   AtC thermal control flux \n \n
     # specify the initial values for the empirical field "temperature"  \n
     #                                 field       node_group  value \n
     fix_modify   AtC initial temperature all         30.\n \n
     # create an output stream for nodal fields \n
     #                                filename      output_frequency  \n
     fix_modify   AtC output atc_fe_output 100\n \n

     run             1000 \n
</TT>

     likewise for this post-processing example: \n

<TT>
     # ... commands to create and initialize the MD system \n

     # initial fix to designate post-processing and the group to apply it to \n
     # no material file is allowed nor required \n
     fix         AtC internal atc hardy \n \n
     # for hardy fix, specific kernel function (function type and range) to
     # be used as a localization function \n
     fix         AtC kernel quartic_sphere 10.0 \n \n
     # create a uniform 1 x 1 x 1 mesh that covers region contain the group \n
     # with periodicity this effectively creats a system average \n
     fix_modify  AtC mesh create 1 1 1 box p p p \n\n
     # change from default lagrangian map to eulerian \n
     #   refreshed every 100 steps \n
     fix_modify  AtC atom_element_map eulerian 100 \n \n
     # start with no field defined \n
     # add mass density, potential energy density, stress and temperature \n
     fix_modify  AtC fields add density energy stress temperature \n\n
     # create an output stream for nodal fields \n
     #                                filename      output_frequency  \n
     fix_modify  AtC output nvtFE 100 text \n

     run             1000 \n
</TT>

    the mesh's linear interpolation functions can be used as the localization function \n
    by using the field option: \n

<TT>
     fix         AtC internal atc field \n \n
     fix_modify  AtC mesh create 1 1 1 box p p p \n\n
     ... \n\n
</TT>

    Note coupling and post-processing can be combined in the same simulations
    using separate fixes.
    \n
    For detailed exposition of the theory and algorithms please see:\n
    - Wagner, GJ; Jones, RE; Templeton, JA; Parks, MA,  <VAR> An
      atomistic-to-continuum coupling method for heat transfer in solids. </VAR>
      Special Issue of Computer Methods and Applied Mechanics (2008) 197:3351. \n
    - Zimmerman, JA; Webb, EB; Hoyt, JJ;. Jones, RE; Klein, PA; Bammann, DJ,
      <VAR> Calculation of stress in atomistic simulation. </VAR>
      Special Issue of Modelling and Simulation in Materials Science and
      Engineering (2004), 12:S319. \n
    - Zimmerman, JA; Jones, RE; Templeton, JA,
      <VAR> A material frame approach for evaluating continuum variables in
       atomistic simulations. </VAR>
      Journal of Computational Physics (2010), 229:2364. \n
    - Templeton, JA; Jones, RE; Wagner, GJ,  <VAR> Application of a field-based method
      to spatially varying thermal transport problems in molecular dynamics. </VAR>
      Modelling and Simulation in Materials Science and Engineering (2010), 18:085007. \n
    - Jones, RE; Templeton, JA; Wagner, GJ; Olmsted, D; Modine, JA, <VAR>
      Electron transport enhanced molecular dynamics for metals and semi-metals.  </VAR>
      International Journal for Numerical Methods in Engineering (2010), 83:940. \n
    - Templeton, JA; Jones, RE; Lee, JW; Zimmerman, JA; Wong, BM,
      <VAR> A long-range electric field solver for molecular dynamics based on
      atomistic-to-continuum modeling.  </VAR>
      Journal of Chemical Theory and Computation (2011), 7:1736. \n
    - Mandadapu, KK; Templeton, JA; Lee, JW, <VAR> Polarization as a field variable
      from molecular dynamics simulations. </VAR>
      Journal of Chemical Physics (2013), 139:054115. \n

    Please refer to the standard
    finite element (FE) texts, e.g. T.J.R Hughes <VAR> The finite element
    method </VAR>, Dover 2003, for the basics of FE simulation.

    \section restrictions
    Thermal and two_temperature (coupling) types use a Verlet time-integration
    algorithm.
    The hardy type does not contain its own time-integrator and must be used
    with a separate fix that does contain one, e.g. nve, nvt, etc.

    Currently,
    - the coupling is restricted to thermal physics
    - the FE computations are done in serial on each processor.

    \section related
    fix_modify commands for setup: \n
    - \ref man_mesh_create
    - \ref man_mesh_quadrature
    - \ref man_mesh_read
    - \ref man_mesh_write
    - \ref man_mesh_create_nodeset
    - \ref man_mesh_add_to_nodeset
    - \ref man_mesh_create_faceset_box
    - \ref man_mesh_create_faceset_plane
    - \ref man_mesh_create_elementset
    - \ref man_mesh_delete_elements
    - \ref man_mesh_nodeset_to_elementset
    - \ref man_boundary
    - \ref man_internal_quadrature
    - \ref man_thermal_time_integration
    - \ref man_momentum_time_integration
    - \ref man_electron_integration
    - \ref man_internal_element_set
    - \ref man_decomposition

    fix_modify commands for boundary and initial conditions:\n
    - \ref man_initial
    - \ref man_fix_nodes
    - \ref man_unfix_nodes
    - \ref man_fix_flux
    - \ref man_unfix_flux
    - \ref man_source
    - \ref man_remove_source

    fix_modify commands for control and filtering: \n
    - \ref man_control
    - \ref man_control_thermal
    - \ref man_control_thermal_correction_max_iterations
    - \ref man_control_momentum
    - \ref man_localized_lambda
    - \ref man_lumped_lambda_solve
    - \ref man_mask_direction
    - \ref man_time_filter
    - \ref man_filter_scale
    - \ref man_filter_type
    - \ref man_equilibrium_start
    - \ref man_extrinsic_exchange
    - \ref man_poisson_solver

    fix_modify commands for output: \n
    - \ref man_output
    - \ref man_output_nodeset
    - \ref man_output_elementset
    - \ref man_boundary_integral
    - \ref man_contour_integral
    - \ref man_mesh_output
    - \ref man_write_restart
    - \ref man_read_restart

    fix_modify commands for post-processing: \n
    - \ref man_hardy_kernel
    - \ref man_hardy_fields
    - \ref man_hardy_gradients
    - \ref man_hardy_rates
    - \ref man_hardy_computes
    - \ref man_hardy_on_the_fly
    - \ref man_pair_interactions
    - \ref man_sample_frequency
    - \ref man_set

    miscellaneous fix_modify commands: \n
    - \ref man_atom_element_map
    - \ref man_atom_weight
    - \ref man_write_atom_weights
    - \ref man_reset_time
    - \ref man_reset_atomic_reference_positions
    - \ref man_fe_md_boundary
    - \ref man_boundary_faceset
    - \ref man_consistent_fe_initialization
    - \ref man_mass_matrix
    - \ref man_material
    - \ref man_atomic_charge
    - \ref man_source_integration
    - \ref man_temperature_definition
    - \ref man_track_displacement
    - \ref man_boundary_dynamics
    - \ref man_add_species
    - \ref man_add_molecule
    - \ref man_remove_species
    - \ref man_remove_molecule

    Note: a set of example input files with the attendant material files are
    included with this package
    \section default
    none
  */


  // Construct new ATC_Method object
  // note use "unfix" to destroy

  int me = ATC::LammpsInterface::instance()->comm_rank();

  string groupName(arg[1]);
  int igroup = group->find(groupName.c_str());
  int atomCount = group->count(igroup);

  try {
    // Postprocessing
    if (strcmp(arg[3],"field")==0)
    {
      if (atomCount == 0) {
        if (me==0) printf("ATC: can't construct transfer, no atoms in group \n");
        throw;
      }
      if (narg < 5) {
        if (me==0) printf("ATC: constructing shape function field estimate\n");
        atc_ = new ATC::ATC_TransferPartitionOfUnity(groupName,
                                                     array_atom,
                                                     this);
      }
      else {
        if (me==0) printf("ATC: constructing shape function field estimate with parameter file %s\n",arg[4]);
        string matParamFile = arg[4];
        atc_ = new ATC::ATC_TransferPartitionOfUnity(groupName,
                                                     array_atom, this,
                                                     matParamFile);
      }
    }
    else if (strcmp(arg[3],"hardy")==0)
    {
      if (atomCount == 0) {
        if (me==0) printf("ATC: Can't construct transfer, no atoms in group \n");
        throw;
      }
      if (narg < 5) {
        if (me==0) printf("ATC: constructing kernel field estimate\n");
        atc_ = new ATC::ATC_TransferKernel(groupName,
                                           array_atom,
                                           this);
      }
      else {
        if (me==0) printf("ATC: constructing kernel field estimate with parameter file %s\n",arg[4]);
        string matParamFile = arg[4];
        atc_ = new ATC::ATC_TransferKernel(groupName,
                                           array_atom, this,
                                           matParamFile);
      }
    }
    // PhysicsTypes
    else if (strcmp(arg[3],"thermal")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing thermal coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingEnergy(groupName,
                                         array_atom, this,
                                         matParamFile);
    }
    else if (strcmp(arg[3],"two_temperature")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing two_temperature coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingEnergy(groupName,
                                         array_atom, this,
                                         matParamFile, ATC::TWO_TEMPERATURE);
    }
    else if (strcmp(arg[3],"drift_diffusion")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing drift_diffusion coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingEnergy(groupName,
                                         array_atom, this,
                                         matParamFile, ATC::DRIFT_DIFFUSION);
    }
    else if (strcmp(arg[3],"drift_diffusion-equilibrium")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing drift_diffusion-equilibrium coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingEnergy(groupName,
                                         array_atom, this,
                                         matParamFile, ATC::DRIFT_DIFFUSION_EQUILIBRIUM);
    }
    else if (strcmp(arg[3],"drift_diffusion-schrodinger")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing drift_diffusion-schrodinger coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingEnergy(groupName,
                                         array_atom, this,
                                         matParamFile, ATC::DRIFT_DIFFUSION_SCHRODINGER);
    }
    else if (strcmp(arg[3],"drift_diffusion-schrodinger-slice")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("Constructing ATC transfer (drift_diffusion-schrodinger-slice) with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingEnergy(groupName,
                                         array_atom, this,
                                         matParamFile, ATC::DRIFT_DIFFUSION_SCHRODINGER_SLICE);
    }
    else if (strcmp(arg[3],"convective_drift_diffusion")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing convective_drift_diffusion coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingEnergy(groupName,
                                         array_atom, this,
                                         matParamFile, ATC::CONVECTIVE_DRIFT_DIFFUSION);
    }
    else if (strcmp(arg[3],"convective_drift_diffusion-equilibrium")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing convective_drift_diffusion-equilibrium coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingEnergy(groupName,
                                         array_atom, this,
                                         matParamFile, ATC::CONVECTIVE_DRIFT_DIFFUSION_EQUILIBRIUM);
    }
    else if (strcmp(arg[3],"convective_drift_diffusion-schrodinger")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing convective_drift_diffusion-schrodinger coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingEnergy(groupName,
                                         array_atom, this,
                                         matParamFile, ATC::CONVECTIVE_DRIFT_DIFFUSION_SCHRODINGER);
    }
    else if (strcmp(arg[3],"elastic")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing elastic coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingMomentum(groupName,
                                           array_atom, this,
                                           matParamFile,
                                           ATC::ELASTIC);
    }
    else if (strcmp(arg[3],"electrostatic")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing electrostatic mechanical coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingMomentum(groupName,
                                           array_atom, this,
                                           matParamFile,
                                           ATC::ELASTIC,
                                           ATC::ELECTROSTATIC);
    }
    else if (strcmp(arg[3],"electrostatic-equilibrium")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing equilibrium electrostatic coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingMomentum(groupName,
                                           array_atom, this,
                                           matParamFile,
                                           ATC::ELASTIC,
                                           ATC::ELECTROSTATIC_EQUILIBRIUM);
    }
    else if (strcmp(arg[3],"shear")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing viscous/shear coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingMomentum(groupName,
                                           array_atom, this,
                                           matParamFile,
                                           ATC::SHEAR);
    }
    else if (strcmp(arg[3],"species")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing species diffusion coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingMass(groupName,
                                       array_atom, this,
                                       matParamFile);
    }
    else if (strcmp(arg[3],"species_electrostatic")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing electrostatic species coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingMass(groupName,
                                       array_atom, this,
                                       matParamFile, ATC::FEM_EFIELD);
    }
    else if (strcmp(arg[3],"thermo_elastic")==0)
    {
      string matParamFile = arg[4];
      if (me==0) printf("ATC: constructing thermo-mechanical coupling with parameter file %s\n",arg[4]);
      atc_ = new ATC::ATC_CouplingMomentumEnergy(groupName,
                                                 array_atom, this,
                                                 matParamFile);
    }
    else
    {
      lmp->error->all(FLERR,"Unknown physics type in ATC");
    }
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }

  lmp->atom->add_callback(0);

  // we write our own restart file
  restart_global = 0;



  // Set output computation data based on transfer info
  scalar_flag = atc_->scalar_flag();
  vector_flag = atc_->vector_flag();
  size_vector = atc_->size_vector();
  global_freq = atc_->global_freq();
  extscalar = atc_->extscalar();
  extvector = atc_->extvector();
  extlist = atc_->extlist();
  thermo_energy = atc_->thermo_energy_flag();

  // set pointer for output
  peratom_flag = atc_->peratom_flag();
  size_peratom_cols = atc_->size_peratom_cols();
  peratom_freq = atc_->peratom_freq();


  // set comm size needed by this fix
  comm_forward = atc_->comm_forward();

  // call this fix every step
  nevery = 1;
}

/*----------------------------------------------------------------------- */
FixATC::~FixATC()
{
  if (lmp->atom) lmp->atom->delete_callback(id,0);
  if (atc_) delete atc_;
}

/* ---------------------------------------------------------------------- */

int FixATC::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_EXCHANGE;
  mask |= PRE_NEIGHBOR;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= MIN_PRE_EXCHANGE;
  mask |= MIN_PRE_NEIGHBOR;
  mask |= MIN_PRE_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_RUN;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixATC::modify_param(int narg, char** arg)
{
  bool match;

  // pass on to transfer layer
  try {
    match = atc_->modify(narg,arg);
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }

  if (!match) return 0;
  return narg;
}

/* ----------------------------------------------------------------------
   create initial list of neighbor partners via call to neighbor->build()
   must be done in setup (not init) since fix init comes before neigh init
   ------------------------------------------------------------------------- */

void FixATC::init()
{
  // Guarantee construction of full neighborlist
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // create computes, if necessary
  atc_->init_computes();
}

void FixATC::min_setup(int vflag)
{
  setup(vflag);
}

void FixATC::setup(int /* vflag */)
{
  comm->forward_comm_fix(this);

  try {
    atc_->initialize();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}

/* ----------------------------------------------------------------------
   pass throughs to atc functions to handle swapping atom data on
   when they move processors
   ------------------------------------------------------------------------- */
void FixATC::pre_exchange()
{
  try {
    atc_->pre_exchange();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}
void FixATC::setup_pre_exchange()
{
  if (atc_->is_initialized()) {
    try {
      atc_->setup_pre_exchange();
    }
    catch (ATC::ATC_Error& atcError) {
      ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
      throw;
    }
  }
}
void FixATC::min_pre_exchange()
{
  try {
    atc_->pre_exchange();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}

double FixATC::memory_usage()
{
  double bytes = (double) atc_->memory_usage() * sizeof(double);
  return bytes;
}

void FixATC::grow_arrays(int nmax)
{
  atc_->grow_arrays(nmax);
}

void FixATC::copy_arrays(int i, int j, int /* delflag */)
{
  atc_->copy_arrays(i,j);
}

int FixATC::pack_exchange(int i, double * buf)
{
  int num = atc_->pack_exchange(i,buf);
  return num;
}

int FixATC::unpack_exchange(int nlocal, double * buf)
{
  int num = atc_->unpack_exchange(nlocal,buf);
  return num;
}

int FixATC::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int num = atc_->pack_comm(n, list, buf, pbc_flag, pbc);
  return num;
}

void FixATC::unpack_forward_comm(int n, int first, double *buf)
{
  atc_->unpack_comm(n, first, buf);
}


/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
   ------------------------------------------------------------------------- */

int FixATC::pack_restart(int /* i */, double * /* buf */){
  return 0;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
   ------------------------------------------------------------------------- */

void FixATC::unpack_restart(int /* nlocal */, int /* nth */){
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
   ------------------------------------------------------------------------- */

int FixATC::maxsize_restart(){
  return 0;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
   ------------------------------------------------------------------------- */

int FixATC::size_restart(int /* nlocal */){
  return 0;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
   ------------------------------------------------------------------------- */

void FixATC::write_restart(FILE * /* fp */){

  char ** args = new char*[2];
  args[0] = new char[50];
  args[1] = new char[50];
  sprintf(args[0],"write_restart");
  sprintf(args[1],"ATC.restart");

  // Then call all objects I own to write their data
  if (comm->me == 0) {
    atc_->modify(2,args);
  }

  delete [] args[0];
  delete [] args[1];
  delete [] args;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
   ------------------------------------------------------------------------- */

void FixATC::restart(char * /* buf */){

  char ** args = new char*[2];
  args[0] = new char[50];
  args[1] = new char[50];
  sprintf(args[0],"read_restart");
  sprintf(args[1],"ATC.restart");

  // Then call all objects I own to write their data
  if (comm->me == 0) {
    atc_->modify(2,args);
  }

  delete [] args[0];
  delete [] args[1];
  delete [] args;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
   ------------------------------------------------------------------------- */

void FixATC::initial_integrate(int /* vflag */)
{
  try {
    atc_->pre_init_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
  // integration of atoms, if desired
  try {
    atc_->init_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}

void FixATC::post_integrate()
{
  try {
    atc_->post_init_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}

/* ---------------------------------------------------------------------- */

void FixATC::final_integrate()
{
  try {
    atc_->pre_final_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }

  try {
    atc_->final_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}

void FixATC::end_of_step()
{
  try {
    atc_->post_final_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
  try {
    atc_->end_of_step();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}

/* ---------------------------------------------------------------------- */
void FixATC::init_list(int id, NeighList *ptr) {
  ATC::LammpsInterface::instance()->set_list(id,ptr);
}
/* ---------------------------------------------------------------------- */
void FixATC::pre_neighbor()
{
  try {
    atc_->pre_neighbor();
    comm->forward_comm_fix(this);
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}
/* ---------------------------------------------------------------------- */
void FixATC::pre_force(int /* vflag */)
{

  try {
    atc_->pre_force();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}
/* ---------------------------------------------------------------------- */
void FixATC::post_force(int /* vflag */)
{

  try {
    atc_->post_force();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}
/* ---------------------------------------------------------------------- */
void FixATC::post_run()
{
  try {
    atc_->finish();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}
/* ---------------------------------------------------------------------- */
void FixATC::setup_pre_neighbor()
{
  if (atc_->is_initialized()) {
    try {
      atc_->pre_neighbor();
    }
    catch (ATC::ATC_Error& atcError) {
      ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
      throw;
    }
  }
}
/* ---------------------------------------------------------------------- */
void FixATC::min_pre_force(int /* vflag */)
{
  try {
    atc_->min_pre_force();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}

/* ---------------------------------------------------------------------- */
void FixATC::min_post_force(int /* vflag */)
{
  try {
    atc_->min_post_force();
  }
  catch (ATC::ATC_Error& atcError) {
    ATC::LammpsInterface::instance()->print_msg(atcError.error_description());
    throw;
  }
}

/* ---------------------------------------------------------------------- */
double FixATC::compute_scalar()
{
  return atc_->compute_scalar();
}
/* ---------------------------------------------------------------------- */
double FixATC::compute_vector(int n)
{
  return atc_->compute_vector(n);
}
/* ---------------------------------------------------------------------- */
double FixATC::compute_array(int irow, int icol)
{
  return atc_->compute_array(irow,icol);
}

