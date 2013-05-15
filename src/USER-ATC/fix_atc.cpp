/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Reese Jones, Jon Zimmerman, Jeremy Templeton (SNL)
------------------------------------------------------------------------- */

#include "stdio.h"
#include "string.h"
#include "fix_atc.h"
#include "ATC_TransferHardy.h"
#include "ATC_TransferThermal.h"
#include "ExtrinsicModel.h"
#include "LammpsInterface.h"
#include "fix_nve.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "pointers.h"
#include "comm.h"
#include "error.h"

using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;

// main page of doxygen documentation
/*! \mainpage AtC : Atom-to-Continuum methods
    fix commands:
    - \ref man_fix_atc (links to all related commands)
*/

/* ------------------------------------------------------------------------- */

FixATC::FixATC(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  // ID GROUP atc PHYSICSTYPE [PARAMETERFILE]
  // can have either 4 or 5 args, only arg[3] and arg[4] are used by this class
  if (narg > 5 || narg < 4) lmp->error->all(FLERR,"Illegal fix atc command");

  // Set LAMMPS pointer on LammpsInterface
  ATC::LammpsInterface::instance()->set_lammps(lmp);

  /*! \page man_fix_atc fix atc command
    \section syntax
    fix AtC transfer <type> <parameter_file>
    - type\n
    = thermal : thermal coupling with fields: temperature  \n
    = two_temperature : electron-phonon coupling with field: temperature and electron_temperature  \n
    = hardy : Hardy on-the-fly post-processing (see "related" section for possible fields) \n
    - parameter_file = name of the file with material parameters. \n
    note: hardy does not require a parameter file
    \section examples
    <TT> fix_modify AtC transfer thermal Ar_thermal.dat </TT> \n
    <TT> fix_modify AtC transfer hardy  </TT>
    \section description
    This fix is the beginning to creating a coupled FE/MD simulation and/or
    an on-the-fly estimation of continuum fields. The coupled versions of this
    fix do Verlet integration and the Hardy/post-processing does not.
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
     fix_modify   AtC fem  create mesh 12 2  2  mdRegion f p p\n \n
     # specify the control method for the type of coupling \n
     #                         physics         control_type \n
     fix_modify   AtC transfer thermal control flux \n \n
     # specify the initial values for the empirical field "temperature"  \n
     #                                 field       node_group  value \n
     fix_modify   AtC transfer initial temperature all         30.\n \n
     # create an output stream for nodal fields \n
     #                                filename      output_frequency  \n
     fix_modify   AtC transfer output atc_fe_output 100\n \n

     run             1000 \n
</TT>

     likewise for this post-processing example: \n

<TT>
     # ... commands to create and initialize the MD system \n

     # initial fix to designate post-processing and the group to apply it to \n
     # no material file is allowed nor required \n
     fix             AtC internal atc hardy \n \n
     # create a uniform 1 x 1 x 1 mesh that covers region contain the group \n
     # with periodicity this effectively creats a system average \n
     fix_modify  AtC fem create mesh 1 1 1 box p p p \n\n
     # change from default lagrangian map to eulerian \n
     #   refreshed every 100 steps \n
     fix_modify  AtC atom_element_map eulerian 100 \n \n
     # start with no field defined \n
     fix_modify  AtC transfer fields none \n \n
     # add mass density, potential energy density, stress and temperature \n
     fix_modify  AtC transfer fields add density energy stress temperature \n\n
     # create an output stream for nodal fields \n
     #                                filename      output_frequency  \n
     fix_modify  AtC transfer output nvtFE 100 text \n

     run             1000 \n
</TT>

    Note coupling and post-processing can be combined in the same simulations
    using separate fixes.
    \n
    For detailed exposition of the theory and algorithms please see:\n
    - Wagner, GJ; Jones, RE; Templeton, JA; Parks, MA.  <VAR> An
      atomistic-to-continuum coupling method for heat transfer in solids. </VAR>
      Special Issue of Computer Methods and Applied Mechanics (2008) 197:3351.
    - Zimmerman, JA; Webb, EB; Hoyt, JJ;. Jones, RE; Klein, PA; Bammann, DJ,
      <VAR> Calculation of stress in atomistic simulation </VAR>
      Special Issue of Modelling and Simulation in Materials Science and
      Engineering (2004), 12:S319

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
    - \ref man_fem_mesh
    - \ref man_mesh_nodeset
    - \ref man_mesh_faceset
    - \ref man_mesh_elemset
    - \ref man_transfer_internal
    - \ref man_transfer_boundary
    - \ref man_internal_quadrature
    - \ref man_time_integration
    - \ref man_electron_integration

    fix_modify commands for boundary and initial conditions:\n
    - \ref man_initial
    - \ref man_fix_nodes
    - \ref man_unfix_nodes
    - \ref man_fix_flux
    - \ref man_unfix_flux
    - \ref man_source
    - \ref man_remove_source

    fix_modify commands for control and filtering: \n
    - \ref man_thermal_control
    - \ref man_time_filter
    - \ref man_filter_scale
    - \ref man_equilibrium_start
    - \ref man_extrinsic_exchange

    fix_modify commands for output: \n
    - \ref man_transfer_output
    - \ref man_transfer_atomic_output
    - \ref man_mesh_output
    - \ref man_write_restart
    - \ref man_read_restart

    fix_modify commands for post-processing: \n
    - \ref man_hardy_fields
    - \ref man_hardy_gradients
    - \ref man_hardy_rates
    - \ref man_hardy_computes
    - \ref man_hardy_set
    - \ref man_hardy_on_the_fly
    - \ref man_boundary_integral
    - \ref man_contour_integral

    miscellaneous fix_modify commands: \n
    - \ref man_atom_element_map
    - \ref man_neighbor_reset_frequency

    Note: a set of example input files with the attendant material files are
    included with this package
    \section default
    none
  */

  // Construct new ATC_Transfer object
  // note use "unfix" to destroy

  int me = ATC::LammpsInterface::instance()->comm_rank();

  string groupName(arg[1]);

  // Postprocessing
  try {
    if (strcmp(arg[3],"hardy")==0)
    {
      if (narg < 5) {
        if (me==0) printf("Constructing ATC transfer (hardy)\n");
        atcTransfer_ = new ATC::ATC_TransferHardy(groupName);
      }
      else {
        if (me==0) printf("Constructing ATC transfer (hardy) with parameter file %s\n",arg[4]);
        std::string matParamFile = arg[4];
        atcTransfer_ = new ATC::ATC_TransferHardy(groupName,matParamFile);
      }
    }
    // PhysicsTypes
    else if (strcmp(arg[3],"thermal")==0)
    {
      std::string matParamFile = arg[4];
      if (me==0) printf("Constructing ATC transfer (thermal) with parameter file %s\n",arg[4]);
      atcTransfer_ = new ATC::ATC_TransferThermal(groupName,matParamFile);
      lmp->atom->add_callback(0); // NOTE what is this?
    }
    else if (strcmp(arg[3],"two_temperature")==0)
    {
      std::string matParamFile = arg[4];
      if (me==0) printf("Constructing ATC transfer (two_temperature) with parameter file %s\n",arg[4]);
      atcTransfer_ = new ATC::ATC_TransferThermal(groupName,matParamFile,
                         ATC::TWO_TEMPERATURE);
    }
    else
    {
      lmp->error->all(FLERR,"Unknown physics type in ATC");
    }
  }
  catch (ATC::ATC_Error& atcError) {
    cout << "ATC ERROR: " << atcError.get_error_description() << " processor: " << atcError.get_id() <<  endl;
    throw;
  }

  lmp->atom->add_callback(0);

  // Tell LAMMPS we want to write_restart to be called (writes fix state)
  //restart_global = 1;
  restart_global = 0; // NOTE : turning off ATC restart
  // Tell LAMMPS we want to pack_restart to be called (writes per-atom data)
  //restart_peratom = 1;

  // Set output computation data based on transfer info
  scalar_flag = atcTransfer_->scalar_flag();
  vector_flag = atcTransfer_->vector_flag();
  size_vector = atcTransfer_->size_vector();
  global_freq = atcTransfer_->global_freq();
  extscalar = atcTransfer_->extscalar();
  extvector = atcTransfer_->extvector();
  extlist = atcTransfer_->extlist();

  // set comm size needed by this fix
  comm_forward = 3;
}

/*----------------------------------------------------------------------- */
FixATC::~FixATC()
{
  if (lmp->atom) lmp->atom->delete_callback(id,0);
  if (atcTransfer_) delete atcTransfer_;
}

/* ---------------------------------------------------------------------- */
void FixATC::init_list(int id, NeighList *ptr) {
  ATC::LammpsInterface::instance()->init_list(id,ptr);
}

/* ---------------------------------------------------------------------- */

int FixATC::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_EXCHANGE;
  mask |= MIN_PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixATC::modify_param(int narg, char** arg)
{
  bool match;

  // pass on to transfer layer
  try {
    match = atcTransfer_->modify(narg,arg);
  }
  catch (ATC::ATC_Error& atcError) {
    cout << "ATC ERROR: " << atcError.get_error_description() << " processor: " << atcError.get_id() <<  endl;
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
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

void FixATC::setup(int vflag)
{
  try {
    atcTransfer_->initialize();
  }
  catch (ATC::ATC_Error& atcError) {
    cout << "ATC ERROR: " << atcError.get_error_description() << " processor: " << atcError.get_id() <<  endl;
    throw;
  }
  groupbit = atcTransfer_->get_groupbit();
}

/* ---------------------------------------------------------------------- */
double FixATC::compute_scalar() {return atcTransfer_->compute_scalar();}

/* ---------------------------------------------------------------------- */

double FixATC::compute_vector(int n) {return atcTransfer_->compute_vector(n);}

/* ----------------------------------------------------------------------
   pass throughs to atc functions to handle swapping atom data on
   when they move processors
   ------------------------------------------------------------------------- */
void FixATC::pre_exchange()
{
  atcTransfer_->pre_exchange();
}
void FixATC::min_pre_exchange()
{
  atcTransfer_->pre_exchange();
}

double FixATC::memory_usage()
{
  double bytes = (double) atcTransfer_->memory_usage() * sizeof(double);
  return bytes;
}

void FixATC::grow_arrays(int nmax)
{
  atcTransfer_->grow_arrays(nmax);
}

void FixATC::copy_arrays(int i, int j, int delflag)
{
  atcTransfer_->copy_arrays(i,j);
}

int FixATC::pack_exchange(int i, double * buf)
{
  int num = atcTransfer_->pack_exchange(i,buf);
  return num;
}

int FixATC::unpack_exchange(int nlocal, double * buf)
{
  int num = atcTransfer_->unpack_exchange(nlocal,buf);
  return num;
}

int FixATC::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int num = atcTransfer_->pack_comm(n, list, buf, pbc_flag, pbc);
  return num;
}

void FixATC::unpack_comm(int n, int first, double *buf)
{
  atcTransfer_->unpack_comm(n, first, buf);
}


/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
   ------------------------------------------------------------------------- */

int FixATC::pack_restart(int i, double *buf){
  return 0;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
   ------------------------------------------------------------------------- */

void FixATC::unpack_restart(int nlocal, int nth){
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

int FixATC::size_restart(int nlocal){
  return 0;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
   ------------------------------------------------------------------------- */

void FixATC::write_restart(FILE *fp){
  // hardcode filename for now
  char ** args = new char*[2];
  args[0] = new char[50];
  args[1] = new char[50];
  sprintf(args[0],"write_restart");
  sprintf(args[1],"ATC.restart");

  // Then call all objects I own to write their data
  if (comm->me == 0) {
    atcTransfer_->modify(2,args);
  }

  delete [] args[0];
  delete [] args[1];
  delete [] args;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
   ------------------------------------------------------------------------- */

void FixATC::restart(char *buf){
  // hardcode filename for now
  char ** args = new char*[2];
  args[0] = new char[50];
  args[1] = new char[50];
  sprintf(args[0],"read_restart");
  sprintf(args[1],"ATC.restart");

  // Then call all objects I own to write their data
  if (comm->me == 0) {
    atcTransfer_->modify(2,args);
  }

  delete [] args[0];
  delete [] args[1];
  delete [] args;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
   ------------------------------------------------------------------------- */

void FixATC::initial_integrate(int vflag)
{
  try {
    atcTransfer_->pre_init_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    cout << "ATC ERROR: " << atcError.get_error_description() << " processor: " << atcError.get_id() <<  endl;
    throw;
  }

  try {
    atcTransfer_->init_integrate_velocity();
  }
  catch (ATC::ATC_Error& atcError) {
    cout << "ATC ERROR: " << atcError.get_error_description() << " processor: " << atcError.get_id() <<  endl;
    throw;
  }

  try {
    atcTransfer_->mid_init_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    cout << "ATC ERROR: " << atcError.get_error_description() << " processor: " << atcError.get_id() <<  endl;
    throw;
  }

  try {
    atcTransfer_->init_integrate_position();
  }
  catch (ATC::ATC_Error& atcError) {
    cout << "ATC ERROR: " << atcError.get_error_description() << " processor: " << atcError.get_id() <<  endl;
    throw;
  }

  try {
    atcTransfer_->post_init_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    cout << "ATC ERROR: " << atcError.get_error_description() << " processor: " << atcError.get_id() <<  endl;
    throw;
  }
}

/* ---------------------------------------------------------------------- */

void FixATC::final_integrate()
{
  // need updated ghost atom positions
  comm->forward_comm_fix(this);

  try {
    atcTransfer_->pre_final_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    cout << "ATC ERROR: " << atcError.get_error_description() << " processor: " << atcError.get_id() <<  endl;
    throw;
  }

  try {
    atcTransfer_->final_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    cout << "ATC ERROR: " << atcError.get_error_description() << " processor: " << atcError.get_id() <<  endl;
    throw;
  }

  try {
    atcTransfer_->post_final_integrate();
  }
  catch (ATC::ATC_Error& atcError) {
    cout << "ATC ERROR: " << atcError.get_error_description() << " processor: " << atcError.get_id() <<  endl;
    throw;
  }
}
