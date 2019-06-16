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
   Contributing authors: Gareth Tribello (Queens U, Belfast)
                         Pablo Piaggi (EPFL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "atom.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "fix_plumed.h"
#include "universe.h"
#include "compute.h"
#include "modify.h"
#include "pair.h"
#include "utils.h"

#include "plumed/wrapper/Plumed.h"

#if defined(__PLUMED_DEFAULT_KERNEL)
#define PLUMED_QUOTE_DIRECT(name) #name
#define PLUMED_QUOTE(macro) PLUMED_QUOTE_DIRECT(macro)
static char plumed_default_kernel[] = "PLUMED_KERNEL=" PLUMED_QUOTE(__PLUMED_DEFAULT_KERNEL);
#endif

/* -------------------------------------------------------------------- */

using namespace LAMMPS_NS;
using namespace FixConst;

#define INVOKED_SCALAR 1

FixPlumed::FixPlumed(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  p(NULL), nlocal(0), gatindex(NULL), masses(NULL), charges(NULL),
  id_pe(NULL), id_press(NULL)
{

  if (!atom->tag_enable)
    error->all(FLERR,"Fix plumed requires atom tags");

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Fix plumed requires consecutive atom IDs");

  if (igroup != 0 && comm->me == 0)
    error->warning(FLERR,"Fix group for fix plumed is not 'all'. "
                   "Group will be ignored.");

#if defined(__PLUMED_DEFAULT_KERNEL)
  if (getenv("PLUMED_KERNEL") == NULL)
    putenv(plumed_default_kernel);
#endif

  p=new PLMD::Plumed;

  // Check API version

  int api_version;
  p->cmd("getApiVersion",&api_version);
  if (api_version > 6)
    error->all(FLERR,"Incompatible API version for PLUMED in fix plumed");

  // If the -partition option is activated then enable
  // inter-partition communication

  if (universe->existflag == 1) {
    int me;
    MPI_Comm inter_comm;
    MPI_Comm_rank(world,&me);

    // Change MPI_COMM_WORLD to universe->uworld which seems more appropriate

    MPI_Comm_split(universe->uworld,me,0,&inter_comm);
    p->cmd("GREX setMPIIntracomm",&world);
    if (me == 0) {
      // The inter-partition communicator is only defined for the root in
      //    each partition (a.k.a. world). This is due to the way in which
      //    it is defined inside plumed.
      p->cmd("GREX setMPIIntercomm",&inter_comm);
    }
    p->cmd("GREX init",NULL);
  }

  // The general communicator is independent of the existence of partitions,
  // if there are partitions, world is defined within each partition,
  // whereas if partitions are not defined then world is equal to
  // MPI_COMM_WORLD.

  p->cmd("setMPIComm",&world);

  // Set up units
  // LAMMPS units wrt kj/mol - nm - ps
  // Set up units

  if(strcmp(update->unit_style,"lj") == 0) {
    // LAMMPS units lj
    p->cmd("setNaturalUnits");
  } else {

    // Conversion factor from LAMMPS energy units to kJ/mol (units of PLUMED)

    double energyUnits=1.0;

    // LAMMPS units real :: kcal/mol;

    if (strcmp(update->unit_style,"real") == 0) {
      energyUnits=4.184;

      // LAMMPS units metal :: eV;

    } else if (strcmp(update->unit_style,"metal") == 0) {
      energyUnits=96.48530749925792;

      // LAMMPS units si :: Joule;

    } else if (strcmp(update->unit_style,"si") == 0) {
      energyUnits=0.001;

      // LAMMPS units cgs :: erg;

    } else if (strcmp(update->unit_style,"cgs") == 0) {
      energyUnits=6.0221418e13;

      // LAMMPS units electron :: Hartree;

    } else if (strcmp(update->unit_style,"electron") == 0) {
      energyUnits=2625.5257;

    } else error->all(FLERR,"Fix plumed cannot handle your choice of units");

    // Conversion factor from LAMMPS length units to nm (units of PLUMED)

    double lengthUnits=0.1/force->angstrom;

    // Conversion factor from LAMMPS time unit to ps (units of PLUMED)

    double timeUnits=0.001/force->femtosecond;

    p->cmd("setMDEnergyUnits",&energyUnits);
    p->cmd("setMDLengthUnits",&lengthUnits);
    p->cmd("setMDTimeUnits",&timeUnits);
  }

  // Read fix parameters:

  int next=0;
  for (int i=3;i<narg;++i) {
    if (!strcmp(arg[i],"outfile")) {
      next=1;
    } else if (next==1) {
      if (universe->existflag == 1) {
        // Each replica writes an independent log file
        //  with suffix equal to the replica id
        char str_num[32], logFile[1024];
        sprintf(str_num,".%d",universe->iworld);
        strncpy(logFile,arg[i],1024-32);
        strcat(logFile,str_num);
        p->cmd("setLogFile",logFile);
        next=0;
      } else {
        // partition option not used
        p->cmd("setLogFile",arg[i]);
        next=0;
      }
    } else if (!strcmp(arg[i],"plumedfile")) {
      next=2;
    } else if (next==2) {
      p->cmd("setPlumedDat",arg[i]);
      next=0;
    } else error->all(FLERR,"Syntax error - use 'fix <fix-ID> plumed "
                      "plumedfile plumed.dat outfile plumed.out' ");
  }
  if (next==1) error->all(FLERR,"missing argument for outfile option");
  if (next==2) error->all(FLERR,"missing argument for plumedfile option");

  p->cmd("setMDEngine","LAMMPS");

  if (atom->natoms > MAXSMALLINT)
    error->all(FLERR,"Fix plumed can only handle up to 2.1 billion atoms");

  natoms=int(atom->natoms);
  p->cmd("setNatoms",&natoms);

  double dt=update->dt;
  p->cmd("setTimestep",&dt);

  virial_flag=1;
  thermo_virial=1;
  scalar_flag = 1;

  // This is the real initialization:

  p->cmd("init");

  // Define compute to calculate potential energy

  id_pe = new char[8];
  strcpy(id_pe,"plmd_pe");
  char **newarg = new char*[3];
  newarg[0] = id_pe;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg);
  delete [] newarg;
  int ipe = modify->find_compute(id_pe);
  c_pe = modify->compute[ipe];

  // Define compute to calculate pressure tensor

  id_press = new char[11];
  strcpy(id_press,"plmd_press");
  newarg = new char*[5];
  newarg[0] = id_press;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = (char *) "NULL";
  newarg[4] = (char *) "virial";
  modify->add_compute(5,newarg);
  delete [] newarg;
  int ipress = modify->find_compute(id_press);
  c_press = modify->compute[ipress];

  for (int i = 0; i < modify->nfix; i++) {
    const char * const check_style = modify->fix[i]->style;

    // There must be only one

    if (strcmp(check_style,"plumed") == 0)
      error->all(FLERR,"There must be only one instance of fix plumed");

    // Avoid conflict with fixes that define internal pressure computes.
    // See comment in the setup method

    if (utils::strmatch(check_style,"^nph") ||
        utils::strmatch(check_style,"^npt") ||
        utils::strmatch(check_style,"^rigid/nph") ||
        utils::strmatch(check_style,"^rigid/npt") ||
        utils::strmatch(check_style,"^msst") ||
        utils::strmatch(check_style,"^nphug") ||
        utils::strmatch(check_style,"^ipi") ||
        utils::strmatch(check_style,"^press/berendsen") ||
        utils::strmatch(check_style,"^qbmsst"))
      error->all(FLERR,"Fix plumed must be defined before any other fixes, "
                 "that compute pressure internally");
  }
}

FixPlumed::~FixPlumed()
{
  delete p;
  modify->delete_compute(id_pe);
  modify->delete_compute(id_press);
  delete[] id_pe;
  delete[] id_press;
  delete[] masses;
  delete[] charges;
  delete[] gatindex;
}

int FixPlumed::setmask()
{
  // set with a bitmask how and when apply the force from plumed
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

void FixPlumed::init()
{
  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // This avoids nan pressure if compute_pressure is called
  // in a setup method

  for (int i=0;i<6;i++) virial[i] = 0.;
}

void FixPlumed::setup(int vflag)
{
  // Here there is a crucial issue connected to constant pressure
  // simulations. The fix_nh will call the compute_pressure inside
  // the setup method, that is executed once and for all at the
  // beginning of the simulation. Since our fix has a contribution
  // to the virial, when this happens the variable virial must have
  // been calculated. In other words, the setup method of fix_plumed
  // has to be executed first. This creates a race condition with the
  // setup method of fix_nh. This is why in the constructor I check if
  // nh fixes have already been called.
  if (utils::strmatch(update->integrate_style,"^respa")) {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  } else {
    post_force(vflag);
  }
}

void FixPlumed::min_setup(int vflag)
{
  // This has to be checked.
  // For instance it might have problems with fix_box_relax
  post_force(vflag);
}

void FixPlumed::post_force(int /* vflag */)
{

  int update_gatindex=0;

  if (natoms != int(atom->natoms))
    error->all(FLERR,"Fix plumed does not support simulations with varying "
               "numbers of atoms");

  // Try to find out if the domain decomposition has been updated:

  if (nlocal != atom->nlocal) {

    if (charges) delete [] charges;
    if (masses) delete [] masses;
    if (gatindex) delete [] gatindex;

    nlocal=atom->nlocal;
    gatindex=new int [nlocal];
    masses=new double [nlocal];
    charges=new double [nlocal];
    update_gatindex=1;

  } else {

    for (int i=0;i<nlocal;i++) {
      if (gatindex[i]!=atom->tag[i]-1) {
        update_gatindex=1;
        break;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,&update_gatindex,1,MPI_INT,MPI_SUM,world);

  // In case it has been updated, rebuild the local mass/charges array
  // and tell plumed about the change:

  if (update_gatindex) {
    for (int i=0;i<nlocal;i++) gatindex[i]=atom->tag[i]-1;
    // Get masses
    if (atom->rmass_flag) {
      for (int i=0;i<nlocal;i++) masses[i]=atom->rmass[i];
    } else {
      for (int i=0;i<nlocal;i++) masses[i]=atom->mass[atom->type[i]];
    }
    // Get charges
    if (atom->q_flag) {
      for (int i=0;i<nlocal;i++) charges[i]=atom->q[i];
    } else {
      for (int i=0;i<nlocal;i++) charges[i]=0.0;
    }
    p->cmd("setAtomsNlocal",&nlocal);
    p->cmd("setAtomsGatindex",gatindex);
  }


  // set up local virial/box. plumed uses full 3x3 matrices
  double plmd_virial[3][3];
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) plmd_virial[i][j]=0.0;
  double box[3][3];
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) box[i][j]=0.0;
  box[0][0]=domain->h[0];
  box[1][1]=domain->h[1];
  box[2][2]=domain->h[2];
  box[2][1]=domain->h[3];
  box[2][0]=domain->h[4];
  box[1][0]=domain->h[5];

  // Make initial of virial of this fix zero
  // The following line is very important, otherwise
  // the compute pressure will include
  for (int i=0;i<6;++i) virial[i] = 0.;

  // local variable with timestep:
  if (update->ntimestep > MAXSMALLINT)
    error->all(FLERR,"Fix plumed can only handle up to 2.1 billion timesteps");
  int step=int(update->ntimestep);

  // pass all pointers to plumed:
  p->cmd("setStep",&step);
  int plumedStopCondition=0; 
  p->cmd("setStopFlag",&plumedStopCondition);
  p->cmd("setPositions",&atom->x[0][0]);
  p->cmd("setBox",&box[0][0]);
  p->cmd("setForces",&atom->f[0][0]);
  p->cmd("setMasses",&masses[0]);
  p->cmd("setCharges",&charges[0]);
  p->cmd("getBias",&bias);

  // Pass virial to plumed
  // If energy is needed plmd_virial is equal to Lammps' virial
  // If energy is not needed plmd_virial is initialized to zero
  // In the first case the virial will be rescaled and an extra term will be added
  // In the latter case only an extra term will be added
  p->cmd("setVirial",&plmd_virial[0][0]);
  p->cmd("prepareCalc");

  plumedNeedsEnergy=0;
  p->cmd("isEnergyNeeded",&plumedNeedsEnergy);

  // Pass potential energy and virial if needed
  double *virial_lmp;
  if (plumedNeedsEnergy) {
    // Error if tail corrections are included
    if (force->pair && force->pair->tail_flag && comm->me == 0)
      error->warning(FLERR,"Tail corrections to the pair potential included."
                     " The energy cannot be biased correctly in this case."
                     " Remove the tail corrections by removing the"
                     " command: pair_modify tail yes");

    // compute the potential energy
    double pot_energy = 0.;
    c_pe->compute_scalar();
    pot_energy = c_pe->scalar;

    // Divide energy by number of processes
    // Plumed wants it this way
    int nprocs;
    MPI_Comm_size(world,&nprocs);
    pot_energy /= nprocs;
    p->cmd("setEnergy",&pot_energy);

    // Compute pressure due to the virial (no kinetic energy term!)
    c_press->compute_vector();
    virial_lmp = c_press->vector;

    // Check if pressure is finite
    if (!std::isfinite(virial_lmp[0]) || !std::isfinite(virial_lmp[1])
        || !std::isfinite(virial_lmp[2]) || !std::isfinite(virial_lmp[3])
        || !std::isfinite(virial_lmp[4]) || !std::isfinite(virial_lmp[5]))
      error->all(FLERR,"Non-numeric virial - Plumed cannot work with that");

    // Convert pressure to virial per number of MPI processes
    // From now on all virials are divided by the number of MPI processes

    double nktv2p = force->nktv2p;
    double inv_volume;
    if (domain->dimension == 3) {
      inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    } else {
      inv_volume = 1.0 / (domain->xprd * domain->yprd);
    }
    for (int i=0;i<6;i++) virial_lmp[i] /= (inv_volume * nktv2p * nprocs);
    // Convert virial from lammps to plumed representation
    plmd_virial[0][0]=-virial_lmp[0];
    plmd_virial[1][1]=-virial_lmp[1];
    plmd_virial[2][2]=-virial_lmp[2];
    plmd_virial[0][1]=-virial_lmp[3];
    plmd_virial[0][2]=-virial_lmp[4];
    plmd_virial[1][2]=-virial_lmp[5];
  }
  // do the real calculation:
  p->cmd("performCalc");

  if(plumedStopCondition) error->all(FLERR,"received instruction from PLUMED to stop");

  // retransform virial to lammps representation and assign it to this
  // fix's virial. If the energy is biased, Plumed is giving back the full
  // virial and therefore we have to subtract the initial virial i.e. virial_lmp.
  // The vector virial contains only the contribution added by plumed.
  // The calculation of the pressure will be done by a compute pressure
  // and will include this contribution.
  if (plumedNeedsEnergy) {
    virial[0] = -plmd_virial[0][0]-virial_lmp[0];
    virial[1] = -plmd_virial[1][1]-virial_lmp[1];
    virial[2] = -plmd_virial[2][2]-virial_lmp[2];
    virial[3] = -plmd_virial[0][1]-virial_lmp[3];
    virial[4] = -plmd_virial[0][2]-virial_lmp[4];
    virial[5] = -plmd_virial[1][2]-virial_lmp[5];
  } else {
    virial[0] = -plmd_virial[0][0];
    virial[1] = -plmd_virial[1][1];
    virial[2] = -plmd_virial[2][2];
    virial[3] = -plmd_virial[0][1];
    virial[4] = -plmd_virial[0][2];
    virial[5] = -plmd_virial[1][2];
  }

  // Ask for the computes in the next time step
  // such that the virial and energy are tallied.
  // This should be changed to something that triggers the
  // calculation only if plumed needs it.
  c_pe->addstep(update->ntimestep+1);
  c_press->addstep(update->ntimestep+1);
}

void FixPlumed::post_force_respa(int vflag, int ilevel, int /* iloop */)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

void FixPlumed::min_post_force(int vflag)
{
  post_force(vflag);
}

void FixPlumed::reset_dt()
{
  error->all(FLERR,"Cannot change the time step when fix plumed is active");
}

double FixPlumed::compute_scalar()
{
  return bias;
}

int FixPlumed::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"pe") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    modify->delete_compute(id_pe);
    delete[] id_pe;
    int n = strlen(arg[1]) + 1;
    id_pe = new char[n];
    strcpy(id_pe,arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify potential energy ID");
    c_pe = modify->compute[icompute];

    if (c_pe->peflag == 0)
      error->all(FLERR,"Fix_modify plmd_pe ID does not compute potential energy");
    if (c_pe->igroup != 0 && comm->me == 0)
      error->warning(FLERR,"Potential for fix PLUMED is not for group all");

    return 2;

  } else if (strcmp(arg[0],"press") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    modify->delete_compute(id_press);
    delete[] id_press;
    int n = strlen(arg[1]) + 1;
    id_press = new char[n];
    strcpy(id_press,arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify pressure ID");
    c_press = modify->compute[icompute];

    if (c_press->pressflag == 0)
      error->all(FLERR,"Fix_modify pressure ID does not compute pressure");
    if (c_press->igroup != 0 && comm->me == 0)
      error->warning(FLERR,"Virial for fix PLUMED is not for group all");

    return 2;
  }
  return 0;
}

double FixPlumed::memory_usage()
{
  return double((8+8+4)*atom->nlocal);
}
