/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------

QMMM with LAMMPS as MDI driver for MM and a quantum code as MDI engine for QM
specified group = QM atoms, remaining atoms are MM atoms
two modes for QMMM coupling: DIRECT and POTENTIAL

-- DIRECT mode:

charge on all QM atoms must be zero
QM atoms have no bonds between them

post_force:
  MDI send QM coords via >COORDS
  MDI send MM coords, MM charges via >CLATTICE, >LATTICE
  invoke the QM code
  MDI recv QM energy via <PE
  MDI recv QM forces via <FORCES ?
  MDI recv MM forces via <FORCES ?
  total force on all atoms = MM forces + received QM/MM forces
  total energy = MM energy + received QM energy

-- POTENTIAL mode:

charge on QM atoms is set by QM code each timestep
QM atoms have no bonds between them

pre_force:
  compute Coulombic portion of FF, pair_style and optional kspace
  ecoul[] = per-atom Coulomb energy for each atom
    may require comm from ghost atoms for pairwise terms
  calculate Coulomb potential for each QM atom
    qpotential[i] = ecoul[i] - double sum over QM atoms (qI qJ / rij)
    divide by Q of QM atom
  MDI send QM coords, QM Coulomb pnential via >COORDS, >POTENTIAL_AT_NUCLEI
  invoke the QM code
  MDI recv QM energy, forces, charges via <PE, <FORCES, <CHARGES
  reset charges on QM atoms, comm to ghost atoms
  store QM energy and forces

post_force:
  QM engine should subtract QM/QM contributions to energy and forces
    else LAMMPS will need to do it
  total force on all atoms = MM forces + stored QM forces
  total energy = MM energy + stored QM energy

------------------------------------------------------------------------- */

#include "fix_mdi_qmmm.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "integrate.h"
#include "kspace.h"
#include "memory.h"
#include "min.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NATIVE, REAL, METAL };    // LAMMPS units which MDI supports
enum { DIRECT, POTENTIAL };      // mode of QMMM coupling

#define MAXELEMENT 103    // used elsewhere in MDI package

// prototype for non-class compare function for sorting QM IDs

static int compare_IDs(const int, const int, void *);

/* ---------------------------------------------------------------------- */

FixMDIQMMM::FixMDIQMMM(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  // check requirements for LAMMPS to work with MDI for a QMMM engine

  if (domain->dimension == 2)
    error->all(FLERR,"Fix mdi/qmmm requires 3d simulation");

  if (!atom->tag_enable) 
    error->all(FLERR,"Fix mdi/qmmm requires atom IDs be defined");

  if (!atom->tag_consecutive()) 
    error->all(FLERR,"Fix mdi/qmmm requires atom IDs be consecutive");

  if (atom->map_style == Atom::MAP_NONE) 
    error->all(FLERR,"Fix mdi/qmmm requires an atom map be defined");

  // confirm LAMMPS is being run as a driver

  int role;
  MDI_Get_role(&role);
  if (role != MDI_DRIVER)
    error->all(FLERR, "Must invoke LAMMPS as an MDI driver to use fix mdi/qmmm");

  // mode arg

  if (strcmp(arg[3],"direct") == 0) mode = DIRECT;
  else if (strcmp(arg[3],"potential") == 0) mode = POTENTIAL;
  else error->all(FLERR,"Illegal fix mdi/qmmm command");
  
  // optional args

  virialflag = 0;
  connectflag = 1;
  elements = nullptr;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "virial") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix mdi/qmmm command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        virialflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        virialflag = 0;
      else
        error->all(FLERR, "Illegal fix mdi/qmmm command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "connect") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix mdi/qmmm command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        connectflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        connectflag = 0;
      else
        error->all(FLERR, "Illegal fix mdi/qmmm command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "elements") == 0) {
      int ntypes = atom->ntypes;
      if (iarg + ntypes + 1 > narg) error->all(FLERR, "Illegal fix mdi/qmmm command");
      delete[] elements;
      elements = new int[ntypes + 1];
      for (int i = 1; i <= ntypes; i++) {
        elements[i] = utils::inumeric(FLERR, arg[iarg + i], false, lmp);
        if (elements[i] < 1 || elements[i] > MAXELEMENT)
          error->all(FLERR, "Illegal fix mdi/qmmm command");
      }
      iarg += ntypes + 1;
    } else
      error->all(FLERR, "Illegal fix mdi/qmmm command");
  }

  // fix output settings are based on optional keywords

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;

  if (virialflag) {
    vector_flag = 1;
    size_vector = 6;
  }

  energy_global_flag = 1;
  thermo_energy = thermo_virial = 1;
  comm_forward = 1;
  comm_reverse = 1;

  // nqm = size of fix group = total # of QM atoms
  //   error if nqm == natoms, should use fix mdi/qm instead
  //   require 3*nqm be a small INT, so can MPI_Allreduce xqm
  // nmm = total # of MM atoms = all non-QM atoms

  nprocs = comm->nprocs;
  bigint ngroup = group->count(igroup);
  
  if (ngroup == 0) error->all(FLERR,"Fix mdi/qmmm has no atoms in quantum group");
  if (ngroup == atom->natoms) error->all(FLERR,"Fix mdi/qmmm has all atoms in quantum group");
  if (3*ngroup > MAXSMALLINT) 
    error->all(FLERR,"Fix mdi/qmmm quantum group has too many atoms");
  nqm = ngroup;
  nmm = atom->natoms = nqm;
  
  // QM atom memory
  
  memory->create(qmIDs,nqm,"mdi/qmmm:qmIDs");
  memory->create(xqm,nqm,3,"mdi/qmmm:xqm");
  memory->create(fqm,nqm,3,"mdi/qmmm:fqm");
  memory->create(qqm,nqm,"mdi/qmmm:qqm");
  memory->create(eqm,nqm,"mdi/qmmm:eqm");
  memory->create(qpotential,nqm,"mdi/qmmm:qpotential");
  memory->create(xqm_mine,nqm,3,"mdi/qmmm:xqm_mine");
  memory->create(qqm_mine,nqm,"mdi/qmmm:qqm_mine");
  memory->create(eqm_mine,nqm,"mdi/qmmm:eqm_mine");
  memory->create(qpotential_mine,nqm,"mdi/qmmm:qpotential_mine");
  memory->create(qm2owned,nqm,"mdi/qmmm:qm2owned");

  // MM atom memory

  memory->create(mmIDs,nmm,"mdi/qmmm:mmIDs");
  memory->create(xmm,nmm,3,"mdi/qmmm:xmm");
  memory->create(fmm,nmm,3,"mdi/qmmm:fmm");
  memory->create(qmm,nmm,"mdi/qmmm:qmm");
  memory->create(emm,nmm,"mdi/qmmm:emm");
  memory->create(xmm_mine,nmm,3,"mdi/qmmm:xmm_mine");
  memory->create(qmm_mine,nmm,"mdi/qmmm:qmm_mine");
  memory->create(emm_mine,nmm,"mdi/qmmm:emm_mine");
  memory->create(mm2owned,nmm,"mdi/qmmm:mm2owned");
  
  // create ordered list of QM and MM atoms IDs

  create_qm_list();
  create_mm_list();
  
  // peratom Coulombic energy

  ecoul = nullptr;
  ncoulmax = 0;
  
  // mdicomm will be initialized in init()
  // cannot do here for a plugin library, b/c mdi plugin command comes later

  mdicomm = MDI_COMM_NULL;
  first_send = 1;

  // set MDI unit conversion factors

  if (strcmp(update->unit_style, "real") == 0)
    lmpunits = REAL;
  else if (strcmp(update->unit_style, "metal") == 0)
    lmpunits = METAL;
  else
    lmpunits = NATIVE;

  unit_conversions();

  // initialize outputs

  qm_energy = 0.0;
  if (virialflag) {
    for (int i = 0; i < 6; i++) {
      qm_virial[i] = 0.0;
      virial[i] = 0.0;
    }
    sumflag = 0;
  }
}

/* ---------------------------------------------------------------------- */

FixMDIQMMM::~FixMDIQMMM()
{
  // send exit command to stand-alone engine code
  // for connnectflag = 0, this is done via "mdi exit" command
  // for plugin, this is done in MDIPlugin::plugin_wrapper()

  if (mdicomm != MDI_COMM_NULL && connectflag && !plugin) {
    int ierr = MDI_Send_command("EXIT", mdicomm);
    if (ierr) error->all(FLERR, "MDI: EXIT command");
  }

  // clean up

  delete [] elements;

  memory->destroy(qmIDs);
  memory->destroy(xqm);
  memory->destroy(fqm);
  memory->destroy(qqm);
  memory->destroy(eqm);
  memory->destroy(qpotential);
  memory->destroy(xqm_mine);
  memory->destroy(qqm_mine);
  memory->destroy(eqm_mine);
  memory->destroy(qpotential_mine);
  memory->destroy(qm2owned);

  memory->destroy(mmIDs);
  memory->destroy(xmm);
  memory->destroy(fmm);
  memory->destroy(qmm);
  memory->destroy(emm);
  memory->destroy(xmm_mine);
  memory->destroy(qmm_mine);
  memory->destroy(emm_mine);
  memory->destroy(mm2owned);

  memory->destroy(ecoul);
}

/* ---------------------------------------------------------------------- */

int FixMDIQMMM::setmask()
{
  int mask = 0;
  mask |= POST_NEIGHBOR;
  mask |= MIN_POST_NEIGHBOR;
  if (mode == POTENTIAL) {
    mask |= PRE_FORCE;
    mask |= MIN_PRE_FORCE;
  }
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::init()
{
  // set local mdicomm one-time only
  // also set plugin = 0/1 for engine = stand-alone code vs plugin library

  if (mdicomm == MDI_COMM_NULL) {

    // this fix makes one-time connection to engine

    if (connectflag) {

      // if MDI's mdicomm not set, need to Accept_comm() with stand-alone engine
      // othewise are already connected to plugin engine

      MDI_Get_communicator(&mdicomm, 0);

      if (mdicomm == MDI_COMM_NULL) {
        plugin = 0;
        MDI_Accept_communicator(&mdicomm);
        if (mdicomm == MDI_COMM_NULL)
          error->all(FLERR, "MDI unable to connect to stand-alone engine");

      } else {
        plugin = 1;
        int method;
        MDI_Get_method(&method, mdicomm);
        if (method != MDI_PLUGIN) error->all(FLERR, "MDI internal error for plugin engine");
      }

      // connection should have been already made by "mdi connect" command
      // only works for stand-alone engines

    } else {
      plugin = 0;

      if (lmp->mdicomm == nullptr)
        error->all(FLERR, "Fix mdi/qmmm is not connected to engine via mdi connect");

      int nbytes = sizeof(MDI_Comm);
      char *ptrcomm = (char *) lmp->mdicomm;
      memcpy(&mdicomm, ptrcomm, nbytes);
    }
  }

  // require per-atom charge for mode = either DIRECT or POTENTIAL
  // POTENTIAL mode requires a pair style that calculates only Coulombic interactions
  //   can also be in conjunction with KSpace solver

  if (!atom->q_flag) error->all(FLERR,"Fix mdi/qmmm requires per-atom charge");

  if (mode == POTENTIAL) {
    if (!force->pair) error->all(FLERR,"Fix mdi/qmmm potential requires a pair style");
    pair_coul = force->pair_match("coul/cut",1,0);
    if (!pair_coul) pair_coul = force->pair_match("coul/long",1,0);
    if (!pair_coul) pair_coul = force->pair_match("coul/msm",1,0);
    if (!pair_coul) error->all(FLERR,"Fix mdi/qmmm potential requires Coulomb-only pair sub-style");
  }

  // send info to QM code which is (possibly) only needed once
  // box size (can be resent)
  // DIRECT and POTENTIAL: number of QM atoms, QM elements
  // DIRECT only: number of MM atoms, MM elements, MM charges
  // NOTE: test if QM code alreads has NATOMS, which implies it set up itself ?
  
  if (first_send) {
    first_send = 0;

    send_box();

    int ierr;

    ierr = MDI_Send_command(">NATOMS", mdicomm);
    if (ierr) error->all(FLERR, "MDI: >NATOMS command");
    ierr = MDI_Send(&nqm, 1, MDI_INT, mdicomm);
    if (ierr) error->all(FLERR, "MDI: >NATOMS data");

    set_qm2owned();
    set_eqm();
    
    ierr = MDI_Send_command(">ELEMENTS", mdicomm);
    if (ierr) error->all(FLERR, "MDI: >ELEMENTS command");
    ierr = MDI_Send(eqm, nqm, MDI_INT, mdicomm);
    if (ierr) error->all(FLERR, "MDI: >ELEMENTS data");

    if (mode == DIRECT) {
      ierr = MDI_Send_command(">NLATTICE", mdicomm);
      if (ierr) error->all(FLERR, "MDI: >NLATTICE command");
      ierr = MDI_Send(&nmm, 1, MDI_DOUBLE, mdicomm);
      if (ierr) error->all(FLERR, "MDI: >NLATTICE data");

      set_mm2owned();
      set_emm();
      
      ierr = MDI_Send_command(">LATTICE_ELEMENTS", mdicomm);
      if (ierr) error->all(FLERR, "MDI: >LATTICE_ELEMENTS command");
      ierr = MDI_Send(&emm, nmm, MDI_INT, mdicomm);
      if (ierr) error->all(FLERR, "MDI: >LATTICE_ELEMENTS data");

      set_qmm();
      ierr = MDI_Send_command(">LATTICE", mdicomm);
      if (ierr) error->all(FLERR, "MDI: >LATTICE command");
      ierr = MDI_Send(&qmm, nmm, MDI_DOUBLE, mdicomm);
      if (ierr) error->all(FLERR, "MDI: >LATTICE data");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::setup_post_neighbor()
{
  post_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::post_neighbor()
{
  set_qm2owned();
  set_mm2owned();
}

/* ----------------------------------------------------------------------
   only called in POTENTIAL mode
   (1) calculate Coulomb potential for each QM atom
   (2) send info on QM atoms to QM code
   (3) invoke the QM solver
   (4) receive results from QM code
---------------------------------------------------------------------- */

void FixMDIQMMM::pre_force(int vflag)
{
  int ilocal,jlocal;
  double rsq;
  double delta[3];

  // invoke pair hybrid sub-style pair coulomb and Kspace directly
  // set eflag = 2 so they calculate per-atom energy

  pair_coul->compute(2,0);
  double *eatom_pair = pair_coul->eatom;

  double *eatom_kspace = nullptr;
  if (force->kspace) {
    force->kspace->compute(2,0);
    eatom_kspace = force->kspace->eatom;
  }

  // allocate ecoul for owned + ghost atoms
  // ghost atoms values only used if newton_pair is set

  if (atom->nmax > ncoulmax) {
    memory->destroy(ecoul);
    ncoulmax = atom->nmax;
    memory->create(ecoul,ncoulmax,"mdi/qmmm:ecoul");
  }

  // ecoul = per-atom energy for my owned atoms
  // if newton_pair, also do reverse_comm for ghost atom contribution

  int nlocal = atom->nlocal;
  int ntotal = nlocal;
  if (force->newton_pair) ntotal += atom->nghost;

  for (int i = 0; i < ntotal; i++)
    ecoul[i] = eatom_pair[i];

  if (force->newton_pair) comm->reverse_comm(this);

  if (force->kspace) {
    for (int i = 0; i < nlocal; i++)
      ecoul[i] += eatom_kspace[i];
  }

  // setup QM inputs: xqm and qpotential
  // xqm = atom coords, mapped into periodic box
  // qpotential[i] = Coulomb potential for each atom
  //   2 * (eatom[i] from pair_coul + kspace) / Qi
  //   factor of 2 comes from need to double count energy for each atom
  //   set for owned atoms, then MPI_Allreduce
  // subtract Qj/Rij energy for QM I interacting with all other QM J atoms
  //   use xqm_mine and qqm_mine for all QM atoms

  set_xqm();

  for (int i = 0; i < nqm; i++) qpotential_mine[i] = 0.0;

  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      if (q[ilocal] == 0.0) qpotential_mine[i] = 0.0;
      else qpotential_mine[i] = 2.0 * ecoul[ilocal] / q[ilocal];
    }
  }

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      for (int j = 0; j < nqm; j++) {
        if (j == i) continue;
        delta[0] = xqm[i][0] - xqm[j][0];
        delta[1] = xqm[i][1] - xqm[j][1];
        delta[2] = xqm[i][2] - xqm[j][2];
        domain->minimum_image_once(delta);
        rsq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
        qpotential_mine[i] -= qqrd2e * qqm[j] / sqrt(rsq);
      }
    }
  }

  MPI_Allreduce(qpotential_mine,qpotential,nqm,MPI_DOUBLE,MPI_SUM,world);

  // unit conversion from LAMMPS to MDI

  for (int i = 0; i < nqm; i++) {
    xqm[i][0] *= lmp2mdi_length;
    xqm[i][1] *= lmp2mdi_length;
    xqm[i][2] *= lmp2mdi_length;
    qpotential[i] *= lmp2mdi_energy;
  }

  // send info to MDI engine with QM atom info
  // first request for results triggers QM calculation
  // QM atoms must be in order of ascending atom ID
  // inputs:
  //   xqm = atom coords
  //   qpotential = vector of zeroes for AIMD
  // outputs:
  //   fqm = forces on QM atoms
  //   qqm = new charges on QM atoms
  //   qm_energy = QM contribution to energy of entire system
  
  if (comm->me == 0) utils::logmesg(lmp, "Invoking QM code ...\n");

  MPI_Barrier(world);
  double tstart = platform::walltime();

  int ierr;
  
  // send current coords of QM atoms to MDI engine

  ierr = MDI_Send_command(">COORDS", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >COORDS command");
  ierr = MDI_Send(&xqm[0][0], 3 * nqm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >COORDS data");
  
  // send Coulomb potential of QM atoms to MDI engine

  ierr = MDI_Send_command(">POTENTIAL_AT_NUCLEI", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >POTENTIAL_AT_NUCLEI command");
  ierr = MDI_Send(qpotential, nqm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >POTENTIAL_AT_NUCLEI data");

  // request QM potential energy from MDI engine
  // this triggers engine to perform QM calculation
  // qm_energy = fix output for global QM energy

  ierr = MDI_Send_command("<PE", mdicomm);
  if (ierr) error->all(FLERR, "MDI: <PE command");
  ierr = MDI_Recv(&qm_energy, 1, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <PE data");
  MPI_Bcast(&qm_energy, 1, MPI_DOUBLE, 0, world);

  // request forces on QM atoms from MDI engine

  ierr = MDI_Send_command("<FORCES", mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES command");
  ierr = MDI_Recv(&fqm[0][0], 3 * nqm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES data");
  MPI_Bcast(&fqm[0][0], 3 * nqm, MPI_DOUBLE, 0, world);

  // request charges on QM atoms from MDI engine

  ierr = MDI_Send_command("<CHARGES", mdicomm);
  if (ierr) error->all(FLERR, "MDI: <CHARGES command");
  ierr = MDI_Recv(qqm, nqm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <CHARGES data");
  MPI_Bcast(qqm, nqm, MPI_DOUBLE, 0, world);

  // end of MDI calls
  
  MPI_Barrier(world);
  if (comm->me == 0) 
    utils::logmesg(lmp, "  time = {:.3f} seconds\n",
                   platform::walltime() - tstart);

  // unit conversion from MDI to LAMMPS

  qm_energy *= mdi2lmp_energy;

  for (int i = 0; i < nqm; i++) {
    fqm[i][0] *= mdi2lmp_force;
    fqm[i][1] *= mdi2lmp_force;
    fqm[i][2] *= mdi2lmp_force;
  }

  // reset owned charges to QM values
  // communicate changes to ghost atoms

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) q[ilocal] = qqm[i];
  }

  comm->forward_comm(this);

  // reset LAMMPS forces to zero
  // NOTE: what about check in force_clear() for external_force_clear = OPENMP ?
  // NOTE: what will whichflag be for single snapshot compute of QMMM forces ?

  if (update->whichflag == 1) 
    update->integrate->force_clear();  
  else if (update->whichflag == 2) 
    update->minimize->force_clear();  
}

/* ----------------------------------------------------------------------
   different methods needed for DIRECT vs POTENTIAL mode
---------------------------------------------------------------------- */

void FixMDIQMMM::post_force(int vflag)
{
  if (mode == DIRECT) post_force_direct(vflag);
  else if (mode == POTENTIAL) post_force_potential(vflag);
}

/* ----------------------------------------------------------------------
   only called in DIRECT mode
   (1) send info on QM atoms to QM code
   (2) invoke the QM solver
   (3) receive results from QM code
---------------------------------------------------------------------- */

void FixMDIQMMM::post_force_direct(int vflag)
{
  // setup QM inputs:
  //   xqm = atom coords
  // setup MM inputs:
  //   xmm = atom coords
  //   qmm = charges on MM atoms
  
  set_xqm();
  set_xmm();
  set_qmm();

  // send info to MDI engine with QM and MM atom info
  // first request for results triggers QM calculation
  // QM and MM atoms must be in order of ascending atom ID
  // inputs:
  //   xqm = QM atom coords, mapped into periodic box
  //   xmm = MM atom coords, mapped into periodic box
  //   qmm = MM atom charges
  // outputs:
  //   qm_energy = QM contribution to energy of entire system
  //   fqm = forces on QM atoms
  //   fmm = forces on MM atoms
  
  if (comm->me == 0) utils::logmesg(lmp, "Invoking QM code ...\n");

  MPI_Barrier(world);
  double tstart = platform::walltime();

  // MDI calls
  // NOTE: add send_box() call when needed for NPT

  int ierr;
  
  // send current coords of QM atoms to MDI engine
  
  ierr = MDI_Send_command(">COORDS", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >COORDS command");
  ierr = MDI_Send(&xqm[0][0], 3 * nqm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >COORDS data");

  // send current coords of MM atoms to MDI engine
  
  ierr = MDI_Send_command(">CLATTICE", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CLATTICE command");
  ierr = MDI_Send(&xmm[0][0], 3 * nmm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CLATTICE data");

  // request QM potential energy from MDI engine
  // this triggers engine to perform QM calculation
  // qm_energy = fix output for global QM energy

  ierr = MDI_Send_command("<PE", mdicomm);
  if (ierr) error->all(FLERR, "MDI: <PE command");
  ierr = MDI_Recv(&qm_energy, 1, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <PE data");
  MPI_Bcast(&qm_energy, 1, MPI_DOUBLE, 0, world);

  // request forces on QM atoms from MDI engine
  
  ierr = MDI_Send_command("<FORCES", mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES command");
  ierr = MDI_Recv(&fqm[0][0], 3 * nqm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES data");
  MPI_Bcast(&fqm[0][0], 3 * nqm, MPI_DOUBLE, 0, world);

  // request forces on MM atoms from MDI engine
  
  ierr = MDI_Send_command("<FORCES_LATTICE", mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES_LATTICE command");
  ierr = MDI_Recv(&fmm[0][0], 3 * nmm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES_LATTICE data");
  MPI_Bcast(&fmm[0][0], 3 * nmm, MPI_DOUBLE, 0, world);

  // end of MDI calls
  
  MPI_Barrier(world);
  if (comm->me == 0) 
    utils::logmesg(lmp, "  time = {:.3f} seconds\n",
                   platform::walltime() - tstart);

  // unit conversion from MDI to LAMMPS

  qm_energy *= mdi2lmp_energy;

  for (int i = 0; i < nqm; i++) {
    fqm[i][0] *= mdi2lmp_force;
    fqm[i][1] *= mdi2lmp_force;
    fqm[i][2] *= mdi2lmp_force;
  }

  for (int i = 0; i < nmm; i++) {
    fmm[i][0] *= mdi2lmp_force;
    fmm[i][1] *= mdi2lmp_force;
    fmm[i][2] *= mdi2lmp_force;
  }

  // add fqm and fmm to LAMMPS forces on respective atoms

  double **f = atom->f;
  int ilocal;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      f[ilocal][0] += fqm[i][0];
      f[ilocal][1] += fqm[i][1];
      f[ilocal][2] += fqm[i][2];
    }
  }

  for (int i = 0; i < nmm; i++) {
    ilocal = mm2owned[i];
    if (ilocal >= 0) {
      f[ilocal][0] += fmm[i][0];
      f[ilocal][1] += fmm[i][1];
      f[ilocal][2] += fmm[i][2];
    }
  }
}

/* ----------------------------------------------------------------------
   only called in POTENTIAL mode
   add QM forces to QM atoms
   called after LAMMPS re-computes all MM forces with new QM charges
---------------------------------------------------------------------- */

void FixMDIQMMM::post_force_potential(int vflag)
{
  // int ilocal,jlocal;
  // double rsq,r2inv,rinv,fpair;
  // double delta[3];

  // subtract pairwise QM energy and forces from energy and forces
  //   LAMMPS just computed for all atoms
  // double loop over all QM pairs, 2nd loop starts at i+1 to include pair once
  // if I own either or both atoms in pair, compute pairwise Coulomb term
  //   subtract force only from owned atoms
  //   subtract half or all energy from qmenergy
  //   this effectively subtract energy from total = pair + kspace + fix

  // NOTE: some codes like NWChem may perform this operation themselves
  //       need to have a fix mdi/qmmm option for this, or different mode ?
  
  /*
  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int nlocal = atom->nlocal;
  double qqrd2e = force->qqrd2e;

  double eqm_mine = 0.0;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    for (int j = i+1; j < nqm; j++) {
      jlocal = qm2owned[j];

      // skip if neither atom is owned

      if (ilocal < 0 && jlocal < 0) continue;
      
      delta[0] = xqm[i][0] - xqm[j][0];
      delta[1] = xqm[i][1] - xqm[j][1];
      delta[2] = xqm[i][2] - xqm[j][2];
      domain->minimum_image_once(delta);
      rsq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
      r2inv = 1.0/rsq;
      rinv = sqrt(r2inv);
      fpair = qqrd2e * qqm[i]*qqm[j]*rinv*r2inv;

      // adjust forces on ilocal and/or jlocal only if they are owned atoms

      if (ilocal >= 0) {
        f[ilocal][0] -= delta[0]*fpair;
        f[ilocal][1] -= delta[1]*fpair;
        f[ilocal][2] -= delta[2]*fpair;
      }
      if (jlocal >= 0) {
        f[jlocal][0] += delta[0]*fpair;
        f[jlocal][1] += delta[1]*fpair;
        f[jlocal][2] += delta[2]*fpair;
      }

      // adjust energy using efactor
      // efactor = 1.0 if both are owned atoms, 0.5 if only one is owned

      double efactor = 0.5;
      if (ilocal >= 0 && jlocal >= 0.0) efactor = 1.0;
      eqm_mine += efactor * qqrd2e * qqm[i]*qqm[j]*rinv;
    }
  }

  // sum eqm_mine across procs, use it to adjust qmenergy
  
  double eqm;
  MPI_Allreduce(&eqm_mine,&eqm,1,MPI_DOUBLE,MPI_SUM,world);

  qmenergy -= eqm;
  */

  // add previously requested QM forces to owned QM atoms
  // do this now, after LAMMPS forces have been re-computed with new QM charges

  double **f = atom->f;

  int ilocal;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      f[ilocal][0] += fqm[i][0];
      f[ilocal][1] += fqm[i][1];
      f[ilocal][2] += fqm[i][2];
    }
  }

  // trigger per-atom energy computation on next step by pair/kspace
  // NOTE: is this needed ?
  //       only if needed for this fix to calc per-atom forces
  //       or needed for this fix to output global (or per-atom) energy

  // c_pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::min_post_neighbor()
{
  post_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

int FixMDIQMMM::pack_forward_comm(int n, int *list, double *buf,
                                  int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  double *q = atom->q;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = q[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  double *q = atom->q;

  for (i = first; i < last; i++) q[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int FixMDIQMMM::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = ecoul[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    ecoul[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   energy from MDI engine
------------------------------------------------------------------------- */

double FixMDIQMMM::compute_scalar()
{
  return qm_energy;
}

/* ----------------------------------------------------------------------
   virial from MDI engine
------------------------------------------------------------------------- */

double FixMDIQMMM::compute_vector(int n)
{
  return qm_virial_symmetric[n];
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double FixMDIQMMM::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)ncoulmax * sizeof(double);
  bytes += (double)(3*3 + 4) * nqm * sizeof(double);  // fpoint QM arrays/vecs
  bytes += nqm * sizeof(tagint);    // qmIDs
  bytes += nqm * sizeof(int);       // qm2owned
  return bytes;
}

// ----------------------------------------------------------------------
// private methods for this fix
// ----------------------------------------------------------------------

void FixMDIQMMM::create_qm_list()
{
  // qmIDs_mine = list of nqm_mine QM atom IDs I own
  // qmIDs = IDs of all QM atoms in ascending order
  // qmIDs created by allgather of qmIDs_mine

  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int nqm_mine = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) nqm_mine++;
  
  tagint *qmIDs_mine;
  memory->create(qmIDs_mine,nqm_mine,"mdi/qmmm:qmIDs_mine");
  
  nqm_mine = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) qmIDs_mine[nqm_mine++] = tag[i];

  int *recvcounts,*displs,*listall;
  memory->create(recvcounts,nprocs,"mdi/qmmm:recvcounts");
  memory->create(displs,nprocs,"mdi/qmmm:displs");

  MPI_Allgather(&nqm_mine,1,MPI_INT,recvcounts,1,MPI_INT,world);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++)
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];
  
  MPI_Allgatherv(qmIDs_mine,nqm_mine,MPI_LMP_TAGINT,qmIDs,recvcounts,displs,
                 MPI_LMP_TAGINT,world);
  
  memory->destroy(qmIDs_mine);
  memory->destroy(recvcounts);
  memory->destroy(displs);

  // sort qmIDs via merge sort

  int *order;
  tagint *qmIDs_sort;

  memory->create(order,nqm,"mdi/qmmm:order");
  memory->create(qmIDs_sort,nqm,"mdi/qmmm:qmIDs_sort");

  for (int i = 0; i < nqm; i++) {
    qmIDs_sort[i] = qmIDs[i];
    order[i] = i;
  }

  utils::merge_sort(order,nqm,(void *) qmIDs_sort,compare_IDs);

  int j;
  for (int i = 0; i < nqm; i++) {
    j = order[i];
    qmIDs_sort[i] = qmIDs[j];
  }

  memcpy(qmIDs,qmIDs_sort,nqm*sizeof(tagint));

  memory->destroy(order);
  memory->destroy(qmIDs_sort);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::create_mm_list()
{
  // mmIDs_mine = list of nmm_mine MM atom IDs I own
  // mmIDs = IDs of all MM atoms in ascending order
  // mmIDs created by allgather of mmIDs_mine

  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int nmm_mine = 0;
  for (int i = 0; i < nlocal; i++)
  if (!(mask[i] & groupbit)) nmm_mine++;
  
  tagint *mmIDs_mine;
  memory->create(mmIDs_mine,nmm_mine,"mdi/qmmm:mmIDs_mine");
  
  nmm_mine = 0;
  for (int i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) mmIDs_mine[nmm_mine++] = tag[i];

  int *recvcounts,*displs,*listall;
  memory->create(recvcounts,nprocs,"mdi/qmmm:recvcounts");
  memory->create(displs,nprocs,"mdi/qmmm:displs");

  MPI_Allgather(&nmm_mine,1,MPI_INT,recvcounts,1,MPI_INT,world);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++)
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];
  
  MPI_Allgatherv(mmIDs_mine,nmm_mine,MPI_LMP_TAGINT,mmIDs,recvcounts,displs,
                 MPI_LMP_TAGINT,world);
  
  memory->destroy(mmIDs_mine);
  memory->destroy(recvcounts);
  memory->destroy(displs);

  // sort mmIDs via merge sort

  int *order;
  tagint *mmIDs_sort;

  memory->create(order,nmm,"mdi/qmmm:order");
  memory->create(mmIDs_sort,nmm,"mdi/qmmm:mmIDs_sort");

  for (int i = 0; i < nmm; i++) {
    mmIDs_sort[i] = mmIDs[i];
    order[i] = i;
  }

  utils::merge_sort(order,nmm,(void *) mmIDs_sort,compare_IDs);

  int j;
  for (int i = 0; i < nmm; i++) {
    j = order[i];
    mmIDs_sort[i] = mmIDs[j];
  }

  memcpy(mmIDs,mmIDs_sort,nmm*sizeof(tagint));

  memory->destroy(order);
  memory->destroy(mmIDs_sort);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_qm2owned()
{
  // qm2owned[i] = index of local atom for each of nqm QM atoms
  // IDs of QM atoms are stored in qmIDs
  // index = -1 if this proc does not own the atom

  int nlocal = atom->nlocal;
  int index;

  for (int i = 0; i < nqm; i++) {
    index = atom->map(qmIDs[i]);
    if (index >= nlocal) qm2owned[i] = -1;
    else qm2owned[i] = index;
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_mm2owned()
{
  // mm2owned[i] = index of local atom for each of nmm MM atoms
  // IDs of MM atoms are stored in mmIDs
  // index = -1 if this proc does not own the atom

  int nlocal = atom->nlocal;
  int index;

  for (int i = 0; i < nmm; i++) {
    index = atom->map(mmIDs[i]);
    if (index >= nlocal) mm2owned[i] = -1;
    else mm2owned[i] = index;
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_qqm()
{
  for (int i = 0; i < nqm; i++) qqm_mine[i] = 0.0;

  double *q = atom->q;
  int ilocal;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) qqm_mine[i] = q[ilocal];
  }
  
  MPI_Allreduce(qqm_mine,qqm,nqm,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_eqm()
{
  for (int i = 0; i < nqm; i++) eqm_mine[i] = 0;

  int *type = atom->type;
  int ilocal;
  
  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) eqm_mine[i] = elements[type[ilocal]];
  }
  
  MPI_Allreduce(eqm_mine,eqm,nqm,MPI_INT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_xqm()
{
  for (int i = 0; i < nqm; i++) {
    xqm_mine[i][0] = 0.0;
    xqm_mine[i][1] = 0.0;
    xqm_mine[i][2] = 0.0;
  }

  double **x = atom->x;
  int ilocal;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      xqm_mine[i][0] = x[ilocal][0];
      xqm_mine[i][1] = x[ilocal][1];
      xqm_mine[i][2] = x[ilocal][2];

      domain->remap(xqm_mine[i]);
    }
  }

  MPI_Allreduce(&xqm_mine[0][0],&xqm[0][0],3*nqm,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_emm()
{
  for (int i = 0; i < nmm; i++) emm_mine[i] = 0;

  int *type = atom->type;
  int ilocal;
  
  for (int i = 0; i < nmm; i++) {
    ilocal = mm2owned[i];
    if (ilocal >= 0) emm_mine[i] = elements[type[ilocal]];
  }
  
  MPI_Allreduce(emm_mine,emm,nmm,MPI_INT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_qmm()
{
  for (int i = 0; i < nmm; i++) qmm_mine[i] = 0.0;

  double *q = atom->q;
  int ilocal;

  for (int i = 0; i < nmm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) qmm_mine[i] = q[ilocal];
  }
  
  MPI_Allreduce(qmm_mine,qmm,nmm,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_xmm()
{
  for (int i = 0; i < nmm; i++) {
    xmm_mine[i][0] = 0.0;
    xmm_mine[i][1] = 0.0;
    xmm_mine[i][2] = 0.0;
  }

  double **x = atom->x;
  int ilocal;

  for (int i = 0; i < nmm; i++) {
    ilocal = mm2owned[i];
    if (ilocal >= 0) {
      xmm_mine[i][0] = x[ilocal][0];
      xmm_mine[i][1] = x[ilocal][1];
      xmm_mine[i][2] = x[ilocal][2];

      domain->remap(xmm_mine[i]);
    }
  }

  MPI_Allreduce(&xmm_mine[0][0],&xmm[0][0],3*nmm,MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   send simulation box size and shape to MDI engine
------------------------------------------------------------------------- */

void FixMDIQMMM::send_box()
{
  double cell[9];

  int cell_displ_exists;
  int ierr = MDI_Check_command_exists("@DEFAULT", ">CELL_DISPL", mdicomm, &cell_displ_exists);
  if (ierr) error->all(FLERR, "MDI: >CELL_DISPL command check");
  MPI_Bcast(&cell_displ_exists, 1, MPI_INT, 0, world);

  if (cell_displ_exists) {
    ierr = MDI_Send_command(">CELL_DISPL", mdicomm);
    if (ierr) error->all(FLERR, "MDI: >CELL_DISPL command");
    cell[0] = domain->boxlo[0] * lmp2mdi_length;
    cell[1] = domain->boxlo[1] * lmp2mdi_length;
    cell[2] = domain->boxlo[2] * lmp2mdi_length;
    ierr = MDI_Send(cell, 3, MDI_DOUBLE, mdicomm);
    if (ierr) error->all(FLERR, "MDI: >CELL_DISPL data");
  }

  ierr = MDI_Send_command(">CELL", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CELL command");
  cell[0] = domain->boxhi[0] - domain->boxlo[0];
  cell[1] = 0.0;
  cell[2] = 0.0;
  cell[3] = domain->xy;
  cell[4] = domain->boxhi[1] - domain->boxlo[1];
  cell[5] = 0.0;
  cell[6] = domain->xz;
  cell[7] = domain->yz;
  cell[8] = domain->boxhi[2] - domain->boxlo[2];

  // convert from LAMMPS to MDI units
  
  for (int icell = 0; icell < 9; icell++) cell[icell] *= lmp2mdi_length;

  ierr = MDI_Send(cell, 9, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CELL data");
}

/* ----------------------------------------------------------------------
   MDI to/from LAMMPS conversion factors
------------------------------------------------------------------------- */

void FixMDIQMMM::unit_conversions()
{
  double angstrom_to_bohr, kelvin_to_hartree, ev_to_hartree, second_to_aut;

  MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
  MDI_Conversion_factor("kelvin_energy", "hartree", &kelvin_to_hartree);
  MDI_Conversion_factor("electron_volt", "hartree", &ev_to_hartree);
  MDI_Conversion_Factor("second", "atomic_unit_of_time", &second_to_aut);

  // length units

  mdi2lmp_length = 1.0;
  lmp2mdi_length = 1.0;

  if (lmpunits == REAL || lmpunits == METAL) {
    lmp2mdi_length = angstrom_to_bohr;
    mdi2lmp_length = 1.0 / angstrom_to_bohr;
  }

  // energy units

  mdi2lmp_energy = 1.0;
  lmp2mdi_energy = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_energy = kelvin_to_hartree / force->boltz;
    mdi2lmp_energy = force->boltz / kelvin_to_hartree;
  } else if (lmpunits == METAL) {
    lmp2mdi_energy = ev_to_hartree;
    mdi2lmp_energy = 1.0 / ev_to_hartree;
  }

  // force units = energy/length

  mdi2lmp_force = 1.0;
  lmp2mdi_force = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_force = (kelvin_to_hartree / force->boltz) / angstrom_to_bohr;
    mdi2lmp_force = 1.0 / lmp2mdi_force;
  } else if (lmpunits == METAL) {
    lmp2mdi_force = ev_to_hartree / angstrom_to_bohr;
    mdi2lmp_force = angstrom_to_bohr / ev_to_hartree;
  }

  // stress units = force/area = energy/volume

  mdi2lmp_pressure = 1.0;
  lmp2mdi_pressure = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_pressure = (kelvin_to_hartree / force->boltz) /
        (angstrom_to_bohr * angstrom_to_bohr * angstrom_to_bohr);
    mdi2lmp_pressure = 1.0 / lmp2mdi_pressure;
  } else if (lmpunits == METAL) {
    lmp2mdi_pressure = ev_to_hartree / (angstrom_to_bohr * angstrom_to_bohr * angstrom_to_bohr);
    mdi2lmp_pressure = 1.0 / lmp2mdi_pressure;
  }
}

/* ----------------------------------------------------------------------
   comparison function invoked by merge_sort()
   void pointer contains list of atom IDs
------------------------------------------------------------------------- */

int compare_IDs(const int i, const int j, void *ptr)
{
  tagint *ids = (int *) ptr;
  if (ids[i] < ids[j]) return -1;
  if (ids[i] > ids[j]) return 1;
  return 0;
}
