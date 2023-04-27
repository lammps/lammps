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

#include "fix_mdi_qm.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NATIVE, REAL, METAL };    // LAMMPS units which MDI supports

#define MAXELEMENT 118

// prototype for non-class compare function for sorting QM IDs

static int compare_IDs(const int, const int, void *);

/* ---------------------------------------------------------------------- */

FixMDIQM::FixMDIQM(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  // check requirements for LAMMPS to work with MDI as an engine
  // atom IDs do not need to be consecutive

  if (atom->tag_enable == 0) error->all(FLERR, "Cannot use MDI engine without atom IDs");

  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Fix mdi/qm requires an atom map be defined");

  // confirm LAMMPS is being run as a driver

  int role;
  MDI_Get_role(&role);
  if (role != MDI_DRIVER)
    error->all(FLERR, "Must invoke LAMMPS as an MDI driver to use fix mdi/qm");

  // optional args

  virialflag = 0;
  addflag = 1;
  every = 1;
  connectflag = 1;
  elements = nullptr;
  mcflag = 0;
  id_mcfix = nullptr;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "virial") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix mdi/qm command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        virialflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        virialflag = 0;
      else
        error->all(FLERR, "Illegal fix mdi/qm command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "add") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix mdi/qm command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        addflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        addflag = 0;
      else
        error->all(FLERR, "Illegal fix mdi/qm command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "every") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix mdi/qm command");
      every = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (every <= 0) error->all(FLERR, "Illegal fix mdi/qm command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "connect") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix mdi/qm command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        connectflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        connectflag = 0;
      else
        error->all(FLERR, "Illegal fix mdi/qm command");
      iarg += 2;

    } else if (strcmp(arg[iarg], "elements") == 0) {
      const char *symbols[] = {
          "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",
          "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
          "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
          "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
          "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re",
          "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
          "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
          "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
      };

      int ntypes = atom->ntypes;
      if (iarg + ntypes + 1 > narg) error->all(FLERR, "Illegal fix mdi/qm command");
      delete[] elements;
      elements = new int[ntypes + 1];
      for (int i = 1; i <= ntypes; i++) {
        int anum;
        for (anum = 0; anum < MAXELEMENT; anum++)
          if (strcmp(arg[iarg + i], symbols[anum]) == 0) break;
        if (anum == MAXELEMENT) error->all(FLERR, "Invalid chemical element in fix mdi/qm command");
        elements[i] = anum + 1;
      }
      iarg += ntypes + 1;

    } else if (strcmp(arg[iarg], "mc") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix mdi/qm command");
      mcflag = 1;
      delete[] id_mcfix;
      id_mcfix = utils::strdup(arg[iarg + 1]);
      iarg += 2;

    } else
      error->all(FLERR, "Illegal fix mdi/qm command");
  }

  // fix output settings are based on optional keywords

  scalar_flag = 1;
  global_freq = every;
  extscalar = 1;

  if (virialflag) {
    vector_flag = 1;
    size_vector = 6;
    extvector = 0;
  }

  peratom_flag = 1;
  size_peratom_cols = 3;
  peratom_freq = every;

  if (addflag) {
    energy_global_flag = 1;
    virial_global_flag = 1;
    thermo_energy = thermo_virial = 1;
  }

  // mdicomm will be initialized in init()
  // cannot do here for a plugin library, b/c mdi plugin command comes later

  mdicomm = MDI_COMM_NULL;

  // set unit conversion factors

  if (strcmp(update->unit_style, "real") == 0)
    lmpunits = REAL;
  else if (strcmp(update->unit_style, "metal") == 0)
    lmpunits = METAL;
  else
    lmpunits = NATIVE;

  unit_conversions();

  nprocs = comm->nprocs;

  // QM atom data

  nqm = nqm_last = max_nqm = 0;
  nexclude = 0;

  qmIDs = nullptr;
  qm2owned = nullptr;

  eqm = nullptr;
  tqm = nullptr;
  xqm = nullptr;
  fqm = nullptr;

  eqm_mine = nullptr;
  tqm_mine = nullptr;
  xqm_mine = nullptr;

  // per-atom data

  nmax = atom->nmax;
  memory->create(array_atom, nmax, 3, "mdi/qm:array_atom");

  // initialize outputs

  qm_energy = 0.0;
  for (int i = 0; i < 6; i++) {
    qm_virial[i] = 0.0;
    virial[i] = 0.0;
  }
  sumflag = 0;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) array_atom[i][0] = array_atom[i][1] = array_atom[i][2] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixMDIQM::~FixMDIQM()
{
  // send exit command to stand-alone engine code
  // for connnectflag = 0, this is done via "mdi exit" command
  // for plugin, this is done in MDIPlugin::plugin_wrapper()

  if (mdicomm != MDI_COMM_NULL && connectflag && !plugin) {
    int ierr = MDI_Send_command("EXIT", mdicomm);
    if (ierr) error->all(FLERR, "MDI: EXIT command");
  }

  // clean up

  delete[] id_mcfix;
  delete[] elements;

  memory->destroy(qmIDs);
  memory->destroy(qm2owned);

  memory->destroy(eqm);
  memory->destroy(tqm);
  memory->destroy(xqm);
  memory->destroy(fqm);

  memory->destroy(eqm_mine);
  memory->destroy(tqm_mine);
  memory->destroy(xqm_mine);

  memory->destroy(array_atom);
}

/* ---------------------------------------------------------------------- */

int FixMDIQM::setmask()
{
  int mask = 0;
  mask |= POST_NEIGHBOR;
  mask |= MIN_POST_NEIGHBOR;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::init()
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
        error->all(FLERR, "Fix mdi/qm is not connected to engine via mdi connect");

      int nbytes = sizeof(MDI_Comm);
      char *ptrcomm = (char *) lmp->mdicomm;
      memcpy(&mdicomm, ptrcomm, nbytes);
    }

    // check which MDI commands engine supports

    int ierr = MDI_Check_command_exists("@DEFAULT", ">NATOMS", mdicomm, &natoms_exists);
    if (ierr) error->all(FLERR, "MDI: >NATOMS command check");
    MPI_Bcast(&natoms_exists, 1, MPI_INT, 0, world);

    ierr = MDI_Check_command_exists("@DEFAULT", ">CELL_DISPL", mdicomm, &celldispl_exists);
    if (ierr) error->all(FLERR, "MDI: >CELL_DISPL command check");
    MPI_Bcast(&celldispl_exists, 1, MPI_INT, 0, world);

    ierr = MDI_Check_command_exists("@DEFAULT", ">ELEMENTS", mdicomm, &elements_exists);
    if (ierr) error->all(FLERR, "MDI: >ELEMENTS command check");
    MPI_Bcast(&elements_exists, 1, MPI_INT, 0, world);

    ierr = MDI_Check_command_exists("@DEFAULT", ">TYPES", mdicomm, &types_exists);
    if (ierr) error->all(FLERR, "MDI: >TYPES command check");
    MPI_Bcast(&types_exists, 1, MPI_INT, 0, world);

    ierr = MDI_Check_command_exists("@DEFAULT", "<STRESS", mdicomm, &stress_exists);
    if (ierr) error->all(FLERR, "MDI: <STRESS command check");
    MPI_Bcast(&stress_exists, 1, MPI_INT, 0, world);

    ierr = MDI_Check_command_exists("@DEFAULT", "<PE", mdicomm, &pe_exists);
    if (ierr) error->all(FLERR, "MDI: <PE command check");
    MPI_Bcast(&pe_exists, 1, MPI_INT, 0, world);

    ierr = MDI_Check_command_exists("@DEFAULT", "<KE_ELEC", mdicomm, &keelec_exists);
    if (ierr) error->all(FLERR, "MDI: <KE_ELEC command check");
    MPI_Bcast(&keelec_exists, 1, MPI_INT, 0, world);
  }

  // extract pointers to MC variables from id_mcfix
  // mc_active = 1/0 for MC moves on/off (required)
  // exclusion_group = index of a group the Fix defines to exclude atoms (optional)

  if (mcflag) {
    Fix *f_mc = modify->get_fix_by_id(id_mcfix);
    if (!f_mc) error->all(FLERR, "Fix mdi/qm could not find Monte Carlo fix ID {}", id_mcfix);
    int dim;
    mc_active_ptr = (int *) f_mc->extract("mc_active", dim);
    if (!mc_active_ptr || dim != 0)
      error->all(FLERR, "Fix mdi/qm could not query mc_active from Monte Carlo fix ID {}",
                 id_mcfix);
    exclusion_group_ptr = (int *) f_mc->extract("exclusion_group", dim);
  }

  // determine whether a new vs incremental QM calc is needed
  // new if first run or if
  //   atom count or elements/types or box has changed between runs
  // otherwise incremental = subsequent run of same system

  int new_system = 0;

  // check if count of QM atoms has changed
  // on first run, old count is 0

  int nqm_old = nqm;
  nqm = set_nqm();

  if (nqm != nqm_old) {
    if (nqm > max_nqm) reallocate();
    create_qm_list();
    set_qm2owned();
    new_system = 1;
  }

  // check if box has changed

  if (new_system)
    set_box();
  else {
    double old_cell[9], old_cell_displ[3];
    memcpy(old_cell, qm_cell, 9 * sizeof(double));
    memcpy(old_cell_displ, qm_cell_displ, 3 * sizeof(double));
    set_box();
    for (int i = 0; i < 9; i++)
      if (qm_cell[i] != old_cell[i]) new_system = 1;
    for (int i = 0; i < 3; i++)
      if (qm_cell_displ[i] != old_cell_displ[i]) new_system = 1;
  }

  // check if atom elements or types have changed

  if (elements && elements_exists) {
    if (new_system)
      set_eqm();
    else {
      int *eqm_old;
      memory->create(eqm_old, nqm, "mdi/qm:eqm_old");
      memcpy(eqm_old, eqm, nqm * sizeof(int));
      set_eqm();
      for (int i = 0; i < nqm; i++)
        if (eqm[i] != eqm_old[i]) new_system = 1;
      memory->destroy(eqm_old);
    }

  } else if (types_exists) {
    if (new_system)
      set_tqm();
    else {
      int *tqm_old;
      memory->create(tqm_old, nqm, "mdi/qm:tqm_old");
      memcpy(tqm_old, tqm, nqm * sizeof(int));
      set_tqm();
      for (int i = 0; i < nqm; i++)
        if (tqm[i] != tqm_old[i]) new_system = 1;
      memory->destroy(tqm_old);
    }
  }

  // if new system, send setup info to MDI engine
  // values that often won't change for AIMD simulations
  // if not sending elements or types, assume engine initialized itself

  if (new_system) {
    send_natoms();
    send_box();
    if (elements && elements_exists)
      send_elements();
    else if (types_exists)
      send_types();
    nqm_last = nqm;
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::setup_post_neighbor()
{
  post_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::post_neighbor()
{
  set_qm2owned();
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::post_force(int vflag)
{
  int ierr;

  // skip if timestep is not a multiple of every

  if (update->ntimestep % every) return;

  // reallocate array_atom if needed

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(array_atom);
    memory->create(array_atom, nmax, 3, "mdi/qm:array_atom");
  }

  // determine whether a new vs incremental QM calc is needed
  // new when atom count has changed (deposit, evaporate)
  //   or when MC fix is active (insertion, deletion, large moves)
  // incremental when a system is slowly evolving (AIMD)

  int new_system = 0;
  if (nqm != nqm_last)
    new_system = 1;
  else if (mcflag && *mc_active_ptr)
    new_system = 1;

  // send new system info to MDI engine: atom count and elements/types
  // reset QM data structures if atom count has changed

  if (new_system) {
    nqm = set_nqm();
    if (nqm > max_nqm) reallocate();
    create_qm_list();
    set_qm2owned();

    send_natoms();
    set_box();
    send_box();
    if (elements && elements_exists) {
      set_eqm();
      send_elements();
    } else if (types_exists) {
      set_tqm();
      send_types();
    }

    nqm_last = nqm;

    // incremental system
    // if simulation box dynamically changes, send current box to MDI engine

  } else if (domain->box_change_size || domain->box_change_shape) {
    set_box();
    send_box();
  }

  // send current coords of QM atoms to MDI engine

  set_xqm();

  ierr = MDI_Send_command(">COORDS", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >COORDS command");
  ierr = MDI_Send(&xqm[0][0], 3 * nqm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >COORDS data");

  // request QM energy from MDI engine
  // this triggers engine to perform QM calculation
  // sets qm_energy = fix output for global QM energy

  request_qm_energy();

  // request forces from MDI engine

  ierr = MDI_Send_command("<FORCES", mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES command");
  ierr = MDI_Recv(&fqm[0][0], 3 * nqm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES data");
  MPI_Bcast(&fqm[0][0], 3 * nqm, MPI_DOUBLE, 0, world);

  // request stress if needed and supported

  if (vflag && virialflag && stress_exists) {
    ierr = MDI_Send_command("<STRESS", mdicomm);
    if (ierr) error->all(FLERR, "MDI: <STRESS command");
    ierr = MDI_Recv(qm_virial, 9, MDI_DOUBLE, mdicomm);
    if (ierr) error->all(FLERR, "MDI: <STRESS data");
    MPI_Bcast(qm_virial, 9, MPI_DOUBLE, 0, world);
  }

  // optionally add QM forces to owned atoms

  if (addflag) {
    double **f = atom->f;
    int ilocal;
    for (int i = 0; i < nqm; i++) {
      ilocal = qm2owned[i];
      if (ilocal >= 0) {
        f[ilocal][0] += fqm[i][0] * mdi2lmp_force;
        f[ilocal][1] += fqm[i][1] * mdi2lmp_force;
        f[ilocal][2] += fqm[i][2] * mdi2lmp_force;
      }
    }
  }

  // array_atom = fix output for peratom QM forces
  // if nexclude, some atoms are not QM atoms, zero array_atom first

  if (nexclude) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) {
      array_atom[i][0] = 0.0;
      array_atom[i][1] = 0.0;
      array_atom[i][2] = 0.0;
    }
  }

  int ilocal;
  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      array_atom[ilocal][0] = fqm[i][0] * mdi2lmp_force;
      array_atom[ilocal][1] = fqm[i][1] * mdi2lmp_force;
      array_atom[ilocal][2] = fqm[i][2] * mdi2lmp_force;
    }
  }

  // qm_virial_symmetric = fix output for global QM virial
  // MDI defines virial tensor as intensive (divided by volume), LAMMPS does not

  if (vflag && virialflag && stress_exists) {
    qm_virial_symmetric[0] = qm_virial[0] * mdi2lmp_pressure;
    qm_virial_symmetric[1] = qm_virial[4] * mdi2lmp_pressure;
    qm_virial_symmetric[2] = qm_virial[8] * mdi2lmp_pressure;
    qm_virial_symmetric[3] = 0.5 * (qm_virial[1] + qm_virial[3]) * mdi2lmp_pressure;
    qm_virial_symmetric[4] = 0.5 * (qm_virial[2] + qm_virial[6]) * mdi2lmp_pressure;
    qm_virial_symmetric[5] = 0.5 * (qm_virial[5] + qm_virial[7]) * mdi2lmp_pressure;

    // optionally set fix->virial
    // multiply by volume to make it extensive
    //   divide by nprocs so each proc stores a portion
    // this is b/c ComputePressure expects this as input from a fix
    //   it will do an MPI_Allreduce and divide by volume

    if (addflag) {
      double volume;
      if (domain->dimension == 2)
        volume = domain->xprd * domain->yprd;
      else if (domain->dimension == 3)
        volume = domain->xprd * domain->yprd * domain->zprd;
      for (int i = 0; i < 6; i++) virial[i] = qm_virial_symmetric[i] * volume / nprocs;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::min_post_neighbor()
{
  post_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy from MDI engine
------------------------------------------------------------------------- */

double FixMDIQM::compute_scalar()
{
  return qm_energy;
}

/* ----------------------------------------------------------------------
   virial from MDI engine
------------------------------------------------------------------------- */

double FixMDIQM::compute_vector(int n)
{
  return qm_virial_symmetric[n];
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double FixMDIQM::memory_usage()
{
  double bytes = 0.0;
  bytes += nqm * sizeof(tagint);            // qmIDs
  bytes += nqm * sizeof(int);               // qm2owned
  bytes += 4 * nqm * sizeof(int);           // int QM arrays/vecs
  bytes += 3 * 3 * nqm * sizeof(double);    // double QM arrays/vecs
  return bytes;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// private methods for this fix
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   reallocate storage for QM atoms
------------------------------------------------------------------------- */

void FixMDIQM::reallocate()
{
  max_nqm = nqm;

  memory->destroy(qmIDs);
  memory->destroy(qm2owned);

  memory->destroy(eqm);
  memory->destroy(tqm);
  memory->destroy(xqm);
  memory->destroy(fqm);

  memory->destroy(eqm_mine);
  memory->destroy(tqm_mine);
  memory->destroy(xqm_mine);

  memory->create(qmIDs, max_nqm, "mdi/qm:qmIDs");
  memory->create(qm2owned, max_nqm, "mdi/qm:qm2owned");

  memory->create(eqm, max_nqm, "mdi/qm:eqm");
  memory->create(tqm, max_nqm, "mdi/qm:tqm");
  memory->create(xqm, max_nqm, 3, "mdi/qm:xqm");
  memory->create(fqm, max_nqm, 3, "mdi/qm:fqm");

  memory->create(eqm_mine, max_nqm, "mdi/qm:eqm_mine");
  memory->create(tqm_mine, max_nqm, "mdi/qm:tqm_mine");
  memory->create(xqm_mine, max_nqm, 3, "mdi/qm:xqm_mine");
}

/* ----------------------------------------------------------------------
   ncount = # of QM atoms
   can be less than all atoms if MC flag is set
   return ncount to set nqm
------------------------------------------------------------------------- */

int FixMDIQM::set_nqm()
{
  // require 3*nqm be a small INT, so can MPI_Allreduce xqm

  if (3 * atom->natoms > MAXSMALLINT) error->all(FLERR, "Fix mdi/qm has too many atoms");

  int ncount = atom->natoms;
  nexclude = 0;

  if (mcflag && exclusion_group_ptr) {
    int nexclude_mine = 0;

    int exclusion_group = *exclusion_group_ptr;
    if (exclusion_group) {
      int excludebit = group->bitmask[exclusion_group];

      int *mask = atom->mask;
      int nlocal = atom->nlocal;

      for (int i = 0; i < nlocal; i++)
        if (mask[i] & excludebit) nexclude_mine = 1;
    }

    MPI_Allreduce(&nexclude_mine, &nexclude, 1, MPI_INT, MPI_SUM, world);
  }

  ncount -= nexclude;
  return ncount;
}

/* ----------------------------------------------------------------------
   create sorted list of QM atom IDs
   ignore excluded atoms, e.g. by fix GCMC
------------------------------------------------------------------------- */

void FixMDIQM::create_qm_list()
{
  int excludebit;
  if (nexclude) {
    int exclusion_group = *exclusion_group_ptr;
    excludebit = group->bitmask[exclusion_group];
  }

  // qmIDs_mine = list of nqm_mine QM atom IDs I own
  // qmIDs = IDs of all QM atoms in ascending order
  // qmIDs created by allgather of qmIDs_mine

  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int nqm_mine = 0;
  for (int i = 0; i < nlocal; i++) {
    if (!nexclude)
      nqm_mine++;
    else if (!(mask[i] & excludebit))
      nqm_mine++;
  }

  tagint *qmIDs_mine;
  memory->create(qmIDs_mine, nqm_mine, "mdi/qm:qmIDs_mine");

  nqm_mine = 0;
  for (int i = 0; i < nlocal; i++) {
    if (!nexclude)
      qmIDs_mine[nqm_mine++] = tag[i];
    else if (!(mask[i] & excludebit))
      qmIDs_mine[nqm_mine++] = tag[i];
  }

  int *recvcounts, *displs;
  memory->create(recvcounts, nprocs, "mdi/qm:recvcounts");
  memory->create(displs, nprocs, "mdi/qm:displs");

  MPI_Allgather(&nqm_mine, 1, MPI_INT, recvcounts, 1, MPI_INT, world);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++)
    displs[iproc] = displs[iproc - 1] + recvcounts[iproc - 1];

  MPI_Allgatherv(qmIDs_mine, nqm_mine, MPI_LMP_TAGINT, qmIDs, recvcounts, displs, MPI_LMP_TAGINT,
                 world);

  memory->destroy(qmIDs_mine);
  memory->destroy(recvcounts);
  memory->destroy(displs);

  // sort qmIDs via merge sort

  int *order;
  tagint *qmIDs_sort;

  memory->create(order, nqm, "mdi/qm:order");
  memory->create(qmIDs_sort, nqm, "mdi/qm:qmIDs_sort");

  for (int i = 0; i < nqm; i++) {
    qmIDs_sort[i] = qmIDs[i];
    order[i] = i;
  }

  utils::merge_sort(order, nqm, (void *) qmIDs_sort, compare_IDs);

  int j;
  for (int i = 0; i < nqm; i++) {
    j = order[i];
    qmIDs_sort[i] = qmIDs[j];
  }

  memcpy(qmIDs, qmIDs_sort, nqm * sizeof(tagint));

  memory->destroy(order);
  memory->destroy(qmIDs_sort);
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::set_qm2owned()
{
  // qm2owned[i] = index of local atom for each of nqm QM atoms
  // IDs of QM atoms are stored in qmIDs
  // index = -1 if this proc does not own the atom

  int nlocal = atom->nlocal;
  int index;

  for (int i = 0; i < nqm; i++) {
    index = atom->map(qmIDs[i]);
    if (index >= nlocal)
      qm2owned[i] = -1;
    else
      qm2owned[i] = index;
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::set_box()
{
  qm_cell_displ[0] = domain->boxlo[0] * lmp2mdi_length;
  qm_cell_displ[1] = domain->boxlo[1] * lmp2mdi_length;
  qm_cell_displ[2] = domain->boxlo[2] * lmp2mdi_length;

  qm_cell[0] = domain->boxhi[0] - domain->boxlo[0];
  qm_cell[1] = 0.0;
  qm_cell[2] = 0.0;
  qm_cell[3] = domain->xy;
  qm_cell[4] = domain->boxhi[1] - domain->boxlo[1];
  qm_cell[5] = 0.0;
  qm_cell[6] = domain->xz;
  qm_cell[7] = domain->yz;
  qm_cell[8] = domain->boxhi[2] - domain->boxlo[2];

  // convert cell units to bohr

  for (int icell = 0; icell < 9; icell++) qm_cell[icell] *= lmp2mdi_length;
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::set_xqm()
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

      xqm_mine[i][0] *= lmp2mdi_length;
      xqm_mine[i][1] *= lmp2mdi_length;
      xqm_mine[i][2] *= lmp2mdi_length;
    }
  }

  MPI_Allreduce(&xqm_mine[0][0], &xqm[0][0], 3 * nqm, MPI_DOUBLE, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::set_eqm()
{
  for (int i = 0; i < nqm; i++) eqm_mine[i] = 0;

  int *type = atom->type;
  int ilocal;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) eqm_mine[i] = elements[type[ilocal]];
  }

  MPI_Allreduce(eqm_mine, eqm, nqm, MPI_INT, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void FixMDIQM::set_tqm()
{
  for (int i = 0; i < nqm; i++) tqm_mine[i] = 0.0;

  int *type = atom->type;
  int ilocal;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) tqm_mine[i] = type[ilocal];
  }

  MPI_Allreduce(tqm_mine, tqm, nqm, MPI_INT, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   send LAMMPS QM atom count to MDI engine
------------------------------------------------------------------------- */

void FixMDIQM::send_natoms()
{
  int ierr;

  // if engine suppports >NATOMS, send it
  // if not, require that engine be consistent with LAMMPS

  if (natoms_exists) {
    ierr = MDI_Send_command(">NATOMS", mdicomm);
    if (ierr) error->all(FLERR, "MDI: >NATOMS command");
    ierr = MDI_Send(&nqm, 1, MDI_INT, mdicomm);
    if (ierr) error->all(FLERR, "MDI: >NATOMS data");

  } else {
    ierr = MDI_Send_command("<NATOMS", mdicomm);
    if (ierr) error->all(FLERR, "MDI: <NATOMS command");
    int n;
    ierr = MDI_Recv(&n, 1, MDI_INT, mdicomm);
    if (ierr) error->all(FLERR, "MDI: <NATOMS data");
    MPI_Bcast(&n, 1, MPI_INT, 0, world);

    if (n != nqm)
      error->all(FLERR, "MDI: Engine has wrong atom count and does not support >NATOMS command");
  }
}

/* ----------------------------------------------------------------------
   send LAMMPS atom types to MDI engine
------------------------------------------------------------------------- */

void FixMDIQM::send_types()
{
  int ierr = MDI_Send_command(">TYPES", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >TYPES command");
  ierr = MDI_Send(tqm, nqm, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >TYPES data");
}

/* ----------------------------------------------------------------------
   send elements to MDI engine = atomic numbers for each type
------------------------------------------------------------------------- */

void FixMDIQM::send_elements()
{
  int ierr = MDI_Send_command(">ELEMENTS", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >ELEMENTS command");
  ierr = MDI_Send(eqm, nqm, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >ELEMENTS data");
}

/* ----------------------------------------------------------------------
   send simulation box size and shape to MDI engine
   only send CELL_DISPL if engine supports it
------------------------------------------------------------------------- */

void FixMDIQM::send_box()
{
  int ierr;

  // only send cell dimensions if fully periodic simulation

  if (domain->nonperiodic == 0) {
    if (celldispl_exists) {
      ierr = MDI_Send_command(">CELL_DISPL", mdicomm);
      if (ierr) error->all(FLERR, "MDI: >CELL_DISPL command");
      ierr = MDI_Send(qm_cell_displ, 3, MDI_DOUBLE, mdicomm);
      if (ierr) error->all(FLERR, "MDI: >CELL_DISPL data");
    }

    ierr = MDI_Send_command(">CELL", mdicomm);
    if (ierr) error->all(FLERR, "MDI: >CELL command");
    ierr = MDI_Send(qm_cell, 9, MDI_DOUBLE, mdicomm);
    if (ierr) error->all(FLERR, "MDI: >CELL data");

  } else if (domain->xperiodic == 1 || domain->yperiodic == 1 ||
             domain->zperiodic == 1) {
    error->all(FLERR,"MDI requires fully periodic or fully non-periodic system");
  }
}

/* ----------------------------------------------------------------------
   request QM energy from MDI engine
   set qm_energy = fix output for global QM energy
------------------------------------------------------------------------- */

void FixMDIQM::request_qm_energy()
{
  int ierr;

  // QM energy = <PE + <KE_ELEC or <ENERGY, depending on engine options

  if (pe_exists && keelec_exists) {
    int pe_energy, keelec_energy;

    ierr = MDI_Send_command("<PE", mdicomm);
    if (ierr) error->all(FLERR, "MDI: <PE command");
    ierr = MDI_Recv(&pe_energy, 1, MDI_DOUBLE, mdicomm);
    if (ierr) error->all(FLERR, "MDI: <PE data");

    ierr = MDI_Send_command("<KE_ELEC", mdicomm);
    if (ierr) error->all(FLERR, "MDI: <KE_ELEC command");
    ierr = MDI_Recv(&keelec_energy, 1, MDI_DOUBLE, mdicomm);
    if (ierr) error->all(FLERR, "MDI: <KE_ELEC data");

    qm_energy = pe_energy + keelec_energy;

  } else {
    ierr = MDI_Send_command("<ENERGY", mdicomm);
    if (ierr) error->all(FLERR, "MDI: <ENERGY command");
    ierr = MDI_Recv(&qm_energy, 1, MDI_DOUBLE, mdicomm);
    if (ierr) error->all(FLERR, "MDI: <ENERGY data");
  }

  MPI_Bcast(&qm_energy, 1, MPI_DOUBLE, 0, world);
  qm_energy *= mdi2lmp_energy;
}

/* ----------------------------------------------------------------------
   MDI to/from LAMMPS conversion factors
------------------------------------------------------------------------- */

void FixMDIQM::unit_conversions()
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
