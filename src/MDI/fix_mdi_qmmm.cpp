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

static constexpr int MAXELEMENT = 118;

// prototype for non-class compare function for sorting QM IDs

static int compare_IDs(const int, const int, void *);

/* ---------------------------------------------------------------------- */

FixMDIQMMM::FixMDIQMMM(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), elements(nullptr), pair_coul(nullptr), qmIDs(nullptr), qm2owned(nullptr),
    eqm(nullptr), eqm_mine(nullptr), tqm(nullptr), tqm_mine(nullptr), xqm(nullptr),
    xqm_mine(nullptr), qqm(nullptr), qqm_mine(nullptr), qpotential(nullptr),
    qpotential_mine(nullptr), fqm(nullptr), ecoul(nullptr)
{
  // check requirements for LAMMPS to work with MDI for a QMMM engine
  // atom IDs do not need to be consecutive

  if (!atom->tag_enable) error->all(FLERR, "Fix mdi/qmmm requires atom IDs be defined");

  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Fix mdi/qmmm requires an atom map be defined");

  // initialize class members

  plugin = 0;
  maxlocal = 0;
  natoms_exists = 0;
  celldispl_exists = 0;
  elements_exists = 0;
  types_exists = 0;
  stress_exists = 0;
  pe_exists = 0;
  keelec_exists = 0;

  // confirm LAMMPS is being run as a driver

  int role;
  MDI_Get_role(&role);
  if (role != MDI_DRIVER)
    error->all(FLERR, "Must invoke LAMMPS as an MDI driver to use fix mdi/qmmm");

  // mode arg

  if (strcmp(arg[3], "direct") == 0)
    mode = DIRECT;
  else if (strcmp(arg[3], "potential") == 0)
    mode = POTENTIAL;
  else
    error->all(FLERR, "Illegal fix mdi/qmmm command");

  // optional args

  virialflag = 0;
  connectflag = 1;

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
      if (iarg + ntypes + 1 > narg) error->all(FLERR, "Illegal fix mdi/qmmm command");
      delete[] elements;
      elements = new int[ntypes + 1];
      for (int i = 1; i <= ntypes; i++) {
        int anum;
        for (anum = 0; anum < MAXELEMENT; anum++)
          if (strcmp(arg[iarg + i], symbols[anum]) == 0) break;
        if (anum == MAXELEMENT)
          error->all(FLERR, "Invalid chemical element in fix mdi/qmmm command");
        elements[i] = anum + 1;
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
    extvector = 0;
  }

  peratom_flag = 1;
  size_peratom_cols = 3;
  peratom_freq = 1;

  energy_global_flag = 1;
  virial_global_flag = 1;
  thermo_energy = thermo_virial = 1;

  comm_forward = 1;
  comm_reverse = 1;

  // mdicomm will be initialized in init()
  // cannot do here for a plugin library, b/c mdi plugin command comes later

  mdicomm = MDI_COMM_NULL;

  // set MDI unit conversion factors

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

  qmIDs = nullptr;
  qm2owned = nullptr;

  eqm = nullptr;
  tqm = nullptr;
  xqm = nullptr;
  fqm = nullptr;
  qqm = nullptr;
  qpotential = nullptr;

  eqm_mine = nullptr;
  tqm_mine = nullptr;
  xqm_mine = nullptr;
  qqm_mine = nullptr;
  qpotential_mine = nullptr;

  // MM atom data

  nmm = nmm_last = max_nmm = 0;

  mmIDs = nullptr;
  mm2owned = nullptr;

  emm = nullptr;
  tmm = nullptr;
  xmm = nullptr;
  fmm = nullptr;
  qmm = nullptr;

  emm_mine = nullptr;
  tmm_mine = nullptr;
  xmm_mine = nullptr;
  qmm_mine = nullptr;

  // peratom Coulombic energy

  ecoul = nullptr;
  ncoulmax = 0;

  // per-atom data

  nmax = atom->nmax;
  memory->create(array_atom, nmax, 3, "mdi/qmmm:array_atom");

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

  delete[] elements;

  memory->destroy(qmIDs);
  memory->destroy(qm2owned);

  memory->destroy(eqm);
  memory->destroy(tqm);
  memory->destroy(xqm);
  memory->destroy(fqm);
  memory->destroy(qqm);
  memory->destroy(qpotential);

  memory->destroy(eqm_mine);
  memory->destroy(tqm_mine);
  memory->destroy(xqm_mine);
  memory->destroy(qqm_mine);
  memory->destroy(qpotential_mine);

  memory->destroy(mmIDs);
  memory->destroy(mm2owned);

  memory->destroy(emm);
  memory->destroy(tmm);
  memory->destroy(xmm);
  memory->destroy(fmm);
  memory->destroy(qmm);

  memory->destroy(emm_mine);
  memory->destroy(tmm_mine);
  memory->destroy(xmm_mine);
  memory->destroy(qmm_mine);

  memory->destroy(ecoul);
  memory->destroy(array_atom);
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

    // check if MDI engine supports commands that match QMMM mode

    if (mode == POTENTIAL) {
      int check1, check2;

      ierr = MDI_Check_command_exists("@DEFAULT", ">POTENTIAL_AT_NUCLEI", mdicomm, &check1);
      if (ierr) error->all(FLERR, "MDI: >POTENTIAL_AT_NUCLEI command check");
      MPI_Bcast(&check1, 1, MPI_INT, 0, world);

      ierr = MDI_Check_command_exists("@DEFAULT", "<CHARGES", mdicomm, &check2);
      if (ierr) error->all(FLERR, "MDI: <CHARGES command check");
      MPI_Bcast(&check2, 1, MPI_INT, 0, world);

      if (!check1 || !check2)
        error->all(FLERR, "Fix mdi/qmmm potential mode not supported by MDI engine");
    }

    if (mode == DIRECT) {
      int check1, check2, check3, check4;

      ierr = MDI_Check_command_exists("@DEFAULT", ">NLATTICE", mdicomm, &check1);
      if (ierr) error->all(FLERR, "MDI: >NLATTICE command check");
      MPI_Bcast(&check1, 1, MPI_INT, 0, world);

      ierr = MDI_Check_command_exists("@DEFAULT", ">CLATTICE", mdicomm, &check2);
      if (ierr) error->all(FLERR, "MDI: >CLATTICE command check");
      MPI_Bcast(&check2, 1, MPI_INT, 0, world);

      ierr = MDI_Check_command_exists("@DEFAULT", ">LATTICE", mdicomm, &check3);
      if (ierr) error->all(FLERR, "MDI: >LATTICE command check");
      MPI_Bcast(&check3, 1, MPI_INT, 0, world);

      ierr = MDI_Check_command_exists("@DEFAULT", "<LATTICE_FORCES", mdicomm, &check4);
      if (ierr) error->all(FLERR, "MDI: <LATTICE_FORCES command check");
      MPI_Bcast(&check4, 1, MPI_INT, 0, world);

      if (!check1 || !check2 || !check3 || !check4)
        error->all(FLERR, "Fix mdi/qmmm direct mode not supported by MDI engine");

      ierr = MDI_Check_command_exists("@DEFAULT", ">LATTICE_ELEMENTS", mdicomm, &check1);
      if (ierr) error->all(FLERR, "MDI: >LATTICE_ELEMENTS command check");
      MPI_Bcast(&check1, 1, MPI_INT, 0, world);

      ierr = MDI_Check_command_exists("@DEFAULT", ">LATTICE_TYPES", mdicomm, &check2);
      if (ierr) error->all(FLERR, "MDI: >LATTICE_TYPES command check");
      MPI_Bcast(&check2, 1, MPI_INT, 0, world);

      if (elements_exists && !check1)
        error->all(FLERR, "Fix mdi/qmmm direct mode elements not supported by MDI engine");
      if (types_exists && !check2)
        error->all(FLERR, "Fix mdi/qmmm direct mode types not supported by MDI engine");
    }
  }

  // require per-atom charge
  // POTENTIAL mode requires pair style to calculate only Coulombic interactions
  //   can be in conjunction with KSpace solver

  if (!atom->q_flag) error->all(FLERR, "Fix mdi/qmmm requires per-atom charge");

  if (mode == POTENTIAL) {
    if (!force->pair) error->all(FLERR, "Fix mdi/qmmm potential requires a pair style");
    pair_coul = force->pair_match("coul/cut", 1, 0);
    if (!pair_coul) pair_coul = force->pair_match("coul/long", 1, 0);
    if (!pair_coul) pair_coul = force->pair_match("coul/msm", 1, 0);
    if (!pair_coul)
      error->all(FLERR, "Fix mdi/qmmm potential requires Coulomb-only pair sub-style");
  }

  // determine whether a new vs incremental QMMM calc is needed
  // new if first run or if
  //   QM/MM atom counts or QM/MM elements/types or box has changed between runs
  // otherwise incremental = subsequent run of same system

  int new_system = 0;

  // check if count of QM or MM atoms has changed
  // on first run, old counts are 0

  int nqm_old = nqm;
  nqm = set_nqm();

  if (nqm != nqm_old) {
    if (nqm > max_nqm) reallocate_qm();
    create_qm_list();
    set_qm2owned();
    new_system = 1;
  }

  int nmm_old = nmm;
  nmm = set_nmm();

  if (nmm != nmm_old) {
    if (nmm > max_nmm) reallocate_mm();
    create_mm_list();
    set_mm2owned();
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

  // check if QM or MM atom elements or types have changed

  if (elements && elements_exists) {
    if (new_system) {
      set_eqm();
      if (mode == DIRECT) set_emm();
    } else {
      int *eqm_old;
      memory->create(eqm_old, nqm, "mdi/qmmm:eqm_old");
      memcpy(eqm_old, eqm, nqm * sizeof(int));
      set_eqm();
      for (int i = 0; i < nqm; i++)
        if (eqm[i] != eqm_old[i]) new_system = 1;
      memory->destroy(eqm_old);

      if (mode == DIRECT) {
        int *emm_old;
        memory->create(emm_old, nmm, "mdi/qmmm:emm_old");
        memcpy(emm_old, emm, nmm * sizeof(int));
        set_emm();
        for (int i = 0; i < nmm; i++)
          if (emm[i] != emm_old[i]) new_system = 1;
        memory->destroy(emm_old);
      }
    }

  } else if (types_exists) {
    if (new_system) {
      set_tqm();
      if (mode == DIRECT) set_tmm();
    } else {
      int *tqm_old;
      memory->create(tqm_old, nqm, "mdi/qmmm:tqm_old");
      memcpy(tqm_old, tqm, nqm * sizeof(int));
      set_tqm();
      for (int i = 0; i < nqm; i++)
        if (tqm[i] != tqm_old[i]) new_system = 1;
      memory->destroy(tqm_old);

      if (mode == DIRECT) {
        int *tmm_old;
        memory->create(tmm_old, nmm, "mdi/qmmm:tmm_old");
        memcpy(tmm_old, tmm, nmm * sizeof(int));
        set_tmm();
        for (int i = 0; i < nmm; i++)
          if (tmm[i] != tmm_old[i]) new_system = 1;
        memory->destroy(tmm_old);
      }
    }
  }

  // if new system, send setup info to MDI engine
  // values that often won't change for QMMM simulations
  // if not sending elements or types, assume engine initialized itself

  if (new_system) {
    send_natoms_qm();
    send_box();
    if (elements && elements_exists)
      send_elements_qm();
    else if (types_exists)
      send_types_qm();
    nqm_last = nqm;

    if (mode == DIRECT) {
      send_natoms_mm();
      if (elements && elements_exists)
        send_elements_mm();
      else if (types_exists)
        send_types_mm();
      set_qmm();
      send_charges_mm();
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
  int ilocal;
  double rsq;

  // invoke pair hybrid sub-style pair coulomb and Kspace directly
  // set eflag = 2 so they calculate per-atom energy

  pair_coul->compute(2, 0);
  double *eatom_pair = pair_coul->eatom;

  double *eatom_kspace = nullptr;
  if (force->kspace) {
    force->kspace->compute(2, 0);
    eatom_kspace = force->kspace->eatom;
  }

  // allocate ecoul for owned + ghost atoms
  // ghost atoms values only used if newton_pair is set

  if (atom->nmax > ncoulmax) {
    memory->destroy(ecoul);
    ncoulmax = atom->nmax;
    memory->create(ecoul, ncoulmax, "mdi/qmmm:ecoul");
  }

  // ecoul = per-atom energy for my owned atoms
  // if newton_pair, also do reverse_comm for ghost atom contribution

  int nlocal = atom->nlocal;
  int ntotal = nlocal;
  if (force->newton_pair) ntotal += atom->nghost;

  for (int i = 0; i < ntotal; i++) ecoul[i] = eatom_pair[i];

  if (force->newton_pair) comm->reverse_comm(this);

  if (force->kspace) {
    for (int i = 0; i < nlocal; i++) ecoul[i] += eatom_kspace[i];
  }

  // setup QM inputs: xqm and qpotential
  // xqm = atom coords, mapped into periodic box
  // qpotential[i] = Coulomb potential for each atom
  //   2 * (eatom[i] from pair_coul + kspace) / Qi
  //   factor of 2 comes from need to double count energy for each atom
  //   set for owned atoms, then MPI_Allreduce
  // subtract Qj/Rij energy for QM I interacting with all other QM J atoms
  //   use xqm_mine and qqm_mine for all QM atoms

  set_xqm(0);
  set_qqm();

  for (int i = 0; i < nqm; i++) qpotential_mine[i] = 0.0;

  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      if (q[ilocal] == 0.0)
        qpotential_mine[i] = 0.0;
      else
        qpotential_mine[i] = 2.0 * ecoul[ilocal] / q[ilocal];
    }
  }

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      for (int j = 0; j < nqm; j++) {
        if (j == i) continue;
        double delx = xqm[i][0] - xqm[j][0];
        double dely = xqm[i][1] - xqm[j][1];
        double delz = xqm[i][2] - xqm[j][2];
        domain->minimum_image(delx, dely, delz);
        rsq = delx * delx + dely * dely + delz * delz;
        qpotential_mine[i] -= qqrd2e * qqm[j] / sqrt(rsq);
      }
    }
  }

  MPI_Allreduce(qpotential_mine, qpotential, nqm, MPI_DOUBLE, MPI_SUM, world);

  // unit conversion from LAMMPS to MDI
  // must be done here, rather than in set_xqm()

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

  //if (comm->me == 0) utils::logmesg(lmp, "Invoking QM code ...\n");

  //MPI_Barrier(world);
  //double tstart = platform::walltime();

  int ierr;

  // if simulation box dynamically changes, send current box to MDI engine

  if (domain->box_change_size || domain->box_change_shape) {
    set_box();
    send_box();
  }

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

  // request QM energy from MDI engine
  // this triggers engine to perform QM calculation
  // sets qm_energy = fix output for global QM energy

  request_qm_energy();

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

  // request stress if needed and supported

  if (vflag && virialflag && stress_exists) {
    ierr = MDI_Send_command("<STRESS", mdicomm);
    if (ierr) error->all(FLERR, "MDI: <STRESS command");
    ierr = MDI_Recv(qm_virial, 9, MDI_DOUBLE, mdicomm);
    if (ierr) error->all(FLERR, "MDI: <STRESS data");
    MPI_Bcast(qm_virial, 9, MPI_DOUBLE, 0, world);
  }

  //MPI_Barrier(world);

  //if (comm->me == 0)
  //  utils::logmesg(lmp, "  time = {:.3f} seconds\n",
  //                 platform::walltime() - tstart);

  // reset owned charges to QM values
  // communicate changes to ghost atoms

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) q[ilocal] = qqm[i];
  }

  comm->forward_comm(this);

  // array_atom = fix output for peratom QM forces
  // zero array_atom first for MM atoms

  for (int i = 0; i < nlocal; i++) {
    array_atom[i][0] = 0.0;
    array_atom[i][1] = 0.0;
    array_atom[i][2] = 0.0;
  }

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

    // multiply by volume to make it extensive
    //   divide by nprocs so each proc stores a portion
    // this is b/c ComputePressure expects this as input from a fix
    //   it will do an MPI_Allreduce and divide by volume

    double volume;
    if (domain->dimension == 2)
      volume = domain->xprd * domain->yprd;
    else if (domain->dimension == 3)
      volume = domain->xprd * domain->yprd * domain->zprd;
    for (int i = 0; i < 6; i++) virial[i] = qm_virial_symmetric[i] * volume / nprocs;
  }

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
  if (mode == DIRECT)
    post_force_direct(vflag);
  else if (mode == POTENTIAL)
    post_force_potential(vflag);
}

/* ----------------------------------------------------------------------
   only called in DIRECT mode
   (1) send info on QM atoms to QM code
   (2) invoke the QM solver
   (3) receive results from QM code
---------------------------------------------------------------------- */

void FixMDIQMMM::post_force_direct(int vflag)
{
  // setup QM inputs: xqm = atom coords
  // setup MM inputs: xmm = atom coords

  set_xqm(1);
  set_xmm();

  // MDI communication with engine
  // first request for results triggers QM calculation
  // inputs: xqm, xmm
  // outputs: qm_energy, fqm, fmm

  //if (comm->me == 0) utils::logmesg(lmp, "Invoking QM code ...\n");

  //MPI_Barrier(world);
  //double tstart = platform::walltime();

  int ierr;

  // if simulation box dynamically changes, send current box to MDI engine

  if (domain->box_change_size || domain->box_change_shape) {
    set_box();
    send_box();
  }

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

  // request QM energy from MDI engine
  // this triggers engine to perform QM calculation
  // sets qm_energy = fix output for global QM energy

  request_qm_energy();

  // request forces on QM atoms from MDI engine

  ierr = MDI_Send_command("<FORCES", mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES command");
  ierr = MDI_Recv(&fqm[0][0], 3 * nqm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <FORCES data");
  MPI_Bcast(&fqm[0][0], 3 * nqm, MPI_DOUBLE, 0, world);

  // request forces on MM atoms from MDI engine

  ierr = MDI_Send_command("<LATTICE_FORCES", mdicomm);
  if (ierr) error->all(FLERR, "MDI: <LATTICE_FORCES command");
  ierr = MDI_Recv(&fmm[0][0], 3 * nmm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <LATTICE_FORCES data");
  MPI_Bcast(&fmm[0][0], 3 * nmm, MPI_DOUBLE, 0, world);

  // request stress if needed and supported

  if (vflag && virialflag && stress_exists) {
    ierr = MDI_Send_command("<STRESS", mdicomm);
    if (ierr) error->all(FLERR, "MDI: <STRESS command");
    ierr = MDI_Recv(qm_virial, 9, MDI_DOUBLE, mdicomm);
    if (ierr) error->all(FLERR, "MDI: <STRESS data");
    MPI_Bcast(qm_virial, 9, MPI_DOUBLE, 0, world);
  }

  //MPI_Barrier(world);
  //if (comm->me == 0)
  //  utils::logmesg(lmp, "  time = {:.3f} seconds\n",
  //                 platform::walltime() - tstart);

  // array_atom = fix output for peratom QM and MM forces

  int ilocal;
  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      array_atom[ilocal][0] = fqm[i][0] * mdi2lmp_force;
      array_atom[ilocal][1] = fqm[i][1] * mdi2lmp_force;
      array_atom[ilocal][2] = fqm[i][2] * mdi2lmp_force;
    }
  }

  for (int i = 0; i < nmm; i++) {
    ilocal = mm2owned[i];
    if (ilocal >= 0) {
      array_atom[ilocal][0] = fmm[i][0] * mdi2lmp_force;
      array_atom[ilocal][1] = fmm[i][1] * mdi2lmp_force;
      array_atom[ilocal][2] = fmm[i][2] * mdi2lmp_force;
    }
  }

  // add fqm and fmm to LAMMPS forces on respective atoms

  double **f = atom->f;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      f[ilocal][0] += fqm[i][0] * mdi2lmp_force;
      f[ilocal][1] += fqm[i][1] * mdi2lmp_force;
      f[ilocal][2] += fqm[i][2] * mdi2lmp_force;
    }
  }

  for (int i = 0; i < nmm; i++) {
    ilocal = mm2owned[i];
    if (ilocal >= 0) {
      f[ilocal][0] += fmm[i][0] * mdi2lmp_force;
      f[ilocal][1] += fmm[i][1] * mdi2lmp_force;
      f[ilocal][2] += fmm[i][2] * mdi2lmp_force;
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

    // set fix->virial
    // multiply by volume to make it extensive
    //   divide by nprocs so each proc stores a portion
    // this is b/c ComputePressure expects this as input from a fix
    //   it will do an MPI_Allreduce and divide by volume

    double volume;
    if (domain->dimension == 2)
      volume = domain->xprd * domain->yprd;
    else if (domain->dimension == 3)
      volume = domain->xprd * domain->yprd * domain->zprd;
    for (int i = 0; i < 6; i++) virial[i] = qm_virial_symmetric[i] * volume / nprocs;
  }
}

/* ----------------------------------------------------------------------
   only called in POTENTIAL mode
   add QM forces to QM atoms
   called after LAMMPS re-computes all MM forces with new QM charges
---------------------------------------------------------------------- */

void FixMDIQMMM::post_force_potential(int /*vflag*/)
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
      f[ilocal][0] += fqm[i][0] * mdi2lmp_force;
      f[ilocal][1] += fqm[i][1] * mdi2lmp_force;
      f[ilocal][2] += fqm[i][2] * mdi2lmp_force;
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

int FixMDIQMMM::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

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
  int i, m, last;

  m = 0;
  last = first + n;

  double *q = atom->q;

  for (i = first; i < last; i++) q[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int FixMDIQMMM::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = ecoul[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

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
  bytes += (double) ncoulmax * sizeof(double);
  bytes += (double) (3 * 3 + 4) * nqm * sizeof(double);    // fpoint QM arrays/vecs
  bytes += nqm * sizeof(tagint);                           // qmIDs
  bytes += nqm * sizeof(int);                              // qm2owned
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

void FixMDIQMMM::reallocate_qm()
{
  max_nqm = nqm;

  memory->destroy(qmIDs);
  memory->destroy(qm2owned);

  memory->destroy(eqm);
  memory->destroy(tqm);
  memory->destroy(xqm);
  memory->destroy(fqm);
  memory->destroy(qqm);
  memory->destroy(qpotential);

  memory->destroy(eqm_mine);
  memory->destroy(tqm_mine);
  memory->destroy(xqm_mine);
  memory->destroy(qqm_mine);
  memory->destroy(qpotential_mine);

  memory->create(qmIDs, max_nqm, "mdi/qmmm:qmIDs");
  memory->create(qm2owned, max_nqm, "mdi/qmmm:qm2owned");

  memory->create(eqm, max_nqm, "mdi/qmmm:eqm");
  memory->create(tqm, max_nqm, "mdi/qmmm:tqm");
  memory->create(xqm, max_nqm, 3, "mdi/qmmm:xqm");
  memory->create(fqm, max_nqm, 3, "mdi/qmmm:fqm");
  memory->create(qqm, max_nqm, "mdi/qmmm:qqm");
  memory->create(qpotential, max_nqm, "mdi/qmmm:qpotential");

  memory->create(eqm_mine, max_nqm, "mdi/qmmm:eqm_mine");
  memory->create(tqm_mine, max_nqm, "mdi/qmmm:tqm_mine");
  memory->create(xqm_mine, max_nqm, 3, "mdi/qmmm:xqm_mine");
  memory->create(qqm_mine, max_nqm, "mdi/qmmm:qqm_mine");
  memory->create(qpotential_mine, max_nqm, "mdi/qmmm:qpotential_mine");
}

/* ----------------------------------------------------------------------
   reallocate storage for MM atoms
------------------------------------------------------------------------- */

void FixMDIQMMM::reallocate_mm()
{
  max_nmm = nmm;

  memory->destroy(mmIDs);
  memory->destroy(mm2owned);

  memory->destroy(emm);
  memory->destroy(tmm);
  memory->destroy(xmm);
  memory->destroy(qmm);
  memory->destroy(fmm);

  memory->destroy(emm_mine);
  memory->destroy(tmm_mine);
  memory->destroy(xmm_mine);
  memory->destroy(qmm_mine);

  memory->create(mmIDs, max_nmm, "mdi/qmmm:mmIDs");
  memory->create(mm2owned, max_nmm, "mdi/qmmm:mm2owned");

  memory->create(emm, max_nmm, "mdi/qmmm:emm");
  memory->create(tmm, max_nmm, "mdi/qmmm:tmm");
  memory->create(xmm, max_nmm, 3, "mdi/qmmm:xmm");
  memory->create(qmm, max_nmm, "mdi/qmmm:qmm");
  memory->create(fmm, max_nmm, 3, "mdi/qmmm:fmm");

  memory->create(emm_mine, max_nmm, "mdi/qmmm:emm_mine");
  memory->create(tmm_mine, max_nmm, "mdi/qmmm:tmm_mine");
  memory->create(xmm_mine, max_nmm, 3, "mdi/qmmm:xmm_mine");
  memory->create(qmm_mine, max_nmm, "mdi/qmmm:qmm_mine");
}

/* ----------------------------------------------------------------------
   ncount = # of QM atoms = # in fix group
   return ncount to set nqm
------------------------------------------------------------------------- */

int FixMDIQMMM::set_nqm()
{
  bigint ngroup = group->count(igroup);

  // require 3*nqm be a small INT, so can MPI_Allreduce QM values

  if (3 * ngroup > MAXSMALLINT) error->all(FLERR, "Fix mdi/qmmm has too many quantum atoms");

  // error if nqm = 0
  // error if nqm = natoms, should use fix mdi/qm instead

  if (ngroup == 0) error->all(FLERR, "Fix mdi/qmmm has no atoms in quantum group");
  if (ngroup == atom->natoms) error->all(FLERR, "Fix mdi/qmmm has all atoms in quantum group");

  int ncount = ngroup;
  return ncount;
}

/* ----------------------------------------------------------------------
   ncount = # of MM atoms = all non-QM atoms
   return ncount to set nmm
------------------------------------------------------------------------- */

int FixMDIQMMM::set_nmm()
{
  // require 3*nmm be a small INT, so can MPI_Allreduce xmm

  if (3 * (atom->natoms - nqm) > MAXSMALLINT)
    error->all(FLERR, "Fix mdi/qmmm has too many classical atoms");

  int ncount = atom->natoms - nqm;
  return ncount;
}

/* ----------------------------------------------------------------------
   create sorted list of QM atom IDs
   ignore excluded atoms if exclude flag if set
------------------------------------------------------------------------- */

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
  memory->create(qmIDs_mine, nqm_mine, "mdi/qmmm:qmIDs_mine");

  nqm_mine = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) qmIDs_mine[nqm_mine++] = tag[i];

  int *recvcounts, *displs;
  memory->create(recvcounts, nprocs, "mdi/qmmm:recvcounts");
  memory->create(displs, nprocs, "mdi/qmmm:displs");

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

  memory->create(order, nqm, "mdi/qmmm:order");
  memory->create(qmIDs_sort, nqm, "mdi/qmmm:qmIDs_sort");

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
  memory->create(mmIDs_mine, nmm_mine, "mdi/qmmm:mmIDs_mine");

  nmm_mine = 0;
  for (int i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) mmIDs_mine[nmm_mine++] = tag[i];

  int *recvcounts, *displs;
  memory->create(recvcounts, nprocs, "mdi/qmmm:recvcounts");
  memory->create(displs, nprocs, "mdi/qmmm:displs");

  MPI_Allgather(&nmm_mine, 1, MPI_INT, recvcounts, 1, MPI_INT, world);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++)
    displs[iproc] = displs[iproc - 1] + recvcounts[iproc - 1];

  MPI_Allgatherv(mmIDs_mine, nmm_mine, MPI_LMP_TAGINT, mmIDs, recvcounts, displs, MPI_LMP_TAGINT,
                 world);

  memory->destroy(mmIDs_mine);
  memory->destroy(recvcounts);
  memory->destroy(displs);

  // sort mmIDs via merge sort

  int *order;
  tagint *mmIDs_sort;

  memory->create(order, nmm, "mdi/qmmm:order");
  memory->create(mmIDs_sort, nmm, "mdi/qmmm:mmIDs_sort");

  for (int i = 0; i < nmm; i++) {
    mmIDs_sort[i] = mmIDs[i];
    order[i] = i;
  }

  utils::merge_sort(order, nmm, (void *) mmIDs_sort, compare_IDs);

  int j;
  for (int i = 0; i < nmm; i++) {
    j = order[i];
    mmIDs_sort[i] = mmIDs[j];
  }

  memcpy(mmIDs, mmIDs_sort, nmm * sizeof(tagint));

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
    if (index >= nlocal)
      qm2owned[i] = -1;
    else
      qm2owned[i] = index;
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
    if (index >= nlocal)
      mm2owned[i] = -1;
    else
      mm2owned[i] = index;
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_box()
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

/* ----------------------------------------------------------------------
   fill xqm with QM atom coords
   if convert, perform LAMMPS to MDI unit conversion
   else no conver (used by POTENTIAL mode to subtract QM/QM Coulombics)
---------------------------------------------------------------------- */

void FixMDIQMMM::set_xqm(int convert)
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

      if (convert) {
        xqm_mine[i][0] *= lmp2mdi_length;
        xqm_mine[i][1] *= lmp2mdi_length;
        xqm_mine[i][2] *= lmp2mdi_length;
      }
    }
  }

  MPI_Allreduce(&xqm_mine[0][0], &xqm[0][0], 3 * nqm, MPI_DOUBLE, MPI_SUM, world);
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

  MPI_Allreduce(eqm_mine, eqm, nqm, MPI_INT, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_tqm()
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

  MPI_Allreduce(qqm_mine, qqm, nqm, MPI_DOUBLE, MPI_SUM, world);
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

      xmm_mine[i][0] *= lmp2mdi_length;
      xmm_mine[i][1] *= lmp2mdi_length;
      xmm_mine[i][2] *= lmp2mdi_length;
    }
  }

  MPI_Allreduce(&xmm_mine[0][0], &xmm[0][0], 3 * nmm, MPI_DOUBLE, MPI_SUM, world);
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

  MPI_Allreduce(emm_mine, emm, nmm, MPI_INT, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_tmm()
{
  for (int i = 0; i < nmm; i++) tmm_mine[i] = 0;

  int *type = atom->type;
  int ilocal;

  for (int i = 0; i < nmm; i++) {
    ilocal = mm2owned[i];
    if (ilocal >= 0) tmm_mine[i] = type[ilocal];
  }

  MPI_Allreduce(tmm_mine, tmm, nmm, MPI_INT, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void FixMDIQMMM::set_qmm()
{
  for (int i = 0; i < nmm; i++) qmm_mine[i] = 0.0;

  double *q = atom->q;
  int ilocal;

  for (int i = 0; i < nmm; i++) {
    ilocal = mm2owned[i];
    if (ilocal >= 0) qmm_mine[i] = q[ilocal];
  }

  MPI_Allreduce(qmm_mine, qmm, nmm, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   send LAMMPS QM atom count to MDI engine
------------------------------------------------------------------------- */

void FixMDIQMMM::send_natoms_qm()
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
   send QM atom types to MDI engine
------------------------------------------------------------------------- */

void FixMDIQMMM::send_types_qm()
{
  int ierr = MDI_Send_command(">TYPES", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >TYPES command");
  ierr = MDI_Send(tqm, nqm, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >TYPES data");
}

/* ----------------------------------------------------------------------
   send QM elements to MDI engine = atomic numbers for each type
------------------------------------------------------------------------- */

void FixMDIQMMM::send_elements_qm()
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

void FixMDIQMMM::send_box()
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

  } else if (domain->xperiodic == 1 || domain->yperiodic == 1 || domain->zperiodic == 1) {
    error->all(FLERR, "MDI requires fully periodic or fully non-periodic system");
  }
}

/* ----------------------------------------------------------------------
   send LAMMPS MM atom count to MDI engine
------------------------------------------------------------------------- */

void FixMDIQMMM::send_natoms_mm()
{
  int ierr = MDI_Send_command(">NLATTICE", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >NLATTICE command");
  ierr = MDI_Send(&nmm, 1, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >NLATTICE data");
}

/* ----------------------------------------------------------------------
   send MM atom types to MDI engine
------------------------------------------------------------------------- */

void FixMDIQMMM::send_types_mm()
{
  int ierr = MDI_Send_command(">LATTICE_TYPES", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >LATTICE_TYPES command");
  ierr = MDI_Send(tqm, nqm, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >LATTICE_TYPES data");
}

/* ----------------------------------------------------------------------
   send MM elements to MDI engine = atomic numbers for each type
------------------------------------------------------------------------- */

void FixMDIQMMM::send_elements_mm()
{
  int ierr = MDI_Send_command(">LATTICE_ELEMENTS", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >LATTICE_ELEMENTS command");
  ierr = MDI_Send(emm, nmm, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >LATTICE_ELEMENTS data");
}

/* ----------------------------------------------------------------------
   send MM charges to MDI engine
------------------------------------------------------------------------- */

void FixMDIQMMM::send_charges_mm()
{
  int ierr = MDI_Send_command(">LATTICE", mdicomm);
  if (ierr) error->all(FLERR, "MDI: >LATTICE command");
  ierr = MDI_Send(qmm, nmm, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >LATTICE data");
}

/* ----------------------------------------------------------------------
   request QM energy from MDI engine
   set qm_energy = fix output for global QM energy
------------------------------------------------------------------------- */

void FixMDIQMMM::request_qm_energy()
{
  int ierr;

  // QM energy = <PE + <KE_ELEC or <ENERGY, depending on engine options

  if (pe_exists && keelec_exists) {
    double pe_energy, keelec_energy;

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
