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
   Contributing authors: Charlles R. A. Abreu (abreu@eq.ufrj.br)
                         Ana J. Silveira (asilveira@plapiqui.edu.ar)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "math_extra.h"
#include "fix_hmc.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "random_park.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "compute.h"
#include "output.h"
#include "neighbor.h"
#include "fix_rigid_nve_small.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 1.0e-7

enum{ ATOMS, VCM_OMEGA, XCM, ITENSOR, ROTATION, FORCE_TORQUE };

/* ---------------------------------------------------------------------- */

FixHMC::FixHMC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),random_equal(NULL)
{
  if (narg < 7) error->all(FLERR, "Illegal fix hmc command");

  // Retrieve user-defined options:
  nevery = utils::numeric(FLERR, arg[3], false, lmp);    // Number of MD steps per MC step
  int seed = utils::numeric(FLERR, arg[4], false, lmp);    // Seed for random number generation
  double temp = utils::numeric(FLERR, arg[5], false, lmp);    // System temperature

  // Retrieve the molecular dynamics integrator type:
  mdi = arg[6];
  if ( strcmp(mdi,"rigid") != 0 && strcmp(mdi,"flexible") != 0 )
    error->all(FLERR, "Illegal fix hmc command");

  KT = force->boltz * temp / force->mvv2e;      // K*T in mvv units
  mbeta = -1.0/(force->boltz * temp);           // -1/(K*T) in energy units

  // Check keywords:
  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"adjust") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix hmc command");
      if (strcmp(arg[iarg+1],"yes") == 0) tune_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) tune_flag = 0;
      else error->all(FLERR,"Illegal fix hmc command");
      iarg += 2;
    }
    else error->all(FLERR,"Illegal fix hmc command");
  }

  // Initialize RNG with a different seed for each process:
  random = new RanPark(lmp,seed + comm->me);
  for (int i = 0; i < 100; i++) random->gaussian();
  random_equal = new RanPark(lmp,seed);
  // Perform initialization of per-atom arrays:
  xu = NULL;
  deltax = NULL;
  scal = NULL;
  vec = NULL;
  eatom = NULL;
  vatom = NULL;

  // Register callback:
  atom->add_callback(0);

  // Initialize arrays and pointers for saving/restoration of states:
  setup_arrays_and_pointers();

  // Add new computes for global and per-atom properties:
  add_new_computes();

  // Define non-default fix attributes:
  global_freq = 1;
  scalar_flag = 1;
  extscalar = 0;
  vector_flag = 1;
  extvector = 0;
  size_vector = 5;
}

/* ---------------------------------------------------------------------- */

FixHMC::~FixHMC()
{
  atom->delete_callback(id,0);
  memory->destroy(xu);
  memory->destroy(deltax);
  memory->destroy(scal);
  memory->destroy(vec);
  memory->destroy(eglobal);
  memory->destroy(vglobal);
  memory->destroy(eatom);
  memory->destroy(vatom);
  delete [] scalptr;
  delete [] vecptr;
  delete [] eglobalptr;
  delete [] vglobalptr;
  delete [] eatomptr;
  delete [] vatomptr;
  delete [] rev_comm;
  delete random_equal;
  modify->delete_compute("hmc_ke");
  modify->delete_compute("hmc_pe");
  modify->delete_compute("hmc_peatom");
  modify->delete_compute("hmc_press");
  modify->delete_compute("hmc_pressatom");

  for (LAMMPS_NS::Atom::PerAtom &stored_peratom_member : stored_peratom) {
    free(stored_peratom_member.address);
  }
}

/* ---------------------------------------------------------------------- */

void FixHMC::post_constructor()
{
  char **newarg = new char*[4];
  newarg[0] = (char *) "hmc_mdi";
  newarg[1] = group->names[igroup];
  if (strcmp(mdi,"flexible") == 0) {
    newarg[2] = (char *) "nve";
    modify->add_fix(3,newarg);
  }
  else {
    newarg[2] = (char *) "rigid/nve/small";
    newarg[3] = (char *) "molecule";
    modify->add_fix(4,newarg);
  }
  class Fix* mdfix = modify->fix[modify->find_fix("hmc_mdi")];
  rigid_flag = mdfix->rigid_flag;
  if (rigid_flag)
    fix_rigid = (class FixRigidSmall*) mdfix;
}

/* ---------------------------------------------------------------------- */

template <typename T>
void store_peratom_member(LAMMPS_NS::Atom::PerAtom &stored_peratom_member,
                          LAMMPS_NS::Atom::PerAtom current_peratom_member, int nlocal)
{
  if (stored_peratom_member.name.compare(current_peratom_member.name)) {
    printf(
        "NONONONONONO\n");    //error->all(FLERR, "fix hmc tried to store incorrect peratom data");
  }
  size_t offset;
  int cols;
  // free old memory if stored_peratom_member isn't a copy of current_peratom_member
  if (stored_peratom_member.address != current_peratom_member.address) {
    free(stored_peratom_member.address);
    stored_peratom_member.address = NULL;
  }
  if (stored_peratom_member.address_maxcols != current_peratom_member.address_maxcols) {
    free(stored_peratom_member.address_maxcols);
    stored_peratom_member.address_maxcols = NULL;
  }
  // peratom scalers
  if (current_peratom_member.cols == 0) {
    stored_peratom_member.address = malloc(sizeof(T) * nlocal);
    memcpy(stored_peratom_member.address, current_peratom_member.address, nlocal * sizeof(T));
  } else {
    // peratom vectors
    if (current_peratom_member.cols < 0) {
      // variable column case
      cols = *(current_peratom_member.address_maxcols);
      stored_peratom_member.address_maxcols = (int *) malloc(sizeof(int));
      *(stored_peratom_member.address_maxcols) = *(current_peratom_member.address_maxcols);
    } else {
      // non-variable column case
      cols = current_peratom_member.cols;
    }
    stored_peratom_member.address = malloc(sizeof(T) * nlocal * cols);
    for (int i = 0; i < nlocal; i++) {
      memcpy((T *) stored_peratom_member.address + i * cols,
             (T *) current_peratom_member.address + i * cols, sizeof(T) * cols);
    }
  }
  stored_peratom_member.cols = current_peratom_member.cols;
  stored_peratom_member.collength = current_peratom_member.collength;
}


template <typename T>
void restore_peratom_member(LAMMPS_NS::Atom::PerAtom stored_peratom_member,
                          LAMMPS_NS::Atom::PerAtom &current_peratom_member, int nlocal)
{
  if (stored_peratom_member.name.compare(current_peratom_member.name)) {
    printf("NONONONONONO\n");    //error->all(FLERR, "fix hmc tried to store incorrect peratom data");
  }
  size_t offset;
  int cols;
  // peratom scalers
  if (stored_peratom_member.cols == 0) {
    memcpy(current_peratom_member.address, stored_peratom_member.address, nlocal * sizeof(T));
  } else {
    // peratom vectors
    if (stored_peratom_member.cols < 0) {
      // variable column case
      cols = *(stored_peratom_member.address_maxcols);
      *(current_peratom_member.address_maxcols) = *(stored_peratom_member.address_maxcols);
    } else {
      // non-variable column case
      cols = stored_peratom_member.cols;
    }
    for (int i = 0; i < nlocal; i++) {
      offset = i * cols * sizeof(T);
      memcpy((T *) current_peratom_member.address + i * cols,
             (T *) stored_peratom_member.address + i * cols,
             sizeof(T) * cols);
    }
  }
  current_peratom_member.cols = stored_peratom_member.cols;
  current_peratom_member.collength = stored_peratom_member.collength;
}


void FixHMC::setup_arrays_and_pointers()
{
  int i, j, m;
  int pair_flag;
  int bond_flag;
  int angle_flag;
  int dihedral_flag;
  int improper_flag;
  int kspace_flag;

  // Per-atom scalar properties to be saved and restored:
  nscal = 0;
  if (atom->erforce_flag) nscal++;
  //if (atom->e_flag) nscal++;
  if (atom->rho_flag) nscal++;
  scalptr = new double**[nscal];
  m = 0;
  if (atom->erforce_flag) scalptr[m++] = &atom->erforce;
  //if (atom->e_flag) scalptr[m++] = &atom->de;
  if (atom->rho_flag) scalptr[m++] = &atom->drho;
  
  current_peratom = atom->peratom;
  stored_nlocal = atom->nlocal;

    // make a copy of all peratom data
  for (LAMMPS_NS::Atom::PerAtom &current_peratom_member : current_peratom) {
    LAMMPS_NS::Atom::PerAtom stored_peratom_member = current_peratom_member;
    switch (current_peratom_member.datatype) {
      case (Atom::INT):
        store_peratom_member<int>(stored_peratom_member, current_peratom_member, stored_nlocal);
        break;
      case (Atom::DOUBLE):
        store_peratom_member<double>(stored_peratom_member, current_peratom_member, stored_nlocal);
        break;
      case (Atom::BIGINT):
        store_peratom_member<bigint>(stored_peratom_member, current_peratom_member, stored_nlocal);
        break;
    }
    stored_peratom.push_back(stored_peratom_member);
  }

  // Per-atom vector properties to be saved and restored:
  nvec = 2;
  if (atom->torque_flag) nvec++;
  vecptr = new double***[nvec];
  m = 0;
  vecptr[m++] = &xu;
  vecptr[m++] = &atom->f;
  if (atom->torque_flag) vecptr[m++] = &atom->torque;

  // Determine which energy contributions must be computed:
  ne = 0;
  if (force->pair) { pair_flag = 1; ne++; } else pair_flag = 0;
  if (force->bond) { bond_flag = 1; ne++; } else bond_flag = 0;
  if (force->angle) { angle_flag = 1; ne++; } else angle_flag = 0;
  if (force->dihedral) { dihedral_flag = 1; ne++; } else dihedral_flag = 0;
  if (force->improper) { improper_flag = 1; ne++; } else improper_flag = 0;
  if (force->kspace) { kspace_flag = 1; ne++; } else kspace_flag = 0;

  // Initialize arrays for managing global energy terms:
  neg = pair_flag ? ne + 1 : ne;
  memory->create(eglobal,neg,"fix_hmc:eglobal");
  eglobalptr = new double*[neg];
  m = 0;
  if (pair_flag) {
    eglobalptr[m++] = &force->pair->eng_vdwl;
    eglobalptr[m++] = &force->pair->eng_coul;
  }
  if (bond_flag) eglobalptr[m++] = &force->bond->energy;
  if (angle_flag) eglobalptr[m++] = &force->angle->energy;
  if (dihedral_flag) eglobalptr[m++] = &force->dihedral->energy;
  if (improper_flag) eglobalptr[m++] = &force->improper->energy;
  if (kspace_flag) eglobalptr[m++] = &force->kspace->energy;

  // Initialize arrays for managing global virial terms:
  nv = ne;
  for (j = 0; j < modify->nfix; j++)
    if (modify->fix[j]->virial_global_flag) nv++;
  memory->create(vglobal,nv,6,"fix_hmc:vglobal");
  vglobalptr = new double**[nv];
  for (m = 0; m < nv; m++) vglobalptr[m] = new double*[6];
  for (i = 0; i < 6; i++) {
    m = 0;
    if (pair_flag) vglobalptr[m++][i] = &force->pair->virial[i];
    if (bond_flag) vglobalptr[m++][i] = &force->bond->virial[i];
    if (angle_flag) vglobalptr[m++][i] = &force->angle->virial[i];
    if (dihedral_flag) vglobalptr[m++][i] = &force->dihedral->virial[i];
    if (improper_flag) vglobalptr[m++][i] = &force->improper->virial[i];
    if (kspace_flag) vglobalptr[m++][i] = &force->kspace->virial[i];
    for (j = 0; j < modify->nfix; j++)
      if (modify->fix[j]->virial_global_flag)
        vglobalptr[m++][i] = &modify->fix[j]->virial[i];
  }

  // Determine which per-atom energy terms require reverse communication:
  rev_comm = new int[nv];
  m = 0;
  if (pair_flag) rev_comm[m++] = force->newton;
  if (bond_flag) rev_comm[m++] = force->newton_bond;
  if (angle_flag) rev_comm[m++] = force->newton_bond;
  if (dihedral_flag) rev_comm[m++] = force->newton_bond;
  if (improper_flag) rev_comm[m++] = force->newton_bond;
  if (kspace_flag) rev_comm[m++] = force->kspace->tip4pflag;
  for (i = ne; i < nv; i++) rev_comm[m++] = 0;

  // Initialize array of pointers to manage per-atom energies:
  eatomptr = new double**[ne];
  m = 0;
  if (pair_flag) eatomptr[m++] = &force->pair->eatom;
  if (bond_flag) eatomptr[m++] = &force->bond->eatom;
  if (angle_flag) eatomptr[m++] = &force->angle->eatom;
  if (dihedral_flag) eatomptr[m++] = &force->dihedral->eatom;
  if (improper_flag) eatomptr[m++] = &force->improper->eatom;
  if (kspace_flag) eatomptr[m++] = &force->kspace->eatom;

  // Initialize array of pointers to manage per-atom virials:
  vatomptr = new double***[nv];
  m = 0;
  if (pair_flag) vatomptr[m++] = &force->pair->vatom;
  if (bond_flag) vatomptr[m++] = &force->bond->vatom;
  if (angle_flag) vatomptr[m++] = &force->angle->vatom;
  if (dihedral_flag) vatomptr[m++] = &force->dihedral->vatom;
  if (improper_flag) vatomptr[m++] = &force->improper->vatom;
  if (kspace_flag) vatomptr[m++] = &force->kspace->vatom;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->virial_peratom_flag)
      vatomptr[m++] = &modify->fix[i]->vatom;

  // Determine the maximum and the actual numbers of per-atom variables in reverse
  // communications:
  // Note: fix-related virials do not communicate (thus 'ne' used instead of 'nv')
  comm_reverse = 0;
  ncommrev = 0;
  for (m = 0; m < ne; m++)
    if (rev_comm[m]) {
      comm_reverse += 7;  // 1 energy + 6 virials
      if (peatom_flag) ncommrev++;
      if (pressatom_flag) ncommrev += 6;
    }

  // Determine maximum number of per-atom variables in forward and reverse
  // communications when dealing with rigid bodies:
  if (rigid_flag) {
    comm_reverse = MAX(comm_reverse,6);
    comm_forward = 12;
  }
}

/* ---------------------------------------------------------------------- */

void FixHMC::add_new_computes()
{
  char **newarg = new char*[5];

  // Add all new computes for group "all":
  newarg[1] = (char *) "all";

  // Kinetic energy:
  newarg[0] = (char *) "hmc_ke";
  newarg[2] = (char *) "ke";
  modify->add_compute(3,newarg);
  ke = modify->compute[modify->ncompute-1];
  //ke = modify->compute[modify->find_compute("thermo_ke")];

  // Potential energy:
  newarg[0] = (char *) "hmc_pe";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg);
  pe = modify->compute[modify->ncompute-1];
  //pe = modify->compute[modify->find_compute("thermo_pe")];

  // Per-atom potential energy:
  newarg[0] = (char *) "hmc_peatom";
  newarg[2] = (char *) "pe/atom";
  modify->add_compute(3,newarg);
  peatom = modify->compute[modify->ncompute-1];

  // System pressure:
  newarg[0] = (char *) "hmc_press";
  newarg[2] = (char *) "pressure";
  newarg[3] = (char *) "NULL";
  newarg[4] = (char *) "virial";
  modify->add_compute(5,newarg);
  press = modify->compute[modify->ncompute-1];

  // Per-atom stress tensor:
  newarg[0] = (char *) "hmc_pressatom";
  newarg[2] = (char *) "stress/atom";
  modify->add_compute(5,newarg);
  pressatom = modify->compute[modify->ncompute-1];

  delete [] newarg;
}

/* ---------------------------------------------------------------------- */

int FixHMC::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHMC::tune_parameter(int *parameter, const char *name)
{
  if (*parameter % nevery != 0) {
    if (tune_flag) {
      *parameter = (*parameter / nevery + 1)*nevery;
      if (comm->me == 0) printf("Fix HMC: adjusting %s to %d\n",name,*parameter);
    }
    else {
      if (comm->me == 0) printf("Fix HMC: %s is not a multiple of nevery\n",name);
      error->all(FLERR,"illegal parameter value");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixHMC::init()
{
  int ntimestep = update->ntimestep;

  // Check conformance of run length:
  tune_parameter( &update->nsteps, "number of steps of the run" );
  update->laststep = ntimestep + update->nsteps;

  // Check conformance of thermo interval:
  if (output->var_thermo)
    error->all(FLERR,"fix HMC does not allow a variable thermo interval");
  else if (output->thermo_every)
    tune_parameter( &output->thermo_every, "thermo interval" );

  // Check conformance of restart interval:
  if (output->restart_flag) {
    if (output->restart_flag_single)
      tune_parameter( &output->restart_every_single, "restart interval" );
    if (output->restart_flag_double)
      tune_parameter( &output->restart_every_double, "restart interval" );
  }

  // Check conformance of dump interval:
  for (int i = 0; i < output->ndump; i++)
    tune_parameter( &output->every_dump[i], "dump interval" );

  // Check whether there is any fixes that change box size and/or shape:
  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i]->box_change)
      error->all(FLERR,"fix hmc is incompatible with fixes that change box size or shape");
  }

  // Check whether there are subsequent fixes with active virial_flag:
  int first = modify->find_fix(this->id) + 1;
  if (rigid_flag) first++;
  for (int i = first; i < modify->nfix; i++)
    if (modify->fix[i]->virial_peratom_flag || modify->fix[i]->virial_global_flag) {
      if (comm->me == 0)
        printf("Fix %s defined after fix hmc.\n",modify->fix[i]->style);
      error->all(FLERR,"fix hmc cannot precede fixes that modify the system pressure");
    }

  // Look for computes with active peatomflag, press_flag, or pressatomflag:
  peatom_flag = 0;
  press_flag = 0;
  pressatom_flag = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strncmp(modify->compute[i]->id,"hmc_",4) != 0) {
      peatom_flag = peatom_flag | modify->compute[i]->peatomflag;
      press_flag = press_flag | modify->compute[i]->pressflag;
      pressatom_flag = pressatom_flag | modify->compute[i]->pressatomflag;
    }

  // Count per-atom properties to be exchanged:
  nvalues = nscal + 3*nvec;
  if (peatom_flag) nvalues += ne;
  if (pressatom_flag) nvalues += 6*nv;

  // (Re)allocate array of per-atom properties:
  grow_arrays(atom->nmax);

  // Activate potential energy and other necessary calculations at setup:
  //pe->addstep(ntimestep);
  if (peatom_flag) peatom->addstep(ntimestep);
  if (press_flag) press->addstep(ntimestep);
  if (pressatom_flag) pressatom->addstep(ntimestep);
}

/* ----------------------------------------------------------------------
   Initialize MC move, save current state, and activate computes
------------------------------------------------------------------------- */

void FixHMC::setup(int vflag)
{
  if (rigid_flag)
      error->all(FLERR, "fix hmc does not yet support rigid");
      //if (fix_rigid->extended)
      //error->all(FLERR,"fix hmc does not support extended particles");

  // Compute properties of the initial state:
  nattempts = 0;
  naccepts = 0;
  DeltaPE = 0;
  DeltaKE = 0;
  if (rigid_flag) {
    //rigid_body_atom_positions(xu);
    atom_positions(xu);
    rigid_body_random_velocities();
  }
  else {
    atom_positions(xu);
    random_velocities();
  }
  for (int i = 0; i < atom->nlocal; i++)
    deltax[i][0] = deltax[i][1] = deltax[i][2] = 0.0;
  update->eflag_global = update->ntimestep;
  PE = pe->compute_scalar();
  KE = ke->compute_scalar();
  save_current_state();

  // Activate potential energy and other necessary calculations:
  int nextstep = update->ntimestep + nevery;
  pe->addstep(nextstep);
  if (peatom_flag) peatom->addstep(nextstep);
  if (press_flag) press->addstep(nextstep);
  if (pressatom_flag) pressatom->addstep(nextstep);
}

/* ----------------------------------------------------------------------
   Apply the Metropolis acceptance criterion
   Restore saved system state if move is rejected
   Activate computes for the next MC step
------------------------------------------------------------------------- */

void FixHMC::end_of_step()
{
  nattempts++;

  // Compute proposed unwrapped positions and displacements:
  if (rigid_flag)
    atom_positions(xu); 
  //  rigid_body_atom_positions(xu);
  else
    atom_positions(xu);
  int nlocal = atom->nlocal;
  double **old_xu = vec[0];
  for (int i = 0; i < nlocal; i++) {
    deltax[i][0] = xu[i][0] - old_xu[i][0];
    deltax[i][1] = xu[i][1] - old_xu[i][1];
    deltax[i][2] = xu[i][2] - old_xu[i][2];
  }

  // Compute potential and kinetic energy variations:
  update->eflag_global = update->ntimestep;
  double newPE = pe->compute_scalar();
  double newKE = ke->compute_scalar();
  DeltaPE = newPE - PE;
  DeltaKE = newKE - KE;

  // Apply the Metropolis criterion:
  double DeltaE = DeltaPE + DeltaKE;
  int accept = DeltaE < 0.0;
  if (~accept) {
    accept = random_equal->uniform() <= exp(mbeta*DeltaE);
    MPI_Bcast(&accept,1,MPI_INT,0,world);
  }
  if (accept) {
    // Update potential energy and save the current state:
    naccepts++;
    PE = newPE;
    save_current_state();
  }
  else {
    // Restore saved state and enforce check_distance/reneighboring in the next step:
    restore_saved_state();
    neighbor->ago = (neighbor->delay/neighbor->every + 1)*neighbor->every;
  }

  // Choose new velocities and compute kinetic energy:
  if (rigid_flag)
    rigid_body_random_velocities();
  else
    random_velocities();
  KE = ke->compute_scalar();

  // Activate potential energy and other necessary calculations:
  int nextstep = update->ntimestep + nevery;
  if (nextstep <= update->laststep) {
    pe->addstep(nextstep);
    if (peatom_flag) peatom->addstep(nextstep);
    if (press_flag) press->addstep(nextstep);
    if (pressatom_flag) pressatom->addstep(nextstep);
  }
}

/* ----------------------------------------------------------------------
   Return the acceptance fraction of proposed MC moves
------------------------------------------------------------------------- */

double FixHMC::compute_scalar()
{
  double acc_frac = naccepts;
  acc_frac /= MAX(1,nattempts);
  return acc_frac;
}

/* ----------------------------------------------------------------------
   Return the acceptance fraction of proposed MC moves, or
   the total energy variation of the last proposed MC move, or
   the mean-square atom displacement in the last proposed MC move
------------------------------------------------------------------------- */

double FixHMC::compute_vector(int item)
{
  int n = item + 1;
  if (n == 1) {
    double acc_frac = naccepts;
    acc_frac /= MAX(1,nattempts);
    return acc_frac;
  }
  else if (n == 2)
    return DeltaPE;
  else if (n == 3)
    return DeltaKE;
  else if (n == 4)
    return DeltaPE + DeltaKE;
  else if (n == 5) {
    double *dx, local_msd = 0.0;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) {
      dx = deltax[i];
      local_msd += dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    }
    double msd;
    MPI_Allreduce(&local_msd,&msd,1,MPI_DOUBLE,MPI_SUM,world);
    msd /= atom->natoms;
    return msd;
  }
  else
    return 0.0;
}

/* ----------------------------------------------------------------------
   Save the system state for eventual restoration if move is rejected
------------------------------------------------------------------------- */

void FixHMC::save_current_state()
{
  int i, m, n;
  int nlocal = atom->nlocal;
  int ntotal = nlocal + atom->nghost;
  double *scalar, **vector, *energy, **stress;

  // Save per-atom scalar properties:
  for (m = 0; m < nscal; m++) {
    scalar = *scalptr[m];
    for (i = 0; i < nlocal; i++)
      scal[m][i] = scalar[i];
  }

  // Save per-atom vector properties:
  for (m = 0; m < nvec; m++) {
    vector = *vecptr[m];
    for (i = 0; i < nlocal; i++)
      memcpy( vec[m][i], vector[i], three );
  }

  // Save global energy terms:
  for (m = 0; m < neg; m++)
    eglobal[m] = *eglobalptr[m];

  // Save global virial terms:
  if (press_flag)
    for (m = 0; m < nv; m++)
      memcpy( vglobal[m], *vglobalptr[m], six );

  // Save per-atom energy terms for all local atoms,
  // and also for ghost atoms when needed:
  if (peatom_flag)
    for (m = 0; m < ne; m++) {
      n = rev_comm[m] ? ntotal : nlocal;
      energy = *eatomptr[m];
      for (i = 0; i < n; i++)
        eatom[m][i] = energy[i];
    }

  // Save per-atom virial terms for all local atoms,
  // and also for ghost atoms when needed:
  if (pressatom_flag)
    for (m = 0; m < nv; m++) {
      n = rev_comm[m] ? ntotal : nlocal;
      stress = *vatomptr[m];
      for (i = 0; i < n; i++)
        memcpy( vatom[m][i], stress[i], six );
    }

  // Perform reverse communication to incorporate ghost atoms info:
  if (comm_reverse && (peatom_flag || pressatom_flag)) {
    comm_flag = ATOMS;
    comm->reverse_comm(this, ncommrev);
  }
}

/* ----------------------------------------------------------------------
   Restore the system state saved at the beginning of the MC step
------------------------------------------------------------------------- */

void FixHMC::restore_saved_state()
{
  int i, m;
  int nlocal = atom->nlocal;
  int ntotal = nlocal + atom->nghost;
  double **x = atom->x;
  double *scalar, **vector, *energy, **stress;

  if (nlocal != stored_nlocal) error->all(FLERR, "bonk");
  for (LAMMPS_NS::Atom::PerAtom &stored_peratom_member : stored_peratom) {
    for (LAMMPS_NS::Atom::PerAtom &current_peratom_member : current_peratom) {
      if (stored_peratom_member.name.compare(current_peratom_member.name)) {
        continue;
      } else {
        switch (current_peratom_member.datatype) {
          case (Atom::INT):
            restore_peratom_member<int>(stored_peratom_member, current_peratom_member, stored_nlocal);
            break;
          case (Atom::DOUBLE):
            restore_peratom_member<double>(stored_peratom_member, current_peratom_member,
                                         stored_nlocal);
            break;
          case (Atom::BIGINT):
            restore_peratom_member<bigint>(stored_peratom_member, current_peratom_member,
                                         stored_nlocal);
            break;
        }
        break;
      }
    }
  }
  //// Restore scalar properties:
  //for (m = 0; m < nscal; m++) {
  //  scalar = *scalptr[m];
  //  for (i = 0; i < nlocal; i++)
  //    scalar[i] = scal[m][i];
  //}

  //// Restore vector properties:
  //for (m = 0; m < nvec; m++) {
  //  vector = *vecptr[m];
  //  for (i = 0; i < nlocal; i++)
  //    memcpy( vector[i], vec[m][i], three );
  //}

  //// Relocate atoms:
  //for (i = 0; i < nlocal; i++) {
  //  x[i][0] -= deltax[i][0];
  //  x[i][1] -= deltax[i][1];
  //  x[i][2] -= deltax[i][2];
  //}

  //// Finish with relocation of rigid bodies:
  //if (rigid_flag) {
  //  rigid_body_restore_positions(deltax);
  //  rigid_body_restore_orientations();
  //  rigid_body_restore_forces();
  //}

  // Restore global energy terms:
  for (i = 0; i < neg; i++)
    *eglobalptr[i] = eglobal[i];

  // Restore global virial terms:
  if (press_flag)
    for (i = 0; i < nv; i++)
      memcpy( *vglobalptr[i], vglobal[i], six );

  //// Restore per-atom energy terms for all local atoms,
  //// and zero those for ghost atoms when needed:
  //if (peatom_flag)
  //  for (m = 0; m < ne; m++) {
  //    energy = *eatomptr[m];
  //    for (i = 0; i < nlocal; i++)
  //      energy[i] = eatom[m][i];
  //    if (rev_comm[m])
  //      for (i = nlocal; i < ntotal; i++)
  //        energy[i] = 0.0;
  //  }

  //// Restore per-atom virial terms for all local atoms,
  //// and zero those for ghost atoms when needed:
  //if (pressatom_flag)
  //  for (m = 0; m < nv; m++) {
  //    stress = *vatomptr[m];
  //    for (i = 0; i < nlocal; i++)
  //      memcpy( stress[i], vatom[m][i], six );
  //    if (rev_comm[m])
  //      for (i = nlocal; i < ntotal; i++)
  //        memset( stress[i], 0, six );
  //  }
}


/* ----------------------------------------------------------------------
  Compute the unwrapped position of an atom present in a rigid body
------------------------------------------------------------------------- */

void FixHMC::atom_positions(double **xu)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  tagint *image = atom->image;
  for (int i = 0; i < nlocal; i++)
      domain->unmap(x[i],image[i],xu[i]);
}

/* ----------------------------------------------------------------------
   Randomly choose velocities from a Maxwell-Boltzmann distribution
------------------------------------------------------------------------- */

void FixHMC::random_velocities()
{
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;

  double stdev;
  int nlocal;

  double *rmass = atom->rmass;
  double *mass = atom->mass;
  nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) stdev = sqrt(KT/rmass[i]);
      else stdev = sqrt(KT/mass[type[i]]);
      v[i][0] = stdev*random->gaussian();
      v[i][1] = stdev*random->gaussian();
      v[i][2] = stdev*random->gaussian();
    }
}

/* ----------------------------------------------------------------------
  Compute the unwrapped positions of atoms in a system with rigid bodies
------------------------------------------------------------------------- */

//void FixHMC::rigid_body_atom_positions(double **xu)
//{
//  int nlocal = atom->nlocal;
//  double **x = atom->x;
//  tagint *image = atom->image;
//  FixRigidSmall::Body *body = fix_rigid->body;
//  int *atom2body = fix_rigid->atom2body;
//  double **displace = fix_rigid->displace;
//
//  int ibody;
//  double delta[3];
//  FixRigidSmall::Body *b;
//  for (int i = 0; i < nlocal; i++) {
//    ibody = atom2body[i];
//    if (ibody < 0)
//      domain->unmap( x[i], image[i], xu[i] );
//    else {
//      b = &body[ibody];
//      domain->unmap(b->xcm,b->image,xu[i]);
//      MathExtra::matvec(b->ex_space,b->ey_space,b->ez_space,displace[i],delta);
//      xu[i][0] += delta[0];
//      xu[i][1] += delta[1];
//      xu[i][2] += delta[2];
//    }
//  }
//}

/* ----------------------------------------------------------------------
  Restore the positions of rigid bodies whose atoms have 
------------------------------------------------------------------------- */

void FixHMC::rigid_body_restore_positions(double **deltax)
{
//  int i,  ibody;
//  int nlocal = atom->nlocal;
//  double *rmass = atom->rmass;
//  double *mass = atom->mass;
//  int *type = atom->type;
//  double massone;
//
//  int nlocal_body = fix_rigid->nlocal_body;
//  int ntotal_body = nlocal_body + fix_rigid->nghost_body;
//  int *atom2body = fix_rigid->atom2body;
//  FixRigidSmall::Body *body = fix_rigid->body;
//
//  FixRigidSmall::Body *b;
//
//  // Multiply xcm by mass for each local body and zero xcm of each ghost body:
//  for (ibody = 0; ibody < ntotal_body; ibody++) {
//    b = &body[ibody];
//    if (ibody < nlocal_body) {
//      b->xcm[0] *= b->mass;
//      b->xcm[1] *= b->mass;
//      b->xcm[2] *= b->mass;
//    }
//    else
//      b->xcm[0] = b->xcm[1] = b->xcm[2] = 0.0;
//  }
//
//  // Each atom's displacement contributes to the displacement of its containing body:
//  for (i = 0; i < nlocal; i++) {
//    ibody = atom2body[i];
//    if (ibody >= 0) {
//      b = &body[ibody];
//      massone = rmass ? rmass[i] : mass[type[i]];
//      b->xcm[0] -= massone * deltax[i][0];
//      b->xcm[1] -= massone * deltax[i][1];
//      b->xcm[2] -= massone * deltax[i][2];
//    }
//  }
//
//  // Reverse communicate sum(mass*xcm) from ghost to local bodies:
//  comm_flag = XCM;
//  comm->reverse_comm_fix(this,3);
//
//  // For each rigid body, divide xcm by mass:
//  for (ibody = 0; ibody < nlocal_body; ibody++) {
//    b = &body[ibody];
//    b->xcm[0] /= b->mass;
//    b->xcm[1] /= b->mass;
//    b->xcm[2] /= b->mass;
//  }
//
//  // Forward communicate xcm from local to ghost bodies:
//  comm->forward_comm_fix(this,3);
}

/* ----------------------------------------------------------------------
  
------------------------------------------------------------------------- */

void FixHMC::rigid_body_restore_orientations()
{
//  double **x = atom->x;
//  double *rmass = atom->rmass;
//  double *mass = atom->mass;
//  int *type = atom->type;
//  imageint *image = atom->image;
//  int nlocal = atom->nlocal;
//
//  int nlocal_body = fix_rigid->nlocal_body;
//  int ntotal_body = nlocal_body + fix_rigid->nghost_body;
//  int *atom2body = fix_rigid->atom2body;
//  double **displace = fix_rigid->displace;
//  FixRigidSmall::Body *body = fix_rigid->body;
//  imageint *xcmimage = fix_rigid->xcmimage;
//  int ibody;
//  double *it;
//  itensor = new double[ntotal_body][6];
//  for (ibody = 0; ibody < ntotal_body; ibody++) {
//    it = itensor[ibody];
//    it[0] = it[1] = it[2] = it[3] = it[4] = it[5] = 0.0;
//  }
//
//  // Compute moments of inertia in Cartesian reference frame:
//  int i;
//  double massone;
//  double unwrap[3];
//  double dx, dy, dz;
//  double *inertia;
//  FixRigidSmall::Body *b;
//  for (i = 0; i < nlocal; i++) {
//    ibody = atom2body[i];
//    if (ibody >= 0) {
//      b = &body[ibody];
//      domain->unmap(x[i],xcmimage[i],unwrap);
//      dx = unwrap[0] - b->xcm[0];
//      dy = unwrap[1] - b->xcm[1];
//      dz = unwrap[2] - b->xcm[2];
//      massone = rmass ? rmass[i] : mass[type[i]];
//      inertia = itensor[ibody];
//      inertia[0] += massone * (dy*dy + dz*dz);
//      inertia[1] += massone * (dx*dx + dz*dz);
//      inertia[2] += massone * (dx*dx + dy*dy);
//      inertia[3] -= massone * dy*dz;
//      inertia[4] -= massone * dx*dz;
//      inertia[5] -= massone * dx*dy;
//    }
//  }
//
//  // Reverse communicate inertia tensors:
//  comm_flag = ITENSOR;
//  comm->reverse_comm_fix(this,6);
//
//  // Diagonalize inertia tensor for each body via Jacobi rotations:
//  int ierror;
//  double tol;
//  double cross[3], tensor[3][3], evectors[3][3];
//  for (ibody = 0; ibody < nlocal_body; ibody++) {
//    inertia = itensor[ibody];
//    tensor[0][0] = inertia[0];
//    tensor[1][1] = inertia[1];
//    tensor[2][2] = inertia[2];
//    tensor[1][2] = tensor[2][1] = inertia[3];
//    tensor[0][2] = tensor[2][0] = inertia[4];
//    tensor[0][1] = tensor[1][0] = inertia[5];
//    b = &body[ibody];
//    ierror = MathExtra::jacobi(tensor,b->inertia,evectors);
//    if (ierror) error->all(FLERR,"Insufficient Jacobi rotations for rigid body");
//    tol = EPSILON*MAX(MAX(b->inertia[0],b->inertia[1]),b->inertia[2]);
//    for (int j = 0; j < 3; j++) {
//      b->ex_space[j] = evectors[j][0];
//      b->ey_space[j] = evectors[j][1];
//      b->ez_space[j] = evectors[j][2];
//      if (b->inertia[j] < tol) b->inertia[j] = 0.0;
//    }
//    MathExtra::cross3(b->ex_space,b->ey_space,cross);
//    if (MathExtra::dot3(cross,b->ez_space) < 0.0)
//      MathExtra::negate3(b->ez_space);
//    MathExtra::exyz_to_q(b->ex_space,b->ey_space,b->ez_space,b->quat);
//  }
//
//  delete [] itensor;
//
//  // Forward communicate body orientations:
//  comm_flag = ROTATION;
//  comm->forward_comm_fix(this,12);
//
//  // Compute atom coordinates in internal body reference frames:
//  double delta[3];
//  for (i = 0; i < nlocal; i++)
//    if (atom2body[i] < 0)
//      displace[i][0] = displace[i][1] = displace[i][2] = 0.0;
//    else {
//      b = &body[atom2body[i]];
//      domain->unmap(x[i],xcmimage[i],unwrap);
//      delta[0] = unwrap[0] - b->xcm[0];
//      delta[1] = unwrap[1] - b->xcm[1];
//      delta[2] = unwrap[2] - b->xcm[2];
//      MathExtra::transpose_matvec(b->ex_space,b->ey_space,b->ez_space,
//                                  delta,displace[i]);
//    }
}

/* ----------------------------------------------------------------------
  
------------------------------------------------------------------------- */

void FixHMC::rigid_body_restore_forces()
{
//  int nlocal = atom->nlocal;
//  double **x = atom->x;
//  double **f = atom->f;
//  tagint *image = atom->image;
//  int ntotal_body = fix_rigid->nlocal_body + fix_rigid->nghost_body;
//  int *atom2body = fix_rigid->atom2body;
//  FixRigidSmall::Body *body = fix_rigid->body;
//  imageint *xcmimage = fix_rigid->xcmimage;
//  int i, ibody;
//  double lx, ly, lz;
//  double fx, fy, fz;
//  double unwrap[3];
//  FixRigidSmall::Body *b;
//
//  // Zero forces and torques on local and ghost bodies:
//  for (ibody = 0; ibody < ntotal_body; ibody++) {
//    b = &body[ibody];
//    b->fcm[0] = b->fcm[1] = b->fcm[2] = 0.0;
//    b->torque[0] = b->torque[1] = b->torque[2] = 0.0;
//  }
//
//  // The force on each atom contributes to the force and torque on a body:
//  for (i = 0; i < nlocal; i++) {
//    ibody = atom2body[i];
//    if (ibody >= 0) {
//      b = &body[atom2body[i]];
//
//      fx = f[i][0];
//      fy = f[i][1];
//      fz = f[i][2];
//
//      b->fcm[0] += fx;
//      b->fcm[1] += fy;
//      b->fcm[2] += fz;
//
//      domain->unmap(x[i],xcmimage[i],unwrap);
//      lx = unwrap[0] - b->xcm[0];
//      ly = unwrap[1] - b->xcm[1];
//      lz = unwrap[2] - b->xcm[2];
//
//      b->torque[0] += ly*fz - lz*fy;
//      b->torque[1] += lz*fx - lx*fz;
//      b->torque[2] += lx*fy - ly*fx;
//    }
//  }
//
//  // Reverse communicate fcm and torque of all bodies:
//  comm_flag = FORCE_TORQUE;
//  comm->reverse_comm_fix(this,6);
}

/* ----------------------------------------------------------------------
  
------------------------------------------------------------------------- */

void FixHMC::rigid_body_random_velocities()
{
//  FixRigidSmall::Body *body = fix_rigid->body;
//  int nlocal = fix_rigid->nlocal_body;
//  int ntotal = nlocal + fix_rigid->nghost_body;
//
//  double stdev, wbody[3], mbody[3];
//  FixRigidSmall::Body *b;
//
//  for (int ibody = 0; ibody < nlocal; ibody++) {
//    b = &body[ibody];
//    stdev = sqrt(KT/b->mass);
//    for (int j = 0; j < 3; j++) {
//      b->vcm[j] = stdev*random->gaussian();
//      if (b->inertia[j] > 0.0)
//        wbody[j] = sqrt(KT/b->inertia[j])*random->gaussian();
//      else
//        wbody[j] = 0.0;
//    }
//    MathExtra::matvec(b->ex_space,b->ey_space,b->ez_space,wbody,b->omega);
//  }
//
//  // Forward communicate vcm and omega to ghost bodies:
//  comm_flag = VCM_OMEGA;
//  comm->forward_comm_fix(this,6);
//
//  // Compute angular momenta of rigid bodies:
//  for (int ibody = 0; ibody < ntotal; ibody++) {
//    b = &body[ibody];
//    MathExtra::omega_to_angmom(b->omega,b->ex_space,b->ey_space,b->ez_space,
//                               b->inertia,b->angmom);
//    MathExtra::transpose_matvec(b->ex_space,b->ey_space,b->ez_space,
//                                b->angmom,mbody);
//    MathExtra::quatvec(b->quat,mbody,b->conjqm);
//    b->conjqm[0] *= 2.0;
//    b->conjqm[1] *= 2.0;
//    b->conjqm[2] *= 2.0;
//    b->conjqm[3] *= 2.0;
//  }
//
//  // Compute velocities of individual atoms:
//  fix_rigid->set_v();
}

/* ----------------------------------------------------------------------
   Pack rigid body info for forward communication
------------------------------------------------------------------------- */

int FixHMC::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  //int *bodyown = fix_rigid->bodyown;
  //FixRigidSmall::Body *body = fix_rigid->body;

  int i, m, ibody;
  m = 0;
  //FixRigidSmall::Body *b;

  //m = 0;
  //for (i = 0; i < n; i++) {
  //  ibody = bodyown[list[i]];
  //  if (ibody >= 0) {
  //    b = &body[ibody];
  //    if (comm_flag == VCM_OMEGA) {
  //      memcpy( &buf[m], b->vcm, three );
  //      memcpy( &buf[m+3], b->omega, three );
  //      m += 6;
  //    }
  //    else if (comm_flag == XCM) {
  //      memcpy( &buf[m], b->xcm, three );
  //      m += 3;
  //    }
  //    else if (comm_flag == ROTATION) {
  //      memcpy( &buf[m], b->ex_space, three );
  //      memcpy( &buf[m+3], b->ey_space, three );
  //      memcpy( &buf[m+6], b->ez_space, three );
  //      memcpy( &buf[m+9], b->quat, four );
  //      m += 12;
  //    }
  //  }
  //}
  return m;
}

/* ----------------------------------------------------------------------
   Unpack rigid body info from forward communication
------------------------------------------------------------------------- */

void FixHMC::unpack_forward_comm(int n, int first, double *buf)
{
  //int *bodyown = fix_rigid->bodyown;
  //FixRigidSmall::Body *body = fix_rigid->body;

  int i, m, last;
  //FixRigidSmall::Body *b;

  //m = 0;
  //last = first + n;
  //for (i = first; i < last; i++)
  //  if (bodyown[i] >= 0) {
  //    b = &body[bodyown[i]];
  //    if (comm_flag == VCM_OMEGA) {
  //      memcpy( b->vcm, &buf[m], three );
  //      memcpy( b->omega, &buf[m+3], three );
  //      m += 6;
  //    }
  //    else if (comm_flag == XCM) {
  //      memcpy( b->xcm, &buf[m], three );
  //      m += 3;
  //    }
  //    else if (comm_flag == ROTATION) {
  //      memcpy( b->ex_space, &buf[m], three );
  //      memcpy( b->ey_space, &buf[m+3], three );
  //      memcpy( b->ez_space, &buf[m+6], three );
  //      memcpy( b->quat, &buf[m+9], four );
  //      m += 12;
  //    }
  //  }
}

/* ----------------------------------------------------------------------
   Pack per-atom energies and/or virials or
   rigid body info for reverse communication
------------------------------------------------------------------------- */

int FixHMC::pack_reverse_comm(int n, int first, double *buf)
{
  int last = first + n;

  int i, k, m;

  m = 0;
  if (comm_flag == ATOMS) {
    for (i = first; i < last; i++)
      for (k = 0; k < ne; k++)
        if (rev_comm[k]) {
          if (peatom_flag)
            buf[m++] = eatom[k][i];
          if (pressatom_flag) {
            memcpy( &buf[m], vatom[k][i], six );
            m += 6;
          }
        }
  } else {
    //int *bodyown = fix_rigid->bodyown;
    //FixRigidSmall::Body *body = fix_rigid->body;

    //FixRigidSmall::Body *b;

    // if (comm_flag == ITENSOR) {
    //    for (i = first; i < last; i++)
    //      if (bodyown[i] >= 0) {
    //        memcpy( &buf[m], itensor[bodyown[i]], six );
    //        m += 6;
    //      }
    //  }
    //  else if (comm_flag == XCM) {
    //    for (i = first; i < last; i++)
    //      if (bodyown[i] >= 0) {
    //        b = &body[bodyown[i]];
    //        memcpy( &buf[m], b->xcm, three );
    //        m += 3;
    //      }
    //  }
    //  else if (comm_flag == FORCE_TORQUE) {
    //    for (i = first; i < last; i++)
    //      if (bodyown[i] >= 0) {
    //        b = &body[bodyown[i]];
    //        memcpy( &buf[m], b->fcm, three );
    //        memcpy( &buf[m+3], b->torque, three );
    //        m += 6;
    //      }
    //  }
    //}
    //return m;
  }
    return 0;
  }

/* ----------------------------------------------------------------------
   Unpack per-atom energies and/or virials or
   rigid body info for reverse communication
------------------------------------------------------------------------- */

void FixHMC::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, k, m;

  m = 0;
  if (comm_flag == ATOMS) {
    for (j = 0; j < n; j++) {
      i = list[j];
      for (k = 0; k < ne; k++)
        if (rev_comm[k]) {
          if (peatom_flag)
            eatom[k][i] += buf[m++];
          if (pressatom_flag) 
            vatom[k][i][0] += buf[m++];
            vatom[k][i][1] += buf[m++];
            vatom[k][i][2] += buf[m++];
            vatom[k][i][3] += buf[m++];
            vatom[k][i][4] += buf[m++];
            vatom[k][i][5] += buf[m++];
        }
    }
  }
  else {
    //int *bodyown = fix_rigid->bodyown;
    //FixRigidSmall::Body *body = fix_rigid->body;

    //int ibody;
    //FixRigidSmall::Body *b;

  //if (comm_flag == ITENSOR) {
  //    for (i = 0; i < n; i++) {
  //      ibody = bodyown[list[i]];
  //      if (ibody >= 0) {
  //        itensor[ibody][0] += buf[m++];
  //        itensor[ibody][1] += buf[m++];
  //        itensor[ibody][2] += buf[m++];
  //        itensor[ibody][3] += buf[m++];
  //        itensor[ibody][4] += buf[m++];
  //        itensor[ibody][5] += buf[m++];
  //      }
  //    }
  //  }
  //  else if (comm_flag == XCM) {
  //    for (i = 0; i < n; i++) {
  //      ibody = bodyown[list[i]];
  //      if (ibody >= 0) {
  //        b = &body[ibody];
  //        b->xcm[0] += buf[m++];
  //        b->xcm[1] += buf[m++];
  //        b->xcm[2] += buf[m++];
  //      }
  //    }
  //  }
  //  else if (comm_flag == FORCE_TORQUE) {
  //    for (i = 0; i < n; i++) {
  //      ibody = bodyown[list[i]];
  //      if (ibody >= 0) {
  //        b = &body[ibody];
  //        b->fcm[0] += buf[m++];
  //        b->fcm[1] += buf[m++];
  //        b->fcm[2] += buf[m++];
  //        b->torque[0] += buf[m++];
  //        b->torque[1] += buf[m++];
  //        b->torque[2] += buf[m++];
  //      }
  //    }
  //  }
  }
}

/* ----------------------------------------------------------------------
   allocate atom-based arrays
------------------------------------------------------------------------- */

void FixHMC::grow_arrays(int nmax)
{
  memory->grow(xu,nmax,3,"fix_hmc:xu");
  memory->grow(deltax,nmax,3,"fix_hmc:deltax");
  memory->grow(scal,nscal,nmax,"fix_hmc:scal");
  memory->grow(vec,nvec,nmax,3,"fix_hmc:vec");
  memory->grow(eatom,ne,nmax,"fix_hmc:eatom");
  memory->grow(vatom,nv,nmax,6,"fix_hmc:vatom");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixHMC::copy_arrays(int i, int j, int delflag)
{
  int m;

  for (m = 0; m < nscal; m++)
    scal[m][j] = scal[m][i];

  for (m = 0; m < nvec; m++)
    memcpy( vec[m][j], vec[m][i], three );

  if (peatom_flag)
    for (m = 0; m < ne; m++)
      eatom[m][j] = eatom[m][i];

  if (pressatom_flag)
    for (m = 0; m < nv; m++)
      memcpy( vatom[m][j], vatom[m][i], six );

    // make a copy of all peratom data

  for (LAMMPS_NS::Atom::PerAtom &stored_peratom_member : stored_peratom) {
    free(stored_peratom_member.address);
    if (stored_peratom_member.address_maxcols != NULL) free(stored_peratom_member.address_maxcols);
  }
  stored_peratom.clear();
  for (LAMMPS_NS::Atom::PerAtom &current_peratom_member : current_peratom) {
    LAMMPS_NS::Atom::PerAtom stored_peratom_member = current_peratom_member;
    switch (current_peratom_member.datatype) {
      case (Atom::INT):
        store_peratom_member<int>(stored_peratom_member, current_peratom_member, stored_nlocal);
        break;
      case (Atom::DOUBLE):
        store_peratom_member<double>(stored_peratom_member, current_peratom_member, stored_nlocal);
        break;
      case (Atom::BIGINT):
        store_peratom_member<bigint>(stored_peratom_member, current_peratom_member, stored_nlocal);
        break;
    }
    stored_peratom.push_back(stored_peratom_member);
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixHMC::pack_exchange(int i, double *buf)
{
  int k, m = 0;

  for (k = 0; k < nscal; k++)
    buf[m++] = scal[k][i];

  for (k = 0; k < nvec; k++) {
    memcpy( &buf[m], vec[k][i], three );
    m += 3;
  }

  if (peatom_flag)
    for (k = 0; k < ne; k++)
      buf[m++] = eatom[k][i];

  if (pressatom_flag)
    for (k = 0; k < nv; k++) {
      memcpy( &buf[m], vatom[k][i], six );
      m += 6;
    }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixHMC::unpack_exchange(int i, double *buf)
{
  int k, m = 0;

  for (k = 0; k < nscal; k++)
    scal[k][i] = buf[m++];

  for (k = 0; k < nvec; k++) {
    memcpy( vec[k][i], &buf[m], three );
    m += 3;
  }

  if (peatom_flag)
    for (k = 0; k < ne; k++) 
      eatom[k][i] = buf[m++];

  if (pressatom_flag)
    for (k = 0; k < nv; k++) {
      memcpy( vatom[k][i], &buf[m], six );
      m += 6;
    }

  return m;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixHMC::memory_usage()
{
  double bytes = nvalues * atom->nmax * sizeof(double);
  return bytes;
}

