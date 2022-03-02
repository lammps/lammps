// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Eric Bylaska (PNNL)
------------------------------------------------------------------------- */

#include "fix_nwchem.h"
#include <cstdio>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "group.h"
#include "integrate.h"
#include "min.h"
#include "neighbor.h"
#include "pair.h"
#include "kspace.h"
#include "modify.h"
#include "compute.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

namespace pwdft {
using namespace pwdft;
  extern void lammps_pspw_input(MPI_Comm, std::string &);
  extern int lammps_pspw_minimizer(MPI_Comm, double *, double *, 
                                   double *, double *, double *);
  //extern int nwchem_abiversion();
}

// the ABIVERSION number here must be kept consistent
// with its counterpart in the NWChem library and the
// prototype above. We want to catch mismatches with
// a meaningful error messages, as they can cause
// difficult to debug crashes or memory corruption.

#define NWCHEM_ABIVERSION 20220101

// unit conversion factors
// NOTE: set all to 1.0 for dummy test

//#define ANGSTROM_TO_BOHR 1.88973
//#define HARTREE_TO_EV 27.2114
//#define HARTREE_TO_KCAL_MOLE 627.5

#define ANGSTROM_TO_BOHR 1.0
#define HARTREE_TO_EV 1.0
#define HARTREE_TO_KCAL_MOLE 1.0

// mode of operation

enum{AIMD,QMMM};

// prototype for non-class compare function for sorting QM IDs

static int compare_IDs(const int, const int, void *);

/* ---------------------------------------------------------------------- */

FixNWChem::FixNWChem(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // error checks

  //if (NWCHEM_ABIVERSION != pwdft::nwchem_abiversion())
  //  error->all(FLERR,"LAMMPS is linked against incompatible NWChem library");

  if ((strcmp(update->unit_style,"metal") != 0) &&
      (strcmp(update->unit_style,"real") != 0))
    error->all(FLERR,"Fix nwchem requires metal or real units");

  if (domain->dimension == 2)
    error->all(FLERR,"Fix nwchem requires 3d simulation");

  if (!atom->tag_enable) 
    error->all(FLERR,"Fix nwchem requires atom IDs be defined");

  if (!atom->tag_consecutive()) 
    error->all(FLERR,"Fix nwchem requires atom IDs be consecutive");

  if (atom->map_style == Atom::MAP_NONE) 
    error->all(FLERR,"Fix nwchem requires an atom map be defined");

  // NOTE: relax this for slab geometries at some point

  if (domain->nonperiodic == 0) pbcflag = 1;
  else if (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic)
    pbcflag = 0;
  else error->all(FLERR,"Fix nwchem requires fully periodic or "
                  "fully non-periodic system");

  // process command args

  if (narg != 4) error->all(FLERR,"Illegal fix nwchem command");

  int n = strlen(arg[3]) + 1;
  nwfile = new char[n];
  strcpy(nwfile,arg[3]);

  // settings for this fix

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  thermo_energy = 1;

  // setup unit conversion factors

  if (strcmp(update->unit_style,"metal") != 0) { 
    lmp2qm_distance = ANGSTROM_TO_BOHR;
    lmp2qm_energy = 1.0 / HARTREE_TO_EV;
    qm2lmp_force = HARTREE_TO_EV * ANGSTROM_TO_BOHR;
    qm2lmp_energy = HARTREE_TO_EV;
  } else if (strcmp(update->unit_style,"real") != 0) {
    lmp2qm_distance = ANGSTROM_TO_BOHR;
    lmp2qm_energy = 1.0 / HARTREE_TO_KCAL_MOLE;
    qm2lmp_force = HARTREE_TO_KCAL_MOLE * ANGSTROM_TO_BOHR;
    qm2lmp_energy = HARTREE_TO_KCAL_MOLE;
  }

  // qflag = 1 if system stores charge/atom, else 0

  qflag = atom->q_flag;

  // nqm = size of fix group = total # of QM atoms
  // if all atoms are QM, mode = AIMD
  // if only some atoms are QM, mode = QMMM
  // require 3*nqm be a small INT, so can MPI_Allreduce xqm

  bigint ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix nwchem has no atoms in quantum group");
  if (3*ngroup > MAXSMALLINT) 
    error->all(FLERR,"Fix nwchem quantum group has too many atoms");
  nqm = ngroup;

  if (nqm == atom->natoms) mode = AIMD;
  else mode = QMMM;

  memory->create(qmIDs,nqm,"nwchem:qmIDs");
  memory->create(xqm,nqm,3,"nwchem:xqm");
  memory->create(fqm,nqm,3,"nwchem:fqm");
  memory->create(qqm,nqm,"nwchem:qqm");
  memory->create(qpotential,nqm,"nwchem:qpotential");
  memory->create(xqm_mine,nqm,3,"nwchem:xqm_mine");
  memory->create(qqm_mine,nqm,"nwchem:qqm_mine");
  memory->create(qpotential_mine,nqm,"nwchem:qpotential_mine");
  memory->create(qm2owned,nqm,"nwchem:qm2owned");
  
  // for QMMM, set qmIDs = IDs of QM atoms in ascending order

  if (mode == QMMM) {
    tagint *tag = atom->tag;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    // qmIDs_mine = list of nqm_mine QM atom IDs I own

    int nqm_mine = 0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) nqm_mine++;

    tagint *qmIDs_mine;
    memory->create(qmIDs_mine,nqm_mine,"nwchem:qmIDs_mine");

    nqm_mine = 0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) qmIDs_mine[nqm_mine++] = tag[i];

    // allgather of qmIDs_mine into qmIDs

    int nprocs = comm->nprocs;

    int *recvcounts,*displs,*listall;
    memory->create(recvcounts,nprocs,"nwchem:recvcounts");
    memory->create(displs,nprocs,"nwchem:displs");

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

    memory->create(order,nqm,"nwchem:order");
    memory->create(qmIDs_sort,nqm,"nwchem:qmIDs_sort");

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

  // nwchem_setup() perform one-time call to lammps_pspw_input() with nwfile

  nwchem_input();

  // flag for one-time init of qqm = charge on QM atoms

  qqm_init = 0;

  // peratom Coulombic energy

  ecoul = nullptr;
  ncoulmax = 0;

  // set QM energy = 0.0 in case accessed on step 0

  qmenergy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixNWChem::~FixNWChem()
{
  delete [] nwfile;

  memory->destroy(qmIDs);
  memory->destroy(xqm);
  memory->destroy(fqm);
  memory->destroy(qqm);
  memory->destroy(qpotential);
  memory->destroy(xqm_mine);
  memory->destroy(qqm_mine);
  memory->destroy(qpotential_mine);
  memory->destroy(qm2owned);
  memory->destroy(ecoul);
}

/* ---------------------------------------------------------------------- */

int FixNWChem::setmask()
{
  int mask = 0;
  mask |= POST_NEIGHBOR;
  mask |= MIN_POST_NEIGHBOR;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNWChem::nwchem_input()
{
  std::string filename = nwfile;

  // read nwfile
  // replace simulation box with new box
  // replace list of atoms with new list of atoms and element names
  //   specify on command line, to match LAMMPS types as for pair styles
  // write out w2.AUTO.nw
  // pass that filename to NWChem

  dummy_pspw_input(world,filename);

  //pwdft::lammps_pspw_input(world,filename);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::init()
{
  // QMMM requires long-range Coulombics for Coulomb potential input to NWChem

  if (mode == QMMM) {
    if (!qflag) 
      error->all(FLERR,"Fix nwchem QMMM mode requires per-atom charge");
    
    if (!force->pair) 
      error->all(FLERR,"Fix nwchem QMMM mode requires a pair style");

    // must be a pair style that calculates only Coulombic interactions
    // can be in conjunction with KSpace solver as well

    pair_coul = force->pair_match("coul/cut",1,0);
    if (!pair_coul) pair_coul = force->pair_match("coul/long",1,0);
    if (!pair_coul) pair_coul = force->pair_match("coul/msm",1,0);
    if (!pair_coul) 
      error->all(FLERR,
                 "Fix nwchem QMMM mode requires Coulomb-only pair sub-style");
  }

  // NOTE: is this needed ?
  // c_pe = instance of compute pe/atom
  
  //int ipe = modify->find_compute(id_pe);
  //if (ipe < 0) error->all(FLERR,"Could not find fix nwchem compute ID");
  //c_pe = modify->compute[ipe];
}

/* ---------------------------------------------------------------------- */

void FixNWChem::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::setup_post_neighbor()
{
  post_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixNWChem::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::post_neighbor()
{
  // qm2owned[i] = index of local atom for each of nqm QM atoms
  // for AIMD, IDs of QM atoms = 1 to Natoms
  // for QMMM, IDs of QM atoms stored in qmIDs
  // index = -1 if this proc does not own the atom
  // NOTE: could store nqm_owned with different qm2owned indexing

  int nlocal = atom->nlocal;
  int index;

  for (int i = 0; i < nqm; i++) {
    if (mode == AIMD) index = atom->map(i+1);
    else index = atom->map(qmIDs[i]);
    if (index >= nlocal) qm2owned[i] = -1;
    else qm2owned[i] = index;
  }

  // onetime setup of qqm = atom charges for all QM atoms
  // calls to NWChem change it for all QM atoms
  // changes only get copied to LAMMPS atoms, never changed by LAMMPS

  if (!qqm_init) {
    for (int i = 0; i < nqm; i++) qqm_mine[i] = 0.0;

    double *q = atom->q;
    int ilocal;

    for (int i = 0; i < nqm; i++) {
      ilocal = qm2owned[i];
      if (ilocal >= 0) qqm_mine[i] = q[ilocal];
    }

    MPI_Allreduce(qqm_mine,qqm,nqm,MPI_DOUBLE,MPI_SUM,world);
    qqm_init = 1;
  }
}

/* ---------------------------------------------------------------------- */

void FixNWChem::pre_force(int vflag)
{
  if (mode == QMMM) pre_force_qmmm(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::pre_force_qmmm(int vflag)
{
  int ilocal,jlocal;
  double rsq;
  double delta[3];

  // invoke pair hybrid sub-style pair coul/long and Kspace directly
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
    memory->create(ecoul,ncoulmax,"nwchem:ecoul");
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

  // setup NWChem inputs
  // xqm = atom coords, mapped into periodic box
  //   set for owned atoms, then MPI_Allreduce

  for (int i = 0; i < nqm; i++) {
    xqm_mine[i][0] = 0.0;
    xqm_mine[i][1] = 0.0;
    xqm_mine[i][2] = 0.0;
  }

  double **x = atom->x;

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

  // qpotential[i] = (eatom[i] from pair_coul + kspace) / Qi
  //   set for owned atoms, then MPI_Allreduce
  // subtract Qj/Rij energy for QM I interacting with all other QM J atoms
  //   use xqm_mine and qqm_mine for all QM atoms

  for (int i = 0; i < nqm; i++) qpotential_mine[i] = 0.0;

  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      if (q[ilocal] == 0.0) qpotential_mine[i] = 0.0;
      else qpotential_mine[i] = ecoul[ilocal] / q[ilocal];
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

  // unit conversion from LAMMPS to NWChem

  for (int i = 0; i < nqm; i++) {
    xqm[i][0] *= lmp2qm_distance;
    xqm[i][1] *= lmp2qm_distance;
    xqm[i][2] *= lmp2qm_distance;
    qpotential[i] *= lmp2qm_energy;
  }
  // call to NWChem with only QM atom info
  // QM atoms must be in order of ascending atom ID
  // inputs:
  //   xqm = atom coords
  //   qpotential = vector of zeroes for AIMD
  // outputs:
  //   fqm,qqm = forces & charges
  //   qmenergy = QM energy of entire system

  int nwerr = dummy_pspw_minimizer(world,
                                   &xqm[0][0],qpotential,
                                   &fqm[0][0],qqm,&qmenergy);

  if (nwerr) error->all(FLERR,"Internal NWChem error");

  // unit conversion from NWChem to LAMMPS

  qmenergy *= qm2lmp_energy;

  for (int i = 0; i < nqm; i++) {
    fqm[i][0] *= qm2lmp_force;
    fqm[i][1] *= qm2lmp_force;
    fqm[i][2] *= qm2lmp_force;
  }

  // reset owned charges to QM values

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    q[ilocal] = qqm[i];
  }

  // reset LAMMPS forces to zero
  // NOTE: what about check in force_clear() for external_force_clear = OPENMP ?
  // NOTE: what will whichflag be for single snapshot compute of QM forces?

  if (update->whichflag == 1) 
    update->integrate->force_clear();  
  else if (update->whichflag == 2) 
    update->minimize->force_clear();  
}

/* ---------------------------------------------------------------------- */

void FixNWChem::post_force(int vflag)
{
  if (mode == AIMD) post_force_aimd(vflag);
  else post_force_qmmm(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::post_force_aimd(int vflag)
{
  int ilocal;

  // setup 2 NWChem inputs
  // xqm = atom coords, mapped into periodic box
  //   set for owned atoms, then MPI_Allreduce
  // qpotential = Coulomb potential = 0.0 for AIMD

  for (int i = 0; i < nqm; i++) {
    xqm_mine[i][0] = 0.0;
    xqm_mine[i][1] = 0.0;
    xqm_mine[i][2] = 0.0;
    qpotential[i] = 0.0;
  }

  double **x = atom->x;

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

  // unit conversion from LAMMPS to NWChem

  for (int i = 0; i < nqm; i++) {
    xqm[i][0] *= lmp2qm_distance;
    xqm[i][1] *= lmp2qm_distance;
    xqm[i][2] *= lmp2qm_distance;
  }

  // call to NWChem with only QM atom info
  // QM atoms must be in order of ascending atom ID
  // inputs:
  //   xqm = atom coords
  //   qpotential = vector of zeroes for AIMD
  // outputs:
  //   fqm,qqm = forces & charges
  //   qmenergy = QM energy of entire system

  int nwerr = dummy_pspw_minimizer(world,
                                   &xqm[0][0],qpotential,
                                   &fqm[0][0],qqm,&qmenergy);

  if (nwerr) error->all(FLERR,"Internal NWChem error");

  // unit conversion from NWChem to LAMMPS

  qmenergy *= qm2lmp_energy;

  for (int i = 0; i < nqm; i++) {
    fqm[i][0] *= qm2lmp_force;
    fqm[i][1] *= qm2lmp_force;
    fqm[i][2] *= qm2lmp_force;
  }

  // add QM forces to owned atoms
  // if qflag: reset owned charges to QM values
  //   allows user to output charges if desired

  double **f = atom->f;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      f[ilocal][0] += fqm[i][0];
      f[ilocal][1] += fqm[i][1];
      f[ilocal][2] += fqm[i][2];
      if (qflag) q[ilocal] = qqm[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNWChem::post_force_qmmm(int vflag)
{
  int ilocal,jlocal;
  double rsq,r2inv,rinv,fpair;
  double delta[3];

  // adjust LAMMPS energy and forces for pairwise interactions between QM atoms
  // double loop over all QM pairs
  // if I own either or both atoms in pair, compute pairwise term
  // subtract half or all energy from qmenergy
  // subtract force from only owned atoms
  // uses xqm and qqm set or computed in pre_force_qmmm() for all QM atoms

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int nlocal = atom->nlocal;
  double qqrd2e = force->qqrd2e;

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
      qmenergy -= efactor * qqrd2e * qqm[i]*qqm[j]*rinv;

    }
  }

  // add NWChem QM forces to owned QM atoms

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      f[ilocal][0] += fqm[i][0];
      f[ilocal][1] += fqm[i][1];
      f[ilocal][2] += fqm[i][2];
    }
  }

  // add QM energy from NWChem
  // NOTE: how to calculate this quantity
  // NOTE: for now, just letting dummy method return it as qmenergy

  // trigger per-atom energy computation on next step by pair/kspace
  // NOTE: is this needed ?
  //       only if needed for this fix to calc per-atom forces
  //       or needed for this fix to output global (or per-atom) energy

  //c_pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::min_post_neighbor()
{
  post_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixNWChem::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

int FixNWChem::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = ecoul[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void FixNWChem::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    ecoul[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   QM energy from NWChem
------------------------------------------------------------------------- */

double FixNWChem::compute_scalar()
{
  return qmenergy;
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double FixNWChem::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)ncoulmax * sizeof(double);
  bytes += (double)(3*3 + 4) * nqm * sizeof(double);  // fpoint QM arrays/vecs
  bytes += nqm * sizeof(tagint);    // qmIDs
  bytes += nqm * sizeof(int);       // qm2owned
  return bytes;
}

/* ----------------------------------------------------------------------
   dummy version of NWChem pwdft_minimizer()
   for dummy test, quantities are in LAMMPS units, not NWChem units

------------------------------------------------------------------------- */

int FixNWChem::dummy_pspw_minimizer(MPI_Comm nwworld, 
                                    double *x, double *qp, 
                                    double *f, double *q, double *eqm)
{
  double rsq,r2inv,rinv,fpair;
  double delta[3];

  // no change to charge in dummy test

  int ilocal;

  if (qflag) {
    for (int i = 0; i < nqm; i++) {
      ilocal = qm2owned[i];
      if (ilocal >= 0) q[i] = atom->q[ilocal];
    }
  }

  // all procs compute entire QM energy/forces

  double eng = 0.0;
  for (int i = 0; i < 3*nqm; i++) f[i] = 0.0;

  double qqrd2e = force->qqrd2e;

  for (int i = 0; i < nqm; i++) {
    for (int j = i+1; j < nqm; j++) {
      delta[0] = x[3*i+0] - x[3*j+0];
      delta[1] = x[3*i+1] - x[3*j+1];
      delta[2] = x[3*i+2] - x[3*j+2];
      domain->minimum_image_once(delta);
      rsq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
      r2inv = 1.0/rsq;
      rinv = sqrt(r2inv);

      eng += qqrd2e * qqm[i]*qqm[j]*rinv;
      fpair = qqrd2e * qqm[i]*qqm[j]*rinv*r2inv;
      
      f[3*i+0] += delta[0]*fpair;
      f[3*i+1] += delta[1]*fpair;
      f[3*i+2] += delta[2]*fpair;
      f[3*j+0] -= delta[0]*fpair;
      f[3*j+1] -= delta[1]*fpair;
      f[3*j+2] -= delta[2]*fpair;
    }
  }

  // return energy as last function argument

  *eqm = eng;

  return 0;
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
