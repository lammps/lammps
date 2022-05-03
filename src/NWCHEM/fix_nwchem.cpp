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

extern void lammps_pspw_input(MPI_Comm, std::string &);
extern int lammps_pspw_aimd_minimizer(MPI_Comm, double *, double *, double *);
extern int lammps_pspw_qmmm_minimizer(MPI_Comm, double *, double *, 
                                      double *, double *, double *);

//extern int nwchem_abiversion();

// the ABIVERSION number here must be kept consistent
// with its counterpart in the NWChem library and the
// prototype above. We want to catch mismatches with
// a meaningful error messages, as they can cause
// difficult to debug crashes or memory corruption.

#define NWCHEM_ABIVERSION 20220101

// unit conversion factors
// NOTE: set all to 1.0 for dummy test

#define ANGSTROM_TO_BOHR 1.88973
#define HARTREE_TO_EV 27.2114
#define HARTREE_TO_KCAL_MOLE 627.5

//#define ANGSTROM_TO_BOHR 1.0
//#define HARTREE_TO_EV 1.0
//#define HARTREE_TO_KCAL_MOLE 1.0

#define MAXLINE 1024

// mode of operation

enum{AIMD,QMMM};

// prototype for non-class compare function for sorting QM IDs

static int compare_IDs(const int, const int, void *);

/* ---------------------------------------------------------------------- */

FixNWChem::FixNWChem(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // error checks

  //if (NWCHEM_ABIVERSION != nwchem_abiversion())
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

  // NOTE: allow slab geometries at some point

  if (domain->nonperiodic == 0) pbcflag = 1;
  else if (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic)
    pbcflag = 0;
  else error->all(FLERR,"Fix nwchem requires fully periodic or "
                  "fully non-periodic system");

  // process command args
  // trailing args map atom types to element strings

  if (narg != 5 + atom->ntypes) 
    error->all(FLERR,"Illegal fix nwchem command");

  int n = strlen(arg[3]) + 1;
  nw_template = new char[n];
  strcpy(nw_template,arg[3]);

  n = strlen(arg[4]) + 1;
  nw_input = new char[n];
  strcpy(nw_input,arg[4]);

  elements = new char*[narg-5];

  for (int i = 5; i < narg; i++) {
    n = strlen(arg[i]) + 1;
    elements[i-5] = new char[n];
    strcpy(elements[i-5],arg[i]);
  }

  // settings for this fix

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  thermo_energy = 1;
  comm_reverse = 1;

  // setup unit conversion factors

  if (strcmp(update->unit_style,"metal") == 0) { 
    lmp2qm_distance = ANGSTROM_TO_BOHR;
    lmp2qm_energy = 1.0 / HARTREE_TO_EV;
    qm2lmp_force = HARTREE_TO_EV * ANGSTROM_TO_BOHR;
    qm2lmp_energy = HARTREE_TO_EV;
  } else if (strcmp(update->unit_style,"real") == 0) {
    lmp2qm_distance = ANGSTROM_TO_BOHR;
    lmp2qm_energy = 1.0 / HARTREE_TO_KCAL_MOLE;
    qm2lmp_force = HARTREE_TO_KCAL_MOLE * ANGSTROM_TO_BOHR;
    qm2lmp_energy = HARTREE_TO_KCAL_MOLE;
  }

  // qflag = 1 if system stores per-atom charge, else 0

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
  else {
    mode = QMMM;
    error->all(FLERR,"Fix nwchem for a non-all group (QMMM) is not yet enabled");
  }

  memory->create(qmIDs,nqm,"nwchem:qmIDs");
  memory->create(xqm,nqm,3,"nwchem:xqm");
  memory->create(fqm,nqm,3,"nwchem:fqm");
  memory->create(qqm,nqm,"nwchem:qqm");
  memory->create(tqm,nqm,"nwchem:tqm");
  memory->create(qpotential,nqm,"nwchem:qpotential");
  memory->create(xqm_mine,nqm,3,"nwchem:xqm_mine");
  memory->create(qqm_mine,nqm,"nwchem:qqm_mine");
  memory->create(tqm_mine,nqm,"nwchem:tqm_mine");
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

  // flag for one-time init of NWChem and qqm

  qm_init = 0;

  // peratom Coulombic energy

  ecoul = nullptr;
  ncoulmax = 0;

  // set QM energy = 0.0 in case accessed on step 0

  qmenergy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixNWChem::~FixNWChem()
{
  delete [] nw_template;
  delete [] nw_input;

  for (int i = 0; i < atom->ntypes; i++)
    delete [] elements[i];
  delete [] elements;

  memory->destroy(qmIDs);
  memory->destroy(xqm);
  memory->destroy(fqm);
  memory->destroy(qqm);
  memory->destroy(tqm);
  memory->destroy(qpotential);
  memory->destroy(xqm_mine);
  memory->destroy(qqm_mine);
  memory->destroy(tqm_mine);
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
  if (mode == QMMM) {
    mask |= PRE_FORCE;
    mask |= MIN_PRE_FORCE;
  }
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
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

  // one-time initialization of NWChem with its input file
  // also one-time setup of qqm = atom charges for all QM atoms
  //   later calls to NWChem change it for all QM atoms
  //   changes get copied to LAMMPS per-atom q which is not changed by LAMMPS
  // setup of xqm needed for nwchem_initialize();

  if (!qm_init) {
    qm_init = 1;
    set_qm2owned();
    set_qqm();
    set_tqm();
    set_xqm();
    if (comm->me == 0) nwchem_initialize();
    MPI_Barrier(world);

    std::string filename = nw_input;
    lammps_pspw_input(world,filename);
    //dummy_pspw_input(world,filename);
  }
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
  set_qm2owned();
}

/* ---------------------------------------------------------------------- */

void FixNWChem::pre_force(int vflag)
{
  pre_force_qmmm(vflag);
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

  // setup 2 NWChem inputs: xqm and qpotential
  // xqm = atom coords, mapped into periodic box
  // qpotential[i] = Coulomb potential for each atom
  //   (eatom[i] from pair_coul + kspace) / Qi
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

  int nwerr = lammps_pspw_qmmm_minimizer(world,&xqm[0][0],qpotential,
                                         &fqm[0][0],qqm,&qmenergy);
  //int nwerr = dummy_pspw_qmmm_minimizer(world,&xqm[0][0],qpotential,
  //                                      &fqm[0][0],qqm,&qmenergy);

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
    if (ilocal >= 0) q[ilocal] = qqm[i];
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

  // setup NWChem input: xqm
  // xqm = atom coords, mapped into periodic box

  set_xqm();

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
  // outputs:
  //   fqm = forces
  //   qmenergy = QM energy of entire system

  int nwerr = lammps_pspw_aimd_minimizer(world,&xqm[0][0],&fqm[0][0],&qmenergy);

  if (nwerr) error->all(FLERR,"Internal NWChem error");

  // unit conversion from NWChem to LAMMPS

  qmenergy *= qm2lmp_energy;

  for (int i = 0; i < nqm; i++) {
    fqm[i][0] *= qm2lmp_force;
    fqm[i][1] *= qm2lmp_force;
    fqm[i][2] *= qm2lmp_force;
  }

  // add QM forces to owned atoms

  double **f = atom->f;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) {
      f[ilocal][0] += fqm[i][0];
      f[ilocal][1] += fqm[i][1];
      f[ilocal][2] += fqm[i][2];
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

// ----------------------------------------------------------------------
// private methods for this fix
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   one-time intialization of NWChem
   only called by proc 0
------------------------------------------------------------------------- */

void FixNWChem::nwchem_initialize()
{
  char *eof;
  char line[MAXLINE];
  FILE *nwfile,*infile;

  // open nw_template file read from, nw_input to write to

  nwfile = fopen(nw_template,"r");
  if (!nwfile) 
    error->one(FLERR,"Cannot open fix nwchem template file {}: {}",
               nw_template,utils::getsyserror());

  infile = fopen(nw_input,"w");
  if (!infile) 
    error->one(FLERR,"Cannot open fix nwchem NWChem input file {}: {}",
               nw_input,utils::getsyserror());
  
  // read nw_template file line by line and echo to nw_input file
  // when FILEINSERT line is read:
  //   replace it with line for filename prefix for NWChem output
  // when GEOMINSERT line is read:
  //   replace it with lines that specify the simulation box and list of atoms
  //   each atom's element ID and coords are written to file

  while (true) {
    eof = fgets(line,MAXLINE,nwfile);
    if (eof == nullptr) break;

    if (strcmp(line,"FILEINSERT\n") == 0) {
      int n = strlen(nw_input);
      char *fileprefix = new char[n];
      strncpy(fileprefix,nw_input,n-3);
      fileprefix[n-3] = '\0';
      fprintf(infile,"start %s\n",fileprefix);
      delete [] fileprefix;
      continue;
    }
    
    if (strcmp(line,"GEOMINSERT\n") == 0) {
      fprintf(infile,"geometry noautosym noautoz nocenter\n");
      fprintf(infile,"system crystal cartesian\n");
      fprintf(infile,"lattice_vectors\n");

      int *type = atom->type;

      // orthogonal box
    
      if (!domain->triclinic) {
        fprintf(infile,"%g %g %g\n",domain->xprd,0.0,0.0);
        fprintf(infile,"%g %g %g\n",0.0,domain->yprd,0.0);
        fprintf(infile,"%g %g %g\n",0.0,0.0,domain->zprd);
        fprintf(infile,"end\n\n");
        
        for (int i = 0; i < nqm; i++) {
          fprintf(infile,"%s %g %g %g\n",
                  elements[tqm[i]-1],xqm[i][0],xqm[i][1],xqm[i][2]);
        }
        fprintf(infile,"end\n");
        
        // triclinic box
        
      } else {
        fprintf(infile,"%g %g %g\n",domain->xprd,0.0,0.0);
        fprintf(infile,"%g %g %g\n",domain->xy,domain->yprd,0.0);
        fprintf(infile,"%g %g %g\n",domain->xz,domain->yz,domain->zprd);
        fprintf(infile,"end\n\n");
        
        for (int i = 0; i < nqm; i++)
          fprintf(infile,"%s %g %g %g\n",
                  elements[tqm[i]-1],xqm[i][0],xqm[i][1],xqm[i][2]);
        fprintf(infile,"end\n");
      }

      continue;
    }

    // just echo all other lines

    fprintf(infile,line);
  }

  // close 2 files

  fclose(nwfile);
  fclose(infile);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::set_qm2owned()
{
  // qm2owned[i] = index of local atom for each of nqm QM atoms
  // for AIMD, IDs of QM atoms = 1 to Natoms
  // for QMMM, IDs of QM atoms stored in qmIDs
  // index = -1 if this proc does not own the atom

  int nlocal = atom->nlocal;
  int index;

  for (int i = 0; i < nqm; i++) {
    if (mode == AIMD) index = atom->map(i+1);
    else index = atom->map(qmIDs[i]);
    if (index >= nlocal) qm2owned[i] = -1;
    else qm2owned[i] = index;
  }
}

/* ---------------------------------------------------------------------- */

void FixNWChem::set_qqm()
{
  for (int i = 0; i < nqm; i++) qqm_mine[i] = 0.0;

  if (qflag) {
    double *q = atom->q;
    int ilocal;

    for (int i = 0; i < nqm; i++) {
      ilocal = qm2owned[i];
      if (ilocal >= 0) qqm_mine[i] = q[ilocal];
    }
  }
  
  MPI_Allreduce(qqm_mine,qqm,nqm,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::set_tqm()
{
  for (int i = 0; i < nqm; i++) tqm_mine[i] = 0;

  int*type = atom->type;
  int ilocal;
  
  for (int i = 0; i < nqm; i++) {
    ilocal = qm2owned[i];
    if (ilocal >= 0) tqm_mine[i] = type[ilocal];
  }
  
  MPI_Allreduce(tqm_mine,tqm,nqm,MPI_INT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::set_xqm()
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

/* ----------------------------------------------------------------------
   dummy version of NWChem pwdft_minimizer()
   for dummy test, quantities are in LAMMPS units, not NWChem units
------------------------------------------------------------------------- */

int FixNWChem::dummy_pspw_qmmm_minimizer(MPI_Comm nwworld, 
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
