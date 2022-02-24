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
  //extern char *util_date();
  //extern void seconds(double *);
  //extern int pspw_geovib(MPI_Comm, std::string&);
  //extern int pspw_minimizer(MPI_Comm, double *, double *, double *, double *);
  //extern int nwchem_abiversion();
}

// the ABIVERSION number here must be kept consistent
// with its counterpart in the NWChem library and the
// prototype above. We want to catch mismatches with
// a meaningful error messages, as they can cause
// difficult to debug crashes or memory corruption.

#define NWCHEM_ABIVERSION 20220101

// unit conversion factors

#define ANGSTROM_TO_BOHR 1.88973
#define HARTREE_TO_EV 27.2114
#define HARTREE_TO_KCAL_MOLE 627.5

/* ---------------------------------------------------------------------- */

FixNWChem::FixNWChem(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if ((strcmp(update->unit_style,"metal") != 0) &&
      (strcmp(update->unit_style,"real") != 0))
    error->all(FLERR,"Must use units metal or real with fix nwchem command");

  // NOTE: change 1-proc restriction later

  if (comm->nprocs != 1)
    error->all(FLERR,"Fix nwchem currently runs only in serial");

  //if (NWCHEM_ABIVERSION != pwdft::nwchem_abiversion())
  //  error->all(FLERR,"LAMMPS is linked against incompatible NWChem library");

  if (narg != 3) error->all(FLERR,"Illegal fix nwchem command");

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
  } else if (strcmp(update->unit_style,"real") != 0) {
    lmp2qm_distance = ANGSTROM_TO_BOHR;
    lmp2qm_energy = 1.0 / HARTREE_TO_KCAL_MOLE;
    qm2lmp_force = HARTREE_TO_KCAL_MOLE * ANGSTROM_TO_BOHR;
  }

  // initializations

  nqm = 0;
  qmIDs = nullptr;
  xqm = fqm = nullptr;
  qpotential = qqm = nullptr;
  qm2lmp = nullptr;

  qmenergy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixNWChem::~FixNWChem()
{
  memory->destroy(qmIDs);
  memory->destroy(xqm);
  memory->destroy(qpotential);
  memory->destroy(fqm);
  memory->destroy(qqm);
  memory->destroy(qm2lmp);
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

void FixNWChem::init()
{
  // error checks

  if (domain->dimension == 2)
    error->all(FLERR,"Fix nwchem requires 3d simulation");

  if (domain->nonperiodic == 0) pbcflag = 1;
  else if (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic)
    pbcflag = 0;
  else error->all(FLERR,"Fix nwchem requires fully periodic or "
                  "fully non-periodic system");
  
  if (atom->q_flag == 0 || force->pair == nullptr || force->kspace == nullptr)
    error->all(FLERR,"Fix nwchem cannot compute Coulomb potential");

  // NOTE: is this needed ?
  // c_pe = instance of compute pe/atom
  
  //int ipe = modify->find_compute(id_pe);
  //if (ipe < 0) error->all(FLERR,"Could not find fix nwchem compute ID");
  //c_pe = modify->compute[ipe];

  // NOTE: could be another coul pair style like coul/msm ?
  
  if (!force->pair) error->all(FLERR,"Fix nwchem requires a pair style");

  pair_coul = force->pair_match("coul/long",1,0);
  if (!pair_coul) 
    error->all(FLERR,"Fix nwchem requires pair sub-style coul/long");

  // fix group = QM atoms
  // one-time initialization of qmIDs
  // NOTE: need to sort in ascending order to match NWChem ?
  
  if (nqm == 0) {
    bigint ngroup = group->count(igroup);
    if (ngroup == 0) error->all(FLERR,"Fix nwchem has no atoms in quantum group");
    if (ngroup > MAXSMALLINT) 
      error->all(FLERR,"Fix nwchem quantum group has too many atoms");
    nqm = ngroup;
    
    memory->create(qmIDs,nqm,"nwchem:qmIDs");
    memory->create(xqm,nqm,3,"nwchem:xqm");
    memory->create(qpotential,nqm,"nwchem:qpotential");
    memory->create(fqm,nqm,3,"nwchem:fqm");
    memory->create(qqm,nqm,"nwchem:qqm");
    memory->create(qm2lmp,nqm,"nwchem:qm2lmp");

    tagint *tag = atom->tag;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    nqm = 0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) qmIDs[nqm++] = tag[i];
    }
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
  for (int i = 0; i < nqm; i++)
    qm2lmp[i] = atom->map(qmIDs[i]);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::pre_force(int vflag)
{
  int ilocal,jlocal;
  double rsq;
  double delta[3];

  // invoke pair hybrid sub-style pair coul/long and Kspace directly
  // set eflag = 2 so they calculate per-atom energy
  // NOTE: need to comm ghost per-atom energy (serial or parallel)

  pair_coul->compute(2,0);
  double *eatom_pair = pair_coul->eatom;

  double *eatom_kspace = nullptr;
  if (force->kspace) {
    force->kspace->compute(2,0);
    eatom_kspace = force->kspace->eatom;
  }

  // create 2 NWChem inputs: xqm and qpotential (Coulomb potential)
  // qpotential[i] = (eatom[i] from pair_coul + kspace) / Qi
  // subtract out Qj/Rij energy for QM I interacting with all other QM J atoms
 
  double **x = atom->x;
  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

  // NOTE: this I,J loop will not work in parallel

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2lmp[i];

    // NOTE: what if LAMMPS atom moves via PBC to other side of LAMMPS box ?
    //       will that mess up NWChem ?

    xqm[i][0] = x[ilocal][0];
    xqm[i][1] = x[ilocal][1];
    xqm[i][2] = x[ilocal][2];

    if (q[ilocal] == 0.0) qpotential[i] = 0.0;
    else if (!force->kspace)
      qpotential[i] = eatom_pair[ilocal] / q[ilocal];
    else
      qpotential[i] = (eatom_pair[ilocal] + eatom_kspace[ilocal]) / q[ilocal];

    for (int j = 0; j < nqm; j++) {
      if (j == i) continue;
      jlocal = qm2lmp[j];
      delta[0] = x[i][0] - x[j][0];
      delta[1] = x[i][1] - x[j][1];
      delta[2] = x[i][2] - x[j][2];
      domain->minimum_image_once(delta);
      rsq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
      qpotential[i] -= qqrd2e * q[jlocal] / sqrt(rsq);
    }
  }

  // unit conversion from LAMMPS to NWChem
  // NOTE: comment out for dummy test

  /*
  for (int i = 0; i < nqm; i++) {
    xqm[i][0] *= lmp2qm_distance;
    xqm[i][1] *= lmp2qm_distance;
    xqm[i][2] *= lmp2qm_distance;
    qpotential[i] *= lmp2qm_energy;
  }
  */

  // call to NWChem
  // input/outputs are only for QM atoms
  // inputs = xqm and qpotential
  // outputs = fqm and qqm
  // NOTE: extra args for dummy test = nqm,qmenergy
  // NOTE: is QM eng returned by NWChem, or needs to be computed by LAMMPS?

  int nwerr = pspw_minimizer(world,nqm,
                             &xqm[0][0],qpotential,&fqm[0][0],qqm,qmenergy);
  if (nwerr) error->all(FLERR,"Internal NWChem error");

  //latte(flags,&natoms,coords,type,&ntypes,mass,boxlo,boxhi,&domain->xy,
  //      &domain->xz,&domain->yz,forces,&maxiter,&latte_energy,
  //      &atom->v[0][0],&update->dt,virial,&newsystem,&latteerror);

  // unit conversion from NWChem to LAMMPS
  // NOTE: comment out for dummy test

  /*
  for (int i = 0; i < nqm; i++) {
    fqm[i][0] *= qm2lmp_force;
    fqm[i][1] *= qm2lmp_force;
    fqm[i][2] *= qm2lmp_force;
  }
  */

  // reset Q of QM atoms

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2lmp[i];
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
  int ilocal,jlocal;
  double rsq,r2inv,rinv,fpair;
  double delta[3];

  // subtract QiQj/Rij energy for pairwise interactions between QM atoms
  // subtract QiQj/Rij^2 force for pairwise interactions between QM atoms

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2lmp[i];
    xqm[i][0] = x[ilocal][0];
    xqm[i][1] = x[ilocal][1];
    xqm[i][2] = x[ilocal][2];
    qqm[i] = q[ilocal];
  }

  qmenergy = 0.0;
  
  for (int i = 0; i < nqm; i++) {
    ilocal = qm2lmp[i];
    for (int j = i+1; j < nqm; j++) {
      jlocal = qm2lmp[j];
      
      delta[0] = xqm[i][0] - xqm[j][0];
      delta[1] = xqm[i][1] - xqm[j][1];
      delta[2] = xqm[i][2] - xqm[j][2];
      domain->minimum_image_once(delta);
      rsq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
      r2inv = 1.0/rsq;
      rinv = sqrt(r2inv);

      qmenergy -= qqrd2e * qqm[i]*qqm[j]*rinv;
      fpair = qqrd2e * qqm[i]*qqm[j]*rinv*r2inv;

      f[ilocal][0] -= delta[0]*fpair;
      f[ilocal][1] -= delta[1]*fpair;
      f[ilocal][2] -= delta[2]*fpair;
      f[jlocal][0] += delta[0]*fpair;
      f[jlocal][1] += delta[1]*fpair;
      f[jlocal][2] += delta[2]*fpair;
    }
  }

  // add QM forces from NWChem to QM atoms

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2lmp[i];
    f[ilocal][0] += fqm[i][0];
    f[ilocal][1] += fqm[i][1];
    f[ilocal][2] += fqm[i][2];
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

/* ----------------------------------------------------------------------
   QM energy from NWChem
------------------------------------------------------------------------- */

double FixNWChem::compute_scalar()
{
  return qmenergy;
}

/* ----------------------------------------------------------------------
   dummy version of NWChem pwdft_minimizer()
   NOTE: for dummy test, quantities are in LAMMPS units, not NWChem units
------------------------------------------------------------------------- */

int FixNWChem::pspw_minimizer(MPI_Comm nwworld, 
                              int n, double *x, double *qp, 
                              double *f, double *q, double &eqm)
{
  double rsq,r2inv,rinv,fpair;
  double delta[3];

  // no change to charge in dummy test
  // NOTE: could test changing of charges

  for (int i = 0; i < nqm; i++) {
    int ilocal = qm2lmp[i];
    q[i] = atom->q[ilocal];
  }

  double qqrd2e = force->qqrd2e;

  eqm = 0.0;
  for (int i = 0; i < 3*nqm; i++) f[i] = 0.0;
  
  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      delta[0] = x[3*i+0] - x[3*j+0];
      delta[1] = x[3*i+1] - x[3*j+1];
      delta[2] = x[3*i+2] - x[3*j+2];
      domain->minimum_image_once(delta);
      rsq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
      r2inv = 1.0/rsq;
      rinv = sqrt(r2inv);

      eqm += qqrd2e * qqm[i]*qqm[j]*rinv;
      fpair = qqrd2e * qqm[i]*qqm[j]*rinv*r2inv;
      
      f[3*i+0] += delta[0]*fpair;
      f[3*i+1] += delta[1]*fpair;
      f[3*i+2] += delta[2]*fpair;
      f[3*j+0] -= delta[0]*fpair;
      f[3*j+1] -= delta[1]*fpair;
      f[3*j+2] -= delta[2]*fpair;
    }
  }
  
  return 0;
}
