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
#include "group.h"
#include "domain.h"
#include "neighbor.h"
#include "force.h"
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
extern char *util_date();
extern void seconds(double *);
extern int pspw_geovib(MPI_Comm, std::string&);
extern int pspw_minimizer(MPI_Comm, double *, double *, double *, double *);
extern int nwchem_abiversion();
}

// the ABIVERSION number here must be kept consistent
// with its counterpart in the NWChem library and the
// prototype above. We want to catch mismatches with
// a meaningful error messages, as they can cause
// difficult to debug crashes or memory corruption.

#define NWCHEM_ABIVERSION 20180622

// NOTE: change 1-proc restriction later

/* ---------------------------------------------------------------------- */

FixNWChem::FixNWChem(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Must use units metal with fix nwchem command");

  if (comm->nprocs != 1)
    error->all(FLERR,"Fix nwchem currently runs only in serial");

  if (NWCHEM_ABIVERSION != pwdft::nwchem_abiversion())
    error->all(FLERR,"LAMMPS is linked against incompatible NWChem library");

  if (narg != 4) error->all(FLERR,"Illegal fix nwchem command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  thermo_energy = 1;

  // store ID of compute pe/atom used to trigger per-atom energy computations

  id_pe = utils::strdup(arg[3]);
  int ipe = modify->find_compute(id_pe);
  if (ipe < 0) error->all(FLERR,"Could not find fix nwchem compute ID");
  if (modify->compute[ipe]->peatomflag == 0)
    error->all(FLERR,"Fix nwchem compute ID does not compute pe/atom");

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
  delete [] id_pe;
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
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNWChem::init()
{
  // error checks

  if (domain->dimension == 2)
    error->all(FLERR,"Fix nwchem requires 3d problem");

  if (domain->nonperiodic == 0) pbcflag = 1;
  else if (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic)
    pbcflag = 0;
  else error->all(FLERR,"Fix nwchem requires fully periodic or "
                  "fully non-periodic system");
  
  if (atom->q_flag == 0 || force->pair == nullptr || force->kspace == nullptr)
    error->all(FLERR,"Fix nwchem cannot compute Coulomb potential");

  // c_pe = instance of compute pe/atom

  int ipe = modify->find_compute(id_pe);
  if (ipe < 0) error->all(FLERR,"Could not find fix nwchem compute ID");
  c_pe = modify->compute[ipe];

  // fix group = QM atoms
  // one-time initialization of qmIDs
  // NOTE: need to sort in ascending order to match NWChem ?
  // NOTE: make nqm an int, check that group count is not 0 or > MAXBIGINT
  
  if (nqm == 0) {
    nqm = group->count(igroup);
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
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) qmIDs[nqm++] = tag[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixNWChem::setup(int vflag)
{
  //newsystem = 1;
  post_force(vflag);
  //newsystem = 0;
}

/* ---------------------------------------------------------------------- */

void FixNWChem::min_setup(int vflag)
{
  //newsystem = 1;
  post_force(vflag);
  //newsystem = 0;
}

/* ---------------------------------------------------------------------- */

void FixNWChem::post_force(int vflag)
{
  int ilocal,jlocal;
  double delx,dely,delz,rsq;

  // on reneigh step, reset qm2lmp vector for indices of QM atoms
  // NOTE: if in parallel, do this differently

  if (neighbor->ago == 0)
    for (int i = 0; i < nqm; i++)
      qm2lmp[i] = atom->map(qmIDs[i]);

  // create 2 NWChem inputs: xqm and qpotential (Coulomb potential)
  // qpotential[i] = (eatom[i] from pair_coul + kspace) / Qi
  // subtract off Qj/Rij for QM I interacting with all other QM J atoms
 
  double **x = atom->x;
  double *q = atom->q;
  double *eatom_pair = pair_coul->eatom;
  double *eatom_kspace = force->kspace->eatom;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2lmp[i];

    // NOTE: what if LAMMPS atom moves via PBC to other side of LAMMPS box ?
    xqm[i][0] = x[ilocal][0];
    xqm[i][1] = x[ilocal][1];
    xqm[i][2] = x[ilocal][2];

    if (q[i] == 0.0) qpotential[i] = 0.0;
    else qpotential[i] = (eatom_pair[ilocal] + eatom_kspace[ilocal]) / q[ilocal];

    for (int j = 0; j < nqm; j++) {
      if (j == i) continue;
      jlocal = qm2lmp[j];
      // NOTE: apply PBC
      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      qpotential[i] -= q[j] / sqrt(rsq);
    }
  }

  // call to NWChem
  // input/outputs are only for QM atoms
  // inputs = xqm and qpotential
  // outputs = fqm and qqm

  int nwerr = pwdft::pspw_minimizer(world,&xqm[0][0],qpotential,&fqm[0][0],qqm);
  if (nwerr) error->all(FLERR,"Internal NWChem error");

  //latte(flags,&natoms,coords,type,&ntypes,mass,boxlo,boxhi,&domain->xy,
  //      &domain->xz,&domain->yz,forces,&maxiter,&latte_energy,
  //      &atom->v[0][0],&update->dt,virial,&newsystem,&latteerror);

  // reset Q of QM atoms

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2lmp[i];
    q[ilocal] = qqm[i];
  }

  // calculate eatom[i] again with new Q values via pair_coul and KSpace
  // NOTE: how to zero f array w/out wiping out forces I need to keep
  // NOTE: need to comm ghost forces (serial or parallel)
  // NOTE: how to have final f array = total force, including LJ
  
  //fhold = f - pair_coul(f) - kspace(f);
  //f = 0;

  pair_coul->compute(2,0);
  force->kspace->compute(2,0);

  //f += fhold;

  // add NWChem QM forces to QM atoms

  double **f = atom->f;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2lmp[i];
    f[ilocal][0] += fqm[i][0];
    f[ilocal][1] += fqm[i][1];
    f[ilocal][2] += fqm[i][2];
  }

  // trigger per-atom energy computation on next timestep by pair/kspace

  c_pe->addstep(update->ntimestep+1);
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
