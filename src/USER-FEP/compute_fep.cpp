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
   Contributing author: Agilio Padua (Univ Blaise Pascal & CNRS)
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <mpi.h>
#include "comm.h"
#include "update.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "input.h"
#include "fix.h"
#include "modify.h"
#include "variable.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "compute_fep.h"

using namespace LAMMPS_NS;

enum{PAIR,ATOM};
enum{CHARGE};

/* ---------------------------------------------------------------------- */

ComputeFEP::ComputeFEP(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal number of arguments in compute fep");

  scalar_flag = 0;
  vector_flag = 1;
  size_vector = 3;
  extvector = 0;

  vector = new double[3];

  fepinitflag = 0;    // avoid init to run entirely when called by write_data

  temp_fep = force->numeric(FLERR,arg[3]);

  // count # of perturbations

  npert = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+6 > narg) error->all(FLERR,
                                    "Illegal pair attribute in compute fep");
      npert++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+4 > narg) error->all(FLERR,
                                    "Illegal atom attribute in compute fep");
      npert++;
      iarg += 4;
    } else break;
  }

  if (npert == 0) error->all(FLERR,"Illegal syntax in compute fep");
  perturb = new Perturb[npert];

  // parse keywords

  npert = 0;
  chgflag = 0;

  iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      perturb[npert].which = PAIR;
      int n = strlen(arg[iarg+1]) + 1;
      perturb[npert].pstyle = new char[n];
      strcpy(perturb[npert].pstyle,arg[iarg+1]);
      n = strlen(arg[iarg+2]) + 1;
      perturb[npert].pparam = new char[n];
      strcpy(perturb[npert].pparam,arg[iarg+2]);
      force->bounds(FLERR,arg[iarg+3],atom->ntypes,
                    perturb[npert].ilo,perturb[npert].ihi);
      force->bounds(FLERR,arg[iarg+4],atom->ntypes,
                    perturb[npert].jlo,perturb[npert].jhi);
      if (strstr(arg[iarg+5],"v_") == arg[iarg+5]) {
        n = strlen(&arg[iarg+5][2]) + 1;
        perturb[npert].var = new char[n];
        strcpy(perturb[npert].var,&arg[iarg+5][2]);
      } else error->all(FLERR,"Illegal variable in compute fep");
      npert++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      perturb[npert].which = ATOM;
      if (strcmp(arg[iarg+1],"charge") == 0) {
        perturb[npert].aparam = CHARGE;
        chgflag = 1;
      } else error->all(FLERR,"Illegal atom argument in compute fep");
      force->bounds(FLERR,arg[iarg+2],atom->ntypes,
                    perturb[npert].ilo,perturb[npert].ihi);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        perturb[npert].var = new char[n];
        strcpy(perturb[npert].var,&arg[iarg+3][2]);
      } else error->all(FLERR,"Illegal variable in compute fep");
      npert++;
      iarg += 4;
    } else break;
  }

  // optional keywords

  tailflag = 0;
  volumeflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"tail") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal optional keyword "
                                    "in compute fep");
      if (strcmp(arg[iarg+1],"no") == 0) tailflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) tailflag = 1;
      else error->all(FLERR,"Illegal optional keyword in compute fep");
      iarg += 2;
    } else if (strcmp(arg[iarg],"volume") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal optional keyword "
                                    "in compute fep");
      if (strcmp(arg[iarg+1],"no") == 0) volumeflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) volumeflag = 1;
      else error->all(FLERR,"Illegal optional keyword in compute fep");
      iarg += 2;
    } else
      error->all(FLERR,"Illegal optional keyword in compute fep");
  }

  // allocate pair style arrays

  int ntype = atom->ntypes;
  for (int m = 0; m < npert; m++) {
    if (perturb[m].which == PAIR)
      memory->create(perturb[m].array_orig,ntype+1,ntype+1,"fep:array_orig");
  }

  // allocate space for charge, force, energy, virial arrays

  f_orig = NULL;
  q_orig = NULL;
  peatom_orig = keatom_orig = NULL;
  pvatom_orig = kvatom_orig = NULL;

  allocate_storage();

  fixgpu = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeFEP::~ComputeFEP()
{
  delete [] vector;

  for (int m = 0; m < npert; m++) {
    delete [] perturb[m].var;
    if (perturb[m].which == PAIR) {
      delete [] perturb[m].pstyle;
      delete [] perturb[m].pparam;
      memory->destroy(perturb[m].array_orig);
    }
  }
  delete [] perturb;

  deallocate_storage();
}

/* ---------------------------------------------------------------------- */

void ComputeFEP::init()
{
  int i,j;

  if (!fepinitflag)    // avoid init to run entirely when called by write_data
      fepinitflag = 1;
  else return;

  // setup and error checks

  pairflag = 0;

  for (int m = 0; m < npert; m++) {
    Perturb *pert = &perturb[m];

    pert->ivar = input->variable->find(pert->var);
    if (pert->ivar < 0)
      error->all(FLERR,"Variable name for compute fep does not exist");
    if (!input->variable->equalstyle(pert->ivar))
      error->all(FLERR,"Variable for compute fep is of invalid style");

    if (force->pair == NULL)
      error->all(FLERR,"compute fep pair requires pair interactions");

    if (pert->which == PAIR) {
      pairflag = 1;

      Pair *pair = force->pair_match(pert->pstyle,1);
      if (pair == NULL) error->all(FLERR,"compute fep pair style "
                                   "does not exist");
      void *ptr = pair->extract(pert->pparam,pert->pdim);
      if (ptr == NULL)
        error->all(FLERR,"compute fep pair style param not supported");

      pert->array = (double **) ptr;

      // if pair hybrid, test that ilo,ihi,jlo,jhi are valid for sub-style

      if ((strcmp(force->pair_style,"hybrid") == 0 ||
           strcmp(force->pair_style,"hybrid/overlay") == 0)) {
        PairHybrid *pair = (PairHybrid *) force->pair;
        for (i = pert->ilo; i <= pert->ihi; i++)
          for (j = MAX(pert->jlo,i); j <= pert->jhi; j++)
            if (!pair->check_ijtype(i,j,pert->pstyle))
              error->all(FLERR,"compute fep type pair range is not valid for "
                         "pair hybrid sub-style");
      }

    } else if (pert->which == ATOM) {
      if (pert->aparam == CHARGE) {
        if (!atom->q_flag)
          error->all(FLERR,"compute fep requires atom attribute charge");
      }
    }
  }

  if (tailflag) {
    if (force->pair->tail_flag == 0)
      error->all(FLERR,"Compute fep tail when pair style does not "
                 "compute tail corrections");
  }

  // detect if package gpu is present

  int ifixgpu = modify->find_fix("package_gpu");
  if (ifixgpu >= 0) fixgpu = modify->fix[ifixgpu];

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen, "FEP settings ...\n");
      fprintf(screen, "  temperature = %f\n", temp_fep);
      fprintf(screen, "  tail %s\n", (tailflag ? "yes":"no"));
      for (int m = 0; m < npert; m++) {
        Perturb *pert = &perturb[m];
        if (pert->which == PAIR)
          fprintf(screen, "  %s %s %d-%d %d-%d\n", pert->pstyle, pert->pparam,
                  pert->ilo, pert->ihi, pert->jlo, pert->jhi);
        else if (pert->which == ATOM)
          fprintf(screen, "  %d-%d charge\n", pert->ilo, pert->ihi);
      }
    }
    if (logfile) {
      fprintf(logfile, "FEP settings ...\n");
      fprintf(logfile, "  temperature = %f\n", temp_fep);
      fprintf(logfile, "  tail %s\n", (tailflag ? "yes":"no"));
      for (int m = 0; m < npert; m++) {
        Perturb *pert = &perturb[m];
        if (pert->which == PAIR)
          fprintf(logfile, "  %s %s %d-%d %d-%d\n", pert->pstyle, pert->pparam,
                  pert->ilo, pert->ihi, pert->jlo, pert->jhi);
        else if (pert->which == ATOM)
          fprintf(logfile, "  %d-%d charge\n", pert->ilo, pert->ihi);
      }
    }
  }

}

/* ---------------------------------------------------------------------- */


void ComputeFEP::compute_vector()
{
  double pe0,pe1;

  eflag = 1;
  vflag = 0;

  invoked_vector = update->ntimestep;

  if (atom->nmax > nmax) {  // reallocate working arrays if necessary
    deallocate_storage();
    allocate_storage();
  }

  backup_qfev();   // backup charge, force, energy, virial array values
  backup_params(); // backup pair parameters

  timer->stamp();
  if (force->pair && force->pair->compute_flag) {
    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);
  }
  if (chgflag && force->kspace && force->kspace->compute_flag) {
    force->kspace->compute(eflag,vflag);
    timer->stamp(Timer::KSPACE);
  }

  // accumulate force/energy/virial from /gpu pair styles
  if (fixgpu) fixgpu->post_force(vflag);

  pe0 = compute_epair();

  perturb_params();

  timer->stamp();
  if (force->pair && force->pair->compute_flag) {
    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);
  }
  if (chgflag && force->kspace && force->kspace->compute_flag) {
    force->kspace->compute(eflag,vflag);
    timer->stamp(Timer::KSPACE);
  }

  // accumulate force/energy/virial from /gpu pair styles
  // this is required as to empty the answer queue,
  // otherwise the force compute on the GPU in the next step would be incorrect
  if (fixgpu) fixgpu->post_force(vflag);

  pe1 = compute_epair();

  restore_qfev();   // restore charge, force, energy, virial array values
  restore_params(); // restore pair parameters

  vector[0] = pe1-pe0;
  vector[1] = exp(-(pe1-pe0)/(force->boltz*temp_fep));
  vector[2] = domain->xprd * domain->yprd * domain->zprd;
  if (volumeflag)
    vector[1] *= vector[2];
}


/* ----------------------------------------------------------------------
   obtain pair energy from lammps accumulators
------------------------------------------------------------------------- */

double ComputeFEP::compute_epair()
{
  double eng, eng_pair;

  eng = 0.0;
  if (force->pair)
    eng = force->pair->eng_vdwl + force->pair->eng_coul;
  MPI_Allreduce(&eng,&eng_pair,1,MPI_DOUBLE,MPI_SUM,world);

  if (tailflag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    eng_pair += force->pair->etail / volume;
  }

  if (chgflag && force->kspace) eng_pair += force->kspace->energy;

  return eng_pair;
}


/* ----------------------------------------------------------------------
   apply perturbation to pair, atom parameters based on variable evaluation
------------------------------------------------------------------------- */

void ComputeFEP::perturb_params()
{
  int i,j;

  for (int m = 0; m < npert; m++) {
    Perturb *pert = &perturb[m];

    double delta = input->variable->compute_equal(pert->ivar);

    if (pert->which == PAIR) {      // modify pair parameters
      for (i = pert->ilo; i <= pert->ihi; i++)
        for (j = MAX(pert->jlo,i); j <= pert->jhi; j++)
          pert->array[i][j] = pert->array_orig[i][j] + delta;

    } else if (pert->which == ATOM) {

      if (pert->aparam == CHARGE) {      // modify charges
        int *atype = atom->type;
        double *q = atom->q;
        int *mask = atom->mask;
        int natom = atom->nlocal + atom->nghost;

        for (i = 0; i < natom; i++)
          if (atype[i] >= pert->ilo && atype[i] <= pert->ihi)
            if (mask[i] & groupbit)
              q[i] += delta;

      }
    }
  }

  // re-initialize pair styles if any PAIR settings were changed
  // this resets other coeffs that may depend on changed values,
  // and also offset and tail corrections

  if (pairflag) force->pair->reinit();

  // reset KSpace charges if charges have changed

  if (chgflag && force->kspace) force->kspace->qsum_qsq();
}


/* ----------------------------------------------------------------------
   backup pair parameters
------------------------------------------------------------------------- */

void ComputeFEP::backup_params()
{
  int i,j;

  for (int m = 0; m < npert; m++) {
    Perturb *pert = &perturb[m];
    if (pert->which == PAIR) {
      for (i = pert->ilo; i <= pert->ihi; i++)
        for (j = MAX(pert->jlo,i); j <= pert->jhi; j++)
          pert->array_orig[i][j] = pert->array[i][j];
    }
  }
}


/* ----------------------------------------------------------------------
   restore pair parameters to original values
------------------------------------------------------------------------- */

void ComputeFEP::restore_params()
{
  int i,j;

  for (int m = 0; m < npert; m++) {
    Perturb *pert = &perturb[m];
    if (pert->which == PAIR) {
      for (i = pert->ilo; i <= pert->ihi; i++)
        for (j = MAX(pert->jlo,i); j <= pert->jhi; j++)
          pert->array[i][j] = pert->array_orig[i][j];
    }
  }

  if (pairflag) force->pair->reinit();

  // reset KSpace charges if charges have changed

  if (chgflag && force->kspace) force->kspace->qsum_qsq();
}


/* ----------------------------------------------------------------------
   manage storage for charge, force, energy, virial arrays
------------------------------------------------------------------------- */

void ComputeFEP::allocate_storage()
{
  nmax = atom->nmax;
  memory->create(f_orig,nmax,3,"fep:f_orig");
  memory->create(peatom_orig,nmax,"fep:peatom_orig");
  memory->create(pvatom_orig,nmax,6,"fep:pvatom_orig");
  if (chgflag) {
    memory->create(q_orig,nmax,"fep:q_orig");
    if (force->kspace) {
      memory->create(keatom_orig,nmax,"fep:keatom_orig");
      memory->create(kvatom_orig,nmax,6,"fep:kvatom_orig");
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeFEP::deallocate_storage()
{
  memory->destroy(f_orig);
  memory->destroy(peatom_orig);
  memory->destroy(pvatom_orig);
  memory->destroy(q_orig);
  memory->destroy(keatom_orig);
  memory->destroy(kvatom_orig);

  f_orig = NULL;
  q_orig = NULL;
  peatom_orig = keatom_orig = NULL;
  pvatom_orig = kvatom_orig = NULL;
}


/* ----------------------------------------------------------------------
   backup and restore arrays with charge, force, energy, virial
------------------------------------------------------------------------- */

void ComputeFEP::backup_qfev()
{
  int i;

  int nall = atom->nlocal + atom->nghost;
  int natom = atom->nlocal;
  if (force->newton || force->kspace->tip4pflag)
      natom += atom->nghost;

  double **f = atom->f;
  for (i = 0; i < natom; i++) {
    f_orig[i][0] = f[i][0];
    f_orig[i][1] = f[i][1];
    f_orig[i][2] = f[i][2];
  }

  eng_vdwl_orig = force->pair->eng_vdwl;
  eng_coul_orig = force->pair->eng_coul;

  pvirial_orig[0] = force->pair->virial[0];
  pvirial_orig[1] = force->pair->virial[1];
  pvirial_orig[2] = force->pair->virial[2];
  pvirial_orig[3] = force->pair->virial[3];
  pvirial_orig[4] = force->pair->virial[4];
  pvirial_orig[5] = force->pair->virial[5];

  if (update->eflag_atom) {
    double *peatom = force->pair->eatom;
    for (i = 0; i < natom; i++)
      peatom_orig[i] = peatom[i];
  }
  if (update->vflag_atom) {
    double **pvatom = force->pair->vatom;
    for (i = 0; i < natom; i++) {
      pvatom_orig[i][0] = pvatom[i][0];
      pvatom_orig[i][1] = pvatom[i][1];
      pvatom_orig[i][2] = pvatom[i][2];
      pvatom_orig[i][3] = pvatom[i][3];
      pvatom_orig[i][4] = pvatom[i][4];
      pvatom_orig[i][5] = pvatom[i][5];
    }
  }

  if (chgflag) {
    double *q = atom->q;
    for (i = 0; i < nall; i++)
      q_orig[i] = q[i];

    if (force->kspace) {
      energy_orig = force->kspace->energy;
      kvirial_orig[0] = force->kspace->virial[0];
      kvirial_orig[1] = force->kspace->virial[1];
      kvirial_orig[2] = force->kspace->virial[2];
      kvirial_orig[3] = force->kspace->virial[3];
      kvirial_orig[4] = force->kspace->virial[4];
      kvirial_orig[5] = force->kspace->virial[5];

      if (update->eflag_atom) {
        double *keatom = force->kspace->eatom;
        for (i = 0; i < natom; i++)
          keatom_orig[i] = keatom[i];
      }
      if (update->vflag_atom) {
        double **kvatom = force->kspace->vatom;
        for (i = 0; i < natom; i++) {
          kvatom_orig[i][0] = kvatom[i][0];
          kvatom_orig[i][1] = kvatom[i][1];
          kvatom_orig[i][2] = kvatom[i][2];
          kvatom_orig[i][3] = kvatom[i][3];
          kvatom_orig[i][4] = kvatom[i][4];
          kvatom_orig[i][5] = kvatom[i][5];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeFEP::restore_qfev()
{
  int i;

  int nall = atom->nlocal + atom->nghost;
  int natom = atom->nlocal;
  if (force->newton || force->kspace->tip4pflag)
      natom += atom->nghost;

  double **f = atom->f;
  for (i = 0; i < natom; i++) {
    f[i][0] = f_orig[i][0];
    f[i][1] = f_orig[i][1];
    f[i][2] = f_orig[i][2];
  }

  force->pair->eng_vdwl = eng_vdwl_orig;
  force->pair->eng_coul = eng_coul_orig;

  force->pair->virial[0] = pvirial_orig[0];
  force->pair->virial[1] = pvirial_orig[1];
  force->pair->virial[2] = pvirial_orig[2];
  force->pair->virial[3] = pvirial_orig[3];
  force->pair->virial[4] = pvirial_orig[4];
  force->pair->virial[5] = pvirial_orig[5];

  if (update->eflag_atom) {
    double *peatom = force->pair->eatom;
    for (i = 0; i < natom; i++)
      peatom[i] = peatom_orig[i];
  }
  if (update->vflag_atom) {
    double **pvatom = force->pair->vatom;
    for (i = 0; i < natom; i++) {
      pvatom[i][0] = pvatom_orig[i][0];
      pvatom[i][1] = pvatom_orig[i][1];
      pvatom[i][2] = pvatom_orig[i][2];
      pvatom[i][3] = pvatom_orig[i][3];
      pvatom[i][4] = pvatom_orig[i][4];
      pvatom[i][5] = pvatom_orig[i][5];
    }
  }

  if (chgflag) {
    double *q = atom->q;
    for (i = 0; i < nall; i++)
      q[i] = q_orig[i];

    if (force->kspace) {
      force->kspace->energy = energy_orig;
      force->kspace->virial[0] = kvirial_orig[0];
      force->kspace->virial[1] = kvirial_orig[1];
      force->kspace->virial[2] = kvirial_orig[2];
      force->kspace->virial[3] = kvirial_orig[3];
      force->kspace->virial[4] = kvirial_orig[4];
      force->kspace->virial[5] = kvirial_orig[5];

      if (update->eflag_atom) {
        double *keatom = force->kspace->eatom;
        for (i = 0; i < natom; i++)
          keatom[i] = keatom_orig[i];
      }
      if (update->vflag_atom) {
        double **kvatom = force->kspace->vatom;
        for (i = 0; i < natom; i++) {
          kvatom[i][0] = kvatom_orig[i][0];
          kvatom[i][1] = kvatom_orig[i][1];
          kvatom[i][2] = kvatom_orig[i][2];
          kvatom[i][3] = kvatom_orig[i][3];
          kvatom[i][4] = kvatom_orig[i][4];
          kvatom[i][5] = kvatom_orig[i][5];
        }
      }
    }
  }
}

