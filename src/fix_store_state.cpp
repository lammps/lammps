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

#include "stdlib.h"
#include "string.h"
#include "fix_store_state.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{KEYWORD,COMPUTE,FIX,VARIABLE};

#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

FixStoreState::FixStoreState(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix store/state command");

  restart_peratom = 1;
  peratom_freq = 1;

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix store/state command");

  // parse values until one isn't recognized
  // customize a new keyword by adding to if statement

  pack_choice = new FnPtrPack[narg-4];
  which = new int[narg-4];
  argindex = new int[narg-4];
  ids = new char*[narg-4];
  value2index = new int[narg-4];
  nvalues = 0;
  cfv_any = 0;

  int iarg = 4;
  while (iarg < narg) {
    which[nvalues] = KEYWORD;
    ids[nvalues] = NULL;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_id;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (!atom->molecule_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_molecule;
    } else if (strcmp(arg[iarg],"type") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_type;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_mass;

    } else if (strcmp(arg[iarg],"x") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_x;
    } else if (strcmp(arg[iarg],"y") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_y;
    } else if (strcmp(arg[iarg],"z") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_z;
    } else if (strcmp(arg[iarg],"xs") == 0) {
      if (domain->triclinic)
        pack_choice[nvalues++] = &FixStoreState::pack_xs_triclinic;
      else pack_choice[nvalues++] = &FixStoreState::pack_xs;
    } else if (strcmp(arg[iarg],"ys") == 0) {
      if (domain->triclinic)
        pack_choice[nvalues++] = &FixStoreState::pack_ys_triclinic;
      else pack_choice[nvalues++] = &FixStoreState::pack_ys;
    } else if (strcmp(arg[iarg],"zs") == 0) {
      if (domain->triclinic)
        pack_choice[nvalues++] = &FixStoreState::pack_zs_triclinic;
      else pack_choice[nvalues++] = &FixStoreState::pack_zs;
    } else if (strcmp(arg[iarg],"xu") == 0) {
      if (domain->triclinic)
        pack_choice[nvalues++] = &FixStoreState::pack_xu_triclinic;
      else pack_choice[nvalues++] = &FixStoreState::pack_xu;
    } else if (strcmp(arg[iarg],"yu") == 0) {
      if (domain->triclinic)
        pack_choice[nvalues++] = &FixStoreState::pack_yu_triclinic;
      else pack_choice[nvalues++] = &FixStoreState::pack_yu;
    } else if (strcmp(arg[iarg],"zu") == 0) {
      if (domain->triclinic)
        pack_choice[nvalues++] = &FixStoreState::pack_zu_triclinic;
      else pack_choice[nvalues++] = &FixStoreState::pack_zu;
    } else if (strcmp(arg[iarg],"ix") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_ix;
    } else if (strcmp(arg[iarg],"iy") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_iy;
    } else if (strcmp(arg[iarg],"iz") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_iz;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_vx;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_vy;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_vz;
    } else if (strcmp(arg[iarg],"fx") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_fx;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_fy;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      pack_choice[nvalues++] = &FixStoreState::pack_fz;

    } else if (strcmp(arg[iarg],"q") == 0) {
      if (!atom->q_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_q;
    } else if (strcmp(arg[iarg],"mux") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_mux;
    } else if (strcmp(arg[iarg],"muy") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_muy;
    } else if (strcmp(arg[iarg],"muz") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_muz;

    } else if (strcmp(arg[iarg],"radius") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_radius;
    } else if (strcmp(arg[iarg],"omegax") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_omegax;
    } else if (strcmp(arg[iarg],"omegay") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_omegay;
    } else if (strcmp(arg[iarg],"omegaz") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_omegaz;
    } else if (strcmp(arg[iarg],"angmomx") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_angmomx;
    } else if (strcmp(arg[iarg],"angmomy") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_angmomy;
    } else if (strcmp(arg[iarg],"angmomz") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_angmomz;
    } else if (strcmp(arg[iarg],"tqx") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_tqx;
    } else if (strcmp(arg[iarg],"tqy") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_tqy;
    } else if (strcmp(arg[iarg],"tqz") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,
                   "Fix store/state for atom property that isn't allocated");
      pack_choice[nvalues++] = &FixStoreState::pack_tqz;

    } else if (strncmp(arg[iarg],"c_",2) == 0 ||
               strncmp(arg[iarg],"f_",2) == 0 ||
               strncmp(arg[iarg],"v_",2) == 0) {
      cfv_any = 1;
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if (arg[iarg][0] == 'v') which[nvalues] = VARIABLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal fix store/state command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else break;

    iarg++;
  }

  // optional args

  comflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"com") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix store/state command");
      if (strcmp(arg[iarg+1],"no") == 0) comflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) comflag = 1;
      else error->all(FLERR,"Illegal fix store/state command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix store/state command");
  }

  // error check

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix store/state does not exist");
      if (modify->compute[icompute]->peratom_flag == 0)
        error->all(FLERR,"Fix store/state compute "
                   "does not calculate per-atom values");
      if (argindex[i] == 0 &&
          modify->compute[icompute]->size_peratom_cols != 0)
        error->all(FLERR,"Fix store/state compute does not "
                   "calculate a per-atom vector");
      if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
        error->all(FLERR,
                   "Fix store/state compute does not "
                   "calculate a per-atom array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_peratom_cols)
        error->all(FLERR,
                   "Fix store/state compute array is accessed out-of-range");

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,
                   "Fix ID for fix store/state does not exist");
      if (modify->fix[ifix]->peratom_flag == 0)
        error->all(FLERR,
                   "Fix store/state fix does not calculate per-atom values");
      if (argindex[i] == 0 && modify->fix[ifix]->size_peratom_cols != 0)
        error->all(FLERR,
                   "Fix store/state fix does not calculate a per-atom vector");
      if (argindex[i] && modify->fix[ifix]->size_peratom_cols == 0)
        error->all(FLERR,
                   "Fix store/state fix does not calculate a per-atom array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_peratom_cols)
        error->all(FLERR,
                   "Fix store/state fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->peratom_freq)
        error->all(FLERR,
                   "Fix for fix store/state not computed at compatible time");

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix store/state does not exist");
      if (input->variable->atomstyle(ivariable) == 0)
        error->all(FLERR,"Fix store/state variable is not atom-style variable");
    }
  }

  // this fix produces either a per-atom vector or array

  peratom_flag = 1;
  if (nvalues == 1) size_peratom_cols = 0;
  else size_peratom_cols = nvalues;

  // perform initial allocation of atom-based array
  // register with Atom class

  values = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // zero the array since dump may access it on timestep 0
  // zero the array since a variable may access it before first run

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    for (int m = 0; m < nvalues; m++)
      values[i][m] = 0.0;

  // store current values for keywords but not for compute, fix, variable

  kflag = 1;
  cfv_flag = 0;
  end_of_step();
  firstflag = 1;
}

/* ---------------------------------------------------------------------- */

FixStoreState::~FixStoreState()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  delete [] which;
  delete [] argindex;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;
  delete [] pack_choice;

  memory->destroy(values);
}

/* ---------------------------------------------------------------------- */

int FixStoreState::setmask()
{
  int mask = 0;
  if (nevery) mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStoreState::init()
{
  // set indices and check validity of all computes,fixes,variables
  // no error check if end_of_step() will not be called

  if (!firstflag && nevery == 0) return;

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix store/state does not exist");
      value2index[m] = icompute;

    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix store/state does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix store/state does not exist");
      value2index[m] = ivariable;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::setup(int vflag)
{
  // if first invocation, store current values for compute, fix, variable

  if (firstflag) {
    kflag = 0;
    cfv_flag = 1;
    end_of_step();
    firstflag = 0;
    kflag = cfv_flag = 1;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::end_of_step()
{
  int i,j,n;

  // compute com if comflag set

  if (comflag) {
    double masstotal = group->mass(igroup);
    group->xcm(igroup,masstotal,cm);
  }

  // if any compute/fix/variable and nevery, wrap with clear/add

  if (cfv_any && nevery) modify->clearstep_compute();

  // fill vector or array with per-atom values

  if (values) vbuf = &values[0][0];
  else vbuf = NULL;

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == KEYWORD && kflag) (this->*pack_choice[m])(m);

    else if (cfv_flag) {
      n = value2index[m];
      j = argindex[m];

      int *mask = atom->mask;
      int nlocal = atom->nlocal;

      // invoke compute if not previously invoked

      if (which[m] == COMPUTE) {
        Compute *compute = modify->compute[n];
        if (!(compute->invoked_flag & INVOKED_PERATOM)) {
          compute->compute_peratom();
          compute->invoked_flag |= INVOKED_PERATOM;
        }

        if (j == 0) {
          double *compute_vector = compute->vector_atom;
          for (i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) values[i][m] = compute_vector[i];
        } else {
          int jm1 = j - 1;
          double **compute_array = compute->array_atom;
          for (i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) values[i][m] = compute_array[i][jm1];
        }

      // access fix fields, guaranteed to be ready

      } else if (which[m] == FIX) {
        if (j == 0) {
          double *fix_vector = modify->fix[n]->vector_atom;
          for (i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) values[i][m] = fix_vector[i];
        } else {
          int jm1 = j - 1;
          double **fix_array = modify->fix[n]->array_atom;
          for (i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) values[i][m] = fix_array[i][jm1];
        }

      // evaluate atom-style variable

      } else if (which[m] == VARIABLE) {
        input->variable->compute_atom(n,igroup,&values[0][m],nvalues,0);
      }
    }
  }

  // if any compute/fix/variable and nevery, wrap with clear/add

  if (cfv_any && nevery) {
    int nextstep = (update->ntimestep/nevery)*nevery + nevery;
    modify->addstep_compute(nextstep);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixStoreState::memory_usage()
{
  double bytes = atom->nmax*nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixStoreState::grow_arrays(int nmax)
{
  memory->grow(values,nmax,nvalues,"store/state:values");
  if (nvalues == 1) {
    if (nmax) vector_atom = &values[0][0];
    else vector_atom = NULL;
  } else array_atom = values;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixStoreState::copy_arrays(int i, int j, int delflag)
{
  for (int m = 0; m < nvalues; m++) values[j][m] = values[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixStoreState::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nvalues; m++) buf[m] = values[i][m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixStoreState::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nvalues; m++) values[nlocal][m] = buf[m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixStoreState::pack_restart(int i, double *buf)
{
  buf[0] = nvalues+1;
  for (int m = 0; m < nvalues; m++) buf[m+1] = values[i][m];
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixStoreState::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  for (int i = 0; i < nvalues; i++) values[nlocal][i] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixStoreState::maxsize_restart()
{
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixStoreState::size_restart(int nlocal)
{
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   one method for every keyword fix store/state can archive
   the atom property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_id(int n)
{
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = tag[i];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_molecule(int n)
{
  tagint *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = molecule[i];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_type(int n)
{
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = type[i];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_mass(int n)
{
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) vbuf[n] = rmass[i];
      else vbuf[n] = 0.0;
      n += nvalues;
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) vbuf[n] = mass[type[i]];
      else vbuf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_x(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = x[i][0];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_y(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = x[i][1];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_z(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = x[i][2];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_xs(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double boxxlo = domain->boxlo[0];
  double invxprd = 1.0/domain->xprd;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = (x[i][0] - boxxlo) * invxprd;
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_ys(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double boxylo = domain->boxlo[1];
  double invyprd = 1.0/domain->yprd;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = (x[i][1] - boxylo) * invyprd;
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_zs(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double boxzlo = domain->boxlo[2];
  double invzprd = 1.0/domain->zprd;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = (x[i][2] - boxzlo) * invzprd;
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_xs_triclinic(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      vbuf[n] = h_inv[0]*(x[i][0]-boxlo[0]) +
        h_inv[5]*(x[i][1]-boxlo[1]) + h_inv[4]*(x[i][2]-boxlo[2]);
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_ys_triclinic(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      vbuf[n] = h_inv[1]*(x[i][1]-boxlo[1]) + h_inv[3]*(x[i][2]-boxlo[2]);
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_zs_triclinic(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      vbuf[n] = h_inv[2]*(x[i][2]-boxlo[2]);
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_xu(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vbuf[n] = x[i][0] + ((image[i] & IMGMASK) - IMGMAX) * xprd;
      if (comflag) vbuf[n] -= cm[0];
    } else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_yu(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double yprd = domain->yprd;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vbuf[n] = x[i][1] + ((image[i] >> IMGBITS & IMGMASK) - IMGMAX) * yprd;
      if (comflag) vbuf[n] -= cm[1];
    } else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_zu(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double zprd = domain->zprd;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vbuf[n] = x[i][2] + ((image[i] >> IMG2BITS) - IMGMAX) * zprd;
      if (comflag) vbuf[n] -= cm[2];
    } else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_xu_triclinic(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  int xbox,ybox,zbox;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xbox = (image[i] & IMGMASK) - IMGMAX;
      ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (image[i] >> IMG2BITS) - IMGMAX;
      vbuf[n] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
      if (comflag) vbuf[n] -= cm[0];
    } else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_yu_triclinic(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  int ybox,zbox;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (image[i] >> IMG2BITS) - IMGMAX;
      vbuf[n] = x[i][1] + h[1]*ybox + h[3]*zbox;
      if (comflag) vbuf[n] -= cm[1];
    } else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_zu_triclinic(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  int zbox;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      zbox = (image[i] >> IMG2BITS) - IMGMAX;
      vbuf[n] = x[i][2] + h[2]*zbox;
      if (comflag) vbuf[n] -= cm[2];
    } else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_ix(int n)
{
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = (image[i] & IMGMASK) - IMGMAX;
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_iy(int n)
{
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_iz(int n)
{
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = (image[i] >> IMG2BITS) - IMGMAX;
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_vx(int n)
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = v[i][0];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_vy(int n)
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = v[i][1];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_vz(int n)
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = v[i][2];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_fx(int n)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = f[i][0];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_fy(int n)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = f[i][1];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_fz(int n)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = f[i][2];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_q(int n)
{
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = q[i];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_mux(int n)
{
  double **mu = atom->mu;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = mu[i][0];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_muy(int n)
{
  double **mu = atom->mu;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = mu[i][1];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_muz(int n)
{
  double **mu = atom->mu;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = mu[i][2];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_radius(int n)
{
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = radius[i];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_omegax(int n)
{
  double **omega = atom->omega;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = omega[i][0];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_omegay(int n)
{
  double **omega = atom->omega;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = omega[i][1];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_omegaz(int n)
{
  double **omega = atom->omega;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = omega[i][2];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_angmomx(int n)
{
  double **angmom = atom->angmom;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = angmom[i][0];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_angmomy(int n)
{
  double **angmom = atom->angmom;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = angmom[i][1];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_angmomz(int n)
{
  double **angmom = atom->angmom;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = angmom[i][2];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_tqx(int n)
{
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = torque[i][0];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_tqy(int n)
{
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = torque[i][1];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_tqz(int n)
{
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = torque[i][2];
    else vbuf[n] = 0.0;
    n += nvalues;
  }
}
