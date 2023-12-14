// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_store_state.h"

#include "arg_info.h"
#include "atom.h"
#include "compute.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixStoreState::FixStoreState(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), avalues(nullptr), vbuf(nullptr)
{
  if (narg < 5) utils::missing_cmd_args(FLERR,"fix store/state", error);

  restart_peratom = 1;
  peratom_freq = 1;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 0) error->all(FLERR,"Invalid fix store/state never value {}", nevery);

  // parse values
  // customize a new keyword by adding to if statement

  values.clear();
  cfv_any = 0;

  int iarg = 4;
  while (iarg < narg) {

    value_t val;
    val.which = ArgInfo::KEYWORD;
    val.argindex = -1;
    val.id = "";
    val.val.c = nullptr;
    val.pack_choice = nullptr;

    if (strcmp(arg[iarg],"id") == 0) {
      val.pack_choice = &FixStoreState::pack_id;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (!atom->molecule_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_molecule;
    } else if (strcmp(arg[iarg],"type") == 0) {
      val.pack_choice = &FixStoreState::pack_type;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      val.pack_choice = &FixStoreState::pack_mass;

    } else if (strcmp(arg[iarg],"x") == 0) {
      val.pack_choice = &FixStoreState::pack_x;
    } else if (strcmp(arg[iarg],"y") == 0) {
      val.pack_choice = &FixStoreState::pack_y;
    } else if (strcmp(arg[iarg],"z") == 0) {
      val.pack_choice = &FixStoreState::pack_z;
    } else if (strcmp(arg[iarg],"xs") == 0) {
      if (domain->triclinic)
        val.pack_choice = &FixStoreState::pack_xs_triclinic;
      else val.pack_choice = &FixStoreState::pack_xs;
    } else if (strcmp(arg[iarg],"ys") == 0) {
      if (domain->triclinic)
        val.pack_choice = &FixStoreState::pack_ys_triclinic;
      else val.pack_choice = &FixStoreState::pack_ys;
    } else if (strcmp(arg[iarg],"zs") == 0) {
      if (domain->triclinic)
        val.pack_choice = &FixStoreState::pack_zs_triclinic;
      else val.pack_choice = &FixStoreState::pack_zs;
    } else if (strcmp(arg[iarg],"xu") == 0) {
      if (domain->triclinic)
        val.pack_choice = &FixStoreState::pack_xu_triclinic;
      else val.pack_choice = &FixStoreState::pack_xu;
    } else if (strcmp(arg[iarg],"yu") == 0) {
      if (domain->triclinic)
        val.pack_choice = &FixStoreState::pack_yu_triclinic;
      else val.pack_choice = &FixStoreState::pack_yu;
    } else if (strcmp(arg[iarg],"zu") == 0) {
      if (domain->triclinic)
        val.pack_choice = &FixStoreState::pack_zu_triclinic;
      else val.pack_choice = &FixStoreState::pack_zu;
    } else if (strcmp(arg[iarg],"xsu") == 0) {
      if (domain->triclinic)
        val.pack_choice = &FixStoreState::pack_xsu_triclinic;
      else val.pack_choice = &FixStoreState::pack_xsu;
    } else if (strcmp(arg[iarg],"ysu") == 0) {
      if (domain->triclinic)
        val.pack_choice = &FixStoreState::pack_ysu_triclinic;
      else val.pack_choice = &FixStoreState::pack_ysu;
    } else if (strcmp(arg[iarg],"zsu") == 0) {
      if (domain->triclinic)
        val.pack_choice = &FixStoreState::pack_zsu_triclinic;
      else val.pack_choice = &FixStoreState::pack_zsu;

    } else if (strcmp(arg[iarg],"ix") == 0) {
      val.pack_choice = &FixStoreState::pack_ix;
    } else if (strcmp(arg[iarg],"iy") == 0) {
      val.pack_choice = &FixStoreState::pack_iy;
    } else if (strcmp(arg[iarg],"iz") == 0) {
      val.pack_choice = &FixStoreState::pack_iz;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      val.pack_choice = &FixStoreState::pack_vx;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      val.pack_choice = &FixStoreState::pack_vy;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      val.pack_choice = &FixStoreState::pack_vz;
    } else if (strcmp(arg[iarg],"fx") == 0) {
      val.pack_choice = &FixStoreState::pack_fx;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      val.pack_choice = &FixStoreState::pack_fy;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      val.pack_choice = &FixStoreState::pack_fz;

    } else if (strcmp(arg[iarg],"q") == 0) {
      if (!atom->q_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_q;
    } else if (strcmp(arg[iarg],"mux") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_mux;
    } else if (strcmp(arg[iarg],"muy") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_muy;
    } else if (strcmp(arg[iarg],"muz") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_muz;
    } else if (strcmp(arg[iarg],"mu") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_mu;

    } else if (strcmp(arg[iarg],"radius") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_radius;
    } else if (strcmp(arg[iarg],"diameter") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_diameter;
    } else if (strcmp(arg[iarg],"omegax") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_omegax;
    } else if (strcmp(arg[iarg],"omegay") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_omegay;
    } else if (strcmp(arg[iarg],"omegaz") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_omegaz;
    } else if (strcmp(arg[iarg],"angmomx") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_angmomx;
    } else if (strcmp(arg[iarg],"angmomy") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_angmomy;
    } else if (strcmp(arg[iarg],"angmomz") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_angmomz;
    } else if (strcmp(arg[iarg],"tqx") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_tqx;
    } else if (strcmp(arg[iarg],"tqy") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_tqy;
    } else if (strcmp(arg[iarg],"tqz") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR, "Cannot use fix store/state {} for atom style {}",
                   arg[iarg], atom->get_style());
      val.pack_choice = &FixStoreState::pack_tqz;

    // compute or fix or variable or custom per-atom vector or array

    } else {
      ArgInfo argi(arg[iarg],ArgInfo::COMPUTE|ArgInfo::FIX|ArgInfo::VARIABLE
                   |ArgInfo::DNAME|ArgInfo::INAME);

      val.which = argi.get_type();
      val.argindex = argi.get_index1();
      val.id = argi.get_name();

      if (val.which == ArgInfo::NONE) break;
      if ((val.which == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
        error->all(FLERR,"Illegal fix store/state argument: {}", arg[iarg]);
    }
    values.push_back(val);
    iarg++;
  }

  // optional args

  comflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"com") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"fix store/state com", error);
      comflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Unknown fix store/state keyword: {}", arg[iarg]);
  }

  // error check

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR,"Compute ID {} for fix store/state does not exist", val.id);
      if (val.val.c->peratom_flag == 0)
        error->all(FLERR,"Fix store/state compute {} does not calculate per-atom values", val.id);
      if (val.argindex == 0 &&
          val.val.c->size_peratom_cols != 0)
        error->all(FLERR,"Fix store/state compute {} does not calculate per-atom vector", val.id);
      if (val.argindex && val.val.c->size_peratom_cols == 0)
        error->all(FLERR, "Fix store/state compute {} does not calculate per-atom array", val.id);
      if (val.argindex && (val.argindex > val.val.c->size_peratom_cols))
        error->all(FLERR, "Fix store/state compute array {} is accessed out-of-range", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f)
        error->all(FLERR, "Fix ID {} for fix store/state does not exist", val.id);
      if (val.val.f->peratom_flag == 0)
        error->all(FLERR, "Fix store/state fix {} does not calculate per-atom values", val.id);
      if (val.argindex == 0 && val.val.f->size_peratom_cols != 0)
        error->all(FLERR, "Fix store/state fix {} does not calculate per-atom vector", val.id);
      if (val.argindex && val.val.f->size_peratom_cols == 0)
        error->all(FLERR, "Fix store/state fix {} does not calculate per-atom array", val.id);
      if (val.argindex && (val.argindex > val.val.f->size_peratom_cols))
        error->all(FLERR, "Fix store/state fix {} array is accessed out-of-range", val.id);
      if (nevery % val.val.f->peratom_freq)
        error->all(FLERR, "Fix {} for fix store/state not computed at compatible time", val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR, "Variable name {} for fix store/state does not exist", val.id);
      if (input->variable->atomstyle(val.val.v) == 0)
        error->all(FLERR,"Fix store/state variable {} is not atom-style variable", val.id);

    } else if (val.which == ArgInfo::DNAME) {
      int iflag,icol;
      val.val.d = atom->find_custom(val.id.c_str(),iflag,icol);
      if (val.val.d < 0)
        error->all(FLERR,"Custom vector/array {} for fix store/state does not exist", val.id);
      if (val.argindex == 0) {
        if (!iflag || icol)
          error->all(FLERR, "Custom property {} for fix store/state is not double vector", val.id);
      } else {
        if (!iflag || !icol)
          error->all(FLERR, "Custom property {} for fix store/state is not double array", val.id);
        if (val.argindex > atom->dcols[val.val.d])
          error->all(FLERR, "Fix store/state custom array {} is accessed out-of-range", val.id);
      }

    } else if (val.which == ArgInfo::INAME) {
      int iflag,icol;
      val.val.i = atom->find_custom(val.id.c_str(),iflag,icol);
      if (val.val.i < 0)
        error->all(FLERR, "Custom vector/array {} for fix store/state does not exist", val.id);
      if (val.argindex == 0) {
        if (iflag || icol)
          error->all(FLERR, "Custom property {} for fix store/state is not integer vector", val.id);
      } else {
        if (iflag || !icol)
          error->all(FLERR, "Custom property {} for fix store/state is not integer array", val.id);
        if (val.argindex > atom->icols[val.val.i])
          error->all(FLERR, "Fix store/state custom array {} is accessed out-of-range", val.id);
      }
    }
  }

  // this fix produces either a per-atom vector or array

  peratom_flag = 1;
  if (values.size() == 1) size_peratom_cols = 0;
  else size_peratom_cols = values.size();

  // perform initial allocation of atom-based array
  // register with Atom class

  avalues = nullptr;
  FixStoreState::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // zero the array since dump may access it on timestep 0
  // zero the array since a variable may access it before first run

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    for (std::size_t m = 0; m < values.size(); m++)
      avalues[i][m] = 0.0;

  // store current values for keywords but not for compute, fix, variable

  kflag = 1;
  cfv_flag = 0;
  FixStoreState::end_of_step();
  firstflag = 1;
}

/* ---------------------------------------------------------------------- */

FixStoreState::~FixStoreState()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,Atom::GROW);
  atom->delete_callback(id,Atom::RESTART);

  memory->destroy(avalues);
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

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR,"Compute ID {} for fix store/state does not exist", val.id);

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f)
        error->all(FLERR,"Fix ID {} for fix store/state does not exist", val.id);

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR,"Variable name {} for fix store/state does not exist", val.id);

    } else if (val.which == ArgInfo::DNAME) {
      int iflag,cols;
      val.val.d = atom->find_custom(val.id.c_str(), iflag, cols);
      if (val.val.d < 0)
        error->all(FLERR,"Custom vector/array {} for fix store/state does not exist", val.id);

    } else if (val.which == ArgInfo::INAME) {
      int iflag,cols;
      val.val.i = atom->find_custom(val.id.c_str(), iflag, cols);
      if (val.val.i < 0)
        error->all(FLERR,"Custom vector/array {} for fix store/state does not exist", val.id);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::setup(int /*vflag*/)
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
  // compute com if comflag set

  if (comflag) {
    double masstotal = group->mass(igroup);
    group->xcm(igroup,masstotal,cm);
  }

  // if any compute/fix/variable and nevery, wrap with clear/add

  if (cfv_any && nevery) modify->clearstep_compute();

  // fill vector or array with per-atom values

  if (avalues) vbuf = &avalues[0][0];
  else vbuf = nullptr;

  int m = 0;
  for (auto &val : values) {
    if (val.which == ArgInfo::KEYWORD && kflag)
      (this->*val.pack_choice)(m);

    else if (cfv_flag) {

      int *mask = atom->mask;
      int nlocal = atom->nlocal;

      // invoke compute if not previously invoked, then access fields

      if (val.which == ArgInfo::COMPUTE) {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_PERATOM)) {
          val.val.c->compute_peratom();
          val.val.c->invoked_flag |= Compute::INVOKED_PERATOM;
        }

        if (val.argindex == 0) {
          double *compute_vector = val.val.c->vector_atom;
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) avalues[i][m] = compute_vector[i];
        } else {
          int jm1 = val.argindex - 1;
          double **compute_array = val.val.c->array_atom;
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) avalues[i][m] = compute_array[i][jm1];
        }

      // access fix fields, guaranteed to be ready

      } else if (val.which == ArgInfo::FIX) {
        if (val.argindex == 0) {
          double *fix_vector = val.val.f->vector_atom;
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) avalues[i][m] = fix_vector[i];
        } else {
          int jm1 = val.argindex - 1;
          double **fix_array = val.val.f->array_atom;
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) avalues[i][m] = fix_array[i][jm1];
        }

      // evaluate atom-style variable

      } else if (val.which == ArgInfo::VARIABLE) {
        input->variable->compute_atom(val.val.v, igroup, &avalues[0][m], values.size(),0);

      // access custom atom vector/array fields

      } else if (val.which == ArgInfo::DNAME) {
        if (val.argindex == 0) {
          double *dvector = atom->dvector[val.val.d];
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) avalues[i][m] = dvector[i];
        } else {
          double **darray = atom->darray[val.val.d];
          int jm1 = val.argindex - 1;
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) avalues[i][m] = darray[i][jm1];
        }

      } else if (val.which == ArgInfo::INAME) {
        if (val.argindex == 0) {
          int *ivector = atom->ivector[val.val.i];
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) avalues[i][m] = ivector[i];
        } else {
          int **iarray = atom->iarray[val.val.i];
          int jm1 = val.argindex - 1;
          for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) avalues[i][m] = iarray[i][jm1];
        }
      }
    }
    ++m;
  }

  // if any compute/fix/variable and nevery, wrap with clear/add

  if (cfv_any && nevery) {
    const bigint nextstep = (update->ntimestep/nevery)*nevery + nevery;
    modify->addstep_compute(nextstep);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixStoreState::memory_usage()
{
  double bytes = (double)atom->nmax*values.size() * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixStoreState::grow_arrays(int nmax)
{
  memory->grow(avalues,nmax,values.size(),"store/state:avalues");
  if (values.size() == 1) {
    if (nmax) vector_atom = &avalues[0][0];
    else vector_atom = nullptr;
  } else array_atom = avalues;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixStoreState::copy_arrays(int i, int j, int /*delflag*/)
{
  for (std::size_t m = 0; m < values.size(); m++) avalues[j][m] = avalues[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixStoreState::pack_exchange(int i, double *buf)
{
  for (std::size_t m = 0; m < values.size(); m++) buf[m] = avalues[i][m];
  return values.size();
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixStoreState::unpack_exchange(int nlocal, double *buf)
{
  for (std::size_t m = 0; m < values.size(); m++) avalues[nlocal][m] = buf[m];
  return values.size();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixStoreState::pack_restart(int i, double *buf)
{
  // pack buf[0] this way because other fixes unpack it
  buf[0] = values.size()+1;
  for (std::size_t m = 0; m < values.size(); m++) buf[m+1] = avalues[i][m];
  return values.size()+1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixStoreState::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int>(extra[nlocal][m]);
  m++;

  for (std::size_t i = 0; i < values.size(); i++) avalues[nlocal][i] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixStoreState::maxsize_restart()
{
  return values.size()+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixStoreState::size_restart(int /*nlocal*/)
{
  return values.size()+1;
}

/* ----------------------------------------------------------------------
   one method for every keyword fix store/state can archive
   the atom property is packed into buf starting at n with stride values.size()
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
      n += values.size();
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) vbuf[n] = mass[type[i]];
      else vbuf[n] = 0.0;
      n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_xsu(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double boxxlo = domain->boxlo[0];
  double invxprd = 1.0/domain->xprd;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      vbuf[n] = (x[i][0]-boxxlo)*invxprd + ((image[i] & IMGMASK) - IMGMAX);
    else vbuf[n] = 0.0;
    n += values.size();
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_ysu(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double boxylo = domain->boxlo[1];
  double invyprd = 1.0/domain->yprd;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      vbuf[n] = (x[i][1]-boxylo)*invyprd +
        (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    else vbuf[n] = 0.0;
    n += values.size();
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_zsu(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double boxzlo = domain->boxlo[2];
  double invzprd = 1.0/domain->zprd;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      vbuf[n] = (x[i][2]-boxzlo)*invzprd + (image[i] >> IMG2BITS) - IMGMAX;
    else vbuf[n] = 0.0;
    n += values.size();
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_xsu_triclinic(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      vbuf[n] = h_inv[0]*(x[i][0]-boxlo[0]) + h_inv[5]*(x[i][1]-boxlo[1]) +
        h_inv[4]*(x[i][2]-boxlo[2]) + (image[i] & IMGMASK) - IMGMAX;
    else vbuf[n] = 0.0;
    n += values.size();
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_ysu_triclinic(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      vbuf[n] = h_inv[1]*(x[i][1]-boxlo[1]) + h_inv[3]*(x[i][2]-boxlo[2]) +
        (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    else vbuf[n] = 0.0;
    n += values.size();
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_zsu_triclinic(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      vbuf[n] = h_inv[2]*(x[i][2]-boxlo[2]) + (image[i] >> IMG2BITS) - IMGMAX;
    else vbuf[n] = 0.0;
    n += values.size();
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
    n += values.size();
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_iy(int n)
{
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      vbuf[n] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    else vbuf[n] = 0.0;
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_mu(int n)
{
  double **mu = atom->mu;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = mu[i][3];
    else vbuf[n] = 0.0;
    n += values.size();
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
    n += values.size();
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreState::pack_diameter(int n)
{
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) vbuf[n] = 2.0*radius[i];
    else vbuf[n] = 0.0;
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
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
    n += values.size();
  }
}
