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

#include "fix.h"
#include <cstring>
#include <cctype>
#include "atom.h"
#include "group.h"
#include "force.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

// allocate space for static class instance variable and initialize it

int Fix::instance_total = 0;

/* ---------------------------------------------------------------------- */

Fix::Fix(LAMMPS *lmp, int /*narg*/, char **arg) :
  Pointers(lmp),
  id(NULL), style(NULL), extlist(NULL), vector_atom(NULL), array_atom(NULL),
  vector_local(NULL), array_local(NULL), eatom(NULL), vatom(NULL)
{
  instance_me = instance_total++;

  // fix ID, group, and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,"Fix ID must be alphanumeric or underscore characters");

  igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find fix group ID");
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  restart_global = restart_peratom = restart_file = 0;
  force_reneighbor = 0;
  box_change_size = box_change_shape = box_change_domain = 0;
  thermo_energy = 0;
  thermo_virial = 0;
  rigid_flag = 0;
  peatom_flag = 0;
  virial_flag = 0;
  no_change_box = 0;
  time_integrate = 0;
  time_depend = 0;
  create_attribute = 0;
  restart_pbc = 0;
  wd_header = wd_section = 0;
  dynamic_group_allow = 0;
  dynamic = 0;
  dof_flag = 0;
  special_alter_flag = 0;
  enforce2d_flag = 0;
  respa_level_support = 0;
  respa_level = -1;

  scalar_flag = vector_flag = array_flag = 0;
  peratom_flag = local_flag = 0;
  global_freq = local_freq = peratom_freq = -1;
  size_vector_variable = size_array_rows_variable = 0;

  comm_forward = comm_reverse = comm_border = 0;
  restart_reset = 0;

  // reasonable defaults
  // however, each fix that uses these values should explicitly set them

  nevery = 1;
  global_freq = 1;

  // per-atom virial
  // set vflag_atom = 0 b/c some fixes grow vatom in grow_arrays()
  //   which may occur outside of timestepping

  maxeatom = maxvatom = 0;
  vflag_atom = 0;

  // KOKKOS per-fix data masks

  execution_space = Host;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;

  kokkosable = 0;
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

Fix::~Fix()
{
  if (copymode) return;

  delete [] id;
  delete [] style;
  memory->destroy(eatom);
  memory->destroy(vatom);
}

/* ----------------------------------------------------------------------
   process params common to all fixes here
   if unknown param, call modify_param specific to the fix
------------------------------------------------------------------------- */

void Fix::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal fix_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"dynamic/dof") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) dynamic = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) dynamic = 1;
      else error->all(FLERR,"Illegal fix_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) thermo_energy = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) {
        if (!(THERMO_ENERGY & setmask()))
          error->all(FLERR,"Illegal fix_modify command");
        thermo_energy = 1;
      } else error->all(FLERR,"Illegal fix_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"virial") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) thermo_virial = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) {
        if (virial_flag == 0)
          error->all(FLERR,"Illegal fix_modify command");
        thermo_virial = 1;
      } else error->all(FLERR,"Illegal fix_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"respa") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (!respa_level_support) error->all(FLERR,"Illegal fix_modify command");
      int lvl = force->inumeric(FLERR,arg[iarg+1]);
      if (lvl < 0) error->all(FLERR,"Illegal fix_modify command");
      respa_level = lvl-1;
      iarg += 2;
    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0) error->all(FLERR,"Illegal fix_modify command");
      iarg += n;
    }
  }
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for values of eflag (0-3) and vflag (0-6)
   fixes call this if they use ev_tally()
------------------------------------------------------------------------- */

void Fix::ev_setup(int eflag, int vflag)
{
  int i,n;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag % 2;
  eflag_atom = eflag / 2;

  vflag_either = vflag;
  vflag_global = vflag % 4;
  vflag_atom = vflag / 4;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nlocal > maxeatom) {
    maxeatom = atom->nmax;
    memory->destroy(eatom);
    memory->create(eatom,maxeatom,"fix:eatom");
  }
  if (vflag_atom && atom->nlocal > maxvatom) {
    maxvatom = atom->nmax;
    memory->destroy(vatom);
    memory->create(vatom,maxvatom,6,"fix:vatom");
  }

  // zero accumulators
  // no global energy variable to zero (unlike pair,bond,angle,etc)
  // fixes tally it individually via fix_modify energy yes and compute_scalar()

  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++) eatom[i] = 0.0;
  }
  if (vflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   if thermo_virial is on:
     setup for virial computation
     see integrate::ev_set() for values of vflag (0-6)
     fixes call this if use v_tally()
   else: set evflag=0
------------------------------------------------------------------------- */

void Fix::v_setup(int vflag)
{
  int i,n;

  if (!thermo_virial) {
    evflag = 0;
    return;
  }

  evflag = 1;

  vflag_global = vflag % 4;
  vflag_atom = vflag / 4;

  // reallocate per-atom array if necessary

  if (vflag_atom && atom->nlocal > maxvatom) {
    maxvatom = atom->nmax;
    memory->destroy(vatom);
    memory->create(vatom,maxvatom,6,"fix:vatom");
  }

  // zero accumulators

  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (vflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   tally per-atom energy and global/per-atom virial into accumulators
   n = # of local owned atoms involved, with local indices in list
   eng = total energy for the interaction involving total atoms
   v = total virial for the interaction involving total atoms
   increment per-atom energy of each atom in list by 1/total fraction
   v_tally tallies virial
   this method can be used when fix computes energy/forces in post_force()
     e.g. fix cmap: compute energy and virial only on owned atoms
       whether newton_bond is on or off
     other procs will tally left-over fractions for atoms they own
------------------------------------------------------------------------- */

void Fix::ev_tally(int n, int *list, double total, double eng, double *v)
{
  if (eflag_atom) {
    double fraction = eng/total;
    for (int i = 0; i < n; i++)
      eatom[list[i]] += fraction;
  }

  v_tally(n,list,total,v);
}


/* ----------------------------------------------------------------------
   tally virial into global and per-atom accumulators
   n = # of local owned atoms involved, with local indices in list
   v = total virial for the interaction involving total atoms
   increment global virial by n/total fraction
   increment per-atom virial of each atom in list by 1/total fraction
   this method can be used when fix computes forces in post_force()
     e.g. fix shake, fix rigid: compute virial only on owned atoms
       whether newton_bond is on or off
     other procs will tally left-over fractions for atoms they own
------------------------------------------------------------------------- */

void Fix::v_tally(int n, int *list, double total, double *v)
{
  int m;

  if (vflag_global) {
    double fraction = n/total;
    virial[0] += fraction*v[0];
    virial[1] += fraction*v[1];
    virial[2] += fraction*v[2];
    virial[3] += fraction*v[3];
    virial[4] += fraction*v[4];
    virial[5] += fraction*v[5];
  }

  if (vflag_atom) {
    double fraction = 1.0/total;
    for (int i = 0; i < n; i++) {
      m = list[i];
      vatom[m][0] += fraction*v[0];
      vatom[m][1] += fraction*v[1];
      vatom[m][2] += fraction*v[2];
      vatom[m][3] += fraction*v[3];
      vatom[m][4] += fraction*v[4];
      vatom[m][5] += fraction*v[5];
    }
  }
}

/* ----------------------------------------------------------------------
   tally virial into global and per-atom accumulators
   i = local index of atom
   v = total virial for the interaction
   increment global virial by v
   increment per-atom virial by v
   this method can be used when fix computes forces in post_force()
   and the force depends on a distance to some external object
     e.g. fix wall/lj93: compute virial only on owned atoms
------------------------------------------------------------------------- */

void Fix::v_tally(int i, double *v)
{
  if (vflag_global) {
    virial[0] += v[0];
    virial[1] += v[1];
    virial[2] += v[2];
    virial[3] += v[3];
    virial[4] += v[4];
    virial[5] += v[5];
  }

  if (vflag_atom) {
    vatom[i][0] += v[0];
    vatom[i][1] += v[1];
    vatom[i][2] += v[2];
    vatom[i][3] += v[3];
    vatom[i][4] += v[4];
    vatom[i][5] += v[5];
  }
}

/* ----------------------------------------------------------------------
   tally virial component into global and per-atom accumulators
   n = index of virial component (0-5)
   i = local index of atom
   vn = nth component of virial for the interaction
   increment nth component of global virial by vn
   increment nth component of per-atom virial by vn
   this method can be used when fix computes forces in post_force()
   and the force depends on a distance to some external object
     e.g. fix wall/lj93: compute virial only on owned atoms
------------------------------------------------------------------------- */

void Fix::v_tally(int n, int i, double vn)
{
  if (vflag_global)
    virial[n] += vn;

  if (vflag_atom)
    vatom[i][n] += vn;
}
