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

#include "fix.h"

#include "atom.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

// allocate space for static class instance variable and initialize it

int Fix::instance_total = 0;

/* ---------------------------------------------------------------------- */

Fix::Fix(LAMMPS *lmp, int /*narg*/, char **arg) :
  Pointers(lmp),
  id(nullptr), style(nullptr), extlist(nullptr), vector_atom(nullptr), array_atom(nullptr),
  vector_local(nullptr), array_local(nullptr), eatom(nullptr), vatom(nullptr),
  cvatom(nullptr)
{
  instance_me = instance_total++;

  // fix ID, group, and style
  // ID must be all alphanumeric chars or underscores

  id = utils::strdup(arg[0]);
  if (!utils::is_id(id))
    error->all(FLERR,"Fix ID must be alphanumeric or underscore characters");

  igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find fix group ID");
  groupbit = group->bitmask[igroup];

  style = utils::strdup(arg[2]);

  restart_global = restart_peratom = restart_file = 0;
  force_reneighbor = 0;
  box_change = NO_BOX_CHANGE;
  thermo_energy = 0;
  thermo_virial = 0;
  energy_global_flag = energy_peratom_flag = 0;
  virial_global_flag = virial_peratom_flag = 0;
  ecouple_flag = 0;
  rigid_flag = 0;
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
  maxexchange = 0;
  maxexchange_dynamic = 0;
  pre_exchange_migrate = 0;
  stores_ids = 0;

  scalar_flag = vector_flag = array_flag = 0;
  peratom_flag = local_flag = pergrid_flag = 0;
  global_freq = local_freq = peratom_freq = pergrid_freq = -1;
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

  maxeatom = maxvatom = maxcvatom = 0;
  vflag_atom = cvflag_atom = 0;
  centroidstressflag = CENTROID_SAME;

  // KOKKOS package

  execution_space = Host;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;

  kokkosable = copymode = 0;
  forward_comm_device = exchange_comm_device = sort_device = 0;
  fuse_integrate_flag = 0;
}

/* ---------------------------------------------------------------------- */

Fix::~Fix()
{
  if (copymode) return;

  delete [] id;
  delete [] style;
  memory->destroy(eatom);
  memory->destroy(vatom);
  memory->destroy(cvatom);
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
      dynamic = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      thermo_energy = utils::logical(FLERR,arg[iarg+1],false,lmp);
      if (thermo_energy && !energy_global_flag && !energy_peratom_flag)
          error->all(FLERR,"Illegal fix_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"virial") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      thermo_virial = utils::logical(FLERR,arg[iarg+1],false,lmp);
      if (thermo_virial && !virial_global_flag && !virial_peratom_flag)
        error->all(FLERR,"Illegal fix_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"respa") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (!respa_level_support) error->all(FLERR,"Illegal fix_modify command");
      int lvl = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
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

void::Fix::set_molecule(int, tagint, int, double *, double *, double *)
{
  error->all(FLERR,"Molecule update not implemented for fix {}", style);
}

/* ----------------------------------------------------------------------
   setup for peratom energy and global/peratom virial computation
   see integrate::ev_set() for values of eflag (0-3) and vflag (0-6)
   fixes call Fix::ev_init() if tally energy and virial values
   if thermo_energy is not set, energy tallying is disabled
   if thermo_virial is not set, virial tallying is disabled
   global energy is tallied separately, output by compute_scalar() method
------------------------------------------------------------------------- */

void Fix::ev_setup(int eflag, int vflag)
{
  int i,n;

  evflag = 1;

  if (!thermo_energy) eflag_either = eflag_global = eflag_atom = 0;
  else {
    eflag_either = eflag;
    eflag_global = eflag & ENERGY_GLOBAL;
    eflag_atom = eflag & ENERGY_ATOM;
  }

  if (!thermo_virial) vflag_either = vflag_global = vflag_atom = 0;
  else {
    vflag_either = vflag;
    vflag_global = vflag & (VIRIAL_PAIR | VIRIAL_FDOTR);
    if (centroidstressflag != CENTROID_AVAIL) {
      vflag_atom = vflag & (VIRIAL_ATOM | VIRIAL_CENTROID);
      cvflag_atom = 0;
    } else {
      vflag_atom = vflag & VIRIAL_ATOM;
      cvflag_atom = vflag & VIRIAL_CENTROID;
    }
  }

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
  if (cvflag_atom && atom->nlocal > maxcvatom) {
    maxcvatom = atom->nmax;
    memory->destroy(cvatom);
    memory->create(cvatom,maxcvatom,9,"fix:cvatom");
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
  if (cvflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++) {
      cvatom[i][0] = 0.0;
      cvatom[i][1] = 0.0;
      cvatom[i][2] = 0.0;
      cvatom[i][3] = 0.0;
      cvatom[i][4] = 0.0;
      cvatom[i][5] = 0.0;
      cvatom[i][6] = 0.0;
      cvatom[i][7] = 0.0;
      cvatom[i][8] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   setup for global/peratom virial computation
   see integrate::ev_set() for values of vflag (0-6)
   fixes call Fix::v_init() if tally virial values but not energy
   if thermo_virial is not set, virial tallying is disabled
------------------------------------------------------------------------- */

void Fix::v_setup(int vflag)
{
  int i,n;

  evflag = 1;
  vflag_global = vflag & (VIRIAL_PAIR | VIRIAL_FDOTR);
  if (centroidstressflag != CENTROID_AVAIL) {
    vflag_atom = vflag & (VIRIAL_ATOM | VIRIAL_CENTROID);
    cvflag_atom = 0;
  } else {
    vflag_atom = vflag & VIRIAL_ATOM;
    cvflag_atom = vflag & VIRIAL_CENTROID;
  }

  // reallocate per-atom array if necessary

  if (vflag_atom && atom->nlocal > maxvatom) {
    maxvatom = atom->nmax;
    memory->destroy(vatom);
    memory->create(vatom,maxvatom,6,"fix:vatom");
  }
  if (cvflag_atom && atom->nlocal > maxcvatom) {
    maxcvatom = atom->nmax;
    memory->destroy(cvatom);
    memory->create(cvatom,maxcvatom,9,"fix:cvatom");
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
  if (cvflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++) {
      cvatom[i][0] = 0.0;
      cvatom[i][1] = 0.0;
      cvatom[i][2] = 0.0;
      cvatom[i][3] = 0.0;
      cvatom[i][4] = 0.0;
      cvatom[i][5] = 0.0;
      cvatom[i][6] = 0.0;
      cvatom[i][7] = 0.0;
      cvatom[i][8] = 0.0;
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
   n = # of local owned atoms involved, with local indices in list
   vtot = total virial for the interaction involving total atoms
   rlist = list of positional vectors
   flist = list of force vectors
   center = centroid coordinate
   increment global virial by n/total fraction
   increment per-atom virial of each atom in list by 1/total fraction
   add centroid form atomic virial contribution for each atom if available
   this method can be used when fix computes forces in post_force()
   and only total forces on each atom in group are easily available
     e.g. fix rigid/small: compute virial only on owned atoms
       whether newton_bond is on or off
     other procs will tally left-over fractions for atoms they own
------------------------------------------------------------------------- */

void Fix::v_tally(int n, int *list, double total, double *vtot,
    double rlist[][3], double flist[][3], double center[])
{

  v_tally(n, list, total, vtot);

  if (cvflag_atom) {
    for (int i = 0; i< n; i++) {
      const double ri0[3] = {
        rlist[i][0]-center[0],
        rlist[i][1]-center[1],
        rlist[i][2]-center[2],
      };
      cvatom[list[i]][0] += ri0[0]*flist[i][0];
      cvatom[list[i]][1] += ri0[1]*flist[i][1];
      cvatom[list[i]][2] += ri0[2]*flist[i][2];
      cvatom[list[i]][3] += ri0[0]*flist[i][1];
      cvatom[list[i]][4] += ri0[0]*flist[i][2];
      cvatom[list[i]][5] += ri0[1]*flist[i][2];
      cvatom[list[i]][6] += ri0[1]*flist[i][0];
      cvatom[list[i]][7] += ri0[2]*flist[i][0];
      cvatom[list[i]][8] += ri0[2]*flist[i][1];
    }
  }

}

/* ----------------------------------------------------------------------
   tally virial into global and per-atom accumulators
   n = # of local owned atoms involved, with local indices in list
   vtot = total virial for the interaction involving total atoms
   npair = # of atom pairs with forces beween them
   pairlist = indice list of pairs
   fpairlist = forces between pairs
   dellist = displacement vectors between pairs
   increment global virial by n/total fraction
   increment per-atom virial of each atom in list by 1/total fraction
   add centroid form atomic virial contribution for each atom if available
   this method can be used when fix computes forces in post_force()
     e.g. fix shake, fix rigid: compute virial only on owned atoms
       whether newton_bond is on or off
     other procs will tally left-over fractions for atoms they own
------------------------------------------------------------------------- */

void Fix::v_tally(int n, int *list, double total, double *vtot, int nlocal,
    int npair, int pairlist[][2], double *fpairlist, double dellist[][3])
{

  v_tally(n, list, total, vtot);

  if (cvflag_atom) {
    double v[6];
    for (int i = 0; i < npair; i++) {
      v[0] = 0.5*dellist[i][0]*dellist[i][0]*fpairlist[i];
      v[1] = 0.5*dellist[i][1]*dellist[i][1]*fpairlist[i];
      v[2] = 0.5*dellist[i][2]*dellist[i][2]*fpairlist[i];
      v[3] = 0.5*dellist[i][0]*dellist[i][1]*fpairlist[i];
      v[4] = 0.5*dellist[i][0]*dellist[i][2]*fpairlist[i];
      v[5] = 0.5*dellist[i][1]*dellist[i][2]*fpairlist[i];
      const int i0 = pairlist[i][0];
      const int i1 = pairlist[i][1];
      if (i0 < nlocal) {
        cvatom[i0][0] += v[0];
        cvatom[i0][1] += v[1];
        cvatom[i0][2] += v[2];
        cvatom[i0][3] += v[3];
        cvatom[i0][4] += v[4];
        cvatom[i0][5] += v[5];
        cvatom[i0][6] += v[3];
        cvatom[i0][7] += v[4];
        cvatom[i0][8] += v[5];
      }
      if (i1 < nlocal) {
        cvatom[i1][0] += v[0];
        cvatom[i1][1] += v[1];
        cvatom[i1][2] += v[2];
        cvatom[i1][3] += v[3];
        cvatom[i1][4] += v[4];
        cvatom[i1][5] += v[5];
        cvatom[i1][6] += v[3];
        cvatom[i1][7] += v[4];
        cvatom[i1][8] += v[5];
      }
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
