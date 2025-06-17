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

#include "set.h"

#include "arg_info.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_body.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "random_park.h"
#include "region.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

enum{SETCOMMAND,FIXSET};         // also used in FixSet class

enum{ATOM_SELECT,MOL_SELECT,TYPE_SELECT,GROUP_SELECT,REGION_SELECT};

enum{ANGLE,ANGMOM,BOND,CC,CHARGE,DENSITY,DIAMETER,DIHEDRAL,DIPOLE,
  DIPOLE_RANDOM,DPD_THETA,EDPD_CV,EDPD_TEMP,EPSILON,IMAGE,IMPROPER,LENGTH,
  MASS,MOLECULE,OMEGA,QUAT,QUAT_RANDOM,RADIUS_ELECTRON,SHAPE,
  SMD_CONTACT_RADIUS,SMD_MASS_DENSITY,SPH_CV,SPH_E,SPH_RHO,
  SPIN_ATOM,SPIN_ATOM_RANDOM,SPIN_ELECTRON,TEMPERATURE,THETA,THETA_RANDOM,
  TRI,TYPE,TYPE_FRACTION,TYPE_RATIO,TYPE_SUBSET,VOLUME,VX,VY,VZ,X,Y,Z,
  IVEC,DVEC,IARRAY,DARRAY};

#define DELTA 4

/* ---------------------------------------------------------------------- */

Set::Set(class LAMMPS *lmp) :
    Command(lmp), id(nullptr), region(nullptr), actions(nullptr), invoke_choice(nullptr),
    vec1(nullptr), vec2(nullptr), vec3(nullptr), vec4(nullptr), select(nullptr)
{
  maxselect = maxvariable = 0;
}

/* ---------------------------------------------------------------------- */

Set::~Set()
{
  memory->sfree(actions);
  memory->sfree(invoke_choice);

  memory->destroy(select);
  memory->destroy(vec1);
  memory->destroy(vec2);
  memory->destroy(vec3);
  memory->destroy(vec4);
}

/* ---------------------------------------------------------------------- */

void Set::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, Error::NOLASTLINE,
               "Set command before simulation box is defined" + utils::errorurl(0));
  if (atom->natoms == 0)
    error->all(FLERR, Error::NOLASTLINE, "Set command on system without atoms");
  if (narg < 4) error->all(FLERR, 1, "Illegal set command: need at least four arguments");

  process_args(SETCOMMAND, narg, arg);

  if (comm->me == 0) utils::logmesg(lmp, "Setting atom values ...\n");

  // select which atoms to act on

  selection(atom->nlocal);

  // loop over list of actions to reset atom attributes

  invoke_actions();

  // print stats for each action
  // for CC option, include species index

  bigint bcount, allcount;

  for (int i = 0; i < naction; i++) {
    Action *action = &actions[i];
    int iarg = action->argindex;

    if (action->count_action < 0) {
      bcount = action->count_select;
    } else {
      bcount = action->count_action;
    }
    MPI_Allreduce(&bcount, &allcount, 1, MPI_LMP_BIGINT, MPI_SUM, world);

    if (comm->me == 0) {
      if (strcmp(arg[iarg], "cc") == 0) {
        utils::logmesg(lmp, "  {} settings made for {} index {}\n", allcount, arg[iarg],
                       arg[iarg + 1]);
      } else {
        utils::logmesg(lmp, "  {} settings made for {}\n", allcount, arg[iarg]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   process args of set command or fix set command
------------------------------------------------------------------------- */

void Set::process_args(int caller_flag, int narg, char **arg)
{
  caller = caller_flag;

  // style and ID info

  id = utils::strdup(arg[1]);

  if (strcmp(arg[0], "atom") == 0) {
    style = ATOM_SELECT;
    if (atom->tag_enable == 0) error->all(FLERR, "Cannot use set atom with no atom IDs defined");
    utils::bounds(FLERR, id, 1, MAXTAGINT, nlobig, nhibig, error);

  } else if (strcmp(arg[0], "mol") == 0) {
    style = MOL_SELECT;
    if (atom->molecule_flag == 0)
      error->all(FLERR, "Cannot use set mol with no molecule IDs defined");
    utils::bounds(FLERR, id, 1, MAXTAGINT, nlobig, nhibig, error);

  } else if (strcmp(arg[0], "type") == 0) {
    style = TYPE_SELECT;
    if (char *typestr = utils::expand_type(FLERR, id, Atom::ATOM, lmp)) {
      delete[] id;
      id = typestr;
    }
    utils::bounds(FLERR, id, 1, atom->ntypes, nlo, nhi, error);

  } else if (strcmp(arg[0], "group") == 0) {
    style = GROUP_SELECT;
    int igroup = group->find(id);
    if (igroup == -1) error->all(FLERR, "Could not find set group ID {}", id);
    groupbit = group->bitmask[igroup];

  } else if (strcmp(arg[0], "region") == 0) {
    style = REGION_SELECT;
    region = domain->get_region_by_id(id);
    if (!region) error->all(FLERR, "Set region {} does not exist", id);

  } else
    error->all(FLERR, "Unknown set or fix set command style: {}", arg[0]);

  delete[] id;

  // loop over remaining keyword/value pairs to create list of actions
  // one action = keyword/value pair

  naction = maxaction = 0;
  actions = nullptr;
  Action *action;

  int iarg = 2;
  while (iarg < narg) {

    // grow list of actions if needed

    if (naction == maxaction) {
      maxaction += DELTA;
      actions = (Action *) memory->srealloc(actions,maxaction*sizeof(Action),"set:actions");
      invoke_choice = (FnPtrPack *)
        memory->srealloc(invoke_choice,maxaction*sizeof(FnPtrPack),"set:invoke_choice");
    }

    action = &actions[naction];
    action->argindex = iarg;
    action->varflag = 0;
    action->varflag1 = action->varflag2 = action->varflag3 = action->varflag4 = 0;

    if (strcmp(arg[iarg],"angle") == 0) {
      action->keyword = ANGLE;
      process_angle(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_angle;
    } else if (strcmp(arg[iarg],"angmom") == 0) {
      action->keyword = ANGMOM;
      process_angmom(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_angmom;
    } else if (strcmp(arg[iarg],"bond") == 0) {
      action->keyword = BOND;
      process_bond(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_bond;
    } else if (strcmp(arg[iarg],"cc") == 0) {
      action->keyword = CC;
      process_cc(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_cc;
    } else if (strcmp(arg[iarg],"charge") == 0) {
      action->keyword = CHARGE;
      process_charge(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_charge;
    } else if (strcmp(arg[iarg],"density") == 0 ||(strcmp(arg[iarg],"density/disc") == 0)) {
      action->keyword = DENSITY;
      process_density(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_density;
    } else if (strcmp(arg[iarg],"diameter") == 0) {
      action->keyword = DIAMETER;
      process_diameter(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_diameter;
    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      action->keyword = DIHEDRAL;
      process_dihedral(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_dihedral;
    } else if (strcmp(arg[iarg],"dipole") == 0) {
      action->keyword = DIPOLE;
      process_dipole(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_dipole;
    } else if (strcmp(arg[iarg],"dipole/random") == 0) {
      action->keyword = DIPOLE_RANDOM;
      process_dipole_random(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_dipole_random;
    } else if (strcmp(arg[iarg],"dpd/theta") == 0) {
      action->keyword = DPD_THETA;
      process_dpd_theta(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_dpd_theta;
    } else if (strcmp(arg[iarg],"edpd/cv") == 0) {
      action->keyword = EDPD_CV;
      process_edpd_cv(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_edpd_cv;
    } else if (strcmp(arg[iarg],"edpd/temp") == 0) {
      action->keyword = EDPD_TEMP;
      process_edpd_temp(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_edpd_temp;
    } else if (strcmp(arg[iarg],"epsilon") == 0) {
      action->keyword = EPSILON;
      process_epsilon(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_epsilon;
    } else if (strcmp(arg[iarg],"image") == 0) {
      action->keyword = IMAGE;
      process_image(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_image;
    } else if (strcmp(arg[iarg],"improper") == 0) {
      action->keyword = IMPROPER;
      process_improper(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_improper;
    } else if (strcmp(arg[iarg],"length") == 0) {
      action->keyword = LENGTH;
      process_length(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_length;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      action->keyword = MASS;
      process_mass(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_mass;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      action->keyword = MOLECULE;
      process_mol(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_mol;
    } else if (strcmp(arg[iarg],"omega") == 0) {
      action->keyword = OMEGA;
      process_omega(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_omega;
    } else if (strcmp(arg[iarg],"quat") == 0) {
      action->keyword = QUAT;
      process_quat(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_quat;
    } else if (strcmp(arg[iarg],"quat/random") == 0) {
      action->keyword = QUAT_RANDOM;
      process_quat_random(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_quat_random;
    } else if (strcmp(arg[iarg],"radius/electron") == 0) {
      action->keyword = RADIUS_ELECTRON;
      process_radius_election(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_radius_election;
    } else if (strcmp(arg[iarg],"shape") == 0) {
      action->keyword = SHAPE;
      process_shape(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_shape;
    } else if (strcmp(arg[iarg],"smd/contact/radius") == 0) {
      action->keyword = SMD_CONTACT_RADIUS;
      process_smd_contact_radius(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_smd_contact_radius;
    } else if (strcmp(arg[iarg],"smd/mass/density") == 0) {
      action->keyword = SMD_MASS_DENSITY;
      process_smd_mass_density(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_smd_mass_density;
    } else if (strcmp(arg[iarg],"sph/cv") == 0) {
      action->keyword = SPH_CV;
      process_sph_cv(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_sph_cv;
    } else if (strcmp(arg[iarg],"sph/e") == 0) {
      action->keyword = SPH_E;
      process_sph_e(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_sph_e;
    } else if (strcmp(arg[iarg],"sph/rho") == 0) {
      action->keyword = SPH_RHO;
      process_sph_rho(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_sph_rho;
    } else if ((strcmp(arg[iarg],"spin/atom") == 0) || (strcmp(arg[iarg],"spin") == 0)) {
      action->keyword = SPIN_ATOM;
      process_spin_atom(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_spin_atom;
    } else if ((strcmp(arg[iarg],"spin/atom/random") == 0) || (strcmp(arg[iarg],"spin/random") == 0)) {
      action->keyword = SPIN_ATOM_RANDOM;
      process_spin_atom_random(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_spin_atom_random;
    } else if (strcmp(arg[iarg],"spin/electron") == 0) {
      action->keyword = SPIN_ELECTRON;
      process_spin_electron(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_spin_electron;
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      action->keyword = TEMPERATURE;
      process_temperature(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_temperature;
    } else if (strcmp(arg[iarg],"theta") == 0) {
      action->keyword = THETA;
      process_theta(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_theta;
    } else if (strcmp(arg[iarg],"theta/random") == 0) {
      action->keyword = THETA_RANDOM;
      process_theta_random(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_theta_random;
    } else if (strcmp(arg[iarg],"tri") == 0) {
      action->keyword = TRI;
      process_tri(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_tri;
    } else if (strcmp(arg[iarg],"type") == 0) {
      action->keyword = TYPE;
      process_type(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_type;
    } else if (strcmp(arg[iarg],"type/fraction") == 0) {
      action->keyword = TYPE_FRACTION;
      process_type_fraction(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_type_fraction;
    } else if (strcmp(arg[iarg],"type/ratio") == 0) {
      action->keyword = TYPE_RATIO;
      process_type_ratio(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_type_ratio;
    } else if (strcmp(arg[iarg],"type/subset") == 0) {
      action->keyword = TYPE_SUBSET;
      process_type_subset(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_type_subset;
    } else if (strcmp(arg[iarg],"volume") == 0) {
      action->keyword = VOLUME;
      process_volume(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_volume;
    } else if (strcmp(arg[iarg],"vx") == 0) {
      action->keyword = VX;
      process_vx(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_vx;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      action->keyword = VY;
      process_vy(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_vy;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      action->keyword = VZ;
      process_vz(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_vz;
    } else if (strcmp(arg[iarg],"x") == 0) {
      action->keyword = X;
      process_x(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_x;
    } else if (strcmp(arg[iarg],"y") == 0) {
      action->keyword = Y;
      process_y(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_y;
    } else if (strcmp(arg[iarg],"z") == 0) {
      action->keyword = Z;
      process_z(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_z;

    } else if (utils::strmatch(arg[iarg],"^[id]2?_")) {
      process_custom(iarg,narg,arg,action);
      invoke_choice[naction++] = &Set::invoke_custom;

    // unrecognized keyword

    } else error->all(FLERR,"Unrecognized set or fix set command keyword {}",arg[iarg]);
  }

  // varflag = 1 if any action uses a per-atom variable
  // varflag1234 = 1 if any action uses a specific per-atom variable

  varflag = 0;
  varflag1 = varflag2 = varflag3 = varflag4 = 0;
  for (int i = 0; i < naction; i++) {
    action = &actions[i];
    if (action->varflag) varflag = 1;
    if (action->varflag1) varflag1 = 1;
    if (action->varflag2) varflag2 = 1;
    if (action->varflag3) varflag3 = 1;
    if (action->varflag4) varflag4 = 1;
  }

  // error if any action of fix set command does not use a per-atom variable
  // b/c fix set is then effectivly a no-op

  if (caller == FIXSET) {
    for (int i = 0; i < naction; i++) {
      action = &actions[i];
      if (!action->varflag)
        error->all(FLERR,"Fix set command keyword {} does not invoke a per-atom variable",
                   arg[action->argindex]);
    }
  }

  // warn if a keyword sets properties for atoms in rigid bodies
  // assume no conflict for properties not in this list of cases

  for (int i = 0; i < naction; i++) {
    switch (actions[i].keyword) {
    case X:
    case Y:
    case Z:
    case MOLECULE:
    case MASS:
    case ANGMOM:
    case SHAPE:
    case DIAMETER:
    case DENSITY:
    case TEMPERATURE:
    case QUAT:
    case IMAGE:
      if (modify->check_rigid_list_overlap(select))
        error->warning(FLERR,"Setting a property of atoms in rigid bodies "
                       "that has no effect unless rigid bodies are re-initialized");
      break;
    default:
      break;
    }
  }

  // in future, could return index to next arg to process
  // i.e. if fix set command appends its own options

  // return iarg;
}

/* ----------------------------------------------------------------------
   select atoms according to ATOM, MOLECULE, TYPE, GROUP, REGION style
   n = nlocal or nlocal+nghost depending on keyword
------------------------------------------------------------------------- */

void Set::selection(int n)
{
  // reallocate select vector if needed
  // this method could be called many times by fix set command

  if (n > maxselect) {
    memory->destroy(select);
    memory->create(select,n,"set:select");
    maxselect = n;
  }

  if (style == ATOM_SELECT) {
    tagint *tag = atom->tag;
    for (int i = 0; i < n; i++)
      if (tag[i] >= nlobig && tag[i] <= nhibig) select[i] = 1;
      else select[i] = 0;

  } else if (style == MOL_SELECT) {
    tagint *molecule = atom->molecule;
    for (int i = 0; i < n; i++)
      if (molecule[i] >= nlobig && molecule[i] <= nhibig) select[i] = 1;
      else select[i] = 0;

  } else if (style == TYPE_SELECT) {
    int *type = atom->type;
    for (int i = 0; i < n; i++)
      if (type[i] >= nlo && type[i] <= nhi) select[i] = 1;
      else select[i] = 0;

  } else if (style == GROUP_SELECT) {
    int *mask = atom->mask;
    for (int i = 0; i < n; i++)
      if (mask[i] & groupbit) select[i] = 1;
      else select[i] = 0;

  } else if (style == REGION_SELECT) {
    region->prematch();
    double **x = atom->x;
    for (int i = 0; i < n; i++)
      if (region->match(x[i][0],x[i][1],x[i][2])) select[i] = 1;
      else select[i] = 0;
  }

  // count_select = count of selected owned atoms

  count_select = 0;
  for (int i = 0; i < n; i++)
    if (select[i]) count_select++;
}

/* ----------------------------------------------------------------------
   loop over list of actions
   perform each action on all selected atoms via call to invoke_choice() method
------------------------------------------------------------------------- */

void Set::invoke_actions()
{
  // reallocate per-atom variable storage if needed

  if (varflag && atom->nlocal > maxvariable) {
    maxvariable = atom->nlocal;
    if (varflag1) {
      memory->destroy(vec1);
      memory->create(vec1,maxvariable,"set:var1");
    }
    if (varflag2) {
      memory->destroy(vec2);
      memory->create(vec2,maxvariable,"set:var2");
    }
    if (varflag3) {
      memory->destroy(vec3);
      memory->create(vec3,maxvariable,"set:var3");
    }
    if (varflag4) {
      memory->destroy(vec4);
      memory->create(vec4,maxvariable,"set:var4");
    }
  }

  // loop over actions

  for (int i = 0; i < naction; i++) {

    Action *action = &actions[i];

    // use count_action to optionally override count_select
    // if stays -1, count_select is used by caller
    // if overwritten by an invoke method, count_action is used
    // only a handful of invoke methods tally their own count

    count_action = -1;

    // evaluate atom-style variable(s) if necessary

    if (action->varflag) {
      if (action->varflag1) {
        input->variable->compute_atom(action->ivar1,0,vec1,1,0);
      }
      if (action->varflag2) {
        input->variable->compute_atom(action->ivar2,0,vec2,1,0);
      }
      if (action->varflag3) {
        input->variable->compute_atom(action->ivar3,0,vec3,1,0);
      }
      if (action->varflag4) {
        input->variable->compute_atom(action->ivar4,0,vec4,1,0);
      }
    }

    // invoke the action to reset per-atom or per-topology values

    (this->*invoke_choice[i])(action);

    action->count_select = count_select;
    action->count_action = count_action;
  }
}

/* ---------------------------------------------------------------------- */

void Set::varparse(const char *name, int m, Action *action)
{
  int ivar = input->variable->find(name+2);
  if (ivar < 0)
    error->all(FLERR,"Variable name {} for set command does not exist", name);
  if (!input->variable->atomstyle(ivar))
    error->all(FLERR,"Variable {} for set command is invalid style", name);

  action->varflag = 1;

  if (m == 1) {
    action->varflag1 = 1; action->ivar1 = ivar;
  } else if (m == 2) {
    action->varflag2 = 1; action->ivar2 = ivar;
  } else if (m == 3) {
    action->varflag3 = 1; action->ivar3 = ivar;
  } else if (m == 4) {
    action->varflag4 = 1; action->ivar4 = ivar;
  }
}

/* ----------------------------------------------------------------------
   set an owned atom property randomly
   set seed based on atom coordinates
   make atom result independent of what proc owns it
------------------------------------------------------------------------- */

void Set::setrandom(int keyword, Action *action)
{
  int i;

  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  auto avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
  auto avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
  auto avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));

  double **x = atom->x;

  // seed is always set to ivalue1 in process() methods

  int seed = action->ivalue1;

  auto ranpark = new RanPark(lmp,1);
  auto ranmars = new RanMars(lmp,seed + comm->me);

  // set approx fraction of atom types to newtype

  if (keyword == TYPE_FRACTION) {
    int nlocal = atom->nlocal;
    double fraction = action->dvalue1;
    int newtype = action->ivalue2;
    int count = 0;

    for (i = 0; i < nlocal; i++)
      if (select[i]) {
        ranpark->reset(seed,x[i]);
        if (ranpark->uniform() > fraction) continue;
        atom->type[i] = newtype;
        count++;
      }

    count_action = count;

  // set exact count of atom types to newtype
  // for TYPE_RATIO, exact = fraction out of total eligible
  // for TYPE_SUBSET, exact = nsubset out of total eligible

  } else if (keyword == TYPE_RATIO || keyword == TYPE_SUBSET) {
    int nlocal = atom->nlocal;
    int newtype = action->ivalue2;

    // convert specified fraction to nsubset of all selected atoms

    bigint bcount = count_select;
    bigint allcount;
    MPI_Allreduce(&bcount,&allcount,1,MPI_LMP_BIGINT,MPI_SUM,world);

    bigint nsubset;
    if (keyword == TYPE_RATIO) {
      double fraction = action->dvalue1;
      nsubset = static_cast<bigint> (fraction * allcount);
    } else if (keyword == TYPE_SUBSET) {
      nsubset = action->bvalue1;
      if (nsubset > allcount)
        error->all(FLERR,"Set type/subset value exceeds eligible atoms");
    }

    // make selection

    int *flag = memory->create(flag,count_select,"set:flag");
    int *work = memory->create(work,count_select,"set:work");

    ranmars->select_subset(nsubset,count_select,flag,work);

    // change types of selected atoms
    // flag vector from select_subset() is only for eligible atoms

    int count = 0;
    int eligible = 0;

    for (i = 0; i < nlocal; i++) {
      if (!select[i]) continue;
      if (flag[eligible]) {
        atom->type[i] = newtype;
        count++;
      }
      eligible++;
    }

    count_action = count;

    // clean up

    memory->destroy(flag);
    memory->destroy(work);

  // set dipole moments to random orientations in 3d or 2d
  // dipole length is determined by dipole type array

  } else if (keyword == DIPOLE_RANDOM) {
    double **mu = atom->mu;
    int nlocal = atom->nlocal;
    double dmag = action->dvalue1;
    double msq,scale;

    if (domain->dimension == 3) {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          ranpark->reset(seed,x[i]);
          mu[i][0] = ranpark->uniform() - 0.5;
          mu[i][1] = ranpark->uniform() - 0.5;
          mu[i][2] = ranpark->uniform() - 0.5;
          msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2];
          scale = dmag/sqrt(msq);
          mu[i][0] *= scale;
          mu[i][1] *= scale;
          mu[i][2] *= scale;
          mu[i][3] = dmag;
        }

    } else {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          ranpark->reset(seed,x[i]);
          mu[i][0] = ranpark->uniform() - 0.5;
          mu[i][1] = ranpark->uniform() - 0.5;
          mu[i][2] = 0.0;
          msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1];
          scale = dmag/sqrt(msq);
          mu[i][0] *= scale;
          mu[i][1] *= scale;
          mu[i][3] = dmag;
        }
    }

  // set spin moments to random orientations in 3d or 2d
  // spin length is fixed to unity

  } else if (keyword == SPIN_ATOM_RANDOM) {
    double **sp = atom->sp;
    int nlocal = atom->nlocal;
    double dlen = action->dvalue1;
    double sp_sq,scale;

    if (domain->dimension == 3) {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          ranpark->reset(seed,x[i]);
          sp[i][0] = ranpark->uniform() - 0.5;
          sp[i][1] = ranpark->uniform() - 0.5;
          sp[i][2] = ranpark->uniform() - 0.5;
          sp_sq = sp[i][0]*sp[i][0] + sp[i][1]*sp[i][1] + sp[i][2]*sp[i][2];
          scale = 1.0/sqrt(sp_sq);
          sp[i][0] *= scale;
          sp[i][1] *= scale;
          sp[i][2] *= scale;
          sp[i][3] = dlen;
        }

    } else {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          ranpark->reset(seed,x[i]);
          sp[i][0] = ranpark->uniform() - 0.5;
          sp[i][1] = ranpark->uniform() - 0.5;
          sp[i][2] = 0.0;
          sp_sq = sp[i][0]*sp[i][0] + sp[i][1]*sp[i][1];
          scale = 1.0/sqrt(sp_sq);
          sp[i][0] *= scale;
          sp[i][1] *= scale;
          sp[i][3] = dlen;
        }
    }

  // set quaternions to random orientations in 3d and 2d

  } else if (keyword == QUAT_RANDOM) {
    int *ellipsoid = atom->ellipsoid;
    int *tri = atom->tri;
    int *body = atom->body;
    double **quat = atom->quat;
    int nlocal = atom->nlocal;
    int quat_flag = atom->quat_flag;
    double *quat_one;

    if (domain->dimension == 3) {
      double s,t1,t2,theta1,theta2;
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          if (avec_ellipsoid && ellipsoid[i] >= 0)
            quat_one = avec_ellipsoid->bonus[ellipsoid[i]].quat;
          else if (avec_tri && tri[i] >= 0)
            quat_one = avec_tri->bonus[tri[i]].quat;
          else if (avec_body && body[i] >= 0)
            quat_one = avec_body->bonus[body[i]].quat;
          else if (quat_flag)
            quat_one = quat[i];
          else
            error->one(FLERR,"Cannot set quaternion for atom that has none");

          ranpark->reset(seed,x[i]);
          s = ranpark->uniform();
          t1 = sqrt(1.0-s);
          t2 = sqrt(s);
          theta1 = 2.0*MY_PI*ranpark->uniform();
          theta2 = 2.0*MY_PI*ranpark->uniform();
          quat_one[0] = cos(theta2)*t2;
          quat_one[1] = sin(theta1)*t1;
          quat_one[2] = cos(theta1)*t1;
          quat_one[3] = sin(theta2)*t2;
        }

    } else {
      double theta2;
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          if (avec_ellipsoid && ellipsoid[i] >= 0)
            quat_one = avec_ellipsoid->bonus[ellipsoid[i]].quat;
          else if (avec_body && body[i] >= 0)
            quat_one = avec_body->bonus[body[i]].quat;
          else if (quat_flag)
            quat_one = quat[i];
          else
            error->one(FLERR,"Cannot set quaternion for atom that has none");

          ranpark->reset(seed,x[i]);
          theta2 = MY_PI*ranpark->uniform();
          quat_one[0] = cos(theta2);
          quat_one[1] = 0.0;
          quat_one[2] = 0.0;
          quat_one[3] = sin(theta2);
        }
    }

  // set theta to random orientation in 2d

  } else if (keyword == THETA_RANDOM) {
    int *line = atom->line;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (select[i]) {
        if (line[i] < 0)
          error->one(FLERR,"Cannot set theta for atom that is not a line");
        ranpark->reset(seed,x[i]);
        avec_line->bonus[line[i]].theta = MY_2PI*ranpark->uniform();
      }
    }
  }

  delete ranpark;
  delete ranmars;
}

/* ---------------------------------------------------------------------- */

void Set::topology(int keyword, Action *action)
{
  int m,atom1,atom2,atom3,atom4;

  // error check

  if (atom->molecular == Atom::TEMPLATE)
    error->all(FLERR,"Cannot set bond topology types for atom style template");

  // border swap to acquire ghost atom info
  // enforce PBC before in case atoms are outside box
  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0) utils::logmesg(lmp,"  system init for set ...\n");
  lmp->init();

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  // select both owned and ghost atoms

  selection(atom->nlocal + atom->nghost);

  int count = 0;

  // for BOND, each of 2 atoms must be in group

  if (keyword == BOND) {
    int *num_bond = atom->num_bond;
    int **bond_type = atom->bond_type;
    tagint **bond_atom = atom->bond_atom;
    int nlocal = atom->nlocal;

    int itype = action->ivalue1;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < num_bond[i]; m++) {
        atom1 = atom->map(bond_atom[i][m]);
        if (atom1 == -1) error->one(FLERR,"Bond atom missing in set command");
        if (select[i] && select[atom1]) {
          bond_type[i][m] = itype;
          count++;
        }
      }
  }

  // for ANGLE, each of 3 atoms must be in group

  if (keyword == ANGLE) {
    int *num_angle = atom->num_angle;
    int **angle_type = atom->angle_type;
    tagint **angle_atom1 = atom->angle_atom1;
    tagint **angle_atom2 = atom->angle_atom2;
    tagint **angle_atom3 = atom->angle_atom3;
    int nlocal = atom->nlocal;

    int itype = action->ivalue1;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < num_angle[i]; m++) {
        atom1 = atom->map(angle_atom1[i][m]);
        atom2 = atom->map(angle_atom2[i][m]);
        atom3 = atom->map(angle_atom3[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1)
          error->one(FLERR,"Angle atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3]) {
          angle_type[i][m] = itype;
          count++;
        }
      }
  }

  // for DIHEDRAL, each of 4 atoms must be in group

  if (keyword == DIHEDRAL) {
    int *num_dihedral = atom->num_dihedral;
    int **dihedral_type = atom->dihedral_type;
    tagint **dihedral_atom1 = atom->dihedral_atom1;
    tagint **dihedral_atom2 = atom->dihedral_atom2;
    tagint **dihedral_atom3 = atom->dihedral_atom3;
    tagint **dihedral_atom4 = atom->dihedral_atom4;
    int nlocal = atom->nlocal;

    int itype = action->ivalue1;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < num_dihedral[i]; m++) {
        atom1 = atom->map(dihedral_atom1[i][m]);
        atom2 = atom->map(dihedral_atom2[i][m]);
        atom3 = atom->map(dihedral_atom3[i][m]);
        atom4 = atom->map(dihedral_atom4[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Dihedral atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3] && select[atom4]) {
          dihedral_type[i][m] = itype;
          count++;
        }
      }
  }

  // for IMPROPER, each of 4 atoms must be in group

  if (keyword == IMPROPER) {
    int *num_improper = atom->num_improper;
    int **improper_type = atom->improper_type;
    tagint **improper_atom1 = atom->improper_atom1;
    tagint **improper_atom2 = atom->improper_atom2;
    tagint **improper_atom3 = atom->improper_atom3;
    tagint **improper_atom4 = atom->improper_atom4;
    int nlocal = atom->nlocal;

    int itype = action->ivalue1;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < num_improper[i]; m++) {
        atom1 = atom->map(improper_atom1[i][m]);
        atom2 = atom->map(improper_atom2[i][m]);
        atom3 = atom->map(improper_atom3[i][m]);
        atom4 = atom->map(improper_atom4[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Improper atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3] && select[atom4]) {
          improper_type[i][m] = itype;
          count++;
        }
      }
  }

  // set count_action for all topology actions

  count_action = count;
}

// ----------------------------------------------------------------------
// pairs of process/invoke methods for each keyword
// process method reads args, stores parameters in Action instance
// invoke method resets atoms properties using Action instance
// separate two operations so can be called by either set or fix set command
// ----------------------------------------------------------------------

/* ---------------------------------------------------------------------- */

void Set::process_angle(int &iarg, int narg, char **arg, Action *action)
{
  if (atom->avec->angles_allow == 0)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set angle", error);

  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::ANGLE,lmp);
  action->ivalue1 = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  if (action->ivalue1 <= 0 || action->ivalue1 > atom->nangletypes)
    error->all(FLERR,"Invalid angle type in set command");
  iarg += 2;
}

void Set::invoke_angle(Action *action)
{
  topology(ANGLE,action);
}

/* ---------------------------------------------------------------------- */

void Set::process_angmom(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->angmom_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set angmom", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2,action);
  else action->dvalue2 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3,action);
  else action->dvalue3 = utils::numeric(FLERR,arg[iarg+3],false,lmp);

  iarg += 4;
}

void Set::invoke_angmom(Action *action)
{
  int nlocal = atom->nlocal;
  double **angmom = atom->angmom;

  int varflag = action->varflag;
  double xvalue,yvalue,zvalue;
  if (!action->varflag1) xvalue = action->dvalue1;
  if (!action->varflag2) yvalue = action->dvalue2;
  if (!action->varflag3) zvalue = action->dvalue3;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      if (action->varflag1) xvalue = vec1[i];
      if (action->varflag1) yvalue = vec2[i];
      if (action->varflag1) zvalue = vec3[i];
    }

    angmom[i][0] = xvalue;
    angmom[i][1] = yvalue;
    angmom[i][2] = zvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_bond(int &iarg, int narg, char **arg, Action *action)
{
  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set bond", error);

  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::BOND,lmp);
  action->ivalue1 = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  if (action->ivalue1 <= 0 || action->ivalue1 > atom->nbondtypes)
    error->all(FLERR,"Invalid bond type in set command");

  iarg += 2;
}

void Set::invoke_bond(Action *action)
{
  topology(BOND,action);
}

/* ---------------------------------------------------------------------- */

void Set::process_cc(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->tdpd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "set cc", error);

  action->ivalue1 = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
  if (action->ivalue1 < 1) error->all(FLERR,"Invalid cc index in set command");

  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    if (action->dvalue1 < 0.0) error->all(FLERR,"Invalid cc value in set command");
  }

  iarg += 3;
}

void Set::invoke_cc(Action *action)
{
  int nlocal = atom->nlocal;
  double **cc = atom->cc;

  int cc_index = action->ivalue1 - 1;
  // NOTE: need to check if cc_index exceeds cc array allocation

  int varflag = action->varflag;
  double ccvalue;
  if (!action->varflag1) ccvalue = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      ccvalue = vec1[i];
      if (ccvalue < 0.0) error->all(FLERR,"Invalid cc value in set command");
    }

    cc[i][cc_index] = ccvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_charge(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->q_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set charge", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);

  iarg += 2;
}

void Set::invoke_charge(Action *action)
{
  int nlocal = atom->nlocal;
  double *q = atom->q;
  double *q_scaled = atom->q_scaled;
  double *epsilon = atom->epsilon;

  int varflag = action->varflag;
  double qvalue;
  if (!action->varflag1) qvalue = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) qvalue = vec1[i];
    q[i] = qvalue;

    // ensure scaled charges are consistent with new charge value

    if (epsilon) q_scaled[i] = qvalue / epsilon[i];
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_density(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->rmass_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set density", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 <= 0.0) error->all(FLERR,"Invalid density in set command");
  }

  action->ivalue1 = 0;
  if (strcmp(arg[iarg],"density/disc") == 0) {
    action->ivalue1 = 1;
    if (domain->dimension != 2) error->all(FLERR,"Set density/disc requires 2d simulation");
  }

  iarg += 2;
}

void Set::invoke_density(Action *action)
{
  int nlocal = atom->nlocal;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;

  int radius_flag = atom->radius_flag;
  int ellipsoid_flag = atom->ellipsoid_flag;
  int line_flag = atom->line_flag;
  int tri_flag = atom->tri_flag;

  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  auto avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
  auto avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));

  int varflag = action->varflag;
  double density;
  if (!action->varflag1) density = action->dvalue1;
  int discflag = action->ivalue1;

  // set rmass via density
  // if radius > 0.0, treat as sphere or disc
  // if shape > 0.0, treat as ellipsoid (or ellipse, when uncomment below)
  // if length > 0.0, treat as line
  // if area > 0.0, treat as tri
  // else set rmass to density directly

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      density = vec1[i];
      if (density <= 0.0) error->one(FLERR,"Invalid density in set command");
    }

    if (radius_flag && radius[i] > 0.0)
      if (discflag) rmass[i] = MY_PI*radius[i]*radius[i] * density;
      else rmass[i] = 4.0*MY_PI/3.0 * radius[i]*radius[i]*radius[i] * density;

    else if (ellipsoid_flag && ellipsoid[i] >= 0) {
      double *shape = avec_ellipsoid->bonus[ellipsoid[i]].shape;
      // could enable 2d ellipse (versus 3d ellipsoid) when time integration
      //   options (fix nve/asphere, fix nh/asphere) are also implemented
      // if (discflag)
      // atom->rmass[i] = MY_PI*shape[0]*shape[1] * dvalue;
      // else
      rmass[i] = 4.0*MY_PI/3.0 * shape[0]*shape[1]*shape[2] * density;

    } else if (line_flag && line[i] >= 0) {
      double length = avec_line->bonus[line[i]].length;
      rmass[i] = length * density;

    } else if (tri_flag && tri[i] >= 0) {
      double *c1 = avec_tri->bonus[tri[i]].c1;
      double *c2 = avec_tri->bonus[tri[i]].c2;
      double *c3 = avec_tri->bonus[tri[i]].c3;
      double c2mc1[3],c3mc1[3];
      MathExtra::sub3(c2,c1,c2mc1);
      MathExtra::sub3(c3,c1,c3mc1);
      double norm[3];
      MathExtra::cross3(c2mc1,c3mc1,norm);
      double area = 0.5 * MathExtra::len3(norm);
      rmass[i] = area * density;

    } else rmass[i] = density;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_diameter(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->radius_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set diameter", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 < 0.0) error->one(FLERR,"Invalid diameter in set command");
  }

  iarg += 2;
}

void Set::invoke_diameter(Action *action)
{
  int nlocal = atom->nlocal;
  double *radius = atom->radius;

  int varflag = action->varflag;
  double diam;
  if (!action->varflag1) diam = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      diam = vec1[i];
      if (diam < 0.0) error->one(FLERR,"Invalid diameter in set command");
    }

    radius[i] = 0.5 * diam;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_dihedral(int &iarg, int narg, char **arg, Action *action)
{
  if (atom->avec->dihedrals_allow == 0)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set dihedral", error);

  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::DIHEDRAL,lmp);
  action->ivalue1 = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  if (action->ivalue1 <= 0 || action->ivalue1 > atom->ndihedraltypes)
    error->all(FLERR,"Invalid dihedral type in set command");

  iarg += 2;
}

void Set::invoke_dihedral(Action *action)
{
  topology(DIHEDRAL,action);
}

/* ---------------------------------------------------------------------- */

void Set::process_dipole(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->mu_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set dipole", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2,action);
  else action->dvalue2 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3,action);
  else action->dvalue3 = utils::numeric(FLERR,arg[iarg+3],false,lmp);

  iarg += 4;
}

void Set::invoke_dipole(Action *action)
{
  int nlocal = atom->nlocal;
  double **mu = atom->mu;

  int varflag = action->varflag;
  double xvalue,yvalue,zvalue;
  if (!action->varflag1) xvalue = action->dvalue1;
  if (!action->varflag2) yvalue = action->dvalue2;
  if (!action->varflag3) zvalue = action->dvalue3;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      if (action->varflag1) xvalue = vec1[i];
      if (action->varflag1) yvalue = vec2[i];
      if (action->varflag1) zvalue = vec3[i];
    }

    mu[i][0] = xvalue;
    mu[i][1] = yvalue;
    mu[i][2] = zvalue;
    mu[i][3] = sqrt(mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2]);
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_dipole_random(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->mu_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "set dipole/random", error);

  action->ivalue1 = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
  action->dvalue1 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (action->ivalue1 <= 0)
    error->all(FLERR,"Invalid random number seed in set command");
  if (action->dvalue1 <= 0.0)
    error->all(FLERR,"Invalid dipole length in set command");

  iarg += 3;
}

void Set::invoke_dipole_random(Action *action)
{
  setrandom(DIPOLE_RANDOM,action);
}

/* ---------------------------------------------------------------------- */

void Set::process_dpd_theta(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->dpd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set dpd/theta", error);

  if (strcmp(arg[iarg+1],"NULL") == 0) action->dvalue1 = -1.0;
  else if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 < 0.0) error->all(FLERR,"Invalid dpd/theta value in set command");
  }

  iarg += 2;
}

void Set::invoke_dpd_theta(Action *action)
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *dpdTheta = atom->dpdTheta;

  double tfactor = force->mvv2e / (domain->dimension * force->boltz);
  double onemass;
  double vx,vy,vz;

  int varflag = action->varflag;
  double theta;
  if (!action->varflag1) theta = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      theta = vec1[i];
      if (theta < 0.0) error->one(FLERR,"Invalid dpd/theta value in set command");
    }

    // if theta is negative, NULL was used, set dpdTheta to KE of particle

    if (theta >= 0.0) dpdTheta[i] = theta;
    else {
      if (rmass) onemass = rmass[i];
      else onemass = mass[type[i]];
      vx = v[i][0];
      vy = v[i][1];
      vz = v[i][2];
      dpdTheta[i] = tfactor * onemass * (vx*vx + vy*vy + vz*vz);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_edpd_cv(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->edpd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set edpd/cv", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 < 0.0) error->all(FLERR,"Invalid edpd/cv value in set command");
  }

  iarg += 2;
}

void Set::invoke_edpd_cv(Action *action)
{
  int nlocal = atom->nlocal;
  double *edpd_cv = atom->edpd_cv;

  int varflag = action->varflag;
  double cv;
  if (!action->varflag1) cv = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      cv = vec1[i];
      if (cv < 0.0) error->one(FLERR,"Invalid edpd/cv value in set command");
    }

    edpd_cv[i] = cv;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_edpd_temp(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->edpd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set edpd/temp", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 < 0.0) error->all(FLERR,"Invalid edpd/temp value in set command");
  }
  iarg += 2;
}

void Set::invoke_edpd_temp(Action *action)
{
  int nlocal = atom->nlocal;
  double *edpd_temp = atom->edpd_temp;

  int varflag = action->varflag;
  double temp;
  if (!action->varflag1) temp = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      temp = vec1[i];
      if (temp < 0.0) error->one(FLERR,"Invalid edpd/temp value in set command");
    }

    edpd_temp[i] = temp;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_epsilon(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->dielectric_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set epsilon", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 <= 0.0) error->all(FLERR,"Invalid epsilon in set command");
  }

  iarg += 2;
}

void Set::invoke_epsilon(Action *action)
{
  int nlocal = atom->nlocal;
  double *epsilon = atom->epsilon;
  double *q = atom->q;
  double *q_scaled = atom->q_scaled;

  int varflag = action->varflag;
  double eps;
  if (!action->varflag1) eps = action->dvalue1;

  // assign local dielectric constant
  // also update scaled charge value

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      eps = vec1[i];
      if (eps <= 0.0) error->one(FLERR,"Invalid epsilon in set command");
    }

    epsilon[i] = eps;
    q_scaled[i] = q[i] / eps;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_image(int &iarg, int narg, char **arg, Action *action)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set image", error);

  if (strcmp(arg[iarg+1],"NULL") == 0) action->ivalue4 = 0;
  else {
    action->ivalue4 = 1;
    if (utils::strmatch(arg[iarg+1],"^v_")) {
      if (!domain->xperiodic)
        error->all(FLERR,"Cannot set variable image flag for non-periodic dimension");
      varparse(arg[iarg+1],1,action);
    } else {
      action->ivalue1 = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (action->ivalue1 && !domain->xperiodic)
        error->all(FLERR,"Cannot set non-zero image flag for non-periodic dimension");
    }
  }

  if (strcmp(arg[iarg+2],"NULL") == 0) action->ivalue5 = 0;
  else {
    action->ivalue5 = 1;
    if (utils::strmatch(arg[iarg+2],"^v_")) {
      if (!domain->yperiodic)
        error->all(FLERR,"Cannot set variable image flag for non-periodic dimension");
      varparse(arg[iarg+2],2,action);
    } else {
      action->ivalue2 = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      if (action->ivalue2 && !domain->yperiodic)
        error->all(FLERR,"Cannot set non-zero image flag for non-periodic dimension");
    }
  }

  if (strcmp(arg[iarg+1],"NULL") == 0) action->ivalue6 = 0;
  else {
    action->ivalue6 = 1;
    if (utils::strmatch(arg[iarg+3],"^v_")) {
      if (!domain->zperiodic)
        error->all(FLERR,"Cannot set variable image flag for non-periodic dimension");
      varparse(arg[iarg+3],3,action);
    } else {
      action->ivalue3 = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      if (action->ivalue3 && !domain->zperiodic)
        error->all(FLERR,"Cannot set non-zero image flag for non-periodic dimension");
    }
  }

  iarg += 4;
}

void Set::invoke_image(Action *action)
{
  int nlocal = atom->nlocal;
  imageint *image = atom->image;
  int xbox,ybox,zbox;

  int ximageflag = action->ivalue4;
  int yimageflag = action->ivalue5;
  int zimageflag = action->ivalue6;

  int varflag = action->varflag;
  int ximage,yimage,zimage;
  if (!action->varflag1) ximage = action->ivalue1;
  if (!action->varflag2) yimage = action->ivalue2;
  if (!action->varflag3) zimage = action->ivalue3;

  // reset any or all of 3 image flags

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      if (action->varflag1) ximage = static_cast<int> (vec1[i]);
      if (action->varflag2) yimage = static_cast<int> (vec2[i]);
      if (action->varflag3) zimage = static_cast<int> (vec3[i]);
    }

    xbox = (image[i] & IMGMASK) - IMGMAX;
    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[i] >> IMG2BITS) - IMGMAX;
    if (ximageflag) xbox = ximage;
    if (yimageflag) ybox = yimage;
    if (zimageflag) zbox = zimage;
    image[i] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
      (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
      (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_improper(int &iarg, int narg, char **arg, Action *action)
{
  if (atom->avec->impropers_allow == 0)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set improper", error);

  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::IMPROPER,lmp);
  action->ivalue1 = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  if (action->ivalue1 <= 0 || action->ivalue1 > atom->nimpropertypes)
    error->all(FLERR,"Invalid value in set command");

  iarg += 2;
}

void Set::invoke_improper(Action *action)
{
  topology(IMPROPER,action);
}

/* ---------------------------------------------------------------------- */

void Set::process_length(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->line_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set length", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 < 0.0) error->one(FLERR,"Invalid length in set command");
  }

  iarg += 2;
}

void Set::invoke_length(Action *action)
{
  int nlocal = atom->nlocal;
  auto avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));

  int varflag = action->varflag;
  double length;
  if (!action->varflag1) length = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      length = vec1[i];
      if (length < 0.0) error->one(FLERR,"Invalid length in set command");
    }

    avec_line->set_length(i,length);
  }

  // update global line count

  bigint nlocal_bonus = avec_line->nlocal_bonus;
  MPI_Allreduce(&nlocal_bonus,&atom->nlines,1,MPI_LMP_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Set::process_mass(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->rmass_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set mass", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 <= 0.0) error->one(FLERR,"Invalid mass in set command");
  }

  iarg += 2;
}

void Set::invoke_mass(Action *action)
{
  int nlocal = atom->nlocal;
  double *rmass = atom->rmass;

  int varflag = action->varflag;
  double mass_one;
  if (!action->varflag1) mass_one = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      mass_one = vec1[i];
      if (mass_one < 0.0) error->one(FLERR,"Invalid mass in set command");
    }

    rmass[i] = mass_one;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_mol(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->molecule_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set mol", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->tvalue1 = utils::tnumeric(FLERR,arg[iarg+1],false,lmp);
    if (action->tvalue1 < 0) error->one(FLERR,"Invalid molecule ID in set command");
  }

  iarg += 2;
}

void Set::invoke_mol(Action *action)
{
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;

  int varflag = action->varflag;
  tagint molID;
  if (!action->varflag1) molID = action->tvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      molID = vec1[i];
      if (molID < 0) error->one(FLERR,"Invalid molecule ID in set command");
    }

    molecule[i] = molID;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_omega(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->omega_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set omega", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2,action);
  else action->dvalue2 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3,action);
  else action->dvalue3 = utils::numeric(FLERR,arg[iarg+3],false,lmp);

  iarg += 4;
}

void Set::invoke_omega(Action *action)
{
  int nlocal = atom->nlocal;
  double **omega = atom->omega;

  int varflag = action->varflag;
  double xvalue,yvalue,zvalue;
  if (!action->varflag1) xvalue = action->dvalue1;
  if (!action->varflag2) yvalue = action->dvalue2;
  if (!action->varflag3) zvalue = action->dvalue3;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      if (action->varflag1) xvalue = vec1[i];
      if (action->varflag2) yvalue = vec2[i];
      if (action->varflag3) zvalue = vec3[i];
    }

    omega[i][0] = xvalue;
    omega[i][1] = yvalue;
    omega[i][2] = zvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_quat(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->ellipsoid_flag && !atom->tri_flag && !atom->body_flag && !atom->quat_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+5 > narg) utils::missing_cmd_args(FLERR, "set quat", error);
  int dimension = domain->dimension;

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (dimension == 2 && action->dvalue1 != 0.0)
      error->one(FLERR,"Cannot set quaternion with xy components for 2d system");
  }
  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2,action);
  else {
    action->dvalue2 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    if (dimension == 2 && action->dvalue2 != 0.0)
      error->one(FLERR,"Cannot set quaternion with xy components for 2d system");
  }
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3,action);
  else action->dvalue3 = utils::numeric(FLERR,arg[iarg+3],false,lmp);
  if (utils::strmatch(arg[iarg+4],"^v_")) varparse(arg[iarg+4],4,action);
  else action->dvalue4 = utils::numeric(FLERR,arg[iarg+4],false,lmp);

  iarg += 5;
}

void Set::invoke_quat(Action *action)
{
  int nlocal = atom->nlocal;
  int *ellipsoid = atom->ellipsoid;
  int *tri = atom->tri;
  int *body = atom->body;
  double **quat = atom->quat;
  int quat_flag = atom->quat_flag;

  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  auto avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
  auto avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));

  int dimension = domain->dimension;
  double radians,sintheta;
  double *quat_one;

  int varflag = action->varflag;
  double xvalue,yvalue,zvalue,theta;
  if (!action->varflag1) xvalue = action->dvalue1;
  if (!action->varflag2) yvalue = action->dvalue2;
  if (!action->varflag3) zvalue = action->dvalue3;
  if (!action->varflag4) theta = action->dvalue4;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (avec_ellipsoid && ellipsoid[i] >= 0)
      quat_one = avec_ellipsoid->bonus[ellipsoid[i]].quat;
    else if (avec_tri && tri[i] >= 0)
      quat_one = avec_tri->bonus[tri[i]].quat;
    else if (avec_body && body[i] >= 0)
      quat_one = avec_body->bonus[body[i]].quat;
    else if (quat_flag)
      quat_one = quat[i];
    else
      error->one(FLERR,"Cannot set quaternion for atom that has none");

    if (varflag) {
      if (action->varflag1) xvalue = vec1[i];
      if (action->varflag2) yvalue = vec2[i];
      if (action->varflag3) zvalue = vec3[i];
      if (action->varflag4) theta = vec4[i];
      if (dimension == 2 && (xvalue != 0.0 || yvalue != 0.0))
        error->one(FLERR,"Cannot set quaternion with xy components for 2d system");
    }

    radians = MY_PI2 * theta/180.0;
    sintheta = sin(radians);
    quat_one[0] = cos(radians);
    quat_one[1] = xvalue * sintheta;
    quat_one[2] = yvalue * sintheta;
    quat_one[3] = zvalue * sintheta;
    MathExtra::qnormalize(quat_one);
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_quat_random(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->ellipsoid_flag && !atom->tri_flag && !atom->body_flag && !atom->quat_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set quat/random", error);

  action->ivalue1 = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
  if (action->ivalue1 <= 0) error->all(FLERR,"Invalid random number seed in set command");

  iarg += 2;
}

void Set::invoke_quat_random(Action *action)
{
  setrandom(QUAT_RANDOM,action);
}

/* ---------------------------------------------------------------------- */

void Set::process_radius_election(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->eradius_flag)
    error->all(FLERR, "Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "set radius/electron", error);

  if (utils::strmatch(arg[iarg + 1], "^v_"))
    varparse(arg[iarg + 1], 1, action);
  else {
    action->dvalue1 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    if (action->dvalue1 < 0.0)
      error->one(FLERR, iarg + 1, "Invalid electron radius {} in set radius/electron command",
                 action->dvalue1);
  }

  iarg += 2;
}

void Set::invoke_radius_election(Action *action)
{
  int nlocal = atom->nlocal;
  double *eradius = atom->eradius;

  int varflag = action->varflag;
  double radius;
  if (!action->varflag1) radius = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      radius = vec1[i];
      if (radius < 0.0) error->one(FLERR,"Invalid electron radius in set command");
    }

    eradius[i] = radius;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_shape(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->ellipsoid_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set shape", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 < 0.0) error->one(FLERR,"Invalid shape in set command");
  }
  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2,action);
  else {
    action->dvalue2 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    if (action->dvalue2 < 0.0) error->one(FLERR,"Invalid shape in set command");
  }
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3,action);
  else {
    action->dvalue3 = utils::numeric(FLERR,arg[iarg+3],false,lmp);
    if (action->dvalue3 < 0.0) error->one(FLERR,"Invalid shape in set command");
  }

  iarg += 4;
}

void Set::invoke_shape(Action *action)
{
  int nlocal = atom->nlocal;
  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));

  int varflag = action->varflag;
  double xvalue,yvalue,zvalue;
  if (!action->varflag1) xvalue = action->dvalue1;
  if (!action->varflag2) yvalue = action->dvalue2;
  if (!action->varflag3) zvalue = action->dvalue3;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      if (action->varflag1) xvalue = vec1[i];
      if (action->varflag1) yvalue = vec2[i];
      if (action->varflag1) zvalue = vec3[i];
      if (xvalue < 0.0 || yvalue < 0.0 || zvalue < 0.0)
        error->one(FLERR,"Invalid shape in set command");
    }

    if (xvalue > 0.0 || yvalue > 0.0 || zvalue > 0.0)
      if (xvalue == 0.0 || yvalue == 0.0 || zvalue == 0.0)
        error->one(FLERR,"Invalid shape in set command");

    avec_ellipsoid->set_shape(i,0.5*xvalue,0.5*yvalue,0.5*zvalue);
  }

  // update global ellipsoid count

  bigint nlocal_bonus = avec_ellipsoid->nlocal_bonus;
  MPI_Allreduce(&nlocal_bonus,&atom->nellipsoids,1,MPI_LMP_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Set::process_smd_contact_radius(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->smd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set smd/contact/radius", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 < 0.0) error->one(FLERR,"Invalid smd/contact/radius in set command");
  }

  iarg += 2;
}

void Set::invoke_smd_contact_radius(Action *action)
{
  int nlocal = atom->nlocal;
  double *contact_radius = atom->contact_radius;

  int varflag = action->varflag;
  double radius;
  if (!action->varflag1) radius = action->dvalue1;


  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      radius = vec1[i];
      if (radius < 0.0) error->one(FLERR,"Invalid smd/contact/radius in set command");
    }

    contact_radius[i] = radius;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_smd_mass_density(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->smd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set smd/mass/density", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 <= 0.0) error->one(FLERR,"Invalid smd/mass/density in set command");
  }

  iarg += 2;
}

void Set::invoke_smd_mass_density(Action *action)
{
  int nlocal = atom->nlocal;
  double *rmass = atom->rmass;
  double *vfrac = atom->vfrac;

  int varflag = action->varflag;
  double density;
  if (!action->varflag1) density = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      density = vec1[i];
      if (density < 0.0) error->one(FLERR,"Invalid smd/mass/density in set command");
    }

    rmass[i] = vfrac[i] * density;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_sph_cv(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->cv_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set sph/cv", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);

  iarg += 2;
}

void Set::invoke_sph_cv(Action *action)
{
  int nlocal = atom->nlocal;
  double *cv = atom->cv;

  int varflag = action->varflag;
  double sph_cv;
  if (!action->varflag1) sph_cv = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (varflag) sph_cv = vec1[i];
    cv[i] = sph_cv;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_sph_e(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->esph_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set sph/e", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);

  iarg += 2;
}

void Set::invoke_sph_e(Action *action)
{
  int nlocal = atom->nlocal;
  double *esph = atom->esph;

  int varflag = action->varflag;
  double sph_e;
  if (!action->varflag1) sph_e = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (varflag) sph_e = vec1[i];
    esph[i] = sph_e;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_sph_rho(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->rho_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set sph/rho", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);

  iarg += 2;
}

void Set::invoke_sph_rho(Action *action)
{
  int nlocal = atom->nlocal;
  double *rho = atom->rho;

  int varflag = action->varflag;
  double sph_rho;
  if (!action->varflag1) sph_rho = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (varflag) sph_rho = vec1[i];
    rho[i] = sph_rho;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_spin_atom(int &iarg, int narg, char **arg, Action *action)
{
  if ((strcmp(arg[iarg],"spin") == 0) && (comm->me == 0))
    error->warning(FLERR, "Set attribute spin is deprecated -- use spin/atom instead");
  if (!atom->sp_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+5 > narg) utils::missing_cmd_args(FLERR, "set spin/atom", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 <= 0.0)
      error->all(FLERR,"Invalid spin magnitude {} in set {} command", action->dvalue1, arg[iarg]);
  }

  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2,action);
  else action->dvalue2 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3,action);
  else action->dvalue3 = utils::numeric(FLERR,arg[iarg+3],false,lmp);
  if (utils::strmatch(arg[iarg+4],"^v_")) varparse(arg[iarg+4],4,action);
  else action->dvalue4 = utils::numeric(FLERR,arg[iarg+4],false,lmp);

  iarg += 5;
}

void Set::invoke_spin_atom(Action *action)
{
  int nlocal = atom->nlocal;
  double **sp = atom->sp;
  double norm;

  int varflag = action->varflag;
  double magnitude,xvalue,yvalue,zvalue;
  if (!action->varflag1) magnitude = action->dvalue1;
  if (!action->varflag2) xvalue = action->dvalue2;
  if (!action->varflag3) yvalue = action->dvalue3;
  if (!action->varflag4) zvalue = action->dvalue4;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      if (action->varflag1) magnitude = vec1[i];
      if (magnitude < 0.0)
        error->one(FLERR,"Invalid spin magnitude in set command");
      if (action->varflag2) xvalue = vec2[i];
      if (action->varflag3) yvalue = vec3[i];
      if (action->varflag4) zvalue = vec4[i];
    }

    if ((xvalue == 0.0) && (yvalue == 0.0) && (zvalue == 0.0))
      error->all(FLERR,"At least one spin vector component must be non-zero");

    norm = 1.0/sqrt(xvalue*xvalue+yvalue*yvalue+zvalue*zvalue);
    sp[i][0] = norm*xvalue;
    sp[i][1] = norm*yvalue;
    sp[i][2] = norm*zvalue;
    sp[i][3] = magnitude;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_spin_atom_random(int &iarg, int narg, char **arg, Action *action)
{
  if ((strcmp(arg[iarg], "spin/random") == 0) && (comm->me == 0))
    error->warning(FLERR,
                   "Set attribute spin/random is deprecated -- use spin/atom/random instead");
  if (!atom->sp_flag)
    error->all(FLERR, "Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "set spin/atom/random", error);

  action->ivalue1 = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
  action->dvalue1 = utils::numeric(FLERR, arg[iarg + 2], false, lmp);

  if (action->ivalue1 <= 0)
    error->all(FLERR, iarg + 1, "Invalid random number seed {} in set {} command", action->ivalue1,
               arg[iarg]);
  if (action->dvalue1 <= 0.0)
    error->all(FLERR, iarg + 2, "Invalid spin magnitude {} in set {} command", action->dvalue1,
               arg[iarg]);
  iarg += 3;
}

void Set::invoke_spin_atom_random(Action *action)
{
  setrandom(SPIN_ATOM_RANDOM, action);
}

/* ---------------------------------------------------------------------- */

void Set::process_spin_electron(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->spin_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set spin/electron", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->ivalue1 = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
    if (action->ivalue1 < -1 || action->ivalue1 > 3)
      error->one(FLERR,"Invalid electron spin {} in set command", action->ivalue1);
  }

  iarg += 2;
}

void Set::invoke_spin_electron(Action *action)
{
  int nlocal = atom->nlocal;
  int *spin = atom->spin;

  int varflag = action->varflag;
  int ispin;
  if (!action->varflag1) ispin = action->ivalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      ispin = static_cast<int> (vec1[i]);
      if (ispin < -1 || ispin > 3)
        error->one(FLERR,"Invalid electron spin in set command");
    }

    atom->spin[i] = ispin;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_temperature(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->temperature_flag)
    error->all(FLERR,"Cannot set this attribute for this atom style");
  if (iarg+2 > narg) error->all(FLERR,"Illegal set command");

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 < 0.0) error->one(FLERR,"Invalid temperature in set command");
  }

  iarg += 2;
}

void Set::invoke_temperature(Action *action)
{
  int nlocal = atom->nlocal;
  double *temperature = atom->temperature;

  int varflag = action->varflag;
  double temp;
  if (!action->varflag1) temp = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      temp = vec1[i];
      if (temp < 0.0) error->one(FLERR,"Invalid temperature in set command");
    }

    temperature[i] = temp;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_theta(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->line_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set theta", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = DEG2RAD * utils::numeric(FLERR,arg[iarg+1],false,lmp);

  iarg += 2;
}

void Set::invoke_theta(Action *action)
{
  int nlocal = atom->nlocal;
  int *line = atom->line;

  auto avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));

  int varflag = action->varflag;
  double theta;
  if (!action->varflag1) theta = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (line[i] < 0) error->one(FLERR,"Cannot set theta for atom which is not a line");
    if (varflag) theta = vec1[i];
    avec_line->bonus[atom->line[i]].theta = theta;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_theta_random(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->line_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set theta/random", error);

  action->ivalue1 = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
  if (action->ivalue1 <= 0) error->all(FLERR,"Invalid random number seed in set command");

  iarg += 2;
}

void Set::invoke_theta_random(Action *action)
{
  setrandom(THETA_RANDOM,action);
}

/* ---------------------------------------------------------------------- */

void Set::process_tri(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->tri_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set tri", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 < 0.0) error->one(FLERR,"Invalid tri size in set command");
  }

  iarg += 2;
}

void Set::invoke_tri(Action *action)
{
  int nlocal = atom->nlocal;
  int *tri = atom->tri;

  auto avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));

  int varflag = action->varflag;
  double trisize;
  if (!action->varflag1) trisize = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
#if 0
    // AK: this seems wrong. Isn't the set command supposed *make* this a triangle?
    if (tri[i] < 0) error->one(FLERR,"Cannot set tri for atom which is not a triangle");
#endif

    if (varflag) {
      trisize = vec1[i];
      if (trisize < 0.0) error->one(FLERR,"Invalid tri size in set command");
    }

    avec_tri->set_equilateral(i,trisize);
  }

  // update bonus tri count

  bigint nlocal_bonus = avec_tri->nlocal_bonus;
  MPI_Allreduce(&nlocal_bonus,&atom->ntris,1,MPI_LMP_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Set::process_type(int &iarg, int narg, char **arg, Action *action)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set type", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::ATOM,lmp);
    action->ivalue1 = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
    delete[] typestr;
    if (action->ivalue1 <= 0 || action->ivalue1 > atom->ntypes)
      error->one(FLERR,"Invalid atom type in set command");
  }

  iarg += 2;
}

void Set::invoke_type(Action *action)
{
  int nlocal = atom->nlocal;
  int *type = atom->type;

  int varflag = action->varflag;
  int itype;
  if (!action->varflag1) itype = action->ivalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (action->varflag) {
      itype = static_cast<int> (vec1[i]);
      if (itype <= 0 || itype > atom->ntypes)
        error->one(FLERR,"Invalid atom type in set command");
    }

    type[i] = itype;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_type_fraction(int &iarg, int narg, char **arg, Action *action)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set type/fraction", error);

  // random seed must be ivalue1 for use in setrandom()

  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::ATOM,lmp);
  action->ivalue2 = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  if (action->ivalue2 <= 0 || action->ivalue2 > atom->ntypes)
      error->one(FLERR,"Invalid atom type in set command");

  action->dvalue1 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (action->dvalue1 < 0.0 || action->dvalue1 > 1.0)
    error->all(FLERR,"Invalid fraction in set command");

  action->ivalue1 = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
  if (action->ivalue1 <= 0)
    error->all(FLERR,"Invalid random number seed in set command");

  iarg += 4;
}

void Set::invoke_type_fraction(Action *action)
{
  setrandom(TYPE_FRACTION,action);
}

/* ---------------------------------------------------------------------- */

void Set::process_type_ratio(int &iarg, int narg, char **arg, Action *action)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set type/ratio", error);

  // random seed must be ivalue1 for use in setrandom()

  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::ATOM,lmp);
  action->ivalue2 = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  if (action->ivalue2 <= 0 || action->ivalue2 > atom->ntypes)
    error->all(FLERR,"Invalid atom type in set command");

  action->dvalue1 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (action->dvalue1 < 0.0 || action->dvalue1 > 1.0)
    error->all(FLERR,"Invalid fraction in set command");

  action->ivalue1 = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
  if (action->ivalue1 <= 0)
    error->all(FLERR,"Invalid random number seed in set command");

  iarg += 4;
}

void Set::invoke_type_ratio(Action *action)
{
  setrandom(TYPE_RATIO,action);
}

/* ---------------------------------------------------------------------- */

void Set::process_type_subset(int &iarg, int narg, char **arg, Action *action)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set type/subset", error);

  // random seed must be ivalue1 for use in setrandom()

  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::ATOM,lmp);
  action->ivalue2 = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  if (action->ivalue2 <= 0 || action->ivalue2 > atom->ntypes)
    error->all(FLERR,"Invalid atom type in set command");

  action->bvalue1 = utils::bnumeric(FLERR,arg[iarg+2],false,lmp);
  if (action->bvalue1 < 0)
    error->all(FLERR,"Invalid subset size in set command");

  action->ivalue1 = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
  if (action->ivalue1 <= 0)
    error->all(FLERR,"Invalid random number seed in set command");

  iarg += 4;
}

void Set::invoke_type_subset(Action *action)
{
  setrandom(TYPE_SUBSET,action);
}

/* ---------------------------------------------------------------------- */

void Set::process_volume(int &iarg, int narg, char **arg, Action *action)
{
  if (!atom->vfrac_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set volume", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else {
    action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (action->dvalue1 <= 0.0) error->all(FLERR,"Invalid volume in set command");
  }

  iarg += 2;
}

void Set::invoke_volume(Action *action)
{
  int nlocal = atom->nlocal;
  double *vfrac = atom->vfrac;

  int varflag = action->varflag;
  double vol;
  if (!action->varflag1) vol = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (varflag) {
      vol = vec1[i];
      if (vol < 0.0) error->one(FLERR,"Invalid volume in set command");
    }

    vfrac[i] = vol;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_vx(int &iarg, int narg, char **arg, Action *action)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set vx", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);

  iarg += 2;
}

void Set::invoke_vx(Action *action)
{
  int nlocal = atom->nlocal;
  double **v = atom->v;

  int varflag = action->varflag;
  double vx;
  if (!action->varflag1) vx = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (varflag) vx = vec1[i];
    v[i][0] = vx;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_vy(int &iarg, int narg, char **arg, Action *action)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set vy", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);

  iarg += 2;
}

void Set::invoke_vy(Action *action)
{
  int nlocal = atom->nlocal;
  double **v = atom->v;

  int varflag = action->varflag;
  double vy;
  if (!action->varflag1) vy = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (varflag) vy = vec1[i];
    v[i][1] = vy;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_vz(int &iarg, int narg, char **arg, Action *action)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set vz", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);

  iarg += 2;
}

void Set::invoke_vz(Action *action)
{
  int nlocal = atom->nlocal;
  double **v = atom->v;

  int varflag = action->varflag;
  double vz;
  if (!action->varflag1) vz = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (varflag) vz = vec1[i];
    v[i][2] = vz;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_x(int &iarg, int narg, char **arg, Action *action)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set x", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);

  iarg += 2;
}

void Set::invoke_x(Action *action)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;

  int varflag = action->varflag;
  double coord;
  if (!action->varflag1) coord = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (varflag) coord = vec1[i];
    x[i][0] = coord;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_y(int &iarg, int narg, char **arg, Action *action)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set y", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);

  iarg += 2;
}

void Set::invoke_y(Action *action)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;

  int varflag = action->varflag;
  double coord;
  if (!action->varflag1) coord = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (varflag) coord = vec1[i];
    x[i][1] = coord;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_z(int &iarg, int narg, char **arg, Action *action)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set z", error);

  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
  else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);

  iarg += 2;
}

void Set::invoke_z(Action *action)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;

  int varflag = action->varflag;
  double coord;
  if (!action->varflag1) coord = action->dvalue1;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (varflag) coord = vec1[i];
    x[i][2] = coord;
  }
}

/* ---------------------------------------------------------------------- */

void Set::process_custom(int &iarg, int narg, char **arg, Action *action)
{
  int flag,cols;
  ArgInfo argi(arg[iarg],ArgInfo::DNAME|ArgInfo::INAME);
  const char *pname = argi.get_name();

  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set", error);
  int index_custom = atom->find_custom(argi.get_name(),flag,cols);
  if (index_custom < 0)
    error->all(FLERR,"Set keyword or custom property {} does not exist",pname);
  action->ivalue2 = index_custom;

  switch (argi.get_type()) {

  case ArgInfo::INAME:
    if (flag != 0) error->all(FLERR,"Set command custom property {} is not integer",pname);
    if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
    else action->ivalue1 = utils::inumeric(FLERR,arg[iarg+1],false,lmp);

    if (argi.get_dim() == 0) {
      if (cols > 0)
        error->all(FLERR,"Set command custom integer property {} is not a vector",pname);
      action->keyword = IVEC;
    } else if (argi.get_dim() == 1) {
      if (cols == 0)
        error->all(FLERR,"Set command custom integer property {} is not an array",pname);
      int icol_custom = argi.get_index1();
      if (icol_custom <= 0 || icol_custom > cols)
        error->all(FLERR,"Set command per-atom custom integer array {} is accessed "
                   "out-of-range",pname);
      action->ivalue3 = icol_custom;
      action->keyword = IARRAY;
    } else error->all(FLERR,"Illegal set command");
    break;

  case ArgInfo::DNAME:
    if (flag != 1) error->all(FLERR,"Custom property {} is not floating-point",argi.get_name());
    if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1,action);
    else action->dvalue1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);

    if (argi.get_dim() == 0) {
      if (cols > 0)
        error->all(FLERR,"Set command custom double property {} is not a vector",pname);
      action->keyword = DVEC;
    } else if (argi.get_dim() == 1) {
      if (cols == 0)
        error->all(FLERR,"Set command custom double property {} is not an array",pname);
      int icol_custom = argi.get_index1();
      if (icol_custom <= 0 || icol_custom > cols)
        error->all(FLERR,"Set command per-atom custom double array {} is "
                   "accessed out-of-range",pname);
      action->ivalue3 = icol_custom;
      action->keyword = DARRAY;
    } else error->all(FLERR,"Illegal set command");
    break;

  default:
    error->all(FLERR,"Illegal set command");
    break;
  }

  iarg += 2;
}

void Set::invoke_custom(Action *action)
{
  int nlocal = atom->nlocal;
  int ivalue;
  double dvalue;

  int varflag = action->varflag;
  int index_custom = action->ivalue2;

  if (action->keyword == IVEC) {
    if (!varflag) ivalue = action->ivalue1;
    int *ivector = atom->ivector[index_custom];
    for (int i = 0; i < nlocal; i++) {
      if (!select[i]) continue;
      if (varflag) ivalue = static_cast<int> (vec1[i]);
      ivector[i] = ivalue;
    }

  } else if (action->keyword == DVEC) {
    if (!varflag) dvalue = action->dvalue1;
    double *dvector = atom->dvector[index_custom];
    for (int i = 0; i < nlocal; i++) {
      if (!select[i]) continue;
      if (varflag) dvalue = vec1[i];
      dvector[i] = dvalue;
    }

  } else if (action->keyword == IARRAY) {
    if (!varflag) ivalue = action->ivalue1;
    int **iarray = atom->iarray[index_custom];
    int icol_custom = action->ivalue3 - 1;
    for (int i = 0; i < nlocal; i++) {
      if (!select[i]) continue;
      if (varflag) ivalue = static_cast<int> (vec1[i]);
      iarray[i][icol_custom] = ivalue;
    }

  } else if (action->keyword == DARRAY) {
    if (!varflag) dvalue = action->dvalue1;
    double **darray = atom->darray[index_custom];
    int icol_custom = action->ivalue3 - 1;
    for (int i = 0; i < nlocal; i++) {
      if (!select[i]) continue;
      if (varflag) dvalue = vec1[i];
      darray[i][icol_custom] = dvalue;
    }
  }
}
