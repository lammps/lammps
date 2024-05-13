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

#include "fix_shake.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_respa.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "respa.h"
#include "update.h"

#include <cctype>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

static constexpr int RVOUS = 1;    // 0 for irregular, 1 for all2all

static constexpr double BIG = 1.0e20;
static constexpr double MASSDELTA = 0.1;

/* ---------------------------------------------------------------------- */

FixShake::FixShake(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), bond_flag(nullptr), angle_flag(nullptr), type_flag(nullptr),
    mass_list(nullptr), bond_distance(nullptr), angle_distance(nullptr), loop_respa(nullptr),
    step_respa(nullptr), x(nullptr), v(nullptr), f(nullptr), ftmp(nullptr), vtmp(nullptr),
    mass(nullptr), rmass(nullptr), type(nullptr), shake_flag(nullptr), shake_atom(nullptr),
    shake_type(nullptr), xshake(nullptr), nshake(nullptr), list(nullptr), closest_list(nullptr),
    b_count(nullptr), b_count_all(nullptr), b_ave(nullptr), b_max(nullptr), b_min(nullptr),
    b_ave_all(nullptr), b_max_all(nullptr), b_min_all(nullptr), a_count(nullptr),
    a_count_all(nullptr), a_ave(nullptr), a_max(nullptr), a_min(nullptr), a_ave_all(nullptr),
    a_max_all(nullptr), a_min_all(nullptr), atommols(nullptr), onemols(nullptr)
{
  energy_global_flag = energy_peratom_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  thermo_energy = thermo_virial = 1;
  create_attribute = 1;
  dof_flag = 1;
  scalar_flag = 1;
  extscalar = 1;
  stores_ids = 1;
  centroidstressflag = CENTROID_AVAIL;
  next_output = -1;

  // to avoid uninitialized access
  vflag_post_force = 0;
  eflag_pre_reverse = 0;
  ebond = 0.0;

  // error check

  molecular = atom->molecular;
  if (molecular == Atom::ATOMIC)
    error->all(FLERR, "Cannot use fix {} with non-molecular system", style);

  // perform initial allocation of atom-based arrays
  // register with Atom class

  shake_flag = nullptr;
  shake_atom = nullptr;
  shake_type = nullptr;
  xshake = nullptr;

  ftmp = nullptr;
  vtmp = nullptr;

  FixShake::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  // set comm size needed by this fix

  comm_forward = 3;

  // parse SHAKE args
  auto mystyle = fmt::format("fix {}", style);

  if (narg < 8) utils::missing_cmd_args(FLERR, mystyle, error);

  tolerance = utils::numeric(FLERR, arg[3], false, lmp);
  max_iter = utils::inumeric(FLERR, arg[4], false, lmp);
  output_every = utils::inumeric(FLERR, arg[5], false, lmp);

  // parse SHAKE args for bond and angle types
  // will be used by find_clusters
  // store args for "b" "a" "t" as flags in (1:n) list for fast access
  // store args for "m" in list of length nmass for looping over
  // for "m" verify that atom masses have been set

  bond_flag = new int[atom->nbondtypes + 1];
  for (int i = 1; i <= atom->nbondtypes; i++) bond_flag[i] = 0;
  angle_flag = new int[atom->nangletypes + 1];
  for (int i = 1; i <= atom->nangletypes; i++) angle_flag[i] = 0;
  type_flag = new int[atom->ntypes + 1];
  for (int i = 1; i <= atom->ntypes; i++) type_flag[i] = 0;
  mass_list = new double[atom->ntypes];
  nmass = 0;

  char mode = '\0';
  int next = 6;
  while (next < narg) {
    if (strcmp(arg[next],"b") == 0) mode = 'b';
    else if (strcmp(arg[next],"a") == 0) mode = 'a';
    else if (strcmp(arg[next],"t") == 0) mode = 't';
    else if (strcmp(arg[next],"m") == 0) {
      mode = 'm';
      atom->check_mass(FLERR);

    // break if keyword that is not b,a,t,m

    } else if (isalpha(arg[next][0])) break;

    // read numeric args of b,a,t,m

    else if (mode == 'b') {
      int i = utils::inumeric(FLERR,arg[next],false,lmp);
      if (i < 1 || i > atom->nbondtypes)
        error->all(FLERR,"Invalid bond type index for {}", mystyle);
      bond_flag[i] = 1;

    } else if (mode == 'a') {
      int i = utils::inumeric(FLERR,arg[next],false,lmp);
      if (i < 1 || i > atom->nangletypes)
        error->all(FLERR,"Invalid angle type index for {}", mystyle);
      angle_flag[i] = 1;

    } else if (mode == 't') {
      int i = utils::inumeric(FLERR,arg[next],false,lmp);
      if (i < 1 || i > atom->ntypes)
        error->all(FLERR,"Invalid atom type index for {}", mystyle);
      type_flag[i] = 1;

    } else if (mode == 'm') {
      double massone = utils::numeric(FLERR,arg[next],false,lmp);
      if (massone == 0.0) error->all(FLERR,"Invalid atom mass for {}", mystyle);
      if (nmass == atom->ntypes)
        error->all(FLERR,"Too many masses for {}", mystyle);
      mass_list[nmass++] = massone;

    } else error->all(FLERR,"Unknown {} command option: {}", mystyle, arg[next]);
    next++;
  }

  // parse optional args

  onemols = nullptr;
  kbond = 1.0e6*force->boltz;

  int iarg = next;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,mystyle+" mol",error);
      int imol = atom->find_molecule(arg[iarg+1]);
      if (imol == -1)
        error->all(FLERR,"Molecule template ID {} for {} does not exist", mystyle, arg[iarg+1]);
      if ((atom->molecules[imol]->nset > 1) && (comm->me == 0))
        error->warning(FLERR,"Molecule template for {} has multiple molecules", mystyle);
      onemols = &atom->molecules[imol];
      nmol = onemols[0]->nset;
      iarg += 2;
    } else if (strcmp(arg[iarg],"kbond") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,mystyle+" kbond",error);
      kbond = utils::numeric(FLERR, arg[iarg+1], false, lmp);
      if (kbond < 0) error->all(FLERR,"Illegal {} kbond value {}. Must be >= 0.0", mystyle, kbond);
      iarg += 2;
    } else error->all(FLERR,"Unknown {} command option: {}", mystyle, arg[iarg]);
  }

  // error check for Molecule template

  if (onemols) {
    for (int i = 0; i < nmol; i++)
      if (onemols[i]->shakeflag == 0)
        error->all(FLERR,"Fix {} molecule template must have shake info", style);
  }

  // allocate bond and angle distance arrays, indexed from 1 to n

  bond_distance = new double[atom->nbondtypes+1];
  angle_distance = new double[atom->nangletypes+1];

  // allocate statistics arrays

  if (output_every) {
    int nb = atom->nbondtypes + 1;
    b_count = new bigint[nb];
    b_count_all = new bigint[nb];
    b_ave = new double[nb];
    b_ave_all = new double[nb];
    b_max = new double[nb];
    b_max_all = new double[nb];
    b_min = new double[nb];
    b_min_all = new double[nb];

    int na = atom->nangletypes + 1;
    a_count = new bigint[na];
    a_count_all = new bigint[na];
    a_ave = new double[na];
    a_ave_all = new double[na];
    a_max = new double[na];
    a_max_all = new double[na];
    a_min = new double[na];
    a_min_all = new double[na];
  }

  // SHAKE vs RATTLE

  rattle = 0;

  // identify all SHAKE clusters

  double time1 = platform::walltime();

  find_clusters();

  if (comm->me == 0)
    utils::logmesg(lmp,"  find clusters CPU = {:.3f} seconds\n",platform::walltime()-time1);

  // initialize list of SHAKE clusters to constrain

  maxlist = 0;
  list = nullptr;
  closest_list = nullptr;
}

/* ---------------------------------------------------------------------- */

FixShake::~FixShake()
{
  if (copymode) return;

  // unregister callbacks to this fix from Atom class

  if (modify->get_fix_by_id(id)) atom->delete_callback(id,Atom::GROW);

  // set bond_type and angle_type back to positive for SHAKE clusters
  // must set for all SHAKE bonds and angles stored by each atom

  int nlocal = atom->nlocal;

  if (shake_flag)
    for (int i = 0; i < nlocal; i++) {
      if (shake_flag[i] == 0) continue;
      else if (shake_flag[i] == 1) {
        bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],1);
        bondtype_findset(i,shake_atom[i][0],shake_atom[i][2],1);
        angletype_findset(i,shake_atom[i][1],shake_atom[i][2],1);
      } else if (shake_flag[i] == 2) {
        bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],1);
      } else if (shake_flag[i] == 3) {
        bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],1);
        bondtype_findset(i,shake_atom[i][0],shake_atom[i][2],1);
      } else if (shake_flag[i] == 4) {
        bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],1);
        bondtype_findset(i,shake_atom[i][0],shake_atom[i][2],1);
        bondtype_findset(i,shake_atom[i][0],shake_atom[i][3],1);
      }
    }

  // delete locally stored arrays

  memory->destroy(shake_flag);
  memory->destroy(shake_atom);
  memory->destroy(shake_type);
  memory->destroy(xshake);
  memory->destroy(ftmp);
  memory->destroy(vtmp);


  delete[] bond_flag;
  delete[] angle_flag;
  delete[] type_flag;
  delete[] mass_list;

  delete[] bond_distance;
  delete[] angle_distance;

  if (output_every) {
    delete[] b_count;
    delete[] b_count_all;
    delete[] b_ave;
    delete[] b_ave_all;
    delete[] b_max;
    delete[] b_max_all;
    delete[] b_min;
    delete[] b_min_all;

    delete[] a_count;
    delete[] a_count_all;
    delete[] a_ave;
    delete[] a_ave_all;
    delete[] a_max;
    delete[] a_max_all;
    delete[] a_min;
    delete[] a_min_all;
  }

  memory->destroy(list);
  memory->destroy(closest_list);
}

/* ---------------------------------------------------------------------- */

int FixShake::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_PRE_REVERSE;
  mask |= MIN_POST_FORCE;
  mask |= POST_RUN;
  return mask;
}

/* ----------------------------------------------------------------------
   set bond and angle distances
   this init must happen after force->bond and force->angle inits
------------------------------------------------------------------------- */

void FixShake::init()
{
  int i,m,flag,flag_all,type1,type2,bond1_type,bond2_type;
  double rsq,angle;

  // error if more than one shake fix
  auto pattern = fmt::format("^{}",style);

  if (modify->get_fix_by_style(pattern).size() > 1)
    error->all(FLERR,"More than one fix {} instance",style);

  // cannot use with minimization since SHAKE turns off bonds
  // that should contribute to potential energy

  if ((comm->me == 0) && (update->whichflag == 2))
    error->warning(FLERR,"Using fix {} with minimization.\n  Substituting constraints with "
                   "harmonic restraint forces using kbond={:.4g}", style, kbond);

  // error if a fix changing the box comes before shake fix
  bool boxflag = false;
  for (auto &ifix : modify->get_fix_list()) {
   if (boxflag && utils::strmatch(ifix->style,pattern))
     error->all(FLERR,"Fix {} must come before any box changing fix", style);
    if (ifix->box_change) boxflag = true;
  }

  // if rRESPA, find associated fix that must exist
  // could have changed locations in fix list since created
  // set ptrs to rRESPA variables

  fix_respa = nullptr;
  if (utils::strmatch(update->integrate_style,"^respa")) {
    if (update->whichflag > 0) {
      auto fixes = modify->get_fix_by_style("^RESPA");
      if (fixes.size() > 0) fix_respa = dynamic_cast<FixRespa *>(fixes.front());
      else error->all(FLERR,"Run style respa did not create fix RESPA");
    }
    auto respa_style = dynamic_cast<Respa *>(update->integrate);
    nlevels_respa = respa_style->nlevels;
    loop_respa = respa_style->loop;
    step_respa = respa_style->step;
  }

  // set equilibrium bond distances

  if (force->bond == nullptr)
    error->all(FLERR,"Bond style must be defined for fix {}",style);
  for (i = 1; i <= atom->nbondtypes; i++)
    bond_distance[i] = force->bond->equilibrium_distance(i);

  // set equilibrium angle distances

  int nlocal = atom->nlocal;

  for (i = 1; i <= atom->nangletypes; i++) {
    if (angle_flag[i] == 0) continue;
    if (force->angle == nullptr)
      error->all(FLERR,"Angle style must be defined for fix {}",style);

    // scan all atoms for a SHAKE angle cluster
    // extract bond types for the 2 bonds in the cluster
    // bond types must be same in all clusters of this angle type,
    //   else set error flag

    flag = 0;
    bond1_type = bond2_type = 0;
    for (m = 0; m < nlocal; m++) {
      if (shake_flag[m] != 1) continue;
      if (shake_type[m][2] != i) continue;
      type1 = MIN(shake_type[m][0],shake_type[m][1]);
      type2 = MAX(shake_type[m][0],shake_type[m][1]);
      if (bond1_type > 0) {
        if (type1 != bond1_type || type2 != bond2_type) {
          flag = 1;
          break;
        }
      }
      bond1_type = type1;
      bond2_type = type2;
    }

    // error check for any bond types that are not the same

    MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_MAX,world);
    if (flag_all) error->all(FLERR,"Fix {} angles have different bond types", style);

    // ensure all procs have bond types

    MPI_Allreduce(&bond1_type,&flag_all,1,MPI_INT,MPI_MAX,world);
    bond1_type = flag_all;
    MPI_Allreduce(&bond2_type,&flag_all,1,MPI_INT,MPI_MAX,world);
    bond2_type = flag_all;

    // if bond types are 0, no SHAKE angles of this type exist
    // just skip this angle

    if (bond1_type == 0) {
      angle_distance[i] = 0.0;
      continue;
    }

    // compute the angle distance as a function of 2 bond distances
    // formula is now correct for bonds of same or different lengths (Oct15)

    angle = force->angle->equilibrium_angle(i);
    const double b1 = bond_distance[bond1_type];
    const double b2 = bond_distance[bond2_type];
    rsq = b1*b1 + b2*b2 - 2.0*b1*b2*cos(angle);
    angle_distance[i] = sqrt(rsq);
  }
}

/* ----------------------------------------------------------------------
   SHAKE as pre-integrator constraint
------------------------------------------------------------------------- */

void FixShake::setup(int vflag)
{
  pre_neighbor();

  if (output_every) stats();

  // setup SHAKE output

  bigint ntimestep = update->ntimestep;
  if (output_every) {
    next_output = ntimestep + output_every;
    if (ntimestep % output_every != 0)
      next_output = (ntimestep/output_every)*output_every + output_every;
  } else next_output = -1;

  // set respa to 0 if verlet is used and to 1 otherwise

  if (utils::strmatch(update->integrate_style,"^verlet"))
    respa = 0;
  else
    respa = 1;

  if (!respa) {
    dtv     = update->dt;
    dtfsq   = 0.5 * update->dt * update->dt * force->ftm2v;
    if (!rattle) dtfsq = update->dt * update->dt * force->ftm2v;
  } else {
    dtv = step_respa[0];
    dtf_innerhalf = 0.5 * step_respa[0] * force->ftm2v;
    dtf_inner = dtf_innerhalf;
  }

  // correct geometry of cluster if necessary

  correct_coordinates(vflag);

  // remove velocities along any bonds

  correct_velocities();

  // precalculate constraining forces for first integration step

  shake_end_of_step(vflag);
}

/* ----------------------------------------------------------------------
   during minimization fix SHAKE adds strong bond forces
------------------------------------------------------------------------- */

void FixShake::min_setup(int vflag)
{
  pre_neighbor();
  min_post_force(vflag);
}

/* --------------------------------------------------------------------- */

void FixShake::setup_pre_reverse(int eflag, int vflag)
{
  min_pre_reverse(eflag,vflag);
}

/* ----------------------------------------------------------------------
   build list of SHAKE clusters to constrain
   if one or more atoms in cluster are on this proc,
     this proc lists the cluster exactly once
------------------------------------------------------------------------- */

void FixShake::pre_neighbor()
{
  ebond = 0.0;
  int atom1,atom2,atom3,atom4;

  // local copies of atom quantities
  // used by SHAKE until next re-neighboring

  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass = atom->mass;
  rmass = atom->rmass;
  type = atom->type;
  nlocal = atom->nlocal;

  // extend size of SHAKE list if necessary

  if (nlocal > maxlist) {
    maxlist = nlocal;
    memory->destroy(list);
    memory->create(list,maxlist,"shake:list");
    memory->destroy(closest_list);
    memory->create(closest_list,maxlist,4,"shake:closest_list");
  }

  // build list of SHAKE clusters I compute

  nlist = 0;

  for (int i = 0; i < nlocal; i++)
    if (shake_flag[i]) {
      if (shake_flag[i] == 2) {
        atom1 = atom->map(shake_atom[i][0]);
        atom2 = atom->map(shake_atom[i][1]);
        if (atom1 == -1 || atom2 == -1)
          error->one(FLERR,"Shake atoms {} {} missing on proc {} at step {}",shake_atom[i][0],
                     shake_atom[i][1],comm->me,update->ntimestep);
        atom1 = domain->closest_image(i, atom1);
        atom2 = domain->closest_image(i, atom2);
        if (i <= atom1 && i <= atom2) {
          list[nlist] = i;
          closest_list[nlist][0] = atom1;
          closest_list[nlist][1] = atom2;
          nlist++;
        }
      } else if (shake_flag[i] % 2 == 1) {
        atom1 = atom->map(shake_atom[i][0]);
        atom2 = atom->map(shake_atom[i][1]);
        atom3 = atom->map(shake_atom[i][2]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1)
          error->one(FLERR,"Shake atoms {} {} {} missing on proc {} at step {}",shake_atom[i][0],
                                       shake_atom[i][1],shake_atom[i][2],
                                       comm->me,update->ntimestep);
        atom1 = domain->closest_image(i, atom1);
        atom2 = domain->closest_image(i, atom2);
        atom3 = domain->closest_image(i, atom3);
        if (i <= atom1 && i <= atom2 && i <= atom3) {
          list[nlist] = i;
          closest_list[nlist][0] = atom1;
          closest_list[nlist][1] = atom2;
          closest_list[nlist][2] = atom3;
          nlist++;
        }
      } else {
        atom1 = atom->map(shake_atom[i][0]);
        atom2 = atom->map(shake_atom[i][1]);
        atom3 = atom->map(shake_atom[i][2]);
        atom4 = atom->map(shake_atom[i][3]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Shake atoms {} {} {} {} missing on proc {} at step {}",shake_atom[i][0],
                                       shake_atom[i][1],shake_atom[i][2],
                                       shake_atom[i][3],comm->me,update->ntimestep);
        atom1 = domain->closest_image(i, atom1);
        atom2 = domain->closest_image(i, atom2);
        atom3 = domain->closest_image(i, atom3);
        atom4 = domain->closest_image(i, atom4);
        if (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4) {
          list[nlist] = i;
          closest_list[nlist][0] = atom1;
          closest_list[nlist][1] = atom2;
          closest_list[nlist][2] = atom3;
          closest_list[nlist][3] = atom4;
          nlist++;
        }
      }
    }
}

/* ----------------------------------------------------------------------
   compute the force adjustment for SHAKE constraint
------------------------------------------------------------------------- */

void FixShake::post_force(int vflag)
{
  if (update->ntimestep == next_output) stats();

  // xshake = unconstrained move with current v,f
  // communicate results if necessary

  unconstrained_update();
  comm->forward_comm(this);

  // virial setup

  int eflag = eflag_pre_reverse;
  ev_init(eflag, vflag);
  ebond = 0.0;

  // loop over clusters to add constraint forces

  int m;
  for (int i = 0; i < nlist; i++) {
    m = list[i];
    if (shake_flag[m] == 2) shake(i);
    else if (shake_flag[m] == 3) shake3(i);
    else if (shake_flag[m] == 4) shake4(i);
    else shake3angle(i);
  }

  // store vflag for coordinate_constraints_end_of_step()

  vflag_post_force = vflag;
}

/* ----------------------------------------------------------------------
   enforce SHAKE constraints from rRESPA
   xshake prediction portion is different than Verlet
------------------------------------------------------------------------- */

void FixShake::post_force_respa(int vflag, int ilevel, int iloop)
{
  // call stats only on outermost level

  if (ilevel == nlevels_respa-1 && update->ntimestep == next_output) stats();

  // might be OK to skip enforcing SHAKE constraings
  // on last iteration of inner levels if pressure not requested
  // however, leads to slightly different trajectories

  //if (ilevel < nlevels_respa-1 && iloop == loop_respa[ilevel]-1 && !vflag)
  //  return;

  // xshake = unconstrained move with current v,f as function of level
  // communicate results if necessary

  unconstrained_update_respa(ilevel);
  comm->forward_comm(this);

  // virial setup only needed on last iteration of innermost level
  //   and if pressure is requested
  // virial accumulation happens via evflag at last iteration of each level

  if (ilevel == 0 && iloop == loop_respa[ilevel]-1 && vflag) v_init(vflag);
  if (iloop == loop_respa[ilevel]-1) evflag = 1;
  else evflag = 0;

  // loop over clusters to add constraint forces

  int m;
  for (int i = 0; i < nlist; i++) {
    m = list[i];
    if (shake_flag[m] == 2) shake(i);
    else if (shake_flag[m] == 3) shake3(i);
    else if (shake_flag[m] == 4) shake4(i);
    else shake3angle(i);
  }

  // store vflag for coordinate_constraints_end_of_step()
  vflag_post_force = vflag;
}

/* ----------------------------------------------------------------------
   store eflag so it can be used in min_post_force
------------------------------------------------------------------------- */

void FixShake::min_pre_reverse(int eflag, int /*vflag*/)
{
  eflag_pre_reverse = eflag;
}

/* ----------------------------------------------------------------------
   substitute shake constraints with very strong bonds
------------------------------------------------------------------------- */

void FixShake::min_post_force(int vflag)
{
  if (output_every) {
    bigint ntimestep = update->ntimestep;
    if (next_output == ntimestep) stats();

    next_output = ntimestep + output_every;
    if (ntimestep % output_every != 0)
      next_output = (ntimestep/output_every)*output_every + output_every;
  } else next_output = -1;

  int eflag = eflag_pre_reverse;
  ev_init(eflag, vflag);

  x = atom->x;
  f = atom->f;
  nlocal = atom->nlocal;
  ebond = 0.0;

  // loop over shake clusters to add restraint forces

  for (int i = 0; i < nlist; i++) {
    int m = list[i];
    if (shake_flag[m] == 2) {
      bond_force(closest_list[i][0], closest_list[i][1], bond_distance[shake_type[m][0]]);
    } else if (shake_flag[m] == 3) {
      bond_force(closest_list[i][0], closest_list[i][1], bond_distance[shake_type[m][0]]);
      bond_force(closest_list[i][0], closest_list[i][2], bond_distance[shake_type[m][1]]);
    } else if (shake_flag[m] == 4) {
      bond_force(closest_list[i][0], closest_list[i][1], bond_distance[shake_type[m][0]]);
      bond_force(closest_list[i][0], closest_list[i][2], bond_distance[shake_type[m][1]]);
      bond_force(closest_list[i][0], closest_list[i][3], bond_distance[shake_type[m][2]]);
    } else {
      bond_force(closest_list[i][0], closest_list[i][1], bond_distance[shake_type[m][0]]);
      bond_force(closest_list[i][0], closest_list[i][2], bond_distance[shake_type[m][1]]);
      bond_force(closest_list[i][1], closest_list[i][2], angle_distance[shake_type[m][2]]);
    }
  }
}

/* ----------------------------------------------------------------------
   count # of degrees-of-freedom removed by SHAKE for atoms in igroup
------------------------------------------------------------------------- */

bigint FixShake::dof(int igroup)
{
  int groupbit = group->bitmask[igroup];

  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  // count dof in a cluster if and only if
  // the central atom is in group and atom i is the central atom

  bigint n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (shake_flag[i] == 0) continue;
    if (shake_atom[i][0] != tag[i]) continue;
    if (shake_flag[i] == 1) n += 3;
    else if (shake_flag[i] == 2) n += 1;
    else if (shake_flag[i] == 3) n += 2;
    else if (shake_flag[i] == 4) n += 3;
  }

  bigint nall;
  MPI_Allreduce(&n,&nall,1,MPI_LMP_BIGINT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   identify whether each atom is in a SHAKE cluster
   only include atoms in fix group and those bonds/angles specified in input
   test whether all clusters are valid
   set shake_flag, shake_atom, shake_type values
   set bond,angle types negative so will be ignored in neighbor lists
------------------------------------------------------------------------- */

void FixShake::find_clusters()
{
  int i,j,m,n,imol,iatom;
  int flag,flag_all;
  tagint tagprev;
  double massone;

  if (comm->me == 0) utils::logmesg(lmp, "Finding {} clusters ...\n",utils::uppercase(style));
  atommols = atom->avec->onemols;

  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;

  int nlocal = atom->nlocal;
  int angles_allow = atom->avec->angles_allow;

  // -----------------------------------------------------
  // allocate arrays for self (1d) and bond partners (2d)
  // max = max # of bond partners for owned atoms = 2nd dim of partner arrays
  // npartner[i] = # of bonds attached to atom i
  // nshake[i] = # of SHAKE bonds attached to atom i
  // partner_tag[i][] = global IDs of each partner
  // partner_mask[i][] = mask of each partner
  // partner_type[i][] = type of each partner
  // partner_massflag[i][] = 1 if partner meets mass criterion, 0 if not
  // partner_bondtype[i][] = type of bond attached to each partner
  // partner_shake[i][] = 1 if SHAKE bonded to partner, 0 if not
  // partner_nshake[i][] = nshake value for each partner
  // -----------------------------------------------------

  int max = 0;
  if (molecular == Atom::MOLECULAR) {
    for (i = 0; i < nlocal; i++) max = MAX(max,nspecial[i][0]);
  } else {
    for (i = 0; i < nlocal; i++) {
      imol = molindex[i];
      if (imol < 0) continue;
      iatom = molatom[i];
      max = MAX(max,atommols[imol]->nspecial[iatom][0]);
    }
  }

  int *npartner;
  memory->create(npartner,nlocal,"shake:npartner");
  memory->create(nshake,nlocal,"shake:nshake");

  tagint **partner_tag;
  int **partner_mask,**partner_type,**partner_massflag;
  int **partner_bondtype,**partner_shake,**partner_nshake;
  memory->create(partner_tag,nlocal,max,"shake:partner_tag");
  memory->create(partner_mask,nlocal,max,"shake:partner_mask");
  memory->create(partner_type,nlocal,max,"shake:partner_type");
  memory->create(partner_massflag,nlocal,max,"shake:partner_massflag");
  memory->create(partner_bondtype,nlocal,max,"shake:partner_bondtype");
  memory->create(partner_shake,nlocal,max,"shake:partner_shake");
  memory->create(partner_nshake,nlocal,max,"shake:partner_nshake");

  // setup atomIDs and procowner vectors in rendezvous decomposition

  atom_owners();

  // -----------------------------------------------------
  // set npartner and partner_tag from special arrays
  // -----------------------------------------------------

  if (molecular == Atom::MOLECULAR) {
    for (i = 0; i < nlocal; i++) {
      npartner[i] = nspecial[i][0];
      for (j = 0; j < npartner[i]; j++)
        partner_tag[i][j] = special[i][j];
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      imol = molindex[i];
      if (imol < 0) continue;
      iatom = molatom[i];
      tagprev = tag[i] - iatom - 1;
      npartner[i] = atommols[imol]->nspecial[iatom][0];
      for (j = 0; j < npartner[i]; j++)
        partner_tag[i][j] = atommols[imol]->special[iatom][j] + tagprev;
    }
  }

  // -----------------------------------------------------
  // set partner_mask, partner_type, partner_massflag,
  //   partner_bondtype for all my bonded partners
  // requires rendezvous communication for off-proc partners
  // -----------------------------------------------------

  partner_info(npartner,partner_tag,partner_mask,partner_type,
               partner_massflag,partner_bondtype);

  // error check for unfilled partner info
  // if partner_type not set, is an error
  // partner_bondtype may not be set if special list is not consistent
  //   with bondatom (e.g. due to delete_bonds command)
  // this is OK if one or both atoms are not in fix group, since
  //   bond won't be SHAKEn anyway
  // else it's an error

  flag = 0;
  for (i = 0; i < nlocal; i++)
    for (j = 0; j < npartner[i]; j++) {
      if (partner_type[i][j] == 0) flag++;
      if (!(mask[i] & groupbit)) continue;
      if (!(partner_mask[i][j] & groupbit)) continue;
    }

  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Did not find fix {} partner info", style);

  // -----------------------------------------------------
  // identify SHAKEable bonds
  // set nshake[i] = # of SHAKE bonds attached to atom i
  // set partner_shake[i][] = 1 if SHAKE bonded to partner, 0 if not
  // both atoms must be in group, bondtype must be > 0
  // check if bondtype is in input bond_flag
  // check if type of either atom is in input type_flag
  // check if mass of either atom is in input mass_list
  // -----------------------------------------------------

  int np;

  for (i = 0; i < nlocal; i++) {
    nshake[i] = 0;
    np = npartner[i];
    for (j = 0; j < np; j++) {
      partner_shake[i][j] = 0;

      if (!(mask[i] & groupbit)) continue;
      if (!(partner_mask[i][j] & groupbit)) continue;
      if (partner_bondtype[i][j] <= 0) continue;

      if (bond_flag[partner_bondtype[i][j]]) {
        partner_shake[i][j] = 1;
        nshake[i]++;
        continue;
      }
      if (type_flag[type[i]] || type_flag[partner_type[i][j]]) {
        partner_shake[i][j] = 1;
        nshake[i]++;
        continue;
      }
      if (nmass) {
        if (partner_massflag[i][j]) {
          partner_shake[i][j] = 1;
          nshake[i]++;
          continue;
        } else {
          if (rmass) massone = rmass[i];
          else massone = mass[type[i]];
          if (masscheck(massone)) {
            partner_shake[i][j] = 1;
            nshake[i]++;
            continue;
          }
        }
      }
    }
  }

  // -----------------------------------------------------
  // set partner_nshake for bonded partners
  // requires rendezvous communication for off-proc partners
  // -----------------------------------------------------

  nshake_info(npartner,partner_tag,partner_nshake);

  // -----------------------------------------------------
  // error checks
  // no atom with nshake > 3
  // no connected atoms which both have nshake > 1
  // -----------------------------------------------------

  flag = 0;
  for (i = 0; i < nlocal; i++) if (nshake[i] > 3) flag++;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Shake cluster of more than 4 atoms");

  flag = 0;
  for (i = 0; i < nlocal; i++) {
    if (nshake[i] <= 1) continue;
    for (j = 0; j < npartner[i]; j++)
      if (partner_shake[i][j] && partner_nshake[i][j] > 1) flag++;
  }
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Shake clusters are connected");

  // -----------------------------------------------------
  // set SHAKE arrays that are stored with atoms & add angle constraints
  // zero shake arrays for all owned atoms
  // if I am central atom set shake_flag & shake_atom & shake_type
  // for 2-atom clusters, I am central atom if my atom ID < partner ID
  // for 3-atom clusters, test for angle constraint
  //   angle will be stored by this atom if it exists
  //   if angle type matches angle_flag, then it is angle-constrained
  // shake_flag[] = 0 if atom not in SHAKE cluster
  //                2,3,4 = size of bond-only cluster
  //                1 = 3-atom angle cluster
  // shake_atom[][] = global IDs of 2,3,4 atoms in cluster
  //                  central atom is 1st
  //                  for 2-atom cluster, lowest ID is 1st
  // shake_type[][] = bondtype of each bond in cluster
  //                  for 3-atom angle cluster, 3rd value is angletype
  // -----------------------------------------------------

  for (i = 0; i < nlocal; i++) {
    shake_flag[i] = 0;
    shake_atom[i][0] = 0;
    shake_atom[i][1] = 0;
    shake_atom[i][2] = 0;
    shake_atom[i][3] = 0;
    shake_type[i][0] = 0;
    shake_type[i][1] = 0;
    shake_type[i][2] = 0;

    if (nshake[i] == 1) {
      for (j = 0; j < npartner[i]; j++)
        if (partner_shake[i][j]) break;
      if (partner_nshake[i][j] == 1 && tag[i] < partner_tag[i][j]) {
        shake_flag[i] = 2;
        shake_atom[i][0] = tag[i];
        shake_atom[i][1] = partner_tag[i][j];
        shake_type[i][0] = partner_bondtype[i][j];
      }
    }

    if (nshake[i] > 1) {
      shake_flag[i] = 1;
      shake_atom[i][0] = tag[i];
      for (j = 0; j < npartner[i]; j++)
        if (partner_shake[i][j]) {
          m = shake_flag[i];
          shake_atom[i][m] = partner_tag[i][j];
          shake_type[i][m-1] = partner_bondtype[i][j];
          shake_flag[i]++;
        }
    }

    if (nshake[i] == 2 && angles_allow) {
      n = angletype_findset(i,shake_atom[i][1],shake_atom[i][2],0);
      if (n <= 0) continue;
      if (angle_flag[n]) {
        shake_flag[i] = 1;
        shake_type[i][2] = n;
      }
    }
  }

  // -----------------------------------------------------
  // set shake_flag,shake_atom,shake_type for non-central atoms
  // requires rendezvous communication for off-proc atoms
  // -----------------------------------------------------

  shake_info(npartner,partner_tag,partner_shake);

  // -----------------------------------------------------
  // free local memory
  // -----------------------------------------------------

  memory->destroy(atomIDs);
  memory->destroy(procowner);

  memory->destroy(npartner);
  memory->destroy(nshake);
  memory->destroy(partner_tag);
  memory->destroy(partner_mask);
  memory->destroy(partner_type);
  memory->destroy(partner_massflag);
  memory->destroy(partner_bondtype);
  memory->destroy(partner_shake);
  memory->destroy(partner_nshake);

  // -----------------------------------------------------
  // set bond_type and angle_type negative for SHAKE clusters
  // must set for all SHAKE bonds and angles stored by each atom
  // -----------------------------------------------------

  for (i = 0; i < nlocal; i++) {
    if (shake_flag[i] == 0) continue;
    else if (shake_flag[i] == 1) {
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],-1);
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][2],-1);
      angletype_findset(i,shake_atom[i][1],shake_atom[i][2],-1);
    } else if (shake_flag[i] == 2) {
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],-1);
    } else if (shake_flag[i] == 3) {
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],-1);
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][2],-1);
    } else if (shake_flag[i] == 4) {
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],-1);
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][2],-1);
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][3],-1);
    }
  }

  // -----------------------------------------------------
  // print info on SHAKE clusters
  // -----------------------------------------------------

  bigint count1,count2,count3,count4;
  count1 = count2 = count3 = count4 = 0;
  for (i = 0; i < nlocal; i++) {
    if (shake_flag[i] == 1) count1++;
    else if (shake_flag[i] == 2) count2++;
    else if (shake_flag[i] == 3) count3++;
    else if (shake_flag[i] == 4) count4++;
  }

  bigint tmp;
  tmp = count1;
  MPI_Allreduce(&tmp,&count1,1,MPI_LMP_BIGINT,MPI_SUM,world);
  tmp = count2;
  MPI_Allreduce(&tmp,&count2,1,MPI_LMP_BIGINT,MPI_SUM,world);
  tmp = count3;
  MPI_Allreduce(&tmp,&count3,1,MPI_LMP_BIGINT,MPI_SUM,world);
  tmp = count4;
  MPI_Allreduce(&tmp,&count4,1,MPI_LMP_BIGINT,MPI_SUM,world);

  if (comm->me == 0) {
    utils::logmesg(lmp,"{:>8} = # of size 2 clusters\n"
                   "{:>8} = # of size 3 clusters\n"
                   "{:>8} = # of size 4 clusters\n"
                   "{:>8} = # of frozen angles\n",
                   count2/2,count3/3,count4/4,count1/3);
  }
}

/* ----------------------------------------------------------------------
   setup atomIDs and procowner
------------------------------------------------------------------------- */

void FixShake::atom_owners()
{
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int *proclist;
  memory->create(proclist,nlocal,"shake:proclist");
  auto idbuf = (IDRvous *) memory->smalloc((bigint) nlocal*sizeof(IDRvous),"shake:idbuf");

  // setup input buf to rendezvous comm
  // input datums = pairs of bonded atoms
  // owning proc for each datum = random hash of atomID
  // one datum for each owned atom: datum = owning proc, atomID

  for (int i = 0; i < nlocal; i++) {
    proclist[i] = tag[i] % comm->nprocs;
    idbuf[i].me = comm->me;
    idbuf[i].atomID = tag[i];
  }

  // perform rendezvous operation
  // each proc assigned every 1/Pth atom

  char *buf;
  comm->rendezvous(RVOUS,nlocal,(char *) idbuf,sizeof(IDRvous),
                   0,proclist,
                   rendezvous_ids,0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->sfree(idbuf);
}

/* ----------------------------------------------------------------------
   setup partner_mask, partner_type, partner_massflag, partner_bondtype
------------------------------------------------------------------------- */

void FixShake::partner_info(int *npartner, tagint **partner_tag,
                            int **partner_mask, int **partner_type,
                            int **partner_massflag, int **partner_bondtype)
{
  int i,j,m,n;
  int nlocal = atom->nlocal;

  // nsend = # of my datums to send
  // one datum for every off-processor partner

  int nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < npartner[i]; j++) {
      m = atom->map(partner_tag[i][j]);
      if (m < 0 || m >= nlocal) nsend++;
    }
  }

  int *proclist;
  memory->create(proclist,nsend,"special:proclist");
  auto inbuf = (PartnerInfo *) memory->smalloc((bigint) nsend*sizeof(PartnerInfo),"special:inbuf");

  // set values in 4 partner arrays for all partner atoms I own
  // also setup input buf to rendezvous comm
  // input datums = pair of bonded atoms where I do not own partner
  // owning proc for each datum = partner_tag % comm->nprocs
  // datum: atomID = partner_tag (off-proc), partnerID = tag (on-proc)
  //        4 values for my owned atom

  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;

  double massone;

  nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < npartner[i]; j++) {
      partner_mask[i][j] = 0;
      partner_type[i][j] = 0;
      partner_massflag[i][j] = 0;
      partner_bondtype[i][j] = 0;

      m = atom->map(partner_tag[i][j]);

      if (m >= 0 && m < nlocal) {
        partner_mask[i][j] = mask[m];
        partner_type[i][j] = type[m];
        if (nmass) {
          if (rmass) massone = rmass[m];
          else massone = mass[type[m]];
          partner_massflag[i][j] = masscheck(massone);
        }
        n = bondtype_findset(i,tag[i],partner_tag[i][j],0);
        if (n) partner_bondtype[i][j] = n;
        else {
          n = bondtype_findset(m,tag[i],partner_tag[i][j],0);
          if (n) partner_bondtype[i][j] = n;
        }

      } else {
        proclist[nsend] = partner_tag[i][j] % comm->nprocs;
        inbuf[nsend].atomID = partner_tag[i][j];
        inbuf[nsend].partnerID = tag[i];
        inbuf[nsend].mask = mask[i];
        inbuf[nsend].type = type[i];
        if (nmass) {
          if (rmass) massone = rmass[i];
          else massone = mass[type[i]];
          inbuf[nsend].massflag = masscheck(massone);
        } else inbuf[nsend].massflag = 0;

        // my atom may own bond, in which case set partner_bondtype
        // else receiver of this datum will own the bond and return the value

        n = bondtype_findset(i,tag[i],partner_tag[i][j],0);
        if (n) {
          partner_bondtype[i][j] = n;
          inbuf[nsend].bondtype = n;
        } else inbuf[nsend].bondtype = 0;

        nsend++;
      }
    }
  }

  // perform rendezvous operation
  // each proc owns random subset of atoms
  // receives all data needed to populate un-owned partner 4 values

  char *buf;
  int nreturn = comm->rendezvous(RVOUS,nsend,(char *) inbuf,sizeof(PartnerInfo),
                                 0,proclist,
                                 rendezvous_partners_info,
                                 0,buf,sizeof(PartnerInfo),
                                 (void *) this);
  auto outbuf = (PartnerInfo *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set partner 4 values for un-onwed partners based on output info
  // outbuf.atomID = my owned atom, outbuf.partnerID = partner the info is for

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    for (j = 0; j < npartner[i]; j++)
      if (partner_tag[i][j] == outbuf[m].partnerID) break;
    partner_mask[i][j] = outbuf[m].mask;
    partner_type[i][j] = outbuf[m].type;
    partner_massflag[i][j] = outbuf[m].massflag;

    // only set partner_bondtype if my atom did not set it
    //   when setting up rendezvous
    // if this proc set it, then sender of this datum set outbuf.bondtype = 0

    if (partner_bondtype[i][j] == 0)
      partner_bondtype[i][j] = outbuf[m].bondtype;
  }

  memory->sfree(outbuf);
}

/* ----------------------------------------------------------------------
   setup partner_nshake
------------------------------------------------------------------------- */

void FixShake::nshake_info(int *npartner, tagint **partner_tag,
                           int **partner_nshake)
{
  int i,j,m;
  int nlocal = atom->nlocal;

  // nsend = # of my datums to send
  // one datum for every off-processor partner

  int nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < npartner[i]; j++) {
      m = atom->map(partner_tag[i][j]);
      if (m < 0 || m >= nlocal) nsend++;
    }
  }

  int *proclist;
  memory->create(proclist,nsend,"special:proclist");
  auto inbuf = (NShakeInfo *) memory->smalloc((bigint) nsend*sizeof(NShakeInfo),"special:inbuf");

  // set partner_nshake for all partner atoms I own
  // also setup input buf to rendezvous comm
  // input datums = pair of bonded atoms where I do not own partner
  // owning proc for each datum = partner_tag % comm->nprocs
  // datum: atomID = partner_tag (off-proc), partnerID = tag (on-proc)
  //        nshake value for my owned atom

  tagint *tag = atom->tag;

  nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < npartner[i]; j++) {
      partner_nshake[i][j] = 0;
      m = atom->map(partner_tag[i][j]);
      if (m >= 0 && m < nlocal) {
        partner_nshake[i][j] = nshake[m];
      } else {
        proclist[nsend] = partner_tag[i][j] % comm->nprocs;
        inbuf[nsend].atomID = partner_tag[i][j];
        inbuf[nsend].partnerID = tag[i];
        inbuf[nsend].nshake = nshake[i];
        nsend++;
      }
    }
  }

  // perform rendezvous operation
  // each proc owns random subset of atoms
  // receives all data needed to populate un-owned partner nshake

  char *buf;
  int nreturn = comm->rendezvous(RVOUS,nsend,(char *) inbuf,sizeof(NShakeInfo),
                                 0,proclist,
                                 rendezvous_nshake,0,buf,sizeof(NShakeInfo),
                                 (void *) this);
  auto outbuf = (NShakeInfo *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set partner nshake for un-onwed partners based on output info
  // outbuf.atomID = my owned atom, outbuf.partnerID = partner the info is for

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    for (j = 0; j < npartner[i]; j++)
      if (partner_tag[i][j] == outbuf[m].partnerID) break;
    partner_nshake[i][j] = outbuf[m].nshake;
  }

  memory->sfree(outbuf);
}

/* ----------------------------------------------------------------------
   setup shake_flag, shake_atom, shake_type
------------------------------------------------------------------------- */

void FixShake::shake_info(int *npartner, tagint **partner_tag,
                          int **partner_shake)
{
  int i,j,m;
  int nlocal = atom->nlocal;

  // nsend = # of my datums to send
  // one datum for every off-processor partner

  int nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < npartner[i]; j++) {
      m = atom->map(partner_tag[i][j]);
      if (m < 0 || m >= nlocal) nsend++;
    }
  }

  int *proclist;
  memory->create(proclist,nsend,"special:proclist");
  auto inbuf = (ShakeInfo *) memory->smalloc((bigint) nsend*sizeof(ShakeInfo),"special:inbuf");

  // set 3 shake arrays for all partner atoms I own
  // also setup input buf to rendezvous comm
  // input datums = partner atom where I do not own partner
  // owning proc for each datum = partner_tag % comm->nprocs
  // datum: atomID = partner_tag (off-proc)
  //        values in 3 shake arrays

  nsend = 0;
  for (i = 0; i < nlocal; i++) {
    if (shake_flag[i] == 0) continue;
    for (j = 0; j < npartner[i]; j++) {
      if (partner_shake[i][j] == 0) continue;
      m = atom->map(partner_tag[i][j]);

      if (m >= 0 && m < nlocal) {
        shake_flag[m] = shake_flag[i];
        shake_atom[m][0] = shake_atom[i][0];
        shake_atom[m][1] = shake_atom[i][1];
        shake_atom[m][2] = shake_atom[i][2];
        shake_atom[m][3] = shake_atom[i][3];
        shake_type[m][0] = shake_type[i][0];
        shake_type[m][1] = shake_type[i][1];
        shake_type[m][2] = shake_type[i][2];

      } else {
        proclist[nsend] = partner_tag[i][j] % comm->nprocs;
        inbuf[nsend].atomID = partner_tag[i][j];
        inbuf[nsend].shake_flag = shake_flag[i];
        inbuf[nsend].shake_atom[0] = shake_atom[i][0];
        inbuf[nsend].shake_atom[1] = shake_atom[i][1];
        inbuf[nsend].shake_atom[2] = shake_atom[i][2];
        inbuf[nsend].shake_atom[3] = shake_atom[i][3];
        inbuf[nsend].shake_type[0] = shake_type[i][0];
        inbuf[nsend].shake_type[1] = shake_type[i][1];
        inbuf[nsend].shake_type[2] = shake_type[i][2];
        nsend++;
      }
    }
  }

  // perform rendezvous operation
  // each proc owns random subset of atoms
  // receives all data needed to populate un-owned shake info

  char *buf;
  int nreturn = comm->rendezvous(RVOUS,nsend,(char *) inbuf,sizeof(ShakeInfo),
                                 0,proclist,
                                 rendezvous_shake,0,buf,sizeof(ShakeInfo),
                                 (void *) this);
  auto outbuf = (ShakeInfo *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set shake info for un-onwed partners based on output info

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    shake_flag[i] = outbuf[m].shake_flag;
    shake_atom[i][0] = outbuf[m].shake_atom[0];
    shake_atom[i][1] = outbuf[m].shake_atom[1];
    shake_atom[i][2] = outbuf[m].shake_atom[2];
    shake_atom[i][3] = outbuf[m].shake_atom[3];
    shake_type[i][0] = outbuf[m].shake_type[0];
    shake_type[i][1] = outbuf[m].shake_type[1];
    shake_type[i][2] = outbuf[m].shake_type[2];
  }

  memory->sfree(outbuf);
}

/* ----------------------------------------------------------------------
   process data for atoms assigned to me in rendezvous decomposition
   inbuf = list of N IDRvous datums
   no outbuf
------------------------------------------------------------------------- */

int FixShake::rendezvous_ids(int n, char *inbuf,
                             int &flag, int *& /*proclist*/, char *& /*outbuf*/,
                             void *ptr)
{
  auto fsptr = (FixShake *) ptr;
  Memory *memory = fsptr->memory;

  tagint *atomIDs;
  int *procowner;

  memory->create(atomIDs,n,"special:atomIDs");
  memory->create(procowner,n,"special:procowner");

  auto in = (IDRvous *) inbuf;

  for (int i = 0; i < n; i++) {
    atomIDs[i] = in[i].atomID;
    procowner[i] = in[i].me;
  }

  // store rendezvous data in FixShake class

  fsptr->nrvous = n;
  fsptr->atomIDs = atomIDs;
  fsptr->procowner = procowner;

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   process data for atoms assigned to me in rendezvous decomposition
   inbuf = list of N PairRvous datums
   outbuf = same list of N PairRvous datums, routed to different procs
------------------------------------------------------------------------- */

int FixShake::rendezvous_partners_info(int n, char *inbuf,
                                       int &flag, int *&proclist, char *&outbuf,
                                       void *ptr)
{
  int i,m;

  auto fsptr = (FixShake *) ptr;
  Atom *atom = fsptr->atom;
  Memory *memory = fsptr->memory;

  // clear atom map so it can be here as a hash table
  // faster than an STL map for large atom counts

  atom->map_clear();

  // hash atom IDs stored in rendezvous decomposition

  int nrvous = fsptr->nrvous;
  tagint *atomIDs = fsptr->atomIDs;

  for (i = 0; i < nrvous; i++)
    atom->map_one(atomIDs[i],i);

  // proclist = owner of atomID in caller decomposition
  // outbuf = info about owned atomID = 4 values

  auto in = (PartnerInfo *) inbuf;
  int *procowner = fsptr->procowner;
  memory->create(proclist,n,"shake:proclist");

  for (i = 0; i < n; i++) {
    m = atom->map(in[i].atomID);
    proclist[i] = procowner[m];
  }

  outbuf = inbuf;

  // re-create atom map

  atom->map_init(0);
  atom->nghost = 0;
  atom->map_set();

  // flag = 1: outbuf = inbuf

  flag = 1;
  return n;
}

/* ----------------------------------------------------------------------
   process data for atoms assigned to me in rendezvous decomposition
   inbuf = list of N NShakeInfo datums
   outbuf = same list of N NShakeInfo datums, routed to different procs
------------------------------------------------------------------------- */

int FixShake::rendezvous_nshake(int n, char *inbuf,
                                int &flag, int *&proclist, char *&outbuf,
                                void *ptr)
{
  int i,m;

  auto fsptr = (FixShake *) ptr;
  Atom *atom = fsptr->atom;
  Memory *memory = fsptr->memory;

  // clear atom map so it can be here as a hash table
  // faster than an STL map for large atom counts

  atom->map_clear();

  // hash atom IDs stored in rendezvous decomposition

  int nrvous = fsptr->nrvous;
  tagint *atomIDs = fsptr->atomIDs;

  for (i = 0; i < nrvous; i++)
    atom->map_one(atomIDs[i],i);

  // proclist = owner of atomID in caller decomposition
  // outbuf = info about owned atomID

  auto in = (NShakeInfo *) inbuf;
  int *procowner = fsptr->procowner;
  memory->create(proclist,n,"shake:proclist");

  for (i = 0; i < n; i++) {
    m = atom->map(in[i].atomID);
    proclist[i] = procowner[m];
  }

  outbuf = inbuf;

  // re-create atom map

  atom->map_init(0);
  atom->nghost = 0;
  atom->map_set();

  // flag = 1: outbuf = inbuf

  flag = 1;
  return n;
}
/* ----------------------------------------------------------------------
   process data for atoms assigned to me in rendezvous decomposition
   inbuf = list of N PairRvous datums
   outbuf = same list of N PairRvous datums, routed to different procs
------------------------------------------------------------------------- */

int FixShake::rendezvous_shake(int n, char *inbuf,
                               int &flag, int *&proclist, char *&outbuf,
                               void *ptr)
{
  int i,m;

  auto fsptr = (FixShake *) ptr;
  Atom *atom = fsptr->atom;
  Memory *memory = fsptr->memory;

  // clear atom map so it can be here as a hash table
  // faster than an STL map for large atom counts

  atom->map_clear();

  // hash atom IDs stored in rendezvous decomposition

  int nrvous = fsptr->nrvous;
  tagint *atomIDs = fsptr->atomIDs;

  for (i = 0; i < nrvous; i++)
    atom->map_one(atomIDs[i],i);

  // proclist = owner of atomID in caller decomposition
  // outbuf = info about owned atomID

  auto in = (ShakeInfo *) inbuf;
  int *procowner = fsptr->procowner;
  memory->create(proclist,n,"shake:proclist");

  for (i = 0; i < n; i++) {
    m = atom->map(in[i].atomID);
    proclist[i] = procowner[m];
  }

  outbuf = inbuf;

  // re-create atom map

  atom->map_init(0);
  atom->nghost = 0;
  atom->map_set();

  // flag = 1: outbuf = inbuf;

  flag = 1;
  return n;
}

/* ----------------------------------------------------------------------
   check if massone is within MASSDELTA of any mass in mass_list
   return 1 if yes, 0 if not
------------------------------------------------------------------------- */

int FixShake::masscheck(double massone)
{
  for (int i = 0; i < nmass; i++)
    if (fabs(mass_list[i]-massone) <= MASSDELTA) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   update the unconstrained position of each atom
   only for SHAKE clusters, else set to 0.0
   assumes NVE update, seems to be accurate enough for NVT,NPT,NPH as well
------------------------------------------------------------------------- */

void FixShake::unconstrained_update()
{
  double dtfmsq;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (shake_flag[i]) {
        dtfmsq = dtfsq / rmass[i];
        xshake[i][0] = x[i][0] + dtv*v[i][0] + dtfmsq*f[i][0];
        xshake[i][1] = x[i][1] + dtv*v[i][1] + dtfmsq*f[i][1];
        xshake[i][2] = x[i][2] + dtv*v[i][2] + dtfmsq*f[i][2];
      } else xshake[i][2] = xshake[i][1] = xshake[i][0] = 0.0;
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (shake_flag[i]) {
        dtfmsq = dtfsq / mass[type[i]];
        xshake[i][0] = x[i][0] + dtv*v[i][0] + dtfmsq*f[i][0];
        xshake[i][1] = x[i][1] + dtv*v[i][1] + dtfmsq*f[i][1];
        xshake[i][2] = x[i][2] + dtv*v[i][2] + dtfmsq*f[i][2];
      } else xshake[i][2] = xshake[i][1] = xshake[i][0] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   update the unconstrained position of each atom in a rRESPA step
   only for SHAKE clusters, else set to 0.0
   assumes NVE update, seems to be accurate enough for NVT,NPT,NPH as well
------------------------------------------------------------------------- */

void FixShake::unconstrained_update_respa(int ilevel)
{
  // xshake = atom coords after next x update in innermost loop
  // depends on rRESPA level
  // for levels > 0 this includes more than one velocity update
  // xshake = predicted position from call to this routine at level N =
  // x + dt0 (v + dtN/m fN + 1/2 dt(N-1)/m f(N-1) + ... + 1/2 dt0/m f0)
  // also set dtfsq = dt0*dtN so that shake,shake3,etc can use it

  double ***f_level = fix_respa->f_level;
  dtfsq = dtf_inner * step_respa[ilevel];

  double invmass,dtfmsq;
  int jlevel;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (shake_flag[i]) {
        invmass = 1.0 / rmass[i];
        dtfmsq = dtfsq * invmass;
        xshake[i][0] = x[i][0] + dtv*v[i][0] + dtfmsq*f[i][0];
        xshake[i][1] = x[i][1] + dtv*v[i][1] + dtfmsq*f[i][1];
        xshake[i][2] = x[i][2] + dtv*v[i][2] + dtfmsq*f[i][2];
        for (jlevel = 0; jlevel < ilevel; jlevel++) {
          dtfmsq = dtf_innerhalf * step_respa[jlevel] * invmass;
          xshake[i][0] += dtfmsq*f_level[i][jlevel][0];
          xshake[i][1] += dtfmsq*f_level[i][jlevel][1];
          xshake[i][2] += dtfmsq*f_level[i][jlevel][2];
        }
      } else xshake[i][2] = xshake[i][1] = xshake[i][0] = 0.0;
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      if (shake_flag[i]) {
        invmass = 1.0 / mass[type[i]];
        dtfmsq = dtfsq * invmass;
        xshake[i][0] = x[i][0] + dtv*v[i][0] + dtfmsq*f[i][0];
        xshake[i][1] = x[i][1] + dtv*v[i][1] + dtfmsq*f[i][1];
        xshake[i][2] = x[i][2] + dtv*v[i][2] + dtfmsq*f[i][2];
        for (jlevel = 0; jlevel < ilevel; jlevel++) {
          dtfmsq = dtf_innerhalf * step_respa[jlevel] * invmass;
          xshake[i][0] += dtfmsq*f_level[i][jlevel][0];
          xshake[i][1] += dtfmsq*f_level[i][jlevel][1];
          xshake[i][2] += dtfmsq*f_level[i][jlevel][2];
        }
      } else xshake[i][2] = xshake[i][1] = xshake[i][0] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   calculate SHAKE constraint forces for size 2 cluster = single bond
------------------------------------------------------------------------- */

void FixShake::shake(int ilist)
{
  int atomlist[2];
  double v[6];
  double invmass0,invmass1;

  // local atom IDs and constraint distances

  int m = list[ilist];
  int i0 = closest_list[ilist][0];
  int i1 = closest_list[ilist][1];
  double bond1 = bond_distance[shake_type[m][0]];

  // r01 = distance vec between atoms

  double r01[3];
  r01[0] = x[i0][0] - x[i1][0];
  r01[1] = x[i0][1] - x[i1][1];
  r01[2] = x[i0][2] - x[i1][2];

  // s01 = distance vec after unconstrained update

  double s01[3];
  s01[0] = xshake[i0][0] - xshake[i1][0];
  s01[1] = xshake[i0][1] - xshake[i1][1];
  s01[2] = xshake[i0][2] - xshake[i1][2];

  // scalar distances between atoms

  double r01sq = r01[0]*r01[0] + r01[1]*r01[1] + r01[2]*r01[2];
  double s01sq = s01[0]*s01[0] + s01[1]*s01[1] + s01[2]*s01[2];

  // a,b,c = coeffs in quadratic equation for lamda

  if (rmass) {
    invmass0 = 1.0/rmass[i0];
    invmass1 = 1.0/rmass[i1];
  } else {
    invmass0 = 1.0/mass[type[i0]];
    invmass1 = 1.0/mass[type[i1]];
  }

  double a = (invmass0+invmass1)*(invmass0+invmass1) * r01sq;
  double b = 2.0 * (invmass0+invmass1) *
    (s01[0]*r01[0] + s01[1]*r01[1] + s01[2]*r01[2]);
  double c = s01sq - bond1*bond1;

  // error check

  double determ = b*b - 4.0*a*c;
  if (determ < 0.0) {
    error->warning(FLERR,"Shake determinant < 0.0");
    determ = 0.0;
  }

  // exact quadratic solution for lamda

  double lamda,lamda1,lamda2;
  lamda1 = (-b+sqrt(determ)) / (2.0*a);
  lamda2 = (-b-sqrt(determ)) / (2.0*a);

  if (fabs(lamda1) <= fabs(lamda2)) lamda = lamda1;
  else lamda = lamda2;

  // update forces if atom is owned by this processor

  lamda /= dtfsq;

  if (i0 < nlocal) {
    f[i0][0] += lamda*r01[0];
    f[i0][1] += lamda*r01[1];
    f[i0][2] += lamda*r01[2];
  }

  if (i1 < nlocal) {
    f[i1][0] -= lamda*r01[0];
    f[i1][1] -= lamda*r01[1];
    f[i1][2] -= lamda*r01[2];
  }

  if (evflag) {
    int count = 0;
    if (i0 < nlocal) atomlist[count++] = i0;
    if (i1 < nlocal) atomlist[count++] = i1;

    v[0] = lamda*r01[0]*r01[0];
    v[1] = lamda*r01[1]*r01[1];
    v[2] = lamda*r01[2]*r01[2];
    v[3] = lamda*r01[0]*r01[1];
    v[4] = lamda*r01[0]*r01[2];
    v[5] = lamda*r01[1]*r01[2];

    double fpairlist[] = {lamda};
    double dellist[][3]  = {{r01[0], r01[1], r01[2]}};
    int pairlist[][2] = {{i0,i1}};
    v_tally(count,atomlist,2.0,v,nlocal,1,pairlist,fpairlist,dellist);
  }
}

/* ----------------------------------------------------------------------
   calculate SHAKE constraint forces for size 3 cluster = two bonds
------------------------------------------------------------------------- */

void FixShake::shake3(int ilist)
{
  int atomlist[3];
  double v[6];
  double invmass0,invmass1,invmass2;

  // local atom IDs and constraint distances

  int m = list[ilist];
  int i0 = closest_list[ilist][0];
  int i1 = closest_list[ilist][1];
  int i2 = closest_list[ilist][2];
  double bond1 = bond_distance[shake_type[m][0]];
  double bond2 = bond_distance[shake_type[m][1]];

  // r01,r02 = distance vec between atoms

  double r01[3];
  r01[0] = x[i0][0] - x[i1][0];
  r01[1] = x[i0][1] - x[i1][1];
  r01[2] = x[i0][2] - x[i1][2];

  double r02[3];
  r02[0] = x[i0][0] - x[i2][0];
  r02[1] = x[i0][1] - x[i2][1];
  r02[2] = x[i0][2] - x[i2][2];

  // s01,s02 = distance vec after unconstrained update

  double s01[3];
  s01[0] = xshake[i0][0] - xshake[i1][0];
  s01[1] = xshake[i0][1] - xshake[i1][1];
  s01[2] = xshake[i0][2] - xshake[i1][2];

  double s02[3];
  s02[0] = xshake[i0][0] - xshake[i2][0];
  s02[1] = xshake[i0][1] - xshake[i2][1];
  s02[2] = xshake[i0][2] - xshake[i2][2];

  // scalar distances between atoms

  double r01sq = r01[0]*r01[0] + r01[1]*r01[1] + r01[2]*r01[2];
  double r02sq = r02[0]*r02[0] + r02[1]*r02[1] + r02[2]*r02[2];
  double s01sq = s01[0]*s01[0] + s01[1]*s01[1] + s01[2]*s01[2];
  double s02sq = s02[0]*s02[0] + s02[1]*s02[1] + s02[2]*s02[2];

  // matrix coeffs and rhs for lamda equations

  if (rmass) {
    invmass0 = 1.0/rmass[i0];
    invmass1 = 1.0/rmass[i1];
    invmass2 = 1.0/rmass[i2];
  } else {
    invmass0 = 1.0/mass[type[i0]];
    invmass1 = 1.0/mass[type[i1]];
    invmass2 = 1.0/mass[type[i2]];
  }

  double a11 = 2.0 * (invmass0+invmass1) *
    (s01[0]*r01[0] + s01[1]*r01[1] + s01[2]*r01[2]);
  double a12 = 2.0 * invmass0 *
    (s01[0]*r02[0] + s01[1]*r02[1] + s01[2]*r02[2]);
  double a21 = 2.0 * invmass0 *
    (s02[0]*r01[0] + s02[1]*r01[1] + s02[2]*r01[2]);
  double a22 = 2.0 * (invmass0+invmass2) *
    (s02[0]*r02[0] + s02[1]*r02[1] + s02[2]*r02[2]);

  // inverse of matrix

  double determ = a11*a22 - a12*a21;
  if (determ == 0.0) error->one(FLERR,"Shake determinant = 0.0");
  double determinv = 1.0/determ;

  double a11inv = a22*determinv;
  double a12inv = -a12*determinv;
  double a21inv = -a21*determinv;
  double a22inv = a11*determinv;

  // quadratic correction coeffs

  double r0102 = (r01[0]*r02[0] + r01[1]*r02[1] + r01[2]*r02[2]);

  double quad1_0101 = (invmass0+invmass1)*(invmass0+invmass1) * r01sq;
  double quad1_0202 = invmass0*invmass0 * r02sq;
  double quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102;

  double quad2_0202 = (invmass0+invmass2)*(invmass0+invmass2) * r02sq;
  double quad2_0101 = invmass0*invmass0 * r01sq;
  double quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102;

  // iterate until converged

  double lamda01 = 0.0;
  double lamda02 = 0.0;
  int niter = 0;
  int done = 0;

  double quad1,quad2,b1,b2,lamda01_new,lamda02_new;

  while (!done && niter < max_iter) {
    quad1 = quad1_0101 * lamda01*lamda01 + quad1_0202 * lamda02*lamda02 +
      quad1_0102 * lamda01*lamda02;
    quad2 = quad2_0101 * lamda01*lamda01 + quad2_0202 * lamda02*lamda02 +
      quad2_0102 * lamda01*lamda02;

    b1 = bond1*bond1 - s01sq - quad1;
    b2 = bond2*bond2 - s02sq - quad2;

    lamda01_new = a11inv*b1 + a12inv*b2;
    lamda02_new = a21inv*b1 + a22inv*b2;

    done = 1;
    if (fabs(lamda01_new-lamda01) > tolerance) done = 0;
    if (fabs(lamda02_new-lamda02) > tolerance) done = 0;

    lamda01 = lamda01_new;
    lamda02 = lamda02_new;

    // stop iterations before we have a floating point overflow
    // max double is < 1.0e308, so 1e150 is a reasonable cutoff

    if (fabs(lamda01) > 1e150 || fabs(lamda02) > 1e150) done = 1;

    niter++;
  }

  // update forces if atom is owned by this processor

  lamda01 = lamda01/dtfsq;
  lamda02 = lamda02/dtfsq;

  if (i0 < nlocal) {
    f[i0][0] += lamda01*r01[0] + lamda02*r02[0];
    f[i0][1] += lamda01*r01[1] + lamda02*r02[1];
    f[i0][2] += lamda01*r01[2] + lamda02*r02[2];
  }

  if (i1 < nlocal) {
    f[i1][0] -= lamda01*r01[0];
    f[i1][1] -= lamda01*r01[1];
    f[i1][2] -= lamda01*r01[2];
  }

  if (i2 < nlocal) {
    f[i2][0] -= lamda02*r02[0];
    f[i2][1] -= lamda02*r02[1];
    f[i2][2] -= lamda02*r02[2];
  }

  if (evflag) {
    int count = 0;
    if (i0 < nlocal) atomlist[count++] = i0;
    if (i1 < nlocal) atomlist[count++] = i1;
    if (i2 < nlocal) atomlist[count++] = i2;

    v[0] = lamda01*r01[0]*r01[0] + lamda02*r02[0]*r02[0];
    v[1] = lamda01*r01[1]*r01[1] + lamda02*r02[1]*r02[1];
    v[2] = lamda01*r01[2]*r01[2] + lamda02*r02[2]*r02[2];
    v[3] = lamda01*r01[0]*r01[1] + lamda02*r02[0]*r02[1];
    v[4] = lamda01*r01[0]*r01[2] + lamda02*r02[0]*r02[2];
    v[5] = lamda01*r01[1]*r01[2] + lamda02*r02[1]*r02[2];

    double fpairlist[] = {lamda01, lamda02};
    double dellist[][3]  = {{r01[0], r01[1], r01[2]},
                            {r02[0], r02[1], r02[2]}};
    int pairlist[][2] = {{i0,i1}, {i0,i2}};
    v_tally(count,atomlist,3.0,v,nlocal,2,pairlist,fpairlist,dellist);
  }
}

/* ----------------------------------------------------------------------
   calculate SHAKE constraint forces for size 4 cluster = three bonds
------------------------------------------------------------------------- */

void FixShake::shake4(int ilist)
{
 int atomlist[4];
  double v[6];
  double invmass0,invmass1,invmass2,invmass3;

  // local atom IDs and constraint distances

  int m = list[ilist];
  int i0 = closest_list[ilist][0];
  int i1 = closest_list[ilist][1];
  int i2 = closest_list[ilist][2];
  int i3 = closest_list[ilist][3];
  double bond1 = bond_distance[shake_type[m][0]];
  double bond2 = bond_distance[shake_type[m][1]];
  double bond3 = bond_distance[shake_type[m][2]];

  // r01,r02,r03 = distance vec between atoms

  double r01[3];
  r01[0] = x[i0][0] - x[i1][0];
  r01[1] = x[i0][1] - x[i1][1];
  r01[2] = x[i0][2] - x[i1][2];

  double r02[3];
  r02[0] = x[i0][0] - x[i2][0];
  r02[1] = x[i0][1] - x[i2][1];
  r02[2] = x[i0][2] - x[i2][2];

  double r03[3];
  r03[0] = x[i0][0] - x[i3][0];
  r03[1] = x[i0][1] - x[i3][1];
  r03[2] = x[i0][2] - x[i3][2];

  // s01,s02,s03 = distance vec after unconstrained update

  double s01[3];
  s01[0] = xshake[i0][0] - xshake[i1][0];
  s01[1] = xshake[i0][1] - xshake[i1][1];
  s01[2] = xshake[i0][2] - xshake[i1][2];

  double s02[3];
  s02[0] = xshake[i0][0] - xshake[i2][0];
  s02[1] = xshake[i0][1] - xshake[i2][1];
  s02[2] = xshake[i0][2] - xshake[i2][2];

  double s03[3];
  s03[0] = xshake[i0][0] - xshake[i3][0];
  s03[1] = xshake[i0][1] - xshake[i3][1];
  s03[2] = xshake[i0][2] - xshake[i3][2];

  // scalar distances between atoms

  double r01sq = r01[0]*r01[0] + r01[1]*r01[1] + r01[2]*r01[2];
  double r02sq = r02[0]*r02[0] + r02[1]*r02[1] + r02[2]*r02[2];
  double r03sq = r03[0]*r03[0] + r03[1]*r03[1] + r03[2]*r03[2];
  double s01sq = s01[0]*s01[0] + s01[1]*s01[1] + s01[2]*s01[2];
  double s02sq = s02[0]*s02[0] + s02[1]*s02[1] + s02[2]*s02[2];
  double s03sq = s03[0]*s03[0] + s03[1]*s03[1] + s03[2]*s03[2];

  // matrix coeffs and rhs for lamda equations

  if (rmass) {
    invmass0 = 1.0/rmass[i0];
    invmass1 = 1.0/rmass[i1];
    invmass2 = 1.0/rmass[i2];
    invmass3 = 1.0/rmass[i3];
  } else {
    invmass0 = 1.0/mass[type[i0]];
    invmass1 = 1.0/mass[type[i1]];
    invmass2 = 1.0/mass[type[i2]];
    invmass3 = 1.0/mass[type[i3]];
  }

  double a11 = 2.0 * (invmass0+invmass1) *
    (s01[0]*r01[0] + s01[1]*r01[1] + s01[2]*r01[2]);
  double a12 = 2.0 * invmass0 *
    (s01[0]*r02[0] + s01[1]*r02[1] + s01[2]*r02[2]);
  double a13 = 2.0 * invmass0 *
    (s01[0]*r03[0] + s01[1]*r03[1] + s01[2]*r03[2]);
  double a21 = 2.0 * invmass0 *
    (s02[0]*r01[0] + s02[1]*r01[1] + s02[2]*r01[2]);
  double a22 = 2.0 * (invmass0+invmass2) *
    (s02[0]*r02[0] + s02[1]*r02[1] + s02[2]*r02[2]);
  double a23 = 2.0 * invmass0 *
    (s02[0]*r03[0] + s02[1]*r03[1] + s02[2]*r03[2]);
  double a31 = 2.0 * invmass0 *
    (s03[0]*r01[0] + s03[1]*r01[1] + s03[2]*r01[2]);
  double a32 = 2.0 * invmass0 *
    (s03[0]*r02[0] + s03[1]*r02[1] + s03[2]*r02[2]);
  double a33 = 2.0 * (invmass0+invmass3) *
    (s03[0]*r03[0] + s03[1]*r03[1] + s03[2]*r03[2]);

  // inverse of matrix;

  double determ = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 -
    a11*a23*a32 - a12*a21*a33 - a13*a22*a31;
  if (determ == 0.0) error->one(FLERR,"Shake determinant = 0.0");
  double determinv = 1.0/determ;

  double a11inv = determinv * (a22*a33 - a23*a32);
  double a12inv = -determinv * (a12*a33 - a13*a32);
  double a13inv = determinv * (a12*a23 - a13*a22);
  double a21inv = -determinv * (a21*a33 - a23*a31);
  double a22inv = determinv * (a11*a33 - a13*a31);
  double a23inv = -determinv * (a11*a23 - a13*a21);
  double a31inv = determinv * (a21*a32 - a22*a31);
  double a32inv = -determinv * (a11*a32 - a12*a31);
  double a33inv = determinv * (a11*a22 - a12*a21);

  // quadratic correction coeffs

  double r0102 = (r01[0]*r02[0] + r01[1]*r02[1] + r01[2]*r02[2]);
  double r0103 = (r01[0]*r03[0] + r01[1]*r03[1] + r01[2]*r03[2]);
  double r0203 = (r02[0]*r03[0] + r02[1]*r03[1] + r02[2]*r03[2]);

  double quad1_0101 = (invmass0+invmass1)*(invmass0+invmass1) * r01sq;
  double quad1_0202 = invmass0*invmass0 * r02sq;
  double quad1_0303 = invmass0*invmass0 * r03sq;
  double quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102;
  double quad1_0103 = 2.0 * (invmass0+invmass1)*invmass0 * r0103;
  double quad1_0203 = 2.0 * invmass0*invmass0 * r0203;

  double quad2_0101 = invmass0*invmass0 * r01sq;
  double quad2_0202 = (invmass0+invmass2)*(invmass0+invmass2) * r02sq;
  double quad2_0303 = invmass0*invmass0 * r03sq;
  double quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102;
  double quad2_0103 = 2.0 * invmass0*invmass0 * r0103;
  double quad2_0203 = 2.0 * (invmass0+invmass2)*invmass0 * r0203;

  double quad3_0101 = invmass0*invmass0 * r01sq;
  double quad3_0202 = invmass0*invmass0 * r02sq;
  double quad3_0303 = (invmass0+invmass3)*(invmass0+invmass3) * r03sq;
  double quad3_0102 = 2.0 * invmass0*invmass0 * r0102;
  double quad3_0103 = 2.0 * (invmass0+invmass3)*invmass0 * r0103;
  double quad3_0203 = 2.0 * (invmass0+invmass3)*invmass0 * r0203;

  // iterate until converged

  double lamda01 = 0.0;
  double lamda02 = 0.0;
  double lamda03 = 0.0;
  int niter = 0;
  int done = 0;

  double quad1,quad2,quad3,b1,b2,b3,lamda01_new,lamda02_new,lamda03_new;

  while (!done && niter < max_iter) {
    quad1 = quad1_0101 * lamda01*lamda01 +
      quad1_0202 * lamda02*lamda02 +
      quad1_0303 * lamda03*lamda03 +
      quad1_0102 * lamda01*lamda02 +
      quad1_0103 * lamda01*lamda03 +
      quad1_0203 * lamda02*lamda03;

    quad2 = quad2_0101 * lamda01*lamda01 +
      quad2_0202 * lamda02*lamda02 +
      quad2_0303 * lamda03*lamda03 +
      quad2_0102 * lamda01*lamda02 +
      quad2_0103 * lamda01*lamda03 +
      quad2_0203 * lamda02*lamda03;

    quad3 = quad3_0101 * lamda01*lamda01 +
      quad3_0202 * lamda02*lamda02 +
      quad3_0303 * lamda03*lamda03 +
      quad3_0102 * lamda01*lamda02 +
      quad3_0103 * lamda01*lamda03 +
      quad3_0203 * lamda02*lamda03;

    b1 = bond1*bond1 - s01sq - quad1;
    b2 = bond2*bond2 - s02sq - quad2;
    b3 = bond3*bond3 - s03sq - quad3;

    lamda01_new = a11inv*b1 + a12inv*b2 + a13inv*b3;
    lamda02_new = a21inv*b1 + a22inv*b2 + a23inv*b3;
    lamda03_new = a31inv*b1 + a32inv*b2 + a33inv*b3;

    done = 1;
    if (fabs(lamda01_new-lamda01) > tolerance) done = 0;
    if (fabs(lamda02_new-lamda02) > tolerance) done = 0;
    if (fabs(lamda03_new-lamda03) > tolerance) done = 0;

    lamda01 = lamda01_new;
    lamda02 = lamda02_new;
    lamda03 = lamda03_new;

    // stop iterations before we have a floating point overflow
    // max double is < 1.0e308, so 1e150 is a reasonable cutoff

    if (fabs(lamda01) > 1e150 || fabs(lamda02) > 1e150
        || fabs(lamda03) > 1e150) done = 1;

    niter++;
  }

  // update forces if atom is owned by this processor

  lamda01 = lamda01/dtfsq;
  lamda02 = lamda02/dtfsq;
  lamda03 = lamda03/dtfsq;

  if (i0 < nlocal) {
    f[i0][0] += lamda01*r01[0] + lamda02*r02[0] + lamda03*r03[0];
    f[i0][1] += lamda01*r01[1] + lamda02*r02[1] + lamda03*r03[1];
    f[i0][2] += lamda01*r01[2] + lamda02*r02[2] + lamda03*r03[2];
  }

  if (i1 < nlocal) {
    f[i1][0] -= lamda01*r01[0];
    f[i1][1] -= lamda01*r01[1];
    f[i1][2] -= lamda01*r01[2];
  }

  if (i2 < nlocal) {
    f[i2][0] -= lamda02*r02[0];
    f[i2][1] -= lamda02*r02[1];
    f[i2][2] -= lamda02*r02[2];
  }

  if (i3 < nlocal) {
    f[i3][0] -= lamda03*r03[0];
    f[i3][1] -= lamda03*r03[1];
    f[i3][2] -= lamda03*r03[2];
  }

  if (evflag) {
    int count = 0;
    if (i0 < nlocal) atomlist[count++] = i0;
    if (i1 < nlocal) atomlist[count++] = i1;
    if (i2 < nlocal) atomlist[count++] = i2;
    if (i3 < nlocal) atomlist[count++] = i3;

    v[0] = lamda01*r01[0]*r01[0]+lamda02*r02[0]*r02[0]+lamda03*r03[0]*r03[0];
    v[1] = lamda01*r01[1]*r01[1]+lamda02*r02[1]*r02[1]+lamda03*r03[1]*r03[1];
    v[2] = lamda01*r01[2]*r01[2]+lamda02*r02[2]*r02[2]+lamda03*r03[2]*r03[2];
    v[3] = lamda01*r01[0]*r01[1]+lamda02*r02[0]*r02[1]+lamda03*r03[0]*r03[1];
    v[4] = lamda01*r01[0]*r01[2]+lamda02*r02[0]*r02[2]+lamda03*r03[0]*r03[2];
    v[5] = lamda01*r01[1]*r01[2]+lamda02*r02[1]*r02[2]+lamda03*r03[1]*r03[2];

    double fpairlist[] = {lamda01, lamda02, lamda03};
    double dellist[][3]  = {{r01[0], r01[1], r01[2]},
                            {r02[0], r02[1], r02[2]},
                            {r03[0], r03[1], r03[2]}};
    int pairlist[][2] = {{i0,i1}, {i0,i2}, {i0,i3}};
    v_tally(count,atomlist,4.0,v,nlocal,3,pairlist,fpairlist,dellist);
  }
}

/* ----------------------------------------------------------------------
   calculate SHAKE constraint forces for size 3 cluster = two bonds + angle
------------------------------------------------------------------------- */

void FixShake::shake3angle(int ilist)
{
  int atomlist[3];
  double v[6];
  double invmass0,invmass1,invmass2;

  // local atom IDs and constraint distances

  int m = list[ilist];
  int i0 = closest_list[ilist][0];
  int i1 = closest_list[ilist][1];
  int i2 = closest_list[ilist][2];
  double bond1 = bond_distance[shake_type[m][0]];
  double bond2 = bond_distance[shake_type[m][1]];
  double bond12 = angle_distance[shake_type[m][2]];

  // r01,r02,r12 = distance vec between atoms

  double r01[3];
  r01[0] = x[i0][0] - x[i1][0];
  r01[1] = x[i0][1] - x[i1][1];
  r01[2] = x[i0][2] - x[i1][2];

  double r02[3];
  r02[0] = x[i0][0] - x[i2][0];
  r02[1] = x[i0][1] - x[i2][1];
  r02[2] = x[i0][2] - x[i2][2];

  double r12[3];
  r12[0] = x[i1][0] - x[i2][0];
  r12[1] = x[i1][1] - x[i2][1];
  r12[2] = x[i1][2] - x[i2][2];

  // s01,s02,s12 = distance vec after unconstrained update

  double s01[3];
  s01[0] = xshake[i0][0] - xshake[i1][0];
  s01[1] = xshake[i0][1] - xshake[i1][1];
  s01[2] = xshake[i0][2] - xshake[i1][2];

  double s02[3];
  s02[0] = xshake[i0][0] - xshake[i2][0];
  s02[1] = xshake[i0][1] - xshake[i2][1];
  s02[2] = xshake[i0][2] - xshake[i2][2];

  double s12[3];
  s12[0] = xshake[i1][0] - xshake[i2][0];
  s12[1] = xshake[i1][1] - xshake[i2][1];
  s12[2] = xshake[i1][2] - xshake[i2][2];

  // scalar distances between atoms

  double r01sq = r01[0]*r01[0] + r01[1]*r01[1] + r01[2]*r01[2];
  double r02sq = r02[0]*r02[0] + r02[1]*r02[1] + r02[2]*r02[2];
  double r12sq = r12[0]*r12[0] + r12[1]*r12[1] + r12[2]*r12[2];
  double s01sq = s01[0]*s01[0] + s01[1]*s01[1] + s01[2]*s01[2];
  double s02sq = s02[0]*s02[0] + s02[1]*s02[1] + s02[2]*s02[2];
  double s12sq = s12[0]*s12[0] + s12[1]*s12[1] + s12[2]*s12[2];

  // matrix coeffs and rhs for lamda equations

  if (rmass) {
    invmass0 = 1.0/rmass[i0];
    invmass1 = 1.0/rmass[i1];
    invmass2 = 1.0/rmass[i2];
  } else {
    invmass0 = 1.0/mass[type[i0]];
    invmass1 = 1.0/mass[type[i1]];
    invmass2 = 1.0/mass[type[i2]];
  }

  double a11 = 2.0 * (invmass0+invmass1) *
    (s01[0]*r01[0] + s01[1]*r01[1] + s01[2]*r01[2]);
  double a12 = 2.0 * invmass0 *
    (s01[0]*r02[0] + s01[1]*r02[1] + s01[2]*r02[2]);
  double a13 = - 2.0 * invmass1 *
    (s01[0]*r12[0] + s01[1]*r12[1] + s01[2]*r12[2]);
  double a21 = 2.0 * invmass0 *
    (s02[0]*r01[0] + s02[1]*r01[1] + s02[2]*r01[2]);
  double a22 = 2.0 * (invmass0+invmass2) *
    (s02[0]*r02[0] + s02[1]*r02[1] + s02[2]*r02[2]);
  double a23 = 2.0 * invmass2 *
    (s02[0]*r12[0] + s02[1]*r12[1] + s02[2]*r12[2]);
  double a31 = - 2.0 * invmass1 *
    (s12[0]*r01[0] + s12[1]*r01[1] + s12[2]*r01[2]);
  double a32 = 2.0 * invmass2 *
    (s12[0]*r02[0] + s12[1]*r02[1] + s12[2]*r02[2]);
  double a33 = 2.0 * (invmass1+invmass2) *
    (s12[0]*r12[0] + s12[1]*r12[1] + s12[2]*r12[2]);

  // inverse of matrix

  double determ = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 -
    a11*a23*a32 - a12*a21*a33 - a13*a22*a31;
  if (determ == 0.0) error->one(FLERR,"Shake determinant = 0.0");
  double determinv = 1.0/determ;

  double a11inv = determinv * (a22*a33 - a23*a32);
  double a12inv = -determinv * (a12*a33 - a13*a32);
  double a13inv = determinv * (a12*a23 - a13*a22);
  double a21inv = -determinv * (a21*a33 - a23*a31);
  double a22inv = determinv * (a11*a33 - a13*a31);
  double a23inv = -determinv * (a11*a23 - a13*a21);
  double a31inv = determinv * (a21*a32 - a22*a31);
  double a32inv = -determinv * (a11*a32 - a12*a31);
  double a33inv = determinv * (a11*a22 - a12*a21);

  // quadratic correction coeffs

  double r0102 = (r01[0]*r02[0] + r01[1]*r02[1] + r01[2]*r02[2]);
  double r0112 = (r01[0]*r12[0] + r01[1]*r12[1] + r01[2]*r12[2]);
  double r0212 = (r02[0]*r12[0] + r02[1]*r12[1] + r02[2]*r12[2]);

  double quad1_0101 = (invmass0+invmass1)*(invmass0+invmass1) * r01sq;
  double quad1_0202 = invmass0*invmass0 * r02sq;
  double quad1_1212 = invmass1*invmass1 * r12sq;
  double quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102;
  double quad1_0112 = - 2.0 * (invmass0+invmass1)*invmass1 * r0112;
  double quad1_0212 = - 2.0 * invmass0*invmass1 * r0212;

  double quad2_0101 = invmass0*invmass0 * r01sq;
  double quad2_0202 = (invmass0+invmass2)*(invmass0+invmass2) * r02sq;
  double quad2_1212 = invmass2*invmass2 * r12sq;
  double quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102;
  double quad2_0112 = 2.0 * invmass0*invmass2 * r0112;
  double quad2_0212 = 2.0 * (invmass0+invmass2)*invmass2 * r0212;

  double quad3_0101 = invmass1*invmass1 * r01sq;
  double quad3_0202 = invmass2*invmass2 * r02sq;
  double quad3_1212 = (invmass1+invmass2)*(invmass1+invmass2) * r12sq;
  double quad3_0102 = - 2.0 * invmass1*invmass2 * r0102;
  double quad3_0112 = - 2.0 * (invmass1+invmass2)*invmass1 * r0112;
  double quad3_0212 = 2.0 * (invmass1+invmass2)*invmass2 * r0212;

  // iterate until converged

  double lamda01 = 0.0;
  double lamda02 = 0.0;
  double lamda12 = 0.0;
  int niter = 0;
  int done = 0;

  double quad1,quad2,quad3,b1,b2,b3,lamda01_new,lamda02_new,lamda12_new;

  while (!done && niter < max_iter) {

    quad1 = quad1_0101 * lamda01*lamda01 +
      quad1_0202 * lamda02*lamda02 +
      quad1_1212 * lamda12*lamda12 +
      quad1_0102 * lamda01*lamda02 +
      quad1_0112 * lamda01*lamda12 +
      quad1_0212 * lamda02*lamda12;

    quad2 = quad2_0101 * lamda01*lamda01 +
      quad2_0202 * lamda02*lamda02 +
      quad2_1212 * lamda12*lamda12 +
      quad2_0102 * lamda01*lamda02 +
      quad2_0112 * lamda01*lamda12 +
      quad2_0212 * lamda02*lamda12;

    quad3 = quad3_0101 * lamda01*lamda01 +
      quad3_0202 * lamda02*lamda02 +
      quad3_1212 * lamda12*lamda12 +
      quad3_0102 * lamda01*lamda02 +
      quad3_0112 * lamda01*lamda12 +
      quad3_0212 * lamda02*lamda12;

    b1 = bond1*bond1 - s01sq - quad1;
    b2 = bond2*bond2 - s02sq - quad2;
    b3 = bond12*bond12 - s12sq - quad3;

    lamda01_new = a11inv*b1 + a12inv*b2 + a13inv*b3;
    lamda02_new = a21inv*b1 + a22inv*b2 + a23inv*b3;
    lamda12_new = a31inv*b1 + a32inv*b2 + a33inv*b3;

    done = 1;
    if (fabs(lamda01_new-lamda01) > tolerance) done = 0;
    if (fabs(lamda02_new-lamda02) > tolerance) done = 0;
    if (fabs(lamda12_new-lamda12) > tolerance) done = 0;

    lamda01 = lamda01_new;
    lamda02 = lamda02_new;
    lamda12 = lamda12_new;

    // stop iterations before we have a floating point overflow
    // max double is < 1.0e308, so 1e150 is a reasonable cutoff

    if (fabs(lamda01) > 1e150 || fabs(lamda02) > 1e150
        || fabs(lamda12) > 1e150) done = 1;

    niter++;
  }

  // update forces if atom is owned by this processor

  lamda01 = lamda01/dtfsq;
  lamda02 = lamda02/dtfsq;
  lamda12 = lamda12/dtfsq;

  if (i0 < nlocal) {
    f[i0][0] += lamda01*r01[0] + lamda02*r02[0];
    f[i0][1] += lamda01*r01[1] + lamda02*r02[1];
    f[i0][2] += lamda01*r01[2] + lamda02*r02[2];
  }

  if (i1 < nlocal) {
    f[i1][0] -= lamda01*r01[0] - lamda12*r12[0];
    f[i1][1] -= lamda01*r01[1] - lamda12*r12[1];
    f[i1][2] -= lamda01*r01[2] - lamda12*r12[2];
  }

  if (i2 < nlocal) {
    f[i2][0] -= lamda02*r02[0] + lamda12*r12[0];
    f[i2][1] -= lamda02*r02[1] + lamda12*r12[1];
    f[i2][2] -= lamda02*r02[2] + lamda12*r12[2];
  }

  if (evflag) {
    int count = 0;
    if (i0 < nlocal) atomlist[count++] = i0;
    if (i1 < nlocal) atomlist[count++] = i1;
    if (i2 < nlocal) atomlist[count++] = i2;

    v[0] = lamda01*r01[0]*r01[0]+lamda02*r02[0]*r02[0]+lamda12*r12[0]*r12[0];
    v[1] = lamda01*r01[1]*r01[1]+lamda02*r02[1]*r02[1]+lamda12*r12[1]*r12[1];
    v[2] = lamda01*r01[2]*r01[2]+lamda02*r02[2]*r02[2]+lamda12*r12[2]*r12[2];
    v[3] = lamda01*r01[0]*r01[1]+lamda02*r02[0]*r02[1]+lamda12*r12[0]*r12[1];
    v[4] = lamda01*r01[0]*r01[2]+lamda02*r02[0]*r02[2]+lamda12*r12[0]*r12[2];
    v[5] = lamda01*r01[1]*r01[2]+lamda02*r02[1]*r02[2]+lamda12*r12[1]*r12[2];

    double fpairlist[] = {lamda01, lamda02, lamda12};
    double dellist[][3]  = {{r01[0], r01[1], r01[2]},
                            {r02[0], r02[1], r02[2]},
                            {r12[0], r12[1], r12[2]}};
    int pairlist[][2] = {{i0,i1}, {i0,i2}, {i1,i2}};
    v_tally(count,atomlist,3.0,v,nlocal,3,pairlist,fpairlist,dellist);
  }
}

/* ----------------------------------------------------------------------
   apply bond force for minimization between atom indices i1 and i2
------------------------------------------------------------------------- */

void FixShake::bond_force(int i1, int i2, double length)
{
  if ((i1 < 0) || (i2 < 0)) return;

  // distance vec between atoms, with PBC

  double delx = x[i1][0] - x[i2][0];
  double dely = x[i1][1] - x[i2][1];
  double delz = x[i1][2] - x[i2][2];

  // compute and apply force

  const double r = sqrt(delx * delx + dely * dely + delz * delz);
  const double dr = r - length;
  const double rk = kbond * dr;
  const double fbond = (r > 0.0) ? -2.0 * rk / r : 0.0;
  const double eb = rk*dr;
  int atomlist[2];
  int count = 0;

  if (i1 < nlocal) {
    f[i1][0] += delx * fbond;
    f[i1][1] += dely * fbond;
    f[i1][2] += delz * fbond;
    atomlist[count++] = i1;
    ebond += 0.5*eb;
  }
  if (i2 < nlocal) {
    f[i2][0] -= delx * fbond;
    f[i2][1] -= dely * fbond;
    f[i2][2] -= delz * fbond;
    atomlist[count++] = i2;
    ebond += 0.5*eb;
  }
  if (evflag) {
    double v[6];
    v[0] = 0.5 * delx * delx * fbond;
    v[1] = 0.5 * dely * dely * fbond;
    v[2] = 0.5 * delz * delz * fbond;
    v[3] = 0.5 * delx * dely * fbond;
    v[4] = 0.5 * delx * delz * fbond;
    v[5] = 0.5 * dely * delz * fbond;
    ev_tally(count, atomlist, 2.0, eb, v);
  }
}

/* ----------------------------------------------------------------------
   print-out bond & angle statistics
------------------------------------------------------------------------- */

void FixShake::stats()
{
  double delx,dely,delz;
  double r,r1,r2,r3,angle;

  // zero out accumulators

  int nb = atom->nbondtypes + 1;
  int na = atom->nangletypes + 1;

  for (int i = 0; i < nb; i++) {
    b_count[i] = 0;
    b_ave[i] = b_max[i] = 0.0;
    b_min[i] = BIG;
  }
  for (int i = 0; i < na; i++) {
    a_count[i] = 0;
    a_ave[i] = a_max[i] = 0.0;
    a_min[i] = BIG;
  }

  // log stats for each bond & angle
  // OK to double count since are just averaging

  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int ii = 0; ii < nlist; ++ii) {
    int i = list[ii];
    int n = shake_flag[i];
    if (n == 0) continue;

    // bond stats

    if (n == 1) n = 3;
    int iatom = closest_list[ii][0];
    for (int j = 1; j < n; j++) {
      int jatom = closest_list[ii][j];
      if (jatom >= nlocal) continue;
      delx = x[iatom][0] - x[jatom][0];
      dely = x[iatom][1] - x[jatom][1];
      delz = x[iatom][2] - x[jatom][2];

      r = sqrt(delx*delx + dely*dely + delz*delz);
      int m = shake_type[i][j-1];
      b_count[m]++;
      b_ave[m] += r;
      b_max[m] = MAX(b_max[m],r);
      b_min[m] = MIN(b_min[m],r);
    }

    // angle stats

    if (shake_flag[i] == 1) {
      int iatom = closest_list[ii][0];
      int jatom = closest_list[ii][1];
      int katom = closest_list[ii][2];
      int n = 0;
      if (iatom < nlocal) ++n;
      if (jatom < nlocal) ++n;
      if (katom < nlocal) ++n;

      delx = x[iatom][0] - x[jatom][0];
      dely = x[iatom][1] - x[jatom][1];
      delz = x[iatom][2] - x[jatom][2];
      r1 = sqrt(delx*delx + dely*dely + delz*delz);

      delx = x[iatom][0] - x[katom][0];
      dely = x[iatom][1] - x[katom][1];
      delz = x[iatom][2] - x[katom][2];
      r2 = sqrt(delx*delx + dely*dely + delz*delz);

      delx = x[jatom][0] - x[katom][0];
      dely = x[jatom][1] - x[katom][1];
      delz = x[jatom][2] - x[katom][2];
      r3 = sqrt(delx*delx + dely*dely + delz*delz);

      angle = acos((r1*r1 + r2*r2 - r3*r3) / (2.0*r1*r2));
      angle *= 180.0/MY_PI;
      int m = shake_type[i][2];
      a_count[m] += n;
      a_ave[m] += n*angle;
      a_max[m] = MAX(a_max[m],angle);
      a_min[m] = MIN(a_min[m],angle);
    }
  }

  // sum across all procs

  MPI_Allreduce(b_count,b_count_all,nb,MPI_LMP_BIGINT,MPI_SUM,world);
  MPI_Allreduce(b_ave,b_ave_all,nb,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(b_max,b_max_all,nb,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(b_min,b_min_all,nb,MPI_DOUBLE,MPI_MIN,world);

  MPI_Allreduce(a_count,a_count_all,na,MPI_LMP_BIGINT,MPI_SUM,world);
  MPI_Allreduce(a_ave,a_ave_all,na,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(a_max,a_max_all,na,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(a_min,a_min_all,na,MPI_DOUBLE,MPI_MIN,world);

  // print stats only for non-zero counts

  if (comm->me == 0) {
    const int width = log10((double)(MAX(MAX(1,nb),na)))+2;
    auto mesg = fmt::format("{} stats (type/ave/delta/count) on step {}\n",
                            utils::uppercase(style), update->ntimestep);
    for (int i = 1; i < nb; i++) {
      const auto bcnt = b_count_all[i];
      if (bcnt)
        mesg += fmt::format("Bond:  {:>{}d}   {:<9.6} {:<11.6} {:>8d}\n",i,width,
                            b_ave_all[i]/bcnt,b_max_all[i]-b_min_all[i],bcnt);
    }
    for (int i = 1; i < na; i++) {
      const auto acnt = a_count_all[i];
      if (acnt)
        mesg += fmt::format("Angle: {:>{}d}   {:<9.6} {:<11.6} {:>8d}\n",i,width,
                            a_ave_all[i]/acnt,a_max_all[i]-a_min_all[i],acnt/3);
    }
    utils::logmesg(lmp,mesg);
  }

  // next timestep for stats

  next_output += output_every;
}

/* ----------------------------------------------------------------------
   find a bond between global atom IDs n1 and n2 stored with local atom i
   if find it:
     if setflag = 0, return bond type
     if setflag = -1/1, set bond type to negative/positive and return 0
   if do not find it, return 0
------------------------------------------------------------------------- */

int FixShake::bondtype_findset(int i, tagint n1, tagint n2, int setflag)
{
  int m,nbonds;
  int *btype;

  if (molecular == Atom::MOLECULAR) {
    tagint *tag = atom->tag;
    tagint **bond_atom = atom->bond_atom;
    nbonds = atom->num_bond[i];

    for (m = 0; m < nbonds; m++) {
      if (n1 == tag[i] && n2 == bond_atom[i][m]) break;
      if (n1 == bond_atom[i][m] && n2 == tag[i]) break;
    }

  } else {
    int imol = atom->molindex[i];
    int iatom = atom->molatom[i];
    tagint *tag = atom->tag;
    tagint tagprev = tag[i] - iatom - 1;
    tagint *batom = atommols[imol]->bond_atom[iatom];
    btype = atommols[imol]->bond_type[iatom];
    nbonds = atommols[imol]->num_bond[iatom];

    for (m = 0; m < nbonds; m++) {
      if (n1 == tag[i] && n2 == batom[m]+tagprev) break;
      if (n1 == batom[m]+tagprev && n2 == tag[i]) break;
    }
  }

  if (m < nbonds) {
    if (setflag == 0) {
      if (molecular == Atom::MOLECULAR) return atom->bond_type[i][m];
      else return btype[m];
    }
    if (molecular == Atom::MOLECULAR) {
      if ((setflag < 0 && atom->bond_type[i][m] > 0) ||
          (setflag > 0 && atom->bond_type[i][m] < 0))
        atom->bond_type[i][m] = -atom->bond_type[i][m];
    } else {
      if ((setflag < 0 && btype[m] > 0) ||
          (setflag > 0 && btype[m] < 0)) btype[m] = -btype[m];
    }
  }

  return 0;
}

/* ----------------------------------------------------------------------
   find an angle with global end atom IDs n1 and n2 stored with local atom i
   if find it:
     if setflag = 0, return angle type
     if setflag = -1/1, set angle type to negative/positive and return 0
   if do not find it, return 0
------------------------------------------------------------------------- */

int FixShake::angletype_findset(int i, tagint n1, tagint n2, int setflag)
{
  int m,nangles;
  int *atype;

  if (molecular == Atom::MOLECULAR) {
    tagint **angle_atom1 = atom->angle_atom1;
    tagint **angle_atom3 = atom->angle_atom3;
    nangles = atom->num_angle[i];

    for (m = 0; m < nangles; m++) {
      if (n1 == angle_atom1[i][m] && n2 == angle_atom3[i][m]) break;
      if (n1 == angle_atom3[i][m] && n2 == angle_atom1[i][m]) break;
    }

  } else {
    int imol = atom->molindex[i];
    int iatom = atom->molatom[i];
    tagint *tag = atom->tag;
    tagint tagprev = tag[i] - iatom - 1;
    tagint *aatom1 = atommols[imol]->angle_atom1[iatom];
    tagint *aatom3 = atommols[imol]->angle_atom3[iatom];
    atype = atommols[imol]->angle_type[iatom];
    nangles = atommols[imol]->num_angle[iatom];

    for (m = 0; m < nangles; m++) {
      if (n1 == aatom1[m]+tagprev && n2 == aatom3[m]+tagprev) break;
      if (n1 == aatom3[m]+tagprev && n2 == aatom1[m]+tagprev) break;
    }
  }

  if (m < nangles) {
    if (setflag == 0) {
      if (molecular == Atom::MOLECULAR) return atom->angle_type[i][m];
      else return atype[m];
    }
    if (molecular == Atom::MOLECULAR) {
      if ((setflag < 0 && atom->angle_type[i][m] > 0) ||
          (setflag > 0 && atom->angle_type[i][m] < 0))
        atom->angle_type[i][m] = -atom->angle_type[i][m];
    } else {
      if ((setflag < 0 && atype[m] > 0) ||
          (setflag > 0 && atype[m] < 0)) atype[m] = -atype[m];
    }
  }

  return 0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixShake::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = (double)nmax * sizeof(int);
  bytes += (double)nmax*4 * sizeof(int);
  bytes += (double)nmax*3 * sizeof(int);
  bytes += (double)nmax*3 * sizeof(double);
  bytes += (double)maxvatom*6 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixShake::grow_arrays(int nmax)
{
  memory->grow(shake_flag,nmax,"shake:shake_flag");
  memory->grow(shake_atom,nmax,4,"shake:shake_atom");
  memory->grow(shake_type,nmax,3,"shake:shake_type");
  memory->destroy(xshake);
  memory->create(xshake,nmax,3,"shake:xshake");
  memory->destroy(ftmp);
  memory->create(ftmp,nmax,3,"shake:ftmp");
  memory->destroy(vtmp);
  memory->create(vtmp,nmax,3,"shake:vtmp");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixShake::copy_arrays(int i, int j, int /*delflag*/)
{
  int flag = shake_flag[j] = shake_flag[i];
  if (flag == 1) {
    shake_atom[j][0] = shake_atom[i][0];
    shake_atom[j][1] = shake_atom[i][1];
    shake_atom[j][2] = shake_atom[i][2];
    shake_type[j][0] = shake_type[i][0];
    shake_type[j][1] = shake_type[i][1];
    shake_type[j][2] = shake_type[i][2];
  } else if (flag == 2) {
    shake_atom[j][0] = shake_atom[i][0];
    shake_atom[j][1] = shake_atom[i][1];
    shake_type[j][0] = shake_type[i][0];
  } else if (flag == 3) {
    shake_atom[j][0] = shake_atom[i][0];
    shake_atom[j][1] = shake_atom[i][1];
    shake_atom[j][2] = shake_atom[i][2];
    shake_type[j][0] = shake_type[i][0];
    shake_type[j][1] = shake_type[i][1];
  } else if (flag == 4) {
    shake_atom[j][0] = shake_atom[i][0];
    shake_atom[j][1] = shake_atom[i][1];
    shake_atom[j][2] = shake_atom[i][2];
    shake_atom[j][3] = shake_atom[i][3];
    shake_type[j][0] = shake_type[i][0];
    shake_type[j][1] = shake_type[i][1];
    shake_type[j][2] = shake_type[i][2];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixShake::set_arrays(int i)
{
  shake_flag[i] = 0;
}

/* ----------------------------------------------------------------------
   update one atom's array values
   called when molecule is created from fix gcmc
------------------------------------------------------------------------- */

void FixShake::update_arrays(int i, int atom_offset)
{
  int flag = shake_flag[i];

  if (flag == 1) {
    shake_atom[i][0] += atom_offset;
    shake_atom[i][1] += atom_offset;
    shake_atom[i][2] += atom_offset;
  } else if (flag == 2) {
    shake_atom[i][0] += atom_offset;
    shake_atom[i][1] += atom_offset;
  } else if (flag == 3) {
    shake_atom[i][0] += atom_offset;
    shake_atom[i][1] += atom_offset;
    shake_atom[i][2] += atom_offset;
  } else if (flag == 4) {
    shake_atom[i][0] += atom_offset;
    shake_atom[i][1] += atom_offset;
    shake_atom[i][2] += atom_offset;
    shake_atom[i][3] += atom_offset;
  }
}

/* ----------------------------------------------------------------------
   initialize a molecule inserted by another fix, e.g. deposit or pour
   called when molecule is created
   nlocalprev = # of atoms on this proc before molecule inserted
   tagprev = atom ID previous to new atoms in the molecule
   xgeom,vcm,quat ignored
------------------------------------------------------------------------- */

void FixShake::set_molecule(int nlocalprev, tagint tagprev, int imol,
                            double * /*xgeom*/, double * /*vcm*/, double * /*quat*/)
{
  int m,flag;

  int nlocal = atom->nlocal;
  if (nlocalprev == nlocal) return;

  tagint *tag = atom->tag;
  tagint **mol_shake_atom = onemols[imol]->shake_atom;
  int **mol_shake_type = onemols[imol]->shake_type;

  for (int i = nlocalprev; i < nlocal; i++) {
    m = tag[i] - tagprev-1;

    flag = shake_flag[i] = onemols[imol]->shake_flag[m];

    if (flag == 1) {
      shake_atom[i][0] = mol_shake_atom[m][0] + tagprev;
      shake_atom[i][1] = mol_shake_atom[m][1] + tagprev;
      shake_atom[i][2] = mol_shake_atom[m][2] + tagprev;
      shake_type[i][0] = mol_shake_type[m][0];
      shake_type[i][1] = mol_shake_type[m][1];
      shake_type[i][2] = mol_shake_type[m][2];
    } else if (flag == 2) {
      shake_atom[i][0] = mol_shake_atom[m][0] + tagprev;
      shake_atom[i][1] = mol_shake_atom[m][1] + tagprev;
      shake_type[i][0] = mol_shake_type[m][0];
    } else if (flag == 3) {
      shake_atom[i][0] = mol_shake_atom[m][0] + tagprev;
      shake_atom[i][1] = mol_shake_atom[m][1] + tagprev;
      shake_atom[i][2] = mol_shake_atom[m][2] + tagprev;
      shake_type[i][0] = mol_shake_type[m][0];
      shake_type[i][1] = mol_shake_type[m][1];
    } else if (flag == 4) {
      shake_atom[i][0] = mol_shake_atom[m][0] + tagprev;
      shake_atom[i][1] = mol_shake_atom[m][1] + tagprev;
      shake_atom[i][2] = mol_shake_atom[m][2] + tagprev;
      shake_atom[i][3] = mol_shake_atom[m][3] + tagprev;
      shake_type[i][0] = mol_shake_type[m][0];
      shake_type[i][1] = mol_shake_type[m][1];
      shake_type[i][2] = mol_shake_type[m][2];
    }
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixShake::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = shake_flag[i];
  int flag = shake_flag[i];
  if (flag == 1) {
    buf[m++] = shake_atom[i][0];
    buf[m++] = shake_atom[i][1];
    buf[m++] = shake_atom[i][2];
    buf[m++] = shake_type[i][0];
    buf[m++] = shake_type[i][1];
    buf[m++] = shake_type[i][2];
  } else if (flag == 2) {
    buf[m++] = shake_atom[i][0];
    buf[m++] = shake_atom[i][1];
    buf[m++] = shake_type[i][0];
  } else if (flag == 3) {
    buf[m++] = shake_atom[i][0];
    buf[m++] = shake_atom[i][1];
    buf[m++] = shake_atom[i][2];
    buf[m++] = shake_type[i][0];
    buf[m++] = shake_type[i][1];
  } else if (flag == 4) {
    buf[m++] = shake_atom[i][0];
    buf[m++] = shake_atom[i][1];
    buf[m++] = shake_atom[i][2];
    buf[m++] = shake_atom[i][3];
    buf[m++] = shake_type[i][0];
    buf[m++] = shake_type[i][1];
    buf[m++] = shake_type[i][2];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixShake::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  int flag = shake_flag[nlocal] = static_cast<int> (buf[m++]);
  if (flag == 1) {
    shake_atom[nlocal][0] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][1] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][2] = static_cast<tagint> (buf[m++]);
    shake_type[nlocal][0] = static_cast<int> (buf[m++]);
    shake_type[nlocal][1] = static_cast<int> (buf[m++]);
    shake_type[nlocal][2] = static_cast<int> (buf[m++]);
  } else if (flag == 2) {
    shake_atom[nlocal][0] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][1] = static_cast<tagint> (buf[m++]);
    shake_type[nlocal][0] = static_cast<int> (buf[m++]);
  } else if (flag == 3) {
    shake_atom[nlocal][0] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][1] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][2] = static_cast<tagint> (buf[m++]);
    shake_type[nlocal][0] = static_cast<int> (buf[m++]);
    shake_type[nlocal][1] = static_cast<int> (buf[m++]);
  } else if (flag == 4) {
    shake_atom[nlocal][0] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][1] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][2] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][3] = static_cast<tagint> (buf[m++]);
    shake_type[nlocal][0] = static_cast<int> (buf[m++]);
    shake_type[nlocal][1] = static_cast<int> (buf[m++]);
    shake_type[nlocal][2] = static_cast<int> (buf[m++]);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int FixShake::pack_forward_comm(int n, int *list, double *buf,
                                int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = xshake[j][0];
      buf[m++] = xshake[j][1];
      buf[m++] = xshake[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = xshake[j][0] + dx;
      buf[m++] = xshake[j][1] + dy;
      buf[m++] = xshake[j][2] + dz;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixShake::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    xshake[i][0] = buf[m++];
    xshake[i][1] = buf[m++];
    xshake[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void FixShake::reset_dt()
{
  if (utils::strmatch(update->integrate_style,"^verlet")) {
    dtv = update->dt;
    if (rattle) dtfsq   = 0.5 * update->dt * update->dt * force->ftm2v;
    else dtfsq = update->dt * update->dt * force->ftm2v;
  } else {
    dtv = step_respa[0];
    dtf_innerhalf = 0.5 * step_respa[0] * force->ftm2v;
    if (rattle) dtf_inner = dtf_innerhalf;
    else dtf_inner = step_respa[0] * force->ftm2v;
  }
}

/* ----------------------------------------------------------------------
   extract Molecule ptr
------------------------------------------------------------------------- */

void *FixShake::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"onemol") == 0) return onemols;
  return nullptr;
}

/* ----------------------------------------------------------------------
   energy due to restraint forces
------------------------------------------------------------------------- */

double FixShake::compute_scalar()
{
  double all;
  MPI_Allreduce(&ebond, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}

/*  ----------------------------------------------------------------------
   print shake stats at the end of a minimization
------------------------------------------------------------------------- */
void FixShake::post_run()
{
  if ((update->whichflag == 2) && (output_every > 0)) stats();
}


/* ----------------------------------------------------------------------
   add coordinate constraining forces
   this method is called at the end of a timestep
------------------------------------------------------------------------- */

void FixShake::shake_end_of_step(int vflag) {

  if (!respa) {
    dtv     = update->dt;
    dtfsq   = 0.5 * update->dt * update->dt * force->ftm2v;
    FixShake::post_force(vflag);
    if (!rattle) dtfsq = update->dt * update->dt * force->ftm2v;

  } else {
    dtv = step_respa[0];
    dtf_innerhalf = 0.5 * step_respa[0] * force->ftm2v;
    dtf_inner = dtf_innerhalf;

    // apply correction to all rRESPA levels

    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel);
      FixShake::post_force_respa(vflag,ilevel,loop_respa[ilevel]-1);
      (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel);
    }
    if (!rattle) dtf_inner = step_respa[0] * force->ftm2v;
  }
}

/* ----------------------------------------------------------------------
   wrapper method for end_of_step fixes which modify velocities
------------------------------------------------------------------------- */

void FixShake::correct_velocities() {}

/* ----------------------------------------------------------------------
   calculate constraining forces based on the current configuration
   change coordinates
------------------------------------------------------------------------- */

void FixShake::correct_coordinates(int vflag) {

  // save current forces and velocities so that you
  // initialize them to zero such that FixShake::unconstrained_coordinate_update has no effect

  for (int j=0; j<nlocal; j++) {
    for (int k=0; k<3; k++) {

      // store current value of forces and velocities

      ftmp[j][k] = f[j][k];
      vtmp[j][k] = v[j][k];

      // set f and v to zero for SHAKE

      v[j][k] = 0;
      f[j][k] = 0;
    }
  }

  // call SHAKE to correct the coordinates which were updated without constraints
  // IMPORTANT: use 1 as argument and thereby enforce velocity Verlet

  dtfsq   = 0.5 * update->dt * update->dt * force->ftm2v;
  FixShake::post_force(vflag);

  // integrate coordinates: x' = xnp1 + dt^2/2m_i * f, where f is the constraining force
  // NOTE: After this command, the coordinates geometry of the molecules will be correct!

  double dtfmsq;
  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      dtfmsq = dtfsq/ rmass[i];
      x[i][0] = x[i][0] + dtfmsq*f[i][0];
      x[i][1] = x[i][1] + dtfmsq*f[i][1];
      x[i][2] = x[i][2] + dtfmsq*f[i][2];
    }
  }
  else {
    for (int i = 0; i < nlocal; i++) {
      dtfmsq = dtfsq / mass[type[i]];
      x[i][0] = x[i][0] + dtfmsq*f[i][0];
      x[i][1] = x[i][1] + dtfmsq*f[i][1];
      x[i][2] = x[i][2] + dtfmsq*f[i][2];
    }
  }

  // copy forces and velocities back

  for (int j=0; j<nlocal; j++) {
    for (int k=0; k<3; k++) {
      f[j][k] = ftmp[j][k];
      v[j][k] = vtmp[j][k];
    }
  }

  if (!rattle) dtfsq = update->dt * update->dt * force->ftm2v;

  // communicate changes
  // NOTE: for compatibility xshake is temporarily set to x, such that pack/unpack_forward
  //       can be used for communicating the coordinates.

  double **xtmp = xshake;
  xshake = x;
  comm->forward_comm(this);
  xshake = xtmp;
}
