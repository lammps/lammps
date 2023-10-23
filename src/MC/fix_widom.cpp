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

/* ----------------------------------------------------------------------
   Contributing author: Evangelos Voyiatzis (Royal DSM)
------------------------------------------------------------------------- */

#include "fix_widom.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "region.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;

#define MAXENERGYTEST 1.0e50
enum { EXCHATOM, EXCHMOL };    // exchmode

/* ---------------------------------------------------------------------- */

FixWidom::FixWidom(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), region(nullptr), idregion(nullptr), full_flag(false), molcoords(nullptr),
    molq(nullptr), molimage(nullptr), random_equal(nullptr), c_pe(nullptr)
{
  if (narg < 8) utils::missing_cmd_args(FLERR, "fix widom", error);

  if (atom->molecular == Atom::TEMPLATE)
    error->all(FLERR, "Fix widom does not (yet) work with atom_style template");

  dynamic_group_allow = 1;

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // required args

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  ninsertions = utils::inumeric(FLERR,arg[4],false,lmp);
  nwidom_type = utils::inumeric(FLERR,arg[5],false,lmp);
  seed = utils::inumeric(FLERR,arg[6],false,lmp);
  insertion_temperature = utils::numeric(FLERR,arg[7],false,lmp);

  if (nevery <= 0) error->all(FLERR,"Invalid fix widom every argument: {}", nevery);
  if (ninsertions < 0) error->all(FLERR,"Invalid fix widom insertions argument: {}", ninsertions);
  if (seed <= 0) error->all(FLERR,"Invalid fix widom seed argument: {}", seed);
  if (insertion_temperature < 0.0)
    error->all(FLERR,"Invalid fix widom temperature argument: {}", insertion_temperature);

  // read options from end of input line

  options(narg-8,&arg[8]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp,seed);

  // error checks on region and its extent being inside simulation box

  region_xlo = region_xhi = region_ylo = region_yhi = region_zlo = region_zhi = 0.0;
  if (region) {
    if (region->bboxflag == 0)
      error->all(FLERR,"Fix widom region {} does not support a bounding box", region->id);
    if (region->dynamic_check())
      error->all(FLERR,"Fix widom region {} cannot be dynamic", region->id);

    region_xlo = region->extent_xlo;
    region_xhi = region->extent_xhi;
    region_ylo = region->extent_ylo;
    region_yhi = region->extent_yhi;
    region_zlo = region->extent_zlo;
    region_zhi = region->extent_zhi;

    if (triclinic) {
      if ((region_xlo < domain->boxlo_bound[0]) || (region_xhi > domain->boxhi_bound[0]) ||
          (region_ylo < domain->boxlo_bound[1]) || (region_yhi > domain->boxhi_bound[1]) ||
          (region_zlo < domain->boxlo_bound[2]) || (region_zhi > domain->boxhi_bound[2]))
        error->all(FLERR,"Fix widom region {} extends outside simulation box", region->id);
    } else {
      if ((region_xlo < domain->boxlo[0]) || (region_xhi > domain->boxhi[0]) ||
          (region_ylo < domain->boxlo[1]) || (region_yhi > domain->boxhi[1]) ||
          (region_zlo < domain->boxlo[2]) || (region_zhi > domain->boxhi[2]))
        error->all(FLERR,"Fix widom region {} extends outside simulation box", region->id);
    }

    // estimate region volume using MC trials

    double coord[3];
    int inside = 0;
    int attempts = 10000000;
    for (int i = 0; i < attempts; i++) {
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      if (region->match(coord[0],coord[1],coord[2]) != 0)
        inside++;
    }

    double max_region_volume = (region_xhi - region_xlo) * (region_yhi - region_ylo)
      * (region_zhi - region_zlo);

    region_volume = max_region_volume * static_cast<double>(inside) / static_cast<double>(attempts);
  }

  // error check and further setup for exchmode = EXCHMOL

  if (exchmode == EXCHMOL) {
    if (onemol->xflag == 0)
      error->all(FLERR,"Fix widom molecule {} must have coordinates", onemol->id);
    if (onemol->typeflag == 0)
      error->all(FLERR,"Fix widom molecule {} must have atom types", onemol->id);
    if (nwidom_type != 0)
      error->all(FLERR,"Atom type must be zero in fix widom mol command");
    if (onemol->qflag == 1 && atom->q == nullptr)
      error->all(FLERR,"Fix widom molecule {} has charges, but atom style does not", onemol->id);

    onemol->check_attributes();
  }

  if (charge_flag && atom->q == nullptr)
    error->all(FLERR,"Fix widom atom has charge, but atom style does not");

  // setup of array of coordinates for molecule insertion

  if (exchmode == EXCHATOM) natoms_per_molecule = 1;
  else natoms_per_molecule = onemol->natoms;
  nmaxmolatoms = natoms_per_molecule;
  grow_molecule_arrays(nmaxmolatoms);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters
  widom_nmax = 0;
  ave_widom_chemical_potential = 0.0;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixWidom::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix widom command");

  // defaults

  exchmode = EXCHATOM;
  region_volume = 0;
  max_region_attempts = 1000;
  molecule_group = 0;
  molecule_group_bit = 0;
  molecule_group_inversebit = 0;
  exclusion_group = 0;
  exclusion_group_bit = 0;
  charge = 0.0;
  charge_flag = false;
  full_flag = false;
  energy_intra = 0.0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix widom mol", error);
      auto onemols = atom->get_molecule_by_id(arg[iarg+1]);
      if (onemols.size() == 0)
        error->all(FLERR,"Molecule template ID {} for fix widom does not exist", arg[iarg+1]);
      if (onemols.size() > 1 && comm->me == 0)
        error->warning(FLERR,"Molecule template {} for fix widom has multiple molecules; "
                       "will use only the first molecule", arg[iarg+1]);
      exchmode = EXCHMOL;
      onemol = onemols[0];
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix widom command");
      region = domain->get_region_by_id(arg[iarg+1]);
      if (!region)
        error->all(FLERR,"Region {} for fix widom does not exist",arg[iarg+1]);
      idregion = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"charge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix widom command");
      charge = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      charge_flag = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"full_energy") == 0) {
      full_flag = true;
      iarg += 1;
    } else if (strcmp(arg[iarg],"intra_energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix widom command");
      energy_intra = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix widom command");
  }
}

/* ---------------------------------------------------------------------- */

FixWidom::~FixWidom()
{
  delete[] idregion;
  delete random_equal;

  memory->destroy(molcoords);
  memory->destroy(molq);
  memory->destroy(molimage);

  // delete exclusion group created in init()
  // delete molecule group created in init()
  // unset neighbor exclusion settings made in init()
  // not necessary if group and neighbor classes already destroyed
  //   when LAMMPS exits

  if (exclusion_group_bit && group) {
    auto group_id = std::string("FixWidom:widom_exclusion_group:") + id;
    try {
      group->assign(group_id + " delete");
    } catch (std::exception &e) {
      if (comm->me == 0)
        fprintf(stderr, "Error deleting group %s: %s\n", group_id.c_str(), e.what());
    }
  }

  if (molecule_group_bit && group) {
    auto group_id = std::string("FixWidom:rotation_gas_atoms:") + id;
    try {
      group->assign(group_id + " delete");
    } catch (std::exception &e) {
      if (comm->me == 0)
        fprintf(stderr, "Error deleting group %s: %s\n", group_id.c_str(), e.what());
    }
  }

  if (full_flag && group && neighbor) {
    int igroupall = group->find("all");
    neighbor->exclusion_group_group_delete(exclusion_group,igroupall);
  }
}

/* ---------------------------------------------------------------------- */

int FixWidom::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWidom::init()
{

  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix widom does not exist", idregion);
  }

  triclinic = domain->triclinic;

  ave_widom_chemical_potential = 0.0;

  if (region) volume = region_volume;
  else volume = domain->xprd * domain->yprd * domain->zprd;

  // decide whether to switch to the full_energy option
  if (!full_flag) {
    if ((force->kspace) ||
        (force->pair == nullptr) ||
        (force->pair->single_enable == 0) ||
        (force->pair_match("^hybrid",0)) ||
        (force->pair_match("^eam",0)) ||
        (force->pair->tail_flag)) {
      full_flag = true;
      if (comm->me == 0)
        error->warning(FLERR,"Fix widom using full_energy option");
    }
  }

  if (full_flag) c_pe = modify->get_compute_by_id("thermo_pe");

  if (exchmode == EXCHATOM) {
    if (nwidom_type <= 0 || nwidom_type > atom->ntypes)
      error->all(FLERR,"Invalid atom type in fix widom command");
  }

  // if molecules are exchanged or moved, check for unset mol IDs
  if (exchmode == EXCHMOL) {
    tagint *molecule = atom->molecule;
    int *mask = atom->mask;
    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if (mask[i] == groupbit)
        if (molecule[i] == 0) flag = 1;
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall && comm->me == 0)
      error->all(FLERR, "All mol IDs should be set for fix widom group atoms");
  }

  if (exchmode == EXCHMOL)
    if (atom->molecule_flag == 0 || !atom->tag_enable
        || (atom->map_style == Atom::MAP_NONE))
      error->all(FLERR, "Fix widom molecule command requires that atoms have molecule attributes");

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix widom in a 2d simulation");

  // create a new group for interaction exclusions
  // used for attempted atom or molecule deletions
  // skip if already exists from previous init()

  if (full_flag && !exclusion_group_bit) {

    // create unique group name for atoms to be excluded

    auto group_id = std::string("FixWidom:widom_exclusion_group:") + id;
    group->assign(group_id + " subtract all all");
    exclusion_group = group->find(group_id);
    if (exclusion_group == -1)
      error->all(FLERR,"Could not find fix widom exclusion group ID");
    exclusion_group_bit = group->bitmask[exclusion_group];

    // neighbor list exclusion setup
    // turn off interactions between group all and the exclusion group

    neighbor->modify_params(fmt::format("exclude group {} all",group_id));
  }

  // create a new group for temporary use with selected molecules

  if (exchmode == EXCHMOL) {

    auto group_id = std::string("FixWidom:rotation_gas_atoms:") + id;
    group->assign(group_id + " molecule -1");
    molecule_group = group->find(group_id);
    if (molecule_group == -1)
      error->all(FLERR,"Could not find fix widom rotation group ID");
    molecule_group_bit = group->bitmask[molecule_group];
    molecule_group_inversebit = molecule_group_bit ^ ~0;
  }

  // get all of the needed molecule data if exchanging
  // or moving molecules, otherwise just get the gas mass

  if (exchmode == EXCHMOL) {

    onemol->compute_mass();
    onemol->compute_com();
    gas_mass = onemol->masstotal;
    for (int i = 0; i < onemol->natoms; i++) {
      onemol->x[i][0] -= onemol->com[0];
      onemol->x[i][1] -= onemol->com[1];
      onemol->x[i][2] -= onemol->com[2];
    }
    onemol->com[0] = 0;
    onemol->com[1] = 0;
    onemol->com[2] = 0;

  } else gas_mass = atom->mass[nwidom_type];

  if (gas_mass <= 0.0) error->all(FLERR,"Illegal fix widom gas mass <= 0");

  // check that no deletable atoms are in atom->firstgroup
  // deleting such an atom would not leave firstgroup atoms first

  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if ((mask[i] == groupbit) && (mask[i] && firstgroupbit)) flag = 1;

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

    if (flagall)
      error->all(FLERR,"Cannot use fix widom on atoms in atom_modify first group");
  }

  // compute beta
  beta = 1.0/(force->boltz*insertion_temperature);

  imagezero = ((imageint) IMGMAX << IMG2BITS) |
             ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  // Current implementation is broken using
  // full_flag on molecules on more than one processor.
  // Print error if this is the current mode
  if (full_flag && (exchmode == EXCHMOL) && comm->nprocs > 1)
    error->all(FLERR,"fix widom does currently not support full_energy option with "
               "molecules on more than 1 MPI process.");

}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixWidom::pre_exchange()
{
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  ave_widom_chemical_potential = 0.0;

  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];
  if (triclinic) {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  } else {
    sublo = domain->sublo;
    subhi = domain->subhi;
  }

  if (region) volume = region_volume;
  else volume = domain->xprd * domain->yprd * domain->zprd;

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  if (full_flag) {
    energy_stored = energy_full();

    if (exchmode == EXCHATOM) {
      attempt_atomic_insertion_full();
    } else {
      attempt_molecule_insertion_full();
    }

    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    atom->nghost = 0;
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  } else {

    if (exchmode == EXCHATOM) {
      attempt_atomic_insertion();
    } else {
      attempt_molecule_insertion();
    }

  }
  next_reneighbor = update->ntimestep + nevery;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixWidom::attempt_atomic_insertion()
{

  double lamda[3];
  double coord[3];

  for (int imove = 0; imove < ninsertions; imove++) {

    // pick coordinates for insertion point

    if (region) {
      int region_attempt = 0;
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      while (region->match(coord[0],coord[1],coord[2]) == 0) {
        coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
        coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
        coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
        region_attempt++;
        if (region_attempt >= max_region_attempts) return;
      }
      if (triclinic) domain->x2lamda(coord,lamda);
    } else {
      if (triclinic == 0) {
        coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
        coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
        coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
      } else {
        lamda[0] = random_equal->uniform();
        lamda[1] = random_equal->uniform();
        lamda[2] = random_equal->uniform();

        // wasteful, but necessary

        if (lamda[0] == 1.0) lamda[0] = 0.0;
        if (lamda[1] == 1.0) lamda[1] = 0.0;
        if (lamda[2] == 1.0) lamda[2] = 0.0;

        domain->lamda2x(lamda,coord);
      }
    }

    int proc_flag = 0;
    if (triclinic == 0) {
      domain->remap(coord);
      if (!domain->inside(coord))
        error->one(FLERR,"Fix widom put atom outside box");
      if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
          coord[1] >= sublo[1] && coord[1] < subhi[1] &&
          coord[2] >= sublo[2] && coord[2] < subhi[2]) proc_flag = 1;
    } else {
      if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
          lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
          lamda[2] >= sublo[2] && lamda[2] < subhi[2]) proc_flag = 1;
    }

    if (proc_flag) {
      int ii = -1;
      if (charge_flag) {
        ii = atom->nlocal + atom->nghost;
        if (ii >= atom->nmax) atom->avec->grow(0);
        atom->q[ii] = charge;
      }
      double insertion_energy = energy(ii,nwidom_type,-1,coord);
      double inst_chem_pot = exp(-insertion_energy*beta);
      double incr_chem_pot = (inst_chem_pot - ave_widom_chemical_potential);
      ave_widom_chemical_potential += incr_chem_pot / (imove + 1);
    }
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixWidom::attempt_molecule_insertion()
{

  double lamda[3];
  double com_coord[3];
  double r[3],rotmat[3][3],quat[4];

  for (int imove = 0; imove < ninsertions; imove++) {

    if (region) {
      int region_attempt = 0;
      com_coord[0] = region_xlo + random_equal->uniform() *
        (region_xhi-region_xlo);
      com_coord[1] = region_ylo + random_equal->uniform() *
        (region_yhi-region_ylo);
      com_coord[2] = region_zlo + random_equal->uniform() *
        (region_zhi-region_zlo);
      while (region->match(com_coord[0],com_coord[1],
                           com_coord[2]) == 0) {
        com_coord[0] = region_xlo + random_equal->uniform() *
          (region_xhi-region_xlo);
        com_coord[1] = region_ylo + random_equal->uniform() *
          (region_yhi-region_ylo);
        com_coord[2] = region_zlo + random_equal->uniform() *
          (region_zhi-region_zlo);
        region_attempt++;
        if (region_attempt >= max_region_attempts) return;
      }
      if (triclinic) domain->x2lamda(com_coord,lamda);
    } else {
      if (triclinic == 0) {
        com_coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
        com_coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
        com_coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
      } else {
        lamda[0] = random_equal->uniform();
        lamda[1] = random_equal->uniform();
        lamda[2] = random_equal->uniform();

        // wasteful, but necessary

        if (lamda[0] == 1.0) lamda[0] = 0.0;
        if (lamda[1] == 1.0) lamda[1] = 0.0;
        if (lamda[2] == 1.0) lamda[2] = 0.0;

        domain->lamda2x(lamda,com_coord);
      }
    }

    // generate point in unit cube
    // then restrict to unit sphere

    double rsq = 1.1;
    while (rsq > 1.0) {
      r[0] = 2.0*random_equal->uniform() - 1.0;
      r[1] = 2.0*random_equal->uniform() - 1.0;
      r[2] = 2.0*random_equal->uniform() - 1.0;
      rsq = MathExtra::dot3(r, r);
    }

    double theta = random_equal->uniform() * MY_2PI;
    MathExtra::norm3(r);
    MathExtra::axisangle_to_quat(r,theta,quat);
    MathExtra::quat_to_mat(quat,rotmat);

    double insertion_energy = 0.0;
    auto procflag = new bool[natoms_per_molecule];

    for (int i = 0; i < natoms_per_molecule; i++) {
      MathExtra::matvec(rotmat,onemol->x[i],molcoords[i]);
      molcoords[i][0] += com_coord[0];
      molcoords[i][1] += com_coord[1];
      molcoords[i][2] += com_coord[2];

      // use temporary variable for remapped position
      // so unmapped position is preserved in molcoords

      double xtmp[3];
      xtmp[0] = molcoords[i][0];
      xtmp[1] = molcoords[i][1];
      xtmp[2] = molcoords[i][2];
      domain->remap(xtmp);
      if (!domain->inside(xtmp))
        error->one(FLERR,"Fix widom put atom outside box");

      procflag[i] = false;
      if (triclinic == 0) {
        if (xtmp[0] >= sublo[0] && xtmp[0] < subhi[0] &&
            xtmp[1] >= sublo[1] && xtmp[1] < subhi[1] &&
            xtmp[2] >= sublo[2] && xtmp[2] < subhi[2]) procflag[i] = true;
      } else {
        domain->x2lamda(xtmp,lamda);
        if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
            lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
            lamda[2] >= sublo[2] && lamda[2] < subhi[2]) procflag[i] = true;
      }

      if (procflag[i]) {
        int ii = -1;
        if (onemol->qflag == 1) {
          ii = atom->nlocal + atom->nghost;
          if (ii >= atom->nmax) atom->avec->grow(0);
          atom->q[ii] = onemol->q[i];
        }
        insertion_energy += energy(ii,onemol->type[i],-1,xtmp);
      }
    }

    double insertion_energy_sum = 0.0;
    MPI_Allreduce(&insertion_energy,&insertion_energy_sum,1,
                  MPI_DOUBLE,MPI_SUM,world);

    // the insertion_energy_sum is the variable with the energy of inserting one molecule
    double inst_chem_pot = exp(-insertion_energy_sum*beta);
    double incr_chem_pot = (inst_chem_pot - ave_widom_chemical_potential);
    ave_widom_chemical_potential += incr_chem_pot / (imove + 1);

    delete[] procflag;
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixWidom::attempt_atomic_insertion_full()
{

  double lamda[3];
  double coord[3];

  for (int imove = 0; imove < ninsertions; imove++) {

    if (region) {
      int region_attempt = 0;
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      while (region->match(coord[0],coord[1],coord[2]) == 0) {
        coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
        coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
        coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
        region_attempt++;
        if (region_attempt >= max_region_attempts) return;
      }
      if (triclinic) domain->x2lamda(coord,lamda);
    } else {
      if (triclinic == 0) {
        coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
        coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
        coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
      } else {
        lamda[0] = random_equal->uniform();
        lamda[1] = random_equal->uniform();
        lamda[2] = random_equal->uniform();

        // wasteful, but necessary

        if (lamda[0] == 1.0) lamda[0] = 0.0;
        if (lamda[1] == 1.0) lamda[1] = 0.0;
        if (lamda[2] == 1.0) lamda[2] = 0.0;

        domain->lamda2x(lamda,coord);
      }
    }

    int proc_flag = 0;
    if (triclinic == 0) {
      domain->remap(coord);
      if (!domain->inside(coord))
        error->one(FLERR,"Fix widom put atom outside box");
      if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
          coord[1] >= sublo[1] && coord[1] < subhi[1] &&
          coord[2] >= sublo[2] && coord[2] < subhi[2]) proc_flag = 1;
    } else {
      if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
          lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
          lamda[2] >= sublo[2] && lamda[2] < subhi[2]) proc_flag = 1;
    }

    if (proc_flag) {
      atom->avec->create_atom(nwidom_type,coord);
      int m = atom->nlocal - 1;

      // add to groups
      // optionally add to type-based groups

      atom->v[m][0] = 0;
      atom->v[m][1] = 0;
      atom->v[m][2] = 0;
      if (charge_flag) atom->q[m] = charge;
      modify->create_attribute(m);
    }

    atom->natoms++;
    if (atom->tag_enable) {
      atom->tag_extend();
      if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    }
    atom->nghost = 0;
    if (triclinic) domain->x2lamda(atom->nlocal);
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();

    double insertion_energy = energy_full() - energy_stored;
    double inst_chem_pot = exp(-insertion_energy*beta);
    double incr_chem_pot = (inst_chem_pot - ave_widom_chemical_potential);
    ave_widom_chemical_potential += incr_chem_pot / (imove + 1);

    atom->natoms--;
    if (proc_flag) atom->nlocal--;
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixWidom::attempt_molecule_insertion_full()
{

  double lamda[3];

  tagint maxmol = 0;
  for (int i = 0; i < atom->nlocal; i++) maxmol = MAX(maxmol,atom->molecule[i]);
  tagint maxmol_all;
  MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
  maxmol_all++;
  if (maxmol_all >= MAXTAGINT)
    error->all(FLERR,"Fix widom ran out of available molecule IDs");
  int insertion_molecule = maxmol_all;

  tagint maxtag = 0;
  for (int i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag,atom->tag[i]);
  tagint maxtag_all;
  MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

  for (int imove = 0; imove < ninsertions; imove++) {
    double com_coord[3];
    if (region) {
      int region_attempt = 0;
      com_coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      com_coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      com_coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      while (region->match(com_coord[0],com_coord[1], com_coord[2]) == 0) {
        com_coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
        com_coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
        com_coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
        region_attempt++;
        if (region_attempt >= max_region_attempts) return;
      }
      if (triclinic) domain->x2lamda(com_coord,lamda);
    } else {
      if (triclinic == 0) {
        com_coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
        com_coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
        com_coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
      } else {
        lamda[0] = random_equal->uniform();
        lamda[1] = random_equal->uniform();
        lamda[2] = random_equal->uniform();

        // wasteful, but necessary

        if (lamda[0] == 1.0) lamda[0] = 0.0;
        if (lamda[1] == 1.0) lamda[1] = 0.0;
        if (lamda[2] == 1.0) lamda[2] = 0.0;

        domain->lamda2x(lamda,com_coord);
      }

    }

    // generate point in unit cube
    // then restrict to unit sphere

    double r[3],rotmat[3][3],quat[4];
    double rsq = 1.1;
    while (rsq > 1.0) {
      r[0] = 2.0*random_equal->uniform() - 1.0;
      r[1] = 2.0*random_equal->uniform() - 1.0;
      r[2] = 2.0*random_equal->uniform() - 1.0;
      rsq = MathExtra::dot3(r, r);
    }

    double theta = random_equal->uniform() * MY_2PI;
    MathExtra::norm3(r);
    MathExtra::axisangle_to_quat(r,theta,quat);
    MathExtra::quat_to_mat(quat,rotmat);

    for (int i = 0; i < natoms_per_molecule; i++) {
      double xtmp[3];
      MathExtra::matvec(rotmat,onemol->x[i],xtmp);
      xtmp[0] += com_coord[0];
      xtmp[1] += com_coord[1];
      xtmp[2] += com_coord[2];

      // need to adjust image flags in remap()

      imageint imagetmp = imagezero;
      domain->remap(xtmp,imagetmp);
      if (!domain->inside(xtmp))
        error->one(FLERR,"Fix widom put atom outside box");

      int proc_flag = 0;
      if (triclinic == 0) {
        if (xtmp[0] >= sublo[0] && xtmp[0] < subhi[0] &&
            xtmp[1] >= sublo[1] && xtmp[1] < subhi[1] &&
            xtmp[2] >= sublo[2] && xtmp[2] < subhi[2]) proc_flag = 1;
      } else {
        domain->x2lamda(xtmp,lamda);
        if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
            lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
            lamda[2] >= sublo[2] && lamda[2] < subhi[2]) proc_flag = 1;
      }

      if (proc_flag) {
        atom->avec->create_atom(onemol->type[i],xtmp);
        int m = atom->nlocal - 1;

        atom->image[m] = imagetmp;
        atom->molecule[m] = insertion_molecule;
        if (maxtag_all+i+1 >= MAXTAGINT)
          error->all(FLERR,"Fix widom ran out of available atom IDs");
        atom->tag[m] = maxtag_all + i + 1;
        atom->v[m][0] = 0;
        atom->v[m][1] = 0;
        atom->v[m][2] = 0;

        atom->add_molecule_atom(onemol,i,m,maxtag_all);
        modify->create_attribute(m);
      }
    }

    atom->natoms += natoms_per_molecule;
    if (atom->natoms < 0) error->all(FLERR,"Too many total atoms");
    atom->nbonds += onemol->nbonds;
    atom->nangles += onemol->nangles;
    atom->ndihedrals += onemol->ndihedrals;
    atom->nimpropers += onemol->nimpropers;
    if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    atom->nghost = 0;
    if (triclinic) domain->x2lamda(atom->nlocal);
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();

    // energy_after corrected by energy_intra
    double insertion_energy = (energy_full() -energy_intra) - energy_stored;
    double inst_chem_pot = exp(-insertion_energy*beta);
    double incr_chem_pot = (inst_chem_pot - ave_widom_chemical_potential);
    ave_widom_chemical_potential += incr_chem_pot / (imove + 1);

    atom->nbonds -= onemol->nbonds;
    atom->nangles -= onemol->nangles;
    atom->ndihedrals -= onemol->ndihedrals;
    atom->nimpropers -= onemol->nimpropers;
    atom->natoms -= natoms_per_molecule;

    int i = 0;
    while (i < atom->nlocal) {
      if (atom->molecule[i] == insertion_molecule) {
        atom->avec->copy(atom->nlocal-1,i,1);
        atom->nlocal--;
      } else i++;
    }
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
  }

}

/* ----------------------------------------------------------------------
   compute particle's interaction energy with the rest of the system
------------------------------------------------------------------------- */

double FixWidom::energy(int i, int itype, tagint imolecule, double *coord)
{
  double delx,dely,delz,rsq;

  double **x = atom->x;
  int *type = atom->type;
  tagint *molecule = atom->molecule;
  int nall = atom->nlocal + atom->nghost;
  pair = force->pair;
  cutsq = force->pair->cutsq;

  double fpair = 0.0;
  double factor_coul = 1.0;
  double factor_lj = 1.0;

  double total_energy = 0.0;

  for (int j = 0; j < nall; j++) {

    if (i == j) continue;
    if (exchmode == EXCHMOL)
      if (imolecule == molecule[j]) continue;

    delx = coord[0] - x[j][0];
    dely = coord[1] - x[j][1];
    delz = coord[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    int jtype = type[j];

    if (rsq < cutsq[itype][jtype])
      total_energy +=
        pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
  }

  return total_energy;
}

/* ----------------------------------------------------------------------
   compute the energy of the given gas molecule in its current position
   sum across all procs that own atoms of the given molecule
------------------------------------------------------------------------- */

double FixWidom::molecule_energy(tagint gas_molecule_id)
{
  double mol_energy = 0.0;
  for (int i = 0; i < atom->nlocal; i++)
    if (atom->molecule[i] == gas_molecule_id) {
      mol_energy += energy(i,atom->type[i],gas_molecule_id,atom->x[i]);
    }

  double mol_energy_sum = 0.0;
  MPI_Allreduce(&mol_energy,&mol_energy_sum,1,MPI_DOUBLE,MPI_SUM,world);

  return mol_energy_sum;
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixWidom::energy_full()
{
  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build(1);
  int eflag = 1;
  int vflag = 0;

  // clear forces so they don't accumulate over multiple
  // calls within fix widom timestep

  size_t nbytes = sizeof(double) * (atom->nlocal + atom->nghost);
  if (nbytes) memset(&atom->f[0][0],0,3*nbytes);

  if (modify->n_pre_force) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) force->kspace->compute(eflag,vflag);

  if (modify->n_post_force_any) modify->post_force(vflag);

  // NOTE: all fixes with energy_global_flag set and which
  //   operate at pre_force() or post_force()
  //   and which user has enabled via fix_modify energy yes,
  //   will contribute to total MC energy via pe->compute_scalar()

  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}

/* ----------------------------------------------------------------------
  return acceptance ratios
------------------------------------------------------------------------- */

double FixWidom::compute_vector(int n)
{
  const double pot = ave_widom_chemical_potential;
  if (n == 0) return (pot > 0) ? -log(pot)/beta : 0.0;
  if (n == 1) return pot;
  if (n == 2) return volume;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixWidom::memory_usage()
{
  double bytes = (double)widom_nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixWidom::write_restart(FILE *fp)
{
  int n = 0;
  double list[3];
  list[n++] = random_equal->state();
  list[n++] = next_reneighbor;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixWidom::restart(char *buf)
{
  int n = 0;
  auto list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  random_equal->reset(seed);

  seed = static_cast<int> (list[n++]);

  next_reneighbor = static_cast<int> (list[n++]);
}

void FixWidom::grow_molecule_arrays(int nmolatoms) {
    nmaxmolatoms = nmolatoms;
    molcoords = memory->grow(molcoords,nmaxmolatoms,3,"widom:molcoords");
    molq = memory->grow(molq,nmaxmolatoms,"widom:molq");
    molimage = memory->grow(molimage,nmaxmolatoms,"widom:molimage");
}
