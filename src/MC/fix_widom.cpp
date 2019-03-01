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
   Contributing author: Paul Crozier, Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_widom.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_hybrid.h"
#include "molecule.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "compute.h"
#include "group.h"
#include "domain.h"
#include "region.h"
#include "random_park.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "thermo.h"
#include "output.h"
#include "neighbor.h"
#include <iostream>

using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define MAXENERGYTEST 1.0e50

enum{EXCHATOM,EXCHMOL}; // exchmode

/* ---------------------------------------------------------------------- */

FixWidom::FixWidom(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  idregion(NULL), full_flag(0),  
  local_gas_list(NULL), molcoords(NULL), molq(NULL), molimage(NULL),
  random_equal(NULL)
{
  if (narg < 8) error->all(FLERR,"Illegal fix Widom command");

  if (atom->molecular == 2)
    error->all(FLERR,"Fix Widom does not (yet) work with atom_style template");

  dynamic_group_allow = 1;

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  //ave_widom_chemical_potential = 0;

  // required args

  nevery = force->inumeric(FLERR,arg[3]);
  ninsertions = force->inumeric(FLERR,arg[4]);
  nwidom_type = force->inumeric(FLERR,arg[5]);
  seed = force->inumeric(FLERR,arg[6]);
  insertion_temperature = force->numeric(FLERR,arg[7]);

  if (nevery <= 0) error->all(FLERR,"Illegal fix Widom command");
  if (ninsertions < 0) error->all(FLERR,"Illegal fix Widom command");
  if (seed <= 0) error->all(FLERR,"Illegal fix Widom command");
  if (insertion_temperature < 0.0)
    error->all(FLERR,"Illegal fix Widom command");

  // read options from end of input line

  options(narg-8,&arg[8]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp,seed);

  // error checks on region and its extent being inside simulation box

  region_xlo = region_xhi = region_ylo = region_yhi =
    region_zlo = region_zhi = 0.0;
  if (regionflag) {
    if (domain->regions[iregion]->bboxflag == 0)
      error->all(FLERR,"Fix Widom region does not support a bounding box");
    if (domain->regions[iregion]->dynamic_check())
      error->all(FLERR,"Fix Widom region cannot be dynamic");

    region_xlo = domain->regions[iregion]->extent_xlo;
    region_xhi = domain->regions[iregion]->extent_xhi;
    region_ylo = domain->regions[iregion]->extent_ylo;
    region_yhi = domain->regions[iregion]->extent_yhi;
    region_zlo = domain->regions[iregion]->extent_zlo;
    region_zhi = domain->regions[iregion]->extent_zhi;

    if (region_xlo < domain->boxlo[0] || region_xhi > domain->boxhi[0] ||
        region_ylo < domain->boxlo[1] || region_yhi > domain->boxhi[1] ||
        region_zlo < domain->boxlo[2] || region_zhi > domain->boxhi[2])
      error->all(FLERR,"Fix Widom region extends outside simulation box");

    // estimate region volume using MC trials

    double coord[3];
    int inside = 0;
    int attempts = 10000000;
    for (int i = 0; i < attempts; i++) {
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      if (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) != 0)
        inside++;
    }

    double max_region_volume = (region_xhi - region_xlo)*
     (region_yhi - region_ylo)*(region_zhi - region_zlo);

    region_volume = max_region_volume*static_cast<double> (inside)/
     static_cast<double> (attempts);
  }

  // error check and further setup for exchmode = EXCHMOL

  if (exchmode == EXCHMOL) {
    if (onemols[imol]->xflag == 0)
      error->all(FLERR,"Fix Widom molecule must have coordinates");
    if (onemols[imol]->typeflag == 0)
      error->all(FLERR,"Fix Widom molecule must have atom types");
    if (nwidom_type != 0)
      error->all(FLERR,"Atom type must be zero in fix Widom mol command");
    if (onemols[imol]->qflag == 1 && atom->q == NULL)
      error->all(FLERR,"Fix Widom molecule has charges, but atom style does not");

    if (atom->molecular == 2 && onemols != atom->avec->onemols)
      error->all(FLERR,"Fix Widom molecule template ID must be same "
                 "as atom_style template ID");
    onemols[imol]->check_attributes(0);
  }

  if (charge_flag && atom->q == NULL)
    error->all(FLERR,"Fix Widom atom has charge, but atom style does not");

  // setup of array of coordinates for molecule insertion

  if (exchmode == EXCHATOM) natoms_per_molecule = 1;
  else natoms_per_molecule = onemols[imol]->natoms;
  nmaxmolatoms = natoms_per_molecule;
  grow_molecule_arrays(nmaxmolatoms);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters
  widom_nmax = 0;
  local_gas_list = NULL;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixWidom::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix Widom command");

  // defaults

  exchmode = EXCHATOM;
  regionflag = 0;
  iregion = -1;
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
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix Widom command");
      imol = atom->find_molecule(arg[iarg+1]);
      if (imol == -1)
        error->all(FLERR,"Molecule template ID for fix Widom does not exist");
      if (atom->molecules[imol]->nset > 1 && comm->me == 0)
        error->warning(FLERR,"Molecule template for "
                       "fix Widom has multiple molecules");
      exchmode = EXCHMOL;
      onemols = atom->molecules;
      nmol = onemols[imol]->nset;
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix Widom command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix Widom does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      regionflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"charge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix Widom command");
      charge = force->numeric(FLERR,arg[iarg+1]);
      charge_flag = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"full_energy") == 0) {
      full_flag = true;
      iarg += 1;
    } else if (strcmp(arg[iarg],"intra_energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix Widom command");
      energy_intra = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix Widom command");
  }
}

/* ---------------------------------------------------------------------- */

FixWidom::~FixWidom()
{
  if (regionflag) delete [] idregion;
  delete random_equal;

  memory->destroy(local_gas_list);
  memory->destroy(molcoords);
  memory->destroy(molq);
  memory->destroy(molimage);

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

  triclinic = domain->triclinic;

  ave_widom_chemical_potential = 0;

  // decide whether to switch to the full_energy option
  if (!full_flag) {
    if ((force->kspace) ||
        (force->pair == NULL) ||
        (force->pair->single_enable == 0) ||
        (force->pair_match("hybrid",0)) ||
        (force->pair_match("eam",0)) ||
        (force->pair->tail_flag)
        ) {
      full_flag = true;
      if (comm->me == 0)
        error->warning(FLERR,"Fix Widom using full_energy option");
    }
  }

  if (full_flag) {
    char *id_pe = (char *) "thermo_pe";
    int ipe = modify->find_compute(id_pe);
    c_pe = modify->compute[ipe];
  }

  int *type = atom->type;

  if (exchmode == EXCHATOM) {
    if (nwidom_type <= 0 || nwidom_type > atom->ntypes)
      error->all(FLERR,"Invalid atom type in fix Widom command");
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
      error->all(FLERR,
       "All mol IDs should be set for fix Widom group atoms");
  }

  if (exchmode == EXCHMOL)
    if (atom->molecule_flag == 0 || !atom->tag_enable || !atom->map_style)
      error->all(FLERR,
       "Fix Widom molecule command requires that "
       "atoms have molecule attributes");

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix Widom in a 2d simulation");

  // create a new group for interaction exclusions
  // used for attempted atom or molecule deletions
  // skip if already exists from previous init()

  if (full_flag && !exclusion_group_bit) {
    char **group_arg = new char*[4];

    // create unique group name for atoms to be excluded

    int len = strlen(id) + 30;
    group_arg[0] = new char[len];
    group_arg[1] = (char *) "subtract";
    group_arg[2] = (char *) "all";
    group_arg[3] = (char *) "all";
    group->assign(4,group_arg);
    exclusion_group = group->find(group_arg[0]);
    if (exclusion_group == -1)
      error->all(FLERR,"Could not find fix Widom exclusion group ID");
    exclusion_group_bit = group->bitmask[exclusion_group];

    // neighbor list exclusion setup
    // turn off interactions between group all and the exclusion group

    int narg = 4;
    char **arg = new char*[narg];;
    arg[0] = (char *) "exclude";
    arg[1] = (char *) "group";
    arg[2] = group_arg[0];
    arg[3] = (char *) "all";
    neighbor->modify_params(narg,arg);
    delete [] group_arg[0];
    delete [] group_arg;
    delete [] arg;
  }

  // create a new group for temporary use with selected molecules

  if (exchmode == EXCHMOL) {
    char **group_arg = new char*[3];
    // create unique group name for atoms to be rotated
    int len = strlen(id) + 30;
    group_arg[0] = new char[len];
    group_arg[1] = (char *) "molecule";
    char digits[12];
    sprintf(digits,"%d",-1);
    group_arg[2] = digits;
    group->assign(3,group_arg);
    molecule_group = group->find(group_arg[0]);
    if (molecule_group == -1)
      error->all(FLERR,"Could not find fix Widom rotation group ID");
    molecule_group_bit = group->bitmask[molecule_group];
    molecule_group_inversebit = molecule_group_bit ^ ~0;
    delete [] group_arg[0];
    delete [] group_arg;
  }

  // get all of the needed molecule data if exchanging
  // or moving molecules, otherwise just get the gas mass

  if (exchmode == EXCHMOL) {

    onemols[imol]->compute_mass();
    onemols[imol]->compute_com();
    gas_mass = onemols[imol]->masstotal;
    for (int i = 0; i < onemols[imol]->natoms; i++) {
      onemols[imol]->x[i][0] -= onemols[imol]->com[0];
      onemols[imol]->x[i][1] -= onemols[imol]->com[1];
      onemols[imol]->x[i][2] -= onemols[imol]->com[2];
    }
    onemols[imol]->com[0] = 0;
    onemols[imol]->com[1] = 0;
    onemols[imol]->com[2] = 0;

  } else gas_mass = atom->mass[nwidom_type];

  if (gas_mass <= 0.0)
    error->all(FLERR,"Illegal fix Widom gas mass <= 0");

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
      error->all(FLERR,"Cannot do Widom on atoms in atom_modify first group");
  }

  // compute beta
  beta = 1.0/(force->boltz*insertion_temperature);

  imagezero = ((imageint) IMGMAX << IMG2BITS) |
             ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  // Current implementation is broken using
  // full_flag on molecules on more than one processor.
  // Print error if this is the current mode
  if (full_flag && (exchmode == EXCHMOL) && comm->nprocs > 1)
    error->all(FLERR,"fix Widom does currently not support full_energy option with molecules on more than 1 MPI process.");

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

  if (regionflag) volume = region_volume;
  else volume = domain->xprd * domain->yprd * domain->zprd;

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  update_gas_atoms_list();

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

    if (regionflag) {
      int region_attempt = 0;
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      while (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) == 0) {
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
        error->one(FLERR,"Fix Widom put atom outside box");
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
      double inst_chem_pot = exp(-insertion_energy/beta);
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

    if (regionflag) {
      int region_attempt = 0;
      com_coord[0] = region_xlo + random_equal->uniform() *
        (region_xhi-region_xlo);
      com_coord[1] = region_ylo + random_equal->uniform() *
        (region_yhi-region_ylo);
      com_coord[2] = region_zlo + random_equal->uniform() *
        (region_zhi-region_zlo);
      while (domain->regions[iregion]->match(com_coord[0],com_coord[1],
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
    bool *procflag = new bool[natoms_per_molecule];

    for (int i = 0; i < natoms_per_molecule; i++) {
      MathExtra::matvec(rotmat,onemols[imol]->x[i],molcoords[i]);
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
        error->one(FLERR,"Fix Widom put atom outside box");

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
        if (onemols[imol]->qflag == 1) {
    ii = atom->nlocal + atom->nghost;
    if (ii >= atom->nmax) atom->avec->grow(0);
    atom->q[ii] = onemols[imol]->q[i];
        }
        insertion_energy += energy(ii,onemols[imol]->type[i],-1,xtmp);
      }
    }

    double insertion_energy_sum = 0.0;
    MPI_Allreduce(&insertion_energy,&insertion_energy_sum,1,
      MPI_DOUBLE,MPI_SUM,world);

    // the insertion_energy_sum is the variable with the energy of inserting one molecule
    double inst_chem_pot = exp(-insertion_energy_sum/beta);
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

    if (regionflag) {
      int region_attempt = 0;
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      while (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) == 0) {
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
        error->one(FLERR,"Fix Widom put atom outside box");
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
      if (atom->map_style) atom->map_init();
    }
    atom->nghost = 0;
    if (triclinic) domain->x2lamda(atom->nlocal);
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();

    double insertion_energy = energy_full() - energy_stored;
    double inst_chem_pot = exp(-insertion_energy/beta);
    double incr_chem_pot = (inst_chem_pot - ave_widom_chemical_potential);
    ave_widom_chemical_potential += incr_chem_pot / (imove + 1);

    atom->natoms--;
    if (proc_flag) atom->nlocal--;
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
    
    update_gas_atoms_list();

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
    error->all(FLERR,"Fix Widom ran out of available molecule IDs");
  int insertion_molecule = maxmol_all;

  tagint maxtag = 0;
  for (int i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag,atom->tag[i]);
  tagint maxtag_all;
  MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

  for (int imove = 0; imove < ninsertions; imove++) {

    int nlocalprev = atom->nlocal;

    double com_coord[3];
    if (regionflag) {
      int region_attempt = 0;
      com_coord[0] = region_xlo + random_equal->uniform() *
        (region_xhi-region_xlo);
      com_coord[1] = region_ylo + random_equal->uniform() *
        (region_yhi-region_ylo);
      com_coord[2] = region_zlo + random_equal->uniform() *
        (region_zhi-region_zlo);
      while (domain->regions[iregion]->match(com_coord[0],com_coord[1],
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
      MathExtra::matvec(rotmat,onemols[imol]->x[i],xtmp);
      xtmp[0] += com_coord[0];
      xtmp[1] += com_coord[1];
      xtmp[2] += com_coord[2];

      // need to adjust image flags in remap()

      imageint imagetmp = imagezero;
      domain->remap(xtmp,imagetmp);
      if (!domain->inside(xtmp))
        error->one(FLERR,"Fix Widom put atom outside box");

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
        atom->avec->create_atom(onemols[imol]->type[i],xtmp);
        int m = atom->nlocal - 1;

        atom->image[m] = imagetmp;
        atom->molecule[m] = insertion_molecule;
        if (maxtag_all+i+1 >= MAXTAGINT)
    error->all(FLERR,"Fix Widom ran out of available atom IDs");
        atom->tag[m] = maxtag_all + i + 1;
        atom->v[m][0] = 0;
        atom->v[m][1] = 0;
        atom->v[m][2] = 0;

        atom->add_molecule_atom(onemols[imol],i,m,maxtag_all);
        modify->create_attribute(m);
      }
    }

    atom->natoms += natoms_per_molecule;
    if (atom->natoms < 0)
      error->all(FLERR,"Too many total atoms");
    atom->nbonds += onemols[imol]->nbonds;
    atom->nangles += onemols[imol]->nangles;
    atom->ndihedrals += onemols[imol]->ndihedrals;
    atom->nimpropers += onemols[imol]->nimpropers;
    if (atom->map_style) atom->map_init();
    atom->nghost = 0;
    if (triclinic) domain->x2lamda(atom->nlocal);
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();

    // energy_after corrected by energy_intra
    double insertion_energy = (energy_full() -energy_intra) - energy_stored;
    double inst_chem_pot = exp(-insertion_energy/beta);
    double incr_chem_pot = (inst_chem_pot - ave_widom_chemical_potential);
    ave_widom_chemical_potential += incr_chem_pot / (imove + 1);

    atom->nbonds -= onemols[imol]->nbonds;
    atom->nangles -= onemols[imol]->nangles;
    atom->ndihedrals -= onemols[imol]->ndihedrals;
    atom->nimpropers -= onemols[imol]->nimpropers;
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
    
    update_gas_atoms_list();

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
  int imolecule;

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

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) force->kspace->compute(eflag,vflag);

  // unlike Verlet, not performing a reverse_comm() or forces here
  // b/c Widom does not care about forces
  // don't think it will mess up energy due to any post_force() fixes

  if (modify->n_post_force) modify->post_force(vflag);
  if (modify->n_end_of_step) modify->end_of_step();

  // NOTE: all fixes with THERMO_ENERGY mask set and which
  //   operate at pre_force() or post_force() or end_of_step()
  //   and which user has enable via fix_modify thermo yes,
  //   will contribute to total MC energy via pe->compute_scalar()

  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}

/* ----------------------------------------------------------------------
   update the list of gas atoms
------------------------------------------------------------------------- */

void FixWidom::update_gas_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

  if (atom->nmax > widom_nmax) {
    memory->sfree(local_gas_list);
    widom_nmax = atom->nmax;
    local_gas_list = (int *) memory->smalloc(widom_nmax*sizeof(int),
     "Widom:local_gas_list");
  }

  ngas_local = 0;

  if (regionflag) {

    if (exchmode == EXCHMOL) {

      tagint maxmol = 0;
      for (int i = 0; i < nlocal; i++) maxmol = MAX(maxmol,molecule[i]);
      tagint maxmol_all;
      MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
      double *comx = new double[maxmol_all];
      double *comy = new double[maxmol_all];
      double *comz = new double[maxmol_all];
      for (int imolecule = 0; imolecule < maxmol_all; imolecule++) {
        for (int i = 0; i < nlocal; i++) {
          if (molecule[i] == imolecule) {
            mask[i] |= molecule_group_bit;
          } else {
            mask[i] &= molecule_group_inversebit;
          }
        }
        double com[3];
        com[0] = com[1] = com[2] = 0.0;
        group->xcm(molecule_group,gas_mass,com);

        // remap unwrapped com into periodic box

        domain->remap(com);
        comx[imolecule] = com[0];
        comy[imolecule] = com[1];
        comz[imolecule] = com[2];
      }

      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          if (domain->regions[iregion]->match(comx[molecule[i]],
             comy[molecule[i]],comz[molecule[i]]) == 1) {
            local_gas_list[ngas_local] = i;
            ngas_local++;
          }
        }
      }
      delete[] comx;
      delete[] comy;
      delete[] comz;
    } else {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]) == 1) {
            local_gas_list[ngas_local] = i;
            ngas_local++;
          }
        }
      }
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        local_gas_list[ngas_local] = i;
        ngas_local++;
      }
    }
  }

  MPI_Allreduce(&ngas_local,&ngas,1,MPI_INT,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
  return acceptance ratios
------------------------------------------------------------------------- */

double FixWidom::compute_vector(int n)
{
  
  if (n == 0) return -beta*log(ave_widom_chemical_potential);
  if (n == 1) return ave_widom_chemical_potential;
  if (n == 2) return volume;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixWidom::memory_usage()
{
  double bytes = widom_nmax * sizeof(int);
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
  double *list = (double *) buf;

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
