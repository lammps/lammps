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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_gcmc.h"
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

enum{ATOM,MOLECULE};

/* ---------------------------------------------------------------------- */

FixGCMC::FixGCMC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 11) error->all(FLERR,"Illegal fix gcmc command");

  if (atom->molecular == 2) 
    error->all(FLERR,"Fix gcmc does not (yet) work with atom_style template");

  dynamic_group_allow = 1;
    
  vector_flag = 1;
  size_vector = 8;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // required args

  nevery = force->inumeric(FLERR,arg[3]);
  nexchanges = force->inumeric(FLERR,arg[4]);
  nmcmoves = force->inumeric(FLERR,arg[5]);
  ngcmc_type = force->inumeric(FLERR,arg[6]);
  seed = force->inumeric(FLERR,arg[7]);
  reservoir_temperature = force->numeric(FLERR,arg[8]);
  chemical_potential = force->numeric(FLERR,arg[9]);
  displace = force->numeric(FLERR,arg[10]);

  if (nexchanges < 0) error->all(FLERR,"Illegal fix gcmc command");
  if (nmcmoves < 0) error->all(FLERR,"Illegal fix gcmc command");
  if (seed <= 0) error->all(FLERR,"Illegal fix gcmc command");
  if (reservoir_temperature < 0.0)
    error->all(FLERR,"Illegal fix gcmc command");
  if (displace < 0.0) error->all(FLERR,"Illegal fix gcmc command");

  // read options from end of input line

  options(narg-11,&arg[11]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp,seed);

  // random number generator, not the same for all procs

  random_unequal = new RanPark(lmp,seed);
  
  // error checks on region and its extent being inside simulation box

  region_xlo = region_xhi = region_ylo = region_yhi = 
    region_zlo = region_zhi = 0.0;
  if (regionflag) {
    if (domain->regions[iregion]->bboxflag == 0)
      error->all(FLERR,"Fix gcmc region does not support a bounding box");
    if (domain->regions[iregion]->dynamic_check())
      error->all(FLERR,"Fix gcmc region cannot be dynamic");
    
    region_xlo = domain->regions[iregion]->extent_xlo;
    region_xhi = domain->regions[iregion]->extent_xhi;
    region_ylo = domain->regions[iregion]->extent_ylo;
    region_yhi = domain->regions[iregion]->extent_yhi;
    region_zlo = domain->regions[iregion]->extent_zlo;
    region_zhi = domain->regions[iregion]->extent_zhi;

    if (region_xlo < domain->boxlo[0] || region_xhi > domain->boxhi[0] ||
        region_ylo < domain->boxlo[1] || region_yhi > domain->boxhi[1] ||
        region_zlo < domain->boxlo[2] || region_zhi > domain->boxhi[2])
      error->all(FLERR,"Fix gcmc region extends outside simulation box");

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

  // error check and further setup for mode = MOLECULE

  if (mode == MOLECULE) {
    if (onemols[imol]->xflag == 0)
      error->all(FLERR,"Fix gcmc molecule must have coordinates");
    if (onemols[imol]->typeflag == 0)
      error->all(FLERR,"Fix gcmc molecule must have atom types");
    if (ngcmc_type+onemols[imol]->ntypes <= 0 || ngcmc_type+onemols[imol]->ntypes > atom->ntypes)
      error->all(FLERR,"Invalid atom type in fix gcmc mol command");

    if (atom->molecular == 2 && onemols != atom->avec->onemols)
      error->all(FLERR,"Fix gcmc molecule template ID must be same "
                 "as atom_style template ID");
    onemols[imol]->check_attributes(0);
  }

  if (shakeflag && mode == ATOM)
    error->all(FLERR,"Cannot use fix gcmc shake and not molecule");

  // setup of coords and imageflags array

  if (mode == ATOM) natoms_per_molecule = 1;
  else natoms_per_molecule = onemols[imol]->natoms;
  memory->create(coords,natoms_per_molecule,3,"gcmc:coords");
  memory->create(imageflags,natoms_per_molecule,"gcmc:imageflags");
  memory->create(atom_coord,natoms_per_molecule,3,"gcmc:atom_coord");

  // compute the number of MC cycles that occur nevery timesteps

  ncycles = nexchanges + nmcmoves;

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters

  ntranslation_attempts = 0.0;
  ntranslation_successes = 0.0;
  nrotation_attempts = 0.0;
  nrotation_successes = 0.0;
  ndeletion_attempts = 0.0;
  ndeletion_successes = 0.0;
  ninsertion_attempts = 0.0;
  ninsertion_successes = 0.0;

  gcmc_nmax = 0;
  local_gas_list = NULL;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixGCMC::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix gcmc command");

  // defaults

  mode = ATOM;
  max_rotation_angle = 10*MY_PI/180;
  regionflag = 0; 
  iregion = -1; 
  region_volume = 0;
  max_region_attempts = 1000; 
  molecule_group = 0;
  molecule_group_bit = 0;
  molecule_group_inversebit = 0;
  exclusion_group = 0;
  exclusion_group_bit = 0;
  exclusion_group_inversebit = 0;
  pressure_flag = false;
  pressure = 0.0;
  fugacity_coeff = 1.0;
  shakeflag = 0;
  charge = 0.0;
  charge_flag = false;
  full_flag = false;
  idshake = NULL;
  
  int iarg = 0;
  while (iarg < narg) {
  if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      imol = atom->find_molecule(arg[iarg+1]);
      if (imol == -1)
        error->all(FLERR,"Molecule template ID for fix gcmc does not exist");
      if (atom->molecules[imol]->nset > 1 && comm->me == 0)
        error->warning(FLERR,"Molecule template for "
                       "fix gcmc has multiple molecules");
      mode = MOLECULE;
      onemols = &atom->molecules[imol];
      nmol = onemols[0]->nset;
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix gcmc does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      regionflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"maxangle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      max_rotation_angle = force->numeric(FLERR,arg[iarg+1]);
      max_rotation_angle *= MY_PI/180;
      iarg += 2;
    } else if (strcmp(arg[iarg],"pressure") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      pressure = force->numeric(FLERR,arg[iarg+1]);
      pressure_flag = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"fugacity_coeff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      fugacity_coeff = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"charge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      charge = force->numeric(FLERR,arg[iarg+1]);
      charge_flag = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"shake") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      int n = strlen(arg[iarg+1]) + 1;
      delete [] idshake;
      idshake = new char[n];
      strcpy(idshake,arg[iarg+1]);
      shakeflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"full_energy") == 0) {
      full_flag = true;
      iarg += 1;
    } else error->all(FLERR,"Illegal fix gcmc command");
  }
}

/* ---------------------------------------------------------------------- */

FixGCMC::~FixGCMC()
{
  if (regionflag) delete [] idregion;
  delete random_equal;
  delete random_unequal;

  // remove groups this fix defined
  // only do it if the groups exists and group itself exists

  if (molecule_group && (strcmp(group->names[0],"all") == 0)) {
    char **group_arg = new char*[2];
    group_arg[0] = group->names[molecule_group];
    group_arg[1] = (char *) "delete";
    group->assign(2,group_arg);
    delete [] group_arg;
  } 
  
  if (exclusion_group && (strcmp(group->names[0],"all") == 0)) {
    char **group_arg = new char*[2];
    group_arg[0] = group->names[exclusion_group];
    group_arg[1] = (char *) "delete";
    group->assign(2,group_arg);
    delete [] group_arg;
  } 

  memory->destroy(local_gas_list);
  memory->destroy(atom_coord);
  memory->destroy(coords);
  memory->destroy(imageflags);
  
  delete [] idshake;
}

/* ---------------------------------------------------------------------- */

int FixGCMC::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGCMC::init()
{ 
  // decide whether to switch to the full_energy option

  if (!full_flag) {
    if ((force->kspace) || 
        (force->pair == NULL) ||
        (force->pair->single_enable == 0) ||
        (force->pair_match("hybrid",0)) ||
        (force->pair_match("eam",0)) ||
        (domain->triclinic == 1)) {
      full_flag = true;
      if (comm->me == 0) 
        error->warning(FLERR,"fix gcmc using full_energy option");
    }
  }
  
  if (full_flag) {
    char *id_pe = (char *) "thermo_pe";
    int ipe = modify->find_compute(id_pe);
    c_pe = modify->compute[ipe];
  }

  int *type = atom->type;

  if (mode == ATOM) {
    if (ngcmc_type <= 0 || ngcmc_type > atom->ntypes)
      error->all(FLERR,"Invalid atom type in fix gcmc command");
  }

  // if mode == ATOM, warn if any deletable atom has a mol ID

  if ((mode == ATOM) && atom->molecule_flag) {
    tagint *molecule = atom->molecule;
    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if (type[i] == ngcmc_type)
        if (molecule[i]) flag = 1;
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall && comm->me == 0)
      error->all(FLERR,
       "Fix gcmc cannot exchange individual atoms belonging to a molecule");
  }

  // if mode == MOLECULE, check for unset mol IDs

  if (mode == MOLECULE) {
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
       "All mol IDs should be set for fix gcmc group atoms");
  }

  if (((mode == MOLECULE) && (atom->molecule_flag == 0)) || 
      ((mode == MOLECULE) && (!atom->tag_enable || !atom->map_style)))
    error->all(FLERR,
               "Fix gcmc molecule command requires that "
               "atoms have molecule attributes");

  // if shakeflag defined, check for SHAKE fix
  // its molecule template must be same as this one

  fixshake = NULL;
  if (shakeflag) {
    int ifix = modify->find_fix(idshake);
    if (ifix < 0) error->all(FLERR,"Fix gcmc shake fix does not exist");
    fixshake = modify->fix[ifix];
    int tmp;
    if (onemols != (Molecule **) fixshake->extract("onemol",tmp))
      error->all(FLERR,"Fix gcmc and fix shake not using "
                 "same molecule template ID");
  }

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix gcmc in a 2d simulation");

  // create a new group for interaction exclusions

  if (full_flag) {
    char **group_arg = new char*[4];
    // create unique group name for atoms to be excluded
    int len = strlen(id) + 30;
    group_arg[0] = new char[len];
    sprintf(group_arg[0],"FixGCMC:gcmc_exclusion_group:%s",id);
    group_arg[1] = (char *) "subtract";
    group_arg[2] = (char *) "all";
    group_arg[3] = (char *) "all";
    group->assign(4,group_arg);
    exclusion_group = group->find(group_arg[0]);
    if (exclusion_group == -1) 
      error->all(FLERR,"Could not find fix gcmc exclusion group ID");
    exclusion_group_bit = group->bitmask[exclusion_group];
    exclusion_group_inversebit = exclusion_group_bit ^ ~0;
    
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
  }
    
  // create a new group for temporary use with selected molecules

  if (mode == MOLECULE) {
    char **group_arg = new char*[3];
    // create unique group name for atoms to be rotated
    int len = strlen(id) + 30;
    group_arg[0] = new char[len];
    sprintf(group_arg[0],"FixGCMC:rotation_gas_atoms:%s",id);
    group_arg[1] = (char *) "molecule";
    char digits[12];
    sprintf(digits,"%d",ngcmc_type);
    group_arg[2] = digits;
    group->assign(3,group_arg);
    molecule_group = group->find(group_arg[0]);
    if (molecule_group == -1) 
      error->all(FLERR,"Could not find fix gcmc rotation group ID");
    molecule_group_bit = group->bitmask[molecule_group];
    molecule_group_inversebit = molecule_group_bit ^ ~0;
    delete [] group_arg[0];
  }
    
  // get all of the needed molecule data if mode == MOLECULE, 
  // otherwise just get the gas mass
  
  if (mode == MOLECULE) {
    onemols[imol]->compute_mass();
    onemols[imol]->compute_com();
    gas_mass = onemols[imol]->masstotal;

    for (int i = 0; i < onemols[imol]->natoms; i++) {
      onemols[imol]->x[i][0] -= onemols[imol]->com[0];
      onemols[imol]->x[i][1] -= onemols[imol]->com[1];
      onemols[imol]->x[i][2] -= onemols[imol]->com[2];
    }
    
  } else gas_mass = atom->mass[ngcmc_type];
  
  if (gas_mass <= 0.0)
    error->all(FLERR,"Illegal fix gcmc gas mass <= 0");
  
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
      error->all(FLERR,"Cannot do GCMC on atoms in atom_modify first group");
  }
  
  // compute beta, lambda, sigma, and the zz factor

  beta = 1.0/(force->boltz*reservoir_temperature);
  double lambda = sqrt(force->hplanck*force->hplanck/
                       (2.0*MY_PI*gas_mass*force->mvv2e*
                        force->boltz*reservoir_temperature));
  sigma = sqrt(force->boltz*reservoir_temperature/gas_mass/force->mvv2e);
  zz = exp(beta*chemical_potential)/(pow(lambda,3.0));
  if (pressure_flag) zz = pressure*fugacity_coeff*beta/force->nktv2p;
  
  imagetmp = ((imageint) IMGMAX << IMG2BITS) | 
             ((imageint) IMGMAX << IMGBITS) | IMGMAX;
}

/* ----------------------------------------------------------------------
   attempt Monte Carlo translations, rotations, insertions, and deletions
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixGCMC::pre_exchange()
{
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];
  sublo = domain->sublo;
  subhi = domain->subhi;

  if (regionflag) volume = region_volume;
  else volume = domain->xprd * domain->yprd * domain->zprd;

  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  update_gas_atoms_list();

  if (full_flag) {
    energy_stored = energy_full();

    if (mode == MOLECULE) {
      for (int i = 0; i < ncycles; i++) {
        int random_int_fraction =
          static_cast<int>(random_equal->uniform()*ncycles) + 1;
        if (random_int_fraction <= nmcmoves) {
          if (random_equal->uniform() < 0.5) attempt_molecule_translation_full();
          else attempt_molecule_rotation_full();
        } else {
          if (random_equal->uniform() < 0.5) attempt_molecule_deletion_full();
          else attempt_molecule_insertion_full();
        }
      }
    } else {
      for (int i = 0; i < ncycles; i++) {
        int random_int_fraction =
          static_cast<int>(random_equal->uniform()*ncycles) + 1;
        if (random_int_fraction <= nmcmoves) {
          attempt_atomic_translation_full();
        } else {
          if (random_equal->uniform() < 0.5) attempt_atomic_deletion_full();
          else attempt_atomic_insertion_full();
        }
      }
    }
    domain->pbc();
    comm->exchange();
    atom->nghost = 0;
    comm->borders();
    
  } else {
    
    if (mode == MOLECULE) {
      for (int i = 0; i < ncycles; i++) {
        int random_int_fraction =
          static_cast<int>(random_equal->uniform()*ncycles) + 1;
        if (random_int_fraction <= nmcmoves) {
          if (random_equal->uniform() < 0.5) attempt_molecule_translation();
          else attempt_molecule_rotation();
        } else {
          if (random_equal->uniform() < 0.5) attempt_molecule_deletion();
          else attempt_molecule_insertion();
        }
      }
    } else {
      for (int i = 0; i < ncycles; i++) {
        int random_int_fraction =
          static_cast<int>(random_equal->uniform()*ncycles) + 1;
        if (random_int_fraction <= nmcmoves) {
          attempt_atomic_translation();
        } else {
          if (random_equal->uniform() < 0.5) attempt_atomic_deletion();
          else attempt_atomic_insertion();
        }
      }
    }
  }
  next_reneighbor = update->ntimestep + nevery;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_atomic_translation()
{
  ntranslation_attempts += 1.0;
  
  if (ngas == 0) return;

  int i = pick_random_gas_atom();
  
  int success = 0;
  if (i >= 0) {
    double **x = atom->x;
    double energy_before = energy(i,ngcmc_type,-1,x[i]);
    double rsq = 1.1;
    double rx,ry,rz;
    rx = ry = rz = 0.0;
    double coord[3];
    while (rsq > 1.0) {
      rx = 2*random_unequal->uniform() - 1.0;
      ry = 2*random_unequal->uniform() - 1.0;
      rz = 2*random_unequal->uniform() - 1.0;
      rsq = rx*rx + ry*ry + rz*rz;
    }    
    coord[0] = x[i][0] + displace*rx;
    coord[1] = x[i][1] + displace*ry;
    coord[2] = x[i][2] + displace*rz;
    if (regionflag) {
      while (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) == 0) {
        while (rsq > 1.0) {
          rx = 2*random_unequal->uniform() - 1.0;
          ry = 2*random_unequal->uniform() - 1.0;
          rz = 2*random_unequal->uniform() - 1.0;
          rsq = rx*rx + ry*ry + rz*rz;
        }    
        coord[0] = x[i][0] + displace*rx;
        coord[1] = x[i][1] + displace*ry;
        coord[2] = x[i][2] + displace*rz;
      }
    }

    double energy_after = energy(i,ngcmc_type,-1,coord);
    if (random_unequal->uniform() < 
        exp(beta*(energy_before - energy_after))) {
      x[i][0] = coord[0];
      x[i][1] = coord[1];
      x[i][2] = coord[2];
      success = 1;
    }
  }

  int success_all = 0;
  MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);

  if (success_all) {
    domain->pbc();
    comm->exchange();
    atom->nghost = 0;
    comm->borders();
    update_gas_atoms_list();
    ntranslation_successes += 1.0;
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_atomic_deletion()
{
  ndeletion_attempts += 1.0;

  if (ngas == 0) return;
  
  int i = pick_random_gas_atom();

  int success = 0;
  if (i >= 0) {
    double deletion_energy = energy(i,ngcmc_type,-1,atom->x[i]);
    if (random_unequal->uniform() < 
        ngas*exp(beta*deletion_energy)/(zz*volume)) {
      atom->avec->copy(atom->nlocal-1,i,1);
      atom->nlocal--;
      success = 1;
    }
  }

  int success_all = 0;
  MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);

  if (success_all) {
    atom->natoms--;
    if (atom->tag_enable) {
      if (atom->map_style) atom->map_init();
    }
    atom->nghost = 0;
    comm->borders();
    update_gas_atoms_list();
    ndeletion_successes += 1.0;
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_atomic_insertion()
{
  ninsertion_attempts += 1.0;

  double coord[3];
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
  } else {
    coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
    coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
    coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
  }

  int proc_flag = 0;
  if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
      coord[1] >= sublo[1] && coord[1] < subhi[1] &&
      coord[2] >= sublo[2] && coord[2] < subhi[2]) proc_flag = 1;

  int success = 0;
  if (proc_flag) {
    double insertion_energy = energy(-1,ngcmc_type,-1,coord);
    if (random_unequal->uniform() <
        zz*volume*exp(-beta*insertion_energy)/(ngas+1)) {
      atom->avec->create_atom(ngcmc_type,coord);
      int m = atom->nlocal - 1;
      atom->mask[m] = 1 | groupbit;
      atom->v[m][0] = random_unequal->gaussian()*sigma;
      atom->v[m][1] = random_unequal->gaussian()*sigma;
      atom->v[m][2] = random_unequal->gaussian()*sigma;
      modify->create_attribute(m);

      success = 1;
    }
  }

  int success_all = 0;
  MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);

  if (success_all) {
    atom->natoms++;
    if (atom->tag_enable) {
      atom->tag_extend();
      if (atom->map_style) atom->map_init();
    }
    atom->nghost = 0;
    comm->borders();
    update_gas_atoms_list();
    ninsertion_successes += 1.0;
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_molecule_translation()
{
  ntranslation_attempts += 1.0;

  if (ngas == 0) return;

  tagint translation_molecule = pick_random_gas_molecule();
  if (translation_molecule == -1) return;

  double energy_before_sum = molecule_energy(translation_molecule);
  
  double **x = atom->x;
  double rx,ry,rz;
  double com_displace[3],coord[3];
  double rsq = 1.1;
  while (rsq > 1.0) {
    rx = 2*random_equal->uniform() - 1.0;
    ry = 2*random_equal->uniform() - 1.0;
    rz = 2*random_equal->uniform() - 1.0;
    rsq = rx*rx + ry*ry + rz*rz;
  }
  com_displace[0] = displace*rx;
  com_displace[1] = displace*ry;
  com_displace[2] = displace*rz;
  
  int nlocal = atom->nlocal;
  if (regionflag) {
    int *mask = atom->mask;
    for (int i = 0; i < nlocal; i++) {
      if (atom->molecule[i] == translation_molecule) {
        mask[i] |= molecule_group_bit;
      } else {
        mask[i] &= molecule_group_inversebit;
      }
    }
    double com[3];
    com[0] = com[1] = com[2] = 0.0;
    group->xcm(molecule_group,gas_mass,com);
    coord[0] = com[0] + displace*rx;
    coord[1] = com[1] + displace*ry;
    coord[2] = com[2] + displace*rz;
    while (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) == 0) {
      while (rsq > 1.0) {
        rx = 2*random_equal->uniform() - 1.0;
        ry = 2*random_equal->uniform() - 1.0;
        rz = 2*random_equal->uniform() - 1.0;
        rsq = rx*rx + ry*ry + rz*rz;
      }    
      coord[0] = com[0] + displace*rx;
      coord[1] = com[1] + displace*ry;
      coord[2] = com[2] + displace*rz;
    }
    com_displace[0] = displace*rx;
    com_displace[1] = displace*ry;
    com_displace[2] = displace*rz;
  }

  double energy_after = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (atom->molecule[i] == translation_molecule) {
      coord[0] = x[i][0] + com_displace[0];
      coord[1] = x[i][1] + com_displace[1];
      coord[2] = x[i][2] + com_displace[2];
      energy_after += energy(i,atom->type[i],translation_molecule,coord);
    }
  }

  double energy_after_sum = 0.0;
  MPI_Allreduce(&energy_after,&energy_after_sum,1,MPI_DOUBLE,MPI_SUM,world);

  if (random_equal->uniform() < 
      exp(beta*(energy_before_sum - energy_after_sum))) {
    for (int i = 0; i < nlocal; i++) {
      if (atom->molecule[i] == translation_molecule) {
        x[i][0] += com_displace[0];
        x[i][1] += com_displace[1];
        x[i][2] += com_displace[2];
      }
    }
    domain->pbc();
    comm->exchange();
    atom->nghost = 0;
    comm->borders();
    update_gas_atoms_list();
    ntranslation_successes += 1.0;
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_molecule_rotation()
{
  nrotation_attempts += 1.0;

  if (ngas == 0) return;

  tagint rotation_molecule = pick_random_gas_molecule();
  if (rotation_molecule == -1) return;
  

  double energy_before_sum = molecule_energy(rotation_molecule);

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  for (int i = 0; i < nlocal; i++) {
    if (atom->molecule[i] == rotation_molecule) {
      mask[i] |= molecule_group_bit;
    } else {
      mask[i] &= molecule_group_inversebit;
    }
  }

  double com[3];
  com[0] = com[1] = com[2] = 0.0;
  group->xcm(molecule_group,gas_mass,com);

  double r[3],rotmat[3][3],quat[4];
  r[0] = random_equal->uniform() - 0.5;
  r[1] = random_equal->uniform() - 0.5;
  r[2] = random_equal->uniform() - 0.5;

  double theta = random_equal->uniform() * max_rotation_angle;
  MathExtra::norm3(r);
  MathExtra::axisangle_to_quat(r,theta,quat);
  MathExtra::quat_to_mat(quat,rotmat);

  double **x = atom->x;
  imageint *image = atom->image;
  double energy_after = 0.0;
  int n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & molecule_group_bit) {
      double xtmp[3];
      domain->unmap(x[i],image[i],xtmp);
      xtmp[0] -= com[0];
      xtmp[1] -= com[1];
      xtmp[2] -= com[2];
      MathExtra::matvec(rotmat,xtmp,atom_coord[n]);    
      atom_coord[n][0] += com[0];
      atom_coord[n][1] += com[1];
      atom_coord[n][2] += com[2];
      xtmp[0] = atom_coord[n][0];
      xtmp[1] = atom_coord[n][1];
      xtmp[2] = atom_coord[n][2];
      domain->remap(xtmp);
      energy_after += energy(i,atom->type[i],rotation_molecule,xtmp);
      n++;
    }
  }

  double energy_after_sum = 0.0;
  MPI_Allreduce(&energy_after,&energy_after_sum,1,MPI_DOUBLE,MPI_SUM,world);

  if (random_equal->uniform() < 
      exp(beta*(energy_before_sum - energy_after_sum))) {
    int n = 0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & molecule_group_bit) {
        image[i] = imagetmp;
        x[i][0] = atom_coord[n][0];
        x[i][1] = atom_coord[n][1];
        x[i][2] = atom_coord[n][2];
        domain->remap(x[i],image[i]);
        n++;
      }
    }
    domain->pbc();
    comm->exchange();
    atom->nghost = 0;
    comm->borders();
    update_gas_atoms_list();
    nrotation_successes += 1.0;
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_molecule_deletion()
{
  ndeletion_attempts += 1.0;

  if (ngas == 0) return;
  
  tagint deletion_molecule = pick_random_gas_molecule();
  if (deletion_molecule == -1) return;

  double deletion_energy_sum = molecule_energy(deletion_molecule);

  if (random_equal->uniform() < 
      ngas*exp(beta*deletion_energy_sum)/(zz*volume*natoms_per_molecule)) {
    int i = 0;
    while (i < atom->nlocal) {
      if (atom->molecule[i] == deletion_molecule) {
        atom->avec->copy(atom->nlocal-1,i,1);
        atom->nlocal--;
      } else i++;
    }
    atom->natoms -= natoms_per_molecule;
    if (atom->map_style) atom->map_init();
    atom->nghost = 0;
    comm->borders();
    update_gas_atoms_list();
    ndeletion_successes += 1.0;
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_molecule_insertion()
{
  ninsertion_attempts += 1.0;

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
  } else {
    com_coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
    com_coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
    com_coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
  }

  double r[3],rotmat[3][3],quat[4];
  r[0] = random_equal->uniform() - 0.5;
  r[1] = random_equal->uniform() - 0.5;
  r[2] = random_equal->uniform() - 0.5;

  double theta = random_equal->uniform() * MY_2PI;
  MathExtra::norm3(r);
  MathExtra::axisangle_to_quat(r,theta,quat);
  MathExtra::quat_to_mat(quat,rotmat);
    
  double insertion_energy = 0.0;
  bool procflag[natoms_per_molecule];

  for (int i = 0; i < natoms_per_molecule; i++) {
    MathExtra::matvec(rotmat,onemols[imol]->x[i],atom_coord[i]);
    atom_coord[i][0] += com_coord[0];
    atom_coord[i][1] += com_coord[1];
    atom_coord[i][2] += com_coord[2];

    double xtmp[3];
    xtmp[0] = atom_coord[i][0];
    xtmp[1] = atom_coord[i][1];
    xtmp[2] = atom_coord[i][2];
    domain->remap(xtmp);

    procflag[i] = false;
    if (xtmp[0] >= sublo[0] && xtmp[0] < subhi[0] &&
        xtmp[1] >= sublo[1] && xtmp[1] < subhi[1] &&
        xtmp[2] >= sublo[2] && xtmp[2] < subhi[2]) {
      procflag[i] = true;
      insertion_energy += energy(-1,onemols[imol]->type[i],-1,xtmp);
    }
  }

  double insertion_energy_sum = 0.0;
  MPI_Allreduce(&insertion_energy,&insertion_energy_sum,1,
                MPI_DOUBLE,MPI_SUM,world);

  if (random_equal->uniform() < zz*volume*natoms_per_molecule*
      exp(-beta*insertion_energy_sum)/(ngas + natoms_per_molecule)) {  

    tagint maxmol = 0;
    for (int i = 0; i < atom->nlocal; i++) maxmol = MAX(maxmol,atom->molecule[i]);
    tagint maxmol_all;
    MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
    maxmol_all++;
    if (maxmol_all >= MAXTAGINT) 
      error->all(FLERR,"Fix gcmc ran out of available molecule IDs");

    tagint maxtag = 0;
    for (int i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag,atom->tag[i]);
    tagint maxtag_all;
    MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
    
    int nfix = modify->nfix;
    Fix **fix = modify->fix;

    int nlocalprev = atom->nlocal;
    
    double vnew[3];
    vnew[0] = random_unequal->gaussian()*sigma;
    vnew[1] = random_unequal->gaussian()*sigma;
    vnew[2] = random_unequal->gaussian()*sigma;
    
    for (int i = 0; i < natoms_per_molecule; i++) {
      if (procflag[i]) {
        atom->avec->create_atom(ngcmc_type+onemols[imol]->type[i],atom_coord[i]);
        int m = atom->nlocal - 1;
        atom->mask[m] = 1 | groupbit;
        atom->image[m] = imagetmp;
        domain->remap(atom->x[m],atom->image[m]);
        atom->molecule[m] = maxmol_all;
        if (maxtag_all+i+1 >= MAXTAGINT)
          error->all(FLERR,"Fix gcmc ran out of available atom IDs");
        atom->tag[m] = maxtag_all + i + 1;
        atom->v[m][0] = vnew[0];
        atom->v[m][1] = vnew[1];
        atom->v[m][2] = vnew[2];
        
        atom->add_molecule_atom(onemols[imol],i,m,maxtag_all);
        modify->create_attribute(m);
      }
    }

    if (shakeflag) 
      fixshake->set_molecule(nlocalprev,maxtag_all,imol,com_coord,vnew,quat);

    atom->natoms += natoms_per_molecule;
    if (atom->natoms < 0 || atom->natoms > MAXBIGINT)
      error->all(FLERR,"Too many total atoms");
    atom->nbonds += onemols[imol]->nbonds;
    atom->nangles += onemols[imol]->nangles;
    atom->ndihedrals += onemols[imol]->ndihedrals;
    atom->nimpropers += onemols[imol]->nimpropers;
    if (atom->map_style) atom->map_init();
    atom->nghost = 0;
    comm->borders();
    update_gas_atoms_list();
    ninsertion_successes += 1.0;
  }
}

/* ----------------------------------------------------------------------
   compute particle's interaction energy with the rest of the system
------------------------------------------------------------------------- */

double FixGCMC::energy(int i, int itype, tagint imolecule, double *coord)
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
    if (mode == MOLECULE)
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

double FixGCMC::molecule_energy(tagint gas_molecule_id)
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
------------------------------------------------------------------------- */

void FixGCMC::attempt_atomic_translation_full()
{
  ntranslation_attempts += 1.0;
  
  if (ngas == 0) return;

  double energy_before = energy_stored;
  
  int i = pick_random_gas_atom();

  double **x = atom->x;
  double xtmp[3];
  
  tagint tmptag = -1;
    
  if (i >= 0) {
  
    double rsq = 1.1;
    double rx,ry,rz;
    rx = ry = rz = 0.0;
    double coord[3];
    while (rsq > 1.0) {
      rx = 2*random_unequal->uniform() - 1.0;
      ry = 2*random_unequal->uniform() - 1.0;
      rz = 2*random_unequal->uniform() - 1.0;
      rsq = rx*rx + ry*ry + rz*rz;
    }    
    coord[0] = x[i][0] + displace*rx;
    coord[1] = x[i][1] + displace*ry;
    coord[2] = x[i][2] + displace*rz;
    if (regionflag) {
      while (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) == 0) {
        while (rsq > 1.0) {
          rx = 2*random_unequal->uniform() - 1.0;
          ry = 2*random_unequal->uniform() - 1.0;
          rz = 2*random_unequal->uniform() - 1.0;
          rsq = rx*rx + ry*ry + rz*rz;
        }    
        coord[0] = x[i][0] + displace*rx;
        coord[1] = x[i][1] + displace*ry;
        coord[2] = x[i][2] + displace*rz;
      }
    }
    xtmp[0] = x[i][0];
    xtmp[1] = x[i][1];
    xtmp[2] = x[i][2];
    x[i][0] = coord[0];
    x[i][1] = coord[1];
    x[i][2] = coord[2];
    
    tmptag = atom->tag[i];
  }
  
  double energy_after = energy_full();
  
  if (random_equal->uniform() < 
      exp(beta*(energy_before - energy_after))) {
    energy_stored = energy_after;
    ntranslation_successes += 1.0;
  } else {
    for (int i = 0; i < atom->nlocal; i++) {
      if (tmptag == atom->tag[i]) { 
        x[i][0] = xtmp[0];
        x[i][1] = xtmp[1];
        x[i][2] = xtmp[2];
      }
    }
    energy_stored = energy_before;
  } 
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_atomic_deletion_full()
{
  double q_tmp;
  
  ndeletion_attempts += 1.0;

  if (ngas == 0) return;
  
  double energy_before = energy_stored;

  int i = pick_random_gas_atom();

  if (i >= 0) {
    atom->mask[i] |= exclusion_group_bit;
    if (atom->q_flag) {
      q_tmp = atom->q[i];
      atom->q[i] = 0.0;
    }
  }
  double energy_after = energy_full();

  if (random_equal->uniform() < 
      ngas*exp(beta*(energy_before - energy_after))/(zz*volume)) {
    if (i >= 0) {
      atom->avec->copy(atom->nlocal-1,i,1);
      atom->nlocal--;  
    }
    atom->natoms--;
    if (atom->map_style) atom->map_init();    
    ndeletion_successes += 1.0;   
    energy_stored = energy_after;    
  } else {
    if (i >= 0) {
      atom->mask[i] &= exclusion_group_inversebit;
      if (atom->q_flag) atom->q[i] = q_tmp;
    }
    energy_stored = energy_before;
  }
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_atomic_insertion_full()
{
  ninsertion_attempts += 1.0;

  double energy_before = energy_stored;
  
  double coord[3];
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
  } else {
    coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
    coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
    coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
  }
  
  int proc_flag = 0;
  if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
      coord[1] >= sublo[1] && coord[1] < subhi[1] &&
      coord[2] >= sublo[2] && coord[2] < subhi[2]) {
    proc_flag = 1;
    atom->avec->create_atom(ngcmc_type,coord);
    int m = atom->nlocal - 1;
    atom->mask[m] = 1 | groupbit;
    atom->v[m][0] = random_unequal->gaussian()*sigma;
    atom->v[m][1] = random_unequal->gaussian()*sigma;
    atom->v[m][2] = random_unequal->gaussian()*sigma;
    if (charge_flag) atom->q[m] = charge;
    modify->create_attribute(m);
  }

  atom->natoms++;
  if (atom->tag_enable) {
    atom->tag_extend();
    if (atom->map_style) atom->map_init();
  }
  atom->nghost = 0;
  comm->borders();
    
  double energy_after = energy_full();
  
  if (random_equal->uniform() <
      zz*volume*exp(beta*(energy_before - energy_after))/(ngas+1)) {

    ninsertion_successes += 1.0;
    energy_stored = energy_after;
  } else {
    atom->natoms--;
    if (proc_flag) atom->nlocal--;
    energy_stored = energy_before;
  }
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_molecule_translation_full()
{
  ntranslation_attempts += 1.0;

  if (ngas == 0) return;

  tagint translation_molecule = pick_random_gas_molecule();
  if (translation_molecule == -1) return;

  double energy_before = energy_stored;
  
  double **x = atom->x;
  double rx,ry,rz;
  double com_displace[3],coord[3];
  double rsq = 1.1;
  while (rsq > 1.0) {
    rx = 2*random_equal->uniform() - 1.0;
    ry = 2*random_equal->uniform() - 1.0;
    rz = 2*random_equal->uniform() - 1.0;
    rsq = rx*rx + ry*ry + rz*rz;
  }
  com_displace[0] = displace*rx;
  com_displace[1] = displace*ry;
  com_displace[2] = displace*rz;
  
  int nlocal = atom->nlocal;
  if (regionflag) {
    int *mask = atom->mask;
    for (int i = 0; i < nlocal; i++) {
      if (atom->molecule[i] == translation_molecule) {
        mask[i] |= molecule_group_bit;
      } else {
        mask[i] &= molecule_group_inversebit;
      }
    }
    double com[3];
    com[0] = com[1] = com[2] = 0.0;
    group->xcm(molecule_group,gas_mass,com);
    coord[0] = com[0] + displace*rx;
    coord[1] = com[1] + displace*ry;
    coord[2] = com[2] + displace*rz;
    while (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) == 0) {
      while (rsq > 1.0) {
        rx = 2*random_equal->uniform() - 1.0;
        ry = 2*random_equal->uniform() - 1.0;
        rz = 2*random_equal->uniform() - 1.0;
        rsq = rx*rx + ry*ry + rz*rz;
      }    
      coord[0] = com[0] + displace*rx;
      coord[1] = com[1] + displace*ry;
      coord[2] = com[2] + displace*rz;
    }
    com_displace[0] = displace*rx;
    com_displace[1] = displace*ry;
    com_displace[2] = displace*rz;
  }
  
  for (int i = 0; i < nlocal; i++) {
    if (atom->molecule[i] == translation_molecule) {
      x[i][0] += com_displace[0];
      x[i][1] += com_displace[1];
      x[i][2] += com_displace[2];
    }
  }

  double energy_after = energy_full();

  if (random_equal->uniform() < 
      exp(beta*(energy_before - energy_after))) {
    ntranslation_successes += 1.0;
    energy_stored = energy_after;
  } else {
    energy_stored = energy_before;
    for (int i = 0; i < nlocal; i++) {
      if (atom->molecule[i] == translation_molecule) {
        x[i][0] -= com_displace[0];
        x[i][1] -= com_displace[1];
        x[i][2] -= com_displace[2];
      }
    }
  }
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_molecule_rotation_full()
{
  nrotation_attempts += 1.0;

  if (ngas == 0) return;
  
  tagint rotation_molecule = pick_random_gas_molecule();
  if (rotation_molecule == -1) return;
  
  double energy_before = energy_stored;

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  for (int i = 0; i < nlocal; i++) {
    if (atom->molecule[i] == rotation_molecule) {
      mask[i] |= molecule_group_bit;
    } else {
      mask[i] &= molecule_group_inversebit;
    }
  }

  double com[3];
  com[0] = com[1] = com[2] = 0.0;
  group->xcm(molecule_group,gas_mass,com);

  double r[3],rotmat[3][3],quat[4];
  r[0] = random_equal->uniform() - 0.5;
  r[1] = random_equal->uniform() - 0.5;
  r[2] = random_equal->uniform() - 0.5;

  double theta = random_equal->uniform() * max_rotation_angle;
  MathExtra::norm3(r);
  MathExtra::axisangle_to_quat(r,theta,quat);
  MathExtra::quat_to_mat(quat,rotmat);

  double **x = atom->x;
  imageint *image = atom->image;
  imageint image_orig[natoms_per_molecule];
  int n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & molecule_group_bit) {
      atom_coord[n][0] = x[i][0];
      atom_coord[n][1] = x[i][1];
      atom_coord[n][2] = x[i][2]; 
      image_orig[n] = image[i];      
      double xtmp[3];
      domain->unmap(x[i],image[i],xtmp);
      xtmp[0] -= com[0];
      xtmp[1] -= com[1];
      xtmp[2] -= com[2];
      MathExtra::matvec(rotmat,xtmp,x[i]); 
      x[i][0] += com[0];
      x[i][1] += com[1];
      x[i][2] += com[2];
      image[i] = imagetmp;
      domain->remap(x[i],image[i]);
      n++;
    }
  }
  
  double energy_after = energy_full();

  if (random_equal->uniform() < 
      exp(beta*(energy_before - energy_after))) {
    nrotation_successes += 1.0;
    energy_stored = energy_after;    
  } else {
    energy_stored = energy_before;
    int n = 0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & molecule_group_bit) {
        x[i][0] = atom_coord[n][0];
        x[i][1] = atom_coord[n][1];
        x[i][2] = atom_coord[n][2];
        image[i] = image_orig[n];
        n++;
      }
    }
  }
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_molecule_deletion_full()
{
  ndeletion_attempts += 1.0;

  if (ngas == 0) return;
  
  tagint deletion_molecule = pick_random_gas_molecule();
  if (deletion_molecule == -1) return;

  double energy_before = energy_stored;
  
  int m = 0;
  double q_tmp[natoms_per_molecule];
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->molecule[i] == deletion_molecule) {
      atom->mask[i] |= exclusion_group_bit;
      if (atom->q_flag) {
        q_tmp[m] = atom->q[i];
        m++;
        atom->q[i] = 0.0;
      }
      
    }
  }
  
  double energy_after = energy_full();  

  if (random_equal->uniform() < 
      ngas*exp(beta*(energy_before - energy_after))/(zz*volume*natoms_per_molecule)) {
    int i = 0;
    while (i < atom->nlocal) {
      if (atom->molecule[i] == deletion_molecule) {
        atom->avec->copy(atom->nlocal-1,i,1);
        atom->nlocal--;
      } else i++;
    }
    atom->natoms -= natoms_per_molecule;
    if (atom->map_style) atom->map_init();
    ndeletion_successes += 1.0;
    energy_stored = energy_after;
  } else {
    energy_stored = energy_before;
    int m = 0;
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->molecule[i] == deletion_molecule) {
        atom->mask[i] &= exclusion_group_inversebit;
        if (atom->q_flag) {
          atom->q[i] = q_tmp[m];
          m++;
        }
      }    
    }
  }
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCMC::attempt_molecule_insertion_full()
{
  ninsertion_attempts += 1.0;
  
  double energy_before = energy_stored;
  
  tagint maxmol = 0;
  for (int i = 0; i < atom->nlocal; i++) maxmol = MAX(maxmol,atom->molecule[i]);
  tagint maxmol_all;
  MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
  maxmol_all++;
  if (maxmol_all >= MAXTAGINT) 
    error->all(FLERR,"Fix gcmc ran out of available molecule IDs");
  int insertion_molecule = maxmol_all;
    
  tagint maxtag = 0;
  for (int i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag,atom->tag[i]);
  tagint maxtag_all;
  MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

  int nfix = modify->nfix;
  Fix **fix = modify->fix;

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
  } else {
    com_coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
    com_coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
    com_coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
  }
  
  double r[3],rotmat[3][3],quat[4];
  r[0] = random_equal->uniform() - 0.5;
  r[1] = random_equal->uniform() - 0.5;
  r[2] = random_equal->uniform() - 0.5;

  double theta = random_equal->uniform() * MY_2PI;
  MathExtra::norm3(r);
  MathExtra::axisangle_to_quat(r,theta,quat);
  MathExtra::quat_to_mat(quat,rotmat);

  double vnew[3];
  vnew[0] = random_unequal->gaussian()*sigma;
  vnew[1] = random_unequal->gaussian()*sigma;
  vnew[2] = random_unequal->gaussian()*sigma;
    
  for (int i = 0; i < natoms_per_molecule; i++) {
    double xtmp[3];
    MathExtra::matvec(rotmat,onemols[imol]->x[i],xtmp);
    xtmp[0] += com_coord[0];
    xtmp[1] += com_coord[1];
    xtmp[2] += com_coord[2];
    
    domain->remap(xtmp);

    if (xtmp[0] >= sublo[0] && xtmp[0] < subhi[0] &&
        xtmp[1] >= sublo[1] && xtmp[1] < subhi[1] &&
        xtmp[2] >= sublo[2] && xtmp[2] < subhi[2]) {

      atom->avec->create_atom(ngcmc_type+onemols[imol]->type[i],xtmp);
      int m = atom->nlocal - 1;
      atom->mask[m] = 1 | groupbit;
      atom->image[m] = imagetmp;
      domain->remap(atom->x[m],atom->image[m]);
      atom->molecule[m] = insertion_molecule;
      if (maxtag_all+i+1 >= MAXTAGINT)
        error->all(FLERR,"Fix gcmc ran out of available atom IDs");
      atom->tag[m] = maxtag_all + i + 1;
      atom->v[m][0] = vnew[0];
      atom->v[m][1] = vnew[1];
      atom->v[m][2] = vnew[2];

      atom->add_molecule_atom(onemols[imol],i,m,maxtag_all);
      modify->create_attribute(m);
    }
  }

  if (shakeflag) 
    fixshake->set_molecule(nlocalprev,maxtag_all,imol,com_coord,vnew,quat);

  atom->natoms += natoms_per_molecule;
  if (atom->natoms < 0 || atom->natoms > MAXBIGINT)
    error->all(FLERR,"Too many total atoms");
  atom->nbonds += onemols[imol]->nbonds;
  atom->nangles += onemols[imol]->nangles;
  atom->ndihedrals += onemols[imol]->ndihedrals;
  atom->nimpropers += onemols[imol]->nimpropers;
  if (atom->map_style) atom->map_init();
  atom->nghost = 0;
  comm->borders();
  double energy_after = energy_full();

  if (random_equal->uniform() < zz*volume*natoms_per_molecule*
      exp(beta*(energy_before - energy_after)/(ngas + natoms_per_molecule))) {  
   
    ninsertion_successes += 1.0;
    energy_stored = energy_after;

  } else {

    atom->nbonds -= onemols[imol]->nbonds;
    atom->nangles -= onemols[imol]->nangles;
    atom->ndihedrals -= onemols[imol]->ndihedrals;
    atom->nimpropers -= onemols[imol]->nimpropers;
    atom->natoms -= natoms_per_molecule;
  
    energy_stored = energy_before;
    int i = 0;
    while (i < atom->nlocal) {
      if (atom->molecule[i] == insertion_molecule) {
        atom->avec->copy(atom->nlocal-1,i,1);
        atom->nlocal--;
      } else i++;
    }
  }
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixGCMC::energy_full()
{ 
  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build();
  int eflag = 1;
  int vflag = 0;
  
  if (modify->n_pre_force) modify->pre_force(vflag);
  
  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) force->kspace->compute(eflag,vflag);
 
  if (modify->n_post_force) modify->post_force(vflag);
  if (modify->n_end_of_step) modify->end_of_step();
 
  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();
  
  return total_energy;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixGCMC::pick_random_gas_atom()
{
  int i = -1;
  int iwhichglobal = static_cast<int> (ngas*random_equal->uniform());
  if ((iwhichglobal >= ngas_before) &&
      (iwhichglobal < ngas_before + ngas_local)) {
    int iwhichlocal = iwhichglobal - ngas_before;
    i = local_gas_list[iwhichlocal];
  }

  return i;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

tagint FixGCMC::pick_random_gas_molecule()
{
  int iwhichglobal = static_cast<int> (ngas*random_equal->uniform());
  tagint gas_molecule_id = 0;
  if ((iwhichglobal >= ngas_before) &&
      (iwhichglobal < ngas_before + ngas_local)) {
    int iwhichlocal = iwhichglobal - ngas_before;
    int i = local_gas_list[iwhichlocal];
    gas_molecule_id = atom->molecule[i];
  }

  tagint gas_molecule_id_all = 0;
  MPI_Allreduce(&gas_molecule_id,&gas_molecule_id_all,1,
                MPI_LMP_TAGINT,MPI_MAX,world);
  
  return gas_molecule_id_all;
}

/* ----------------------------------------------------------------------
   update the list of gas atoms
------------------------------------------------------------------------- */

void FixGCMC::update_gas_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;
      
  if (nlocal > gcmc_nmax) {
    memory->sfree(local_gas_list);
    gcmc_nmax = atom->nmax;
    local_gas_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "GCMC:local_gas_list");
  }

  ngas_local = 0;
  
  if (regionflag) {
  
    if (mode == MOLECULE) {
    
      tagint maxmol = 0;
      for (int i = 0; i < nlocal; i++) maxmol = MAX(maxmol,molecule[i]);
      tagint maxmol_all;
      MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
      double comx[maxmol_all];
      double comy[maxmol_all];
      double comz[maxmol_all];
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
  MPI_Scan(&ngas_local,&ngas_before,1,MPI_INT,MPI_SUM,world);
  ngas_before -= ngas_local;
}

/* ----------------------------------------------------------------------
  return acceptance ratios
------------------------------------------------------------------------- */

double FixGCMC::compute_vector(int n)
{
  if (n == 0) return ntranslation_attempts;
  if (n == 1) return ntranslation_successes;
  if (n == 2) return ninsertion_attempts;
  if (n == 3) return ninsertion_successes;
  if (n == 4) return ndeletion_attempts;
  if (n == 5) return ndeletion_successes;
  if (n == 6) return nrotation_attempts;
  if (n == 7) return nrotation_successes;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixGCMC::memory_usage()
{
  double bytes = gcmc_nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixGCMC::write_restart(FILE *fp)
{
  int n = 0;
  double list[4];
  list[n++] = random_equal->state();
  list[n++] = random_unequal->state();
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

void FixGCMC::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  random_equal->reset(seed);

  seed = static_cast<int> (list[n++]);
  random_unequal->reset(seed);

  next_reneighbor = static_cast<int> (list[n++]);
}
