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
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "group.h"
#include "domain.h"
#include "region.h"
#include "random_park.h"
#include "force.h"
#include "pair.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include <iostream>

using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixGCMC::FixGCMC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 11) error->all(FLERR,"Illegal fix gcmc command");

  vector_flag = 1;
  size_vector = 8;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // required args

  nevery = atoi(arg[3]);
  nexchanges = atoi(arg[4]);
  nmcmoves = atoi(arg[5]);
  ngcmc_type = atoi(arg[6]);
  seed = atoi(arg[7]);
  reservoir_temperature = atof(arg[8]);
  chemical_potential = atof(arg[9]);
  displace = atof(arg[10]);

  if (nexchanges < 0) error->all(FLERR,"Illegal fix gcmc command");
  if (nmcmoves < 0) error->all(FLERR,"Illegal fix gcmc command");
  if (seed <= 0) error->all(FLERR,"Illegal fix gcmc command");
  if (reservoir_temperature < 0.0)
    error->all(FLERR,"Illegal fix gcmc command");
  if (displace < 0.0) error->all(FLERR,"Illegal fix gcmc command");

  // set defaults

  molflag = 0;
  max_rotation_angle = 10*MY_PI/180;
  regionflag = 0; 
  iregion = -1; 
  region_volume = 0;
  max_region_attempts = 1000; 
  rotation_group = 0;
  rotation_groupbit = 0;
  rotation_inversegroupbit = 0;
  pressure_flag = false;
  pressure = 0.0;
  fugacity_coeff = 1.0;

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
      if (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) != 0) inside++;
    }

    double max_region_volume = (region_xhi - region_xlo)*
     (region_yhi - region_ylo)*(region_zhi - region_zlo);

    region_volume = max_region_volume*static_cast<double> (inside)/
     static_cast<double> (attempts);
  }

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

   atom_coord = NULL;
   model_atom = NULL;
   model_atom_buf = NULL;

}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixGCMC::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix gcmc command");

  int iarg = 0;
  while (iarg < narg) {
  if (strcmp(arg[iarg],"molecule") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      if (strcmp(arg[iarg+1],"no") == 0) molflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) molflag = 1;
      else error->all(FLERR,"Illegal fix gcmc command");
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
      max_rotation_angle = atof(arg[iarg+1]);
      max_rotation_angle *= MY_PI/180;
      iarg += 2;
    } else if (strcmp(arg[iarg],"pressure") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      pressure = atof(arg[iarg+1]);
      pressure_flag = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"fugacity_coeff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      fugacity_coeff = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix gcmc command");
  }
}

/* ---------------------------------------------------------------------- */

FixGCMC::~FixGCMC()
{
  if (regionflag) delete [] idregion;
  delete random_equal;
  delete random_unequal;
  memory->destroy(local_gas_list);
  memory->destroy(atom_coord);
  memory->destroy(model_atom_buf);
  delete model_atom;
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
  int *type = atom->type;

  if (molflag == 0) {
    if (ngcmc_type <= 0 || ngcmc_type > atom->ntypes)
      error->all(FLERR,"Invalid atom type in fix gcmc command");
  }

  // if molflag not set, warn if any deletable atom has a mol ID

  if (molflag == 0 && atom->molecule_flag) {
    int *molecule = atom->molecule;
    int *mask = atom->mask;
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

  // if molflag set, check for unset mol IDs

  if (molflag == 1) {
    int *molecule = atom->molecule;
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

  if ((molflag && (atom->molecule_flag == 0)) || 
      (molflag && ((!atom->tag_enable) || (!atom->map_style))))
    error->all(FLERR,
     "Fix gcmc molecule command requires that atoms have molecule attributes");

  if (force->pair->single_enable == 0)
    error->all(FLERR,"Fix gcmc incompatible with given pair_style");

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix gcmc in a 2d simulation");

  if (domain->triclinic == 1)
    error->all(FLERR,"Cannot use fix gcmc with a triclinic box");

  // create a new group for rotation molecules

  if (molflag) {
    char **group_arg = new char*[3];
    group_arg[0] = (char *) "rotation_gas_atoms";
    group_arg[1] = (char *) "molecule";
    char digits[12];
    sprintf(digits,"%d",ngcmc_type);
    group_arg[2] = digits;
    group->assign(3,group_arg);
    rotation_group = group->find(group_arg[0]);
    if (rotation_group == -1) 
      error->all(FLERR,"Could not find fix gcmc rotation group ID");
    rotation_groupbit = group->bitmask[rotation_group];
    rotation_inversegroupbit = rotation_groupbit ^ ~0;
    delete [] group_arg;
  }
    
  // get all of the needed molecule data if molflag, 
  // otherwise just get the gas mass
  
  if (molflag) get_model_molecule();
  else gas_mass = atom->mass[ngcmc_type];
  
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
  
  imagetmp = ((tagint) IMGMAX << IMG2BITS) | 
             ((tagint) IMGMAX << IMGBITS) | IMGMAX;
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

  if (molflag) {
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
    while (rsq > 1.0) {
      rx = 2*random_unequal->uniform() - 1.0;
      ry = 2*random_unequal->uniform() - 1.0;
      rz = 2*random_unequal->uniform() - 1.0;
      rsq = rx*rx + ry*ry + rz*rz;
    }
    double coord[3];
    coord[0] = x[i][0] + displace*rx;
    coord[1] = x[i][1] + displace*ry;
    coord[2] = x[i][2] + displace*rz;
    double energy_after = energy(i,ngcmc_type,-1,coord);
    if (random_unequal->uniform() < exp(-beta*(energy_after - energy_before))) {
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
    if (random_unequal->uniform() < ngas*exp(beta*deletion_energy)/(zz*volume)) {
      atom->avec->copy(atom->nlocal-1,i,1);
      atom->nlocal--;
      success = 1;
    }
  }

  int success_all = 0;
  MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);

  if (success_all) {
    if (atom->tag_enable) {
      atom->natoms--;
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

      int nfix = modify->nfix;
      Fix **fix = modify->fix;
      for (int j = 0; j < nfix; j++)
        if (fix[j]->create_attribute) fix[j]->set_arrays(m);

      success = 1;
    }
  }

  int success_all = 0;
  MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);

  if (success_all) {
    if (atom->tag_enable) {
      atom->natoms++;
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

  int translation_molecule = pick_random_gas_molecule();
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

  double energy_after = 0.0;
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->molecule[i] == translation_molecule) {
      coord[0] = x[i][0] + com_displace[0];
      coord[1] = x[i][1] + com_displace[1];
      coord[2] = x[i][2] + com_displace[2];
      energy_after += energy(i,atom->type[i],translation_molecule,coord);
    }
  }

  double energy_after_sum = 0.0;
  MPI_Allreduce(&energy_after,&energy_after_sum,1,MPI_DOUBLE,MPI_SUM,world);

  if (random_equal->uniform() < exp(-beta*(energy_after_sum - energy_before_sum))) {
    for (int i = 0; i < atom->nlocal; i++) {
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

  int rotation_molecule = pick_random_gas_molecule();
  if (rotation_molecule == -1) return;
  
  double energy_before_sum = molecule_energy(rotation_molecule);

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  for (int i = 0; i < nlocal; i++) {
    if (atom->molecule[i] == rotation_molecule) {
      mask[i] |= rotation_groupbit;
    } else {
      mask[i] &= rotation_inversegroupbit;
    }
  }

  double com[3];
  com[0] = com[1] = com[2] = 0.0;
  group->xcm(rotation_group,gas_mass,com);

  double rot[9];
  get_rotation_matrix(max_rotation_angle,&rot[0]);

  double **x = atom->x;
  tagint *image = atom->image;
  double energy_after = 0.0;
  int n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & rotation_groupbit) {
      double xtmp[3];
      domain->unmap(x[i],image[i],xtmp);
      xtmp[0] -= com[0];
      xtmp[1] -= com[1];
      xtmp[2] -= com[2];
      atom_coord[n][0] = rot[0]*xtmp[0] + rot[1]*xtmp[1] + rot[2]*xtmp[2] + com[0];
      atom_coord[n][1] = rot[3]*xtmp[0] + rot[4]*xtmp[1] + rot[5]*xtmp[2] + com[1];
      atom_coord[n][2] = rot[6]*xtmp[0] + rot[7]*xtmp[1] + rot[8]*xtmp[2] + com[2];
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

  if (random_equal->uniform() < exp(-beta*(energy_after_sum - energy_before_sum))) {
    int n = 0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & rotation_groupbit) {
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
  
  int deletion_molecule = pick_random_gas_molecule();
  if (deletion_molecule == -1) return;

  double deletion_energy_sum = molecule_energy(deletion_molecule);

  if (random_equal->uniform() < ngas*exp(beta*deletion_energy_sum)/(zz*volume*natoms_per_molecule)) {
    int i = 0;
    while (i < atom->nlocal) {
      if (atom->molecule[i] == deletion_molecule) {
        atom->avec->copy(atom->nlocal-1,i,1);
        atom->nlocal--;
      } else i++;
    }
    atom->natoms -= natoms_per_molecule;
    atom->map_init();
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

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double com_coord[3];
  if (regionflag) {
    int region_attempt = 0;
    com_coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
    com_coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
    com_coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
    while (domain->regions[iregion]->match(com_coord[0],com_coord[1],com_coord[2]) == 0) {
      com_coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      com_coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      com_coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      region_attempt++;
      if (region_attempt >= max_region_attempts) return;
    }
  } else {
    com_coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
    com_coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
    com_coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
  }

  double rot[9];
  get_rotation_matrix(MY_2PI,&rot[0]);

  double **model_x = model_atom->x;
  double insertion_energy = 0.0;
  bool procflag[natoms_per_molecule];
  for (int i = 0; i < natoms_per_molecule; i++) {
    atom_coord[i][0] = rot[0]*model_x[i][0] + rot[1]*model_x[i][1] + rot[2]*model_x[i][2] + com_coord[0];
    atom_coord[i][1] = rot[3]*model_x[i][0] + rot[4]*model_x[i][1] + rot[5]*model_x[i][2] + com_coord[1];
    atom_coord[i][2] = rot[6]*model_x[i][0] + rot[7]*model_x[i][1] + rot[8]*model_x[i][2] + com_coord[2];

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
      insertion_energy += energy(-1,model_atom->type[i],-1,xtmp);
    }
  }

  double insertion_energy_sum = 0.0;
  MPI_Allreduce(&insertion_energy,&insertion_energy_sum,1,MPI_DOUBLE,MPI_SUM,world);

  if (random_equal->uniform() < zz*volume*natoms_per_molecule*exp(-beta*insertion_energy_sum)/(ngas+1)) {  
    maxmol++;
    if (maxmol >= MAXSMALLINT) 
      error->all(FLERR,"Fix gcmc ran out of available molecule IDs");

    int maxtag = 0;
    for (int i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag,atom->tag[i]);
    int maxtag_all;
    MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_INT,MPI_MAX,world);
    int atom_offset = maxtag_all;

    int k = 0;
    double **x = atom->x;
    double **v = atom->v;
    tagint *image = atom->image;
    int *molecule = atom->molecule;
    int *tag = atom->tag;
    for (int i = 0; i < natoms_per_molecule; i++) {
      k += atom->avec->unpack_exchange(&model_atom_buf[k]);
      if (procflag[i]) {
        int m = atom->nlocal - 1;
        image[m] = imagetmp;
        x[m][0] = atom_coord[i][0];
        x[m][1] = atom_coord[i][1];
        x[m][2] = atom_coord[i][2];
        domain->remap(x[m],image[m]);
        atom->molecule[m] = maxmol;
        tag[m] += atom_offset;
        v[m][0] = random_unequal->gaussian()*sigma;
        v[m][1] = random_unequal->gaussian()*sigma;
        v[m][2] = random_unequal->gaussian()*sigma;
        
        if (atom->avec->bonds_allow)
          for (int j = 0; j < atom->num_bond[m]; j++)
            atom->bond_atom[m][j] += atom_offset;
        if (atom->avec->angles_allow)
          for (int j = 0; j < atom->num_angle[m]; j++) {
            atom->angle_atom1[m][j] += atom_offset;
            atom->angle_atom2[m][j] += atom_offset;
            atom->angle_atom3[m][j] += atom_offset;
          }
        if (atom->avec->dihedrals_allow)
          for (int j = 0; j < atom->num_dihedral[m]; j++) {
            atom->dihedral_atom1[m][j] += atom_offset;
            atom->dihedral_atom2[m][j] += atom_offset;
            atom->dihedral_atom3[m][j] += atom_offset;
            atom->dihedral_atom4[m][j] += atom_offset;
          }
        if (atom->avec->impropers_allow)
          for (int j = 0; j < atom->num_improper[m]; j++) {
            atom->improper_atom1[m][j] += atom_offset;
            atom->improper_atom2[m][j] += atom_offset;
            atom->improper_atom3[m][j] += atom_offset;
            atom->improper_atom4[m][j] += atom_offset;
          }

        for (int j = 0; j < atom->nspecial[m][2]; j++)
          atom->special[m][j] += atom_offset;

        int nfix = modify->nfix;
        Fix **fix = modify->fix;
        for (int j = 0; j < nfix; j++) {
          if (strcmp(modify->fix[j]->style,"shake") == 0) {
            fix[j]->update_arrays(m,atom_offset);
          } else if (fix[j]->create_attribute) fix[j]->set_arrays(m);
        }

      } else atom->nlocal--;
    }
    atom->natoms += natoms_per_molecule;
    atom->map_init();
    atom->nghost = 0;
    comm->borders();
    update_gas_atoms_list();
    ninsertion_successes += 1.0;
  }
}

/* ----------------------------------------------------------------------
   compute particle's interaction energy with the rest of the system
------------------------------------------------------------------------- */

double FixGCMC::energy(int i, int itype, int imolecule, double *coord)
{
  double delx,dely,delz,rsq;

  double **x = atom->x;
  int *type = atom->type;
  int *molecule = atom->molecule;
  int nall = atom->nlocal + atom->nghost;
  pair = force->pair;
  cutsq = force->pair->cutsq;

  double fpair = 0.0;
  double factor_coul = 1.0;
  double factor_lj = 1.0;

  double total_energy = 0.0;
  for (int j = 0; j < nall; j++) {

    if (i == j) continue;
    if (molflag)
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

int FixGCMC::pick_random_gas_molecule()
{
  int iwhichglobal = static_cast<int> (ngas*random_equal->uniform());
  int gas_molecule_id = 0;
  if ((iwhichglobal >= ngas_before) &&
      (iwhichglobal < ngas_before + ngas_local)) {
    int iwhichlocal = iwhichglobal - ngas_before;
    int i = local_gas_list[iwhichlocal];
    gas_molecule_id = atom->molecule[i];
  }

  int gas_molecule_id_all = 0;
  MPI_Allreduce(&gas_molecule_id,&gas_molecule_id_all,1,MPI_INT,MPI_MAX,world);
  
  return gas_molecule_id_all;
}

/* ----------------------------------------------------------------------
   compute the energy of the given gas molecule in its current position 
   sum across all procs that own atoms of the given molecule
------------------------------------------------------------------------- */

double FixGCMC::molecule_energy(int gas_molecule_id)
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
   compute a 3x3 rotation matrix using 3 random Euler angles, 
   each with a random maximum value supplied by the caller
------------------------------------------------------------------------- */

void FixGCMC::get_rotation_matrix(double max_angle, double *rot)
{
  double angle_x = max_angle*random_equal->uniform();
  double angle_y = max_angle*random_equal->uniform();
  double angle_z = max_angle*random_equal->uniform();
  
  double a = cos(angle_x);
  double b = sin(angle_x);
  double c = cos(angle_y);
  double d = sin(angle_y);
  double e = cos(angle_z);
  double f = sin(angle_z);
  double ad = a*d;
  double bd = b*d;

  rot[0] = c*e;
  rot[1] = -c*f;
  rot[2] = -d;
  rot[3] = -bd*e + a*f;
  rot[4] = bd*f + a*e;
  rot[5] = -b*c;
  rot[6] = ad*e + b*f;
  rot[7] = -ad*f + b*e;
  rot[8] = a*c;
}

/* ----------------------------------------------------------------------
   when using the molecule capability, populate model atom arrays from
   the model molecule provided by the user that will then be used to build 
   inserted molecules
------------------------------------------------------------------------- */

void FixGCMC::get_model_molecule()
{
  // find out how many atoms are in the model molecule
  // just loop through all of the atoms I own, then sum up across procs
  
  int model_molecule_number = ngcmc_type;
  int natoms_per_molecule_local = 0;
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->molecule[i] == model_molecule_number) {
      natoms_per_molecule_local++;
    }
  }

  natoms_per_molecule = 0;
  MPI_Allreduce(&natoms_per_molecule_local,&natoms_per_molecule,1,MPI_INT,MPI_SUM,world);

  if (natoms_per_molecule == 0)
    error->all(FLERR,"Fix gcmc could not find any atoms in the user-supplied template molecule");
  
  memory->create(atom_coord,natoms_per_molecule,3,"fixGCMC:atom_coord");

  // maxmol = largest molecule tag across all existing atoms

  maxmol = 0;
  if (atom->molecular) {
    for (int i = 0; i < atom->nlocal; i++) maxmol = MAX(atom->molecule[i],maxmol);
    int maxmol_all;
    MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_INT,MPI_MAX,world);
    maxmol = maxmol_all;
  }

  // communication buffer for model atom's info
  // max_size = largest buffer needed by any proc
  // must do before new Atom class created,
  //   since size_restart() uses atom->nlocal

  int max_size;
  int buf_send_size = atom->avec->size_restart(); 
  
  MPI_Allreduce(&buf_send_size,&max_size,1,MPI_INT,MPI_MAX,world);
  max_size *= 2;
  double *buf;
  memory->create(buf,max_size,"fixGCMC:buf");
  
  // create storage space for the model molecule's atoms
  // create a new atom object called atom to store the data
  
  // old_atom = original atom class
  // atom = new model atom class
  // if old_atom style was hybrid, pass sub-style names to create_avec

  Atom *old_atom = atom;
  atom = new Atom(lmp);
  atom->settings(old_atom);
  
  int nstyles = 0;
  char **keywords = NULL;
  if (strcmp(old_atom->atom_style,"hybrid") == 0) {
    AtomVecHybrid *avec_hybrid = (AtomVecHybrid *) old_atom->avec;
    nstyles = avec_hybrid->nstyles;
    keywords = avec_hybrid->keywords;
  }
  
  atom->create_avec(old_atom->atom_style,nstyles,keywords);
  
  // assign atom and topology counts in model atom class from old_atom

  atom->ntypes = old_atom->ntypes;
  atom->nbondtypes = old_atom->nbondtypes;
  atom->nangletypes = old_atom->nangletypes;
  atom->ndihedraltypes = old_atom->ndihedraltypes;
  atom->nimpropertypes = old_atom->nimpropertypes;
  atom->bond_per_atom = old_atom->bond_per_atom;
  atom->angle_per_atom = old_atom->angle_per_atom;
  atom->dihedral_per_atom = old_atom->dihedral_per_atom;
  atom->improper_per_atom = old_atom->improper_per_atom;
  atom->maxspecial = old_atom->maxspecial;
  atom->nextra_grow = old_atom->nextra_grow;

  if (atom->nextra_grow) {
    memory->grow(atom->extra_grow,old_atom->nextra_grow_max,"fixGCMC:extra_grow");
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      atom->extra_grow[iextra] = old_atom->extra_grow[iextra];
  }

  atom->extra_bond_per_atom = old_atom->extra_bond_per_atom;
  atom->allocate_type_arrays();
  atom->avec->grow(natoms_per_molecule + old_atom->nlocal);

  // copy type arrays to model atom class
  
  if (atom->mass) {
    for (int itype = 1; itype <= atom->ntypes; itype++) {
      atom->mass_setflag[itype] = old_atom->mass_setflag[itype];
      if (atom->mass_setflag[itype]) atom->mass[itype] = old_atom->mass[itype];
    }
  }
  // loop over all procs
  // if this iteration of loop is me:
  //   pack my atom data into buf
  //   bcast it to all other procs

  AtomVec *old_avec = old_atom->avec;
  AtomVec *model_avec = atom->avec;

  int model_buf_size = 0;
  for (int iproc = 0; iproc < comm->nprocs; iproc++) {
    int nbuf_iproc = 0;
    if (comm->me == iproc) {
      for (int i = 0; i < old_atom->nlocal; i++) {
        if (old_atom->molecule[i] == model_molecule_number) {
          nbuf_iproc += old_avec->pack_exchange(i,&buf[nbuf_iproc]);
        }
      }
    }
    MPI_Bcast(&nbuf_iproc,1,MPI_INT,iproc,world);
    MPI_Bcast(buf,nbuf_iproc,MPI_DOUBLE,iproc,world);

    model_buf_size += nbuf_iproc;
     
    int m = 0;
    while (m < nbuf_iproc)
      m += model_avec->unpack_exchange(&buf[m]);
  }

  // free communication buffer 

  memory->destroy(buf);
  
  // make sure that the number of model atoms is equal to the number of atoms per gas molecule
  
  int nlocal = atom->nlocal;
  if (nlocal != natoms_per_molecule)
    error->all(FLERR,"Fix gcmc incorrect number of atoms per molecule");
  
  // compute the model molecule's mass and center-of-mass
  // then recenter model molecule on the origin

  double com[3]; 
  gas_mass = group->mass(0);
  group->xcm(0,gas_mass,com);

  double **x = atom->x;  
  for (int i = 0; i < nlocal; i++) {
    domain->unmap(x[i],atom->image[i]);
    x[i][0] -= com[0];
    x[i][1] -= com[1];
    x[i][2] -= com[2];
  }

  int mintag = atom->tag[0];
  for (int i = 0; i < atom->nlocal; i++) mintag = MIN(mintag,atom->tag[i]);
  int atom_offset = mintag - 1;
  
  for (int i = 0; i < nlocal; i++) {
    atom->mask[i] = 1 | groupbit;
    atom->tag[i] -= atom_offset;
    if (atom->avec->bonds_allow)
      for (int j = 0; j < atom->num_bond[i]; j++)
        atom->bond_atom[i][j] -= atom_offset;
    if (atom->avec->angles_allow)
      for (int j = 0; j < atom->num_angle[i]; j++) {
        atom->angle_atom1[i][j] -= atom_offset;
        atom->angle_atom2[i][j] -= atom_offset;
        atom->angle_atom3[i][j] -= atom_offset;
      }
    if (atom->avec->dihedrals_allow)
      for (int j = 0; j < atom->num_dihedral[i]; j++) {
        atom->dihedral_atom1[i][j] -= atom_offset;
        atom->dihedral_atom2[i][j] -= atom_offset;
        atom->dihedral_atom3[i][j] -= atom_offset;
        atom->dihedral_atom4[i][j] -= atom_offset;
      }
    if (atom->avec->impropers_allow)
      for (int j = 0; j < atom->num_improper[i]; j++) {
        atom->improper_atom1[i][j] -= atom_offset;
        atom->improper_atom2[i][j] -= atom_offset;
        atom->improper_atom3[i][j] -= atom_offset;
        atom->improper_atom4[i][j] -= atom_offset;
      }
    for (int j = 0; j < atom->nspecial[i][2]; j++)
      atom->special[i][j] -= atom_offset;
  }

  // pack model atoms into a buffer for use during molecule insertions
  
  memory->create(model_atom_buf,model_buf_size,"fixGCMC:model_atom_buf");
  int n = 0;
  for (int i = 0; i < nlocal; i++) 
    n += model_avec->pack_exchange(i,&model_atom_buf[n]);

  // move atom to model_atom and restore old_atom class pointer back to atom

  model_atom = atom;
  atom = old_atom;
}

/* ----------------------------------------------------------------------
   update the list of gas atoms
------------------------------------------------------------------------- */

void FixGCMC::update_gas_atoms_list()
{
  if (atom->nlocal > gcmc_nmax) {
    memory->sfree(local_gas_list);
    gcmc_nmax = atom->nmax;
    local_gas_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "GCMC:local_gas_list");
  }

  ngas_local = 0;
  if (regionflag) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->mask[i] & groupbit) {
        double **x = atom->x;
        if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]) == 1) {
          local_gas_list[ngas_local] = i;
          ngas_local++;
        }
      }
    }
  } else {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->mask[i] & groupbit) {
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
