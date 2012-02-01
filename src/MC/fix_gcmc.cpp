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
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "group.h"
#include "domain.h"
#include "random_park.h"
#include "force.h"
#include "pair.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixGCMC::FixGCMC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 11) error->all(FLERR,"Illegal fix GCMC command");

  vector_flag = 1;
  size_vector = 6;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // required args
  
  nevery = atoi(arg[3]);
  nexchanges = atoi(arg[4]);
  nmcmoves = atoi(arg[5]);
  ntype = atoi(arg[6]);
  seed = atoi(arg[7]);
  reservoir_temperature = atof(arg[8]);
  chemical_potential = atof(arg[9]);
  displace = atof(arg[10]);

  if (ntype <= 0 || ntype > atom->ntypes) 
    error->all(FLERR,"Invalid atom type in fix GCMC command");
  if (nexchanges < 0) error->all(FLERR,"Illegal fix GCMC command");
  if (nmcmoves < 0) error->all(FLERR,"Illegal fix GCMC command");
  if (seed <= 0) error->all(FLERR,"Illegal fix GCMC command");
  if (reservoir_temperature < 0.0) error->all(FLERR,"Illegal fix GCMC command");  
  if (displace < 0.0) error->all(FLERR,"Illegal fix GCMC command"); 

  // compute beta, lambda, sigma, and the zz factor
  
  beta = 1.0/(force->boltz*reservoir_temperature);
  double gas_mass = atom->mass[ntype];
  double lambda = sqrt(force->hplanck*force->hplanck/
		       (2.0*MY_PI*gas_mass*force->mvv2e*
			force->boltz*reservoir_temperature));
  sigma = sqrt(force->boltz*reservoir_temperature/gas_mass/force->mvv2e);
  zz = exp(beta*chemical_potential)/(pow(lambda,3));

  // set defaults

  molflag = 0;

  // read options from end of input line

  options(narg-11,&arg[11]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp,seed);

  // random number generator, not the same for all procs

  random_unequal = new RanPark(lmp,seed);
  
  // compute the number of MC cycles that occur nevery timesteps
  
  ncycles = nexchanges + nmcmoves;
  
  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters
  
  nmove_attempts = 0.0;   
  nmove_successes = 0.0;  
  ndel_attempts = 0.0;    
  ndel_successes = 0.0;   
  ninsert_attempts = 0.0; 
  ninsert_successes = 0.0;

  nmax = 0;
  local_gas_list = NULL;
}

/* ---------------------------------------------------------------------- */

FixGCMC::~FixGCMC()
{
  delete random_equal;
  delete random_unequal;
  memory->sfree(local_gas_list);
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
  // check that no deletable atoms are in atom->firstgroup
  // deleting such an atom would not leave firstgroup atoms first

  int *type = atom->type;
  
  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < nlocal; i++)
      if ((type[i] == ntype) && (mask[i] && firstgroupbit)) flag = 1;

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

    if (flagall)
      error->all(FLERR,"Cannot do GCMC on atoms in atom_modify first group");
  }

  // if molflag not set, warn if any deletable atom has a mol ID

  if (molflag == 0 && atom->molecule_flag) {
    int *molecule = atom->molecule;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int flag = 0;
    for (int i = 0; i < nlocal; i++)
      if (type[i] == ntype)
        if (molecule[i]) flag = 1;
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall && comm->me == 0)
      error->warning(FLERR,"Fix GCMC may delete atom with non-zero molecule ID");
  }

  if (molflag && atom->molecule_flag == 0)
      error->all(FLERR,"Fix GCMC molecule command requires atom attribute molecule");
      
  if (molflag != 0) error->all(FLERR,"Fix GCMC molecule feature does not yet work"); 
  
  if (force->pair->single_enable == 0) 
    error->all(FLERR,"Fix GCMC incompatible with given pair_style");
}

/* ----------------------------------------------------------------------
   attempt particle insertions and deletions
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixGCMC::pre_exchange()
{
  // just return if should not be called on this timestep
  
  if (next_reneighbor != update->ntimestep) return;

  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  } 
  
  volume = domain->xprd * domain->yprd * domain->zprd;
  
  // grow local_gas_list array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(local_gas_list);
    nmax = atom->nmax;
    local_gas_list = (int *) memory->smalloc(nmax*sizeof(int),
					     "GCMC:local_gas_list");
  }
  
  int *type = atom->type;
  ngas_local = 0;
  for (int i = 0; i < atom->nlocal; i++)
    if (type[i] == ntype)
      local_gas_list[ngas_local++] = i;
                                           
  MPI_Allreduce(&ngas_local,&ngas,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&ngas_local,&ngas_before,1,MPI_INT,MPI_SUM,world);
  ngas_before -= ngas_local;
      
  // perform ncycles MC cycles
  
  for (int i = 0; i < ncycles; i++) {
    int random_int_fraction = 
      static_cast<int>(random_equal->uniform()*ncycles) + 1;	
    if (random_int_fraction <= nmcmoves) {
      attempt_move();
    } else {
      if (random_equal->uniform() < 0.5) attempt_deletion();
      else attempt_insertion();
    }
  }
  
  next_reneighbor = update->ntimestep + nevery;
} 

/* ----------------------------------------------------------------------
   choose particle randomly across all procs and attempt displacement
------------------------------------------------------------------------- */

void FixGCMC::attempt_move()
{  
  int i,iwhichglobal,iwhichlocal;
  double rx,ry,rz;
  double coord[3];
  double **x = atom->x;

  nmove_attempts += 1.0;
  
  int success = 0;
  iwhichglobal = static_cast<int> (ngas*random_equal->uniform());
  if ((iwhichglobal >= ngas_before) && 
      (iwhichglobal < ngas_before + ngas_local)) {
    iwhichlocal = iwhichglobal - ngas_before;
    i = local_gas_list[iwhichlocal];
    double energy_before = energy(i,x[i]);
    double rsq = 1.1;
    while (rsq > 1.0) {
      rx = 2*random_unequal->uniform() - 1.0;
      ry = 2*random_unequal->uniform() - 1.0;
      if (domain->dimension == 3) rz = 2*random_unequal->uniform() - 1.0;
      else rz = 0.0;
      rsq = rx*rx + ry*ry + rz*rz;
    }
    coord[0] = x[i][0] + displace*rx;
    coord[1] = x[i][1] + displace*ry;
    coord[2] = x[i][2] + displace*rz;
    double energy_after = energy(i,coord);
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
    nmove_successes += 1.0;  
    comm->borders(); 
  }
  
}

/* ----------------------------------------------------------------------
   attempt particle deletion
------------------------------------------------------------------------- */

void FixGCMC::attempt_deletion()
{ 
  ndel_attempts += 1.0;

  if (ngas == 0) return;

  int i,iwhichglobal,iwhichlocal;
  AtomVec *avec = atom->avec;

  // choose particle randomly across all procs and delete it
  // keep ngas, ngas_local, ngas_before, and local_gas_list current
  // after each deletion
  
  int success = 0;
  iwhichglobal = static_cast<int> (ngas*random_equal->uniform());
  if ((iwhichglobal >= ngas_before) && 
      (iwhichglobal < ngas_before + ngas_local)) {
    iwhichlocal = iwhichglobal - ngas_before;
    i = local_gas_list[iwhichlocal];
    double deletion_energy = energy(i,atom->x[i]);      
    if (random_unequal->uniform() < 
	ngas*exp(beta*deletion_energy)/(zz*volume)) {
      avec->copy(atom->nlocal-1,i,1);
      atom->nlocal--;
      local_gas_list[iwhichlocal] = local_gas_list[ngas_local-1];
      ngas_local--;
      success = 1; 
    }
  }
 
  int success_all = 0;
  MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);
  
  if (success_all) {
    ngas--;
    ndel_successes += 1.0;
    atom->natoms--;
    if (iwhichglobal < ngas_before) ngas_before--;      
    comm->borders();
    if (atom->tag_enable) {
      atom->tag_extend();
      if (atom->map_style) {
        atom->map_init();
        atom->map_set();
      }
    }
  }
}

/* ----------------------------------------------------------------------
   attempt particle insertion
------------------------------------------------------------------------- */

void FixGCMC::attempt_insertion()
{  
  int flag,success;
  double coord[3],lamda[3];
  double *newcoord; 
  
  ninsert_attempts += 1.0;
  
  // choose random position for new atom within box

  coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
  coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
  coord[2] = zlo + random_equal->uniform() * (zhi-zlo);

  // check if new atom is in my sub-box or above it if I'm highest proc
  // if so, add to my list via create_atom()
  // initialize info about the atoms
  // set group mask to "all" plus fix group
  
  if (domain->triclinic) {
    domain->x2lamda(coord,lamda);
    newcoord = lamda;
  } else newcoord = coord;

  flag = 0;
  if (newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
    newcoord[1] >= sublo[1] && newcoord[1] < subhi[1] &&
    newcoord[2] >= sublo[2] && newcoord[2] < subhi[2]) flag = 1;
    
  success = 0;
  if (flag) {
    int nall = atom->nlocal + atom->nghost;
    double insertion_energy = energy(nall,coord);
    if (random_unequal->uniform() < 
	zz*volume*exp(-beta*insertion_energy)/(ngas+1)) {
      atom->avec->create_atom(ntype,coord);    
      int m = atom->nlocal - 1;
      atom->type[m] = ntype;
      atom->mask[m] = 1 | groupbit;
      atom->v[m][0] = random_unequal->gaussian()*sigma;
      atom->v[m][1] = random_unequal->gaussian()*sigma;
      atom->v[m][2] = random_unequal->gaussian()*sigma;
      
      int nfix = modify->nfix;
      Fix **fix = modify->fix;
      for (int j = 0; j < nfix; j++)
        if (fix[j]->create_attribute) fix[j]->set_arrays(m);

      if (atom->nlocal > nmax) {
        nmax = atom->nmax;
        local_gas_list = (int *) 
	  memory->srealloc(local_gas_list,nmax*sizeof(int),
			   "GCMC:local_gas_list");
      }
      
      local_gas_list[ngas_local] = atom->nlocal;
      ngas_local++; 
      success = 1; 
    }
  }
  
  int success_all = 0;
  MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);
  
  if (success_all) {
    ngas++;
    ninsert_successes += 1.0;
    MPI_Scan(&ngas_local,&ngas_before,1,MPI_INT,MPI_SUM,world);
    ngas_before -= ngas_local;
    comm->borders();
    if (atom->tag_enable) {
      atom->natoms++;
      atom->tag_extend();
      if (atom->map_style) {
        atom->map_init();
        atom->map_set();
      }
    }   
  } 
}

/* ----------------------------------------------------------------------
   compute particle's interaction energy with the rest of the system
------------------------------------------------------------------------- */

double FixGCMC::energy(int i, double *coord) 
{
  double delx,dely,delz,rsq;

  double **x = atom->x;
  int *type = atom->type;
  int nall = atom->nlocal + atom->nghost;
  pair = force->pair;
  cutsq = force->pair->cutsq;

  double fpair = 0.0;
  double factor_coul = 1.0;
  double factor_lj = 1.0;
  
  double total_energy = 0.0;
  for (int j = 0; j < nall; j++) {

    if (i == j) continue;

    delx = coord[0] - x[j][0];
    dely = coord[1] - x[j][1];
    delz = coord[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    int jtype = type[j];

    if (rsq < cutsq[ntype][jtype])
      total_energy += 
	pair->single(i,j,ntype,jtype,rsq,factor_coul,factor_lj,fpair);
  }

  return total_energy;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line 
------------------------------------------------------------------------- */

void FixGCMC::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix GCMC command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"molecule") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix GCMC command");
      if (strcmp(arg[iarg+1],"no") == 0) molflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) molflag = 1;
      else error->all(FLERR,"Illegal fix evaporate command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix GCMC command");
  }
}

/* ----------------------------------------------------------------------
  return acceptance ratios
------------------------------------------------------------------------- */

double FixGCMC::compute_vector(int n)
{
  if (n == 0) return nmove_attempts;
  if (n == 1) return nmove_successes;
  if (n == 2) return ndel_attempts;
  if (n == 3) return ndel_successes;
  if (n == 4) return ninsert_attempts;
  if (n == 5) return ninsert_successes;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixGCMC::memory_usage()
{
  double bytes = nmax * sizeof(int);
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
