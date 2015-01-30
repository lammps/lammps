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
   Contributing authors: Paul Crozier (SNL)
                         Alexander Stukowski
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_atom_swap.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_hybrid.h"
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

/* ---------------------------------------------------------------------- */

FixAtomSwap::FixAtomSwap(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 10) error->all(FLERR,"Illegal fix atom/swap command");

  dynamic_group_allow = 1;
    
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // required args

  nevery = force->inumeric(FLERR,arg[3]);
  ncycles = force->inumeric(FLERR,arg[4]);
  seed = force->inumeric(FLERR,arg[5]); 
  double temperature = force->numeric(FLERR,arg[6]);
  beta = 1.0/(force->boltz*temperature);   

  if (ncycles < 0) error->all(FLERR,"Illegal fix atom/swap command");
  if (seed <= 0) error->all(FLERR,"Illegal fix atom/swap command");

  memory->create(delta_mu,atom->ntypes+1,"atom/swap:delta_mu");
  for (int i = 1; i <= atom->ntypes; i++) delta_mu[i] = 0.0;
  
  // read options from end of input line

  options(narg-7,&arg[7]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp,seed);
  
  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters

  nswap_attempts = 0.0;
  nswap_successes = 0.0;

  atom_swap_nmax = 0;
  local_swap_atom_list = NULL;

  // set comm size needed by this Fix

  if (atom->q_flag) comm_forward = 2;
  else comm_forward = 1;

}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixAtomSwap::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix atom/swap command");

  regionflag = 0; 
  conserve_ke_flag = 1;
  semi_grand_flag = 0;
  iregion = -1; 
  
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix atom/swap command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix atom/swap does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      regionflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"ke") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix atom/swap command");
      if (strcmp(arg[iarg+1],"no") == 0) conserve_ke_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) conserve_ke_flag = 1;
      else error->all(FLERR,"Illegal fix atom/swap command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"semi-grand") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix atom/swap command");
      if (strcmp(arg[iarg+1],"no") == 0) semi_grand_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) semi_grand_flag = 1;
      else error->all(FLERR,"Illegal fix atom/swap command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"types") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix atom/swap command");
      iarg++;
      memory->create(type_list,atom->ntypes,"atom/swap:type_list");
      nswaptypes = 0;
      while (iarg < narg) {
        if (isalpha(arg[iarg][0])) break;
        type_list[nswaptypes++] = force->numeric(FLERR,arg[iarg]);
        iarg++;
      }
    } else if (strcmp(arg[iarg],"delta_mu") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix atom/swap command");
      iarg++;
      ndeltamutypes = 0;
      while (iarg < narg) {
        if (isalpha(arg[iarg][0])) break;
        delta_mu[ndeltamutypes+2] = force->numeric(FLERR,arg[iarg]);
        ndeltamutypes++;
        iarg++;
      }
    } else error->all(FLERR,"Illegal fix atom/swap command");
  }
}

/* ---------------------------------------------------------------------- */

FixAtomSwap::~FixAtomSwap()
{
  memory->destroy(type_list);
  memory->destroy(delta_mu);
  memory->destroy(qtype);
  memory->destroy(sqrt_mass_ratio);
  if (regionflag) delete [] idregion;
  delete random_equal;
}

/* ---------------------------------------------------------------------- */

int FixAtomSwap::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAtomSwap::init()
{ 
  char *id_pe = (char *) "thermo_pe";
  int ipe = modify->find_compute(id_pe);
  c_pe = modify->compute[ipe];

  int *type = atom->type;

  if (nswaptypes < 2)
    error->all(FLERR,"Must specify at least 2 types in fix atom/swap command");
    
  if (semi_grand_flag) {
    if (atom->ntypes-1 != ndeltamutypes)
      error->all(FLERR,"Need ntypes-1 delta_mu values in fix atom/swap command");
  } else {
    if (nswaptypes != 2)
      error->all(FLERR,"Only 2 types allowed when not using semi-grand in fix atom/swap command");
    if (ndeltamutypes != 0)
      error->all(FLERR,"Delta_mu not allowed when not using semi-grand in fix atom/swap command");
  }
  
  for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++)
    if (type_list[iswaptype] <= 0 || type_list[iswaptype] > atom->ntypes)
      error->all(FLERR,"Invalid atom type in fix atom/swap command");

  memory->create(qtype,nswaptypes,"atom/swap:qtype");
  for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++) {
    if (atom->q_flag) {
      bool first = true;
      for (int i = 0; i < atom->nlocal; i++) {
        if (type[i] == type_list[iswaptype]) {
          if (first) qtype[iswaptype] = atom->q[i];
          first = false;
          if (qtype[iswaptype] != atom->q[i])
            error->all(FLERR,"All atoms of a swapped type must have the same charge.");
        }
      }
    }
  }
  
  memory->create(sqrt_mass_ratio,atom->ntypes+1,atom->ntypes+1,"atom/swap:sqrt_mass_ratio");
  for (int itype = 1; itype <= atom->ntypes; itype++)
    for (int jtype = 1; jtype <= atom->ntypes; jtype++)
      sqrt_mass_ratio[itype][jtype] = sqrt(atom->mass[itype]/atom->mass[jtype]);
  
  // check to see if itype and jtype cutoffs are the same
  // if not, reneighboring will be needed between swaps
  
  double **cutsq = force->pair->cutsq;
  unequal_cutoffs = false;
  for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++)
    for (int jswaptype = 0; jswaptype < nswaptypes; jswaptype++)
      for (int ktype = 1; ktype <= atom->ntypes; ktype++)
        if (cutsq[type_list[iswaptype]][ktype] != cutsq[type_list[jswaptype]][ktype])
          unequal_cutoffs = true;
  
  // check that no swappable atoms are in atom->firstgroup
  // swapping such an atom might not leave firstgroup atoms first
  
  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if ((mask[i] == groupbit) && (mask[i] && firstgroupbit)) flag = 1;
    
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

    if (flagall)
      error->all(FLERR,"Cannot do atom/swap on atoms in atom_modify first group");
  }
}

/* ----------------------------------------------------------------------
   attempt Monte Carlo swaps
------------------------------------------------------------------------- */

void FixAtomSwap::pre_exchange()
{
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;
  
  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build();
  
  energy_stored = energy_full();
  
  int nsuccess = 0;
  if (semi_grand_flag) {
    update_semi_grand_atoms_list();
    for (int i = 0; i < ncycles; i++) nsuccess += attempt_semi_grand();
  } else {
    update_swap_atoms_list();
    for (int i = 0; i < ncycles; i++) nsuccess += attempt_swap();
  }
  
  nswap_attempts += ncycles;
  nswap_successes += nsuccess;

  energy_full();
  next_reneighbor = update->ntimestep + nevery;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixAtomSwap::attempt_semi_grand()
{
  if (nswap == 0) return 0;

  double energy_before = energy_stored;

  int itype,jtype,jswaptype;
  double qtmp;

  int i = pick_semi_grand_atom();
  if (i >= 0) {
    jswaptype = static_cast<int> (nswaptypes*random_equal->uniform());
    jtype = type_list[jswaptype];
    itype = atom->type[i];
    while (itype == jtype) {
      jswaptype = static_cast<int> (nswaptypes*random_equal->uniform());
      jtype = type_list[jswaptype];
    }
    atom->type[i] = jtype;
    if (atom->q_flag) { 
      qtmp = atom->q[i];
      atom->q[i] = qtype[jswaptype];
    }
  }
  
  if (unequal_cutoffs) {
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    if (modify->n_pre_neighbor) modify->pre_neighbor();
    neighbor->build();
  } else {
    comm->forward_comm_fix(this);
  }
  
  double energy_after = energy_full();

  if (random_equal->uniform() < 
      exp(-beta*(energy_after - energy_before + 
              delta_mu[jtype] - delta_mu[itype]))) {
    update_semi_grand_atoms_list();
    energy_stored = energy_after;
    if (conserve_ke_flag) {
      if (i >= 0) {
        atom->v[i][0] *= sqrt_mass_ratio[itype][jtype];
        atom->v[i][1] *= sqrt_mass_ratio[itype][jtype];
        atom->v[i][2] *= sqrt_mass_ratio[itype][jtype];
      }
    }
    return 1;
  } else {
    if (i >= 0) {
      atom->type[i] = itype;
      if (atom->q_flag) atom->q[i] = qtmp;
    }
    energy_stored = energy_before;
    
    if (unequal_cutoffs) {
      if (domain->triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      comm->exchange();
      comm->borders();
      if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
      if (modify->n_pre_neighbor) modify->pre_neighbor();
      neighbor->build();
    } else {
      comm->forward_comm_fix(this);
    }
  } 
  return 0;
}


/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixAtomSwap::attempt_swap()
{
  if ((niswap == 0) || (njswap == 0)) return 0;
  
  double energy_before = energy_stored;

  int i = pick_i_swap_atom();
  int j = pick_j_swap_atom();
  int itype = type_list[0];
  int jtype = type_list[1];
  
  if (i >= 0) {
    atom->type[i] = jtype;
    if (atom->q_flag) atom->q[i] = qtype[1];
  }
  if (j >= 0) {
    atom->type[j] = itype;
    if (atom->q_flag) atom->q[j] = qtype[0];
  }
  
  if (unequal_cutoffs) {
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    if (modify->n_pre_neighbor) modify->pre_neighbor();
    neighbor->build();
  } else {
    comm->forward_comm_fix(this);
  }
  
  double energy_after = energy_full();

  if (random_equal->uniform() < 
      exp(beta*(energy_before - energy_after))) {
    update_swap_atoms_list();
    energy_stored = energy_after;
    if (conserve_ke_flag) {
      if (i >= 0) {
        atom->v[i][0] *= sqrt_mass_ratio[itype][jtype];
        atom->v[i][1] *= sqrt_mass_ratio[itype][jtype];
        atom->v[i][2] *= sqrt_mass_ratio[itype][jtype];
      }
      if (j >= 0) {
        atom->v[j][0] *= sqrt_mass_ratio[jtype][itype];
        atom->v[j][1] *= sqrt_mass_ratio[jtype][itype];
        atom->v[j][2] *= sqrt_mass_ratio[jtype][itype];
      }
    }
    return 1;
  } else {
    if (i >= 0) {
      atom->type[i] =  type_list[0];
      if (atom->q_flag) atom->q[i] = qtype[0];
    }
    if (j >= 0) {
      atom->type[j] =  type_list[1];
      if (atom->q_flag) atom->q[j] = qtype[1];
    }
    energy_stored = energy_before;
    
    if (unequal_cutoffs) {
      if (domain->triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      comm->exchange();
      comm->borders();
      if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
      if (modify->n_pre_neighbor) modify->pre_neighbor();
      neighbor->build();
    } else {
      comm->forward_comm_fix(this);
    }
  }
  return 0;
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixAtomSwap::energy_full()
{ 
  int eflag = 1;
  int vflag = 0;

  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) force->kspace->compute(eflag,vflag);
 
  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixAtomSwap::pick_semi_grand_atom()
{
  int i = -1;
  int iwhichglobal = static_cast<int> (nswap*random_equal->uniform());
  if ((iwhichglobal >= nswap_before) &&
      (iwhichglobal < nswap_before + nswap_local)) {
    int iwhichlocal = iwhichglobal - nswap_before;
    i = local_swap_atom_list[iwhichlocal];
  }

  return i;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixAtomSwap::pick_i_swap_atom()
{
  int i = -1;
  int iwhichglobal = static_cast<int> (niswap*random_equal->uniform());
  if ((iwhichglobal >= niswap_before) &&
      (iwhichglobal < niswap_before + niswap_local)) {
    int iwhichlocal = iwhichglobal - niswap_before;
    i = local_swap_iatom_list[iwhichlocal];
  }

  return i;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixAtomSwap::pick_j_swap_atom()
{
  int j = -1;
  int jwhichglobal = static_cast<int> (njswap*random_equal->uniform());
  if ((jwhichglobal >= njswap_before) &&
      (jwhichglobal < njswap_before + njswap_local)) {
    int jwhichlocal = jwhichglobal - njswap_before;
    j = local_swap_jatom_list[jwhichlocal];
  }

  return j;
}

/* ----------------------------------------------------------------------
   update the list of gas atoms
------------------------------------------------------------------------- */

void FixAtomSwap::update_semi_grand_atoms_list()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **x = atom->x;
      
  if (nlocal > atom_swap_nmax) {
    memory->sfree(local_swap_atom_list);
    atom_swap_nmax = atom->nmax;
    local_swap_atom_list = (int *) memory->smalloc(atom_swap_nmax*sizeof(int),
     "MCSWAP:local_swap_atom_list");
  }

  nswap_local = 0;
  
  if (regionflag) {
  
    for (int i = 0; i < nlocal; i++) {
      if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]) == 1) {
        if (atom->mask[i] & groupbit) {
          local_swap_atom_list[nswap_local] = i;
          nswap_local++;
        }
      }
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      if (atom->mask[i] & groupbit) {
        local_swap_atom_list[nswap_local] = i;
        nswap_local++;
      }
    }
  }

  MPI_Allreduce(&nswap_local,&nswap,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&nswap_local,&nswap_before,1,MPI_INT,MPI_SUM,world);
  nswap_before -= nswap_local;
}


/* ----------------------------------------------------------------------
   update the list of gas atoms
------------------------------------------------------------------------- */

void FixAtomSwap::update_swap_atoms_list()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **x = atom->x;
      
  if (nlocal > atom_swap_nmax) {
    memory->sfree(local_swap_iatom_list);
    memory->sfree(local_swap_jatom_list);
    atom_swap_nmax = atom->nmax;
    local_swap_iatom_list = (int *) memory->smalloc(atom_swap_nmax*sizeof(int),
     "MCSWAP:local_swap_iatom_list");
    local_swap_jatom_list = (int *) memory->smalloc(atom_swap_nmax*sizeof(int),
     "MCSWAP:local_swap_jatom_list");
  }

  niswap_local = 0;
  njswap_local = 0;
  
  if (regionflag) {
  
    for (int i = 0; i < nlocal; i++) {
      if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]) == 1) {
        if (atom->mask[i] & groupbit) {
          if (type[i] ==  type_list[0]) {
            local_swap_iatom_list[niswap_local] = i;
            niswap_local++;
          } else if (type[i] ==  type_list[1]) {
            local_swap_jatom_list[njswap_local] = i;
            njswap_local++;
          }
        }
      }
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      if (atom->mask[i] & groupbit) {
        if (type[i] ==  type_list[0]) {
          local_swap_iatom_list[niswap_local] = i;
          niswap_local++;
        } else if (type[i] ==  type_list[1]) {
          local_swap_jatom_list[njswap_local] = i;
          njswap_local++;
        }
      }
    }
  }

  MPI_Allreduce(&niswap_local,&niswap,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&niswap_local,&niswap_before,1,MPI_INT,MPI_SUM,world);
  niswap_before -= niswap_local;
  
  MPI_Allreduce(&njswap_local,&njswap,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&njswap_local,&njswap_before,1,MPI_INT,MPI_SUM,world);
  njswap_before -= njswap_local;
}

/* ---------------------------------------------------------------------- */

int FixAtomSwap::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  int *type = atom->type;
  double *q = atom->q;
  
  m = 0;

  if (atom->q_flag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
      buf[m++] = q[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixAtomSwap::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  int *type = atom->type;
  double *q = atom->q;
  
  m = 0;
  last = first + n;

  if (atom->q_flag) {
    for (i = first; i < last; i++) {
      type[i] = static_cast<int> (buf[m++]);
      q[i] = buf[m++];
    }
  } else {
    for (i = first; i < last; i++)
      type[i] = static_cast<int> (buf[m++]);
  }
}

/* ----------------------------------------------------------------------
  return acceptance ratio
------------------------------------------------------------------------- */

double FixAtomSwap::compute_vector(int n)
{
  if (n == 0) return nswap_attempts;
  if (n == 1) return nswap_successes;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixAtomSwap::memory_usage()
{
  double bytes = atom_swap_nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixAtomSwap::write_restart(FILE *fp)
{
  int n = 0;
  double list[4];
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

void FixAtomSwap::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  random_equal->reset(seed);

  next_reneighbor = static_cast<int> (list[n++]);
}
