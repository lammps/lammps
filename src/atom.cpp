/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "output.h"
#include "thermo.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "memory.h"
#include "error.h"

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define DELTA 10000
#define DELTA_CALLBACK 1

/* ---------------------------------------------------------------------- */

Atom::Atom(int narg, char **arg)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);
  
  // arg[0] sets one style flag to 1
  
  style_angle = style_atomic = style_bond = style_charge = style_dipole =
    style_dpd = style_full = style_granular = style_molecular = 
    style_peri = style_hybrid = 0;
  
  set_style(arg[0]);

  // if hybrid style, set additional style for each additional arg
  
  if (style_hybrid) {
    if (narg < 2) error->all("Illegal atom style hybrid command");
    for (int i = 1; i < narg; i++) {
      if (strcmp(arg[i],"hybrid") == 0) 
	error->all("Atom style hybrid cannot have hybrid as an argument");
      set_style(arg[i]);
    }
  }

  // set low-level flags from style flags

  mass_require = 0;
  if (style_angle || style_atomic || style_bond || style_charge ||
      style_dipole || style_dpd || style_full ||
      style_molecular || style_hybrid) mass_require = 1;

  mass_allow = 0;
  if (style_granular || style_peri) mass_allow = 1;

  charge_allow = 0;
  if (style_charge || style_full || style_dipole) charge_allow = 1;

  dipole_require = 0;
  if (style_dipole) dipole_require = 1;

  molecular = 0;
  if (style_angle || style_bond || style_full || style_molecular)
    molecular = 1;

  bonds_allow = 0;
  if (style_angle || style_bond || style_full || style_molecular)
    bonds_allow = 1;

  angles_allow = 0;
  if (style_angle || style_full || style_molecular) angles_allow = 1;

  dihedrals_allow = 0;
  if (style_full || style_molecular) dihedrals_allow = 1;

  impropers_allow = 0;
  if (style_full || style_molecular) impropers_allow = 1;

  // set size variables from styles
  // if size_comm and size_reverse are changed,
  //   must also change direct_flag in comm::init()

  size_comm = 3;
  if (style_dpd) size_comm += 3;                  // v
  if (style_dipole) size_comm += 3;               // mu
  if (style_granular) size_comm += 6;             // v,phiv

  size_reverse = 3;
  if (style_dipole) size_reverse += 3;            // torque
  if (style_granular) size_reverse += 3;          // phia

  size_border = 6;
  if (charge_allow) size_border += 1;             // q
  if (style_dpd) size_border += 3;                // v
  if (style_dipole) size_border += 3;             // mu
  if (style_granular) size_border += 8;           // v,phiv,radius,rmass
  if (molecular) size_border += 1;                // molecule

  size_atom_valid = 5;
  if (molecular) size_atom_valid += 1;            // molecule
  if (charge_allow) size_atom_valid += 1;         // q
  if (style_granular) size_atom_valid += 2;       // radius,density
  if (style_peri) size_atom_valid += 2;           // vfrac,rmass
  if (style_dipole) size_atom_valid += 3;         // mu

  // initialize atom arrays to empty

  tag = type = mask = image = NULL;
  x = v = f = NULL;

  q = NULL;
  mu = omega = torque = NULL;
  phix = phiv = phia = NULL;
  radius = density = rmass = vfrac = NULL;

  molecule = NULL;

  maxspecial = 1;
  nspecial = NULL;
  special = NULL;

  num_bond = NULL;
  bond_type = bond_atom = NULL;

  num_angle = NULL;
  angle_type = angle_atom1 = angle_atom2 = angle_atom3 = NULL;

  num_dihedral = NULL;
  dihedral_type = dihedral_atom1 = dihedral_atom2 = NULL;
  dihedral_atom3 = dihedral_atom4 = NULL;

  num_improper = NULL;
  improper_type = improper_atom1 = improper_atom2 = NULL;
  improper_atom3 = improper_atom4 = NULL;

  // set data and restart file header values to defaults

  natoms = 0;
  nbonds = nangles = ndihedrals = nimpropers = 0;
  ntypes = 0;
  nbondtypes = nangletypes = ndihedraltypes = nimpropertypes = 0;
  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;

  // no atoms initially

  nlocal = 0;
  nghost = 0;
  nmax = 0;

  // ntype length arrays

  mass = NULL;
  mass_setflag = NULL;
  dipole = NULL;
  dipole_setflag = NULL;

  // callback lists & extra restart info

  nextra_grow = nextra_restart = 0;
  extra_grow = extra_restart = NULL;
  nextra_grow_max = nextra_restart_max = 0;
  nextra_store = 0;
  extra = NULL;

  // set default mapping values and hash table primes

  tag_enable = 1;
  if (molecular) map_style = 1;
  else map_style = 0;
  map_tag_max = 0;
  map_nhash = 0;
  
  nprimes = 38;
  primes = new int[nprimes];
  int plist[] = {5041,10007,20011,30011,40009,50021,60013,70001,80021,
		 90001,100003,110017,120011,130003,140009,150001,160001,
		 170003,180001,190027,200003,210011,220009,230003,240007,
		 250007,260003,270001,280001,290011,300007,310019,320009,
		 330017,340007,350003,362881,3628801};
  for (int i = 0; i < nprimes; i++) primes[i] = plist[i];

  // constant

  PI = 4.0*atan(1.0);
}

/* ---------------------------------------------------------------------- */

Atom::~Atom()
{
  // delete atom arrays

  memory->sfree(tag);
  memory->sfree(type);
  memory->sfree(mask);
  memory->sfree(image);
  memory->destroy_2d_double_array(x);
  memory->destroy_2d_double_array(v);
  memory->destroy_2d_double_array(f);

  memory->sfree(q);
  memory->destroy_2d_double_array(mu);
  memory->destroy_2d_double_array(omega);
  memory->destroy_2d_double_array(torque);
  memory->destroy_2d_double_array(phix);
  memory->destroy_2d_double_array(phiv);
  memory->destroy_2d_double_array(phia);
  memory->sfree(radius);
  memory->sfree(density);
  memory->sfree(rmass);
  memory->sfree(vfrac);

  memory->sfree(molecule);

  memory->destroy_2d_int_array(nspecial);
  memory->destroy_2d_int_array(special);

  memory->sfree(num_bond);
  memory->destroy_2d_int_array(bond_type);
  memory->destroy_2d_int_array(bond_atom);

  memory->sfree(num_angle);
  memory->destroy_2d_int_array(angle_type);
  memory->destroy_2d_int_array(angle_atom1);
  memory->destroy_2d_int_array(angle_atom2);
  memory->destroy_2d_int_array(angle_atom3);

  memory->sfree(num_dihedral);
  memory->destroy_2d_int_array(dihedral_type);
  memory->destroy_2d_int_array(dihedral_atom1);
  memory->destroy_2d_int_array(dihedral_atom2);
  memory->destroy_2d_int_array(dihedral_atom3);
  memory->destroy_2d_int_array(dihedral_atom4);

  memory->sfree(num_improper);
  memory->destroy_2d_int_array(improper_type);
  memory->destroy_2d_int_array(improper_atom1);
  memory->destroy_2d_int_array(improper_atom2);
  memory->destroy_2d_int_array(improper_atom3);
  memory->destroy_2d_int_array(improper_atom4);

  // delete auxiliary arrays

  delete [] style;

  delete [] mass;
  delete [] mass_setflag;
  delete [] dipole;
  delete [] dipole_setflag;

  memory->sfree(extra_grow);
  memory->sfree(extra_restart);
  memory->destroy_2d_double_array(extra);

  map_delete();
  delete [] primes;
}

/* ---------------------------------------------------------------------- */

void Atom::set_style(char *name)
{
  if (strcmp(name,"angle") == 0) style_angle = 1;
  else if (strcmp(name,"atomic") == 0) style_atomic = 1;
  else if (strcmp(name,"bond") == 0) style_bond = 1;
  else if (strcmp(name,"charge") == 0) style_charge = 1;
  else if (strcmp(name,"dipole") == 0) style_dipole = 1;
  else if (strcmp(name,"dpd") == 0) style_dpd = 1;
  else if (strcmp(name,"full") == 0) style_full = 1;
  else if (strcmp(name,"granular") == 0) style_granular = 1;
  else if (strcmp(name,"molecular") == 0) style_molecular = 1;
  else if (strcmp(name,"peri") == 0) style_peri = 1;
  else if (strcmp(name,"hybrid") == 0) style_hybrid = 1;
  else error->all("Illegal atom_style command");
}

/* ----------------------------------------------------------------------
   return 1 if named atom style is set as pure style or sub-class of hybrid
   else return 0
------------------------------------------------------------------------- */

int Atom::check_style(char *name)
{
  if (!strcmp(name,"angle") && style_angle) return 1;
  else if (!strcmp(name,"atomic") && style_atomic) return 1;
  else if (!strcmp(name,"bond") && style_bond) return 1;
  else if (!strcmp(name,"charge") && style_charge) return 1;
  else if (!strcmp(name,"dipole") && style_dipole) return 1;
  else if (!strcmp(name,"dpd") && style_dpd) return 1;
  else if (!strcmp(name,"full") && style_full) return 1;
  else if (!strcmp(name,"granular") && style_granular) return 1;
  else if (!strcmp(name,"molecular") && style_molecular) return 1;
  else if (!strcmp(name,"peri") && style_peri) return 1;
  else if (!strcmp(name,"hybrid") && style_hybrid) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   convert style settings to list of keywords
   put "hybrid" first if it is set
   return list and number of words in list
------------------------------------------------------------------------- */

int Atom::style2arg(char **&name)
{
  int n = style_angle + style_atomic + style_bond + style_charge +
    style_dipole + style_dpd + style_full + 
    style_granular + style_molecular + style_peri + style_hybrid;

  name = new char*[n];

  n = 0;
  if (style_hybrid) name[n++] = style2word("hybrid");
  if (style_angle) name[n++] = style2word("angle");
  if (style_atomic) name[n++] = style2word("atomic");
  if (style_bond) name[n++] = style2word("bond");
  if (style_charge) name[n++] = style2word("charge");
  if (style_dipole) name[n++] = style2word("dipole");
  if (style_dpd) name[n++] = style2word("dpd");
  if (style_full) name[n++] = style2word("full");
  if (style_granular) name[n++] = style2word("granular");
  if (style_molecular) name[n++] = style2word("molecular");
  if (style_peri) name[n++] = style2word("peri");

  return n;
}

/* ----------------------------------------------------------------------
   copy modify settings from old Atom class to current Atom class
------------------------------------------------------------------------- */

void Atom::settings(Atom *old)
{
  map_style = old->map_style;
}

/* ----------------------------------------------------------------------
   copy name into new word string
------------------------------------------------------------------------- */

char *Atom::style2word(char *name)
{
  int n = strlen(name) + 1;
  char *word = new char[n];
  strcpy(word,name);
  return word;
}

/* ---------------------------------------------------------------------- */

void Atom::init()
{
  // delete extra array since it doesn't persist past first run

  if (nextra_store) {
    memory->destroy_2d_double_array(extra);
    extra = NULL;
    nextra_store = 0;
  }

  // check arrays that are atom type in length

  check_mass();
  check_dipole();

  // for dipole style:
  // for dipole type atoms, check that dipole moment is set, normalize it
  // for non-dipole type atoms, check that dipole moment is 0.0

  if (style_dipole) {
    double msq,scale;
    int flag = 0;
    for (int i = 0; i < nlocal; i++) {
      msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2];
      if (dipole[type[i]] > 0.0 && msq == 0.0) flag++;
      else if (dipole[type[i]] > 0.0) {
	scale = dipole[type[i]]/sqrt(msq);
	mu[i][0] *= scale;
	mu[i][1] *= scale;
	mu[i][2] *= scale;
      }
      else if (dipole[type[i]] == 0.0 && msq > 0.0) flag++;
    }
    int flag_all;
    MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
    if (flag) error->all("Inconsistent dipole settings for some atoms");
  }

  // for granular style:
  // insure LJ units and thermo style granular or custom is used

  if (style_granular) {
    if (strcmp(update->unit_style,"lj") != 0)
      error->all("Must use lj units with atom style granular");
    if ((strcmp(output->thermo->style,"granular") != 0) &&
	(strcmp(output->thermo->style,"custom") != 0))
      error->all("Must use atom style granular with granular thermo output");
  }

  // don't allow granular and dpd styles together
  // only reason is that they both communicate v in pack_comm
  // ok if use better logic for setting size_comm and atom_hybrid::pack_comm

  if (style_granular && style_dpd)
    error->all("Atom style granular and dpd cannot be used together");
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n 
------------------------------------------------------------------------- */

void Atom::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;

  tag = (int *) memory->srealloc(tag,nmax*sizeof(int),"atom:tag");
  type = (int *) memory->srealloc(type,nmax*sizeof(int),"atom:type");
  mask = (int *) memory->srealloc(mask,nmax*sizeof(int),"atom:mask");
  image = (int *) memory->srealloc(image,nmax*sizeof(int),"atom:image");
  x = memory->grow_2d_double_array(x,nmax,3,"atom:x");
  v = memory->grow_2d_double_array(v,nmax,3,"atom:v");
  f = memory->grow_2d_double_array(f,nmax,3,"atom:f");

  if (charge_allow)
    q = (double *) memory->srealloc(q,nmax*sizeof(double),"atom:q");
  if (mass_allow)
    rmass = (double *) memory->srealloc(rmass,nmax*sizeof(double),
					"atom:rmass");

  if (style_dipole) {
    mu = memory->grow_2d_double_array(mu,nmax,3,"atom:mu");
    omega = memory->grow_2d_double_array(omega,nmax,3,"atom:omega");
    torque = memory->grow_2d_double_array(torque,nmax,3,"atom:torque");
  }

  if (style_granular) {
    phix = memory->grow_2d_double_array(phix,nmax,3,"atom:phix");
    phiv = memory->grow_2d_double_array(phiv,nmax,3,"atom:phiv");
    phia = memory->grow_2d_double_array(phia,nmax,3,"atom:phia");
    radius = (double *) memory->srealloc(radius,nmax*sizeof(double),
					 "atom:radius");
    density = (double *) memory->srealloc(density,nmax*sizeof(double),
					  "atom:density");
  }

  if (style_peri)
    vfrac = (double *) memory->srealloc(vfrac,nmax*sizeof(double),
					"atom:vfrac");

  if (molecular) {
    molecule = (int *) 
      memory->srealloc(molecule,nmax*sizeof(int),"atom:molecule");

    nspecial = memory->grow_2d_int_array(nspecial,nmax,3,"atom:nspecial");
    special = 
      memory->grow_2d_int_array(special,nmax,maxspecial,"atom:special");

    if (bonds_allow) {
      num_bond = (int *) 
	memory->srealloc(num_bond,nmax*sizeof(int),"atom:num_bond");
      bond_type = memory->grow_2d_int_array(bond_type,nmax,bond_per_atom,
					    "atom:bond_type");
      bond_atom = memory->grow_2d_int_array(bond_atom,nmax,bond_per_atom,
					    "atom:bond_atom");
    }

    if (angles_allow) {
      num_angle = (int *) 
	memory->srealloc(num_angle,nmax*sizeof(int),"atom:num_angle");
      angle_type = memory->grow_2d_int_array(angle_type,nmax,angle_per_atom,
					     "atom:angle_type");
      angle_atom1 = memory->grow_2d_int_array(angle_atom1,nmax,angle_per_atom,
					      "atom:angle_atom1");
      angle_atom2 = memory->grow_2d_int_array(angle_atom2,nmax,angle_per_atom,
					      "atom:angle_atom2");
      angle_atom3 = memory->grow_2d_int_array(angle_atom3,nmax,angle_per_atom,
					      "atom:angle_atom3");
    }

    if (dihedrals_allow) {
      num_dihedral = (int *) 
	memory->srealloc(num_dihedral,nmax*sizeof(int),"atom:num_dihedral");
      dihedral_type = 
	memory->grow_2d_int_array(dihedral_type,nmax,dihedral_per_atom,
				  "atom:dihedral_type");
      dihedral_atom1 = 
	memory->grow_2d_int_array(dihedral_atom1,nmax,dihedral_per_atom,
				  "atom:dihedral_atom1");
      dihedral_atom2 = 
	memory->grow_2d_int_array(dihedral_atom2,nmax,dihedral_per_atom,
				  "atom:dihedral_atom2");
      dihedral_atom3 = 
	memory->grow_2d_int_array(dihedral_atom3,nmax,dihedral_per_atom,
				  "atom:dihedral_atom3");
      dihedral_atom4 = 
	memory->grow_2d_int_array(dihedral_atom4,nmax,dihedral_per_atom,
				  "atom:dihedral_atom4");
    }
    
    if (impropers_allow) {
      num_improper = (int *) 
	memory->srealloc(num_improper,nmax*sizeof(int),"atom:num_improper");
      improper_type = 
	memory->grow_2d_int_array(improper_type,nmax,improper_per_atom,
				  "atom:improper_type");
      improper_atom1 = 
	memory->grow_2d_int_array(improper_atom1,nmax,improper_per_atom,
				  "atom:improper_atom1");
      improper_atom2 = 
	memory->grow_2d_int_array(improper_atom2,nmax,improper_per_atom,
				  "atom:improper_atom2");
      improper_atom3 = 
	memory->grow_2d_int_array(improper_atom3,nmax,improper_per_atom,
				  "atom:improper_atom3");
      improper_atom4 = 
	memory->grow_2d_int_array(improper_atom4,nmax,improper_per_atom,
				  "atom:improper_atom4");
    }
  }

  if (nextra_grow)
    for (int iextra = 0; iextra < nextra_grow; iextra++) 
      modify->fix[extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   modify parameters of the atom style
------------------------------------------------------------------------- */

void Atom::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal atom_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"map") == 0) {
      if (iarg+2 > narg) error->all("Illegal atom_modify command");
      if (strcmp(arg[iarg+1],"array") == 0) map_style = 1;
      else if (strcmp(arg[iarg+1],"hash") == 0) map_style = 2;
      else error->all("Illegal atom_modify command");
      iarg += 2;
    } else error->all("Illegal atom_modify command");
  }
}

/* ----------------------------------------------------------------------
   allocate and initialize array or hash table for global -> local map
   set map_tag_max = largest atom ID (may be larger than natoms)
   for array option:
     array length = 1 to largest tag of any atom
     set entire array to -1 as initial values
   for hash option:
     map_nhash = length of hash table
     map_nbucket = # of hash buckets, prime larger than map_nhash
       so buckets will only be filled with 0 or 1 atoms on average
------------------------------------------------------------------------- */

void Atom::map_init()
{
  map_delete();

  if (tag_enable == 0)
    error->all("Cannot create an atom map unless atoms have IDs");

  int max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&map_tag_max,1,MPI_INT,MPI_MAX,world);

  if (map_style == 1) {
    map_array = (int *) 
      memory->smalloc((map_tag_max+1)*sizeof(int),"atom:map_array");
    for (int i = 0; i <= map_tag_max; i++) map_array[i] = -1;

  } else {

    // map_nhash = max of atoms/proc or total atoms, times 2, at least 1000

    int nper = static_cast<int> (natoms/comm->nprocs);
    map_nhash = MAX(nper,nmax);
    if (map_nhash > natoms) map_nhash = static_cast<int> (natoms);
    if (comm->nprocs > 1) map_nhash *= 2;
    map_nhash = MAX(map_nhash,1000);

    // map_nbucket = prime just larger than map_nhash

    int n = map_nhash/10000;
    n = MIN(n,nprimes-1);
    map_nbucket = primes[n];
    if (map_nbucket < map_nhash && n < nprimes-1) map_nbucket = primes[n+1];

    // set all buckets to empty
    // set hash to map_nhash in length
    // put all hash entries in free list and point them to each other

    map_bucket = new int[map_nbucket];
    for (int i = 0; i < map_nbucket; i++) map_bucket[i] = -1;

    map_hash = new HashElem[map_nhash];
    map_nused = 0;
    map_free = 0;
    for (int i = 0; i < map_nhash; i++) map_hash[i].next = i+1;
    map_hash[map_nhash-1].next = -1;
  }
}

/* ----------------------------------------------------------------------
   clear global -> local map for all of my own and ghost atoms
   for hash table option:
     global ID may not be in table if image atom was already cleared
------------------------------------------------------------------------- */

void Atom::map_clear()
{
  if (map_style == 1) {
    int nall = nlocal + nghost;
    for (int i = 0; i < nall; i++) map_array[tag[i]] = -1;

  } else {
    int previous,global,ibucket,index;
    int nall = nlocal + nghost;
    for (int i = 0; i < nall; i++) {

      // search for key
      // if don't find it, done

      previous = -1;
      global = tag[i];
      ibucket = global % map_nbucket;
      index = map_bucket[ibucket];
      while (index > -1) {
	if (map_hash[index].global == global) break;
	previous = index;
	index = map_hash[index].next;
      }
      if (index == -1) continue;

      // delete the hash entry and add it to free list
      // special logic if entry is 1st in the bucket
    
      if (previous == -1) map_bucket[ibucket] = map_hash[index].next;
      else map_hash[previous].next = map_hash[index].next;
    
      map_hash[index].next = map_free;
      map_free = index;
      map_nused--;
    }
  }
}

/* ----------------------------------------------------------------------
   set global -> local map for all of my own and ghost atoms
   loop in reverse order so that nearby images take precedence over far ones
     and owned atoms take precedence over images
   this enables valid lookups of bond topology atoms 
   for hash table option:
     if hash table too small, re-init
     global ID may already be in table if image atom was set
------------------------------------------------------------------------- */

void Atom::map_set()
{
  if (map_style == 1) {
    int nall = nlocal + nghost;
    for (int i = nall-1; i >= 0 ; i--) map_array[tag[i]] = i;

  } else {
    int previous,global,ibucket,index;
    int nall = nlocal + nghost;
    if (nall > map_nhash) map_init();

    for (int i = nall-1; i >= 0 ; i--) {

      // search for key
      // if found it, just overwrite local value with index

      previous = -1;
      global = tag[i];
      ibucket = global % map_nbucket;
      index = map_bucket[ibucket];
      while (index > -1) {
	if (map_hash[index].global == global) break;
	previous = index;
	index = map_hash[index].next;
      }
      if (index > -1) {
	map_hash[index].local = i;
	continue;
      }

      // take one entry from free list
      // add the new global/local pair as entry at end of bucket list
      // special logic if this entry is 1st in bucket

      index = map_free;
      map_free = map_hash[map_free].next;
      if (previous == -1) map_bucket[ibucket] = index;
      else map_hash[previous].next = index;
      map_hash[index].global = global;
      map_hash[index].local = i;
      map_hash[index].next = -1;
      map_nused++;
    }
  }
}

/* ----------------------------------------------------------------------
   set global to local map for one atom
   for hash table option:
     global ID may already be in table if atom was already set
------------------------------------------------------------------------- */

void Atom::map_one(int global, int local)
{
  if (map_style == 1) map_array[global] = local;

  else {
    // search for key
    // if found it, just overwrite local value with index
    
    int previous = -1;
    int ibucket = global % map_nbucket;
    int index = map_bucket[ibucket];
    while (index > -1) {
      if (map_hash[index].global == global) break;
      previous = index;
      index = map_hash[index].next;
    }
    if (index > -1) {
      map_hash[index].local = local;
      return;
    }
    
    // take one entry from free list
    // add the new global/local pair as entry at end of bucket list
    // special logic if this entry is 1st in bucket
    
    index = map_free;
    map_free = map_hash[map_free].next;
    if (previous == -1) map_bucket[ibucket] = index;
    else map_hash[previous].next = index;
    map_hash[index].global = global;
    map_hash[index].local = local;
    map_hash[index].next = -1;
    map_nused++;
  }
}

/* ----------------------------------------------------------------------
   free the array or hash table for global to local mapping
------------------------------------------------------------------------- */

void Atom::map_delete()
{
  if (map_style == 1) {
    if (map_tag_max) memory->sfree(map_array);
  } else {
    if (map_nhash) {
      delete [] map_bucket;
      delete [] map_hash;
    }
    map_nhash = 0;
  }
  map_tag_max = 0;
}

/* ----------------------------------------------------------------------
   lookup global ID in hash table, return local index
------------------------------------------------------------------------- */

int Atom::map_find_hash(int global)
{
  int local = -1;
  int index = map_bucket[global % map_nbucket];
  while (index > -1) {
    if (map_hash[index].global == global) {
      local = map_hash[index].local;
      break;
    }
    index = map_hash[index].next;
  }
  return local;
}

/* ----------------------------------------------------------------------
   add unique tags to any atoms with tag = 0
   new tags are grouped by proc and start after max current tag
   called after creating new atoms 
------------------------------------------------------------------------- */

void Atom::tag_extend()
{
  // maxtag_all = max tag for all atoms

  int maxtag = 0;
  for (int i = 0; i < nlocal; i++) maxtag = MAX(maxtag,tag[i]);
  int maxtag_all;
  MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_INT,MPI_MAX,world);

  // notag = # of atoms I own with no tag (tag = 0)
  // notag_sum = # of total atoms on procs <= me with no tag

  int notag = 0;
  for (int i = 0; i < nlocal; i++) if (tag[i] == 0) notag++;
  int notag_sum;
  MPI_Scan(&notag,&notag_sum,1,MPI_INT,MPI_SUM,world);

  // itag = 1st new tag that my untagged atoms should use

  int itag = maxtag_all + notag_sum - notag + 1;
  for (int i = 0; i < nlocal; i++) if (tag[i] == 0) tag[i] = itag++;
}

/* ----------------------------------------------------------------------
   check if atom tags are consecutive from 1 to Natoms
   return 0 if any tag <= 0 or maxtag > Natoms
   return 1 if OK (doesn't actually check if all tag values are used)
------------------------------------------------------------------------- */

int Atom::tag_consecutive()
{
  // check[0] = flagged if any tag <= 0
  // check[1] = max tag of any atom

  int check[2],check_all[2];

  check[0] = check[1] = 0;
  for (int i = 0; i < nlocal; i++) {
    if (tag[i] <= 0) check[0] = 1;
    if (tag[i] > check[1]) check[1] = tag[i];
  }
  MPI_Allreduce(check,check_all,2,MPI_INT,MPI_MAX,world);

  if (check_all[0] || check_all[1] > natoms) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   parse one atom line from data file, check if valid for atom style 
------------------------------------------------------------------------- */

int Atom::parse_data(char *line)
{
  size_atom_actual = count_words(line);
  if (size_atom_actual == size_atom_valid || 
      size_atom_actual == size_atom_valid + 3) return 0;
  else return 1;
}

/* ----------------------------------------------------------------------
   count and return words in a single line
   make copy of line before using strtok so as not to change line
   trim anything from '#' onward
------------------------------------------------------------------------- */

int Atom::count_words(char *line)
{
  int n = strlen(line) + 1;
  char *copy = (char *) memory->smalloc(n*sizeof(char),"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->sfree(copy);
  return n;
}

/* ----------------------------------------------------------------------
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::unpack_vels(int n, char *buf)
{
  int m,tagtmp;
  double vxtmp,vytmp,vztmp;
  char *next;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    sscanf(buf,"%d %lg %lg %lg",&tagtmp,&vxtmp,&vytmp,&vztmp);
    if (tagtmp <= 0 || tagtmp > map_tag_max)
      error->one("Invalid atom ID in Velocities section of data file");
    if ((m = map(tagtmp)) >= 0) {
      v[m][0] = vxtmp;
      v[m][1] = vytmp;
      v[m][2] = vztmp;
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::unpack_bonds(int n, char *buf)
{
  int m,tmp,itype,atom1,atom2;
  char *next;
  int newton_bond = force->newton_bond;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    sscanf(buf,"%d %d %d %d",&tmp,&itype,&atom1,&atom2);
    if (atom1 <= 0 || atom1 > map_tag_max || 
	atom2 <= 0 || atom2 > map_tag_max)
      error->one("Invalid atom ID in Bonds section of data file");
    if (itype <= 0 || itype > nbondtypes)
      error->one("Invalid bond type in Bonds section of data file");
    if ((m = map(atom1)) >= 0) {
      bond_type[m][num_bond[m]] = itype;
      bond_atom[m][num_bond[m]] = atom2;
      num_bond[m]++;
    }
    if (newton_bond == 0) {
      if ((m = map(atom2)) >= 0) {
	bond_type[m][num_bond[m]] = itype;
	bond_atom[m][num_bond[m]] = atom1;
	num_bond[m]++;
      }
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::unpack_angles(int n, char *buf)
{
  int m,tmp,itype,atom1,atom2,atom3;
  char *next;
  int newton_bond = force->newton_bond;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    sscanf(buf,"%d %d %d %d %d",&tmp,&itype,&atom1,&atom2,&atom3);
    if (atom1 <= 0 || atom1 > map_tag_max || 
	atom2 <= 0 || atom2 > map_tag_max || 
	atom3 <= 0 || atom3 > map_tag_max)
      error->one("Invalid atom ID in Angles section of data file");
    if (itype <= 0 || itype > nangletypes)
      error->one("Invalid angle type in Angles section of data file");
    if ((m = map(atom2)) >= 0) {
      angle_type[m][num_angle[m]] = itype;
      angle_atom1[m][num_angle[m]] = atom1;
      angle_atom2[m][num_angle[m]] = atom2;
      angle_atom3[m][num_angle[m]] = atom3;
      num_angle[m]++;
    }
    if (newton_bond == 0) {
      if ((m = map(atom1)) >= 0) {
	angle_type[m][num_angle[m]] = itype;
	angle_atom1[m][num_angle[m]] = atom1;
	angle_atom2[m][num_angle[m]] = atom2;
	angle_atom3[m][num_angle[m]] = atom3;
	num_angle[m]++;
      }
      if ((m = map(atom3)) >= 0) {
	angle_type[m][num_angle[m]] = itype;
	angle_atom1[m][num_angle[m]] = atom1;
	angle_atom2[m][num_angle[m]] = atom2;
	angle_atom3[m][num_angle[m]] = atom3;
	num_angle[m]++;
      }
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::unpack_dihedrals(int n, char *buf)
{
  int m,tmp,itype,atom1,atom2,atom3,atom4;
  char *next;
  int newton_bond = force->newton_bond;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    sscanf(buf,"%d %d %d %d %d %d",&tmp,&itype,&atom1,&atom2,&atom3,&atom4);
    if (atom1 <= 0 || atom1 > map_tag_max || 
	atom2 <= 0 || atom2 > map_tag_max || 
	atom3 <= 0 || atom3 > map_tag_max || 
	atom4 <= 0 || atom4 > map_tag_max)
      error->one("Invalid atom ID in Dihedrals section of data file");
    if (itype <= 0 || itype > ndihedraltypes)
      error->one("Invalid dihedral type in Dihedrals section of data file");
    if ((m = map(atom2)) >= 0) {
      dihedral_type[m][num_dihedral[m]] = itype;
      dihedral_atom1[m][num_dihedral[m]] = atom1;
      dihedral_atom2[m][num_dihedral[m]] = atom2;
      dihedral_atom3[m][num_dihedral[m]] = atom3;
      dihedral_atom4[m][num_dihedral[m]] = atom4;
      num_dihedral[m]++;
    }
    if (newton_bond == 0) {
      if ((m = map(atom1)) >= 0) {
	dihedral_type[m][num_dihedral[m]] = itype;
	dihedral_atom1[m][num_dihedral[m]] = atom1;
	dihedral_atom2[m][num_dihedral[m]] = atom2;
	dihedral_atom3[m][num_dihedral[m]] = atom3;
	dihedral_atom4[m][num_dihedral[m]] = atom4;
	num_dihedral[m]++;
      }
      if ((m = map(atom3)) >= 0) {
	dihedral_type[m][num_dihedral[m]] = itype;
	dihedral_atom1[m][num_dihedral[m]] = atom1;
	dihedral_atom2[m][num_dihedral[m]] = atom2;
	dihedral_atom3[m][num_dihedral[m]] = atom3;
	dihedral_atom4[m][num_dihedral[m]] = atom4;
	num_dihedral[m]++;
      }
      if ((m = map(atom4)) >= 0) {
	dihedral_type[m][num_dihedral[m]] = itype;
	dihedral_atom1[m][num_dihedral[m]] = atom1;
	dihedral_atom2[m][num_dihedral[m]] = atom2;
	dihedral_atom3[m][num_dihedral[m]] = atom3;
	dihedral_atom4[m][num_dihedral[m]] = atom4;
	num_dihedral[m]++;
      }
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::unpack_impropers(int n, char *buf)
{
  int m,tmp,itype,atom1,atom2,atom3,atom4;
  char *next;
  int newton_bond = force->newton_bond;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    sscanf(buf,"%d %d %d %d %d %d",&tmp,&itype,&atom1,&atom2,&atom3,&atom4);
    if (atom1 <= 0 || atom1 > map_tag_max || 
	atom2 <= 0 || atom2 > map_tag_max || 
	atom3 <= 0 || atom3 > map_tag_max || 
	atom4 <= 0 || atom4 > map_tag_max)
      error->one("Invalid atom ID in Impropers section of data file");
    if (itype <= 0 || itype > nimpropertypes)
      error->one("Invalid improper type in Impropers section of data file");
    if ((m = map(atom2)) >= 0) {
      improper_type[m][num_improper[m]] = itype;
      improper_atom1[m][num_improper[m]] = atom1;
      improper_atom2[m][num_improper[m]] = atom2;
      improper_atom3[m][num_improper[m]] = atom3;
      improper_atom4[m][num_improper[m]] = atom4;
      num_improper[m]++;
    }
    if (newton_bond == 0) {
      if ((m = map(atom1)) >= 0) {
	improper_type[m][num_improper[m]] = itype;
	improper_atom1[m][num_improper[m]] = atom1;
	improper_atom2[m][num_improper[m]] = atom2;
	improper_atom3[m][num_improper[m]] = atom3;
	improper_atom4[m][num_improper[m]] = atom4;
	num_improper[m]++;
      }
      if ((m = map(atom3)) >= 0) {
	improper_type[m][num_improper[m]] = itype;
	improper_atom1[m][num_improper[m]] = atom1;
	improper_atom2[m][num_improper[m]] = atom2;
	improper_atom3[m][num_improper[m]] = atom3;
	improper_atom4[m][num_improper[m]] = atom4;
	num_improper[m]++;
      }
      if ((m = map(atom4)) >= 0) {
	improper_type[m][num_improper[m]] = itype;
	improper_atom1[m][num_improper[m]] = atom1;
	improper_atom2[m][num_improper[m]] = atom2;
	improper_atom3[m][num_improper[m]] = atom3;
	improper_atom4[m][num_improper[m]] = atom4;
	num_improper[m]++;
      }
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   allocate arrays of length ntypes
   only done after ntypes is set
------------------------------------------------------------------------- */

void Atom::allocate_type_arrays()
{
  if (mass_require) {
    mass = new double[ntypes+1];
    mass_setflag = new int[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) mass_setflag[itype] = 0;
  }
  if (dipole_require) {
    dipole = new double[ntypes+1];
    dipole_setflag = new int[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) dipole_setflag[itype] = 0;
  }
}

/* ----------------------------------------------------------------------
   set a mass and flag it as set
   called from reading of data file
------------------------------------------------------------------------- */

void Atom::set_mass(char *str)
{
  if (mass_require == 0) error->all("Cannot set mass for this atom style");

  int itype;
  double mass_one;
  sscanf(str,"%d %lg",&itype,&mass_one);

  if (itype < 1 || itype > ntypes) error->all("Invalid type for mass set");

  mass[itype] = mass_one;
  mass_setflag[itype] = 1;
}

/* ----------------------------------------------------------------------
   set a mass and flag it as set
   called from EAM pair routine
------------------------------------------------------------------------- */

void Atom::set_mass(int itype, double value)
{
  if (mass_require == 0) error->all("Cannot set mass for this atom style");
  if (itype < 1 || itype > ntypes) error->all("Invalid type for mass set");

  mass[itype] = value;
  mass_setflag[itype] = 1;
}

/* ----------------------------------------------------------------------
   set one or more masses and flag them as set
   called from reading of input script
------------------------------------------------------------------------- */

void Atom::set_mass(int narg, char **arg)
{
  if (mass_require == 0) error->all("Cannot set mass for this atom style");

  int lo,hi;
  force->bounds(arg[0],ntypes,lo,hi);
  if (lo < 1 || hi > ntypes) error->all("Invalid type for mass set");

  for (int itype = lo; itype <= hi; itype++) {
    mass[itype] = atof(arg[1]);
    mass_setflag[itype] = 1;
  }
}

/* ----------------------------------------------------------------------
   set all masses as read in from restart file
------------------------------------------------------------------------- */

void Atom::set_mass(double *values)
{
  for (int itype = 1; itype <= ntypes; itype++) {
    mass[itype] = values[itype];
    mass_setflag[itype] = 1;
  }
}

/* ----------------------------------------------------------------------
   check that all masses have been set
------------------------------------------------------------------------- */

void Atom::check_mass()
{
  if (!mass_require) return;
  for (int itype = 1; itype <= ntypes; itype++)
    if (mass_setflag[itype] == 0) error->all("All masses are not set");
}

/* ----------------------------------------------------------------------
   set a dipole moment and flag it as set
   called from reading of data file
------------------------------------------------------------------------- */

void Atom::set_dipole(char *str)
{
  if (dipole_require == 0) 
    error->all("Cannot set dipole for this atom style");

  int i;
  double dipole_one;
  sscanf(str,"%d %lg",&i,&dipole_one);

  dipole[i] = dipole_one;
  dipole_setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   set one or more dipole moments and flag them as set
   called from reading of input script
------------------------------------------------------------------------- */

void Atom::set_dipole(int narg, char **arg)
{
  if (dipole_require == 0) 
    error->all("Cannot set dipole for this atom style");

  int lo,hi;
  force->bounds(arg[0],ntypes,lo,hi);
  if (lo < 1 || hi > ntypes) error->all("Invalid type for dipole set");

  for (int itype = lo; itype <= hi; itype++) {
    dipole[itype] = atof(arg[1]);
    dipole_setflag[itype] = 1;
  }
}

/* ----------------------------------------------------------------------
   set all dipole moments as read in from restart file
------------------------------------------------------------------------- */

void Atom::set_dipole(double *values)
{
  for (int itype = 1; itype <= ntypes; itype++) {
    dipole[itype] = values[itype];
    dipole_setflag[itype] = 1;
  }
}

/* ----------------------------------------------------------------------
   check that all dipole moments have been set
------------------------------------------------------------------------- */

void Atom::check_dipole()
{
  if (!dipole_require) return;
  for (int itype = 1; itype <= ntypes; itype++)
    if (dipole_setflag[itype] == 0)
      error->all("All dipole moments are not set");
}

/* ----------------------------------------------------------------------
   register a callback to a fix so it can manage atom-based arrays
   happens when fix is created
   flag = 0 for grow, 1 for restart
------------------------------------------------------------------------- */

void Atom::add_callback(int flag)
{
  int ifix;

  // find the fix
  // if find NULL ptr:
  //   it's this one, since it is deleted at this point in creation
  // if don't find NULL ptr:
  //   i is set to nfix = new one currently being added at end of list

  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (modify->fix[ifix] == NULL) break;

  // add callback to lists, reallocating if necessary

  if (flag == 0) {
    if (nextra_grow == nextra_grow_max) {
      nextra_grow_max += DELTA_CALLBACK;
      extra_grow = (int *) 
	memory->srealloc(extra_grow,nextra_grow_max*sizeof(int),
			 "atom:extra_grow");
    }
    extra_grow[nextra_grow] = ifix;
    nextra_grow++;
  } else if (flag == 1) {
    if (nextra_restart == nextra_restart_max) {
      nextra_restart_max += DELTA_CALLBACK;
      extra_restart = (int *) 
	memory->srealloc(extra_restart,nextra_restart_max*sizeof(int),
			 "atom:extra_restart");
    }
    extra_restart[nextra_restart] = ifix;
    nextra_restart++;
  }
}

/* ----------------------------------------------------------------------
   unregister a callback to a fix
   happens when fix is deleted
   flag = 0 for grow, 1 for restart
------------------------------------------------------------------------- */

void Atom::delete_callback(char *id, int flag)
{
  int ifix;

  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(id,modify->fix[ifix]->id) == 0) break;

  // compact the list of callbacks

  if (flag == 0) {
    for (int i = ifix; i < nextra_grow-1; i++)
      extra_grow[i] = extra_grow[i+1];
    nextra_grow--;
  } else if (flag == 1) {
    for (int i = ifix; i < nextra_restart-1; i++)
      extra_restart[i] = extra_restart[i+1];
    nextra_restart--;
  }
}

/* ----------------------------------------------------------------------
   decrement ptrs in callback lists to fixes beyond the deleted ifix
   happens after fix is deleted
------------------------------------------------------------------------- */

void Atom::update_callback(int ifix)
{
  for (int i = 0; i < nextra_grow; i++)
    if (extra_grow[i] > ifix) extra_grow[i]--;
  for (int i = 0; i < nextra_restart; i++)
    if (extra_restart[i] > ifix) extra_restart[i]--;
}

/* ----------------------------------------------------------------------
   unpack n lines of atoms from buf
   set all atom values and defaults
------------------------------------------------------------------------- */

void Atom::unpack_data(int n, char *buf)
{
  int m,imagetmp,xptr,iptr;
  double xtmp,ytmp,ztmp;
  char *next;

  double subxlo = domain->subxlo;
  double subxhi = domain->subxhi;
  double subylo = domain->subylo;
  double subyhi = domain->subyhi;
  double subzlo = domain->subzlo;
  double subzhi = domain->subzhi;

  char **values = new char*[size_atom_actual];

  // xptr = which word in line is start of xyz coords
  // iptr = which word in line is start of ix,iy,iz image flags

  xptr = 2;
  if (molecular) xptr++;
  if (charge_allow) xptr++;
  if (style_granular) xptr += 2;
  if (style_peri) xptr += 2;

  int imageflag = 0;
  if (size_atom_actual > size_atom_valid) imageflag = 1;
  if (imageflag) iptr = size_atom_actual - 3;

  // loop over lines of atom data
  // tokenize the line into values
  // extract xyz coords and image flags
  // if atom is in my sub-domain, set its values
  // scan thru values in ascending order, storing ones appropriate to style

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    for (m = 1; m < size_atom_actual; m++)
      values[m] = strtok(NULL," \t\n\r\f");

    xtmp = atof(values[xptr]);
    ytmp = atof(values[xptr+1]);
    ztmp = atof(values[xptr+2]);
    if (imageflag)
      imagetmp = ((atoi(values[iptr+2]) + 512 & 1023) << 20) |
	((atoi(values[iptr+1]) + 512 & 1023) << 10) |
	(atoi(values[iptr]) + 512 & 1023);
    else imagetmp = (512 << 20) | (512 << 10) | 512;

    domain->remap(xtmp,ytmp,ztmp,imagetmp);
    if (xtmp >= subxlo && xtmp < subxhi &&
	ytmp >= subylo && ytmp < subyhi &&
	ztmp >= subzlo && ztmp < subzhi) {

      if (nlocal == nmax) grow(0);

      // parse quantities in prescribed order depending on atom style(s)

      m = 0;
      tag[nlocal] = atoi(values[m++]);
      if (tag[nlocal] <= 0)
	error->one("Invalid atom ID in Atoms section of data file");

      if (molecular) molecule[nlocal] = atoi(values[m++]);

      type[nlocal] = atoi(values[m++]);
      if (type[nlocal] <= 0 || type[nlocal] > ntypes)
	error->one("Invalid atom type in Atoms section of data file");

      if (charge_allow) q[nlocal] = atof(values[m++]);
      if (style_granular) {
	radius[nlocal] = 0.5 * atof(values[m++]);
	density[nlocal] = atof(values[m++]);
	if (force->dimension == 3)
	  rmass[nlocal] = 4.0*PI/3.0 *
	    radius[nlocal]*radius[nlocal]*radius[nlocal] * density[nlocal];
	else
	  rmass[nlocal] = PI * radius[nlocal]*radius[nlocal] * density[nlocal];
      }
      if (style_peri) {
	vfrac[nlocal] = atof(values[m++]);
	rmass[nlocal] = atof(values[m++]);
      }

      x[nlocal][0] = xtmp;
      x[nlocal][1] = ytmp;
      x[nlocal][2] = ztmp;
      m += 3;

      if (style_dipole) {
	mu[nlocal][0] = atof(values[m++]);
	mu[nlocal][1] = atof(values[m++]);
	mu[nlocal][2] = atof(values[m++]);
      }
      image[nlocal] = imagetmp;

      // initialize quantities not included in Atom section of data file

      mask[nlocal] = 1;
      v[nlocal][0] = 0.0;
      v[nlocal][1] = 0.0;
      v[nlocal][2] = 0.0;

      if (style_dipole) {
	omega[nlocal][0] = 0.0;
	omega[nlocal][1] = 0.0;
	omega[nlocal][2] = 0.0;
      }

      if (style_granular) {
	phix[nlocal][0] = 0.0;
	phix[nlocal][1] = 0.0;
	phix[nlocal][2] = 0.0;
	phiv[nlocal][0] = 0.0;
	phiv[nlocal][1] = 0.0;
	phiv[nlocal][2] = 0.0;
      }

      if (bonds_allow) num_bond[nlocal] = 0;
      if (angles_allow) num_angle[nlocal] = 0;
      if (dihedrals_allow) num_dihedral[nlocal] = 0;
      if (impropers_allow) num_improper[nlocal] = 0;
      
      nlocal++;
    }

    buf = next + 1;
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   create one atom of type itype at x0,y0,z0
   set all other values to defaults
------------------------------------------------------------------------- */

void Atom::create_one(int itype, double x0, double y0, double z0)
{
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = x0;
  x[nlocal][1] = y0;
  x[nlocal][2] = z0;
  mask[nlocal] = 1;
  image[nlocal] = (512 << 20) | (512 << 10) | 512;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  if (charge_allow) q[nlocal] = 0.0;

  if (style_dipole) {
    mu[nlocal][0] = 0.0;
    mu[nlocal][1] = 0.0;
    mu[nlocal][2] = 0.0;
    omega[nlocal][0] = 0.0;
    omega[nlocal][1] = 0.0;
    omega[nlocal][2] = 0.0;
  }

  if (style_granular) {
    radius[nlocal] = 0.5;
    density[nlocal] = 1.0;
    if (force->dimension == 3)
      rmass[nlocal] = 4.0*PI/3.0 *
	radius[nlocal]*radius[nlocal]*radius[nlocal] * density[nlocal];
    else
      rmass[nlocal] = PI * radius[nlocal]*radius[nlocal] * density[nlocal];
    phix[nlocal][0] = 0.0;
    phix[nlocal][1] = 0.0;
    phix[nlocal][2] = 0.0;
    phiv[nlocal][0] = 0.0;
    phiv[nlocal][1] = 0.0;
    phiv[nlocal][2] = 0.0;
  }

  if (style_peri) {
    vfrac[nlocal] = 1.0;
    rmass[nlocal] = 1.0;
  }

  if (molecular) {
    molecule[nlocal] = 0;
    if (bonds_allow) num_bond[nlocal] = 0;
    if (angles_allow) num_angle[nlocal] = 0;
    if (dihedrals_allow) num_dihedral[nlocal] = 0;
    if (impropers_allow) num_improper[nlocal] = 0;
  }

  nlocal++;
}

/* ---------------------------------------------------------------------- */

int Atom::size_restart()
{
  int i;

  int n = 0;
  for (i = 0; i < nlocal; i++) {
    n += 11;
    if (charge_allow) n++;
    if (style_dipole) n += 6;
    if (style_granular) n += 8;
    if (style_peri) n += 5;
    if (molecular) {
      n++;
      if (bonds_allow) n += 1 + 2*num_bond[i];
      if (angles_allow) n += 1 + 4*num_angle[i];
      if (dihedrals_allow) n += 1 + 5*num_dihedral[i];
      if (impropers_allow) n += 1 + 5*num_improper[i];
    }
  }

  if (nextra_restart)
    for (int iextra = 0; iextra < nextra_restart; iextra++) 
      for (i = 0; i < nlocal; i++)
	n += modify->fix[extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack all atom quantities for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   interaction types may be negative, but write as positive   
------------------------------------------------------------------------- */

int Atom::pack_restart(int i, double *buf)
{
  int k;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m++] = image[i];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  if (charge_allow) buf[m++] = q[i];

  if (style_dipole) {
    buf[m++] = mu[i][0];
    buf[m++] = mu[i][1];
    buf[m++] = mu[i][2];
    buf[m++] = omega[i][0];
    buf[m++] = omega[i][1];
    buf[m++] = omega[i][2];
  }

  if (style_granular) {
    buf[m++] = phix[i][0];
    buf[m++] = phix[i][1];
    buf[m++] = phix[i][2];
    buf[m++] = phiv[i][0];
    buf[m++] = phiv[i][1];
    buf[m++] = phiv[i][2];
    buf[m++] = radius[i];
    buf[m++] = density[i];
  }

  if (style_peri) {
    buf[m++] = vfrac[i];
    buf[m++] = rmass[i];
  }

  if (molecular) {
    buf[m++] = molecule[i];

    if (bonds_allow) {
      buf[m++] = num_bond[i];
      for (k = 0; k < num_bond[i]; k++) {
	buf[m++] = MAX(bond_type[i][k],-bond_type[i][k]);
	buf[m++] = bond_atom[i][k];
      }
    }

    if (angles_allow) {
      buf[m++] = num_angle[i];
      for (k = 0; k < num_angle[i]; k++) {
	buf[m++] = MAX(angle_type[i][k],-angle_type[i][k]);
	buf[m++] = angle_atom1[i][k];
	buf[m++] = angle_atom2[i][k];
	buf[m++] = angle_atom3[i][k];
      }
    }

    if (dihedrals_allow) {
      buf[m++] = num_dihedral[i];
      for (k = 0; k < num_dihedral[i]; k++) {
	buf[m++] = MAX(dihedral_type[i][k],-dihedral_type[i][k]);
	buf[m++] = dihedral_atom1[i][k];
	buf[m++] = dihedral_atom2[i][k];
	buf[m++] = dihedral_atom3[i][k];
	buf[m++] = dihedral_atom4[i][k];
      }
    }

    if (impropers_allow) {
      buf[m++] = num_improper[i];
      for (k = 0; k < num_improper[i]; k++) {
	buf[m++] = MAX(improper_type[i][k],-improper_type[i][k]);
	buf[m++] = improper_atom1[i][k];
	buf[m++] = improper_atom2[i][k];
	buf[m++] = improper_atom3[i][k];
	buf[m++] = improper_atom4[i][k];
      }
    }
  }

  if (nextra_restart)
    for (int iextra = 0; iextra < nextra_restart; iextra++) 
      m += modify->fix[extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack all atom quantities from restart file including extra quantities
------------------------------------------------------------------------- */

int Atom::unpack_restart(double *buf)
{
  int k;

  if (nlocal == nmax) {
    grow(0);
    if (nextra_store)
      extra = memory->grow_2d_double_array(extra,nmax,nextra_store,
					   "atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = static_cast<int> (buf[m++]);
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  if (charge_allow) q[nlocal] = buf[m++];

  if (style_dipole) {
    mu[nlocal][0] = buf[m++];
    mu[nlocal][1] = buf[m++];
    mu[nlocal][2] = buf[m++];
    omega[nlocal][0] = buf[m++];
    omega[nlocal][1] = buf[m++];
    omega[nlocal][2] = buf[m++];
  }

  if (style_granular) {
    phix[nlocal][0] = buf[m++];
    phix[nlocal][1] = buf[m++];
    phix[nlocal][2] = buf[m++];
    phiv[nlocal][0] = buf[m++];
    phiv[nlocal][1] = buf[m++];
    phiv[nlocal][2] = buf[m++];
    radius[nlocal] = buf[m++];
    density[nlocal] = buf[m++];
    if (force->dimension == 3) 
      rmass[nlocal] = 4.0*PI/3.0 * 
	radius[nlocal]*radius[nlocal]*radius[nlocal] * density[nlocal];
    else
      rmass[nlocal] = PI * radius[nlocal]*radius[nlocal] * density[nlocal];
  }

  if (style_peri) {
    vfrac[nlocal] = buf[m++];
    rmass[nlocal] = buf[m++];
  }

  if (molecular) {
    molecule[nlocal] = static_cast<int> (buf[m++]);

    if (bonds_allow) {
      num_bond[nlocal] = static_cast<int> (buf[m++]);
      for (k = 0; k < num_bond[nlocal]; k++) {
	bond_type[nlocal][k] = static_cast<int> (buf[m++]);
	bond_atom[nlocal][k] = static_cast<int> (buf[m++]);
      }
    }

    if (angles_allow) {
      num_angle[nlocal] = static_cast<int> (buf[m++]);
      for (k = 0; k < num_angle[nlocal]; k++) {
	angle_type[nlocal][k] = static_cast<int> (buf[m++]);
	angle_atom1[nlocal][k] = static_cast<int> (buf[m++]);
	angle_atom2[nlocal][k] = static_cast<int> (buf[m++]);
	angle_atom3[nlocal][k] = static_cast<int> (buf[m++]);
      }
    }

    if (dihedrals_allow) {
      num_dihedral[nlocal] = static_cast<int> (buf[m++]);
      for (k = 0; k < num_dihedral[nlocal]; k++) {
	dihedral_type[nlocal][k] = static_cast<int> (buf[m++]);
	dihedral_atom1[nlocal][k] = static_cast<int> (buf[m++]);
	dihedral_atom2[nlocal][k] = static_cast<int> (buf[m++]);
	dihedral_atom3[nlocal][k] = static_cast<int> (buf[m++]);
	dihedral_atom4[nlocal][k] = static_cast<int> (buf[m++]);
      }
    }

    if (impropers_allow) {
      num_improper[nlocal] = static_cast<int> (buf[m++]);
      for (k = 0; k < num_improper[nlocal]; k++) {
	improper_type[nlocal][k] = static_cast<int> (buf[m++]);
	improper_atom1[nlocal][k] = static_cast<int> (buf[m++]);
	improper_atom2[nlocal][k] = static_cast<int> (buf[m++]);
	improper_atom3[nlocal][k] = static_cast<int> (buf[m++]);
	improper_atom4[nlocal][k] = static_cast<int> (buf[m++]);
      }
    }
  }

  if (nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory 
------------------------------------------------------------------------- */

int Atom::memory_usage()
{
  int bytes = 0;

  bytes += 4 * nmax * sizeof(int);
  bytes += 3 * nmax*3 * sizeof(double);

  if (charge_allow) bytes += 1 * nmax * sizeof(double);     // q
  if (mass_allow) bytes += 1 * nmax * sizeof(double);       // rmass
  if (style_dipole) bytes += 3 * nmax*3 * sizeof(double);   // mu,omega,torque
  if (style_granular) {
    bytes += 3 * nmax*3 * sizeof(double);                   // phix,phiv,phia
    bytes += 2 * nmax * sizeof(double);                     // radius,density
  }
  if (style_peri) bytes += nmax * sizeof(double);           // vfrac

  if (map_style == 1)
    bytes += map_tag_max * sizeof(int);
  else if (map_style == 2) {
    bytes += map_nbucket*sizeof(int);
    bytes += map_nhash*sizeof(HashElem);
  }

  if (molecular) {
    bytes += nmax * sizeof(int);
    bytes += nmax*3 * sizeof(int);
    bytes += nmax*maxspecial * sizeof(int);

    if (bonds_allow) {
      bytes += nmax * sizeof(int);
      bytes += 2 * nmax*bond_per_atom * sizeof(int);
    }

    if (angles_allow) {
      bytes += nmax * sizeof(int);
      bytes += 4 * nmax*angle_per_atom * sizeof(int);
    }

    if (dihedrals_allow) {
      bytes += nmax * sizeof(int);
      bytes += 5 * nmax*dihedral_per_atom * sizeof(int);
    }

    if (impropers_allow) {
      bytes += nmax * sizeof(int);
      bytes += 5 * nmax*improper_per_atom * sizeof(int);
    }
  }

  return bytes;
}
