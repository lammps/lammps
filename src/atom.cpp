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

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "atom_vec.h"
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

#define AtomInclude
#include "style.h"
#undef AtomInclude

using namespace LAMMPS_NS;

#define DELTA 1
#define DELTA_MEMSTR 1024
#define EPSILON 1.0e-6

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

Atom::Atom(LAMMPS *lmp) : Pointers(lmp)
{
  natoms = nlocal = nghost = nmax = 0;
  ntypes = 0;
  nbondtypes = nangletypes = ndihedraltypes = nimpropertypes = 0;
  nbonds = nangles = ndihedrals = nimpropers = 0;
  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;
  extra_bond_per_atom = 0;

  firstgroupname = NULL;

  // initialize atom arrays
  // customize by adding new array

  tag = type = mask = image = NULL;
  x = v = f = NULL;

  molecule = NULL;
  q = NULL;
  mu = NULL;
  quat = omega = angmom = torque = NULL;
  radius = density = rmass = NULL;
  vfrac = s0 = NULL;
  x0 = NULL;

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

  // initialize atom array existence flags
  // customize by adding new flag

  molecule_flag = 0;
  q_flag = mu_flag = 0;
  quat_flag = omega_flag = angmom_flag = torque_flag = 0;
  radius_flag = density_flag = rmass_flag = vfrac_flag = 0;

  // ntype-length arrays

  mass = NULL;
  mass_setflag = NULL;
  shape = NULL;
  shape_setflag = NULL;
  dipole = NULL;
  dipole_setflag = NULL;

  // callback lists & extra restart info

  nextra_grow = nextra_restart = 0;
  extra_grow = extra_restart = NULL;
  nextra_grow_max = nextra_restart_max = 0;
  nextra_store = 0;
  extra = NULL;

  // default mapping values and hash table primes

  tag_enable = 1;
  map_style = 0;
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

  // default atom style = atomic

  atom_style = NULL;
  avec = NULL;
  create_avec("atomic",0,NULL);
}

/* ---------------------------------------------------------------------- */

Atom::~Atom()
{
  delete [] atom_style;
  delete avec;
  delete [] firstgroupname;

  // delete atom arrays
  // customize by adding new array

  memory->sfree(tag);
  memory->sfree(type);
  memory->sfree(mask);
  memory->sfree(image);
  memory->destroy_2d_double_array(x);
  memory->destroy_2d_double_array(v);
  memory->destroy_2d_double_array(f);

  memory->sfree(q);
  memory->destroy_2d_double_array(mu);
  memory->destroy_2d_double_array(quat);
  memory->destroy_2d_double_array(omega);
  memory->destroy_2d_double_array(angmom);
  memory->destroy_2d_double_array(torque);

  memory->sfree(radius);
  memory->sfree(density);
  memory->sfree(rmass);
  memory->sfree(vfrac);
  memory->sfree(s0);
  memory->destroy_2d_double_array(x0);

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

  // delete per-type arrays

  delete [] mass;
  delete [] mass_setflag;
  memory->destroy_2d_double_array(shape);
  delete [] shape_setflag;
  delete [] dipole;
  delete [] dipole_setflag;

  memory->sfree(extra_grow);
  memory->sfree(extra_restart);
  memory->destroy_2d_double_array(extra);

  map_delete();
  delete [] primes;
}

/* ----------------------------------------------------------------------
   copy modify settings from old Atom class to current Atom class
------------------------------------------------------------------------- */

void Atom::settings(Atom *old)
{
  map_style = old->map_style;
}

/* ----------------------------------------------------------------------
   create an AtomVec style
   called from input script, restart file, replicate
------------------------------------------------------------------------- */

void Atom::create_avec(const char *style, int narg, char **arg)
{
  delete [] atom_style;
  if (avec) delete avec;

  avec = new_avec(style,narg,arg);
  int n = strlen(style) + 1;
  atom_style = new char[n];
  strcpy(atom_style,style);

  // if molecular system, default is to have array map

  molecular = avec->molecular;
  if (map_style == 0 && molecular) map_style = 1;
}

/* ----------------------------------------------------------------------
   generate an AtomVec class
------------------------------------------------------------------------- */

AtomVec *Atom::new_avec(const char *style, int narg, char **arg)
{
  if (0) return NULL;

#define AtomClass
#define AtomStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(lmp,narg,arg);
#include "style.h"
#undef AtomClass

  else error->all("Invalid atom style");
  return NULL;
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
  check_shape();
  check_dipole();

  // setup of firstgroup

  if (firstgroupname) {
    firstgroup = group->find(firstgroupname);
    if (firstgroup < 0)
      error->all("Could not find atom_modify first group ID");
  } else firstgroup = -1;

  // init sub-style

  avec->init();
}

/* ----------------------------------------------------------------------
   return 1 if style matches atom style or hybrid sub-style
   else return 0
------------------------------------------------------------------------- */

int Atom::style_match(const char *style)
{
  if (strcmp(atom_style,style) == 0) return 1;
  else if (strcmp(atom_style,"hybrid") == 0) {
    AtomVecHybrid *avec_hybrid = (AtomVecHybrid *) avec;
    for (int i = 0; i < avec_hybrid->nstyles; i++)
      if (strcmp(avec_hybrid->keywords[i],style) == 0) return 1;
  }
  return 0;
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
    } else if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) error->all("Illegal atom_modify command");
      if (strcmp(arg[iarg+1],"all") == 0) delete [] firstgroupname;
      else {
	int n = strlen(arg[iarg+1]) + 1;
	firstgroupname = new char[n];
	strcpy(firstgroupname,arg[iarg+1]);
      }
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
   check that atom IDs span range from 1 to Natoms
   return 0 if mintag != 1 or maxtag != Natoms
   return 1 if OK
   doesn't actually check if all tag values are used
------------------------------------------------------------------------- */

int Atom::tag_consecutive()
{
  int idmin = static_cast<int> (natoms);
  int idmax = 0;
  
  for (int i = 0; i < nlocal; i++) {
    idmin = MIN(idmin,tag[i]);
    idmax = MAX(idmax,tag[i]);
  }
  int idminall,idmaxall;
  MPI_Allreduce(&idmin,&idminall,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&idmax,&idmaxall,1,MPI_INT,MPI_MAX,world);

  if (idminall != 1 || idmaxall != static_cast<int> (natoms)) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   count and return words in a single line
   make copy of line before using strtok so as not to change line
   trim anything from '#' onward
------------------------------------------------------------------------- */

int Atom::count_words(const char *line)
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
   unpack n lines from Atom section of data file
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_atoms(int n, char *buf)
{
  int m,imagedata,xptr,iptr;
  double xdata[3],lamda[3],sublo[3],subhi[3];
  double *coord;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = count_words(buf);
  *next = '\n';

  if (nwords != avec->size_data_atom && nwords != avec->size_data_atom + 3)
    error->all("Incorrect atom format in data file");

  char **values = new char*[nwords];

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all data atoms will be owned even with round-off

  int triclinic = domain->triclinic;
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
  }

  if (domain->xperiodic) {
    if (comm->myloc[0] == 0) sublo[0] -= EPSILON;
    if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] += EPSILON;
  }
  if (domain->yperiodic) {
    if (comm->myloc[1] == 0) sublo[1] -= EPSILON;
    if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] += EPSILON;
  }
  if (domain->zperiodic) {
    if (comm->myloc[2] == 0) sublo[2] -= EPSILON;
    if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] += EPSILON;
  }

  // xptr = which word in line starts xyz coords
  // iptr = which word in line starts ix,iy,iz image flags

  xptr = avec->xcol_data - 1;
  int imageflag = 0;
  if (nwords > avec->size_data_atom) imageflag = 1;
  if (imageflag) iptr = nwords - 3;

  // loop over lines of atom data
  // tokenize the line into values
  // extract xyz coords and image flags
  // remap atom into simulation box
  // if atom is in my sub-domain, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    for (m = 1; m < nwords; m++)
      values[m] = strtok(NULL," \t\n\r\f");

    if (imageflag)
      imagedata = ((atoi(values[iptr+2]) + 512 & 1023) << 20) |
	((atoi(values[iptr+1]) + 512 & 1023) << 10) |
	(atoi(values[iptr]) + 512 & 1023);
    else imagedata = (512 << 20) | (512 << 10) | 512;

    xdata[0] = atof(values[xptr]);
    xdata[1] = atof(values[xptr+1]);
    xdata[2] = atof(values[xptr+2]);
    domain->remap(xdata,imagedata);
    if (triclinic) {
      domain->x2lamda(xdata,lamda);
      coord = lamda;
    } else coord = xdata;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
	coord[1] >= sublo[1] && coord[1] < subhi[1] &&
	coord[2] >= sublo[2] && coord[2] < subhi[2])
      avec->data_atom(xdata,imagedata,values);

    buf = next + 1;
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   unpack n lines from Velocity section of data file
   check that atom IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_vels(int n, char *buf)
{
  int j,m,tagdata;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = count_words(buf);
  *next = '\n';

  if (nwords != avec->size_data_vel)
    error->all("Incorrect velocity format in data file");

  char **values = new char*[nwords];

  // loop over lines of atom velocities
  // tokenize the line into values
  // if I own atom tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL," \t\n\r\f");

    tagdata = atoi(values[0]);
    if (tagdata <= 0 || tagdata > map_tag_max)
      error->one("Invalid atom ID in Velocities section of data file");
    if ((m = map(tagdata)) >= 0) avec->data_vel(m,&values[1]);

    buf = next + 1;
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::data_bonds(int n, char *buf)
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

void Atom::data_angles(int n, char *buf)
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

void Atom::data_dihedrals(int n, char *buf)
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

void Atom::data_impropers(int n, char *buf)
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
  if (avec->mass_type) {
    mass = new double[ntypes+1];
    mass_setflag = new int[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) mass_setflag[itype] = 0;
  }
  if (avec->shape_type) {
    shape = memory->create_2d_double_array(ntypes+1,3,"atom:shape");
    shape_setflag = new int[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) shape_setflag[itype] = 0;
  }    
  if (avec->dipole_type) {
    dipole = new double[ntypes+1];
    dipole_setflag = new int[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) dipole_setflag[itype] = 0;
  }
}

/* ----------------------------------------------------------------------
   set a mass and flag it as set
   called from reading of data file
------------------------------------------------------------------------- */

void Atom::set_mass(const char *str)
{
  if (mass == NULL) error->all("Cannot set mass for this atom style");

  int itype;
  double mass_one;
  int n = sscanf(str,"%d %lg",&itype,&mass_one);
  if (n != 2) error->all("Invalid mass line in data file");

  if (itype < 1 || itype > ntypes) error->all("Invalid type for mass set");

  mass[itype] = mass_one;
  mass_setflag[itype] = 1;

  if (mass[itype] <= 0.0) error->all("Invalid mass value");
}

/* ----------------------------------------------------------------------
   set a mass and flag it as set
   called from EAM pair routine
------------------------------------------------------------------------- */

void Atom::set_mass(int itype, double value)
{
  if (mass == NULL) error->all("Cannot set mass for this atom style");
  if (itype < 1 || itype > ntypes) error->all("Invalid type for mass set");

  mass[itype] = value;
  mass_setflag[itype] = 1;

  if (mass[itype] <= 0.0) error->all("Invalid mass value");
}

/* ----------------------------------------------------------------------
   set one or more masses and flag them as set
   called from reading of input script
------------------------------------------------------------------------- */

void Atom::set_mass(int narg, char **arg)
{
  if (mass == NULL) error->all("Cannot set mass for this atom style");

  int lo,hi;
  force->bounds(arg[0],ntypes,lo,hi);
  if (lo < 1 || hi > ntypes) error->all("Invalid type for mass set");

  for (int itype = lo; itype <= hi; itype++) {
    mass[itype] = atof(arg[1]);
    mass_setflag[itype] = 1;

    if (mass[itype] <= 0.0) error->all("Invalid mass value");
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
  if (mass == NULL) return;
  for (int itype = 1; itype <= ntypes; itype++)
    if (mass_setflag[itype] == 0) error->all("All masses are not set");
}

/* ----------------------------------------------------------------------
   set particle shape and flag it as set
   called from reading of data file
------------------------------------------------------------------------- */

void Atom::set_shape(const char *str)
{
  if (shape == NULL) error->all("Cannot set shape for this atom style");

  int itype;
  double a,b,c;
  int n = sscanf(str,"%d %lg %lg %lg",&itype,&a,&b,&c);
  if (n != 4) error->all("Invalid shape line in data file");

  if (itype < 1 || itype > ntypes) error->all("Invalid type for shape set");

  // store shape as radius, though specified as diameter

  shape[itype][0] = 0.5*a;
  shape[itype][1] = 0.5*b;
  shape[itype][2] = 0.5*c;
  shape_setflag[itype] = 1;

  if (shape[itype][0] < 0.0 || shape[itype][1] < 0.0 || 
      shape[itype][2] < 0.0)
    error->all("Invalid shape value");
  if (shape[itype][0] > 0.0 || shape[itype][1] > 0.0 || 
      shape[itype][2] > 0.0) {
    if (shape[itype][0] == 0.0 || shape[itype][1] == 0.0 || 
	shape[itype][2] == 0.0)
      error->all("Invalid shape value");
  }
}

/* ----------------------------------------------------------------------
   set one or more particle shapes and flag them as set
   called from reading of input script
------------------------------------------------------------------------- */

void Atom::set_shape(int narg, char **arg)
{
  if (shape == NULL) error->all("Cannot set shape for this atom style");

  int lo,hi;
  force->bounds(arg[0],ntypes,lo,hi);
  if (lo < 1 || hi > ntypes) 
	error->all("Invalid type for shape set");

  // store shape as radius, though specified as diameter

  for (int itype = lo; itype <= hi; itype++) {
    shape[itype][0] = 0.5*atof(arg[1]);
    shape[itype][1] = 0.5*atof(arg[2]);
    shape[itype][2] = 0.5*atof(arg[3]);
    shape_setflag[itype] = 1;

    if (shape[itype][0] < 0.0 || shape[itype][1] < 0.0 || 
	shape[itype][2] < 0.0)
      error->all("Invalid shape value");
    if (shape[itype][0] > 0.0 || shape[itype][1] > 0.0 || 
	shape[itype][2] > 0.0) {
      if (shape[itype][0] == 0.0 || shape[itype][1] == 0.0 || 
	  shape[itype][2] == 0.0)
	error->all("Invalid shape value");
    }
  }
}

/* ----------------------------------------------------------------------
   set all particle shapes as read in from restart file
------------------------------------------------------------------------- */

void Atom::set_shape(double **values)
{
  for (int itype = 1; itype <= ntypes; itype++) {
    shape[itype][0] = values[itype][0];
    shape[itype][1] = values[itype][1];
    shape[itype][2] = values[itype][2];
    shape_setflag[itype] = 1;
  }
}

/* ----------------------------------------------------------------------
   check that all particle shapes have been set
------------------------------------------------------------------------- */

void Atom::check_shape()
{
  if (shape == NULL) return;
  for (int itype = 1; itype <= ntypes; itype++)
    if (shape_setflag[itype] == 0) error->all("All shapes are not set");
}

/* ----------------------------------------------------------------------
   set a dipole moment and flag it as set
   called from reading of data file
------------------------------------------------------------------------- */

void Atom::set_dipole(const char *str)
{
  if (dipole == NULL) error->all("Cannot set dipole for this atom style");

  int itype;
  double dipole_one;
  int n = sscanf(str,"%d %lg",&itype,&dipole_one);
  if (n != 2) error->all("Invalid dipole line in data file");

  dipole[itype] = dipole_one;
  dipole_setflag[itype] = 1;

  if (dipole[itype] < 0.0) error->all("Invalid dipole value");
}

/* ----------------------------------------------------------------------
   set one or more dipole moments and flag them as set
   called from reading of input script
------------------------------------------------------------------------- */

void Atom::set_dipole(int narg, char **arg)
{
  if (dipole == NULL) error->all("Cannot set dipole for this atom style");

  int lo,hi;
  force->bounds(arg[0],ntypes,lo,hi);
  if (lo < 1 || hi > ntypes) error->all("Invalid type for dipole set");

  for (int itype = lo; itype <= hi; itype++) {
    dipole[itype] = atof(arg[1]);
    dipole_setflag[itype] = 1;

    if (dipole[itype] < 0.0) error->all("Invalid dipole value");
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
  if (dipole == NULL) return;
  for (int itype = 1; itype <= ntypes; itype++)
    if (dipole_setflag[itype] == 0)
      error->all("All dipole moments are not set");
}

/* ----------------------------------------------------------------------
   reorder owned atoms so those in firstgroup appear first
   called by comm->exchange() if atom_modify first group is set
   only owned atoms exist at this point, no ghost atoms
------------------------------------------------------------------------- */

void Atom::first_reorder()
{
  // insure there is one extra atom location at end of arrays for swaps

  if (nlocal == nmax) avec->grow(0);

  // loop over owned atoms
  // nfirst = index of first atom not in firstgroup
  // when find firstgroup atom out of place, swap it with atom nfirst

  int bitmask = group->bitmask[firstgroup];
  nfirst = 0;
  while (nfirst < nlocal && mask[nfirst] & bitmask) nfirst++;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & bitmask && i > nfirst) {
      avec->copy(i,nlocal);
      avec->copy(nfirst,i);
      avec->copy(nlocal,nfirst);
      while (nfirst < nlocal && mask[nfirst] & bitmask) nfirst++;
    }
  }
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
      nextra_grow_max += DELTA;
      extra_grow = (int *) 
	memory->srealloc(extra_grow,nextra_grow_max*sizeof(int),
			 "atom:extra_grow");
    }
    extra_grow[nextra_grow] = ifix;
    nextra_grow++;
  } else if (flag == 1) {
    if (nextra_restart == nextra_restart_max) {
      nextra_restart_max += DELTA;
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

void Atom::delete_callback(const char *id, int flag)
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
   return a pointer to a named internal variable
   if don't recognize name, return NULL
------------------------------------------------------------------------- */

void *Atom::extract(char *name)
{
  if (strcmp(name,"natoms") == 0) return (void *) &natoms;
  if (strcmp(name,"nlocal") == 0) return (void *) &nlocal;

  if (strcmp(name,"id") == 0) return (void *) tag;
  if (strcmp(name,"type") == 0) return (void *) type;
  if (strcmp(name,"x") == 0) return (void *) x;
  if (strcmp(name,"v") == 0) return (void *) v;
  if (strcmp(name,"f") == 0) return (void *) f;
  if (strcmp(name,"mass") == 0) return (void *) mass;
  if (strcmp(name,"rmass") == 0) return (void *) rmass;

  return NULL;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   call to avec sums per-atom vectors
   add in global to local mapping storage
------------------------------------------------------------------------- */

double Atom::memory_usage()
{
  memlength = DELTA_MEMSTR;
  memstr = (char *) memory->smalloc(memlength*sizeof(char),"atom:memstr");
  memstr[0] = '\0';

  double bytes = avec->memory_usage();
  if (map_style == 1)
    bytes += map_tag_max * sizeof(int);
  else if (map_style == 2) {
    bytes += map_nbucket*sizeof(int);
    bytes += map_nhash*sizeof(HashElem);
  }

  memory->sfree(memstr);
  return bytes;
}

/* ----------------------------------------------------------------------
   accumulate per-atom vec names in memstr, padded by spaces
   return 1 if padded str is not already in memlist, else 0
------------------------------------------------------------------------- */

int Atom::memcheck(const char *str)
{
  int n = strlen(str) + 3;
  char *padded = new char[n];
  strcpy(padded," ");
  strcat(padded,str);
  strcat(padded," ");
  
  if (strstr(memstr,padded)) {
    delete [] padded;
    return 0;
  }

  if (strlen(memstr) + n >= memlength) {
    memlength += DELTA_MEMSTR;
    memstr = (char *) memory->srealloc(memstr,memlength*sizeof(char),
				       "atom:memstr");
  }

  strcat(memstr,padded);
  delete [] padded;
  return 1;
}
