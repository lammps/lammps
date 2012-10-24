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

#include "math.h"
#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EXTRA 1000

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
    error->all(FLERR,"Cannot create an atom map unless atoms have IDs");

  int max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&map_tag_max,1,MPI_INT,MPI_MAX,world);

  memory->destroy(sametag);
  smax = nlocal + nghost + EXTRA;
  memory->create(sametag,smax,"atom:sametag");

  if (map_style == 1) {
    memory->create(map_array,map_tag_max+1,"atom:map_array");
    for (int i = 0; i <= map_tag_max; i++) map_array[i] = -1;

  } else {

    // map_nhash = max # of atoms that can be hashed on this proc
    // set to max of ave atoms/proc or atoms I can store
    // multiply by 2, require at least 1000
    // doubling means hash table will be re-init only rarely

    int nper = static_cast<int> (natoms/comm->nprocs);
    map_nhash = MAX(nper,nmax);
    map_nhash *= 2;
    map_nhash = MAX(map_nhash,1000);

    // map_nbucket = prime just larger than map_nhash

    map_nbucket = next_prime(map_nhash);

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
    for (int i = 0; i < nall; i++) {
      sametag[i] = -1;
      map_array[tag[i]] = -1;
    }

  } else {
    int previous,global,ibucket,index;
    int nall = nlocal + nghost;
    for (int i = 0; i < nall; i++) {
      sametag[i] = -1;

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
  int nall = nlocal + nghost;
  if (nall > smax) {
    smax = nall + EXTRA;
    memory->destroy(sametag);
    memory->create(sametag,smax,"atom:sametag");
  }

  if (map_style == 1) {
    for (int i = nall-1; i >= 0 ; i--) {
      sametag[i] = map_array[tag[i]];
      map_array[tag[i]] = i;
    }

  } else {
    int previous,global,ibucket,index;
    if (nall > map_nhash) map_init();

    for (int i = nall-1; i >= 0 ; i--) {
      sametag[i] = map_find_hash(tag[i]);

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
   called by Special class
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
  memory->destroy(sametag);
  sametag = NULL;

  if (map_style == 1) {
    memory->destroy(map_array);
    map_array = NULL;
  } else {
    if (map_nhash) {
      delete [] map_bucket;
      delete [] map_hash;
      map_bucket = NULL;
      map_hash = NULL;
    }
    map_nhash = 0;
  }
  map_tag_max = 0;
}

/* ----------------------------------------------------------------------
   lookup global ID in hash table, return local index
   called by map() in atom.h
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
   return next prime larger than n
------------------------------------------------------------------------- */

int Atom::next_prime(int n)
{
  int factor;

  int nprime = n+1;
  if (nprime % 2 == 0) nprime++;
  int root = static_cast<int> (sqrt(1.0*n)) + 2;

  while (nprime <= MAXSMALLINT) {
    for (factor = 3; factor < root; factor++)
      if (nprime % factor == 0) break;
    if (factor == root) return nprime;
    nprime += 2;
  }

  return MAXSMALLINT;
}
