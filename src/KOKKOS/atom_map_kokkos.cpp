/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_kokkos.h"

#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

#define EXTRA 1000

/* ----------------------------------------------------------------------
   allocate and initialize array or hash table for global -> local map
   for array option:
     array length = 1 to map_tag_max
     set entire array to -1 as initial values
   for hash option:
     map_nhash = length of hash table
     map_nbucket = # of hash buckets, prime larger than map_nhash * 2
       so buckets will only be filled with 0 or 1 atoms on average
------------------------------------------------------------------------- */

void AtomKokkos::map_init(int check)
{
  // check for new map style if max atomID changed (check = 1 = default)
  // recreate = 1 if must delete old map and create new map
  // recreate = 0 if can re-use old map w/out realloc and just adjust settings
  // map_maxarray/map_nhash initially -1, to force recreate even when no atoms

  int recreate = 0;
  if (check) recreate = map_style_set();

  if (map_style == MAP_ARRAY && map_tag_max > map_maxarray)
    recreate = 1;
  else if (map_style == MAP_HASH && nlocal + nghost > map_nhash)
    recreate = 1;

  // if not recreating:
  // for array, initialize current map_tag_max values
  // for hash, set all buckets to empty, put all entries in free list

  if (!recreate) {
    if (map_style == MAP_ARRAY) {
      for (int i = 0; i <= map_tag_max; i++) map_array[i] = -1;
    } else {
      for (int i = 0; i < map_nbucket; i++) map_bucket[i] = -1;
      map_nused = 0;
      map_free = 0;
      for (int i = 0; i < map_nhash; i++) map_hash[i].next = i + 1;
      if (map_nhash > 0) map_hash[map_nhash - 1].next = -1;
    }

    // recreating: delete old map and create new one for array or hash

  } else {
    map_delete();

    if (map_style == MAP_ARRAY) {
      map_maxarray = map_tag_max;
      memoryKK->create_kokkos(k_map_array, map_array, map_maxarray + 1, "atom:map_array");
      for (int i = 0; i <= map_tag_max; i++) map_array[i] = -1;

    } else {

      // map_nhash = max # of atoms that can be hashed on this proc
      // set to max of ave atoms/proc or atoms I can store
      // multiply by 2, require at least 1000
      // doubling means hash table will need to be re-init only rarely

      int nper = static_cast<int>(natoms / comm->nprocs);
      map_nhash = MAX(nper, nmax);
      map_nhash *= 2;
      map_nhash = MAX(map_nhash, 1000);

      // map_nbucket = prime just larger than map_nhash
      // next_prime() should be fast enough,
      //   about 10% of odd integers are prime above 1M

      map_nbucket = next_prime(map_nhash);

      // set all buckets to empty
      // set hash to map_nhash in length
      // put all hash entries in free list and point them to each other

      map_bucket = new int[map_nbucket];
      for (int i = 0; i < map_nbucket; i++) map_bucket[i] = -1;

      map_hash = new HashElem[map_nhash];
      map_nused = 0;
      map_free = 0;
      for (int i = 0; i < map_nhash; i++) map_hash[i].next = i + 1;
      map_hash[map_nhash - 1].next = -1;

      // use "view" template method to avoid unnecessary deep_copy

      auto h_map_hash = k_map_hash.view<LMPHostType>();    // get type
      h_map_hash = decltype(h_map_hash)(map_nhash);
      k_map_hash.view<LMPHostType>() = h_map_hash;
    }
  }

  k_sametag.modify_host();
  if (map_style == Atom::MAP_ARRAY) k_map_array.modify_host();
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

void AtomKokkos::map_set()
{
  int nall = nlocal + nghost;

  atomKK->sync(Host, TAG_MASK);

  k_sametag.sync_host();
  if (map_style == Atom::MAP_ARRAY) k_map_array.sync_host();

  if (map_style == MAP_ARRAY) {

    // possible reallocation of sametag must come before loop over atoms
    // since loop sets sametag

    if (nall > max_same) {
      max_same = nall + EXTRA;
      memoryKK->destroy_kokkos(k_sametag, sametag);
      memoryKK->create_kokkos(k_sametag, sametag, max_same, "atom:sametag");
    }

    for (int i = nall - 1; i >= 0; i--) {
      sametag[i] = map_array[tag[i]];
      map_array[tag[i]] = i;
    }

  } else {

    // if this proc has more atoms than hash table size, call map_init()
    //   call with 0 since max atomID in system has not changed
    // possible reallocation of sametag must come after map_init(),
    //   b/c map_init() may invoke map_delete(), whacking sametag

    if (nall > map_nhash) map_init(0);
    if (nall > max_same) {
      max_same = nall + EXTRA;
      memoryKK->destroy_kokkos(k_sametag, sametag);
      memoryKK->create_kokkos(k_sametag, sametag, max_same, "atom:sametag");
    }

    int previous, ibucket, index;
    tagint global;

    for (int i = nall - 1; i >= 0; i--) {
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
      if (previous == -1)
        map_bucket[ibucket] = index;
      else
        map_hash[previous].next = index;
      map_hash[index].global = global;
      map_hash[index].local = i;
      map_hash[index].next = -1;
      map_nused++;
    }

    // Copy to Kokkos hash

    // use "view" template method to avoid unnecessary deep_copy

    auto h_map_hash = k_map_hash.view<LMPHostType>();
    h_map_hash.clear();

    for (int i = nall - 1; i >= 0; i--) {

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

      int local = map_hash[index].local;

      auto insert_result = h_map_hash.insert(global, local);
      if (insert_result.failed()) error->one(FLERR, "Kokkos::UnorderedMap insertion failed");
    }
  }

  k_sametag.modify_host();
  if (map_style == Atom::MAP_ARRAY)
    k_map_array.modify_host();
  else if (map_style == Atom::MAP_HASH) {

    // use "view" template method to avoid unnecessary deep_copy

    auto h_map_hash = k_map_hash.view<LMPHostType>();
    auto d_map_hash = k_map_hash.view<LMPDeviceType>();

    // check if fix shake or neigh bond needs a device hash

    int device_hash_flag = 0;

    auto neighborKK = (NeighborKokkos *) neighbor;
    if (neighborKK->device_flag) device_hash_flag = 1;

    for (int n = 0; n < modify->nfix; n++)
      if (utils::strmatch(modify->fix[n]->style, "^shake"))
        if (modify->fix[n]->execution_space == Device) device_hash_flag = 1;

    if (device_hash_flag) {
      Kokkos::deep_copy(d_map_hash, h_map_hash);
      k_map_hash.view<LMPDeviceType>() = d_map_hash;
    }
  }
}

/* ----------------------------------------------------------------------
   free the array or hash table for global to local mapping
------------------------------------------------------------------------- */

void AtomKokkos::map_delete()
{
  memoryKK->destroy_kokkos(k_sametag, sametag);
  sametag = nullptr;

  if (map_style == MAP_ARRAY) {
    memoryKK->destroy_kokkos(k_map_array, map_array);
    map_array = nullptr;
  } else {
    k_map_hash.h_view = host_hash_type();
    k_map_hash.d_view = hash_type();
  }

  Atom::map_delete();
}
