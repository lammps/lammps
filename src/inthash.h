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

#ifndef INTHASH_INCLUDED
#define INTHASH_INCLUDED

/* re-usable integer hash table code with static linkage. */

/** hash table top level data structure */
typedef struct inthash_t {
  struct inthash_node_t **bucket;        /* array of hash nodes */
  int size;                           /* size of the array */
  int entries;                        /* number of entries in table */
  int downshift;                      /* shift cound, used in hash function */
  int mask;                           /* used to select bits for hashing */
} inthash_t;

/** hash table node data structure */
typedef struct inthash_node_t {
  int data;                           /* data in hash node */
  int key;                            /* key for hash lookup */
  struct inthash_node_t *next;        /* next node in hash chain */
} inthash_node_t;

#define HASH_FAIL  -1
#define HASH_LIMIT  0.5

/* initialize new hash table  */
static void inthash_init(inthash_t *tptr, int buckets);
/* lookup entry in hash table */
static int inthash_lookup(const inthash_t *tptr, int key);
/* generate list of keys for reverse lookups. */
static int *inthash_keys(inthash_t *tptr);
/* insert an entry into hash table. */
static int inthash_insert(inthash_t *tptr, int key, int data);
/* delete the hash table */
static void inthash_destroy(inthash_t *tptr);
/* adapted sort for in-place sorting of map indices. */
static void id_sort(int *idmap, int left, int right);

/************************************************************************
 * integer hash code:
 ************************************************************************/

/* inthash() - Hash function returns a hash number for a given key.
 * tptr: Pointer to a hash table, key: The key to create a hash number for */
static int inthash(const inthash_t *tptr, int key) {
  int hashvalue;

  hashvalue = (((key*1103515249)>>tptr->downshift) & tptr->mask);
  if (hashvalue < 0) {
    hashvalue = 0;
  }    

  return hashvalue;
}

/*
 *  rebuild_table_int() - Create new hash table when old one fills up.
 *
 *  tptr: Pointer to a hash table
 */
static void rebuild_table_int(inthash_t *tptr) {
  inthash_node_t **old_bucket, *old_hash, *tmp;
  int old_size, h, i;

  old_bucket=tptr->bucket;
  old_size=tptr->size;

  /* create a new table and rehash old buckets */
  inthash_init(tptr, old_size<<1);
  for (i=0; i<old_size; i++) {
    old_hash=old_bucket[i];
    while(old_hash) {
      tmp=old_hash;
      old_hash=old_hash->next;
      h=inthash(tptr, tmp->key);
      tmp->next=tptr->bucket[h];
      tptr->bucket[h]=tmp;
      tptr->entries++;
    } /* while */
  } /* for */

  /* free memory used by old table */
  free(old_bucket);

  return;
}

/*
 *  inthash_init() - Initialize a new hash table.
 *
 *  tptr: Pointer to the hash table to initialize
 *  buckets: The number of initial buckets to create
 */
void inthash_init(inthash_t *tptr, int buckets) {

  /* make sure we allocate something */
  if (buckets==0)
    buckets=16;

  /* initialize the table */
  tptr->entries=0;
  tptr->size=2;
  tptr->mask=1;
  tptr->downshift=29;

  /* ensure buckets is a power of 2 */
  while (tptr->size<buckets) {
    tptr->size<<=1;
    tptr->mask=(tptr->mask<<1)+1;
    tptr->downshift--;
  } /* while */

  /* allocate memory for table */
  tptr->bucket=(inthash_node_t **) calloc(tptr->size, sizeof(inthash_node_t *));

  return;
}

/*
 *  inthash_lookup() - Lookup an entry in the hash table and return a pointer to
 *    it or HASH_FAIL if it wasn't found.
 *
 *  tptr: Pointer to the hash table
 *  key: The key to lookup
 */
int inthash_lookup(const inthash_t *tptr, int key) {
  int h;
  inthash_node_t *node;


  /* find the entry in the hash table */
  h=inthash(tptr, key);
  for (node=tptr->bucket[h]; node!=NULL; node=node->next) {
    if (node->key == key)
      break;
  }

  /* return the entry if it exists, or HASH_FAIL */
  return(node ? node->data : HASH_FAIL);
}


/*
 *  inthash_keys() - Return a list of keys.
 *  NOTE: the returned list must be freed with free(3).
 */
int *inthash_keys(inthash_t *tptr) {

  int *keys;
  inthash_node_t *node;
  
  keys = (int *)calloc(tptr->entries, sizeof(int));
  
  for (int i=0; i < tptr->size; ++i) {
    for (node=tptr->bucket[i]; node != NULL; node=node->next) {
      keys[node->data] = node->key;
    }
  }
  
  return keys;
}

/*
 *  inthash_insert() - Insert an entry into the hash table.  If the entry already
 *  exists return a pointer to it, otherwise return HASH_FAIL.
 *
 *  tptr: A pointer to the hash table
 *  key: The key to insert into the hash table
 *  data: A pointer to the data to insert into the hash table
 */
int inthash_insert(inthash_t *tptr, int key, int data) {
  int tmp;
  inthash_node_t *node;
  int h;

  /* check to see if the entry exists */
  if ((tmp=inthash_lookup(tptr, key)) != HASH_FAIL)
    return(tmp);

  /* expand the table if needed */
  while (tptr->entries>=HASH_LIMIT*tptr->size)
    rebuild_table_int(tptr);

  /* insert the new entry */
  h=inthash(tptr, key);
  node=(struct inthash_node_t *) malloc(sizeof(inthash_node_t));
  node->data=data;
  node->key=key;
  node->next=tptr->bucket[h];
  tptr->bucket[h]=node;
  tptr->entries++;

  return HASH_FAIL;
}

/*
 * inthash_destroy() - Delete the entire table, and all remaining entries.
 * 
 */
void inthash_destroy(inthash_t *tptr) {
  inthash_node_t *node, *last;
  int i;

  for (i=0; i<tptr->size; i++) {
    node = tptr->bucket[i];
    while (node != NULL) { 
      last = node;   
      node = node->next;
      free(last);
    }
  }     

  /* free the entire array of buckets */
  if (tptr->bucket != NULL) {
    free(tptr->bucket);
    memset(tptr, 0, sizeof(inthash_t));
  }
}

/************************************************************************
 * integer list sort code:
 ************************************************************************/

/* sort for integer map. initial call  id_sort(idmap, 0, natoms - 1); */
static void id_sort(int *idmap, int left, int right)
{
  int pivot, l_hold, r_hold;

  l_hold = left;
  r_hold = right;
  pivot = idmap[left];

  while (left < right) {
    while ((idmap[right] >= pivot) && (left < right))
      right--;
    if (left != right) {
      idmap[left] = idmap[right];
      left++;
    }
    while ((idmap[left] <= pivot) && (left < right))
      left++;
    if (left != right) {
      idmap[right] = idmap[left];
      right--;
    }
  }
  idmap[left] = pivot;
  pivot = left;
  left = l_hold;
  right = r_hold;

  if (left < pivot)
    id_sort(idmap, left, pivot-1);
  if (right > pivot)
    id_sort(idmap, pivot+1, right);
}

#endif /* INTHASH_INCLUDED */

// Local Variables:
// mode: c++
// End:
