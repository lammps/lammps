/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MERGESORT
#define LMP_MERGESORT

/* ---------------------------------------------------------------------- */

// custom upward merge sort implementation which allows to pass a custom
// pointer to the comparison function for access to class instances.
// this avoids having to use global variables.

// part 1. merge two sublists.

static void do_merge(int *idx, int *buf, int llo, int lhi, int rlo, int rhi,
                     void *ptr, int (*comp)(int, int, void *))
{
  int i = llo;
  int l = llo;
  int r = rlo;
  while ((l < lhi) && (r < rhi)) {
    if ((*comp)(buf[l],buf[r],ptr) < 0)
      idx[i++] = buf[l++];
    else idx[i++] = buf[r++];
  }
    
  while(l < lhi) idx[i++] = buf[l++];
  while(r < rhi) idx[i++] = buf[r++];
}

// part 2: loop over sublists doubling in size with each iteration

static void merge_sort(int *index, int num, void *ptr,
                       int (*comp)(int, int, void *))
{
  if (num < 2) return;

  int *hold = new int[num];
  int i,j,k,m;

  i = 1;
  while (i < num) {
    memcpy(hold,index,sizeof(int)*num);
    for (j=0; j < num-1; j += 2*i) {
      k = j + 2*i;
      if (k > num) k=num;
      m = j+i;
      if (m > num) m=num;
      do_merge(index,hold,j,m,m,k,ptr,comp);
    }
    i *= 2;
  }

  delete[] hold;
}

#endif
