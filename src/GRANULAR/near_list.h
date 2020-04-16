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

#ifndef LMP_NEAR_LIST_H
#define LMP_NEAR_LIST_H

#include "pointers.h"

namespace LAMMPS_NS {

/***************************************************************************
 * Interfaces
 ***************************************************************************/

class INearList {
public:
  virtual ~INearList() = default;
  virtual void insert(double * x, double r) = 0;
  virtual size_t count() const = 0;
  virtual bool has_overlap(double *x, double r) const = 0;
};

class IDistributedNearList : public virtual INearList {
public:
  virtual ~IDistributedNearList() = default;
  virtual void allgather(INearList * local_nlist) = 0;
};

/***************************************************************************
 * Implementations
 ***************************************************************************/

/***************************************************************************
 * Original implementation of near list used by FixPour, refactored into
 * a class
 ***************************************************************************/
class NearList : protected Pointers, public virtual INearList {
  friend class DistributedNearList;

  double ** elements;
  int ncount;
  int nallocated;

public:
  NearList(class LAMMPS *);
  virtual ~NearList();

  void allocate(int nmax);

  // Implementation of INearList interface
  virtual void insert(double * x, double r);
  virtual size_t count() const;
  virtual bool has_overlap(double *x, double r) const;
};


/***************************************************************************
 * Extends NearList with allgather method to collect local lists from other
 * processors
 ***************************************************************************/
class DistributedNearList : public NearList, public virtual IDistributedNearList {
  int nprocs;
  int *recvcounts, *displs;

public:
  DistributedNearList(class LAMMPS * lmp);
  virtual ~DistributedNearList();

  // Implementation of IDistributedNearList
  virtual void allgather(INearList * local_nlist);
};

}

#endif
