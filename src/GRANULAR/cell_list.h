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

#ifndef LMP_CELL_LIST_H
#define LMP_CELL_LIST_H

#include <vector>
#include <mpi.h>
#include "pointers.h"
#include "near_list.h"

namespace LAMMPS_NS {

class CellList : protected Pointers, public virtual INearList {
  friend class DistributedCellList;

  int nbinx, nbiny, nbinz;         // # of global bins
  int mbinx, mbiny, mbinz;
  int mbinxlo, mbinylo, mbinzlo;

  double binsizex, binsizey, binsizez;  // bin sizes and inverse sizes
  double bininvx, bininvy, bininvz;
  double bboxlo[3];
  double bboxhi[3];

protected:
  struct Element {
    double x[3];
    double r;

    Element() {
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 0.0;
        r = 0.0;
    }

    Element(const Element & other) {
        x[0] = other.x[0];
        x[1] = other.x[1];
        x[2] = other.x[2];
        r = other.r;
    }

    Element(double * x, double r) {
        this->x[0] = x[0];
        this->x[1] = x[1];
        this->x[2] = x[2];
        this->r = r;
    }
  };

  std::vector<int> binhead;
  std::vector<int> next;
  std::vector<Element> elements;
  int nbins;

  const static int NSTENCIL = 27;
  int stencil[NSTENCIL];

  int coord2bin(double *x) const;

public:
  CellList(class LAMMPS *);
  virtual ~CellList() = default;

  void setup(double * bboxlo, double * bboxhi, double binsize);
  void clear();

  // Implementation of INearList interface
  virtual void insert(double * x, double r);
  virtual size_t count() const;
  virtual bool has_overlap(double *x, double r) const;
};


class DistributedCellList : public CellList, public IDistributedNearList {
  std::vector<int> recvcounts;
  std::vector<int> displs;
  MPI_Datatype mpi_element_type;

public:
  DistributedCellList(LAMMPS * lmp);
  virtual ~DistributedCellList() = default;

  // Implementation of IDistributedNearList
  virtual void allgather(INearList * local_nlist);
};

}

#endif
