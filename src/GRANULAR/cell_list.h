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
#include "pointers.h"

namespace LAMMPS_NS {

class CellList : protected Pointers {
  struct Element {
    double x[3];
    double r;

    Element(double * x, double r) {
        this->x[0] = x[0]; 
        this->x[1] = x[1]; 
        this->x[2] = x[2]; 
        this->r = r;
    }
  };

  int nbinx, nbiny, nbinz;         // # of global bins
  int mbinx, mbiny, mbinz;
  int mbinxlo, mbinylo, mbinzlo;

  double binsizex, binsizey, binsizez;  // bin sizes and inverse sizes
  double bininvx, bininvy, bininvz;
  double bboxlo[3];
  double bboxhi[3];

  std::vector<int> binhead;
  std::vector<int> next;
  std::vector<Element> elements;

  const static int NSTENCIL = 27;
  int stencil[NSTENCIL];

  int coord2bin(double *x) const;

public:
    CellList(class LAMMPS *);

    void setup(double * bboxlo, double * bboxhi, double binsize);

    void insert(double * x, double r);
    void clear();

    size_t count() const;

    bool has_overlap(double *x, double r) const;
};

}

#endif
