/*
Copyright 2021 Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
^2: University of Cambridge, Cambridge, United Kingdom
^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
^4: University of British Columbia, Vancouver, BC, Canada


    This FILENAME is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


//
// Created by Lysogorskiy Yury on 27.02.20.
//


#ifdef PAIR_CLASS

PairStyle(pace,PairPACE)

#else

#ifndef LMP_PAIR_PACE_H
#define LMP_PAIR_PACE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPACE : public Pair {
 public:
  PairPACE(class LAMMPS *);
  virtual ~PairPACE();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);

  void *extract(const char *, int &);

 protected:
  struct ACEImpl *aceimpl;

  virtual void allocate();

  double **scale;
  bool recursive;             // "recursive" option for ACERecursiveEvaluator
};
}

#endif
#endif
