/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*
  Do not link plumed directly but rather do it at runtime
*/
#define __PLUMED_WRAPPER_LINK_RUNTIME 1

/*
  Make sure the inline C++ interface is not included here.
  Should not be necessary, but it doesn't hurt.
*/
#define __PLUMED_WRAPPER_CXX 0

/*
  Tell Plumed.h to emit the whole implementation
*/
#define __PLUMED_WRAPPER_IMPLEMENTATION 1

/*
  Emit fortran wrappers
*/
#define __PLUMED_WRAPPER_FORTRAN 1

#include "Plumed.h"
