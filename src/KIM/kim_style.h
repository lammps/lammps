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

/* ----------------------------------------------------------------------
   Contributing authors: Axel Kohlmeyer (Temple U),
                         Ryan S. Elliott (UMN)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program; if not, see <https://www.gnu.org/licenses>.

   Linking LAMMPS statically or dynamically with other modules is making a
   combined work based on LAMMPS. Thus, the terms and conditions of the GNU
   General Public License cover the whole combination.

   In addition, as a special exception, the copyright holders of LAMMPS give
   you permission to combine LAMMPS with free software programs or libraries
   that are released under the GNU LGPL and with code included in the standard
   release of the "kim-api" under the CDDL (or modified versions of such code,
   with unchanged license). You may copy and distribute such a system following
   the terms of the GNU GPL for LAMMPS and the licenses of the other code
   concerned, provided that you include the source code of that other code
   when and as the GNU GPL requires distribution of source code.

   Note that people who make modified versions of LAMMPS are not obligated to
   grant this special exception for their modified versions; it is their choice
   whether to do so. The GNU General Public License gives permission to release
   a modified version without this exception; this exception also makes it
   possible to release a modified version which carries forward this exception.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Designed for use with the kim-api-2.0.2 (and newer) package
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(kim_style,KimStyle)

#else

#ifndef LMP_KIM_STYLE_H
#define LMP_KIM_STYLE_H

#include "pointers.h"

namespace LAMMPS_NS {

class KimStyle : protected Pointers {
 public:
  KimStyle(class LAMMPS *lmp) : Pointers(lmp) {};
  void command(int, char **);
 private:
  void do_init(char *);
  void do_defn(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal kim_style command

Incorrect number or kind of arguments to kim_style.

E: Must use 'kim_style init' command before simulation box is defined

Self-explanatory.

E: Must use 'kim_style define' command after simulation box is defined

Self-explanatory.

E: Must use 'kim_style init' command before 'kim_style define'

Self-explanatory.

E: Incompatible KIM Simulator Model

The requested KIM Simulator Model was defined for a different MD code
and thus is not compatible with LAMMPS.

E: Species XXX is not supported by this KIM Simulator Model

The kim_style define command was referencing a species that is not
present in the requested KIM Simulator Model.

E: Incompatible units for KIM Simulator Model

The selected unit style is not compatible with the requested KIM
Simulator Model.

E: KIM Simulator Model has no Model definition

There is no model definition (key: model-defn) in the KIM Simulator
Model.  Please contact the OpenKIM database maintainers to verify
and potentially correct this.

*/
