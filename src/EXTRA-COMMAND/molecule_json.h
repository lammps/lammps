/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(molecule/json,MoleculeJsonCommand);
// clang-format on
#else

#ifndef LMP_MOLECULE_JSON_H
#define LMP_MOLECULE_JSON_H

#include "molecule.h"
#include "command.h"

namespace LAMMPS_NS {

class MoleculeJsonCommand : public Command {
public:
   MoleculeJsonCommand(class LAMMPS * lmp) : Command(lmp) {};;
   void command(int, char **) override;

   class MoleculeJson : public Molecule {
   public:
      MoleculeJson(class LAMMPS *, int, char **, int &);
      ~MoleculeJson();

   private:
      void read(int) override;
      void coords();
      void types();
      void molecules();
      void charges();
      void bonds(int);
      void impropers(int);
      void angles(int);
      void dihedrals(int);

      struct JsonImpl *json_impl;     // pointer to json PIMPL struct
  };
};

}

#endif
#endif
