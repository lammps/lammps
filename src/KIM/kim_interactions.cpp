// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Axel Kohlmeyer  (Temple U),
                         Ryan S. Elliott (UMN),
                         Ellad B. Tadmor (UMN),
                         Ronald Miller   (Carleton U),
                         Yaser Afshar    (UMN)
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
   Designed for use with the kim-api-2.1.0 (and newer) package
------------------------------------------------------------------------- */

#include "kim_interactions.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_kim.h"
#include "input.h"
#include "variable.h"
#include "modify.h"
#include "update.h"

#include <cstring>
#include <vector>

extern "C" {
#include "KIM_SimulatorHeaders.h"
}

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

void KimInteractions::command(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "kim interactions", error);

  if (!domain->box_exist)
    error->all(FLERR, "Use of 'kim interactions' before simulation box is defined");

  do_setup(narg, arg);
}

/* ---------------------------------------------------------------------- */

void KimInteractions::do_setup(int narg, char **arg)
{
  bool fixed_types;
  const std::string arg_str(arg[0]);
  if ((narg == 1) && (arg_str == "fixed_types")) {
    fixed_types = true;
  } else if (narg != atom->ntypes) {
    error->all(FLERR, "Illegal 'kim interactions' command.\nThe LAMMPS simulation has {} atom "
               "type(s), but {} chemical species passed to the 'kim interactions' command",
               atom->ntypes, narg);
  } else {
    fixed_types = false;
  }

  char *model_name = nullptr;
  KIM_SimulatorModel *simulatorModel(nullptr);

  // check if we had a kim init command by finding fix STORE/KIM
  // retrieve model name and pointer to simulator model class instance.
  // validate model name if not given as null pointer.

  auto fix_store = dynamic_cast<FixStoreKIM *>(modify->get_fix_by_id("KIM_MODEL_STORE"));
  if (fix_store) {
    model_name = (char *)fix_store->getptr("model_name");
    simulatorModel = (KIM_SimulatorModel *)fix_store->getptr("simulator_model");
  } else error->all(FLERR, "Must use 'kim init' before 'kim interactions'");

  // Begin output to log file
  input->write_echo("#=== BEGIN kim interactions ==================================\n");

  if (simulatorModel) {
    auto first_visit = input->variable->find("kim_update");
    if (!fixed_types) {
      std::string atom_type_sym_list =
        fmt::format("{}", fmt::join(arg, arg + narg, " "));

      std::string atom_type_num_list =
        fmt::format("{}", species_to_atomic_no(arg[0]));

      for (int i = 1; i < narg; ++i)
        atom_type_num_list += fmt::format(" {}", species_to_atomic_no(arg[i]));

      KIM_SimulatorModel_AddTemplateMap(
          simulatorModel, "atom-type-sym-list", atom_type_sym_list.c_str());
      KIM_SimulatorModel_AddTemplateMap(
          simulatorModel, "atom-type-num-list", atom_type_num_list.c_str());
      KIM_SimulatorModel_CloseTemplateMap(simulatorModel);

      // validate species selection

      int sim_num_species;
      bool species_is_supported;
      char const *sim_species;
      KIM_SimulatorModel_GetNumberOfSupportedSpecies(
          simulatorModel, &sim_num_species);

      for (auto atom_type_sym : utils::split_words(atom_type_sym_list)) {
        species_is_supported = false;
        if (atom_type_sym == "NULL") continue;
        for (int i = 0; i < sim_num_species; ++i) {
          KIM_SimulatorModel_GetSupportedSpecies(
            simulatorModel, i, &sim_species);
          if (atom_type_sym == sim_species) species_is_supported = true;
        }
        if (!species_is_supported) {
          error->all(FLERR, "Species '{}' is not supported by this KIM Simulator Model",
                     atom_type_sym);
        }
      }
    } else {
      KIM_SimulatorModel_CloseTemplateMap(simulatorModel);
    }

    // check if units are unchanged

    int sim_fields, sim_lines;
    const char *sim_field, *sim_value;
    KIM_SimulatorModel_GetNumberOfSimulatorFields(simulatorModel, &sim_fields);
    for (int i = 0; i < sim_fields; ++i) {
      KIM_SimulatorModel_GetSimulatorFieldMetadata(simulatorModel, i, &sim_lines, &sim_field);

      const std::string sim_field_str(sim_field);
      if (sim_field_str == "units") {
        KIM_SimulatorModel_GetSimulatorFieldLine(simulatorModel, i, 0, &sim_value);

        if (strcmp(sim_value, update->unit_style) != 0)
          error->all(FLERR, "Incompatible units for KIM Simulator Model: {} vs {}",
                     sim_value, update->unit_style);
      }
    }

    bool no_model_definition = true;
    for (int i = 0; i < sim_fields; ++i) {
      KIM_SimulatorModel_GetSimulatorFieldMetadata(
        simulatorModel, i, &sim_lines, &sim_field);

      const std::string sim_field_str(sim_field);
      if (sim_field_str == "model-defn") {
        if (first_visit < 0) input->one("variable kim_update equal 0");
        else  input->one("variable kim_update equal 1");
        if (domain->periodicity[0] &&
            domain->periodicity[1] &&
            domain->periodicity[2])
          input->one("variable kim_periodic equal 1");
        else if (!domain->periodicity[0] &&
                 !domain->periodicity[1] &&
                 !domain->periodicity[2])
          input->one("variable kim_periodic equal 0");
        else input->one("variable kim_periodic equal 2");

        // KIM Simulator Model has a Model definition
        no_model_definition = false;

        for (int j = 0; j < sim_lines; ++j) {
          KIM_SimulatorModel_GetSimulatorFieldLine(simulatorModel, i, j, &sim_value);
          if (utils::strmatch(sim_value, "^KIM_SET_TYPE_PARAMETERS")) {
            // Notes regarding the KIM_SET_TYPE_PARAMETERS command
            //  * This is an INTERNAL command.
            //  * It is intended for use only by KIM Simulator Models.
            //  * It is not possible to use this command outside of the context
            //    of the kim interactions command and KIM Simulator Models.
            //  * The command performs a transformation from symbolic
            //    string-based atom types to lammps numeric atom types for
            //    the pair_coeff and charge settings.
            //  * The command is not documented fully as it is expected to be
            //    temporary.  Eventually it should be replaced by a more
            //    comprehensive symbolic types support in lammps.
            KIM_SET_TYPE_PARAMETERS(sim_value);
          } else {
            input->one(sim_value);
          }
        }
      }
    }

    if (no_model_definition)
      error->all(FLERR, "KIM Simulator Model has no Model definition");

    KIM_SimulatorModel_OpenAndInitializeTemplateMap(simulatorModel);

  } else {


    // not a simulator model. issue pair_style and pair_coeff commands.

    if (fixed_types)
      error->all(FLERR, "fixed_types cannot be used with a KIM Portable Model");

    // NOTE: all references to arg must appear before calls to input->one()
    // as that will reset the argument vector.

    auto cmd1 = fmt::format("pair_style kim {}", model_name);
    auto cmd2 =
      fmt::format("pair_coeff * * {}", fmt::join(arg, arg + narg, " "));

    input->one(cmd1);
    input->one(cmd2);
  }

  // End output to log file
  input->write_echo("#=== END kim interactions ====================================\n\n");
}

/* ---------------------------------------------------------------------- */

void KimInteractions::KIM_SET_TYPE_PARAMETERS(const std::string &input_line) const
{
  auto words = utils::split_words(input_line);

  const std::string key = words[1];
  if (key != "pair" && key != "charge")
    error->one(FLERR, "Unrecognized KEY {} for KIM_SET_TYPE_PARAMETERS command", key);

  std::string filename = words[2];
  std::vector<std::string> species(words.begin() + 3, words.end());
  if ((int)species.size() != atom->ntypes)
    error->one(FLERR, "Incorrect args for KIM_SET_TYPE_PARAMETERS command");

  FILE *fp = nullptr;
  if (comm->me == 0) {
    fp = fopen(filename.c_str(), "r");
    if (fp == nullptr) error->one(FLERR, "Parameter file {} not found", filename);
  }

  char line[MAXLINE], *ptr;
  int n, eof = 0;

  while (true) {
    if (comm->me == 0) {
      ptr = fgets(line, MAXLINE,fp);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof, 1, MPI_INT, 0, world);
    if (eof) break;
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    MPI_Bcast(line, n, MPI_CHAR, 0, world);

    auto trimmed = utils::trim_comment(line);
    if (trimmed.find_first_not_of(" \t\n\r") == std::string::npos) continue;

    words = utils::split_words(trimmed);
    if (key == "pair") {
      for (int ia = 0; ia < atom->ntypes; ++ia) {
        for (int ib = ia; ib < atom->ntypes; ++ib)
          if (((species[ia] == words[0]) && (species[ib] == words[1]))
              || ((species[ib] == words[0]) && (species[ia] == words[1])))
            input->one(fmt::format("pair_coeff {} {} {}", ia + 1, ib + 1,
              fmt::join(words.begin() + 2, words.end(), " ")));
      }
    } else {
      for (int ia = 0; ia < atom->ntypes; ++ia)
        if (species[ia] == words[0])
          input->one(fmt::format("set type {} charge {}", ia + 1, words[1]));
    }
  }
}

/* ---------------------------------------------------------------------- */

int KimInteractions::species_to_atomic_no(const std::string &species) const
{
  if (species == "H") return 1;
  else if (species == "He") return 2;
  else if (species == "Li") return 3;
  else if (species == "Be") return 4;
  else if (species == "B") return 5;
  else if (species == "C") return 6;
  else if (species == "N") return 7;
  else if (species == "O") return 8;
  else if (species == "F") return 9;
  else if (species == "Ne") return 10;
  else if (species == "Na") return 11;
  else if (species == "Mg") return 12;
  else if (species == "Al") return 13;
  else if (species == "Si") return 14;
  else if (species == "P") return 15;
  else if (species == "S") return 16;
  else if (species == "Cl") return 17;
  else if (species == "Ar") return 18;
  else if (species == "K") return 19;
  else if (species == "Ca") return 20;
  else if (species == "Sc") return 21;
  else if (species == "Ti") return 22;
  else if (species == "V") return 23;
  else if (species == "Cr") return 24;
  else if (species == "Mn") return 25;
  else if (species == "Fe") return 26;
  else if (species == "Co") return 27;
  else if (species == "Ni") return 28;
  else if (species == "Cu") return 29;
  else if (species == "Zn") return 30;
  else if (species == "Ga") return 31;
  else if (species == "Ge") return 32;
  else if (species == "As") return 33;
  else if (species == "Se") return 34;
  else if (species == "Br") return 35;
  else if (species == "Kr") return 36;
  else if (species == "Rb") return 37;
  else if (species == "Sr") return 38;
  else if (species == "Y") return 39;
  else if (species == "Zr") return 40;
  else if (species == "Nb") return 41;
  else if (species == "Mo") return 42;
  else if (species == "Tc") return 43;
  else if (species == "Ru") return 44;
  else if (species == "Rh") return 45;
  else if (species == "Pd") return 46;
  else if (species == "Ag") return 47;
  else if (species == "Cd") return 48;
  else if (species == "In") return 49;
  else if (species == "Sn") return 50;
  else if (species == "Sb") return 51;
  else if (species == "Te") return 52;
  else if (species == "I") return 53;
  else if (species == "Xe") return 54;
  else if (species == "Cs") return 55;
  else if (species == "Ba") return 56;
  else if (species == "La") return 57;
  else if (species == "Ce") return 58;
  else if (species == "Pr") return 59;
  else if (species == "Nd") return 60;
  else if (species == "Pm") return 61;
  else if (species == "Sm") return 62;
  else if (species == "Eu") return 63;
  else if (species == "Gd") return 64;
  else if (species == "Tb") return 65;
  else if (species == "Dy") return 66;
  else if (species == "Ho") return 67;
  else if (species == "Er") return 68;
  else if (species == "Tm") return 69;
  else if (species == "Yb") return 70;
  else if (species == "Lu") return 71;
  else if (species == "Hf") return 72;
  else if (species == "Ta") return 73;
  else if (species == "W") return 74;
  else if (species == "Re") return 75;
  else if (species == "Os") return 76;
  else if (species == "Ir") return 77;
  else if (species == "Pt") return 78;
  else if (species == "Au") return 79;
  else if (species == "Hg") return 80;
  else if (species == "Tl") return 81;
  else if (species == "Pb") return 82;
  else if (species == "Bi") return 83;
  else if (species == "Po") return 84;
  else if (species == "At") return 85;
  else if (species == "Rn") return 86;
  else if (species == "Fr") return 87;
  else if (species == "Ra") return 88;
  else if (species == "Ac") return 89;
  else if (species == "Th") return 90;
  else if (species == "Pa") return 91;
  else if (species == "U") return 92;
  else if (species == "Np") return 93;
  else if (species == "Pu") return 94;
  else if (species == "Am") return 95;
  else if (species == "Cm") return 96;
  else if (species == "Bk") return 97;
  else if (species == "Cf") return 98;
  else if (species == "Es") return 99;
  else if (species == "Fm") return 100;
  else if (species == "Md") return 101;
  else if (species == "No") return 102;
  else if (species == "Lr") return 103;
  else if (species == "Rf") return 104;
  else if (species == "Db") return 105;
  else if (species == "Sg") return 106;
  else if (species == "Bh") return 107;
  else if (species == "Hs") return 108;
  else if (species == "Mt") return 109;
  else if (species == "Ds") return 110;
  else if (species == "Rg") return 111;
  else if (species == "Cn") return 112;
  else if (species == "Nh") return 113;
  else if (species == "Fl") return 114;
  else if (species == "Mc") return 115;
  else if (species == "Lv") return 116;
  else if (species == "Ts") return 117;
  else if (species == "Og") return 118;
  else return -1;
}
