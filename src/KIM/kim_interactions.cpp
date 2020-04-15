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

/* ----------------------------------------------------------------------
   Contributing authors: Axel Kohlmeyer (Temple U),
                         Ryan S. Elliott (UMN)
                         Ellad B. Tadmor (UMN)
                         Ronald Miller   (Carleton U)
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
#include <cstring>
#include <cstdio>
#include <string>
#include <sstream>
#include <vector>
#include "error.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "update.h"
#include "universe.h"
#include "input.h"
#include "variable.h"
#include "fix_store_kim.h"

extern "C" {
#include "KIM_SimulatorHeaders.h"
}

#define SNUM(x)                                                \
  static_cast<std::ostringstream const &>(std::ostringstream() \
                                          << std::dec << x).str()

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

void KimInteractions::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal kim_interactions command");

  if (!domain->box_exist)
    error->all(FLERR,"Must use 'kim_interactions' command after "
                     "simulation box is defined");
  do_setup(narg,arg);
}

/* ---------------------------------------------------------------------- */

void KimInteractions::kim_interactions_log_delimiter(
    std::string const begin_end) const
{
  if (comm->me == 0) {
    std::string mesg;
    if (begin_end == "begin")
      mesg =
          "#=== BEGIN kim_interactions ==================================\n";
    else if (begin_end == "end")
      mesg =
          "#=== END kim_interactions ====================================\n\n";

    input->write_echo(mesg.c_str());
  }
}

/* ---------------------------------------------------------------------- */

void KimInteractions::do_setup(int narg, char **arg)
{
  bool fixed_types;
  if ((narg == 1) && (0 == strcmp("fixed_types",arg[0]))) {
    fixed_types = true;
  }
  else if (narg != atom->ntypes) {
    error->all(FLERR,"Illegal kim_interactions command");
  }
  else {
    fixed_types = false;
  }

  char *model_name = NULL;
  KIM_SimulatorModel *simulatorModel(NULL);

  // check if we had a kim_init command by finding fix STORE/KIM
  // retrieve model name and pointer to simulator model class instance.
  // validate model name if not given as NULL.

  int ifix = modify->find_fix("KIM_MODEL_STORE");
  if (ifix >= 0) {
    FixStoreKIM *fix_store = (FixStoreKIM *) modify->fix[ifix];
    model_name = (char *)fix_store->getptr("model_name");
    simulatorModel = (KIM_SimulatorModel *)fix_store->getptr("simulator_model");
  } else error->all(FLERR,"Must use 'kim_init' before 'kim_interactions'");

  // Begin output to log file
  kim_interactions_log_delimiter("begin");

  if (simulatorModel) {

    if (!fixed_types) {
      std::string delimiter("");
      std::string atom_type_sym_list;
      std::string atom_type_num_list;

      for (int i = 0; i < narg; i++)
      {
        atom_type_sym_list += delimiter + arg[i];
        atom_type_num_list += delimiter + SNUM(species_to_atomic_no(arg[i]));
        delimiter = " ";
      }

      KIM_SimulatorModel_AddTemplateMap(
          simulatorModel,"atom-type-sym-list",atom_type_sym_list.c_str());
      KIM_SimulatorModel_AddTemplateMap(
          simulatorModel,"atom-type-num-list",atom_type_num_list.c_str());
      KIM_SimulatorModel_CloseTemplateMap(simulatorModel);

      int len = strlen(atom_type_sym_list.c_str())+1;
      char *strbuf = new char[len];
      char *strword;

      // validate species selection

      int sim_num_species;
      bool species_is_supported;
      char const *sim_species;
      KIM_SimulatorModel_GetNumberOfSupportedSpecies(
          simulatorModel,&sim_num_species);
      strcpy(strbuf,atom_type_sym_list.c_str());
      strword = strtok(strbuf," \t");
      while (strword) {
        species_is_supported = false;
        if (strcmp(strword,"NULL") == 0) continue;
        for (int i=0; i < sim_num_species; ++i) {
          KIM_SimulatorModel_GetSupportedSpecies(simulatorModel,i,&sim_species);
          if (strcmp(sim_species,strword) == 0)
            species_is_supported = true;
        }
        if (!species_is_supported) {
          std::string msg("Species '");
          msg += strword;
          msg += "' is not supported by this KIM Simulator Model";
          error->all(FLERR,msg.c_str());
        }
        strword = strtok(NULL," \t");
      }
      delete[] strbuf;
    }
    else
    {
      KIM_SimulatorModel_CloseTemplateMap(simulatorModel);
    }

    // check if units are unchanged

    int sim_fields, sim_lines;
    const char *sim_field, *sim_value;
    KIM_SimulatorModel_GetNumberOfSimulatorFields(simulatorModel, &sim_fields);
    for (int i=0; i < sim_fields; ++i) {
      KIM_SimulatorModel_GetSimulatorFieldMetadata(
          simulatorModel,i,&sim_lines,&sim_field);

      if (0 == strcmp(sim_field,"units")) {
        KIM_SimulatorModel_GetSimulatorFieldLine(simulatorModel,i,0,&sim_value);
        if (0 != strcmp(sim_value,update->unit_style))
          error->all(FLERR,"Incompatible units for KIM Simulator Model");
      }
    }

    int sim_model_idx=-1;
    for (int i=0; i < sim_fields; ++i) {
      KIM_SimulatorModel_GetSimulatorFieldMetadata(
          simulatorModel,i,&sim_lines,&sim_field);
      if (0 == strcmp(sim_field,"model-defn")) {
        sim_model_idx = i;
        for (int j=0; j < sim_lines; ++j) {
          KIM_SimulatorModel_GetSimulatorFieldLine(
              simulatorModel,sim_model_idx,j,&sim_value);
	  char strbuf[MAXLINE];
	  char * strword;
	  strcpy(strbuf,sim_value);
	  strword = strtok(strbuf," \t");
	  if (0==strcmp(strword,"KIM_MATCH_PAIRS")) {
            // Notes regarding the KIM_MATCH_PAIRS command
            //  * This is an INTERNAL command.
            //  * It is intended for use only by KIM Simulator Models.
            //  * It is not possible to use this command outside of the context
            //    of the kim_interactions command and KIM Simulator Models.
            //  * The command performs a transformation from symbolic
            //    string-based atom types to lammps numeric atom types for
            //    the pair_coeff settings.
            //  * The command is not documented fully as it is expected to be
            //    temporary.  Eventually it should be replaced by a more
            //    comprehensive symbolic types support in lammps.
	    KIM_MATCH_PAIRS(sim_value);
	  } else {
            input->one(sim_value);
          }
	}
      }
    }

    if (sim_model_idx < 0)
      error->all(FLERR,"KIM Simulator Model has no Model definition");

    KIM_SimulatorModel_OpenAndInitializeTemplateMap(simulatorModel);

  } else {

    // not a simulator model. issue pair_style and pair_coeff commands.

    if (fixed_types)
      error->all(FLERR,"fixed_types cannot be used with a KIM Portable Model");

    // NOTE: all references to arg must appear before calls to input->one()
    // as that will reset the argument vector.

    std::string cmd1("pair_style kim ");
    cmd1 += model_name;

    std::string cmd2("pair_coeff * * ");
    for (int i=0; i < narg; ++i) {
      cmd2 += arg[i];
      cmd2 += " ";
    }

    input->one(cmd1.c_str());
    input->one(cmd2.c_str());
  }

  // End output to log file
  kim_interactions_log_delimiter("end");

}

/* ---------------------------------------------------------------------- */

void KimInteractions::KIM_MATCH_PAIRS(char const *const input_line) const
{
  char strbuf[MAXLINE];
  strcpy(strbuf,input_line);
  char *cmd, *filename;
  cmd = strtok(strbuf," \t");
  filename = strtok(NULL," \t");

  FILE *fp;
  fp = fopen(filename,"r");
  if (fp == NULL) {
    error->one(FLERR,"Parameter file not found");
  }

  std::vector<char *> species;
  for (int i = 0; i < atom->ntypes; ++i)
  {
    char *str;
    str = strtok(NULL," \t");
    if (str == NULL)
      error->one(FLERR,"Incorrect args for pair_species_coeff command");
    species.push_back(str);
  }

  char line[MAXLINE],*ptr;
  int n, eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    char *species1, *species2, *the_rest;
    ptr = line;
    species1 = strtok(ptr," \t");
    species2 = strtok(NULL," \t");
    the_rest = strtok(NULL,"\n");

    for (int type_a = 0; type_a < atom->ntypes; ++type_a) {
      for (int type_b = type_a; type_b < atom->ntypes; ++type_b) {
	if(((strcmp(species[type_a],species1) == 0) &&
            (strcmp(species[type_b],species2) == 0))
           ||
	   ((strcmp(species[type_b],species1) == 0) &&
            (strcmp(species[type_a],species2) == 0))
          ) {
          char pair_command[MAXLINE];
	  sprintf(pair_command,"pair_coeff %i %i %s",type_a+1,type_b+1,
                  the_rest);
	  input->one(pair_command);
	}
      }
    }
  }
  fclose(fp);
}

/* ---------------------------------------------------------------------- */

int KimInteractions::species_to_atomic_no(std::string const species) const
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
