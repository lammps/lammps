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
   Contributing authors: Yaser Afshar (UMN),
                         Ryan S. Elliott (UMN),
                         Ellad B. Tadmor (UMN)
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

#include "kim_param.h"
#include <mpi.h>
#include <cstring>
#include <string>
#include <sstream>
#include "comm.h"
#include "error.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "force.h"
#include "fix_store_kim.h"
#include "pair_kim.h"

extern "C"
{
#include "KIM_SimulatorHeaders.h"
}

#ifdef SNUM
#undef SNUM
#endif

#define SNUM(x)                                                \
  static_cast<std::ostringstream const &>(std::ostringstream() \
                                          << std::dec << x).str()

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

namespace
{
void get_kim_unit_names(
    char const *const system,
    KIM_LengthUnit &lengthUnit,
    KIM_EnergyUnit &energyUnit,
    KIM_ChargeUnit &chargeUnit,
    KIM_TemperatureUnit &temperatureUnit,
    KIM_TimeUnit &timeUnit,
    Error *error)
{
  if ((strcmp(system, "real") == 0)) {
    lengthUnit = KIM_LENGTH_UNIT_A;
    energyUnit = KIM_ENERGY_UNIT_kcal_mol;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_fs;
  } else if ((strcmp(system, "metal") == 0)) {
    lengthUnit = KIM_LENGTH_UNIT_A;
    energyUnit = KIM_ENERGY_UNIT_eV;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_ps;
  } else if ((strcmp(system, "si") == 0)) {
    lengthUnit = KIM_LENGTH_UNIT_m;
    energyUnit = KIM_ENERGY_UNIT_J;
    chargeUnit = KIM_CHARGE_UNIT_C;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_s;
  } else if ((strcmp(system, "cgs") == 0)) {
    lengthUnit = KIM_LENGTH_UNIT_cm;
    energyUnit = KIM_ENERGY_UNIT_erg;
    chargeUnit = KIM_CHARGE_UNIT_statC;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_s;
  } else if ((strcmp(system, "electron") == 0)) {
    lengthUnit = KIM_LENGTH_UNIT_Bohr;
    energyUnit = KIM_ENERGY_UNIT_Hartree;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_fs;
  } else if ((strcmp(system, "lj") == 0)) {
    error->all(FLERR, "LAMMPS unit_style lj not supported by KIM models");
  } else
    error->all(FLERR, "Unknown unit_style");
}
} // namespace

/* ---------------------------------------------------------------------- */

KimParam::KimParam(LAMMPS *lmp) : Pointers(lmp) {}

KimParam::~KimParam() {}

void KimParam::command(int narg, char **arg)
{
  // kim_param is a command for
  // getting/setting the value of a %KIM PM parameter
  //
  // kim_param get param_name index_range variables formatarg
  // kim_param set param_name index_range values

  // kim_param   get paramname 1 varname
  // kim_param   get paramname index_range varname_1, ..., varname_N
  // kim_param   get paramname index_range varname_base split
  // kim_param   get paramname index_range varname_base list
  // kim_param   set paramname index_range values

  if (narg < 4)
    error->all(FLERR, "Illegal kim_param command");

  kim_param_get = (strcmp(arg[0], "get") == 0);
  kim_param_set = (strcmp(arg[0], "set") == 0);

  if (!kim_param_get && !kim_param_set) {
    std::string msg("Incorrect arguments in kim_param command.\n");
    msg += "'kim_param get/set' is mandatory.";
    error->all(FLERR, msg.c_str());
  }

  // Check if we called a kim_init command
  // by finding fix STORE/KIM
  // retrieve model name and model units.

  char *model_name;
  char *model_units;

  bool isPortableModel(false);

  int const ifix = modify->find_fix("KIM_MODEL_STORE");
  if (ifix >= 0) {
    FixStoreKIM *fix_store = reinterpret_cast<FixStoreKIM *>(modify->fix[ifix]);

    KIM_SimulatorModel *simulatorModel =
        reinterpret_cast<KIM_SimulatorModel *>(
            fix_store->getptr("simulator_model"));

    isPortableModel = simulatorModel ? false : true;
    if (!isPortableModel)
      error->all(FLERR, "kim_param can only be used with a KIM Portable Model");

    model_name = (char *)fix_store->getptr("model_name");
    model_units = (char *)fix_store->getptr("model_units");
  }
  else
    error->all(FLERR, "Must use 'kim_init' before 'kim_param'");

  kim_param_log_delimiter("begin");

  KIM_Model *pkim = NULL;

  std::string atom_type_list;

  int kim_error;

  bool isPairStyleAssigned = force->pair ? true : false;
  if (isPairStyleAssigned) {
    Pair *pair = force->pair_match("kim", 1, 0);
    if (pair) {
      PairKIM *pairKIM = reinterpret_cast<PairKIM *>(pair);

      pkim = pairKIM->get_KIM_Model();
      if (!pkim)
        error->all(FLERR, "Unable to get the KIM Portable Model.");

      if (kim_param_set) {
        atom_type_list = pairKIM->get_atom_type_list();
        if (atom_type_list.empty())
          error->all(FLERR, "The requested atom type list is empty.");
      }
    } else
      error->all(FLERR, "Pair style is defined,"
                        " but there is no match for kim style in lammps.");
  } else {
    if (kim_param_set) {
      std::string msg("Wrong kim_param set command.\n");
      msg += "To set the new parameter values, pair style must be assigned.\n";
      msg += "Must use 'kim_interactions' or";
      msg += "'pair_style kim ' before 'kim_param set'";
      error->all(FLERR, msg.c_str());
    } else {
      KIM_LengthUnit lengthUnit;
      KIM_EnergyUnit energyUnit;
      KIM_ChargeUnit chargeUnit;
      KIM_TemperatureUnit temperatureUnit;
      KIM_TimeUnit timeUnit;

      get_kim_unit_names(model_units, lengthUnit, energyUnit,
                         chargeUnit, temperatureUnit, timeUnit,
                         error);

      int units_accepted;

      kim_error = KIM_Model_Create(KIM_NUMBERING_zeroBased,
                                   lengthUnit,
                                   energyUnit,
                                   chargeUnit,
                                   temperatureUnit,
                                   timeUnit,
                                   model_name,
                                   &units_accepted,
                                   &pkim);
      if (kim_error)
        error->all(FLERR, "Unable to create KIM Portable Model.");
    }
  }

  // Get the number of mutable parameters in the kim model
  int numberOfParameters(0);

  KIM_Model_GetNumberOfParameters(pkim, &numberOfParameters);
  if (numberOfParameters) {
    // Get the parameters
    if (kim_param_get) {
      // Parameter name
      char *paramname = NULL;
      // Variable name
      char *varname = NULL;

      // Loop over all the arguments
      for (int i = 1; i < narg;) {
        // Parameter name
        if (i < narg)
          paramname = arg[i++];
        else
          break;

        // Find the requested parameter within the model parameters
        int param_index;
        KIM_DataType kim_DataType;
        int extent;
        char const *str_name = NULL;
        char const *str_desc = NULL;

        for (param_index = 0; param_index < numberOfParameters; ++param_index) {
          kim_error = KIM_Model_GetParameterMetadata(pkim, param_index,
                                                     &kim_DataType, &extent,
                                                     &str_name, &str_desc);
          if (kim_error)
            error->all(FLERR, "KIM GetParameterMetadata returned error.");

          if (strcmp(paramname, str_name) == 0)
            break;
        }

        if (param_index >= numberOfParameters) {
          std::string msg("Wrong argument in kim_param get command.\n");
          msg += "This Model does not have the requested '";
          msg += paramname;
          msg += "' parameter.";
          error->all(FLERR, msg.c_str());
        }

        // Get the index_range for the requested parameter
        int nlbound(0);
        int nubound(0);

        if (i < narg) {
          std::string argtostr(arg[i++]);

          // Check to see if the indices range contains
          // only integer numbers and/or range :
          if (argtostr.find_first_not_of("0123456789:") != std::string::npos) {
            std::string msg("Illegal index_range.\n");
            msg += "Expected integer parameter(s) instead of '";
            msg += argtostr;
            msg += "' in index_range.";
            error->all(FLERR, msg.c_str());
          }

          std::string::size_type npos = argtostr.find(':');
          if (npos != std::string::npos) {
            argtostr[npos] = ' ';
            std::stringstream str(argtostr);
            str >> nlbound >> nubound;
            if (nubound < 1 || nubound > extent ||
                nlbound < 1 || nlbound > nubound) {
              std::string msg("Illegal index_range '");
              msg += SNUM(nlbound) + "-" + SNUM(nubound);
              msg += "' for '";
              msg += paramname;
              msg += "' parameter with extent of '";
              msg += SNUM(extent);
              msg += "' .";
              error->all(FLERR, msg.c_str());
            }
          } else {
            std::stringstream str(argtostr);
            str >> nlbound;
            if (nlbound < 1 || nlbound > extent) {
              std::string msg("Illegal index '");
              msg += SNUM(nlbound) + "' for parameter '";
              msg += paramname;
              msg += "' with the extent of '";
              msg += SNUM(extent);
              msg += "' .";
              error->all(FLERR, msg.c_str());
            }
            nubound = nlbound;
          }
        } else {
          std::string msg("Wrong number of arguments in ");
          msg += "kim_param get command.\n";
          msg += "Index range after parameter name is mandatory.";
          error->all(FLERR, msg.c_str());
        }

        int const nvars = nubound - nlbound + 1;
        char **varsname = NULL;

        if (i < narg) {
          // Get the variable/variable_base name
          varname = arg[i++];
        } else {
          std::string msg("Wrong number of arguments in ");
          msg += "kim_param get command.\n";
          msg += "The LAMMPS variable name is mandatory.";
          error->all(FLERR, msg.c_str());
        }

        // indicator flag for list request
        bool list_requested(false);

        if (nvars > 1) {
          if (i < narg) {
            if (strcmp(arg[i], "split") == 0) {
              varsname = new char *[nvars];
              for (int j = 0, k = nlbound; j < nvars; ++j, ++k) {
                std::stringstream str;
                str << varname << "_" << k;
                varsname[j] = const_cast<char *>(str.str().c_str());
              }
            } else if (strcmp(arg[i], "list") == 0) {
              list_requested = true;
              varsname = new char *[1];
              varsname[0] = varname;
            // Default explicit (optional) formatarg
            } else if (i - 1 + nvars < narg) {
              varsname = new char *[nvars];
              --i;
              for (int j = 0; j < nvars; ++j, ++i)
                varsname[j] = arg[i];
              if (i < narg) {
                if (strcmp(arg[i], "explicit") == 0)
                  ++i;
              }
            } else {
              std::string msg("Wrong number of arguments in ");
              msg += "kim_param get command.\n";
              msg += "The LAMMPS '";
              msg += SNUM(nvars);
              msg += "' variable names or '";
              msg += varname;
              msg += " split' is mandatory.";
              error->all(FLERR, msg.c_str());
            }
          } else {
            std::string msg("Wrong number of arguments in ");
            msg += "kim_param get command.\n";
            msg += "The LAMMPS '";
            msg += SNUM(nvars);
            msg += "' variable names or '";
            msg += varname;
            msg += " split/list' is mandatory.";
            error->all(FLERR, msg.c_str());
          }
        } else {
          varsname = new char *[1];
          if (i < narg)
          {
            if (strcmp(arg[i], "split") == 0)
            {
              std::stringstream str;
              str << varname << "_" << nlbound;
              varsname[0] = const_cast<char *>(str.str().c_str());
              ++i;
            } else {
              if ((strcmp(arg[i], "list") == 0) ||
                  (strcmp(arg[i], "explicit") == 0))
                ++i;
              varsname[0] = varname;
            }
          } else {
            varsname[0] = varname;
          }
        }

        char **varcmd = new char *[3];
        varcmd[1] = const_cast<char *>("string");

        if (KIM_DataType_Equal(kim_DataType, KIM_DATA_TYPE_Double)) {
          if (list_requested) {
            std::stringstream str;
            double V;
            {
              kim_error = KIM_Model_GetParameterDouble(pkim, param_index,
                                                       nlbound - 1, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterDouble returned error.");
              str << V;
            }
            for (int j = 1; j < nvars; ++j) {
              kim_error = KIM_Model_GetParameterDouble(pkim, param_index,
                                                       nlbound - 1 + j, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterDouble returned error.");
              str << " " << V;
            }
            varcmd[0] = varsname[0];
            varcmd[2] = const_cast<char *>(str.str().c_str());
            input->variable->set(3, varcmd);
            echo_var_assign(varcmd[0], varcmd[2]);
          } else {
            for (int j = 0; j < nvars; ++j) {
              varcmd[0] = varsname[j];
              double V;
              kim_error = KIM_Model_GetParameterDouble(pkim, param_index,
                                                       nlbound - 1 + j, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterDouble returned error.");
              std::stringstream str;
              str << V;
              varcmd[2] = const_cast<char *>(str.str().c_str());
              input->variable->set(3, varcmd);
              echo_var_assign(varcmd[0], varcmd[2]);
            }
          }
        } else if (KIM_DataType_Equal(kim_DataType, KIM_DATA_TYPE_Integer)) {
          if (list_requested) {
            std::stringstream str;
            int V;
            {
              kim_error = KIM_Model_GetParameterInteger(pkim, param_index,
                                                        nlbound - 1, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterInteger returned error.");
              str << V;
            }
            for (int j = 1; j < nvars; ++j) {
              kim_error = KIM_Model_GetParameterInteger(pkim, param_index,
                                                        nlbound - 1 + j, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterInteger returned error.");
              str << " " << V;
            }
            varcmd[0] = varsname[0];
            varcmd[2] = const_cast<char *>(str.str().c_str());
            input->variable->set(3, varcmd);
            echo_var_assign(varcmd[0], varcmd[2]);
          } else {
            for (int j = 0; j < nvars; ++j) {
              varcmd[0] = varsname[j];
              int V;
              kim_error = KIM_Model_GetParameterInteger(pkim, param_index,
                                                        nlbound - 1 + j, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterInteger returned error.");
              std::stringstream str;
              str << V;
              varcmd[2] = const_cast<char *>(str.str().c_str());
              input->variable->set(3, varcmd);
              echo_var_assign(varcmd[0], varcmd[2]);
            }
          }
        } else
          error->all(FLERR, "Wrong parameter type.");

        delete[] varcmd;
        delete[] varsname;
      } // End of loop over all the arguments
    // Set the parameters
    } else {
      std::string set_cmd("pair_coeff * * ");
      set_cmd += atom_type_list;
      for (int i = 1; i < narg; ++i) {
        set_cmd += " ";
        set_cmd += arg[i];
      }
      input->one(set_cmd.c_str());
    }
  } else
    error->all(FLERR, "This model has No mutable parameters.");

  if (!isPairStyleAssigned)
    KIM_Model_Destroy(&pkim);

  kim_param_log_delimiter("end");
}

/* ---------------------------------------------------------------------- */

void KimParam::kim_param_log_delimiter(std::string const &begin_end) const
{
  if (comm->me == 0) {
    std::string msg;
    if (begin_end == "begin") {
      msg = "#=== BEGIN kim-param ";
      msg += kim_param_get ? "get " : "set ";
      msg += "=====================================\n";
    } else if (begin_end == "end") {
      msg = "#=== END kim-param ";
      msg += kim_param_get ? "get " : "set ";
      msg += "=======================================\n\n";
    }
    input->write_echo(msg.c_str());
  }
}

/* ---------------------------------------------------------------------- */

void KimParam::echo_var_assign(std::string const &name,
                               std::string const &value) const
{
  if (comm->me == 0) {
    std::string msg;
    msg += "variable " + name + " string " + value + "\n";
    input->write_echo(msg.c_str());
  }
}

#undef SNUM
