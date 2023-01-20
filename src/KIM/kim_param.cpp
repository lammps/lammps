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

#include "comm.h"
#include "error.h"
#include "fix_store_kim.h"
#include "force.h"
#include "input.h"
#include "modify.h"
#include "pair_kim.h"
#include "variable.h"

#include <cstdlib>
#include <cstring>
#include <vector>

extern "C"
{
#include "KIM_SimulatorHeaders.h"
}

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
  const std::string system_str(system);
  if (system_str == "real") {
    lengthUnit = KIM_LENGTH_UNIT_A;
    energyUnit = KIM_ENERGY_UNIT_kcal_mol;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_fs;
  } else if (system_str == "metal") {
    lengthUnit = KIM_LENGTH_UNIT_A;
    energyUnit = KIM_ENERGY_UNIT_eV;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_ps;
  } else if (system_str == "si") {
    lengthUnit = KIM_LENGTH_UNIT_m;
    energyUnit = KIM_ENERGY_UNIT_J;
    chargeUnit = KIM_CHARGE_UNIT_C;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_s;
  } else if (system_str == "cgs") {
    lengthUnit = KIM_LENGTH_UNIT_cm;
    energyUnit = KIM_ENERGY_UNIT_erg;
    chargeUnit = KIM_CHARGE_UNIT_statC;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_s;
  } else if (system_str == "electron") {
    lengthUnit = KIM_LENGTH_UNIT_Bohr;
    energyUnit = KIM_ENERGY_UNIT_Hartree;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_fs;
  } else if ((system_str == "lj") ||
             (system_str == "micro") ||
             (system_str == "nano")) {
    error->all(FLERR, "LAMMPS unit_style {} not supported "
                                  "by KIM models", system_str);
  } else {
    error->all(FLERR, "Unknown unit_style");
  }
}
} // namespace

/* ---------------------------------------------------------------------- */

KimParam::KimParam(LAMMPS *lmp) : Pointers(lmp) {}

void KimParam::command(int narg, char **arg)
{
  // kim param is a command for
  // getting/setting the value of a %KIM PM parameter
  //
  // kim param get param_name index_range variables formatarg
  // kim param set param_name index_range values
  //
  // kim param get paramname 1 varname
  // kim param get paramname index_range varname_1, ..., varname_N
  // kim param get paramname index_range varname_base split
  // kim param get paramname index_range varname_base list
  // kim param set paramname index_range values

  if (narg < 4) error->all(FLERR, "Illegal 'kim param' command");

  std::string kim_param_get_set(arg[0]);

  if ((kim_param_get_set != "get") && (kim_param_get_set != "set"))
    error->all(FLERR, "Incorrect arguments in 'kim param' command.\n"
               "'kim param get/set' is mandatory");

  int const ifix = modify->find_fix("KIM_MODEL_STORE");
  if (ifix >= 0) {
    auto fix_store = reinterpret_cast<FixStoreKIM *>(modify->fix[ifix]);

    KIM_SimulatorModel *simulatorModel =
        reinterpret_cast<KIM_SimulatorModel *>(
            fix_store->getptr("simulator_model"));

    if (simulatorModel)
      error->all(FLERR, "'kim param' can only be used with a KIM Portable Model");
  }

  input->write_echo(fmt::format("#=== BEGIN kim param {} ==================="
                                "==================\n", kim_param_get_set));

  KIM_Model *pkim = nullptr;

  std::string atom_type_list;

  if (force->pair) {
    Pair *pair = force->pair_match("kim", 1, 0);
    if (pair) {
      auto pairKIM = reinterpret_cast<PairKIM *>(pair);

      pkim = pairKIM->get_kim_model();
      if (!pkim) error->all(FLERR, "Unable to get the KIM Portable Model");

      if (kim_param_get_set == "set") {
        atom_type_list = pairKIM->get_atom_type_list();
        if (atom_type_list.empty()) error->all(FLERR, "The requested atom type list is empty");
      }
    } else
      error->all(FLERR, "Pair style is defined, but there is "
                        "no match for kim style in lammps");
  } else
    error->all(FLERR, "Illegal 'kim param {0}' command.\nTo {0} the new parameter values, "
               "pair style must be assigned.\nMust use 'kim interactions' or 'pair_style kim' "
               "before 'kim param {0}'", kim_param_get_set);

  // Get the number of mutable parameters in the kim model
  int numberOfParameters(0);

  KIM_Model_GetNumberOfParameters(pkim, &numberOfParameters);
  if (numberOfParameters) {
    // Get the parameters
    if (kim_param_get_set == "get") {
      int kim_error;
      // Parameter name
      std::string paramname;
      // Variable name
      std::string varname;

      // Loop over all the arguments
      for (int i = 1; i < narg;) {
        // Parameter name
        if (i < narg)
          paramname = std::string(arg[i++]);
        else
          break;

        // Find the requested parameter within the model parameters
        int param_index;
        KIM_DataType kim_DataType;
        int extent;
        char const *str_name = nullptr;
        char const *str_desc = nullptr;

        for (param_index = 0; param_index < numberOfParameters; ++param_index) {
          kim_error = KIM_Model_GetParameterMetadata(pkim, param_index,
                                                     &kim_DataType, &extent,
                                                     &str_name, &str_desc);
          if (kim_error)
            error->all(FLERR, "KIM GetParameterMetadata returned error");

          const std::string str_name_str(str_name);
          if (paramname == str_name_str) break;
        }

        if (param_index >= numberOfParameters)
          error->all(FLERR, "Wrong argument in 'kim param get' command.\n"
                     "This Model does not have the requested '{}' parameter", paramname);

        // Get the index_range for the requested parameter
        int nlbound(0);
        int nubound(0);

        if (i < narg) {
          std::string argtostr(arg[i++]);

          // Check to see if the indices range contains
          // only integer numbers and/or range :
          if (argtostr.find_first_not_of("0123456789:") != std::string::npos)
            error->all(FLERR, "Illegal index_range.\nExpected integer parameter(s) instead "
                       "of '{}' in index_range", argtostr);

          std::string::size_type npos = argtostr.find(':');
          if (npos != std::string::npos) {
            argtostr[npos] = ' ';
            auto words = utils::split_words(argtostr);
            nlbound = atoi(words[0].c_str());
            nubound = atoi(words[1].c_str());

            if ((nubound < 1) || (nubound > extent) || (nlbound < 1) || (nlbound > nubound))
              error->all(FLERR, "Illegal index_range '{}-{}' for '{}' parameter with the "
                         "extent of '{}'",nlbound, nubound, paramname, extent);
          } else {
            nlbound = atoi(argtostr.c_str());

            if (nlbound < 1 || nlbound > extent)
              error->all(FLERR, "Illegal index '{}' for '{}' parameter with the extent of '{}'",
                         nlbound, paramname, extent);
            nubound = nlbound;
          }
        } else
          error->all(FLERR, "Wrong number of arguments in 'kim param get' command.\n"
                     "Index range after parameter name is mandatory");

        int const nvars = nubound - nlbound + 1;
        std::vector<std::string> varsname;

        if (i < narg) {
          // Get the variable/variable_base name
          varname = std::string(arg[i++]);
          if ((varname == "split") || (varname == "list") || (varname == "explicit"))
            error->all(FLERR, "Illegal variable name '{}' in 'kim param get'", varname);
        } else
          error->all(FLERR, "Wrong number of arguments in 'kim param get' command.\n"
                     "The LAMMPS variable name is mandatory");

        // indicator flag for list request
        bool list_requested(false);

        if (nvars > 1) {
          if (i < narg) {
            std::string formatarg(arg[i]);
            if (formatarg == "split") {
              varsname.resize(nvars);
              for (int j = 0, k = nlbound; j < nvars; ++j, ++k) {
                varsname[j] = fmt::format("{}_{}", varname, k);
              }
              ++i;
            } else if (formatarg == "list") {
              list_requested = true;
              varsname.resize(1);
              varsname[0] = varname;
              ++i;
            // Default explicit (optional) formatarg
            } else if (i - 1 + nvars - 1 < narg) {
              varsname.resize(nvars);
              --i;
              for (int j = 0; j < nvars; ++j, ++i) {
                varsname[j] = std::string(arg[i]);
                if (varsname[j] == "split" || varsname[j] == "list" || varsname[j] == "explicit")
                  error->all(FLERR, "Illegal variable name '{}' in 'kim param get'", varsname[j]);
              }
              if (i < narg) {
                formatarg = std::string(arg[i]);
                if (formatarg == "explicit") ++i;
              }
            } else {
              error->all(FLERR, "Wrong number of arguments in 'kim param get' command.\n"
                         "The LAMMPS '{}' variable names or '{} split' is mandatory",
                         nvars, varname);
            }
          } else
            error->all(FLERR, "Wrong number of arguments in 'kim param get' command.\n"
                       "The LAMMPS '{}' variable names or '{} split/list' is mandatory",
                       nvars, varname);
        } else {
          varsname.resize(1);
          if (i < narg) {
            const std::string formatarg(arg[i]);
            if (formatarg == "split") {
              varsname[0] = fmt::format("{}_{}", varname, nlbound);
              ++i;
            } else {
              if (formatarg == "list" || formatarg == "explicit") ++i;
              varsname[0] = varname;
            }
          } else {
            varsname[0] = varname;
          }
        }

        if (KIM_DataType_Equal(kim_DataType, KIM_DATA_TYPE_Double)) {
          if (list_requested) {
            std::string str;
            double V;
            {
              kim_error = KIM_Model_GetParameterDouble(pkim, param_index,
                                                       nlbound - 1, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterDouble returned error");

              str = fmt::format("{}", V);
            }
            for (int j = 1; j < nvars; ++j) {
              kim_error = KIM_Model_GetParameterDouble(pkim, param_index,
                                                       nlbound - 1 + j, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterDouble returned error");

              str += fmt::format(" {}", V);
            }

            auto setcmd = fmt::format("{} string \"{}\"", varsname[0], str);
            input->variable->set(setcmd);
            input->write_echo(fmt::format("variable {}\n", setcmd));

          } else {
            double V;
            for (int j = 0; j < nvars; ++j) {
              kim_error = KIM_Model_GetParameterDouble(pkim, param_index,
                                                       nlbound - 1 + j, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterDouble returned error");

              auto setcmd = fmt::format("{} string {}", varsname[j], V);
              input->variable->set(setcmd);
              input->write_echo(fmt::format("variable {}\n", setcmd));
            }
          }
        } else if (KIM_DataType_Equal(kim_DataType, KIM_DATA_TYPE_Integer)) {
          if (list_requested) {
            std::string str;
            int V;
            {
              kim_error = KIM_Model_GetParameterInteger(pkim, param_index,
                                                        nlbound - 1, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterInteger returned error");

              str = fmt::format("{}", V);
            }
            for (int j = 1; j < nvars; ++j) {
              kim_error = KIM_Model_GetParameterInteger(pkim, param_index,
                                                        nlbound - 1 + j, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterInteger returned error");

              str += fmt::format(" {}", V);
            }

            auto setcmd = fmt::format("{} string \"{}\"", varsname[0], str);
            input->variable->set(setcmd);
            input->write_echo(fmt::format("variable {}\n", setcmd));

          } else {
            int V;
            for (int j = 0; j < nvars; ++j) {
              kim_error = KIM_Model_GetParameterInteger(pkim, param_index,
                                                        nlbound - 1 + j, &V);
              if (kim_error)
                error->all(FLERR, "KIM GetParameterInteger returned error");

              auto setcmd = fmt::format("{} string {}", varsname[j], V);
              input->variable->set(setcmd);
              input->write_echo(fmt::format("variable {}\n", setcmd));
            }
          }
        } else
          error->all(FLERR, "Wrong parameter type");
      } // End of loop over all the arguments
    // Set the parameters
    } else {
      auto setcmd = fmt::format("pair_coeff * * {} {}", atom_type_list,
                                fmt::join(arg + 1, arg + narg, " "));
      input->one(setcmd);
    }
  } else
    error->all(FLERR, "This model has No mutable parameters");

  input->write_echo(fmt::format("#=== END kim param {} ====================="
                                "==================\n", kim_param_get_set));
}
