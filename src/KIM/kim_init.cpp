/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Axel Kohlmeyer (Temple U),
                         Ryan S. Elliott (UMN),
                         Ellad B. Tadmor (UMN),
                         Yaser Afshar (UMN)
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

#include "kim_init.h"

#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_kim.h"
#include "input.h"
#include "kim_units.h"
#include "modify.h"
#include "universe.h"
#include "variable.h"

#include <cstring>

extern "C" {
#include "KIM_SimulatorHeaders.h"
}

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void KimInit::command(int narg, char **arg)
{
  if ((narg < 2) || (narg > 3)) error->all(FLERR, "Illegal 'kim init' command");

  if (domain->box_exist)
    error->all(FLERR, "Must use 'kim init' command before simulation box is defined");

  char *model_name = utils::strdup(arg[0]);
  char *user_units = utils::strdup(arg[1]);
  if (narg == 3) {
    auto arg_str = std::string(arg[2]);
    if (arg_str == "unit_conversion_mode")
      unit_conversion_mode = true;
    else {
      error->all(FLERR,
                 "Illegal 'kim init' command.\n"
                 "The argument followed by unit_style {} is an optional argument and when "
                 "is used must be unit_conversion_mode",
                 user_units);
    }
  } else
    unit_conversion_mode = false;

  char *model_units;
  KIM_Model *pkim = nullptr;

  if (universe->me == 0) std::remove("kim.log");
  if (universe->nprocs > 1) MPI_Barrier(universe->uworld);

  determine_model_type_and_units(model_name, user_units, &model_units, pkim);

  write_log_cite(lmp, model_type, model_name);

  do_init(model_name, user_units, model_units, pkim);
}

/* ---------------------------------------------------------------------- */

namespace {
void get_kim_unit_names(char const *const system, KIM_LengthUnit &lengthUnit,
                        KIM_EnergyUnit &energyUnit, KIM_ChargeUnit &chargeUnit,
                        KIM_TemperatureUnit &temperatureUnit, KIM_TimeUnit &timeUnit, Error *error)
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
  } else if ((system_str == "lj") || (system_str == "micro") || (system_str == "nano")) {
    error->all(FLERR, "LAMMPS unit_style {} not supported by KIM models", system_str);
  } else {
    error->all(FLERR, "Unknown unit_style");
  }
}
}    // namespace

void KimInit::print_dirs(struct KIM_Collections *const collections) const
{
  int kim_error = 0;
  int dirListExtent = 0;
  int dirCounter = 0;

  std::string mesg = "#=== KIM is looking for 'Portable Models' in these directories ===\n";
  std::vector<struct KIM_Collection> collection_list;
  collection_list.push_back(KIM_COLLECTION_currentWorkingDirectory);
  collection_list.push_back(KIM_COLLECTION_environmentVariable);
  collection_list.push_back(KIM_COLLECTION_user);
  collection_list.push_back(KIM_COLLECTION_system);

  for (auto col : collection_list) {
    kim_error = KIM_Collections_CacheListOfDirectoryNames(
        collections, col, KIM_COLLECTION_ITEM_TYPE_portableModel, &dirListExtent);
    if (!kim_error) {
      for (int i = 0; i < dirListExtent; ++i) {
        char const *name;
        kim_error = KIM_Collections_GetDirectoryName(collections, i, &name);
        // Don't check for error due to bug in kim-api-2.2.1 and below.
#if ((KIM_VERSION_MAJOR * 1000 + KIM_VERSION_MINOR) * 1000 + KIM_VERSION_PATCH) <= 2002001
        kim_error = 0;
#endif
        if (!kim_error) mesg += fmt::format("# {:2}: {}\n", ++dirCounter, name);
      }
    }
  }

  dirCounter = 0;
  mesg += "#=== KIM is looking for 'Simulator Models' in these directories ===\n";
  for (auto col : collection_list) {
    kim_error = KIM_Collections_CacheListOfDirectoryNames(
        collections, col, KIM_COLLECTION_ITEM_TYPE_simulatorModel, &dirListExtent);
    if (!kim_error) {
      for (int i = 0; i < dirListExtent; ++i) {
        char const *name;
        kim_error = KIM_Collections_GetDirectoryName(collections, i, &name);
        // Don't check for error due to bug in kim-api-2.2.1 and below.
#if ((KIM_VERSION_MAJOR * 1000 + KIM_VERSION_MINOR) * 1000 + KIM_VERSION_PATCH) <= 2002001
        kim_error = 0;
#endif
        if (!kim_error) mesg += fmt::format("# {:2}: {}\n", ++dirCounter, name);
      }
    }
  }
  input->write_echo(mesg);
}

void KimInit::determine_model_type_and_units(char *model_name, char *user_units, char **model_units,
                                             KIM_Model *&pkim)
{
  KIM_LengthUnit lengthUnit;
  KIM_EnergyUnit energyUnit;
  KIM_ChargeUnit chargeUnit;
  KIM_TemperatureUnit temperatureUnit;
  KIM_TimeUnit timeUnit;
  int units_accepted;
  KIM_Collections *collections;
  KIM_CollectionItemType itemType;

  int kim_error = KIM_Collections_Create(&collections);
  if (kim_error) error->all(FLERR, "Unable to access KIM Collections to find Model");

  auto logID = fmt::format("{}_Collections", comm->me);
  KIM_Collections_SetLogID(collections, logID.c_str());

  print_dirs(collections);

  kim_error = KIM_Collections_GetItemType(collections, model_name, &itemType);
  if (kim_error) error->all(FLERR, "KIM Model name not found");
  KIM_Collections_Destroy(&collections);

  if (KIM_CollectionItemType_Equal(itemType, KIM_COLLECTION_ITEM_TYPE_portableModel)) {
    get_kim_unit_names(user_units, lengthUnit, energyUnit, chargeUnit, temperatureUnit, timeUnit,
                       error);
    int kim_error = KIM_Model_Create(KIM_NUMBERING_zeroBased, lengthUnit, energyUnit, chargeUnit,
                                     temperatureUnit, timeUnit, model_name, &units_accepted, &pkim);

    if (kim_error) error->all(FLERR, "Unable to load KIM Simulator Model");

    model_type = MO;

    if (units_accepted) {
      logID = fmt::format("{}_Model", comm->me);
      KIM_Model_SetLogID(pkim, logID.c_str());
      *model_units = utils::strdup(user_units);
      return;
    } else if (unit_conversion_mode) {
      KIM_Model_Destroy(&pkim);
      const char *unit_systems[] = {"metal", "real", "si", "cgs", "electron"};
      for (auto units : unit_systems) {
        get_kim_unit_names(units, lengthUnit, energyUnit, chargeUnit, temperatureUnit, timeUnit,
                           error);
        kim_error = KIM_Model_Create(KIM_NUMBERING_zeroBased, lengthUnit, energyUnit, chargeUnit,
                                     temperatureUnit, timeUnit, model_name, &units_accepted, &pkim);
        if (units_accepted) {
          logID = fmt::format("{}_Model", comm->me);
          KIM_Model_SetLogID(pkim, logID.c_str());
          *model_units = utils::strdup(units);
          return;
        }
        KIM_Model_Destroy(&pkim);
      }
      error->all(FLERR, "KIM Model does not support any lammps unit system");
    } else {
      KIM_Model_Destroy(&pkim);
      error->all(FLERR, "KIM Model does not support the requested unit system");
    }
  } else if (KIM_CollectionItemType_Equal(itemType, KIM_COLLECTION_ITEM_TYPE_simulatorModel)) {
    KIM_SimulatorModel *simulatorModel;
    kim_error = KIM_SimulatorModel_Create(model_name, &simulatorModel);
    if (kim_error) error->all(FLERR, "Unable to load KIM Simulator Model");
    model_type = SM;

    logID = fmt::format("{}_SimulatorModel", comm->me);
    KIM_SimulatorModel_SetLogID(simulatorModel, logID.c_str());

    int sim_fields;
    int sim_lines;
    char const *sim_field;
    char const *sim_value;
    KIM_SimulatorModel_GetNumberOfSimulatorFields(simulatorModel, &sim_fields);
    KIM_SimulatorModel_CloseTemplateMap(simulatorModel);
    for (int i = 0; i < sim_fields; ++i) {
      KIM_SimulatorModel_GetSimulatorFieldMetadata(simulatorModel, i, &sim_lines, &sim_field);

      const std::string sim_field_str(sim_field);
      if (sim_field_str == "units") {
        KIM_SimulatorModel_GetSimulatorFieldLine(simulatorModel, i, 0, &sim_value);
        *model_units = utils::strdup(sim_value);
        break;
      }
    }
    KIM_SimulatorModel_Destroy(&simulatorModel);

    const std::string model_units_str(*model_units);
    const std::string user_units_str(user_units);
    if ((!unit_conversion_mode) && (model_units_str != user_units_str)) {
      error->all(FLERR, "Incompatible units for KIM Simulator Model, required units = {}",
                 model_units_str);
    }
  }
}

/* ---------------------------------------------------------------------- */

void KimInit::do_init(char *model_name, char *user_units, char *model_units, KIM_Model *&pkim)
{
  // create storage proxy fix. delete existing fix, if needed.

  int ifix = modify->find_fix("KIM_MODEL_STORE");
  if (ifix >= 0) modify->delete_fix(ifix);
  modify->add_fix("KIM_MODEL_STORE all STORE/KIM");
  ifix = modify->find_fix("KIM_MODEL_STORE");

  auto fix_store = dynamic_cast<FixStoreKIM *>(modify->fix[ifix]);
  fix_store->setptr("model_name", (void *) model_name);
  fix_store->setptr("user_units", (void *) user_units);
  fix_store->setptr("model_units", (void *) model_units);

  // Begin output to log file
  input->write_echo("#=== BEGIN kim init ==========================================\n");

  KIM_SimulatorModel *simulatorModel;
  if (model_type == SM) {
    int kim_error = KIM_SimulatorModel_Create(model_name, &simulatorModel);
    if (kim_error) error->all(FLERR, "Unable to load KIM Simulator Model");

    auto logID = fmt::format("{}_SimulatorModel", comm->me);
    KIM_SimulatorModel_SetLogID(simulatorModel, logID.c_str());

    char const *sim_name, *sim_version;
    KIM_SimulatorModel_GetSimulatorNameAndVersion(simulatorModel, &sim_name, &sim_version);

    const std::string sim_name_str(sim_name);
    if (sim_name_str != "LAMMPS") error->all(FLERR, "Incompatible KIM Simulator Model");

    if (comm->me == 0) {
      auto mesg = fmt::format("# Using KIM Simulator Model : {}\n"
                              "# For Simulator             : {} {}\n"
                              "# Running on                : LAMMPS {}\n#\n",
                              model_name, sim_name_str, sim_version, lmp->version);
      utils::logmesg(lmp, mesg);
    }

    fix_store->setptr("simulator_model", (void *) simulatorModel);

    // need to call this to have access to (some) simulator model init data.

    KIM_SimulatorModel_CloseTemplateMap(simulatorModel);
  }

  // Define unit conversion factor variables and print to log
  if (unit_conversion_mode) do_variables(user_units, model_units);

  // set units

  const std::string model_units_str(model_units);
  auto cmd = fmt::format("units {}", model_units_str);
  input->one(cmd);

  // Set the skin and timestep default values as
  // 2.0 Angstroms and 1.0 femtosecond

  const std::string skin_cmd = (model_units_str == "real") ? "neighbor 2.0 bin   # Angstroms"
      : (model_units_str == "metal")                       ? "neighbor 2.0 bin   # Angstroms"
      : (model_units_str == "si")                          ? "neighbor 2e-10 bin   # meters"
      : (model_units_str == "cgs")                         ? "neighbor 2e-8 bin   # centimeters"
                                                           : "neighbor 3.77945224 bin   # Bohr";
  const std::string step_cmd = (model_units_str == "real") ? "timestep 1.0       # femtoseconds"
      : (model_units_str == "metal")                       ? "timestep 1.0e-3    # picoseconds"
      : (model_units_str == "si")                          ? "timestep 1e-15       # seconds"
      : (model_units_str == "cgs")                         ? "timestep 1e-15      # seconds"
                                   : "timestep 1.0              # femtoseconds";
  input->one(skin_cmd);
  input->one(step_cmd);

  if (model_type == SM) {
    int sim_fields, sim_lines;
    char const *sim_field, *sim_value;
    KIM_SimulatorModel_GetNumberOfSimulatorFields(simulatorModel, &sim_fields);

    // init model

    for (int i = 0; i < sim_fields; ++i) {
      KIM_SimulatorModel_GetSimulatorFieldMetadata(simulatorModel, i, &sim_lines, &sim_field);

      const std::string sim_field_str(sim_field);
      if (sim_field_str == "model-init") {
        for (int j = 0; j < sim_lines; ++j) {
          KIM_SimulatorModel_GetSimulatorFieldLine(simulatorModel, i, j, &sim_value);
          input->one(sim_value);
        }
        break;
      }
    }

    // reset template map.
    KIM_SimulatorModel_OpenAndInitializeTemplateMap(simulatorModel);
  } else if (model_type == MO) {
    int numberOfParameters;
    KIM_Model_GetNumberOfParameters(pkim, &numberOfParameters);

    std::string mesg = "\nThis model has ";
    if (numberOfParameters) {
      KIM_DataType kim_DataType;
      int extent;
      char const *str_name = nullptr;
      char const *str_desc = nullptr;

      mesg += std::to_string(numberOfParameters) + " mutable parameters. \n";

      int max_len(0);
      for (int i = 0; i < numberOfParameters; ++i) {
        KIM_Model_GetParameterMetadata(pkim, i, &kim_DataType, &extent, &str_name, &str_desc);
        max_len = MAX(max_len, (int) strlen(str_name));
      }
      max_len = MAX(18, max_len + 1);
      mesg += fmt::format(" No.      | {:<{}} | data type  | extent\n", "Parameter name", max_len);
      mesg += fmt::format("{:-<{}}\n", "-", max_len + 35);
      for (int i = 0; i < numberOfParameters; ++i) {
        KIM_Model_GetParameterMetadata(pkim, i, &kim_DataType, &extent, &str_name, &str_desc);
        auto data_type = std::string("\"");
        data_type += KIM_DataType_ToString(kim_DataType) + std::string("\"");
        mesg += fmt::format(" {:<8} | {:<{}} | {:<10} | {}\n", i + 1, str_name, max_len, data_type,
                            extent);
      }
    } else
      mesg += "No mutable parameters.\n";

    KIM_Model_Destroy(&pkim);
    input->write_echo(mesg);
  }

  // End output to log file
  input->write_echo("#=== END kim init ============================================\n\n");
}

/* ---------------------------------------------------------------------- */

void KimInit::do_variables(const std::string &from, const std::string &to)
{
  // refuse conversion from or to reduced units

  if ((from == "lj") || (to == "lj"))
    error->all(FLERR, "Cannot set up conversion variables for 'lj' units");

  // get index to internal style variables. create, if needed.
  // set conversion factors for newly created variables.
  double conversion_factor;
  int ier;
  std::string var_str;
  int v_unit;
  const char *units[] = {"mass",   "distance", "time",        "energy",   "velocity",
                         "force",  "torque",   "temperature", "pressure", "viscosity",
                         "charge", "dipole",   "efield",      "density",  nullptr};

  input->write_echo(fmt::format("# Conversion factors from {} to {}:\n", from, to));

  auto variable = input->variable;
  for (int i = 0; units[i] != nullptr; ++i) {
    var_str = std::string("_u_") + units[i];
    v_unit = variable->find(var_str.c_str());
    if (v_unit < 0) {
      variable->set(var_str + " internal 1.0");
      v_unit = variable->find(var_str.c_str());
    }
    ier = lammps_unit_conversion(units[i], from, to, conversion_factor);
    if (ier != 0)
      error->all(FLERR,
                 "Unable to obtain conversion factor: "
                 "unit = {}; from = {}; to = {}",
                 units[i], from, to);

    variable->internal_set(v_unit, conversion_factor);
    input->write_echo(
        fmt::format("variable {:<15s} internal {:<15.12e}\n", var_str, conversion_factor));
  }
  input->write_echo("#\n");
}

/* ---------------------------------------------------------------------- */

void KimInit::write_log_cite(class LAMMPS *lmp, KimInit::model_type_enum model_type,
                             char *model_name)
{
  if (!lmp->citeme) return;

  std::string model_name_str(model_name);
  std::string re = "[MS][OM]_\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d_\\d\\d\\d";
  std::string kim_id = utils::strfind(model_name_str, re);

  std::string cite_id;
  if (kim_id.empty()) {
    cite_id = fmt::format("KIM potential: unpublished, \"{}\"\n", model_name_str);
  } else {
    KIM_Collections *collections;
    int err = KIM_Collections_Create(&collections);
    if (err) return;

    auto logID = fmt::format("{}_Collections", lmp->comm->me);
    KIM_Collections_SetLogID(collections, logID.c_str());

    int extent;
    if (model_type == MO) {
      err = KIM_Collections_CacheListOfItemMetadataFiles(
          collections, KIM_COLLECTION_ITEM_TYPE_portableModel, model_name, &extent);
    } else if (model_type == SM) {
      err = KIM_Collections_CacheListOfItemMetadataFiles(
          collections, KIM_COLLECTION_ITEM_TYPE_simulatorModel, model_name, &extent);
    } else {
      lmp->error->all(FLERR, "Unknown model type");
    }

    if (err) {
      KIM_Collections_Destroy(&collections);
      return;
    }

    cite_id = fmt::format("OpenKIM potential: https://openkim.org/cite/"
                          "{}#item-citation\n\n",
                          kim_id);

    for (int i = 0; i < extent; ++i) {
      char const *fileName;
      int availableAsString;
      char const *fileString;
      err = KIM_Collections_GetItemMetadataFile(collections, i, &fileName, nullptr, nullptr,
                                                &availableAsString, &fileString);
      if (err) continue;

      if (utils::strmatch(fileName, "^kimcite") && availableAsString) cite_id += fileString;
    }
    KIM_Collections_Destroy(&collections);
  }

  lmp->citeme->add(cite_id);
}
