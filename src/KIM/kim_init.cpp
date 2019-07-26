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
#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>
#include "error.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "update.h"
#include "universe.h"
#include "input.h"
#include "variable.h"
#include "citeme.h"
#include "fix_store_kim.h"
#include "kim_units.h"

extern "C" {
#include "KIM_SimulatorHeaders.h"
}

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void KimInit::command(int narg, char **arg)
{
  if ((narg < 2) || (narg > 3)) error->all(FLERR,"Illegal kim_init command");

  if (domain->box_exist)
    error->all(FLERR,"Must use 'kim_init' command before "
                     "simulation box is defined");
  char *model_name = new char[strlen(arg[0])+1];
  strcpy(model_name,arg[0]);
  char *user_units = new char[strlen(arg[1])+1];
  strcpy(user_units,arg[1]);
  if (narg == 3) {
    if (strcmp(arg[2],"unit_conversion_mode")==0) unit_conversion_mode = true;
    else { error->all(FLERR,"Illegal kim_init command"); }
  } else unit_conversion_mode = false;

  char *model_units;
  determine_model_type_and_units(model_name, user_units, &model_units);

  write_log_cite(model_name);

  do_init(model_name, user_units, model_units);
}


/* ---------------------------------------------------------------------- */
namespace {
void get_kim_unit_names(
    char const * const system,
    KIM_LengthUnit & lengthUnit,
    KIM_EnergyUnit & energyUnit,
    KIM_ChargeUnit & chargeUnit,
    KIM_TemperatureUnit & temperatureUnit,
    KIM_TimeUnit & timeUnit,
    Error * error)
{
  if ((strcmp(system,"real")==0)) {
    lengthUnit = KIM_LENGTH_UNIT_A;
    energyUnit = KIM_ENERGY_UNIT_kcal_mol;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_fs;
  } else if ((strcmp(system,"metal")==0)) {
    lengthUnit = KIM_LENGTH_UNIT_A;
    energyUnit = KIM_ENERGY_UNIT_eV;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_ps;
  } else if ((strcmp(system,"si")==0)) {
    lengthUnit = KIM_LENGTH_UNIT_m;
    energyUnit = KIM_ENERGY_UNIT_J;
    chargeUnit = KIM_CHARGE_UNIT_C;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_s;
  } else if ((strcmp(system,"cgs")==0)) {
    lengthUnit = KIM_LENGTH_UNIT_cm;
    energyUnit = KIM_ENERGY_UNIT_erg;
    chargeUnit = KIM_CHARGE_UNIT_statC;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_s;
  } else if ((strcmp(system,"electron")==0)) {
    lengthUnit = KIM_LENGTH_UNIT_Bohr;
    energyUnit = KIM_ENERGY_UNIT_Hartree;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_fs;
  } else if ((strcmp(system,"lj")==0)) {
    error->all(FLERR,"LAMMPS unit_style lj not supported by KIM models");
  } else {
    error->all(FLERR,"Unknown unit_style");
  }
}
}  // namespace
void KimInit::determine_model_type_and_units(char * model_name,
                                             char * user_units,
                                             char ** model_units)
{
  KIM_LengthUnit lengthUnit;
  KIM_EnergyUnit energyUnit;
  KIM_ChargeUnit chargeUnit;
  KIM_TemperatureUnit temperatureUnit;
  KIM_TimeUnit timeUnit;
  int units_accepted;
  KIM_Collections * kim_Coll;
  KIM_CollectionItemType itemType;

  int kim_error = KIM_Collections_Create(&kim_Coll);
  if (kim_error) {
    error->all(FLERR,"Unable to access KIM Collections to find Model.");
  }

  kim_error = KIM_Collections_GetItemType(kim_Coll, model_name, &itemType);
  if (kim_error) {
    error->all(FLERR,"KIM Model name not found.");
  }
  KIM_Collections_Destroy(&kim_Coll);

  if (KIM_CollectionItemType_Equal(itemType,
                                   KIM_COLLECTION_ITEM_TYPE_portableModel))
  {
    get_kim_unit_names(user_units, lengthUnit, energyUnit,
                       chargeUnit, temperatureUnit, timeUnit, error);
    KIM_Model * kim_MO;
    int kim_error = KIM_Model_Create(KIM_NUMBERING_zeroBased,
                                     lengthUnit,
                                     energyUnit,
                                     chargeUnit,
                                     temperatureUnit,
                                     timeUnit,
                                     model_name,
                                     &units_accepted,
                                     &kim_MO);

    if (kim_error)
      error->all(FLERR,"Unable to load KIM Simulator Model.");

    model_type = MO;
    KIM_Model_Destroy(&kim_MO);

    if (units_accepted) {
      *model_units = new char[strlen(user_units)+1];
      strcpy(*model_units,user_units);
      return;
    } else if (unit_conversion_mode) {
      int const num_systems = 5;
      char const * const systems[num_systems]
          = {"metal", "real", "si", "cgs", "electron"};
      for (int i=0; i < num_systems; ++i) {
        get_kim_unit_names(systems[i], lengthUnit, energyUnit,
                           chargeUnit, temperatureUnit, timeUnit, error);
        kim_error = KIM_Model_Create(KIM_NUMBERING_zeroBased,
                                     lengthUnit,
                                     energyUnit,
                                     chargeUnit,
                                     temperatureUnit,
                                     timeUnit,
                                     model_name,
                                     &units_accepted,
                                     &kim_MO);
        KIM_Model_Destroy(&kim_MO);
        if (units_accepted) {
          *model_units = new char[strlen(systems[i])+1];
          strcpy(*model_units,systems[i]);
          return;
        }
      } error->all(FLERR,"KIM Model does not support any lammps unit system");
    } else {
      error->all(FLERR,"KIM Model does not support the requested unit system");
    }
  }
  else if (KIM_CollectionItemType_Equal(
               itemType, KIM_COLLECTION_ITEM_TYPE_simulatorModel)) {
    KIM_SimulatorModel * kim_SM;
    kim_error = KIM_SimulatorModel_Create(model_name, &kim_SM);
    if (kim_error)
      error->all(FLERR,"Unable to load KIM Simulator Model.");
    model_type = SM;

    int sim_fields;
    int sim_lines;
    char const * sim_field;
    char const * sim_value;
    KIM_SimulatorModel_GetNumberOfSimulatorFields(kim_SM, &sim_fields);
    KIM_SimulatorModel_CloseTemplateMap(kim_SM);
    for (int i=0; i < sim_fields; ++i) {
      KIM_SimulatorModel_GetSimulatorFieldMetadata(
          kim_SM,i,&sim_lines,&sim_field);

      if (0 == strcmp(sim_field,"units")) {
        KIM_SimulatorModel_GetSimulatorFieldLine(kim_SM,i,0,&sim_value);
        int len=strlen(sim_value)+1;
        *model_units = new char[len]; strcpy(*model_units,sim_value);
        break;
      }
    }
    KIM_SimulatorModel_Destroy(&kim_SM);

    if ((! unit_conversion_mode) && (strcmp(*model_units, user_units)!=0)) {
      std::string mesg("Incompatible units for KIM Simulator Model, "
                       "required units = ");
      mesg += *model_units;
      error->all(FLERR,mesg.c_str());
    }
  }
}


/* ---------------------------------------------------------------------- */

void KimInit::do_init(char *model_name, char *user_units, char *model_units)
{
  // create storage proxy fix. delete existing fix, if needed.

  int ifix = modify->find_fix("KIM_MODEL_STORE");
  if (ifix >= 0) modify->delete_fix(ifix);
  char *fixarg[3];
  fixarg[0] = (char *)"KIM_MODEL_STORE";
  fixarg[1] = (char *)"all";
  fixarg[2] = (char *)"STORE/KIM";
  modify->add_fix(3,fixarg);
  ifix = modify->find_fix("KIM_MODEL_STORE");

  FixStoreKIM *fix_store = (FixStoreKIM *) modify->fix[ifix];
  fix_store->setptr("model_name", (void *) model_name);
  fix_store->setptr("user_units", (void *) user_units);
  fix_store->setptr("model_units", (void *) model_units);

  // Begin output to log file
  kim_init_log_delimiter("begin");

  int kimerror;
  KIM_SimulatorModel * simulatorModel;
  if (model_type == SM)
  {
    kimerror = KIM_SimulatorModel_Create(model_name,&simulatorModel);

    char const *sim_name, *sim_version;
    KIM_SimulatorModel_GetSimulatorNameAndVersion(
        simulatorModel,&sim_name, &sim_version);

    if (0 != strcmp(sim_name,"LAMMPS"))
      error->all(FLERR,"Incompatible KIM Simulator Model");

    if (comm->me == 0) {
      std::string mesg("# Using KIM Simulator Model : ");
      mesg += model_name;
      mesg += "\n";
      mesg += "# For Simulator             : ";
      mesg += std::string(sim_name) + " " + sim_version + "\n";
      mesg += "# Running on                : LAMMPS ";
      mesg += universe->version;
      mesg += "\n";
      mesg += "#\n";

      if (screen) fputs(mesg.c_str(),screen);
      if (logfile) fputs(mesg.c_str(),logfile);
    }

    fix_store->setptr("simulator_model", (void *) simulatorModel);

    // need to call this to have access to (some) simulator model init data.

    KIM_SimulatorModel_CloseTemplateMap(simulatorModel);
  }

  // Define unit conversion factor variables and print to log
  if (unit_conversion_mode) do_variables(user_units, model_units);

  // set units

  std::string cmd("units ");
  cmd += model_units;
  input->one(cmd.c_str());

  if (model_type == SM) {
    int sim_fields, sim_lines;
    char const *sim_field, *sim_value;
    KIM_SimulatorModel_GetNumberOfSimulatorFields(simulatorModel, &sim_fields);

    // init model

    for (int i=0; i < sim_fields; ++i) {
      KIM_SimulatorModel_GetSimulatorFieldMetadata(
          simulatorModel,i,&sim_lines,&sim_field);
      if (0 == strcmp(sim_field,"model-init")) {
        for (int j=0; j < sim_lines; ++j) {
          KIM_SimulatorModel_GetSimulatorFieldLine(
              simulatorModel,i,j,&sim_value);
          input->one(sim_value);
        }
        break;
      }
    }

    // reset template map.
    KIM_SimulatorModel_OpenAndInitializeTemplateMap(simulatorModel);
  }

  // End output to log file
  kim_init_log_delimiter("end");

}

/* ---------------------------------------------------------------------- */

void KimInit::kim_init_log_delimiter(std::string const begin_end) const
{
  if (comm->me == 0) {
    std::string mesg;
    if (begin_end == "begin")
      mesg =
          "#=== BEGIN kim-init ==========================================\n";
    else if (begin_end == "end")
      mesg =
          "#=== END kim-init ============================================\n\n";

    if ((screen) && (input->echo_screen)) fputs(mesg.c_str(),screen);
    if ((logfile) && (input->echo_log)) fputs(mesg.c_str(),logfile);
  }
}

/* ---------------------------------------------------------------------- */

void KimInit::do_variables(char *user_units, char *model_units)
{
  char *from = user_units, *to = model_units;
  Variable *variable = input->variable;

  // refuse conversion from or to reduced units

  if ((strcmp(from,"lj") == 0) || (strcmp(to,"lj") == 0))
    error->all(FLERR,"Cannot set up conversion variables for 'lj' units");

  // get index to internal style variables. create, if needed.
  // set conversion factors for newly created variables.
  double conversion_factor;
  int ier;
  char *args[3];
  std::string var_str;
  args[1] = (char *)"internal";
  args[2] = (char *)"1.0";
  int v_unit;
  int const nunits = 14;
  char *units[nunits] = {(char *)"mass",
                         (char *)"distance",
                         (char *)"time",
                         (char *)"energy",
                         (char *)"velocity",
                         (char *)"force",
                         (char *)"torque",
                         (char *)"temperature",
                         (char *)"pressure",
                         (char *)"viscosity",
                         (char *)"charge",
                         (char *)"dipole",
                         (char *)"efield",
                         (char *)"density"};

  if (comm->me == 0) {
    std::string mesg("# Conversion factors from ");
    mesg += from;
    mesg += " to ";
    mesg += to;
    mesg += ":\n";
    if ((screen) && (input->echo_screen)) fputs(mesg.c_str(),screen);
    if ((logfile) && (input->echo_log)) fputs(mesg.c_str(),logfile);
  }

  for (int i = 0; i < nunits; i++) {
    var_str = std::string("_u_") + std::string(units[i]);
    args[0] = (char *)var_str.c_str();
    v_unit = variable->find(args[0]);
    if (v_unit < 0) {
      variable->set(3,args);
      v_unit = variable->find(args[0]);
    }
    ier = lammps_unit_conversion(units[i],
                                 from,
                                 to,
                                 conversion_factor);
    if (ier != 0) {
      std::string err = std::string("Unable to obtain conversion factor: ") +
                        "unit = " + units[i] + "; "
                        "from = " + from + "; "
                        "to = " + to + ".";
      error->all(FLERR,err.c_str());
    }
    variable->internal_set(v_unit,conversion_factor);
    if (comm->me == 0) {
      std::stringstream mesg;
      mesg << "variable " << std::setw(15) << std::left << var_str
           << " internal "
           << std::setprecision(12) << std::scientific << conversion_factor
           << std::endl;
      if ((screen) && (input->echo_screen)) fputs(mesg.str().c_str(),screen);
      if ((logfile) && (input->echo_log)) fputs(mesg.str().c_str(),logfile);
    }
  }
  if (comm->me == 0) {
    if ((screen) && (input->echo_screen)) fputs("#\n",screen);
    if ((logfile) && (input->echo_log)) fputs("#\n",logfile);
  }
}

/* ---------------------------------------------------------------------- */

void KimInit::write_log_cite(char * model_name)
{
  KIM_Collections * coll;
  int err = KIM_Collections_Create(&coll);
  if (err) return;

  int extent;
  if (model_type == MO)
  {
    err = KIM_Collections_CacheListOfItemMetadataFiles(
        coll,KIM_COLLECTION_ITEM_TYPE_portableModel,model_name,&extent);
  }
  else if (model_type == SM)
  {
    err = KIM_Collections_CacheListOfItemMetadataFiles(
        coll,KIM_COLLECTION_ITEM_TYPE_simulatorModel,model_name,&extent);
  }
  else
  {
    error->all(FLERR,"Unknown model type.");
  }

  if (err)
  {
    KIM_Collections_Destroy(&coll);
    return;
  }

  for (int i = 0; i < extent;++i)
  {
    char const * fileName;
    int availableAsString;
    char const * fileString;
    err = KIM_Collections_GetItemMetadataFile(
        coll,i,&fileName,NULL,NULL,&availableAsString,&fileString);
    if (err) continue;

    if (0 == strncmp("kimcite",fileName,7))
    {
      if ((lmp->citeme) && (availableAsString)) lmp->citeme->add(fileString);
    }
  }

  KIM_Collections_Destroy(&coll);
}
