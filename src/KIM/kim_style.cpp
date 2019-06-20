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

#include <cstring>
#include <string>
#include "kim_style.h"
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
//@@@@@ Need to switch to c-bindings when they are available.
#include "KIM_SimulatorModel.hpp"
//@@@@@

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void KimStyle::command(int narg, char **arg)
{
  if ((narg < 2) || (narg > 4)) error->all(FLERR,"Illegal kim_style command");

  if (strcmp(arg[0],"model") == 0) {
    if (domain->box_exist)
      error->all(FLERR,"Must use 'kim_style model' command before "
                 "simulation box is defined");
    int len1 = strlen(arg[1])+1;
    int len2 = strlen(arg[2])+1;
    char *model_name = new char[len1]; strcpy(model_name,arg[1]);
    char *user_units = new char[len2]; strcpy(user_units,arg[2]);
    if (narg == 4) {
      if (strcmp(arg[3],"unit_conversion_mode")==0) unit_conversion_mode = true;
      else { error->all(FLERR,"Illegal kim_style command"); }
    } else unit_conversion_mode = false;

    char *model_units;
    determine_model_type_and_units(model_name, user_units, &model_units);

    do_init(model_name, user_units, model_units);
    if (unit_conversion_mode) do_variables(user_units, model_units);
  } else if (strcmp(arg[0],"setup") == 0) {
    if (!domain->box_exist)
      error->all(FLERR,"Must use 'kim_style setup' command after "
                 "simulation box is defined");
    do_setup(narg-1,++arg);
  } else error->all(FLERR,"Illegal kim_style command");
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
void KimStyle::determine_model_type_and_units(char * model_name,
                                              char * user_units,
                                              char ** model_units)
{
  KIM_LengthUnit lengthUnit;
  KIM_EnergyUnit energyUnit;
  KIM_ChargeUnit chargeUnit;
  KIM_TemperatureUnit temperatureUnit;
  KIM_TimeUnit timeUnit;
  int units_accepted;
  KIM_Model * kim_MO;

  get_kim_unit_names(user_units, lengthUnit, energyUnit,
                     chargeUnit, temperatureUnit, timeUnit, error);
  int kim_error = KIM_Model_Create(KIM_NUMBERING_zeroBased,
                                   lengthUnit,
                                   energyUnit,
                                   chargeUnit,
                                   temperatureUnit,
                                   timeUnit,
                                   model_name,
                                   &units_accepted,
                                   &kim_MO);

  if (!kim_error)  // model is an MO
  {
    model_type = MO;
    KIM_Model_Destroy(&kim_MO);

    if (units_accepted)
    {
      int len=strlen(user_units);
      *model_units = new char[len]; strcpy(*model_units,user_units);
      return;
    }
    else if (unit_conversion_mode)
    {
      int const num_systems = 5;
      char const * const systems[num_systems]
          = {"metal", "real", "si", "cgs", "electron"};
      for (int i=0; i < num_systems; ++i)
      {
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
        if (units_accepted)
        {
          int len=strlen(systems[i]);
          *model_units = new char[len]; strcpy(*model_units,systems[i]);
          return;
        }
      }
      error->all(FLERR,"KIM Model does not support any lammps unit system");
    }
    else
    {
      error->all(FLERR,"KIM Model does not support the requested unit system");
    }
  }

  KIM::SimulatorModel * kim_SM;
  kim_error = KIM::SimulatorModel::Create(model_name, &kim_SM);
  if (kim_error)
  {
    error->all(FLERR,"KIM model name not found");
  }
  model_type = SM;

  int sim_fields;
  int sim_lines;
  std::string const * sim_field;
  std::string const * sim_value;
  kim_SM->GetNumberOfSimulatorFields(&sim_fields);
  kim_SM->CloseTemplateMap();
  for (int i=0; i < sim_fields; ++i) {
    kim_SM->GetSimulatorFieldMetadata(i,&sim_lines,&sim_field);

    if (*sim_field == "units") {
      kim_SM->GetSimulatorFieldLine(i,0,&sim_value);
      int len=(*sim_value).length();
      *model_units = new char[len]; strcpy(*model_units,sim_value->c_str());
      break;
    }
  }
  KIM::SimulatorModel::Destroy(&kim_SM);

  if ((! unit_conversion_mode) && (strcmp(*model_units, user_units)!=0))
  {
    error->all(FLERR,"Incompatible units for KIM Simulator Model");
  }
}


/* ---------------------------------------------------------------------- */

void KimStyle::do_init(char *model_name, char *user_units, char* model_units)
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

  // set units

  std::string cmd("units ");
  cmd += model_units;
  input->one(cmd.c_str());

  int kimerror;
  // @@@@@ switch to c-bindings when they are available
  KIM::SimulatorModel * simulatorModel;
  kimerror = KIM::SimulatorModel::Create(model_name,&simulatorModel);

  // not a Kim Simulator Model; nothing else to do here.

  if (kimerror) return;

  fix_store->setptr("simulator_model", (void *) simulatorModel);

  // need to call this to have access to (some) simulator model init data.

  simulatorModel->CloseTemplateMap();

  int sim_fields, sim_lines;
  const std::string *sim_field, *sim_value;
  simulatorModel->GetNumberOfSimulatorFields(&sim_fields);

  // init model

  for (int i=0; i < sim_fields; ++i) {
    simulatorModel->GetSimulatorFieldMetadata(i,&sim_lines,&sim_field);
    if (*sim_field == "model-init") {
      for (int j=0; j < sim_lines; ++j) {
        simulatorModel->GetSimulatorFieldLine(i,j,&sim_value);
        input->one(sim_value->c_str());
      }
      break;
    }
  }

  // reset template map.

  simulatorModel->OpenAndInitializeTemplateMap();
}

/* ---------------------------------------------------------------------- */

void KimStyle::do_setup(int narg, char **arg)
{
  if (narg != atom->ntypes)
    error->all(FLERR,"Illegal kim_style command");

  char *model = NULL;
  KIM::SimulatorModel *simulatorModel(NULL);

  // check if we had a kim_style init command by finding fix STORE/KIM
  // retrieve model name and pointer to simulator model class instance.
  // validate model name if not given as NULL.
  // if kim_style init wasn't run try to initialize simulator model now.

  int ifix = modify->find_fix("KIM_MODEL_STORE");
  if (ifix >= 0) {
    FixStoreKIM *fix_store = (FixStoreKIM *) modify->fix[ifix];
    model = (char *)fix_store->getptr("model_name");
    simulatorModel = (KIM::SimulatorModel *)fix_store->getptr("simulator_model");
  } else error->all(FLERR,"Must use 'kim_style model' before 'kim_style setup'");

  if (simulatorModel) {

    const std::string *sim_name, *sim_version;
    std::string atom_type_sym_list;

    simulatorModel->GetSimulatorNameAndVersion(&sim_name, &sim_version);

    if (comm->me == 0) {
      std::string mesg("Using KIM Simulator Model : ");
      mesg += model;
      mesg += "\n";
      mesg += "For Simulator             : ";
      mesg += *sim_name + " " + *sim_version + "\n";
      mesg += "Running on                : LAMMPS ";
      mesg += universe->version;
      mesg += "\n";

      if (screen) fputs(mesg.c_str(),screen);
      if (logfile) fputs(mesg.c_str(),logfile);
    }

    if (*sim_name != "LAMMPS")
      error->all(FLERR,"Incompatible KIM Simulator Model");

    for (int i = 0; i < narg; i++)
      atom_type_sym_list += std::string(" ") + arg[i];

    simulatorModel->AddTemplateMap("atom-type-sym-list",atom_type_sym_list);
    simulatorModel->CloseTemplateMap();

    int len = strlen(atom_type_sym_list.c_str())+1;
    char *strbuf = new char[len];
    char *strword;

    // validate species selection

    int sim_num_species;
    bool species_is_supported;
    const std::string *sim_species;
    simulatorModel->GetNumberOfSupportedSpecies(&sim_num_species);
    strcpy(strbuf,atom_type_sym_list.c_str());
    strword = strtok(strbuf," \t");
    while (strword) {
      species_is_supported = false;
      if (strcmp(strword,"NULL") == 0) continue;
      for (int i=0; i < sim_num_species; ++i) {
        simulatorModel->GetSupportedSpecies(i, &sim_species);
        if (strcmp(sim_species->c_str(),strword) == 0)
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

    // check if units are unchanged, and if kim_style init was required

    int sim_fields, sim_lines;
    const std::string *sim_field, *sim_value;
    simulatorModel->GetNumberOfSimulatorFields(&sim_fields);
    for (int i=0; i < sim_fields; ++i) {
      simulatorModel->GetSimulatorFieldMetadata(i,&sim_lines,&sim_field);

      if (*sim_field == "units") {
        simulatorModel->GetSimulatorFieldLine(i,0,&sim_value);
        if (*sim_value != update->unit_style)
          error->all(FLERR,"Incompatible units for KIM Simulator Model");
      }
    }

    int sim_model_idx=-1;
    for (int i=0; i < sim_fields; ++i) {
      simulatorModel->GetSimulatorFieldMetadata(i,&sim_lines,&sim_field);
      if (*sim_field == "model-defn") {
        sim_model_idx = i;
        for (int j=0; j < sim_lines; ++j) {
          simulatorModel->GetSimulatorFieldLine(sim_model_idx,j,&sim_value);
          input->one(sim_value->c_str());
        }
      }
    }

    if (sim_model_idx < 0)
      error->all(FLERR,"KIM Simulator Model has no Model definition");

    simulatorModel->OpenAndInitializeTemplateMap();

  } else {

    // not a simulator model. issue pair_style and pair_coeff commands.
    // NOTE: all references to arg must appear before calls to input->one()
    // as that will reset the argument vector.

    std::string cmd1("pair_style kim ");
    cmd1 += model;

    std::string cmd2("pair_coeff * * ");
    for (int i=0; i < narg; ++i) {
      cmd2 += arg[i];
      cmd2 += " ";
    }

    input->one(cmd1.c_str());
    input->one(cmd2.c_str());
  }
}

/* ---------------------------------------------------------------------- */

void KimStyle::do_variables(char *user_units, char *model_units)
{
  char *from = user_units, *to = model_units;
  Variable *variable = input->variable;

  // refuse convertion from or to reduced units

  if ((strcmp(from,"lj") == 0) || (strcmp(to,"lj") == 0))
    error->all(FLERR,"Cannot set up conversion variables for 'lj' units");

  // get index to internal style variables. create, if needed.
  // default to conversion factor 1.0 for newly created variables

  // @@@@@@ below needs to be updated to use Ellad's luc.
  int v_length, v_mass, v_time;
  char *args[3];
  args[1] = (char *)"internal";
  args[2] = (char *)"1.0";

  args[0] = (char *)"_u_length";
  v_length = variable->find(args[0]);
  if (v_length < 0) {
    variable->set(3,args);
    v_length = variable->find(args[0]);
  }

  args[0] = (char *)"_u_mass";
  v_mass = variable->find(args[0]);
  if (v_mass < 0) {
    variable->set(3,args);
    v_mass = variable->find(args[0]);
  }

  args[0] = (char *)"_u_time";
  v_time = variable->find(args[0]);
  if (v_time < 0) {
    variable->set(3,args);
    v_time = variable->find(args[0]);
  }

  // special case: both unit styles are the same => conversion factor 1.0

  if (strcmp(from,to) == 0) {
    variable->internal_set(v_length,1.0);
    variable->internal_set(v_mass,1.0);
    variable->internal_set(v_time,1.0);
    return;
  }

  if (strcmp(from,"real") == 0) {
    if (strcmp(to,"metal") == 0) {
      variable->internal_set(v_length,1.0);
      variable->internal_set(v_mass,1.0);
      variable->internal_set(v_time,0.001);
    } else {
      std::string err("Do not know how to set up conversion variables ");
      err += "between '";
      err += from;
      err += "' and '";
      err += to;
      err += "' units";
      error->all(FLERR,err.c_str());
    }
  } else if (strcmp(from,"metal") == 0) {
    if (strcmp(to,"real") == 0) {
      variable->internal_set(v_length,1.0);
      variable->internal_set(v_mass,1.0);
      variable->internal_set(v_time,1000.0);
    } else {
      std::string err("Do not know how to set up conversion variables ");
      err += "between '";
      err += from;
      err += "' and '";
      err += to;
      err += "' units";
      error->all(FLERR,err.c_str());
    }
  } else {
    std::string err("Do not know how to set up conversion variables ");
    err += "between '";
    err += from;
    err += "' and '";
    err += to;
    err += "' units";
    error->all(FLERR,err.c_str());
  }
}
