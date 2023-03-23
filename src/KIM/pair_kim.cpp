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
   Contributing authors: Ryan S. Elliott (UMinn)
                         Axel Kohlmeyer (Temple U)
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
   Designed for use with the kim-api-2.0.2 (and newer) package
------------------------------------------------------------------------- */
#include "pair_kim.h"
#include "kim_init.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cstdlib>
#include <cstring>
#include <vector>

using namespace LAMMPS_NS;

static constexpr const char *const cite_openkim =
  "OpenKIM Project: doi:10.1007/s11837-011-0102-6\n\n"
  "@Article{tadmor:elliott:2011,\n"
  " author = {E. B. Tadmor and R. S. Elliott and J. P. Sethna and R. E. Miller "
  "and C. A. Becker},\n"
  " title = {The potential of atomistic simulations and the {K}nowledgebase of "
  "{I}nteratomic {M}odels},\n"
  " journal = {{JOM}},\n"
  " year =    2011,\n"
  " volume =  63,\n"
  " number =  17,\n"
  " pages =   {17},\n"
  " doi =     {10.1007/s11837-011-0102-6}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairKIM::PairKIM(LAMMPS *lmp) :
  Pair(lmp),
  settings_call_count(0),
  init_style_call_count(0),
  kim_modelname(nullptr),
  lmps_map_species_to_unique(nullptr),
  lmps_unique_elements(nullptr),
  lmps_num_unique_elements(0),
  lmps_units(METAL),
  lengthUnit(KIM_LENGTH_UNIT_unused),
  energyUnit(KIM_ENERGY_UNIT_unused),
  chargeUnit(KIM_CHARGE_UNIT_unused),
  temperatureUnit(KIM_TEMPERATURE_UNIT_unused),
  timeUnit(KIM_TIME_UNIT_unused),
  pkim(nullptr),
  pargs(nullptr),
  kim_model_support_for_energy(KIM_SUPPORT_STATUS_notSupported),
  kim_model_support_for_forces(KIM_SUPPORT_STATUS_notSupported),
  kim_model_support_for_particleEnergy(KIM_SUPPORT_STATUS_notSupported),
  kim_model_support_for_particleVirial(KIM_SUPPORT_STATUS_notSupported),
  lmps_local_tot_num_atoms(0),
  kim_global_influence_distance(0.0),
  kim_number_of_neighbor_lists(0),
  kim_cutoff_values(nullptr),
  modelWillNotRequestNeighborsOfNoncontributingParticles(nullptr),
  neighborLists(nullptr),
  kim_particle_codes(nullptr),
  lmps_maxalloc(0),
  kim_particleSpecies(nullptr),
  kim_particleContributing(nullptr),
  lmps_stripped_neigh_list(nullptr),
  lmps_stripped_neigh_ptr(nullptr)
{
  // Initialize Pair data members to appropriate values
  single_enable = 0;  // We do not provide the Single() function
  restartinfo = 0;    // We do not write any restart info
  one_coeff = 1;      // We only allow one coeff * * call
  // set to 1, regardless use of fdotr, to avoid ev_set()'s futzing with
  // vflag_global
  no_virial_fdotr_compute = 1;

  // initial values that determine the KIM state (used by kim_free(), etc.)
  kim_init_ok = false;
  kim_particle_codes_ok = false;

  if (lmp->citeme) lmp->citeme->add(cite_openkim);
}

/* ---------------------------------------------------------------------- */

PairKIM::~PairKIM()
{
  // clean up kim_modelname
  if (kim_modelname != nullptr) delete[] kim_modelname;

  // clean up lammps atom species number to unique particle names mapping
  if (lmps_unique_elements)
    for (int i = 0; i < lmps_num_unique_elements; i++)
      delete[] lmps_unique_elements[i];
  delete[] lmps_unique_elements;

  if (kim_particle_codes_ok) {
    delete[] kim_particle_codes;
    kim_particle_codes = nullptr;
    kim_particle_codes_ok = false;
  }

  // clean up local memory used to support KIM interface
  memory->destroy(kim_particleSpecies);
  memory->destroy(kim_particleContributing);
  memory->destroy(lmps_stripped_neigh_list);
  // clean up lmps_stripped_neigh_ptr
  if (lmps_stripped_neigh_ptr) {
    delete[] lmps_stripped_neigh_ptr;
    lmps_stripped_neigh_ptr = nullptr;
  }

  // clean up allocated memory for standard Pair class usage
  // also, we allocate lmps_map_species_to_uniuqe in the allocate() function
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete[] lmps_map_species_to_unique;
    lmps_map_species_to_unique = nullptr;
  }

  // clean up neighborlist pointers
  if (neighborLists) {
    delete[] neighborLists;
    neighborLists = nullptr;
  }

  // clean up KIM interface (if necessary)
  kim_free();
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_contributing()
{
  int const nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; ++i)
    kim_particleContributing[i] = ( (i < atom->nlocal) ? 1 : 0 );
}

/* ---------------------------------------------------------------------- */

void PairKIM::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  // grow kim_particleSpecies and kim_particleContributing array if necessary
  // needs to be atom->nmax in length
  if (atom->nmax > lmps_maxalloc) {
    memory->destroy(kim_particleSpecies);
    memory->destroy(kim_particleContributing);

    lmps_maxalloc = atom->nmax;
    memory->create(kim_particleSpecies,lmps_maxalloc,
                   "pair:kim_particleSpecies");
    int kimerror = KIM_ComputeArguments_SetArgumentPointerInteger(pargs,
                     KIM_COMPUTE_ARGUMENT_NAME_particleSpeciesCodes,
                     kim_particleSpecies);
    memory->create(kim_particleContributing,lmps_maxalloc,
                   "pair:kim_particleContributing");
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerInteger(pargs,
                             KIM_COMPUTE_ARGUMENT_NAME_particleContributing,
                             kim_particleContributing);
    if (kimerror)
      error->all(FLERR,"Unable to set KIM particle species codes and/or contributing");
  }

  // kim_particleSpecies = KIM atom species for each LAMMPS atom

  int *species = atom->type;
  int nall = atom->nlocal + atom->nghost;
  int ielement;

  for (int i = 0; i < nall; i++) {
    ielement = lmps_map_species_to_unique[species[i]];
    kim_particleSpecies[i] = kim_particle_codes[ielement];
  }

  // Set kim contributing flags
  set_contributing();

  // pass current atom pointers to KIM
  set_argument_pointers();

  // set number of particles
  lmps_local_tot_num_atoms = (int) nall;

  // compute via KIM model
  int kimerror = KIM_Model_Compute(pkim, pargs);
  if (kimerror) error->all(FLERR,"KIM Compute returned error {}", kimerror);

  // compute virial before reverse comm!
  if (vflag_global)
    virial_fdotr_compute();

  // if newton is off, perform reverse comm
  if (!lmps_using_newton) {
    comm->reverse_comm(this);
  }

  if ((vflag_atom != 0) &&
      KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                 KIM_SUPPORT_STATUS_notSupported)) {
    // flip sign and order of virial if KIM is computing it
    double tmp;
    for (int i = 0; i < nall; ++i) {
      for (int j = 0; j < 3; ++j) vatom[i][j] = -1.0*vatom[i][j];
      tmp = vatom[i][3];
      vatom[i][3] = -vatom[i][5];
      vatom[i][4] = -vatom[i][4];
      vatom[i][5] = -tmp;
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairKIM::allocate()
{
  int n = atom->ntypes;

  // allocate standard Pair class arrays
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  // allocate mapping array
  lmps_map_species_to_unique = new int[n+1];

  allocated = 1;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairKIM::settings(int narg, char **arg)
{
  // This is called when "pair_style kim ..." is read from input
  // may be called multiple times
  ++settings_call_count;
  init_style_call_count = 0;

  // arg[0] is the KIM Model name
  if (narg == 0) utils::missing_cmd_args(FLERR, "pair_style kim", error);
  if (narg > 1) {
    const std::string arg_str(arg[1]);
    if ((arg_str == "KIMvirial") || (arg_str == "LAMMPSvirial")) {
      error->all(FLERR,"'KIMvirial' or 'LAMMPSvirial' not supported with kim-api");
    } else error->all(FLERR,"Unknown pair_style kim keyword: {}", arg_str);
  }

  lmps_using_molecular = (atom->molecular > 0);

  // ensure we are in a clean state for KIM (needed on repeated call)
  // first time called will do nothing...
  kim_free();

  // set lmps_* bool flags
  set_lmps_flags();

  // set KIM Model name
  if (kim_modelname != nullptr) {
    delete[] kim_modelname;
    kim_modelname = nullptr;
  }
  kim_modelname = utils::strdup(arg[0]);

  // initialize KIM Model
  kim_init();

  // add citation
  KimInit::write_log_cite(lmp, KimInit::MO, kim_modelname);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairKIM::coeff(int narg, char **arg)
{
  // This is called when "pair_coeff ..." is read from input
  // may be called multiple times

  int i,j;

  if (!allocated) allocate();

  if (narg < 2 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom species to KIM elements
  // lmps_map_species_to_unique[i] =
  // which element the Ith atom type is
  // lmps_num_unique_elements = # of unique elements
  // lmps_unique_elements = list of element names

  // if called multiple times: update lmps_unique_elements
  if (lmps_unique_elements) {
    for (i = 0; i < lmps_num_unique_elements; i++)
      delete[] lmps_unique_elements[i];
    delete[] lmps_unique_elements;
  }
  lmps_unique_elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) lmps_unique_elements[i] = nullptr;

  // Assume all species arguments are valid
  // errors will be detected by below
  atom_type_list.clear();
  lmps_num_unique_elements = 0;
  for (i = 2; i < 2 + atom->ntypes; i++) {
    atom_type_list += std::string(" ") + arg[i];
    for (j = 0; j < lmps_num_unique_elements; j++)
      if (strcmp(arg[i],lmps_unique_elements[j]) == 0) break;
    lmps_map_species_to_unique[i-1] = j;
    if (j == lmps_num_unique_elements) {
      lmps_unique_elements[j] = utils::strdup(arg[i]);
      lmps_num_unique_elements++;
    }
  }

  int count = 0;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (lmps_map_species_to_unique[i] >= 0 &&
          lmps_map_species_to_unique[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  // setup mapping between LAMMPS unique elements and KIM species codes
  if (kim_particle_codes_ok) {
    delete[] kim_particle_codes;
    kim_particle_codes = nullptr;
    kim_particle_codes_ok = false;
  }
  kim_particle_codes = new int[lmps_num_unique_elements];
  kim_particle_codes_ok = true;

  for (int i = 0; i < lmps_num_unique_elements; i++) {
    int supported;
    int code;
    KIM_Model_GetSpeciesSupportAndCode(
      pkim,
      KIM_SpeciesName_FromString(lmps_unique_elements[i]),
      &supported,
      &code);
    if (supported) {
      kim_particle_codes[i] = code;
    } else {
      error->all(FLERR,"GetSpeciesSupportAndCode: symbol not found: {}",lmps_unique_elements[i]);
    }
  }
  // Set the new values for PM parameters
  if (narg > 2 + atom->ntypes) {
    // Get the number of mutable parameters in the kim model
    int numberOfParameters(0);
    KIM_Model_GetNumberOfParameters(pkim, &numberOfParameters);

    if (!numberOfParameters)
      error->all(FLERR,"Incorrect args for pair coefficients\n"
                 "This model has No mutable parameters");

    int kimerror;

    // Parameter name
    std::string paramname;

    for (int i = 2 + atom->ntypes; i < narg;) {
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
        kimerror = KIM_Model_GetParameterMetadata(pkim, param_index, &kim_DataType,
                                                  &extent, &str_name, &str_desc);
        if (kimerror)
          error->all(FLERR,"KIM GetParameterMetadata returned error {}", kimerror);

        const std::string str_name_str(str_name);
        if (paramname == str_name_str) break;
      }

      if (param_index >= numberOfParameters)
        error->all(FLERR,"Wrong argument for pair coefficients.\n"
                   "This Model does not have the requested '{}' parameter", paramname);

      // Get the index_range for the requested parameter
      int nlbound(0);
      int nubound(0);

      if (i < narg) {
        std::string argtostr(arg[i++]);

        // Check to see if the indices range contains only integer numbers & :
        if (argtostr.find_first_not_of("0123456789:") != std::string::npos)
          error->all(FLERR,"Illegal index_range.\nExpected integer parameter(s) instead "
                     "of '{}' in index_range", argtostr);

        std::string::size_type npos = argtostr.find(':');
        if (npos != std::string::npos) {
          argtostr[npos] = ' ';
          auto words = utils::split_words(argtostr);
          nlbound = atoi(words[0].c_str());
          nubound = atoi(words[1].c_str());

          if ((nubound < 1) || (nubound > extent) || (nlbound < 1) || (nlbound > nubound))
            error->all(FLERR,"Illegal index_range '{}-{}' for '{}' parameter with the extent "
                       "of '{}'", nlbound, nubound, paramname, extent);
        } else {
          nlbound = atoi(argtostr.c_str());

          if ((nlbound < 1) || (nlbound > extent))
            error->all(FLERR,"Illegal index '{}' for '{}' parameter with the extent of '{}'",
                       nlbound, paramname, extent);

          nubound = nlbound;
        }
      } else {
        error->all(FLERR,"Wrong number of arguments for pair coefficients.\n"
                   "Index range after parameter name is mandatory");
      }

      // Parameter values
      if (i + nubound - nlbound < narg) {
        if (KIM_DataType_Equal(kim_DataType, KIM_DATA_TYPE_Double)) {
          for (int j = 0; j < nubound - nlbound + 1; ++j) {
            double const V = utils::numeric(FLERR, arg[i++], true, lmp);
            kimerror = KIM_Model_SetParameterDouble(pkim, param_index, nlbound - 1 + j, V);
            if (kimerror)
              error->all(FLERR,"KIM SetParameterDouble returned error");
          }
        } else if (KIM_DataType_Equal(kim_DataType, KIM_DATA_TYPE_Integer)) {
          for (int j = 0; j < nubound - nlbound + 1; ++j) {
            int const V = utils::inumeric(FLERR, arg[i++], true, lmp);
            kimerror = KIM_Model_SetParameterInteger(pkim, param_index, nlbound - 1 + j, V);
            if (kimerror)
              error->all(FLERR,"KIM SetParameterInteger returned error");
          }
        } else
          error->all(FLERR,"Wrong parameter type to update");
      } else {
        error->all(FLERR,"Wrong number of variable values for pair coefficients.\n"
                   " '{}' values are requested for '{}' parameter",
                   nubound - nlbound + 1, paramname);
      }
    }

    kimerror = KIM_Model_ClearThenRefresh(pkim);
    if (kimerror)
      error->all(FLERR,"KIM KIM_Model_ClearThenRefresh returned error");

    // Update cached quantities that may have changed due to Refresh
    KIM_Model_GetInfluenceDistance(pkim, &kim_global_influence_distance);
    KIM_Model_GetNeighborListPointers(
      pkim,
      &kim_number_of_neighbor_lists,
      &kim_cutoff_values,
      &modelWillNotRequestNeighborsOfNoncontributingParticles);
    if (neighborLists) {
      delete[] neighborLists;
      neighborLists = nullptr;
    }
    neighborLists = new NeighList*[kim_number_of_neighbor_lists];
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairKIM::init_style()
{
  // This is called for each "run ...", "minimize ...", etc. read from input
  ++init_style_call_count;

  if (domain->dimension != 3)
    error->all(FLERR,"PairKIM only works with 3D problems");

  // setup lmps_stripped_neigh_list for neighbors of one atom, if needed
  if (lmps_using_molecular) {
    memory->destroy(lmps_stripped_neigh_list);
    memory->create(lmps_stripped_neigh_list,
                   kim_number_of_neighbor_lists*neighbor->oneatom,
                   "pair:lmps_stripped_neigh_list");
    delete[] lmps_stripped_neigh_ptr;
    lmps_stripped_neigh_ptr = new int*[kim_number_of_neighbor_lists];
    for (int i = 0; i < kim_number_of_neighbor_lists; ++i)
      lmps_stripped_neigh_ptr[i]
        = &(lmps_stripped_neigh_list[i*(neighbor->oneatom)]);
  }

  // make sure comm_reverse expects (at most) 9 values when newton is off
  if (!lmps_using_newton) comm_reverse_off = 9;

  // request full neighbor list
  for (int i = 0; i < kim_number_of_neighbor_lists; ++i) {
    int neighflags = NeighConst::REQ_FULL | NeighConst::REQ_NEWTON_OFF;
    if (!modelWillNotRequestNeighborsOfNoncontributingParticles[i])
      neighflags |= NeighConst::REQ_GHOST;
    auto req = neighbor->add_request(this, neighflags);
    req->set_id(i);

    // set cutoff
    if (kim_cutoff_values[i] <= neighbor->skin)
      error->all(FLERR,"Illegal neighbor request (force cutoff {:.3} <= skin {:.3})",
                 kim_cutoff_values[i], neighbor->skin);
    req->set_cutoff(kim_cutoff_values[i] + neighbor->skin);
  }
  // increment instance_me in case of need to change the neighbor list
  // request settings
  instance_me += 1;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   half or full
------------------------------------------------------------------------- */

void PairKIM::init_list(int id, NeighList *ptr)
{
  neighborLists[id] = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairKIM::init_one(int i, int j)
{
  // This is called once of each (unordered) i,j pair for each
  // "run ...", "minimize ...", etc. read from input

  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return kim_global_influence_distance;
}

/* ---------------------------------------------------------------------- */

int PairKIM::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int const last = first + n;
  if (KIM_SupportStatus_NotEqual(kim_model_support_for_forces,
                                 KIM_SUPPORT_STATUS_notSupported) &&
      (KIM_SupportStatus_Equal(kim_model_support_for_particleVirial,
                               KIM_SUPPORT_STATUS_notSupported) ||
       (vflag_atom == 0))) {
    const double *const fp = &(atom->f[0][0]);
    for (int i = first; i < last; i++) {
      buf[m++] = fp[3*i+0];
      buf[m++] = fp[3*i+1];
      buf[m++] = fp[3*i+2];
    }
    return m;
  }
  // ----------------------------------------------------------------------
  // see Pair::ev_setup & Integrate::ev_set()
  // for values of eflag (0-3) and vflag (0-14)
  // -------------------------------------------------------------------------
  if ((vflag_atom != 0) &&
      KIM_SupportStatus_NotEqual(kim_model_support_for_forces,
                                 KIM_SUPPORT_STATUS_notSupported) &&
      KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                 KIM_SUPPORT_STATUS_notSupported)) {
    const double *const fp = &(atom->f[0][0]);
    const double *const va = &(vatom[0][0]);
    for (int i = first; i < last; i++) {
      buf[m++] = fp[3*i+0];
      buf[m++] = fp[3*i+1];
      buf[m++] = fp[3*i+2];

      buf[m++] = va[6*i+0];
      buf[m++] = va[6*i+1];
      buf[m++] = va[6*i+2];
      buf[m++] = va[6*i+3];
      buf[m++] = va[6*i+4];
      buf[m++] = va[6*i+5];
    }
    return m;
  }
  if ((vflag_atom != 0) &&
      KIM_SupportStatus_Equal(kim_model_support_for_forces,
                              KIM_SUPPORT_STATUS_notSupported) &&
      KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                 KIM_SUPPORT_STATUS_notSupported)) {
    const double *const va = &(vatom[0][0]);
    for (int i = first; i < last; i++) {
      buf[m++] = va[6*i+0];
      buf[m++] = va[6*i+1];
      buf[m++] = va[6*i+2];
      buf[m++] = va[6*i+3];
      buf[m++] = va[6*i+4];
      buf[m++] = va[6*i+5];
    }
    return m;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void PairKIM::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  if (KIM_SupportStatus_NotEqual(kim_model_support_for_forces,
                                 KIM_SUPPORT_STATUS_notSupported) &&
      (KIM_SupportStatus_Equal(kim_model_support_for_particleVirial,
                               KIM_SUPPORT_STATUS_notSupported) ||
       (vflag_atom == 0))) {
    double *const fp = &(atom->f[0][0]);
    for (int i = 0; i < n; i++) {
      const int j = list[i];
      fp[3*j+0]+= buf[m++];
      fp[3*j+1]+= buf[m++];
      fp[3*j+2]+= buf[m++];
    }
    return;
  }
  if ((vflag_atom != 0) &&
      KIM_SupportStatus_NotEqual(kim_model_support_for_forces,
                                 KIM_SUPPORT_STATUS_notSupported) &&
      KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                 KIM_SUPPORT_STATUS_notSupported)) {
    double *const fp = &(atom->f[0][0]);
    double *const va = &(vatom[0][0]);
    for (int i = 0; i < n; i++) {
      const int j = list[i];
      fp[3*j+0]+= buf[m++];
      fp[3*j+1]+= buf[m++];
      fp[3*j+2]+= buf[m++];

      va[j*6+0]+=buf[m++];
      va[j*6+1]+=buf[m++];
      va[j*6+2]+=buf[m++];
      va[j*6+3]+=buf[m++];
      va[j*6+4]+=buf[m++];
      va[j*6+5]+=buf[m++];
    }
    return;
  }
  if ((vflag_atom != 0) &&
      KIM_SupportStatus_Equal(kim_model_support_for_forces,
                              KIM_SUPPORT_STATUS_notSupported) &&
      KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                 KIM_SUPPORT_STATUS_notSupported)) {
    double *const va=&(vatom[0][0]);
    for (int i = 0; i < n; i++) {
      const int j = list[i];
      va[j*6+0]+=buf[m++];
      va[j*6+1]+=buf[m++];
      va[j*6+2]+=buf[m++];
      va[j*6+3]+=buf[m++];
      va[j*6+4]+=buf[m++];
      va[j*6+5]+=buf[m++];
    }
  }
  // do nothing
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairKIM::memory_usage()
{
  double bytes = 2 * lmps_maxalloc * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   KIM-specific interface
------------------------------------------------------------------------- */

int PairKIM::get_neigh(void const * const dataObject,
                       int const numberOfNeighborLists,
                       double const * const cutoffs,
                       int const neighborListIndex,
                       int const particleNumber,
                       int * const numberOfNeighbors,
                       int const ** const neighborsOfParticle)
{
  auto  const Model = reinterpret_cast<PairKIM const *>(dataObject);

  if (numberOfNeighborLists != Model->kim_number_of_neighbor_lists)
    return true;
  for (int i = 0; i < numberOfNeighborLists; ++i) {
    if (Model->kim_cutoff_values[i] < cutoffs[i]) return true;
  }

  // neighborListIndex and particleNumber are validated by KIM API

  // initialize numNeigh
  *numberOfNeighbors = 0;

  NeighList *neiobj = Model->neighborLists[neighborListIndex];

  int *numneigh, **firstneigh;
  numneigh = neiobj->numneigh;     // # of J neighbors for each I atom
  firstneigh = neiobj->firstneigh; // ptr to 1st J int value of each I atom

  *numberOfNeighbors = numneigh[particleNumber];

  // strip off neighbor mask for molecular systems
  if (!Model->lmps_using_molecular)
    *neighborsOfParticle = firstneigh[particleNumber];
  else {
    int n = *numberOfNeighbors;
    int *ptr = firstneigh[particleNumber];
    int *lmps_stripped_neigh_list
      = Model->lmps_stripped_neigh_ptr[neighborListIndex];
    for (int i = 0; i < n; i++)
      lmps_stripped_neigh_list[i] = *(ptr++) & NEIGHMASK;
    *neighborsOfParticle = lmps_stripped_neigh_list;
  }
  return false;
}

/* ---------------------------------------------------------------------- */

void PairKIM::kim_free()
{
  if (kim_init_ok) {
    int kimerror = KIM_Model_ComputeArgumentsDestroy(pkim, &pargs);
    if (kimerror)
      error->all(FLERR,"Unable to destroy Compute Arguments Object");

    KIM_Model_Destroy(&pkim);

    lmps_maxalloc = 0;  // reinitialize member variable
  }
  kim_init_ok = false;
}

/* ---------------------------------------------------------------------- */

void PairKIM::kim_init()
{
  int kimerror;

  // initialize KIM model
  int requestedUnitsAccepted;
  kimerror = KIM_Model_Create(
    KIM_NUMBERING_zeroBased,
    lengthUnit, energyUnit, chargeUnit, temperatureUnit, timeUnit,
    kim_modelname,
    &requestedUnitsAccepted,
    &pkim);
  if (kimerror) error->all(FLERR,"KIM ModelCreate failed");
  else if (!requestedUnitsAccepted)
    error->all(FLERR,"KIM Model did not accept the requested unit system");

  auto logID = fmt::format("{}_Model", comm->me);
  KIM_Model_SetLogID(pkim, logID.c_str());

  // check that the model does not require unknown capabilities
  kimerror = check_for_routine_compatibility();
  if (kimerror)
    error->all(FLERR,"KIM Model requires unknown Routines. Unable to proceed");

  kimerror = KIM_Model_ComputeArgumentsCreate(pkim, &pargs);
  if (kimerror) {
    KIM_Model_Destroy(&pkim);
    error->all(FLERR,"KIM ComputeArgumentsCreate failed");
  } else kim_init_ok = true;

  logID = fmt::format("{}_ComputeArguments", comm->me);
  KIM_ComputeArguments_SetLogID(pargs, logID.c_str());

  // determine KIM Model capabilities (used in this function below)
  set_kim_model_has_flags();

  KIM_Model_GetInfluenceDistance(pkim, &kim_global_influence_distance);
  KIM_Model_GetNeighborListPointers(
    pkim,
    &kim_number_of_neighbor_lists,
    &kim_cutoff_values,
    &modelWillNotRequestNeighborsOfNoncontributingParticles);
  if (neighborLists) {
    delete[] neighborLists;
    neighborLists = nullptr;
  }
  neighborLists = new NeighList*[kim_number_of_neighbor_lists];

  kimerror = KIM_ComputeArguments_SetArgumentPointerInteger(pargs,
    KIM_COMPUTE_ARGUMENT_NAME_numberOfParticles,
    &lmps_local_tot_num_atoms);
  if (kimerror) error->all(FLERR,"Unable to set KIM argument pointer");

  kimerror = KIM_ComputeArguments_SetCallbackPointer(pargs,
    KIM_COMPUTE_CALLBACK_NAME_GetNeighborList,
    KIM_LANGUAGE_NAME_cpp,
    reinterpret_cast<KIM_Function *>(get_neigh),
    reinterpret_cast<void *>(this));

  if (kimerror) error->all(FLERR,"Unable to set KIM call back pointer");
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_argument_pointers()
{
  int kimerror;
  kimerror = KIM_ComputeArguments_SetArgumentPointerDouble(
    pargs, KIM_COMPUTE_ARGUMENT_NAME_coordinates, &(atom->x[0][0]));

  // Set KIM pointer appropriately for Energy
  if (KIM_SupportStatus_NotEqual(kim_model_support_for_energy,
                                 KIM_SUPPORT_STATUS_notSupported)) {
      if (KIM_SupportStatus_Equal(kim_model_support_for_energy,
                                  KIM_SUPPORT_STATUS_required)
        || (eflag_global != 0)) {
        kimerror = kimerror ||
        KIM_ComputeArguments_SetArgumentPointerDouble(
          pargs,KIM_COMPUTE_ARGUMENT_NAME_partialEnergy,&(eng_vdwl));
      } else {
        kimerror = kimerror ||
        KIM_ComputeArguments_SetArgumentPointerDouble(
          pargs,KIM_COMPUTE_ARGUMENT_NAME_partialEnergy,
          static_cast<double *>(nullptr));
      }
  }

  // Set KIM pointer appropriately for particalEnergy
  if (KIM_SupportStatus_Equal(kim_model_support_for_particleEnergy,
                              KIM_SUPPORT_STATUS_required)
      && (eflag_atom == 0)) {
    // reallocate per-atom energy array if necessary
    if (atom->nmax > maxeatom) {
      maxeatom = atom->nmax;
      memory->destroy(eatom);
      memory->create(eatom,comm->nthreads*maxeatom,"pair:eatom");
    }
  }

  if (KIM_SupportStatus_Equal(kim_model_support_for_particleEnergy,
                              KIM_SUPPORT_STATUS_optional)
      && (eflag_atom == 0)) {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
      pargs,
      KIM_COMPUTE_ARGUMENT_NAME_partialParticleEnergy,
      static_cast<double *>(nullptr));
  } else if (KIM_SupportStatus_NotEqual(kim_model_support_for_particleEnergy,
                                        KIM_SUPPORT_STATUS_notSupported)) {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
        pargs, KIM_COMPUTE_ARGUMENT_NAME_partialParticleEnergy, eatom);
  }

  // Set KIM pointer appropriately for forces
  if (KIM_SupportStatus_Equal(kim_model_support_for_forces,
                              KIM_SUPPORT_STATUS_notSupported)) {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
      pargs,
      KIM_COMPUTE_ARGUMENT_NAME_partialForces,
      static_cast<double *>(nullptr));
  } else {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
        pargs, KIM_COMPUTE_ARGUMENT_NAME_partialForces, &(atom->f[0][0]));
  }

  // Set KIM pointer appropriately for particleVirial
  if (KIM_SupportStatus_Equal(kim_model_support_for_particleVirial,
                              KIM_SUPPORT_STATUS_required)
      && (vflag_atom == 0)) {
    // reallocate per-atom virial array if necessary
    if (atom->nmax > maxvatom) {
      maxvatom = atom->nmax;
      memory->destroy(vatom);
      memory->create(vatom,comm->nthreads*maxvatom,6,"pair:vatom");
    }
  }

  if (KIM_SupportStatus_Equal(kim_model_support_for_particleVirial,
                              KIM_SUPPORT_STATUS_optional)
      && (vflag_atom == 0)) {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
      pargs,
      KIM_COMPUTE_ARGUMENT_NAME_partialParticleVirial,
      static_cast<double *>(nullptr));
  } else if (KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                        KIM_SUPPORT_STATUS_notSupported)) {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
      pargs, KIM_COMPUTE_ARGUMENT_NAME_partialParticleVirial, &(vatom[0][0]));
  }

  if (kimerror) error->all(FLERR,"Unable to set KIM argument pointers");
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_lmps_flags()
{
  // determint if newton is on or off
  lmps_using_newton = (force->newton_pair == 1);

  // determine if running with pair hybrid
  if (force->pair_match("hybrid",0))
    error->all(FLERR,"Pair style must not be used as a hybrid sub-style");

  const std::string unit_style_str(update->unit_style);

  // determine unit system and set lmps_units flag
  if (unit_style_str == "real") {
    lmps_units = REAL;
    lengthUnit = KIM_LENGTH_UNIT_A;
    energyUnit = KIM_ENERGY_UNIT_kcal_mol;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_fs;
  } else if (unit_style_str == "metal") {
    lmps_units = METAL;
    lengthUnit = KIM_LENGTH_UNIT_A;
    energyUnit = KIM_ENERGY_UNIT_eV;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_ps;
  } else if (unit_style_str == "si") {
    lmps_units = SI;
    lengthUnit = KIM_LENGTH_UNIT_m;
    energyUnit = KIM_ENERGY_UNIT_J;
    chargeUnit = KIM_CHARGE_UNIT_C;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_s;
  } else if (unit_style_str == "cgs") {
    lmps_units = CGS;
    lengthUnit = KIM_LENGTH_UNIT_cm;
    energyUnit = KIM_ENERGY_UNIT_erg;
    chargeUnit = KIM_CHARGE_UNIT_statC;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_s;
  } else if (unit_style_str == "electron") {
    lmps_units = ELECTRON;
    lengthUnit = KIM_LENGTH_UNIT_Bohr;
    energyUnit = KIM_ENERGY_UNIT_Hartree;
    chargeUnit = KIM_CHARGE_UNIT_e;
    temperatureUnit = KIM_TEMPERATURE_UNIT_K;
    timeUnit = KIM_TIME_UNIT_fs;
  } else if ((unit_style_str == "lj") ||
             (unit_style_str == "micro") ||
             (unit_style_str == "nano")) {
    error->all(FLERR,"LAMMPS unit_style {} not supported by KIM models", unit_style_str);
  } else {
    error->all(FLERR,"Unknown unit_style {}", unit_style_str);
  }
}

/* ---------------------------------------------------------------------- */

int PairKIM::check_for_routine_compatibility()
{
  /* Check that we know about all required routines */
  int numberOfModelRoutineNames;
  KIM_MODEL_ROUTINE_NAME_GetNumberOfModelRoutineNames(
      &numberOfModelRoutineNames);
  for (int i = 0; i < numberOfModelRoutineNames; ++i) {
    KIM_ModelRoutineName modelRoutineName;
    KIM_MODEL_ROUTINE_NAME_GetModelRoutineName(i, &modelRoutineName);

    int present;
    int required;
    int error = KIM_Model_IsRoutinePresent(
        pkim, modelRoutineName, &present, &required);
    if (error) return true;

    if (present && required) {
      if (!(KIM_ModelRoutineName_Equal(modelRoutineName, KIM_MODEL_ROUTINE_NAME_Create)
            || KIM_ModelRoutineName_Equal(modelRoutineName,
                                          KIM_MODEL_ROUTINE_NAME_ComputeArgumentsCreate)
            || KIM_ModelRoutineName_Equal(modelRoutineName,
                                          KIM_MODEL_ROUTINE_NAME_Compute)
            || KIM_ModelRoutineName_Equal(modelRoutineName,
                                          KIM_MODEL_ROUTINE_NAME_Refresh)
            || KIM_ModelRoutineName_Equal(modelRoutineName,
                                          KIM_MODEL_ROUTINE_NAME_ComputeArgumentsDestroy)
            || KIM_ModelRoutineName_Equal(modelRoutineName,
                                          KIM_MODEL_ROUTINE_NAME_Destroy))) {
        return true;
      }
    }
  }

  /* everything is good */
  return false;
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_kim_model_has_flags()
{
  int numberOfComputeArgumentNames;
  KIM_COMPUTE_ARGUMENT_NAME_GetNumberOfComputeArgumentNames(
    &numberOfComputeArgumentNames);
  for (int i = 0; i < numberOfComputeArgumentNames; ++i) {
    KIM_ComputeArgumentName computeArgumentName;
    KIM_COMPUTE_ARGUMENT_NAME_GetComputeArgumentName(
      i, &computeArgumentName);
    KIM_SupportStatus supportStatus;
    KIM_ComputeArguments_GetArgumentSupportStatus(
      pargs, computeArgumentName, &supportStatus);

    if (KIM_ComputeArgumentName_Equal(computeArgumentName, KIM_COMPUTE_ARGUMENT_NAME_partialEnergy))
      kim_model_support_for_energy = supportStatus;
    else if (KIM_ComputeArgumentName_Equal(
               computeArgumentName, KIM_COMPUTE_ARGUMENT_NAME_partialForces))
      kim_model_support_for_forces = supportStatus;
    else if (KIM_ComputeArgumentName_Equal(
               computeArgumentName,
               KIM_COMPUTE_ARGUMENT_NAME_partialParticleEnergy))
      kim_model_support_for_particleEnergy = supportStatus;
    else if (KIM_ComputeArgumentName_Equal(
               computeArgumentName,
               KIM_COMPUTE_ARGUMENT_NAME_partialParticleVirial))
      kim_model_support_for_particleVirial = supportStatus;
    else if (KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_required)) {
      error->all(FLERR, "KIM Model requires unsupported compute argument: {}",
                 KIM_ComputeArgumentName_ToString(computeArgumentName));
    }
  }

  if (comm->me == 0) {
    if (KIM_SupportStatus_Equal(kim_model_support_for_energy, KIM_SUPPORT_STATUS_notSupported))
      error->warning(FLERR,"KIM Model does not provide 'partialEnergy'; "
                     "Potential energy will be zero");

    if (KIM_SupportStatus_Equal(kim_model_support_for_forces,
                                KIM_SUPPORT_STATUS_notSupported))
      error->warning(FLERR,"KIM Model does not provide 'partialForce'; Forces will be zero");

    if (KIM_SupportStatus_Equal(kim_model_support_for_particleEnergy,
                                KIM_SUPPORT_STATUS_notSupported))
      error->warning(FLERR,"KIM Model does not provide 'partialParticleEnergy'; "
                     "energy per atom will be zero");

    if (KIM_SupportStatus_Equal(kim_model_support_for_particleVirial,
                                KIM_SUPPORT_STATUS_notSupported))
      error->warning(FLERR,"KIM Model does not provide 'partialParticleVirial'; "
                     "virial per atom will be zero");
  }

  int numberOfComputeCallbackNames;
  KIM_COMPUTE_CALLBACK_NAME_GetNumberOfComputeCallbackNames(
    &numberOfComputeCallbackNames);
  for (int i = 0; i < numberOfComputeCallbackNames; ++i) {
    KIM_ComputeCallbackName computeCallbackName;
    KIM_COMPUTE_CALLBACK_NAME_GetComputeCallbackName(i, &computeCallbackName);
    KIM_SupportStatus supportStatus;
    KIM_ComputeArguments_GetCallbackSupportStatus(pargs,
                                                  computeCallbackName,
                                                  &supportStatus);

    if (KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_required))
      error->all(FLERR,"KIM Model requires unsupported compute callback");
  }
}

KIM_Model *PairKIM::get_kim_model() { return pkim; }

std::string PairKIM::get_atom_type_list() { return atom_type_list; }
