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
   Contributing authors: Ryan S. Elliott (UMinn)
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

#include <cstring>
#include <cstdlib>

// includes from LAMMPS
#include "pair_kim.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "memory.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairKIM::PairKIM(LAMMPS *lmp) :
   Pair(lmp),
   settings_call_count(0),
   init_style_call_count(0),
   kim_modelname(NULL),
   lmps_map_species_to_unique(NULL),
   lmps_unique_elements(NULL),
   lmps_num_unique_elements(0),
   lmps_units(METAL),
   lengthUnit(KIM_LENGTH_UNIT_unused),
   energyUnit(KIM_ENERGY_UNIT_unused),
   chargeUnit(KIM_CHARGE_UNIT_unused),
   temperatureUnit(KIM_TEMPERATURE_UNIT_unused),
   timeUnit(KIM_TIME_UNIT_unused),
   pkim(NULL),
   pargs(NULL),
   kim_model_support_for_energy(KIM_SUPPORT_STATUS_notSupported),
   kim_model_support_for_forces(KIM_SUPPORT_STATUS_notSupported),
   kim_model_support_for_particleEnergy(KIM_SUPPORT_STATUS_notSupported),
   kim_model_support_for_particleVirial(KIM_SUPPORT_STATUS_notSupported),
   lmps_local_tot_num_atoms(0),
   kim_global_influence_distance(0.0),
   kim_number_of_neighbor_lists(0),
   kim_cutoff_values(NULL),
   modelWillNotRequestNeighborsOfNoncontributingParticles(NULL),
   neighborLists(NULL),
   kim_particle_codes(NULL),
   lmps_maxalloc(0),
   kim_particleSpecies(NULL),
   kim_particleContributing(NULL),
   lmps_stripped_neigh_list(NULL),
   lmps_stripped_neigh_ptr(NULL)
{
   // Initialize Pair data members to appropriate values
   single_enable = 0;  // We do not provide the Single() function
   restartinfo = 0;    // We do not write any restart info
   one_coeff = 1;      // We only allow one coeff * * call
   // set to 1, regardless use of fdotr, to avoid ev_set()'s futzing with
   // vflag_global
   no_virial_fdotr_compute = 1;

   // BEGIN: initial values that determine the KIM state
   // (used by kim_free(), etc.)
   kim_init_ok = false;
   kim_particle_codes_ok = false;
   // END

   return;
}

/* ---------------------------------------------------------------------- */

PairKIM::~PairKIM()
{
   // clean up kim_modelname
   if (kim_modelname != 0) delete [] kim_modelname;

   // clean up lammps atom species number to unique particle names mapping
   if (lmps_unique_elements)
      for (int i = 0; i < lmps_num_unique_elements; i++)
        delete [] lmps_unique_elements[i];
   delete [] lmps_unique_elements;

   if (kim_particle_codes_ok)
   {
      delete [] kim_particle_codes;
      kim_particle_codes = NULL;
      kim_particle_codes_ok = false;
   }

   // clean up local memory used to support KIM interface
   memory->destroy(kim_particleSpecies);
   memory->destroy(kim_particleContributing);
   memory->destroy(lmps_stripped_neigh_list);
   // clean up lmps_stripped_neigh_ptr
   if (lmps_stripped_neigh_ptr)
   {
     delete [] lmps_stripped_neigh_ptr;
     lmps_stripped_neigh_ptr = 0;
   }

   // clean up allocated memory for standard Pair class usage
   // also, we allocate lmps_map_species_to_uniuqe in the allocate() function
   if (allocated) {
      memory->destroy(setflag);
      memory->destroy(cutsq);
      delete [] lmps_map_species_to_unique;
   }

   // clean up neighborlist pointers
   if (neighborLists)
   {
     delete [] neighborLists;
     neighborLists = 0;
   }

   // clean up KIM interface (if necessary)
   kim_free();

   return;
}

/* ---------------------------------------------------------------------- */
void PairKIM::set_contributing()
{
  int const nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; ++i)
  {
    kim_particleContributing[i] = ( (i < atom->nlocal) ? 1 : 0 );
  }
}

/* ---------------------------------------------------------------------- */

void PairKIM::compute(int eflag , int vflag)
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
      kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerInteger(
          pargs,
          KIM_COMPUTE_ARGUMENT_NAME_particleContributing,
          kim_particleContributing);
      if (kimerror)
        error->all(
            FLERR,
            "Unable to set KIM particle species codes and/or contributing");
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
   if (kimerror) error->all(FLERR,"KIM Compute returned error");

   // compute virial before reverse comm!
   if (vflag_global)
   {
     virial_fdotr_compute();
   }

   // if newton is off, perform reverse comm
   if (!lmps_using_newton)
   {
     comm->reverse_comm_pair(this);
   }

   if ((vflag_atom) &&
       KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                  KIM_SUPPORT_STATUS_notSupported)
      )
   {  // flip sign and order of virial if KIM is computing it
      double tmp;
      for (int i = 0; i < nall; ++i)
      {
         for (int j = 0; j < 3; ++j) vatom[i][j] = -1.0*vatom[i][j];
         tmp = vatom[i][3];
         vatom[i][3] = -vatom[i][5];
         vatom[i][4] = -vatom[i][4];
         vatom[i][5] = -tmp;
      }
   }

   return;
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

   return;
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

   if (narg != 1)
   {
     if ((narg > 0) && ((0 == strcmp("KIMvirial", arg[0])) ||
                        (0 == strcmp("LAMMPSvirial", arg[0]))))
     {
       error->all(FLERR,"'KIMvirial' or 'LAMMPSvirial' not supported with "
                  "kim-api.");
     }
     else
       error->all(FLERR,"Illegal pair_style command");
   }
   // arg[0] is the KIM Model name

   lmps_using_molecular = (atom->molecular > 0);

   // ensure we are in a clean state for KIM (needed on repeated call)
   // first time called will do nothing...
   kim_free();

   // make sure things are allocated
   if (allocated != 1) allocate();

   // clear setflag to ensure coeff() is called after settings()
   int n = atom->ntypes;
   for (int i = 1; i <= n; i++)
      for (int j = i; j <= n; j++)
         setflag[i][j] = 0;

    // set lmps_* bool flags
   set_lmps_flags();

   // set KIM Model name
   int nmlen = strlen(arg[0]);
   if (kim_modelname != 0)
   {
      delete [] kim_modelname;
      kim_modelname = 0;
   }
   kim_modelname = new char[nmlen+1];
   strcpy(kim_modelname, arg[0]);

   // initialize KIM Model
   kim_init();

   return;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairKIM::coeff(int narg, char **arg)
{
   // This is called when "pair_coeff ..." is read from input
   // may be called multiple times

   int i,j,n;

   if (!allocated) allocate();

   if (narg != 2 + atom->ntypes)
      error->all(FLERR,"Incorrect args for pair coefficients");

   // insure I,J args are * *

   if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
     error->all(FLERR,"Incorrect args for pair coefficients");


   int ilo,ihi,jlo,jhi;
   force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
   force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

   // read args that map atom species to KIM elements
   // lmps_map_species_to_unique[i] =
   // which element the Ith atom type is
   // lmps_num_unique_elements = # of unique elements
   // lmps_unique_elements = list of element names

   // if called multiple times: update lmps_unique_elements
   if (lmps_unique_elements) {
      for (i = 0; i < lmps_num_unique_elements; i++)
        delete [] lmps_unique_elements[i];
      delete [] lmps_unique_elements;
   }
   lmps_unique_elements = new char*[atom->ntypes];
   for (i = 0; i < atom->ntypes; i++) lmps_unique_elements[i] = 0;


   // Assume all species arguments are valid
   // errors will be detected by below
   lmps_num_unique_elements = 0;
   for (i = 2; i < narg; i++) {
     for (j = 0; j < lmps_num_unique_elements; j++)
       if (strcmp(arg[i],lmps_unique_elements[j]) == 0) break;
     lmps_map_species_to_unique[i-1] = j;
     if (j == lmps_num_unique_elements) {
       n = strlen(arg[i]) + 1;
       lmps_unique_elements[j] = new char[n];
       strcpy(lmps_unique_elements[j],arg[i]);
       lmps_num_unique_elements++;
     }
   }

   int count = 0;
   for (int i = ilo; i <= ihi; i++) {
     for (int j = MAX(jlo,i); j <= jhi; j++) {
       if (lmps_map_species_to_unique[i] >= 0 &&
           lmps_map_species_to_unique[j] >= 0) {
         setflag[i][j] = 1;
         count++;
       }
     }
   }

   if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

   // setup mapping between LAMMPS unique elements and KIM species codes
   if (kim_particle_codes_ok)
   {
     delete [] kim_particle_codes;
     kim_particle_codes = NULL;
     kim_particle_codes_ok = false;
   }
   kim_particle_codes = new int[lmps_num_unique_elements];
   kim_particle_codes_ok = true;
   for(int i = 0; i < lmps_num_unique_elements; i++){
      int supported;
      int code;
      KIM_Model_GetSpeciesSupportAndCode(
          pkim,
          KIM_SpeciesName_FromString(lmps_unique_elements[i]),
          &supported,
          &code);
      if (supported)
        kim_particle_codes[i] = code;
      else
      {
        std::stringstream msg;
        msg << "create_kim_particle_codes: symbol not found: "
            << lmps_unique_elements[i];
        error->all(FLERR, msg.str().c_str());
      }
   }

   return;
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
      delete [] lmps_stripped_neigh_ptr;
      lmps_stripped_neigh_ptr = new int*[kim_number_of_neighbor_lists];
      for (int i = 0; i < kim_number_of_neighbor_lists; ++i)
      {
        lmps_stripped_neigh_ptr[i]
            = &(lmps_stripped_neigh_list[i*(neighbor->oneatom)]);
      }

   }

   // make sure comm_reverse expects (at most) 9 values when newton is off
   if (!lmps_using_newton) comm_reverse_off = 9;

   // request full neighbor
   for (int i = 0; i < kim_number_of_neighbor_lists; ++i)
   {
     int irequest = neighbor->request(this,instance_me);
     neighbor->requests[irequest]->id = i;
     neighbor->requests[irequest]->half = 0;
     neighbor->requests[irequest]->full = 1;

     if (modelWillNotRequestNeighborsOfNoncontributingParticles[i])
     {
       neighbor->requests[irequest]->ghost = 0;
     }
     else
     {
       neighbor->requests[irequest]->ghost = 1;

     }
     // always want all owned/ghost pairs
     neighbor->requests[irequest]->newton = 2;
     // set cutoff
     neighbor->requests[irequest]->cut = 1;
     neighbor->requests[irequest]->cutoff
         = kim_cutoff_values[i] + neighbor->skin;
   }

   return;
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
   int i,m,last;
   double *fp;
   fp = &(atom->f[0][0]);

   m = 0;
   last = first + n;
   if (KIM_SupportStatus_NotEqual(kim_model_support_for_forces,
                                  KIM_SUPPORT_STATUS_notSupported)
       &&
       ((vflag_atom == 0) ||
        KIM_SupportStatus_Equal(kim_model_support_for_particleVirial,
                                KIM_SUPPORT_STATUS_notSupported)))
   {
      for (i = first; i < last; i++)
      {
         buf[m++] = fp[3*i+0];
         buf[m++] = fp[3*i+1];
         buf[m++] = fp[3*i+2];
      }
      return m;
   }
   else if (KIM_SupportStatus_NotEqual(kim_model_support_for_forces,
                                       KIM_SUPPORT_STATUS_notSupported) &&
            (vflag_atom == 1) &&
            KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                       KIM_SUPPORT_STATUS_notSupported))
   {
      double *va=&(vatom[0][0]);
      for (i = first; i < last; i++)
      {
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
   else if (KIM_SupportStatus_Equal(kim_model_support_for_forces,
                                    KIM_SUPPORT_STATUS_notSupported)
            &&
            (vflag_atom == 1) &&
            KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                       KIM_SUPPORT_STATUS_notSupported))
   {
      double *va=&(vatom[0][0]);
      for (i = first; i < last; i++)
      {
         buf[m++] = va[6*i+0];
         buf[m++] = va[6*i+1];
         buf[m++] = va[6*i+2];
         buf[m++] = va[6*i+3];
         buf[m++] = va[6*i+4];
         buf[m++] = va[6*i+5];
      }
      return m;
   }
   else
      return 0;
}

/* ---------------------------------------------------------------------- */

void PairKIM::unpack_reverse_comm(int n, int *list, double *buf)
{
   int i,j,m;
   double *fp;
   fp = &(atom->f[0][0]);

   m = 0;
   if (KIM_SupportStatus_NotEqual(kim_model_support_for_forces,
                                  KIM_SUPPORT_STATUS_notSupported)
       &&
       ((vflag_atom == 0) ||
        KIM_SupportStatus_Equal(kim_model_support_for_particleVirial,
                                KIM_SUPPORT_STATUS_notSupported)))
   {
      for (i = 0; i < n; i++)
      {
         j = list[i];
         fp[3*j+0]+= buf[m++];
         fp[3*j+1]+= buf[m++];
         fp[3*j+2]+= buf[m++];
      }
   }
   else if (KIM_SupportStatus_NotEqual(kim_model_support_for_forces,
                                       KIM_SUPPORT_STATUS_notSupported)
            &&
            (vflag_atom == 1) &&
            KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                       KIM_SUPPORT_STATUS_notSupported))
   {
      double *va=&(vatom[0][0]);
      for (i = 0; i < n; i++)
      {
         j = list[i];
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
   }
   else if (KIM_SupportStatus_Equal(kim_model_support_for_forces,
                                    KIM_SUPPORT_STATUS_notSupported)
            &&
            (vflag_atom == 1) &&
            KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                       KIM_SUPPORT_STATUS_notSupported))
   {
      double *va=&(vatom[0][0]);
      for (i = 0; i < n; i++)
      {
         j = list[i];
         va[j*6+0]+=buf[m++];
         va[j*6+1]+=buf[m++];
         va[j*6+2]+=buf[m++];
         va[j*6+3]+=buf[m++];
         va[j*6+4]+=buf[m++];
         va[j*6+5]+=buf[m++];
      }
   } else {
     ; // do nothing
   }

   return;
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
                       int const neighborListIndex, int const particleNumber,
                       int * const numberOfNeighbors,
                       int const ** const neighborsOfParticle)
{
   PairKIM const * const Model
       = reinterpret_cast<PairKIM const * const>(dataObject);

   if (numberOfNeighborLists != Model->kim_number_of_neighbor_lists)
     return true;
   for (int i = 0; i < numberOfNeighborLists; ++i)
   {
     if (Model->kim_cutoff_values[i] < cutoffs[i]) return true;
   }

   // neighborListIndex and particleNumber are validated by KIM API

   // initialize numNeigh
   *numberOfNeighbors = 0;

   NeighList * neiobj = Model->neighborLists[neighborListIndex];

   int *numneigh, **firstneigh;
   numneigh = neiobj->numneigh;     // # of J neighbors for each I atom
   firstneigh = neiobj->firstneigh; // ptr to 1st J int value of each I atom

   *numberOfNeighbors = numneigh[particleNumber];

   // strip off neighbor mask for molecular systems
   if (!Model->lmps_using_molecular)
     *neighborsOfParticle = firstneigh[particleNumber];
   else
   {
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
   if (kim_init_ok)
   {
     int kimerror = KIM_Model_ComputeArgumentsDestroy(pkim, &pargs);
     if (kimerror)
       error->all(FLERR,"Unable to destroy Compute Arguments Object");

     KIM_Model_Destroy(&pkim);
   }
   kim_init_ok = false;

   return;
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
   if (kimerror)
     error->all(FLERR,"KIM ModelCreate failed");
   else {
     if (!requestedUnitsAccepted) {
       error->all(FLERR,"KIM Model did not accept the requested unit system");
     }

     // check that the model does not require unknown capabilities
     kimerror = check_for_routine_compatibility();
     if (kimerror)
     {
       error->all(FLERR,
                  "KIM Model requires unknown Routines.  Unable to proceed.");
     }

     kimerror = KIM_Model_ComputeArgumentsCreate(pkim, &pargs);
     if (kimerror)
     {
       KIM_Model_Destroy(&pkim);
       error->all(FLERR,"KIM ComputeArgumentsCreate failed");
     }
     else
     {
       kim_init_ok = true;
     }
   }

   // determine KIM Model capabilities (used in this function below)
   set_kim_model_has_flags();

   KIM_Model_GetInfluenceDistance(pkim, &kim_global_influence_distance);
   KIM_Model_GetNeighborListPointers(
       pkim,
       &kim_number_of_neighbor_lists,
       &kim_cutoff_values,
       &modelWillNotRequestNeighborsOfNoncontributingParticles);
   if (neighborLists)
   {
     delete [] neighborLists;
     neighborLists = 0;
   }
   neighborLists = new NeighList*[kim_number_of_neighbor_lists];

   kimerror = KIM_ComputeArguments_SetArgumentPointerInteger(pargs,
       KIM_COMPUTE_ARGUMENT_NAME_numberOfParticles,
       &lmps_local_tot_num_atoms);
   if (KIM_SupportStatus_NotEqual(kim_model_support_for_energy,
                                  KIM_SUPPORT_STATUS_notSupported))
     kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(pargs,
         KIM_COMPUTE_ARGUMENT_NAME_partialEnergy,
         &(eng_vdwl));

   kimerror = KIM_ComputeArguments_SetCallbackPointer(pargs,
       KIM_COMPUTE_CALLBACK_NAME_GetNeighborList,
       KIM_LANGUAGE_NAME_cpp,
       reinterpret_cast<KIM_Function *>(get_neigh),
       reinterpret_cast<void *>(this));

   if (kimerror)
     error->all(FLERR,"Unable to register KIM pointers");

   return;
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_argument_pointers()
{
  int kimerror;
  kimerror = KIM_ComputeArguments_SetArgumentPointerDouble(
      pargs, KIM_COMPUTE_ARGUMENT_NAME_coordinates, &(atom->x[0][0]));

  // Set KIM pointer appropriately for particalEnergy
  if (KIM_SupportStatus_Equal(kim_model_support_for_particleEnergy,
                              KIM_SUPPORT_STATUS_required)
      && (eflag_atom != 1))
  {
    // reallocate per-atom energy array if necessary
    if (atom->nmax > maxeatom)
    {
      maxeatom = atom->nmax;
      memory->destroy(eatom);
      memory->create(eatom,comm->nthreads*maxeatom,"pair:eatom");
    }
  }
  if (KIM_SupportStatus_Equal(kim_model_support_for_particleEnergy,
                              KIM_SUPPORT_STATUS_optional)
      && (eflag_atom != 1))
  {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
        pargs,
        KIM_COMPUTE_ARGUMENT_NAME_partialParticleEnergy,
        reinterpret_cast<double * const>(NULL));
  }
  else if (KIM_SupportStatus_NotEqual(kim_model_support_for_particleEnergy,
                                      KIM_SUPPORT_STATUS_notSupported))
  {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
        pargs, KIM_COMPUTE_ARGUMENT_NAME_partialParticleEnergy, eatom);
  }

  // Set KIM pointer appropriately for forces
  if (KIM_SupportStatus_Equal(kim_model_support_for_forces,
                              KIM_SUPPORT_STATUS_notSupported))
  {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
        pargs,
        KIM_COMPUTE_ARGUMENT_NAME_partialForces,
        reinterpret_cast<double * const>(NULL));
  }
  else
  {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
        pargs, KIM_COMPUTE_ARGUMENT_NAME_partialForces, &(atom->f[0][0]));
  }

  // Set KIM pointer appropriately for particleVirial
  if (KIM_SupportStatus_Equal(kim_model_support_for_particleVirial,
                              KIM_SUPPORT_STATUS_required)
      && (vflag_atom != 1))
  {
    // reallocate per-atom virial array if necessary
    if (atom->nmax > maxeatom)
    {
      maxvatom = atom->nmax;
      memory->destroy(vatom);
      memory->create(vatom,comm->nthreads*maxvatom,6,"pair:vatom");
    }
  }
  if (KIM_SupportStatus_Equal(kim_model_support_for_particleVirial,
                              KIM_SUPPORT_STATUS_optional)
      && (vflag_atom != 1))
  {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
        pargs,
        KIM_COMPUTE_ARGUMENT_NAME_partialParticleVirial,
        reinterpret_cast<double * const>(NULL));
  }
  else if (KIM_SupportStatus_NotEqual(kim_model_support_for_particleVirial,
                                      KIM_SUPPORT_STATUS_notSupported))
  {
    kimerror = kimerror || KIM_ComputeArguments_SetArgumentPointerDouble(
        pargs, KIM_COMPUTE_ARGUMENT_NAME_partialParticleEnergy, &(vatom[0][0]));
  }

  if (kimerror)
  {
    error->all(FLERR,"Unable to set KIM argument pointers");
  }

  return;
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_lmps_flags()
{
   // determint if newton is on or off
   lmps_using_newton = (force->newton_pair == 1);

   // determine if running with pair hybrid
   if (force->pair_match("hybrid",0))
   {
     error->all(FLERR,"pair_kim does not support hybrid");
   }

   // determine unit system and set lmps_units flag
   if ((strcmp(update->unit_style,"real")==0)) {
     lmps_units = REAL;
     lengthUnit = KIM_LENGTH_UNIT_A;
     energyUnit = KIM_ENERGY_UNIT_kcal_mol;
     chargeUnit = KIM_CHARGE_UNIT_e;
     temperatureUnit = KIM_TEMPERATURE_UNIT_K;
     timeUnit = KIM_TIME_UNIT_fs;
   } else if ((strcmp(update->unit_style,"metal")==0)) {
     lmps_units = METAL;
     lengthUnit = KIM_LENGTH_UNIT_A;
     energyUnit = KIM_ENERGY_UNIT_eV;
     chargeUnit = KIM_CHARGE_UNIT_e;
     temperatureUnit = KIM_TEMPERATURE_UNIT_K;
     timeUnit = KIM_TIME_UNIT_ps;
   } else if ((strcmp(update->unit_style,"si")==0)) {
     lmps_units = SI;
     lengthUnit = KIM_LENGTH_UNIT_m;
     energyUnit = KIM_ENERGY_UNIT_J;
     chargeUnit = KIM_CHARGE_UNIT_C;
     temperatureUnit = KIM_TEMPERATURE_UNIT_K;
     timeUnit = KIM_TIME_UNIT_s;
   } else if ((strcmp(update->unit_style,"cgs")==0)) {
     lmps_units = CGS;
     lengthUnit = KIM_LENGTH_UNIT_cm;
     energyUnit = KIM_ENERGY_UNIT_erg;
     chargeUnit = KIM_CHARGE_UNIT_statC;
     temperatureUnit = KIM_TEMPERATURE_UNIT_K;
     timeUnit = KIM_TIME_UNIT_s;
   } else if ((strcmp(update->unit_style,"electron")==0)) {
     lmps_units = ELECTRON;
     lengthUnit = KIM_LENGTH_UNIT_Bohr;
     energyUnit = KIM_ENERGY_UNIT_Hartree;
     chargeUnit = KIM_CHARGE_UNIT_e;
     temperatureUnit = KIM_TEMPERATURE_UNIT_K;
     timeUnit = KIM_TIME_UNIT_fs;
   } else if ((strcmp(update->unit_style,"lj")==0)) {
     error->all(FLERR,"LAMMPS unit_style lj not supported by KIM models");
   } else {
     error->all(FLERR,"Unknown unit_style");
   }

   return;
}

/* ---------------------------------------------------------------------- */

int PairKIM::check_for_routine_compatibility()
{
  /* Check that we know about all required routines */
  int numberOfModelRoutineNames;
  KIM_MODEL_ROUTINE_NAME_GetNumberOfModelRoutineNames(
      &numberOfModelRoutineNames);
  for (int i = 0; i < numberOfModelRoutineNames; ++i)
  {
    KIM_ModelRoutineName modelRoutineName;
    KIM_MODEL_ROUTINE_NAME_GetModelRoutineName(i, &modelRoutineName);

    int present;
    int required;
    int error = KIM_Model_IsRoutinePresent(
        pkim, modelRoutineName, &present, &required);
    if (error) { return true; }

    if ((present == true) && (required == true))
    {
      if (!(KIM_ModelRoutineName_Equal(modelRoutineName,
                                       KIM_MODEL_ROUTINE_NAME_Create)
            || KIM_ModelRoutineName_Equal(
                   modelRoutineName,
                   KIM_MODEL_ROUTINE_NAME_ComputeArgumentsCreate)
            || KIM_ModelRoutineName_Equal(modelRoutineName,
                                          KIM_MODEL_ROUTINE_NAME_Compute)
            || KIM_ModelRoutineName_Equal(modelRoutineName,
                                          KIM_MODEL_ROUTINE_NAME_Refresh)
            || KIM_ModelRoutineName_Equal(
                   modelRoutineName,
                   KIM_MODEL_ROUTINE_NAME_ComputeArgumentsDestroy)
            || KIM_ModelRoutineName_Equal(modelRoutineName,
                                          KIM_MODEL_ROUTINE_NAME_Destroy)))
      { return true; }
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
  for (int i = 0; i < numberOfComputeArgumentNames; ++i)
  {
    KIM_ComputeArgumentName computeArgumentName;
    KIM_COMPUTE_ARGUMENT_NAME_GetComputeArgumentName(
        i, &computeArgumentName);
    KIM_SupportStatus supportStatus;
    KIM_ComputeArguments_GetArgumentSupportStatus(
        pargs, computeArgumentName, &supportStatus);

    if (KIM_ComputeArgumentName_Equal(computeArgumentName,
                                      KIM_COMPUTE_ARGUMENT_NAME_partialEnergy)
       )
      kim_model_support_for_energy = supportStatus;
    else if (KIM_ComputeArgumentName_Equal(
                 computeArgumentName, KIM_COMPUTE_ARGUMENT_NAME_partialForces)
            )
      kim_model_support_for_forces = supportStatus;
    else if
        (KIM_ComputeArgumentName_Equal(
            computeArgumentName,
            KIM_COMPUTE_ARGUMENT_NAME_partialParticleEnergy)\
        )
      kim_model_support_for_particleEnergy = supportStatus;
    else if
        (KIM_ComputeArgumentName_Equal(
            computeArgumentName,
            KIM_COMPUTE_ARGUMENT_NAME_partialParticleVirial)
        )
      kim_model_support_for_particleVirial = supportStatus;
    else if (KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_required)
            )
    {
      std::stringstream msg;
      msg << "KIM Model requires unsupported compute argument: "
          << KIM_ComputeArgumentName_ToString(computeArgumentName);
      error->all(FLERR, msg.str().c_str());
    }
  }

  if (KIM_SupportStatus_Equal(kim_model_support_for_energy,
                              KIM_SUPPORT_STATUS_notSupported))
    error->warning(FLERR,"KIM Model does not provide `partialEnergy'; "
                   "Potential energy will be zero");

  if (KIM_SupportStatus_Equal(kim_model_support_for_forces,
                              KIM_SUPPORT_STATUS_notSupported))
    error->warning(FLERR,"KIM Model does not provide `partialForce'; "
                   "Forces will be zero");

  if (KIM_SupportStatus_Equal(kim_model_support_for_particleEnergy,
                              KIM_SUPPORT_STATUS_notSupported))
    error->warning(FLERR,"KIM Model does not provide `partialParticleEnergy'; "
                   "energy per atom will be zero");

  if (KIM_SupportStatus_Equal(kim_model_support_for_particleVirial,
                              KIM_SUPPORT_STATUS_notSupported))
    error->warning(FLERR,"KIM Model does not provide `partialParticleVirial'; "
                   "virial per atom will be zero");

  int numberOfComputeCallbackNames;
  KIM_COMPUTE_CALLBACK_NAME_GetNumberOfComputeCallbackNames(
      &numberOfComputeCallbackNames);
  for (int i = 0; i < numberOfComputeCallbackNames; ++i)
  {
    KIM_ComputeCallbackName computeCallbackName;
    KIM_COMPUTE_CALLBACK_NAME_GetComputeCallbackName(
        i, &computeCallbackName);
    KIM_SupportStatus supportStatus;
    KIM_ComputeArguments_GetCallbackSupportStatus(
        pargs, computeCallbackName, &supportStatus);

    if (KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_required))
    {
      error->all(FLERR,"KIM Model requires unsupported compute callback");
    }
  }

  return;
}
