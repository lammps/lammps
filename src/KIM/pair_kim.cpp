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
   Contributing authors: Ryan S. Elliott,
                         Valeriu Smirichinski,
                         Ellad Tadmor
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Designed for use with the kim-api-v1.6.0 (and newer) package
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
   kim_modelname(0),
   lmps_map_species_to_unique(0),
   lmps_unique_elements(0),
   lmps_num_unique_elements(0),
   lmps_units(METAL),
   pkim(0),
   pargs(0),
   kim_particle_codes(0),
   lmps_local_tot_num_atoms(0),
   kim_global_influence_distance(0.0),
   kim_number_of_cutoffs(0),
   kim_cutoff_values(0),
   lmps_maxalloc(0),
   kim_particleSpecies(0),
   kim_particleContributing(0),
   lmps_stripped_neigh_list(0)
{
   // Initialize Pair data members to appropriate values
   single_enable = 0;  // We do not provide the Single() function
   restartinfo = 0;    // We do not write any restart info
   one_coeff = 1;      // We only allow one coeff * * call

   // BEGIN: initial values that determine the KIM state
   // (used by kim_free(), etc.)
   kim_model_init_ok = false;
   kim_init_ok = false;
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

   // clean up local memory used to support KIM interface
   memory->destroy(kim_particleSpecies);
   memory->destroy(kim_particleContributing);
   memory->destroy(lmps_stripped_neigh_list);

   // clean up allocated memory for standard Pair class usage
   // also, we allocate lmps_map_species_to_uniuqe in the allocate() function
   if (allocated) {
      memory->destroy(setflag);
      memory->destroy(cutsq);
      delete [] lmps_map_species_to_unique;
   }

   // clean up KIM interface (if necessary)
   kim_free();

   return;
}

/* ---------------------------------------------------------------------- */

void PairKIM::compute(int eflag , int vflag)
{
   int kimerror;

   if (eflag || vflag)
      ev_setup(eflag,vflag);
   else
      ev_unset();

   // grow kim_particleSpecies and kim_particleContributing array if necessary
   // needs to be atom->nmax in length
   if (atom->nmax > lmps_maxalloc) {
      memory->destroy(kim_particleSpecies);
      memory->destroy(kim_particleContributing);

      lmps_maxalloc = atom->nmax;
      memory->create(kim_particleSpecies,lmps_maxalloc,
                     "pair:kim_particleSpecies");
      int kimerror = pargs->SetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes,
          kim_particleSpecies);
      memory->create(kim_particleContributing,lmps_maxalloc,
                     "pair:kim_particleContributing");
      kimerror = kimerror || pargs->SetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::particleContributing,
          kim_particleContributing);
      if (kimerror)
        error->all(
            FLERR,
            "Unable to set KIM particle species codes and/or contributing");
   }

   // kim_particleSpecies = KIM atom species for each LAMMPS atom
   // set ielement to valid 0 if lmps_map_species_to_unique[] stores an un-used -1

   int *species = atom->type;
   int nall = atom->nlocal + atom->nghost;
   int ielement;

   for (int i = 0; i < nall; i++) {
      ielement = lmps_map_species_to_unique[species[i]];
      ielement = MAX(ielement,0);
      kim_particleSpecies[i] = kim_particle_codes[ielement];

      kim_particleContributing[i] = ( (i<atom->nlocal) ? 1 : 0 );
   }

   // pass current atom pointers to KIM
   set_volatiles();

   // compute via KIM model
   kimerror = pkim->Compute(pargs);
   if (kimerror) error->all(FLERR,"KIM Compute returned error");

   // assemble force and particleVirial if needed
   if (!lmps_using_newton) comm->reverse_comm_pair(this);

   if ((no_virial_fdotr_compute == 1) && (vflag_global))
   {  // flip sign and order of virial if KIM is computing it
      for (int i = 0; i < 3; ++i) virial[i] = -1.0*virial[i];
      double tmp = virial[3];
      virial[3] = -virial[5];
      virial[4] = -virial[4];
      virial[5] = -tmp;
   }
   else
   {  // compute virial via LAMMPS fdotr mechanism
      if (vflag_fdotr) virial_fdotr_compute();
   }

   if ((kim_model_has_particleVirial) && (vflag_atom))
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

   if (narg < 2) error->all(FLERR,"Illegal pair_style command");
   // arg[0] is the virial handling option: "LAMMPSvirial" or "KIMvirial"
   // arg[1] is the KIM Model name
   // arg[2] is the print-kim-file flag: 0/1 do-not/do print (default 0)

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

   // set virial handling
   if (strcmp(arg[0],"LAMMPSvirial") == 0)
   {
      no_virial_fdotr_compute = 0;
   }
   else if (strcmp(arg[0],"KIMvirial") == 0)
   {
      no_virial_fdotr_compute = 1;
   }
   else
   {
      error->all(FLERR,"Unrecognized virial argument in pair_style command");
   }

   // set KIM Model name
   int nmlen = strlen(arg[1]);
   if (kim_modelname != 0)
   {
      delete [] kim_modelname;
      kim_modelname = 0;
   }
   kim_modelname = new char[nmlen+1];
   strcpy(kim_modelname, arg[1]);

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
   // errors will be detected by kim_api_init() matching
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

   int kimerror;
   // KIM and Model initialization (only once)
   // also sets kim_ind_* and kim_* bool flags
   if (!kim_init_ok)
   {
      kim_init();
   }

   // make sure comm_reverse expects (at most) 9 values when newton is off
   if (!lmps_using_newton) comm_reverse_off = 9;

   // request full neighbor list
   int irequest = neighbor->request(this,instance_me);
   neighbor->requests[irequest]->half = 0;
   neighbor->requests[irequest]->full = 1;
   // make sure full lists also include local-ghost pairs
   if (lmps_using_newton) neighbor->requests[irequest]->newton = 0;

   return;
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
   if ((kim_model_has_forces) && ((vflag_atom == 0) ||
                                  (!kim_model_has_particleVirial)))
   {
      for (i = first; i < last; i++)
      {
         buf[m++] = fp[3*i+0];
         buf[m++] = fp[3*i+1];
         buf[m++] = fp[3*i+2];
      }
      return m;
   }
   else if ((kim_model_has_forces) && (vflag_atom == 1) &&
            (kim_model_has_particleVirial))
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
   else if ((!kim_model_has_forces) && (vflag_atom == 1) &&
            (kim_model_has_particleVirial))
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
   if ((kim_model_has_forces) && ((vflag_atom == 0) ||
                                  (!kim_model_has_particleVirial)))
   {
      for (i = 0; i < n; i++)
      {
         j = list[i];
         fp[3*j+0]+= buf[m++];
         fp[3*j+1]+= buf[m++];
         fp[3*j+2]+= buf[m++];
      }
   }
   else if ((kim_model_has_forces) && (vflag_atom == 1) &&
            (kim_model_has_particleVirial))
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
   else if ((!kim_model_has_forces) && (vflag_atom == 1) &&
            (kim_model_has_particleVirial))
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
   }
   else
      ;// do nothing

   return;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairKIM::memory_usage()
{
   double bytes = lmps_maxalloc * sizeof(int);
   return bytes;
}

/* ----------------------------------------------------------------------
   KIM-specific interface
------------------------------------------------------------------------- */

int PairKIM::get_neigh(void const * const dataObject,
                       int const numberOfCutoffs, double const * const cutoffs,
                       int const neighborListIndex, int const particleNumber,
                       int * const numberOfNeighbors,
                       int const ** const neighborsOfParticle)
{
   PairKIM const * const Model
       = reinterpret_cast<PairKIM const * const>(dataObject);

   if ((numberOfCutoffs != 1) || (cutoffs[0] > Model->kim_cutoff_values[0]))
     return true;

   if (neighborListIndex != 0) return true;

   // initialize numNeigh
   *numberOfNeighbors = 0;

   if ((particleNumber >= Model->lmps_local_tot_num_atoms) ||
       (particleNumber < 0)) /* invalid id */
   {
     return true;
   }

   NeighList * neiobj = Model->list;
   int nAtoms = Model->lmps_local_tot_num_atoms;

   int j, jj, inum, *ilist, *numneigh, **firstneigh;
   inum = neiobj->inum;             //# of I atoms neighbors are stored for
   ilist = neiobj->ilist;           //local indices of I atoms
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
     int *lmps_stripped_neigh_list = Model->lmps_stripped_neigh_list;
     for (int i = 0; i < n; i++)
       lmps_stripped_neigh_list[i] = *(ptr++) & NEIGHMASK;
     *neighborsOfParticle = lmps_stripped_neigh_list;
   }
   return false;
}

/* ---------------------------------------------------------------------- */

void PairKIM::kim_free()
{
   int kimerror;

   if (kim_model_init_ok)
   {
     int kimerror = pkim->ComputeArgumentsDestroy(&pargs);
     if (kimerror)
       error->all(FLERR,"Unable to destroy Compute Arguments Object");

     KIM::Model::Destroy(&pkim);
     kim_model_init_ok = false;
   }
   kim_init_ok = false;

   if (kim_particle_codes_ok)
   {
      delete [] kim_particle_codes;
      kim_particle_codes = 0;
      kim_particle_codes_ok = false;
   }

   return;
}

/* ---------------------------------------------------------------------- */

void PairKIM::kim_init()
{
   int kimerror;

   // initialize KIM model
   int requestedUnitsAccepted;
   kimerror = KIM::Model::Create(
       KIM::NUMBERING::zeroBased,
       lengthUnit, energyUnit, chargeUnit, temperatureUnit, timeUnit,
       kim_modelname,
       &requestedUnitsAccepted,
       &pkim);
   if (kimerror)
     error->all(FLERR,"KIM ModelCreate failed");
   else {
     if (!requestedUnitsAccepted) {
       // @@@ error for now.  Fix as needed
       error->all(FLERR,"KIM Model did not accept the requested unit system");
     }

     kimerror = pkim->ComputeArgumentsCreate(&pargs);
     if (kimerror)
       error->all(FLERR,"KIM ComputeArgumentsCreate failed");
     else
       kim_init_ok = true;
   }

   // determine KIM Model capabilities (used in this function below)
   set_kim_model_has_flags();

   // setup mapping between LAMMPS unique elements and KIM species codes
   kim_particle_codes = new int[lmps_num_unique_elements];
   kim_particle_codes_ok = true;
   for(int i = 0; i < lmps_num_unique_elements; i++){
      int kimerror;
      int supported;
      int code;
      kimerror = pkim->GetSpeciesSupportAndCode(
          KIM::SpeciesName(lmps_unique_elements[i]),
          &supported,
          &code);
      if (supported)
        kim_particle_codes[i] = code;
      else
        error->all(FLERR,"create_kim_particle_codes: symbol not found ");
   }

   // set pointer values in KIM API object that will not change during run
   set_statics();

   return;
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_statics()
{
   // set total number of atoms
   lmps_local_tot_num_atoms = (int) (atom->nghost + atom->nlocal);

   pkim->GetInfluenceDistance(&kim_global_influence_distance);
   pkim->GetNeighborListCutoffsPointer(&kim_number_of_cutoffs,
                                       &kim_cutoff_values);

   int kimerror = pargs->SetArgumentPointer(
       KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles,
       &lmps_local_tot_num_atoms);
   if (kim_model_has_energy)
     kimerror = kimerror || pargs->SetArgumentPointer(
         KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
         &(eng_vdwl));

   kimerror = pargs->SetCallbackPointer(
       KIM::COMPUTE_CALLBACK_NAME::GetNeighborList,
       KIM::LANGUAGE_NAME::cpp,
       reinterpret_cast<KIM::func *>(get_neigh),
       reinterpret_cast<void *>(this));

   if (kimerror)
     error->all(FLERR,"Unable to register KIM static pointers");

   return;
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_volatiles()
{
   int kimerror;
   lmps_local_tot_num_atoms = (int) (atom->nghost + atom->nlocal);

   kimerror = pargs->SetArgumentPointer(
       KIM::COMPUTE_ARGUMENT_NAME::coordinates,
       &(atom->x[0][0]));

   if (kim_model_has_particleEnergy && (eflag_atom == 1))
   {
     kimerror = kimerror || pargs->SetArgumentPointer(
         KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
         eatom);
   }

   if (kim_model_has_forces)
   {
     kimerror = kimerror || pargs->SetArgumentPointer(
         KIM::COMPUTE_ARGUMENT_NAME::partialForces,
         &(atom->f[0][0]));
   }

   if (kim_model_has_particleVirial)
   {
     if(vflag_atom != 1) {
       kimerror = kimerror || pargs->SetArgumentPointer(
           KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
           &(vatom[0][0]));
     } else {
       kimerror = kimerror || pargs->SetArgumentPointer(
           KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
           reinterpret_cast<double const * const>(NULL));
     }
   }

   if (no_virial_fdotr_compute == 1)
   {
     if (kim_model_has_virial)
     {
       if (vflag_global == 1)
         kimerror = kimerror || pargs->SetArgumentPointer(
             KIM::COMPUTE_ARGUMENT_NAME::partialVirial,
             &(virial[0]));
       else
         kimerror = kimerror || pargs->SetArgumentPointer(
             KIM::COMPUTE_ARGUMENT_NAME::partialVirial,
             reinterpret_cast<double const * const>(NULL));
     }
   }

   if (kimerror)
   {
     error->all(FLERR,"Unable to set KIM volatile pointers");
   }

   return;
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_lmps_flags()
{
   // determint if newton is on or off
   lmps_using_newton = (force->newton_pair == 1);

   // setup lmps_stripped_neigh_list for neighbors of one atom, if needed
   lmps_using_molecular = (atom->molecular > 0);
   if (lmps_using_molecular) {
      memory->destroy(lmps_stripped_neigh_list);
      memory->create(lmps_stripped_neigh_list,neighbor->oneatom,
                     "pair:lmps_stripped_neigh_list");
   }

   // determine if running with pair hybrid
   if (force->pair_match("hybrid",0))
   {
     error->all(FLERR,"pair_kim does not support hybrid.");
   }

   // determine unit system and set lmps_units flag
   if ((strcmp(update->unit_style,"real")==0)) {
     lmps_units = REAL;
     lengthUnit = KIM::LENGTH_UNIT::A;
     energyUnit = KIM::ENERGY_UNIT::kcal_mol;
     chargeUnit = KIM::CHARGE_UNIT::e;
     temperatureUnit = KIM::TEMPERATURE_UNIT::K;
     timeUnit = KIM::TIME_UNIT::fs;
   } else if ((strcmp(update->unit_style,"metal")==0)) {
     lmps_units = METAL;
     lengthUnit = KIM::LENGTH_UNIT::A;
     energyUnit = KIM::ENERGY_UNIT::eV;
     chargeUnit = KIM::CHARGE_UNIT::e;
     temperatureUnit = KIM::TEMPERATURE_UNIT::K;
     timeUnit = KIM::TIME_UNIT::ps;
   } else if ((strcmp(update->unit_style,"si")==0)) {
     lmps_units = SI;
     lengthUnit = KIM::LENGTH_UNIT::m;
     energyUnit = KIM::ENERGY_UNIT::J;
     chargeUnit = KIM::CHARGE_UNIT::C;
     temperatureUnit = KIM::TEMPERATURE_UNIT::K;
     timeUnit = KIM::TIME_UNIT::s;
   } else if ((strcmp(update->unit_style,"cgs")==0)) {
     lmps_units = CGS;
     lengthUnit = KIM::LENGTH_UNIT::cm;
     energyUnit = KIM::ENERGY_UNIT::erg;
     chargeUnit = KIM::CHARGE_UNIT::statC;
     temperatureUnit = KIM::TEMPERATURE_UNIT::K;
     timeUnit = KIM::TIME_UNIT::s;
   } else if ((strcmp(update->unit_style,"electron")==0)) {
     lmps_units = ELECTRON;
     lengthUnit = KIM::LENGTH_UNIT::Bohr;
     energyUnit = KIM::ENERGY_UNIT::Hartree;
     chargeUnit = KIM::CHARGE_UNIT::e;
     temperatureUnit = KIM::TEMPERATURE_UNIT::K;
     timeUnit = KIM::TIME_UNIT::fs;
   } else if ((strcmp(update->unit_style,"lj")==0)) {
     error->all(FLERR,"LAMMPS unit_style lj not supported by KIM models");
   } else {
     error->all(FLERR,"Unknown unit_style");
   }

   return;
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_kim_model_has_flags()
{
  // @@ the procedure below should be improved to be more comprehensive
  // @@ and ensure that there are no additions/changes to the kim-api
  // @@ that could cause a problem.  This should be done using the
  // @@ "discoverability" features of the kim-api
   int kimerror;
   KIM::SupportStatus supportStatus;

   // determine if the KIM Model can compute the total partialEnergy

   // determine if the KIM Model can compute the energy
   kimerror = pargs->GetArgumentSupportStatus(
       KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
       &supportStatus);
   if (kimerror)
     error->all(FLERR,"Unable to get KIM Support Status");
   if (KIM::SUPPORT_STATUS::notSupported == supportStatus) {
     kim_model_has_energy = false;
     error->warning(FLERR,"KIM Model does not provide `partialEnergy'; "
                    "Potential energy will be zero");
   } else {
     kim_model_has_energy = true;
   }

   // determine if the KIM Model can compute the partialForces
   kimerror = pargs->GetArgumentSupportStatus(
       KIM::COMPUTE_ARGUMENT_NAME::partialForces,
       &supportStatus);
   if (kimerror)
     error->all(FLERR,"Unable to get KIM Support Status");
   if (KIM::SUPPORT_STATUS::notSupported == supportStatus) {
     kim_model_has_forces = false;
     error->warning(FLERR,"KIM Model does not provide `partialForce'; "
                    "Forces will be zero");
   } else {
     kim_model_has_forces = true;
   }

   // determine if the KIM Model can compute the partialVirial
   kimerror = pargs->GetArgumentSupportStatus(
       KIM::COMPUTE_ARGUMENT_NAME::partialVirial,
       &supportStatus);
   if (kimerror)
     error->all(FLERR,"Unable to get KIM Support Status");
   if (KIM::SUPPORT_STATUS::notSupported == supportStatus) {
     kim_model_has_virial = false;
     error->warning(FLERR,"KIM Model does not provide `partialVirial'. "
                    "pair_kim now using `LAMMPSvirial' option.");
     no_virial_fdotr_compute = 0;
   } else {
     kim_model_has_virial = true;
   }

   // determine if the KIM Model can compute the partialParticleEnergy
   kimerror = pargs->GetArgumentSupportStatus(
       KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
       &supportStatus);
   if (kimerror)
     error->all(FLERR,"Unable to get KIM Support Status");
   if (KIM::SUPPORT_STATUS::notSupported == supportStatus) {
     kim_model_has_particleEnergy = false;
     error->warning(FLERR,"KIM Model does not provide `partialParticleEnergy'; "
                    "energy per atom will be zero");
   } else {
     kim_model_has_particleEnergy = true;
   }

   // determine if the KIM Model can compute the partialParticleVirial
   kimerror = pargs->GetArgumentSupportStatus(
       KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
       &supportStatus);
   if (kimerror)
     error->all(FLERR,"Unable to get KIM Support Status");
   if (KIM::SUPPORT_STATUS::notSupported == supportStatus) {
     kim_model_has_particleVirial = false;
     error->warning(FLERR,"KIM Model does not provide `partialParticleVirial'; "
                    "virial per atom will be zero");
   } else {
     kim_model_has_particleVirial = true;
   }

   return;
}
