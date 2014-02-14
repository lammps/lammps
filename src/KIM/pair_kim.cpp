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
   Designed for use with the openkim-api-v1.2.0 (and newer) package
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

// includes from KIM
#include "KIM_API.h"
#include "KIM_API_status.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairKIM::PairKIM(LAMMPS *lmp) :
   Pair(lmp),
   kim_modelname(0),
   lmps_map_types_to_unique(0),
   lmps_unique_elements(0),
   lmps_num_unique_elements(0),
   lmps_units(METAL),
   pkim(0),
   kim_ind_coordinates(-1),
   kim_ind_numberOfParticles(-1),
   kim_ind_numberContributingParticles(-1),
   kim_ind_numberParticleTypes(-1),
   kim_ind_particleTypes(-1),
   kim_ind_get_neigh(-1),
   kim_ind_neighObject(-1),
   kim_ind_cutoff(-1),
   kim_ind_energy(-1),
   kim_ind_particleEnergy(-1),
   kim_ind_forces(-1),
   kim_ind_virial(-1),
   kim_ind_particleVirial(-1),
   kim_particle_codes(0),
   lmps_local_tot_num_atoms(0),   
   kim_global_cutoff(0.0),
   lmps_maxalloc(0),
   kim_particleTypes(0),
   lmps_force_tmp(0),
   lmps_stripped_neigh_list(0),
   kim_iterator_position(0)
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

   // allocate enough memory to ensure we are safe (by using neighbor->oneatom)
   memory->create(Rij,3*(neighbor->oneatom),"pair:Rij");

   return;
}

/* ---------------------------------------------------------------------- */

PairKIM::~PairKIM()
{
   // clean up kim_modelname
   if (kim_modelname != 0) delete [] kim_modelname;

   // clean up lammps atom type number to unique particle names mapping
   if (lmps_unique_elements)
      for (int i = 0; i < lmps_num_unique_elements; i++) 
        delete [] lmps_unique_elements[i];
   delete [] lmps_unique_elements;

   // clean up local memory used to support KIM interface
   memory->destroy(kim_particleTypes);
   memory->destroy(lmps_force_tmp);
   memory->destroy(lmps_stripped_neigh_list);

   // clean up allocated memory for standard Pair class usage
   // also, we allocate lmps_map_types_to_uniuqe in the allocate() function
   if (allocated) {
      memory->destroy(setflag);
      memory->destroy(cutsq);
      delete [] lmps_map_types_to_unique;
   }

   // clean up Rij array
   memory->destroy(Rij);

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

   // grow kim_particleTypes array if necessary
   // needs to be atom->nmax in length
   if (atom->nmax > lmps_maxalloc) {
      memory->destroy(kim_particleTypes);
      memory->destroy(lmps_force_tmp);
      
      lmps_maxalloc = atom->nmax;
      memory->create(kim_particleTypes,lmps_maxalloc,"pair:kim_particleTypes");
      memory->create(lmps_force_tmp,lmps_maxalloc,3,"pair:lmps_force_tmp");
   }

   // kim_particleTypes = KIM atom type for each LAMMPS atom
   // set ielement to valid 0 if lmps_map_types_to_unique[] stores an un-used -1

   int *type = atom->type;
   int nall = atom->nlocal + atom->nghost;
   int ielement;

   for (int i = 0; i < nall; i++) {
      ielement = lmps_map_types_to_unique[type[i]];
      ielement = MAX(ielement,0);
      // @@ this (above line) provides bogus info 
      // @@ (when lmps_map_types_to_unique[type[i]]==-1) to KIM, but I guess
      // @@ this only happens when lmps_hybrid==true, 
      // @@ and we are sure that iterator mode will
      // @@ not use these atoms.... (?)
      kim_particleTypes[i] = kim_particle_codes[ielement];
   }

   // pass current atom pointers to KIM
   set_volatiles();

   pkim->setm_compute_by_index(&kimerror,3*3,
                               kim_ind_particleEnergy, eflag_atom,
                                (int) kim_model_has_particleEnergy,
                               kim_ind_particleVirial, vflag_atom,
                                (int) kim_model_has_particleVirial,
                               kim_ind_virial, vflag_global!=0, 
                                no_virial_fdotr_compute);
   kim_error(__LINE__,"setm_compute_by_index",kimerror);

   // compute via KIM model
   kimerror = pkim->model_compute();
   kim_error(__LINE__,"PairKIM::pkim->model_compute() error",kimerror);
   // assemble force and particleVirial if needed
   if (!lmps_using_newton) comm->reverse_comm_pair(this);

   // sum lmps_force_tmp to f if running in hybrid mode
   if (lmps_hybrid) {
      double **f = atom->f;
      for (int i = 0; i < nall; i++) {
         f[i][0] += lmps_force_tmp[i][0];
         f[i][1] += lmps_force_tmp[i][1];
         f[i][2] += lmps_force_tmp[i][2];
      }
   }

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
   memory->create(cutsq,n+1,n+1,"pair:cutsq");

   // allocate mapping array
   lmps_map_types_to_unique = new int[n+1];

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

   if (narg != 2) error->all(FLERR,"Illegal pair_style command");
   // arg[0] is the virial handling option: "LAMMPSvirial" or "KIMvirial"
   // arg[1] is the KIM Model name

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

   // ensure I,J args are * *

   if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
      error->all(FLERR,"Incorrect args for pair coefficients");

   // read args that map atom types to KIM elements
   // lmps_map_types_to_unique[i] = 
   // which element the Ith atom type is, -1 if NULL
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

   lmps_num_unique_elements = 0;
   for (i = 2; i < narg; i++) {
      if (strcmp(arg[i],"NULL") == 0) {
         if (!lmps_hybrid) 
           error->all(FLERR,"Invalid args for non-hybrid pair coefficients");
         lmps_map_types_to_unique[i-1] = -1;
         continue;
      }
      for (j = 0; j < lmps_num_unique_elements; j++)
         if (strcmp(arg[i],lmps_unique_elements[j]) == 0) break;
      lmps_map_types_to_unique[i-1] = j;
      if (j == lmps_num_unique_elements) {
         n = strlen(arg[i]) + 1;
         lmps_unique_elements[j] = new char[n];
         strcpy(lmps_unique_elements[j],arg[i]);
         lmps_num_unique_elements++;
      }
   }

   // clear setflag since coeff() called once with I,J = * *
   n = atom->ntypes;
   for (int i = 1; i <= n; i++)
      for (int j = i; j <= n; j++)
         setflag[i][j] = 0;

   // set setflag i,j for type pairs where both are mapped to elements
   int count = 0;
   for (int i = 1; i <= n; i++)
      for (int j = i; j <= n; j++)
         if (lmps_map_types_to_unique[i] >= 0 && 
             lmps_map_types_to_unique[j] >= 0) {
            setflag[i][j] = 1;
            count++;
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

   if (domain->dimension != 3)
      error->all(FLERR,"PairKIM only works with 3D problems");

   // set lmps_* bool flags
   set_lmps_flags();
   
   int kimerror;
   // KIM and Model initialization (only once)
   // also sets kim_ind_* and kim_* bool flags
   if (!kim_init_ok)
   {
      kim_init();
      kimerror = pkim->model_init();
      if (kimerror != KIM_STATUS_OK)
         kim_error(__LINE__, "KIM API:model_init() failed", kimerror);
      else
         kim_model_init_ok = true;
   }

   // request none, half, or full neighbor list
   // depending on KIM model requirement

   int irequest = neighbor->request(this);
   if (kim_model_using_cluster)
   {
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->full = 0;
   }
   else
   {
      // make sure comm_reverse expects (at most) 9 values when newton is off
      if (!lmps_using_newton) comm_reverse_off = 9;

      if (kim_model_using_half)
      {
         neighbor->requests[irequest]->half = 1;
         neighbor->requests[irequest]->full = 0;
         // make sure half lists also include local-ghost pairs
         if (lmps_using_newton) neighbor->requests[irequest]->newton = 2;
      }
      else
      {
         neighbor->requests[irequest]->half = 0;
         neighbor->requests[irequest]->full = 1;
         // make sure full lists also include local-ghost pairs
         if (lmps_using_newton) neighbor->requests[irequest]->newton = 0;
      }
   }

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

   return kim_global_cutoff;
}

/* ---------------------------------------------------------------------- */

int PairKIM::pack_reverse_comm(int n, int first, double *buf)
{
   int i,m,last;
   double *fp;
   if (lmps_hybrid) fp = &(lmps_force_tmp[0][0]);
   else fp = &(atom->f[0][0]);

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
      return 3;
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
      return 9;
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
      return 6;
   }
   else
      return 0;
}

/* ---------------------------------------------------------------------- */

void PairKIM::unpack_reverse_comm(int n, int *list, double *buf)
{
   int i,j,m;
   double *fp;
   if (lmps_hybrid) fp = &(lmps_force_tmp[0][0]);
   else fp = &(atom->f[0][0]);

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

void PairKIM::kim_error(int ln, const char* msg, int errcode)
{
   if (errcode == KIM_STATUS_OK) return;
   KIM_API_model::report_error(ln,(char *) __FILE__, (char *) msg,errcode);
   error->all(__FILE__,ln,"Internal KIM error");

   return;
}

/* ---------------------------------------------------------------------- */

int PairKIM::get_neigh(void **kimmdl,int *mode,int *request,
                       int *atom, int *numnei, int **nei1atom, double **pRij)
{
   KIM_API_model *pkim = (KIM_API_model *) *kimmdl;

   int kimerror;
   PairKIM *self = (PairKIM *) pkim->get_test_buffer(&kimerror);

   if (self->kim_model_using_Rij) {
      *pRij = &(self->Rij[0]);
   } else {
      *pRij = 0;
   }
   

   // subvert KIM api by using direct access to self->list
   //
   // get neighObj from KIM API obj
   // NeighList * neiobj = (NeighList * ) 
   // (*pkim).get_data_by_index(self->kim_ind_neighObject, &kimerror);
   NeighList * neiobj = self->list;

   // subvert KIM api by using direct acces to self->lmps_local_tot_num_atoms
   //
   //int * pnAtoms = (int *)
   // (*pkim).get_data_by_index(self->kim_ind_numberOfParticles, &kimerror);
   //int nAtoms = *pnAtoms;
   int nAtoms = self->lmps_local_tot_num_atoms;

   int j, jj, inum, *ilist, *numneigh, **firstneigh;
   inum = neiobj->inum;             //# of I atoms neighbors are stored for
   ilist = neiobj->ilist;           //local indices of I atoms
   numneigh = neiobj->numneigh;     // # of J neighbors for each I atom
   firstneigh = neiobj->firstneigh; // ptr to 1st J int value of each I atom

   if (*mode==0){ //iterator mode
      if (*request==1) { //increment iterator
         if (self->kim_iterator_position < inum) {
            *atom = ilist[self->kim_iterator_position];
            *numnei = numneigh[*atom];

            // strip off neighbor mask for molecular systems
            if (!self->lmps_using_molecular)
               *nei1atom = firstneigh[*atom];
            else
            {
               int n = *numnei;
               int *ptr = firstneigh[*atom];
               int *lmps_stripped_neigh_list = self->lmps_stripped_neigh_list;
               for (int i = 0; i < n; i++)
                  lmps_stripped_neigh_list[i] = *(ptr++) & NEIGHMASK;
               *nei1atom = lmps_stripped_neigh_list;
            }

            // set Rij if needed
            if (self->kim_model_using_Rij) {
               double* x = (double *) 
                 (*pkim).get_data_by_index(self->kim_ind_coordinates, 
                                           &kimerror);
               for (jj=0; jj < *numnei; jj++) {
                  int i = *atom;
                  j = (*nei1atom)[jj];
                  self->Rij[jj*3 +0] = -x[i*3+0] + x[j*3+0];
                  self->Rij[jj*3 +1] = -x[i*3+1] + x[j*3+1];
                  self->Rij[jj*3 +2] = -x[i*3+2] + x[j*3+2];
               }
            }

            // increment iterator
            self->kim_iterator_position++;

            return KIM_STATUS_OK; //successful increment
         } else if (self->kim_iterator_position == inum) {
            *numnei = 0;
            return KIM_STATUS_NEIGH_ITER_PAST_END; //reached end by iterator
         } else if (self->kim_iterator_position > inum || inum < 0){
            self->error->one(FLERR, "KIM neighbor iterator exceeded range");
         }
      } else if (*request == 0){ //restart iterator
         self->kim_iterator_position = 0;
         *numnei = 0;
         return KIM_STATUS_NEIGH_ITER_INIT_OK; //succsesful restart
      }
   } else if (*mode == 1){//locator mode
      //...
      if (*request < inum) {
         *atom = *request;
         *numnei = numneigh[*atom];

         // strip off neighbor mask for molecular systems
         if (!self->lmps_using_molecular)
            *nei1atom = firstneigh[*atom];
         else
         {
            int n = *numnei;
            int *ptr = firstneigh[*atom];
            int *lmps_stripped_neigh_list = self->lmps_stripped_neigh_list;
            for (int i = 0; i < n; i++)
               lmps_stripped_neigh_list[i] = *(ptr++) & NEIGHMASK;
            *nei1atom = lmps_stripped_neigh_list;
         }

         // set Rij if needed
         if (self->kim_model_using_Rij){
            double* x = (double *) 
              (*pkim).get_data_by_index(self->kim_ind_coordinates, &kimerror);
            for(int jj=0; jj < *numnei; jj++){
               int i = *atom;
               int j = (*nei1atom)[jj];
               self->Rij[jj*3 +0] = -x[i*3+0] + x[j*3+0];
               self->Rij[jj*3 +1] = -x[i*3+1] + x[j*3+1];
               self->Rij[jj*3 +2] = -x[i*3+2] + x[j*3+2];
            }
         }
         return KIM_STATUS_OK; //successful end
      }
      else if (*request >= nAtoms || inum < 0)
         return KIM_STATUS_NEIGH_INVALID_REQUEST;
      else if (*request >= inum) {
         *atom = *request;
         *numnei = 0;
         return KIM_STATUS_OK; //successfull but no neighbors in the list
      }
   } else return KIM_STATUS_NEIGH_INVALID_MODE; //invalid mode

   return -16; //should not get here: unspecified error
}

/* ---------------------------------------------------------------------- */

void PairKIM::kim_free()
{
   int kimerror;

   if (kim_model_init_ok)
   {
      kimerror = pkim->model_destroy();
      kim_model_init_ok = false;
   }
   if (kim_init_ok)
   {
      pkim->free(&kimerror);
      kim_init_ok = false;
   }
   if (pkim != 0)
   {
      delete pkim;
      pkim = 0;
   }
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

   // determine KIM Model capabilities (used in this function below)
   set_kim_model_has_flags();

   // create appropriate KIM descriptor file
   char* test_descriptor_string = 0;
   // allocate memory for test_descriptor_string and write descriptor file
   write_descriptor(&test_descriptor_string);

   // initialize KIM model
   pkim = new KIM_API_model();
   kimerror = pkim->string_init(test_descriptor_string, kim_modelname);
   if (kimerror != KIM_STATUS_OK)
      kim_error(__LINE__,"KIM initialization failed", kimerror);
   else
   {
      kim_init_ok = true;
      delete [] test_descriptor_string;
      test_descriptor_string = 0;
   }

   // determine kim_model_using_* true/false values
   //
   // check for half or full list
   kim_model_using_half = (pkim->is_half_neighbors(&kimerror));
   //
   char * NBC_method =(char *) pkim->get_NBC_method(&kimerror);
   kim_error(__LINE__,"NBC method not set",kimerror);
   // check for CLUSTER mode
   kim_model_using_cluster = (strcmp(NBC_method,"CLUSTER")==0);
   // check if Rij needed for get_neigh
   kim_model_using_Rij = ((strcmp(NBC_method,"NEIGH_RVEC_H")==0) ||
                          (strcmp(NBC_method,"NEIGH_RVEC_F")==0));
   free((void*)NBC_method);

   // get correct index of each variable in kim_api object
   pkim->getm_index(&kimerror, 3*13,
    "coordinates", &kim_ind_coordinates, 1,
    "cutoff", &kim_ind_cutoff, 1,
    "numberOfParticles", &kim_ind_numberOfParticles, 1,
    "numberParticleTypes", &kim_ind_numberParticleTypes, 1,
    "particleTypes", &kim_ind_particleTypes, 1,
    "numberContributingParticles", &kim_ind_numberContributingParticles,
                                   kim_model_using_half,
    "particleEnergy", &kim_ind_particleEnergy,
                      (int) kim_model_has_particleEnergy,
    "energy", &kim_ind_energy, (int) kim_model_has_energy,
    "forces", &kim_ind_forces, (int) kim_model_has_forces,
    "neighObject", &kim_ind_neighObject, (int) !kim_model_using_cluster,
    "get_neigh", &kim_ind_get_neigh, (int) !kim_model_using_cluster,
    "particleVirial", &kim_ind_particleVirial,
                      (int) kim_model_has_particleVirial,
    "virial", &kim_ind_virial, no_virial_fdotr_compute);
   kim_error(__LINE__,"getm_index",kimerror);

   // setup mapping between LAMMPS unique elements and KIM particle type codes
   kim_particle_codes = new int[lmps_num_unique_elements];
   kim_particle_codes_ok = true;
   for(int i = 0; i < lmps_num_unique_elements; i++){
      int kimerror;
      kim_particle_codes[i]
         = pkim->get_partcl_type_code(lmps_unique_elements[i], &kimerror);
      kim_error(__LINE__, "create_kim_particle_codes: symbol not found ",
                kimerror);
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

   int kimerror;
   pkim->setm_data_by_index(&kimerror, 4*6,
    kim_ind_numberParticleTypes, 1, (void *) &(atom->ntypes), 1,
    kim_ind_cutoff, 1, (void *) &(kim_global_cutoff), 1,
    kim_ind_numberOfParticles, 1, (void *) &lmps_local_tot_num_atoms,  1,
    kim_ind_numberContributingParticles, 1, (void *) &(atom->nlocal),
                                         (int) kim_model_using_half,
    kim_ind_energy, 1, (void *) &(eng_vdwl), (int) kim_model_has_energy,
    kim_ind_virial, 1, (void *) &(virial[0]), no_virial_fdotr_compute);
   kim_error(__LINE__, "setm_data_by_index", kimerror);
   if (!kim_model_using_cluster)
   {
      kimerror = pkim->set_method_by_index(kim_ind_get_neigh, 1,
                                           (func_ptr) &get_neigh);
      kim_error(__LINE__, "set_method_by_index", kimerror);
   }

   pkim->set_test_buffer((void *)this, &kimerror);
   kim_error(__LINE__, "set_test_buffer", kimerror);

   return;
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_volatiles()
{
   int kimerror;
   lmps_local_tot_num_atoms = (int) (atom->nghost + atom->nlocal);
   intptr_t nall = (intptr_t) lmps_local_tot_num_atoms;

   pkim->setm_data_by_index(&kimerror, 4*2,
    kim_ind_coordinates, 3*nall, (void*) &(atom->x[0][0]), 1,
    kim_ind_particleTypes, nall, (void*) kim_particleTypes, 1);
   kim_error(__LINE__, "setm_data_by_index", kimerror);

   if (kim_model_has_particleEnergy && (eflag_atom == 1))
   {
      kimerror = pkim->set_data_by_index(kim_ind_particleEnergy, nall,
                                         (void*) eatom);
      kim_error(__LINE__, "set_data_by_index", kimerror);
   }

   if (kim_model_has_particleVirial && (vflag_atom == 1))
   {
      kimerror = pkim->set_data_by_index(kim_ind_particleVirial, 6*nall,
                                         (void*) &(vatom[0][0]));
      kim_error(__LINE__, "set_data_by_index", kimerror);
   }

   if (kim_model_has_forces)
   {
      if (lmps_hybrid)
         kimerror = pkim->set_data_by_index(kim_ind_forces, nall*3,
                                            (void*) &(lmps_force_tmp[0][0]));
      else
         kimerror = pkim->set_data_by_index(kim_ind_forces, nall*3,
                                            (void*) &(atom->f[0][0]));
      kim_error(__LINE__, "setm_data_by_index", kimerror);
   }

   // subvert the KIM api by direct access to this->list in get_neigh
   //
   //if (!kim_model_using_cluster)
   //   kimerror = pkim->set_data_by_index(kim_ind_neighObject, 1,
   //                                      (void*) this->list);

   if (kim_model_has_particleVirial)
   {
      if(vflag_atom != 1) {
         pkim->set_compute_by_index(kim_ind_particleVirial, KIM_COMPUTE_FALSE,
                                    &kimerror);
      } else {
         pkim->set_compute_by_index(kim_ind_particleVirial, KIM_COMPUTE_TRUE,
                                    &kimerror);
      }
   }

   if (no_virial_fdotr_compute == 1)
   {
      pkim->set_compute_by_index(kim_ind_virial,
       ((vflag_global != 1) ? KIM_COMPUTE_FALSE : KIM_COMPUTE_TRUE),
       &kimerror);
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
   lmps_hybrid = (force->pair_match("hybrid",0));

   // support cluster mode if everything is just right
   lmps_support_cluster = ((domain->xperiodic == 0 &&
                            domain->yperiodic == 0 &&
                            domain->zperiodic == 0
                           )
                           &&
                           (comm->nprocs == 1)
                          );

   // determine unit system and set lmps_units flag
   if ((strcmp(update->unit_style,"real")==0))
      lmps_units = REAL;
   else if ((strcmp(update->unit_style,"metal")==0))
      lmps_units = METAL;
   else if ((strcmp(update->unit_style,"si")==0))
      lmps_units = SI;
   else if ((strcmp(update->unit_style,"cgs")==0))
      lmps_units = CGS;
   else if ((strcmp(update->unit_style,"electron")==0))
      lmps_units = ELECTRON;
   else if ((strcmp(update->unit_style,"lj")==0))
      error->all(FLERR,"LAMMPS unit_style lj not supported by KIM models");
   else
      error->all(FLERR,"Unknown unit_style");

   return;
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_kim_model_has_flags()
{
   KIM_API_model mdl;

   int kimerror;

   // get KIM API object representing the KIM Model only
   kimerror = mdl.model_info(kim_modelname);
   kim_error(__LINE__,"KIM initialization failed", kimerror);

   // determine if the KIM Model can compute the total energy
   mdl.get_index((char*) "energy", &kimerror);
   kim_model_has_energy = (kimerror == KIM_STATUS_OK);
   if (!kim_model_has_energy) 
     error->warning(FLERR,"KIM Model does not provide `energy'; "
                    "Potential energy will be zero");

   // determine if the KIM Model can compute the forces
   mdl.get_index((char*) "forces", &kimerror);
   kim_model_has_forces = (kimerror == KIM_STATUS_OK);
   if (!kim_model_has_forces) 
     error->warning(FLERR,"KIM Model does not provide `forces'; "
                    "Forces will be zero");

   // determine if the KIM Model can compute the particleEnergy
   mdl.get_index((char*) "particleEnergy", &kimerror);
   kim_model_has_particleEnergy = (kimerror == KIM_STATUS_OK);
   if (!kim_model_has_particleEnergy) 
     error->warning(FLERR,"KIM Model does not provide `particleEnergy'; "
                    "energy per atom will be zero");

   // determine if the KIM Model can compute the particleVerial
   mdl.get_index((char*) "particleVirial", &kimerror);
   kim_model_has_particleVirial = (kimerror == KIM_STATUS_OK);
   mdl.get_index((char*) "process_dEdr", &kimerror);
   kim_model_has_particleVirial = kim_model_has_particleVirial || 
     (kimerror == KIM_STATUS_OK);
   if (!kim_model_has_particleVirial) 
     error->warning(FLERR,"KIM Model does not provide `particleVirial'; "
                    "virial per atom will be zero");

   // tear down KIM API object
   mdl.free(&kimerror);
   // now destructor will do the remaining tear down for mdl

   return;
}

/* ---------------------------------------------------------------------- */

void PairKIM::write_descriptor(char** test_descriptor_string)
{
   // allocate memory
   if (*test_descriptor_string != 0) 
     error->all(FLERR, "Test_descriptor_string already allocated");
   // assuming 75 lines at 100 characters each (should be plenty)
   *test_descriptor_string = new char[100*75]; 
   // initialize
   strcpy(*test_descriptor_string, "");

   // Write Test name and units
   strcat(*test_descriptor_string,
      "# This file is automatically generated from LAMMPS pair_style "
      "PairKIM command\n"
      "\n"
      "# Base units\n");
   switch (lmps_units)
   {
      case REAL:
         strcat(*test_descriptor_string,
      "Unit_length      := A\n"
      "Unit_energy      := kcal/mol\n"
      "Unit_charge      := e\n"
      "Unit_temperature := K\n"
      "Unit_time        := fs\n\n");
      break;
      case METAL:
         strcat(*test_descriptor_string,
      "Unit_length      := A\n"
      "Unit_energy      := eV\n"
      "Unit_charge      := e\n"
      "Unit_temperature := K\n"
      "Unit_time        := ps\n\n");
      break;
      case SI:
         strcat(*test_descriptor_string,
      "Unit_length      := m\n"
      "Unit_energy      := J\n"
      "Unit_charge      := C\n"
      "Unit_temperature := K\n"
      "Unit_time        := s\n\n");
      break;
      case CGS:
         strcat(*test_descriptor_string,
      "Unit_length      := cm\n"
      "Unit_energy      := erg\n"
      "Unit_charge      := statC\n"
      "Unit_temperature := K\n"
      "Unit_time        := s\n\n");
      break;
      case ELECTRON:
         strcat(*test_descriptor_string,
      "Unit_length      := Bohr\n"
      "Unit_energy      := Hartree\n"
      "Unit_charge      := e\n"
      "Unit_temperature := K\n"
      "Unit_time        := fs\n\n");
      break;
   }

   // Write Supported types section
   strcat(*test_descriptor_string,
      "\n"
      "SUPPORTED_ATOM/PARTICLES_TYPES:\n"
      "# Symbol/name           Type            code\n\n");
   int code=1;
   char* tmp_line = 0;
   tmp_line = new char[100];
   for (int i=0; i < lmps_num_unique_elements; i++){
      sprintf(tmp_line, "%-24s%-16s%-3i\n", lmps_unique_elements[i], 
              "spec", code++);
      strcat(*test_descriptor_string, tmp_line);
   }
   delete [] tmp_line;
   tmp_line = 0;
   strcat(*test_descriptor_string, "\n");

   // Write conventions section
   strcat(*test_descriptor_string,
      "\n"
      "CONVENTIONS:\n"
      "# Name                  Type\n\n"
      "ZeroBasedLists          flag\n");
   // can use iterator or locator neighbor mode, unless in hybrid mode
   if (lmps_hybrid)
      strcat(*test_descriptor_string,
      "Neigh_IterAccess        flag\n");
   else
      strcat(*test_descriptor_string,
      "Neigh_BothAccess        flag\n\n");

   strcat(*test_descriptor_string,
      "NEIGH_PURE_H            flag\n"
      "NEIGH_PURE_F            flag\n"
      "NEIGH_RVEC_H            flag\n"
      "NEIGH_RVEC_F            flag\n");
   // @@ add code for MI_OPBC_? support ????
   if (lmps_support_cluster)
   {
      strcat(*test_descriptor_string,
      "CLUSTER                 flag\n\n");
   }
   else
   {
      strcat(*test_descriptor_string, "\n");
   }

   // Write input section
   strcat(*test_descriptor_string,
      "\n"
      "MODEL_INPUT:\n"
      "# Name                         Type         Unit    Shape\n\n"
      "numberOfParticles              integer      none    []\n\n"
      "numberContributingParticles    integer      none    []\n\n"
      "numberParticleTypes            integer      none    []\n\n"
      "particleTypes                  integer      none    "
       "[numberOfParticles]\n\n"
      "coordinates                    double       length  "
       "[numberOfParticles,3]\n\n"
      "neighObject                    pointer      none    []\n\n"
      "get_neigh                      method       none    []\n\n");

   // Write output section
   strcat(*test_descriptor_string,
      "\n"
      "MODEL_OUPUT:\n"
      "# Name                         Type         Unit    Shape\n\n"
      "compute                        method       none    []\n\n"
      "destroy                        method       none    []\n\n"
      "cutoff                         double       length  []\n\n");
   if (kim_model_has_energy) strcat(*test_descriptor_string,
      "energy                         double       energy  []\n\n");
   if (kim_model_has_forces) strcat(*test_descriptor_string,
      "forces                         double       force   "
       "[numberOfParticles,3]\n\n");
   if (kim_model_has_particleEnergy) strcat(*test_descriptor_string,
      "particleEnergy                 double       energy  "
       "[numberOfParticles]\n\n");
   if (no_virial_fdotr_compute == 1) strcat(*test_descriptor_string,
      "virial                         double       energy  [6] \n\n");
   if (kim_model_has_particleVirial) strcat(*test_descriptor_string,
      "particleVirial                 double       energy  "
       "[numberOfParticles,6] \n\n");

   return;
}
