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

#include <mpi.h>
#include <cstring>
#include "server_md.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "integrate.h"
#include "kspace.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

// CSlib interface

#include "cslib.h"

using namespace LAMMPS_NS;
using namespace CSLIB_NS;

enum{OTHER,REAL,METAL};
enum{SETUP=1,STEP};
enum{DIM=1,PERIODICITY,ORIGIN,BOX,NATOMS,NTYPES,TYPES,COORDS,UNITS,CHARGE};
enum{FORCES=1,ENERGY,PRESSURE,ERROR};

/* ---------------------------------------------------------------------- */

ServerMD::ServerMD(LAMMPS *lmp) : Pointers(lmp)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Server command before simulation box is defined");

  if (!atom->map_style) error->all(FLERR,"Server md requires atom map");
  if (atom->tag_enable == 0) error->all(FLERR,"Server md requires atom IDs");
  if (sizeof(tagint) != 4) error->all(FLERR,"Server md requires 4-byte atom IDs");

  if (strcmp(update->unit_style,"real") == 0) units = REAL;
  else if (strcmp(update->unit_style,"metal") == 0) units = METAL;
  else units = OTHER;

  // unit conversion factors for REAL
  // otherwise not needed
  // local computation in REAL units, send message in METAL units

  fconvert = econvert = pconvert = 1.0;
  if (units == REAL) {
    fconvert = econvert = 1.0 / 23.06035;    // Kcal/mole -> eV
    pconvert = 1.0 / 0.986923;               // atmospheres -> bars
  }

  fcopy = NULL;
}

/* ---------------------------------------------------------------------- */

ServerMD::~ServerMD()
{
  memory->destroy(fcopy);
}

/* ---------------------------------------------------------------------- */

void ServerMD::loop()
{
  int i,j,m;

  // cs = instance of CSlib

  CSlib *cs = (CSlib *) lmp->cslib;

  // counters

  int forcecalls = 0;
  int neighcalls = 0;

  // loop on messages
  // receive a message, process it, send return message

  int msgID,nfield;
  int *fieldID,*fieldtype,*fieldlen;

  while (1) {
    msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);
    if (msgID < 0) break;

    // SETUP receive at beginning of each run
    // required fields: DIM, PERIODICTY, ORIGIN, BOX,
    //                  NATOMS, NTYPES, TYPES, COORDS
    // optional fields: others in enum above

    if (msgID == SETUP) {

      int dim = 0;
      int *periodicity = NULL;
      int natoms = -1;
      int ntypes = -1;
      double *origin = NULL;
      double *box = NULL;
      int *types = NULL;
      double *coords = NULL;
      char *unit_style = NULL;
      double *charge = NULL;

      for (int ifield = 0; ifield < nfield; ifield++) {
        if (fieldID[ifield] == DIM) {
          dim = cs->unpack_int(DIM);
          if (dim != domain->dimension)
            error->all(FLERR,"Server md dimension mis-match with client");
        } else if (fieldID[ifield] == PERIODICITY) {
          periodicity = (int *) cs->unpack(PERIODICITY);
          if (periodicity[0] != domain->periodicity[0] ||
              periodicity[1] != domain->periodicity[1] ||
              periodicity[2] != domain->periodicity[2])
            error->all(FLERR,"Server md periodicity mis-match with client");
        } else if (fieldID[ifield] == ORIGIN) {
          origin = (double *) cs->unpack(ORIGIN);
        } else if (fieldID[ifield] == BOX) {
          box = (double *) cs->unpack(BOX);
        } else if (fieldID[ifield] == NATOMS) {
          natoms = cs->unpack_int(NATOMS);
          if (3*natoms > INT_MAX)
            error->all(FLERR,"Server md max atoms is 1/3 of 2^31");
        } else if (fieldID[ifield] == NTYPES) {
          ntypes = cs->unpack_int(NTYPES);
          if (ntypes != atom->ntypes)
            error->all(FLERR,"Server md ntypes mis-match with client");
        } else if (fieldID[ifield] == TYPES) {
          types = (int *) cs->unpack(TYPES);
        } else if (fieldID[ifield] == COORDS) {
          coords = (double *) cs->unpack(COORDS);

        } else if (fieldID[ifield] == UNITS) {
          unit_style = (char *) cs->unpack(UNITS);
        } else if (fieldID[ifield] == CHARGE) {
          charge = (double *) cs->unpack(CHARGE);
        } else error->all(FLERR,"Server md setup field unknown");
      }

      if (dim == 0 || !periodicity || !origin || !box ||
          natoms < 0 || ntypes < 0 || !types || !coords)
        error->all(FLERR,"Required server md setup field not received");

      if (unit_style && strcmp(unit_style,update->unit_style) != 0)
        error->all(FLERR,"Server md does not match client units");

      if (charge && atom->q_flag == 0)
        error->all(FLERR,"Server md does not define atom attribute q");

      // reset box, global and local
      // reset proc decomposition

      if ((box[3] != 0.0 || box[6] != 0.0 || box[7] != 0.0) &&
          domain->triclinic == 0)
        error->all(FLERR,"Server md is not initialized for a triclinic box");

      box_change(origin,box);

      domain->set_initial_box();
      domain->set_global_box();
      comm->set_proc_grid();
      domain->set_local_box();

      // clear all atoms

      atom->nlocal = 0;
      atom->nghost = 0;

      // add atoms that are in my sub-box

      int nlocal = 0;
      for (int i = 0; i < natoms; i++) {
        if (!domain->ownatom(i+1,&coords[3*i],NULL,0)) continue;
        atom->avec->create_atom(types[i],&coords[3*i]);
        atom->tag[nlocal] = i+1;
        if (charge) atom->q[nlocal] = charge[i];
        nlocal++;
      }

      int ntotal;
      MPI_Allreduce(&atom->nlocal,&ntotal,1,MPI_INT,MPI_SUM,world);
      if (ntotal != natoms)
        error->all(FLERR,"Server md atom count does not match client");

      atom->map_init();
      atom->map_set();
      atom->natoms = natoms;

      // allocate fcopy if needed

      if (units == REAL) {
        memory->destroy(fcopy);
        memory->create(fcopy,atom->nlocal,3,"server/md:fcopy");
      }

      // perform system setup() for dynamics
      // also OK for minimization, since client runs minimizer
      // return forces, energy, virial to client

      update->whichflag = 1;
      lmp->init();
      update->integrate->setup_minimal(1);
      neighcalls++;
      forcecalls++;

      send_fev(msgID);

    // STEP receive at each timestep of run or minimization
    // required fields: COORDS
    // optional fields: ORIGIN, BOX

    } else if (msgID == STEP) {

      double *coords = NULL;
      double *origin = NULL;
      double *box = NULL;

      for (int ifield = 0; ifield < nfield; ifield++) {
        if (fieldID[ifield] == COORDS) {
          coords = (double *) cs->unpack(COORDS);
        } else if (fieldID[ifield] == ORIGIN) {
          origin = (double *) cs->unpack(ORIGIN);
        } else if (fieldID[ifield] == BOX) {
          box = (double *) cs->unpack(BOX);
        } else error->all(FLERR,"Server md step field unknown");
      }

      if (!coords)
        error->all(FLERR,"Required server md step field not received");

      // change box size/shape, only if origin and box received
      // reset global/local box like FixDeform at end_of_step()

      if (origin && box) {
        if ((box[3] != 0.0 || box[6] != 0.0 || box[7] != 0.0) &&
            domain->triclinic == 0)
          error->all(FLERR,"Server md is not initialized for a triclinic box");
        box_change(origin,box);
        domain->set_global_box();
        domain->set_local_box();
      }

      // assign received coords to owned atoms
      // closest_image() insures xyz matches current server PBC

      double **x = atom->x;
      int nlocal = atom->nlocal;
      int nall = atom->natoms;

      j = 0;
      for (tagint id = 1; id <= nall; id++) {
        m = atom->map(id);
        if (m < 0 || m >= nlocal) j += 3;
        else {
          domain->closest_image(x[m],&coords[j],x[m]);
          j += 3;
        }
      }

      // if no need to reneighbor:
      //   ghost comm
      //   setup_minimal(0) which just computes forces
      // else:
      //   setup_minimal(1) for pbc, reset_box, reneigh, forces

      int nflag = neighbor->decide();
      if (nflag == 0) {
        comm->forward_comm();
        update->integrate->setup_minimal(0);
      } else {
        update->integrate->setup_minimal(1);
        neighcalls++;
      }

      forcecalls++;

      send_fev(msgID);

    } else error->all(FLERR,"Server MD received unrecognized message");
  }

  // final reply to client

  cs->send(0,0);

  // stats

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Server MD calls = %d\n",forcecalls);
      fprintf(screen,"Server MD reneighborings = %d\n",neighcalls);
    }
    if (logfile) {
      fprintf(logfile,"Server MD calls = %d\n",forcecalls);
      fprintf(logfile,"Server MD reneighborings %d\n",neighcalls);
    }
  }

  // clean up

  delete cs;
  lmp->cslib = NULL;
}

/* ----------------------------------------------------------------------
   box change due to received message
------------------------------------------------------------------------- */

void ServerMD::box_change(double *origin, double *box)
{
  domain->boxlo[0] = origin[0];
  domain->boxlo[1] = origin[1];
  domain->boxlo[2] = origin[2];

  domain->boxhi[0] = origin[0] + box[0];
  domain->boxhi[1] = origin[1] + box[4];
  domain->boxhi[2] = origin[2] + box[8];

  domain->xy = box[3];
  domain->xz = box[6];
  domain->yz = box[7];
}

/* ----------------------------------------------------------------------
   return message with forces, energy, pressure tensor
   pressure tensor should be just pair and KSpace contributions
   required fields: FORCES, ENERGY, PRESSURE
   optional field: ERROR (not ever sending)
------------------------------------------------------------------------- */

void ServerMD::send_fev(int msgID)
{
  CSlib *cs = (CSlib *) lmp->cslib;

  cs->send(msgID,3);

  double *forces = NULL;
  if (atom->nlocal) {
    if (units != REAL) forces = &atom->f[0][0];
    else {
      double **f = atom->f;
      int nlocal = atom->nlocal;
      for (int i = 0; i < nlocal; i++) {
        fcopy[i][0] = fconvert*f[i][0];
        fcopy[i][1] = fconvert*f[i][1];
        fcopy[i][2] = fconvert*f[i][2];
      }
      forces = &fcopy[0][0];
    }
  }
  cs->pack_parallel(FORCES,4,atom->nlocal,atom->tag,3,forces);

  double eng = force->pair->eng_vdwl + force->pair->eng_coul;
  double engall;
  MPI_Allreduce(&eng,&engall,1,MPI_DOUBLE,MPI_SUM,world);
  engall *= econvert;
  cs->pack_double(ENERGY,engall);

  double v[6],vall[6];
  for (int i = 0; i < 6; i++)
    v[i] = force->pair->virial[i];
  MPI_Allreduce(&v,&vall,6,MPI_DOUBLE,MPI_SUM,world);

  if (force->kspace)
    for (int i = 0; i < 6; i++)
      vall[i] += force->kspace->virial[i];

  double nktv2p = force->nktv2p;
  double volume = domain->xprd * domain->yprd * domain->zprd;
  double factor = pconvert / volume * nktv2p;
  for (int i = 0; i < 6; i++) vall[i] *= factor;

  cs->pack(PRESSURE,4,6,vall);
}
