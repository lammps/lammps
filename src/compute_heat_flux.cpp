// clang-format off
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
   Contributing authors: German Samolyuk (ORNL) and
                         Mario Pinto (Computational Research Lab, Pune, India)
------------------------------------------------------------------------- */

#include "compute_heat_flux.h"

#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

ComputeHeatFlux::ComputeHeatFlux(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  id_ke(nullptr), id_pe(nullptr), id_stress(nullptr)
{
  if (narg != 6) error->all(FLERR,"Illegal compute heat/flux command");

  vector_flag = 1;
  size_vector = 6;
  extvector = 1;

  // store ke/atom, pe/atom, stress/atom IDs used by heat flux computation
  // insure they are valid for these computations

  id_ke = utils::strdup(arg[3]);
  id_pe = utils::strdup(arg[4]);
  id_stress = utils::strdup(arg[5]);

  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  int istress = modify->find_compute(id_stress);
  if (ike < 0 || ipe < 0 || istress < 0)
    error->all(FLERR,"Could not find compute heat/flux compute ID");
  if (strcmp(modify->compute[ike]->style,"ke/atom") != 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute ke/atom");
  if (modify->compute[ipe]->peatomflag == 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute pe/atom");
  if (modify->compute[istress]->pressatomflag != 1
      && modify->compute[istress]->pressatomflag != 2)
    error->all(FLERR,
               "Compute heat/flux compute ID does not compute stress/atom or centroid/stress/atom");

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeHeatFlux::~ComputeHeatFlux()
{
  delete [] id_ke;
  delete [] id_pe;
  delete [] id_stress;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFlux::init()
{
  // error checks

  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  int istress = modify->find_compute(id_stress);
  if (ike < 0 || ipe < 0 || istress < 0)
    error->all(FLERR,"Could not find compute heat/flux compute ID");

  c_ke = modify->compute[ike];
  c_pe = modify->compute[ipe];
  c_stress = modify->compute[istress];
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFlux::compute_vector()
{
  invoked_vector = update->ntimestep;

  // invoke 3 computes if they haven't been already

  if (!(c_ke->invoked_flag & Compute::INVOKED_PERATOM)) {
    c_ke->compute_peratom();
    c_ke->invoked_flag |= Compute::INVOKED_PERATOM;
  }
  if (!(c_pe->invoked_flag & Compute::INVOKED_PERATOM)) {
    c_pe->compute_peratom();
    c_pe->invoked_flag |= Compute::INVOKED_PERATOM;
  }
  if (!(c_stress->invoked_flag & Compute::INVOKED_PERATOM)) {
    c_stress->compute_peratom();
    c_stress->invoked_flag |= Compute::INVOKED_PERATOM;
  }

  // heat flux vector = jc[3] + jv[3]
  // jc[3] = convective portion of heat flux = sum_i (ke_i + pe_i) v_i[3]
  // jv[3] = virial portion of heat flux = sum_i (stress_tensor_i . v_i[3])
  // normalization by volume is not included

  double *ke = c_ke->vector_atom;
  double *pe = c_pe->vector_atom;
  double **stress = c_stress->array_atom;

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double jc[3] = {0.0,0.0,0.0};
  double jv[3] = {0.0,0.0,0.0};
  double eng;

  // heat flux via centroid atomic stress
  if (c_stress->pressatomflag == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        eng = pe[i] + ke[i];
        jc[0] += eng*v[i][0];
        jc[1] += eng*v[i][1];
        jc[2] += eng*v[i][2];
        // stress[0]: rijx*fijx
        // stress[1]: rijy*fijy
        // stress[2]: rijz*fijz
        // stress[3]: rijx*fijy
        // stress[4]: rijx*fijz
        // stress[5]: rijy*fijz
        // stress[6]: rijy*fijx
        // stress[7]: rijz*fijx
        // stress[8]: rijz*fijy
        // jv[0]  = rijx fijx vjx + rijx fijy vjy + rijx fijz vjz
        jv[0] -= stress[i][0]*v[i][0] + stress[i][3]*v[i][1] +
          stress[i][4]*v[i][2];
        // jv[1]  = rijy fijx vjx + rijy fijy vjy + rijy fijz vjz
        jv[1] -= stress[i][6]*v[i][0] + stress[i][1]*v[i][1] +
          stress[i][5]*v[i][2];
        // jv[2]  = rijz fijx vjx + rijz fijy vjy + rijz fijz vjz
        jv[2] -= stress[i][7]*v[i][0] + stress[i][8]*v[i][1] +
          stress[i][2]*v[i][2];
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        eng = pe[i] + ke[i];
        jc[0] += eng*v[i][0];
        jc[1] += eng*v[i][1];
        jc[2] += eng*v[i][2];
        jv[0] -= stress[i][0]*v[i][0] + stress[i][3]*v[i][1] +
          stress[i][4]*v[i][2];
        jv[1] -= stress[i][3]*v[i][0] + stress[i][1]*v[i][1] +
          stress[i][5]*v[i][2];
        jv[2] -= stress[i][4]*v[i][0] + stress[i][5]*v[i][1] +
          stress[i][2]*v[i][2];
      }
    }
  }

  // convert jv from stress*volume to energy units via nktv2p factor

  double nktv2p = force->nktv2p;
  jv[0] /= nktv2p;
  jv[1] /= nktv2p;
  jv[2] /= nktv2p;

  // sum across all procs
  // 1st 3 terms are total heat flux
  // 2nd 3 terms are just conductive portion

  double data[6] = {jc[0]+jv[0],jc[1]+jv[1],jc[2]+jv[2],jc[0],jc[1],jc[2]};
  MPI_Allreduce(data,vector,6,MPI_DOUBLE,MPI_SUM,world);
}
