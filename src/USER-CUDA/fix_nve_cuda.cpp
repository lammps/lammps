/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

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

#include <cstdio>
#include <cstring>
#include "fix_nve_cuda.h"
#include "fix_nve_cuda_cu.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "user_cuda.h"
#include "cuda_modify_flags.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

/* ---------------------------------------------------------------------- */

FixNVECuda::FixNVECuda(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  cuda = lmp->cuda;

  if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  if (strcmp(style,"nve/sphere") != 0 && narg < 3)
                error->all(FLERR,"Illegal fix nve command");

        time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVECuda::setmask()
{
        int mask = 0;
        mask |= INITIAL_INTEGRATE_CUDA;
        mask |= FINAL_INTEGRATE_CUDA;
        // mask |= INITIAL_INTEGRATE_RESPA_CUDA;
        // mask |= FINAL_INTEGRATE_RESPA_CUDA;
        return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVECuda::init()
{
        dtv = update->dt;
        dtf = 0.5 * update->dt * force->ftm2v;

        if (strstr(update->integrate_style,"respa"))
                step_respa = ((Respa *) update->integrate)->step;

        triggerneighsq= cuda->shared_data.atom.triggerneighsq;
    cuda->neighbor_decide_by_integrator=1;
    Cuda_FixNVECuda_Init(&cuda->shared_data,dtv,dtf);

}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVECuda::initial_integrate(int vflag)
{
        if(triggerneighsq!=cuda->shared_data.atom.triggerneighsq)
        {
                triggerneighsq= cuda->shared_data.atom.triggerneighsq;
                Cuda_FixNVECuda_Init(&cuda->shared_data,dtv,dtf);
        }
        int nlocal = atom->nlocal;
        if(igroup == atom->firstgroup) nlocal = atom->nfirst;

    Cuda_FixNVECuda_InitialIntegrate(& cuda->shared_data, groupbit,nlocal);
}

/* ---------------------------------------------------------------------- */

void FixNVECuda::final_integrate()
{
        int nlocal = atom->nlocal;
        if(igroup == atom->firstgroup) nlocal = atom->nfirst;

        Cuda_FixNVECuda_FinalIntegrate(& cuda->shared_data, groupbit,nlocal);
}

/* ---------------------------------------------------------------------- */

void FixNVECuda::initial_integrate_respa(int vflag, int ilevel, int flag)
{
        //this point should not be reached yet since RESPA is not supported
        if (flag) return;             // only used by NPT,NPH

        dtv = step_respa[ilevel];
        dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

        // innermost level - NVE update of v and x
        // all other levels - NVE update of v

        if(ilevel == 0) initial_integrate(vflag);
        else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVECuda::final_integrate_respa(int ilevel, int iloop)
{
        //this point should not be reached yet since RESPA is not supported
        dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
        final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVECuda::reset_dt()
{
        dtv = update->dt;
        dtf = 0.5 * update->dt * force->ftm2v;
        Cuda_FixNVECuda_Init(&cuda->shared_data,dtv,dtf);
}
