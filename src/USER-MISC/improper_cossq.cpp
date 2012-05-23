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
   Contributing author: Georgios G. Vogiatzis (CoMSE, NTU Athens),
     gvog@chemeng.ntua.gr
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "improper_cossq.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

ImproperCossq::ImproperCossq(LAMMPS *lmp) : Improper(lmp) {}

/* ---------------------------------------------------------------------- */

ImproperCossq::~ImproperCossq()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(chi);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperCossq::compute(int eflag, int vflag)
{
   int i1,i2,i3,i4,m,n,type;
   double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z ;
   double eimproper,f1[3],f2[3],f3[3],f4[3];
   double rjisq, rji, rlksq, rlk, cosphi, angfac;
   double cjiji, clkji, clklk, cfact1, cfact2, cfact3;


   eimproper = 0.0;
   if (eflag || vflag) ev_setup(eflag,vflag);
   else evflag = 0;

   double **x = atom->x;
   double **f = atom->f;
   int **improperlist = neighbor->improperlist;
   int nimproperlist = neighbor->nimproperlist;
   int nlocal = atom->nlocal;
   int newton_bond = force->newton_bond;

   for (n = 0; n < nimproperlist; n++) {
      /* Ask the improper list for the atom types. */
      i1 = improperlist[n][0]; 
      i2 = improperlist[n][1];
      i3 = improperlist[n][2];
      i4 = improperlist[n][3];
      type = improperlist[n][4];

      /* separation vector between i1 and i2, (i2-i1) */
      vb1x = x[i2][0] - x[i1][0];
      vb1y = x[i2][1] - x[i1][1];
      vb1z = x[i2][2] - x[i1][2];
      domain->minimum_image(vb1x,vb1y,vb1z);
      rjisq = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z ;
      rji = sqrt(rjisq);

      /* separation vector between i2 and i3 (i3-i2) */
      vb2x = x[i3][0] - x[i2][0];
      vb2y = x[i3][1] - x[i2][1];
      vb2z = x[i3][2] - x[i2][2];
      domain->minimum_image(vb2x,vb2y,vb2z);

      /* separation vector between i3 and i4, (i4-i3) */
      vb3x = x[i4][0] - x[i3][0];
      vb3y = x[i4][1] - x[i3][1];
      vb3z = x[i4][2] - x[i3][2];
      domain->minimum_image(vb3x,vb3y,vb3z);
      rlksq = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z ;
      rlk = sqrt(rlksq);

      cosphi = (vb3x*vb1x + vb3y*vb1y + vb3z*vb1z)/(rji * rlk);
     
      /* Check that cos(phi) is in the correct limits. */     
      if (cosphi > 1.0 + TOLERANCE || cosphi < (-1.0 - TOLERANCE)) 
      {
         int me;
         MPI_Comm_rank(world,&me);
         if (screen) {
            char str[128];
            sprintf(str,"Improper problem: %d " BIGINT_FORMAT " %d %d %d %d",
               me,update->ntimestep,atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
            error->warning(FLERR,str,0);
            fprintf(screen,"  1st atom: %d %g %g %g\n",me,x[i1][0],x[i1][1],x[i1][2]);
            fprintf(screen,"  2nd atom: %d %g %g %g\n",me,x[i2][0],x[i2][1],x[i2][2]);
            fprintf(screen,"  3rd atom: %d %g %g %g\n",me,x[i3][0],x[i3][1],x[i3][2]);
            fprintf(screen,"  4th atom: %d %g %g %g\n",me,x[i4][0],x[i4][1],x[i4][2]);
            }
      }

      
      /* Apply corrections to round-off errors. */
      if (cosphi > 1.0)  cosphi -= SMALL;
      if (cosphi < -1.0) cosphi += SMALL;
      
      /* Calculate the angle: */
      double torangle = acos(cosphi);
      cosphi = cos(torangle - chi[type]);

      if (eflag) eimproper = 0.5 * k[type] * cosphi * cosphi;
     
      /*
      printf("The tags: %d-%d-%d-%d, of type %d .\n",atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4],type);   
      printf("The ji vector: %f, %f, %f.\nThe lk vector: %f, %f, %f.\n", vb1x,vb1y,vb1z,vb3x,vb3y,vb3z);
      printf("The cosine of the angle: %-1.16e.\n", cosphi);
      printf("The energy of the improper: %-1.16e with prefactor %-1.16e.\n", eimproper, 0.5*k[type]);
      */

      /* Work out forces. */
      angfac = - k[type] * cosphi;

      cjiji = rjisq;
      clklk = rlksq;
      /*CLKJI = RXLK * RXJI + RYLK * RYJI + RZLK * RZJI */
      clkji = vb3x*vb1x + vb3y*vb1y + vb3z*vb1z; 

      /*CFACT1 = CLKLK * CJIJI
        CFACT1 = SQRT(CFACT1)
        CFACT1 = ANGFAC / CFACT1*/
      cfact1 = angfac / sqrt(clklk * cjiji);
      /*CFACT2 = CLKJI / CLKLK*/
      cfact2 = clkji / clklk;
      /*CFACT3 = CLKJI / CJIJI*/
      cfact3 = clkji / cjiji;

      /*FIX = -RXLK + CFACT3 * RXJI
        FIY = -RYLK + CFACT3 * RYJI
        FIZ = -RZLK + CFACT3 * RZJI*/
      f1[0] = - vb3x + cfact3 * vb1x;
      f1[1] = - vb3y + cfact3 * vb1y;
      f1[2] = - vb3z + cfact3 * vb1z;

      /*FJX = -FIX
        FJY = -FIY
        FJZ = -FIZ*/
      f2[0] = - f1[0];
      f2[1] = - f1[1];
      f2[2] = - f1[2];
      
      /*FKX = CFACT2 * RXLK - RXJI
        FKY = CFACT2 * RYLK - RYJI
        FKZ = CFACT2 * RZLK - RZJI*/
      f3[0] = cfact2 * vb3x - vb1x;
      f3[1] = cfact2 * vb3y - vb1y;
      f3[2] = cfact2 * vb3z - vb1z;

      /*FLX = -FKX
        FLY = -FKY
        FLZ = -FKZ*/
      f4[0] = - f3[0];
      f4[1] = - f3[1];
      f4[2] = - f3[2];

      /*FIX = FIX * CFACT1
        FIY = FIY * CFACT1
        FIZ = FIZ * CFACT1*/
      f1[0] *= cfact1;
      f1[1] *= cfact1; 
      f1[2] *= cfact1;

      /*FJX = FJX * CFACT1
        FJY = FJY * CFACT1
        FJZ = FJZ * CFACT1*/
      f2[0] *= cfact1;
      f2[1] *= cfact1;
      f2[2] *= cfact1;

      /*FKX = FKX * CFACT1
        FKY = FKY * CFACT1
        FKZ = FKZ * CFACT1*/
      f3[0] *= cfact1; 
      f3[1] *= cfact1; 
      f3[2] *= cfact1;

      /*FLX = FLX * CFACT1
        FLY = FLY * CFACT1
        FLZ = FLZ * CFACT1*/
      f4[0] *= cfact1;
      f4[1] *= cfact1; 
      f4[2] *= cfact1;

      /* Apply force to each of 4 atoms */
      if (newton_bond || i1 < nlocal) {
         f[i1][0] += f1[0];
         f[i1][1] += f1[1];
         f[i1][2] += f1[2];
      }

      if (newton_bond || i2 < nlocal) {
         f[i2][0] += f2[0];
         f[i2][1] += f2[1];
         f[i2][2] += f2[2];
      }

      if (newton_bond || i3 < nlocal) {
         f[i3][0] += f3[0];
         f[i3][1] += f3[1];
         f[i3][2] += f3[2];
      }

      if (newton_bond || i4 < nlocal) {
         f[i4][0] += f4[0];
         f[i4][1] += f4[1];
         f[i4][2] += f4[2];
      }

      if (evflag)
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,eimproper,f1,f3,f4,
	       -vb1x,-vb1y,-vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
   }
}

/* ---------------------------------------------------------------------- */

void ImproperCossq::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(k,n+1,"improper:k");
  memory->create(chi,n+1,"improper:chi");

  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void ImproperCossq::coeff(int narg, char **arg)
{
   /* Check whether there exist sufficient number of arguments. 
      0: type of improper to be applied to
      1: energetic constant
      2: equilibrium angle in degrees */
   if (narg != 3) error->all(FLERR,"Incorrect args for cossq improper coefficients");
   if (!allocated) allocate();

   int ilo,ihi;
   force->bounds(arg[0],atom->nimpropertypes,ilo,ihi);

   double k_one = force->numeric(arg[1]);
   double chi_one = force->numeric(arg[2]);

   int count = 0;
   for (int i = ilo; i <= ihi; i++) {
      k[i] = k_one;
      chi[i] = ((chi_one * MY_PI)/180.0);
      setflag[i] = 1;
      count++;
   }

   if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */
void ImproperCossq::write_restart(FILE *fp)
{
   fwrite(&k[1],sizeof(double),atom->nimpropertypes,fp);
   fwrite(&chi[1],sizeof(double),atom->nimpropertypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */
void ImproperCossq::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&chi[1],sizeof(double),atom->nimpropertypes,fp);
  }
  MPI_Bcast(&k[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&chi[1],atom->nimpropertypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}
