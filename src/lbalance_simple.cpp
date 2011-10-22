/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com


LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov


Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.


See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "lbalance_simple.h"
#include "math.h"
#include "domain.h"
#include "mpi.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "neighbor.h"
#include "update.h"
#include "stdlib.h"
#include "string.h"

using namespace LAMMPS_NS;

/*NL*/ #define LMP_DEBUGMODE_LBALANCE_SIMPLE false
/*NL*/ #define LMP_DEBUGMODE_LBALANCE_SIMPLE_RESULTS false //(idim==2)//(update->ntimestep>7000)
/*NL*/ #define LMP_DEBUG_OUT_LBALANCE_SIMPLE screen

/* ---------------------------------------------------------------------- */

LbalanceSimple::LbalanceSimple(class LAMMPS *lmp,int narg, char **arg): Lbalance(lmp,narg,arg)
{
    //NP do not parse here for derived classes such as hybrid
    if(strncmp(style,"nlocal/",6)) return;

    if(narg < iarg+2)  error->all(FLERR,"Loadbalance simple: Not enough arguments");

    if(strcmp(arg[iarg++],"ntry")) error->all(FLERR,"Loadbalance nlocal/simple: Expecting 'ntry'");
    ntry_simple = atoi(arg[iarg++]);
    if(ntry_simple < 1) error->all(FLERR,"Loadbalance max: ntry too small");
    if(ntry_simple > 10) error->warning(FLERR,"Loadbalance max: ntry >10 might result in high comm cost");
}

/* ---------------------------------------------------------------------- */

LbalanceSimple::~LbalanceSimple() {}

/*NP ----------------------------------------------------------------------
   reset_box() has been called before, global box is set
   this function sets the local boxes and should yield an equal number of
   particles in each domain for low # of processes

   function is called each reneighboring if lbalance == 1

   for each dim, borders are calculated so atom distribution is equal

   possible problem: if atoms jump from proc i to i+2 in a dim
     b/c communication works on a stencil (domain->procneigh)
     atoms could jump if half skin is larger than smallest subhi-sublo
     or if neigh list build is off for a while

------------------------------------------------------------------------- */
void LbalanceSimple::loadbalance_local_boxes()
{
    // do not do anything if there is no reasonable # of particles in the system
    if(atom->natoms < 10 * comm->nprocs) return;

    if(domain->triclinic) error->all(FLERR,"Load balancing not implemented for triclinic boxes");

    /*NL*/ //if(comm->me == 0) fprintf(screen,"Loadbalancing: %f particles\n",atom->natoms);

    procgrid = comm->procgrid;
    myloc = comm->myloc;

    lodim[0] = lodim[1] = lodim[2] = 0;
    hidim[0] = procgrid[0]-1;
    hidim[1] = procgrid[1]-1;
    hidim[2] = procgrid[2]-1;

    loadbalance_local_boxes_simple(lodim,hidim);
}

void LbalanceSimple::loadbalance_local_boxes_simple(int *lodim, int *hidim)
{
  // new borders for load-balanced system
  double bal_sublo[3],bal_subhi[3];

  double border[3];
  int idim,idim_proc,ncount[3],ncount_ideal,nproc_grid;
  subhi = domain->subhi;
  sublo = domain->sublo;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;

  //NP minimum box extent equivalent to max cutoff
  minextent = 1.05 * cutneighmax();

  //NP count total particles in proc box to handle
  //NP equal to total particle number if proc stencil extends whole box
  int natoms = count_particles(0,boxhi[0]);

  //NP error if domain is too small to be loadbalanced
  //NP should not occur since minextent is accounted for in apply_border()
  for (int i = 0; i < 3; i++) 
      if(boxhi[i] - boxlo[i] < procgrid[i] * minextent) 
         error->all(FLERR,"Domain too small for this processor grid and this cutoff size: Enlarge domain, reduce # prcessors or choose smaller cutoff");

  /*NL*/ if(LMP_DEBUGMODE_LBALANCE_SIMPLE) fprintf(LMP_DEBUG_OUT_LBALANCE_SIMPLE,"minextent %f\n",minextent);

  for(idim = 0; idim < 3; idim++)
  {
        nproc_grid = hidim[idim] - lodim[idim] + 1;
        /*NL*/ //if(LMP_DEBUGMODE_LBALANCE_SIMPLE && idim == 2) fprintf(LMP_DEBUG_OUT_LBALANCE_SIMPLE,"Processor %d:, nlocal %d\n",comm->me,atom->nlocal);

        //NP loop proc stencil
        for(idim_proc = lodim[idim]; idim_proc <= hidim[idim]; idim_proc++)
        {
            //NP first calculate the maximum allowable border shift
            //NP this ensures that no proc is skipped - i.e. particles are always communicated to a neigh proc
            calc_max_shift(idim, idim_proc);

           //--------------------------
           //step 1- calculate lo bound
           //--------------------------

           //NP for the first proc - take boxlo
           if(idim_proc == 0) bal_sublo[idim] = boxlo[idim];
           //otherwise last hi limit as lo limit
           else bal_sublo[idim] = bal_subhi[idim];

           //--------------------------
           //step2 - calculate hi bound
           //--------------------------

           //NP for the last proc - take boxhi
           if(idim_proc == procgrid[idim]-1) bal_subhi[idim] = boxhi[idim];
           else
           {
              //NP ideal particle count for the slice
              ncount_ideal = (idim_proc - lodim[idim] + 1) * static_cast<int>(natoms / nproc_grid);

              /*NL*/ //if(LMP_DEBUGMODE_LBALANCE_SIMPLE) fprintf(LMP_DEBUG_OUT_LBALANCE_SIMPLE,"ncount_ideal %d\n",ncount_ideal);

              //NP init and perform recursive search
              border[0] = boxlo[idim];
              border[2] = boxhi[idim];
              ncount[0] = 0;
              ncount[2] = natoms;
              bal_subhi[idim] = calc_border(ntry_simple,idim,ncount_ideal,border,ncount);
           }

           //--------------------------
           //step3 - apply bounds
           //--------------------------
           // this may change the value for bal_subhi if necesary
           apply_border(bal_sublo,bal_subhi,idim,idim_proc);
        }
  }
  //NP error->all(FLERR,"loadbalance finished");
}

/*NP ----------------------------------------------------------------------
   recursive function to calc optimal domain decomposition
------------------------------------------------------------------------- */

double LbalanceSimple::calc_border(int ntry,int dim,int ncount_ideal,double *border,int *ncount)
{
  double btemp;
  int ntemp;

  border[1] = 0.5 * (border[0] + border[2]);
  ncount[1] = count_particles(dim,border[1]);

  /*NL*/ //if(LMP_DEBUGMODE_LBALANCE_SIMPLE) fprintf(LMP_DEBUG_OUT_LBALANCE_SIMPLE, "proc %d, iteration %d, border %f\n",comm->me,ntry,border[1]);

  //NP recursive binary search
  if(ntry > 0)
  {
      if(ncount[0] > ncount_ideal || ncount[2] < ncount_ideal) error->all(FLERR,"Illegal situation in LbalanceSimple::calc_border");
      if(ncount[1] == ncount_ideal) return border[1];

      //NP adjust hi and lo border accordingly
      else if(ncount[1] < ncount_ideal)
      {
         border[0] = border[1];
         ncount[0] = ncount[1];
      }
      else if(ncount[1] > ncount_ideal)
      {
         border[2] = border[1];
         ncount[2] = ncount[1];
      }
      return calc_border(ntry-1,dim,ncount_ideal,border,ncount);
  }

  //NP calc relative deviation for the three results
  double rel_dev[3];
  for(int i = 0; i < 3; i++)
     rel_dev[i] = fabs(static_cast<double>(ncount[i]-ncount_ideal) / static_cast<double>(ncount_ideal)) ;

  //NP return the border with minimum deviation
  if(rel_dev[0] < rel_dev[1] && rel_dev[0] < rel_dev[2]) return border[0];
  else if(rel_dev[1] < rel_dev[2]) return border[1];
  else return border[2];
}

/*NP ----------------------------------------------------------------------
   function that counts number of particles
   only count particles in the proc stencil
------------------------------------------------------------------------- */
inline int LbalanceSimple::count_particles(int dim,double border)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;

  int count = 0, count_all;

  //NP only take the right processors
  if( myloc[dim] < lodim[dim] || myloc[dim] > hidim[dim])
     count = 0;
  else
  {
    for(int i = 0; i < nlocal; i++)
      if(x[i][dim] < border) count++;
  }

  MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);
  /*NL*/ //if(LMP_DEBUGMODE_LBALANCE_SIMPLE) fprintf(LMP_DEBUG_OUT_LBALANCE_SIMPLE,"  proc %d counted %d particles for this border \n",comm->me,count_all);
  return count_all;
}

/*NP ----------------------------------------------------------------------
   calculate the maximum allowable shift
------------------------------------------------------------------------- */

void LbalanceSimple::calc_max_shift(int idim,int idim_proc)
{
            //NP get max shift distance for negative and positive direction = subbox extent in this dim for
            //NP this ensures that no proc is skipped - i.e. particles are always communicated to a neigh proc
            double max_shift_all[2];
            for(int i = 0; i < 2; i++)
            {
                if (myloc[idim] == idim_proc+i) max_shift[i] = subhi[idim] - sublo[idim];
                else max_shift[i] = BIG;
                MPI_Allreduce(&(max_shift[i]),&(max_shift_all[i]),1,MPI_DOUBLE,MPI_MIN,world);
                max_shift[i] = 0.90*max_shift_all[i];
            }
            //NP handle special cases - may not shift last border
            if(idim_proc == comm->procgrid[idim]-1) max_shift[1] = 0.;
}


/*NP ----------------------------------------------------------------------
   apply that has been calculated border
------------------------------------------------------------------------- */

void LbalanceSimple::apply_border(double *bal_sublo, double *bal_subhi,int idim,int idim_proc)
{
           double subhi_final,subhi_final_all, boxhi_stencil,boxhi_stencil_all,bal_subhi_max;

           //NP get boxhi for the stencil
           boxhi_stencil = - BIG;
           if(myloc[idim] <= hidim[idim]) boxhi_stencil = subhi[idim];
           MPI_Allreduce(&boxhi_stencil,&boxhi_stencil_all,1,MPI_DOUBLE,MPI_MAX,world);
           boxhi_stencil = boxhi_stencil_all;

           //NP prevent the case that choosing a border does not leave enough place
           //NP for the other procs to fulfil the minimum extent
           bal_subhi_max = boxhi_stencil - (hidim[idim] - idim_proc) * minextent;
           if(bal_subhi[idim] > bal_subhi_max) bal_subhi[idim] = bal_subhi_max;

           //NP only take the calculated subbox size if relevant for me
           subhi_final = BIG;
           if(myloc[idim] == idim_proc)
           {
              //NP can take sublo directly
              sublo[idim] = bal_sublo[idim];

              /*NL*/ if(LMP_DEBUGMODE_LBALANCE_SIMPLE) fprintf(LMP_DEBUG_OUT_LBALANCE_SIMPLE,"original border for proc %d of dimension %d: %f, new border %f, maxshift: %f %f\n",idim_proc,idim,subhi[idim],bal_subhi[idim],-max_shift[0],max_shift[1]);

              //NP make sure the minumim extent is obeyed
              if(bal_subhi[idim] - bal_sublo[idim] < minextent) bal_subhi[idim] = bal_sublo[idim] + minextent;

              /*NL*/ if(LMP_DEBUGMODE_LBALANCE_SIMPLE) fprintf(LMP_DEBUG_OUT_LBALANCE_SIMPLE,"    border result after min extent: %f\n",bal_subhi[idim]);

              //NP make sure border is not shifted more than allowed
              double diff = bal_subhi[idim] - subhi[idim];

              //if(idim == 2) fprintf(screen,"Processor %d: boundaries initally calculated for dim %d: %f / %f\n",comm->me,idim,bal_sublo[idim],bal_subhi[idim]);

              if     (diff < 0 && diff < -max_shift[0])
                  subhi[idim] = subhi[idim] - max_shift[0];
              else if(diff > 0 && diff >  max_shift[1])
                  subhi[idim] = subhi[idim] + max_shift[1];
              else subhi[idim] = bal_subhi[idim];

              subhi_final = subhi[idim];

              /*NL*/ if(LMP_DEBUGMODE_LBALANCE_SIMPLE) fprintf(LMP_DEBUG_OUT_LBALANCE_SIMPLE,"    border result after max shift extent: %f\n",subhi[idim]);

              //NP warn if maximum shift distance did override
              if(subhi[idim] - sublo[idim] < minextent ) error->warning(FLERR,"Minimum sub-domain extent could not be obeyed because particles would have been lost (did you insert large particles?). Inaccuracies may result.");

              //if(idim ==2) fprintf(screen,"Processor %d: boundaries finally calculated for dim %d: %f / %f\n",comm->me,idim,sublo[idim],subhi[idim]);
           }

           //NP comunicate subhi that came out of the calculation to all procs
           MPI_Allreduce(&subhi_final,&subhi_final_all,1,MPI_DOUBLE,MPI_MIN,world);
           bal_subhi[idim] = subhi_final_all;

           /*NL*/ if(LMP_DEBUGMODE_LBALANCE_SIMPLE && comm->me == 0) fprintf(LMP_DEBUG_OUT_LBALANCE_SIMPLE,"FINAL result for proc %d of dimension %d: %f\n",idim_proc,idim,bal_subhi[idim]);
           /*NL*/ if(LMP_DEBUGMODE_LBALANCE_SIMPLE_RESULTS && comm->me == 0) fprintf(LMP_DEBUG_OUT_LBALANCE_SIMPLE,"FINAL result for proc %d of dimension %d: %f\n",idim_proc,idim,bal_subhi[idim]);
}
