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
   Contributing authors: Aidan Thompson (SNL), Hansohl Cho (MIT)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "pair_reax.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

#include "../REAX/reax_fortran.h"
#include "../REAX/reax_params.h"
#include "../REAX/reax_cbkc.h"
#include "../REAX/reax_cbkd.h"
#include "../REAX/reax_cbkch.h"
#include "../REAX/reax_cbkabo.h"
#include "../REAX/reax_cbkia.h"
#include "../REAX/reax_cbkpairs.h"
#include "../REAX/reax_energies.h"
#include "../REAX/reax_small.h"
#include "../REAX/reax_functions.h"


using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SMALL 0.0001

/* ---------------------------------------------------------------------- */

PairREAX::PairREAX(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  one_coeff = 1;
  
  Lmidpoint = false;
  cutmax=0.0;
  hbcut=6.0;
  iprune=4;
  ihb=1;

  chpot = 0;

  // Set communication size needed by ReaxFF
  comm_forward = 2;
  comm_reverse = 1;

  precision = 1.0e-6;

}

/* ----------------------------------------------------------------------
   free all arrays
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairREAX::~PairREAX()
{
  // memory->sfree(w);
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

  }
}

/* ---------------------------------------------------------------------- */

void PairREAX::compute(int eflag, int vflag)
{
  
  double evdwl, ecoul;
  double reax_energy_pieces[13];

  double energy_charge_equilibration;

  evdwl = 0.0;
  ecoul = 0.0;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int newton_pair = force->newton_pair;

  if (eflag)
    {
      for (int k = 0; k < 14; k++)
	{
	  reax_energy_pieces[k]=0;
	}
    }

  if (vflag)
    {
      FORTRAN(cbkvirial, CBKVIRIAL).Lvirial = 1;
    }
  else
    {
      FORTRAN(cbkvirial, CBKVIRIAL).Lvirial = 0;
    }

  if (evflag)
    {
      FORTRAN(cbkvirial, CBKVIRIAL).Lvirial = 1;
      FORTRAN(cbkvirial, CBKVIRIAL).Latomvirial = 1;
    }
  else
    {
      FORTRAN(cbkvirial, CBKVIRIAL).Latomvirial = 0;
    }

  // Call compute_charge() to calculate the atomic charge distribution
  compute_charge(energy_charge_equilibration);

  // Call write_reax methods to transfer LAMMPS positions and neighbor lists
  // into REAX Fortran positions and Verlet lists
  write_reax_positions();
  write_reax_vlist();

  // Determine whether this bond is owned by the processor or not. 
  FORTRAN(srtbon1, SRTBON1)(&iprune, &ihb, &hbcut);

  // Need to communicate with other processors for the atomic bond order calculations
  FORTRAN(cbkabo, CBKABO).abo;
  // Communicate the local atomic bond order in this processor into the ghost atomic bond order in other processors
  comm -> comm_pair(this);

  FORTRAN(molec, MOLEC)();
  FORTRAN(encalc, ENCALC)();

  int node = comm -> me;
  FORTRAN(mdsav, MDSAV)(&node);

  // Read in the forces from ReaxFF Fortran
  read_reax_forces();

  // Read in the reax energy pieces from ReaxFF Fortran
  if (eflag)
    {
      reax_energy_pieces[0] = FORTRAN(cbkenergies, CBKENERGIES).eb;
      reax_energy_pieces[1] = FORTRAN(cbkenergies, CBKENERGIES).ea;
      reax_energy_pieces[2] = FORTRAN(cbkenergies, CBKENERGIES).elp;
      reax_energy_pieces[3] = FORTRAN(cbkenergies, CBKENERGIES).emol;
      reax_energy_pieces[4] = FORTRAN(cbkenergies, CBKENERGIES).ev;
      reax_energy_pieces[5] = FORTRAN(cbkenergies, CBKENERGIES).epen;
      reax_energy_pieces[6] = FORTRAN(cbkenergies, CBKENERGIES).ecoa;
      reax_energy_pieces[7] = FORTRAN(cbkenergies, CBKENERGIES).ehb;
      reax_energy_pieces[8] = FORTRAN(cbkenergies, CBKENERGIES).et;
      reax_energy_pieces[9] = FORTRAN(cbkenergies, CBKENERGIES).eco;
      reax_energy_pieces[10] = FORTRAN(cbkenergies, CBKENERGIES).ew;
      reax_energy_pieces[11] = FORTRAN(cbkenergies, CBKENERGIES).ep;
      reax_energy_pieces[12] = FORTRAN(cbkenergies, CBKENERGIES).efi;

      for (int k = 0; k < 13; k++)
        {
          evdwl += reax_energy_pieces[k];
        }

      // eVDWL energy in LAMMPS does not include the Coulomb energy in ReaxFF
      evdwl = evdwl - reax_energy_pieces[11];
      // eCOUL energy in LAMMPS does include the Coulomb energy and charge equilibation energy based on the calculated charge distribution in ReaxFF
      ecoul = reax_energy_pieces[11]+energy_charge_equilibration;

      // Call the global energy pieces at this step
      eng_vdwl += evdwl;
      eng_coul += ecoul;
    }

  int nval, ntor, nhb, nbonall, nbon, na, na_local, nvpair;
  int reaxsize[8];

  na_local = FORTRAN(rsmall, RSMALL).na_local;
  na = FORTRAN(rsmall, RSMALL).na;
  FORTRAN(getnbonall, GETNBONALL)(&nbonall);
  nbon = FORTRAN(rsmall, RSMALL).nbon;
  FORTRAN(getnval, GETNVAL)(&nval);
  FORTRAN(getntor, GETNTOR)(&ntor);
  FORTRAN(getnhb, GETNHB)(&nhb);
  nvpair = FORTRAN(cbkpairs, CBKPAIRS).nvpair;

  reaxsize[0] = na_local;
  reaxsize[1] = na;
  reaxsize[2] = nbonall;
  reaxsize[3] = nbon;
  reaxsize[4] = nval;
  reaxsize[5] = ntor;
  reaxsize[6] = nhb;
  reaxsize[7] = nvpair;
 
  // Need to call ev_tally to update the global energy and virial

  if (vflag)
    {
      virial[0] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[0];
      virial[1] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[1];
      virial[2] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[2];
      virial[3] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[3];
      virial[4] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[4];
      virial[5] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[5];
    }

  if (vflag_atom)
    {
      read_reax_atom_virial();
    }
}

/* ---------------------------------------------------------------------- */

void PairREAX::write_reax_positions()
{
  // Write atomic positions used in ReaxFF Fortran
  // Copy from atomic position data in LAMMPS into ReaxFF Fortran
  double **x = atom->x;
  // Copy calculated charge distribution into ReaxFF Fortran
  double *q = atom->q;
  int *type = atom->type;
  int *tag = atom->tag;
  double xtmp, ytmp, ztmp;
  int itype;
  int itag; 
  int j, jx, jy, jz, jia;
  int nlocal, nghost;
  nlocal = atom->nlocal;
  nghost = atom->nghost;

  FORTRAN(rsmall, RSMALL).na = nlocal+nghost;
  FORTRAN(rsmall, RSMALL).na_local = nlocal;

  j = 0;
  jx = 0;
  jy = ReaxParams::nat;
  jz = 2*ReaxParams::nat;
  jia = 0;

  for (int i = 0; i < nlocal+nghost; i++)
    {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      itag = tag[i];
      FORTRAN(cbkc, CBKC).c[j+jx] = xtmp;
      FORTRAN(cbkc, CBKC).c[j+jy] = ytmp;
      FORTRAN(cbkc, CBKC).c[j+jz] = ztmp;
      FORTRAN(cbkch, CBKCH).ch[j] = q[i];
      FORTRAN(cbkia, CBKIA).ia[j+jia] = itype;
      FORTRAN(cbkia, CBKIA).iag[j+jia] = itype;
      FORTRAN(cbkc, CBKC).itag[j]= itag;
      j++;
    }

}

void PairREAX::write_reax_vlist()
{
  double **x = atom->x;
  int *type = atom->type;
  int *tag = atom->tag;
 
  int nlocal, nghost;
  nlocal = atom->nlocal;
  nghost = atom->nghost;

  int inum, jnum;
  int *ilist;
  int *jlist;
  int *numneigh, **firstneigh;

  double xitmp, yitmp, zitmp;
  double xjtmp, yjtmp, zjtmp;

  int itype, jtype, itag, jtag;

  int ii, jj, i, j, iii, jjj;

  int nvpair, nvlself, nvpairmax;
  int nbond;
  double delr2;
  double delx, dely, delz;
  double rtmp[3];

  nvpairmax = ReaxParams::nneighmax * ReaxParams::nat;

  nvpair = 0;
  nvlself =0;
  nbond = 0;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

   for (ii = 0; ii < inum; ii++)
    {
      i = ilist[ii]; // the local index for the ii th atom having neighbors
      xitmp = x[i][0];
      yitmp = x[i][1];
      zitmp = x[i][2];
      itype = type[i];
      itag = tag[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++)
	{
	  j = jlist[jj];
	  xjtmp = x[j][0];
	  yjtmp = x[j][1];
	  zjtmp = x[j][2];
	  jtype = type[j];
	  jtag = tag[j];

	  delx = xitmp - xjtmp;
	  dely = yitmp - yjtmp;
	  delz = zitmp - zjtmp;

	  delr2 = delx*delx+dely*dely+delz*delz;

	  if (delr2 <= rcutvsq)
	    {
	      // Avoid the double check for the vpairs
	      if (i < j)
		{
		  // Index for the FORTRAN array
		  iii = i+1;
		  jjj = j+1;
		}
	      else
		{
		  iii = j+1;
		  jjj = i+1;
		}
	      if (nvpair >= nvpairmax)
		{
		  break;
		}
	      FORTRAN(cbkpairs, CBKPAIRS).nvl1[nvpair] = iii;
	      FORTRAN(cbkpairs, CBKPAIRS).nvl2[nvpair] = jjj;
	      FORTRAN(cbknvlbo, CBKNVLBO).nvlbo[nvpair] = 0;
	      if (delr2 <= rcutbsq)
		{
		  FORTRAN(cbknvlbo, CBKNVLBO).nvlbo[nvpair] = 1;
		  nbond++;
		}

	      FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 0;

	      // Midpoint check
	      if (Lmidpoint)
		{
		/*	{
		  if (jjj <= nlocal)
		    {
		      FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
		    }
		  else
		    {
		      rtmp[0] = 0.5*(xitmp + xjtmp);
		      rtmp[1] = 0.5*(yitmp + yjtmp);
		      rtmp[2] = 0.5*(zitmp + zjtmp);
		      if (sub_check_strict(rtmp))
			{
			  FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] =1;
			}
		    }
		    }*/
		}
	      else
		{    
		  if (i < nlocal)
		    {
		      if (j < nlocal)
			{
			  FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
			}
		      else if (itag < jtag)
			{
			  FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
			}
		      else if (itag == jtag)
			{
			  if (delz > SMALL)
			    {
			      FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
			    }
			  else if (fabs(delz) < SMALL)
			    {
			      if (dely > SMALL)
				{
				  FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
				}
			      else if (fabs(dely) < SMALL)
				{
				  if (delx > SMALL)
				    {
				      FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
				    }
				}
			    }
			}
		    }
		}
	      nvpair++;
	    }
	}
    }

  for (int i = nlocal; i < nlocal + nghost; i++)
    {
      xitmp = x[i][0];
      yitmp = x[i][1];
      zitmp = x[i][2];
      itype = type[i];
      itag = tag[i];
      for (int j = i+1; j < nlocal + nghost; j++)
	{
	  xjtmp = x[j][0];
	  yjtmp = x[j][1];
	  zjtmp = x[j][2];
	  jtype = type[j];
	  jtag = tag[j];

	  delx = xitmp - xjtmp;
	  dely = yitmp - yjtmp;
	  delz = zitmp - zjtmp;

	  delr2 = delx*delx+dely*dely+delz*delz;

	  // Don't need to check the double count since i < j in the ghost region
	  if (delr2 <= rcutvsq)
	    {
	      {
		iii = i+1;
		jjj = j+1;
	      }

	      if (nvpair >= nvpairmax)
		{
		  break;
		}
	      FORTRAN(cbkpairs, CBKPAIRS).nvl1[nvpair] = iii;
	      FORTRAN(cbkpairs, CBKPAIRS).nvl2[nvpair] = jjj;
	      FORTRAN(cbknvlbo, CBKNVLBO).nvlbo[nvpair] = 0;
	      if (delr2 <= rcutbsq)
		{
		  FORTRAN(cbknvlbo, CBKNVLBO).nvlbo[nvpair] = 1;
		  nbond++;
		}

	      FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 0;

	      if (Lmidpoint)
		{
		  /*		  if (jjj <= nlocal)
		    {
		      FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
		    }
		  else
		    {
		      rtmp[0] = 0.5*(xitmp + xjtmp);
		      rtmp[1] = 0.5*(yitmp + yjtmp);
		      rtmp[2] = 0.5*(zitmp + zjtmp);
		      if (sub_check_strict(rtmp))
			{
			  FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] =1;
			}
			}*/
		}
	      
	      else
		{
		  if (i < nlocal)
		    {
		      if (j < nlocal)
			{
			  FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
			}
		      else if (itag < jtag)
			{
			  FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
			}
		      else if (itag == jtag)
			{
			  if (delz > SMALL)
			    {
			      FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
			    }
			  else if (fabs(delz) < SMALL)
			    {
			      if (dely > SMALL)
				{
				  FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
				}
			      else if (fabs(dely) < SMALL)
				{
				  if (delx > SMALL)
				    {
				      FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
				    }
				}
			    }
			}
		    }
		}
	      nvpair++;
	    }
	}
    }

  FORTRAN(cbkpairs, CBKPAIRS).nvpair = nvpair;
  FORTRAN(cbkpairs, CBKPAIRS).nvlself = nvlself;
}
						      



void PairREAX::read_reax_forces()
{
  double **f = atom->f;
  double ftmp[3];
  int j, jx, jy, jz;
  int nlocal, nghost;
  nlocal = atom->nlocal;
  nghost = atom->nghost;

  j = 0;
  jx = 0;
  jy = 1;
  jz = 2;

  for (int i = 0; i < nlocal + nghost; i++)
    {
      ftmp[0] = -FORTRAN(cbkd, CBKD).d[j+jx];
      ftmp[1] = -FORTRAN(cbkd, CBKD).d[j+jy];
      ftmp[2] = -FORTRAN(cbkd, CBKD).d[j+jz];
      f[i][0] = ftmp[0];
      f[i][1] = ftmp[1];
      f[i][2] = ftmp[2];
      j += 3;
    }
}

void PairREAX::read_reax_atom_virial()
{
  error->warning("read_reax_atom_virial");

  double vatomtmp[6];
  int j;
  int nlocal, nghost;
  nlocal = atom->nlocal;
  nghost = atom->nghost;

  j = 0;

  for (int i =0; i < nlocal + nghost; i++)
    {
      vatomtmp[0] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+0];
      vatomtmp[1] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+1];
      vatomtmp[2] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+2];
      vatomtmp[3] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+3];
      vatomtmp[4] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+4];
      vatomtmp[5] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+5];

      vatom[i][0] = vatomtmp[0];
      vatom[i][1] = vatomtmp[1];
      vatom[i][2] = vatomtmp[2];
      vatom[i][3] = vatomtmp[3];
      vatom[i][4] = vatomtmp[4];
      vatom[i][5] = vatomtmp[5];

      j += 6;
    }
}

/* ---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */

void PairREAX::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairREAX::settings(int narg, char **arg)
{

  if (narg != 0 && narg !=2) error->all("Illegal pair_style command");
  
  if (narg == 2) {
    hbcut = atof (arg[0]); // User-specifed hydrogen-bond cutoff
    precision = atof (arg[1]); // User-specified charge equilibration precision
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairREAX::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all("Incorrect args for pair coefficients");

  // Clear setflag since coeff() called once with I,J = * *

  int  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass for i,i in atom class

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      {
	setflag[i][j] = 1;
	count++;
      }

  if (count == 0) error->all("Incorrect args for pair coefficients");
 
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairREAX::init_style()
{
  // Need a half neighbor list for REAX
  // if (force->newton_pair == 0)
  // error->all("Pair interactions in ReaxFF require newton pair on");

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 1;
  neighbor->requests[irequest]->full = 0;

  // Initial setup for ReaxFF
  int node = comm->me;

  FORTRAN(readc, READC)();
  FORTRAN(init, INIT)();
  FORTRAN(ffinpt, FFINPT)();
  FORTRAN(tap7th, TAP7TH)();

  // Need to turn off read_in by fort.3 in REAX Fortran
  int ngeofor_tmp = -1;
  FORTRAN(setngeofor, SETNGEOFOR)(&ngeofor_tmp);
  if (node == 0)
    {
      FORTRAN(readgeo, READGEO)();
    }

  // Initial setup for cutoff radius of VLIST and BLIST in ReaxFF

  double vlbora;
  
  FORTRAN(getswb, GETSWB)(&swb);
  cutmax=MAX(swb, hbcut);
  rcutvsq=cutmax*cutmax;
  FORTRAN(getvlbora, GETVLBORA)(&vlbora);
  rcutbsq=vlbora*vlbora;

  // Parameters for charge equilibration from ReaxFF input, fort.4
  FORTRAN(getswa, GETSWA)(&swa);
  ff_params ff_param_tmp;
  double chi, eta, gamma;
  int nparams;
  nparams = 5;
  // FORTRAN(getnso, GETNSO)(&ntypes);
  param_list = new ff_params[atom->ntypes+1];
  for (int itype = 1; itype <= atom->ntypes; itype++)
    {
      chi = FORTRAN(cbkchb, CBKCHB).chi[itype-1];
      eta = FORTRAN(cbkchb, CBKCHB).eta[itype-1];
      gamma = FORTRAN(cbkchb, CBKCHB).gam[itype-1];
      ff_param_tmp.np = nparams;
      ff_param_tmp.rcutsq = cutmax;
      ff_param_tmp.params = new double[nparams];
      ff_param_tmp.params[0] = chi;
      ff_param_tmp.params[1] = eta;
      ff_param_tmp.params[2] = gamma;
      ff_param_tmp.params[3] = swa;
      ff_param_tmp.params[4] = swb;
      param_list[itype] = ff_param_tmp;
    }

  taper_setup();

  // Need to turn off newton_pair flag for the neighbor list for ReaxFF
  // In ReaxFF, local-ghost pairs should be stored in both processors
  force -> newton_pair = 0;
  force -> newton = 1;
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairREAX::init_one(int i, int j)
{

  return cutmax;
}


/* ---------------------------------------------------------------------- */

int PairREAX::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) 
    {
      j = list[i];
      buf[m++] = FORTRAN(cbkabo, CBKABO).abo[j];

      buf[m++] = w[j];
    }
  return 2;
}

/* ---------------------------------------------------------------------- */

void PairREAX::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) 
     {
       FORTRAN(cbkabo, CBKABO).abo[i] = buf[m++];
       w[i] = buf[m++];
     } 
}

/* ---------------------------------------------------------------------- */

int PairREAX::pack_reverse_comm(int n, int first, double *buf)
{
  int i,k,m,last,size;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) 
    {
      buf[m++] = w[i];
    }

  return 1;
}

/* ---------------------------------------------------------------------- */

void PairREAX::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,k,m;

  m = 0;
  for (i = 0; i < n; i++) 
    {
    j = list[i];
    w[j] += buf[m++];
    }
}
 

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------Charge Equilibation Caculation------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

void PairREAX::taper_setup()
{
  double swb2,swa2,swb3,swa3,d1,d7;

  d1=swb-swa;
  d7=pow(d1,7);
  swa2=swa*swa;
  swa3=swa2*swa;
  swb2=swb*swb;
  swb3=swb2*swb;
 
  swc7=  20.0e0/d7;
  swc6= -70.0e0*(swa+swb)/d7;
  swc5=  84.0e0*(swa2+3.0e0*swa*swb+swb2)/d7;
  swc4= -35.0e0*(swa3+9.0e0*swa2*swb+9.0e0*swa*swb2+swb3)/d7;
  swc3= 140.0e0*(swa3*swb+3.0e0*swa2*swb2+swa*swb3)/d7;
  swc2=-210.0e0*(swa3*swb2+swa2*swb3)/d7;
  swc1= 140.0e0*swa3*swb3/d7;
  swc0=(-35.0e0*swa3*swb2*swb2+21.0e0*swa2*swb3*swb2+
	7.0e0*swa*swb3*swb3+swb3*swb3*swb)/d7;
}

double PairREAX::taper_E(const double &r, const double &r2)
{
  double r3;
  r3=r2*r;
  return swc7*r3*r3*r+swc6*r3*r3+swc5*r3*r2+swc4*r2*r2+swc3*r3+swc2*r2+
     swc1*r+swc0;
}

double PairREAX::taper_F(const double &r, const double &r2)
{
  double r3;
  r3=r2*r;
  return 7.0e0*swc7*r3*r3+6.0e0*swc6*r3*r2+5.0e0*swc5*r2*r2+
    4.0e0*swc4*r3+3.0e0*swc3*r2+2.0e0*swc2*r+swc1;
}

// Compute current charge distributions based on the charge equilibration
void PairREAX::compute_charge(double &energy_charge_equilibration)
{
  
  double **x = atom->x;
  int *type = atom->type;
  int *tag = atom->tag;

  double *q = atom->q;

  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  int inum, jnum;
  int *ilist;
  int *jlist;
  int *numneigh, **firstneigh;

  double xitmp, yitmp, zitmp;
  double xjtmp, yjtmp, zjtmp;

  int itype, jtype, itag, jtag;

  int ii, jj, i, j;
  double delr2, delr_norm, gamt, hulp1, hulp2;
  double delx, dely, delz;
  int node, nprocs;

  double qsum, qi;
  double *ch;
  double *elcvec;
  double *aval;
  int *arow_ptr;
  int *acol_ind;
  int maxmatentries, nmatentries;
  double sw;
  double rtmp[3];

  node = comm -> me;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int numneigh_total = 0;
  for (ii =0; ii < inum; ii++)
    {
      numneigh_total += numneigh[ilist[ii]];
    }

  arow_ptr = new int[nlocal+nghost+1];
  ch = new double[nlocal+nghost+1];
  elcvec = new double[nlocal+nghost+1];
  maxmatentries = numneigh_total+2*nlocal;
  aval = new double[maxmatentries];
  acol_ind = new int[maxmatentries];
  nmatentries = 0;

  // Construct a linear system for this processor
  for (ii =0; ii<inum; ii++)
    {
      i = ilist[ii];
      xitmp = x[i][0];
      yitmp = x[i][1];
      zitmp = x[i][2];
      itype = type[i];
      itag = tag[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // Caution : index of i and ii
      arow_ptr[i] = nmatentries;
      aval[nmatentries] = 2.0*param_list[itype].params[1];
      // Caution : index of i and ii
      acol_ind[nmatentries] = i;

      nmatentries++;

      aval[nmatentries] = 1.0;
      acol_ind[nmatentries] = nlocal + nghost;
      nmatentries++;   
	  
      for (jj = 0; jj < jnum; jj++)
	{
	  j = jlist[jj];
	  xjtmp = x[j][0];
	  yjtmp = x[j][1];
	  zjtmp = x[j][2];
	  jtype = type[j];
	  jtag = tag[j];

	  delx = xitmp - xjtmp;
	  dely = yitmp - yjtmp;
	  delz = zitmp - zjtmp;
	  
	  delr2 = delx*delx+dely*dely+delz*delz;

	  // We did construct half neighbor list with no Newton flag for ReaxFF
	  // However, in Charge Equilibration, local-ghost pair should be stored only in one processor
	  // So, filtering is necessary when local-ghost pair is counted twice
	  neighbor->style;
	  if (j >= nlocal) 
	    {
	      if (neighbor->style==0)
		{  
		  if (itag > jtag) 
		    {
		      if ((itag+jtag) % 2 == 0) continue;
		    } 
		  else if (itag < jtag) 
		    {
		      if ((itag+jtag) % 2 == 1) continue;
		    }
		  else
	      // if (itag == jtag)
		    {
		      if (zjtmp < zitmp) continue;
		      else if (zjtmp == zitmp && yjtmp < yitmp) continue;
		      else if (zjtmp == zitmp && yjtmp == yitmp && xjtmp < xitmp)
			continue;
		    } 
		}
	      else if (neighbor->style==1)
		{
		  if (zjtmp < zitmp) continue;
		  else if (zjtmp == zitmp && yjtmp < yitmp) continue;
		  else if (zjtmp == zitmp && yjtmp == yitmp && xjtmp < xitmp)
		    continue;
		}
	      else
		{
		  error->all("ReaxFF does not support multi style neiborlist");
		}
	    }

 
	  // rcutvsq = cutmax*cutmax, in ReaxFF
	  if (delr2 <= rcutvsq)
	    {
	      gamt = sqrt(param_list[itype].params[2]*
			  param_list[jtype].params[2]);
	      delr_norm = sqrt(delr2);
	      sw = taper_E(delr_norm, delr2);	
	      hulp1=(delr_norm*delr2+(1.0/(gamt*gamt*gamt)));
	      hulp2=sw*14.40/cbrt(hulp1);
	      aval[nmatentries] = hulp2;
	      acol_ind[nmatentries] = j;
	      nmatentries++;
	    }
	}
    }

  // In this case, we don't use Midpoint method
  // So, we don't need to consider ghost-ghost interactions
  // But, need to fill the arow_ptr[] arrays for the ghost atoms
  for (i = nlocal; i < nlocal+nghost; i++)
    {
      arow_ptr[i] = nmatentries;
    }

  arow_ptr[nlocal+nghost] = nmatentries;

  // Add rhs matentries to linear system
 
  for (ii =0; ii<inum; ii++)
    {
      i = ilist[ii];
      itype = type[i];
      elcvec[i] = -param_list[itype].params[0];
    }

  for (i = nlocal; i < nlocal+nghost; i++)
    {
      elcvec[i] = 0.0;
    }

  // Assign current charges to charge vector
  qsum = 0.0;
  for (ii =0; ii<inum; ii++)
    {
      i = ilist[ii];
      qi = q[i];
      ch[i] = qi;
      if (i<nlocal)
	{
	  qsum += qi;
	}
    }

  for (i = nlocal; i < nlocal+nghost; i++)
    {
      qi = q[i];
      ch[i] = qi;
      if (i<nlocal)
	{
	  qsum += qi;
	}
    }

  double qtot;

  MPI_Allreduce(&qsum, &qtot, 1, MPI_DOUBLE, MPI_SUM, world);
  elcvec[nlocal+nghost] = 0.0;
  ch[nlocal+nghost] = chpot;

  // Solve the linear system using CG sover
  charge_reax(nlocal, nghost, ch, aval, acol_ind, arow_ptr, elcvec);

  // Calculate the charge equilibration energy
  energy_charge_equilibration = 0;

  // We already have the updated charge distributions for the current structure
  for (ii = 0; ii < inum; ii++)
    {
      i = ilist[ii];
      itype = type[i];

      // 23.02 is the ReaxFF conversion from eV to kcal/mol
      // Should really use constants.evfactor ~ 23.06, but
      // that would break consistency with serial ReaxFF code.
      qi = 23.02 * (param_list[itype].params[0]*ch[i]+
		    param_list[itype].params[1]*ch[i]*ch[i]);
      energy_charge_equilibration += qi;
    }

  // Copy charge vector back to particles from the calculated values

  for (i = 0; i < nlocal+nghost; i++)
    {
      q[i] = ch[i];
    }
  chpot = ch[nlocal+nghost];

  delete []aval;
  delete []arow_ptr;
  delete []elcvec;
  delete []ch;
  delete []acol_ind;

}

// Call the CG solver
void PairREAX::charge_reax(const int & nlocal, const int & nghost,
			   double ch[], double aval[], int acol_ind[],
			   int arow_ptr[], double elcvec[])
{
  double chpottmp, suma;
  double sumtmp;

  cg_solve(nlocal, nghost, aval, acol_ind, arow_ptr, ch, elcvec);
}

// Conjugate gradient solver for linear systems
void PairREAX::cg_solve(const int & nlocal, const int & nghost, 
			double aval[], int acol_ind[], int arow_ptr[],
			double x[], double b[])

{
  double one, zero, rho, rho_old, alpha, beta, gamma;
  int iter, maxiter;
  int n;
  double sumtmp;
  double *r;
  double *p; 
  //  double *w;
  double *p_old;
  double *q;

  //  Sketch of parallel CG method by A. P. Thompson
  //
  //  Distributed (partial) vectors: b, r, q, A
  //  Accumulated (full) vectors: x, w, p
  //
  //  r = b-A.x
  //  w = r            /* (ReverseComm + Comm) */
  // 
  
  r = new double[nlocal+nghost+1];
  p = new double[nlocal+nghost+1];
  w = new double[nlocal+nghost+1];
  p_old = new double[nlocal+nghost+1];
  q = new double[nlocal+nghost+1];

  n = nlocal+nghost+1;

  one = 1.0;
  zero = 0.0;
  maxiter = 100;

  // For w1, need to zero
  for (int i = 0; i < n; i++)
    {
      w[i] = 0;
    }
  
  // Construct r = b-Ax
  sparse_product(n, nlocal, nghost, aval, acol_ind, arow_ptr, x, r);
  
  // We will not use BLAS library

  for (int i=0; i<n; i++)
    {
      r[i] = b[i] - r[i];
      w[i] = r[i];
    }
  

  comm -> reverse_comm_pair(this);
  comm -> comm_pair(this);
  
  MPI_Allreduce(&w[n-1], &sumtmp, 1, MPI_DOUBLE, MPI_SUM, world);
  w[n-1] = sumtmp;
  rho_old = one;

   for (iter = 1; iter < maxiter; iter++)
    {
      rho = 0.0;
      for (int i=0; i<nlocal; i++)
	{
	  rho += w[i]*w[i];
	}

      MPI_Allreduce(&rho, &sumtmp, 1, MPI_DOUBLE, MPI_SUM, world);
      rho = sumtmp + w[n-1]*w[n-1];

      if (rho < precision) break;

      for (int i = 0; i<n; i++)
	{
	  p[i] = w[i];
	}

      if (iter > 1)
	{
	  beta = rho/rho_old;
	  for (int i = 0; i<n; i++)
	    {
	      p[i] += beta*p_old[i];
	    }
	}

      sparse_product(n, nlocal, nghost, aval, acol_ind, arow_ptr, p, q);
    
      gamma = 0.0;
      for (int i=0; i<n; i++)
	{
	  gamma += p[i]*q[i];

	}

      MPI_Allreduce(&gamma, &sumtmp, 1, MPI_DOUBLE, MPI_SUM, world);

      gamma = sumtmp;
      alpha = rho/gamma;

      for (int i=0; i<n; i++)
	{
	  x[i] += alpha*p[i];
	  r[i] -= alpha*q[i];
	  w[i] = r[i];
	}
      
      comm -> reverse_comm_pair(this);
      comm -> comm_pair(this);

      MPI_Allreduce(&w[n-1], &sumtmp, 1, MPI_DOUBLE, MPI_SUM, world);
      w[n-1] = sumtmp;

      for (int i=0; i<n; i++)
	{
	  p_old[i] = p[i];
	}

      rho_old = rho;
    }
  
  delete []r;
  delete []p;
  delete []p_old;
  delete []q; 
}

// Sparse maxtrix operations
void PairREAX::sparse_product(const int & n, const int & nlocal,
			      const int & nghost,
			      double aval[], int acol_ind[], int arow_ptr[],
			      double x[], double r[])
{
  int jj;
  for (int i=0; i<n; i++)
    {
      r[i] = 0.0;
    }

  // Loop over local particle matentries
  for (int i=0; i<nlocal; i++)
    {
      // Diagonal terms
      r[i]+=aval[arow_ptr[i]]*x[i];
      // Loop over remaining matrix entries and transposes
      for (int j=arow_ptr[i]+1; j<arow_ptr[i+1]; j++)
	{
	  jj = acol_ind[j];
	  r[i] += aval[j]*x[jj];
	  r[jj] += aval[j]*x[i];
	}
    }

  // Loop over ghost particle matentries
  for (int i=nlocal; i<nlocal+nghost; i++)
    {
      for (int j=arow_ptr[i]; j<arow_ptr[i+1]; j++)
	{
	  jj = acol_ind[j];
	  r[i] += aval[j]*x[jj];
	  r[jj] += aval[j]*x[i];
	}
    }

}
