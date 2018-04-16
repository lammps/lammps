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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)
   
   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018). 
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics 
   and molecular dynamics. arXiv preprint arXiv:1801.10233.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_nve_spin.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "pair_hybrid_overlay.h"
#include "pair_spin.h"
#include "pair_spin_dmi.h"
#include "pair_spin_exchange.h"
#include "pair_spin_neel.h"
#include "pair_spin_me.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSpin::PairSpin(LAMMPS *lmp) : Pair(lmp),
lockfixnvespin(NULL), pair_keyword(NULL), pair_spin_keywords(NULL),
exchange_spin_styles(NULL), dmi_spin_styles(NULL), 
neel_spin_styles(NULL), me_spin_styles(NULL)
{
  hbar = force->hplanck/MY_2PI;
  single_enable = 0;
  no_virial_fdotr_compute = 1;
  lattice_flag = 0;

  // init # of Pair/Spin styles

  nspinstyle = 0;		// # of PairSpin styles
  nexchangespinstyle = 0;	// # of PairSpinExchange styles
  ndmispinstyle = 0;		// # of PairSpinDmi styles
  nneelspinstyle = 0;		// # of PairSpinNeel styles
  nmespinstyle = 0;		// # of PairSpinMe styles

  // init larger Pair/Spin style cutoff

  larger_cutoff = 0.0;
}

/* ---------------------------------------------------------------------- */

PairSpin::~PairSpin()
{
  
  if (nspinstyle) {
    for (int m = 0; m < nspinstyle; m++) {
      delete [] pair_spin_keywords[m];
    }
  }
  delete [] pair_keyword;
  delete [] exchange_spin_styles;
  delete [] dmi_spin_styles;
  delete [] neel_spin_styles;
  delete [] me_spin_styles;
  delete [] pair_spin_keywords;

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpin::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect number of args in pair_style pair/spin command");

  // pair/spin need the metal unit style

  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"pair/spin style requires metal units");

}

/* ----------------------------------------------------------------------
   global compute, defined in Pair/Spin subclasses
------------------------------------------------------------------------- */

void PairSpin::compute(int eflag, int vflag) {}

/* ----------------------------------------------------------------------
   compute all Pair/Spin interactions for atom ii
------------------------------------------------------------------------- */

void PairSpin::compute_pair_single_spin(int ii, double fmi[3])
{
  const int nlocal = atom->nlocal;
  int *type = atom->type;
  double **x = atom->x;
  double **sp = atom->sp;

  double xi[3], rij[3], eij[3];
  double spi[3], spj[3];

  int iexchange, idmi, ineel, ime;
  int i,j,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double rsq, rd, inorm;

  iexchange = idmi = ineel = ime = 0;
  
  for (int ipair=0; ipair < nspinstyle; ipair++) {

    if (strstr(pair_spin_keywords[ipair],"pair/spin/exchange")) {
      inum = exchange_spin_styles[iexchange]->list->inum;
      ilist = exchange_spin_styles[iexchange]->list->ilist;
      numneigh = exchange_spin_styles[iexchange]->list->numneigh;
      firstneigh = exchange_spin_styles[iexchange]->list->firstneigh;

      i = ilist[ii];
      
      spi[0] = sp[i][0];
      spi[1] = sp[i][1];
      spi[2] = sp[i][2];
 
      xi[0] = x[i][0];
      xi[1] = x[i][1];
      xi[2] = x[i][2];
  
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {

	j = jlist[jj];
	j &= NEIGHMASK;
	itype = type[ii];
	jtype = type[j];

	spj[0] = sp[j][0];
	spj[1] = sp[j][1];
	spj[2] = sp[j][2];

	rij[0] = x[j][0] - xi[0];
	rij[1] = x[j][1] - xi[1];
	rij[2] = x[j][2] - xi[2];
	rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
	inorm = 1.0/sqrt(rsq);

	exchange_spin_styles[iexchange]->compute_exchange(i,j,rsq,fmi,spi,spj);
      }
      iexchange++;
    }

    if (strstr(pair_spin_keywords[ipair],"pair/spin/dmi")) {
      inum = dmi_spin_styles[idmi]->list->inum;
      ilist = dmi_spin_styles[idmi]->list->ilist;
      numneigh = dmi_spin_styles[idmi]->list->numneigh;
      firstneigh = dmi_spin_styles[idmi]->list->firstneigh;

      i = ilist[ii];
      
      spi[0] = sp[i][0];
      spi[1] = sp[i][1];
      spi[2] = sp[i][2];
 
      xi[0] = x[i][0];
      xi[1] = x[i][1];
      xi[2] = x[i][2];
  
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {

	j = jlist[jj];
	j &= NEIGHMASK;
	itype = type[ii];
	jtype = type[j];

	spj[0] = sp[j][0];
	spj[1] = sp[j][1];
	spj[2] = sp[j][2];

	rij[0] = x[j][0] - xi[0];
	rij[1] = x[j][1] - xi[1];
	rij[2] = x[j][2] - xi[2];
	rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
	inorm = 1.0/sqrt(rsq);

	dmi_spin_styles[idmi]->compute_dmi(i,j,rsq,fmi,spi,spj);
      }
      idmi++;
    }
    
    if (strstr(pair_spin_keywords[ipair],"pair/spin/neel")) {
      inum = neel_spin_styles[ineel]->list->inum;
      ilist = neel_spin_styles[ineel]->list->ilist;
      numneigh = neel_spin_styles[ineel]->list->numneigh;
      firstneigh = neel_spin_styles[ineel]->list->firstneigh;

      i = ilist[ii];
      
      spi[0] = sp[i][0];
      spi[1] = sp[i][1];
      spi[2] = sp[i][2];
 
      xi[0] = x[i][0];
      xi[1] = x[i][1];
      xi[2] = x[i][2];
  
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {

	j = jlist[jj];
	j &= NEIGHMASK;
	itype = type[ii];
	jtype = type[j];

	spj[0] = sp[j][0];
	spj[1] = sp[j][1];
	spj[2] = sp[j][2];

	rij[0] = x[j][0] - xi[0];
	rij[1] = x[j][1] - xi[1];
	rij[2] = x[j][2] - xi[2];
	rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
	inorm = 1.0/sqrt(rsq);
	eij[0] = rij[0]*inorm;
	eij[1] = rij[1]*inorm;
	eij[2] = rij[2]*inorm;

	neel_spin_styles[ineel]->compute_neel(i,j,rsq,eij,fmi,spi,spj);
      }
      ineel++;
    }
    
    if (strstr(pair_spin_keywords[ipair],"pair/spin/me")) {
      inum = me_spin_styles[ime]->list->inum;
      ilist = me_spin_styles[ime]->list->ilist;
      numneigh = me_spin_styles[ime]->list->numneigh;
      firstneigh = me_spin_styles[ime]->list->firstneigh;

      i = ilist[ii];
      
      spi[0] = sp[i][0];
      spi[1] = sp[i][1];
      spi[2] = sp[i][2];
 
      xi[0] = x[i][0];
      xi[1] = x[i][1];
      xi[2] = x[i][2];
  
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {

	j = jlist[jj];
	j &= NEIGHMASK;
	itype = type[ii];
	jtype = type[j];

	spj[0] = sp[j][0];
	spj[1] = sp[j][1];
	spj[2] = sp[j][2];

	rij[0] = x[j][0] - xi[0];
	rij[1] = x[j][1] - xi[1];
	rij[2] = x[j][2] - xi[2];
	rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
	inorm = 1.0/sqrt(rsq);
	eij[0] = rij[0]*inorm;
	eij[1] = rij[1]*inorm;
	eij[2] = rij[2]*inorm;

	me_spin_styles[ime]->compute_me(i,j,rsq,eij,fmi,spi,spj);
      }
      ime++;
    }

  }

}

/* ----------------------------------------------------------------------
   called from Fix/NVE/Spin
   initialize the # and the lists of Pair/Spin styles
------------------------------------------------------------------------- */

int PairSpin::init_pair()
{
  int nsub = 0;

  // getting the pair style

  pair_keyword = new char[strlen(force->pair_style) + 1];
  strcpy(pair_keyword,force->pair_style);


  // searching for number of PairSpin styles

  int temp_npair = 0;
  nspinstyle = 0;
  pair = pair_spin_match("spin",0,nsub);

  // init lists of PairSpin styles

  exchange_spin_styles = new PairSpinExchange*[nspinstyle];
  dmi_spin_styles = new PairSpinDmi*[nspinstyle];
  neel_spin_styles = new PairSpinNeel*[nspinstyle];
  me_spin_styles = new PairSpinMe*[nspinstyle];
  
  // init lists of PairSpin names

  pair_spin_keywords = new char*[nspinstyle];

  nexchangespinstyle = 0;
  ndmispinstyle = 0;
  nneelspinstyle = 0;
  nmespinstyle = 0;

  // loop to define lists of Pair/Spin styles

  int ispin = 0;
  if (strstr(pair_keyword,"spin")) {
    
    // Pair/Spin/Exchange style

    if (strstr(pair_keyword,"pair/spin/exchange")) {
      int n = strlen(pair_keyword) + 1;
      pair_spin_keywords[ispin] = new char[n];
      strcpy(pair_spin_keywords[ispin],pair_keyword);
      exchange_spin_styles[nexchangespinstyle] = (PairSpinExchange *) force->pair;
      nexchangespinstyle++;
      ispin++;
    }

    // Pair/Spin/Dmi style

    if (strstr(pair_keyword,"pair/spin/dmi")) {
      int n = strlen(pair_keyword) + 1;
      pair_spin_keywords[ispin] = new char[n];
      strcpy(pair_spin_keywords[ispin],pair_keyword);
      dmi_spin_styles[ndmispinstyle] = (PairSpinDmi *) force->pair;
      ndmispinstyle++;
      ispin++;
    }

    // Pair/Spin/Neel style

    if (strstr(pair_keyword,"pair/spin/neel")) {
      int n = strlen(pair_keyword) + 1;
      pair_spin_keywords[ispin] = new char[n];
      strcpy(pair_spin_keywords[ispin],pair_keyword);
      neel_spin_styles[nneelspinstyle] = (PairSpinNeel *) force->pair;
      nneelspinstyle++;
      ispin++;
    }
      
    // list Pair/Spin/Me styles

    if (strstr(pair_keyword,"pair/spin/me")) {
      int n = strlen(pair_keyword) + 1;
      pair_spin_keywords[ispin] = new char[n];
      strcpy(pair_spin_keywords[ispin],pair_keyword);
      me_spin_styles[nmespinstyle] = (PairSpinMe *) force->pair;
      nmespinstyle++;
      ispin++;
    }

  } else if (strstr(pair_keyword,"hybrid/overlay")) {		// if hybrid/overlay
    PairHybrid *lockhybrid = (PairHybrid *) force->pair;
    for (int i =0; i < lockhybrid->nstyles; i++) {
      
      // error checks

      if (strcmp(lockhybrid->keywords[i],"hybrid") == 0)
	error->all(FLERR,"Pair style hybrid cannot have hybrid as an argument");
      if (strcmp(lockhybrid->keywords[i],"none") == 0)
	error->all(FLERR,"Pair style hybrid cannot have none as an argument");

      // list Pair/Spin/Exchange styles

      if (strstr(lockhybrid->keywords[i],"pair/spin/exchange")) {
	int n = strlen(lockhybrid->keywords[i]) + 1;
	pair_spin_keywords[ispin] = new char[n];
	strcpy(pair_spin_keywords[ispin],lockhybrid->keywords[i]);
	exchange_spin_styles[nexchangespinstyle] = (PairSpinExchange *) lockhybrid->styles[i];
	nexchangespinstyle++;
	ispin++;
      }

      // list Pair/Spin/Dmi styles

      if (strstr(lockhybrid->keywords[i],"pair/spin/dmi")) {
	int n = strlen(lockhybrid->keywords[i]) + 1;
	pair_spin_keywords[ispin] = new char[n];
	strcpy(pair_spin_keywords[ispin],lockhybrid->keywords[i]);
	dmi_spin_styles[ndmispinstyle] = (PairSpinDmi *) lockhybrid->styles[i];
	ndmispinstyle++;
	ispin++;
      }
    
      // list Pair/Spin/Neel styles

      if (strstr(lockhybrid->keywords[i],"pair/spin/neel")) {
	int n = strlen(lockhybrid->keywords[i]) + 1;
	pair_spin_keywords[ispin] = new char[n];
	strcpy(pair_spin_keywords[ispin],lockhybrid->keywords[i]);
	neel_spin_styles[nneelspinstyle] = (PairSpinNeel *) lockhybrid->styles[i];
	nneelspinstyle++;
	ispin++;
      }
    
      // list Pair/Spin/Me styles

      if (strstr(lockhybrid->keywords[i],"pair/spin/me")) {
	int n = strlen(lockhybrid->keywords[i]) + 1;
	pair_spin_keywords[ispin] = new char[n];
	strcpy(pair_spin_keywords[ispin],lockhybrid->keywords[i]);
	me_spin_styles[nmespinstyle] = (PairSpinMe *) lockhybrid->styles[i];
	nmespinstyle++;
	ispin++;
      }

    }
  } else if (strstr(pair_keyword,"hybrid")) { 		// no hybrid style with PairSpin
    error->all(FLERR,"Pair/Spin styles need hybrid/overlay Pair style");
  } else error->all(FLERR,"Wrong arguments in PairSpin style"); 
 
  if (ispin != nspinstyle) 
    error->all(FLERR,"Wrong number of PairSpin styles");

  if ((nexchangespinstyle + ndmispinstyle + nneelspinstyle + nmespinstyle) != nspinstyle)
    error->all(FLERR,"Wrong number of PairSpin styles");

  if (strstr(pair_keyword,"spin") && nspinstyle != 1)
    error->all(FLERR,"Wrong number of PairSpin styles");

  int mag_pair_flag = 0; 
  if (nspinstyle >= 1) mag_pair_flag = 1;


  // get larger_cutoff

  larger_cutoff = larger_spin_cutoff();

  if (mag_pair_flag == 1 && larger_cutoff == 0.0)
    error->all(FLERR,"Wrong arguments for PairSpin styles");

  return mag_pair_flag;
}

/* ----------------------------------------------------------------------
   get the larger Pair/Spin style cutoff for the sectoring operation
------------------------------------------------------------------------- */

double PairSpin::larger_spin_cutoff()
{
  int iexchange, idmi, ineel, ime;
  double local_cutoff = 0.0;
 
  iexchange = idmi = ineel = ime = 0;

  for (int ipair=0; ipair < nspinstyle; ipair++) {

    if (strstr(pair_spin_keywords[ipair],"pair/spin/exchange")) {
      if (local_cutoff < exchange_spin_styles[iexchange]->cut_spin_exchange_global)
	local_cutoff = exchange_spin_styles[iexchange]->cut_spin_exchange_global;
      iexchange++;
    }

    if (strstr(pair_spin_keywords[ipair],"pair/spin/dmi")) {
      if (local_cutoff < dmi_spin_styles[idmi]->cut_spin_dmi_global)
	local_cutoff = dmi_spin_styles[idmi]->cut_spin_dmi_global;
      idmi++;
    }

    if (strstr(pair_spin_keywords[ipair],"pair/spin/neel")) {
      if (local_cutoff < neel_spin_styles[ineel]->cut_spin_neel_global)
	local_cutoff = neel_spin_styles[ineel]->cut_spin_neel_global;
      ineel++;
    }
    
    if (strstr(pair_spin_keywords[ipair],"pair/spin/me")) {
      if (local_cutoff < me_spin_styles[ime]->cut_spin_me_global)
	local_cutoff = me_spin_styles[ime]->cut_spin_me_global;
      ime++;
    }

  }

  if ((iexchange + idmi + ineel + ime) != nspinstyle)
    error->all(FLERR,"Wrong number of PairSpin styles");

  return local_cutoff;
}


/* ---------------------------------------------------------------------- */

Pair *PairSpin::pair_spin_match(const char *word, int exact, int nsub)
{
  int iwhich,count;


  if (exact && strcmp(pair_keyword,word) == 0) return pair;
  else if (!exact && strstr(pair_keyword,word)) return pair;

  else if (strstr(pair_keyword,"hybrid/overlay")) {
    PairHybrid *lockhybrid = (PairHybrid *) force->pair;
    count = 0;
    for (int i = 0; i < lockhybrid->nstyles; i++)
      if ((exact && strcmp(lockhybrid->keywords[i],word) == 0) ||
          (!exact && strstr(lockhybrid->keywords[i],word))) {
        iwhich = i;
        count++;
	nspinstyle = count;
        if (nsub == count) return lockhybrid->styles[iwhich];
      }
    if (count == 1) return lockhybrid->styles[iwhich];

  } else if (strstr(pair_keyword,"hybrid")) {
    PairHybrid *lockhybrid = (PairHybrid *) force->pair;
    count = 0;
    for (int i = 0; i < lockhybrid->nstyles; i++)
      if ((exact && strcmp(lockhybrid->keywords[i],word) == 0) ||
          (!exact && strstr(lockhybrid->keywords[i],word))) {
        iwhich = i;
        count++;
	nspinstyle = count;
        if (nsub == count) return lockhybrid->styles[iwhich];
      }
    if (count == 1) return lockhybrid->styles[iwhich];
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSpin::allocate()
{
  allocated = 1;
  char *pair_keyword = new char[strlen(force->pair_style) + 1];
  
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type spin pairs (only one for now)
------------------------------------------------------------------------- */

void PairSpin::coeff(int narg, char **arg)
{
  const double hbar = force->hplanck/MY_2PI;

  if (!allocated) allocate();
 
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpin::init_style()
{
  if (!atom->sp_flag)
    error->all(FLERR,"Pair spin requires atom/spin style");

  // checking if nve/spin is a listed fix

  int ifix = 0;
  while (ifix < modify->nfix) {
    if (strcmp(modify->fix[ifix]->style,"nve/spin") == 0) break;
    ifix++;
  }
  if (ifix == modify->nfix)
    error->all(FLERR,"pair/spin style requires nve/spin");

  // get the lattice_flag from nve/spin

  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"nve/spin") == 0) {
      lockfixnvespin = (FixNVESpin *) modify->fix[i];
      lattice_flag = lockfixnvespin->lattice_flag;
    }
  }
}
