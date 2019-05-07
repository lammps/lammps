/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   contributing authors: Steven E Strong and Nicholas J Hestand
   contact: stevene.strong at gmail dot com
------------------------------------------------------------------------- */

#include "pair_e3b.h"

#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "domain.h"
#include "citeme.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

//these are defined here to avoid confusing hardcoded indicies, but
//they do not allow flexibility. If they are changed the code will break
#define DIM 3
#define NUMH 2  //number of hydrogen atoms per water molecule
#define NUMO 2  //number of oxygen atoms per pair of water molecules
#define BOND_DELTA 1.01   //buffer for OH bonds that aren't perfectly rigid

using namespace LAMMPS_NS;
/* ---------------------------------------------------------------------- */

PairE3B::PairE3B(LAMMPS *lmp) : Pair(lmp),pairPerAtom(10)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;

  allocatedE3B = false;
  pairO  = NULL;
  pairH  = NULL;
  exps   = NULL;
  del3   = NULL;
  fpair3 = NULL;
  sumExp = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairE3B::~PairE3B()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
  if (allocatedE3B) {
    memory->destroy(pairO);
    memory->destroy(pairH);
    memory->destroy(exps);
    memory->destroy(del3);
    memory->destroy(fpair3);
    memory->destroy(sumExp);
  }
}

/* ---------------------------------------------------------------------- */

void PairE3B::compute(int eflag, int vflag)
{
  int i,j,k,h,ii,jj,hh,kk,inum,jnum,otherO;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair,rsq,tmpexp;
  double fxtmp,fytmp,fztmp,fix,fiy,fiz;
  double delxh,delyh,delzh,rsqh,tmpr;
  double scFact1,scFact2,scEng,scDer;
  int *ilist,*jlist,*numneigh,**firstneigh;
  bool addedH;

  if (natoms != atom->natoms)
    error->all(FLERR,"pair E3B requires a fixed number of atoms");

  //clear sumExp array
  memset(sumExp,0.0,nbytes);

  evdwl = 0.0;
  etot[0]=etot[1]=etot[2]=etot[3]=0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  int npair = 0;
  // loop over half neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (type[i]!=typeO)
      continue;

    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fix = fiy = fiz = 0.0;

    // two-body interactions
    jlist = firstneigh[i];
    jnum  = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      //skip unless O-O interaction
      if (type[j]!=typeO)
        continue;

      jtag = tag[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq  = delx*delx + dely*dely + delz*delz; //OO distance

      //two body interaction
      //not shifted b/c k2=4.87/A, so at cutoff (5.2A) e^(-kr) = 1e-11
      if (rsq < rc2sq) {
	tmpr    = sqrt(rsq);
	tmpexp  = e2 * exp(-k2*tmpr);
	fpair   = k2 * tmpexp / tmpr;

	fxtmp    = delx*fpair;
	fytmp    = dely*fpair;
	fztmp    = delz*fpair;
	fix     += fxtmp;
	fiy     += fytmp;
	fiz     += fztmp;
	f[j][0] -= fxtmp;
	f[j][1] -= fytmp;
	f[j][2] -= fztmp;

	if (evflag) {
	  ev_tally(i,j,nlocal,newton_pair,tmpexp,0.0,fpair,delx,dely,delz);
	  etot[3] += tmpexp;
	}
      } //end if rsq<rc2sq

      //accumulate info about each OH pair for later 3body stuff
      //test OO distance with augmented cutoff to account for dangling Hs
      if (rsq < rc3deltaSq) {
	//pairO and pairH are set here even if no Hs are within the cutoff
	//in that case, npair is not incremented and they will be overwritten
	pairO[npair][0] = i;
	pairO[npair][1] = j;
	addedH = false;

	for (kk=0; kk<NUMO; kk++) {
	  k  = pairO[npair][kk];
	  otherO = pairO[npair][(kk+1)%2];
	  for (hh=0; hh<NUMH; hh++) {
	    h=atom->map(tag[otherO]+hh+1);
	    //if hydrogen atom is missing, bond potential or shake will
	    //catch this, so don't need to check here
	    //if (h<0)
	    //  error->one(FLERR,"hydrogen atom missing");
	    h = domain->closest_image(otherO,h);
	    pairH[npair][kk][hh] = h;

	    delxh = x[k][0] - x[h][0];
	    delyh = x[k][1] - x[h][1];
	    delzh = x[k][2] - x[h][2];
	    rsqh = delxh*delxh + delyh*delyh + delzh*delzh;

	    if (rsqh < rc3sq) {

	      tmpr         = sqrt(rsqh);
	      tmpexp       = exp(-k3*tmpr);
	      if (tmpr > rs) {
		scFact1    = rc3-tmpr;
		scFact2    = sc_num + 2*tmpr;
		scEng      = scFact1*scFact1*scFact2*sc_denom;
		scDer      = k3*scEng - 6*scFact1*(rs-tmpr)*sc_denom;
	      } else  {
		scDer = k3;
		scEng = 1.0;
	      }

	      //need to keep fpair3 separate from del3 for virial
	      fpair3[npair][kk][hh]  = scDer*tmpexp/tmpr;
	      tmpexp                *= scEng;
	      exps[npair][kk][hh]    = tmpexp;
	      del3[npair][kk][hh][0] = delxh;
	      del3[npair][kk][hh][1] = delyh;
	      del3[npair][kk][hh][2] = delzh;

	      //accumulate global vector of sum(e^kr)
	      //tags start at 1, so subtract one to index sumExp
	      sumExp[tag[k]-1] += tmpexp;
	      sumExp[tag[h]-1] += tmpexp;

	      addedH = true;
	    } else {
	      exps  [npair][kk][hh] = 0.0;
	      fpair3[npair][kk][hh] = 0.0;
	    } //if < rc3sq
	  } //end loop through 2 Hs
	} //end for kk in NUMO
	//if added a pair, check if array is too small and grow
	if (addedH) {
	  npair++;
	  if (npair >= pairmax)
	    error->one(FLERR,"neigh is too small");
	}
      } //end if < rc3deltaSq
    } //end for jj neigh

    //add 2-body forces on i
    f[i][0] += fix;
    f[i][1] += fiy;
    f[i][2] += fiz;
  } //end for ii

  //communicate sumExp array
  //tested that no change in speed with MPI_IN_PLACE
  MPI_Allreduce(MPI_IN_PLACE,sumExp,maxID,MPI_DOUBLE,MPI_SUM,world);

  //now loop through list of pairs, calculating 3body forces
  int j2,otherH;
  double partA,partB,partC;
  for (ii = 0; ii < npair; ii++) {

    for (kk=0; kk<NUMO; kk++) {
      i      = pairO[ii][kk];
      otherO = (kk+1)%2;
      partB = eb*(  sumExp[tag[pairO[ii][otherO]   ]-1]
		    + sumExp[tag[pairH[ii][otherO][0]]-1]
		    + sumExp[tag[pairH[ii][otherO][1]]-1]
		    - 2*(exps[ii][otherO][0] + exps[ii][otherO][1]));
      partC = ec*(sumExp[tag[i]-1] - exps[ii][kk][0] - exps[ii][kk][1]);

      for (hh=0; hh<NUMH; hh++) {
	j  = pairH[ii][kk][hh];

	//type A
	otherH = (hh+1)%2;
	j2 = pairH[ii][kk][otherH];

	partA = ea*(sumExp[tag[j2]-1] - exps[ii][kk][otherH]); //not full energy yet
	fpair = partA*fpair3[ii][kk][hh];
	fxtmp = fpair*del3[ii][kk][hh][0];
	fytmp = fpair*del3[ii][kk][hh][1];
	fztmp = fpair*del3[ii][kk][hh][2];

	f[i][0] += fxtmp;
	f[i][1] += fytmp;
	f[i][2] += fztmp;
	f[j][0] -= fxtmp;
	f[j][1] -= fytmp;
	f[j][2] -= fztmp;

	if (evflag) {
	  evdwl = partA*exps[ii][kk][hh]*0.5; //mult by exp on this H
	  ev_tally(i,j,nlocal,newton_pair,evdwl, 0.0,fpair,
		   del3[ii][kk][hh][0],del3[ii][kk][hh][1],del3[ii][kk][hh][2]);
	  etot[0] += evdwl;
	}

	//type B
	fpair = partB*fpair3[ii][kk][hh];
	fxtmp = fpair*del3[ii][kk][hh][0];
	fytmp = fpair*del3[ii][kk][hh][1];
	fztmp = fpair*del3[ii][kk][hh][2];

	f[i][0] += fxtmp;
	f[i][1] += fytmp;
	f[i][2] += fztmp;
	f[j][0] -= fxtmp;
	f[j][1] -= fytmp;
	f[j][2] -= fztmp;

	if (evflag) {
	  evdwl = partB*exps[ii][kk][hh]*0.5; //mult by exp on this H
	  ev_tally(i,j,nlocal,newton_pair,evdwl, 0.0,fpair,
		   del3[ii][kk][hh][0],del3[ii][kk][hh][1],del3[ii][kk][hh][2]);
	  etot[1] += evdwl;
	}

	//type C
	fpair = partC*fpair3[ii][kk][hh];
	fxtmp = fpair*del3[ii][kk][hh][0];
	fytmp = fpair*del3[ii][kk][hh][1];
	fztmp = fpair*del3[ii][kk][hh][2];

	f[i][0] += fxtmp;
	f[i][1] += fytmp;
	f[i][2] += fztmp;
	f[j][0] -= fxtmp;
	f[j][1] -= fytmp;
	f[j][2] -= fztmp;

	if (evflag) {
	  evdwl = partC*exps[ii][kk][hh]*0.5; //mult by exp on this H
	  ev_tally(i,j,nlocal,newton_pair,evdwl, 0.0,fpair,
		   del3[ii][kk][hh][0],del3[ii][kk][hh][1],del3[ii][kk][hh][2]);
	  etot[2] += evdwl;
	}
      } //end for hh in NUMH
    } //end for kk in NUMO
  } //end for ii in npairs

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairE3B::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

void PairE3B::allocateE3B()
{
  allocatedE3B = true;

  //TODO: get memory->grow working for 4d arrays
  pairmax = atom->nlocal*pairPerAtom;    //initial guess for size of pair lists
  memory->create(pairO ,pairmax,NUMO         ,"pair:pairO");
  memory->create(pairH ,pairmax,NUMO,NUMH    ,"pair:pairH");
  memory->create(exps  ,pairmax,NUMO,NUMH    ,"pair:exps");
  memory->create(fpair3,pairmax,NUMO,NUMH    ,"pair:fpair3");
  memory->create(del3  ,pairmax,NUMO,NUMH,DIM,"pair:del3");

  //set del3 to zero to silence valgrind memcheck errors
  //don't need to do this in every call to compute() because we set
  //exps and fpair3 to zero, and all uses of del3 are multiplied by one of those
  int ii,jj,kk,ll;
  for (ii=0; ii<pairmax; ii++)
    for (jj=0; jj<NUMO; jj++)
      for (kk=0; kk<NUMH; kk++)
	for (ll=0; ll<DIM; ll++)
	  del3[ii][jj][kk][ll] = 0.0;

  natoms=atom->natoms;
  maxID=find_maxID();
  if (!natoms)
    error->all(FLERR,"No atoms found");
  memory->create(sumExp,maxID,"pair:sumExp");
  nbytes = sizeof(double) * maxID;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairE3B::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");
  typeO=force->inumeric(FLERR,arg[0]);
  if (typeO<1 || typeO>atom->ntypes)
    error->all(FLERR,"Invalid Otype: out of bounds");
}

/* ----------------------------------------------------------------------
   coeffs must be * * keyword/value
   keyword/values set the potential parameters
------------------------------------------------------------------------- */
void PairE3B::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  //1=* 2=* 3/4=1st keyword/value
  if (narg < 4)
    error->all(FLERR,"There must be at least one keyword given to pair_coeff");

  // ensure I,J args are * *
  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

 // clear setflag since coeff() called once with I,J = * *
  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  setflag[typeO][typeO]=1;

  //parse keyword/value pairs
  double bondL=0.0; //OH bond length
  bool repeatFlag=false;
  int presetFlag;

  //clear parameters
  e2=ea=eb=ec=k3=k2=NAN;
  rs=rc3=rc2=0.0;

  int iarg = 2; //beginning of keyword/value pairs
  while(iarg < narg) {
    char *keyword = arg[iarg++];
    if (checkKeyword(keyword,"Ea",1,narg-iarg))
      ea=force->numeric(FLERR,arg[iarg++]);
    else if (checkKeyword(keyword,"Eb",1,narg-iarg))
      eb=force->numeric(FLERR,arg[iarg++]);
    else if (checkKeyword(keyword,"Ec",1,narg-iarg))
      ec=force->numeric(FLERR,arg[iarg++]);
    else if (checkKeyword(keyword,"K3",1,narg-iarg))
      k3=force->numeric(FLERR,arg[iarg++]);
    else if (checkKeyword(keyword,"Rs",1,narg-iarg))
      rs=force->numeric(FLERR,arg[iarg++]);
    else if (checkKeyword(keyword,"Rc3",1,narg-iarg))
      rc3=force->numeric(FLERR,arg[iarg++]);
    else if (checkKeyword(keyword,"Rc2",1,narg-iarg))
      rc2=force->numeric(FLERR,arg[iarg++]);
    else if (checkKeyword(keyword,"bondL",1,narg-iarg))
      bondL=force->numeric(FLERR,arg[iarg++]);
    else if (checkKeyword(keyword,"E2",1,narg-iarg))
      e2=force->numeric(FLERR,arg[iarg++]);
    else if (checkKeyword(keyword,"K2",1,narg-iarg))
      k2=force->numeric(FLERR,arg[iarg++]);
    else if (checkKeyword(keyword,"neigh",1,narg-iarg))
      pairPerAtom=force->inumeric(FLERR,arg[iarg++]);
    else if (checkKeyword(keyword,"preset",1,narg-iarg)) {
      presetFlag=force->inumeric(FLERR,arg[iarg++]);
      presetParam(presetFlag,repeatFlag,bondL);
    } else {
      char str[256];
      snprintf(str,256,"Keyword %s is unknown",keyword);
      error->all(FLERR,str);
    }
  }

  checkInputs(bondL);

  //cutmax for neighbor listing
  cutmax = std::max(rc2,rc3);
  rc2sq  = rc2*rc2;
  rc3sq  = rc3*rc3;
  rc3deltaSq = (rc3+bondL)*(rc3+bondL);

  double tmpfact=1.0/(rc3-rs);
  sc_denom=tmpfact*tmpfact*tmpfact;
  sc_num=rc3-3*rs;
}

/* ----------------------------------------------------------------------
   init specific to this pair styles
------------------------------------------------------------------------- */

void PairE3B::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style E3B requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style E3B requires newton pair on");

  // need a half neighbor list
  int irequest = neighbor->request(this,instance_me);
  //don't need this, half is default
  //neighbor->requests[irequest]->half = 0;

  if (!force->pair_match("tip4p",false,0))
    if (comm->me==0) error->warning(FLERR,"E3B pair_style is designed for use with hybrid/overlay tip4p style");

  if (!allocatedE3B) allocateE3B();
}


static const char cite_E3B1[] =
  "Explicit Three-Body (E3B) potential for water:\n\n"
  "@article{kumar_water_2008,\n"
  "title = {Water Simulation Model with Explicit Three-Molecule Interactions},\n"
  "volume = {112},\n"
  "doi = {10.1021/jp8009468},\n"
  "number = {28},\n"
  "journal = {J Phys. Chem. B},\n"
  "author = {Kumar, R. and Skinner, J. L.},\n"
  "year = {2008},\n"
  "pages = {8311--8318}\n"
  "}\n\n";

static const char cite_E3B2[] =
  "Explicit Three-Body (E3B) potential for water:\n\n"
  "@article{tainter_robust_2011,\n"
  "title = {Robust three-body water simulation model},\n"
  "volume = {134},\n"
  "doi = {10.1063/1.3587053},\n"
  "number = {18},\n"
  "journal = {J. Chem. Phys},\n"
  "author = {Tainter, C. J. and Pieniazek, P. A. and Lin, Y.-S. and Skinner, J. L.},\n"
  "year = {2011},\n"
  "pages = {184501}\n"
  "}\n\n";

static const char cite_E3B3[] =
  "Explicit Three-Body (E3B) potential for water:\n\n"
  "@article{tainter_reparametrized_2015,\n"
  "title = {Reparametrized {E3B} (Explicit Three-Body) Water Model Using the {TIP4P/2005} Model as a Reference},\n"
  "volume = {11},\n"
  "doi = {10.1021/acs.jctc.5b00117},\n"
  "number = {5},\n"
  "journal = {J. Chem. Theory Comput.},\n"
  "author = {Tainter, Craig J. and Shi, Liang and Skinner, James L.},\n"
  "year = {2015},\n"
  "pages = {2268--2277}\n"
  "}\n\n";

void PairE3B::presetParam(const int flag,bool &repeatFlag,double &bondL) {
  if (repeatFlag) {
    error->all(FLERR,
	       "Cannot request two different sets of preset parameters");
  }
  repeatFlag=true;

  if (!std::isnan(ea) || !std::isnan(eb) || !std::isnan(ec) ||
      !std::isnan(e2) || !std::isnan(k3) || !std::isnan(k2) ||
      bondL!=0.0 || rs!=0.0 || rc3!=0.0 || rc2!=0.0 )
    error->all(FLERR,"Preset keyword will overwrite another keyword setting");

  double econv,lconv;
  if (strcmp(update->unit_style,"real") == 0) {
    econv=1.0/4.184;
    lconv=1.0;
  } else if (strcmp(update->unit_style,"metal") == 0) {
    econv=0.103653271;
    lconv=1.0;
  } else if (strcmp(update->unit_style,"si") == 0) {
    econv=1.660578e-21;
    lconv=1e-10;
  } else if (strcmp(update->unit_style,"cgs") == 0) {
    econv=1.660578e-14;
    lconv=1e-8;
  } else {
    char str[256];
    snprintf(str,256,
	     "Pre-defined E3B parameters have not been set for %s units.",
	     update->unit_style);
      error->all(FLERR,str);
  }

  //here parameters are defined in kJ/mol and A
  //they will be converted to the lammps units after
  if (flag==2008) {
    error->all(FLERR,"\"preset 2008\" is not yet supported, because this would require distinct k3 coefficients, use \"preset 2011\" or \"preset 2015\"");
    if (lmp->citeme) lmp->citeme->add(cite_E3B1);
    ea  = 4699.6;
    eb  =-2152.9;
    ec  = 1312.7;
    //ka  = 1.0/1.88;
    //kb  = 1.0/1.71;
    //kc  = 1.0/1.56;
    e2  = 1.925e6;
    k2  = 4.67;
    rs  = 5.0;
    rc3 = 5.2;
    rc2 = 5.2;
    bondL = 0.9572;
  } else if (flag==2011) {
    if (lmp->citeme) lmp->citeme->add(cite_E3B2);
    ea  = 1745.7;
    eb  =-4565.0;
    ec  = 7606.8;
    k3  = 1.907;
    e2  = 2.349e6;
    k2  = 4.872;
    rs  = 5.0;
    rc3 = 5.2;
    rc2 = 5.2;
    bondL = 0.9572;
  } else if (flag==2015) {
    if (lmp->citeme) lmp->citeme->add(cite_E3B3);
    ea  = 150.0;
    eb  =-1005.0;
    ec  = 1880.0;
    k3  = 1.907;
    e2  = 0.453e6;
    k2  = 4.872;
    rs  = 5.0;
    rc3 = 5.2;
    rc2 = 5.2;
    bondL = 0.9572;
  } else
    error->all(FLERR,"Unknown argument: preset only takes 2011 or 2015 as arguments");

  //convert units
  ea *= econv;
  eb *= econv;
  ec *= econv;
  e2 *= econv;
  k3 /= lconv;
  k2 /= lconv;
  rs *= lconv;
  rc2 *= lconv;
  rc3 *= lconv;
  bondL *= lconv*BOND_DELTA;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
//pair.cpp::init uses this to set cutsq array, used for neighboring, etc
double PairE3B::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

bool PairE3B::checkKeyword(const char *thiskey,const char *test,const int nVal, const int nRem) {
  if (strcmp(thiskey,test) == 0) {
    if(nRem<nVal) {
      char str[256];
      snprintf(str,256,"Too few arguments to \"%s\" keyword.",test);
      error->all(FLERR,str);
    }
    return true;
  }
  return false;
}

//find max atom ID for all atoms
//from fix_deposit.cpp
tagint PairE3B::find_maxID()
{
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  tagint max = 0;
  tagint maxID;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&maxID,1,MPI_LMP_TAGINT,MPI_MAX,world);

  return maxID;
}

void PairE3B::checkInputs(const double &bondL) {
  //first check that all necessary values were set
  if (rc2==0.0)
    error->all(FLERR,"rc2 keyword missing");
  if (rs==0.0)
    error->all(FLERR,"Rs keyword missing");
  if (rc3==0.0)
    error->all(FLERR,"Rc3 keyword missing");
  if (bondL==0.0)
    error->all(FLERR,"bondL keyword missing");
  if (std::isnan(ea))
    error->all(FLERR,"Ea keyword missing");
  if (std::isnan(eb))
    error->all(FLERR,"Eb keyword missing");
  if (std::isnan(ec))
    error->all(FLERR,"Ec keyword missing");
  if (std::isnan(k3))
    error->all(FLERR,"K3 keyword missing");
  if (std::isnan(e2))
    error->all(FLERR,"E2 keyword missing");
  if (std::isnan(k2))
    error->all(FLERR,"K2 keyword missing");

  //now test that values are within acceptable ranges
  if (k2 < 0.0 or k3 < 0.0)
    error->all(FLERR,"exponential decay is negative");
  if (bondL<0.0)
    error->all(FLERR,"OH bond length is negative");
  if (rc2 < 0.0 || rc3 < 0.0 || rs < 0.0)
    error->all(FLERR,"potential cutoff is negative");
  if (rs > rc3)
    error->all(FLERR,"potential switching distance is larger than cutoff");
  if (rs==rc3)
    error->warning(FLERR,"potential switching distance is equal to cutoff: this is untested and not conserve energy");
  if (pairPerAtom<0)
    error->all(FLERR,"neigh is negative");
}
