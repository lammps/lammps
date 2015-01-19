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
   Contributing authors: Aidan Thompson (Sandia, athomps@sandia.gov)
                         Hansohl Cho (MIT, hansohl@mit.edu)
   LAMMPS implementation of the Reactive Force Field (ReaxFF) is based on
     Aidan Thompson's GRASP code
       (General Reactive Atomistic Simulation Program)
     and Ardi Van Duin's original ReaxFF code
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_reax.h"
#include "pair_reax_fortran.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define SMALL 0.0001

/* ---------------------------------------------------------------------- */

PairREAX::PairREAX(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  no_virial_fdotr_compute = 1;

  nextra = 14;
  pvector = new double[nextra];

  cutmax = 0.0;
  hbcut = 6.0;
  ihbnew = 1;
  itripstaball = 1;
  iprune = 4;
  ihb = 1;
  chpot = 0;

  nmax = 0;
  arow_ptr = NULL;
  ch = NULL;
  elcvec = NULL;
  rcg = NULL;
  wcg = NULL;
  pcg = NULL;
  poldcg = NULL;
  qcg = NULL;

  matmax = 0;
  aval = NULL;
  acol_ind = NULL;

  comm_forward = 1;
  comm_reverse = 1;

  precision = 1.0e-6;
}

/* ----------------------------------------------------------------------
   free all arrays
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairREAX::~PairREAX()
{
  delete [] pvector;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    for (int i = 1; i <= atom->ntypes; i++)
      delete [] param_list[i].params;
    delete [] param_list;

    delete [] map;
  }

  memory->destroy(arow_ptr);
  memory->destroy(ch);
  memory->destroy(elcvec);
  memory->destroy(rcg);
  memory->destroy(wcg);
  memory->destroy(pcg);
  memory->destroy(poldcg);
  memory->destroy(qcg);

  memory->destroy(aval);
  memory->destroy(acol_ind);
}

/* ---------------------------------------------------------------------- */

void PairREAX::compute(int eflag, int vflag)
{
  int i,j;
  double evdwl,ecoul;
  double energy_charge_equilibration;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else ev_unset();

  if (vflag_global) FORTRAN(cbkvirial, CBKVIRIAL).Lvirial = 1;
  else FORTRAN(cbkvirial, CBKVIRIAL).Lvirial = 0;

  if (vflag_atom) FORTRAN(cbkvirial, CBKVIRIAL).Latomvirial = 1;
  else FORTRAN(cbkvirial, CBKVIRIAL).Latomvirial = 0;

  // reallocate charge equilibration and CG arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(rcg);
    memory->destroy(wcg);
    memory->destroy(pcg);
    memory->destroy(poldcg);
    memory->destroy(qcg);

    nmax = atom->nmax;
    int n = nmax+1;

    memory->create(arow_ptr,n,"reax:arow_ptr");
    memory->create(ch,n,"reax:ch");
    memory->create(elcvec,n,"reax:elcvec");
    memory->create(rcg,n,"reax:rcg");
    memory->create(wcg,n,"reax:wcg");
    memory->create(pcg,n,"reax:pcg");
    memory->create(poldcg,n,"reax:poldcg");
    memory->create(qcg,n,"reax:qcg");
  }

  // calculate the atomic charge distribution

  compute_charge(energy_charge_equilibration);

  // transfer LAMMPS positions and neighbor lists to REAX

  write_reax_positions();
  write_reax_vlist();

  // determine whether this bond is owned by the processor or not

  FORTRAN(srtbon1, SRTBON1)(&iprune, &ihb, &hbcut, &ihbnew, &itripstaball);

  // communicate with other processors for the atomic bond order calculations

  FORTRAN(cbkabo, CBKABO).abo;

  // communicate local atomic bond order to ghost atomic bond order

  packflag = 0;
  comm->forward_comm_pair(this);

  FORTRAN(molec, MOLEC)();
  FORTRAN(encalc, ENCALC)();
  FORTRAN(mdsav, MDSAV)(&comm->me);

  // read forces from ReaxFF Fortran

  read_reax_forces();

  // extract global and per-atom energy from ReaxFF Fortran
  // compute_charge already contributed to eatom

  if (eflag_global) {
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).eb;
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).ea;
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).elp;
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).emol;
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).ev;
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).epen;
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).ecoa;
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).ehb;
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).et;
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).eco;
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).ew;
    evdwl += FORTRAN(cbkenergies, CBKENERGIES).efi;

    ecoul += FORTRAN(cbkenergies, CBKENERGIES).ep;
    ecoul += energy_charge_equilibration;

    eng_vdwl += evdwl;
    eng_coul += ecoul;

    // Store the different parts of the energy
    // in a list for output by compute pair command

    pvector[0] = FORTRAN(cbkenergies, CBKENERGIES).eb;
    pvector[1] = FORTRAN(cbkenergies, CBKENERGIES).ea;
    pvector[2] = FORTRAN(cbkenergies, CBKENERGIES).elp;
    pvector[3] = FORTRAN(cbkenergies, CBKENERGIES).emol;
    pvector[4] = FORTRAN(cbkenergies, CBKENERGIES).ev;
    pvector[5] = FORTRAN(cbkenergies, CBKENERGIES).epen;
    pvector[6] = FORTRAN(cbkenergies, CBKENERGIES).ecoa;
    pvector[7] = FORTRAN(cbkenergies, CBKENERGIES).ehb;
    pvector[8] = FORTRAN(cbkenergies, CBKENERGIES).et;
    pvector[9] = FORTRAN(cbkenergies, CBKENERGIES).eco;
    pvector[10] = FORTRAN(cbkenergies, CBKENERGIES).ew;
    pvector[11] = FORTRAN(cbkenergies, CBKENERGIES).ep;
    pvector[12] = FORTRAN(cbkenergies, CBKENERGIES).efi;
    pvector[13] = energy_charge_equilibration;

  }

  if (eflag_atom) {
    int ntotal = atom->nlocal + atom->nghost;
    for (i = 0; i < ntotal; i++)
      eatom[i] += FORTRAN(cbkd,CBKD).estrain[i];
  }

  // extract global and per-atom virial from ReaxFF Fortran

  if (vflag_global) {
    virial[0] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[0];
    virial[1] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[1];
    virial[2] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[2];
    virial[3] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[3];
    virial[4] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[4];
    virial[5] = -FORTRAN(cbkvirial, CBKVIRIAL).virial[5];
  }

  if (vflag_atom) {
    int ntotal = atom->nlocal + atom->nghost;
    j = 0;
    for (i = 0; i < ntotal; i++) {
      vatom[i][0] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+0];
      vatom[i][1] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+1];
      vatom[i][2] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+2];
      vatom[i][3] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+3];
      vatom[i][4] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+4];
      vatom[i][5] = -FORTRAN(cbkvirial, CBKVIRIAL).atomvirial[j+5];
      j += 6;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairREAX::write_reax_positions()
{
  int j, jx, jy, jz, jia;

  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  FORTRAN(rsmall, RSMALL).na = nlocal+nghost;
  FORTRAN(rsmall, RSMALL).na_local = nlocal;

  if (nlocal+nghost > ReaxParams::nat)
    error->one(FLERR,"Reax_defs.h setting for NATDEF is too small");

  jx = 0;
  jy = ReaxParams::nat;
  jz = 2*ReaxParams::nat;
  jia = 0;

  j = 0;
  for (int i = 0; i < nlocal+nghost; i++, j++) {
    FORTRAN(cbkc, CBKC).c[j+jx] = x[i][0];
    FORTRAN(cbkc, CBKC).c[j+jy] = x[i][1];
    FORTRAN(cbkc, CBKC).c[j+jz] = x[i][2];
    FORTRAN(cbkch, CBKCH).ch[j] = q[i];
    FORTRAN(cbkia, CBKIA).ia[j+jia] = map[type[i]];
    FORTRAN(cbkia, CBKIA).iag[j+jia] = map[type[i]];
    FORTRAN(cbkc, CBKC).itag[j] = tag[i];
  }
}

/* ---------------------------------------------------------------------- */

void PairREAX::write_reax_vlist()
{
  int ii, jj, i, j, iii, jjj;
  double xitmp, yitmp, zitmp;
  double xjtmp, yjtmp, zjtmp;
  int itag,jtag;
  int nvpair, nvlself, nvpairmax;
  int nbond;
  int inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double delr2;
  double delx, dely, delz;

  double **x = atom->x;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  nvpairmax = ReaxParams::nneighmax * ReaxParams::nat;

  nvpair = 0;
  nvlself =0;
  nbond = 0;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xitmp = x[i][0];
    yitmp = x[i][1];
    zitmp = x[i][2];
    itag = tag[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      xjtmp = x[j][0];
      yjtmp = x[j][1];
      zjtmp = x[j][2];
      jtag = tag[j];

      delx = xitmp - xjtmp;
      dely = yitmp - yjtmp;
      delz = zitmp - zjtmp;

      delr2 = delx*delx+dely*dely+delz*delz;

      if (delr2 <= rcutvsq) {
        if (i < j) {
          iii = i+1;
          jjj = j+1;
        } else {
          iii = j+1;
          jjj = i+1;
        }
        if (nvpair >= nvpairmax)
          error->one(FLERR,"Reax_defs.h setting for NNEIGHMAXDEF is too small");

        FORTRAN(cbkpairs, CBKPAIRS).nvl1[nvpair] = iii;
        FORTRAN(cbkpairs, CBKPAIRS).nvl2[nvpair] = jjj;
        FORTRAN(cbknvlbo, CBKNVLBO).nvlbo[nvpair] = 0;

        if (delr2 <= rcutbsq) {
          FORTRAN(cbknvlbo, CBKNVLBO).nvlbo[nvpair] = 1;
          nbond++;
        }

        FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 0;

        if (j < nlocal)
          FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
        else if (itag < jtag)
          FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
        else if (itag == jtag) {
          if (delz > SMALL)
            FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
          else if (fabs(delz) < SMALL) {
            if (dely > SMALL)
              FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
            else if (fabs(dely) < SMALL && delx > SMALL)
              FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 1;
          }
        }
        nvpair++;
      }
    }
  }

  int ntotal = nlocal + nghost;

  for (int i = nlocal; i < ntotal; i++) {
    xitmp = x[i][0];
    yitmp = x[i][1];
    zitmp = x[i][2];
    itag = tag[i];

    for (int j = i+1; j < ntotal; j++) {
      xjtmp = x[j][0];
      yjtmp = x[j][1];
      zjtmp = x[j][2];
      jtag = tag[j];

      delx = xitmp - xjtmp;
      dely = yitmp - yjtmp;
      delz = zitmp - zjtmp;

      delr2 = delx*delx+dely*dely+delz*delz;

      // don't need to check the double count since i < j in the ghost region

      if (delr2 <= rcutvsq) {
        iii = i+1;
        jjj = j+1;

        if (nvpair >= nvpairmax)
          error->one(FLERR,"Reax_defs.h setting for NNEIGHMAXDEF is too small");

        FORTRAN(cbkpairs, CBKPAIRS).nvl1[nvpair] = iii;
        FORTRAN(cbkpairs, CBKPAIRS).nvl2[nvpair] = jjj;
        FORTRAN(cbknvlbo, CBKNVLBO).nvlbo[nvpair] = 0;

        if (delr2 <= rcutbsq) {
          FORTRAN(cbknvlbo, CBKNVLBO).nvlbo[nvpair] = 1;
          nbond++;
        }

        FORTRAN(cbknvlown, CBKNVLOWN).nvlown[nvpair] = 0;
        nvpair++;
      }
    }
  }

  FORTRAN(cbkpairs, CBKPAIRS).nvpair = nvpair;
  FORTRAN(cbkpairs, CBKPAIRS).nvlself = nvlself;
}

/* ---------------------------------------------------------------------- */

void PairREAX::read_reax_forces()
{
  double ftmp[3];

  double **f = atom->f;
  int ntotal = atom->nlocal + atom->nghost;

  int j = 0;
  for (int i = 0; i < ntotal; i++) {
    ftmp[0] = -FORTRAN(cbkd, CBKD).d[j];
    ftmp[1] = -FORTRAN(cbkd, CBKD).d[j+1];
    ftmp[2] = -FORTRAN(cbkd, CBKD).d[j+2];
    f[i][0] = ftmp[0];
    f[i][1] = ftmp[1];
    f[i][2] = ftmp[2];
    j += 3;
  }
}

/* ---------------------------------------------------------------------- */

void PairREAX::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  param_list = new ff_params[n+1];
  for (int i = 1; i <= n; i++)
    param_list[i].params = new double[5];

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairREAX::settings(int narg, char **arg)
{
  if (narg != 0 && narg !=4) error->all(FLERR,"Illegal pair_style command");

  if (narg == 4) {
    hbcut = force->numeric(FLERR,arg[0]);
    ihbnew = static_cast<int> (force->numeric(FLERR,arg[1]));
    itripstaball = static_cast<int> (force->numeric(FLERR,arg[2]));
    precision = force->numeric(FLERR,arg[3]);

    if (hbcut <= 0.0 ||
        (ihbnew != 0 && ihbnew != 1) ||
        (itripstaball != 0 && itripstaball != 1) ||
        precision <= 0.0)
      error->all(FLERR,"Illegal pair_style command");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairREAX::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure filename is ffield.reax

  if (strcmp(arg[2],"ffield.reax") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // NOTE: for now throw an error if NULL is used to disallow use with hybrid
  //       qEq matrix solver needs to be modified to exclude atoms

  for (int i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      error->all(FLERR,"Cannot currently use pair reax with pair hybrid");
      continue;
    }
    map[i-2] = force->inumeric(FLERR,arg[i]);
  }

  int n = atom->ntypes;

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 1;
      count++;
    }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairREAX::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style reax requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style reax requires newton pair on");
  if (!atom->q_flag)
    error->all(FLERR,"Pair style reax requires atom attribute q");
  if (strcmp(update->unit_style,"real") != 0 && comm->me == 0)
    error->warning(FLERR,"Not using real units with pair reax");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->newton = 2;

  FORTRAN(readc, READC)();
  FORTRAN(reaxinit, REAXINIT)();
  FORTRAN(ffinpt, FFINPT)();
  FORTRAN(tap7th, TAP7TH)();

  // turn off read_in by fort.3 in REAX Fortran

  int ngeofor_tmp = -1;
  FORTRAN(setngeofor, SETNGEOFOR)(&ngeofor_tmp);
  if (comm->me == 0) FORTRAN(readgeo, READGEO)();

  // initial setup for cutoff radius of VLIST and BLIST in ReaxFF

  double vlbora;

  FORTRAN(getswb, GETSWB)(&swb);
  cutmax=MAX(swb, hbcut);
  rcutvsq=cutmax*cutmax;
  FORTRAN(getvlbora, GETVLBORA)(&vlbora);
  rcutbsq=vlbora*vlbora;

  // parameters for charge equilibration from ReaxFF input, fort.4
  // verify that no LAMMPS type to REAX type mapping was invalid

  int nelements;
  FORTRAN(getnso, GETNSO)(&nelements);

  FORTRAN(getswa, GETSWA)(&swa);
  double chi, eta, gamma;
  for (int itype = 1; itype <= atom->ntypes; itype++) {
    if (map[itype] < 1 || map[itype] > nelements)
      error->all(FLERR,"Invalid REAX atom type");
    chi = FORTRAN(cbkchb, CBKCHB).chi[map[itype]-1];
    eta = FORTRAN(cbkchb, CBKCHB).eta[map[itype]-1];
    gamma = FORTRAN(cbkchb, CBKCHB).gam[map[itype]-1];
    param_list[itype].np = 5;
    param_list[itype].rcutsq = cutmax;
    param_list[itype].params[0] = chi;
    param_list[itype].params[1] = eta;
    param_list[itype].params[2] = gamma;
    param_list[itype].params[3] = swa;
    param_list[itype].params[4] = swb;
  }

  taper_setup();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairREAX::init_one(int i, int j)
{
  return cutmax;
}

/* ---------------------------------------------------------------------- */

int PairREAX::pack_forward_comm(int n, int *list, double *buf, 
                                int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;

  if (packflag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = FORTRAN(cbkabo, CBKABO).abo[j];
    }

  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = wcg[j];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairREAX::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (packflag == 0) {
    for (i = first; i < last; i++)
      FORTRAN(cbkabo, CBKABO).abo[i] = buf[m++];

  } else {
    for (i = first; i < last; i++)
      wcg[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairREAX::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    buf[m++] = wcg[i];

  return m;
}

/* ---------------------------------------------------------------------- */

void PairREAX::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    wcg[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   charge equilibration routines
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void PairREAX::taper_setup()
{
  double swb2,swa2,swb3,swa3,d1,d7;

  d1=swb-swa;
  d7=pow(d1,7.0);
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

/* ---------------------------------------------------------------------- */

double PairREAX::taper_E(const double &r, const double &r2)
{
  double r3=r2*r;
  return swc7*r3*r3*r+swc6*r3*r3+swc5*r3*r2+swc4*r2*r2+swc3*r3+swc2*r2+
     swc1*r+swc0;
}

/* ---------------------------------------------------------------------- */

double PairREAX::taper_F(const double &r, const double &r2)
{
  double r3=r2*r;
  return 7.0e0*swc7*r3*r3+6.0e0*swc6*r3*r2+5.0e0*swc5*r2*r2+
    4.0e0*swc4*r3+3.0e0*swc3*r2+2.0e0*swc2*r+swc1;
}

/* ----------------------------------------------------------------------
   compute current charge distributions based on the charge equilibration
------------------------------------------------------------------------- */

void PairREAX::compute_charge(double &energy_charge_equilibration)
{
  double xitmp, yitmp, zitmp;
  double xjtmp, yjtmp, zjtmp;
  int itype, jtype, itag, jtag;
  int ii, jj, i, j;
  double delr2, delr_norm, gamt, hulp1, hulp2;
  double delx, dely, delz;
  double qsum,qi;
  int nmatentries;
  double sw;
  int inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int *tag = atom->tag;

  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // realloc neighbor based arrays if necessary

  int numneigh_total = 0;
  for (ii = 0; ii < inum; ii++)
    numneigh_total += numneigh[ilist[ii]];

  if (numneigh_total + 2*nlocal > matmax) {
    memory->destroy(aval);
    memory->destroy(acol_ind);
    matmax = numneigh_total + 2*nlocal;
    memory->create(aval,matmax,"reax:aval");
    memory->create(acol_ind,matmax,"reax:acol_ind");
  }

  // build linear system

  nmatentries = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xitmp = x[i][0];
    yitmp = x[i][1];
    zitmp = x[i][2];
    itype = type[i];
    itag = tag[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    arow_ptr[i] = nmatentries;
    aval[nmatentries] = 2.0*param_list[itype].params[1];
    acol_ind[nmatentries] = i;
    nmatentries++;

    aval[nmatentries] = 1.0;
    acol_ind[nmatentries] = nlocal + nghost;
    nmatentries++;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      xjtmp = x[j][0];
      yjtmp = x[j][1];
      zjtmp = x[j][2];
      jtype = type[j];
      jtag = tag[j];

      delx = xitmp - xjtmp;
      dely = yitmp - yjtmp;
      delz = zitmp - zjtmp;

      delr2 = delx*delx+dely*dely+delz*delz;

      // avoid counting local-ghost pair twice since
      // ReaxFF uses half neigh list with newton off

      if (j >= nlocal) {
        if (itag > jtag) {
          if ((itag+jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag+jtag) % 2 == 1) continue;
        } else {
          if (zjtmp < zitmp) continue;
          if (zjtmp == zitmp && yjtmp < yitmp) continue;
          if (zjtmp == zitmp && yjtmp == yitmp && xjtmp < xitmp) continue;
        }
      }

      // rcutvsq = cutmax*cutmax, in ReaxFF

      if (delr2 <= rcutvsq) {
        gamt = sqrt(param_list[itype].params[2]*param_list[jtype].params[2]);
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

  // in this case, we don't use Midpoint method
  // so, we don't need to consider ghost-ghost interactions
  // but, need to fill the arow_ptr[] arrays for the ghost atoms

  for (i = nlocal; i < nlocal+nghost; i++)
    arow_ptr[i] = nmatentries;
  arow_ptr[nlocal+nghost] = nmatentries;

  // add rhs matentries to linear system

  for (ii =0; ii<inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    elcvec[i] = -param_list[itype].params[0];
  }

  for (i = nlocal; i < nlocal+nghost; i++) elcvec[i] = 0.0;

  // assign current charges to charge vector

  qsum = 0.0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qi = q[i];
    ch[i] = qi;
    if (i < nlocal) qsum += qi;
  }

  for (i = nlocal; i < nlocal+nghost; i++) {
    qi = q[i];
    ch[i] = qi;
  }

  double qtot;
  MPI_Allreduce(&qsum,&qtot,1,MPI_DOUBLE,MPI_SUM,world);
  elcvec[nlocal+nghost] = 0.0;
  ch[nlocal+nghost] = chpot;

  // solve the linear system using CG sover

  charge_reax(nlocal,nghost,ch,aval,acol_ind,arow_ptr,elcvec);

  // calculate the charge equilibration energy

  energy_charge_equilibration = 0;

  // have already updated charge distributions for the current structure

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    // 23.02 is the ReaxFF conversion from eV to kcal/mol
    // should really use constants.evfactor ~23.06
    // but that would break consistency with serial ReaxFF code
    // NOTE: this hard-wired real units
    //       if want other units would have to change params[] in file

    qi = 23.02 * (param_list[itype].params[0]*ch[i]+
                  param_list[itype].params[1]*ch[i]*ch[i]);
    energy_charge_equilibration += qi;
    if (eflag_atom) eatom[i] += qi;
  }

  // copy charge vector back to particles from the calculated values

  for (i = 0; i < nlocal+nghost; i++) q[i] = ch[i];
  chpot = ch[nlocal+nghost];
}

/* ---------------------------------------------------------------------- */

void PairREAX::charge_reax(const int & nlocal, const int & nghost,
                           double ch[], double aval[], int acol_ind[],
                           int arow_ptr[], double elcvec[])
{
  cg_solve(nlocal,nghost,aval,acol_ind,arow_ptr,ch,elcvec);
}

/* ----------------------------------------------------------------------
   CG solver for linear systems
------------------------------------------------------------------------- */

void PairREAX::cg_solve(const int & nlocal, const int & nghost,
                        double aval[], int acol_ind[], int arow_ptr[],
                        double x[], double b[])
{
  double one, zero, rho, rho_old, alpha, beta, gamma;
  int iter, maxiter;
  int n;
  double sumtmp;

  // parallel CG method by A. P. Thompson
  // distributed (partial) vectors: b, r, q, A
  // accumulated (full) vectors: x, w, p
  // r = b-A.x
  // w = r            (ReverseComm + Comm)

  double *r = rcg;
  double *w = wcg;
  double *p = pcg;
  double *p_old = poldcg;
  double *q = qcg;

  n = nlocal+nghost+1;

  one = 1.0;
  zero = 0.0;
  maxiter = 100;

  for (int i = 0; i < n; i++) w[i] = 0;

  // construct r = b-Ax

  sparse_product(n, nlocal, nghost, aval, acol_ind, arow_ptr, x, r);

  // not using BLAS library

  for (int i=0; i<n; i++) {
    r[i] = b[i] - r[i];
    w[i] = r[i];
  }

  packflag = 1;
  comm->reverse_comm_pair(this);
  comm->forward_comm_pair(this);

  MPI_Allreduce(&w[n-1], &sumtmp, 1, MPI_DOUBLE, MPI_SUM, world);
  w[n-1] = sumtmp;
  rho_old = one;

  for (iter = 1; iter < maxiter; iter++) {
    rho = 0.0;
    for (int i=0; i<nlocal; i++) rho += w[i]*w[i];

    MPI_Allreduce(&rho, &sumtmp, 1, MPI_DOUBLE, MPI_SUM, world);
    rho = sumtmp + w[n-1]*w[n-1];
    if (rho < precision) break;

    for (int i = 0; i<n; i++) p[i] = w[i];

    if (iter > 1) {
      beta = rho/rho_old;
      for (int i = 0; i<n; i++) p[i] += beta*p_old[i];
    }

    sparse_product(n, nlocal, nghost, aval, acol_ind, arow_ptr, p, q);

    gamma = 0.0;
    for (int i=0; i<n; i++) gamma += p[i]*q[i];
    MPI_Allreduce(&gamma, &sumtmp, 1, MPI_DOUBLE, MPI_SUM, world);

    gamma = sumtmp;
    alpha = rho/gamma;

    for (int i=0; i<n; i++) {
      x[i] += alpha*p[i];
      r[i] -= alpha*q[i];
      w[i] = r[i];
    }

    comm->reverse_comm_pair(this);
    comm->forward_comm_pair(this);

    MPI_Allreduce(&w[n-1], &sumtmp, 1, MPI_DOUBLE, MPI_SUM, world);
    w[n-1] = sumtmp;

    for (int i=0; i<n; i++) p_old[i] = p[i];
    rho_old = rho;
  }
}

/* ----------------------------------------------------------------------
   sparse maxtrix operations
------------------------------------------------------------------------- */

void PairREAX::sparse_product(const int &n, const int &nlocal,
                              const int &nghost,
                              double aval[], int acol_ind[], int arow_ptr[],
                              double *x, double *r)
{
  int i,j,jj;

  for (i=0; i<n; i++) r[i] = 0.0;

  for (i=0; i<nlocal; i++) {
    r[i] += aval[arow_ptr[i]]*x[i];
    for (j=arow_ptr[i]+1; j<arow_ptr[i+1]; j++) {
      jj = acol_ind[j];
      r[i] += aval[j]*x[jj];
      r[jj] += aval[j]*x[i];
    }
  }

  for (i=nlocal; i<nlocal+nghost; i++)
    for (j=arow_ptr[i]; j<arow_ptr[i+1]; j++) {
      jj = acol_ind[j];
      r[i] += aval[j]*x[jj];
      r[jj] += aval[j]*x[i];
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairREAX::memory_usage()
{
  double bytes = nmax * sizeof(int);
  bytes += 7 * nmax * sizeof(double);
  bytes += matmax * sizeof(int);
  bytes += matmax * sizeof(double);
  return bytes;
}
