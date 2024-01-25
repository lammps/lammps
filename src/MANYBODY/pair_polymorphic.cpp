// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Reese Jones, Xiaowang Zhou (SNL)
   This modifies from pair_tersoff.cpp by Aidan Thompson (SNL)
------------------------------------------------------------------------- */

// uncomment define to enable writing table files for debugging
// #define LMP_POLYMORPHIC_WRITE_TABLES 1

#include "pair_polymorphic.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "tabular_function.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathExtra;

static constexpr int MAXLINE = 1024;
static constexpr int DELTA = 4;


/* ---------------------------------------------------------------------- */

PairPolymorphic::PairParameters::PairParameters()
{
  cut = 0.0;
  cutsq = 0.0;
  xi = 0.0;
  U = nullptr;
  V = nullptr;
  W = nullptr;
  F = nullptr;
}

PairPolymorphic::PairParameters::~PairParameters()
{
  delete U;
  delete V;
  delete W;
  delete F;
}

/* ---------------------------------------------------------------------- */

PairPolymorphic::TripletParameters::TripletParameters()
{
  P = nullptr;
  G = nullptr;
}

PairPolymorphic::TripletParameters::~TripletParameters()
{
  delete P;
  delete G;
}

/* ---------------------------------------------------------------------- */

PairPolymorphic::PairPolymorphic(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;

  match = nullptr;
  pairParameters = nullptr;
  tripletParameters = nullptr;
  elem2param = nullptr;
  elem3param = nullptr;
  epsilon = 0.0;
  neighsize = 0;
  firstneighV = nullptr;
  firstneighW = nullptr;
  firstneighW1 = nullptr;
  delxV = nullptr;
  delyV = nullptr;
  delzV = nullptr;
  drV = nullptr;
  delxW = nullptr;
  delyW = nullptr;
  delzW = nullptr;
  drW = nullptr;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairPolymorphic::~PairPolymorphic()
{
  delete [] match;

  delete [] pairParameters;
  delete [] tripletParameters;

  memory->destroy(elem2param);
  memory->destroy(elem3param);
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] firstneighV;
    delete [] firstneighW;
    delete [] firstneighW1;
    delete [] delxV;
    delete [] delyV;
    delete [] delzV;
    delete [] drV;
    delete [] delxW;
    delete [] delyW;
    delete [] delzW;
    delete [] drW;
  }
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::compute(int eflag, int vflag)
{
  tagint itag,jtag;
  int i,j,k,ii,jj,kk,kk1,inum,jnum;
  int itype,jtype,ktype;
  int iparam_ii,iparam_jj,iparam_kk,iparam_ij,iparam_ik,iparam_ijk;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r0,r1,r2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor,wfac,pfac,gfac,fa,fa_d,bij,bij_d;
  double costheta;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double emb;

  evdwl = 0.0;
  ev_init(eflag,vflag);

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

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (neighsize < jnum) {
      delete [] firstneighV;
      delete [] delxV;
      delete [] delyV;
      delete [] delzV;
      delete [] drV;
      delete [] firstneighW;
      delete [] delxW;
      delete [] delyW;
      delete [] delzW;
      delete [] drW;
      delete [] firstneighW1;
      neighsize = jnum + 20;
      firstneighV = new int[neighsize];
      delxV = new double[neighsize];
      delyV = new double[neighsize];
      delzV = new double[neighsize];
      drV = new double[neighsize];
      firstneighW = new int[neighsize];
      delxW = new double[neighsize];
      delyW = new double[neighsize];
      delzW = new double[neighsize];
      drW = new double[neighsize];
      firstneighW1 = new int[neighsize];
    }

    if (eta == 1) {
      iparam_ii = elem2param[itype][itype];
      PairParameters &p = pairParameters[iparam_ii];
      emb = (p.F)->get_vmax();
    }

    numneighV = -1;
    numneighW = -1;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= cutmaxsq) continue;
      r0 = sqrt(rsq);

      iparam_ij = elem2param[itype][jtype];
      PairParameters &p = pairParameters[iparam_ij];

// do not include the neighbor if get_vmax() <= epsilon because the function is near zero
      if (eta == 1) {
        if (emb > epsilon) {
          iparam_jj = elem2param[jtype][jtype];
          PairParameters &q = pairParameters[iparam_jj];
          if (rsq < (q.W)->get_xmaxsq() && (q.W)->get_vmax() > epsilon) {
            numneighW = numneighW + 1;
            firstneighW[numneighW] = j;
            delxW[numneighW] = delx;
            delyW[numneighW] = dely;
            delzW[numneighW] = delz;
            drW[numneighW] =  r0;
          }
        }
      } else {
        if ((p.F)->get_vmax() > epsilon) {
          if (rsq < (p.V)->get_xmaxsq() && (p.V)->get_vmax() > epsilon) {
            numneighV = numneighV + 1;
            firstneighV[numneighV] = j;
            delxV[numneighV] = delx;
            delyV[numneighV] = dely;
            delzV[numneighV] = delz;
            drV[numneighV] =  r0;
          }
          if (rsq < (p.W)->get_xmaxsq() && (p.W)->get_vmax() > epsilon) {
            numneighW = numneighW + 1;
            firstneighW[numneighW] = j;
            delxW[numneighW] = delx;
            delyW[numneighW] = dely;
            delzW[numneighW] = delz;
            drW[numneighW] =  r0;
          }
        }
      }

    // two-body interactions, skip half of them

      jtag = tag[j];
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < x[i][2]) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      if (rsq >= (p.U)->get_xmaxsq() || (p.U)->get_vmax() <= epsilon) continue;
      (p.U)->value(r0,evdwl,eflag,fpair,1);
      fpair = -fpair/r0;

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
    }

    if (eta == 1) {

      if (emb > epsilon) {

        iparam_ii = elem2param[itype][itype];
        PairParameters &p = pairParameters[iparam_ii];

        // accumulate bondorder zeta for each i-j interaction via loop over k

        zeta_ij = 0.0;

        for (kk = 0; kk <= numneighW; kk++) {
          k = firstneighW[kk];
          ktype = map[type[k]];

          iparam_kk = elem2param[ktype][ktype];
          PairParameters &q = pairParameters[iparam_kk];

          (q.W)->value(drW[kk],wfac,1,fpair,0);

          zeta_ij += wfac;
        }

        // pairwise force due to zeta

        (p.F)->value(zeta_ij,bij,1,bij_d,1);

        prefactor = 0.5* bij_d;
        if (eflag) evdwl = -0.5*bij;

        if (evflag) ev_tally(i,i,nlocal,newton_pair,evdwl,0.0,0.0,delx,dely,delz);

        // attractive term via loop over k

        for (kk = 0; kk <= numneighW; kk++) {
          k = firstneighW[kk];
          ktype = map[type[k]];

          delr2[0] = -delxW[kk];
          delr2[1] = -delyW[kk];
          delr2[2] = -delzW[kk];

          iparam_kk = elem2param[ktype][ktype];
          PairParameters &q = pairParameters[iparam_kk];

          (q.W)->value(drW[kk],wfac,0,fpair,1);
          fpair = -prefactor*fpair/drW[kk];

          f[i][0] += delr2[0]*fpair;
          f[i][1] += delr2[1]*fpair;
          f[i][2] += delr2[2]*fpair;
          f[k][0] -= delr2[0]*fpair;
          f[k][1] -= delr2[1]*fpair;
          f[k][2] -= delr2[2]*fpair;

          if (vflag_either) v_tally2(i, k, -fpair, delr2);
        }
      }

    } else {

      for (jj = 0; jj <= numneighV; jj++) {
        j = firstneighV[jj];
        jtype = map[type[j]];

        iparam_ij = elem2param[itype][jtype];
        PairParameters &p = pairParameters[iparam_ij];

        delr1[0] = -delxV[jj];
        delr1[1] = -delyV[jj];
        delr1[2] = -delzV[jj];
        r1 = drV[jj];

        // accumulate bondorder zeta for each i-j interaction via loop over k

        zeta_ij = 0.0;

        numneighW1 = -1;
        for (kk = 0; kk <= numneighW; kk++) {
          k = firstneighW[kk];
          if (j == k) continue;
          ktype = map[type[k]];
          iparam_ijk = elem3param[jtype][itype][ktype];
          TripletParameters &trip = tripletParameters[iparam_ijk];
          if ((trip.G)->get_vmax() <= epsilon) continue;

          numneighW1 = numneighW1 + 1;
          firstneighW1[numneighW1] = kk;

          delr2[0] = -delxW[kk];
          delr2[1] = -delyW[kk];
          delr2[2] = -delzW[kk];
          r2 = drW[kk];

          costheta = (delr1[0]*delr2[0] + delr1[1]*delr2[1] +
                      delr1[2]*delr2[2]) / (r1*r2);

          iparam_ik = elem2param[itype][ktype];
          PairParameters &q = pairParameters[iparam_ik];

          (q.W)->value(r2,wfac,1,fpair,0);
          (trip.P)->value(r1-(p.xi)*r2,pfac,1,fpair,0);
          (trip.G)->value(costheta,gfac,1,fpair,0);

          zeta_ij += wfac*pfac*gfac;
        }

        // pairwise force due to zeta

        (p.V)->value(r1,fa,1,fa_d,1);
        (p.F)->value(zeta_ij,bij,1,bij_d,1);
        fpair = -0.5*bij*fa_d / r1;
        prefactor = 0.5* fa * bij_d;
        if (eflag) evdwl = -0.5*bij*fa;

        f[i][0] += delr1[0]*fpair;
        f[i][1] += delr1[1]*fpair;
        f[i][2] += delr1[2]*fpair;
        f[j][0] -= delr1[0]*fpair;
        f[j][1] -= delr1[1]*fpair;
        f[j][2] -= delr1[2]*fpair;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,
                             -fpair,-delr1[0],-delr1[1],-delr1[2]);

        // attractive term via loop over k

        for (kk1 = 0; kk1 <= numneighW1; kk1++) {
          kk = firstneighW1[kk1];
          k = firstneighW[kk];
          ktype = map[type[k]];
          iparam_ijk = elem3param[jtype][itype][ktype];
          TripletParameters &trip = tripletParameters[iparam_ijk];

          delr2[0] = -delxW[kk];
          delr2[1] = -delyW[kk];
          delr2[2] = -delzW[kk];
          r2 = drW[kk];

          iparam_ik = elem2param[itype][ktype];
          PairParameters &q = pairParameters[iparam_ik];

          attractive(&p,&q,&trip,prefactor,r1,r2,delr1,delr2,fi,fj,fk);

          f[i][0] += fi[0];
          f[i][1] += fi[1];
          f[i][2] += fi[2];
          f[j][0] += fj[0];
          f[j][1] += fj[1];
          f[j][2] += fj[2];
          f[k][0] += fk[0];
          f[k][1] += fk[1];
          f[k][2] += fk[2];

          if (vflag_either) v_tally3(i,j,k,fj,fk,delr1,delr2);
        }
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];

  neighsize = 40;
  firstneighV = new int[neighsize];
  delxV = new double[neighsize];
  delyV = new double[neighsize];
  delzV = new double[neighsize];
  drV = new double[neighsize];
  firstneighW = new int[neighsize];
  delxW = new double[neighsize];
  delyW = new double[neighsize];
  delzW = new double[neighsize];
  drW = new double[neighsize];
  firstneighW1 = new int[neighsize];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPolymorphic::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPolymorphic::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  // parse and remove optional last parameter

  if (narg == 4 + atom->ntypes)
    epsilon = utils::numeric(FLERR,arg[--narg],false,lmp);

  map_element2type(narg-3,arg+3);

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPolymorphic::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style polymorphic requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style polymorphic requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPolymorphic::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::read_file(char *file)
{
  PotentialFileReader * reader = nullptr;

  // read potential
  if (comm->me == 0) {
    try {
      reader = new PotentialFileReader(lmp, file, "polymorphic");

      ValueTokenizer values = reader->next_values(2);

      int ntypes = values.next_int();

      if (ntypes != nelements)
        error->one(FLERR,"Incorrect number of elements (expected: {} found: {}) in potential file",
                   nelements, ntypes);

      eta = values.next_int();

      // map the elements in the potential file to LAMMPS atom types
      delete [] match;
      match = new int[nelements];

      for (int i = 0; i < nelements; i++) {
        values = reader->next_values(3);
        values.next_double(); // atomic number
        values.next_double(); // atomic mass
        std::string name = values.next_string();

        int j;
        for (j = 0; j < nelements; j++) {
          if (name == elements[j]) break;
        }
        if (j == nelements) error->one(FLERR, "Element {} in potential file not used", name);
        match[i] = j;
      }

      // sizes
      // Note: the format of this line has changed between the
      // 2015-06-06 and 2015-12-09 versions of the pair style.
      try {
        values = reader->next_values(4);
        nr = ng = nx = 0;
        nr = values.next_int();
        ng = values.next_int();
        nx = values.next_int();
        maxX = values.next_double();

        if ((ng == 0) || (nr == 0) || (nx == 0))
          error->one(FLERR,"Error reading potential file header");
      } catch (TokenizerException &) {
        error->one(FLERR,"Potential file incompatible with this pair style version");
      }

      // cutoffs
      npair = nelements*(nelements+1)/2;
      ntriple = nelements*nelements*nelements;
      delete [] pairParameters;
      delete [] tripletParameters;
      pairParameters = new PairParameters[npair];
      tripletParameters = new TripletParameters[ntriple];

      for (int i = 0; i < npair; i++) {
        PairParameters &p = pairParameters[i];
        values = reader->next_values(2);
        p.cut = values.next_double();
        p.cutsq = p.cut*p.cut;
        p.xi = values.next_double();
      }
    } catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }
  }

  MPI_Bcast(&nr, 1, MPI_INT, 0, world);
  MPI_Bcast(&ng, 1, MPI_INT, 0, world);
  MPI_Bcast(&nx, 1, MPI_INT, 0, world);
  MPI_Bcast(&eta, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxX, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&npair, 1, MPI_INT, 0, world);
  MPI_Bcast(&ntriple, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    delete [] match;
    match = new int[nelements];
    delete [] pairParameters;
    delete [] tripletParameters;
    pairParameters = new PairParameters[npair];
    tripletParameters = new TripletParameters[ntriple];
  }

  MPI_Bcast(match, nelements, MPI_INT, 0, world);
  MPI_Bcast(pairParameters, npair*sizeof(PairParameters), MPI_BYTE, 0, world);

  // start reading tabular functions
  auto  singletable = new double[nr];
  for (int i = 0; i < npair; i++) { // U
    PairParameters &p = pairParameters[i];
    if (comm->me == 0) reader->next_dvector(singletable, nr);
    MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
    p.U = new TabularFunction;
    (p.U)->set_values(nr,0.0,p.cut,singletable);
  }
  for (int i = 0; i < npair; i++) { // V
    PairParameters &p = pairParameters[i];
    if (comm->me == 0) reader->next_dvector(singletable, nr);
    MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
    p.V = new TabularFunction;
    (p.V)->set_values(nr,0.0,p.cut,singletable);
  }
  for (int i = 0; i < npair; i++) { // W
    PairParameters &p = pairParameters[i];
    if (comm->me == 0) reader->next_dvector(singletable, nr);
    MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
    p.W = new TabularFunction;
    (p.W)->set_values(nr,0.0,p.cut,singletable);
  }

  cutmax = 0.0;
  for (int i = 0; i < npair; i++) {
    PairParameters &p = pairParameters[i];
    if (p.cut > cutmax) cutmax = p.cut;
  }
  cutmaxsq = cutmax*cutmax;

  if (eta != 3) {
    for (int j = 0; j < nelements; j++) { // P
      if (comm->me == 0) reader->next_dvector(singletable, nr);
      MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
      for (int i = 0; i < nelements; i++) {
        TripletParameters &p = tripletParameters[i*nelements*nelements+j*nelements+j];
        p.P = new TabularFunction;
        (p.P)->set_values(nr,-cutmax,cutmax,singletable);
      }
    }
    for (int j = 0; j < nelements-1; j++) { // P
      for (int k = j+1; k < nelements; k++) {
        if (comm->me == 0) reader->next_dvector(singletable, nr);
        MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
        for (int i = 0; i < nelements; i++) {
          TripletParameters &p = tripletParameters[i*nelements*nelements+j*nelements+k];
          p.P = new TabularFunction;
          (p.P)->set_values(nr,-cutmax,cutmax,singletable);
          TripletParameters &q = tripletParameters[i*nelements*nelements+k*nelements+j];
          q.P = new TabularFunction;
          (q.P)->set_values(nr,-cutmax,cutmax,singletable);
        }
      }
    }
  }
  if (eta == 3) {
    for (int i = 0; i < ntriple; i++) { // P
      TripletParameters &p = tripletParameters[i];
      if (comm->me == 0) reader->next_dvector(singletable, nr);
      MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
      p.P = new TabularFunction;
      (p.P)->set_values(nr,-cutmax,cutmax,singletable);
    }
  }
  delete[] singletable;
  singletable = new double[ng];
  for (int i = 0; i < ntriple; i++) { // G
    TripletParameters &p = tripletParameters[i];
    if (comm->me == 0) reader->next_dvector(singletable, ng);
    MPI_Bcast(singletable,ng,MPI_DOUBLE,0,world);
    p.G = new TabularFunction;
    (p.G)->set_values(ng,-1.0,1.0,singletable);
  }
  delete[] singletable;
  singletable = new double[nx];
  for (int i = 0; i < npair; i++) { // F
    PairParameters &p = pairParameters[i];
    if (comm->me == 0) reader->next_dvector(singletable, nx);
    MPI_Bcast(singletable,nx,MPI_DOUBLE,0,world);
    p.F = new TabularFunction;
    (p.F)->set_values(nx,0.0,maxX,singletable);
  }
  delete[] singletable;
  if (comm->me == 0) delete reader;

  // recalculate cutoffs of all params
  for (int i = 0; i < npair; i++) {
    PairParameters &p = pairParameters[i];
    p.cut = (p.U)->get_xmax();
    if (p.cut < (p.V)->get_xmax()) p.cut = (p.V)->get_xmax();
    if (p.cut < (p.W)->get_xmax()) p.cut = (p.W)->get_xmax();
    p.cutsq = p.cut*p.cut;
  }

  // set cutmax to max of all params
  cutmax = 0.0;
  for (int i = 0; i < npair; i++) {
    PairParameters &p = pairParameters[i];
    if (cutmax < p.cut) cutmax = p.cut;
  }
  cutmaxsq = cutmax*cutmax;
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::setup_params()
{
  int i,j,k,n;

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,"pair:elem2param");
  memory->destroy(elem3param);
  memory->create(elem3param,nelements,nelements,nelements,"pair:elem3param");

  // map atom pair to parameter index

  n = 0;
  for (i = 0; i < nelements; i++) {
    elem2param[match[i]][match[i]] = n;
    n++;
  }
  for (i = 0; i < nelements-1; i++) {
    for (j = i+1; j < nelements; j++) {
      elem2param[match[i]][match[j]] = n;
      elem2param[match[j]][match[i]] = n;
      n++;
    }
  }

  // map atom triplet to parameter index

  n = 0;
  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        elem3param[match[i]][match[j]][match[k]] = n;
        n++;
      }

// for debugging, call write_tables() to check the tabular functions
#if defined(LMP_POLYMORPHIC_WRITE_TABLES)
  if (comm->me == 0) write_tables(51);
  error->all(FLERR,"Test potential tables");
#endif
}

/* ----------------------------------------------------------------------
   attractive term
------------------------------------------------------------------------- */

void PairPolymorphic::attractive(PairParameters *p, PairParameters *q,
                                 TripletParameters *trip,
                                 double prefactor, double rij, double rik,
                                 double *delrij, double *delrik,
                                 double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rijinv,rikinv;

  rijinv = 1.0/rij;
  scale3(rijinv,delrij,rij_hat);

  rikinv = 1.0/rik;
  scale3(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,p,q,trip);
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::ters_zetaterm_d(double prefactor,
                                      double *rij_hat, double rij,
                                      double *rik_hat, double rik,
                                      double *dri, double *drj, double *drk,
                                      PairParameters *p, PairParameters *q,
                                      TripletParameters *trip)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  cos_theta = dot3(rij_hat,rik_hat);

  (q->W)->value(rik,fc,1,dfc,1);
  (trip->P)->value(rij-(p->xi)*rik,ex_delr,1,ex_delr_d,1);
  (trip->G)->value(cos_theta,gijk,1,gijk_d,1);

  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri

  scale3(-dfc*gijk*ex_delr,rik_hat,dri);
  scaleadd3(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  scaleadd3(fc*gijk*ex_delr_d*(p->xi),rik_hat,dri,dri);
  scaleadd3(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  scale3(prefactor,dri);

  // compute the derivative wrt Rj

  scale3(fc*gijk_d*ex_delr,dcosdrj,drj);
  scaleadd3(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  scale3(prefactor,drj);

  // compute the derivative wrt Rk

  scale3(dfc*gijk*ex_delr,rik_hat,drk);
  scaleadd3(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  scaleadd3(-fc*gijk*ex_delr_d*(p->xi),rik_hat,drk,drk);
  scale3(prefactor,drk);
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::costheta_d(double *rij_hat, double rij,
                                 double *rik_hat, double rik,
                                 double *dri, double *drj, double *drk)
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = dot3(rij_hat,rik_hat);

  scaleadd3(-cos_theta,rij_hat,rik_hat,drj);
  scale3(1.0/rij,drj);
  scaleadd3(-cos_theta,rik_hat,rij_hat,drk);
  scale3(1.0/rik,drk);
  add3(drj,drk,dri);
  scale3(-1.0,dri);
}

/* ---------------------------------------------------------------------- */
#if defined(LMP_POLYMORPHIC_WRITE_TABLES)
void PairPolymorphic::write_tables(int npts)
{
  FILE* fp =  nullptr;
  double  xmin,xmax,x,uf,vf,wf,pf,gf,ff,ufp,vfp,wfp,pfp,gfp,ffp;
  std::string filename;
  for (int i = 0; i < nelements; i++) {
    for (int j = 0; j < nelements; j++) {
      filename = fmt::format("{}{}_UVW{}",elements[i],
                             elements[j],comm->me);
      fp = fopen(filename.c_str(), "w");
      int iparam_ij = elem2param[i][j];
      PairParameters &pair = pairParameters[iparam_ij];
      xmin = (pair.U)->get_xmin();
      xmax = (pair.U)->get_xmax();
      double xl = xmax - xmin;
      xmin = xmin - 0.5*xl;
      xmax = xmax + 0.5*xl;
      for (int k = 0; k < npts; k++) {
        x = xmin + (xmax-xmin) * k / (npts-1);
        (pair.U)->value(x, uf, 1, ufp, 1);
        (pair.V)->value(x, vf, 1, vfp, 1);
        (pair.W)->value(x, wf, 1, wfp, 1);
        fprintf(fp,"%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f \n",
                x,uf,vf,wf,ufp,vfp,wfp);
      }
      fclose(fp);
    }
  }
  for (int i = 0; i < nelements; i++) {
    for (int j = 0; j < nelements; j++) {
      for (int k = 0; k < nelements; k++) {
        filename = fmt::format("{}{}{}_P{}",elements[i],elements[j],
                               elements[k],comm->me);
        fp = fopen(filename.c_str(), "w");
        int iparam_ij = elem3param[i][j][k];
        TripletParameters &pair = tripletParameters[iparam_ij];
        xmin = (pair.P)->get_xmin();
        xmax = (pair.P)->get_xmax();
        double xl = xmax - xmin;
        xmin = xmin - 0.5*xl;
        xmax = xmax + 0.5*xl;
        for (int n = 0; n < npts; n++) {
          x = xmin + (xmax-xmin) * n / (npts-1);
          (pair.P)->value(x, pf, 1, pfp, 1);
          fprintf(fp,"%12.4f %12.4f %12.4f \n",x,pf,pfp);
        }
        fclose(fp);
      }
    }
  }
  for (int i = 0; i < nelements; i++) {
    for (int j = 0; j < nelements; j++) {
      for (int k = 0; k < nelements; k++) {
        filename = fmt::format("{}{}{}_G{}",elements[i],elements[j],
                               elements[k],comm->me);
        fp = fopen(filename.c_str(), "w");
        int iparam_ij = elem3param[i][j][k];
        TripletParameters &pair = tripletParameters[iparam_ij];
        xmin = (pair.G)->get_xmin();
        xmax = (pair.G)->get_xmax();
        for (int n = 0; n < npts; n++) {
          x = xmin + (xmax-xmin) * n / (npts-1);
          (pair.G)->value(x, gf, 1, gfp, 1);
          fprintf(fp,"%12.4f %12.4f %12.4f \n",x,gf,gfp);
        }
        fclose(fp);
      }
    }
  }
  for (int i = 0; i < nelements; i++) {
    for (int j = 0; j < nelements; j++) {
      filename = fmt::format("{}{}_F{}",elements[i],
                             elements[j],comm->me);
      fp = fopen(filename.c_str(), "w");
      int iparam_ij = elem2param[i][j];
      PairParameters &pair = pairParameters[iparam_ij];
      xmin = (pair.F)->get_xmin();
      xmax = (pair.F)->get_xmax();
      double xl = xmax - xmin;
      xmin = xmin - 0.5*xl;
      xmax = xmax + 0.5*xl;
      for (int k = 0; k < npts; k++) {
        x = xmin + (xmax-xmin) * k / (npts-1);
        (pair.F)->value(x, ff, 1, ffp, 1);
        fprintf(fp,"%12.4f %12.4f %12.4f \n",x,ff,ffp);
      }
      fclose(fp);
    }
  }
}

#endif
