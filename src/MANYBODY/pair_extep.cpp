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
   Contributing author: Jan Los
------------------------------------------------------------------------- */

#include "pair_extep.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathExtra;

static constexpr int DELTA = 4;
static constexpr int PGDELTA = 1;

/* ---------------------------------------------------------------------- */

PairExTeP::PairExTeP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  ghostneigh = 1;

  params = nullptr;

  maxlocal = 0;
  SR_numneigh = nullptr;
  SR_firstneigh = nullptr;
  ipage = nullptr;
  pgsize = oneatom = 0;

  Nt = nullptr;
  Nd = nullptr;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairExTeP::~PairExTeP()
{
  memory->destroy(params);
  memory->destroy(elem3param);

  memory->destroy(SR_numneigh);
  memory->sfree(SR_firstneigh);
  delete[] ipage;
  memory->destroy(Nt);
  memory->destroy(Nd);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);
  }
}

/* ----------------------------------------------------------------------
   create SR neighbor list from main neighbor list
   SR neighbor list stores neighbors of ghost atoms
------------------------------------------------------------------------- */

void PairExTeP::SR_neigh()
{
  int i,j,ii,jj,n,allnum,jnum,itype,jtype,iparam_ij;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;

  if (atom->nmax > maxlocal) {  // ensure there is enough space
    maxlocal = atom->nmax;      // for atoms and ghosts allocated
    memory->destroy(SR_numneigh);
    memory->sfree(SR_firstneigh);
    memory->destroy(Nt);
    memory->destroy(Nd);
    memory->create(SR_numneigh,maxlocal,"ExTeP:numneigh");
    SR_firstneigh = (int **) memory->smalloc(maxlocal*sizeof(int *),
                           "ExTeP:firstneigh");
    memory->create(Nt,maxlocal,"ExTeP:Nt");
    memory->create(Nd,maxlocal,"ExTeP:Nd");
  }

  allnum = list->inum + list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // store all SR neighs of owned and ghost atoms
  // scan full neighbor list of I

  ipage->reset();

  for (ii = 0; ii < allnum; ii++) {
    i = ilist[ii];
    itype=map[type[i]];

    n = 0;
    neighptr = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    Nt[i] = 0.0;
    Nd[i] = 0.0;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      jtype=map[type[j]];
      iparam_ij = elem3param[itype][jtype][jtype];

      if (rsq < params[iparam_ij].cutsq) {
        neighptr[n++] = j;
        double tmp_fc = ters_fc(sqrt(rsq),&params[iparam_ij]);
        Nt[i] += tmp_fc;
        if (itype!=jtype) {
          Nd[i] += tmp_fc;
        }
      }
    }
  //printf("SR_neigh : N[%d] = %f\n",i,N[i]);

    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
}


/* ---------------------------------------------------------------------- */

void PairExTeP::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2,r2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  SR_neigh();

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

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
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

      jtype = map[type[j]];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      iparam_ij = elem3param[itype][jtype][jtype];
      if (rsq > params[iparam_ij].cutsq) continue;

      repulsive(&params[iparam_ij],rsq,fpair,eflag,evdwl);

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }

    // three-body interactions      -(bij + Fcorrection) * fA
    // skip immediately if I-J is not within cutoff

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtag = tag[j];
      jtype = map[type[j]];
      iparam_ij = elem3param[itype][jtype][jtype];

      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 > params[iparam_ij].cutsq) continue;

      // accumulate bondorder zeta for each i-j interaction via loop over k

      zeta_ij = 0.0;

      /* F_IJ (1) */
      // compute correction to energy and forces
      // dE/dr = -Fij(Zi,Zj) dV/dr
      //         - dFij/dZi dZi/dr V
      //         (conjugate term is computed when j is a central atom)

      double FXY, dFXY_dNdij, dFXY_dNdji, fa, fa_d, deng, fpair;
      double Ntij = Nt[i];
      double Ndij = Nd[i];
      double Ntji = Nt[j];
      double Ndji = Nd[j];
      double r = sqrt(rsq1);
      double fc_ij = ters_fc(r,&params[iparam_ij]);

      Ntij -= fc_ij;
      Ntji -= fc_ij;
      if (jtype!=itype) {
        Ndij -= fc_ij;
        Ndji -= fc_ij;
      }
      if (Ntij<0) { Ntij=0.; }
      if (Ndij<0) { Ndij=0.; }
      if (Ntji<0) { Ntji=0.; }
      if (Ndji<0) { Ndji=0.; }
      FXY = F_corr(itype, jtype, Ndij, Ndji, &dFXY_dNdij, &dFXY_dNdji);

      // envelop functions
      double fenv, dfenv_ij;
      fenv = envelop_function(Ntij, Ntji, &dfenv_ij);
      //
      double Fc = fenv * FXY;
      double dFc_dNtij = dfenv_ij * FXY;
      double dFc_dNdij = fenv * dFXY_dNdij;

      fa = ters_fa(r,&params[iparam_ij]);
      fa_d = ters_fa_d(r,&params[iparam_ij]);
      deng = 0.5 * fa * Fc;
      fpair = 0.5 * fa_d * Fc / r;

      f[i][0] += delr1[0]*fpair;
      f[i][1] += delr1[1]*fpair;
      f[i][2] += delr1[2]*fpair;
      f[j][0] -= delr1[0]*fpair;
      f[j][1] -= delr1[1]*fpair;
      f[j][2] -= delr1[2]*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           deng,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);
      /* END F_IJ (1) */

      for (kk = 0; kk < jnum; kk++) {
        if (jj == kk) continue;
        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        iparam_ijk = elem3param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 > params[iparam_ijk].cutsq) continue;

        r2 = sqrt(rsq2);

        zeta_ij += zeta(&params[iparam_ijk],r,r2,delr1,delr2);

        /* F_IJ (2) */
        // compute force components due to spline derivatives
        // uses only the part with FXY_x (FXY_y is done when i and j are inversed)
        int iparam_ik = elem3param[itype][ktype][0];
        double fc_ik_d = ters_fc_d(r2,&params[iparam_ik]);
        double fc_prefac_ik_0 = 1.0 * fc_ik_d * fa / r2;
        double fc_prefac_ik = dFc_dNtij * fc_prefac_ik_0;
        f[i][0] += fc_prefac_ik * delr2[0];
        f[i][1] += fc_prefac_ik * delr2[1];
        f[i][2] += fc_prefac_ik * delr2[2];
        f[k][0] -= fc_prefac_ik * delr2[0];
        f[k][1] -= fc_prefac_ik * delr2[1];
        f[k][2] -= fc_prefac_ik * delr2[2];
        if (vflag_either) v_tally2(i,k,-fc_prefac_ik,delr2);
        if (itype != ktype) {
          fc_prefac_ik = dFc_dNdij * fc_prefac_ik_0;
          f[i][0] += fc_prefac_ik * delr2[0];
          f[i][1] += fc_prefac_ik * delr2[1];
          f[i][2] += fc_prefac_ik * delr2[2];
          f[k][0] -= fc_prefac_ik * delr2[0];
          f[k][1] -= fc_prefac_ik * delr2[1];
          f[k][2] -= fc_prefac_ik * delr2[2];
          if (vflag_either) v_tally2(i,k,-fc_prefac_ik,delr2);
        }
        /* END F_IJ (2) */

      }

      // pairwise force due to zeta

      force_zeta(&params[iparam_ij],r,zeta_ij,fpair,prefactor,eflag,evdwl);

      f[i][0] += delr1[0]*fpair;
      f[i][1] += delr1[1]*fpair;
      f[i][2] += delr1[2]*fpair;
      f[j][0] -= delr1[0]*fpair;
      f[j][1] -= delr1[1]*fpair;
      f[j][2] -= delr1[2]*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);

      // attractive term via loop over k

      for (kk = 0; kk < jnum; kk++) {
        if (jj == kk) continue;
        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        iparam_ijk = elem3param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 > params[iparam_ijk].cutsq) continue;

        attractive(&params[iparam_ijk],prefactor,
                   rsq1,rsq2,delr1,delr2,fi,fj,fk);


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

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairExTeP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutghost,n+1,n+1,"pair:cutghost");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairExTeP::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairExTeP::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg-3,arg+3);

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  spline_init();
  setup();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairExTeP::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style ExTeP requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style ExTeP requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);

  // create pages if first time or if neighbor pgsize/oneatom has changed

  int create = 0;
  if (ipage == nullptr) create = 1;
  if (pgsize != neighbor->pgsize) create = 1;
  if (oneatom != neighbor->oneatom) create = 1;

  if (create) {
    delete[] ipage;
    pgsize = neighbor->pgsize;
    oneatom = neighbor->oneatom;

    int nmypage= comm->nthreads;
    ipage = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      ipage[i].init(oneatom,pgsize,PGDELTA);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairExTeP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  cutghost[i][j] = cutmax ;
  cutghost[j][i] = cutghost[i][j];

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairExTeP::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "ExTeP");
    char *line;

    while ((line = reader.next_line(17))) {
      try {
        ValueTokenizer values(line);

        std::string iname = values.next_string();
        std::string jname = values.next_string();
        std::string kname = values.next_string();

        // ielement,jelement,kelement = 1st args
        // if all 3 args are in element list, then parse this line
        // else skip to next entry in file
        int ielement, jelement, kelement;

        for (ielement = 0; ielement < nelements; ielement++)
          if (iname == elements[ielement]) break;
        if (ielement == nelements) continue;
        for (jelement = 0; jelement < nelements; jelement++)
          if (jname == elements[jelement]) break;
        if (jelement == nelements) continue;
        for (kelement = 0; kelement < nelements; kelement++)
          if (kname == elements[kelement]) break;
        if (kelement == nelements) continue;

        // load up parameter settings and error check their values

        if (nparams == maxparam) {
          maxparam += DELTA;
          params = (Param *) memory->srealloc(params,maxparam*sizeof(Param), "pair:params");

          // make certain all addional allocated storage is initialized
          // to avoid false positives when checking with valgrind

          memset(params + nparams, 0, DELTA*sizeof(Param));
        }

        params[nparams].ielement  = ielement;
        params[nparams].jelement  = jelement;
        params[nparams].kelement  = kelement;
        params[nparams].powerm    = values.next_double();
        params[nparams].gamma     = values.next_double();
        params[nparams].lam3      = values.next_double();
        params[nparams].c         = values.next_double();
        params[nparams].d         = values.next_double();
        params[nparams].h         = values.next_double();
        params[nparams].powern    = values.next_double();
        params[nparams].beta      = values.next_double();
        params[nparams].lam2      = values.next_double();
        params[nparams].bigb      = values.next_double();
        params[nparams].bigr      = values.next_double();
        params[nparams].bigd      = values.next_double();
        params[nparams].lam1      = values.next_double();
        params[nparams].biga      = values.next_double();

        // currently only allow m exponent of 1 or 3

        params[nparams].powermint = int(params[nparams].powerm);

      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }

      if (params[nparams].c < 0.0 ||
          params[nparams].d < 0.0 ||
          params[nparams].powern < 0.0 ||
          params[nparams].beta < 0.0 ||
          params[nparams].lam2 < 0.0 ||
          params[nparams].bigb < 0.0 ||
          params[nparams].bigr < 0.0 ||
          params[nparams].bigd < 0.0 ||
          params[nparams].bigd > params[nparams].bigr ||
          params[nparams].lam1 < 0.0 ||
          params[nparams].biga < 0.0 ||
          params[nparams].powerm - params[nparams].powermint != 0.0 ||
          (params[nparams].powermint != 3 &&
           params[nparams].powermint != 1) ||
          params[nparams].gamma < 0.0)
        error->one(FLERR,"Illegal ExTeP parameter");

      nparams++;
      if (nparams >= pow((double)nelements,3)) break;
    }

    /* F_IJ (3) */
    // initialize F_corr_data to all zeros
    for (int iel=0; iel < nelements; iel++)
      for (int jel=0; jel < nelements; jel++)
        for (int in=0; in < 4; in++)
          for (int jn=0; jn < 4; jn++)
            for (int ivar=0; ivar < 3; ivar++)
              F_corr_data[iel][jel][in][jn][ivar]=0;

    // read the spline coefficients

    while ((line = reader.next_line(8))) {
      try {
        ValueTokenizer values(line);

        std::string iname = values.next_string();
        std::string jname = values.next_string();
        std::string kname = values.next_string();

        // ielement,jelement = 1st args
        // if all 2 args are in element list, then parse this line
        // else skip to next entry in file
        int ielement, jelement;

        for (ielement = 0; ielement < nelements; ielement++)
          if (iname == elements[ielement]) break;
        if (ielement == nelements) continue;
        for (jelement = 0; jelement < nelements; jelement++)
          if (jname == elements[jelement]) break;
        if (jelement == nelements) continue;

        // skip line if it is a leftover from the previous section,
        // which can be identified by having 3 elements (instead of 2)
        // as first words.

        if (!utils::is_integer(kname))
          continue;

        int Ni  = atoi(kname.c_str());
        int Nj  = values.next_int();
        double spline_val = values.next_double();
        double spline_derx = values.next_double();
        double spline_dery = values.next_double();

        // Set value for all pairs of ielement,jelement  (any kelement)
        for (int iparam = 0; iparam < nparams; iparam++) {
          if ( ielement == params[iparam].ielement && jelement == params[iparam].jelement) {
            F_corr_data[ielement][jelement][Ni][Nj][0] = spline_val;
            F_corr_data[ielement][jelement][Ni][Nj][1] = spline_derx;
            F_corr_data[ielement][jelement][Ni][Nj][2] = spline_dery;

            F_corr_data[jelement][ielement][Nj][Ni][0] = spline_val;
            F_corr_data[jelement][ielement][Nj][Ni][1] = spline_dery;
            F_corr_data[jelement][ielement][Nj][Ni][2] = spline_derx;
          }
        }
      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }
    }
  }

  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param), "pair:params");
  }
  MPI_Bcast(params, maxparam*sizeof(Param), MPI_BYTE, 0, world);
  MPI_Bcast(&F_corr_data[0][0][0][0][0], MAXTYPES*MAXTYPES*NSPLINE*NSPLINE*3, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairExTeP::setup()
{
  int i,j,k,m,n;

  // set elem3param for all element triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem3param);
  memory->create(elem3param,nelements,nelements,nelements,"pair:elem3param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has a duplicate entry for: {} {} {}",
                                   elements[i], elements[j], elements[k]);
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry for: {} {} {}",
                              elements[i], elements[j], elements[k]);
        elem3param[i][j][k] = n;
      }

  // compute parameter values derived from inputs

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].bigr + params[m].bigd;
    params[m].cutsq = params[m].cut*params[m].cut;

    params[m].c1 = pow(2.0*params[m].powern*1.0e-16,-1.0/params[m].powern);
    params[m].c2 = pow(2.0*params[m].powern*1.0e-8,-1.0/params[m].powern);
    params[m].c3 = 1.0/params[m].c2;
    params[m].c4 = 1.0/params[m].c1;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].cut > cutmax) cutmax = params[m].cut;
}

/* ---------------------------------------------------------------------- */

void PairExTeP::repulsive(Param *param, double rsq, double &fforce,
                            int eflag, double &eng)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp;

  r = sqrt(rsq);
  tmp_fc = ters_fc(r,param);
  tmp_fc_d = ters_fc_d(r,param);
  tmp_exp = exp(-param->lam1 * r);
  fforce = -param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1) / r;
  if (eflag) eng = tmp_fc * param->biga * tmp_exp;
}

/* ---------------------------------------------------------------------- */

double PairExTeP::zeta(Param *param, double rij, double rik,
                         double *delrij, double *delrik)
{
  double costheta,arg,ex_delr;

  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) / (rij*rik);

  if (param->powermint == 3) arg = pow(param->lam3 * (rij-rik),3.0);
  else arg = param->lam3 * (rij-rik);

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return ters_fc(rik,param) * ters_gijk(costheta,param) * ex_delr;
}

/* ---------------------------------------------------------------------- */

void PairExTeP::force_zeta(Param *param, double r, double zeta_ij,
                             double &fforce, double &prefactor,
                             int eflag, double &eng)
{
  double fa,fa_d,bij;

  fa = ters_fa(r,param);
  fa_d = ters_fa_d(r,param);
  bij = ters_bij(zeta_ij,param);
  fforce = 0.5*bij*fa_d / r;
  prefactor = -0.5*fa * ( ters_bij_d(zeta_ij,param) );
  if (eflag) eng = 0.5*bij*fa;
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

void PairExTeP::attractive(Param *param, double prefactor,
                             double rsqij, double rsqik,
                             double *delrij, double *delrik,
                             double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  scale3(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  scale3(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_fc(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - sin(MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_fc_d(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(MY_PI4/ters_D) * cos(MY_PI2*(r - ters_R)/ters_D);
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_fa(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp(-param->lam2 * r) * ters_fc(r,param);
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_fa_d(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp(-param->lam2 * r) *
    (param->lam2 * ters_fc(r,param) - ters_fc_d(r,param));
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_bij(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return 1.0/sqrt(tmp);
  if (tmp > param->c2)
    return (1.0 - pow(tmp,-param->powern) / (2.0*param->powern))/sqrt(tmp);
  if (tmp < param->c4) return 1.0;
  if (tmp < param->c3)
    return 1.0 - pow(tmp,param->powern)/(2.0*param->powern);
  return pow(1.0 + pow(tmp,param->powern), -1.0/(2.0*param->powern));
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_bij_d(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return param->beta * -0.5*pow(tmp,-1.5);
  if (tmp > param->c2)
    return param->beta * (-0.5*pow(tmp,-1.5) *
                          (1.0 - 0.5*(1.0 +  1.0/(2.0*param->powern)) *
                           pow(tmp,-param->powern)));
  if (tmp < param->c4) return 0.0;
  if (tmp < param->c3)
    return -0.5*param->beta * pow(tmp,param->powern-1.0);

  double tmp_n = pow(tmp,param->powern);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*param->powern)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

void PairExTeP::ters_zetaterm_d(double prefactor,
                                double *rij_hat, double rij,
                                double *rik_hat, double rik,
                                double *dri, double *drj, double *drk,
                                Param *param)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param);
  dfc = ters_fc_d(rik,param);
  if (param->powermint == 3) tmp = pow(param->lam3 * (rij-rik),3.0);
  else tmp = param->lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (param->powermint == 3)
    ex_delr_d = 3.0*pow(param->lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else ex_delr_d = param->lam3 * ex_delr;

  cos_theta = dot3(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta,param);
  gijk_d = ters_gijk_d(cos_theta,param);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  scale3(-dfc*gijk*ex_delr,rik_hat,dri);
  scaleadd3(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  scaleadd3(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  scaleadd3(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  scale3(prefactor,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  scale3(fc*gijk_d*ex_delr,dcosdrj,drj);
  scaleadd3(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  scale3(prefactor,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  scale3(dfc*gijk*ex_delr,rik_hat,drk);
  scaleadd3(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  scaleadd3(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  scale3(prefactor,drk);
}

/* ---------------------------------------------------------------------- */

void PairExTeP::costheta_d(double *rij_hat, double rij,
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

/* F_IJ (4) */
// initialize spline for F_corr (based on PairLCBOP::F_conj)

void PairExTeP::spline_init() {
  for ( int iel=0; iel<nelements; iel++) {
    for ( int jel=0; jel<nelements; jel++) {
      for (int N_ij=0; N_ij<4; N_ij++) {
        for (int N_ji=0; N_ji<4; N_ji++) {
          TF_corr_param &f = F_corr_param[iel][jel][N_ij][N_ji];

          // corner points for each spline function
          f.f_00 = F_corr_data[iel][jel][N_ij  ][N_ji  ][0];
          f.f_01 = F_corr_data[iel][jel][N_ij  ][N_ji+1][0];
          f.f_10 = F_corr_data[iel][jel][N_ij+1][N_ji  ][0];
          f.f_11 = F_corr_data[iel][jel][N_ij+1][N_ji+1][0];

          // compute f tilde according to [Los & Fasolino, PRB 68, 024107 2003]
          f.f_x_00 =   F_corr_data[iel][jel][N_ij  ][N_ji  ][1] - f.f_10 + f.f_00;
          f.f_x_01 =   F_corr_data[iel][jel][N_ij  ][N_ji+1][1] - f.f_11 + f.f_01;
          f.f_x_10 = -(F_corr_data[iel][jel][N_ij+1][N_ji  ][1] - f.f_10 + f.f_00);
          f.f_x_11 = -(F_corr_data[iel][jel][N_ij+1][N_ji+1][1] - f.f_11 + f.f_01);

          f.f_y_00 =   F_corr_data[iel][jel][N_ij  ][N_ji  ][2] - f.f_01 + f.f_00;
          f.f_y_01 = -(F_corr_data[iel][jel][N_ij  ][N_ji+1][2] - f.f_01 + f.f_00);
          f.f_y_10 =   F_corr_data[iel][jel][N_ij+1][N_ji  ][2] - f.f_11 + f.f_10;
          f.f_y_11 = -(F_corr_data[iel][jel][N_ij+1][N_ji+1][2] - f.f_11 + f.f_10);
        }
      }
    }
  }
}

double PairExTeP::envelop_function(double x, double y, double *func_der) {
  double fx,fy,fxy,dfx,dfxydx;
  double del, delsq;

  fxy = 1.0;
  dfxydx = 0.0;

  if (x <= 3.0) {
    fx = 1.0;
    dfx = 0.0;
    if (x < 1.0 && y < 1.0) {
      double gx=(1.0-x);
      double gy=(1.0-y);
      double gxsq=gx*gx;
      double gysq=gy*gy;
      fxy = 1.0 - gxsq*gysq;
      dfxydx = 2.0*gx*gysq;
    }
  } else if (x < 4.0) {
    del = 4.0-x;
    delsq = del*del;
    fx = (3.0-2.0*del)*delsq;
    dfx = - 6.0*del*(1.0-del);
  } else {
    fx = 0.0;
    dfx = 0.0;
  }
  if (y <= 3.0) {
    fy = 1.0;
  } else if (y < 4.0) {
    del = 4.0-y;
    delsq = del*del;
    fy = (3.0-2.0*del)*delsq;
  } else {
    fy = 0.0;
  }

  double func_val = fxy*fx*fy;
  *func_der = dfxydx*fx*fy+fxy*dfx*fy;

  return func_val;
}

double PairExTeP::F_corr(int iel, int jel, double Ndij, double Ndji, double *dFN_x, double *dFN_y) {

  // compute F_XY

  int Ndij_int         = static_cast<int>( floor( Ndij ) );
  int Ndji_int         = static_cast<int>( floor( Ndji ) );
  double x                = Ndij - Ndij_int;
  double y                = Ndji - Ndji_int;
  TF_corr_param &f  = F_corr_param[iel][jel][Ndij_int][Ndji_int];
  double F   = 0;
  double dF_dx = 0, dF_dy = 0;
  double l, r;
  if (Ndij_int < 4 && Ndji_int < 4) {
    l = (1-y)* (1-x);
    r = ( f.f_00 + x*x* f.f_x_10 + y*y* f.f_y_01 );
    F += l*r;
    dF_dx += -(1-y)*r +l*2*x* f.f_x_10;
    dF_dy += -(1-x)*r +l*2*y* f.f_y_01;
    l = (1-y)*x;
    r = ( f.f_10 + (1-x)*(1-x)*f.f_x_00 + y*  y* f.f_y_11 );
    F += l*r;
    dF_dx += (1-y)*r -l*2*(1-x)*f.f_x_00;
    dF_dy += -x*r +l*2*y* f.f_y_11;
    l = y*  (1-x);
    r = ( f.f_01 + x*x* f.f_x_11 + (1-y)*(1-y)*f.f_y_00 );
    F += l*r;
    dF_dx += -y*r +l*2*x* f.f_x_11;
    dF_dy += (1-x)*r -l*2*(1-y)*f.f_y_00;
    l = y*  x;
    r = ( f.f_11 + (1-x)*(1-x)*f.f_x_01 + (1-y)*(1-y)*f.f_y_10 );
    F += l*r;
    dF_dx += y*r -l*2*(1-x)*f.f_x_01;
    dF_dy += x*r -l*2*(1-y)*f.f_y_10;
  }
  double result = F;
  *dFN_x = dF_dx;
  *dFN_y = dF_dy;

  return result;
}
/* F_IJ (4) */
