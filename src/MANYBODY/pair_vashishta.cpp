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
   Contributing author:  Yongnan Xiong (HNU), xyn@hnu.edu.cn
                         Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "pair_vashishta.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define DELTA 4

/* ---------------------------------------------------------------------- */

PairVashishta::PairVashishta(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);

  params = nullptr;

  r0max = 0.0;
  maxshort = 10;
  neighshort = nullptr;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairVashishta::~PairVashishta()
{
  if (copymode) return;

  memory->destroy(params);
  memory->destroy(elem3param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(neighshort);
  }
}

/* ---------------------------------------------------------------------- */

void PairVashishta::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,jnumm1;
  int itype,jtype,ktype,ijparam,ikparam,ijkparam;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fj[3],fk[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  const double cutshortsq = r0max*r0max;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutshortsq) {
        neighshort[numshort++] = j;
        if (numshort >= maxshort) {
          maxshort += maxshort/2;
          memory->grow(neighshort,maxshort,"pair:neighshort");
        }
      }

      jtag = tag[j];
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      if (rsq >= params[ijparam].cutsq) continue;

      twobody(&params[ijparam],rsq,fpair,eflag,evdwl);

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }

    jnumm1 = numshort - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = neighshort[jj];
      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 >= params[ijparam].cutsq2) continue;

      double fjxtmp,fjytmp,fjztmp;
      fjxtmp = fjytmp = fjztmp = 0.0;

      for (kk = jj+1; kk < numshort; kk++) {
        k = neighshort[kk];
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[ikparam].cutsq2) continue;

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
                  rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

        fxtmp -= fj[0] + fk[0];
        fytmp -= fj[1] + fk[1];
        fztmp -= fj[2] + fk[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairVashishta::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(neighshort,maxshort,"pair:neighshort");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairVashishta::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairVashishta::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg-3,arg+3);

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairVashishta::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style Vashishta requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style Vashishta requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairVashishta::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairVashishta::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "vashishta", unit_convert_flag);
    char * line;

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY,
                                                            unit_convert);

    while ((line = reader.next_line(NPARAMS_PER_LINE))) {
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
          params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                              "pair:params");

          // make certain all addional allocated storage is initialized
          // to avoid false positives when checking with valgrind

          memset(params + nparams, 0, DELTA*sizeof(Param));
        }

        params[nparams].ielement = ielement;
        params[nparams].jelement = jelement;
        params[nparams].kelement = kelement;
        params[nparams].bigh     = values.next_double();
        params[nparams].eta      = values.next_double();
        params[nparams].zi       = values.next_double();
        params[nparams].zj       = values.next_double();
        params[nparams].lambda1  = values.next_double();
        params[nparams].bigd     = values.next_double();
        params[nparams].lambda4  = values.next_double();
        params[nparams].bigw     = values.next_double();
        params[nparams].cut      = values.next_double();
        params[nparams].bigb     = values.next_double();
        params[nparams].gamma    = values.next_double();
        params[nparams].r0       = values.next_double();
        params[nparams].bigc     = values.next_double();
        params[nparams].costheta = values.next_double();

        if (unit_convert) {
          params[nparams].bigh *= conversion_factor;
          params[nparams].bigd *= conversion_factor;
          params[nparams].bigw *= conversion_factor;
          params[nparams].bigb *= conversion_factor;
        }

      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }

      if (params[nparams].bigb < 0.0 || params[nparams].gamma < 0.0 ||
          params[nparams].r0 < 0.0 || params[nparams].bigc < 0.0 ||
          params[nparams].bigh < 0.0 || params[nparams].eta < 0.0 ||
          params[nparams].lambda1 < 0.0 || params[nparams].bigd < 0.0 ||
          params[nparams].lambda4 < 0.0 || params[nparams].bigw < 0.0 ||
          params[nparams].cut < 0.0)
        error->one(FLERR,"Illegal Vashishta parameter");

      nparams++;
    }
  }

  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param), "pair:params");
  }

  MPI_Bcast(params, maxparam*sizeof(Param), MPI_BYTE, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairVashishta::setup_params()
{
  int i,j,k,m,n;

  // set elem3param for all triplet combinations
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

  // set cutsq using shortcut to reduce neighbor list for accelerated
  // calculations. cut must remain unchanged as it is a potential parameter
  double tmp_par;

  for (m = 0; m < nparams; m++) {
    params[m].cutsq = params[m].cut * params[m].cut;
    params[m].cutsq2 = params[m].r0 * params[m].r0;

    tmp_par = params[m].lambda1;
    params[m].lam1inv = (tmp_par == 0.0) ? 0.0 : 1.0/tmp_par;
    tmp_par = params[m].lambda4;
    params[m].lam4inv = (tmp_par == 0.0) ? 0.0 : 1.0/tmp_par;
    params[m].zizj = params[m].zi*params[m].zj * force->qqr2e;
    // note that bigd does not have 1/2 factor
    params[m].mbigd = params[m].bigd;
    params[m].heta = params[m].bigh*params[m].eta;
    params[m].big2b = 2.0*params[m].bigb;
    params[m].big6w = 6.0*params[m].bigw;

    tmp_par = params[m].cut;
    params[m].rcinv =  (tmp_par == 0.0) ? 0.0 : 1.0/tmp_par;
    params[m].rc2inv = params[m].rcinv*params[m].rcinv;
    params[m].rc4inv = params[m].rc2inv*params[m].rc2inv;
    params[m].rc6inv = params[m].rc2inv*params[m].rc4inv;
    params[m].rceta = pow(params[m].rcinv,params[m].eta);
    params[m].lam1rc = params[m].cut*params[m].lam1inv;
    params[m].lam4rc = params[m].cut*params[m].lam4inv;
    params[m].vrcc2 = params[m].zizj*params[m].rcinv *
      exp(-params[m].lam1rc);
    params[m].vrcc3 = params[m].mbigd*params[m].rc4inv *
      exp(-params[m].lam4rc);
    params[m].vrc = params[m].bigh*params[m].rceta +
      params[m].vrcc2 - params[m].vrcc3 -
      params[m].bigw*params[m].rc6inv;

    params[m].dvrc =
      params[m].vrcc3 * (4.0*params[m].rcinv+params[m].lam4inv)
      + params[m].big6w * params[m].rc6inv * params[m].rcinv
      - params[m].heta * params[m].rceta*params[m].rcinv
      - params[m].vrcc2 * (params[m].rcinv+params[m].lam1inv);
    params[m].c0 = params[m].cut*params[m].dvrc - params[m].vrc;
  }

  // set cutmax to max of all cutoff params. r0max only for r0

  cutmax = 0.0;
  r0max = 0.0;
  for (m = 0; m < nparams; m++) {
    if (params[m].cut > cutmax) cutmax = params[m].cut;
    if (params[m].r0 > r0max) r0max = params[m].r0;
  }
  if (r0max > cutmax) cutmax = r0max;
}

/* ---------------------------------------------------------------------- */

void PairVashishta::twobody(Param *param, double rsq, double &fforce,
                            int eflag, double &eng)
{
  double r,rinvsq,r4inv,r6inv,reta,lam1r,lam4r,vc2,vc3;

  r = sqrt(rsq);
  rinvsq = 1.0/rsq;
  r4inv = rinvsq*rinvsq;
  r6inv = rinvsq*r4inv;
  reta = pow(r,-param->eta);
  lam1r = r*param->lam1inv;
  lam4r = r*param->lam4inv;
  vc2 = param->zizj * exp(-lam1r)/r;
  vc3 = param->mbigd * r4inv*exp(-lam4r);

  fforce = (param->dvrc*r
            - (4.0*vc3 + lam4r*vc3+param->big6w*r6inv
               - param->heta*reta - vc2 - lam1r*vc2)
            ) * rinvsq;
  if (eflag) eng = param->bigh*reta
               + vc2 - vc3 - param->bigw*r6inv
               - r*param->dvrc + param->c0;
}

/* ---------------------------------------------------------------------- */

void PairVashishta::threebody(Param *paramij, Param *paramik, Param *paramijk,
                       double rsq1, double rsq2,
                       double *delr1, double *delr2,
                       double *fj, double *fk, int eflag, double &eng)
{
  double r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
  double r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
  double rinv12,cs,delcs,delcssq,facexp,facrad,frad1,frad2,pcsinv,pcsinvsq,pcs;
  double facang,facang12,csfacang,csfac1,csfac2;

  r1 = sqrt(rsq1);
  rinvsq1 = 1.0/rsq1;
  rainv1 = 1.0/(r1 - paramij->r0);
  gsrainv1 = paramij->gamma * rainv1;
  gsrainvsq1 = gsrainv1*rainv1/r1;
  expgsrainv1 = exp(gsrainv1);

  r2 = sqrt(rsq2);
  rinvsq2 = 1.0/rsq2;
  rainv2 = 1.0/(r2 - paramik->r0);
  gsrainv2 = paramik->gamma * rainv2;
  gsrainvsq2 = gsrainv2*rainv2/r2;
  expgsrainv2 = exp(gsrainv2);

  rinv12 = 1.0/(r1*r2);
  cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;
  delcs = cs - paramijk->costheta;
  delcssq = delcs*delcs;
  pcsinv = paramijk->bigc*delcssq + 1.0;
  pcsinvsq = pcsinv*pcsinv;
  pcs = delcssq/pcsinv;

  facexp = expgsrainv1*expgsrainv2;

  facrad = paramijk->bigb * facexp * pcs;
  frad1 = facrad*gsrainvsq1;
  frad2 = facrad*gsrainvsq2;
  facang = paramijk->big2b * facexp * delcs/pcsinvsq;
  facang12 = rinv12*facang;
  csfacang = cs*facang;
  csfac1 = rinvsq1*csfacang;

  fj[0] = delr1[0]*(frad1+csfac1)-delr2[0]*facang12;
  fj[1] = delr1[1]*(frad1+csfac1)-delr2[1]*facang12;
  fj[2] = delr1[2]*(frad1+csfac1)-delr2[2]*facang12;

  csfac2 = rinvsq2*csfacang;

  fk[0] = delr2[0]*(frad2+csfac2)-delr1[0]*facang12;
  fk[1] = delr2[1]*(frad2+csfac2)-delr1[1]*facang12;
  fk[2] = delr2[2]*(frad2+csfac2)-delr1[2]*facang12;

  if (eflag) eng = facrad;
}
