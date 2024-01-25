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
   Contributing authors: Axel Kohlmeyer (Temple U), Venkatesh Botu, James Chapman (Lawrence Livermore National Lab)
------------------------------------------------------------------------- */

#include "pair_agni.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "force.h"
#include "error.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathSpecial;

static const char cite_pair_agni[] =
  "pair agni command: doi:10.1021/acs.jpcc.9b04207\n\n"
  "@article{huan2019jpc,\n"
  " author    = {Huan, T. and Batra, R. and Chapman, J. and Kim, C. and Chandrasekaran, A. and Ramprasad, Rampi},\n"
  " journal   = {J.~Phys.\\ Chem.~C},\n"
  " volume    = {123},\n"
  " number    = {34},\n"
  " pages     = {20715--20722},\n"
  " year      = {2019},\n"
  "}\n\n";

static constexpr int MAXLINE = 10240;
static constexpr int MAXWORD = 40;

/* ---------------------------------------------------------------------- */

PairAGNI::PairAGNI(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_agni);

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  atomic_feature_version = 0;

  centroidstressflag = CENTROID_NOTAVAIL;

  no_virial_fdotr_compute = 1;

  params = nullptr;
  cutmax = 0.0;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairAGNI::~PairAGNI()
{
  if (params) {
    for (int i = 0; i < nparams; ++i) {
      memory->destroy(params[i].eta);
      memory->destroy(params[i].alpha);
      memory->destroy(params[i].xU);
    }
    memory->destroy(params);
    params = nullptr;
  }
  memory->destroy(elem1param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairAGNI::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,inum,jnum,itype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double fxtmp,fytmp,fztmp;
  double *Vx, *Vy, *Vz;

  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    const Param &iparam = params[elem1param[itype]];
    Vx = new double[iparam.numeta];
    Vy = new double[iparam.numeta];
    Vz = new double[iparam.numeta];
    memset(Vx,0,iparam.numeta*sizeof(double));
    memset(Vy,0,iparam.numeta*sizeof(double));
    memset(Vz,0,iparam.numeta*sizeof(double));

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if ((rsq > 0.0) && (rsq < iparam.cutsq)) {
        const double r = sqrt(rsq);
        const double cF = 0.5*(cos((MathConst::MY_PI*r)/iparam.cut)+1.0);
        const double wX = cF*delx/r;
        const double wY = cF*dely/r;
        const double wZ = cF*delz/r;

        for (k = 0; k < iparam.numeta; ++k) {
          double e = 0.0;

          if (atomic_feature_version == AGNI_VERSION_1)
            e = exp(-(iparam.eta[k]*rsq));
          else if (atomic_feature_version == AGNI_VERSION_2)
            e = (1.0 / (square(iparam.eta[k]) * iparam.gwidth * sqrt(MathConst::MY_2PI))) * exp(-(square(r - iparam.eta[k])) / (2.0 * square(iparam.gwidth)));

          Vx[k] += wX*e;
          Vy[k] += wY*e;
          Vz[k] += wZ*e;
        }
      }
    }

    for (j = 0; j < iparam.numtrain; ++j) {
      double kx = 0.0;
      double ky = 0.0;
      double kz = 0.0;

      for (int k = 0; k < iparam.numeta; ++k) {
        const double xu = iparam.xU[k][j];

        kx += square(Vx[k] - xu);
        ky += square(Vy[k] - xu);
        kz += square(Vz[k] - xu);
      }
      const double e = -0.5/(square(iparam.sigma));
      fxtmp += iparam.alpha[j]*exp(kx*e);
      fytmp += iparam.alpha[j]*exp(ky*e);
      fztmp += iparam.alpha[j]*exp(kz*e);
    }
    fxtmp += iparam.b;
    fytmp += iparam.b;
    fztmp += iparam.b;
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;

    if (evflag) ev_tally_xyz_full(i,0.0,0.0,fxtmp,fytmp,fztmp,delx,dely,delz);

    delete [] Vx;
    delete [] Vy;
    delete [] Vz;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairAGNI::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairAGNI::settings(int narg, char ** /* arg */)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairAGNI::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg-3,arg+3);

  if (nelements != 1)
    error->all(FLERR,"Cannot handle multi-element systems with this potential");

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAGNI::init_style()
{
  // need a full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAGNI::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairAGNI::read_file(char *filename)
{
  int i,j,k,curparam,wantdata,fp_counter;

  if (params) {
    for (int i = 0; i < nparams; ++i) {
      memory->destroy(params[i].eta);
      memory->destroy(params[i].alpha);
      memory->destroy(params[i].xU);
    }
    memory->destroy(params);
    params = nullptr;
  }
  nparams = 0;

  fp_counter = 0;
  wantdata = -1;

  // read potential file
  if (comm->me == 0) {
    PotentialFileReader reader(lmp, filename, "agni", unit_convert_flag);

    char * line;

    while ((line = reader.next_line())) {
      try {
        ValueTokenizer values(line);
        if (wantdata == -1) {
          std::string tag = values.next_string();

          if (tag == "n_elements") {
            nparams = values.next_int();

            if ((nparams < 1) || params) // sanity check
              error->all(FLERR,"Invalid AGNI potential file");
            params = memory->create(params,nparams,"pair:params");
            memset(params,0,nparams*sizeof(Param));

            curparam = -1;
            wantdata = -1;
          } else if (tag == "interaction") {
            for (j = 0; j < nparams; ++j) {
              std::string element = values.next_string();
              if (element == elements[params[j].ielement])
                curparam = j;
            }
          } else if (tag == "element") {
            for (j = 0; j < nparams; ++j) {
              std::string element = values.next_string();
              for (k = 0; k < nelements; ++k)
                if (element == elements[k]) break;
              if (k == nelements)
                error->all(FLERR,"No suitable parameters for requested element found");
              else params[j].ielement = k;
            }
          } else if (tag == "generation") {
            atomic_feature_version = values.next_int();
            if (atomic_feature_version != AGNI_VERSION_1 && atomic_feature_version != AGNI_VERSION_2)
              error->all(FLERR,"Incompatible AGNI potential file version");
          } else if (tag == "eta") {
            params[curparam].numeta = values.count() - 1;
            memory->create(params[curparam].eta,params[curparam].numeta,"agni:eta");
            for (j = 0; j < params[curparam].numeta; j++)
              params[curparam].eta[j] = values.next_double();
          } else if (tag == "gwidth") {
            params[curparam].gwidth = values.next_double();
          } else if (tag == "Rc") {
            params[curparam].cut = values.next_double();
          } else if (tag == "n_train") {
            params[curparam].numtrain = values.next_int();
            memory->create(params[curparam].alpha,params[curparam].numtrain,"agni:alpha");
            memory->create(params[curparam].xU,params[curparam].numtrain,params[curparam].numtrain,"agni:xU");
          } else if (tag == "sigma") {
            params[curparam].sigma = values.next_double();
          } else if (tag == "b") {
            params[curparam].b = values.next_double();
          } else if (tag == "endVar") {
            if (atomic_feature_version == AGNI_VERSION_1)
              params[curparam].gwidth = 0.0;
            wantdata = curparam;
            curparam = -1;
          } else error->warning(FLERR,"Ignoring unknown tag '{}' in AGNI potential file.",tag);
        } else {
          if (params && wantdata >= 0) {
            if ((int)values.count() == params[wantdata].numeta + 2) {
              for (k = 0; k < params[wantdata].numeta; ++k)
                params[wantdata].xU[k][fp_counter] = values.next_double();
              values.next_double(); // ignore
              params[wantdata].alpha[fp_counter] = values.next_double();
              fp_counter++;
            } else if ((int)values.count() == params[wantdata].numeta + 3) {
              values.next_double(); // ignore
              for (k = 0; k < params[wantdata].numeta; ++k)
                params[wantdata].xU[k][fp_counter] = values.next_double();
              values.next_double(); // ignore
              params[wantdata].alpha[fp_counter] = values.next_double();
              fp_counter++;
            } else error->all(FLERR,"Invalid AGNI potential file");
          }
        }
      }
      catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }
    }
  }
  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&atomic_feature_version, 1, MPI_INT, 0, world);
  if (comm->me != 0) {
    params = memory->create(params,nparams,"pair:params");
    memset(params,0,nparams*sizeof(Param));
  }
  MPI_Bcast(params, nparams*sizeof(Param), MPI_BYTE, 0, world);
  for (i = 0; i < nparams; i++) {
    if (comm->me != 0) {
      memory->create(params[i].alpha,params[i].numtrain,"agni:alpha");
      memory->create(params[i].eta,params[i].numeta,"agni:eta");
      memory->create(params[i].xU,params[i].numeta,params[i].numtrain,"agni:xU");
    }

    MPI_Bcast(params[i].alpha, params[i].numtrain, MPI_DOUBLE, 0, world);
    MPI_Bcast(params[i].eta, params[i].numeta, MPI_DOUBLE, 0, world);
    MPI_Bcast(&params[i].xU[0][0],params[i].numtrain*params[i].numeta,MPI_DOUBLE,0,world);
  }
}

/* ---------------------------------------------------------------------- */

void PairAGNI::setup_params()
{
  int i,m,n;
  double rtmp;

  // set elem1param for all elements

  memory->destroy(elem1param);
  memory->create(elem1param,nelements,"pair:elem1param");

  for (i = 0; i < nelements; i++) {
    n = -1;
    for (m = 0; m < nparams; m++) {
      if (i == params[m].ielement) {
        if (n >= 0) error->all(FLERR,"Potential file has a duplicate entry for: {}", elements[i]);
        n = m;
      }
    }
    if (n < 0) error->all(FLERR,"Potential file is missing an entry for: {}", elements[i]);
    elem1param[i] = n;
  }

  // compute parameter values derived from inputs

  // set cutsq using shortcut to reduce neighbor list for accelerated
  // calculations. cut must remain unchanged as it is a potential parameter
  // (cut = a*sigma)

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = params[m].cut;
    params[m].cutsq = rtmp * rtmp;
    if (rtmp > cutmax) cutmax = rtmp;
  }
}
