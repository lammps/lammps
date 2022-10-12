// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Ray Shan (Sandia, tnshan@sandia.gov)
   Updates and debug: Tao Liang (U Florida, liang75@ufl.edu)
                      Dundar Yilmaz (dundar.yilmaz@zirve.edu.tr)
------------------------------------------------------------------------- */

#include "pair_comb3.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "text_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathExtra;
using namespace MathSpecial;

#define DELTA 4
#define PGDELTA 1
#define MAXNEIGH 24

/* ---------------------------------------------------------------------- */

PairComb3::PairComb3(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  ghostneigh = 1;

  nmax = 0;
  NCo = nullptr;
  bbij = nullptr;
  map = nullptr;
  esm = nullptr;

  params = nullptr;

  intype = nullptr;
  afb = nullptr;
  dafb = nullptr;
  fafb = nullptr;
  dfafb = nullptr;
  ddfafb = nullptr;
  phin = nullptr;
  dphin = nullptr;
  erpaw = nullptr;
  vvdw = nullptr;
  vdvdw = nullptr;
  dpl = nullptr;
  xcctmp = nullptr;
  xchtmp = nullptr;
  xcotmp = nullptr;

  sht_num = nullptr;
  sht_first = nullptr;

  ipage = nullptr;
  pgsize = oneatom = 0;

  cflag = 0;

  // set comm size needed by this Pair
  comm_forward = 1;
  comm_reverse = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairComb3::~PairComb3()
{
  memory->destroy(NCo);

  memory->sfree(params);
  memory->destroy(elem3param);

  memory->destroy(afb);
  memory->destroy(dpl);
  memory->destroy(dafb);
  memory->destroy(fafb);
  memory->destroy(phin);
  memory->destroy(bbij);
  memory->destroy(vvdw);
  memory->destroy(vdvdw);
  memory->destroy(dphin);
  memory->destroy(erpaw);
  memory->destroy(dfafb);
  memory->destroy(ddfafb);
  memory->destroy(xcctmp);
  memory->destroy(xchtmp);
  memory->destroy(xcotmp);
  memory->destroy(intype);
  memory->destroy(sht_num);
  memory->sfree(sht_first);

  delete[] ipage;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);
    delete[] esm;
  }

}

/* ---------------------------------------------------------------------- */

void PairComb3::allocate()
{
 allocated = 1;
 int n = atom->ntypes;

 memory->create(setflag,n+1,n+1,"pair:setflag");
 memory->create(cutsq,n+1,n+1,"pair:cutsq");
 memory->create(cutghost,n+1,n+1,"pair:cutghost");

 map = new int[n+1];
 esm = new double[n];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairComb3::settings(int narg, char **arg)
{

  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  if (strcmp(arg[0],"polar_on") == 0) pol_flag = 1;
  else if (strcmp(arg[0],"polar_off") == 0) pol_flag = 0;
  else error->all(FLERR,"Illegal pair_style command");

  if (comm->me == 0 && screen)
    fmt::print(screen,"   PairComb3: polarization is {} \n",
               pol_flag ? "on" : "off");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairComb3::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  cflag = 0;
  for (int i = 3; i < narg; i++)
    if (strcmp(arg[i],"C") == 0) {
      cflag = 1;
      break;
    }

  if (cflag) {
    if (comm->me == 0 && screen)
      fputs(" PairComb3: Found C: reading additional library file\n",screen);
    read_lib();
  }

  map_element2type(narg-3,arg+3);

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();

  // generate Wolf 1/r energy and van der Waals look-up tables

  tables();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairComb3::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style COMB3 requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style COMB3 requires newton pair on");
  if (!atom->q_flag)
    error->all(FLERR,"Pair style COMB3 requires atom attribute q");

// need a full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);

  // local Comb neighbor list
  // create pages if first time or if neighbor pgsize/oneatom has changed

  int create = 0;
  if (ipage == nullptr) create = 1;
  if (pgsize != neighbor->pgsize) create = 1;
  if (oneatom != neighbor->oneatom) create = 1;

  if (create) {
    delete[] ipage;
    pgsize = neighbor->pgsize;
    oneatom = neighbor->oneatom;

    int nmypage = comm->nthreads;
    ipage = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      ipage[i].init(oneatom,pgsize);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairComb3::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  cutghost[j][i] = cutghost[i][j] = cutmax;
  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairComb3::read_lib()
{
  int i,j,k,l;
  int ii,jj,kk,ll,mm,iii;

  // open library file on proc 0

  if (comm->me == 0) {
    try {
      PotentialFileReader reader(lmp, "lib.comb3", "comb3");
      reader.next_dvector(ccutoff, 6);
      reader.next_dvector(ch_a, 7);

      ValueTokenizer values = reader.next_values(3);
      nsplpcn = values.next_int();
      nsplrad = values.next_int();
      nspltor = values.next_int();

      values = reader.next_values(3);
      maxx = values.next_int();
      maxy = values.next_int();
      maxz = values.next_int();

      values = reader.next_values(3);
      maxxc   = values.next_int();
      maxyc   = values.next_int();
      maxconj = values.next_int();

      for (l = 0; l < nsplpcn; l++) {
        values = reader.next_values(3);
        values.skip(1);
        maxxcn[l]   = values.next_int();
        vmaxxcn[l]  = values.next_double();
        dvmaxxcn[l] = values.next_double();
      }

      ntab = reader.next_int();

      for (i = 0; i < (ntab + 1); i++) {
        values = reader.next_values(4);
        values.skip(1);
        pang[i]   = values.next_double();
        dpang[i]  = values.next_double();
        ddpang[i] = values.next_double();
      }

      for (l = 0; l < nsplpcn; l++)
        for (i = 0; i < (maxx + 1); i++)
          for (j = 0; j < (maxy + 1); j++)
            for (k = 0; k < (maxz + 1); k++) {
              values = reader.next_values(8);

              ll = values.next_int() - 1;
              ii = values.next_int();
              jj = values.next_int();
              kk = values.next_int();
              pcn_grid[ll][ii][jj][kk]  = values.next_double();
              pcn_gridx[ll][ii][jj][kk] = values.next_double();
              pcn_gridy[ll][ii][jj][kk] = values.next_double();
              pcn_gridz[ll][ii][jj][kk] = values.next_double();
            }

      for (l = 0; l < nsplpcn; l++)
        for (i = 0; i < maxx; i++)
          for (j = 0; j < maxy; j++)
            for (k = 0; k < maxz; k++) {
              values = reader.next_values(4);

              ll = values.next_int()-1;
              ii = values.next_int();
              jj = values.next_int();
              kk = values.next_int();

              for (iii = 0; iii < 2; iii++) {
                mm = iii*32;
                reader.next_dvector(&pcn_cubs[ll][ii][jj][kk][mm], 32);
              }
            }

      for (l = 0; l < nsplrad; l++)
        for (i = 0; i < (maxxc + 1); i++)
          for (j = 0; j < (maxyc + 1); j++)
            for (k = 0; k < maxconj; k++) {
              values = reader.next_values(8);

              ll = values.next_int() - 1;
              ii = values.next_int();
              jj = values.next_int();
              kk = values.next_int() - 1;
              rad_grid[ll][ii][jj][kk]  = values.next_double();
              rad_gridx[ll][ii][jj][kk] = values.next_double();
              rad_gridy[ll][ii][jj][kk] = values.next_double();
              rad_gridz[ll][ii][jj][kk] = values.next_double();
            }

      for (l = 0; l < nsplrad; l++)
        for (i = 0; i < maxxc; i++)
          for (j = 0; j < maxyc; j++)
            for (k = 0; k < (maxconj - 1); k++) {
              values = reader.next_values(4);

              ll = values.next_int() - 1;
              ii = values.next_int();
              jj = values.next_int();
              kk = values.next_int() - 1;

              for (iii = 0; iii < 2; iii++) {
                mm = iii * 32;
                reader.next_dvector(&rad_spl[ll][ii][jj][kk][mm], 32);
              }
            }

      for (l=0; l<nspltor; l++)
        for (i=0; i<maxxc+1; i++)
          for (j=0; j<maxyc+1; j++)
            for (k=0; k<maxconj; k++) {
              values = reader.next_values(8);

              ll = values.next_int() - 1;
              ii = values.next_int();
              jj = values.next_int();
              kk = values.next_int() - 1;
              tor_grid[ll][ii][jj][kk]  = values.next_double();
              tor_gridx[ll][ii][jj][kk] = values.next_double();
              tor_gridy[ll][ii][jj][kk] = values.next_double();
              tor_gridz[ll][ii][jj][kk] = values.next_double();
            }

      for (l=0; l<nspltor; l++)
        for (i=0; i<maxxc; i++)
          for (j=0; j<maxyc; j++)
            for (k=0; k<maxconj-1; k++) {
              values = reader.next_values(4);

              ll = values.next_int() - 1;
              ii = values.next_int();
              jj = values.next_int();
              kk = values.next_int() - 1;
              for (iii=0; iii<2; iii++) {
                mm=iii*32;
                reader.next_dvector(&tor_spl[ll][ii][jj][kk][mm], 32);
              }
            }

      k = 0;
      for (i=0; i<4; i++)
        for (j=0; j<4; j++) {
          iin2[k][0] = i;
          iin2[k][1] = j;
          k++;
        }

      l = 0;
      for (i=0; i<4; i++)
        for (j=0; j<4; j++)
          for (k=0; k<4; k++) {
            iin3[l][0] = i;
            iin3[l][1] = j;
            iin3[l][2] = k;
            l++;
          }

    } catch (FileReaderException &fre) {
      error->one(FLERR, fre.what());
    } catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }
  }

  MPI_Bcast(&ccutoff[0],6,MPI_DOUBLE,0,world);
  MPI_Bcast(&ch_a[0],7,MPI_DOUBLE,0,world);
  MPI_Bcast(&nsplpcn,1,MPI_INT,0,world);
  MPI_Bcast(&nsplrad,1,MPI_INT,0,world);
  MPI_Bcast(&nspltor,1,MPI_INT,0,world);
  MPI_Bcast(&maxx,1,MPI_INT,0,world);
  MPI_Bcast(&maxy,1,MPI_INT,0,world);
  MPI_Bcast(&maxz,1,MPI_INT,0,world);
  MPI_Bcast(&maxxc,1,MPI_INT,0,world);
  MPI_Bcast(&maxyc,1,MPI_INT,0,world);
  MPI_Bcast(&maxconj,1,MPI_INT,0,world);
  MPI_Bcast(&maxxcn,4,MPI_INT,0,world);
  MPI_Bcast(&vmaxxcn,4,MPI_DOUBLE,0,world);
  MPI_Bcast(&dvmaxxcn,4,MPI_DOUBLE,0,world);
  MPI_Bcast(&ntab,1,MPI_INT,0,world);
  MPI_Bcast(&pang[0],20001,MPI_DOUBLE,0,world);
  MPI_Bcast(&dpang[0],20001,MPI_DOUBLE,0,world);
  MPI_Bcast(&ddpang[0],20001,MPI_DOUBLE,0,world);
  MPI_Bcast(&pcn_grid[0][0][0][0],500,MPI_DOUBLE,0,world);
  MPI_Bcast(&pcn_gridx[0][0][0][0],500,MPI_DOUBLE,0,world);
  MPI_Bcast(&pcn_gridy[0][0][0][0],500,MPI_DOUBLE,0,world);
  MPI_Bcast(&pcn_gridz[0][0][0][0],500,MPI_DOUBLE,0,world);
  MPI_Bcast(&pcn_cubs[0][0][0][0][0],16384,MPI_DOUBLE,0,world);

  MPI_Bcast(&rad_grid[0][0][0][0],825,MPI_DOUBLE,0,world);
  MPI_Bcast(&rad_gridx[0][0][0][0],825,MPI_DOUBLE,0,world);
  MPI_Bcast(&rad_gridy[0][0][0][0],825,MPI_DOUBLE,0,world);
  MPI_Bcast(&rad_gridz[0][0][0][0],825,MPI_DOUBLE,0,world);
  MPI_Bcast(&rad_spl[0][0][0][0][0],30720,MPI_DOUBLE,0,world);

  MPI_Bcast(&tor_grid[0][0][0][0],275,MPI_DOUBLE,0,world);
  MPI_Bcast(&tor_gridx[0][0][0][0],275,MPI_DOUBLE,0,world);
  MPI_Bcast(&tor_gridy[0][0][0][0],275,MPI_DOUBLE,0,world);
  MPI_Bcast(&tor_gridz[0][0][0][0],275,MPI_DOUBLE,0,world);
  MPI_Bcast(&tor_spl[0][0][0][0][0],10240,MPI_DOUBLE,0,world);

  MPI_Bcast(&iin2[0][0],32,MPI_DOUBLE,0,world);
  MPI_Bcast(&iin3[0][0],192,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

void PairComb3::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = 0;
  maxparam = 0;

  // open file on proc 0
  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "COMB3");
    char * line;

    while ((line = reader.next_line(NPARAMS_PER_LINE))) {
      try {
        ValueTokenizer values(line);

        std::string iname = values.next_string();
        std::string jname = values.next_string();
        std::string kname = values.next_string();

        // ielement,jelement,kelement = 1st args
        // if all 3 args are in element list, then parse this line
        // else skip to next line
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
        params[nparams].ielementgp = values.next_int();
        params[nparams].jelementgp = values.next_int();
        params[nparams].kelementgp = values.next_int();
        params[nparams].ang_flag   = values.next_int();
        params[nparams].pcn_flag   = values.next_int();
        params[nparams].rad_flag   = values.next_int();
        params[nparams].tor_flag   = values.next_int();
        params[nparams].vdwflag    = values.next_double();
        params[nparams].powerm     = values.next_double();
        params[nparams].veps       = values.next_double();
        params[nparams].vsig       = values.next_double();
        params[nparams].paaa       = values.next_double();
        params[nparams].pbbb       = values.next_double();
        params[nparams].lami       = values.next_double();
        params[nparams].alfi       = values.next_double();
        params[nparams].powern     = values.next_double();
        params[nparams].QL         = values.next_double();
        params[nparams].QU         = values.next_double();
        params[nparams].DL         = values.next_double();
        params[nparams].DU         = values.next_double();
        params[nparams].qmin       = values.next_double();
        params[nparams].qmax       = values.next_double();
        params[nparams].chi        = values.next_double();
        params[nparams].dj         = values.next_double();
        params[nparams].dk         = values.next_double();
        params[nparams].dl         = values.next_double();
        params[nparams].esm        = values.next_double();
        params[nparams].cmn1       = values.next_double();
        params[nparams].cmn2       = values.next_double();
        params[nparams].pcmn1      = values.next_double();
        params[nparams].pcmn2      = values.next_double();
        params[nparams].coulcut    = values.next_double();
        params[nparams].polz       = values.next_double();
        params[nparams].curl       = values.next_double();
        params[nparams].curlcut1   = values.next_double();
        params[nparams].curlcut2   = values.next_double();
        params[nparams].curl0      = values.next_double();
        params[nparams].alpha1     = values.next_double();
        params[nparams].bigB1      = values.next_double();
        params[nparams].alpha2     = values.next_double();
        params[nparams].bigB2      = values.next_double();
        params[nparams].alpha3     = values.next_double();
        params[nparams].bigB3      = values.next_double();
        params[nparams].lambda     = values.next_double();
        params[nparams].bigA       = values.next_double();
        params[nparams].beta       = values.next_double();
        params[nparams].bigr       = values.next_double();
        params[nparams].bigd       = values.next_double();
        params[nparams].pcos6      = values.next_double();
        params[nparams].pcos5      = values.next_double();
        params[nparams].pcos4      = values.next_double();
        params[nparams].pcos3      = values.next_double();
        params[nparams].pcos2      = values.next_double();
        params[nparams].pcos1      = values.next_double();
        params[nparams].pcos0      = values.next_double();
        params[nparams].pcna       = values.next_double();
        params[nparams].pcnb       = values.next_double();
        params[nparams].pcnc       = values.next_double();
        params[nparams].pcnd       = values.next_double();
        params[nparams].p6p0       = values.next_double();
        params[nparams].p6p1       = values.next_double();
        params[nparams].p6p2       = values.next_double();
        params[nparams].p6p3       = values.next_double();
        params[nparams].p6p4       = values.next_double();
        params[nparams].p6p5       = values.next_double();
        params[nparams].p6p6       = values.next_double();
        params[nparams].ptork1     = values.next_double();
        params[nparams].ptork2     = values.next_double();
        params[nparams].addrepr    = values.next_double();
        params[nparams].addrep     = values.next_double();
        params[nparams].pcross     = values.next_double();
        params[nparams].powermint = int(params[nparams].powerm);
      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }

      // parameter sanity checks

      if (params[nparams].lambda < 0.0 || params[nparams].powern < 0.0 ||
          params[nparams].beta < 0.0 || params[nparams].alpha1 < 0.0 ||
          params[nparams].bigB1< 0.0 || params[nparams].bigA< 0.0 ||
          params[nparams].bigB2< 0.0 || params[nparams].alpha2 <0.0 ||
          params[nparams].bigB3< 0.0 || params[nparams].alpha3 <0.0 ||
          params[nparams].bigr < 0.0 || params[nparams].bigd < 0.0 ||
          params[nparams].bigd > params[nparams].bigr ||
          params[nparams].powerm - params[nparams].powermint != 0.0 ||
          params[nparams].addrepr < 0.0 || params[nparams].powermint < 1.0 ||
          params[nparams].QL > 0.0 || params[nparams].QU < 0.0 ||
          params[nparams].DL < 0.0 || params[nparams].DU > 0.0 ||
          params[nparams].pcross < 0.0 ||
          params[nparams].esm < 0.0 || params[nparams].veps < 0.0 ||
          params[nparams].vsig < 0.0 || params[nparams].vdwflag < 0.0
          )
        error->one(FLERR,"Illegal COMB3 parameter");

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

void PairComb3::setup_params()
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
    params[m].c2 = pow(2.0*params[m].powern*1.00e-8,-1.0/params[m].powern);
    params[m].c3 = 1.0/params[m].c2;
    params[m].c4 = 1.0/params[m].c1;

    params[m].Qo = (params[m].QU+params[m].QL)/2.0; // (A22)
    params[m].dQ = (params[m].QU-params[m].QL)/2.0; // (A21)
    params[m].aB = 1.0 /
      (1.0-powint(fabs(params[m].Qo/params[m].dQ),10)); // (A20)
    params[m].bB = pow(fabs(params[m].aB),0.1)/params[m].dQ; // (A19)
    params[m].nD = log(params[m].DU/(params[m].DU-params[m].DL))/
                    log(params[m].QU/(params[m].QU-params[m].QL));
    params[m].bD = (pow((params[m].DL-params[m].DU),(1.0/params[m].nD)))/
                    (params[m].QU-params[m].QL);

    params[m].lcut = params[m].coulcut;
    params[m].lcutsq = params[m].lcut*params[m].lcut;
  }

  // set cutmax to max of all params

  cutmin = cutmax = 0.0;
  polar = 0;
  for (m = 0; m < nparams; m++) {
    if (params[m].cutsq > cutmin) cutmin = params[m].cutsq + 2.0;
    if (params[m].lcut > cutmax) cutmax = params[m].lcut;
  }
  chicut1 = 7.0;
  chicut2 = cutmax;
}

/* ---------------------------------------------------------------------- */

void PairComb3::Short_neigh()
{
  int nj,*neighptrj,icontrol;
  int iparam_ij,*ilist,*jlist,*numneigh,**firstneigh;
  int inum,jnum,i,j,ii,jj,itype,jtype;
  double rr1,rsq1,delrj[3];

  double **x = atom->x;
  int *type  = atom->type;

  if (atom->nmax > nmax) {
    memory->sfree(sht_first);
    nmax = atom->nmax;
    sht_first = (int **) memory->smalloc(nmax*sizeof(int *),
                                        "pair:sht_first");
    memory->grow(dpl,nmax,3,"pair:dpl");
    memory->grow(xcctmp,nmax,"pair:xcctmp");
    memory->grow(xchtmp,nmax,"pair:xchtmp");
    memory->grow(xcotmp,nmax,"pair:xcotmp");
    memory->grow(NCo,nmax,"pair:NCo");
    memory->grow(sht_num,nmax,"pair:sht_num");
    memory->grow(bbij,nmax,MAXNEIGH,"pair:bbij");
  }

  inum = list->inum + list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // create COMB neighbor list

  ipage->reset();

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    dpl[i][0] = dpl[i][1] = dpl[i][2] = 0.0;

    nj = 0;
    neighptrj = ipage->vget();

    itype = map[type[i]];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    NCo[i] = 0.0;
    xcctmp[i] = 0.0;
    xchtmp[i] = 0.0;
    xcotmp[i] = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj] & NEIGHMASK;

      delrj[0] = x[i][0] - x[j][0];
      delrj[1] = x[i][1] - x[j][1];
      delrj[2] = x[i][2] - x[j][2];
      rsq1 = dot3(delrj,delrj);
      jtype = map[type[j]];
      iparam_ij = elem3param[itype][jtype][jtype];

      if (rsq1 > cutmin) continue;

      neighptrj[nj++] = j;
      rr1 = sqrt(rsq1);
      NCo[i] += comb_fc(rr1,&params[iparam_ij]) * params[iparam_ij].pcross;

      icontrol = params[iparam_ij].jelementgp;

      if (icontrol == 1)
          xcctmp[i] += comb_fc(rr1,&params[iparam_ij]) * params[iparam_ij].pcross;
      if (icontrol == 2)
          xchtmp[i] += comb_fc(rr1,&params[iparam_ij]) * params[iparam_ij].pcross;
      if (icontrol == 3)
          xcotmp[i] += comb_fc(rr1,&params[iparam_ij]) * params[iparam_ij].pcross;

    }

    sht_first[i] = neighptrj;
    sht_num[i] = nj;
    ipage->vgot(nj);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  // communicating coordination number to all nodes
  pack_flag = 2;
  comm->forward_comm(this);

}

/* ---------------------------------------------------------------------- */

void PairComb3::compute(int eflag, int vflag)
{
  int i,ii,k,kk,j,jj,im,inum,jnum,itype,jtype,ktype;
  int iparam_i,iparam_ij,iparam_ji;
  int iparam_ijk,iparam_jik,iparam_ikj,iparam_jli,iparam_ikl;
  int sht_jnum,*sht_jlist,sht_lnum,*sht_llist;
  int sht_mnum,*sht_mlist,sht_pnum,*sht_plist;
  int *ilist,*jlist,*numneigh,**firstneigh,mr1,mr2,mr3,inty,nj;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2,rsq3,iq,jq,yaself;
  double eng_tmp,vionij,fvionij,sr1,sr2,sr3;
  double zeta_ij,prefac_ij1,prefac_ij2,prefac_ij3,prefac_ij4,prefac_ij5;
  double zeta_ji,prefac_ji1,prefac_ji2,prefac_ji3,prefac_ji4,prefac_ji5;
  double delrj[3],delrk[3],fi[3],fj[3],fk[3],fl[3];
  double ep6p_ij,ep6p_ji,fip6p[3],fjp6p[3],fkp6p[3],flp6p[3];
  double potal,fac11,fac11e;
  tagint itag, jtag;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  tagint *tag = atom->tag;
  int *type = atom->type;
  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;

  // coordination terms
  double xcn, ycn;
  double kcn, lcn;
  int torindx;

  // torsion and radical variables
  int l, ll, ltype, m, mm, mtype, p, pp, ptype;
  int iparam_jil, iparam_ijl, iparam_ki, iparam_lj;
  int iparam_jl, iparam_ik, iparam_km, iparam_lp;
  double kconjug, lconjug, kradtot, lradtot;
  double delrl[3], delrm[3], delrp[3], ddprx[3], srmu;
  double zet_addi,zet_addj;

  evdwl = eng_tmp = 0.0;

  ev_init(eflag,vflag);

  // Build short range neighbor list
  Short_neigh();

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  yaself = vionij = fvionij = fpair = 0.0;

  // self energy correction term: potal
  potal_calc(potal,fac11,fac11e);

  // generate initial dipole tensor
  if (pol_flag )
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      itag = tag[i];
      itype = map[type[i]];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      iq = q[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj] & NEIGHMASK;
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
        jq = q[j];

        delrj[0] = x[j][0] - xtmp;
        delrj[1] = x[j][1] - ytmp;
        delrj[2] = x[j][2] - ztmp;
        rsq = dot3(delrj,delrj);

        iparam_ij = elem3param[itype][jtype][jtype];
        iparam_ji = elem3param[jtype][itype][itype];

        if (rsq > params[iparam_ij].lcutsq) continue;

        tri_point(rsq, mr1, mr2, mr3, sr1, sr2, sr3);

        dipole_init(&params[iparam_ij],&params[iparam_ji],fac11,delrj,
                    rsq,mr1,mr2,mr3,sr1,sr2,sr3,iq,jq,i,j);
      }
    }

  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    iq = q[i];
    nj = 0;
    iparam_i = elem3param[itype][itype][itype];

    // self energy, only on i atom
    yaself = self(&params[iparam_i],iq);

    // dipole self energy
    if (pol_flag)
      yaself += dipole_self(&params[iparam_i],i);

    if (evflag) ev_tally(i,i,nlocal,0,0.0,yaself,0.0,0.0,0.0,0.0);

    // two-body interactions (long:R + A, short: only R)

    jlist = firstneigh[i];
    jnum = numneigh[i];
    sht_jlist = sht_first[i];
    sht_jnum = sht_num[i];

    // long range interactions
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj] & NEIGHMASK;

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
      jq = q[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      iparam_ij = elem3param[itype][jtype][jtype];
      iparam_ji = elem3param[jtype][itype][itype];

      if (rsq > params[iparam_ij].lcutsq) continue;

      inty = intype[itype][jtype];

      // three-point interpolation
      tri_point(rsq, mr1, mr2, mr3, sr1, sr2, sr3);

      // Q-indenpendent: van der Waals
      vdwaals(inty,mr1,mr2,mr3,rsq,sr1,sr2,sr3,eng_tmp,fpair);
      evdwl = eng_tmp;
      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag)
        ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);

      // Q-dependent: Coulombic, field, polarization
      // 1/r energy and forces

      direct(&params[iparam_ij], &params[iparam_ji],
                mr1, mr2, mr3, rsq, sr1, sr2, sr3, iq, jq,
                fac11, fac11e, eng_tmp, fvionij, i, j);

      vionij = eng_tmp;

      // field correction to self energy
      field(&params[iparam_ij], &params[iparam_ji],rsq,iq,jq,
             eng_tmp,fvionij);
      vionij += eng_tmp;

      // sums up long range Q-dependent forces (excluding dipole)
      f[i][0] += delx*fvionij;
      f[i][1] += dely*fvionij;
      f[i][2] += delz*fvionij;
      f[j][0] -= delx*fvionij;
      f[j][1] -= dely*fvionij;
      f[j][2] -= delz*fvionij;

      // sums up long range Q-dependent energies (excluding dipole)
      if (evflag)
        ev_tally(i,j,nlocal,newton_pair,0.0,vionij,fvionij,delx,dely,delz);

      // polarization field
      if (pol_flag) {
        dipole_calc(&params[iparam_ij], &params[iparam_ji],fac11,
                delx,dely,delz,rsq,mr1,mr2,mr3,
                sr1,sr2,sr3,iq,jq,i,j,eng_tmp,fvionij,ddprx);
        vionij = eng_tmp;

        // sums up dipole energies
        if (evflag)
          ev_tally(i,j,nlocal,newton_pair,0.0,vionij,fvionij,delx,dely,delz);

        // sums up dipole forces
        f[i][0] += (ddprx[0] + delx*fvionij);
        f[i][1] += (ddprx[1] + dely*fvionij);
        f[i][2] += (ddprx[2] + delz*fvionij);
        f[j][0] -= (ddprx[0] + delx*fvionij);
        f[j][1] -= (ddprx[1] + dely*fvionij);
        f[j][2] -= (ddprx[2] + delz*fvionij);
      }

      if (rsq > params[iparam_ij].cutsq) continue;

      repulsive(&params[iparam_ij], &params[iparam_ji], rsq,
              fpair, eflag, eng_tmp, iq, jq);

      evdwl = eng_tmp;

      // repulsion is pure two-body, sums up pair repulsive forces
      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag)
        ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
    }

    // many-body interactions: start of short-range
    xcn = NCo[i];
    for (jj = 0; jj < sht_jnum; jj++) {
      j = sht_jlist[jj];
      sht_llist = sht_first[j];
      sht_lnum = sht_num[j];

      jtag = tag[j];
      if (jtag <= itag) continue ;
      ycn = NCo[j];

      jtype = map[type[j]];
      iparam_ij = elem3param[itype][jtype][jtype];
      iparam_ji = elem3param[jtype][itype][itype];

      delrj[0] = x[j][0] - xtmp;
      delrj[1] = x[j][1] - ytmp;
      delrj[2] = x[j][2] - ztmp;
      rsq1 = dot3(delrj,delrj);
      if (rsq1 > params[iparam_ij].cutsq) continue;
      nj ++;

      // this Qj for q-dependent BSi
      jq = q[j];

      // accumulate bondorder zeta for each i-j interaction via k and l loops
      zeta_ij = 0.0;
      bbtor = 0.0;
      kconjug = 0.0;

      for (kk = 0; kk < sht_jnum; kk++) {       // kk is neighbor of ii
        k = sht_jlist[kk];
        if (j == k) continue;

        ktype = map[type[k]];
        iparam_ijk = elem3param[itype][jtype][ktype];
        iparam_ikj = elem3param[itype][ktype][jtype];
        iparam_jik = elem3param[jtype][itype][ktype];
        iparam_ik  = elem3param[itype][ktype][ktype];
        delrk[0] = x[k][0] - xtmp;
        delrk[1] = x[k][1] - ytmp;
        delrk[2] = x[k][2] - ztmp;
        rsq2 = dot3(delrk,delrk);

        if (rsq2 > params[iparam_ik].cutsq) continue;

        // 3-body zeta in bond order
        zeta_ij += zeta(&params[iparam_ijk], &params[iparam_ik],
                        rsq1, rsq2, delrj, delrk, i, xcn);

        // radical initialization: apply only to CC,CO,OC bonds
        if (params[iparam_ij].rad_flag > 0 &&
            params[iparam_ik].ielementgp == 1 &&
            params[iparam_ik].jelementgp == 1) {
          iparam_ki = elem3param[ktype][itype][itype];
          kcn=NCo[k];
          kconjug += rad_init(rsq2,&params[iparam_ki],i,kradtot,kcn);

        }

        // torsion: i-j-k-l: apply to all C-C bonds

        if (params[iparam_ij].tor_flag != 0) {
          srmu = dot3(delrj,delrk)/(sqrt(rsq1*rsq2));
          srmu = sqrt(1.0-srmu*srmu);

          if (srmu > 0.1) {
            for (ll = 0; ll < sht_lnum; ll++) { // ll is neighbor of jj
              l = sht_llist[ll];

              if (l==i || l==j || l==k) continue;

              ltype = map[type[l]];

              delrl[0] = x[l][0] - x[j][0];
              delrl[1] = x[l][1] - x[j][1];
              delrl[2] = x[l][2] - x[j][2];
              rsq3 = dot3(delrl,delrl);
              iparam_jl = elem3param[jtype][ltype][ltype];

              if (rsq3 > params[iparam_jl].cutsq) continue;

              iparam_ikl = elem3param[itype][ktype][ltype];
              torindx = params[iparam_ij].tor_flag;
              bbtor += bbtor1(torindx, &params[iparam_ikl],&params[iparam_jl],
                       rsq1,rsq2,rsq3,delrj,delrk,delrl,srmu);
            }
          }
        }
      }

      zeta_ji = 0.0;
      lconjug = 0.0;

      for (ll = 0; ll < sht_lnum; ll++) {
        l = sht_llist[ll];
        if (l == i) continue;

        ltype = map[type[l]];
        iparam_jil = elem3param[jtype][itype][ltype];
        iparam_ijl = elem3param[itype][jtype][ltype];
        iparam_jl  = elem3param[jtype][ltype][ltype];
        iparam_lj  = elem3param[ltype][jtype][jtype];

        delrk[0] = x[l][0] - x[j][0];
        delrk[1] = x[l][1] - x[j][1];
        delrk[2] = x[l][2] - x[j][2];
        rsq2 = dot3(delrk,delrk);

        delrl[0] = x[l][0] - x[j][0];
        delrl[1] = x[l][1] - x[j][1];
        delrl[2] = x[l][2] - x[j][2];
        rsq2 = dot3(delrl,delrl);

        if (rsq2 > params[iparam_jl].cutsq) continue;

        scale3(-1,delrj,delrl);     // ji_hat is -(ij_hat)

        zeta_ji += zeta(&params[iparam_jil], &params[iparam_jl]
                        , rsq1, rsq2, delrl, delrk, j, ycn);

        // radical initialization: apply only to CC,CO,OC bonds
        if (params[iparam_ji].rad_flag > 0
          && params[iparam_jl].ielementgp == 1
          && params[iparam_jl].jelementgp == 1) {
          iparam_lj = elem3param[ltype][jtype][jtype];
          lcn=NCo[l];
          lconjug += rad_init(rsq2,&params[iparam_lj],j,lradtot,lcn);
        }
      }

      force_zeta(&params[iparam_ij], &params[iparam_ji],
         rsq1, xcn, ycn, zeta_ij, zeta_ji, fpair,
         prefac_ij1, prefac_ij2, prefac_ij3, prefac_ij4, prefac_ij5,
         prefac_ji1, prefac_ji2, prefac_ji3, prefac_ji4, prefac_ji5,
         eflag, eng_tmp, iq, jq, i, j, nj, bbtor, kconjug, lconjug);

      evdwl = eng_tmp;
      selfp6p(&params[iparam_ij],&params[iparam_ji],rsq1,eng_tmp,fpair);

      evdwl += eng_tmp;
      f[i][0] += delrj[0]*fpair;
      f[i][1] += delrj[1]*fpair;
      f[i][2] += delrj[2]*fpair;
      f[j][0] -= delrj[0]*fpair;
      f[j][1] -= delrj[1]*fpair;
      f[j][2] -= delrj[2]*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,-fpair,-delrj[0],-delrj[1],-delrj[2]);

      // attractive term via loop over k (3-body forces: i-j-k)
      zet_addi=0;
      zet_addj=0;

      for (kk = 0; kk < sht_jnum; kk++) {
        k = sht_jlist[kk];
        if (j == k) continue;
        sht_mlist = sht_first[k];
        sht_mnum = sht_num[k];

        ktype = map[type[k]];
        iparam_ijk = elem3param[itype][jtype][ktype];
        iparam_ikj = elem3param[itype][ktype][jtype];
        iparam_jik = elem3param[jtype][itype][ktype];
        iparam_ik  = elem3param[itype][ktype][ktype];
        delrk[0] = x[k][0] - xtmp;
        delrk[1] = x[k][1] - ytmp;
        delrk[2] = x[k][2] - ztmp;
        rsq2 = dot3(delrk,delrk);
        if (rsq2 > params[iparam_ik].cutsq) continue;

        // BO-dependent 3-body E & F
        attractive(&params[iparam_ijk], &params[iparam_jik],&params[iparam_ikj],
                prefac_ij1, prefac_ij2, prefac_ij3, prefac_ij4, prefac_ij5,
                rsq1,rsq2,delrj,delrk,fi,fj,fk,i,xcn);

        ep6p_ij = ep6p(&params[iparam_ijk],&params[iparam_ikj],rsq1,rsq2,delrj,delrk,zet_addi);
        fp6p(&params[iparam_ijk],&params[iparam_ikj],rsq1,rsq2,delrj,delrk,fip6p,fjp6p,fkp6p);

        // Sums up i-j-k forces: LP contribution
        for (im = 0; im < 3; im++) {
          fi[im] += fip6p[im];
          fj[im] += fjp6p[im];
          fk[im] += fkp6p[im];
        }

        // Sums up i-j-k forces: Tallies into global force vector
        for (im = 0; im < 3; im++) {
          f[i][im] += fi[im];
          f[j][im] += fj[im];
          f[k][im] += fk[im];
        }

        // torsion and radical: apply to all C-C bonds
        if (params[iparam_ijk].tor_flag != 0 && fabs(ptorr)>1.0e-8) {
          srmu = dot3(delrj,delrk)/(sqrt(rsq1*rsq2));
          srmu = sqrt(1.0-srmu*srmu);

          if (srmu > 0.1) {
            for (ll = 0; ll < sht_lnum; ll++) { // ll is neighbor of jj
              l = sht_llist[ll];
              if (l==i||l==j||l==k) continue;

              ltype = map[type[l]];

              delrl[0] = x[l][0] - x[j][0];
              delrl[1] = x[l][1] - x[j][1];
              delrl[2] = x[l][2] - x[j][2];
              rsq3 = dot3(delrl,delrl);

              iparam_jl = elem3param[jtype][ltype][ltype];
              if (rsq3 > params[iparam_jl].cutsq) continue;
              iparam_ikl = elem3param[itype][ktype][ltype];
              torindx = params[iparam_ij].tor_flag;
              tor_force(torindx, &params[iparam_ikl], &params[iparam_jl],srmu,
                          rsq1,rsq2,rsq3,delrj,delrk,delrl);

              for (im = 0; im < 3; im++) {
                f[i][im] += fi_tor[im];
                f[j][im] += fj_tor[im];
                f[k][im] += fk_tor[im];
                f[l][im] += fl_tor[im];
              }
            }
          }
        }

        if ( params[iparam_ijk].rad_flag>=1 &&
          params[iparam_ijk].ielementgp==1 &&
          params[iparam_ijk].kelementgp==1) {
          iparam_ki = elem3param[ktype][itype][itype];
          kcn=NCo[k];
          double rik=sqrt(rsq2);
          kradtot = -comb_fc(rik,&params[iparam_ki])*params[iparam_ki].pcross+kcn;

          rad_forceik(&params[iparam_ki],rsq2,delrk,kconjug,kradtot);

          for (im = 0; im < 3; im++) {
            f[i][im] += fi_rad[im];
            f[k][im] += fk_rad[im];
          }

          if (fabs(radtmp) > 1.0e-12) {
            for (mm = 0; mm < sht_mnum; mm++) { // mm is neighbor of kk
               m = sht_mlist[mm];
               if (m == k) continue;

               mtype = map[type[m]];

              delrm[0] = x[m][0] - x[k][0];
              delrm[1] = x[m][1] - x[k][1];
              delrm[2] = x[m][2] - x[k][2];
              rsq3 = dot3(delrm,delrm);

              iparam_km = elem3param[ktype][mtype][mtype];
              iparam_ki = elem3param[ktype][itype][itype];

              if (rsq3 > params[iparam_km].cutsq) continue;

              rad_force(&params[iparam_km],rsq3,delrm,radtmp);

              for (im = 0; im < 3; im++) {
                f[k][im] += fj_rad[im];
                f[m][im] += fk_rad[im];
              }
            }
          }
        }

        if (evflag)
          ev_tally(i,j,nlocal,newton_pair,ep6p_ij,0.0,0.0,0.0,0.0,0.0);
        if (vflag_either)
          v_tally3(i,j,k,fj,fk,delrj,delrk);

      } // k-loop

      // attractive term via loop over l (3-body forces: j-i-l)
      for (ll = 0; ll < sht_lnum; ll++) {
        l = sht_llist[ll];
        if (l == i) continue;

        sht_plist = sht_first[l];
        sht_pnum = sht_num[l];

        ltype = map[type[l]];
        iparam_jil = elem3param[jtype][itype][ltype];
        iparam_jli = elem3param[jtype][ltype][itype];
        iparam_ijl = elem3param[itype][jtype][ltype];
        iparam_jl  = elem3param[jtype][ltype][ltype];
        delrk[0] = x[l][0] - x[j][0];
        delrk[1] = x[l][1] - x[j][1];
        delrk[2] = x[l][2] - x[j][2];

        rsq2 = dot3(delrk,delrk);
        if (rsq2 > params[iparam_jl].cutsq) continue;
        scale3(-1,delrj,delrl);

        attractive(&params[iparam_jil],&params[iparam_ijl],&params[iparam_jli],
                prefac_ji1,prefac_ji2,prefac_ji3,prefac_ji4,prefac_ji5,
                rsq1,rsq2,delrl,delrk,fj,fi,fl,j,ycn);

        // BO-independent 3-body j-i-l LP and BB correction and forces
        ep6p_ji = ep6p(&params[iparam_jil],&params[iparam_jli],rsq1,rsq2,delrl,delrk,zet_addj);
        fp6p(&params[iparam_jil],&params[iparam_jli],rsq1,rsq2,delrl,delrk,fjp6p,fip6p,flp6p);

        if (evflag)
          ev_tally(j,i,nlocal,newton_pair,ep6p_ji,0.0,0.0,0.0,0.0,0.0);

        // BO-dependent 3-body E & F
        for (im = 0; im < 3; im++) {
          fj[im] += fjp6p[im];
          fi[im] += fip6p[im];
          fl[im] += flp6p[im];
        }

        // Sums up j-i-l forces: Tallies into global force vector
        for (im = 0; im < 3; im++) {
          f[j][im] += fj[im];
          f[i][im] += fi[im];
          f[l][im] += fl[im];
        }

        // radical i-j-l-p: apply to all CC,CO,OC bonds
        if ( params[iparam_jil].rad_flag >= 1 &&
          params[iparam_jil].ielementgp == 1 &&
          params[iparam_jil].kelementgp == 1) {
            iparam_lj = elem3param[ltype][jtype][jtype];
            lcn=NCo[l];
            double rjl=sqrt(rsq2);
            lradtot=-comb_fc(rjl,&params[iparam_lj])*params[iparam_lj].pcross +lcn;

            rad_forceik(&params[iparam_lj],rsq2,delrk,lconjug,lradtot);

            for (im = 0; im < 3; im++) {
              f[j][im] += fi_rad[im];
              f[l][im] += fk_rad[im];
            }

            if (fabs(radtmp)>1.0e-12) {
              for (pp = 0; pp < sht_pnum; pp++) {       // pp is neighbor of ll
                p = sht_plist[pp];
                if (p == l) continue;
                ptype = map[type[p]];

                delrp[0] = x[p][0] - x[l][0];
                delrp[1] = x[p][1] - x[l][1];
                delrp[2] = x[p][2] - x[l][2];
                rsq3 = dot3(delrp,delrp);

                iparam_lp = elem3param[ltype][ptype][ptype];

                if (rsq3 > params[iparam_lp].cutsq) continue;

                scale3(-1,delrj);
                rad_force(&params[iparam_lp],rsq3,delrp,radtmp);
                scale3(-1,delrj);
                for (im = 0; im < 3; im++) {
                  f[l][im] += fj_rad[im];
                  f[p][im] += fk_rad[im];
                }
              }
            }
        }

        if (vflag_either)
          v_tally3(j,i,l,fi,fl,delrl,delrk);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

}

/* ---------------------------------------------------------------------- */

void PairComb3::repulsive(Param *parami, Param *paramj, double rsq,
        double &fforce,int /*eflag*/, double &eng, double iq, double jq)
{
  double r,tmp_fc,tmp_fc_d,Di,Dj;
  double caj,vrcs,fvrcs;
  double LamDiLamDj,fcdA,rlm1,bigA;

  double romi = parami->addrep;
  double rrcs = parami->bigr + parami->bigd;
  double addr = parami->addrepr;

  r = sqrt(rsq);
  if (r > rrcs) return ;

  tmp_fc = comb_fc(r,parami);
  tmp_fc_d = comb_fc_d(r,parami);

  Di = parami->DU + pow(fabs(parami->bD*(parami->QU-iq)),parami->nD);
  Dj = paramj->DU + pow(fabs(paramj->bD*(paramj->QU-jq)),paramj->nD);

  bigA = parami->bigA;
  rlm1 = parami->lambda;

  fcdA = tmp_fc_d - tmp_fc * rlm1;
  LamDiLamDj = exp(0.5*(parami->lami*Di+paramj->lami*Dj)-rlm1*r);
  caj = bigA * LamDiLamDj;

  fforce = -caj * fcdA;

  // additional repulsion

  vrcs = 1.0; fvrcs = 0.0;
  if (romi != 0.0 && r < addr) {
    vrcs += romi * square(1.0-r/addr);
    fvrcs = romi * 2.0 * (r/addr-1.0)/addr;
    fforce = fforce*vrcs - caj * tmp_fc * vrcs * fvrcs;
  }
  fforce /= r;

  // eng = repulsive energy
  eng = caj * tmp_fc * vrcs;
}

/* ---------------------------------------------------------------------- */

double PairComb3::zeta(Param *parami, Param *paramj, double rsqij,
        double rsqik, double *delrij, double *delrik, int /*i*/, double xcn)
{
  double rij,rik,costheta,arg,ex_delr,rlm3;

  rij = sqrt(rsqij);
  if (rij > parami->bigr+parami->bigd) return 0.0;
  rik = sqrt(rsqik);
  costheta = dot3(delrij,delrik) / (rij*rik);

  rlm3 = parami->beta;
  arg = powint(rlm3*(rij-rik),int(parami->powermint));
  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return comb_fc(rik,paramj) * comb_gijk(costheta,parami,xcn) * ex_delr;
}

/* ---------------------------------------------------------------------- */

void PairComb3::selfp6p(Param *parami, Param *paramj, double rsq,
        double &eng, double &force)
{
  double r,comtti,comttj,fcj,fcj_d;

  r=sqrt(rsq);
  fcj=comb_fc(r,parami);
  fcj_d=comb_fc_d(r,parami);
  comtti = comttj = 0.0;

    double pilp0 = parami->p6p0;
    double pilp1 = parami->p6p1, pilp2 = parami->p6p2, pilp3 = parami->p6p3;
    double pilp4 = parami->p6p4, pilp5 = parami->p6p5, pilp6 = parami->p6p6;
    comtti = pilp0 + pilp1 + pilp2 + pilp3 + pilp4 + pilp5 + pilp6;

    double pjlp0 = paramj->p6p0;
    double pjlp1 = paramj->p6p1, pjlp2 = paramj->p6p2, pjlp3 = paramj->p6p3;
    double pjlp4 = paramj->p6p4, pjlp5 = paramj->p6p5, pjlp6 = paramj->p6p6;
    comttj = pjlp0 + pjlp1 + pjlp2 + pjlp3 + pjlp4 + pjlp5 + pjlp6;

  eng = 0.5 * fcj * (comtti + comttj);
  force += 0.5 * fcj_d * (comtti + comttj)/r;
}

/* ---------------------------------------------------------------------- */

double PairComb3::ep6p(Param *paramj, Param *paramk, double rsqij, double rsqik,
                     double *delrij, double *delrik , double &/*zet_add*/)
{
  double comtt;
  double pplp0 = paramj->p6p0;
  double pplp1 = paramj->p6p1, pplp2 = paramj->p6p2, pplp3 = paramj->p6p3;
  double pplp4 = paramj->p6p4, pplp5 = paramj->p6p5, pplp6 = paramj->p6p6;
  double rij,rik,costheta,lp0,lp1,lp2,lp3,lp4,lp5,lp6;
  double rmu,rmu2,rmu3,rmu4,rmu5,rmu6,fcj,fck;
    comtt=0.0;
    rij = sqrt(rsqij);
    rik = sqrt(rsqik);
    costheta = dot3(delrij,delrik) / (rij*rik);
    fcj = comb_fc(rij,paramj);
    fck = comb_fc(rik,paramk);

    rmu = costheta;
    rmu2 = rmu *rmu; rmu3 = rmu2*rmu; rmu4 = rmu3*rmu;
    rmu5 = rmu4*rmu; rmu6 = rmu5*rmu;
    lp0 = pplp0;
    lp1 = pplp1*rmu;
    lp2 = pplp2*rmu2;
    lp3 = pplp3*rmu3;
    lp4 = pplp4*rmu4;
    lp5 = pplp5*rmu5;
    lp6 = pplp6*rmu6;
    comtt = lp0 + lp1 + lp2 + lp3 + lp4 + lp5 + lp6;
    return 0.5 * fck * comtt *fcj;
  }

/*---------------------------------------------------------------------- */

void PairComb3::fp6p(Param *paramij,Param *paramik, double rsqij, double rsqik,
                   double *delrij, double *delrik, double *drilp,
                   double *drjlp, double *drklp)
{
  double pplp0 = paramij->p6p0;
  double pplp1 = paramij->p6p1, pplp2 = paramij->p6p2, pplp3 = paramij->p6p3;
  double pplp4 = paramij->p6p4, pplp5 = paramij->p6p5, pplp6 = paramij->p6p6;
  double ffj1,ffj2,ffk1,ffk2;
  double rij,rik,costheta;
  double rmu,comtt,comtt_d,com4k,com5,com5k,fcj,fcj_d,fck,fck_d;
  double lp0,lp1,lp2,lp3,lp4,lp5,lp6;
  double lp1_d,lp2_d,lp3_d,lp4_d,lp5_d,lp6_d;
  double rmu2, rmu3, rmu4, rmu5, rmu6;

    ffj1 = 0.0, ffj2 = 0.0;
    ffk1 = 0.0, ffk2 = 0.0;

    rij = sqrt(rsqij); rik = sqrt(rsqik);
    costheta = dot3(delrij,delrik) / (rij*rik);
    fcj = comb_fc(rij,paramij);
    fck = comb_fc(rik,paramik);
    fcj_d = comb_fc_d(rij,paramij);
    fck_d = comb_fc_d(rik,paramik);
    rmu = costheta;

      rmu2 = rmu *rmu; rmu3 = rmu2*rmu;
      rmu4 = rmu3*rmu; rmu5 = rmu4*rmu; rmu6 = rmu5*rmu;
      lp0 = pplp0;
      lp1 = pplp1*rmu;
      lp2 = pplp2*rmu2;
      lp3 = pplp3*rmu3;
      lp4 = pplp4*rmu4;
      lp5 = pplp5*rmu5;
      lp6 = pplp6*rmu6;
      lp1_d = pplp1;
      lp2_d = pplp2*2.0*rmu;
      lp3_d = pplp3*3.0*rmu2;
      lp4_d = pplp4*4.0*rmu3;
      lp5_d = pplp5*5.0*rmu4;
      lp6_d = pplp6*6.0*rmu5;
      comtt = lp0 + lp1 + lp2 + lp3 + lp4 + lp5 + lp6;
      comtt_d = lp1_d + lp2_d + lp3_d + lp4_d + lp5_d + lp6_d;

      com4k = fcj * fck_d * comtt;
      com5  = fcj * fck * comtt_d;
      com5k = fck * comtt * fcj_d;

      ffj1 = 0.5*(-com5/(rij*rik));
      ffj2 = 0.5*(com5*rmu/rsqij-com5k/rij);
      ffk1 = ffj1;
      ffk2 = 0.5*(-com4k/rik+com5*rmu/rsqik);

  // j-atom
  scale3(ffj1,delrik,drjlp);
  scaleadd3(ffj2,delrij,drjlp,drjlp);

  // k-atom
  scale3(ffk1,delrij,drklp);
  scaleadd3(ffk2,delrik,drklp,drklp);

  // i-atom
  add3(drjlp,drklp,drilp);
  scale3(-1.0,drilp);
}

/* ---------------------------------------------------------------------- */

void PairComb3::force_zeta(Param *parami, Param *paramj, double rsq,
        double xcn, double ycn, double &zeta_ij, double &zeta_ji, double &fforce,
        double &prefac_ij1, double &prefac_ij2, double &prefac_ij3,
        double &prefac_ij4, double &prefac_ij5,
        double &prefac_ji1, double &prefac_ji2, double &prefac_ji3,
        double &prefac_ji4, double &prefac_ji5,
        int eflag, double &eng, double iq, double jq,
        int i, int j, int nj, double bbtor, double kconjug, double lconjug)
{
  double r,att_eng,att_force,bij;  // att_eng is -cbj
  double boij, dbij1, dbij2, dbij3, dbij4, dbij5;
  double boji, dbji1, dbji2, dbji3, dbji4, dbji5;
  double pradx, prady;
  r = sqrt(rsq);

  if (r > parami->bigr + parami->bigd) return;
  comb_fa(r, parami, paramj, iq, jq, att_eng, att_force);
  comb_bij_d(zeta_ij,parami,r,i,boij,dbij1,dbij2,dbij3,dbij4,dbij5,xcn);
  comb_bij_d(zeta_ji,paramj,r,j,boji,dbji1,dbji2,dbji3,dbji4,dbji5,ycn);
  bij = 0.5*(boij + boji);

  // radical energy

  if (parami->rad_flag>0) {
    rad_calc( r, parami, paramj, kconjug, lconjug, i, j, xcn, ycn);
    bij +=  brad[0];
    pradx = brad[1]*att_eng;
    prady = brad[2]*att_eng;
    brad[3] = 1.0 * brad[3]*att_eng;
   }

  // torsion energy
  if ( parami->tor_flag!=0) {
     tor_calc( r, parami, paramj, kconjug, lconjug, i, j, xcn, ycn);
     bij += btor[0] * bbtor;
     ptorr =  att_eng * btor[0];
     pradx  += 1.0 *  btor[1] * bbtor * att_eng;
     prady  += 1.0 *  btor[2] * bbtor * att_eng;
     brad[3]+= 1.0 *  btor[3] * bbtor * att_eng;
  }

  fforce = 1.0*bij*att_force/r; // divide by r will done compute
  bbij[i][nj] = bij;

  prefac_ij1 = -0.5*att_eng*dbij1;      // prefac_ij1 = -pfij
  prefac_ij2 = -0.5*att_eng*dbij2;      // prefac_ij2 = -pfij1
  prefac_ij3 = -0.5*att_eng*dbij3;      // prefac_ij3 = -pfij2
  prefac_ij4 = -0.5*att_eng*dbij4;      // prefac_ij4 = -pfij3
  prefac_ij5 = -0.5*att_eng*dbij5;      // prefac_ij5 = -pfij4

  prefac_ji1 = -0.5*att_eng*dbji1;      // prefac_ji1 = -pfji
  prefac_ji2 = -0.5*att_eng*dbji2;      // prefac_ji2 = -pfji1
  prefac_ji3 = -0.5*att_eng*dbji3;      // prefac_ji3 = -pfji2
  prefac_ji4 = -0.5*att_eng*dbji4;      // prefac_ji4 = -pfji3
  prefac_ji5 = -0.5*att_eng*dbji5;      // prefac_ji5 = -pfji4

  // combines com6 & com7 below
  if (parami->rad_flag>0 || parami->tor_flag!=0) {
    prefac_ij2-=pradx;
    prefac_ji2-=prady;
  }

  // eng = attraction energy
  if (eflag) eng = 1.0*bij*att_eng;
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_fc(double r, Param *param)
{
  double r_inn = param->bigr - param->bigd;
  double r_out = param->bigr + param->bigd;
  if (r <= r_inn) return 1.0;
  if (r >= r_out) return 0.0;
  return 0.5*(1.0 + cos(MY_PI*(r-r_inn)/(r_out-r_inn)));
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_fc_d(double r, Param *param)
{
  double r_inn = param->bigr - param->bigd;
  double r_out = param->bigr + param->bigd;
  if (r <= r_inn) return 0.0;
  if (r >= r_out) return 0.0;
  return -MY_PI2/(r_out-r_inn)*sin(MY_PI*(r-r_inn)/(r_out-r_inn));
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_fccc(double xcn)
{
  double cut1 = ccutoff[0];
  double cut2 = ccutoff[1];

  if (xcn <= cut1) return 1.0;
  if (xcn >= cut2) return 0.0;
  return 0.5*(1.0 + cos(MY_PI*(xcn-cut1)/(cut2-cut1)));
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_fccc_d(double xcn)
{
  double cut1 = ccutoff[0];
  double cut2 = ccutoff[1];

  if (xcn <= cut1) return 0.0;
  if (xcn >= cut2) return 0.0;
  return -MY_PI2/(cut2-cut1)*sin(MY_PI*(xcn-cut1)/(cut2-cut1));
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_fcch(double xcn)
{
  double cut1 = ccutoff[2];
  double cut2 = ccutoff[3];

  if (xcn <= cut1) return 1.0;
  if (xcn >= cut2) return 0.0;
  return 0.5*(1.0 + cos(MY_PI*(xcn-cut1)/(cut2-cut1)));
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_fcch_d(double xcn)
{
  double cut1 = ccutoff[2];
  double cut2 = ccutoff[3];

  if (xcn <= cut1) return 0.0;
  if (xcn >= cut2) return 0.0;
  return -MY_PI2/(cut2-cut1)*sin(MY_PI*(xcn-cut1)/(cut2-cut1));
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_fccch(double xcn)
{
  double cut1 = ccutoff[4];
  double cut2 = ccutoff[5];

  if (xcn <= cut1) return 1.0;
  if (xcn >= cut2) return 0.0;
  return 0.5*(1.0 + cos(MY_PI*(xcn-cut1)/(cut2-cut1)));
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_fccch_d(double xcn)
{
  double cut1 = ccutoff[4];
  double cut2 = ccutoff[5];

  if (xcn <= cut1) return 0.0;
  if (xcn >= cut2) return 0.0;
  return -MY_PI2/(cut2-cut1)*sin(MY_PI*(xcn-cut1)/(cut2-cut1));
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_fcsw(double rsq)
{
  double r = sqrt(rsq);

  if (r <= chicut1) return 1.0;
  if (r >= chicut2) return 0.0;
  return 0.5*(1.0 + cos(MY_PI*(r-chicut1)/(chicut2-chicut1)));
}

/* ---------------------------------------------------------------------- */

double PairComb3::self(Param *param, double qi)
{
 double self_tmp, cmin, cmax, qmin, qmax;
 double s1=param->chi, s2=param->dj, s3=param->dk, s4=param->dl;

 self_tmp = 0.0;

 qmin = param->qmin;
 qmax = param->qmax;
 cmin = cmax = 100.0;

 self_tmp = qi*(s1+qi*(s2+qi*(s3+qi*s4)));

 if (qi < qmin) self_tmp += cmin * powint((qi-qmin),4);
 if (qi > qmax) self_tmp += cmax * powint((qi-qmax),4);

 return self_tmp;
}

/* ---------------------------------------------------------------------- */

void PairComb3::comb_fa(double r, Param *parami, Param *paramj, double iq,
        double jq, double &att_eng, double &att_force)
{
  double Bsi;
  double qi,qj,Di,Dj;
  double AlfDiAlfDj, YYBn, YYBj;
  double alfij1= parami->alpha1;
  double alfij2= parami->alpha2;
  double alfij3= parami->alpha3;
  double pbij1= parami->bigB1;
  double pbij2= parami->bigB2;
  double pbij3= parami->bigB3;
  if (r > parami->bigr + parami->bigd) Bsi = 0.0;

  qi = iq; qj = jq;
  Di = Dj = Bsi = 0.0;
  Di = parami->DU + pow(fabs(parami->bD*(parami->QU-qi)),parami->nD);
  Dj = paramj->DU + pow(fabs(paramj->bD*(paramj->QU-qj)),paramj->nD);
  YYBn = (parami->aB-fabs(powint(parami->bB*(qi-parami->Qo),10)));
  YYBj = (paramj->aB-fabs(powint(paramj->bB*(qj-paramj->Qo),10)));

  if (YYBn*YYBj > 0.0) {
    AlfDiAlfDj = exp(0.5*(parami->alfi*Di+paramj->alfi*Dj));
    Bsi = (pbij1*exp(-alfij1*r)+pbij2*exp(-alfij2*r)+pbij3*exp(-alfij3*r))*
      sqrt(YYBn*YYBj)*AlfDiAlfDj;                               // Bsi is cbj

    att_eng = -Bsi * comb_fc(r,parami);
    att_force = -(Bsi*comb_fc_d(r,parami)-comb_fc(r,parami)*sqrt(YYBn*YYBj)*
        AlfDiAlfDj*(alfij1*pbij1*exp(-alfij1*r)+
        alfij2*pbij2*exp(-alfij2*r)+alfij3*pbij3*exp(-alfij3*r)));

  } else {
    att_eng = 0.0;
    att_force = 0.0;
  }

}

/* ---------------------------------------------------------------------- */

void PairComb3::comb_bij_d(double zet, Param *param, double r, int i,
        double &tbij, double &tbij1, double &tbij2,
        double &tbij3, double &tbij4, double &tbij5, double xcn)
{
  double pcorn,dpcorn,dxccij,dxchij,dxcoij;
  double zeta = zet;
  double zetang,tmp_tbij, pow_n;

  pcorn = dpcorn = dxccij = dxchij = dxcoij = 0.0;
  coord(param,r,i,pcorn,dpcorn,dxccij,dxchij,dxcoij,xcn);       // coordination term

  zetang=zeta;
  pow_n=param->powern;
  zeta = pow(zetang,pow_n)+pcorn;
  tmp_tbij=pow_n*pow(zetang,(pow_n-1.0));

  if ((1.0 + zeta) < 0.1) {
    zeta=0.1-1.0;
    tbij = pow(1.0 + zeta, -0.5/pow_n);
    tbij1=0.0;
   }
   else if (zeta > param->c1) {
    tbij = pow(zeta,-0.5/pow_n);
    tbij1 = -0.5/pow_n*pow(zeta,(-0.5/pow_n-1.0));
   } else if (zeta > param->c2) {
    tbij = pow(zeta,-0.5/pow_n)-0.5/pow_n*pow(zeta,(-0.5/pow_n-1.0));
    tbij1 = -0.5/pow_n/zeta;
   } else if (fabs(zeta) < param->c4) {
    tbij = 1.0;
    tbij1 = 0.0;
   } else if (fabs(zeta) < param->c3) {
    tbij = 1.0 - zeta/(2.0*pow_n);
    tbij1 = -1/(2.0*pow_n);
   } else {
    tbij = pow(1.0 + zeta, -0.5/pow_n);
    tbij1 = -0.5/pow_n * pow(1.0 + zeta,(-1.0-0.5/pow_n));
   }

  tbij2 = tbij1 * dpcorn;
  tbij3 = tbij1 * dxccij;
  tbij4 = tbij1 * dxchij;
  tbij5 = tbij1 * dxcoij;
  tbij1 = tbij1 * tmp_tbij;

}

/* ---------------------------------------------------------------------- */

void PairComb3::coord(Param *param, double r, int i,
        double &pcorn, double &dpcorn, double &dxccij,
        double &dxchij, double &dxcoij, double xcn)
{
  int ixmin,iymin,izmin;
  double xcntot,xcccn,xchcn,xcocn;
  int tri_flag= param-> pcn_flag;
  int jele_gp= param->jelementgp;
  double pan = param->pcna;
  double pbn = param->pcnb;
  double pcn = param->pcnc;
  double pdn = param->pcnd;

  xcccn = xchcn = xcocn = 0.0;

  xcccn = xcctmp[i];
  xchcn = xchtmp[i];
  xcocn = xcotmp[i];
  xcntot = -comb_fc(r,param)*param->pcross + xcn;
  pcorn = dpcorn = dxccij = dxchij = dxcoij = 0.0;
  pcorn = 0.0; dpcorn = 0.0;

  if (xcntot  < 0.0) xcntot  = 0.0;

  if (tri_flag>0) {
    if (jele_gp==1) xcccn = xcccn-comb_fc(r,param)*param->pcross;
    if (jele_gp==2) xchcn = xchcn-comb_fc(r,param)*param->pcross;
    if (jele_gp==3) xcocn = xcocn-comb_fc(r,param)*param->pcross;
    if (xcccn < 0.0) xcccn = 0.0;
    if (xchcn < 0.0) xchcn = 0.0;
    if (xcocn < 0.0) xcocn = 0.0;
    if (xcccn > maxx) xcccn = maxx;
    if (xchcn > maxy) xchcn = maxy;
    if (xcocn > maxz) xcocn = maxz;

    double xcntritot=xcccn+xchcn+xcocn;

    if (xcntritot > maxxcn[tri_flag-1]) {
      pcorn  = vmaxxcn[tri_flag-1]+(xcntot-maxxcn[tri_flag-1])*dvmaxxcn[tri_flag-1];
      dxccij = dxchij = dxcoij = dvmaxxcn[tri_flag-1];
    }
    else {
      ixmin=int(xcccn+1.0e-12);
      iymin=int(xchcn+1.0e-12);
      izmin=int(xcocn+1.0e-12);
      if (fabs(float(ixmin)-xcccn)>1.0e-8 ||
          fabs(float(iymin)-xchcn)>1.0e-8 ||
          fabs(float(izmin)-xcocn)>1.0e-8) {
            cntri_int(tri_flag,xcccn,xchcn,xcocn,ixmin,iymin,izmin,
            pcorn,dxccij,dxchij,dxcoij,param);
      }
      else  {
        pcorn  = pcn_grid[tri_flag-1][ixmin][iymin][izmin];
        dxccij = pcn_gridx[tri_flag-1][ixmin][iymin][izmin];
        dxchij = pcn_gridy[tri_flag-1][ixmin][iymin][izmin];
        dxcoij = pcn_gridz[tri_flag-1][ixmin][iymin][izmin];
      }
    }
  } else {
    pcorn = pan*xcntot+pbn*exp(pcn*xcntot)+pdn;
    dpcorn = pan+pbn*pcn*exp(pcn*xcntot);
  }
}

/* ---------------------------------------------------------------------- */

void PairComb3::cntri_int(int tri_flag, double xval, double yval,
                double zval, int ixmin, int iymin, int izmin, double &vval,
                double &dvalx, double &dvaly, double &dvalz, Param * /*param*/)
{
  double x;
  vval = 0.0; dvalx = 0.0; dvaly = 0.0; dvalz = 0.0;
  if (ixmin >= maxx-1) { ixmin=maxx-1; }
  if (iymin >= maxy-1) { iymin=maxy-1; }
  if (izmin >= maxz-1) { izmin=maxz-1; }
  for (int j=0; j<64; j++) {
      x = pcn_cubs[tri_flag-1][ixmin][iymin][izmin][j]
          *pow(xval,iin3[j][0])*pow(yval,iin3[j][1])
          *pow(zval,iin3[j][2]);
    vval += x;
    if (xval>1.0e-8) {dvalx += x*iin3[j][0]/xval;}
    if (yval>1.0e-8) {dvaly += x*iin3[j][1]/yval;}
    if (zval>1.0e-8) {dvalz += x*iin3[j][2]/zval;}
  }
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_gijk(double costheta, Param *param, double nco_tmp)
{
  double rmu1 = costheta;
  double rmu2 = rmu1*rmu1;
  double rmu3 = rmu2*rmu1;
  double rmu4 = rmu3*rmu1;
  double rmu5 = rmu4*rmu1;
  double rmu6 = rmu5*rmu1;
  double co6 = param->pcos6*rmu6;
  double co5 = param->pcos5*rmu5;
  double co4 = param->pcos4*rmu4;
  double co3 = param->pcos3*rmu3;
  double co2 = param->pcos2*rmu2;
  double co1 = param->pcos1*rmu1;
  double co0 = param->pcos0;
  double pcross = param->pcross;
  double gmu;

  if (param->ang_flag==1) {
    double qtheta, gmu1, gmu2, rrmu, astep;
    int k;

    qtheta = comb_fccc(nco_tmp);
    astep = 2.0/ntab;
    rrmu = (rmu1+1.0)/astep;
    k = int(rrmu);
    gmu1 = co6+co5+co4+co3+co2+co1+co0;
    gmu2 = pang[k]+(pang[k+1]-pang[k])*(rrmu-k);
    gmu = gmu2+qtheta*(gmu1-gmu2);
    return gmu*pcross;

  } else if (param->ang_flag==2) {
    double qtheta, gmu1, gmu2;
    double ch6 = ch_a[6]*rmu6;
    double ch5 = ch_a[5]*rmu5;
    double ch4 = ch_a[4]*rmu4;
    double ch3 = ch_a[3]*rmu3;
    double ch2 = ch_a[2]*rmu2;
    double ch1 = ch_a[1]*rmu1;
    double ch0 = ch_a[0];
    qtheta = comb_fccch(nco_tmp);
    gmu1 = co6+co5+co4+co3+co2+co1+co0;
    gmu2 = ch6+ch5+ch4+ch3+ch2+ch1+ch0;
    gmu = gmu2+qtheta*(gmu1-gmu2);
    return gmu*pcross;
  } else {
    gmu = co6+co5+co4+co3+co2+co1+co0;
    return gmu*pcross;
  }
}

/* ---------------------------------------------------------------------- */

void PairComb3::comb_gijk_d(double costheta, Param *param, double nco_tmp,
                double &gijk_d, double &com3jk)
{
  double rmu1 = costheta;
  double rmu2 = rmu1*rmu1;
  double rmu3 = rmu2*rmu1;
  double rmu4 = rmu3*rmu1;
  double rmu5 = rmu4*rmu1;
  double rmu6 = rmu5*rmu1;
  double co6 = param->pcos6; //*rmu5*6.0;
  double co5 = param->pcos5; //*rmu4*5.0;
  double co4 = param->pcos4; //*rmu3*4.0;
  double co3 = param->pcos3; //*rmu2*3.0;
  double co2 = param->pcos2; //*rmu1*2.0;
  double co1 = param->pcos1;
  double co0 = param->pcos0;
  double pcross = param->pcross;

  gijk_d = com3jk = 0.0;
  if (param->ang_flag==1) {
    double qtheta, dqtheta, gmu1, gmu2, dgmu1,dgmu2, rrmu, astep;
    int k;
    qtheta = comb_fccc(nco_tmp);
    dqtheta = comb_fccc_d(nco_tmp);

    astep = 2.0/ntab;
    rrmu = (rmu1+1.0)/astep;
    k = int(rrmu);

    gmu1 =rmu6*co6+rmu5*co5+rmu4*co4
         +rmu3*co3+rmu2*co2+rmu1*co1+co0;
    dgmu1 =6.0*rmu5*co6+5.0*rmu4*co5+4.0*rmu3*co4
           +3.0*rmu2*co3+2.0*rmu1*co2+co1;
    gmu2 = pang[k]+(pang[k+1]-pang[k])*(rrmu-k);
    dgmu2 = dpang[k]+(dpang[k+1]-dpang[k])*(rrmu-k);
    gijk_d = pcross*(dgmu2+qtheta*(dgmu1-dgmu2));
    com3jk = dqtheta * (gmu1-gmu2);
  } else if (param->ang_flag==2) {
    double qtheta, dqtheta, gmu1, gmu2, dgmu1,dgmu2;
    double ch6 = ch_a[6];
    double ch5 = ch_a[5];
    double ch4 = ch_a[4];
    double ch3 = ch_a[3];
    double ch2 = ch_a[2];
    double ch1 = ch_a[1];
    double ch0 = ch_a[0];
    qtheta = comb_fccch(nco_tmp);
    dqtheta = comb_fccch_d(nco_tmp);

    gmu1 =rmu6*co6+rmu5*co5+rmu4*co4
         +rmu3*co3+rmu2*co2+rmu1*co1+co0;
    dgmu1 =6.0*rmu5*co6+5.0*rmu4*co5+4.0*rmu3*co4
           +3.0*rmu2*co3+2.0*rmu1*co2+co1;
    gmu2 =rmu6*ch6+rmu5*ch5+rmu4*ch4
         +rmu3*ch3+rmu2*ch2+rmu1*ch1+ch0;
    dgmu2 =6.0*rmu5*ch6+5.0*rmu4*ch5+4.0*rmu3*ch4
           +3.0*rmu2*ch3+2.0*rmu1*ch2+ch1;
    gijk_d = pcross*(dgmu2+qtheta*(dgmu1-dgmu2));
    com3jk = dqtheta * (gmu1-gmu2);

  } else {
    gijk_d = pcross*(6.0*rmu5*co6+5.0*rmu4*co5+4.0*rmu3*co4
                    +3.0*rmu2*co3+2.0*rmu1*co2+co1);
    com3jk = 0.0;
  }
}

/*------------------------------------------------------------------------- */

void PairComb3::attractive(Param *parami, Param *paramj , Param *paramk, double prefac_ij1,
        double prefac_ij2, double prefac_ij3, double prefac_ij4,
        double prefac_ij5, double rsqij, double rsqik, double *delrij,
        double *delrik, double *fi, double *fj,double *fk, int /*i*/, double xcn)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  scale3(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  scale3(rikinv,delrik,rik_hat);

  comb_zetaterm_d(prefac_ij1, prefac_ij2, prefac_ij3, prefac_ij4, prefac_ij5,
        rij_hat, rij,rik_hat, rik, fi, fj, fk, parami, paramj, paramk,xcn);

}

/* ---------------------------------------------------------------------- */

void PairComb3::comb_zetaterm_d(double prefac_ij1, double prefac_ij2,
        double prefac_ij3, double prefac_ij4, double prefac_ij5,
        double *rij_hat, double rij, double *rik_hat, double rik, double *dri,
        double *drj, double *drk, Param *parami, Param *paramj, Param *paramk, double xcn)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc_k,cos_theta,tmp,rlm3;
  double dcosdri[3],dcosdrj[3],dcosdrk[3],dfc_i,dfc_k;
  double com6, com3j, com3k, com3jk;

  int mint = int(parami->powermint);
  double pcrossi = parami->pcross;
  double pcrossj = paramj->pcross;
  double pcrossk = paramk->pcross;
  int icontrol = parami->pcn_flag;

  dfc_i = comb_fc_d(rij,parami);
  fc_k = comb_fc(rik,paramk);
  dfc_k = comb_fc_d(rik,paramk);
  rlm3 = parami->beta;
  tmp = powint(rlm3*(rij-rik),mint);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);
  ex_delr *= pcrossi;

  cos_theta = dot3(rij_hat,rik_hat);
  gijk = comb_gijk(cos_theta,parami,xcn);
  comb_gijk_d(cos_theta,parami,xcn,gijk_d,com3jk);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // com6 & com7
  if (icontrol > 0) {
    if (parami->kelementgp==1) {com6 = prefac_ij3*pcrossk*dfc_k;}
    if (parami->kelementgp==2) {com6 = prefac_ij4*pcrossk*dfc_k;}
    if (parami->kelementgp==3) {com6 = prefac_ij5*pcrossk*dfc_k;}
    if (parami->rad_flag>=1 || parami->tor_flag!=0)
            {com6+=prefac_ij2*pcrossk*dfc_k;}
  } else {
    com6 = prefac_ij2*pcrossi*dfc_k;
  }

  if (parami->ang_flag==1 || parami->ang_flag==2) {
    com3j = com3jk*ex_delr*pcrossk*pcrossj*fc_k*dfc_i;
    com3k = com3jk*ex_delr*pcrossk*pcrossk*fc_k*dfc_k;
  } else {
    com3j = 0.0;
    com3k = 0.0;
  }

  ex_delr_d = mint*powint(rlm3,mint)*powint((rij-rik),(mint-1))*ex_delr; // com3
  scale3(-dfc_k*gijk*ex_delr,rik_hat,dri);          // com1
  scaleadd3(fc_k*gijk_d*ex_delr,dcosdri,dri,dri);   // com2
  scaleadd3(fc_k*gijk*ex_delr_d,rik_hat,dri,dri);   // com3 cont'd
  scaleadd3(-fc_k*gijk*ex_delr_d,rij_hat,dri,dri);  // com3 sums j
  scaleadd3(-com3k,rik_hat,dri,dri);                // com3k
  scaleadd3(-com3j,rij_hat,dri,dri);                // com3j
  scale3(prefac_ij1,dri);
  scaleadd3(-com6,rik_hat,dri,dri);                 // com6

  scale3(fc_k*gijk_d*ex_delr,dcosdrj,drj);          // com2
  scaleadd3(fc_k*gijk*ex_delr_d,rij_hat,drj,drj);   // com3 cont'd
  scaleadd3(com3j,rij_hat,drj,drj);                 // com3j
  scale3(prefac_ij1,drj);

  scale3(dfc_k*gijk*ex_delr,rik_hat,drk);           // com1
  scaleadd3(fc_k*gijk_d*ex_delr,dcosdrk,drk,drk);   // com2
  scaleadd3(-fc_k*gijk*ex_delr_d,rik_hat,drk,drk);  // com3 cont'd
  scaleadd3(com3k,rik_hat,drk,drk);                 // com3k
  scale3(prefac_ij1,drk);
  scaleadd3(com6,rik_hat,drk,drk);                  // com6
}

/* ---------------------------------------------------------------------- */

void PairComb3::costheta_d(double *rij_hat, double rij, double *rik_hat,
        double rik, double *dri, double *drj, double *drk)
{
  double cos_theta = dot3(rij_hat,rik_hat);

  scaleadd3(-cos_theta,rij_hat,rik_hat,drj);
  scale3(1.0/rij,drj);
  scaleadd3(-cos_theta,rik_hat,rij_hat,drk);
  scale3(1.0/rik,drk);
  add3(drj,drk,dri);
  scale3(-1.0,dri);
}

/* ---------------------------------------------------------------------- */

void PairComb3::tables()

{
  int i,j,k,m, nntypes, ncoul,nnbuf, ncoul_lim, inty, itype, jtype;
  int iparam_i, iparam_ij, iparam_ji;
  double r,dra,drin,drbuf,rc,z,zr,zrc,ea,eb,ea3,eb3,alf;
  double exp2er,exp2ersh,fafash,dfafash,F1,dF1,ddF1,E1,E2,E3,E4;
  double exp2ear,exp2ebr,exp2earsh,exp2ebrsh,fafbsh,dfafbsh;
  double afbshift, dafbshift, exp2ershift;

  int n = nelements;

  memory->destroy(intype);
  memory->destroy(erpaw);
  memory->destroy(fafb);
  memory->destroy(dfafb);
  memory->destroy(ddfafb);
  memory->destroy(phin);
  memory->destroy(dphin);
  memory->destroy(afb);
  memory->destroy(dafb);
  memory->destroy(vvdw);
  memory->destroy(vdvdw);
  memory->destroy(dpl);
  memory->destroy(bbij);
  memory->destroy(xcctmp);
  memory->destroy(xchtmp);
  memory->destroy(xcotmp);
  memory->destroy(NCo);
  memory->destroy(sht_num);
  memory->sfree(sht_first);

  dra  = 0.001;
  drin = 0.100;
  drbuf = 0.100;
  nnbuf = int(drbuf/dra) +1;
  rc = cutmax;
  alf = 0.20;
  nmax = atom->nmax;

  nntypes = int((n+1)*n/2.0)+1;
  ncoul = int((rc-drin)/dra)+ nnbuf;
  ncoul_lim = int(ncoul * 1.20);

  // allocate arrays
  memory->create(intype,n,n,"pair:intype");
  memory->create(erpaw,ncoul_lim,3,"pair:erpaw");
  memory->create(fafb,ncoul_lim,nntypes,"pair:fafb");
  memory->create(dfafb,ncoul_lim,nntypes,"pair:dfafb");
  memory->create(ddfafb,ncoul_lim,nntypes,"pair:ddfafb");
  memory->create(phin,ncoul_lim,nntypes,"pair:phin");
  memory->create(dphin,ncoul_lim,nntypes,"pair:dphin");
  memory->create(afb,ncoul_lim,nntypes,"pair:afb");
  memory->create(dafb,ncoul_lim,nntypes,"pair:dafb");
  memory->create(vvdw,ncoul,nntypes,"pair:vvdw");
  memory->create(vdvdw,ncoul,nntypes,"pair:vdvdw");
  memory->create(dpl,nmax,3,"pair:dpl");
  memory->create(bbij,nmax,MAXNEIGH,"pair:bbij");
  memory->create(xcctmp,nmax,"pair:xcctmp");
  memory->create(xchtmp,nmax,"pair:xchtmp");
  memory->create(xcotmp,nmax,"pair:xcotmp");
  memory->create(NCo,nmax,"pair:NCo");
  memory->create(sht_num,nmax,"pair:sht_num");
  sht_first = (int **) memory->smalloc(nmax*sizeof(int *),
        "pair:sht_first");

  // set interaction number: 0-0=0, 1-1=1, 0-1=1-0=2

  m = 0; k = n;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j == i) {
        intype[i][j] = m;
        m += 1;
      } else if (j != i && j > i) {
        intype[i][j] = k;
        k += 1;
      } else if (j != i && j < i) {
        intype[i][j] = intype[j][i];
      }
    }
  }

  // default arrays to zero

  for (i = 0; i < ncoul; i ++) {
    for (j = 0; j < nntypes; j ++) {
      fafb[i][j]   = 0.0;
      dfafb[i][j]  = 0.0;
      ddfafb[i][j] = 0.0;
      phin[i][j]   = 0.0;
      dphin[i][j]  = 0.0;
      afb[i][j]    = 0.0;
      dafb[i][j]   = 0.0;
    }
  }

  // direct 1/r energy with Slater 1S orbital overlap

  for (i = 0; i < n; i++) {
    if (map[i+1] < 0) continue;
    r = drin - dra;
    itype = map[i+1];
    iparam_i = elem3param[itype][itype][itype];
    z = params[iparam_i].esm;
    exp2ershift = exp(-2.0*z*rc);
    afbshift = -exp2ershift*(z+1.0/rc);
    dafbshift = exp2ershift*(2.0*z*z+2.0*z/rc+1.0/(rc*rc));

    if (comm->me == 0 && screen)
      fprintf(screen,"  element[%d] = %-2s, z = %g\n",i+1,elements[map[i+1]],z);

    for (j = 0; j < ncoul; j++) {
      exp2er = exp(-2.0 * z * r);
      phin[j][i] = 1.0 - exp2er * (1.0 + 2.0 * z * r * (1.0 + z * r));
      dphin[j][i] = (4.0 * exp2er * z * z * z * r * r);
      afb[j][i] = -exp2er*(z+1.0/r)-afbshift-(r-rc)*dafbshift;
      dafb[j][i] = -(exp2er*(2.0*z*z+2.0*z/r+1.0/(r*r))-dafbshift);
      r += dra;
    }
  }

  for (i = 0; i < n; i ++) {
    if (map[i+1] < 0) continue;
    for (j = 0; j < n; j ++) {
      if (map[j+1] < 0) continue;
      r = drin - dra;
      if (j == i) {
        itype = map[i+1];
        inty = intype[itype][itype];
        iparam_i = elem3param[itype][itype][itype];
        z = params[iparam_i].esm;
        zrc = z * rc;
        exp2ersh = exp(-2.0 * zrc);
        fafash = -exp2ersh * (1.0 / rc +
                              z * (11.0/8.0 + 3.0/4.0*zrc + zrc*zrc/6.0));
        dfafash = exp2ersh * (1.0/(rc*rc) + 2.0*z/rc +
                              z*z*(2.0 + 7.0/6.0*zrc + zrc*zrc/3.0));
        for (k = 0; k < ncoul; k ++) {
          zr = z * r;
          exp2er = exp(-2.0*zr);
          F1 = -exp2er * (1.0 / r +
                          z * (11.0/8.0 + 3.0/4.0*zr + zr*zr/6.0));
          dF1 = exp2er * (1.0/(r*r) + 2.0*z/r +
                          z*z*(2.0 + 7.0/6.0*zr + zr*zr/3.0));
          ddF1 = -exp2er * (2.0/(r*r*r) + 4.0*z/(r*r) + 4.0*z*z/r +
                            z*z*z/3.0*(17.0/2.0 + 5.0*zr + 2.0*zr*zr));
          fafb[k][inty] = F1-fafash-(r-rc)*dfafash;
          dfafb[k][inty] = -(dF1 - dfafash);
          ddfafb[k][inty] = ddF1;
                  r += dra;
        }
      } else if (j != i) {
        itype = map[i+1];
        jtype = map[j+1];
        inty = intype[itype][jtype];
        iparam_ij = elem3param[itype][jtype][jtype];
        ea = params[iparam_ij].esm;
        ea3 = ea*ea*ea;
        iparam_ji = elem3param[jtype][itype][itype];
        eb = params[iparam_ji].esm;
        eb3 = eb*eb*eb;
        E1 = ea*eb3*eb/((ea+eb)*(ea+eb)*(ea-eb)*(ea-eb));
        E2 = eb*ea3*ea/((ea+eb)*(ea+eb)*(eb-ea)*(eb-ea));
        E3 = (3.0*ea*ea*eb3*eb-eb3*eb3) /
          ((ea+eb)*(ea+eb)*(ea+eb)*(ea-eb)*(ea-eb)*(ea-eb));
        E4 = (3.0*eb*eb*ea3*ea-ea3*ea3) /
          ((ea+eb)*(ea+eb)*(ea+eb)*(eb-ea)*(eb-ea)*(eb-ea));
        exp2earsh = exp(-2.0*ea*rc);
        exp2ebrsh = exp(-2.0*eb*rc);
        fafbsh = -exp2earsh*(E1 + E3/rc)-exp2ebrsh*(E2 + E4/rc);
        dfafbsh =
          exp2earsh*(2.0*ea*(E1+E3/rc)+E3/(rc*rc)) +
          exp2ebrsh*(2.0*eb*(E2+E4/rc)+E4/(rc*rc));
        for (k = 0; k < ncoul; k ++) {
          exp2ear = exp(-2.0*ea*r);
          exp2ebr = exp(-2.0*eb*r);
          fafb[k][inty] =
            - exp2ear*(E1+E3/r) - exp2ebr*(E2+E4/r)
            - fafbsh - (r-rc) * dfafbsh;
          dfafb[k][inty] = -(exp2ear*(2.0*ea*(E1+E3/r) + E3/(r*r))
                           + exp2ebr*(2.0*eb*(E2+E4/r) + E4/(r*r))- dfafbsh);
          ddfafb[k][inty] = -exp2ear*(4.0*ea*ea*(E1+E3/r)+4.0*ea*E3/(r*r)
                                +2.0*E3/(r*r*r))
                            -exp2ebr*(4.0*eb*eb*(E2+E4/r)+4.0*eb*E4/(r*r)
                                +2.0*E4/(r*r*r));
          r += dra;
        }
      }
    }
  }

  for (i = 0; i < ncoul_lim; i ++) {
    r = dra * (i-1) + drin;
    erpaw[i][0] = erfc(r*alf);
    erpaw[i][1] = exp(-r*r*alf*alf);
  }
  // end wolf summation

  // van der Waals
  int ii,jj;
  double **rvdw, *cc2, *cc3, *vrc, *rrc;
  double r6, r7, r12, r13, rf6, rf12, drf7, drf13;
  double drcc, temp6, temp7, temp12, temp13;
  double vsigt, vepst, vdwt, dvdwt;

  vrc = new double[13];
  rrc = new double[13];
  cc2 = new double[nntypes];
  cc3 = new double[nntypes];
  memory->create(rvdw,2,nntypes,"pair:rvdw");

  vrc[0] = rc;
  for (i=1; i<13; i++) {
    vrc[i] = vrc[i-1] * vrc[0];
  }

  // generate spline coefficients for CC, CH, HH vdw
  for (ii = 0; ii < n; ii ++) {
    for (jj = ii; jj < n; jj ++) {
      itype = ii;
      jtype = jj;
      inty = intype[itype][jtype];
      iparam_ij = elem3param[itype][jtype][jtype];

      // parameter check: eps > 0
      if (params[iparam_ij].vdwflag > 0) {

        if (params[iparam_ij].vdwflag==1) {
          rvdw[0][inty] = params[iparam_ij].bigr + params[iparam_ij].bigd;
        }
        else {
          rvdw[0][inty] = params[iparam_ij].bigr - params[iparam_ij].bigd;
        }

        rvdw[1][inty] = params[iparam_ij].vsig * 0.950;

        // radius check: outer radius vs. sigma
        if (rvdw[0][inty] > rvdw[1][inty])
          error->all(FLERR,"Error in vdw spline: inner radius > outer radius");

        rrc[0] = rvdw[1][inty];

        for (i=1; i<13; i++)
          rrc[i] = rrc[i-1] * rrc[0];

        drcc = rrc[0] - rvdw[0][inty];
        temp6 = 1.0/rrc[5]-1.0/vrc[5]+6.0*(rrc[0]-vrc[0])/vrc[6];
        temp7 = 6.0*(1.0/vrc[6]-1.0/rrc[6]);
        temp12 = 1.0/rrc[11]-1.0/vrc[11]+(rrc[0]-vrc[0])*12.0/vrc[12];
        temp13 = 12.0*(1.0/vrc[12]-1.0/rrc[12]);

        vsigt = params[iparam_ij].vsig;
        vepst = params[iparam_ij].veps;
        vsigt = vsigt*vsigt*vsigt*vsigt*vsigt*vsigt;

        vdwt = vepst*(vsigt*vsigt*temp12-vsigt*temp6);
        dvdwt = vepst*(vsigt*vsigt*temp13-vsigt*temp7);
        cc2[inty] = (3.0/drcc*vdwt-dvdwt)/drcc;
        cc3[inty] = (vdwt/(drcc*drcc)-cc2[inty] )/drcc;
      }
    }
  }

  // generate vdw look-up table
  for (ii = 0; ii < n; ii ++) {
    for (jj = ii; jj < n; jj ++) {
      itype = ii;
      jtype = jj;
      inty = intype[itype][jtype];
      iparam_ij = elem3param[itype][jtype][jtype];
      r = drin;
      for (k = 0; k < ncoul; k ++) {
        r6 = r*r*r*r*r*r;
        r7 = r6 * r;
        rf6 = 1.0/r6-1.0/vrc[5]+(r-vrc[0])*6.0/vrc[6];
        drf7 = 6.0*(1.0/vrc[6]-1.0/r7);
        vsigt = params[iparam_ij].vsig;
        vepst = params[iparam_ij].veps;
        vsigt = vsigt*vsigt*vsigt*vsigt*vsigt*vsigt;

        if (params[iparam_ij].vdwflag>0) {
          if (r <= rvdw[0][inty]) {
            vvdw[k][inty] = 0.0;
            vdvdw[k][inty] = 0.0;
          }
          else if ( r > rvdw[0][inty] && r <= rvdw[1][inty]) {
            drcc = r-rvdw[0][inty];
            vvdw[k][inty] = drcc*drcc*(drcc*cc3[inty]+cc2[inty]);
            vdvdw[k][inty] = drcc*(3.0*drcc*cc3[inty]+2.0*cc2[inty]);
          } else {
            r12 = r6*r6;
            r13 = r6*r7;
            rf12 = 1.0/r12-1.0/vrc[11]+(r-vrc[0])*12.0/vrc[12];
            drf13= 12.0*(1.0/vrc[12]-1.0/r13);
            vvdw[k][inty] = vepst*(vsigt*vsigt*rf12-vsigt*rf6);
            vdvdw[k][inty] = vepst*(vsigt*vsigt*drf13-vsigt*drf7);
          }
        } else {
          vvdw[k][inty]=0.0;
          vdvdw[k][inty]=0.0;
        }
          r += dra;
      }
    }
  }

  delete[] vrc;
  delete[] rrc;
  delete[] cc2;
  delete[] cc3;
  memory->destroy(rvdw);
}

/* ---------------------------------------------------------------------- */

void PairComb3::potal_calc(double &calc1, double &calc2, double &calc3)
{
  double alf,rcoul,esucon;
  int m;

  rcoul = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].lcut > rcoul) rcoul = params[m].lcut;

  alf = 0.20;
  esucon = force->qqr2e;

  calc2 = (erfc(rcoul*alf)/rcoul/rcoul+2.0*alf/MY_PIS*
           exp(-alf*alf*rcoul*rcoul)/rcoul)*esucon/rcoul;
  calc3 = (erfc(rcoul*alf)/rcoul)*esucon;
  calc1 = -(alf/MY_PIS*esucon+calc3*0.5);
}

/* ---------------------------------------------------------------------- */

void PairComb3::tri_point(double rsq, int &mr1, int &mr2,
                 int &mr3, double &sr1, double &sr2, double &sr3)
{
  double r, rin, dr, dd, rr1, rridr, rridr2;

  rin = 0.1000; dr = 0.0010;
  r = sqrt(rsq);
  if (r < rin + 2.0*dr)    r = rin + 2.0*dr;
  if (r > cutmax - 2.0*dr) r = cutmax - 2.0*dr;
  rridr = (r-rin)/dr;

  mr1 = int(rridr) ;
  dd = rridr - float(mr1);
  if (dd > 0.5) mr1 += 1;

  rr1 = float(mr1)*dr;
  rridr = (r - rin - rr1)/dr;
  rridr2 = rridr * rridr;

  sr1 = (rridr2 - rridr) * 0.50;
  sr2 = 1.0 - rridr2;
  sr3 = (rridr2 + rridr) * 0.50;

  mr2 = mr1 + 1;
  mr3 = mr1 + 2;
}

/* ---------------------------------------------------------------------- */

void PairComb3::vdwaals(int inty, int mr1, int mr2, int mr3, double rsq,
                      double sr1, double sr2, double sr3,
                      double &eng, double &fforce)
{
  double r = sqrt(rsq);

  eng = 1.0*(sr1*vvdw[mr1-1][inty]+sr2*vvdw[mr2-1][inty]+sr3*vvdw[mr3-1][inty]);
  fforce = -1.0/r*(sr1*vdvdw[mr1-1][inty]+sr2*vdvdw[mr2-1][inty]+sr3*vdvdw[mr3-1][inty]);
}

/* ---------------------------------------------------------------------- */

void PairComb3::direct(Param *parami, Param *paramj, int mr1,
        int mr2, int mr3, double rsq, double sr1, double sr2, double sr3,
        double iq, double jq, double fac11, double fac11e,
        double &pot_tmp, double &for_tmp, int i, int j)
{
  double r,erfcc,fafbnl,potij,esucon;
  double r3,erfcd,dfafbnl,smf2,dvdrr,alf,alfdpi;
  double afbn,afbj,sme1n,sme1j,sme1,sme2,dafbn, dafbj,smf1n,smf1j;
  double curli = parami->curl;
  double curlj = paramj->curl;
  int inti = parami->ielement;
  int intj = paramj->ielement;
  int inty = intype[inti][intj];

  double curlij0 = parami->curl0;
  double curlji0 = paramj->curl0;
  double curlij1,curlji1,dcurlij,dcurlji;
  double fcp1j,xcoij,xcoji;
  int icurl, jcurl;
  int ielegp = parami->ielementgp;
  int jelegp = paramj->ielementgp;

  r = sqrt(rsq);
  r3 = r * rsq;
  alf = 0.20;
  alfdpi = 2.0*alf/MY_PIS;
  esucon = force->qqr2e;
  pot_tmp = for_tmp = 0.0;
  icurl=jcurl=0;

  if (ielegp==2 && curli>curlij0) {
    icurl=1;
    curlij1=curli;
  }

  if (jelegp==2 && curlj>curlji0) {
    jcurl=1;
    curlji1=curlj;
  }

  if (icurl==1 || jcurl ==1) {
    xcoij = xcotmp[i];
    xcoji = xcotmp[j];
    fcp1j = comb_fc_d(r,parami);

    if (icurl==1) {
      curli=curlij1+(curlij0-curlij1)*comb_fc_curl(xcoij,parami);
      dcurlij=fcp1j*(curlij0-curlij1)*comb_fc_curl_d(xcoij,parami);
    }

    if (jcurl==1) {
      curlj=curlji1+(curlji0-curlji1)*comb_fc_curl(xcoji,paramj);
      dcurlji=fcp1j*(curlji0-curlji1)*comb_fc_curl_d(xcoji,paramj);
    }
  }

  erfcc = sr1*erpaw[mr1][0] + sr2*erpaw[mr2][0] + sr3*erpaw[mr3][0];
  afbn = sr1*afb[mr1][inti] + sr2*afb[mr2][inti] + sr3*afb[mr3][inti];
  afbj = sr1*afb[mr1][intj] + sr2*afb[mr2][intj] + sr3*afb[mr3][intj];
  fafbnl= sr1*fafb[mr1][inty] + sr2*fafb[mr2][inty] + sr3*fafb[mr3][inty];
  potij = (erfcc/r * esucon - fac11e);

  sme1n = iq*curlj*(afbn-fafbnl)*esucon;
  sme1j = jq*curli*(afbj-fafbnl)*esucon;
  sme1 = sme1n + sme1j;
  sme2 = (potij + fafbnl * esucon) * iq * jq;
  pot_tmp = 1.0 * (sme1+sme2);

  // 1/r force (wrt r)

  erfcd = sr1*erpaw[mr1][1] + sr2*erpaw[mr2][1] + sr3*erpaw[mr3][1];
  dafbn = sr1*dafb[mr1][inti] + sr2*dafb[mr2][inti] + sr3*dafb[mr3][inti];
  dafbj = sr1*dafb[mr1][intj] + sr2*dafb[mr2][intj] + sr3*dafb[mr3][intj];
  dfafbnl= sr1*dfafb[mr1][inty] + sr2*dfafb[mr2][inty] + sr3*dfafb[mr3][inty];

  dvdrr = (erfcc/r3+alfdpi*erfcd/rsq)*esucon-fac11;
  smf1n = iq * curlj * (dafbn-dfafbnl)*esucon/r;
  smf1j = jq * curli * (dafbj-dfafbnl)*esucon/r;

  if (jcurl==1 && ielegp == 3 && dcurlji != 0.0) {
   smf1n += dcurlji*iq*(afbn-fafbnl)*esucon/r;
  }
  if (icurl==1 && jelegp == 3 && dcurlij != 0.0) {
   smf1j += dcurlij*jq*(afbj-fafbnl)*esucon/r;
  }

  smf2 = dvdrr + dfafbnl * esucon/r;
  for_tmp =  1.0 * iq * jq * smf2 + smf1n + smf1j;
}

/* ---------------------------------------------------------------------- */

void PairComb3::field(Param *parami, Param *paramj, double rsq, double iq,
                double jq, double &eng_tmp,double &for_tmp)
{
  double r,r3,r4,r5,rc,rc2,rc3,rc4,rc5;
  double cmi1,cmi2,cmj1,cmj2,pcmi1,pcmi2;
  double rf3i,rcf3i,rf5i,rcf5i;
  double drf3i,drcf3i,drf5i,drcf5i;
  double rf3,rf5,drf4,drf6;
  double smpn,smpl,rfx1,rfx2;

  r = sqrt(rsq);
  r3 = r * r * r;
  r4 = r3 * r;
  r5 = r4 * r;
  rc = parami->lcut;
  rc2 = rc * rc;
  rc3 = rc*rc*rc;
  rc4 = rc3 * rc;
  rc5 = rc4 * rc;
  cmi1 = parami->cmn1;
  cmi2 = parami->cmn2;
  cmj1 = paramj->cmn1;
  cmj2 = paramj->cmn2;
  pcmi1 = parami->pcmn1;
  pcmi2 = parami->pcmn2;

  rf3i = r3/(square(r3)+cube(pcmi1));
  rcf3i = rc3/(square(rc3)+cube(pcmi1));
  rf5i = r5/(square(r5)+powint(pcmi2,5));
  rcf5i = rc5/(square(rc5)+powint(pcmi2,5));

  drf3i = 3/r*rf3i-6*rsq*rf3i*rf3i;
  drcf3i = 3/rc*rcf3i-6*rc2*rcf3i*rcf3i;
  drf5i = 5/r*rf5i-10*r4*rf5i*rf5i;
  drcf5i = 5/rc*rcf5i-10*rc4*rcf5i*rcf5i;

  rf3 = rf3i-rcf3i-(r-rc)*drcf3i;
  rf5 = rf5i-rcf5i-(r-rc)*drcf5i;
  drf4 = drf3i - drcf3i;
  drf6 = drf5i - drcf5i;

 // field correction energy
  smpn = jq*(cmi1*rf3+jq*cmi2*rf5);
  smpl = iq*(cmj1*rf3+iq*cmj2*rf5);
  eng_tmp = 1.0 * (smpn + smpl);

 // field correction force
  rfx1 = jq*(cmi1*drf4+jq*cmi2*drf6)/r;
  rfx2 = iq*(cmj1*drf4+iq*cmj2*drf6)/r;
  for_tmp -= 1.0 * (rfx1 + rfx2);
}

/* ---------------------------------------------------------------------- */

double PairComb3::rad_init(double rsq2,Param *param,int /*i*/,
                double &radtot, double cnconj)
{
  double r, fc1k, radcut;

  r = sqrt(rsq2);
  fc1k = comb_fc(r,param);
  radtot = -fc1k * param->pcross + cnconj;
  radcut = comb_fcch(radtot);
  return fc1k * param->pcross * radcut;
}

/* ---------------------------------------------------------------------- */

void PairComb3::rad_calc(double r, Param *parami, Param *paramj,
        double kconjug, double lconjug, int /*i*/, int /*j*/, double xcn, double ycn)
{
  int ixmin, iymin, izmin;
  int radindx;
  double xrad, yrad, zcon, vrad, pradx, prady, pradz;

  vrad = pradx = prady = pradz = 0.0;
  xrad = -comb_fc(r,parami)*parami->pcross + xcn;
  yrad = -comb_fc(r,paramj)*paramj->pcross + ycn;
  zcon = 1.0 + square(kconjug) + square(lconjug);

  if (xrad < 0.0) xrad = 0.0;
  if (yrad < 0.0) yrad = 0.0;
  if (zcon < 1.0) zcon = 1.0;
  if (xrad > maxxc) xrad = maxxc;
  if (yrad > maxyc) yrad = maxyc;
  if (zcon > maxconj) zcon = maxconj;
  ixmin = int(xrad+1.0e-12);
  iymin = int(yrad+1.0e-12);
  izmin = int(zcon+1.0e-12);
  radindx=parami->rad_flag-1;
  if (fabs(float(ixmin)-xrad)>1.0e-8 ||
      fabs(float(iymin)-yrad)>1.0e-8 ||
      fabs(float(izmin)-zcon)>1.0e-8) {
    rad_int(radindx,xrad,yrad,zcon,ixmin,iymin,izmin,
              vrad,pradx,prady,pradz);
  } else {
    vrad  = rad_grid[radindx][ixmin][iymin][izmin-1];
    pradx = rad_gridx[radindx][ixmin][iymin][izmin-1];
    prady = rad_gridy[radindx][ixmin][iymin][izmin-1];
    pradz = rad_gridz[radindx][ixmin][iymin][izmin-1];
  }

  brad[0] = vrad;
  brad[1] = pradx;
  brad[2] = prady;
  brad[3] = pradz;
}

/* ---------------------------------------------------------------------- */

void PairComb3::rad_int(int radindx,double xrad, double yrad, double zcon, int l,
                int m, int n, double &vrad, double &pradx, double &prady,
                double &pradz)
{
  int j;
  double x;
  vrad = pradx = prady = pradz = 0.0;
  if (l >= maxxc-1) { l=maxxc-1;}
  if (m >= maxyc-1) { m=maxyc-1; }
  if (n >= maxconj-1) { n=maxconj-1;}

  for (j=0; j<64; j++) {
    x = rad_spl[radindx][l][m][n-1][j] * pow(xrad,iin3[j][0])
          * pow(yrad,iin3[j][1]) * pow(zcon,iin3[j][2]);
    vrad  += x;
    if (xrad > 1.0e-8) pradx += x*iin3[j][0]/xrad;
    if (yrad > 1.0e-8) prady += x*iin3[j][1]/yrad;
    if (zcon > 1.0e-8) pradz += x*iin3[j][2]/zcon;
  }
}


/* ---------------------------------------------------------------------- */
void PairComb3::rad_forceik(Param *paramk, double rsq2, double *delrk,
        double conjug, double radtot)
{
  int nm;
  double  rik, fc1k, fcp1k;
  double pradk, ffkk2, fktmp[3];
  double radcut = comb_fcch(radtot);
  double dradcut = comb_fcch_d(radtot);

  for (nm=0; nm<3; nm++) {
    fi_rad[nm] =  fk_rad[nm] = 0.0;
  }
    radtmp =0.0;

  rik = sqrt(rsq2);

  fc1k = comb_fc(rik, paramk);
  fcp1k = comb_fc_d(rik,paramk);

  pradk = brad[3]*fcp1k*radcut*paramk->pcross*2.0*conjug;
  radtmp= brad[3]*fc1k*dradcut*paramk->pcross*2.0*conjug;

  ffkk2 = -pradk/rik;

  for (nm=0; nm<3; nm++) {
    fktmp[nm] = - ffkk2 * delrk[nm];
  }

  for (nm=0; nm<3; nm++) {
    fi_rad[nm] =  fktmp[nm];
    fk_rad[nm] = -fktmp[nm];
  }
}

/* ---------------------------------------------------------------------- */

void PairComb3::rad_force(Param *paramm, double rsq3,
        double *delrm, double dpradk)
{
  int nm;
  double rkm, fcp1m;
  double comkm, ffmm2, fkm[3];

  for (nm=0; nm<3; nm++) {
    fj_rad[nm] = fk_rad[nm] = 0.0;
    fkm[nm]=0.0;
  }

  rkm = sqrt(rsq3);

  fcp1m = comb_fc_d(rkm, paramm);

  comkm = dpradk * fcp1m *  paramm->pcross;
  ffmm2 = -comkm/rkm;

  for (nm=0; nm<3; nm++) {
    fkm[nm] = -ffmm2 * delrm[nm];
  }

  for (nm=0; nm<3; nm++) {
    fj_rad[nm] =  fkm[nm];
    fk_rad[nm] = -fkm[nm];
  }
}

/* ---------------------------------------------------------------------- */

double PairComb3::bbtor1(int torindx, Param *paramk, Param *paraml,
        double rsq1, double rsq2, double rsq3, double *delrj,
        double *delrk, double *delrl, double srmu)
{
  double rmul, rij, rik, rjl;

  rij = sqrt(rsq1);
  rik = sqrt(rsq2);
  rjl = sqrt(rsq3);

  scale3(-1.0,delrl,delrl);
  rmul = dot3(delrj,delrl)/(rij*rjl);
  scale3(-1.0,delrl,delrl);
  rmul = sqrt(1.0-rmul*rmul);

  if (rmul > 0.1) {
    double fc1k, fc1l, TT1, TT2, rmut, btt, tork[3], torl[3];

    fc1k = comb_fc(rik,paramk);
    fc1l = comb_fc(rjl,paraml);

    TT1 = rik*rjl*rij*rij*srmu*rmul;
    tork[0] = delrk[1]*delrj[2] - delrk[2]*delrj[1];
    torl[0] = delrj[1]*delrl[2] - delrj[2]*delrl[1];
    tork[1] = delrk[2]*delrj[0] - delrk[0]*delrj[2];
    torl[1] = delrj[2]*delrl[0] - delrj[0]*delrl[2];
    tork[2] = delrk[0]*delrj[1] - delrk[1]*delrj[0];
    torl[2] = delrj[0]*delrl[1] - delrj[1]*delrl[0];
    TT2 = dot3(tork,torl);
    rmut = square(TT2/TT1);
    if (torindx>=1) {
      btt = 1.0 - rmut;
      return btt * fc1k * fc1l;
    }
    else {
      btt=paramk->ptork1-TT2/TT1;
      btt=paramk->ptork2*square(btt);
      return btt * fc1k * fc1l;
    }

  } else {
    return 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void PairComb3::tor_calc(double r, Param *parami, Param *paramj,
        double kconjug, double lconjug, int /*i*/, int /*j*/, double xcn, double ycn)
{
  int ixmin, iymin, izmin;
  double vtor, dtorx, dtory, dtorz;
  double xtor, ytor, zcon;
  int torindx;

  vtor = dtorx = dtory = dtorz = 0.0;
  torindx=parami->tor_flag;

  if (torindx<0) {
    vtor=1.0;
    dtorx=0.0;
    dtory=0.0;
    dtorz=0.0;
  } else {
    xtor = -comb_fc(r, parami) * parami->pcross + xcn;
    ytor = -comb_fc(r, paramj) * paramj->pcross + ycn;
    zcon = 1.0 + square(kconjug) + square(lconjug);
    if (xtor < 0.0) xtor = 0.0;
    if (ytor < 0.0) ytor = 0.0;
    if (zcon < 1.0) zcon = 1.0;
    if (xtor > maxxc) xtor = maxxc;
    if (ytor > maxyc) ytor = maxyc;
    if (zcon > maxconj) zcon = maxconj;

    ixmin = int(xtor+1.0e-12);
    iymin = int(ytor+1.0e-12);
    izmin = int(zcon+1.0e-12);

    torindx=torindx-1;

    if (fabs(float(ixmin)-xtor)>1.0e-8 ||
      fabs(float(iymin)-ytor)>1.0e-8 ||
      fabs(float(izmin)-zcon)>1.0e-8) {
      tor_int(torindx,xtor,ytor,zcon,ixmin,iymin,izmin,
              vtor,dtorx,dtory,dtorz);
    } else {
      vtor  = tor_grid[torindx][ixmin][iymin][izmin-1];
      dtorx = tor_gridx[torindx][ixmin][iymin][izmin-1];
      dtory = tor_gridy[torindx][ixmin][iymin][izmin-1];
      dtorz = tor_gridz[torindx][ixmin][iymin][izmin-1];
    }
  }

  btor[0] = vtor;
  btor[1] = dtorx;
  btor[2] = dtory;
  btor[3] = dtorz;
}

/* ---------------------------------------------------------------------- */

void PairComb3::tor_int(int torindx,double xtor, double ytor, double zcon, int l,
        int m, int n, double &vtor, double &dtorx, double &dtory, double &dtorz)
{
  int j;
  double x;

  vtor = dtorx = dtory = dtorz = 0.0;
  if (l >= maxxc-1) { l=maxxc-1; }  //boundary condition changed
  if (m >= maxyc-1) { m=maxyc-1; }
  if (n >= maxconj-1) { n=maxconj-1; }

  for (j=0; j<64; j++) {
    x = tor_spl[torindx][l][m][n-1][j] * pow(xtor,iin3[j][0])
          * pow(ytor,iin3[j][1]) * pow(zcon,iin3[j][2]);
    vtor += x;

  if (xtor > 1.0e-8) dtorx += x*iin3[j][0]/xtor;
  if (ytor > 1.0e-8) dtory += x*iin3[j][1]/ytor;
  if (zcon > 1.0e-8) dtorz += x*iin3[j][2]/zcon;
  }
}

/* ---------------------------------------------------------------------- */

void PairComb3::tor_force(int torindx, Param *paramk, Param *paraml,
        double srmu, double rsq1,double rsq2, double rsq3,
        double *delrj, double *delrk, double *delrl)
{
  int nm;
  double rmu, rmul, srmul, rij, rik, rjl;

  for (nm=0; nm<3; nm++) {
    fi_tor[nm] = fj_tor[nm] = fk_tor[nm] = fl_tor[nm] = 0.0;
  }

  rij = sqrt(rsq1);
  rik = sqrt(rsq2);
  rjl = sqrt(rsq3);

  rmu = dot3(delrj,delrk)/(rij*rik);
  scale3(-1.0,delrl,delrl);
  rmul = dot3(delrj,delrl)/(rij*rjl);
  scale3(-1.0,delrl,delrl);
  srmul = sqrt(1.0-rmul*rmul);
  if (acos(rmul) > MY_PI) srmul = -srmul;

  if (srmul > 0.1) {
    double fc1k, fcp1k, fc1l, fcp1l, srmul2, dt1dik, dt1djl;
    double TT1, TT2, rmut, btt, tork[3], torl[3];
    double dt2dik[3], dt2djl[3], dt2dij[3], AA, AA2;
    double tfij[4], tfik[2], tfjl[2], tjx[3], tjy[3], tjz[3];
    double tkx[2], tky[2], tkz[2], tlx[2], tly[2], tlz[2];

    fc1k  = comb_fc(rik,paramk);
    fcp1k = comb_fc_d(rik,paramk);
    fc1l  = comb_fc(rjl,paraml);
    fcp1l = comb_fc_d(rjl,paraml);
    srmul2 = square(srmul);

    TT1 = rik*rjl*rij*rij*srmu*srmul;
    dt1dik = -rmu/square(srmu);
    dt1djl = -rmul/srmul2;
    tork[0] = delrk[1]*delrj[2] - delrk[2]*delrj[1];
    torl[0] = delrj[1]*delrl[2] - delrj[2]*delrl[1];
    tork[1] = delrk[2]*delrj[0] - delrk[0]*delrj[2];
    torl[1] = delrj[2]*delrl[0] - delrj[0]*delrl[2];
    tork[2] = delrk[0]*delrj[1] - delrk[1]*delrj[0];
    torl[2] = delrj[0]*delrl[1] - delrj[1]*delrl[0];
    TT2 = dot3(tork,torl);

    dt2dik[0] = -delrj[1]*torl[2] + delrj[2]*torl[1];
    dt2dik[1] = -delrj[2]*torl[0] + delrj[0]*torl[2];
    dt2dik[2] = -delrj[0]*torl[1] + delrj[1]*torl[0];
    dt2djl[0] =  delrj[1]*tork[2] - delrj[2]*tork[1];
    dt2djl[1] =  delrj[2]*tork[0] - delrj[0]*tork[2];
    dt2djl[2] =  delrj[0]*tork[1] - delrj[1]*tork[0];
    dt2dij[0] = -delrk[2]*torl[1] + delrl[2]*tork[1]
               + delrk[1]*torl[2] - delrl[1]*tork[2];
    dt2dij[1] = -delrk[0]*torl[2] + delrl[0]*tork[2]
               + delrk[2]*torl[0] - delrl[2]*tork[0];
    dt2dij[2] = -delrk[1]*torl[0] + delrl[1]*tork[0]
               + delrk[0]*torl[1] - delrl[0]*tork[1];

    rmut = TT2/TT1;

    if (torindx>=1) {
        btt = 1.0 - square(rmut);
        AA = -2.0 * ptorr * rmut * fc1k * fc1l / TT1;
    }
    else {
        btt=paramk->ptork1-rmut;
        btt=paramk->ptork2*square(btt);
        AA = -2.0 * ptorr * paramk->ptork2 *
          (paramk->ptork1-rmut) * fc1k * fc1l /TT1;
   }

    AA2 = AA * TT2;
    tfij[0] = -(dt1dik*AA2)/rij/rik;
    tfij[1] = AA2/rij/rij - dt1dik*AA2*rmu/rij/rij;
    tfij[2] = -dt1djl*AA2/rij/rjl;
    tfij[3] = AA2/rij/rij - dt1djl*AA2*rmul/rij/rij;
    tfik[0] = tfij[0];
    tfik[1] = (AA2/rik - btt*ptorr*fc1l*fcp1k)/rik -
            dt1dik*AA2*rmu/rik/rik;
    tfjl[0] = tfij[2];
    tfjl[1] = (AA2/rjl - btt*ptorr*fc1k*fcp1l)/rjl -
            dt1djl*AA2*rmul/rjl/rjl;

    tjx[0] = tfij[0]*delrk[0] - tfij[1]*delrj[0];
    tjy[0] = tfij[0]*delrk[1] - tfij[1]*delrj[1];
    tjz[0] = tfij[0]*delrk[2] - tfij[1]*delrj[2];
    tjx[1] = -tfij[2]*delrl[0] - tfij[3]*delrj[0];
    tjy[1] = -tfij[2]*delrl[1] - tfij[3]*delrj[1];
    tjz[1] = -tfij[2]*delrl[2] - tfij[3]*delrj[2];
    tjx[2] = -dt2dij[0] * AA;
    tjy[2] = -dt2dij[1] * AA;
    tjz[2] = -dt2dij[2] * AA;

    tkx[0] = tfik[0]*delrj[0] - tfik[1]*delrk[0];
    tky[0] = tfik[0]*delrj[1] - tfik[1]*delrk[1];
    tkz[0] = tfik[0]*delrj[2] - tfik[1]*delrk[2];
    tkx[1] = -dt2dik[0] * AA;
    tky[1] = -dt2dik[1] * AA;
    tkz[1] = -dt2dik[2] * AA;

    tlx[0] = -tfjl[0]*delrj[0] - tfjl[1]*delrl[0];
    tly[0] = -tfjl[0]*delrj[1] - tfjl[1]*delrl[1];
    tlz[0] = -tfjl[0]*delrj[2] - tfjl[1]*delrl[2];
    tlx[1] = -dt2djl[0] * AA;
    tly[1] = -dt2djl[1] * AA;
    tlz[1] = -dt2djl[2] * AA;

    fi_tor[0] = tjx[0]+tjx[1]+tjx[2]+tkx[0]+tkx[1];
    fi_tor[1] = tjy[0]+tjy[1]+tjy[2]+tky[0]+tky[1];
    fi_tor[2] = tjz[0]+tjz[1]+tjz[2]+tkz[0]+tkz[1];

    fj_tor[0] = -tjx[0]-tjx[1]-tjx[2]+tlx[0]+tlx[1];
    fj_tor[1] = -tjy[0]-tjy[1]-tjy[2]+tly[0]+tly[1];
    fj_tor[2] = -tjz[0]-tjz[1]-tjz[2]+tlz[0]+tlz[1];

    fk_tor[0] = -tkx[0]-tkx[1];
    fk_tor[1] = -tky[0]-tky[1];
    fk_tor[2] = -tkz[0]-tkz[1];

    fl_tor[0] = -tlx[0]-tlx[1];
    fl_tor[1] = -tly[0]-tly[1];
    fl_tor[2] = -tlz[0]-tlz[1];

  }
}

/* ---------------------------------------------------------------------- */

double PairComb3::combqeq(double *qf_fix, int &igroup)
{
  int i,j,ii, jj,itype,jtype,jnum;
  int iparam_i,iparam_ji,iparam_ij;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int mr1,mr2,mr3,inty,nj;
  double xtmp,ytmp,ztmp,rsq1,delrj[3];
  double iq,jq,fqi,fqij,fqji,sr1,sr2,sr3;
  double potal,fac11,fac11e;
  int sht_jnum,*sht_jlist;
  tagint itag, jtag;

  double **x = atom->x;
  double *q = atom->q;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int inum = list->inum;
  int *mask = atom->mask;
  int groupbit = group->bitmask[igroup];

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  qf = qf_fix;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      qf[i] = 0.0;
      dpl[i][0] = dpl[i][1] = dpl[i][2] = 0.0;
    }
  }
  // communicating charge force to all nodes, first forward then reverse

  pack_flag = 1;
  comm->forward_comm(this);

  // self energy correction term: potal

  potal_calc(potal,fac11,fac11e);

  // loop over full neighbor list of my atoms

  fqi = fqij = fqji = 0.0;

  for (ii = 0; ii < inum; ii ++) {
    i = ilist[ii];
    itag = tag[i];
    nj = 0;
    if (mask[i] & groupbit) {
      itype = map[type[i]];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      iq = q[i];
      iparam_i = elem3param[itype][itype][itype];

      // charge force from self energy
      fqi = qfo_self(&params[iparam_i],iq);

      jlist = firstneigh[i];
      jnum = numneigh[i];

      sht_jlist = sht_first[i];
      sht_jnum = sht_num[i];

      // two-body interactions

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj] & NEIGHMASK;

        jtag = tag[j];
        if (itag >= jtag) continue;

        jtype = map[type[j]];
        inty = intype[itype][jtype];
        jq = q[j];

        delrj[0] = xtmp - x[j][0];
        delrj[1] = ytmp - x[j][1];
        delrj[2] = ztmp - x[j][2];
        rsq1 = dot3(delrj,delrj);

        iparam_ij = elem3param[itype][jtype][jtype];
        iparam_ji = elem3param[jtype][itype][itype];

        // long range q-dependent

        if (rsq1 > params[iparam_ij].lcutsq) continue;

        // polynomial three-point interpolation
        tri_point(rsq1,mr1,mr2,mr3,sr1,sr2,sr3);

        // 1/r charge forces
        qfo_direct(&params[iparam_ij],&params[iparam_ji],
                mr1,mr2,mr3,rsq1,sr1,sr2,sr3,fac11e,fqij,fqji,
                iq,jq,i,j);

        fqi += fqij;  qf[j] += fqji;

       // field correction to self energy and charge force
        qfo_field(&params[iparam_ij],&params[iparam_ji],rsq1,
                iq,jq,fqij,fqji);

        fqi += fqij;  qf[j] += fqji;

        // polarization field charge force
        if (pol_flag) {
          qfo_dipole(fac11,mr1,mr2,mr3,inty,rsq1,delrj,sr1,sr2,sr3,
                fqij,fqji,i,j);

          fqi += fqij;  qf[j] += fqji;
        }
      }

      for (jj = 0; jj < sht_jnum; jj++) {
        j = sht_jlist[jj];

        jtag = tag[j];
        if (itag >= jtag) continue;

        jtype = map[type[j]];
        inty = intype[itype][jtype];
        jq = q[j];

        delrj[0] = xtmp - x[j][0];
        delrj[1] = ytmp - x[j][1];
        delrj[2] = ztmp - x[j][2];
        rsq1 = dot3(delrj,delrj);

        iparam_ij = elem3param[itype][jtype][jtype];
        iparam_ji = elem3param[jtype][itype][itype];

        if (rsq1 >= params[iparam_ij].cutsq) continue;
        nj ++;

        // charge force in Aij and Bij
         qfo_short(&params[iparam_ij],&params[iparam_ji],
                rsq1,iq,jq,fqij,fqji,i,j,nj);

        fqi += fqij;  qf[j] += fqji;
      }
      qf[i] += fqi;
    }
  }

  comm->reverse_comm(this);

  // sum charge force on each node and return it

  double eneg = 0.0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      eneg += qf[i];
          itag=tag[i];
    }
  }

  MPI_Allreduce(&eneg,&enegtot,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Bcast(&enegtot,1,MPI_DOUBLE,0,world);
  return enegtot;

}

/* ---------------------------------------------------------------------- */

double PairComb3::qfo_self(Param *param, double qi)
{
  double self_d,cmin,cmax,qmin,qmax;
  double s1 = param->chi;
  double s2 = param->dj;
  double s3 = param->dk;
  double s4 = param->dl;

  self_d = 0.0;

  qmin = param->qmin;
  qmax = param->qmax;
  cmin = cmax = 100.0;
  self_d = s1+qi*(2.0*s2+qi*(3.0*s3+qi*4.0*s4));

  if (qi < qmin) self_d += 4.0 * cmin * cube((qi-qmin));
  if (qi > qmax) self_d += 4.0 * cmax * cube((qi-qmax));

  return self_d;
}

/* ---------------------------------------------------------------------- */

void PairComb3::qfo_direct(Param *parami, Param *paramj, int mr1,
          int mr2, int mr3, double rsq, double sr1, double sr2,
          double sr3, double fac11e, double &fqij, double &fqji,
          double iq, double jq, int i, int j)
{
  double r, erfcc, fafbnl, vm, vmfafb, esucon;
  double afbn, afbj, sme1n, sme1j;
  double curli = parami->curl;
  double curlj = paramj->curl;
  int inti = parami->ielement;
  int intj = paramj->ielement;
  int inty = intype[inti][intj];

  double curlij0 = parami->curl0;
  double curlji0 = paramj->curl0;
  double curlij1,curlji1;
  int icurl, jcurl;
  int ielegp = parami->ielementgp;
  int jelegp = paramj->ielementgp;

  r = sqrt(rsq);
  esucon=force->qqr2e;

  icurl = jcurl = 0;
  if (ielegp==2 && curli>curlij0) {
    icurl=1;
    curlij1=curli;
  }

  if (jelegp==2 && curlj>curlji0) {
    jcurl=1;
    curlji1=curlj;
  }
  if (icurl==1 || jcurl ==1) {
    double xcoij= xcotmp[i];
    double xcoji= xcotmp[j];

    if (icurl==1) {
      curli=curlij1+(curlij0-curlij1)*comb_fc_curl(xcoij,parami);
    }
    if (jcurl==1) {
      curlj=curlji1+(curlji0-curlji1)*comb_fc_curl(xcoji,paramj);
    }
  }

  // 1/r force (wrt q)

  erfcc = sr1*erpaw[mr1][0]   + sr2*erpaw[mr2][0]   + sr3*erpaw[mr3][0];
  fafbnl= sr1*fafb[mr1][inty] + sr2*fafb[mr2][inty] + sr3*fafb[mr3][inty];
  afbn  = sr1*afb[mr1][inti]  + sr2*afb[mr2][inti]  + sr3*afb[mr3][inti];
  afbj  = sr1*afb[mr1][intj]  + sr2*afb[mr2][intj]  + sr3*afb[mr3][intj];
  vm = (erfcc/r * esucon - fac11e);
  vmfafb = vm + esucon * fafbnl;
  sme1n = curlj * (afbn - fafbnl) * esucon;
  sme1j = curli * (afbj - fafbnl) * esucon;
  fqij = 1.0 * (jq * vmfafb + sme1n);
  fqji = 1.0 * (iq * vmfafb + sme1j);
}

/* ---------------------------------------------------------------------- */

void PairComb3::qfo_field(Param *parami, Param *paramj, double rsq,
        double iq,double jq, double &fqij, double &fqji)
{
  double r,r3,r5,rc,rc2,rc3,rc4,rc5;
  double cmi1,cmi2,cmj1,cmj2,pcmi1,pcmi2;
  double rf3i,rcf3i,rf5i,rcf5i;
  double drcf3i,drcf5i,rf3,rf5;

  r  = sqrt(rsq);
  r3 = r * rsq;
  r5 = r3 * rsq;
  rc = parami->lcut;
  rc2=  rc*rc;
  rc3 = rc*rc*rc;
  rc4 = rc3 * rc;
  rc5 = rc4 * rc;
  cmi1 = parami->cmn1;
  cmi2 = parami->cmn2;
  cmj1 = paramj->cmn1;
  cmj2 = paramj->cmn2;
  pcmi1 = parami->pcmn1;
  pcmi2 = parami->pcmn2;

  rf3i = r3/(square(r3)+cube(pcmi1));
  rcf3i = rc3/(square(rc3)+cube(pcmi1));
  rf5i = r5/(square(r5)+powint(pcmi2,5));
  rcf5i = rc5/(square(rc5)+powint(pcmi2,5));

  drcf3i = 3/rc*rcf3i-6*rc2*rcf3i*rcf3i;
  drcf5i = 5/rc*rcf5i-10*rc4*rcf5i*rcf5i;

  rf3 = rf3i-rcf3i-(r-rc)*drcf3i;
  rf5 = rf5i-rcf5i-(r-rc)*drcf5i;

  // field correction charge force
  fqij = 1.0 * cmj1*rf3+2.0*iq*cmj2*rf5;
  fqji = 1.0 * cmi1*rf3+2.0*jq*cmi2*rf5;

}

/* ---------------------------------------------------------------------- */

void PairComb3::qfo_dipole(double fac11, int mr1, int mr2, int mr3,
        int inty, double rsq, double *delrj, double sr1, double sr2,
        double sr3, double &fqij, double &fqji, int i, int j)
{
  double erfcc, erfcd, dvdrr, dfafbnl, smf2;
  double r, r3, alfdpi, esucon;

  r = sqrt(rsq);
  r3 = r * rsq;
  alfdpi = 0.4/MY_PIS;
  esucon = force->qqr2e;

  erfcc = sr1*erpaw[mr1][0] + sr2*erpaw[mr2][0] + sr3*erpaw[mr3][0];
  erfcd = sr1*erpaw[mr1][1] + sr2*erpaw[mr2][1] + sr3*erpaw[mr3][1];
  dvdrr = (erfcc/r3+alfdpi*erfcd/rsq)*esucon-fac11;
  dfafbnl= sr1*dfafb[mr1][inty] + sr2*dfafb[mr2][inty] + sr3*dfafb[mr3][inty];
  smf2 = (dvdrr + dfafbnl*esucon)/r;

  fqij = dpl[i][0]*delrj[0] + dpl[i][1]*delrj[1] +dpl[i][2]*delrj[2];
  fqji = dpl[j][0]*delrj[0] + dpl[j][1]*delrj[1] +dpl[j][2]*delrj[2];
  fqij *= smf2;
  fqji *= smf2;

}

/* ---------------------------------------------------------------------- */

void PairComb3::qfo_short(Param *parami, Param *paramj, double rsq,
        double iq, double jq, double &fqij, double &fqji,
        int i, int /*j*/, int nj)
{
  double r, tmp_fc;
  double Di, Dj, dDi, dDj, Bsi, Bsj, dBsi, dBsj;
  double QUchi, QOchi, QUchj, QOchj;
  double bij, caj, cbj, caqpn, caqpj, cbqpn, cbqpj;
  double LamDiLamDj, AlfDiAlfDj;
  double rlm1 = parami->lambda;
  double alfij1= parami->alpha1;
  double alfij2= parami->alpha2;
  double alfij3= parami->alpha3;
  double pbij1= parami->bigB1;
  double pbij2= parami->bigB2;
  double pbij3= parami->bigB3;

  caj = cbj = caqpn = caqpj = cbqpn = cbqpj = 0.0;
  r = sqrt(rsq);
  tmp_fc = comb_fc(r,parami);
  bij = bbij[i][nj];

  QUchi = (parami->QU - iq) * parami->bD;
  QUchj = (paramj->QU - jq) * paramj->bD;
  QOchi = (iq - parami->Qo) * parami->bB;
  QOchj = (jq - paramj->Qo) * paramj->bB;

  if (iq < parami->QL-0.2) {
    iq = parami->QL-0.2;
    Di = parami->DL;
    dDi = Bsi = dBsi = 0.0;
  } else if (iq > parami->QU+0.2) {
    iq = parami->QU+0.2;
    Di = parami->DU;
    dDi = Bsi = dBsi = 0.0;
  } else {
    Di = parami->DU + pow(QUchi,parami->nD);                            // YYDin
    dDi = -parami->nD * parami->bD * pow(QUchi,(parami->nD-1.0));       // YYDiqp
    Bsi = parami->aB - pow(QOchi,10);                                   // YYBsin
    dBsi = -parami->bB * 10.0 * powint(QOchi,9);                        // YYBsiqp
  }

  if (jq < paramj->QL-0.2) {
    jq = paramj->QL-0.2;
    Dj = paramj->DL;
    dDj = Bsj = dBsj = 0.0;
  } else if (jq > paramj->QU+0.2) {
    jq = paramj->QU+0.2;
    Dj = paramj->DU;
    dDj = Bsj = dBsj = 0.0;
  } else {
    Dj = paramj->DU + pow(QUchj,paramj->nD);                            // YYDij
    dDj = -paramj->nD * paramj->bD * pow(QUchj,(paramj->nD-1.0));       // YYDiqpj
    Bsj = paramj->aB - pow(QOchj,10);                                   // YYBsij
    dBsj = -paramj->bB * 10.0 * powint(QOchj,9);                        // YYBsiqpj
  }

  LamDiLamDj = exp(0.5*(parami->lami*Di+paramj->lami*Dj)-rlm1*r);
  caj = 0.5 * tmp_fc * parami->bigA * LamDiLamDj;

  if (Bsi*Bsj > 0.0) {
    AlfDiAlfDj = exp(0.5*(parami->alfi*Di+paramj->alfi*Dj));
    cbj=-0.5*tmp_fc*bij*sqrt(Bsi*Bsj)*AlfDiAlfDj*
                (pbij1*exp(-alfij1*r)+pbij2*exp(-alfij2*r)+pbij3*exp(-alfij3*r));
    cbqpn = cbj * (parami->alfi * dDi + dBsi/Bsi);
    cbqpj = cbj * (paramj->alfi * dDj + dBsj/Bsj);
  } else {
    cbj = cbqpn = cbqpj = 0.0;
  }

  caqpn = caj * parami->lami * dDi;
  caqpj = caj * paramj->lami * dDj;

  fqij = 1.0 * (caqpn + cbqpn);
  fqji = 1.0 * (caqpj + cbqpj);

}

/* ---------------------------------------------------------------------- */

void PairComb3::dipole_init(Param *parami, Param *paramj, double fac11,
        double *delrj, double rsq, int mr1, int mr2, int mr3, double sr1,
        double sr2, double sr3, double iq, double jq, int i, int j)
{
  double erfcc, erfcd, dvdrr, dfafbnl, smf2, phinn, phinj, efn, efj;
  double r, r3, alfdpi, esucon;
  double rcd, rct, tmurn, tmurj, poln[3], polj[3], Qext[3];
  int nm;
  int inti = parami->ielement;
  int intj = paramj->ielement;
  int inty = intype[inti][intj];

  for (nm=0; nm<3; nm++) Qext[nm] = 0.0;

  r = sqrt(rsq);
  r3 = r * rsq;
  rcd = 1.0/(r3);
  rct = 3.0*rcd/rsq;
  alfdpi = 0.4/MY_PIS;
  esucon = force->qqr2e;

  erfcc = sr1*erpaw[mr1][0] + sr2*erpaw[mr2][0] + sr3*erpaw[mr3][0];
  erfcd = sr1*erpaw[mr1][1] + sr2*erpaw[mr2][1] + sr3*erpaw[mr3][1];
  dvdrr = (erfcc/r3+alfdpi*erfcd/rsq)*esucon-fac11;
  dfafbnl= sr1*dfafb[mr1][inty] + sr2*dfafb[mr2][inty] + sr3*dfafb[mr3][inty];
  smf2 = dvdrr/esucon + dfafbnl/r;
  phinn = sr1*phin[mr1][inti] + sr2*phin[mr2][inti] + sr3*phin[mr3][inti];
  phinj = sr1*phin[mr1][intj] + sr2*phin[mr2][intj] + sr3*phin[mr3][intj];
  efn = jq * smf2;
  efj = iq * smf2;

  tmurn = dpl[i][0]*delrj[0] + dpl[i][1]*delrj[1] + dpl[i][2]*delrj[2];
  tmurj = dpl[j][0]*delrj[0] + dpl[j][1]*delrj[1] + dpl[j][2]*delrj[2];

  for (nm=0; nm<3; nm++) {
    poln[nm] = (tmurj*delrj[nm]*rct - dpl[j][nm]*rcd)*phinj;
    polj[nm] = (tmurn*delrj[nm]*rct - dpl[i][nm]*rcd)*phinn;
  }

  for (nm=0; nm<3; nm++) {
    dpl[i][nm] += (Qext[nm]/esucon + delrj[nm]*efn + poln[nm])*parami->polz*0.50;
    dpl[j][nm] += (Qext[nm]/esucon - delrj[nm]*efj + polj[nm])*paramj->polz*0.50;
  }

}

/* ---------------------------------------------------------------------- */

double PairComb3::dipole_self(Param *parami, int i)
{
  double esucon = force->qqr2e;
  double apn = parami->polz;
  double selfdpV = 0.0;

  if (apn != 0.0) {
      selfdpV= (dpl[i][0]*dpl[i][0]+dpl[i][1]*dpl[i][1]+dpl[i][2]*dpl[i][2])
                *esucon/(2.0*apn); }
  return selfdpV;
}

/* ---------------------------------------------------------------------- */

void PairComb3::dipole_calc(Param *parami, Param *paramj, double fac11,
        double delx, double dely, double delz, double rsq,
        int mr1, int mr2, int mr3, double sr1, double sr2, double sr3,
        double iq, double jq, int i, int j, double &vionij,
        double &fvionij, double *ddprx)
{
  double erfcc, erfcd, dvdrr, dfafbnl, ef, phinn, phinj, efn, efj;
  double r, r3, alf, alfdpi, esucon, dphinn, dphinj, ddfafbnl;
  double def, defn, defj, tmun, tmuj, emuTmu, edqn, edqj, ddvdrr;
  double rcd, rct, tmurn, tmurj, tmumu, poln[3], polj[3], delr1[3];
  double demuTmu, ddpr, dcoef;
  int nm;
  int inti = parami->ielement;
  int intj = paramj->ielement;
  int inty = intype[inti][intj];

  r = sqrt(rsq);
  r3 = r * rsq;
  esucon = force->qqr2e;
  rcd = esucon/r3;
  rct = 3.0*rcd/rsq;
  alf = 0.2;
  alfdpi = 2.0*alf/MY_PIS;
  delr1[0] = delx;
  delr1[1] = dely;
  delr1[2] = delz;

  // generate energy & force information from tables
  erfcc = sr1*erpaw[mr1][0] + sr2*erpaw[mr2][0] + sr3*erpaw[mr3][0];
  erfcd = sr1*erpaw[mr1][1] + sr2*erpaw[mr2][1] + sr3*erpaw[mr3][1];
  dvdrr = (erfcc/r3+alfdpi*erfcd/rsq)*esucon-fac11;
  ddvdrr = (2.0*erfcc/r3 + 2.0*alfdpi*erfcd*(1.0/rsq+alf*alf))*esucon;
  dfafbnl= sr1*dfafb[mr1][inty] + sr2*dfafb[mr2][inty] + sr3*dfafb[mr3][inty];
  phinn = sr1*phin[mr1][inti] + sr2*phin[mr2][inti] + sr3*phin[mr3][inti];
  phinj = sr1*phin[mr1][intj] + sr2*phin[mr2][intj] + sr3*phin[mr3][intj];
  dphinn = sr1*dphin[mr1][inti] + sr2*dphin[mr2][inti] + sr3*dphin[mr3][inti];
  dphinj = sr1*dphin[mr1][intj] + sr2*dphin[mr2][intj] + sr3*dphin[mr3][intj];
  ddfafbnl= sr1*ddfafb[mr1][inty] + sr2*ddfafb[mr2][inty] + sr3*ddfafb[mr3][inty];
  ef = (dvdrr + dfafbnl * esucon)/r;
  efn =  jq * ef;
  efj = -iq * ef;
  def = (ddvdrr + ddfafbnl * esucon)/r;
  defn =  jq * def;
  defj = -iq * def;

  // dipole - dipole field tensor (Tij)
  tmurn = dpl[i][0]*delr1[0] + dpl[i][1]*delr1[1] + dpl[i][2]*delr1[2];
  tmurj = dpl[j][0]*delr1[0] + dpl[j][1]*delr1[1] + dpl[j][2]*delr1[2];
  tmumu = dpl[i][0]*dpl[j][0] + dpl[i][1]*dpl[j][1] + dpl[i][2]*dpl[j][2];

  for (nm=0; nm<3; nm++) {
    poln[nm] = (tmurj*delr1[nm]*rct - dpl[j][nm]*rcd);
    polj[nm] = (tmurn*delr1[nm]*rct - dpl[i][nm]*rcd);
  }
  tmun = dpl[j][0]*polj[0] + dpl[j][1]*polj[1] + dpl[j][2]*polj[2];
  tmuj = dpl[i][0]*poln[0] + dpl[i][1]*poln[1] + dpl[i][2]*poln[2];

  // dipole - dipole energy
  emuTmu = -0.5*(tmun*phinn+tmuj*phinj);

  // dipole - charge energy
  edqn = -0.5 * (tmurn * efn);
  edqj = -0.5 * (tmurj * efj);

  // overall dipole energy
  vionij = emuTmu + edqn + edqj;

  // dipole - dipole force
  demuTmu = (tmun*dphinn + tmuj*dphinj)/r;
  ddpr =  5.0*tmurn*tmurj/rsq - tmumu;
  dcoef = rct * (phinn+phinj);

  for (nm = 0; nm < 3; nm ++) {
    ddprx[nm] = dcoef * (ddpr*delr1[nm] - tmurn*dpl[j][nm] - tmurj*dpl[i][nm])
          + demuTmu * delr1[nm];
  }

  // dipole - charge force
  fvionij = -tmurn*defn - tmurj*defj;
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_fc_curl(double rocn, Param *param)
{
  double r_inn = param->curlcut1;
  double r_out = param->curlcut2;
  if (rocn <= r_inn) return 1.0;
  if (rocn >= r_out) return 0.0;
  return 0.5*(1.0 + cos(MY_PI*(rocn-r_inn)/(r_out-r_inn)));
}

/* ---------------------------------------------------------------------- */

double PairComb3::comb_fc_curl_d(double rocn, Param *param)
{
  double r_inn = param->curlcut1;
  double r_out = param->curlcut2;
  if (rocn <= r_inn) return 0.0;
  if (rocn >= r_out) return 0.0;
  return -MY_PI2/(r_out-r_inn)*sin(MY_PI*(rocn-r_inn)/(r_out-r_inn));
}

/* ---------------------------------------------------------------------- */

int PairComb3::heaviside(double rr)
{
  if (rr <= 0.0) return 0;
  else return 1;
}

/* ---------------------------------------------------------------------- */

double PairComb3::switching(double rr)
{
  if (rr <= 0.0) return 1.0;
  else if (rr >= 1.0) return 0.0;
  else return heaviside(-rr)+heaviside(rr)*heaviside(1.0-rr)
          * (1.0-(3.0-2.0*rr)*rr*rr);
}

/* ---------------------------------------------------------------------- */

double PairComb3::switching_d(double rr)
{
  if (rr <= 0.0) return 0.0;
  else if (rr >= 1.0) return 0.0;
  else return heaviside(rr)*heaviside(1.0-rr)
          * 6.0*rr*(rr-1.0);
}

/* ---------------------------------------------------------------------- */

int PairComb3::pack_forward_comm(int n, int *list, double *buf,
                                 int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  if (pack_flag == 1) {
    for (i = 0; i < n; i ++) {
      j = list[i];
      buf[m++] = qf[j];
    }
  } else if (pack_flag == 2) {
    for (i = 0; i < n; i ++) {
      j = list[i];
      buf[m++] = NCo[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairComb3::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n ;
  if (pack_flag == 1) {
    for (i = first; i < last; i++)
      qf[i] = buf[m++];
  } else if (pack_flag == 2) {
    for (i = first; i < last; i++)
      NCo[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairComb3::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (pack_flag == 1) {
    for (i = first; i < last; i++)
      buf[m++] = qf[i];
  } else if (pack_flag == 2) {
    for (i = first; i < last; i++)
      buf[m++] = NCo[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairComb3::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (pack_flag == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      qf[j] += buf[m++];
    }
  } else if (pack_flag == 2) {
    for (i = 0; i < n; i++) {
      j = list[i];
      NCo[j] += buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairComb3::memory_usage()
{
  double bytes = (double)maxeatom * sizeof(double);
  bytes += (double)maxvatom*6 * sizeof(double);
  bytes += (double)nmax * sizeof(int);
  bytes += (double)nmax * 8.0 * sizeof(double);
  bytes += (double)25000*2*sizeof(double);

  for (int i = 0; i < comm->nthreads; i++)
    bytes += ipage[i].size();

  return bytes;
}
