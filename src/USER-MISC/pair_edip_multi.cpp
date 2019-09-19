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
   Environment Dependent Interatomic Potential

   Contributing author: Chao Jiang
------------------------------------------------------------------------- */

#include "pair_edip_multi.h"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "citeme.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4


static const char cite_pair_edip[] =
  "@article{cjiang2012\n"
  " author    = {Jian, Chao and Morgan, Dane, and Szlufarska, Izabella},\n"
  " title     = {Carbon tri-interstitial defect: A model for DII center},\n"
  " journal   = {Physical Review B},\n"
  " volume    = {86},\n"
  " pages     = {144118},\n"
  " year      = {2012},\n"
  "}\n\n"
  "@article{lpizzagalli2010,\n"
  " author    = {G. Lucas, M. Bertolus, and L. Pizzagalli},\n"
  " journal   = {J. Phys. : Condens. Matter 22},\n"
  " volume    = {22},\n"
  " pages     = {035802},\n"
  " year      = {2010},\n"
  "}\n\n";



/* ---------------------------------------------------------------------- */

PairEDIPMulti::PairEDIPMulti(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_edip);

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairEDIPMulti::~PairEDIPMulti()
{
  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;

//XXX    deallocateGrids();
    deallocatePreLoops();
  }
}

/* ---------------------------------------------------------------------- */

void PairEDIPMulti::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum;
  int itype,jtype,ktype,ijparam,ikparam,ijkparam;
  double xtmp,ytmp,ztmp,evdwl;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int preForceCoord_counter;

  double zeta_i;
  double dzetair;
  double fpair;
  double costheta;
  double dpairZ,dtripleZ;

 // eflag != 0 means compute energy contributions in this step
 // vflag != 0 means compute virial contributions in this step

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;//total number of atoms in the cell
  ilist = list->ilist;//list of atoms
  numneigh = list->numneigh;//number of near neighbors
  firstneigh = list->firstneigh;//list of neighbors

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    zeta_i = 0.0;
    int numForceCoordPairs = 0;

    i = ilist[ii];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // all the neighbors of atom i

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // pre-loop to compute environment coordination f(Z)

    for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        double delx, dely, delz, r_ij;

        delx = x[j][0] - xtmp;
        dely = x[j][1] - ytmp;
        delz = x[j][2] - ztmp;
        r_ij = delx * delx + dely * dely + delz * delz;

        jtype = map[type[j]];
        ijparam = elem2param[itype][jtype][jtype];
        if (r_ij > params[ijparam].cutsq) continue;

        r_ij = sqrt(r_ij);

        // zeta and its derivative dZ/dr

        if (r_ij < params[ijparam].cutoffC) zeta_i += 1.0;
        else {
            double f, fdr;
            edip_fc(r_ij, &params[ijparam], f, fdr);
            zeta_i += f;
            dzetair = -fdr / r_ij;

            preForceCoord_counter=numForceCoordPairs*5;
            preForceCoord[preForceCoord_counter+0]=dzetair;
            preForceCoord[preForceCoord_counter+1]=delx;
            preForceCoord[preForceCoord_counter+2]=dely;
            preForceCoord[preForceCoord_counter+3]=delz;
            preForceCoord[preForceCoord_counter+4]=j;
            numForceCoordPairs++;
        }
    }

    // two-body interactions

    dpairZ=0;
    dtripleZ=0;

    for (jj = 0; jj < jnum; jj++) {
      double dr_ij[3], r_ij, f_ij[3];

      j = jlist[jj];
      j &= NEIGHMASK;

      dr_ij[0] = x[j][0] - xtmp;
      dr_ij[1] = x[j][1] - ytmp;
      dr_ij[2] = x[j][2] - ztmp;
      r_ij = dr_ij[0]*dr_ij[0] + dr_ij[1]*dr_ij[1] + dr_ij[2]*dr_ij[2];

      jtype = map[type[j]];
      ijparam = elem2param[itype][jtype][jtype];
      if (r_ij > params[ijparam].cutsq) continue;

      r_ij = sqrt(r_ij);

      // potential energy and force
      // since pair i-j is different from pair j-i, double counting is
      // already considered in constructing the potential

      double fdr, fdZ;
      edip_pair(r_ij, zeta_i, &params[ijparam], evdwl, fdr, fdZ);
      fpair = -fdr / r_ij;
      dpairZ += fdZ;

      f[i][0] -= fpair * dr_ij[0];
      f[i][1] -= fpair * dr_ij[1];
      f[i][2] -= fpair * dr_ij[2];

      f[j][0] += fpair * dr_ij[0];
      f[j][1] += fpair * dr_ij[1];
      f[j][2] += fpair * dr_ij[2];

      if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, -dr_ij[0], -dr_ij[1], -dr_ij[2]);

      // three-body Forces

      for (kk = jj + 1; kk < jnum; kk++) {
          double dr_ik[3], r_ik, f_ik[3];

          k = jlist[kk];
          k &= NEIGHMASK;
          ktype = map[type[k]];
          ikparam = elem2param[itype][ktype][ktype];
          ijkparam = elem2param[itype][jtype][ktype];

          dr_ik[0] = x[k][0] - xtmp;
          dr_ik[1] = x[k][1] - ytmp;
          dr_ik[2] = x[k][2] - ztmp;
          r_ik = dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2];

          if (r_ik > params[ikparam].cutsq) continue;

          r_ik = sqrt(r_ik);

          costheta=vec3_dot(dr_ij, dr_ik) / r_ij / r_ik;

          double v1, v2, v3, v4, v5, v6, v7;

          edip_fcut3(r_ij, &params[ijparam], v1, v2);
          edip_fcut3(r_ik, &params[ikparam], v3, v4);
          edip_h(costheta, zeta_i, &params[ijkparam], v5, v6, v7);

          // potential energy and forces
          evdwl = v1 * v3 * v5;
          dtripleZ += v1 * v3 * v7;

          double dri[3], drj[3], drk[3];
          double dhl, dfr;

          dhl = v1 * v3 * v6;

          costheta_d(dr_ij, r_ij, dr_ik, r_ik, dri, drj, drk);

          f_ij[0] = -dhl * drj[0];
          f_ij[1] = -dhl * drj[1];
          f_ij[2] = -dhl * drj[2];
          f_ik[0] = -dhl * drk[0];
          f_ik[1] = -dhl * drk[1];
          f_ik[2] = -dhl * drk[2];

          dfr = v2 * v3 * v5;
          fpair = -dfr / r_ij;

          f_ij[0] += fpair * dr_ij[0];
          f_ij[1] += fpair * dr_ij[1];
          f_ij[2] += fpair * dr_ij[2];

          dfr = v1 * v4 * v5;
          fpair = -dfr / r_ik;

          f_ik[0] += fpair * dr_ik[0];
          f_ik[1] += fpair * dr_ik[1];
          f_ik[2] += fpair * dr_ik[2];

          f[j][0] += f_ij[0];
          f[j][1] += f_ij[1];
          f[j][2] += f_ij[2];

          f[k][0] += f_ik[0];
          f[k][1] += f_ik[1];
          f[k][2] += f_ik[2];

          f[i][0] -= f_ij[0] + f_ik[0];
          f[i][1] -= f_ij[1] + f_ik[1];
          f[i][2] -= f_ij[2] + f_ik[2];

          if (evflag) ev_tally3(i,j,k,evdwl,0.0,f_ij,f_ik,dr_ij,dr_ik);
      }
    }

    // forces due to environment coordination f(Z)
    for (int idx = 0; idx < numForceCoordPairs; idx++) {
        double delx, dely, delz;

        preForceCoord_counter = idx * 5;
        dzetair = preForceCoord[preForceCoord_counter+0];
        delx = preForceCoord[preForceCoord_counter+1];
        dely = preForceCoord[preForceCoord_counter+2];
        delz = preForceCoord[preForceCoord_counter+3];
        j = static_cast<int> (preForceCoord[preForceCoord_counter+4]);

        dzetair *= (dpairZ + dtripleZ);

        f[j][0] += dzetair * delx;
        f[j][1] += dzetair * dely;
        f[j][2] += dzetair * delz;

        f[i][0] -= dzetair * delx;
        f[i][1] -= dzetair * dely;
        f[i][2] -= dzetair * delz;

        evdwl = 0.0;
        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, dzetair, -delx, -dely, -delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

double sqr(double x)
{
  return x * x;
}

//pair Vij, partial derivatives dVij(r,Z)/dr and dVij(r,Z)/dZ
void PairEDIPMulti::edip_pair(double r, double z, Param *param, double &eng,
                         double &fdr, double &fZ)
{
  double A = param->A;
  double B = param->B;
  double rho = param->rho;
  double beta = param->beta;
  double v1,v2,v3,v4;

  v1 = pow(B / r, rho);
  v2 = exp(-beta * z * z);
  edip_fcut2(r, param, v3, v4);

  eng = A * (v1 - v2) * v3;
  fdr = A * (v1 - v2) * v4 + A * (-rho * v1 / r) * v3;
  fZ = A * (2 * beta * z * v2) * v3;
}

//function fc(r) in calculating coordination Z and derivative fc'(r)
void PairEDIPMulti::edip_fc(double r, Param *param, double &f, double &fdr)
{
  double a = param->cutoffA;
  double c = param->cutoffC;
  double alpha = param->alpha;
  double x;
  double v1, v2;

  if(r < c + 1E-6)
  {
    f=1.0;
    fdr=0.0;
    return;
  }

  if(r > a - 1E-6)
  {
    f=0.0;
    fdr=0.0;
    return;
  }

  x = (a - c) / (r - c);
  v1 = x * x * x;
  v2 = 1.0 / (1.0 - v1);

  f = exp(alpha * v2);
  fdr = (3.0 * x * v1 / (a - c)) * (-alpha * v2 * v2) * f;
}

//cut-off function for Vij and its derivative fcut2'(r)
void PairEDIPMulti::edip_fcut2(double r, Param *param, double &f, double &fdr)
{
  double sigma = param->sigma;
  double a = param->cutoffA;
  double v1;

  if(r > a - 1E-6)
  {
    f=0.0;
    fdr=0.0;
    return;
  }

  v1 = 1.0 / (r - a);
  f = exp(sigma * v1);
  fdr = (-sigma * v1 * v1) * f;
}

//function tau(Z) and its derivative tau'(Z)
void PairEDIPMulti::edip_tau(double z, Param *param, double &f, double &fdZ)
{
  double u1 = param->u1;
  double u2 = param->u2;
  double u3 = param->u3;
  double u4 = param->u4;
  double v1, v2;

  v1 = exp(-u4 * z);
  v2 = exp(-2.0 * u4 * z);

  f = u1 + u2 * u3 * v1 - u2 * v2;
  fdZ = -u2 * u3 * u4 * v1 + 2.0 * u2 * u4 * v2;
}

//function h(l,Z) and its partial derivatives dh(l,Z)/dl and dh(l,Z)/dZ
void PairEDIPMulti::edip_h(double l, double z, Param *param, double &f,
                      double &fdl, double &fdZ)
{
  double lambda = param->lambda;
  double eta = param->eta;
  double Q0 = param->Q0;
  double mu = param->mu;
  double Q, QdZ, Tau, TaudZ;
  double u2, du2l, du2Z;
  double v1, v2, v3;

  //function Q(Z)
  Q = Q0 * exp(-mu * z);
  //derivative Q'(Z)
  QdZ= -mu * Q;

  edip_tau(z, param, Tau, TaudZ);

  v1 = sqr(l + Tau);
  u2 = Q * v1;
  v2 = exp(-u2);

  f = lambda * (1 - v2 + eta * u2);

  //df/du2
  v3 = lambda * (v2 + eta);

  //du2/dl
  du2l = Q * 2 * (l + Tau);
  fdl = v3 * du2l;

  //du2/dZ
  du2Z = QdZ * v1 + Q * 2 * (l + Tau) * TaudZ;
  fdZ = v3 * du2Z;
}

//cut-off function for Vijk and its derivative fcut3'(r)
void PairEDIPMulti::edip_fcut3(double r, Param *param, double &f, double &fdr)
{
  double gamma = param->gamma;
  double a = param->cutoffA;
  double v1;

  if(r > a - 1E-6)
  {
    f=0.0;
    fdr=0.0;
    return;
  }

  v1 = 1.0 / (r - a);
  f = exp(gamma * v1);
  fdr = (-gamma * v1 * v1) * f;
}

/* ----------------------------------------------------------------------
   pre-calculated structures
------------------------------------------------------------------------- */

void PairEDIPMulti::allocatePreLoops(void)
{
  int nthreads = comm->nthreads;

  memory->create(preForceCoord,5*nthreads*leadDimInteractionList,"edip:preForceCoord");
}

/* ----------------------------------------------------------------------
   deallocate preLoops
------------------------------------------------------------------------- */

void PairEDIPMulti::deallocatePreLoops(void)
{
  memory->destroy(preForceCoord);
}

/* ---------------------------------------------------------------------- */

void PairEDIPMulti::allocate()
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

void PairEDIPMulti::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairEDIPMulti::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup();

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  // allocate tables and internal structures

  allocatePreLoops();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEDIPMulti::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style edip/multi requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style edip/multi requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEDIPMulti::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairEDIPMulti::read_file(char *file)
{
  int params_per_line = 20;
  char **words = new char*[params_per_line+1];

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open EDIP potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each set of params from potential file
  // one set of params can span multiple lines
  // store params if all 3 element tags are in element list

  int n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in EDIP potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next entry in file

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].A = atof(words[3]);
    params[nparams].B = atof(words[4]);
    params[nparams].cutoffA = atof(words[5]);
    params[nparams].cutoffC = atof(words[6]);
    params[nparams].alpha = atof(words[7]);
    params[nparams].beta = atof(words[8]);
    params[nparams].eta = atof(words[9]);
    params[nparams].gamma = atof(words[10]);
    params[nparams].lambda = atof(words[11]);
    params[nparams].mu = atof(words[12]);
    params[nparams].rho = atof(words[13]);
    params[nparams].sigma = atof(words[14]);
    params[nparams].Q0 = atof(words[15]);
    params[nparams].u1 = atof(words[16]);
    params[nparams].u2 = atof(words[17]);
    params[nparams].u3 = atof(words[18]);
    params[nparams].u4 = atof(words[19]);

    if (params[nparams].A < 0.0 || params[nparams].B < 0.0 ||
        params[nparams].cutoffA < 0.0 || params[nparams].cutoffC < 0.0 ||
        params[nparams].alpha < 0.0 || params[nparams].beta < 0.0 ||
        params[nparams].eta < 0.0 || params[nparams].gamma < 0.0 ||
        params[nparams].lambda < 0.0 || params[nparams].mu < 0.0 ||
        params[nparams].rho < 0.0 || params[nparams].sigma < 0.0)
      error->all(FLERR,"Illegal EDIP parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairEDIPMulti::setup()
{
  int i,j,k,m,n;
  double rtmp;

  // set elem2param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem2param[i][j][k] = n;
      }

  // set cutoff square

  for (m = 0; m < nparams; m++) {
    params[m].cutsq = params[m].cutoffA*params[m].cutoffA;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = sqrt(params[m].cutsq);
    if (rtmp > cutmax) cutmax = rtmp;
  }

}
