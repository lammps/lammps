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
   Contributing authors: Aidan P. Thompson (SNL)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "compute_adf.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum{DEGREE, RADIAN, COSINE};

/* ----------------------------------------------------------------------
   compute angular distribution functions for I, J, K atoms
 ---------------------------------------------------------------------- */

ComputeADF::ComputeADF(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  ilo(NULL), ihi(NULL), jlo(NULL), jhi(NULL), klo(NULL), khi(NULL),
  hist(NULL), histall(NULL),
  rcutinnerj(NULL), rcutinnerk(NULL),
  rcutouterj(NULL), rcutouterk(NULL),
  list(NULL),
  iatomcount(NULL), iatomcountall(NULL), iatomflag(NULL),
  maxjatom(NULL), maxkatom(NULL),
  numjatom(NULL), numkatom(NULL),
  neighjatom(NULL),neighkatom(NULL),
  jatomflag(NULL), katomflag(NULL),
  maxjkatom(NULL), numjkatom(NULL),
  neighjkatom(NULL), bothjkatom(NULL), delrjkatom(NULL)
{
  int nargsperadf = 7;

  if (narg < 4 ) error->all(FLERR,"Illegal compute adf command");

  array_flag = 1;
  extarray = 0;

  // default values

  ordinate_style = DEGREE;
  cutflag = 0;

  nbin = force->inumeric(FLERR,arg[3]);
  if (nbin < 1) error->all(FLERR,"Illegal compute adf command");

  // optional args
  // nargtriple = # of triplewise args, starting at iarg = 4

  int iarg;
  for (iarg = 4; iarg < narg; iarg++)
    if (strcmp(arg[iarg],"ordinate") == 0) break;

  int nargtriple = iarg - 4;

  //  loop over keywords

  while (iarg < narg) {
    if (strcmp(arg[iarg],"ordinate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute adf command");
      if (strcmp(arg[iarg+1],"degree") == 0) ordinate_style = DEGREE;
      else if (strcmp(arg[iarg+1],"radian") == 0) ordinate_style = RADIAN;
      else if (strcmp(arg[iarg+1],"cosine") == 0) ordinate_style = COSINE;
      else error->all(FLERR,"Illegal compute adf command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute adf command");
  }

  // triplewise args

  if (!nargtriple) ntriples = 1;
  else {
    if (nargtriple % nargsperadf)
      error->all(FLERR,"Illegal compute adf command");
    ntriples = nargtriple/nargsperadf;
  }

  size_array_rows = nbin;
  size_array_cols = 1 + 2*ntriples;

  int ntypes = atom->ntypes;
  memory->create(iatomflag,ntriples,ntypes+1,
                 "adf:iatomflag");
  memory->create(jatomflag,ntriples,ntypes+1,
                 "adf:jatomflag");
  memory->create(katomflag,ntriples,ntypes+1,
                 "adf:katomflag");

  ilo = new int[ntriples];
  ihi = new int[ntriples];
  jlo = new int[ntriples];
  jhi = new int[ntriples];
  klo = new int[ntriples];
  khi = new int[ntriples];
  iatomcount = new int[ntriples];
  iatomcountall = new int[ntriples];

  rcutinnerj = new double[ntriples];
  rcutouterj = new double[ntriples];
  rcutinnerk = new double[ntriples];
  rcutouterk = new double[ntriples];

  if (!nargtriple) {
    ilo[0] = 1; ihi[0] = ntypes;
    jlo[0] = 1; jhi[0] = ntypes;
    klo[0] = 1; khi[0] = ntypes;
  } else {
    cutflag = 1;
    iarg = 4;
    for (int m = 0; m < ntriples; m++) {
      force->bounds(FLERR,arg[iarg],atom->ntypes,ilo[m],ihi[m]);
      force->bounds(FLERR,arg[iarg+1],atom->ntypes,jlo[m],jhi[m]);
      force->bounds(FLERR,arg[iarg+2],atom->ntypes,klo[m],khi[m]);
      if (ilo[m] > ihi[m] ||
          jlo[m] > jhi[m] ||
          klo[m] > khi[m])
        error->all(FLERR,"Illegal compute adf command");
      rcutinnerj[m] = force->numeric(FLERR,arg[iarg+3]);
      rcutouterj[m] = force->numeric(FLERR,arg[iarg+4]);
      if (rcutinnerj[m] < 0.0 || rcutinnerj[m] >= rcutouterj[m])
        error->all(FLERR,"Illegal compute adf command");
      rcutinnerk[m] = force->numeric(FLERR,arg[iarg+5]);
      rcutouterk[m] = force->numeric(FLERR,arg[iarg+6]);
      if (rcutinnerk[m] < 0.0 || rcutinnerk[m] >= rcutouterk[m])
        error->all(FLERR,"Illegal compute adf command");
      iarg += nargsperadf;
    }
  }

  // identify central atom types

  int i,j,k;

  for (int m = 0; m < ntriples; m++) {
    for (i = 1; i <= ntypes; i++)
      iatomflag[m][i] = 0;
    for (i = ilo[m]; i <= ihi[m]; i++)
      iatomflag[m][i] = 1;
  }

    // identify j atom types

  for (int m = 0; m < ntriples; m++) {
    for (j = 1; j <= ntypes; j++)
      jatomflag[m][j] = 0;
    for (j = jlo[m]; j <= jhi[m]; j++)
      jatomflag[m][j] = 1;
  }

    // identify k atom types

  for (int m = 0; m < ntriples; m++) {
    for (k = 1; k <= ntypes; k++)
      katomflag[m][k] = 0;
    for (k = klo[m]; k <= khi[m]; k++)
      katomflag[m][k] = 1;
  }

  memory->create(hist,ntriples,nbin,"adf:hist");
  memory->create(histall,ntriples,nbin,"adf:histall");
  memory->create(array,nbin,1+2*ntriples,"adf:array");

  // initialize list of jatom neighbor lists

  maxjatom = new int[ntriples];
  numjatom = new int[ntriples];
  neighjatom = new int*[ntriples];
  for (int m = 0; m < ntriples; m++) {
    maxjatom[m] = 10;
    memory->create(neighjatom[m],maxjatom[m],"adf:neighjatom");
  }

  // initialize list of katom neighbor lists

  maxkatom = new int[ntriples];
  numkatom = new int[ntriples];
  neighkatom = new int*[ntriples];
  for (int m = 0; m < ntriples; m++) {
    maxkatom[m] = 10;
    memory->create(neighkatom[m],maxkatom[m],"adf:neighkatom");
  }

  // initialize list of short neighbor lists

  maxjkatom = new int[ntriples];
  numjkatom = new int[ntriples];
  neighjkatom = new int*[ntriples];
  bothjkatom = new int*[ntriples];
  delrjkatom = new double**[ntriples];
  for (int m = 0; m < ntriples; m++) {
    maxjkatom[m] = 10;
    memory->create(neighjkatom[m],maxjkatom[m],"adf:neighjkatom");
    memory->create(bothjkatom[m],maxjkatom[m],"adf:bothjkatom");
    memory->create(delrjkatom[m],maxjkatom[m],4,"adf:delrjkatom");
  }

  rad2deg = 180.0 / MY_PI;
}

/* ---------------------------------------------------------------------- */

ComputeADF::~ComputeADF()
{
  memory->destroy(iatomflag);
  delete [] ilo;
  delete [] ihi;
  delete [] jlo;
  delete [] jhi;
  delete [] klo;
  delete [] khi;
  delete [] iatomcount;
  delete [] iatomcountall;
  memory->destroy(hist);
  memory->destroy(histall);
  memory->destroy(array);

  memory->destroy(jatomflag);
  delete [] rcutinnerj;
  delete [] rcutouterj;
  delete [] maxjatom;
  delete [] numjatom;
  for (int m = 0; m < ntriples; m++)
    memory->destroy(neighjatom[m]);
  delete [] neighjatom;

  memory->destroy(katomflag);
  delete [] rcutinnerk;
  delete [] rcutouterk;
  delete [] maxkatom;
  delete [] numkatom;
  for (int m = 0; m < ntriples; m++)
    memory->destroy(neighkatom[m]);
  delete [] neighkatom;

  delete [] maxjkatom;
  delete [] numjkatom;
  for (int m = 0; m < ntriples; m++)
    memory->destroy(neighjkatom[m]);
  delete [] neighjkatom;
  for (int m = 0; m < ntriples; m++)
    memory->destroy(bothjkatom[m]);
  delete [] bothjkatom;
  for (int m = 0; m < ntriples; m++)
    memory->destroy(delrjkatom[m]);
  delete [] delrjkatom;
}

/* ---------------------------------------------------------------------- */

void ComputeADF::init()
{
  double mycutneigh = 0.0;
  double maxouter = 0.0;

  if (!cutflag) {
    if (!force->pair)
      error->all(FLERR,"Compute adf requires a pair style be defined "
                 "or an outer cutoff specified");
    rcutinnerj[0] = 0.0;
    rcutinnerk[0] = 0.0;
    rcutouterj[0] = force->pair->cutforce;
    rcutouterk[0] = force->pair->cutforce;
    maxouter = force->pair->cutforce;;
  } else {
    for (int m = 0; m < ntriples; m++) {
      maxouter = MAX(rcutouterj[m],maxouter);
      maxouter = MAX(rcutouterk[m],maxouter);
    }
  }

  // specify mycutneigh if force cutoff too small or non-existent

  if (!(force->pair) || maxouter > force->pair->cutforce) {
    double skin = neighbor->skin;
    mycutneigh = maxouter + skin;
    if (mycutneigh > comm->cutghostuser)
      error->all(FLERR,"Compute adf outer cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");
  }

  // assign ordinate values to 1st column of output array

  int x0;
  if (ordinate_style == DEGREE) {
    deltax = MY_PI / nbin * rad2deg;
    deltaxinv = nbin / MY_PI;
    x0 = 0.0;

  } else if (ordinate_style == RADIAN) {
    deltax = MY_PI / nbin;
    deltaxinv = nbin / MY_PI;
    x0 = 0.0;

  } else if (ordinate_style == COSINE) {
    deltax = 2.0 / nbin;
    deltaxinv = 1.0/deltax;
    x0 = -1.0;
  }

  for (int i = 0; i < nbin; i++)
    array[i][0] = x0 + (i+0.5) * deltax;

  // need an occasional full neighbor list
  // if mycutneigh specified, request a cutoff = maxouter + skin
  // skin is included b/c Neighbor uses this value similar
  //   to its cutneighmax = force cutoff + skin
  // also, this NeighList may be used by this compute for multiple steps
  //   (until next reneighbor), so it needs to contain atoms further
  //   than maxouter apart, just like a normal neighbor list does

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
  if (mycutneigh > 0.0) {
    neighbor->requests[irequest]->cut = 1;
    neighbor->requests[irequest]->cutoff = mycutneigh;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeADF::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeADF::compute_array()
{
  int i,j,k,m,ii,jj,jatom,katom,jk,jjj,kkk;
  int inum,jnum,itype,jtype,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double factor_lj,factor_coul;
  double delr1[3],delr2[3],rinv1,rinv2,rinv12,cs,theta;

  invoked_array = update->ntimestep;

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero the histogram counts

  for (m = 0; m < ntriples; m++)
    for (ibin = 0; ibin < nbin; ibin++)
      hist[m][ibin] = 0.0;

  // zero the central atom counts

  for (m = 0; m < ntriples; m++)
    iatomcount[m] = 0;

  // tally the ADFs
  // all three atoms i, j, and k must be in fix group
  // tally I,J,K triple only if I is central atom
  // and J,K matches unordered neighbor types (JJ,KK)

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;

  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // count atom i in each matching ADF
    // zero the jatom, katom neighbor list counts

    for (m = 0; m < ntriples; m++) {
      if (iatomflag[m][itype]) {
        iatomcount[m]++;
      }
      numjatom[m] = 0;
      numkatom[m] = 0;
      numjkatom[m] = 0;
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      // if both weighting factors are 0, skip this pair
      // could be 0 and still be in neigh list for long-range Coulombics
      // want consistency with non-charged triples which wouldn't be in list

      if (factor_lj == 0.0 && factor_coul == 0.0) continue;

      if (!(mask[j] & groupbit)) continue;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      for (m = 0; m < ntriples; m++) {

        // check if itype, jtype, ktype match this ADF definition
        // if yes, add j to jatom, katom, and jkatom lists

        if (!iatomflag[m][itype]) continue;

        int jflag = 0;
        if (jatomflag[m][jtype] &&
            rsq >= rcutinnerj[m]*rcutinnerj[m] &&
            rsq <= rcutouterj[m]*rcutouterj[m]) {
          jflag = 1;
          jatom = numjatom[m]++;
          neighjatom[m][jatom] = numjkatom[m];
          if (numjatom[m] >= maxjatom[m]) {
            maxjatom[m] += maxjatom[m]/2;
            memory->grow(neighjatom[m],maxjatom[m],"adf:neighjatom");
          }
        }

        int kflag = 0;
        if (katomflag[m][jtype] &&
            rsq >= rcutinnerk[m]*rcutinnerk[m] &&
            rsq <= rcutouterk[m]*rcutouterk[m]) {
          kflag = 1;
          katom = numkatom[m]++;
          neighkatom[m][katom] = numjkatom[m];
          if (numkatom[m] >= maxkatom[m]) {
            maxkatom[m] += maxkatom[m]/2;
            memory->grow(neighkatom[m],maxkatom[m],"adf:neighkatom");
          }
        }

        // if atom in either list, add to jklist

        if (jflag || kflag) {
          jk = numjkatom[m]++;
          neighjkatom[m][jk] = j;
          delrjkatom[m][jk][0] = delx;
          delrjkatom[m][jk][1] = dely;
          delrjkatom[m][jk][2] = delz;
          delrjkatom[m][jk][3] = 1.0/sqrt(rsq);
          if (numjkatom[m] >= maxjkatom[m]) {
            maxjkatom[m] += maxjkatom[m]/2;
            memory->grow(neighjkatom[m],maxjkatom[m],"adf:neighjkatom");
            memory->grow(bothjkatom[m],maxjkatom[m],"adf:bothjkatom");
            memory->grow(delrjkatom[m],maxjkatom[m],4,"adf:delrjkatom");
          }

          // indicate if atom in both lists

          if (jflag && kflag)
            bothjkatom[m][jk] = 1;
          else
            bothjkatom[m][jk] = 0;
        }
      }
    }

    // loop over ADFs

    for (m = 0; m < ntriples; m++) {

      // loop over (j,k) pairs

      for (jatom = 0; jatom < numjatom[m]; jatom++) {
        jjj = neighjatom[m][jatom];
        j = neighjkatom[m][jjj];

        delr1[0] = delrjkatom[m][jjj][0];
        delr1[1] = delrjkatom[m][jjj][1];
        delr1[2] = delrjkatom[m][jjj][2];
        rinv1 = delrjkatom[m][jjj][3];

        for (katom = 0; katom < numkatom[m]; katom++) {
          kkk = neighkatom[m][katom];
          k = neighjkatom[m][kkk];

          // skip if j==k, or j > k and both are in both lists

          if (k == j || (j > k && bothjkatom[m][jjj] && bothjkatom[m][kkk])) continue;

          delr2[0] = delrjkatom[m][kkk][0];
          delr2[1] = delrjkatom[m][kkk][1];
          delr2[2] = delrjkatom[m][kkk][2];
          rinv2 = delrjkatom[m][kkk][3];
          rinv12 = rinv1*rinv2;
          cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;

          if (ordinate_style == COSINE) {
            ibin = static_cast<int> ((cs+1.0)*deltaxinv);
            if (ibin >= nbin || ibin < 0) continue;
            hist[m][ibin] += 1.0;
          } else {
            if (cs > 1.0) cs = 1.0;
            if (cs < -1.0) cs = -1.0;
            theta = acos(cs);
            ibin = static_cast<int> (theta*deltaxinv);
            if (ibin >= nbin || ibin < 0) continue;
            hist[m][ibin] += 1.0;
          }
        }
      }
    }
  }

  // sum histograms across procs

  MPI_Allreduce(hist[0],histall[0],ntriples*nbin,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(iatomcount,iatomcountall,ntriples,MPI_INT,MPI_SUM,world);

  // convert counts to pdf(theta) and adf(theta)
  // copy into output array

  for (m = 0; m < ntriples; m++) {

    double count = 0;
    for (ibin = 0; ibin < nbin; ibin++)
      count += histall[m][ibin];

    double normfac1, pdftheta, normfac2, adftheta;
    if (count > 0.0) normfac1 = 1.0/deltax/count;
    else normfac1 = 0.0;
    if (iatomcountall[m] > 0.0) normfac2 = 1.0/iatomcountall[m];
    else normfac2 = 0.0;
    adftheta = 0.0;
    for (ibin = 0; ibin < nbin; ibin++) {
      pdftheta = histall[m][ibin] * normfac1;
      adftheta += histall[m][ibin] * normfac2;
      array[ibin][1+2*m] = pdftheta;
      array[ibin][2+2*m] = adftheta;
    }
  }

}
