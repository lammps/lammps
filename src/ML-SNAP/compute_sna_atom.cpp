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

#include "compute_sna_atom.h"

#include "sna.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

ComputeSNAAtom::ComputeSNAAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), cutsq(nullptr), list(nullptr), sna(nullptr),
  radelem(nullptr), wjelem(nullptr), sinnerelem(nullptr), dinnerelem(nullptr)

{
  // begin code common to all SNAP computes

  double rfac0, rmin0;
  int twojmax, switchflag, bzeroflag, bnormflag, wselfallflag;

  int ntypes = atom->ntypes;
  int nargmin = 6 + 2 * ntypes;

  if (narg < nargmin) error->all(FLERR, "Illegal compute {} command", style);

  // default values

  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;
  chemflag = 0;
  bnormflag = 0;
  wselfallflag = 0;
  switchinnerflag = 0;
  nelements = 1;
  nnn = 12;
  wmode = 0;
  delta = 1.e-3;
  nearest_neighbors_mode = false;
  // process required arguments

  memory->create(radelem, ntypes + 1, "sna/atom:radelem"); // offset by 1 to match up with types
  memory->create(wjelem, ntypes + 1, "sna/atom:wjelem");

  rcutfac = utils::numeric(FLERR, arg[3], false, lmp);
  rfac0 = utils::numeric(FLERR, arg[4], false, lmp);
  twojmax = utils::inumeric(FLERR, arg[5], false, lmp);

  for (int i = 0; i < ntypes; i++)
    radelem[i + 1] =
        utils::numeric(FLERR, arg[6 + i], false, lmp);
  for (int i = 0; i < ntypes; i++)
    wjelem[i + 1] =
        utils::numeric(FLERR, arg[6 + ntypes + i], false, lmp);

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq, ntypes + 1, ntypes + 1, "sna/atom:cutsq");
  for (int i = 1; i <= ntypes; i++) {
    cut = 2.0 * radelem[i] * rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut * cut;
    for (int j = i + 1; j <= ntypes; j++) {
      cut = (radelem[i] + radelem[j]) * rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut * cut;
    }
  }

  // set local input checks

  int sinnerflag = 0;
  int dinnerflag = 0;

  // process optional args

  int iarg = nargmin;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "rmin0") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      rmin0 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "switchflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      switchflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "bzeroflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      bzeroflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "quadraticflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      quadraticflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"nnn") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      nnn = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      nearest_neighbors_mode = true;
      if (nnn <= 0) error->all(FLERR, "Illegal compute compute {} command", style);
      iarg += 2;
    } else if (strcmp(arg[iarg],"wmode") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      wmode = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (wmode < 0) error->all(FLERR, "Illegal compute compute {} command", style);
      iarg += 2;
    } else if (strcmp(arg[iarg],"delta") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      delta = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (delta < 1.0e-3) error->all(FLERR, "Illegal compute compute {} command", style);
      iarg += 2;
    } else if (strcmp(arg[iarg], "chem") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      chemflag = 1;
      memory->create(map, ntypes + 1, "compute_sna_grid:map");
      nelements = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      for (int i = 0; i < ntypes; i++) {
        int jelem = utils::inumeric(FLERR, arg[iarg + 2 + i], false, lmp);
        if (jelem < 0 || jelem >= nelements) error->all(FLERR, "Illegal compute {} command", style);
        map[i + 1] = jelem;
      }
      iarg += 2 + ntypes;
    } else if (strcmp(arg[iarg], "bnormflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      bnormflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "wselfallflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      wselfallflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "switchinnerflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      switchinnerflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "sinner") == 0) {
      iarg++;
      if (iarg + ntypes > narg) error->all(FLERR, "Illegal compute {} command", style);
      memory->create(sinnerelem, ntypes + 1, "snap:sinnerelem");
      for (int i = 0; i < ntypes; i++)
        sinnerelem[i + 1] = utils::numeric(FLERR, arg[iarg + i], false, lmp);
      sinnerflag = 1;
      iarg += ntypes;
    } else if (strcmp(arg[iarg], "dinner") == 0) {
      iarg++;
      if (iarg + ntypes > narg) error->all(FLERR, "Illegal compute {} command", style);
      memory->create(dinnerelem, ntypes + 1, "snap:dinnerelem");
      for (int i = 0; i < ntypes; i++)
        dinnerelem[i + 1] = utils::numeric(FLERR, arg[iarg + i], false, lmp);
      dinnerflag = 1;
      iarg += ntypes;
    } else
      error->all(FLERR, "Illegal compute {} command", style);
  }

  if (switchinnerflag && !(sinnerflag && dinnerflag))
    error->all(
        FLERR,
        "Illegal compute {} command: switchinnerflag = 1, missing sinner/dinner keyword",
        style);

  if (!switchinnerflag && (sinnerflag || dinnerflag))
    error->all(
        FLERR,
        "Illegal compute {} command: switchinnerflag = 0, unexpected sinner/dinner keyword",
        style);

  snaptr = new SNA(lmp, rfac0, twojmax, rmin0, switchflag, bzeroflag, chemflag, bnormflag,
                   wselfallflag, nelements, switchinnerflag);

  ncoeff = snaptr->ncoeff;
  nvalues = ncoeff;
  if (quadraticflag) nvalues += (ncoeff * (ncoeff + 1)) / 2;

  // end code common to all SNAP computes

  size_peratom_cols = nvalues;
  peratom_flag = 1;

  nmax = 0;
  sna = nullptr;

}

/* ---------------------------------------------------------------------- */

ComputeSNAAtom::~ComputeSNAAtom()
{
  memory->destroy(sna);
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  delete snaptr;

  if (chemflag) memory->destroy(map);

  if (switchinnerflag) {
    memory->destroy(sinnerelem);
    memory->destroy(dinnerelem);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSNAAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute sna/atom requires a pair style be defined");
  rcutsq = force->pair->cutforce * force->pair->cutforce;
  if (cutmax > force->pair->cutforce)
    error->all(FLERR,"Compute sna/atom cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style("sna/atom").size() > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute sna/atom");
  snaptr->init();
}

/* ---------------------------------------------------------------------- */

void ComputeSNAAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSNAAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow sna array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(sna);
    nmax = atom->nmax;
    memory->create(sna,nmax,size_peratom_cols,"sna/atom:sna");
    array_atom = sna;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int* const ilist = list->ilist;
  const int* const numneigh = list->numneigh;
  int** const firstneigh = list->firstneigh;
  int * const type = atom->type;

  // compute sna for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double** const x = atom->x;
  const int* const mask = atom->mask;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {

      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const int itype = type[i];
      int ielem = 0;
      if (chemflag)
        ielem = map[itype];
      const double radi = radelem[itype];
      const int* const jlist = firstneigh[i];
      const int jnum = numneigh[i];


      // ############################################################################## //
      // ##### Start of section for computing bispectrum on nnn nearest neighbors ##### //
      // ############################################################################## //
      if (nearest_neighbors_mode == true) {
        // ##### 1) : consider full neighbor list in rlist
        memory->create(distsq, jnum, "snann/atom:distsq");
        memory->create(rlist, jnum, 3, "snann/atom:rlist");

        int ncount = 0;
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= NEIGHMASK;
          int jtype = type[j];

          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx * delx + dely * dely + delz * delz;

          if (rsq < rcutsq) {
            distsq[ncount] = rsq;
            rlist[ncount][0] = delx;
            rlist[ncount][1] = dely;
            rlist[ncount][2] = delz;
            ncount++;
          }
        }

        // ##### 2) : compute optimal cutoff such that sum weights S_target = nnn
        double S_target=1.*nnn;
        double rc_start=0.1;
        double rc_max=sqrt(rcutsq);
        double tol=1.e-8;
        double * sol_dich = dichotomie(S_target, rc_start, rc_max, tol, distsq, ncount, wmode, delta);
        memory->destroy(distsq);

        // ##### 3) : assign that optimal cutoff radius to bispectrum context using rcsol
        double rcsol = (sol_dich[0]+sol_dich[1])/2.;
        memory->destroy(sol_dich);
        snaptr->grow_rij(ncount);

        int ninside = 0;
        for (int jj = 0; jj < ncount; jj++) {
          int j = jlist[jj];
          j &= NEIGHMASK;

          const double rsq = rlist[jj][0]*rlist[jj][0]+rlist[jj][1]*rlist[jj][1]+rlist[jj][2]*rlist[jj][2];
          int jtype = type[j];
          int jelem = 0;
          if (chemflag)
            jelem = map[jtype];

          if (rsq < rcsol*rcsol) {
            snaptr->rij[ninside][0] = rlist[jj][0];//rijmax;
            snaptr->rij[ninside][1] = rlist[jj][1];//rijmax;
            snaptr->rij[ninside][2] = rlist[jj][2];//rijmax;
            snaptr->inside[ninside] = j;
            snaptr->wj[ninside] = 1.;
            snaptr->rcutij[ninside] = rcsol;

            if (switchinnerflag) {
              snaptr->sinnerij[ninside] = 0.5*(sinnerelem[itype]+sinnerelem[jtype]);
              snaptr->dinnerij[ninside] = 0.5*(dinnerelem[itype]+dinnerelem[jtype]);
            }
            if (chemflag) snaptr->element[ninside] = jelem;
            ninside++;
          }
        }

        memory->destroy(rlist);

        // ############################################################################ //
        // ##### End of section for computing bispectrum on nnn nearest neighbors ##### //
        // ############################################################################ //
        snaptr->compute_ui(ninside, ielem);
        snaptr->compute_zi();
        snaptr->compute_bi(ielem);

        for (int icoeff = 0; icoeff < ncoeff; icoeff++)
          sna[i][icoeff] = snaptr->blist[icoeff];
        if (quadraticflag) {
          int ncount = ncoeff;
          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
            double bi = snaptr->blist[icoeff];

            // diagonal element of quadratic matrix

            sna[i][ncount++] = 0.5*bi*bi;

            // upper-triangular elements of quadratic matrix

            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++)
              sna[i][ncount++] = bi*snaptr->blist[jcoeff];
          }
        }

      } else {
        // ensure rij, inside, and typej  are of size jnum

        snaptr->grow_rij(jnum);

        // rij[][3] = displacements between atom I and those neighbors
        // inside = indices of neighbors of I within cutoff
        // typej = types of neighbors of I within cutoff

        int ninside = 0;
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= NEIGHMASK;

          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx*delx + dely*dely + delz*delz;
          int jtype = type[j];
          int jelem = 0;
          if (chemflag)
            jelem = map[jtype];
          if (rsq < cutsq[itype][jtype] && rsq>1e-20) {
            snaptr->rij[ninside][0] = delx;
            snaptr->rij[ninside][1] = dely;
            snaptr->rij[ninside][2] = delz;
            snaptr->inside[ninside] = j;
            snaptr->wj[ninside] = wjelem[jtype];
            snaptr->rcutij[ninside] = (radi+radelem[jtype])*rcutfac;
            if (switchinnerflag) {
              snaptr->sinnerij[ninside] = 0.5*(sinnerelem[itype]+sinnerelem[jtype]);
              snaptr->dinnerij[ninside] = 0.5*(dinnerelem[itype]+dinnerelem[jtype]);
            }
            if (chemflag) snaptr->element[ninside] = jelem;
            ninside++;
          }
        }
        snaptr->compute_ui(ninside, ielem);
        snaptr->compute_zi();
        snaptr->compute_bi(ielem);

        for (int icoeff = 0; icoeff < ncoeff; icoeff++)
          sna[i][icoeff] = snaptr->blist[icoeff];
        if (quadraticflag) {
          int ncount = ncoeff;
          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
            double bi = snaptr->blist[icoeff];

            // diagonal element of quadratic matrix

            sna[i][ncount++] = 0.5*bi*bi;

            // upper-triangular elements of quadratic matrix

            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++)
              sna[i][ncount++] = bi*snaptr->blist[jcoeff];
          }
        }

      }

    } else {
      for (int icoeff = 0; icoeff < size_peratom_cols; icoeff++)
        sna[i][icoeff] = 0.0;
    }
  }

}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeSNAAtom::memory_usage()
{
  double bytes = (double)nmax*size_peratom_cols * sizeof(double); // sna
  bytes += snaptr->memory_usage();                        // SNA object

  return bytes;
}

/* ----------------------------------------------------------------------
   select3 routine from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   sort auxiliary arrays at same time
------------------------------------------------------------------------- */

// Use no-op do while to create single statement

#define SWAP(a, b) \
  do {             \
    tmp = a;       \
    (a) = b;       \
    (b) = tmp;     \
  } while (0)

#define ISWAP(a, b) \
  do {              \
    itmp = a;       \
    (a) = b;        \
    (b) = itmp;     \
  } while (0)

#define SWAP3(a, b)  \
  do {               \
    tmp = (a)[0];    \
    (a)[0] = (b)[0]; \
    (b)[0] = tmp;    \
    tmp = (a)[1];    \
    (a)[1] = (b)[1]; \
    (b)[1] = tmp;    \
    tmp = (a)[2];    \
    (a)[2] = (b)[2]; \
    (b)[2] = tmp;    \
  } while (0)

/* ---------------------------------------------------------------------- */

void ComputeSNAAtom::select3(int k, int n, double *arr, int *iarr, double **arr3)
{
  int i, ir, j, l, mid, ia, itmp;
  double a, tmp, a3[3];

  arr--;
  iarr--;
  arr3--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l + 1) {
      if (ir == l + 1 && arr[ir] < arr[l]) {
        SWAP(arr[l], arr[ir]);
        ISWAP(iarr[l], iarr[ir]);
        SWAP3(arr3[l], arr3[ir]);
      }
      return;
    } else {
      mid = (l + ir) >> 1;
      SWAP(arr[mid], arr[l + 1]);
      ISWAP(iarr[mid], iarr[l + 1]);
      SWAP3(arr3[mid], arr3[l + 1]);
      if (arr[l] > arr[ir]) {
        SWAP(arr[l], arr[ir]);
        ISWAP(iarr[l], iarr[ir]);
        SWAP3(arr3[l], arr3[ir]);
      }
      if (arr[l + 1] > arr[ir]) {
        SWAP(arr[l + 1], arr[ir]);
        ISWAP(iarr[l + 1], iarr[ir]);
        SWAP3(arr3[l + 1], arr3[ir]);
      }
      if (arr[l] > arr[l + 1]) {
        SWAP(arr[l], arr[l + 1]);
        ISWAP(iarr[l], iarr[l + 1]);
        SWAP3(arr3[l], arr3[l + 1]);
      }
      i = l + 1;
      j = ir;
      a = arr[l + 1];
      ia = iarr[l + 1];
      a3[0] = arr3[l + 1][0];
      a3[1] = arr3[l + 1][1];
      a3[2] = arr3[l + 1][2];
      for (;;) {
        do i++;
        while (arr[i] < a);
        do j--;
        while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i], arr[j]);
        ISWAP(iarr[i], iarr[j]);
        SWAP3(arr3[i], arr3[j]);
      }
      arr[l + 1] = arr[j];
      arr[j] = a;
      iarr[l + 1] = iarr[j];
      iarr[j] = ia;
      arr3[l + 1][0] = arr3[j][0];
      arr3[l + 1][1] = arr3[j][1];
      arr3[l + 1][2] = arr3[j][2];
      arr3[j][0] = a3[0];
      arr3[j][1] = a3[1];
      arr3[j][2] = a3[2];
      if (j >= k) ir = j - 1;
      if (j <= k) l = i;
    }
  }
}

double * ComputeSNAAtom::weights(double * rsq, double rcut, int ncounts)
{
  double * w=nullptr;
  memory->destroy(w);
  memory->create(w, ncounts, "snann:gauss_weights");
  double rloc=0.;
  for (int i=0; i<ncounts; i++)
    {
      rloc = sqrt(rsq[i]);
      if (rloc > rcut){
        w[i]=0.;
      } else {
        w[i]=1.;
      }
    }
  return w;
}

double * ComputeSNAAtom::tanh_weights(double * rsq, double rcut, double delta, int ncounts)
{
  double * w=nullptr;
  memory->destroy(w);
  memory->create(w, ncounts, "snann:gauss_weights");
  double rloc=0.;

  for (int i=0; i<ncounts; i++)
    {
      rloc = sqrt(rsq[i]);
      w[i] = 0.5*(1.-tanh((rloc-rcut)/delta));
    }
  return w;
}

double ComputeSNAAtom::sum_weights(double * rsq, double * w, int ncounts)
{
  double S=0.;
  double rloc=0.;
  for (int i=0; i<ncounts; i++)
    {
      S += w[i];
    }
  return S;
}

double ComputeSNAAtom::get_target_rcut(double S_target, double * rsq, double rcut, int ncounts, int weightmode, double delta)
{
  double S_sol;
  if (weightmode == 0) {
    double * www = weights(rsq, rcut, ncounts);
    S_sol = sum_weights(rsq, www, ncounts);
    memory->destroy(www);
  } else if (weightmode == 1) {
    double * www = tanh_weights(rsq, rcut, delta, ncounts);
    S_sol = sum_weights(rsq, www, ncounts);
    memory->destroy(www);
  }
  double err = S_sol - S_target;
  return err;
}

double * ComputeSNAAtom::dichotomie(double S_target, double a, double b, double e, double * rsq, int ncounts, int weightmode, double delta)
{

  double d=b-a;
  double * sol = nullptr;
  memory->destroy(sol);
  memory->create(sol, 2, "snann:sol");
  double m=0.;

  int cnt=0;
  do
    {
      m = ( a + b ) / 2.;
      d = fabs( b - a );
      double f_ra = get_target_rcut(S_target, rsq, a, ncounts, weightmode, delta);
      double f_rm = get_target_rcut(S_target, rsq, m, ncounts, weightmode, delta);
      if (f_rm == 0.)
        {
          sol[0]=m;
          sol[1]=m;
          return sol;
        }
      else if (f_rm*f_ra > 0.)
        {
          a = m;
        }
      else
        {
          b = m;
        }
      cnt+=1;
    } while ( d > e );
  sol[0]=a;
  sol[1]=b;
  return sol;
}
