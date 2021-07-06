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

#include "pair_snap.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "sna.h"
#include "tokenizer.h"

#include <cmath>

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

PairSNAP::PairSNAP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;

  radelem = nullptr;
  wjelem = nullptr;
  coeffelem = nullptr;

  beta_max = 0;
  beta = nullptr;
  bispectrum = nullptr;
  snaptr = nullptr;
}

/* ---------------------------------------------------------------------- */

PairSNAP::~PairSNAP()
{
  if (copymode) return;

  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(coeffelem);

  memory->destroy(beta);
  memory->destroy(bispectrum);

  delete snaptr;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
  }

}

/* ----------------------------------------------------------------------
   This version is a straightforward implementation
   ---------------------------------------------------------------------- */

void PairSNAP::compute(int eflag, int vflag)
{
  int i,j,jnum,ninside;
  double delx,dely,delz,evdwl,rsq;
  double fij[3];
  int *jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  if (beta_max < list->inum) {
    memory->grow(beta,list->inum,ncoeff,"PairSNAP:beta");
    memory->grow(bispectrum,list->inum,ncoeff,"PairSNAP:bispectrum");
    beta_max = list->inum;
  }

  // compute dE_i/dB_i = beta_i for all i in list

  if (quadraticflag || eflag)
    compute_bispectrum();
  compute_beta();

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int ielem = map[itype];
    const double radi = radelem[ielem];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // insure rij, inside, wj, and rcutij are of size jnum

    snaptr->grow_rij(jnum);

    // rij[][3] = displacements between atom I and those neighbors
    // inside = indices of neighbors of I within cutoff
    // wj = weights for neighbors of I within cutoff
    // rcutij = cutoffs for neighbors of I within cutoff
    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      int jelem = map[jtype];

      if (rsq < cutsq[itype][jtype]&&rsq>1e-20) {
        snaptr->rij[ninside][0] = delx;
        snaptr->rij[ninside][1] = dely;
        snaptr->rij[ninside][2] = delz;
        snaptr->inside[ninside] = j;
        snaptr->wj[ninside] = wjelem[jelem];
        snaptr->rcutij[ninside] = (radi + radelem[jelem])*rcutfac;
        snaptr->element[ninside] = jelem;
        ninside++;
      }
    }

    // compute Ui, Yi for atom I

    if (chemflag)
      snaptr->compute_ui(ninside, ielem);
    else
      snaptr->compute_ui(ninside, 0);

    // for neighbors of I within cutoff:
    // compute Fij = dEi/dRj = -dEi/dRi
    // add to Fi, subtract from Fj
    // scaling is that for type I

    snaptr->compute_yi(beta[ii]);

    for (int jj = 0; jj < ninside; jj++) {
      int j = snaptr->inside[jj];
      if (chemflag)
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, snaptr->element[jj]);
      else
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, 0);

      snaptr->compute_deidrj(fij);

      f[i][0] += fij[0]*scale[itype][itype];
      f[i][1] += fij[1]*scale[itype][itype];
      f[i][2] += fij[2]*scale[itype][itype];
      f[j][0] -= fij[0]*scale[itype][itype];
      f[j][1] -= fij[1]*scale[itype][itype];
      f[j][2] -= fij[2]*scale[itype][itype];

      // tally per-atom virial contribution

      if (vflag)
        ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,
                     fij[0],fij[1],fij[2],
                     -snaptr->rij[jj][0],-snaptr->rij[jj][1],
                     -snaptr->rij[jj][2]);
    }

    // tally energy contribution

    if (eflag) {

      // evdwl = energy of atom I, sum over coeffs_k * Bi_k

      double* coeffi = coeffelem[ielem];
      evdwl = coeffi[0];
      // snaptr->copy_bi2bvec();

      // E = beta.B + 0.5*B^t.alpha.B

      // linear contributions

      for (int icoeff = 0; icoeff < ncoeff; icoeff++)
        evdwl += coeffi[icoeff+1]*bispectrum[ii][icoeff];

      // quadratic contributions

      if (quadraticflag) {
        int k = ncoeff+1;
        for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
          double bveci = bispectrum[ii][icoeff];
          evdwl += 0.5*coeffi[k++]*bveci*bveci;
          for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
            double bvecj = bispectrum[ii][jcoeff];
            evdwl += coeffi[k++]*bveci*bvecj;
          }
        }
      }
      evdwl *= scale[itype][itype];
      ev_tally_full(i,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);
    }

  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   compute beta
------------------------------------------------------------------------- */

void PairSNAP::compute_beta()
{
  int i;
  int *type = atom->type;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];
    const int itype = type[i];
    const int ielem = map[itype];
    double* coeffi = coeffelem[ielem];

    for (int icoeff = 0; icoeff < ncoeff; icoeff++)
      beta[ii][icoeff] = coeffi[icoeff+1];

    if (quadraticflag) {
      int k = ncoeff+1;
      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        double bveci = bispectrum[ii][icoeff];
        beta[ii][icoeff] += coeffi[k]*bveci;
        k++;
        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
          double bvecj = bispectrum[ii][jcoeff];
          beta[ii][icoeff] += coeffi[k]*bvecj;
          beta[ii][jcoeff] += coeffi[k]*bveci;
          k++;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute bispectrum
------------------------------------------------------------------------- */

void PairSNAP::compute_bispectrum()
{
  int i,j,jnum,ninside;
  double delx,dely,delz,rsq;
  int *jlist;

  double **x = atom->x;
  int *type = atom->type;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int ielem = map[itype];
    const double radi = radelem[ielem];

    jlist = list->firstneigh[i];
    jnum = list->numneigh[i];

    // insure rij, inside, wj, and rcutij are of size jnum

    snaptr->grow_rij(jnum);

    // rij[][3] = displacements between atom I and those neighbors
    // inside = indices of neighbors of I within cutoff
    // wj = weights for neighbors of I within cutoff
    // rcutij = cutoffs for neighbors of I within cutoff
    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      int jelem = map[jtype];

      if (rsq < cutsq[itype][jtype]&&rsq>1e-20) {
        snaptr->rij[ninside][0] = delx;
        snaptr->rij[ninside][1] = dely;
        snaptr->rij[ninside][2] = delz;
        snaptr->inside[ninside] = j;
        snaptr->wj[ninside] = wjelem[jelem];
        snaptr->rcutij[ninside] = (radi + radelem[jelem])*rcutfac;
        snaptr->element[ninside] = jelem;
        ninside++;
      }
    }

    if (chemflag)
      snaptr->compute_ui(ninside, ielem);
    else
      snaptr->compute_ui(ninside, 0);
    snaptr->compute_zi();
    if (chemflag)
      snaptr->compute_bi(ielem);
    else
      snaptr->compute_bi(0);

    for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
      bispectrum[ii][icoeff] = snaptr->blist[icoeff];
    }
  }

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSNAP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(scale,n+1,n+1,"pair:scale");
  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSNAP::settings(int narg, char ** /* arg */)
{
  if (narg > 0)
    error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSNAP::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  if (narg != 4 + atom->ntypes) error->all(FLERR,"Incorrect args for pair coefficients");

  map_element2type(narg-4,arg+4);

  // read snapcoeff and snapparam files

  read_files(arg[2],arg[3]);

  if (!quadraticflag)
    ncoeff = ncoeffall - 1;
  else {

    // ncoeffall should be (ncoeff+2)*(ncoeff+1)/2
    // so, ncoeff = floor(sqrt(2*ncoeffall))-1

    ncoeff = sqrt(2*ncoeffall)-1;
    ncoeffq = (ncoeff*(ncoeff+1))/2;
    int ntmp = 1+ncoeff+ncoeffq;
    if (ntmp != ncoeffall) {
      error->all(FLERR,"Incorrect SNAP coeff file");
    }
  }

  snaptr = new SNA(lmp, rfac0, twojmax,
                   rmin0, switchflag, bzeroflag,
                   chemflag, bnormflag, wselfallflag, nelements);

  if (ncoeff != snaptr->ncoeff) {
    if (comm->me == 0)
      printf("ncoeff = %d snancoeff = %d \n",ncoeff,snaptr->ncoeff);
    error->all(FLERR,"Incorrect SNAP parameter file");
  }

  // Calculate maximum cutoff for all elements
  rcutmax = 0.0;
  for (int ielem = 0; ielem < nelements; ielem++)
    rcutmax = MAX(2.0*radelem[ielem]*rcutfac,rcutmax);

  // set default scaling
  int n = atom->ntypes;
  for (int ii = 0; ii < n+1; ii++)
    for (int jj = 0; jj < n+1; jj++)
      scale[ii][jj] = 1.0;

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSNAP::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style SNAP requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  snaptr->init();

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSNAP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  scale[j][i] = scale[i][j];
  return (radelem[map[i]] +
          radelem[map[j]])*rcutfac;
}

/* ---------------------------------------------------------------------- */

void PairSNAP::read_files(char *coefffilename, char *paramfilename)
{

  // open SNAP coefficient file on proc 0

  FILE *fpcoeff;
  if (comm->me == 0) {
    fpcoeff = utils::open_potential(coefffilename,lmp,nullptr);
    if (fpcoeff == nullptr)
      error->one(FLERR,"Cannot open SNAP coefficient file {}: ",
                                   coefffilename, utils::getsyserror());
  }

  char line[MAXLINE],*ptr;
  int eof = 0;
  int nwords = 0;
  while (nwords == 0) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpcoeff);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fpcoeff);
      }
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(line,MAXLINE,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    nwords = utils::count_words(utils::trim_comment(line));
  }
  if (nwords != 2)
    error->all(FLERR,"Incorrect format in SNAP coefficient file");

  // strip single and double quotes from words

  int nelemtmp = 0;
  try {
    ValueTokenizer words(utils::trim_comment(line),"\"' \t\n\r\f");
    nelemtmp = words.next_int();
    ncoeffall = words.next_int();
  } catch (TokenizerException &e) {
    error->all(FLERR,"Incorrect format in SNAP coefficient "
                                 "file: {}", e.what());
  }

  // clean out old arrays and set up element lists

  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(coeffelem);
  memory->create(radelem,nelements,"pair:radelem");
  memory->create(wjelem,nelements,"pair:wjelem");
  memory->create(coeffelem,nelements,ncoeffall,"pair:coeffelem");

  // initialize checklist for all required nelements

  int *elementflags = new int[nelements];
  for (int jelem = 0; jelem < nelements; jelem++)
      elementflags[jelem] = 0;

  // loop over nelemtmp blocks in the SNAP coefficient file

  for (int ielem = 0; ielem < nelemtmp; ielem++) {

    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpcoeff);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fpcoeff);
      }
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof)
      error->all(FLERR,"Incorrect format in SNAP coefficient file");
    MPI_Bcast(line,MAXLINE,MPI_CHAR,0,world);

    std::vector<std::string> words;
    try {
      words = Tokenizer(utils::trim_comment(line),"\"' \t\n\r\f").as_vector();
    } catch (TokenizerException &e) {
      // ignore
    }
    if (words.size() != 3)
      error->all(FLERR,"Incorrect format in SNAP coefficient file");

    int jelem;
    for (jelem = 0; jelem < nelements; jelem++)
      if (words[0] == elements[jelem]) break;

    // if this element not needed, skip this block

    if (jelem == nelements) {
      if (comm->me == 0) {
        for (int icoeff = 0; icoeff < ncoeffall; icoeff++) {
          ptr = fgets(line,MAXLINE,fpcoeff);
          if (ptr == nullptr) {
            eof = 1;
            fclose(fpcoeff);
          }
        }
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof)
        error->all(FLERR,"Incorrect format in SNAP coefficient file");
      continue;
    }

    if (elementflags[jelem] == 1)
      error->all(FLERR,"Incorrect format in SNAP coefficient file");
    else
      elementflags[jelem] = 1;

    radelem[jelem] = utils::numeric(FLERR,words[1].c_str(),false,lmp);
    wjelem[jelem] = utils::numeric(FLERR,words[2].c_str(),false,lmp);

    if (comm->me == 0)
      utils::logmesg(lmp,"SNAP Element = {}, Radius {}, Weight {}\n",
                     elements[jelem], radelem[jelem], wjelem[jelem]);

    for (int icoeff = 0; icoeff < ncoeffall; icoeff++) {
      if (comm->me == 0) {
        ptr = fgets(line,MAXLINE,fpcoeff);
        if (ptr == nullptr) {
          eof = 1;
          fclose(fpcoeff);
        }
      }

      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof)
        error->all(FLERR,"Incorrect format in SNAP coefficient file");
      MPI_Bcast(line,MAXLINE,MPI_CHAR,0,world);

      try {
        ValueTokenizer coeff(utils::trim_comment(line));
        if (coeff.count() != 1)
          error->all(FLERR,"Incorrect format in SNAP coefficient file");

        coeffelem[jelem][icoeff] = coeff.next_double();
      } catch (TokenizerException &e) {
        error->all(FLERR,"Incorrect format in SNAP coefficient "
                                     "file: {}", e.what());
      }
    }
  }

  if (comm->me == 0) fclose(fpcoeff);

  for (int jelem = 0; jelem < nelements; jelem++) {
    if (elementflags[jelem] == 0)
      error->all(FLERR,"Element {} not found in SNAP coefficient "
                                   "file", elements[jelem]);
  }
  delete[] elementflags;

  // set flags for required keywords

  rcutfacflag = 0;
  twojmaxflag = 0;

  // Set defaults for optional keywords

  rfac0 = 0.99363;
  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;
  chemflag = 0;
  bnormflag = 0;
  wselfallflag = 0;
  chunksize = 4096;

  // open SNAP parameter file on proc 0

  FILE *fpparam;
  if (comm->me == 0) {
    fpparam = utils::open_potential(paramfilename,lmp,nullptr);
    if (fpparam == nullptr)
      error->one(FLERR,"Cannot open SNAP parameter file {}: {}",
                                   paramfilename, utils::getsyserror());
  }

  eof = 0;
  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpparam);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fpparam);
      }
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(line,MAXLINE,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    std::vector<std::string> words;
    try {
      words = Tokenizer(utils::trim_comment(line),"\"' \t\n\r\f").as_vector();
    } catch (TokenizerException &e) {
      // ignore
    }

    if (words.size() == 0) continue;
    if (words.size() != 2)
      error->all(FLERR,"Incorrect format in SNAP parameter file");

    auto keywd = words[0];
    auto keyval = words[1];

    if (comm->me == 0)
      utils::logmesg(lmp,"SNAP keyword {} {}\n",keywd,keyval);

    if (keywd == "rcutfac") {
      rcutfac = utils::numeric(FLERR,keyval.c_str(),false,lmp);
      rcutfacflag = 1;
    } else if (keywd == "twojmax") {
      twojmax = utils::inumeric(FLERR,keyval.c_str(),false,lmp);
      twojmaxflag = 1;
    } else if (keywd == "rfac0")
      rfac0 = utils::numeric(FLERR,keyval.c_str(),false,lmp);
    else if (keywd == "rmin0")
      rmin0 = utils::numeric(FLERR,keyval.c_str(),false,lmp);
    else if (keywd == "switchflag")
      switchflag = utils::inumeric(FLERR,keyval.c_str(),false,lmp);
    else if (keywd == "bzeroflag")
      bzeroflag = utils::inumeric(FLERR,keyval.c_str(),false,lmp);
    else if (keywd == "quadraticflag")
      quadraticflag = utils::inumeric(FLERR,keyval.c_str(),false,lmp);
    else if (keywd == "chemflag")
      chemflag = utils::inumeric(FLERR,keyval.c_str(),false,lmp);
    else if (keywd == "bnormflag")
      bnormflag = utils::inumeric(FLERR,keyval.c_str(),false,lmp);
    else if (keywd == "wselfallflag")
      wselfallflag = utils::inumeric(FLERR,keyval.c_str(),false,lmp);
    else if (keywd == "chunksize")
      chunksize = utils::inumeric(FLERR,keyval.c_str(),false,lmp);
    else
      error->all(FLERR,"Unknown parameter '{}' in SNAP "
                                   "parameter file", keywd);
  }

  if (rcutfacflag == 0 || twojmaxflag == 0)
    error->all(FLERR,"Incorrect SNAP parameter file");

  if (chemflag && nelemtmp != nelements)
    error->all(FLERR,"Incorrect SNAP parameter file");

}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double PairSNAP::memory_usage()
{
  double bytes = Pair::memory_usage();

  int n = atom->ntypes+1;
  bytes += (double)n*n*sizeof(int);      // setflag
  bytes += (double)n*n*sizeof(double);   // cutsq
  bytes += (double)n*n*sizeof(double);   // scale
  bytes += (double)n*sizeof(int);        // map
  bytes += (double)beta_max*ncoeff*sizeof(double); // bispectrum
  bytes += (double)beta_max*ncoeff*sizeof(double); // beta

  bytes += snaptr->memory_usage(); // SNA object

  return bytes;
}

void *PairSNAP::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  return nullptr;
}
