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

#if defined(__AVX512F__)
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)

#include "pair_snap_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "sna_intel.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace ip_simd;

static constexpr int MAXLINE = 1024;
static constexpr int MAXWORD = 3;

/* ---------------------------------------------------------------------- */

PairSNAPIntel::PairSNAPIntel(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;

  radelem = nullptr;
  wjelem = nullptr;
  coeffelem = nullptr;
  sinnerelem = nullptr;
  dinnerelem = nullptr;

  beta = nullptr;
  bispectrum = nullptr;
  snaptr = nullptr;
}

/* ---------------------------------------------------------------------- */

PairSNAPIntel::~PairSNAPIntel()
{
  if (copymode) return;

  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(coeffelem);
  memory->destroy(sinnerelem);
  memory->destroy(dinnerelem);

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

void PairSNAPIntel::compute(int eflag, int vflag)
{
  SNA_DVEC fij[3];
  int *jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);
  int tally_xyz = 0;
  if (vflag_atom || (vflag && !vflag_fdotr)) tally_xyz = 1;

  double **x = atom->x;
  double *_x = atom->x[0];
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // compute dE_i/dB_i = beta_i for all i in list

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  SNA_DVEC sevdwl(0);

  const int vw = snaptr->vector_width();
  for (int ii = 0; ii < list->inum; ii+=vw) {
    SNA_IVEC i, jnum;
    int max_jnum = 0;
    for (int l = 0; l < vw; l++) {
      if (ii + l < list->inum) {
        i[l] = list->ilist[ii + l];
        jnum[l] = numneigh[i[l]];
      } else {
        i[l] = list->ilist[0];
        jnum[l] = 0;
      }
      if (jnum[l] > max_jnum) max_jnum = jnum[l];
    }

    // ensure rij, inside, wj, and rcutij are of size jnum

    snaptr->grow_rij(max_jnum);

    SNA_IVEC zero_vec(0);

    const SNA_DVEC xtmp = SIMD_gather(_x, i * 3);
    const SNA_DVEC ytmp = SIMD_gather(_x, i * 3 + 1);
    const SNA_DVEC ztmp = SIMD_gather(_x, i * 3 + 2);
    const SNA_IVEC itype = SIMD_gather(type, i);
    const SNA_IVEC ielem = SIMD_gather(map, itype);
    const SNA_DVEC radi = SIMD_gather(radelem, ielem);

    // rij[][3] = displacements between atom I and those neighbors
    // inside = indices of neighbors of I within cutoff
    // wj = weights for neighbors of I within cutoff
    // rcutij = cutoffs for neighbors of I within cutoff
    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    SNA_IVEC ninside(0);
    for (int jj = 0; jj < max_jnum; jj++) {
      SIMD_mask m(SIMD256_set(jj) < jnum);

      SNA_IVEC j;
      SV_for (int l = 0; l < vw; l++) {
        jlist = firstneigh[i[l]];
        if (jj < jnum[l]) j[l] = jlist[jj];
        else j[l] = 0;
      }
      j &= NEIGHMASK;

      const SNA_DVEC delx = SIMD_gather(m, _x, j * 3) - xtmp;
      const SNA_DVEC dely = SIMD_gather(m, _x, j * 3 + 1) - ytmp;
      const SNA_DVEC delz = SIMD_gather(m, _x, j * 3 + 2) - ztmp;
      const SNA_IVEC jtype = SIMD_gather(type, j);
      const SNA_DVEC rsq = delx*delx + dely*dely + delz*delz;
      const SNA_DVEC vcut = SIMD_gather(m, cutsq[0],
                                        itype * (atom->ntypes + 1) + jtype);

      m &= rsq < vcut;
      m &= rsq > SIMD_set(1e-20);
      const SNA_IVEC jelem = SIMD_gather(map, jtype);
      const SNA_IVEC ni3 = ninside * vw * 3 + SIMD256_count();
      SIMD_scatter(m, (double *)(snaptr->rij[0]), ni3, delx);
      SIMD_scatter(m, (double *)(snaptr->rij[0] + 1), ni3, dely);
      SIMD_scatter(m, (double *)(snaptr->rij[0] + 2), ni3, delz);
      const SNA_IVEC ni = ninside * vw + SIMD256_count();
      SIMD_scatter(m, (int *)(snaptr->inside), ni, j);
      SIMD_scatter(m, (double *)(snaptr->wj), ni,
                   SIMD_gather(m, wjelem, jelem));
      SIMD_scatter(m, (double *)(snaptr->rcutij), ni,
                   (radi + SIMD_gather(m, radelem, jelem)) * rcutfac);
      if (switchinnerflag) {
        SIMD_scatter(m, (double *)(snaptr->sinnerij), ni,
                     (SIMD_gather(m, sinnerelem, ielem) +
                      SIMD_gather(m, sinnerelem, jelem)) * 0.5);
        SIMD_scatter(m, (double *)(snaptr->dinnerij), ni,
                     (SIMD_gather(m, dinnerelem, ielem) +
                      SIMD_gather(m, dinnerelem, jelem)) * 0.5);
      }
      if (chemflag)
        SIMD_scatter(m, (int *)(snaptr->element), ni, jelem);
      ninside = SIMD_add(m, ninside, 1);
    } // for jj

    // compute Ui, Yi for atom I

    if (chemflag)
      snaptr->compute_ui(ninside, ielem, max_jnum);
    else
      snaptr->compute_ui(ninside, zero_vec, max_jnum);

    // Compute bispectrum
    if (quadraticflag || eflag) {
      snaptr->compute_zi_or_yi<0>(beta);
      if (chemflag)
        snaptr->compute_bi(ielem);
      else
        snaptr->compute_bi(zero_vec);
      for (int icoeff = 0; icoeff < ncoeff; icoeff++)
        SIMD_store(bispectrum + icoeff, SIMD_load(snaptr->blist + icoeff));
    }

    // Compute beta
    for (int icoeff = 0; icoeff < ncoeff; icoeff++)
      SIMD_store(beta + icoeff, SIMD_gather(coeffelem[0],
                                            ielem * ncoeffall + icoeff + 1));

    if (quadraticflag) {
      int k = ncoeff+1;
      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        SNA_DVEC bveci = SIMD_load(bispectrum + icoeff);
        SNA_DVEC beta_i = SIMD_load(beta + icoeff) +
          SIMD_gather(coeffelem[0], ielem * ncoeffall + k) * bveci;
        k++;
        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
          const SNA_DVEC ci = SIMD_gather(coeffelem[0], ielem * ncoeffall + k);
          beta_i = beta_i + ci * SIMD_load(bispectrum + jcoeff);
          SIMD_store(beta + jcoeff, ci * bveci + SIMD_load(beta + jcoeff));
          k++;
        }
        SIMD_store(beta + icoeff, beta_i);
      }
    }

    // for neighbors of I within cutoff:
    // compute Fij = dEi/dRj = -dEi/dRi
    // add to Fi, subtract from Fj
    // scaling is that for type I

    if (quadraticflag || eflag)
      snaptr->compute_yi_from_zi(beta);
    else
      snaptr->compute_zi_or_yi<1>(beta);

    SNA_DVEC fi_x(0.0), fi_y(0.0), fi_z(0.0);
    SNA_DVEC scalev = SIMD_gather(scale[0], itype * (atom->ntypes+1) + itype);
    for (int jj = 0; jj < max_jnum; jj++) {
      snaptr->compute_duidrj(jj, ninside);
      if (chemflag && nelements > 1)
        snaptr->compute_deidrj_e(jj, ninside, fij);
      else
        snaptr->compute_deidrj(jj, ninside, fij);

      SNA_DVEC fijs_x = fij[0] * scalev;
      SNA_DVEC fijs_y = fij[1] * scalev;
      SNA_DVEC fijs_z = fij[2] * scalev;

      fi_x += fijs_x;
      fi_y += fijs_y;
      fi_z += fijs_z;

      for (int l = 0; l < vw; l++) {
        if (jj < ninside[l]) {
          int j = snaptr->inside[jj][l];
          f[j][0] -= fijs_x[l];
          f[j][1] -= fijs_y[l];
          f[j][2] -= fijs_z[l];

          if (tally_xyz)
            ev_tally_xyz(i[l],j,nlocal,newton_pair,0.0,0.0,
                         fij[0][l],fij[1][l],fij[2][l],
                         -snaptr->rij[jj][0][l],-snaptr->rij[jj][1][l],
                         -snaptr->rij[jj][2][l]);
        }
      } // for l
    } // for jj
    SIMD_mask m((SIMD256_count() + ii) < list->inum);
    SNA_DVEC fix = SIMD_gather(m, f[0], i * 3) +  fi_x;
    SIMD_scatter(m, f[0], i * 3, fix);
    SNA_DVEC fiy = SIMD_gather(m, f[0], i * 3 + 1) +  fi_y;
    SIMD_scatter(m, f[0], i * 3 + 1, fiy);
    SNA_DVEC fiz = SIMD_gather(m, f[0], i * 3 + 2) +  fi_z;
    SIMD_scatter(m, f[0], i * 3 + 2, fiz);

    // tally energy contribution

    if (eflag) {
      SNA_DVEC evdwl = SIMD_gather(coeffelem[0], ielem * ncoeffall);
      for (int icoeff = 0; icoeff < ncoeff; icoeff++)
        evdwl += SIMD_gather(coeffelem[0], ielem * ncoeffall + icoeff +1) *
          bispectrum[icoeff];

      if (quadraticflag) {
        int k = ncoeff+1;
        for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
          SNA_DVEC bveci = SIMD_load(bispectrum + icoeff);
          SNA_DVEC c = SIMD_gather(coeffelem[0], ielem * ncoeffall + k);
          k++;
          evdwl += c * 0.5 * bveci * bveci;
          for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
            SNA_DVEC bvecj = SIMD_load(bispectrum + jcoeff);
            SNA_DVEC cj = SIMD_gather(coeffelem[0], ielem * ncoeffall + k);
            k++;
            evdwl += cj * bveci * bvecj;
          }
        }
      }
      sevdwl += scalev * evdwl;
      if (eatom) {
        SNA_DVEC ea = SIMD_gather(m, eatom, i) + scalev * evdwl;
        SIMD_scatter(m, eatom, i, ea);
      }
    } // if (eflag)
  } // for ii
  if (eflag) eng_vdwl += SIMD_sum(sevdwl);
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSNAPIntel::allocate()
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

void PairSNAPIntel::settings(int narg, char ** /* arg */)
{
  if (narg > 0)
    error->all(FLERR,"Illegal pair_style command");
  if ((comm->me == 0) && (comm->nthreads > 1))
    error->warning(FLERR, "Pair style snap/intel does not use OpenMP threads");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSNAPIntel::coeff(int narg, char **arg)
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

    ncoeff = sqrt(2.0*ncoeffall)-1;
    ncoeffq = (ncoeff*(ncoeff+1))/2;
    int ntmp = 1+ncoeff+ncoeffq;
    if (ntmp != ncoeffall) {
      error->all(FLERR,"Incorrect SNAP coeff file");
    }
  }

  snaptr = new SNAIntel(lmp, rfac0, twojmax,
                        rmin0, switchflag, bzeroflag,
                        chemflag, bnormflag, wselfallflag,
                        nelements, switchinnerflag);

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

void PairSNAPIntel::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style SNAP requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);

  snaptr->init();

  fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");

  fix->pair_init_check();

  memory->create(bispectrum,ncoeff,"PairSNAP:bispectrum");
  memory->create(beta,ncoeff,"PairSNAP:beta");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSNAPIntel::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  scale[j][i] = scale[i][j];
  return (radelem[map[i]] +
          radelem[map[j]])*rcutfac;
}

/* ---------------------------------------------------------------------- */

void PairSNAPIntel::read_files(char *coefffilename, char *paramfilename)
{

  // open SNAP coefficient file on proc 0

  FILE *fpcoeff;
  if (comm->me == 0) {
    fpcoeff = utils::open_potential(coefffilename,lmp,nullptr);
    if (fpcoeff == nullptr)
      error->one(FLERR,"Cannot open SNAP coefficient file {}: ",
                                   coefffilename, utils::getsyserror());
  }

  char line[MAXLINE] = {'\0'};
  char *ptr;
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
    error->all(FLERR,"Incorrect format in SNAP coefficient file: {}", e.what());
  }

  // clean out old arrays and set up element lists

  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(coeffelem);
  memory->destroy(sinnerelem);
  memory->destroy(dinnerelem);
  memory->create(radelem,nelements,"pair:radelem");
  memory->create(wjelem,nelements,"pair:wjelem");
  memory->create(coeffelem,nelements,ncoeffall,"pair:coeffelem");
  memory->create(sinnerelem,nelements,"pair:sinnerelem");
  memory->create(dinnerelem,nelements,"pair:dinnerelem");

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
    } catch (TokenizerException &) {
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

    radelem[jelem] = utils::numeric(FLERR,words[1],false,lmp);
    wjelem[jelem] = utils::numeric(FLERR,words[2],false,lmp);

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
        error->all(FLERR,"Incorrect format in SNAP coefficient file: {}", e.what());
      }
    }
  }

  if (comm->me == 0) fclose(fpcoeff);

  for (int jelem = 0; jelem < nelements; jelem++) {
    if (elementflags[jelem] == 0)
      error->all(FLERR,"Element {} not found in SNAP coefficient file", elements[jelem]);
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
  switchinnerflag = 0;
  chunksize = 32768;
  parallel_thresh = 8192;

  // set local input checks

  int sinnerflag = 0;
  int dinnerflag = 0;

  // open SNAP parameter file on proc 0

  FILE *fpparam;
  if (comm->me == 0) {
    fpparam = utils::open_potential(paramfilename,lmp,nullptr);
    if (fpparam == nullptr)
      error->one(FLERR,"Cannot open SNAP parameter file {}: {}",
                                   paramfilename, utils::getsyserror());
  }

  eof = 0;
  while (true) {
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

    // words = ptrs to all words in line
    // strip single and double quotes from words

    std::vector<std::string> words;
    try {
      words = Tokenizer(utils::trim_comment(line),"\"' \t\n\r\f").as_vector();
    } catch (TokenizerException &) {
      // ignore
    }

    if (words.size() == 0) continue;

    if (words.size() < 2)
      error->all(FLERR,"Incorrect format in SNAP parameter file");

    auto keywd = words[0];
    auto keyval = words[1];

    // check for keywords with more than one value per element

    if (keywd == "sinner" || keywd == "dinner") {

      if ((int)words.size() != nelements+1)
        error->all(FLERR,"Incorrect SNAP parameter file");

      // innerlogstr collects all values of sinner or dinner for log output below

      std::string innerlogstr;

      int iword = 1;

      if (keywd == "sinner") {
        for (int ielem = 0; ielem < nelements; ielem++) {
          keyval = words[iword];
          sinnerelem[ielem] = utils::numeric(FLERR,keyval,false,lmp);
          iword++;
          innerlogstr += keyval + " ";
        }
        sinnerflag = 1;
      } else if (keywd == "dinner") {
        for (int ielem = 0; ielem < nelements; ielem++) {
          keyval = words[iword];
          dinnerelem[ielem] = utils::numeric(FLERR,keyval,false,lmp);
          iword++;
          innerlogstr += keyval + " ";
        }
        dinnerflag = 1;
      }

      if (comm->me == 0)
        utils::logmesg(lmp,"SNAP keyword {} {} ... \n", keywd, innerlogstr);

    } else {

      // all other keywords take one value

      if (nwords != 2)
        error->all(FLERR,"Incorrect SNAP parameter file");

      if (comm->me == 0)
        utils::logmesg(lmp,"SNAP keyword {} {}\n",keywd,keyval);

      if (keywd == "rcutfac") {
        rcutfac = utils::numeric(FLERR,keyval,false,lmp);
        rcutfacflag = 1;
      } else if (keywd == "twojmax") {
        twojmax = utils::inumeric(FLERR,keyval,false,lmp);
        twojmaxflag = 1;
      } else if (keywd == "rfac0")
        rfac0 = utils::numeric(FLERR,keyval,false,lmp);
      else if (keywd == "rmin0")
        rmin0 = utils::numeric(FLERR,keyval,false,lmp);
      else if (keywd == "switchflag")
        switchflag = utils::inumeric(FLERR,keyval,false,lmp);
      else if (keywd == "bzeroflag")
        bzeroflag = utils::inumeric(FLERR,keyval,false,lmp);
      else if (keywd == "quadraticflag")
        quadraticflag = utils::inumeric(FLERR,keyval,false,lmp);
      else if (keywd == "chemflag")
        chemflag = utils::inumeric(FLERR,keyval,false,lmp);
      else if (keywd == "bnormflag")
        bnormflag = utils::inumeric(FLERR,keyval,false,lmp);
      else if (keywd == "wselfallflag")
        wselfallflag = utils::inumeric(FLERR,keyval,false,lmp);
      else if (keywd == "switchinnerflag")
        switchinnerflag = utils::inumeric(FLERR,keyval,false,lmp);
      else if (keywd == "chunksize")
        chunksize = utils::inumeric(FLERR,keyval,false,lmp);
      else if (keywd == "parallelthresh")
        parallel_thresh = utils::inumeric(FLERR,keyval,false,lmp);
      else
        error->all(FLERR,"Unknown parameter '{}' in SNAP parameter file", keywd);
    }
  }

  if (rcutfacflag == 0 || twojmaxflag == 0)
    error->all(FLERR,"Incorrect SNAP parameter file");

  if (chemflag && nelemtmp != nelements)
    error->all(FLERR,"Incorrect SNAP parameter file");

  if (switchinnerflag && !(sinnerflag && dinnerflag))
    error->all(FLERR,"Incorrect SNAP parameter file");

  if (!switchinnerflag && (sinnerflag || dinnerflag))
    error->all(FLERR,"Incorrect SNAP parameter file");
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double PairSNAPIntel::memory_usage()
{
  double bytes = Pair::memory_usage();

  int n = atom->ntypes+1;
  bytes += (double)n*n*sizeof(int);         // setflag
  bytes += (double)n*n*sizeof(double);      // cutsq
  bytes += (double)n*n*sizeof(double);      // scale
  bytes += (double)n*sizeof(int);           // map
  bytes += (double)ncoeff*sizeof(SNA_DVEC); // bispectrum
  bytes += (double)ncoeff*sizeof(SNA_DVEC); // beta

  bytes += snaptr->memory_usage(); // SNA object

  return bytes;
}

/* ---------------------------------------------------------------------- */

void *PairSNAPIntel::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  return nullptr;
}

#endif
#endif
