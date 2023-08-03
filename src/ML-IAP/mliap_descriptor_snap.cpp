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
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "mliap_descriptor_snap.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "mliap_data.h"
#include "pair_mliap.h"
#include "sna.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

MLIAPDescriptorSNAP::MLIAPDescriptorSNAP(LAMMPS *_lmp, char *paramfilename) : Pointers(_lmp),  MLIAPDescriptor(_lmp)
{
  radelem = nullptr;
  wjelem = nullptr;
  snaptr = nullptr;
  sinnerelem = nullptr;
  dinnerelem = nullptr;

  read_paramfile(paramfilename);

  snaptr = new SNA(lmp, rfac0, twojmax, rmin0, switchflag, bzeroflag, chemflag, bnormflag,
                   wselfallflag, nelements, switchinnerflag);

  ndescriptors = snaptr->ncoeff;
}

/* ---------------------------------------------------------------------- */

MLIAPDescriptorSNAP::~MLIAPDescriptorSNAP()
{
  memory->destroy(radelem);
  memory->destroy(wjelem);
  delete snaptr;
  memory->destroy(sinnerelem);
  memory->destroy(dinnerelem);
}

/* ----------------------------------------------------------------------
   compute descriptors for each atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_descriptors(class MLIAPData *data)
{
  int ij = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int ielem = data->ielems[ii];

    // ensure rij, inside, wj, and rcutij are of size jnum

    const int jnum = data->numneighs[ii];
    snaptr->grow_rij(jnum);

    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      const int j = data->jatoms[ij];
      const int jelem = data->jelems[ij];
      const double *delr = data->rij[ij];

      snaptr->rij[ninside][0] = delr[0];
      snaptr->rij[ninside][1] = delr[1];
      snaptr->rij[ninside][2] = delr[2];
      snaptr->inside[ninside] = j;
      snaptr->wj[ninside] = wjelem[jelem];
      snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
      if (switchinnerflag) {
        snaptr->sinnerij[ninside] = 0.5 * (sinnerelem[ielem] + sinnerelem[jelem]);
        snaptr->dinnerij[ninside] = 0.5 * (dinnerelem[ielem] + dinnerelem[jelem]);
      }
      if (chemflag) snaptr->element[ninside] = jelem;
      ninside++;
      ij++;
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

    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
      data->descriptors[ii][icoeff] = snaptr->blist[icoeff];
  }
}

/* ----------------------------------------------------------------------
   compute forces for each atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_forces(class MLIAPData *data)
{
  double fij[3];
  double **f = atom->f;

  int ij = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int i = data->iatoms[ii];
    const int ielem = data->ielems[ii];

    // ensure rij, inside, wj, and rcutij are of size jnum

    const int jnum = data->numneighs[ii];
    snaptr->grow_rij(jnum);

    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      const int j = data->jatoms[ij];
      const int jelem = data->jelems[ij];
      const double *delr = data->rij[ij];

      snaptr->rij[ninside][0] = delr[0];
      snaptr->rij[ninside][1] = delr[1];
      snaptr->rij[ninside][2] = delr[2];
      snaptr->inside[ninside] = j;
      snaptr->wj[ninside] = wjelem[jelem];
      snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
      if (switchinnerflag) {
        snaptr->sinnerij[ninside] = 0.5 * (sinnerelem[ielem] + sinnerelem[jelem]);
        snaptr->dinnerij[ninside] = 0.5 * (dinnerelem[ielem] + dinnerelem[jelem]);
      }
      if (chemflag) snaptr->element[ninside] = jelem;
      ninside++;
      ij++;
    }

    // compute Ui, Yi for atom I

    if (chemflag)
      snaptr->compute_ui(ninside, ielem);
    else
      snaptr->compute_ui(ninside, 0);

    // for neighbors of I within cutoff:
    // compute Fij = dEi/dRj = -dEi/dRi
    // add to Fi, subtract from Fj

    snaptr->compute_yi(data->betas[ii]);

    for (int jj = 0; jj < ninside; jj++) {
      int j = snaptr->inside[jj];
      snaptr->compute_duidrj(jj);

      snaptr->compute_deidrj(fij);

      f[i][0] += fij[0];
      f[i][1] += fij[1];
      f[i][2] += fij[2];
      f[j][0] -= fij[0];
      f[j][1] -= fij[1];
      f[j][2] -= fij[2];

      // add in global and per-atom virial contributions
      // this is optional and has no effect on force calculation

      if (data->vflag) data->pairmliap->v_tally(i, j, fij, snaptr->rij[jj]);
    }
  }
}

/* ----------------------------------------------------------------------
   calculate gradients of forces w.r.t. parameters
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_force_gradients(class MLIAPData *data)
{
  int ij = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int i = data->iatoms[ii];
    const int ielem = data->ielems[ii];

    // ensure rij, inside, wj, and rcutij are of size jnum

    const int jnum = data->numneighs[ii];
    snaptr->grow_rij(jnum);

    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      const int j = data->jatoms[ij];
      const int jelem = data->jelems[ij];

      const double *delr = data->rij[ij];

      snaptr->rij[ninside][0] = delr[0];
      snaptr->rij[ninside][1] = delr[1];
      snaptr->rij[ninside][2] = delr[2];
      snaptr->inside[ninside] = j;
      snaptr->wj[ninside] = wjelem[jelem];
      snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
      if (switchinnerflag) {
        snaptr->sinnerij[ninside] = 0.5 * (sinnerelem[ielem] + sinnerelem[jelem]);
        snaptr->dinnerij[ninside] = 0.5 * (dinnerelem[ielem] + dinnerelem[jelem]);
      }
      if (chemflag) snaptr->element[ninside] = jelem;
      ninside++;
      ij++;
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

    for (int jj = 0; jj < ninside; jj++) {
      const int j = snaptr->inside[jj];

      snaptr->compute_duidrj(jj);
      snaptr->compute_dbidrj();

      // Accumulate gamma_lk*dB_k/dRi, -gamma_lk**dB_k/dRj

      for (int inz = 0; inz < data->gamma_nnz; inz++) {
        const int l = data->gamma_row_index[ii][inz];
        const int k = data->gamma_col_index[ii][inz];
        data->gradforce[i][l] += data->gamma[ii][inz] * snaptr->dblist[k][0];
        data->gradforce[i][l + data->yoffset] += data->gamma[ii][inz] * snaptr->dblist[k][1];
        data->gradforce[i][l + data->zoffset] += data->gamma[ii][inz] * snaptr->dblist[k][2];
        data->gradforce[j][l] -= data->gamma[ii][inz] * snaptr->dblist[k][0];
        data->gradforce[j][l + data->yoffset] -= data->gamma[ii][inz] * snaptr->dblist[k][1];
        data->gradforce[j][l + data->zoffset] -= data->gamma[ii][inz] * snaptr->dblist[k][2];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute descriptor gradients for each neighbor atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_descriptor_gradients(class MLIAPData *data)
{
  int ij = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int ielem = data->ielems[ii];

    // ensure rij, inside, wj, and rcutij are of size jnum

    const int jnum = data->numneighs[ii];
    snaptr->grow_rij(jnum);

    int ij0 = ij;
    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      const int j = data->jatoms[ij];
      const int jelem = data->jelems[ij];
      const double *delr = data->rij[ij];

      snaptr->rij[ninside][0] = delr[0];
      snaptr->rij[ninside][1] = delr[1];
      snaptr->rij[ninside][2] = delr[2];
      snaptr->inside[ninside] = j;
      snaptr->wj[ninside] = wjelem[jelem];
      snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
      if (switchinnerflag) {
        snaptr->sinnerij[ninside] = 0.5 * (sinnerelem[ielem] + sinnerelem[jelem]);
        snaptr->dinnerij[ninside] = 0.5 * (dinnerelem[ielem] + dinnerelem[jelem]);
      }
      if (chemflag) snaptr->element[ninside] = jelem;
      ninside++;
      ij++;
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

    ij = ij0;
    for (int jj = 0; jj < ninside; jj++) {
      snaptr->compute_duidrj(jj);
      snaptr->compute_dbidrj();

      // Accumulate dB_k^i/dRi, dB_k^i/dRj

      for (int k = 0; k < data->ndescriptors; k++) {
        data->graddesc[ij][k][0] = snaptr->dblist[k][0];
        data->graddesc[ij][k][1] = snaptr->dblist[k][1];
        data->graddesc[ij][k][2] = snaptr->dblist[k][2];
      }
      ij++;
    }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::init()
{
  snaptr->init();
}

/* ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::read_paramfile(char *paramfilename)
{

  // set flags for required keywords

  int rcutfacflag = 0;
  int twojmaxflag = 0;
  int nelementsflag = 0;
  int elementsflag = 0;
  int radelemflag = 0;
  int wjelemflag = 0;

  // Set defaults for optional keywords

  rfac0 = 0.99363;
  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  chemflag = 0;
  bnormflag = 0;
  wselfallflag = 0;
  switchinnerflag = 0;

  // set local input checks

  int sinnerflag = 0;
  int dinnerflag = 0;

  for (int i = 0; i < nelements; i++) delete[] elements[i];
  delete[] elements;
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);

  // open SNAP parameter file on proc 0

  FILE *fpparam;
  if (comm->me == 0) {
    fpparam = utils::open_potential(paramfilename, lmp, nullptr);
    if (fpparam == nullptr)
      error->one(FLERR, "Cannot open SNAP parameter file {}: {}", paramfilename,
                 utils::getsyserror());
  }

  char line[MAXLINE], *ptr;
  int eof = 0;
  int n;

  while (true) {
    if (comm->me == 0) {
      ptr = fgets(line, MAXLINE, fpparam);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fpparam);
      } else
        n = strlen(line) + 1;
    }
    MPI_Bcast(&eof, 1, MPI_INT, 0, world);
    if (eof) break;
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    MPI_Bcast(line, n, MPI_CHAR, 0, world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line, '#'))) *ptr = '\0';

    // strip single and double quotes from words
    Tokenizer words(line, "\"' \t\n\t\f");
    if (words.count() == 0) continue;

    auto keywd = words.next();

    // check for keywords with one value per element

    if ((keywd == "elems") || (keywd == "radelems") || (keywd == "welems") ||
        (keywd == "sinnerelems") || (keywd == "dinnerelems")) {

      if ((nelementsflag == 0) || ((int) words.count() != nelements + 1))
        error->all(FLERR, "Incorrect SNAP parameter file");

      if (comm->me == 0) utils::logmesg(lmp, "SNAP keyword {} \n", utils::trim(line));

      if (keywd == "elems") {
        for (int ielem = 0; ielem < nelements; ielem++)
          elements[ielem] = utils::strdup(words.next());
        elementsflag = 1;
      } else if (keywd == "radelems") {
        for (int ielem = 0; ielem < nelements; ielem++)
          radelem[ielem] = utils::numeric(FLERR, words.next(), false, lmp);
        radelemflag = 1;
      } else if (keywd == "welems") {
        for (int ielem = 0; ielem < nelements; ielem++)
          wjelem[ielem] = utils::numeric(FLERR, words.next(), false, lmp);
        wjelemflag = 1;
      } else if (keywd == "sinnerelems") {
        for (int ielem = 0; ielem < nelements; ielem++)
          sinnerelem[ielem] = utils::numeric(FLERR, words.next(), false, lmp);
        sinnerflag = 1;
      } else if (keywd == "dinnerelems") {
        for (int ielem = 0; ielem < nelements; ielem++)
          dinnerelem[ielem] = utils::numeric(FLERR, words.next(), false, lmp);
        dinnerflag = 1;
      }

    } else {

      // all other keywords take one value

      if (words.count() != 2) error->all(FLERR, "Incorrect SNAP parameter file");

      auto keyval = words.next();
      if (comm->me == 0) utils::logmesg(lmp, "SNAP keyword {} {} \n", keywd, keyval);

      if (keywd == "nelems") {
        nelements = utils::inumeric(FLERR, keyval, false, lmp);
        elements = new char *[nelements];
        memory->create(radelem, nelements, "mliap_snap_descriptor:radelem");
        memory->create(wjelem, nelements, "mliap_snap_descriptor:wjelem");
        memory->create(sinnerelem, nelements, "mliap_snap_descriptor:sinner");
        memory->create(dinnerelem, nelements, "mliap_snap_descriptor:dinner");
        nelementsflag = 1;
      } else if (keywd == "rcutfac") {
        rcutfac = utils::numeric(FLERR, keyval, false, lmp);
        rcutfacflag = 1;
      } else if (keywd == "twojmax") {
        twojmax = utils::inumeric(FLERR, keyval, false, lmp);
        twojmaxflag = 1;
      } else if (keywd == "rfac0")
        rfac0 = utils::numeric(FLERR, keyval, false, lmp);
      else if (keywd == "rmin0")
        rmin0 = utils::numeric(FLERR, keyval, false, lmp);
      else if (keywd == "switchflag")
        switchflag = utils::inumeric(FLERR, keyval, false, lmp);
      else if (keywd == "bzeroflag")
        bzeroflag = utils::inumeric(FLERR, keyval, false, lmp);
      else if (keywd == "chemflag")
        chemflag = utils::inumeric(FLERR, keyval, false, lmp);
      else if (keywd == "bnormflag")
        bnormflag = utils::inumeric(FLERR, keyval, false, lmp);
      else if (keywd == "wselfallflag")
        wselfallflag = utils::inumeric(FLERR, keyval, false, lmp);
      else if (keywd == "switchinnerflag")
        switchinnerflag = utils::inumeric(FLERR, keyval, false, lmp);
      else
        error->all(FLERR, "Incorrect SNAP parameter file");
    }
  }

  if (!rcutfacflag || !twojmaxflag || !nelementsflag || !elementsflag || !radelemflag ||
      !wjelemflag)
    error->all(FLERR, "Incorrect SNAP parameter file");

  if (switchinnerflag && !(sinnerflag && dinnerflag))
    error->all(FLERR, "Incorrect SNAP parameter file");

  if (!switchinnerflag && (sinnerflag || dinnerflag))
    error->all(FLERR, "Incorrect SNAP parameter file");

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq, nelements, nelements, "mliap/descriptor/snap:cutsq");
  for (int ielem = 0; ielem < nelements; ielem++) {
    cut = 2.0 * radelem[ielem] * rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[ielem][ielem] = cut * cut;
    for (int jelem = ielem + 1; jelem < nelements; jelem++) {
      cut = (radelem[ielem] + radelem[jelem]) * rcutfac;
      cutsq[ielem][jelem] = cutsq[jelem][ielem] = cut * cut;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double MLIAPDescriptorSNAP::memory_usage()
{
  double bytes = MLIAPDescriptor::memory_usage();
  bytes += snaptr->memory_usage();    // SNA object

  return bytes;
}
