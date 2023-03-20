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
 Contributing authors: Byungkyun Kang (University of Nevada, Las Vegas)
 ------------------------------------------------------------------------- */

#include "mliap_descriptor_so3.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "mliap_data.h"
#include "mliap_so3.h"
#include "pair_mliap.h"
#include "tokenizer.h"

#include <cstring>

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

MLIAPDescriptorSO3::MLIAPDescriptorSO3(LAMMPS *lmp, char *paramfilename) : Pointers(lmp), MLIAPDescriptor(lmp)
{
  radelem = nullptr;
  wjelem = nullptr;
  so3ptr = nullptr;
  read_paramfile(paramfilename);

  so3ptr = new MLIAP_SO3(lmp, rcutfac, lmax, nmax, alpha);

  ndescriptors = so3ptr->ncoeff;
}

/* ---------------------------------------------------------------------- */

MLIAPDescriptorSO3::~MLIAPDescriptorSO3()
{
  memory->destroy(radelem);
  memory->destroy(wjelem);
  delete so3ptr;
}

/* ---------------------------------------------------------------------- */

void MLIAPDescriptorSO3::read_paramfile(char *paramfilename)
{
  int rcutfacflag = 0;
  int nelementsflag = 0;
  int elementsflag = 0;
  int radelemflag = 0;
  int wjelemflag = 0;
  int nmaxflag = 0;
  int lmaxflag = 0;
  int alphaflag = 0;

  // set defaults for optional keywords

  rfac0 = 0.99363;
  rmin0 = 0.0;

  for (int i = 0; i < nelements; i++) delete[] elements[i];
  delete[] elements;
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);

  // open SO3 parameter file on proc 0

  FILE *fpparam;
  if (comm->me == 0) {
    fpparam = utils::open_potential(paramfilename, lmp, nullptr);
    if (fpparam == nullptr)
      error->one(FLERR, "Cannot open SO3 parameter file {}: {}", paramfilename,
                 utils::getsyserror());
  }

  char line[MAXLINE], *ptr;
  int eof = 0;
  int n, nwords;

  while (true) {
    if (comm->me == 0) {
      ptr = utils::fgets_trunc(line, MAXLINE, fpparam);
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
    nwords = utils::count_words(line);
    if (nwords == 0) continue;

    // words = ptrs to all words in line
    // strip single and double quotes from words

    Tokenizer p(line, "' \t\n\r\f");
    std::string skeywd = p.next();
    std::string skeyval = p.next();

    if (comm->me == 0) { utils::logmesg(lmp, "SO3 keyword {} {} \n", skeywd, skeyval); }

    // check for keywords with one value per element

    if ((skeywd == "elems") || (skeywd == "radelems") || (skeywd == "welems")) {

      if (nelementsflag == 0 || nwords != nelements + 1)
        error->all(FLERR, "Incorrect SO3 parameter file");

      if (skeywd == "elems") {
        for (int ielem = 0; ielem < nelements; ielem++) {
          elements[ielem] = utils::strdup(skeyval);
          if (ielem < nelements - 1) skeyval = p.next();
        }

        elementsflag = 1;
      } else if (skeywd == "radelems") {
        for (int ielem = 0; ielem < nelements; ielem++) {
          radelem[ielem] = utils::numeric(FLERR, skeyval, false, lmp);
          if (ielem < nelements - 1) skeyval = p.next();
        }
        radelemflag = 1;
      } else if (skeywd == "welems") {
        for (int ielem = 0; ielem < nelements; ielem++) {
          wjelem[ielem] = utils::numeric(FLERR, skeyval, false, lmp);
          if (ielem < nelements - 1) skeyval = p.next();
        }
        wjelemflag = 1;
      }

    } else {

      // all other keywords take one value

      if (nwords != 2) error->all(FLERR, "Incorrect SO3 parameter file");

      if (skeywd == "nelems") {
        nelements = utils::inumeric(FLERR, skeyval, false, lmp);
        elements = new char *[nelements];
        memory->create(radelem, nelements, "mliap_so3_descriptor:radelem");
        memory->create(wjelem, nelements, "mliap_so3_descriptor:wjelem");
        nelementsflag = 1;
      } else if (skeywd == "rcutfac") {
        rcutfac = utils::numeric(FLERR, skeyval, false, lmp);
        rcutfacflag = 1;
      } else if (skeywd == "nmax") {
        nmax = utils::inumeric(FLERR, skeyval, false, lmp);
        nmaxflag = 1;
      } else if (skeywd == "lmax") {
        lmax = utils::inumeric(FLERR, skeyval, false, lmp);
        lmaxflag = 1;
      } else if (skeywd == "alpha") {
        alpha = utils::numeric(FLERR, skeyval, false, lmp);
        alphaflag = 1;
      } else
        error->all(FLERR, "Incorrect SO3 parameter file");
    }
  }

  if (!rcutfacflag || !nelementsflag || !elementsflag || !radelemflag || !wjelemflag || !nmaxflag ||
      !lmaxflag || !alphaflag)
    error->all(FLERR, "Incorrect SO3 parameter file");

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq, nelements, nelements, "mliap/descriptor/so3:cutsq");
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

/* ---------------------------------------------------------------------- */

void MLIAPDescriptorSO3::compute_descriptors(class MLIAPData *data)
{
  so3ptr->spectrum(data->nlistatoms, data->numneighs, data->jelems, wjelem, data->rij, nmax, lmax,
                   rcutfac, alpha, data->ndescriptors);

  for (int ii = 0; ii < data->nlistatoms; ii++) {
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
      data->descriptors[ii][icoeff] = so3ptr->m_plist_r[ii * (data->ndescriptors) + icoeff];
  }
}

/* ---------------------------------------------------------------------- */

void MLIAPDescriptorSO3::compute_forces(class MLIAPData *data)
{
  bigint npairs = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) npairs += data->numneighs[ii];

  so3ptr->spectrum_dxdr(data->nlistatoms, data->numneighs, data->jelems, wjelem, data->rij, nmax,
                        lmax, rcutfac, alpha, npairs, data->ndescriptors);
  double fij[3];
  double **f = atom->f;
  int ij = 0;

  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int i = data->iatoms[ii];

    // ensure rij, inside, wj, and rcutij are of size jnum

    const int jnum = data->numneighs[ii];

    for (int jj = 0; jj < jnum; jj++) {
      int j = data->jatoms[ij];

      for (int ir = 0; ir < 3; ir++) {
        fij[ir] = 0.0;
        for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
          fij[ir] += data->betas[ii][icoeff] *
              so3ptr->m_dplist_r[(ij * (data->ndescriptors) + icoeff) * 3 + ir];
      }

      f[i][0] += fij[0];
      f[i][1] += fij[1];
      f[i][2] += fij[2];
      f[j][0] -= fij[0];
      f[j][1] -= fij[1];
      f[j][2] -= fij[2];

      // add in global and per-atom virial contributions
      // this is optional and has no effect on force calculation

      if (data->vflag) data->pairmliap->v_tally(i, j, fij, data->rij[ij]);
      ij++;
    }
  }
}

/* ---------------------------------------------------------------------- */

void MLIAPDescriptorSO3::compute_force_gradients(class MLIAPData *data)
{
  bigint npairs = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) npairs += data->numneighs[ii];

  so3ptr->spectrum_dxdr(data->nlistatoms, data->numneighs, data->jelems, wjelem, data->rij, nmax,
                        lmax, rcutfac, alpha, npairs, data->ndescriptors);
  int ij = 0;

  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int i = data->iatoms[ii];

    // ensure rij, inside, wj, and rcutij are of size jnum

    const int jnum = data->numneighs[ii];

    for (int jj = 0; jj < jnum; jj++) {
      int j = data->jatoms[ij];

      for (int inz = 0; inz < data->gamma_nnz; inz++) {
        const int l = data->gamma_row_index[ii][inz];
        const int k = data->gamma_col_index[ii][inz];

        data->gradforce[i][l] +=
            data->gamma[ii][inz] * so3ptr->m_dplist_r[(ij * (data->ndescriptors) + k) * 3];
        data->gradforce[i][l + data->yoffset] +=
            data->gamma[ii][inz] * so3ptr->m_dplist_r[(ij * (data->ndescriptors) + k) * 3 + 1];
        data->gradforce[i][l + data->zoffset] +=
            data->gamma[ii][inz] * so3ptr->m_dplist_r[(ij * (data->ndescriptors) + k) * 3 + 2];
        data->gradforce[j][l] -=
            data->gamma[ii][inz] * so3ptr->m_dplist_r[(ij * (data->ndescriptors) + k) * 3];
        data->gradforce[j][l + data->yoffset] -=
            data->gamma[ii][inz] * so3ptr->m_dplist_r[(ij * (data->ndescriptors) + k) * 3 + 1];
        data->gradforce[j][l + data->zoffset] -=
            data->gamma[ii][inz] * so3ptr->m_dplist_r[(ij * (data->ndescriptors) + k) * 3 + 2];
      }
      ij++;
    }
  }
}

/* ---------------------------------------------------------------------- */

void MLIAPDescriptorSO3::compute_descriptor_gradients(class MLIAPData *data)
{
  bigint npairs = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) npairs += data->numneighs[ii];

  so3ptr->spectrum_dxdr(data->nlistatoms, data->numneighs, data->jelems, wjelem, data->rij, nmax,
                        lmax, rcutfac, alpha, npairs, data->ndescriptors);
  int ij = 0;

  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int jnum = data->numneighs[ii];
    for (int jj = 0; jj < jnum; jj++) {
      for (int k = 0; k < data->ndescriptors; k++) {
        data->graddesc[ij][k][0] = so3ptr->m_dplist_r[(ij * (data->ndescriptors) + k) * 3];
        data->graddesc[ij][k][1] = so3ptr->m_dplist_r[(ij * (data->ndescriptors) + k) * 3 + 1];
        data->graddesc[ij][k][2] = so3ptr->m_dplist_r[(ij * (data->ndescriptors) + k) * 3 + 2];
      }
      ij++;
    }
  }
}

/* ---------------------------------------------------------------------- */

void MLIAPDescriptorSO3::init()
{
  so3ptr->init();
}

/* ---------------------------------------------------------------------- */

double MLIAPDescriptorSO3::memory_usage()
{
  double bytes = MLIAPDescriptor::memory_usage();
  bytes += so3ptr->memory_usage();

  return bytes;
}
