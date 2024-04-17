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
   Contributing author: James Goff (SNL)
------------------------------------------------------------------------- */
#ifdef MLIAP_ACE

#include "mliap_descriptor_ace.h"

#include "ace-evaluator/ace_abstract_basis.h"
#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_types.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "mliap_data.h"
#include "neigh_list.h"
#include "pair_mliap.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

namespace LAMMPS_NS {
struct ACE_ML_impl {
  ACE_ML_impl() : basis_set(nullptr), ace(nullptr) {}
  ~ACE_ML_impl()
  {
    delete basis_set;
    delete ace;
  }
  ACECTildeBasisSet *basis_set;
  ACECTildeEvaluator *ace;
};
}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MLIAPDescriptorACE::MLIAPDescriptorACE(LAMMPS *_lmp, char *yacefilename) :
    Pointers(_lmp), MLIAPDescriptor(_lmp)
{

  acemlimpl = new ACE_ML_impl;
  allocated_elements = 0;
  //read in file with CG coefficients or c_tilde coefficients
  ctilde_file = yacefilename;
  delete acemlimpl->basis_set;
  acemlimpl->basis_set = new ACECTildeBasisSet(ctilde_file);
  nelements = acemlimpl->basis_set->nelements;
  int tot_num = 0;
  for (int mu = 0; mu < nelements; mu++) {
    if (max_num < acemlimpl->basis_set->total_basis_size_rank1[mu] +
            acemlimpl->basis_set->total_basis_size[mu]) {
      max_num = acemlimpl->basis_set->total_basis_size_rank1[mu] +
          acemlimpl->basis_set->total_basis_size[mu];
    }
    tot_num += acemlimpl->basis_set->total_basis_size_rank1[mu] +
        acemlimpl->basis_set->total_basis_size[mu];
  }

  ndescriptors = max_num;
  nelements = acemlimpl->basis_set->nelements;

  memory->destroy(cutsq);

  if (allocated_elements) {
    for (int iielem = 0; iielem < nelements; iielem++) delete[] elements[iielem];
    delete[] elements;
    allocated_elements = 0;
  } else if (not allocated_elements) {
    elements = new char *[nelements];
    for (int iielem = 0; iielem < nelements; iielem++) {
      elements[iielem] = utils::strdup(acemlimpl->basis_set->elements_name[iielem]);
    }
    allocated_elements = 1;
  }

  memory->create(cutsq, nelements + 1, nelements + 1, "mliap/descriptor/ace:cutsq");
  float icmax = 0.0;
  float icuti, icutj;
  for (int mui = 0; mui < acemlimpl->basis_set->nelements; mui++) {
    icuti = acemlimpl->basis_set->radial_functions->cut(mui, mui);
    if (icuti > icmax) icmax = icuti;
    for (int muj = mui + 1; muj < acemlimpl->basis_set->nelements; muj++) {
      icutj = acemlimpl->basis_set->radial_functions->cut(mui, muj);
      if (icutj > icmax) icmax = icutj;
    }
  }
  float cutmax = 0.0;
  float cuti, cutj;
  float cutfac = 1.0;
  for (int mui = 0; mui < acemlimpl->basis_set->nelements; mui++) {
    cuti = acemlimpl->basis_set->radial_functions->cut(mui, mui);
    if (cuti > cutmax) cutmax = cuti;

    cutsq[mui][mui] = ((2 * cuti * cutfac) * (2 * cuti * cutfac));
    for (int muj = mui + 1; muj < nelements; muj++) {
      cutj = acemlimpl->basis_set->radial_functions->cut(mui, muj);
      cutsq[mui][muj] = cutsq[muj][mui] = ((2 * cuti * cutfac) * (2 * cutj * cutfac));
    }
  }
}

void MLIAPDescriptorACE::allocate() {}

/* ---------------------------------------------------------------------- */

MLIAPDescriptorACE::~MLIAPDescriptorACE()
{
  delete acemlimpl;
}

/* ----------------------------------------------------------------------
   compute descriptors for each atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorACE::compute_descriptors(class MLIAPData *data)
{
  int max_jnum = -1;
  int nei = 0;
  int jtmp = 0;
  for (int iitmp = 0; iitmp < data->nlistatoms; iitmp++) {
    int itmp = data->iatoms[iitmp];
    jtmp = data->numneighs[iitmp];
    nei = nei + jtmp;
    if (jtmp > max_jnum) { max_jnum = jtmp; }
  }

  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int i = data->iatoms[ii];
    const int ielemx = data->ielems[ii];
    const int jnum = data->numneighs[ii];

    delete acemlimpl->ace;
    acemlimpl->ace = new ACECTildeEvaluator(*acemlimpl->basis_set);
    acemlimpl->ace->compute_projections = 1;
    acemlimpl->ace->compute_b_grad = 1;

    acemlimpl->ace->element_type_mapping.init(nelements + 1);
    for (int ik = 1; ik <= nelements; ik++) {
      for (int mu = 0; mu < nelements; mu++) {
        if (mu != -1) {
          if (mu == ik - 1) { acemlimpl->ace->element_type_mapping(ik) = mu; }
        }
      }
    }

    acemlimpl->ace->resize_neighbours_cache(jnum);
    acemlimpl->ace->compute_atom(i, atom->x, atom->type, data->numneighs[ii],
                                 data->lmp_firstneigh[ii]);
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
      data->descriptors[ii][icoeff] = acemlimpl->ace->projections(icoeff);
    }
  }
}

/* ----------------------------------------------------------------------
   compute forces for each atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorACE::compute_forces(class MLIAPData *data)
{
  double fij[3];
  double **f = atom->f;
  int ij = 0;

  int max_jnum = -1;
  int nei = 0;
  int jtmp = 0;
  for (int iitmp = 0; iitmp < data->nlistatoms; iitmp++) {
    int itmp = data->iatoms[iitmp];
    jtmp = data->numneighs[iitmp];
    nei = nei + jtmp;
    if (jtmp > max_jnum) { max_jnum = jtmp; }
  }

  // BEGIN force loop
  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int i = data->iatoms[ii];
    const int ielem = data->ielems[ii];
    delete acemlimpl->ace;
    acemlimpl->ace = new ACECTildeEvaluator(*acemlimpl->basis_set);
    acemlimpl->ace->compute_projections = 1;
    acemlimpl->ace->compute_b_grad = 1;
    acemlimpl->ace->element_type_mapping.init(nelements + 1);
    for (int ik = 1; ik <= nelements; ik++) {
      for (int mu = 0; mu < acemlimpl->basis_set->nelements; mu++) {
        if (mu != -1) {
          if (mu == ik - 1) { acemlimpl->ace->element_type_mapping(ik) = mu; }
        }
      }
    }

    const int jnum = data->numneighs[ii];
    acemlimpl->ace->resize_neighbours_cache(jnum);
    acemlimpl->ace->compute_atom(i, atom->x, atom->type, data->numneighs[ii],
                                 data->lmp_firstneigh[ii]);
    int ij0 = ij;
    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      ninside++;
      ij++;
    }

    ij = ij0;
    const int *const jlist = data->lmp_firstneigh[ii];
    double **x = atom->x;
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];

    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj];
      for (int idim = 0; idim < 3; idim++) { fij[idim] = 0.0; }
      for (int iicoeff = 0; iicoeff < ndescriptors; iicoeff++) {
        DOUBLE_TYPE fx_dB =
            acemlimpl->ace->neighbours_dB(iicoeff, jj, 0) * data->betas[ii][iicoeff];
        DOUBLE_TYPE fy_dB =
            acemlimpl->ace->neighbours_dB(iicoeff, jj, 1) * data->betas[ii][iicoeff];
        DOUBLE_TYPE fz_dB =
            acemlimpl->ace->neighbours_dB(iicoeff, jj, 2) * data->betas[ii][iicoeff];
        // add force contribution from each descriptor
        fij[0] += fx_dB;
        fij[1] += fy_dB;
        fij[2] += fz_dB;
        f[i][0] += fx_dB;
        f[i][1] += fy_dB;
        f[i][2] += fz_dB;
        f[j][0] -= fx_dB;
        f[j][1] -= fy_dB;
        f[j][2] -= fz_dB;
      }
      const double delx = x[j][0] - xtmp;
      const double dely = x[j][1] - ytmp;
      const double delz = x[j][2] - ztmp;
      double rij_tmp[3] = {delx, dely, delz};
      if (data->vflag) data->pairmliap->v_tally(i, j, fij, rij_tmp);
      ij++;
    }
  }
}

/* ----------------------------------------------------------------------
   calculate gradients of forces w.r.t. parameters
   ---------------------------------------------------------------------- */

void MLIAPDescriptorACE::compute_force_gradients(class MLIAPData *data)
{
  int ij = 0;

  int max_jnum = -1;
  int nei = 0;
  int jtmp = 0;
  for (int iitmp = 0; iitmp < data->nlistatoms; iitmp++) {
    int itmp = data->iatoms[iitmp];
    jtmp = data->numneighs[iitmp];
    nei = nei + jtmp;
    if (jtmp > max_jnum) { max_jnum = jtmp; }
  }

  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int i = data->iatoms[ii];
    const int ielem = data->ielems[ii];
    delete acemlimpl->ace;
    acemlimpl->ace = new ACECTildeEvaluator(*acemlimpl->basis_set);
    acemlimpl->ace->compute_projections = 1;
    acemlimpl->ace->compute_b_grad = 1;

    acemlimpl->ace->element_type_mapping.init(nelements + 1);
    for (int ik = 1; ik <= nelements; ik++) {
      for (int mu = 0; mu < acemlimpl->basis_set->nelements; mu++) {
        if (mu != -1) {
          if (mu == ik - 1) { acemlimpl->ace->element_type_mapping(ik) = mu; }
        }
      }
    }

    const int jnum = data->numneighs[ii];
    const int *const jlist = data->lmp_firstneigh[ii];
    acemlimpl->ace->resize_neighbours_cache(jnum);
    acemlimpl->ace->compute_atom(i, atom->x, atom->type, data->numneighs[ii],
                                 data->lmp_firstneigh[ii]);
    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj];
      for (int inz = 0; inz < data->gamma_nnz; inz++) {
        const int l = data->gamma_row_index[ii][inz];
        const int k = data->gamma_col_index[ii][inz];
        DOUBLE_TYPE fx_dB = acemlimpl->ace->neighbours_dB(k, jj, 0);
        DOUBLE_TYPE fy_dB = acemlimpl->ace->neighbours_dB(k, jj, 1);
        DOUBLE_TYPE fz_dB = acemlimpl->ace->neighbours_dB(k, jj, 2);
        data->gradforce[i][l] += data->gamma[ii][inz] * fx_dB;
        data->gradforce[i][l + data->yoffset] += data->gamma[ii][inz] * fy_dB;
        data->gradforce[i][l + data->zoffset] += data->gamma[ii][inz] * fz_dB;
        data->gradforce[j][l] -= data->gamma[ii][inz] * fx_dB;
        data->gradforce[j][l + data->yoffset] -= data->gamma[ii][inz] * fy_dB;
        data->gradforce[j][l + data->zoffset] -= data->gamma[ii][inz] * fz_dB;
      }
      ij++;
    }
  }
}

/* ----------------------------------------------------------------------
   compute descriptor gradients for each neighbor atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorACE::compute_descriptor_gradients(class MLIAPData *data)
{
  int ij = 0;
  int max_jnum = -1;
  int nei = 0;
  int jtmp = 0;
  for (int iitmp = 0; iitmp < data->nlistatoms; iitmp++) {
    int itmp = data->iatoms[iitmp];
    jtmp = data->numneighs[iitmp];
    nei = nei + jtmp;
    if (jtmp > max_jnum) { max_jnum = jtmp; }
  }
  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int i = data->iatoms[ii];
    const int ielem = data->ielems[ii];
    delete acemlimpl->ace;
    acemlimpl->ace = new ACECTildeEvaluator(*acemlimpl->basis_set);
    acemlimpl->ace->compute_projections = 1;
    acemlimpl->ace->compute_b_grad = 1;

    acemlimpl->ace->element_type_mapping.init(nelements + 1);
    for (int ik = 1; ik <= nelements; ik++) {
      for (int mu = 0; mu < acemlimpl->basis_set->nelements; mu++) {
        if (mu != -1) {
          if (mu == ik - 1) { acemlimpl->ace->element_type_mapping(ik) = mu; }
        }
      }
    }

    const int *const jlist = data->lmp_firstneigh[ii];
    const int jnum = data->numneighs[ii];
    acemlimpl->ace->resize_neighbours_cache(jnum);
    acemlimpl->ace->compute_atom(i, atom->x, atom->type, data->numneighs[ii],
                                 data->lmp_firstneigh[ii]);
    int ij0 = ij;
    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      ninside++;
      ij++;
    }

    ij = ij0;
    for (int jj = 0; jj < data->numneighs[ii]; jj++) {
      const int jt = data->jatoms[ij];

      const int j = jlist[jj];
      int yoffset = ndescriptors;
      int zoffset = ndescriptors * 2;
      for (int iicoeff = 0; iicoeff < ndescriptors; iicoeff++) {
        DOUBLE_TYPE fx_dB = acemlimpl->ace->neighbours_dB(iicoeff, jj, 0);
        DOUBLE_TYPE fy_dB = acemlimpl->ace->neighbours_dB(iicoeff, jj, 1);
        DOUBLE_TYPE fz_dB = acemlimpl->ace->neighbours_dB(iicoeff, jj, 2);
        // Accumulate dB_k^i/dRi, dB_k^i/dRj
        data->graddesc[ij][iicoeff][0] = fx_dB;
        data->graddesc[ij][iicoeff][1] = fy_dB;
        data->graddesc[ij][iicoeff][2] = fz_dB;
      }
      ij++;
    }
  }
}

void MLIAPDescriptorACE::init() {}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double MLIAPDescriptorACE::memory_usage()
{
  double bytes = MLIAPDescriptor::memory_usage();

  return bytes;
}

#endif
