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
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "reaxff_library.h"

#define LAMMPS_LIB_MPI 1
#include <mpi.h>

#include "error.h"
#include "force.h"
#include "pair_reaxff.h"
#include "reaxff_api.h"

#include <cstdio>

#if defined(LMP_PYTHON)
#include <Python.h>
#endif

/// string buffer for error messages of global errors
static std::string lammps_last_global_errormessage;

using namespace LAMMPS_NS;

// __func__ is supposed to be portable for compilers also supporting C99: MSVC, GCC, Clang.
#if defined(FNERR)
#undef FNERR
#endif
#define FNERR __func__

// ----------------------------------------------------------------------
// utility macros
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   macros for optional code path which captures all exceptions
   and stores the last error message. These assume there is a variable lmp
   which is a pointer to the current LAMMPS instance.

   Usage:

   BEGIN_CAPTURE
   {
     // code paths which might throw exception
     ...
   }
   END_CAPTURE
------------------------------------------------------------------------- */

#define BEGIN_CAPTURE \
  Error *error = lmp->error; \
  try

#define END_CAPTURE \
  catch(LAMMPSAbortException &ae) { \
    int nprocs = 0; \
    MPI_Comm_size(ae.get_universe(), &nprocs ); \
    \
    if (nprocs > 1) { \
      error->set_last_error(ae.what(), ERROR_ABORT); \
    } else { \
      error->set_last_error(ae.what(), ERROR_NORMAL); \
    } \
  } catch(LAMMPSException &e) { \
    error->set_last_error(e.what(), ERROR_NORMAL); \
  }

// ----------------------------------------------------------------------
// Library functions to set reaxff parameters
// ----------------------------------------------------------------------

void lammps_set_reaxff_gen_parameter(void *handle, int parameter_index, double value) {

  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    PairReaxFF *reaxff = static_cast<PairReaxFF *>(lmp->force->pair);
    auto reax = &(reaxff->api->system->reax_param);
    auto &gp = reax->gp;

    // only used for gp.l[34] by fitsnap-reaxff at this time
    // but as general as possible
    gp.l[parameter_index] = value;

    // FIXME: handle gp.l[35] for lg-dispersion later
  }
  END_CAPTURE

}

void lammps_set_reaxff_atm_parameter(void *handle, int type, int parameter_index, double value) {

  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    PairReaxFF *reaxff = static_cast<PairReaxFF *>(lmp->force->pair);
    auto reax = &(reaxff->api->system->reax_param);
    auto &gp = reax->gp;
    auto &sbp = reax->sbp;
    auto &tbp = reax->tbp;
    int i = type - 1;
    int j, k, n = reax->num_atom_types;

    switch(parameter_index) {

      // line one

      case 0:
        sbp[i].r_s = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].r_s=tbp[j][k].r_s=0.5*(sbp[j].r_s+sbp[k].r_s);
        break;

      case 1:
        sbp[i].valency = value;
        sbp[i].nlp_opt = 0.5 * (sbp[i].valency_e-sbp[i].valency);
        break;

      case 2:    sbp[i].mass        = value;      break;

      case 3:
        sbp[i].r_vdw = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].r_vdw=tbp[j][k].r_vdw=2.0*sqrt(sbp[j].r_vdw*sbp[k].r_vdw);
        break;

      case 4:
        sbp[i].epsilon = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].D=tbp[j][k].D=sqrt(sbp[j].epsilon*sbp[k].epsilon);
        break;

      case 5:
        sbp[i].gamma = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].gamma=tbp[j][k].gamma=pow(sbp[j].gamma*sbp[k].gamma,-1.5);
        break;

      case 6:
        sbp[i].r_p = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].r_p=tbp[j][k].r_p=0.5*(sbp[j].r_p+sbp[k].r_p);
        break;

      case 7:
        sbp[i].valency_e   = value;
        sbp[i].nlp_opt = 0.5 * (sbp[i].valency_e-sbp[i].valency);
        break;

      // line two

      case 8:
        sbp[i].alpha = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].alpha=tbp[j][k].alpha=sqrt(sbp[j].alpha*sbp[k].alpha);
        break;

      case 9:
        sbp[i].gamma_w = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].gamma_w=tbp[j][k].gamma_w=sqrt(sbp[j].gamma_w*sbp[k].gamma_w);
        break;

      case 10:
        sbp[i].valency_boc = value;
        if((sbp[i].mass<21)&&(sbp[i].valency_val!=sbp[i].valency_boc)) sbp[i].valency_val=sbp[i].valency_boc;
        break;

      case 11:   sbp[i].p_ovun5     = value;      break;
      //case 12:    values.skip();
      case 13:   sbp[i].chi         = value;      break;
      case 14:   sbp[i].eta         = 2.0*value;  break;
      case 15:   sbp[i].p_hbond     = (int)value; break;

      // line three
      case 16:
        sbp[i].r_pp = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].r_pp=tbp[j][k].r_pp=0.5*(sbp[j].r_pp+sbp[k].r_pp);
        break;

      case 17:   sbp[i].p_lp2       = value;      break;
      //case 18:    values.skip();

      case 19:
        sbp[i].b_o_131 = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].p_boc4=tbp[j][k].p_boc4=sqrt(sbp[j].b_o_131*sbp[k].b_o_131);
        break;

      case 20:
        sbp[i].b_o_132 = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].p_boc3=tbp[j][k].p_boc3=sqrt(sbp[j].b_o_132*sbp[k].b_o_132);
        break;

      case 21:
        sbp[i].b_o_133 = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].p_boc5=tbp[j][k].p_boc5=sqrt(sbp[j].b_o_133*sbp[k].b_o_133);
        break;

      case 22:   sbp[i].bcut_acks2  = value;      break;
      //case 23:    values.skip();

      // line four
      case 24:   sbp[i].p_ovun2     = value;      break;
      case 25:   sbp[i].p_val3      = value;      break;
      //case 26:    values.skip();

      case 27:
        sbp[i].valency_val = value;
        if((sbp[i].mass<21)&&(sbp[i].valency_val!=sbp[i].valency_boc)) sbp[i].valency_val=sbp[i].valency_boc;
        break;

      case 28:   sbp[i].p_val5      = value;      break;

      case 29:
        sbp[i].rcore2 = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].rcore=tbp[j][k].rcore=sqrt(sbp[k].rcore2*sbp[j].rcore2);
        break;

      case 30:
        sbp[i].ecore2 = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].ecore=tbp[j][k].ecore=sqrt(sbp[k].ecore2*sbp[j].ecore2);
        break;

      case 31:
        sbp[i].acore2 = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].acore=tbp[j][k].acore=sqrt(sbp[k].acore2*sbp[j].acore2);
        break;

      case 32:
        sbp[i].lgcij = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].lgcij=tbp[j][k].lgcij=sqrt(sbp[k].lgcij*sbp[j].lgcij);
        break;

      case 33:
        sbp[i].lgre = value;
        for(j=0;j<n;++j) for(k=j;k<n;++k) tbp[k][j].lgre=tbp[j][k].lgre=2.0*gp.l[35]*sqrt(sbp[k].lgre*sbp[j].lgre);
        break;

    }

    // FIXME: van der Waals settings check

  }
  END_CAPTURE

}

void lammps_set_reaxff_bnd_parameter(void *handle, int type1, int type2, int parameter_index, double value) {

  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    PairReaxFF *reaxff = static_cast<PairReaxFF *>(lmp->force->pair);
    auto &tbp = reaxff->api->system->reax_param.tbp;
    int j = type1 - 1;
    int k = type2 - 1;

    switch(parameter_index) {

      // line one
      case 0:   tbp[j][k].De_s    = tbp[k][j].De_s    = value;  break;
      case 1:   tbp[j][k].De_p    = tbp[k][j].De_p    = value;  break;
      case 2:   tbp[j][k].De_pp   = tbp[k][j].De_pp   = value;  break;
      case 3:   tbp[j][k].p_be1   = tbp[k][j].p_be1   = value;  break;
      case 4:   tbp[j][k].p_bo5   = tbp[k][j].p_bo5   = value;  break;
      case 5:   tbp[j][k].v13cor  = tbp[k][j].v13cor  = value;  break;
      case 6:   tbp[j][k].p_bo6   = tbp[k][j].p_bo6   = value;  break;
      case 7:   tbp[j][k].p_ovun1 = tbp[k][j].p_ovun1 = value;  break;

      // line two
      case 8:   tbp[j][k].p_be2   = tbp[k][j].p_be2   = value;  break;
      case 9:   tbp[j][k].p_bo3   = tbp[k][j].p_bo3   = value;  break;
      case 10:  tbp[j][k].p_bo4   = tbp[k][j].p_bo4   = value;  break;
      //case 11:      values.skip();
      case 12:  tbp[j][k].p_bo1   = tbp[k][j].p_bo1   = value;  break;
      case 13:  tbp[j][k].p_bo2   = tbp[k][j].p_bo2   = value;  break;
      case 14:  tbp[j][k].ovc     = tbp[k][j].ovc     = value;  break;

    }



  }
  END_CAPTURE

}

void lammps_set_reaxff_ofd_parameter(void *handle, int type1, int type2, int parameter_index, double value) {

  auto lmp = (LAMMPS *) handle;

  if (value <= 0.0) return;

  BEGIN_CAPTURE
  {
    PairReaxFF *reaxff = static_cast<PairReaxFF *>(lmp->force->pair);
    auto &tbp = reaxff->api->system->reax_param.tbp;
    int j = type1 - 1;
    int k = type2 - 1;

    switch(parameter_index) {
      case 0:  tbp[j][k].D     = tbp[k][j].D     = value;     break;
      case 1:  tbp[j][k].r_vdw = tbp[k][j].r_vdw = 2.0*value; break;
      case 2:  tbp[j][k].alpha = tbp[k][j].alpha = value;     break;
      case 3:  tbp[j][k].r_s   = tbp[k][j].r_s   = value;     break;
      case 4:  tbp[j][k].r_p   = tbp[k][j].r_p   = value;     break;
      case 5:  tbp[j][k].r_pp  = tbp[k][j].r_pp  = value;     break;
      case 6:  tbp[j][k].lgcij = tbp[k][j].lgcij = value;     break;
    }

  }
  END_CAPTURE

}

void lammps_set_reaxff_ang_parameter(void *handle, int type1, int type2, int type3, int parameter_index, double value) {

  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    PairReaxFF *reaxff = static_cast<PairReaxFF *>(lmp->force->pair);
    auto &thbp = reaxff->api->system->reax_param.thbp;
    int j = type1 - 1;
    int k = type2 - 1;
    int l = type3 - 1;

    // check if cnt>1 is possible
    // yes it is, eg. N-N-N gets cnt+=2
    if( thbp[j][k][l].cnt > 2 )
      lmp->error->all(FLERR,
        "lammps_set_reaxff_ang_parameter(): thbp[{}][{}][{}].cnt = {} > 2.", j, k, l, thbp[j][k][l].cnt );

    switch(parameter_index) {
      case 0:  thbp[j][k][l].prm[0].theta_00 = thbp[l][k][j].prm[0].theta_00 = value;  break;
      case 1:  thbp[j][k][l].prm[0].p_val1   = thbp[l][k][j].prm[0].p_val1   = value;  break;
      case 2:  thbp[j][k][l].prm[0].p_val2   = thbp[l][k][j].prm[0].p_val2   = value;  break;
      case 3:  thbp[j][k][l].prm[0].p_coa1   = thbp[l][k][j].prm[0].p_coa1   = value;  break;
      case 4:  thbp[j][k][l].prm[0].p_val7   = thbp[l][k][j].prm[0].p_val7   = value;  break;
      case 5:  thbp[j][k][l].prm[0].p_pen1   = thbp[l][k][j].prm[0].p_pen1   = value;  break;
      case 6:  thbp[j][k][l].prm[0].p_val4   = thbp[l][k][j].prm[0].p_val4   = value;  break;
    }

  }
  END_CAPTURE

}

template <int parameter_index>
void _set_tor_parameter(void *handle, int j, int k, int l, int m, double value) {

  auto lmp = static_cast<LAMMPS *>(handle);
  BEGIN_CAPTURE
  {
    auto *reaxff = static_cast<PairReaxFF *>(lmp->force->pair);
    auto &param = reaxff->api->system->reax_param;
    auto &fbp = param.fbp;

    if (j >= 0 && m >= 0) {
      if constexpr (parameter_index == 0) fbp[j][k][l][m].prm[0].V1     = fbp[m][l][k][j].prm[0].V1     = value;
      if constexpr (parameter_index == 1) fbp[j][k][l][m].prm[0].V2     = fbp[m][l][k][j].prm[0].V2     = value;
      if constexpr (parameter_index == 2) fbp[j][k][l][m].prm[0].V3     = fbp[m][l][k][j].prm[0].V3     = value;
      if constexpr (parameter_index == 3) fbp[j][k][l][m].prm[0].p_tor1 = fbp[m][l][k][j].prm[0].p_tor1 = value;
      if constexpr (parameter_index == 4) fbp[j][k][l][m].prm[0].p_cot1 = fbp[m][l][k][j].prm[0].p_cot1 = value;
    } else {
      auto &tor_flag = param.tor_flag;
      const int ntypes = param.num_atom_types;
      for (int p = 0; p < ntypes; ++p) {
        for (int o = 0; o < ntypes; ++o) {
          if( tor_flag[p][k][l][o]==0 ) {
            if constexpr (parameter_index == 0) fbp[p][k][l][o].prm[0].V1     = value;
            if constexpr (parameter_index == 1) fbp[p][k][l][o].prm[0].V2     = value;
            if constexpr (parameter_index == 2) fbp[p][k][l][o].prm[0].V3     = value;
            if constexpr (parameter_index == 3) fbp[p][k][l][o].prm[0].p_tor1 = value;
            if constexpr (parameter_index == 4) fbp[p][k][l][o].prm[0].p_cot1 = value;
          }
          if( tor_flag[o][l][k][p]==0 ) {
            if constexpr (parameter_index == 0) fbp[o][l][k][p].prm[0].V1     = value;
            if constexpr (parameter_index == 1) fbp[o][l][k][p].prm[0].V2     = value;
            if constexpr (parameter_index == 2) fbp[o][l][k][p].prm[0].V3     = value;
            if constexpr (parameter_index == 3) fbp[o][l][k][p].prm[0].p_tor1 = value;
            if constexpr (parameter_index == 4) fbp[o][l][k][p].prm[0].p_cot1 = value;
          }
        }
      }
    }
  }
  END_CAPTURE
}

void lammps_set_reaxff_tor_parameter(void *handle, int type1, int type2, int type3, int type4, int parameter_index, double value) {

  const int j = type1 - 1;
  const int k = type2 - 1;
  const int l = type3 - 1;
  const int m = type4 - 1;

  switch(parameter_index) {
    case 0: _set_tor_parameter<0>(handle, j, k, l, m, value); break;
    case 1: _set_tor_parameter<1>(handle, j, k, l, m, value); break;
    case 2: _set_tor_parameter<2>(handle, j, k, l, m, value); break;
    case 3: _set_tor_parameter<3>(handle, j, k, l, m, value); break;
    case 4: _set_tor_parameter<4>(handle, j, k, l, m, value); break;
  }
}

void lammps_set_reaxff_hbd_parameter(void *handle, int type1, int type2, int type3, int parameter_index, double value) {

  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    PairReaxFF *reaxff = static_cast<PairReaxFF *>(lmp->force->pair);
    auto &hbp = reaxff->api->system->reax_param.hbp;
    int j = type1 - 1;
    int k = type2 - 1;
    int l = type3 - 1;

    switch(parameter_index) {
      case 0:  hbp[j][k][l].r0_hb = value;  break;
      case 1:  hbp[j][k][l].p_hb1 = value;  break;
      case 2:  hbp[j][k][l].p_hb2 = value;  break;
      case 3:  hbp[j][k][l].p_hb3 = value;  break;
    }

  }
  END_CAPTURE

}

void lammps_reset_box_single(void *handle, double *boxlo, double *boxhi,
                      double xy, double yz, double xz)
{
  auto *lmp = (LAMMPS *) handle;
  if (!lmp || !lmp->error || !lmp->domain || !lmp->comm) {
    lammps_last_global_errormessage = fmt::format("ERROR: {}(): Invalid LAMMPS handle\n", FNERR);
    return;
  }
  Domain *domain = lmp->domain;

  BEGIN_CAPTURE
  {

    if (lmp->comm->nprocs > 1)
      lmp->error->all(FLERR, Error::NOLASTLINE,
                      "Calling lammps_reset_box_single() with more than one processor");

    // warn and do nothing if no box exists
    if (lmp->domain->box_exist == 0) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Call to lammps_reset_box() without a box ignored");
      return;
    }

    domain->boxlo[0] = boxlo[0];
    domain->boxlo[1] = boxlo[1];
    domain->boxlo[2] = boxlo[2];
    domain->boxhi[0] = boxhi[0];
    domain->boxhi[1] = boxhi[1];
    domain->boxhi[2] = boxhi[2];

    domain->xy = xy;
    domain->yz = yz;
    domain->xz = xz;

    domain->set_global_box();
    //lmp->comm->set_proc_grid();
    domain->set_local_box();
  }
  END_CAPTURE
}

