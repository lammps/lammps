/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
------------------------------------------------------------------------- */

#include "pair_mbx.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

#include "fix_mbx.h"
#include "bblock/system.h"

#define TTMNRG

// subject for removal
// Systems::DispersionPME()
// Systems::ElectrostaticsMPI()

namespace LAMMPS_NS{
//PImpl idiom to hide MBX implementation details
struct MBXImpl {
  MBXImpl() : ptr_mbx(nullptr), ptr_mbx_local(nullptr) {}
  ~MBXImpl()
  {
    delete ptr_mbx;
    delete ptr_mbx_local;
  }
  bblock::System *ptr_mbx;
  bblock::System *ptr_mbx_local;
};
} // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairMBX::PairMBX(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  no_virial_fdotr_compute = 1;
  one_coeff = 1;

  mbx_total_energy = 0;

  me = comm->me;

  // energy terms available to pair compute

  nextra = 10;
  pvector = new double[nextra];
}

/* ---------------------------------------------------------------------- */

PairMBX::~PairMBX()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }

  delete[] pvector;
}

/* ---------------------------------------------------------------------- */

void PairMBX::compute(int eflag, int vflag)
{

  ev_init(eflag, vflag);

  // compute energy+gradients in parallel

  bblock::System *ptr_mbx = fix_MBX->mbx_impl->ptr_mbx;              // compute terms in parallel
  bblock::System *ptr_mbx_local =
      fix_MBX->mbx_impl->ptr_mbx_local;    // compute PME terms in parallel w/ sub-domains

  double mbx_e2b_local, mbx_e2b_ghost;
  double mbx_e3b_local, mbx_e3b_ghost;
  double mbx_e4b_local, mbx_e4b_ghost;
  double mbx_disp_real, mbx_disp_pme;

  // compute energy

  mbx_e1b = 0.0;
  mbx_e2b = 0.0;
  mbx_e3b = 0.0;
  mbx_e4b = 0.0;

  mbx_e2b_local = 0.0;
  mbx_e2b_ghost = 0.0;
  mbx_e3b_local = 0.0;
  mbx_e3b_ghost = 0.0;
  mbx_e4b_local = 0.0;
  mbx_e4b_ghost = 0.0;

  mbx_disp_real = 0.0;
  mbx_disp_pme = 0.0;

  mbx_disp = 0.0;
  mbx_buck = 0.0;
  mbx_ele = 0.0;

  for (int i = 0; i < 6; ++i) mbx_virial[i] = 0.0;

  if (fix_MBX->mbx_num_atoms > 0) {

    fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::E1B);
    mbx_e1b = ptr_mbx->OneBodyEnergy(true);
    fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::E1B);
    accumulate_f(false);

    fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::E2B_GHOST);
    mbx_e2b_ghost = ptr_mbx->TwoBodyEnergy(true, true);
    fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::E2B_GHOST);
    accumulate_f_all(false);
    mbx_e2b = mbx_e2b_ghost;

    fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::E3B_GHOST);
    mbx_e3b_ghost = ptr_mbx->ThreeBodyEnergy(true, true);
    fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::E3B_GHOST);
    accumulate_f_all(false);
    mbx_e3b = mbx_e3b_ghost;

    fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::E4B_GHOST);
    mbx_e4b_ghost = ptr_mbx->FourBodyEnergy(true, true);
    fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::E4B_GHOST);
    accumulate_f_all(false);

    mbx_e4b = mbx_e4b_ghost;
  }


  fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::ELE);
  mbx_ele = ptr_mbx_local->ElectrostaticsMPIlocal(true, true);
  fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::ELE);
  accumulate_f_local(true);


  fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::DISP);
  mbx_disp_real = ptr_mbx_local->Dispersion(
      true, true);    // computes real-space with local-local & local-ghost pairs
  fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::DISP);
  accumulate_f_local(false);

  fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::DISP_PME);
  mbx_disp_pme = ptr_mbx_local->DispersionPMElocal(
      true, true);    // computes PME-space with local-local & local-ghost pairs
  fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::DISP_PME);
  accumulate_f_local(false);


  if (fix_MBX->mbx_num_atoms > 0) {

#ifdef TTMNRG
    fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::BUCK);
    mbx_buck = ptr_mbx->Buckingham(true, true);
    fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::BUCK);
    accumulate_f(false);

    fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::BUCK);
    mbx_buck += ptr_mbx_local->LennardJones(true, true);
    fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::BUCK);
    accumulate_f_local(false);
#else
    fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::BUCK);
    mbx_buck = 0.0;
    fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::BUCK);
#endif
  }

  mbx_disp = mbx_disp_real + mbx_disp_pme;

  mbx_total_energy = mbx_e1b + mbx_e2b + mbx_disp + mbx_buck + mbx_e3b + mbx_e4b + mbx_ele;

  for (int i = 0; i < 6; ++i) virial[i] += mbx_virial[i];

  // save total energy from mbx as vdwl

  if (evflag) {
    eng_vdwl = mbx_total_energy;

    // generally useful

    pvector[0] = mbx_e1b;
    pvector[1] = mbx_e2b;
    pvector[2] = mbx_e3b;
    pvector[3] = mbx_e4b;
    pvector[4] = mbx_disp;
    pvector[5] = mbx_buck;
    pvector[6] = ptr_mbx_local->GetPermanentElectrostaticEnergy();
    pvector[7] = ptr_mbx_local->GetInducedElectrostaticEnergy();
    pvector[8] = mbx_ele;
    pvector[9] = mbx_total_energy;

  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMBX::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMBX::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMBX::coeff(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Incorrect num args for pair coefficients");
  if (!allocated) allocate();

  int count = 0;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");

  std::string fix_args = "";
  for (int i = 2; i < narg; ++i) { fix_args += std::string(arg[i]) + " "; }

  fix_args = fmt::format("_FIX_MBX_INTERNAL all MBX {}", fix_args);

  fix_MBX = dynamic_cast<FixMBX *>(modify->add_fix(fix_args));
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMBX::init_one(int i, int j)
{
  return cut_global;
}


/* ----------------------------------------------------------------------
   update forces with MBX contribution
------------------------------------------------------------------------- */

void PairMBX::accumulate_f(bool include_ext)
{

  fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::ACCUMULATE_F);

  bblock::System *ptr_mbx = fix_MBX->mbx_impl->ptr_mbx;

  const int nlocal = atom->nlocal;
  double **f = atom->f;

  const int *const mol_anchor = fix_MBX->mol_anchor;
  const int *const mol_type = fix_MBX->mol_type;
  char **mol_names = fix_MBX->mol_names;

  std::vector<double> grads = ptr_mbx->GetRealGrads();

  std::vector<double> grads_ext;
  if (include_ext)
    grads_ext = ptr_mbx->GetExternalChargesGradients();
  else
    grads_ext = std::vector<double>(fix_MBX->mbx_num_ext * 3, 0.0);

  // accumulate forces on local particles
  // -- forces on ghost particles ignored/not needed
  // -- should use a map created from earlier loop loading particles into mbx

  int indx = 0;
  int indx_ext = 0;

  for (int i = 0; i < nlocal; ++i) {
    if (mol_anchor[i]) {
      const int mtype = mol_type[i];

      // to be replaced with integer comparison

      bool include_monomer = true;
      bool is_ext = false;
      tagint anchor = atom->tag[i];

      int na = fix_MBX->get_include_monomer(mol_names[mtype], anchor, include_monomer, is_ext);

      if (include_monomer) {
        for (int j = 0; j < na; ++j) {
          const int ii = atom->map(anchor + j);
          f[ii][0] -= grads[indx++];
          f[ii][1] -= grads[indx++];
          f[ii][2] -= grads[indx++];
        }
      } else if (is_ext) {
        f[i][0] -= grads_ext[indx_ext++];
        f[i][1] -= grads_ext[indx_ext++];
        f[i][2] -= grads_ext[indx_ext++];
      }

    }    // if(anchor)
  }

  // accumulate virial: only global is supported
  // MBX: xx, xy, xz, yx, yy, yz, zx, zy, zz
  // LAMMPS: xx, yy, zz, xy, xz, yz

  if (vflag_either) {
    std::vector<double> mbx_vir = ptr_mbx->GetVirial();

    mbx_virial[0] += mbx_vir[0];
    mbx_virial[1] += mbx_vir[4];
    mbx_virial[2] += mbx_vir[8];
    mbx_virial[3] += mbx_vir[1];
    mbx_virial[4] += mbx_vir[2];
    mbx_virial[5] += mbx_vir[5];
  }

  fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::ACCUMULATE_F);
}

/* ----------------------------------------------------------------------
   update forces with MBX contribution
------------------------------------------------------------------------- */

void PairMBX::accumulate_f_all(bool include_ext)
{

  fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::ACCUMULATE_F);

  bblock::System *ptr_mbx = fix_MBX->mbx_impl->ptr_mbx;

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  double **f = atom->f;

  const int *const mol_anchor = fix_MBX->mol_anchor;
  const int *const mol_type = fix_MBX->mol_type;
  char **mol_names = fix_MBX->mol_names;

  std::vector<double> grads = ptr_mbx->GetRealGrads();

  std::vector<double> grads_ext;
  if (include_ext)
    grads_ext = ptr_mbx->GetExternalChargesGradients();
  else
    grads_ext = std::vector<double>(fix_MBX->mbx_num_ext * 3, 0.0);

  // accumulate forces on local + ghost particles
  // -- should use a map created from earlier loop loading particles into mbx

  int indx = 0;
  int indx_ext = 0;

  for (int i = 0; i < nall; ++i) {
    if (mol_anchor[i]) {
      const int mtype = mol_type[i];

      // to be replaced with integer comparison

      bool include_monomer = true;
      bool is_ext = false;
      tagint anchor = atom->tag[i];

      int na = fix_MBX->get_include_monomer(mol_names[mtype], anchor, include_monomer, is_ext);

      if (include_monomer) {
        for (int j = 0; j < na; ++j) {
          const int ii = atom->map(anchor + j);
          f[ii][0] -= grads[indx++];
          f[ii][1] -= grads[indx++];
          f[ii][2] -= grads[indx++];
        }
      } else if (is_ext) {
        f[i][0] -= grads_ext[indx_ext++];
        f[i][1] -= grads_ext[indx_ext++];
        f[i][2] -= grads_ext[indx_ext++];
      }

    }    // if(anchor)
  }

  // accumulate virial: only global is supported
  // MBX: xx, xy, xz, yx, yy, yz, zx, zy, zz
  // LAMMPS: xx, yy, zz, xy, xz, yz

  if (vflag_either) {
    std::vector<double> mbx_vir = ptr_mbx->GetVirial();

    mbx_virial[0] += mbx_vir[0];
    mbx_virial[1] += mbx_vir[4];
    mbx_virial[2] += mbx_vir[8];
    mbx_virial[3] += mbx_vir[1];
    mbx_virial[4] += mbx_vir[2];
    mbx_virial[5] += mbx_vir[5];
  }

  fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::ACCUMULATE_F);
}

/* ----------------------------------------------------------------------
   update forces with MBX contribution
------------------------------------------------------------------------- */

void PairMBX::accumulate_f_local(bool include_ext)
{

  fix_MBX->mbxt_start(FixMBX::MBXT_LABELS::ACCUMULATE_F_LOCAL);

  bblock::System *ptr_mbx = fix_MBX->mbx_impl->ptr_mbx_local;

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  double **f = atom->f;

  const int *const mol_anchor = fix_MBX->mol_anchor;
  const int *const mol_local = fix_MBX->mol_local;
  const int *const mol_type = fix_MBX->mol_type;
  char **mol_names = fix_MBX->mol_names;

  std::vector<double> grads = ptr_mbx->GetRealGrads();

  std::vector<double> grads_ext;
  if (include_ext)
    grads_ext = ptr_mbx->GetExternalChargesGradients();
  else
    grads_ext = std::vector<double>(fix_MBX->mbx_num_ext_local * 3, 0.0);

  // accumulate forces on monomers with at least one local particle
  // -- forces on ghost particles ignored/not needed ??
  // -- should use a map created from earlier loop loading particles into mbx

  int indx = 0;
  int indx_ext = 0;

  for (int i = 0; i < nall; ++i) {
    if (mol_anchor[i] && mol_local[i]) {
      const int mtype = mol_type[i];

      // to be replaced with integer comparison

      bool include_monomer = true;
      bool is_ext = false;
      tagint anchor = atom->tag[i];

      int na = fix_MBX->get_include_monomer(mol_names[mtype], anchor, include_monomer, is_ext);

      if (include_monomer) {
        for (int j = 0; j < na; ++j) {
          const int ii = atom->map(anchor + j);
          f[ii][0] -= grads[indx++];
          f[ii][1] -= grads[indx++];
          f[ii][2] -= grads[indx++];
        }
      } else if (is_ext) {
        f[i][0] -= grads_ext[indx_ext++];
        f[i][1] -= grads_ext[indx_ext++];
        f[i][2] -= grads_ext[indx_ext++];
      }

    }    // if(anchor)
  }

  // accumulate virial: only global is supported
  // MBX: xx, xy, xz, yx, yy, yz, zx, zy, zz
  // LAMMPS: xx, yy, zz, xy, xz, yz

  if (vflag_either) {
    std::vector<double> mbx_vir = ptr_mbx->GetVirial();

    mbx_virial[0] += mbx_vir[0];
    mbx_virial[1] += mbx_vir[4];
    mbx_virial[2] += mbx_vir[8];
    mbx_virial[3] += mbx_vir[1];
    mbx_virial[4] += mbx_vir[2];
    mbx_virial[5] += mbx_vir[5];
  }

  fix_MBX->mbxt_stop(FixMBX::MBXT_LABELS::ACCUMULATE_F_LOCAL);
}

