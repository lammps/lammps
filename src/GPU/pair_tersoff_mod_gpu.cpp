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
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
------------------------------------------------------------------------- */

#include "pair_tersoff_mod_gpu.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "gpu_extra.h"
#include "info.h"
#include "lammps_gpu.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;


/* ---------------------------------------------------------------------- */

PairTersoffMODGPU::PairTersoffMODGPU(LAMMPS *lmp) : PairTersoffMOD(lmp), gpu_mode(GPU_FORCE)
{
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);

  cutghost = nullptr;
  ghostneigh = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairTersoffMODGPU::~PairTersoffMODGPU()
{
  tersoff_mod_gpu_clear();
  if (allocated) memory->destroy(cutghost);
}

/* ---------------------------------------------------------------------- */

void PairTersoffMODGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    double sublo[3], subhi[3];
    if (domain->triclinic == 0) {
      sublo[0] = domain->sublo[0];
      sublo[1] = domain->sublo[1];
      sublo[2] = domain->sublo[2];
      subhi[0] = domain->subhi[0];
      subhi[1] = domain->subhi[1];
      subhi[2] = domain->subhi[2];
    } else {
      domain->bbox(domain->sublo_lamda, domain->subhi_lamda, sublo, subhi);
    }
    inum = atom->nlocal;
    firstneigh = tersoff_mod_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo,
                                           subhi, atom->tag, atom->nspecial, atom->special, eflag,
                                           vflag, eflag_atom, vflag_atom, &ilist,
                                           &numneigh, success);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    tersoff_mod_gpu_compute(neighbor->ago, inum, nall, inum + list->gnum, atom->x, atom->type,
                            ilist, numneigh, firstneigh, eflag, vflag, eflag_atom, vflag_atom,
                            success);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");
  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
}

/* ---------------------------------------------------------------------- */

void PairTersoffMODGPU::allocate()
{
  PairTersoffMOD::allocate();
  int n = atom->ntypes;

  memory->create(cutghost, n + 1, n + 1, "pair:cutghost");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairTersoffMODGPU::init_style()
{
  double cell_size = cutmax + neighbor->skin;

  if (atom->tag_enable == 0) error->all(FLERR, "Pair style tersoff/mod/gpu requires atom IDs");

  double *lam1, *lam2, *lam3, *powermint;
  double *biga, *bigb, *bigr, *bigd;
  double *c1, *c2, *c3, *c4, *c5, *h;
  double *beta, *powern, *ca1, *powern_del, *_cutsq;
  lam1 = lam2 = lam3 = powermint = nullptr;
  biga = bigb = bigr = bigd = nullptr;
  powern_del = ca1 = nullptr;
  c1 = c2 = c3 = c4 = c5 = h = nullptr;
  beta = powern = _cutsq = nullptr;

  memory->create(lam1, nparams, "pair:lam1");
  memory->create(lam2, nparams, "pair:lam2");
  memory->create(lam3, nparams, "pair:lam3");
  memory->create(powermint, nparams, "pair:powermint");
  memory->create(biga, nparams, "pair:biga");
  memory->create(bigb, nparams, "pair:bigb");
  memory->create(bigr, nparams, "pair:bigr");
  memory->create(bigd, nparams, "pair:bigd");
  memory->create(c1, nparams, "pair:c1");
  memory->create(c2, nparams, "pair:c2");
  memory->create(c3, nparams, "pair:c3");
  memory->create(c4, nparams, "pair:c4");
  memory->create(c5, nparams, "pair:c5");
  memory->create(h, nparams, "pair:h");
  memory->create(beta, nparams, "pair:beta");
  memory->create(powern, nparams, "pair:powern");
  memory->create(powern_del, nparams, "pair:powern_del");
  memory->create(ca1, nparams, "pair:ca1");
  memory->create(_cutsq, nparams, "pair:_cutsq");

  for (int i = 0; i < nparams; i++) {
    lam1[i] = params[i].lam1;
    lam2[i] = params[i].lam2;
    lam3[i] = params[i].lam3;
    powermint[i] = params[i].powermint;
    biga[i] = params[i].biga;
    bigb[i] = params[i].bigb;
    bigr[i] = params[i].bigr;
    bigd[i] = params[i].bigd;
    c1[i] = params[i].c1;
    c2[i] = params[i].c2;
    c3[i] = params[i].c3;
    c4[i] = params[i].c4;
    c5[i] = params[i].c5;
    h[i] = params[i].h;
    beta[i] = params[i].beta;
    powern[i] = params[i].powern;
    powern_del[i] = params[i].powern_del;
    ca1[i] = params[i].ca1;
    _cutsq[i] = params[i].cutsq;
  }

  int mnf = 5e-2 * neighbor->oneatom;
  int success = tersoff_mod_gpu_init(atom->ntypes + 1, atom->nlocal, atom->nlocal + atom->nghost,
                                     mnf, cell_size, gpu_mode, screen, map, nelements, elem3param,
                                     nparams, lam1, lam2, lam3, powermint, biga, bigb, bigr, bigd,
                                     c1, c2, c3, c4, c5, h, beta, powern, powern_del, ca1, _cutsq);

  memory->destroy(lam1);
  memory->destroy(lam2);
  memory->destroy(lam3);
  memory->destroy(powermint);
  memory->destroy(biga);
  memory->destroy(bigb);
  memory->destroy(bigr);
  memory->destroy(bigd);
  memory->destroy(c1);
  memory->destroy(c2);
  memory->destroy(c3);
  memory->destroy(c4);
  memory->destroy(c5);
  memory->destroy(h);
  memory->destroy(beta);
  memory->destroy(powern);
  memory->destroy(ca1);
  memory->destroy(powern_del);
  memory->destroy(_cutsq);

  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE)
    neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
  if (comm->get_comm_cutoff() < (2.0 * cutmax + neighbor->skin)) {
    comm->cutghostuser = 2.0 * cutmax + neighbor->skin;
    if (comm->me == 0)
      error->warning(FLERR, "Increasing communication cutoff to {:.8} for GPU pair style",
                     comm->cutghostuser);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTersoffMODGPU::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, Error::NOLASTLINE,
               "All pair coeffs are not set. Status:\n" + Info::get_pair_coeff_status(lmp));
  cutghost[i][j] = cutmax;
  cutghost[j][i] = cutmax;

  return cutmax;
}
