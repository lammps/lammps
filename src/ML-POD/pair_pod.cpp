/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ngoc Cuong Nguyen (MIT) and Andrew Rohskopf (SNL)
------------------------------------------------------------------------- */

#include "pair_pod.h"

#include "mlpod.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairPOD::PairPOD(LAMMPS *lmp) :
    Pair(lmp), gd(nullptr), gdall(nullptr), podcoeff(nullptr), newpodcoeff(nullptr),
    energycoeff(nullptr), forcecoeff(nullptr), podptr(nullptr), tmpmem(nullptr), typeai(nullptr),
    numneighsum(nullptr), rij(nullptr), idxi(nullptr), ai(nullptr), aj(nullptr), ti(nullptr),
    tj(nullptr)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  peratom_warn = true;

  dim = 3;
  nablockmax = 0;
  nij = 0;
  nijmax = 0;
  szd = 0;
}

/* ---------------------------------------------------------------------- */

PairPOD::~PairPOD()
{
  free_tempmemory();
  memory->destroy(podcoeff);
  memory->destroy(newpodcoeff);
  memory->destroy(gd);
  memory->destroy(gdall);
  memory->destroy(energycoeff);
  memory->destroy(forcecoeff);

  delete podptr;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

void PairPOD::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  // we must enforce using F dot r, since we have no energy or stress tally calls.
  vflag_fdotr = 1;

  if (peratom_warn && (vflag_atom || eflag_atom)) {
    peratom_warn = false;
    if (comm->me == 0)
      error->warning(FLERR, "Pair style pod does not support per-atom energies or stresses");
  }

  double **x = atom->x;
  double **f = atom->f;
  int **firstneigh = list->firstneigh;
  int *numneigh = list->numneigh;
  int *type = atom->type;
  int *ilist = list->ilist;
  int inum = list->inum;

  // initialize global descriptors to zero

  int nd1234 = podptr->pod.nd1234;
  podptr->podArraySetValue(gd, 0.0, nd1234);

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    int jnum = numneigh[i];

    // allocate temporary memory

    if (nijmax < jnum) {
      nijmax = MAX(nijmax, jnum);
      nablockmax = 1;
      free_tempmemory();
      estimate_tempmemory();
      allocate_tempmemory();
    }

    // get neighbor pairs for atom i

    lammpsNeighPairs(x, firstneigh, type, map, numneigh, i);

    // compute global POD descriptors for atom i

    podptr->linear_descriptors_ij(gd, tmpmem, rij, &tmpmem[nd1234], numneighsum, typeai, idxi, ti,
                                  tj, 1, nij);
  }

  int nd22 = podptr->pod.nd22;
  int nd23 = podptr->pod.nd23;
  int nd24 = podptr->pod.nd24;
  int nd33 = podptr->pod.nd33;
  int nd34 = podptr->pod.nd34;
  int nd44 = podptr->pod.nd44;
  int nd = podptr->pod.nd;
  bigint natom = atom->natoms;

  for (int j = nd1234; j < (nd1234 + nd22 + nd23 + nd24 + nd33 + nd34 + nd44); j++)
    newpodcoeff[j] = podcoeff[j] / (natom);

  for (int j = (nd1234 + nd22 + nd23 + nd24 + nd33 + nd34 + nd44); j < nd; j++)
    newpodcoeff[j] = podcoeff[j] / (natom * natom);

  // compute energy and effective coefficients
  eng_vdwl = podptr->calculate_energy(energycoeff, forcecoeff, gd, gdall, newpodcoeff);

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];

    // get neighbor pairs for atom i

    lammpsNeighPairs(x, firstneigh, type, map, numneigh, i);

    // compute atomic force for atom i

    podptr->calculate_force(f, forcecoeff, rij, tmpmem, numneighsum, typeai, idxi, ai, aj, ti, tj,
                            1, nij);
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPOD::settings(int narg, char ** /* arg */)
{
  if (narg > 0) error->all(FLERR, "Pair style pod accepts no arguments");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPOD::coeff(int narg, char **arg)
{
  const int np1 = atom->ntypes + 1;
  memory->destroy(setflag);
  memory->destroy(cutsq);
  memory->create(setflag, np1, np1, "pair:setflag");
  memory->create(cutsq, np1, np1, "pair:cutsq");
  delete[] map;
  map = new int[np1];
  allocated = 1;

  if (narg < 4) utils::missing_cmd_args(FLERR, "pair_coeff", error);
  map_element2type(narg - 4, arg + 4);

  std::string pod_file = std::string(arg[2]);      // pod input file
  std::string coeff_file = std::string(arg[3]);    // coefficient input file

  delete podptr;
  podptr = new MLPOD(lmp, pod_file, coeff_file);

  if (coeff_file != "") {
    memory->destroy(podcoeff);
    memory->destroy(newpodcoeff);
    memory->destroy(energycoeff);
    memory->destroy(forcecoeff);
    memory->destroy(gd);
    memory->destroy(gdall);
    memory->create(podcoeff, podptr->pod.nd, "pair:podcoeff");
    memory->create(newpodcoeff, podptr->pod.nd, "pair:newpodcoeff");
    memory->create(energycoeff, podptr->pod.nd1234, "pair:energycoeff");
    memory->create(forcecoeff, podptr->pod.nd1234, "pair:forcecoeff");
    memory->create(gd, podptr->pod.nd1234, "pair:gd");
    memory->create(gdall, podptr->pod.nd1234, "pair:gdall");
    podptr->podArrayCopy(podcoeff, podptr->pod.coeff, podptr->pod.nd);
    podptr->podArrayCopy(newpodcoeff, podptr->pod.coeff, podptr->pod.nd);
  }

  for (int ii = 0; ii < np1; ii++)
    for (int jj = 0; jj < np1; jj++) cutsq[ii][jj] = podptr->pod.rcut * podptr->pod.rcut;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPOD::init_style()
{
  if (force->newton_pair == 0) error->all(FLERR, "Pair style pod requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);

  // reset flag to print warning about per-atom energies or stresses
  peratom_warn = true;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPOD::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
  return podptr->pod.rcut;
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double PairPOD::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes;
}

void PairPOD::free_tempmemory()
{
  memory->destroy(rij);
  memory->destroy(idxi);
  memory->destroy(ai);
  memory->destroy(aj);
  memory->destroy(ti);
  memory->destroy(tj);
  memory->destroy(numneighsum);
  memory->destroy(typeai);
  memory->destroy(tmpmem);
}

void PairPOD::allocate_tempmemory()
{
  memory->create(rij, dim * nijmax, "pair:rij");
  memory->create(idxi, nijmax, "pair:idxi");
  memory->create(ai, nijmax, "pair:ai");
  memory->create(aj, nijmax, "pair:aj");
  memory->create(ti, nijmax, "pair:ti");
  memory->create(tj, nijmax, "pair:tj");
  memory->create(numneighsum, nablockmax + 1, "pair:numneighsum");
  memory->create(typeai, nablockmax, "pair:typeai");
  memory->create(tmpmem, szd, "pair:tmpmem");
}

void PairPOD::estimate_tempmemory()
{
  int nrbf2 = podptr->pod.nbf2;
  int nabf3 = podptr->pod.nabf3;
  int nrbf3 = podptr->pod.nrbf3;
  int ns2 = podptr->pod.ns2;
  int ns3 = podptr->pod.ns3;

  szd = dim * nijmax + (1 + dim) * nijmax * MAX(nrbf2 + ns2, nrbf3 + ns3) + (nabf3 + 1) * 7;
  int szsnap = 0;
  if (podptr->sna.twojmax > 0) {
    szsnap += nijmax * dim;
    szsnap += MAX(2 * podptr->sna.idxu_max * nijmax,
                  2 * podptr->sna.idxz_max * podptr->sna.ndoubles *
                      nablockmax);                        // (Ur, Ui) and (Zr, Zi)
    szsnap += 2 * podptr->sna.idxu_max * dim * nijmax;    // dUr, dUi
    szsnap += MAX(podptr->sna.idxb_max * podptr->sna.ntriples * dim * nijmax,
                  2 * podptr->sna.idxu_max * podptr->sna.nelements *
                      nablockmax);    // dblist and (Utotr, Utoti)
  }

  szd = MAX(szsnap, szd);
  szd = nablockmax * (podptr->pod.nd1234) + szd;
}

void PairPOD::lammpsNeighPairs(double **x, int **firstneigh, int *atomtypes, int *map,
                               int *numneigh, int gi)
{

  double rcutsq = podptr->pod.rcut * podptr->pod.rcut;

  nij = 0;
  int itype = map[atomtypes[gi]] + 1;
  int m = numneigh[gi];
  typeai[0] = itype;
  for (int l = 0; l < m; l++) {           // loop over each atom around atom i
    int gj = firstneigh[gi][l];           // atom j
    double delx = x[gj][0] - x[gi][0];    // xj - xi
    double dely = x[gj][1] - x[gi][1];    // xj - xi
    double delz = x[gj][2] - x[gi][2];    // xj - xi
    double rsq = delx * delx + dely * dely + delz * delz;
    if (rsq < rcutsq && rsq > 1e-20) {
      rij[nij * 3 + 0] = delx;
      rij[nij * 3 + 1] = dely;
      rij[nij * 3 + 2] = delz;
      idxi[nij] = 0;
      ai[nij] = gi;
      aj[nij] = gj;
      ti[nij] = itype;
      tj[nij] = map[atomtypes[gj]] + 1;
      nij++;
    }
  }

  numneighsum[0] = 0;
  numneighsum[1] = nij;
}
