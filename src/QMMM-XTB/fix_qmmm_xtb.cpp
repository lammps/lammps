/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This file is distributed under the GNU General Public License.

   GFN-xTB electrostatic-embedding QM/MM with a self-consistent PPPM
   partition.  The implementation supports a compact, non-covalently
   embedded QM region in a fixed 3-D periodic cell.
------------------------------------------------------------------------- */

#include "fix_qmmm_xtb.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "pair_hybrid_overlay.h"
#include "pppm_tip4p_xtb.h"
#include "pppm_xtb.h"
#include "qmmm_xtb_adapter.h"
#include "qmmm_xtb_ewald.h"
#include "update.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <memory>
#include <numeric>
#include <string>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

namespace {

constexpr double BOHR_TO_ANGSTROM = 0.529177210903;
const char *ELEMENT_SYMBOLS[] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",
    "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re",
    "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
    "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

int atomic_number(const char *symbol)
{
  // A zero atomic number represents an extra point-charge site.  It uses
  // hydrogen's hardness when element-dependent MM hardness is requested,
  // but it is not permitted in the QM group.
  if (std::strcmp(symbol, "X") == 0 || std::strcmp(symbol, "NULL") == 0) return 0;
  for (int i = 0; i < static_cast<int>(sizeof(ELEMENT_SYMBOLS) / sizeof(const char *)); ++i)
    if (std::strcmp(symbol, ELEMENT_SYMBOLS[i]) == 0) return i + 1;
  return -1;
}

}    // namespace

FixQMMMXTB::FixQMMMXTB(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), pair_long(nullptr), pair_coulomb_mm_only(false),
    tip4p_qm_group_validated(false), xtb_method(2), qm_charge(0), qm_uhf(0), maxiter(250),
    cutoff(-1.0), accuracy(1.0e-3), electronic_temperature(300.0), mm_hardness(0.0),
    image_alpha(-1.0), image_kmax({8, 8, 8}), image_ksqmax(64), nqm(0), pair_mm_energy(0.0),
    pair_full_energy(0.0), mm_kspace_energy(0.0), qm_energy(0.0), energy_correction(0.0),
    pppm_xtb(nullptr), pppm_tip4p_xtb(nullptr), image_ewald(new QMMMXTBEwald), adapter_active(false)
{
  if (!atom->tag_enable) error->all(FLERR, "Fix qmmm/xtb requires atom IDs");
  if (atom->map_style == Atom::MAP_NONE) error->all(FLERR, "Fix qmmm/xtb requires an atom map");

  elements.assign(atom->ntypes + 1, 0);
  bool have_elements = false;
  bool have_image_alpha = false;

  int iarg = 3;
  while (iarg < narg) {
    if (std::strcmp(arg[iarg], "elements") == 0) {
      if (iarg + atom->ntypes >= narg)
        utils::missing_cmd_args(FLERR, "fix qmmm/xtb elements", error);
      for (int itype = 1; itype <= atom->ntypes; ++itype) {
        elements[itype] = atomic_number(arg[iarg + itype]);
        if (elements[itype] < 0)
          error->all(FLERR, "Invalid chemical element {} in fix qmmm/xtb", arg[iarg + itype]);
      }
      have_elements = true;
      iarg += atom->ntypes + 1;
    } else if (std::strcmp(arg[iarg], "charge") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix qmmm/xtb charge", error);
      qm_charge = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "uhf") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix qmmm/xtb uhf", error);
      qm_uhf = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "method") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix qmmm/xtb method", error);
      if (std::strcmp(arg[iarg + 1], "gfn1") == 0 || std::strcmp(arg[iarg + 1], "gfn1-xtb") == 0)
        xtb_method = 1;
      else if (std::strcmp(arg[iarg + 1], "gfn2") == 0 ||
               std::strcmp(arg[iarg + 1], "gfn2-xtb") == 0)
        xtb_method = 2;
      else
        error->all(FLERR, "Unknown fix qmmm/xtb method {}", arg[iarg + 1]);
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix qmmm/xtb cutoff", error);
      cutoff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "accuracy") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix qmmm/xtb accuracy", error);
      accuracy = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "maxiter") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix qmmm/xtb maxiter", error);
      maxiter = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "etemp") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix qmmm/xtb etemp", error);
      electronic_temperature = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "mmhardness") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix qmmm/xtb mmhardness", error);
      mm_hardness = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "ewald_alpha") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix qmmm/xtb ewald_alpha", error);
      image_alpha = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      have_image_alpha = true;
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "kmax") == 0) {
      if (iarg + 3 >= narg) utils::missing_cmd_args(FLERR, "fix qmmm/xtb kmax", error);
      for (int dim = 0; dim < 3; ++dim)
        image_kmax[dim] = utils::inumeric(FLERR, arg[iarg + 1 + dim], false, lmp);
      iarg += 4;
    } else if (std::strcmp(arg[iarg], "ksqmax") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix qmmm/xtb ksqmax", error);
      image_ksqmax = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown fix qmmm/xtb keyword {}", arg[iarg]);
    }
  }

  if (!have_elements) error->all(FLERR, "Fix qmmm/xtb requires the elements keyword");
  if (cutoff <= 0.0 || accuracy <= 0.0 || maxiter <= 0 || qm_uhf < 0 ||
      electronic_temperature < 0.0 || (have_image_alpha && image_alpha <= 0.0) ||
      image_ksqmax <= 0 ||
      std::any_of(
          image_kmax.begin(), image_kmax.end(),
          [](int value) {
            return value < 0;
          }) ||
      std::all_of(
          image_kmax.begin(), image_kmax.end(),
          [](int value) {
            return value == 0;
          }))
    error->all(FLERR, "Invalid numeric parameter in fix qmmm/xtb command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  virial_global_flag = 1;
  thermo_energy = thermo_virial = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  comm_forward = 1;

  pair_mm_virial.fill(0.0);
  pair_full_virial.fill(0.0);
  qmqm_kspace_virial.fill(0.0);
}

FixQMMMXTB::~FixQMMMXTB()
{
  if (comm->me == 0 && adapter_active) lammps_qmmm_xtb_destroy();
}

FixQMMMXTB::PairCoulombMapping FixQMMMXTB::classify_pair_coulomb_mapping(char *substyle)
{
  auto *hybrid = dynamic_cast<PairHybridOverlay *>(force->pair);
  if (!hybrid) return PairCoulombMapping::FULL;

  std::vector<int> local_qm_type(atom->ntypes + 1, 0);
  std::vector<int> local_mm_type(atom->ntypes + 1, 0);
  for (int i = 0; i < atom->nlocal; ++i) {
    if (atom->mask[i] & groupbit)
      local_qm_type[atom->type[i]] = 1;
    else
      local_mm_type[atom->type[i]] = 1;
  }

  std::vector<int> qm_type(atom->ntypes + 1, 0);
  std::vector<int> mm_type(atom->ntypes + 1, 0);
  MPI_Allreduce(local_qm_type.data(), qm_type.data(), atom->ntypes + 1, MPI_INT, MPI_MAX, world);
  MPI_Allreduce(local_mm_type.data(), mm_type.data(), atom->ntypes + 1, MPI_INT, MPI_MAX, world);

  bool full_mapping = true;
  bool mm_only_mapping = true;
  for (int itype = 1; itype <= atom->ntypes; ++itype) {
    if (qm_type[itype] && mm_type[itype]) mm_only_mapping = false;
    if (!qm_type[itype] && !mm_type[itype]) continue;
    for (int jtype = itype; jtype <= atom->ntypes; ++jtype) {
      if (!qm_type[jtype] && !mm_type[jtype]) continue;
      const bool mapped = hybrid->check_ijtype(itype, jtype, substyle);
      full_mapping &= mapped;
      const bool should_map_mm_only = mm_type[itype] && mm_type[jtype];
      mm_only_mapping &= mapped == should_map_mm_only;
    }
  }

  if (full_mapping) return PairCoulombMapping::FULL;

  // Type-pair routing can replace the reference subtraction only for a
  // Coulomb-only sub-style.  Combined LJ/Coulomb styles would also change
  // the retained Lennard-Jones interactions when their mappings are pruned.
  const bool coulomb_only_style = std::strcmp(substyle, "coul/long") == 0 ||
      std::strcmp(substyle, "coul/long/omp") == 0 || std::strcmp(substyle, "tip4p/long") == 0 ||
      std::strcmp(substyle, "tip4p/long/omp") == 0;
  if (coulomb_only_style && mm_only_mapping) return PairCoulombMapping::MM_ONLY;
  return PairCoulombMapping::INVALID;
}

int FixQMMMXTB::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

void FixQMMMXTB::init()
{
  pair_coulomb_mm_only = false;
  tip4p_qm_group_validated = false;

  if (modify->get_fix_by_style("^qmmm/xtb$").size() > 1)
    error->all(FLERR, "Only one instance of fix qmmm/xtb is supported");
  if (update->integrate_style && utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, "Fix qmmm/xtb does not yet support r-RESPA");
  if (!modify->get_fix_by_style("(npt|nph|press/|deform)").empty() ||
      !modify->get_fix_by_style("^box/relax").empty())
    error->all(FLERR, "Fix qmmm/xtb does not yet support a changing simulation box");
  if (update->whichflag == 2 && !thermo_energy)
    error->all(FLERR, "Fix qmmm/xtb requires fix_modify energy yes during minimization");
  if (!atom->q_flag) error->all(FLERR, "Fix qmmm/xtb requires per-atom charge");
  if (domain->dimension != 3 || !domain->xperiodic || !domain->yperiodic || !domain->zperiodic)
    error->all(FLERR, "Fix qmmm/xtb requires 3-D periodic boundaries");
  if (std::strcmp(update->unit_style, "real") != 0 && std::strcmp(update->unit_style, "metal") != 0)
    error->all(FLERR, "Fix qmmm/xtb currently supports real and metal units");
  if (std::fabs(force->dielectric - 1.0) > 1.0e-12)
    error->all(FLERR, "Fix qmmm/xtb requires dielectric 1.0");
  if (!force->kspace || !force->kspace->pppmflag || !force->kspace->xtbflag ||
      force->kspace->dispersionflag)
    error->all(FLERR, "Fix qmmm/xtb requires kspace style pppm/xtb or pppm/tip4p/xtb");
  if (force->kspace->tip4pflag) {
    pppm_tip4p_xtb = dynamic_cast<PPPMTIP4PXTB *>(force->kspace);
    if (!pppm_tip4p_xtb)
      error->all(FLERR, "Fix qmmm/xtb found an incompatible TIP4P xTB KSpace implementation");
  } else {
    pppm_xtb = dynamic_cast<PPPMXTB *>(force->kspace);
    if (!pppm_xtb)
      error->all(FLERR, "Fix qmmm/xtb found an incompatible xTB KSpace implementation");
  }

  if (force->kspace->tip4pflag) {
    // Accept the standard combined TIP4P/LJ styles as well as the Coulomb-only
    // style.  Non-Coulomb terms are identical in the two reference captures
    // and cancel from the QM/MM correction.
    const char *tip4p_styles[] = {"tip4p/long",
                                  "tip4p/long/omp",
                                  "lj/cut/tip4p/long",
                                  "lj/cut/tip4p/long/omp",
                                  "lj/cut/tip4p/long/opt",
                                  "lj/cut/tip4p/long/gpu"};
    for (const char *style : tip4p_styles) {
      pair_long = force->pair_match(style, 1, 0);
      if (pair_long) break;
    }
  } else {
    pair_long = force->pair_match("coul/long", 1, 0);
  }
  if (!pair_long)
    error->all(FLERR, "Fix qmmm/xtb requires a compatible long-range Coulomb pair style");
  // PairHybrid::check_ijtype() uses the sub-style keyword as its public lookup
  // interface.  Obtain the canonical keyword from the matched Pair instance.
  char *pair_long_style = force->pair_match_ptr(pair_long);
  if (!pair_long_style)
    error->all(FLERR, "Fix qmmm/xtb could not identify the long-range Coulomb pair style");
  int cut_dim = 0;
  auto *coulomb_cutoff = static_cast<double *>(pair_long->extract("cut_coul", cut_dim));
  if (!coulomb_cutoff || std::fabs(*coulomb_cutoff - cutoff) > 1.0e-8 * std::max(1.0, cutoff))
    error->all(FLERR, "Fix qmmm/xtb cutoff must equal the Coulomb pair-style cutoff");

  gather_qm_atoms(true);
  if (nqm == 0) error->all(FLERR, "Fix qmmm/xtb QM group is empty");

  const PairCoulombMapping pair_mapping = classify_pair_coulomb_mapping(pair_long_style);
  if (pair_mapping == PairCoulombMapping::INVALID)
    error->all(FLERR,
               "Fix qmmm/xtb Coulomb pair sub-style must cover all populated type pairs or "
               "exactly the MM-MM type pairs");
  pair_coulomb_mm_only = pair_mapping == PairCoulombMapping::MM_ONLY;
  if (pair_coulomb_mm_only && comm->me == 0)
    utils::logmesg(lmp,
                   "Fix qmmm/xtb detected an MM-only {} type-pair mapping; skipping the "
                   "MM/full pair reference evaluations\n",
                   pair_long_style);

  int status = 0;
  if (comm->me == 0) {
    std::vector<double> qm_bohr(3 * nqm);
    for (int i = 0; i < 3 * nqm; ++i) qm_bohr[i] = qm_x[i] / BOHR_TO_ANGSTROM;
    status = lammps_qmmm_xtb_create(nqm, qm_atomic_numbers.data(), qm_bohr.data(), xtb_method,
                                    qm_charge, qm_uhf, accuracy, maxiter, electronic_temperature);
  }
  MPI_Bcast(&status, 1, MPI_INT, 0, world);
  if (status) error->all(FLERR, "Could not initialize the xTB QM/MM adapter");
  adapter_active = true;
}

void FixQMMMXTB::setup_pre_force(int vflag)
{
  // Verlet and Min invoke setup_pre_force() before the first KSpace::setup().
  // The QM/MM reference solve needs volume-dependent PPPM data (notably
  // volume, influence functions, and the Ewald parameter), so initialize it
  // here.  The immediately following regular setup call is harmless and keeps
  // the normal KSpace lifecycle intact.
  force->kspace->setup();
  // Fix::init() runs before the first border communication.  Defer this check
  // until TIP4P's enlarged communication cutoff has acquired bonded H ghosts.
  if (!tip4p_qm_group_validated) {
    validate_tip4p_qm_group();
    tip4p_qm_group_validated = true;
  }
  pre_force(vflag);
}

void FixQMMMXTB::setup(int vflag)
{ post_force(vflag); }

void FixQMMMXTB::min_setup(int vflag)
{
  // Min::setup() performs the initial pre-force QM/MM solve through
  // setup_pre_force(), followed by the regular pair and KSpace evaluations.
  // Complete that first energy/force evaluation with the same correction used
  // by Verlet before the minimizer records its initial objective and gradient.
  post_force(vflag);
}

void FixQMMMXTB::min_pre_force(int vflag)
{
  // Every line-search trial needs a fresh self-consistent QM/MM state at its
  // trial coordinates before LAMMPS computes the production pair/KSpace terms.
  pre_force(vflag);
}

void FixQMMMXTB::min_post_force(int vflag)
{
  // Assemble the energy correction and replace double-counted classical
  // forces after the minimizer's production pair/KSpace evaluation.
  post_force(vflag);
}

void FixQMMMXTB::save_entry_forces()
{
  const int nall = atom->nlocal + atom->nghost;
  entry_force.resize(static_cast<std::size_t>(3) * nall);
  for (int i = 0; i < nall; ++i)
    for (int dim = 0; dim < 3; ++dim) entry_force[3 * i + dim] = atom->f[i][dim];
}

void FixQMMMXTB::restore_entry_forces()
{
  const int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; ++i)
    for (int dim = 0; dim < 3; ++dim) atom->f[i][dim] = entry_force[3 * i + dim];
}

void FixQMMMXTB::clear_forces()
{
  const int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; ++i) atom->f[i][0] = atom->f[i][1] = atom->f[i][2] = 0.0;
}

void FixQMMMXTB::set_qm_charges(double value)
{
  for (int i = 0; i < atom->nlocal; ++i)
    if (atom->mask[i] & groupbit) atom->q[i] = value;
  comm->forward_comm(this);
  force->kspace->qsum_qsq(0);
}

void FixQMMMXTB::restore_qm_charges()
{
  for (int i = 0; i < atom->nlocal; ++i) {
    if (!(atom->mask[i] & groupbit)) continue;
    const auto iter = qm_index.find(atom->tag[i]);
    if (iter == qm_index.end()) error->one(FLERR, "Fix qmmm/xtb lost a QM atom ID");
    atom->q[i] = qm_charge_scf[iter->second];
  }
  comm->forward_comm(this);
  force->kspace->qsum_qsq(0);
}

int FixQMMMXTB::pack_forward_comm(int n, int *list, double *buffer, int, int *)
{
  // Charges are not part of AtomVec's regular forward communication.  The
  // SCC Mulliken charges must nevertheless be present on ghost atoms before
  // the real-space Coulomb reference calculations are repeated.
  for (int i = 0; i < n; ++i) buffer[i] = atom->q[list[i]];
  return n;
}

void FixQMMMXTB::unpack_forward_comm(int n, int first, double *buffer)
{
  for (int i = 0; i < n; ++i) atom->q[first + i] = buffer[i];
}

int FixQMMMXTB::get_charge_site(int i, double *site, int *indices, double *weights)
{
  if (pppm_tip4p_xtb) return pppm_tip4p_xtb->get_charge_site(i, site, indices, weights);
  return pppm_xtb->get_charge_site(i, site, indices, weights);
}

void FixQMMMXTB::compute_group_potential(double *potential, int sensor_groupbit,
                                         int source_groupbit, bool invert_source)
{
  if (pppm_tip4p_xtb) {
    pppm_tip4p_xtb->compute_group_potential(potential, sensor_groupbit, source_groupbit,
                                            invert_source);
  } else {
    pppm_xtb->compute_group_potential(potential, sensor_groupbit, source_groupbit, invert_source);
  }
}

void FixQMMMXTB::gather_qm_atoms(bool initialize)
{
  int nlocal_qm = 0;
  for (int i = 0; i < atom->nlocal; ++i)
    if (atom->mask[i] & groupbit) ++nlocal_qm;

  std::vector<int> counts(comm->nprocs), displs(comm->nprocs, 0);
  MPI_Allgather(&nlocal_qm, 1, MPI_INT, counts.data(), 1, MPI_INT, world);
  for (int iproc = 1; iproc < comm->nprocs; ++iproc)
    displs[iproc] = displs[iproc - 1] + counts[iproc - 1];
  const int total = std::accumulate(counts.begin(), counts.end(), 0);

  std::vector<tagint> local_tags(nlocal_qm), tags(total);
  std::vector<int> local_types(nlocal_qm), types(total);
  std::vector<double> local_x(3 * nlocal_qm), local_xu(3 * nlocal_qm);
  std::vector<double> all_x(3 * total), all_xu(3 * total);
  int offset = 0;
  for (int i = 0; i < atom->nlocal; ++i) {
    if (!(atom->mask[i] & groupbit)) continue;
    local_tags[offset] = atom->tag[i];
    local_types[offset] = atom->type[i];
    double xu[3];
    domain->unmap(atom->x[i], atom->image[i], xu);
    for (int dim = 0; dim < 3; ++dim) {
      local_x[3 * offset + dim] = atom->x[i][dim];
      local_xu[3 * offset + dim] = xu[dim];
    }
    ++offset;
  }

  std::vector<int> counts3(comm->nprocs), displs3(comm->nprocs);
  for (int iproc = 0; iproc < comm->nprocs; ++iproc) {
    counts3[iproc] = 3 * counts[iproc];
    displs3[iproc] = 3 * displs[iproc];
  }
  MPI_Allgatherv(local_tags.data(), nlocal_qm, MPI_LMP_TAGINT, tags.data(), counts.data(),
                 displs.data(), MPI_LMP_TAGINT, world);
  MPI_Allgatherv(local_types.data(), nlocal_qm, MPI_INT, types.data(), counts.data(), displs.data(),
                 MPI_INT, world);
  MPI_Allgatherv(local_x.data(), 3 * nlocal_qm, MPI_DOUBLE, all_x.data(), counts3.data(),
                 displs3.data(), MPI_DOUBLE, world);
  MPI_Allgatherv(local_xu.data(), 3 * nlocal_qm, MPI_DOUBLE, all_xu.data(), counts3.data(),
                 displs3.data(), MPI_DOUBLE, world);

  std::vector<int> order(total);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return tags[i] < tags[j];
  });

  if (initialize) {
    nqm = total;
    qm_tags.resize(nqm);
    qm_atomic_numbers.resize(nqm);
    qm_index.clear();
    for (int iqm = 0; iqm < nqm; ++iqm) {
      const int source = order[iqm];
      qm_tags[iqm] = tags[source];
      qm_atomic_numbers[iqm] = elements[types[source]];
      if (qm_atomic_numbers[iqm] == 0)
        error->all(FLERR, "Fix qmmm/xtb does not allow an extra point site in the QM group");
      qm_index[qm_tags[iqm]] = iqm;
    }
    qm_charge_scf.assign(nqm, 0.0);
  } else if (total != nqm) {
    error->all(FLERR, "Fix qmmm/xtb does not support a dynamic QM group");
  }

  qm_x.resize(3 * nqm);
  qm_x_wrapped.resize(3 * nqm);
  // Anchor the compact QM cluster to the unwrapped position of its lowest-ID
  // atom, then place every other QM atom in the nearest periodic image.  This
  // preserves continuous motion while also handling an initially split
  // molecule whose image flags were not prepared by the input data file.
  const int anchor_source = order[0];
  const double anchor_unwrapped[3] = {all_xu[3 * anchor_source], all_xu[3 * anchor_source + 1],
                                      all_xu[3 * anchor_source + 2]};
  const double anchor_wrapped[3] = {all_x[3 * anchor_source], all_x[3 * anchor_source + 1],
                                    all_x[3 * anchor_source + 2]};
  for (int iqm = 0; iqm < nqm; ++iqm) {
    const int source = order[iqm];
    if (tags[source] != qm_tags[iqm])
      error->all(FLERR, "Fix qmmm/xtb does not support changing QM atom IDs");
    double dx = all_x[3 * source] - anchor_wrapped[0];
    double dy = all_x[3 * source + 1] - anchor_wrapped[1];
    double dz = all_x[3 * source + 2] - anchor_wrapped[2];
    domain->minimum_image(FLERR, dx, dy, dz);
    qm_x[3 * iqm] = anchor_unwrapped[0] + dx;
    qm_x[3 * iqm + 1] = anchor_unwrapped[1] + dy;
    qm_x[3 * iqm + 2] = anchor_unwrapped[2] + dz;
    for (int dim = 0; dim < 3; ++dim) qm_x_wrapped[3 * iqm + dim] = all_x[3 * source + dim];
  }
}

void FixQMMMXTB::gather_mm_points()
{
  std::vector<MMPoint> local;
  const double cutsq = cutoff * cutoff;
  for (int i = 0; i < atom->nlocal; ++i) {
    if (atom->mask[i] & groupbit) continue;

    double site[3], force_weights[3];
    int force_indices[3];
    const int nforce = get_charge_site(i, site, force_indices, force_weights);
    if (nforce < 1 || nforce > 3)
      error->one(FLERR, "Invalid charge-site force mapping for fix qmmm/xtb");

    bool selected = false;
    for (int iqm = 0; iqm < nqm && !selected; ++iqm) {
      double dx = site[0] - qm_x_wrapped[3 * iqm];
      double dy = site[1] - qm_x_wrapped[3 * iqm + 1];
      double dz = site[2] - qm_x_wrapped[3 * iqm + 2];
      domain->minimum_image(FLERR, dx, dy, dz);
      selected = dx * dx + dy * dy + dz * dz < cutsq;
    }
    if (!selected) continue;

    MMPoint point;
    point.tag = atom->tag[i];
    // A virtual charge site has no chemical element of its own.  This also
    // selects the documented extra-site hardness for implicit TIP4P M sites.
    point.atomic_number = nforce == 1 ? elements[atom->type[i]] : 0;
    double dx = site[0] - qm_x_wrapped[0];
    double dy = site[1] - qm_x_wrapped[1];
    double dz = site[2] - qm_x_wrapped[2];
    domain->minimum_image(FLERR, dx, dy, dz);
    point.x[0] = qm_x[0] + dx;
    point.x[1] = qm_x[1] + dy;
    point.x[2] = qm_x[2] + dz;
    point.charge = atom->q[i];
    point.nforce = nforce;
    for (int iparent = 0; iparent < 3; ++iparent) {
      point.force_tags[iparent] = 0;
      point.force_weights[iparent] = 0.0;
    }
    for (int iparent = 0; iparent < nforce; ++iparent) {
      point.force_tags[iparent] = atom->tag[force_indices[iparent]];
      point.force_weights[iparent] = force_weights[iparent];
    }
    local.push_back(point);
  }

  int nlocal_mm = local.size();
  std::vector<int> counts(comm->nprocs), displs(comm->nprocs, 0);
  MPI_Allgather(&nlocal_mm, 1, MPI_INT, counts.data(), 1, MPI_INT, world);
  for (int iproc = 1; iproc < comm->nprocs; ++iproc)
    displs[iproc] = displs[iproc - 1] + counts[iproc - 1];
  const int total = std::accumulate(counts.begin(), counts.end(), 0);

  std::vector<tagint> local_tags(nlocal_mm), tags(total);
  std::vector<int> local_atomic_numbers(nlocal_mm), atomic_numbers(total);
  std::vector<int> local_nforce(nlocal_mm), nforce(total);
  std::vector<double> local_data(4 * nlocal_mm), data(4 * total);
  std::vector<tagint> local_force_tags(3 * nlocal_mm), force_tags(3 * total);
  std::vector<double> local_force_weights(3 * nlocal_mm), force_weights(3 * total);
  for (int i = 0; i < nlocal_mm; ++i) {
    local_tags[i] = local[i].tag;
    local_atomic_numbers[i] = local[i].atomic_number;
    local_nforce[i] = local[i].nforce;
    local_data[4 * i] = local[i].x[0];
    local_data[4 * i + 1] = local[i].x[1];
    local_data[4 * i + 2] = local[i].x[2];
    local_data[4 * i + 3] = local[i].charge;
    for (int iparent = 0; iparent < 3; ++iparent) {
      local_force_tags[3 * i + iparent] = local[i].force_tags[iparent];
      local_force_weights[3 * i + iparent] = local[i].force_weights[iparent];
    }
  }
  std::vector<int> counts3(comm->nprocs), displs3(comm->nprocs);
  std::vector<int> counts4(comm->nprocs), displs4(comm->nprocs);
  for (int iproc = 0; iproc < comm->nprocs; ++iproc) {
    counts3[iproc] = 3 * counts[iproc];
    displs3[iproc] = 3 * displs[iproc];
    counts4[iproc] = 4 * counts[iproc];
    displs4[iproc] = 4 * displs[iproc];
  }
  MPI_Allgatherv(local_tags.data(), nlocal_mm, MPI_LMP_TAGINT, tags.data(), counts.data(),
                 displs.data(), MPI_LMP_TAGINT, world);
  MPI_Allgatherv(local_atomic_numbers.data(), nlocal_mm, MPI_INT, atomic_numbers.data(),
                 counts.data(), displs.data(), MPI_INT, world);
  MPI_Allgatherv(local_nforce.data(), nlocal_mm, MPI_INT, nforce.data(), counts.data(),
                 displs.data(), MPI_INT, world);
  MPI_Allgatherv(local_data.data(), 4 * nlocal_mm, MPI_DOUBLE, data.data(), counts4.data(),
                 displs4.data(), MPI_DOUBLE, world);
  MPI_Allgatherv(local_force_tags.data(), 3 * nlocal_mm, MPI_LMP_TAGINT, force_tags.data(),
                 counts3.data(), displs3.data(), MPI_LMP_TAGINT, world);
  MPI_Allgatherv(local_force_weights.data(), 3 * nlocal_mm, MPI_DOUBLE, force_weights.data(),
                 counts3.data(), displs3.data(), MPI_DOUBLE, world);

  std::vector<int> order(total);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return tags[i] < tags[j];
  });
  mm_points.resize(total);
  for (int imm = 0; imm < total; ++imm) {
    const int source = order[imm];
    mm_points[imm].tag = tags[source];
    mm_points[imm].atomic_number = atomic_numbers[source];
    for (int dim = 0; dim < 3; ++dim) mm_points[imm].x[dim] = data[4 * source + dim];
    mm_points[imm].charge = data[4 * source + 3];
    mm_points[imm].nforce = nforce[source];
    for (int iparent = 0; iparent < 3; ++iparent) {
      mm_points[imm].force_tags[iparent] = force_tags[3 * source + iparent];
      mm_points[imm].force_weights[iparent] = force_weights[3 * source + iparent];
    }
  }
}

void FixQMMMXTB::validate_tip4p_qm_group()
{
  if (!force->kspace->tip4pflag) return;

  int local_invalid = 0;
  for (int i = 0; i < atom->nlocal && !local_invalid; ++i) {
    double site[3], weights[3];
    int indices[3];
    const int nforce = get_charge_site(i, site, indices, weights);
    if (nforce == 1) continue;
    for (int iparent = 0; iparent < nforce; ++iparent)
      if (atom->mask[indices[iparent]] & groupbit) local_invalid = 1;
  }

  int invalid = 0;
  MPI_Allreduce(&local_invalid, &invalid, 1, MPI_INT, MPI_MAX, world);
  if (invalid)
    error->all(FLERR, "Fix qmmm/xtb does not allow an implicit TIP4P molecule in the QM group");
}

void FixQMMMXTB::capture_pair(std::vector<double> &captured_force, double &captured_energy,
                              std::array<double, 6> &captured_virial)
{
  clear_forces();
  pair_long->compute(1, 1);
  if (force->newton_pair) comm->reverse_comm();
  captured_force.resize(static_cast<std::size_t>(3) * atom->nlocal);
  for (int i = 0; i < atom->nlocal; ++i)
    for (int dim = 0; dim < 3; ++dim) captured_force[3 * i + dim] = atom->f[i][dim];

  double local_energy = pair_long->eng_coul + pair_long->eng_vdwl;
  MPI_Allreduce(&local_energy, &captured_energy, 1, MPI_DOUBLE, MPI_SUM, world);
  for (int i = 0; i < 6; ++i) captured_virial[i] = pair_long->virial[i];
}

void FixQMMMXTB::capture_kspace(std::vector<double> &captured_force, double &captured_energy,
                                std::array<double, 6> &captured_virial)
{
  clear_forces();
  force->kspace->qsum_qsq(0);
  force->kspace->compute(1, 1);
  captured_force.resize(static_cast<std::size_t>(3) * atom->nlocal);
  for (int i = 0; i < atom->nlocal; ++i)
    for (int dim = 0; dim < 3; ++dim) captured_force[3 * i + dim] = atom->f[i][dim];
  captured_energy = force->kspace->energy;
  for (int i = 0; i < 6; ++i) captured_virial[i] = force->kspace->virial[i];
}

void FixQMMMXTB::compute_mm_shift(std::vector<double> &shift)
{
  std::vector<double> potential(atom->nmax, 0.0);
  compute_group_potential(potential.data(), groupbit, groupbit, true);

  shift.assign(nqm, 0.0);
  for (int i = 0; i < atom->nlocal; ++i) {
    if (!(atom->mask[i] & groupbit)) continue;
    shift[qm_index.at(atom->tag[i])] = potential[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, shift.data(), nqm, MPI_DOUBLE, MPI_SUM, world);

  double local_mm_charge = 0.0;
  for (int i = 0; i < atom->nlocal; ++i)
    if (!(atom->mask[i] & groupbit)) local_mm_charge += atom->q[i];
  double total_mm_charge = 0.0;
  MPI_Allreduce(&local_mm_charge, &total_mm_charge, 1, MPI_DOUBLE, MPI_SUM, world);

  const double alpha = force->kspace->g_ewald;
  const double volume = domain->xprd * domain->yprd * domain->zprd;
  const double background = -MY_PI * total_mm_charge / (alpha * alpha * volume);
  for (int iqm = 0; iqm < nqm; ++iqm) {
    shift[iqm] += background;
    for (const auto &point : mm_points) {
      const double dx = qm_x[3 * iqm] - point.x[0];
      const double dy = qm_x[3 * iqm + 1] - point.x[1];
      const double dz = qm_x[3 * iqm + 2] - point.x[2];
      const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      shift[iqm] -= point.charge * std::erf(alpha * r) / r;
    }
    shift[iqm] *= BOHR_TO_ANGSTROM;
  }
}

void FixQMMMXTB::run_xtb(const std::vector<double> &mm_shift,
                         const std::vector<double> &image_response)
{
  int status = 0;
  double energy_hartree = 0.0;
  qm_gradient.assign(3 * nqm, 0.0);
  qm_charge_scf.assign(nqm, 0.0);
  mm_gradient.assign(3 * mm_points.size(), 0.0);

  if (comm->me == 0) {
    std::vector<double> qm_bohr(3 * nqm);
    for (int i = 0; i < 3 * nqm; ++i) qm_bohr[i] = qm_x[i] / BOHR_TO_ANGSTROM;
    std::vector<double> point_bohr(std::max<std::size_t>(1, 3 * mm_points.size()), 0.0);
    std::vector<double> point_charge(std::max<std::size_t>(1, mm_points.size()), 0.0);
    std::vector<int> point_atomic_numbers(std::max<std::size_t>(1, mm_points.size()), 0);
    for (int imm = 0; imm < static_cast<int>(mm_points.size()); ++imm) {
      for (int dim = 0; dim < 3; ++dim)
        point_bohr[3 * imm + dim] = mm_points[imm].x[dim] / BOHR_TO_ANGSTROM;
      point_charge[imm] = mm_points[imm].charge;
      point_atomic_numbers[imm] = mm_points[imm].atomic_number;
    }
    std::vector<double> point_gradient(std::max<std::size_t>(1, 3 * mm_points.size()), 0.0);
    status = lammps_qmmm_xtb_calculate(
        nqm, qm_bohr.data(), mm_points.size(), point_bohr.data(), point_charge.data(),
        point_atomic_numbers.data(), mm_hardness, mm_shift.data(), image_response.data(),
        &energy_hartree, qm_gradient.data(), qm_charge_scf.data(), point_gradient.data());
    if (!mm_points.empty())
      std::copy(point_gradient.begin(), point_gradient.begin() + mm_gradient.size(),
                mm_gradient.begin());
  }

  MPI_Bcast(&status, 1, MPI_INT, 0, world);
  if (status) error->all(FLERR, "xTB SCC calculation failed in fix qmmm/xtb");
  MPI_Bcast(&energy_hartree, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(qm_gradient.data(), 3 * nqm, MPI_DOUBLE, 0, world);
  MPI_Bcast(qm_charge_scf.data(), nqm, MPI_DOUBLE, 0, world);
  if (!mm_gradient.empty()) MPI_Bcast(mm_gradient.data(), mm_gradient.size(), MPI_DOUBLE, 0, world);

  qm_energy = energy_hartree * force->qqrd2e / BOHR_TO_ANGSTROM;
}

void FixQMMMXTB::build_periodic_forces()
{
  const double xtb_force_to_lmp = force->qqrd2e / (BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM);
  qm_force_correction.resize(3 * nqm);
  for (int i = 0; i < 3 * nqm; ++i) qm_force_correction[i] = -qm_gradient[i] * xtb_force_to_lmp;
  mm_force_correction.resize(3 * mm_points.size());
  for (int i = 0; i < static_cast<int>(mm_force_correction.size()); ++i)
    mm_force_correction[i] = -mm_gradient[i] * xtb_force_to_lmp;

  std::vector<double> image_force(3 * nqm, 0.0);
  try {
    image_ewald->add_force(qm_x, qm_charge_scf, image_force);
  } catch (const std::exception &exception) {
    error->all(FLERR, "Fix qmmm/xtb direct-Ewald force failed: {}", exception.what());
  }
  for (int i = 0; i < 3 * nqm; ++i) qm_force_correction[i] += force->qqrd2e * image_force[i];

  const double alpha = force->kspace->g_ewald;
  try {
    for (int iqm = 0; iqm < nqm; ++iqm) {
      for (int imm = 0; imm < static_cast<int>(mm_points.size()); ++imm) {
        double fqm[3] = {0.0, 0.0, 0.0};
        double fmm[3] = {0.0, 0.0, 0.0};
        QMMMXTBEwald::add_erf_pair_force(&qm_x[3 * iqm], mm_points[imm].x, qm_charge_scf[iqm],
                                         mm_points[imm].charge, alpha, fqm, fmm);
        for (int dim = 0; dim < 3; ++dim) {
          qm_force_correction[3 * iqm + dim] += force->qqrd2e * fqm[dim];
          mm_force_correction[3 * imm + dim] += force->qqrd2e * fmm[dim];
        }
      }
    }
  } catch (const std::exception &exception) {
    error->all(FLERR, "Fix qmmm/xtb QM/MM Ewald force failed: {}", exception.what());
  }
}

void FixQMMMXTB::pre_force(int vflag)
{
  v_init(vflag);
  save_entry_forces();
  gather_qm_atoms(false);

  // First PPPM pass: QM charges are zero, so this is the MM-only reference
  // used both for the fixed SCC potential and for energy de-duplication.
  set_qm_charges(0.0);
  gather_mm_points();
  if (!pair_coulomb_mm_only) capture_pair(pair_mm_force, pair_mm_energy, pair_mm_virial);

  std::vector<double> mm_shift;
  compute_mm_shift(mm_shift);

  std::vector<double> unused_force;
  std::array<double, 6> unused_virial;
  capture_kspace(unused_force, mm_kspace_energy, unused_virial);

  const double alpha = image_alpha > 0.0
      ? image_alpha
      : 10.0 / std::cbrt(domain->xprd * domain->yprd * domain->zprd);
  std::vector<double> image_response;
  try {
    // LAMMPS stores a restricted triclinic cell as columns
    // a=(xprd,0,0), b=(xy,yprd,0), c=(xz,yz,zprd).  Passing the full matrix
    // lets the direct QM-image Ewald operator construct 2*pi*H^-T*n.
    const QMMMXTBEwald::CellMatrix cell = {domain->h[0],
                                           domain->triclinic ? domain->h[5] : 0.0,
                                           domain->triclinic ? domain->h[4] : 0.0,
                                           0.0,
                                           domain->h[1],
                                           domain->triclinic ? domain->h[3] : 0.0,
                                           0.0,
                                           0.0,
                                           domain->h[2]};
    image_ewald->setup(cell, alpha, image_kmax, image_ksqmax);
    image_ewald->response(qm_x, image_response);
  } catch (const std::exception &exception) {
    error->all(FLERR, "Fix qmmm/xtb direct-Ewald response failed: {}", exception.what());
  }
  for (double &value : image_response) value *= BOHR_TO_ANGSTROM;

  run_xtb(mm_shift, image_response);
  restore_qm_charges();
  if (!pair_coulomb_mm_only) capture_pair(pair_full_force, pair_full_energy, pair_full_virial);

  // The final production PPPM keeps the cross QM-MM reciprocal force, but the
  // QM-QM reciprocal force is already represented by the SCC image term.
  const int nall = atom->nlocal + atom->nghost;
  std::vector<double> saved_charge(nall);
  for (int i = 0; i < nall; ++i) saved_charge[i] = atom->q[i];
  for (int i = 0; i < atom->nlocal; ++i)
    if (!(atom->mask[i] & groupbit)) atom->q[i] = 0.0;
  comm->forward_comm(this);
  capture_kspace(qmqm_kspace_force, energy_correction, qmqm_kspace_virial);
  for (int i = 0; i < nall; ++i) atom->q[i] = saved_charge[i];
  comm->forward_comm(this);
  force->kspace->qsum_qsq(0);

  build_periodic_forces();
  restore_entry_forces();
}

void FixQMMMXTB::post_force(int vflag)
{
  energy_correction = qm_energy + mm_kspace_energy - force->kspace->energy;
  if (!pair_coulomb_mm_only) energy_correction += pair_mm_energy - pair_full_energy;

  std::vector<double> xtb_local(static_cast<std::size_t>(3) * atom->nlocal, 0.0);
  for (int iqm = 0; iqm < nqm; ++iqm) {
    const int ilocal = atom->map(qm_tags[iqm]);
    if (ilocal < 0 || ilocal >= atom->nlocal) continue;
    for (int dim = 0; dim < 3; ++dim)
      xtb_local[3 * ilocal + dim] += qm_force_correction[3 * iqm + dim];
  }
  for (int imm = 0; imm < static_cast<int>(mm_points.size()); ++imm) {
    // Forces returned at a virtual charge site are distributed with the same
    // linear weights used by the electrostatic model (O/H/H for TIP4P).
    for (int iparent = 0; iparent < mm_points[imm].nforce; ++iparent) {
      const int ilocal = atom->map(mm_points[imm].force_tags[iparent]);
      if (ilocal < 0 || ilocal >= atom->nlocal) continue;
      for (int dim = 0; dim < 3; ++dim)
        xtb_local[3 * ilocal + dim] +=
            mm_points[imm].force_weights[iparent] * mm_force_correction[3 * imm + dim];
    }
  }

  for (int i = 0; i < atom->nlocal; ++i) {
    for (int dim = 0; dim < 3; ++dim) {
      atom->f[i][dim] += -qmqm_kspace_force[3 * i + dim] + xtb_local[3 * i + dim];
      if (!pair_coulomb_mm_only)
        atom->f[i][dim] += pair_mm_force[3 * i + dim] - pair_full_force[3 * i + dim];
    }
  }

  if (vflag_global) {
    for (int i = 0; i < 6; ++i) virial[i] = -qmqm_kspace_virial[i] / comm->nprocs;
    if (!pair_coulomb_mm_only)
      for (int i = 0; i < 6; ++i) virial[i] += pair_mm_virial[i] - pair_full_virial[i];
    for (int i = 0; i < atom->nlocal; ++i) {
      virial[0] += xtb_local[3 * i] * atom->x[i][0];
      virial[1] += xtb_local[3 * i + 1] * atom->x[i][1];
      virial[2] += xtb_local[3 * i + 2] * atom->x[i][2];
      virial[3] += xtb_local[3 * i + 1] * atom->x[i][0];
      virial[4] += xtb_local[3 * i + 2] * atom->x[i][0];
      virial[5] += xtb_local[3 * i + 2] * atom->x[i][1];
    }
  }
}

double FixQMMMXTB::compute_scalar()
{ return energy_correction; }

double FixQMMMXTB::memory_usage()
{
  return sizeof(double) *
      (qm_x.capacity() + qm_x_wrapped.capacity() + qm_charge_scf.capacity() +
       qm_gradient.capacity() + qm_force_correction.capacity() + mm_gradient.capacity() +
       mm_force_correction.capacity() + entry_force.capacity() + pair_mm_force.capacity() +
       pair_full_force.capacity() + qmqm_kspace_force.capacity()) +
      sizeof(MMPoint) * mm_points.capacity();
}
