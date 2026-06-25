/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Knut Nikolas Lausch, Gunnar Schmitz,
                        Alexander Knoll
------------------------------------------------------------------------- */

#include "pair_runner.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "lmptype.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <vector>

// External RuNNer library interface
extern "C" {
int runner_lammps_api_version();
void runner_lammps_interface_init(const char *path, int *npath, double *cutoff, double *cfenergy,
                                  double *cflength, int *nnp_generation, int *num_committee_members,
                                  bool *l_hirshfeld_vdw, bool *ltwo_body, bool *lcheck_extrap,
                                  int *rank, int *size);

void runner_lammps_interface_transfer_atoms_and_neighbor_lists(
    int *nlocal, int *nghost, int *atomic_numbers, int *inum, int *sum_num_neigh, int *ilist,
    int *num_neigh, int *first_neigh, int *neigh, double *lattice, double *xyz, bool *lperiodic,
    int *lstress);

void runner_interface_short_range(int *nlocal, int *nghost, int *inum, int *nmax, int *ilist,
                                  double *energy, double *forces, double *d_energy_d_strain,
                                  double *hirsh_volumes, double *atomic_charges,
                                  double *elec_negativities, double *hardness);

void runner_interface_evaluate_electrostatics_3g_part_1(int *num_atoms, double *xyz,
                                                        double *total_charge, double *lattice,
                                                        bool *lperiodic, double *atomic_charges,
                                                        double *energy, double *force_outer,
                                                        double *d_energy_d_q,
                                                        double *d_energy_d_strain);

void runner_interface_reinitialize_electrostatics(int *num_atoms_global, double *pos_global,
                                                  int *atomic_numbers_global);

void runner_interface_calc_screening(int *nlocal, int *nghost, double *atomic_charges,
                                     double *energy, double *forces, double *de_dq,
                                     double *d_energy_d_strain);

void runner_interface_evaluate_electrostatics_3g_part_2(int *nlocal, int *nghost, int *nglobal,
                                                        int icomm, double *energy, double *forces,
                                                        double *d_energy_d_q,
                                                        double *d_energy_d_q_sum_global,
                                                        double *d_energy_d_strain);

void runner_interface_compute_charges_4g(int *num_atoms, double *total_charge,
                                         double *atomic_electronegativities,
                                         double *atomic_hardness, double *atomic_charges,
                                         bool *luse_prev_q, int icomm);

void runner_interface_finalize_step();

void runner_interface_short_range_4g(int *nlocal, int *nghost, int *inum, int *nmax, int *ilist,
                                     double *, double *energy, double *forces,
                                     double *d_energy_d_strain, double *d_energy_d_q);

void runner_interface_evaluate_electrostatics_4g_part_1(int *nglobal, double *d_energy_dq,
                                                        double *energy, double *forces,
                                                        double *d_energy_d_strain,
                                                        double *lagrange_charges, int icomm);

void runner_interface_evaluate_electrostatics_4g_part_2(int *nlocal, int *nghost, int icomm,
                                                        double *lagrange_charges, double *charges,
                                                        double *forces, double *d_energy_d_strain);

void runner_interface_hirshfeld_vdw(int *nlocal, int *nghost, int *inum, int *ilist, int icomm,
                                    double *hirsh_volumes, double *energy, double *forces,
                                    double *d_energy_d_strain);

void runner_interface_two_body(int *nlocal, int *nghost, double *energy, double *forces,
                               double *d_energy_d_strain);

void runner_interface_extrapolation_warnings(char **c_ptr_extrap_msg, int *len_extrap_msg,
                                             int *global_atom_ids, int *nlocal);

void runner_interface_dealloc_extrapolation_warnings();

void runner_interface_extrapolation_count(int64_t *extraplation_count,
                                          int64_t *total_extrapolation_count, bool *lreset);
}

using namespace LAMMPS_NS;

namespace {
enum {
  COMM_NONE = 0,
  COMM_AT_CHARGE = 1,
  COMM_ELEC_NEGATIVITY = 2,
  COMM_HIRSH_VOLUME = 3,
  COMM_DEDQ = 4,
  COMM_SCREENING_DEDQ = 5,
  COMM_LAMBDA_CHARGE = 6,
  COMM_COMMITTEE_STORAGE = 7
};
}

// Initialize the static instance counter
int PairRuNNer::instances = 0;

PairRuNNer::PairRuNNer(LAMMPS *lmp) :
    Pair(lmp), map(nullptr), atomic_charge(nullptr), hirshfeld_volume(nullptr),
    electronegativity(nullptr), lagrange_charges(nullptr), de_dq(nullptr), screening_de_dq(nullptr),
    committee_storage(nullptr)
{
  // Sanity check: Prevent multiple instances due to static Fortran interface
  if (instances > 0) { error->all(FLERR, "Only one pair runner instance can be active at a time"); }
  instances++;

  // HDNNP is not pairwise additive, due to three body terms
  single_enable = 0;
  // Do not write binary restart files
  restartinfo = 0;
  // Parameters are read from input.nn, therefore pair_coeff only has single
  // command.
  one_coeff = 1;
  // Many-body potential flag.
  manybody_flag = 1;
  // Currently no unit conversion necessary as RuNNer is unit-agnostic.
  unit_convert_flag = 0;
  // We generally calculate the virial ourselves
  // and do not need to call virial_fdotr().
  // In case of 2G HDNNP, we call virial_fdotr()
  // for performance reasons.
  // flag is then overwritten in init_style().
  no_virial_fdotr_compute = 1;

  nmax = 0;

  // Set commsize needed by this pairstyle.
  comm_forward = 1;    // Forward communication (1 double per atom)
  comm_reverse = 1;    // Reverse communication (1 double per atom)
  commstyle = COMM_NONE;

  // Sum of extrapolation is zero upon initialization of this pair style
  // Tracks number of extrapolations over multiple timesteps
  local_extrap_sum = 0;

  // variable for output using pair compute
  nextra = 0;
  pvector = nullptr;

  if (atom->natoms > MAXSMALLINT / 4) error->all(FLERR, "Too many total atoms");
}

PairRuNNer::~PairRuNNer()
{
  // Decrement instance counter
  instances--;

  // Deallocate member variables
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(atomic_charge);
    memory->destroy(hirshfeld_volume);
    memory->destroy(electronegativity);
    memory->destroy(lagrange_charges);
    memory->destroy(de_dq);
    memory->destroy(screening_de_dq);
    memory->destroy(committee_storage);
    delete[] map;
    delete[] directory;
    delete[] pvector;
  }
}

void PairRuNNer::compute(int eflag, int vflag)
{
  int inum, jnum, ii, jj, i, j;
  int *ilist;
  int *jlist;
  int *numneigh, **firstneigh;

  // Initializes flags, which signal if energy and virial need to be tallied.
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  // Number of atoms owned by this process.
  int nlocal = atom->nlocal;
  // Number of ghost atoms on this process. Ghost atoms are atoms owned by another
  // process which are inside the neighbor list cutoff for local atoms on this
  // process. We have to calculate the force contribution from local atoms to
  // these ghost atoms.
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  // Total number of atoms in the simulation box.
  int natoms = static_cast<int>(atom->natoms);
  int *type = atom->type;
  // Global atom-ids. Used for correct mapping of extrapolation warnings and
  // sorting of global properties.
  tagint *tag = atom->tag;
  std::vector<int> global_atom_ids(nlocal);
  for (i = 0; i < nlocal; i++) global_atom_ids[i] = static_cast<int>(tag[i]);

  // Interface variables
  bool lperiodic;
  double lattice[9] = {0.0};

  // MPI
  int size, rank;
  rank = comm->me;
  size = comm->nprocs;

  // Allocate additional per-atom arrays
  if (atom->nmax > nmax) {
    memory->destroy(atomic_charge);
    memory->destroy(hirshfeld_volume);
    memory->destroy(electronegativity);
    memory->destroy(lagrange_charges);
    memory->destroy(de_dq);
    memory->destroy(screening_de_dq);
    memory->destroy(committee_storage);
    nmax = atom->nmax;

    memory->create(atomic_charge, nmax, "pair:atomic_charge");
    memory->create(hirshfeld_volume, nmax, "pair:hirshfeld_volume");
    memory->create(electronegativity, nmax, "pair:electronegativity");
    memory->create(lagrange_charges, nmax, "pair:lagrange_charges");
    memory->create(de_dq, nmax, "pair:de_dq");
    memory->create(screening_de_dq, nmax, "pair:screening_de_dq");
    memory->create(committee_storage, nmax, "pair:committee_storage");
  }

  // Set additional per-atom arrays to zero
  memset(atomic_charge, 0, nmax * (sizeof *atomic_charge));
  memset(hirshfeld_volume, 0, nmax * (sizeof *hirshfeld_volume));
  memset(electronegativity, 0, nmax * (sizeof *electronegativity));
  memset(lagrange_charges, 0, nmax * (sizeof *lagrange_charges));
  memset(de_dq, 0, nmax * (sizeof *de_dq));
  memset(screening_de_dq, 0, nmax * (sizeof *screening_de_dq));
  memset(committee_storage, 0, nmax * (sizeof *committee_storage));

  std::vector<double> committee_energy(num_committee_members, 0.0);
  std::vector<double> committee_force(nall * 3 * num_committee_members, 0.0);
  std::vector<double> committee_d_energy_d_strain(9 * num_committee_members, 0.0);
  std::vector<double> committee_atomic_charge(nmax * num_committee_members, 0.0);
  std::vector<double> committee_hirshfeld_volume(nmax * num_committee_members, 0.0);
  std::vector<double> committee_electronegativity(nmax * num_committee_members, 0.0);
  std::vector<double> committee_hardness(nmax * num_committee_members, 0.0);

  // Neighborlist information
  // Number of local atoms for which a neighbor list has been built
  // This is usually equal to nlocal, but may also be smaller when some
  // atoms have been removed from the neighbor list calculation by the user.
  inum = list->inum;
  // Local index of those inum atoms for which a neighbor list was built.
  ilist = list->ilist;
  // The number of neighbors for each atom in ilist (dimension inum).
  numneigh = list->numneigh;
  // Pointer array to the first neighbor of each atom in ilist (dimension inum).
  firstneigh = list->firstneigh;

  // Calculate total number of neighbors
  int num_neigh_sum = 0;
  for (ii = 0; ii < inum; ii++) num_neigh_sum += numneigh[ii];

  // create arrays for interface
  std::vector<int> runner_num_neigh(inum);
  std::vector<int> runner_first_neighbor(inum);
  std::vector<int> runner_jlist(num_neigh_sum);
  std::vector<int> runner_types(nall);

  // Collect neighbor data
  int irunner = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];            // get local index of atom ii
    jlist = firstneigh[i];    // pointer to first neighbor of atom
    jnum = numneigh[i];       // number of neighbors jj of atom
    runner_num_neigh[ii] = jnum;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];     // index of neighbor jj
      j &= NEIGHMASK;    // masks bits, which encode pair information
      if (jj == 0) runner_first_neighbor[ii] = irunner + 1;    // plus one due to Fortran
      runner_jlist[irunner] = j;
      irunner++;    // runs till num_neigh_sum
    }
  }
  // collect atomic numbers of all atoms by converting element id to atomic
  // number using map
  for (ii = 0; ii < nall; ii++) runner_types[ii] = map[type[ii]];

  // get lattice parameters
  lattice[0] = domain->xprd;
  lattice[1] = 0.0;
  lattice[2] = 0.0;
  lattice[3] = domain->xy;
  lattice[4] = domain->yprd;
  lattice[5] = 0.0;
  lattice[6] = domain->xz;
  lattice[7] = domain->yz;
  lattice[8] = domain->zprd;

  // get periodicity information
  switch (domain->periodicity[0] + domain->periodicity[1] + domain->periodicity[2]) {
    case 3:
      lperiodic = true;
      break;
    case 0:
      lperiodic = false;
      break;
    default:
      if (nnp_generation == 2) {
        lperiodic = false;
        break;
      } else {
        error->all(
            FLERR,
            "Periodicity in only one or two dimensions not supported when using 3G or 4G HDNNPs,");
      }
  }

  if (lperiodic && total_charge != 0.0)
    error->all(FLERR,
               "Periodic systems must be charge neutral (total_charge = 0.0) when using pair_style "
               "runner.");

  runner_lammps_interface_transfer_atoms_and_neighbor_lists(
      &nlocal, &nghost, runner_types.data(), &inum, &num_neigh_sum, ilist, runner_num_neigh.data(),
      runner_first_neighbor.data(), runner_jlist.data(), lattice, &x[0][0], &lperiodic,
      &vflag_global);

  runner_interface_short_range(&nlocal, &nghost, &inum, &nmax, ilist, committee_energy.data(),
                               committee_force.data(), committee_d_energy_d_strain.data(),
                               committee_hirshfeld_volume.data(), committee_atomic_charge.data(),
                               committee_electronegativity.data(), committee_hardness.data());

  if (lhirshfeld_vdw) {

    for (int i = 0; i < num_committee_members; i++) {

      int icomm_fortran = i + 1;

      double vdw_energy = 0.0;
      std::vector<double> vdw_forces(nall * 3, 0.0);
      double vdw_d_energy_d_strain[9] = {0.0};

      for (int ii = 0; ii < nmax; ii++)
        hirshfeld_volume[ii] = committee_hirshfeld_volume[nmax * i + ii];

      // Communicate Hirshfeld volumes from local atoms to ghost atoms.
      commstyle = COMM_HIRSH_VOLUME;
      comm->forward_comm(this);

      // Calculate dispersion energies and forces using Hirshfeld volumes
      // and volume gradients (stored on runner side)
      runner_interface_hirshfeld_vdw(&nlocal, &nghost, &inum, ilist, icomm_fortran,
                                     hirshfeld_volume, &vdw_energy, vdw_forces.data(),
                                     vdw_d_energy_d_strain);

      // Add vdw interactions to short-range results
      committee_energy[i] += vdw_energy;

      for (ii = 0; ii < nall * 3; ii++) committee_force[nall * 3 * i + ii] += vdw_forces[ii];

      for (ii = 0; ii < 9; ii++)
        committee_d_energy_d_strain[i * 9 + ii] += vdw_d_energy_d_strain[ii];
    }
  }

  if (ltwo_body) {

    double two_body_energy = 0.0;
    std::vector<double> two_body_forces(nall * 3, 0.0);
    double two_body_d_energy_d_strain[9] = {0.0};

    // Calculate two-body energies and forces
    runner_interface_two_body(&nlocal, &nghost, &two_body_energy, two_body_forces.data(),
                              two_body_d_energy_d_strain);

    for (int i = 0; i < num_committee_members; i++) {
      // Add two-body interactions to short-range results
      committee_energy[i] += two_body_energy;

      for (ii = 0; ii < nall * 3; ii++) committee_force[nall * 3 * i + ii] += two_body_forces[ii];

      for (ii = 0; ii < 9; ii++)
        committee_d_energy_d_strain[i * 9 + ii] += two_body_d_energy_d_strain[ii];
    }
  }

  if (nnp_generation >= 3) {
    // Collect local atoms into a global structure. This needs to be done only
    // once for all committee members.
    std::vector<double> xyz_global(natoms * 3, 0.0);
    std::vector<int> z_global(natoms, 0);

    pack_structure(rank, size, natoms, inum, ilist, tag, x, runner_types.data(), xyz_global.data(),
                   z_global.data());

    if (rank == 0) {
      // Reinitialize the electrostatics calculator on the root process.
      // It is completely reallocated if the global
      // number of atoms changed. Otherwise, only the ordering
      // of atoms is updated.
      runner_interface_reinitialize_electrostatics(&natoms, &xyz_global[0], z_global.data());
    }

    if (nnp_generation == 3) {

      for (int i = 0; i < num_committee_members; i++) {

        int icomm_fortran = i + 1;

        // Collect the charge predicted by the committee members into one global array.
        std::vector<double> q_global(natoms, 0.0);

        pack_atomic_property(rank, size, natoms, inum, ilist, tag,
                             &committee_atomic_charge[nmax * i], q_global.data());

        // long-range electrostatics variables
        double runner_elec_energy = 0.0;
        std::vector<double> elec_force_global(natoms * 3, 0.0);
        std::vector<double> de_dq_global(natoms, 0.0);
        double runner_elec_d_energy_d_strain[9] = {0.0};

        if (rank == 0) {
          // Calculate long-range electrostatics on root using the global structure.
          runner_interface_evaluate_electrostatics_3g_part_1(
              &natoms, &xyz_global[0], &total_charge, lattice, &lperiodic, &q_global[0],
              &runner_elec_energy, elec_force_global.data(), de_dq_global.data(),
              runner_elec_d_energy_d_strain);
        }

        MPI_Barrier(world);

        // Broadcast and unpack electrostatic results back to each local process.
        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 1, q_global.data(),
                                       atomic_charge);
        for (ii = 0; ii < nmax; ii++) committee_atomic_charge[nmax * i + ii] = atomic_charge[ii];

        memset(de_dq, 0, nall * (sizeof *de_dq));
        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 1, de_dq_global.data(),
                                       de_dq);

        std::vector<double> runner_elec_forces(nall * 3, 0.0);
        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 3,
                                       elec_force_global.data(), runner_elec_forces.data());

        // Communicate scaled atomic charges from local atoms to ghost
        // atoms for screening calculation.
        commstyle = COMM_AT_CHARGE;
        comm->forward_comm(this);

        double screening_energy = 0.0;
        std::vector<double> screening_forces(nall * 3, 0.0);
        memset(screening_de_dq, 0, nall * (sizeof *screening_de_dq));
        double screening_d_energy_d_strain[9] = {0.0};

        // Apply screening
        runner_interface_calc_screening(&nlocal, &nghost, atomic_charge, &screening_energy,
                                        screening_forces.data(), screening_de_dq,
                                        screening_d_energy_d_strain);

        // Communicate screening de_dq from ghost atoms to local atoms
        commstyle = COMM_SCREENING_DEDQ;
        comm->reverse_comm(this);

        // Determine sum of de_dq of local atoms across all procs.
        double de_dq_sum_local = 0.0;
        for (j = 0; j < inum; j++) {
          ii = ilist[j];
          de_dq[ii] -= screening_de_dq[ii];
          de_dq_sum_local += de_dq[ii];
        }
        double de_dq_sum_global = 0.0;
        MPI_Allreduce(&de_dq_sum_local, &de_dq_sum_global, 1, MPI_DOUBLE, MPI_SUM, world);

        runner_interface_evaluate_electrostatics_3g_part_2(
            &nlocal, &nghost, &natoms, icomm_fortran, &runner_elec_energy,
            runner_elec_forces.data(), de_dq, &de_dq_sum_global, runner_elec_d_energy_d_strain);

        // Add electrostatic interactions to short-range results
        committee_energy[i] += runner_elec_energy - screening_energy;

        for (ii = 0; ii < nall * 3; ii++)
          committee_force[nall * 3 * i + ii] += runner_elec_forces[ii] - screening_forces[ii];

        for (ii = 0; ii < 9; ii++)
          committee_d_energy_d_strain[i * 9 + ii] +=
              runner_elec_d_energy_d_strain[ii] - screening_d_energy_d_strain[ii];
      }

    } else if (nnp_generation == 4) {

      // Part 1. For each committee member:
      //   - collect predicted electronegativities and hardness values on root.
      //   - perform QeQ to obtain charges for all atoms
      //   - transfer charges back to each process (local and ghost atoms)
      for (int i = 0; i < num_committee_members; i++) {

        int icomm_fortran = i + 1;

        // Collect electronegativities and hardness values from each process into
        // one structure on the root process for the computation of charges via
        // QeQ.
        std::vector<double> electronegativity_global(natoms, 0.0);
        std::vector<double> hardness_global(natoms, 0.0);
        std::vector<double> q_global(natoms, 0.0);

        pack_atomic_property(rank, size, natoms, inum, ilist, tag,
                             &committee_electronegativity[nmax * i],
                             electronegativity_global.data());

        pack_atomic_property(rank, size, natoms, inum, ilist, tag, &committee_hardness[nmax * i],
                             hardness_global.data());

        if (rank == 0) {
          // compute charges using qeq on root using global structure
          runner_interface_compute_charges_4g(
              &natoms, &total_charge, electronegativity_global.data(), hardness_global.data(),
              q_global.data(), &luse_prev_q, icomm_fortran);
        }
        MPI_Barrier(world);

        // unpack global charges from root to respective processes
        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 1, q_global.data(),
                                       atomic_charge);

        // Communicate qeq charges from local atoms to ghost atoms
        // for screening calculation
        commstyle = COMM_AT_CHARGE;
        comm->forward_comm(this);

        for (ii = 0; ii < nmax; ii++) committee_atomic_charge[nmax * i + ii] = atomic_charge[ii];
      }

      std::vector<double> committee_d_energy_d_q(nall * num_committee_members, 0.0);

      // Perform short-range prediction for all committee members at once.
      runner_interface_short_range_4g(&nlocal, &nghost, &inum, &nmax, ilist,
                                      committee_atomic_charge.data(), committee_energy.data(),
                                      committee_force.data(), committee_d_energy_d_strain.data(),
                                      committee_d_energy_d_q.data());
      MPI_Barrier(world);

      // Part 2. For each committee member:
      //   - Calculate electrostatic screening
      //   - sum up short-range and screening dE/dQ
      //   - determine lagrange charges, evaluate electrostatics, and sum up
      //     contributions
      //   - calculate force contributions of dChi/dxyz and dHardness/dxyz.
      for (int i = 0; i < num_committee_members; i++) {

        int icomm_fortran = i + 1;

        double screening_energy = 0.0;
        std::vector<double> screening_forces(nall * 3, 0.0);
        memset(screening_de_dq, 0, nall * (sizeof *screening_de_dq));
        double screening_d_energy_d_strain[9] = {0.0};

        // Apply screening
        runner_interface_calc_screening(&nlocal, &nghost, &committee_atomic_charge[nmax * i],
                                        &screening_energy, screening_forces.data(), screening_de_dq,
                                        screening_d_energy_d_strain);

        // Communicate screening de_dq from ghost atoms to local atoms
        commstyle = COMM_SCREENING_DEDQ;
        comm->reverse_comm(this);

        // Determine sum of de_dq of local atoms across all procs.
        memset(de_dq, 0, nall * (sizeof *de_dq));
        for (j = 0; j < inum; j++) {
          ii = ilist[j];
          de_dq[ii] += committee_d_energy_d_q[nall * i + ii] - screening_de_dq[ii];
        }

        double runner_elec_energy = 0.0;
        std::vector<double> elec_force_global(natoms * 3, 0.0);
        std::vector<double> de_dq_global(natoms, 0.0);
        double runner_elec_d_energy_d_strain[9] = {0.0};
        std::vector<double> lagrange_global(natoms, 0.0);

        // Communicate de_dq and pack into global structure for
        // determination of lagrange charges.

        // communicate ghost atom contributions to local atoms
        commstyle = COMM_DEDQ;
        comm->reverse_comm(this);

        pack_atomic_property(rank, size, natoms, inum, ilist, tag, de_dq, de_dq_global.data());

        if (rank == 0) {
          // serial step determining lagrange charges and
          // electrostatic contribution from global de_dq
          runner_interface_evaluate_electrostatics_4g_part_1(
              &natoms, de_dq_global.data(), &runner_elec_energy, elec_force_global.data(),
              runner_elec_d_energy_d_strain, lagrange_global.data(), icomm_fortran);
        }

        MPI_Barrier(world);

        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 1,
                                       lagrange_global.data(), lagrange_charges);

        // communicate lagrange charges from local atoms to ghost
        // atoms for calculation of force trick part 2
        commstyle = COMM_LAMBDA_CHARGE;
        comm->forward_comm(this);

        std::vector<double> runner_elec_forces(nall * 3, 0.0);
        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 3,
                                       elec_force_global.data(), runner_elec_forces.data());

        // Apply remaining force contributions from predicted
        // electronegativities and lagrange charges to
        // electrostatic forces.
        runner_interface_evaluate_electrostatics_4g_part_2(
            &nlocal, &nghost, icomm_fortran, lagrange_charges, &committee_atomic_charge[nmax * i],
            runner_elec_forces.data(), runner_elec_d_energy_d_strain);

        // Add electrostatic interactions to short-range results
        committee_energy[i] += runner_elec_energy - screening_energy;

        for (ii = 0; ii < nall * 3; ii++)
          committee_force[nall * 3 * i + ii] += runner_elec_forces[ii] - screening_forces[ii];

        for (ii = 0; ii < 9; ii++)
          committee_d_energy_d_strain[i * 9 + ii] +=
              runner_elec_d_energy_d_strain[ii] - screening_d_energy_d_strain[ii];
      }
    }
  }

  // Copy results from RuNNer back into LAMMPS atom array

  // Forces
  irunner = 0;
  for (i = 0; i < num_committee_members; i++) {
    for (ii = 0; ii < nall; ii++) {
      for (jj = 0; jj < 3; jj++) {
        f[ii][jj] += committee_force[irunner] / cfenergy * cflength /
            num_committee_members;    // runner_force is a vector
        irunner++;
      }
    }
  }

  // write individual committee forces into f_comm array
  if (lwrite_f_comm) {
    int flag, cols, ghost;
    int idx = atom->find_custom_ghost("f_comm", flag, cols, ghost);
    double **f_comm = atom->darray[idx];
    for (jj = 0; jj < 3; jj++) {                       // xyz components
      for (i = 0; i < num_committee_members; i++) {    // committee members
        irunner = jj + i * nall * 3;
        memset(committee_storage, 0, nmax * (sizeof *committee_storage));
        for (ii = 0; ii < nall; ii++) {    // nlocal + nghost atoms
          committee_storage[ii] = committee_force[irunner];
          irunner += 3;
        }
        commstyle = COMM_COMMITTEE_STORAGE;
        comm->reverse_comm(
            this);    // communicate forces so that local atoms have full contribution
        for (ii = 0; ii < nall; ii++)
          f_comm[ii][jj + i * 3] =
              committee_storage[ii] / cfenergy * cflength;    // copy into custom d2_array
      }
    }
  }

  // write individual committee charges into q_comm array
  if (lwrite_q_comm) {
    int flag, cols, ghost;
    int idx = atom->find_custom_ghost("q_comm", flag, cols, ghost);
    double **q_comm = atom->darray[idx];
    for (i = 0; i < num_committee_members; i++) {    // committee members
      for (ii = 0; ii < nall; ii++) q_comm[ii][i] = committee_atomic_charge[ii + i * nmax];
      // copy into custom d2_array
    }
  }

  // Potential energy
  if (eflag_global) eng_vdwl = 0.0;
  for (i = 0; i < num_committee_members; i++) {
    if (eflag_global) eng_vdwl = eng_vdwl + committee_energy[i] / cfenergy / num_committee_members;
  }

  // Write committee energies into pair compute vector
  memset(pvector, 0, num_committee_members * (sizeof *pvector));
  for (i = 0; i < num_committee_members; i++) pvector[i] = committee_energy[i] / cfenergy;

  // Charges if charge atom style is used
  if (q != NULL) {
    // The q array does not seem to get reset to zero every timestep by LAMMPS
    memset(q, 0, nmax * (sizeof *q));
    for (ii = 0; ii < nall; ii++) {
      for (jj = 0; jj < num_committee_members; jj++) {
        q[ii] += committee_atomic_charge[ii + jj * nmax] / num_committee_members;
      }
    }
  }

  // In case of 2G HDNNP, virial could be calculated via F dot r.
  if (vflag_fdotr) virial_fdotr_compute();

  // Virial is -1.0 * d_energy_d_strain
  if (vflag_global) {
    memset(virial, 0, 9 * (sizeof *virial));
    for (i = 0; i < num_committee_members; i++) {
      virial[0] -= committee_d_energy_d_strain[0 + 9 * i] / cfenergy / num_committee_members;
      virial[1] -= committee_d_energy_d_strain[4 + 9 * i] / cfenergy / num_committee_members;
      virial[2] -= committee_d_energy_d_strain[8 + 9 * i] / cfenergy / num_committee_members;
      virial[3] -= committee_d_energy_d_strain[1 + 9 * i] / cfenergy / num_committee_members;
      virial[4] -= committee_d_energy_d_strain[6 + 9 * i] / cfenergy / num_committee_members;
      virial[5] -= committee_d_energy_d_strain[7 + 9 * i] / cfenergy / num_committee_members;
    }
  }

  // Handling of feature extrapolations
  // - Adding extrapolation warnings to log and screen if `show_extrap` is set to yes.
  // - Stops the MD simulation if the number of extrapolations exceeds
  // the threshold defined by `max_extrap`.
  // - Resets the number of total number of recorded extrapolations during a simulation if
  // the timestep is a multiple of `reset_ew_freq` and larger than 0.
  // - Prints a summary of the recorded extrapolations at every interval until the timestep is
  // a multiple of `sum_ew_freq` and larger than 0.
  if (lcheck_extrap) {

    bigint timestep = update->ntimestep;

    // Add extrapolation warnings to log
    if (lshow_ew) {
      char *c_ptr_extrap_msg = nullptr;
      int len_extrap_msg = 0;

      runner_interface_extrapolation_warnings(&c_ptr_extrap_msg, &len_extrap_msg,
                                              global_atom_ids.data(), &nlocal);

      if (rank == 0) {
        // Write extrapolation message owned by rank 0
        if (len_extrap_msg > 0) {
          utils::logmesg(lmp, "RuNNer2: Feature extrapolations on process 0 at timestep {:d}\n",
                         timestep);
          std::string extrap_msg = std::string(c_ptr_extrap_msg);
          utils::logmesg(lmp, extrap_msg);
        }
      }

      // Gather and print warnings from other ranks
      std::vector<int> len_extrap_msg_array;
      if (rank == 0) len_extrap_msg_array.resize(size);

      MPI_Gather(&len_extrap_msg, 1, MPI_INT, len_extrap_msg_array.data(), 1, MPI_INT, 0, world);

      if (rank == 0) {
        for (int irank = 1; irank < size; irank++) {
          if (len_extrap_msg_array[irank] > 0) {
            std::vector<char> c_extrap_msg(len_extrap_msg_array[irank]);
            MPI_Status mpi_stat;
            MPI_Recv(c_extrap_msg.data(), len_extrap_msg_array[irank], MPI_BYTE, irank, 0, world,
                     &mpi_stat);

            utils::logmesg(lmp,
                           "RuNNer2: Feature extrapolations on process {:d} at timestep {:d}\n",
                           irank, timestep);
            std::string extrap_msg = std::string(c_extrap_msg.data());
            utils::logmesg(lmp, extrap_msg);
          }
        }
      } else if (len_extrap_msg > 0) {
        MPI_Send(c_ptr_extrap_msg, len_extrap_msg, MPI_BYTE, 0, 0, world);
      }
      c_ptr_extrap_msg = nullptr;
    }

    // Total number of extrapolation accumulated on this process during the simulation
    bigint local_extrap_count_total = 0;
    // Number of extrapolation accumulated on this process during this this timestep
    bigint extrap_count_timestep = 0;
    // Sets the flag `lreset` to reset the total extrapolation count if the timestep
    // is a multiple of reset_ew_freq and larger than zero
    bool lreset = false;
    if (reset_ew_freq > 0 && timestep % reset_ew_freq == 0 && timestep > 0) { lreset = true; }

    // Retrieve the number of extrapolations during this timestep and during the whole simulation
    // on each process and reset the latter if `lreset` is true.
    runner_interface_extrapolation_count(&extrap_count_timestep, &local_extrap_count_total,
                                         &lreset);

    // Number of extrapolations recorded on this process (reset at every summary)
    local_extrap_sum += extrap_count_timestep;

    // Total number of extrapolation accumulated across all processes
    bigint global_extrap_count_total = 0;
    MPI_Reduce(&local_extrap_count_total, &global_extrap_count_total, 1, MPI_LMP_BIGINT, MPI_SUM, 0,
               world);

    // Abort simulation if the total number of extrapolations across all processes exceeds `max_extrap`
    if (max_extrap > -1 && global_extrap_count_total > max_extrap && rank == 0) {
      error->one(
          FLERR,
          "Maximal number of allowed extrapolations have been exceeded during the simulation!\n"
          "Current extrapolation count: {:10.3e}",
          static_cast<double>(global_extrap_count_total));
    }

    // Prints a summary of the recorded extrapolations at every interval until the timestep is
    // a multiple of `sum_ew_freq` and larger than 0.
    if (sum_ew_freq > 0 && timestep % sum_ew_freq == 0 && timestep > 0) {
      bigint global_extrap_sum = 0;
      MPI_Reduce(&local_extrap_sum, &global_extrap_sum, 1, MPI_LMP_BIGINT, MPI_SUM, 0, world);

      if (rank == 0) {
        utils::logmesg(lmp,
                       "RuNNer2: HDNNP Extrapolation Summary\n"
                       "Timestep: {:10d} Sum: {:10.3e} Frequency: {:10.3e} per timestep\n",
                       timestep, static_cast<double>(global_extrap_sum),
                       static_cast<double>(global_extrap_sum) / static_cast<double>(sum_ew_freq));
      }
      local_extrap_sum = 0;
    }
    // Deallocates the character array containing the extrapolation message on the Fortran side
    // and frees up the internal memory of the `ExtrapolationHandler` (see RuNNer 2 documentation)
    runner_interface_dealloc_extrapolation_warnings();
  }

  MPI_Barrier(world);
  runner_interface_finalize_step();
}

void PairRuNNer::settings(int narg, char **arg)
{
  int iarg = 0;

  // default settings
  directory = utils::strdup("./");
  cflength = 1.0;
  cfenergy = 1.0;
  luse_prev_q = false;
  lwrite_f_comm = false;
  lwrite_q_comm = false;
  total_charge = 0.0;
  lcheck_extrap = false;
  max_extrap = 100;
  lshow_ew = false;
  sum_ew_freq = 0;
  reset_ew_freq = 0;
  nextra = 1;    // defaults to committee size of 1

  while (iarg < narg) {
    // set RuNNer potential directory
    if (strcmp(arg[iarg], "dir") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      delete[] directory;
      directory = utils::strdup(arg[iarg + 1]);
      iarg += 2;
      // length unit conversion factor
    } else if (strcmp(arg[iarg], "cflength") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      cflength = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
      // energy unit conversion factor
    } else if (strcmp(arg[iarg], "cfenergy") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      cfenergy = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "f_comm") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      lwrite_f_comm = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "q_comm") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      lwrite_q_comm = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "committee_size") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      nextra = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "total_charge") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      total_charge = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "use_prev_q") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      luse_prev_q = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "check_extrap") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      lcheck_extrap = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "max_extrap") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      max_extrap = utils::bnumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "show_ew") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      lshow_ew = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "sum_ew_freq") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      sum_ew_freq = utils::bnumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "reset_ew_freq") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      reset_ew_freq = utils::bnumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal pair_style command");
  }

  // check if linked to the correct RuNNer library API version
  if (runner_lammps_api_version() != 2)
    error->all(FLERR,
               "RuNNer LAMMPS wrapper API version is not compatible "
               "with this version of LAMMPS");
}

void PairRuNNer::allocate()
{
  allocated = 1;                 // sets allocated flag, checked in coeff function
  int np1 = atom->ntypes + 1;    // +1 because typeID start at 1 not 0

  setflag = memory->create(setflag, np1, np1, "pair:setflag");
  cutsq = memory->create(cutsq, np1, np1, "pair:cutsq");
  map = new int[np1];
}

void PairRuNNer::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ntypes = atom->ntypes;

  // We only have the mapping of the element to type in the pair_coeff
  // line narg 0 and 1 are two wildcards * * for I,J, so mapping
  // declaration starts at 2
  if (narg != (2 + ntypes))
    error->all(FLERR, "Number of arguments {} is not correct, it should be {}", narg, 2 + ntypes);

  // iarg - 1, because we start at 2, so first entry is 1.
  for (int iarg = 2; iarg < narg; iarg++)
    map[iarg - 1] = utils::inumeric(FLERR, arg[iarg], false, lmp);

  // clear setflag since coeff() might be called once with I,J = * *
  for (int iat = 1; iat <= ntypes; iat++)
    for (int jat = iat; jat <= ntypes; jat++) setflag[iat][jat] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  int count = 0;
  for (int iat = 1; iat <= ntypes; iat++)
    for (int jat = iat; jat <= ntypes; jat++)
      if (map[iat] >= 0 && map[jat] >= 0) {
        setflag[iat][jat] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

void PairRuNNer::init_style()
{
  // Require newton pair on
  // Switches reverse communication on on every time step, which adds
  // data from ghost atoms to corresponding local atoms Required for
  // RuNNer forces
  if (force->newton_pair != 1) error->all(FLERR, "Pair style runner requires newton pair on");

  // request full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);

  // Read model coefficients and do initialization on RuNNer side.
  // Returns max cutoff for LAMMPS neighborlist.
  // Also returns booleans if additional atomic properties are
  // predicted, which need to be communicated between local and ghost
  // atoms.
  int n_directory_len = strlen(directory);
  int rank = comm->me;
  int size = comm->nprocs;

  runner_lammps_interface_init(directory, &n_directory_len, &cutoff, &cfenergy, &cflength,
                               &nnp_generation, &num_committee_members, &lhirshfeld_vdw, &ltwo_body,
                               &lcheck_extrap, &rank, &size);

  // In the 2G case, this has a performance benefit when multiple MPI processes
  // are used.
  if (nnp_generation == 2) no_virial_fdotr_compute = 0;    // Overwrite default flag

  // Error checking for output by compute pair command
  if (nextra == num_committee_members) {
    // array for storing committee energies for output by compute pair command
    if (pvector) delete[] pvector;
    pvector = new double[nextra];
  } else {
    error->all(FLERR,
               "committee_size specification for pair style runner "
               "inconsistent between lammps input and input.nn.");
  }
  // Error checking for custom f_comm per-atom array
  if (lwrite_f_comm) {
    int flag, cols, ghost;
    int idx = atom->find_custom_ghost("f_comm", flag, cols, ghost);
    if (idx < 0 || ghost == 0)
      error->all(FLERR, "f_comm array not correctly specified in fix for pair_style runner.");
    if (flag != 1)
      error->all(FLERR, "f_comm array needs to be of type double for pair_style runner.");
    if (cols != num_committee_members * 3)
      error->all(FLERR,
                 "Not enough columns specified for f_comm array in fix for pair_style runner.");
  }

  // Error checking for custom q_comm per-atom array
  if (lwrite_q_comm) {
    int flag, cols, ghost;
    int idx = atom->find_custom_ghost("q_comm", flag, cols, ghost);
    if (idx < 0 || ghost == 0)
      error->all(FLERR, "q_comm array not correctly specified in fix for pair_style runner.");
    if (flag != 1)
      error->all(FLERR, "q_comm array needs to be of type double for pair_style runner.");
    if (cols != num_committee_members)
      error->all(FLERR,
                 "Not enough columns specified for q_comm array in fix for pair_style runner.");
  }
}

double PairRuNNer::init_one(int /*i*/, int /*j*/)
{
  return cutoff;
}

/* ----------------------------------------------------------------------
communication between local and ghost atoms
------------------------------------------------------------------------- */

int PairRuNNer::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m = 0;

  switch (commstyle) {
    case COMM_HIRSH_VOLUME:
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = hirshfeld_volume[j];
      }
      break;
    case COMM_AT_CHARGE:
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = atomic_charge[j];
      }
      break;
    case COMM_ELEC_NEGATIVITY:
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = electronegativity[j];
      }
      break;
    case COMM_DEDQ:
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = de_dq[j];
      }
      break;
    case COMM_SCREENING_DEDQ:
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = screening_de_dq[j];
      }
      break;
    case COMM_LAMBDA_CHARGE:
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = lagrange_charges[j];
      }
      break;
    case COMM_COMMITTEE_STORAGE:
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = committee_storage[j];
      }
      break;
  }
  return m;
}

void PairRuNNer::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m = 0, last = first + n;

  switch (commstyle) {
    case COMM_HIRSH_VOLUME:
      for (i = first; i < last; i++) hirshfeld_volume[i] = buf[m++];
      break;
    case COMM_AT_CHARGE:
      for (i = first; i < last; i++) atomic_charge[i] = buf[m++];
      break;
    case COMM_ELEC_NEGATIVITY:
      for (i = first; i < last; i++) electronegativity[i] = buf[m++];
      break;
    case COMM_DEDQ:
      for (i = first; i < last; i++) de_dq[i] = buf[m++];
      break;
    case COMM_SCREENING_DEDQ:
      for (i = first; i < last; i++) screening_de_dq[i] = buf[m++];
      break;
    case COMM_LAMBDA_CHARGE:
      for (i = first; i < last; i++) lagrange_charges[i] = buf[m++];
      break;
    case COMM_COMMITTEE_STORAGE:
      for (i = first; i < last; i++) committee_storage[i] = buf[m++];
      break;
  }
}

int PairRuNNer::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m = 0, last = first + n;

  switch (commstyle) {
    case COMM_HIRSH_VOLUME:
      for (i = first; i < last; i++) buf[m++] = hirshfeld_volume[i];
      break;
    case COMM_AT_CHARGE:
      for (i = first; i < last; i++) buf[m++] = atomic_charge[i];
      break;
    case COMM_ELEC_NEGATIVITY:
      for (i = first; i < last; i++) buf[m++] = electronegativity[i];
      break;
    case COMM_DEDQ:
      for (i = first; i < last; i++) buf[m++] = de_dq[i];
      break;
    case COMM_SCREENING_DEDQ:
      for (i = first; i < last; i++) buf[m++] = screening_de_dq[i];
      break;
    case COMM_LAMBDA_CHARGE:
      for (i = first; i < last; i++) buf[m++] = lagrange_charges[i];
      break;
    case COMM_COMMITTEE_STORAGE:
      for (i = first; i < last; i++) buf[m++] = committee_storage[i];
      break;
  }
  return m;
}

void PairRuNNer::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m = 0;

  switch (commstyle) {
    case COMM_HIRSH_VOLUME:
      for (i = 0; i < n; i++) {
        j = list[i];
        hirshfeld_volume[j] += buf[m++];
      }
      break;
    case COMM_AT_CHARGE:
      for (i = 0; i < n; i++) {
        j = list[i];
        atomic_charge[j] += buf[m++];
      }
      break;
    case COMM_ELEC_NEGATIVITY:
      for (i = 0; i < n; i++) {
        j = list[i];
        electronegativity[j] += buf[m++];
      }
      break;
    case COMM_DEDQ:
      for (i = 0; i < n; i++) {
        j = list[i];
        de_dq[j] += buf[m++];
      }
      break;
    case COMM_SCREENING_DEDQ:
      for (i = 0; i < n; i++) {
        j = list[i];
        screening_de_dq[j] += buf[m++];
      }
      break;
    case COMM_LAMBDA_CHARGE:
      for (i = 0; i < n; i++) {
        j = list[i];
        lagrange_charges[j] += buf[m++];
      }
      break;
    case COMM_COMMITTEE_STORAGE:
      for (i = 0; i < n; i++) {
        j = list[i];
        committee_storage[j] += buf[m++];
      }
      break;
  }
}

void PairRuNNer::pack_structure(int rank, int size, int natoms, int inum, int *ilist, tagint *tag,
                                double **x, int *runner_types, double *xyz_global, int *z_global)
{
  int i, ii;
  int start;
  int tag_minus_one;

  // First collect all Cartesian coordinates.
  std::vector<double> xyz_local(natoms * 3, 0.0);

  double xtmp, ytmp, ztmp;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    // positions of atom i
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // get atom ID of atom i
    tag_minus_one = static_cast<int>(tag[i]) - 1;    // tags start at 1
    start = tag_minus_one * 3;

    // store at atom ID position of local array to get consistent order between timesteps
    xyz_local[start] = xtmp;
    xyz_local[start + 1] = ytmp;
    xyz_local[start + 2] = ztmp;
  }
  MPI_Barrier(world);

  // Communicate local positions to root process.
  MPI_Reduce(xyz_local.data(), xyz_global, natoms * 3, MPI_DOUBLE, MPI_SUM, 0, world);

  // Second collect the atomic numbers z.
  std::vector<int> z_local(natoms, 0);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    // get atom ID of atom i
    tag_minus_one = static_cast<int>(tag[i]) - 1;
    z_local[tag_minus_one] = runner_types[i];
  }
  MPI_Barrier(world);

  // Communicate local atomic numbers to root process.
  MPI_Reduce(z_local.data(), z_global, natoms, MPI_INT, MPI_SUM, 0, world);
}

void PairRuNNer::pack_atomic_property(int rank, int size, int natoms, int inum, int *ilist,
                                      tagint *tag, double *local_property, double *global_property)
{
  int i, ii;
  int tag_minus_one;

  // Prepare an array of the local properties that is natoms long,
  // where the local atoms of each process are sorted according to their unique atom ID
  // to achieve a consistent order between timesteps.
  std::vector<double> local_property_sorted(natoms, 0.0);

  // Clear global property as well
  memset(global_property, 0, natoms * (sizeof *global_property));

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    tag_minus_one = static_cast<int>(tag[i]) - 1;
    local_property_sorted[tag_minus_one] = local_property[i];
  }
  MPI_Barrier(world);

  // Communicate local charges to root process.
  MPI_Reduce(local_property_sorted.data(), global_property, natoms, MPI_DOUBLE, MPI_SUM, 0, world);
}

void PairRuNNer::unpack_local_atomic_properties(int rank, int size, int natoms, int inum,
                                                int *ilist, tagint *tag, int nprop,
                                                double *global_properties, double *local_properties)
{
  int i, ii, iprop;
  int start;
  int tag_minus_one;

  // Broadcast global property results
  MPI_Bcast(global_properties, natoms * nprop, MPI_DOUBLE, 0, world);

  // sort back from atom ID position in global array to local order
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    tag_minus_one = static_cast<int>(tag[i]) - 1;
    start = tag_minus_one * nprop;

    for (iprop = 0; iprop < nprop; iprop++) {
      local_properties[i * nprop + iprop] = global_properties[start + iprop];
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairRuNNer::memory_usage()
{
  double bytes = Pair::memory_usage();
  // 7 per-atom double arrays: atomic_charge, hirshfeld_volume, electronegativity,
  // lagrange_charges, de_dq, screening_de_dq, committee_storage
  bytes += (double) nmax * 7 * sizeof(double);
  return bytes;
}
