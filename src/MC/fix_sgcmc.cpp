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
 * Parallel Monte-Carlo code for the semi-grandcanonical ensemble (SGC)
 * and the variance-constrained semi-grandcanonical ensemble (VC-SGC).
 *
 * See Sadigh et al., Phys. Rev. B 85, 184203 (2012) for a
 * description of the algorithm.
 *
 * Code author: Alexander Stukowski (stukowski@mm.tu-darmstadt.de)
 *
 * Updates for integrtion into LAMMPS: Aidan Thompson, SNL and Axel Kohlmeyer, Temple U
------------------------------------------------------------------------- */

#include "fix_sgcmc.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "integrate.h"
#include "kspace.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "universe.h"
#include "update.h"

#include "pair_eam.h"
#include "random_park.h"

#include <algorithm>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/*********************************************************************
 * Constructs the fix object and parses the input parameters
 * that control the Monte Carlo routine.
 *********************************************************************/
FixSemiGrandCanonicalMC::FixSemiGrandCanonicalMC(LAMMPS *_lmp, int narg, char **arg) :
    Fix(_lmp, narg, arg), random(nullptr), localRandom(nullptr), neighborList(nullptr),
    pairEAM(nullptr), compute_pe(nullptr)
{
  scalar_flag = 0;
  vector_flag = 1;
  extvector = 0;
  global_freq = 1;

  // Specifies the number of output fields this fix produces for thermo output.
  // It calculates the
  //    - Number of accepted trial moves
  //    - Number of rejected trial moves
  //    - Atom counts for each species.
  size_vector = 2 + atom->ntypes;

  // Let LAMMPS know the number of data values per atom to transfer in MPI communication.
  comm_forward = 4;
  comm_reverse = 3;

  samplingWindowUserSize = 0;
  samplingWindowPosition = 5;
  nAcceptedSwaps = 0;
  nRejectedSwaps = 0;
  kappa = 0;
  serialMode = false;

  if (domain->triclinic)
    error->all(FLERR, "Fix sgcmc does not support non-orthogonal simulation boxes.");

  if (narg < 6) utils::missing_cmd_args(FLERR, "fix sgcmc", error);

  // Parse the number of MD timesteps to do between MC.
  nevery_mdsteps = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery_mdsteps <= 0) error->all(FLERR, "Invalid number of MD timesteps {}", nevery_mdsteps);
  if (comm->me == 0) utils::logmesg(lmp, "  SGC - Number of MD timesteps: {}\n", nevery_mdsteps);

  // Parse the fraction of atoms swaps attempted during each cycle.
  swap_fraction = utils::numeric(FLERR, arg[4], false, lmp);
  if ((swap_fraction < 0.0) || (swap_fraction > 1.0))
    error->all(FLERR, "Invalid fraction {} of swap atoms", swap_fraction);
  if (comm->me == 0) utils::logmesg(lmp, "  SGC - Fraction of swap atoms: {}\n", swap_fraction);

  // Parse temperature for MC.
  double temperature = utils::numeric(FLERR, arg[5], false, lmp);
  if (temperature <= 0) error->all(FLERR, "Temperature {} invalid", temperature);
  if (comm->me == 0) utils::logmesg(lmp, "  SGC - Temperature: {}\n", temperature);
  beta = 1.0 / (force->boltz * temperature);

  // Parse chemical potentials.
  int iarg = 6;
  deltamu.resize(atom->ntypes + 1);
  deltamu[0] = 0.0;
  deltamu[1] = 0.0;
  if (atom->ntypes < 2)
    error->all(FLERR, "Fix sgcmc can only be used in simulations with at least two atom types.");
  for (int i = 2; i <= atom->ntypes; i++, iarg++) {
    if (iarg >= narg) error->all(FLERR, "Too few chemical potentials specified");
    deltamu[i] = utils::numeric(FLERR, arg[iarg], false, lmp);
    if (comm->me == 0)
      utils::logmesg(lmp, "  SGC - Chemical potential of species {}: {}\n", i, deltamu[i]);
  }

  // Default values for optional parameters (where applicable).
  numSamplingWindowMoves = 8;
  seed = 324234;

  // Parse extra/optional parameters
  while (iarg < narg) {

    if (strcmp(arg[iarg], "randseed") == 0) {
      // Random number seed.
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix sgcmc randseed", error);
      seed = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (seed <= 0) error->all(FLERR, "Random number seed {} must be positive", seed);
      if (comm->me == 0) utils::logmesg(lmp, "  SGC - Random number seed: {}\n", seed);
      iarg += 2;

    } else if (strcmp(arg[iarg], "window_moves") == 0) {
      // Parse number of window moves.
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix sgcmc window_moves", error);
      numSamplingWindowMoves = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (numSamplingWindowMoves <= 0)
        error->all(FLERR, "Invalid number {} of sampling window moves", numSamplingWindowMoves);
      if (comm->me == 0)
        utils::logmesg(lmp, "  SGC - Number of sampling window moves: {}\n",
                       numSamplingWindowMoves);
      iarg += 2;

    } else if (strcmp(arg[iarg], "window_size") == 0) {
      // Parse sampling window size parameter.
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix sgcmc window_moves", error);
      samplingWindowUserSize = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if ((samplingWindowUserSize < 0.5) || (samplingWindowUserSize > 1.0))
        error->all(FLERR, "Sampling window size {} is out of range.", samplingWindowUserSize);
      if (comm->me == 0)
        utils::logmesg(lmp, "  SGC - Sampling window size: {}\n", samplingWindowUserSize);
      iarg += 2;

    } else if (strcmp(arg[iarg], "variance") == 0) {
      // Parse parameters for variance constraint ensemble.
      if (iarg + 1 + atom->ntypes > narg)
        utils::missing_cmd_args(FLERR, "fix sgcmc variance", error);
      iarg++;

      kappa = utils::numeric(FLERR, arg[iarg], false, lmp);
      if (kappa < 0) error->all(FLERR, "Variance constraint parameter must not be negative.");
      if (comm->me == 0) utils::logmesg(lmp, "  SGC - Kappa: {}\n", kappa);
      iarg++;

      targetConcentration.resize(atom->ntypes + 1);
      targetConcentration[0] = 1.0;
      targetConcentration[1] = 1.0;
      for (int i = 2; i <= atom->ntypes; i++, iarg++) {
        targetConcentration[i] = utils::numeric(FLERR, arg[iarg], false, lmp);
        targetConcentration[1] -= targetConcentration[i];
      }
      for (int i = 1; i <= atom->ntypes; i++, iarg++) {
        if ((targetConcentration[i] < 0.0) || (targetConcentration[i] > 1.0))
          error->all(FLERR, "Target concentration {} for species {} is out of range",
                     targetConcentration[i], i);
        if (comm->me == 0)
          utils::logmesg(lmp, "  SGC - Target concentration of species {}: {}\n", i,
                         targetConcentration[i]);
      }

    } else if (strcmp(arg[iarg], "serial") == 0) {
      // Switch off second rejection.
      serialMode = true;
      if (comm->me == 0)
        utils::logmesg(lmp, "  SGC - Using serial MC version without second rejection.\n");
      iarg++;

      if (comm->nprocs != 1)
        error->all(FLERR, "Cannot use serial mode Monte Carlo in a parallel simulation.");
    } else {
      error->all(FLERR, "Unknown fix sgcmc keyword: {}", arg[iarg]);
    }
  }

  // Initialize random number generators.
  random = new RanPark(lmp, seed);
  localRandom = new RanPark(lmp, seed + universe->me);
}

/*********************************************************************
 * Destructor. Cleans up the random number generators.
 *********************************************************************/
FixSemiGrandCanonicalMC::~FixSemiGrandCanonicalMC()
{
  delete random;
  delete localRandom;
}

// clang-format off

/*********************************************************************
 * The return value of this method specifies at which points the
 * fix is invoked during the simulation.
 *********************************************************************/
int FixSemiGrandCanonicalMC::setmask()
{
  // We want the MC routine to be called in between the MD steps.
  // We need the electron densities for each atom, so after the
  // EAM potential has computed them in the force routine is a good
  // time to invoke the MC routine.
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/*********************************************************************
 * This gets called by the system before the simulation starts.
 *********************************************************************/
void FixSemiGrandCanonicalMC::init()
{
  // Make sure the user has defined only one Monte-Carlo fix.
  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"sgcmc") == 0) count++;
  if (count > 1) error->all(FLERR, "More than one fix sgcmc defined.");

  // Save a pointer to the EAM potential.
  pairEAM = dynamic_cast<PairEAM*>(force->pair);
  if (!pairEAM) {
    if (comm->me == 0)
      utils::logmesg(lmp, "  SGC - Using naive total energy calculation for MC -> SLOW!\n");

    if (comm->nprocs > 1)
      error->all(FLERR, "Can not run fix vcsgc with naive total energy calculation and more than one MPI process.");

    // Create a compute that will provide the total energy of the system.
    // This is needed by computeTotalEnergy().
    char* id_pe = (char*)"thermo_pe";
    int ipe = modify->find_compute(id_pe);
    compute_pe = modify->compute[ipe];
  }
  interactionRadius = force->pair->cutforce;
  if (comm->me == 0) utils::logmesg(lmp, "  SGC - Interaction radius: {}\n", interactionRadius);

  // This fix needs a full neighbor list.
  neighbor->add_request(this, NeighConst::REQ_FULL);

  // Count local number of atoms from each species.
  const int *type = atom->type;
  const int *mask = atom->mask;
  std::vector<int> localSpeciesCounts(atom->ntypes+1, 0);
  for (int i = 0; i < atom->nlocal; i++, ++type) {
    if (mask[i] & groupbit)
      localSpeciesCounts[*type]++;
  }

  // MPI sum to get global concentrations.
  speciesCounts.resize(atom->ntypes+1);
  MPI_Allreduce(&localSpeciesCounts.front(), &speciesCounts.front(), localSpeciesCounts.size(),
                MPI_INT, MPI_SUM, world);
}

/*********************************************************************
 * Assigns the requested neighbor list to the fix.
 *********************************************************************/
void FixSemiGrandCanonicalMC::init_list(int /*id*/, NeighList *ptr)
{
  neighborList = ptr;
}

/*********************************************************************
 * Called after the EAM force calculation during each timestep.
 * This method triggers the MC routine from time to time.
 *********************************************************************/
void FixSemiGrandCanonicalMC::post_force(int /*vflag*/)
{
  if ((update->ntimestep % nevery_mdsteps) == 0)
    doMC();
}

/*********************************************************************
 * This routine does one full MC step.
 *********************************************************************/
void FixSemiGrandCanonicalMC::doMC()
{
  // Allocate array memory.
  changedAtoms.resize(atom->nmax);

  // During the last MD timestep the EAM potential routine has computed the
  // electron densities for all atoms that belong to this processor.
  // They are stored in the rho array of the PairEAM class.
  // But computing the energy change caused by flipping one atom of this processor
  // might require the electron densities of atoms that belong to other processors.
  // So we first need to fetch those electron densities for our ghost atoms now.
  fetchGhostAtomElectronDensities();

  const int *mask = atom->mask;

  // Reset counters.
  int nAcceptedSwapsLocal = 0;
  int nRejectedSwapsLocal = 0;

  int oldSpecies, newSpecies;
  std::vector<int> deltaN(atom->ntypes+1, 0);         //< Local change in number of atoms of each species.
  std::vector<int> deltaNGlobal(atom->ntypes+1, 0);   //< Global change in number of atoms of each species.

  for (int i = 0; i < numSamplingWindowMoves; i++) {

    // Reset flag array that keeps track of changed per-atom quantities.
    std::fill(changedAtoms.begin(), changedAtoms.end(), false);

    // Position the sampling window within the node's boundaries.
    // By default the size of the sampling window is the size of the processor bounds minus two cutoff radii.
    // This ensures that changing atoms in the sampling windows of two adjacent processors cannot affect
    // the same atoms in the region between the two sampling windows.
    // For debugging purposes the sampling window can be chosen larger than the default size. Then it is
    // considered an 'oversize' window and we have to exchange atom information after each and
    // and every swap step, which is very slow.
    bool oversizeWindow = placeSamplingWindow();

    /// The number of times we want to swap an atom.
    int nDice = (int)(swap_fraction * numFixAtomsLocal / numSamplingWindowMoves);

    // This number must be synchronized with the other nodes. We take the largest
    // of all nodes and skip trial moves later.
    int largestnDice;
    MPI_Allreduce(&nDice, &largestnDice, 1, MPI_INT, MPI_MAX, world);

    // The probability to do one swap step.
    double diceProbability = (double)nDice / (double)largestnDice;

    // Inner MC loop that swaps atom types.
    for (int j = 0; j < largestnDice; j++) {

      double deltaE = 0;
      std::fill(deltaN.begin(), deltaN.end(), 0);
      int selectedAtom = -1, selectedAtomNL = -1;

      // As already said above, we have to do swap steps only with a certain probability
      // to keep nodes in sync.
      if (localRandom->uniform() <= diceProbability) {

        // Choose a random atom from the pool of atoms that are inside the sampling window.
        int index = (int)(localRandom->uniform() * (double)samplingWindowAtoms.size());
        selectedAtomNL = samplingWindowAtoms[index];

        // Get the real atom index.
        selectedAtom = neighborList->ilist[selectedAtomNL];
        oldSpecies = atom->type[selectedAtom];

        // Choose the new type for the swapping atom by random.
        if (atom->ntypes > 2) {
          // Use a random number to choose the new species if there are three or more atom types.
          newSpecies = (int)(localRandom->uniform() * (atom->ntypes-1)) + 1;
          if (newSpecies >= oldSpecies) newSpecies++;
        }
        else {
          // If there are only two atom types, then the decision is clear.
          newSpecies = (oldSpecies == 1) ? 2 : 1;
        }
        deltaN[oldSpecies] = -1;
        deltaN[newSpecies] = +1;

        // Compute the energy difference that swapping this atom would cost or gain.
        if (pairEAM) {
          deltaE = computeEnergyChangeEAM(selectedAtom, selectedAtomNL, oldSpecies, newSpecies);
        } else {
          // Generic case:
          deltaE = computeEnergyChangeGeneric(selectedAtom, oldSpecies, newSpecies);
        }

        // Perform inner MC acceptance test.
        double dm = 0.0;
        if (serialMode && kappa != 0.0) {
          for (int i = 2; i <= atom->ntypes; i++)
            dm += (deltamu[i] + kappa / atom->natoms * (2.0 * speciesCounts[i] + deltaN[i])) * deltaN[i];
        }
        else {
          for (int i = 2; i <= atom->ntypes; i++)
            dm += deltamu[i] * deltaN[i];
        }
        double deltaB = -(deltaE + dm) * beta;
        if (deltaB < 0.0) {
          if (deltaB < log(localRandom->uniform())) {
            std::fill(deltaN.begin(), deltaN.end(), 0);
            selectedAtom = -1;
            deltaE = 0;
          }
        }
      }

      if (kappa != 0.0 && serialMode == false) {

        // What follows is the second rejection test for the variance-constrained
        // semi-grandcanonical method.

        // MPI sum of total change in number of particles.
        MPI_Allreduce(&deltaN.front(), &deltaNGlobal.front(), deltaN.size(), MPI_INT, MPI_SUM, world);

        // Perform outer MC acceptance test.
        // This is done in sync by all processors.
        double A = 0.0;
        for (int i = 1; i <= atom->ntypes; i++) {
          A += deltaNGlobal[i] * deltaNGlobal[i];
          A += 2.0 * deltaNGlobal[i] * (speciesCounts[i] - (int)(targetConcentration[i] * atom->natoms));
        }
        double deltaB = -(kappa / atom->natoms) * A;
        if (deltaB < 0.0) {
          if (deltaB < log(random->uniform())) {
            std::fill(deltaN.begin(), deltaN.end(), 0);
            std::fill(deltaNGlobal.begin(), deltaNGlobal.end(), 0);
            selectedAtom = -1;
          }
        }

        // Update global species counters.
        for (int i = 1; i <= atom->ntypes; i++)
          speciesCounts[i] += deltaNGlobal[i];
      }
      else if (serialMode) {
        // Update the local species counters.
        for (int i = 1; i <= atom->ntypes; i++)
          speciesCounts[i] += deltaN[i];
      }

      // Make accepted atom swap permanent.
      if (selectedAtom >= 0) {
        if (pairEAM)
          flipAtomEAM(selectedAtom, selectedAtomNL, oldSpecies, newSpecies);
        else
          flipAtomGeneric(selectedAtom, oldSpecies, newSpecies);
        nAcceptedSwapsLocal++;
      }
      else {
        nRejectedSwapsLocal++;
      }

      if (oversizeWindow) {
        // In case of an oversized sampling window we have to exchange the atom types and all other
        // per-atom quantities after each and every swap step. This is very slow and should only be used
        // for debugging purposes.
        communicateRhoAndTypes();
      }
    }

    // Finally the changed electron densities and atom types must be exchanged before
    // the sampling window is moved.
    if (!oversizeWindow)
      communicateRhoAndTypes();
  }

  // MPI sum total number of accepted/rejected swaps.
  MPI_Allreduce(&nAcceptedSwapsLocal, &nAcceptedSwaps, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nRejectedSwapsLocal, &nRejectedSwaps, 1, MPI_INT, MPI_SUM, world);

  // For (parallelized) semi-grandcanonical MC we have to determine the current concentrations now.
  // For the serial version and variance-constrained MC it has already been done in the loop.
  if (kappa == 0.0 && serialMode == false) {
    const int *type = atom->type;
    std::vector<int> localSpeciesCounts(atom->ntypes+1, 0);
    for (int i = 0; i < atom->nlocal; i++, ++type) {
      if (mask[i] & groupbit)
        localSpeciesCounts[*type]++;
    }
    MPI_Allreduce(&localSpeciesCounts.front(), &speciesCounts.front(), localSpeciesCounts.size(), MPI_INT, MPI_SUM, world);
  }
}

/*********************************************************************
 * Fetches the electron densities for the local ghost atoms
 * from the neighbor nodes.
 *********************************************************************/
void FixSemiGrandCanonicalMC::fetchGhostAtomElectronDensities()
{
  if (pairEAM) {
    // Transfer original EAM rho values.
    communicationStage = 1;
    comm->forward_comm(this);
  }
}

/*********************************************************************
 * Transfers the locally changed electron densities and atom
 * types to the neighbors.
 *********************************************************************/
void FixSemiGrandCanonicalMC::communicateRhoAndTypes()
{
  // Electron densities can have changed for real atoms as well as ghost atoms during the last MC step.
  // So we have to perform a forward and a reverse communication to keep everything in sync.
  // In the array changedAtoms we kept track of which rhos have been changed by the MC. This helps us
  // here to not overwrite values when doing the bidirectional exchange.

  if (pairEAM) {
    // Transfer changed electron densities of ghost atoms to the real atoms.
    communicationStage = 2;
    comm->reverse_comm(this);
  }

  // Transfer changed atom types and electron densities of the real atoms to the ghost atoms.
  communicationStage = 3;
  comm->forward_comm(this);
}

/*********************************************************************
 * This is for MPI communication with neighbor nodes.
 *********************************************************************/
int FixSemiGrandCanonicalMC::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                               int * /*pbc*/)
{
  int m = 0;
  if (communicationStage == 1) {
    // Send electron densities of local atoms to neighbors.
    for (int i = 0; i < n; i++) buf[m++] = pairEAM->rho[list[i]];
  } else if (communicationStage == 3) {
    if (pairEAM) {
      // Send types and rhos of real atoms to the ghost atoms of the neighbor proc.
      for (int i = 0; i < n; i++) {
        buf[m++] = atom->type[list[i]];
        buf[m++] = pairEAM->rho[list[i]];
      }
    } else {
      // Generic potential case:
      for (int i = 0; i < n; i++) {
        buf[m++] = atom->type[list[i]];
      }
    }
  }
  return m;
}

/*********************************************************************
 * This is for MPI communication with neighbor nodes.
 *********************************************************************/
void FixSemiGrandCanonicalMC::unpack_forward_comm(int n, int first, double* buf)
{
  if (communicationStage == 1) {
    // Receive electron densities of ghost atoms from neighbors.
    int last = first + n;
    for (int i = first; i < last; i++) pairEAM->rho[i] = *buf++;
  } else if (communicationStage == 3) {
    int last = first + n;
    if (pairEAM) {
      // Receive types and rhos of real atoms of the neighbor proc and assign them
      // to the local ghost atoms.
      for (int i = first; i < last; i++, buf += 2) {
        atom->type[i] = (int)buf[0];
        // We have to make sure that rhos changed locally do not get overridden by the rhos
        // sent by the neighbor procs.
        if (!changedAtoms[i])
          pairEAM->rho[i] = buf[1];
      }
    } else {
      // Generic potential case:
      for (int i = first; i < last; i++, buf += 1) {
        atom->type[i] = (int)buf[0];
      }
    }
  }
}

/*********************************************************************
 * This is for MPI communication with neighbor nodes.
 *********************************************************************/
int FixSemiGrandCanonicalMC::pack_reverse_comm(int n, int first, double* buf)
{
  int m = 0;

  // Send changed electron densities of ghost atoms to the real atoms of neighbor procs.
  int last = first + n;
  for (int i = first; i < last; i++) buf[m++] = pairEAM->rho[i];
  return m;
}

/*********************************************************************
 * This is for MPI communication with neighbor nodes.
 *********************************************************************/
void FixSemiGrandCanonicalMC::unpack_reverse_comm(int n, int *list, double* buf)
{
  // Received changed electron densities of ghost atoms of neighbor procs and assign them to our
  // real atoms.
  for (int i = 0; i < n; i++, buf++) {
    // We have to make sure that rhos changed locally do not get overridden by the rhos
    // sent by the neighbor procs.
    if (!changedAtoms[list[i]])
      pairEAM->rho[list[i]] = *buf;
  }
}

/*********************************************************************
 * Positions the sampling window inside the node's bounding box.
 *********************************************************************/
bool FixSemiGrandCanonicalMC::placeSamplingWindow()
{
  // By default the size of the sampling window is the size of the processor bounds minus two cutoff radii.
  // This ensures that changing atoms in the sampling windows of two adjacent processors cannot affect
  // the same atoms in the region between the two sampling windows.
  // For debugging purposes the sampling window can be chosen larger than the default size. Then it is
  // considered an 'oversize' window.
  bool oversizeWindow = false;

  // Align the sampling window to one of the 8 corners of the processor cell.
  double samplingWindowLo[3];
  double samplingWindowHi[3];
  double margin[3];
  for (int i = 0; i < 3; i++) {

    margin[i] = interactionRadius * 2.0;
    if (samplingWindowUserSize > 0.0) {
      margin[i] = (domain->subhi[i] - domain->sublo[i]) * (1.0 - samplingWindowUserSize);
      if (margin[i] < interactionRadius * 2.0)
        oversizeWindow = true;
    }

    double shift = (double)((samplingWindowPosition >> i) & 1) * margin[i];
    samplingWindowLo[i] = domain->sublo[i] + shift;
    samplingWindowHi[i] = domain->subhi[i] + shift - margin[i];

    // Check if processor cells are large enough.
    // Node bounds must be at least four times as large as the atom interaction radius.
    // Sampling window must be at least half as wise as the processor cell to cover the cell completely.
    if (samplingWindowHi[i] - samplingWindowLo[i] + 1e-6 < (domain->subhi[i] - domain->sublo[i]) * 0.5) {
      error->one(FLERR, "Per-node simulation cell is too small for fix sgcmc. Processor cell size must be at least 4 times cutoff radius.");
    }
  }
  // Increase counter by one.
  // Since we are only using the lower 3 bits of the integer value the alignment will
  // be the same after 8 iterations.
  samplingWindowPosition += 1;

  // Compile a list of atoms that are inside the sampling window.
  samplingWindowAtoms.resize(0);
  samplingWindowAtoms.reserve(atom->nlocal);
  numSamplingWindowAtoms = 0;
  numFixAtomsLocal = 0;

  const int *mask = atom->mask;
  for (int ii = 0; ii < neighborList->inum; ii++) {
    int i = neighborList->ilist[ii];
    if (mask[i] & groupbit) {
      numFixAtomsLocal++;
      const double* x = atom->x[i];
      // Is atom inside window region?
      if (x[0] >= samplingWindowLo[0] && x[0] < samplingWindowHi[0] &&
          x[1] >= samplingWindowLo[1] && x[1] < samplingWindowHi[1] &&
          x[2] >= samplingWindowLo[2] && x[2] < samplingWindowHi[2]) {
          // Atoms within a distance of two times the interaction radius from the cell border
          // are less often inside the sampling window than atoms in the center of the node cell,
          // which are always inside the window.
          // We therefore have to increase their probability here to make them chosen
          // as often as the core atoms.
          int multiplicity = 1;
          for (int k=0; k < 3; k++) {
            if (x[k] < domain->sublo[k] + margin[k] ||
                x[k] > domain->subhi[k] - margin[k])
              multiplicity *= 2;
          }

          for (int m = 0; m < multiplicity; m++)
            samplingWindowAtoms.push_back(ii);

          numSamplingWindowAtoms++;
        }
    }
  }

  return oversizeWindow;
}

/*********************************************************************
 * Calculates the change in energy that swapping the given
 * atom would produce. This routine is for the standard EAM potential.
 *
 * Parameters:
 *
 * flipAtom [in]
 *   This specifies the atom to be swapped. It's an index into the local list of atoms.
 *
 * flipAtomNL [in]
 *   This specifies the atom to be swapped. It's an index into the neighbor list.
 *
 * oldSpecies [in]
 *   The current species of the atom before the routine is called.
 *
 * newSpecies [in]
 *   The new species of the atom. The atom's type is not changed by this routine. It only computes the induced energy change.
 *
 * Return value:
 *   The expected change in total potential energy.
 *********************************************************************/
double FixSemiGrandCanonicalMC::computeEnergyChangeEAM(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies)
{
  double p;
  int m;
  double const* rho = pairEAM->rho;
  double* coeff;
  double new_total_rho_i = 0.0;
  double deltaE = 0.0;

  // Calculate change of electron density at the surrounding
  // sites induced by the swapped atom. Then calculate the change of embedding energy for each neighbor atom.
  // Also recalculate the total electron density at the site of the swapped atom.

  double xi = atom->x[flipAtom][0];
  double yi = atom->x[flipAtom][1];
  double zi = atom->x[flipAtom][2];

  // Loop over all neighbors of the selected atom.
  int* jlist = neighborList->firstneigh[flipAtomNL];
  int jnum = neighborList->numneigh[flipAtomNL];
  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj];

    double delx = xi - atom->x[j][0];
    double dely = yi - atom->x[j][1];
    double delz = zi - atom->x[j][2];
    double rsq = delx*delx + dely*dely + delz*delz;
    if (rsq >= pairEAM->cutforcesq) continue;

    int jtype = atom->type[j];
    double r = sqrt(rsq);

    p = r * pairEAM->rdr + 1.0;
    m = static_cast<int>(p);
    m = MIN(m, pairEAM->nr - 1);
    p -= m;
    p = MIN(p, 1.0);

    // Calculate change of pair energy ij.
    coeff = pairEAM->z2r_spline[pairEAM->type2z2r[oldSpecies][jtype]][m];
    double oldz2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    coeff = pairEAM->z2r_spline[pairEAM->type2z2r[newSpecies][jtype]][m];
    double newz2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    deltaE += (newz2 - oldz2) / r;

    // Calculate change of electron density at site j.
    coeff = pairEAM->rhor_spline[pairEAM->type2rhor[oldSpecies][jtype]][m];
    double oldrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    coeff = pairEAM->rhor_spline[pairEAM->type2rhor[newSpecies][jtype]][m];
    double newrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    double delta_rho = newrho_contr - oldrho_contr;

    // Sum total rho at site of swapped atom.
    coeff = pairEAM->rhor_spline[pairEAM->type2rhor[jtype][newSpecies]][m];
    new_total_rho_i += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

    // Calculate old embedding energy of atom j.
    p = rho[j] * pairEAM->rdrho + 1.0;
    m = static_cast<int>(p);
    m = MAX(1, MIN(m, pairEAM->nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    coeff = pairEAM->frho_spline[pairEAM->type2frho[jtype]][m];
    double oldF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

    // Calculate new embedding energy of atom j.
    p = (rho[j] + delta_rho) * pairEAM->rdrho + 1.0;
    m = static_cast<int>(p);
    m = MAX(1, MIN(m, pairEAM->nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    coeff = pairEAM->frho_spline[pairEAM->type2frho[jtype]][m];
    double newF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

    deltaE += newF - oldF;
  }

  // Compute the change in embedding energy of the changing atom.
  p = rho[flipAtom] * pairEAM->rdrho + 1.0;
  m = static_cast<int>(p);
  m = MAX(1, MIN(m, pairEAM->nrho - 1));
  p -= m;
  p = MIN(p, 1.0);
  coeff = pairEAM->frho_spline[pairEAM->type2frho[oldSpecies]][m];
  double oldF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

  p = new_total_rho_i * pairEAM->rdrho + 1.0;
  m = static_cast<int>(p);
  m = MAX(1, MIN(m, pairEAM->nrho - 1));
  p -= m;
  p = MIN(p, 1.0);
  coeff = pairEAM->frho_spline[pairEAM->type2frho[newSpecies]][m];
  double newF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

  deltaE += newF - oldF;

  return deltaE;
}

/*********************************************************************
 * Calculates the change in energy that swapping the given atom would produce.
 * This routine is for the general case of an arbitrary potential and
 * IS VERY SLOW! It computes the total energies of the system for the unmodified state
 * and for the modified state and then returns the difference of both values.
 * This routine should only be used for debugging purposes.
 *
 * Parameters:
 *
 * flipAtom [in]
 *   This specifies the atom to be swapped. It's an index into the local list of atoms.
 *
 * oldSpecies [in]
 *   The current species of the atom before the routine is called.
 *
 * newSpecies [in]
 *   The new species of the atom. The atom's type is not changed by this method. It only computes the induced energy change.
 *
 * Return value:
 *   The expected change in total potential energy.
 *********************************************************************/
double FixSemiGrandCanonicalMC::computeEnergyChangeGeneric(int flipAtom, int oldSpecies, int newSpecies)
{
  // This routine is called even when no trial move is being performed during the
  // the current iteration to keep the parallel processors in sync. If no trial
  // move is performed then the energy is calculated twice for the same state of the system.
  if (flipAtom >= 0) {
    // Change system. Perform trial move.
    atom->type[flipAtom] = newSpecies;
  }
  // Transfer changed atom types of the real atoms to the ghost atoms.
  communicationStage = 3;
  comm->forward_comm(this);

  // Calculate new total energy.
  double newEnergy = computeTotalEnergy();

  // Undo trial move. Restore old system state.
  if (flipAtom >= 0) {
    atom->type[flipAtom] = oldSpecies;
  }
  // Transfer changed atom types of the real atoms to the ghost atoms.
  communicationStage = 3;
  comm->forward_comm(this);

  // Calculate old total energy.
  double oldEnergy = computeTotalEnergy();

  // Restore the correct electron densities and forces.
  update->integrate->setup_minimal(0);
  fetchGhostAtomElectronDensities();

  return newEnergy - oldEnergy;
}

/*********************************************************************
 * Lets LAMMPS calculate the total potential energy of the system.
 *********************************************************************/
double FixSemiGrandCanonicalMC::computeTotalEnergy()
{
  int eflag = 1;
  int vflag = 0;

  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) force->kspace->compute(eflag,vflag);

  update->eflag_global = update->ntimestep;
  return compute_pe->compute_scalar();
}

/*********************************************************************
 * Flips the type of one atom and changes the electron densities
 * of nearby atoms accordingly.
 * This routine is for the case of a standard EAM potential.
 *
 * Parameters:
 *
 * flipAtom [in]
 *   This specifies the atom to be swapped. It's an index into the local list of atoms.
 *
 * flipAtomNL [in]
 *   This specifies the atom to be swapped. It's an index into the neighbor list.
 *
 * oldSpecies [in]
 *   The current species of the atom before the routine is called.
 *
 * newSpecies [in]
 *   The new type to be assigned to the atom.
 *********************************************************************/
void FixSemiGrandCanonicalMC::flipAtomEAM(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies)
{
  double p;
  int m;
  double* rho = pairEAM->rho;
  double* coeff;
  double new_total_rho_i = 0.0;

  // Change atom's type and mark it for exchange.
  atom->type[flipAtom] = newSpecies;
  changedAtoms[flipAtom] = true;

  // Rescale particle velocity vector to conserve kinetic energy.
  double vScaleFactor = sqrt(atom->mass[oldSpecies] / atom->mass[newSpecies]);
  atom->v[flipAtom][0] *= vScaleFactor;
  atom->v[flipAtom][1] *= vScaleFactor;
  atom->v[flipAtom][2] *= vScaleFactor;

  double xi = atom->x[flipAtom][0];
  double yi = atom->x[flipAtom][1];
  double zi = atom->x[flipAtom][2];

  // Loop over all neighbors of the selected atom.
  int* jlist = neighborList->firstneigh[flipAtomNL];
  int jnum = neighborList->numneigh[flipAtomNL];
  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj];

    double delx = xi - atom->x[j][0];
    double dely = yi - atom->x[j][1];
    double delz = zi - atom->x[j][2];
    double rsq = delx*delx + dely*dely + delz*delz;
    if (rsq >= pairEAM->cutforcesq) continue;

    int jtype = atom->type[j];
    double r = sqrt(rsq);
    p = r * pairEAM->rdr + 1.0;
    m = static_cast<int>(p);
    m = MIN(m, pairEAM->nr - 1);
    p -= m;
    p = MIN(p, 1.0);

    // Calculate change of electron density at site j.
    coeff = pairEAM->rhor_spline[pairEAM->type2rhor[oldSpecies][jtype]][m];
    double oldrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    coeff = pairEAM->rhor_spline[pairEAM->type2rhor[newSpecies][jtype]][m];
    double newrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    double delta_rho = newrho_contr - oldrho_contr;

    rho[j] += delta_rho;

    // Sum total rho at site of swapped atom.
    coeff = pairEAM->rhor_spline[pairEAM->type2rhor[jtype][newSpecies]][m];
    new_total_rho_i += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

    // Set the flag for this atom to indicate that its rho has changed and needs
    // to be transfered at end of MC step.
    changedAtoms[j] = true;
  }

  // Store newly calculated electron density at swapped atom site.
  rho[flipAtom] = new_total_rho_i;
}

/*********************************************************************
 * Flips the type of one atom.
 * This routine is for the generic case.
 *
 * Parameters:
 *
 * flipAtom [in]
 *   This specifies the atom to be swapped. It's an index into the local list of atoms.
 *
 * oldSpecies [in]
 *   The current species of the atom before the routine is called.
 *
 * newSpecies [in]
 *   The new type to be assigned to the atom.
 *********************************************************************/
void FixSemiGrandCanonicalMC::flipAtomGeneric(int flipAtom, int oldSpecies, int newSpecies)
{
  atom->type[flipAtom] = newSpecies;

  // Rescale particle velocity vector to conserve kinetic energy.
  double vScaleFactor = sqrt(atom->mass[oldSpecies] / atom->mass[newSpecies]);
  atom->v[flipAtom][0] *= vScaleFactor;
  atom->v[flipAtom][1] *= vScaleFactor;
  atom->v[flipAtom][2] *= vScaleFactor;

  changedAtoms[flipAtom] = true;
}

/*********************************************************************
 * Lets the fix report one of its internal state variables to LAMMPS.
 *********************************************************************/
double FixSemiGrandCanonicalMC::compute_vector(int index)
{
  if (index == 0) return nAcceptedSwaps;
  if (index == 1) return nRejectedSwaps;
  index -= 1;
  int totalAtoms = 0;
  for (int i = 0; i < (int)speciesCounts.size(); i++)
    totalAtoms += speciesCounts[i];
  if (index <= atom->ntypes)
    return (double)speciesCounts[index] / (totalAtoms > 0 ? totalAtoms : 1);
  return 0.0;
}

/*********************************************************************
 * Reports the memory usage of this fix to LAMMPS.
 *********************************************************************/
double FixSemiGrandCanonicalMC::memory_usage()
{
  return (changedAtoms.size() * sizeof(bool)) +
    (samplingWindowAtoms.size() * sizeof(int));
}

