/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:  Karl D. Hammond <hammondkd@missouri.edu>
                         University of Missouri, Columbia (USA), 2018

   Updated July 8, 2022 by the author.

   Refactored for inclusion into LAMMPS 23 June 2026
   and to use dump image for graphics: by Axel Kohlmeyer, Temple U
------------------------------------------------------------------------- */

#include "compute_frenkel.h"
#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "graphics.h"
#include "group.h"
#include "lattice.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "pair.h"
#include "text_file_reader.h"
#include "update.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>

static constexpr int BIN_GROW_SIZE = 32;    // bins start this size and grow by this
static constexpr int MAX_OCCUPANTS = 8;     // max number of occupants of one lattice site
static constexpr double BIG = 1.0e20;
static constexpr double SMALL = 1.0e-10;

// small, fixed MPI message tags for the ghost-site exchange.  Using the MPI
// rank (or rank + n*nprocs) as a tag can exceed the guaranteed MPI_TAG_UB
// (>= 32767) on large runs; the source rank is already given to MPI_Irecv, so
// one constant tag per message type is sufficient and portable.
static constexpr int TAG_COUNT = 0;
static constexpr int TAG_SITE = 1;
static constexpr int TAG_X = 2;
static constexpr int TAG_OCCUP = 3;
static constexpr int TAG_CLUST = 4;
static constexpr int TAG_RTAG = 5;    // reverse (ghost -> owner) label push: site tags
static constexpr int TAG_RCID = 6;    // reverse (ghost -> owner) label push: cluster IDs

using namespace LAMMPS_NS;
using MathConst::MY_2PI;
using MathConst::MY_PI;

static const char cite_compute_frenkel_c[] =
    "compute_frenkel command: doi:10.1016/j.cpc.2019.106862\n\n"
    "@article{Hammond2020a,\n"
    "  author  = \"Hammond, Karl D.\",\n"
    "  title   = \"Parallel Point Defect Identification in Molecular Dynamics\n"
    "             Simulations Without Post-Processing: A Compute and Dump Style\n"
    "             for {LAMMPS}\",\n"
    "  journal = \"Computer Physics Communications\",\n"
    "  volume  = 247,\n"
    "  pages   = 106862,\n"
    "  doi     = 10.1016/j.cpc.2019.106862,\n"
    "  year    = 2020,\n"
    "}\n\n";

ComputeFrenkel::ComputeFrenkel(class LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), image_objvec(nullptr), image_objarray(nullptr), image_nmax(0),
    region(nullptr), ifgroup(igroup), fgroupbit(groupbit), rescale(false), mindist(nullptr),
    site_mindist(nullptr), noccupants(nullptr), occupant_tag(nullptr), nnormal(0), normal(nullptr),
    nlatsites(0), nlatghosts(0), latsites(nullptr), latsites0(nullptr), site_tag(nullptr),
    first_local_tag(0), nlatbins{0, 0, 0, 0}, latbins(nullptr), clusterID(nullptr),
    cluster_center(nullptr), noccupied(0), occupied_cluster_ID(nullptr), old_boxlo{0.0, 0.0, 0.0},
    old_boxhi{0.0, 0.0, 0.0}, invoked_find_defects(-1), invoked_find_clusters(-1),
    invoked_construct_WS_cell(-1)
{
  if (narg != 3) error->all(FLERR, "Illegal compute frenkel command");

  if (lmp->citeme) lmp->citeme->add(cite_compute_frenkel_c);

  // every atom must have an ID, and we need a lattice to define the reference sites
  if (!atom->tag_enable)
    error->all(FLERR, "Cannot use compute style frenkel unless atoms have IDs");
  if (!domain->lattice || domain->lattice->nbasis == 0)
    error->all(FLERR, "Use of compute style frenkel with undefined lattice");

  // base-class Compute settings
  comm_reverse = 1;
  vector_flag = array_flag = peratom_flag = local_flag = image_flag = 1;
  size_peratom_cols = 0;
  size_local_rows = 0;
  size_local_cols = 5;     // cluster tag, cluster size, x, y, z
  size_vector = 3;         // vacancies (0), interstitials (1), irregulars (2)
  size_array_rows = 2;     // vacancies and interstitials
  size_array_cols = 20;    // clusters larger than this are counted in the last column
  extvector = 0;
  extarray = 0;
  vector_atom = nullptr;
  array_local = nullptr;
  memory->create(vector, size_vector, "ComputeFrenkel:vector");
  memory->create(array, size_array_rows, size_array_cols, "ComputeFrenkel:array");

  // default vacancy/interstitial search radii from the lattice spacing
  double a_max =
      MAX(MAX(domain->lattice->xlattice, domain->lattice->ylattice), domain->lattice->zlattice);
  cut_vac = 1.01 * a_max;
  cut_int = 1.42 * a_max;
  cutoff = MAX(cut_vac, cut_int);
  binwidth = cutoff;
}

/****************************************************************************/

ComputeFrenkel::~ComputeFrenkel()
{
  memory->destroy(vector);
  memory->destroy(array);
  memory->destroy(array_local);

  memory->destroy(latsites);
  memory->destroy(latsites0);
  memory->destroy(site_tag);
  memory->destroy(normal);
  memory->destroy(noccupants);
  memory->destroy(occupant_tag);
  memory->destroy(latbins);
  memory->destroy(mindist);
  memory->destroy(site_mindist);
  memory->destroy(clusterID);
  memory->destroy(cluster_center);
  memory->destroy(occupied_cluster_ID);
  memory->destroy(image_objvec);
  memory->destroy(image_objarray);
}

/****************************************************************************/

// Handle the style-specific compute_modify keywords.  The base class
// Compute::modify_params() handles the generic ones (extra/dof, dynamic/dof)
// and calls this for anything it does not recognize; we return the number of
// arguments consumed, or 0 if the keyword is not one of ours.
int ComputeFrenkel::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "drvac") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "compute_modify drvac", error);
    cut_vac = utils::numeric(FLERR, arg[1], false, lmp);
    if (cut_vac <= 0.0) error->all(FLERR, "compute_modify drvac value must be > 0.0");
    return 2;
  } else if (strcmp(arg[0], "drint") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "compute_modify drint", error);
    cut_int = utils::numeric(FLERR, arg[1], false, lmp);
    if (cut_int <= 0.0) error->all(FLERR, "compute_modify drint value must be > 0.0");
    return 2;
  } else if (strcmp(arg[0], "region") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "compute_modify region", error);
    if (strcmp(arg[1], "none") == 0) {
      region = nullptr;
    } else {
      region = domain->get_region_by_id(arg[1]);
      if (!region) error->all(FLERR, "compute_modify region {} does not exist", arg[1]);
    }
    return 2;
  } else if (strcmp(arg[0], "frenkelgroup") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "compute_modify frenkelgroup", error);
    ifgroup = group->find(arg[1]);
    if (ifgroup == -1) error->all(FLERR, "compute_modify frenkelgroup {} does not exist", arg[1]);
    fgroupbit = group->bitmask[ifgroup];
    return 2;
  } else if (strcmp(arg[0], "rescale") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "compute_modify rescale", error);
    rescale = utils::logical(FLERR, arg[1], false, lmp);
    return 2;
  } else if (strcmp(arg[0], "site_file") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "compute_modify site_file", error);
    sitefile.clear();
    if (strcmp(arg[1], "none") != 0) {
      sitefile = arg[1];
      if (!platform::file_is_readable(sitefile))
        error->one(FLERR, "Compute frenkel site file {} is not readable", sitefile);
    }
    return 2;
  }
  return 0;
}

/****************************************************************************/

void ComputeFrenkel::init()
{
  // Note that invoked_vector and invoked_array are reset inside
  // modify->compute[]->init during the initialization routines.
  invoked_find_defects = -1;
  invoked_find_clusters = -1;
  invoked_construct_WS_cell = -1;

  // Make sure we have a way to generate lattice sites
  if (!domain->lattice)
    error->all(FLERR, Error::NOLASTLINE, "Use of compute style frenkel with undefined lattice");

  // Recompute the derived search/bin cutoffs from the (possibly compute_modify
  // changed) vacancy and interstitial radii, so drvac/drint stay consistent
  // with the values used for binning, ghost exchange, and nearest-site search.
  cutoff = MAX(cut_vac, cut_int);
  binwidth = cutoff;

  // Ghost atoms must reach at least one search cutoff beyond the subdomain or
  // defects near subdomain boundaries can be silently misidentified.  The
  // actual ghost range is not finalized until comm->setup() (after this), so
  // estimate it the way the neighbor/comm code will (pair cutoff + skin, or a
  // user-set comm cutoff), as compute rdf does.
  double ghostcut = comm->cutghostuser;
  if (force->pair) ghostcut = MAX(ghostcut, force->pair->cutforce + neighbor->skin);
  if ((comm->me == 0) && (comm->nprocs > 1) && (cutoff > ghostcut))
    error->warning(FLERR,
                   "Compute frenkel cutoff {} exceeds the communication cutoff {}; "
                   "use the comm_modify cutoff command to increase the ghost range",
                   cutoff, ghostcut);

  create_lattice_sites();
  put_sites_in_bins();

  // Initialize occupancy lists for all (local) sites
  memory->destroy(noccupants);
  memory->destroy(occupant_tag);
  memory->create(noccupants, nlatsites, "ComputeFrenkel:noccupants");
  memory->create(occupant_tag, nlatsites, MAX_OCCUPANTS, "ComputeFrenkel:occupant_tag");
  std::fill_n(noccupants, nlatsites, 0);
  std::fill_n(&occupant_tag[0][0], nlatsites * MAX_OCCUPANTS, tagint(-1));    // contiguous block

  construct_WS_cell();

  // Build neighbor list, which should never change until the sites do.  Note
  // ghost sites are deliberately NOT binned: the atom-to-site assignment in
  // find_defects must only ever pick a local site (occupancy is counted by the
  // owning rank).  Cross-boundary clusters still connect because each ghost
  // site's neighbor list (from its clamped edge bin) points back at the nearby
  // local sites, which find_clusters merges.
  exchange_lattice_ghosts();    // should skip occupancies until second one
  construct_site_nlists();
  // These two are redone every time find_defects() runs, so the initial pass
  // here is partly redundant; it is kept because it is a one-time setup cost
  // (negligible against the per-invocation defect search) and primes the ghost
  // sites and neighbor lists before the first compute call.

  // Store current box boundaries
  old_boxlo[0] = domain->boxlo[0];
  old_boxlo[1] = domain->boxlo[1];
  old_boxlo[2] = domain->boxlo[2];
  old_boxhi[0] = domain->boxhi[0];
  old_boxhi[1] = domain->boxhi[1];
  old_boxhi[2] = domain->boxhi[2];
}

/****************************************************************************/

int ComputeFrenkel::pack_reverse_comm(int n, int first, double *buf)
{
  int m, last;
  m = 0;
  last = first + n;
  for (int i = first; i < last; i++) buf[m++] = mindist[i];
  return m;
}

/****************************************************************************/

void ComputeFrenkel::unpack_reverse_comm(int n, int *list, double *buf)
{
  int j = 0, m = 0;
  double tmp;

  for (int i = 0; i < n; i++) {
    j = list[i];
    tmp = buf[m++];    // So that m doesn't get incremented TWICE
    mindist[j] = MIN(mindist[j], tmp);
  }
}

/****************************************************************************/

void ComputeFrenkel::find_defects()
{
  invoked_find_defects = update->ntimestep;

  rescale_lattice_sites();
  exchange_lattice_ghosts();    // In case sites are NOW close enough to exchange
  construct_site_nlists();

  // Allocate and zero out mindist, site_mindist, and noccupants
  memory->destroy(site_mindist);
  memory->create(site_mindist, nlatsites, "ComputeFrenkel:site_mindist");
  std::fill_n(noccupants, nlatsites, 0);
  std::fill_n(&occupant_tag[0][0], nlatsites * MAX_OCCUPANTS, tagint(-1));    // contiguous block
  std::fill_n(site_mindist, nlatsites, BIG);

  memory->destroy(mindist);
  memory->create(mindist, atom->nmax, "ComputeFrenkel:mindist");
  std::fill_n(mindist, atom->nmax, BIG);

  double cutsq = cutoff * cutoff;

  // Loop over all atoms and ghosts
  for (int n = 0; n < atom->nlocal + atom->nghost; n++) {

    // Find closest lattice bin
    int ii, jj, kk;
    find_closest_bin(atom->x[n], ii, jj, kk);

    // Find closest site within this bin
    double drsq_min = BIG;
    int closest_site = -1;
    double dx, dy, dz, drsq;
    for (int m = 0; m < nlatbins[3]; m++) {
      int s = latbins[ii][jj][kk][m];
      if (s < 0) break;    // Reached end of list
      double *r_s = latsites[s];
      dx = r_s[0] - atom->x[n][0];
      dy = r_s[1] - atom->x[n][1];
      dz = r_s[2] - atom->x[n][2];
      // Apply periodic boundary conditions, if necessary
      domain->minimum_image(FLERR, dx, dy, dz);
      drsq = dx * dx + dy * dy + dz * dz;
      if (drsq < drsq_min && drsq <= cutsq) {
        drsq_min = drsq;
        closest_site = s;
      }
    }

    // Now loop over the neighbor list of that site, to find ghosts that
    // might be even closer (ghosts are not always in the right bin)
    if (closest_site >= 0) {
      int closest_site_before = closest_site;
      do {
        closest_site_before = closest_site;
        for (tagint t : nlist[closest_site]) {
          int s = site_tag2index(t);
          if (s < 0) continue;    // Reached end of list
          double *r_s = latsites[s];
          dx = r_s[0] - atom->x[n][0];
          dy = r_s[1] - atom->x[n][1];
          dz = r_s[2] - atom->x[n][2];
          // Apply periodic boundary conditions, if necessary
          domain->minimum_image(FLERR, dx, dy, dz);
          drsq = dx * dx + dy * dy + dz * dz;
          if (drsq < drsq_min && drsq <= cutsq) {
            drsq_min = drsq;
            closest_site = s;
            break;
          }
        }
      } while (closest_site != closest_site_before);
    }

    // If a nearby site was found on this process, set mindist
    // Otherwise, leave it set to BIG (implies it's on another process)
    if (closest_site < 0) {
      mindist[n] = BIG;
      continue;
    } else
      mindist[n] = sqrt(drsq_min);

    site_mindist[closest_site] = MIN(site_mindist[closest_site], mindist[n]);

    // If inside the Wigner-Seitz cell, add it to the occupancy list (if it's
    // not already in it)
    if (inside_WS_cell(n, closest_site)) {
      // If already in site's occupancy list, we're done
      if (tag_is_already_in_occupancy_list(atom->tag[n], closest_site)) continue;
      // If not, add it to the list and increase site occupancy
      noccupants[closest_site] += 1;
      int s = next_free_occupant_tag_index(closest_site, __LINE__);
      occupant_tag[closest_site][s] = atom->tag[n];
    }
  }

  // Exchange mindist info amongst processes for ghost atoms
  //comm->reverse_comm_compute (this);
  comm->reverse_comm(this);
}

/****************************************************************************/

void ComputeFrenkel::find_clusters()
{
  invoked_find_clusters = update->ntimestep;

  memory->destroy(clusterID);
  memory->create(clusterID, nlatsites, "ComputeFrenkel:clusterID");

  // Exchange "ghost" sites with nearby processes
  // includes site tag, cluster ID, position, and number of occupants
  //  exchange_lattice_ghosts ();

  // Start each defect off in its own cluster with cluster ID = site ID
  for (int k = 0; k < nlatsites; k++) {
    if (noccupants[k] == 1)
      clusterID[k] = 0;    // "regular" sites
    else
      clusterID[k] = site_tag[k];    // interstitials, vacancies, and "irregulars"
  }

  double dx, dy, dz, drsq;
  double cutvacsq = cut_vac * cut_vac;
  double cutintsq = cut_int * cut_int;

  // For each local site, parse the neighbor list to find clusters.  This is a
  // distributed connected-components (union-find) labelling: adjacent defect
  // sites repeatedly adopt the smaller of their two cluster IDs until a global
  // fixpoint is reached.  Each outer pass:
  //   (1) exchange_lattice_ghosts() stamps each owner's current cluster ID onto
  //       the ghost copies held by neighbouring subdomains;
  //   (2) a local relaxation merges neighbouring sites (owned and ghost) by MIN;
  //   (3) push_ghost_labels_to_owners() reverse-communicates any label a ghost
  //       picked up in step (2) back to its owner (also by MIN).
  //
  // Step (3) is essential.  Without it a relabelling that can only be found on
  // the *receiver* side of a boundary lives solely on a ghost copy, is
  // overwritten by the next exchange, and is rediscovered every pass.  For a
  // compact cluster wrapping a corner where >= 2 subdomain cuts meet that made
  // the old "did any (owned-or-ghost) label change?" flag stay true forever -- a
  // livelock.  We instead iterate until no *owned* label changes on any rank,
  // which is a genuine fixpoint because owned labels only ever decrease and are
  // bounded below (so the loop is guaranteed to terminate).
  int global_changes_made = false;
  std::vector<tagint> prev_clusterID(nlatsites);
  do {
    exchange_lattice_ghosts();

    // snapshot the owned labels so we can detect real progress this pass
    for (int k = 0; k < nlatsites; k++) prev_clusterID[k] = clusterID[k];

    // local relaxation over owned sites and ghosts (only ever lowers labels)
    int done;
    do {
      done = true;
      for (int n = 0; n < nlatsites + nlatghosts; n++) {
        if (noccupants[n] == 1) continue;    // Skip "normal" sites
        for (tagint t : nlist[n]) {
          int m = site_tag2index(t);
          if (m != n && m >= 0 && noccupants[m] != 1 && clusterID[m] != clusterID[n]) {
            dx = latsites[n][0] - latsites[m][0];
            dy = latsites[n][1] - latsites[m][1];
            dz = latsites[n][2] - latsites[m][2];
            domain->minimum_image(FLERR, dx, dy, dz);
            drsq = dx * dx + dy * dy + dz * dz;
            if ((noccupants[m] >= 2 && drsq <= cutintsq) ||
                (noccupants[m] == 0 && drsq <= cutvacsq)) {
              clusterID[n] = clusterID[m] = MIN(clusterID[n], clusterID[m]);
              done = false;
            }
          }
        }
      }
    } while (!done);

    // send labels learned on ghost copies back to their owners (MIN)
    push_ghost_labels_to_owners();

    // a pass is productive only if some *owned* label actually decreased
    int changes_made = 0;
    for (int k = 0; k < nlatsites; k++)
      if (clusterID[k] != prev_clusterID[k]) {
        changes_made = 1;
        break;
      }

    // we stop only once no process lowered an owned label this pass
    MPI_Allreduce(&changes_made, &global_changes_made, 1, MPI_INT, MPI_LOR, world);

  } while (global_changes_made);

  // Find the size of each cluster by counting the sites in each cluster;
  // vacancies are counted as -1 and interstitials as +1.

  // This algorithm is O(nlatsites_global), or O(natoms), in memory if we do
  // it with global site IDs.  This works great for small systems, but if
  // you have a billion atoms, it sort of defeats the purpose of
  // parallelization.... Instead, we generate an array containing a list of
  // all cluster IDs that are non-zero and store ONLY those cluster sizes and
  // centers of mass.
  find_occupied_clusters();

  // If we found NO occupied clusters, we're done
  if (noccupied == 0) return;

  // Find the number of sites involved in each cluster (its size), and the
  // center of each cluster (NOT the center of mass, the center)
  std::vector<int> local_cluster_size(noccupied, 0);
  std::vector<int> local_cluster_nsites(noccupied, 0);
  double **local_cluster_xi;
  double **local_cluster_zeta;
  memory->create(local_cluster_xi, noccupied, 3, "ComputeFrenkel::local_cluster_xi");
  memory->create(local_cluster_zeta, noccupied, 3, "ComputeFrenkel::local_cluster_zeta");
  for (int i = 0; i < noccupied * 3; i++) local_cluster_xi[0][i] = local_cluster_zeta[0][i] = 0.0;
  int n = 0;
  double theta;
  for (int k = 0; k < nlatsites; k++) {
    if (clusterID[k] == 0) continue;
    n = clusterID2occupied_index(clusterID[k]);
    if (n < 0) {
      // should NEVER happen...
      error->warning(FLERR, "Did not find cluster index");
      continue;
    }
    if (noccupants[k] == 0)    // vacancy
      local_cluster_size[n] += -1;
    else if (noccupants[k] >= 2)    // interstitial
      local_cluster_size[n] += +1;

    // This calculates the center of mass using the image as determined from
    // the method of Bai and Breen (J. Graph. Tools 13(4): 53-60 (2008).
    // Note that their method does NOT calculate the center of mass;
    // here's what I actually do:
    //  (1) For each spatial dimension, wrap the simulation box around the
    //      complex unit circle:  z_j = exp(2*pi*i*x_j/L) = xi_j + i*zeta_j
    //  (2) What we want is the /geometric/ mean,
    //      <z> = prod_{j=1}^N z_j^(1/N)
    //  but that has N-1 non-degenerate branches, all of which produce a
    //  different value of <z>!  So we need to pick the "right" geometric mean
    //  (3) Calculate the /arithmetic/ mean in the complex plane,
    //      <z>_a = sum_{j=1}^N z_j/N
    //  (4) Back out <x>_a, the arithmetic mean position, via
    //      <x>_a = arg(<z>_a)
    //  (5) Calculate the "correct" center of mass by taking the image closest
    //    to <x>_a for each particle and averaging in the usual way
    if (noccupants[k] != 1) {
      local_cluster_nsites[n] += 1;
      if (domain->xperiodic) {
        theta = (latsites[k][0] - domain->boxlo[0]) / domain->xprd * MY_2PI;
        local_cluster_zeta[n][0] += sin(theta);
        local_cluster_xi[n][0] += cos(theta);
      } else    // If not periodic, xi = x (no projection required)
        local_cluster_xi[n][0] += latsites[k][0];
      if (domain->yperiodic) {
        theta = (latsites[k][1] - domain->boxlo[1]) / domain->yprd * MY_2PI;
        local_cluster_zeta[n][1] += sin(theta);
        local_cluster_xi[n][1] += cos(theta);
      } else
        local_cluster_xi[n][1] += latsites[k][1];
      if (domain->zperiodic) {
        theta = (latsites[k][2] - domain->boxlo[2]) / domain->zprd * MY_2PI;
        local_cluster_zeta[n][2] += sin(theta);
        local_cluster_xi[n][2] += cos(theta);
      } else
        local_cluster_xi[n][2] += latsites[k][2];
    }
  }

  // Find the global cluster size for each cluster across all processes
  cluster_size.assign(noccupied, 0);
  cluster_nsites.assign(noccupied, 0);
  MPI_Allreduce(local_cluster_size.data(), cluster_size.data(), noccupied, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(local_cluster_nsites.data(), cluster_nsites.data(), noccupied, MPI_INT, MPI_SUM,
                world);

  // Now average the center in the complex plane <z>_a = (<xi>,<zeta>) across
  // all processes for all clusters
  double **cluster_xi;
  double **cluster_zeta;
  memory->create(cluster_xi, noccupied, 3, "ComputeFrenkel:cluster_xi");
  memory->create(cluster_zeta, noccupied, 3, "ComputeFrenkel:cluster_zeta");
  for (int i = 0; i < noccupied * 3; i++) cluster_xi[0][i] = cluster_zeta[0][i] = 0.0;
  MPI_Allreduce(local_cluster_xi[0], cluster_xi[0], noccupied * 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(local_cluster_zeta[0], cluster_zeta[0], noccupied * 3, MPI_DOUBLE, MPI_SUM, world);
  memory->destroy(local_cluster_xi);
  memory->destroy(local_cluster_zeta);

  // If not periodic in that direction, then just divide xi by the cluster
  // size to get the center; else back out the angle and the center
  double **cluster_approx_center;
  memory->create(cluster_approx_center, noccupied, 3, "ComputeFrenkel::cluster_approx_center");
  for (int n = 0; n < noccupied; n++) {
    if (cluster_nsites[n] == 0) continue;
    cluster_xi[n][0] /= cluster_nsites[n];
    cluster_xi[n][1] /= cluster_nsites[n];
    cluster_xi[n][2] /= cluster_nsites[n];
    cluster_zeta[n][0] /= cluster_nsites[n];
    cluster_zeta[n][1] /= cluster_nsites[n];
    cluster_zeta[n][2] /= cluster_nsites[n];
    if (domain->xperiodic) {
      theta = atan2(-cluster_zeta[n][0], -cluster_xi[n][0]) + MY_PI;
      cluster_approx_center[n][0] = domain->xprd * theta / MY_2PI + domain->boxlo[0];
    } else
      cluster_approx_center[n][0] = cluster_xi[n][0];
    if (domain->yperiodic) {
      theta = atan2(-cluster_zeta[n][1], -cluster_xi[n][1]) + MY_PI;
      cluster_approx_center[n][1] = domain->yprd * theta / MY_2PI + domain->boxlo[1];
    } else
      cluster_approx_center[n][1] = cluster_xi[n][1];
    if (domain->zperiodic) {
      theta = atan2(-cluster_zeta[n][2], -cluster_xi[n][2]) + MY_PI;
      cluster_approx_center[n][2] = domain->zprd * theta / MY_2PI + domain->boxlo[2];
    } else
      cluster_approx_center[n][2] = cluster_xi[n][2];

    domain->remap(cluster_approx_center[n]);
  }
  memory->destroy(cluster_zeta);
  memory->destroy(cluster_xi);

  // Now find the center using the approximate center to pick the right image
  double **local_cluster_x;
  double nearest_image[3];
  memory->create(local_cluster_x, noccupied, 3, "ComputeFrenkel:local_cluster_x");
  for (int i = 0; i < noccupied * 3; i++) local_cluster_x[0][i] = 0.0;
  memory->destroy(cluster_center);
  memory->create(cluster_center, noccupied, 3, "ComputeFrenkel::cluster_center");
  for (int i = 0; i < noccupied * 3; i++) cluster_center[0][i] = 0.0;
  for (int k = 0; k < nlatsites; k++) {
    if (clusterID[k] == 0) continue;
    n = clusterID2occupied_index(clusterID[k]);
    domain->closest_image(cluster_approx_center[n], latsites[k], nearest_image);
    local_cluster_x[n][0] += nearest_image[0] / cluster_nsites[n];
    local_cluster_x[n][1] += nearest_image[1] / cluster_nsites[n];
    local_cluster_x[n][2] += nearest_image[2] / cluster_nsites[n];
  }
  MPI_Allreduce(local_cluster_x[0], cluster_center[0], noccupied * 3, MPI_DOUBLE, MPI_SUM, world);

  memory->destroy(local_cluster_x);
  memory->destroy(cluster_approx_center);
}

/****************************************************************************/

void ComputeFrenkel::find_occupied_clusters()
{
  tagint *local_occupied_cluster_ID;
  int local_noccupied = 0;
  for (int k = 0; k < nlatsites; k++)
    if (clusterID[k] > 0) local_noccupied += 1;
  memory->create(local_occupied_cluster_ID, local_noccupied, "ComputeFrenkel:occupied_cluster_ID");

  // Assign occupied clusters to an array
  int n = 0;
  for (int k = 0; k < nlatsites; k++)
    if (clusterID[k] > 0) local_occupied_cluster_ID[n++] = clusterID[k];
  // sort, then remove duplicates (defects with more than one site involved);
  // all stored IDs are positive, so sort+unique yields the distinct cluster IDs
  std::sort(local_occupied_cluster_ID, local_occupied_cluster_ID + local_noccupied);
  local_noccupied = static_cast<int>(
      std::unique(local_occupied_cluster_ID, local_occupied_cluster_ID + local_noccupied) -
      local_occupied_cluster_ID);
  memory->grow(local_occupied_cluster_ID, local_noccupied,
               "ComputeFrenkel:reallocate_local_occupied_cluster_ID");

  // Now make the array global (size will be number of clusters across ALL
  // processes)
  MPI_Allreduce(&local_noccupied, &noccupied, 1, MPI_INT, MPI_SUM, world);
  memory->destroy(occupied_cluster_ID);
  memory->create(occupied_cluster_ID, noccupied, "ComputeFrenkel:occupied_cluster_ID");
  std::vector<int> nreceive(comm->nprocs, 0);
  MPI_Allgather(&local_noccupied, 1, MPI_INT, nreceive.data(), 1, MPI_INT, world);
  int displ;
  MPI_Scan(&local_noccupied, &displ, 1, MPI_INT, MPI_SUM, world);
  displ -= local_noccupied;
  std::vector<int> displacement(comm->nprocs, 0);
  MPI_Allgather(&displ, 1, MPI_INT, displacement.data(), 1, MPI_INT, world);
  MPI_Allgatherv(local_occupied_cluster_ID, local_noccupied, MPI_LMP_TAGINT, occupied_cluster_ID,
                 nreceive.data(), displacement.data(), MPI_LMP_TAGINT, world);
  memory->destroy(local_occupied_cluster_ID);

  // sort the gathered array and once again remove duplicates
  std::sort(occupied_cluster_ID, occupied_cluster_ID + noccupied);
  noccupied = static_cast<int>(std::unique(occupied_cluster_ID, occupied_cluster_ID + noccupied) -
                               occupied_cluster_ID);
  memory->grow(occupied_cluster_ID, noccupied, "ComputeFrenkel:reallocate_occupied_cluster_ID");
}

/****************************************************************************/

int ComputeFrenkel::clusterID2occupied_index(tagint id)
{
  for (int i = 0; i < noccupied; i++)
    if (id == occupied_cluster_ID[i]) return i;
  return -1;
}

/****************************************************************************/

void ComputeFrenkel::compute_vector()
{
  if (invoked_vector == update->ntimestep) return;
  invoked_vector = update->ntimestep;

  // If we haven't identified the defects yet, do so now
  if (invoked_find_defects != update->ntimestep) find_defects();

  // Zero out the array
  for (int n = 0; n < size_vector; n++) vector[n] = 0.0;

  // Add up all sites with zero (vacancies), two (interstitials), and more
  // than two occupants (irregular).
  int nints = 0, nvacs = 0, nirreg = 0;
  for (int k = 0; k < nlatsites; k++) {
    switch (noccupants[k]) {
      case 0:
        nvacs += 1;
        break;
      case 1:
        break;
      case 2:
        nints += 1;
        break;
      default:
        nints += 1;
        nirreg += 1;
    }
  }

  int nvacancies = 0, ninterstitials = 0, nirregular = 0;
  MPI_Allreduce(&nvacs, &nvacancies, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nints, &ninterstitials, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nirreg, &nirregular, 1, MPI_INT, MPI_SUM, world);

  vector[0] = static_cast<double>(nvacancies);
  vector[1] = static_cast<double>(ninterstitials);
  vector[2] = static_cast<double>(nirregular);
}

/****************************************************************************/

void ComputeFrenkel::compute_array()
{
  if (invoked_array == update->ntimestep) return;
  invoked_array = update->ntimestep;

  // If we haven't identified the defects or clusters yet, do so now
  if (invoked_find_defects != update->ntimestep) find_defects();
  if (invoked_find_clusters != update->ntimestep) find_clusters();

  // Initialize the array
  for (int i = 0; i < size_array_cols * size_array_rows; i++) array[0][i] = 0;

  // Add up # of -1's, -2's, etc. in the cluster size lists.
  for (int n = 0; n < noccupied; n++) {
    if (cluster_size[n] == 0) continue;
    if (cluster_size[n] > 0 && cluster_size[n] < size_array_cols)
      array[1][cluster_size[n] - 1] += 1;
    else if (cluster_size[n] >= size_array_cols)
      array[1][size_array_cols - 1] += 1;
    else if (cluster_size[n] < 0 && cluster_size[n] > -size_array_cols)
      array[0][-cluster_size[n] - 1] += 1;
    else
      array[0][size_array_cols - 1] += 1;
  }
}

/****************************************************************************/

void ComputeFrenkel::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // If we haven't identified the defects yet, do so now
  if (invoked_find_defects != update->ntimestep) find_defects();

  vector_atom = mindist;
}

/****************************************************************************/

// True if point c lies in this processor's subdomain.  Used to assign each
// defect cluster to exactly one rank (by its center) so global per-cluster
// output is neither duplicated nor dropped in parallel runs.
bool ComputeFrenkel::center_in_subdomain(const double *c) const
{
  return (c[0] >= domain->sublo[0] && c[0] < domain->subhi[0] && c[1] >= domain->sublo[1] &&
          c[1] < domain->subhi[1] && c[2] >= domain->sublo[2] && c[2] < domain->subhi[2]);
}

/****************************************************************************/

void ComputeFrenkel::compute_local()
{
  invoked_local = update->ntimestep;

  // If we haven't identified defects yet, do that, too
  if (invoked_find_defects != update->ntimestep) find_defects();

  // Ditto for clusters
  if (invoked_find_clusters != update->ntimestep) find_clusters();

  //  if ( noccupied == 0 ) return;

  // Find out whether the cluster belongs to this subdomain
  // (this prevents duplicate output)
  std::vector<bool> owned(noccupied);
  int nowned = 0;
  for (int i = 0; i < noccupied; i++) {
    owned[i] = center_in_subdomain(cluster_center[i]);
    if (owned[i]) nowned++;
  }

  size_local_rows = nowned;

  if (nowned == 0) { return; }

  memory->destroy(array_local);
  memory->create(array_local, size_local_rows, size_local_cols, "ComputeFrenkel:array_local");

  // Store the following quantities in the local array:
  // id, size, center_x, center_y, center_z
  int i = 0;
  for (int j = 0; j < noccupied; j++) {
    if (!owned[j]) continue;
    array_local[i][0] = static_cast<double>(occupied_cluster_ID[j]);
    array_local[i][1] = static_cast<double>(cluster_size[j]);
    /*    array_local[i][2] = (cluster_center[j][0] - domain->boxlo[0])/domain->xprd;
    array_local[i][3] = (cluster_center[j][1] - domain->boxlo[1])/domain->yprd;
    array_local[i][4] = (cluster_center[j][2] - domain->boxlo[2])/domain->zprd;
*/
    array_local[i][2] = cluster_center[j][0];
    array_local[i][3] = cluster_center[j][1];
    array_local[i][4] = cluster_center[j][2];
    i++;
  }
}

/****************************************************************************/

// Provide one sphere per defect cluster to "dump image ... compute <ID> ...".
// objarray columns are {color-index, x, y, z, diameter}; color index 1 marks
// vacancy clusters and 2 marks interstitial clusters.
int ComputeFrenkel::compute_image(int *&objs, double **&parms)
{
  if (invoked_image != update->ntimestep) {
    invoked_image = update->ntimestep;
    if (invoked_find_defects != update->ntimestep) find_defects();
    if (invoked_find_clusters != update->ntimestep) find_clusters();

    // only draw clusters whose center is owned by this subdomain, so each
    // defect is rendered exactly once when running in parallel
    int nowned = 0;
    for (int j = 0; j < noccupied; j++)
      if (center_in_subdomain(cluster_center[j])) nowned++;

    memory->destroy(image_objvec);
    memory->destroy(image_objarray);
    int alloc = (nowned > 0) ? nowned : 1;
    memory->create(image_objvec, alloc, "ComputeFrenkel:image_objvec");
    memory->create(image_objarray, alloc, 5, "ComputeFrenkel:image_objarray");

    double a = domain->lattice->xlattice;
    if (domain->lattice->ylattice < a) a = domain->lattice->ylattice;
    if (domain->lattice->zlattice < a) a = domain->lattice->zlattice;
    double diam = 0.6 * a;

    int n = 0;
    for (int j = 0; j < noccupied; j++) {
      double *c = cluster_center[j];
      if (!center_in_subdomain(c)) continue;
      image_objvec[n] = Graphics::SPHERE;
      image_objarray[n][0] = (cluster_size[j] < 0) ? 1.0 : 2.0;
      image_objarray[n][1] = c[0];
      image_objarray[n][2] = c[1];
      image_objarray[n][3] = c[2];
      image_objarray[n][4] = diam;
      n++;
    }
    image_nmax = n;
  }

  objs = image_objvec;
  parms = image_objarray;
  return image_nmax;
}

/****************************************************************************/

void ComputeFrenkel::create_lattice_sites()
{
  if (domain->lattice == nullptr)
    error->all(FLERR, Error::NOLASTLINE, "Use of compute style frenkel with undefined lattice");

  // Generate the reference lattice sites directly from the lattice definition,
  // keeping the ones whose position falls in this processor's subdomain.  This
  // is the same tiling create_atoms does, but a compute must not add atoms,
  // create groups, or run input commands, so it is done inline here.

  static constexpr double EPSILON = 1.0e-6;
  Lattice *const lattice = domain->lattice;
  const int triclinic = domain->triclinic;

  // my subdomain bounds (box coords if orthogonal, lamda coords if triclinic),
  // shrunk by EPSILON at the global periodic boundaries so that a site sitting
  // exactly on a periodic face is generated on exactly one side of the box.
  double sublo[3], subhi[3];
  for (int d = 0; d < 3; ++d) {
    sublo[d] = triclinic ? domain->sublo_lamda[d] : domain->sublo[d];
    subhi[d] = triclinic ? domain->subhi_lamda[d] : domain->subhi[d];
  }
  const int periodic[3] = {domain->xperiodic, domain->yperiodic, domain->zperiodic};
  for (int d = 0; d < 3; ++d) {
    if (!periodic[d]) continue;
    const double eps = triclinic ? EPSILON : domain->prd[d] * EPSILON;
    if (comm->layout != Comm::LAYOUT_TILED) {
      if (comm->myloc[d] == 0) sublo[d] -= eps;
      if (comm->myloc[d] == comm->procgrid[d] - 1) subhi[d] -= 2.0 * eps;
    } else {
      if (comm->mysplit[d][0] == 0.0) sublo[d] -= eps;
      if (comm->mysplit[d][1] == 1.0) subhi[d] -= 2.0 * eps;
    }
  }

  // true if box-coordinate point x lies in my (shrunk) subdomain
  auto in_subdomain = [&](double *x) {
    double lamda[3];
    const double *c = x;
    if (triclinic) {
      domain->x2lamda(x, lamda);
      c = lamda;
    }
    return (c[0] >= sublo[0] && c[0] < subhi[0] && c[1] >= sublo[1] && c[1] < subhi[1] &&
            c[2] >= sublo[2] && c[2] < subhi[2]);
  };

  std::vector<std::array<double, 3>> sites;

  if (!sitefile.empty()) {
    // read explicit site coordinates, one "x y z" per line
    try {
      TextFileReader reader(sitefile, "compute frenkel site");
      char *line;
      while ((line = reader.next_line(3))) {
        ValueTokenizer values(line);
        double x[3];
        x[0] = values.next_double();
        x[1] = values.next_double();
        x[2] = values.next_double();
        if (region && !region->match(x[0], x[1], x[2])) continue;
        if (in_subdomain(x)) sites.push_back({x[0], x[1], x[2]});
      }
    } catch (std::exception &e) {
      error->one(FLERR, Error::NOLASTLINE, e.what());
    }
  } else {
    // lattice-index bounds covering my subdomain (use the unshrunk subbox)
    double bboxlo[3], bboxhi[3];
    if (triclinic) {
      domain->bbox(domain->sublo_lamda, domain->subhi_lamda, bboxlo, bboxhi);
    } else {
      for (int d = 0; d < 3; ++d) {
        bboxlo[d] = domain->sublo[d];
        bboxhi[d] = domain->subhi[d];
      }
    }
    double lo[3] = {BIG, BIG, BIG}, hi[3] = {-BIG, -BIG, -BIG};
    for (int c = 0; c < 8; ++c) {
      double cx = (c & 1) ? bboxhi[0] : bboxlo[0];
      double cy = (c & 2) ? bboxhi[1] : bboxlo[1];
      double cz = (c & 4) ? bboxhi[2] : bboxlo[2];
      lattice->bbox(1, cx, cy, cz, lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]);
    }
    int ilo = static_cast<int>(lo[0]) - 1;
    int jlo = static_cast<int>(lo[1]) - 1;
    int klo = static_cast<int>(lo[2]) - 1;
    int ihi = static_cast<int>(hi[0]) + 1;
    int jhi = static_cast<int>(hi[1]) + 1;
    int khi = static_cast<int>(hi[2]) + 1;
    if (lo[0] < 0.0) --ilo;
    if (lo[1] < 0.0) --jlo;
    if (lo[2] < 0.0) --klo;
    if (domain->dimension == 2) klo = khi = 0;

    const double *const *const basis = lattice->basis;
    const int nbasis = lattice->nbasis;
    for (int k = klo; k <= khi; ++k)
      for (int j = jlo; j <= jhi; ++j)
        for (int i = ilo; i <= ihi; ++i)
          for (int m = 0; m < nbasis; ++m) {
            double x[3] = {i + basis[m][0], j + basis[m][1], k + basis[m][2]};
            lattice->lattice2box(x[0], x[1], x[2]);
            if (region && !region->match(x[0], x[1], x[2])) continue;
            if (in_subdomain(x)) sites.push_back({x[0], x[1], x[2]});
          }
  }

  nlatsites = static_cast<int>(sites.size());
  nlatghosts = 0;
  if (nlatsites == 0) error->warning(FLERR, "compute frenkel generated no lattice sites");

  memory->destroy(latsites);
  memory->create(latsites, MAX(nlatsites, 1), 3, "ComputeFrenkel:latsites");
  for (int n = 0; n < nlatsites; ++n) {
    latsites[n][0] = sites[n][0];
    latsites[n][1] = sites[n][1];
    latsites[n][2] = sites[n][2];
  }

  // give every site a globally unique tag, numbered above the real atom tags
  memory->destroy(site_tag);
  memory->create(site_tag, MAX(nlatsites, 1), "ComputeFrenkel:site_tag");
  tagint nlats = nlatsites;
  first_local_tag = 0;
  MPI_Scan(&nlats, &first_local_tag, 1, MPI_LMP_TAGINT, MPI_SUM, world);
  first_local_tag = first_local_tag - nlatsites + 1 + atom->natoms;
  for (int n = 0; n < nlatsites; ++n) site_tag[n] = first_local_tag + n;

  // keep a pristine copy of the coordinates when rescaling is enabled
  memory->destroy(latsites0);
  if (rescale) {
    memory->create(latsites0, MAX(nlatsites, 1), 3, "ComputeFrenkel:latsites0");
    memcpy(*latsites0, *latsites, nlatsites * 3 * sizeof(double));
  } else
    latsites0 = nullptr;
}

/****************************************************************************/

int ComputeFrenkel::site_tag2index(tagint tag)
{
  int idx = static_cast<int>(tag - this->first_local_tag);
  if (idx >= 0 && idx < nlatsites) return idx;

  // Tag belongs only to a ghost; look it up in the O(1) ghost index
  auto it = ghost_index.find(tag);
  if (it != ghost_index.end()) return it->second;

  // Hmm...this tag doesn't seem to exist!  Why are you looking for it?
  error->warning(FLERR, "(proc {}): Didn't find an index for tag {}\n", comm->me, tag);
  return -1;
}

/****************************************************************************/

// Post nonblocking sends of sbuf(p) and matching receives into rbuf(p) for all
// p != me (message length in datatype units = n_send[p]*stride on send,
// n_recv[p]*stride on receive), then wait for all to complete.  The caller
// handles the self (p == me) copy.  One fixed tag per message type keeps the
// tags well below MPI_TAG_UB; the source rank already disambiguates messages.
void ComputeFrenkel::exchange_one(const std::function<void *(int)> &sbuf,
                                  const std::function<void *(int)> &rbuf,
                                  const std::vector<int> &n_send, const std::vector<int> &n_recv,
                                  int stride, MPI_Datatype type, int tag)
{
  const int me = comm->me;
  const int nprocs = comm->nprocs;
  if (nprocs == 1) return;    // nothing to exchange in serial

  std::vector<MPI_Request> sreq(nprocs, MPI_REQUEST_NULL);
  std::vector<MPI_Request> rreq(nprocs, MPI_REQUEST_NULL);
  for (int p = 0; p < nprocs; p++) {
    if (p == me) continue;
    if (n_send[p] > 0) MPI_Isend(sbuf(p), n_send[p] * stride, type, p, tag, world, &sreq[p]);
    if (n_recv[p] > 0) MPI_Irecv(rbuf(p), n_recv[p] * stride, type, p, tag, world, &rreq[p]);
  }
  MPI_Waitall(nprocs, sreq.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(nprocs, rreq.data(), MPI_STATUSES_IGNORE);
}

/****************************************************************************/

void ComputeFrenkel::exchange_lattice_ghosts()
{
  // Exchange sites within one cutoff distance of the edge of the box with
  // adjacent processes.  Quantities exchanged:  site tag, site coordinates,
  // cluster ID's, and site occupancies.

  const int me = comm->me;
  const int nprocs = comm->nprocs;
  tagint **tag_send;
  tagint **tag_recv;
  tagint **clusterID_send;
  tagint **clusterID_recv;
  double ***x_send;
  double ***x_recv;
  int **occup_send;
  int **occup_recv;
  const int NINCR = 100;
  int max_size = NINCR;

  nlatghosts = 0;
  ghost_index.clear();
  std::vector<int> n_send(nprocs, 0);
  std::vector<int> n_recv(nprocs, 0);
  std::vector<int> idx(nprocs, 0);
  memory->create(tag_send, nprocs, max_size, "ComputeFrenkel:tag_send");
  memory->create(x_send, nprocs, max_size, 3, "ComputeFrenkel:x_send");
  memory->create(occup_send, nprocs, max_size, "ComputeFrenkel:occup_send");
  memory->create(clusterID_send, nprocs, max_size, "ComputeFrenkel:clusterID_send");
  for (int i = 0; i < nprocs * max_size; i++) tag_send[0][i] = -1;
  for (int i = 0; i < nprocs * max_size * 3; i++) x_send[0][0][i] = 0.0;
  for (int i = 0; i < nprocs * max_size; i++) occup_send[0][i] = 0;
  for (int i = 0; i < nprocs * max_size; i++) clusterID_send[0][i] = 0;

  // Find out which sites are near the boundaries
  // and which processes those boundaries correspond to.
  double dr[2][3];
  std::vector<bool> already_sent(nprocs);
  for (int k = 0; k < nlatsites; k++) {
    dr[0][0] = latsites[k][0] - domain->sublo[0];
    dr[0][1] = latsites[k][1] - domain->sublo[1];
    dr[0][2] = latsites[k][2] - domain->sublo[2];
    dr[1][0] = domain->subhi[0] - latsites[k][0];
    dr[1][1] = domain->subhi[1] - latsites[k][1];
    dr[1][2] = domain->subhi[2] - latsites[k][2];

    // Prepare list of processes to which we've sent this site to (to avoid
    // duplication)
    for (int p = 0; p < nprocs; p++) already_sent[p] = false;

    // Test processes in all 27 directions from this one
    for (int ix = -1; ix <= 1; ix++) {
      // If at the end of the domain, ignore those directions unless
      // box is periodic in that direction
      if (ix == -1 && comm->myloc[0] == 0 && !domain->xperiodic) continue;
      if (ix == 1 && comm->myloc[0] == comm->procgrid[0] - 1 && !domain->xperiodic) continue;
      if (ix == -1 && dr[0][0] > cutoff + SMALL) continue;
      if (ix == 1 && dr[1][0] > cutoff + SMALL) continue;
      for (int iy = -1; iy <= 1; iy++) {
        if (iy == -1 && comm->myloc[1] == 0 && !domain->yperiodic) continue;
        if (iy == 1 && comm->myloc[1] == comm->procgrid[1] - 1 && !domain->yperiodic) continue;
        if (iy == -1 && dr[0][1] > cutoff + SMALL) continue;
        if (iy == 1 && dr[1][1] > cutoff + SMALL) continue;
        for (int iz = -1; iz <= 1; iz++) {
          if (ix == 0 && iy == 0 && iz == 0) continue;
          if (iz == -1 && comm->myloc[2] == 0 && !domain->zperiodic) continue;
          if (iz == 1 && comm->myloc[2] == comm->procgrid[2] - 1 && !domain->zperiodic) continue;
          if (iz == -1 && dr[0][2] > cutoff + SMALL) continue;
          if (iz == 1 && dr[1][2] > cutoff + SMALL) continue;
          // Should only get here when exchanging with an adjacent process
          int p = process_neighbor(ix, iy, iz);
          // Don't add site to process multiple times
          if (already_sent[p]) continue;
          n_send[p] += 1;
          if (n_send[p] > max_size) {    // Grow arrays if necessary
            this->reallocate_array(tag_send, nprocs, max_size, nprocs, max_size + NINCR);
            this->reallocate_array(x_send, nprocs, max_size, 3, nprocs, max_size + NINCR, 3);
            this->reallocate_array(occup_send, nprocs, max_size, nprocs, max_size + NINCR);
            this->reallocate_array(clusterID_send, nprocs, max_size, nprocs, max_size + NINCR);
            max_size += NINCR;
          }
          tag_send[p][idx[p]] = site_tag[k];
          occup_send[p][idx[p]] = noccupants[k];
          if (clusterID) clusterID_send[p][idx[p]] = clusterID[k];
          x_send[p][idx[p]][0] = latsites[k][0];
          x_send[p][idx[p]][1] = latsites[k][1];
          x_send[p][idx[p]][2] = latsites[k][2];
          idx[p] += 1;
          already_sent[p] = true;
        }
      }
    }
  }
  // Send your sites to them and get theirs in return
  std::vector<MPI_Request> send_request(nprocs);
  std::vector<MPI_Request> recv_request(nprocs);
  send_request[me] = MPI_REQUEST_NULL;
  recv_request[me] = MPI_REQUEST_NULL;
  for (int p = 0; p < nprocs; p++) {
    if (p == me)
      n_recv[p] = n_send[p];
    else {
      MPI_Isend(&n_send[p], 1, MPI_INT, p, TAG_COUNT, world, &send_request[p]);
      MPI_Irecv(&n_recv[p], 1, MPI_INT, p, TAG_COUNT, world, &recv_request[p]);
    }
  }
  if (nprocs > 1) {
    MPI_Waitall(nprocs, send_request.data(), MPI_STATUSES_IGNORE);
    MPI_Waitall(nprocs, recv_request.data(), MPI_STATUSES_IGNORE);
  }

  // Remember the per-proc message sizes so push_ghost_labels_to_owners() can
  // run the reverse (ghost -> owner) exchange.  Record them now, before the
  // "no ghosts received" early return below: a rank may receive no ghosts of
  // its own yet still own sites that are ghosts elsewhere, and those owners
  // will reverse-send label updates back to it.
  ghost_send_counts = n_send;
  ghost_recv_counts = n_recv;

  max_size = 0;
  for (int p = 0; p < nprocs; p++) max_size = MAX(max_size, n_recv[p]);
  if (max_size == 0) {
    memory->destroy(tag_send);
    memory->destroy(x_send);
    memory->destroy(occup_send);
    memory->destroy(clusterID_send);
    return;    // Should only happen if regions don't cross domain boundaries
  }
  memory->create(tag_recv, nprocs, max_size, "ComputeFrenkel:tag_recv");
  memory->create(x_recv, nprocs, max_size, 3, "ComputeFrenkel:x_recv");
  memory->create(occup_recv, nprocs, max_size, "ComputeFrenkel:occup_recv");
  memory->create(clusterID_recv, nprocs, max_size, "ComputeFrenkel:clusterID_recv");
  for (int i = 0; i < nprocs * max_size; i++) tag_recv[0][i] = 0;
  for (int i = 0; i < nprocs * max_size; i++) occup_recv[0][i] = 0;
  for (int i = 0; i < nprocs * max_size * 3; i++) x_recv[0][0][i] = 0.0;
  for (int i = 0; i < nprocs * max_size; i++) clusterID_recv[0][i] = 0;
  // self (p == me) is just a local copy; the four per-site quantities are then
  // exchanged with the neighbor processes one message type at a time.  A site is
  // sent to itself only as a periodic self-image (procgrid == 1 in a periodic
  // dimension); the ghost occupancy is a verbatim copy of the owner's, exactly
  // like the cross-proc path -- it is only used to classify the ghost as
  // vacancy/normal/interstitial during clustering and never re-counted (real
  // occupancy is tallied on the owning rank in find_defects).  Accumulating
  // (+=) would be equivalent anyway: already_sent[p] sends each owned site to a
  // given proc at most once, so every self slot is written exactly once.
  for (int i = 0; i < n_send[me]; i++) {
    tag_recv[me][i] = tag_send[me][i];
    occup_recv[me][i] = occup_send[me][i];
    clusterID_recv[me][i] = clusterID_send[me][i];
    x_recv[me][i][0] = x_send[me][i][0];
    x_recv[me][i][1] = x_send[me][i][1];
    x_recv[me][i][2] = x_send[me][i][2];
  }
  exchange_one(
      [&](int p) {
        return (void *) tag_send[p];
      },
      [&](int p) {
        return (void *) tag_recv[p];
      },
      n_send, n_recv, 1, MPI_LMP_TAGINT, TAG_SITE);
  exchange_one(
      [&](int p) {
        return (void *) x_send[p][0];
      },
      [&](int p) {
        return (void *) x_recv[p][0];
      },
      n_send, n_recv, 3, MPI_DOUBLE, TAG_X);
  exchange_one(
      [&](int p) {
        return (void *) occup_send[p];
      },
      [&](int p) {
        return (void *) occup_recv[p];
      },
      n_send, n_recv, 1, MPI_INT, TAG_OCCUP);
  exchange_one(
      [&](int p) {
        return (void *) clusterID_send[p];
      },
      [&](int p) {
        return (void *) clusterID_recv[p];
      },
      n_send, n_recv, 1, MPI_LMP_TAGINT, TAG_CLUST);

  // Reallocate the necessary memory
  nlatghosts = 0;
  for (int p = 0; p < nprocs; p++) nlatghosts += n_recv[p];
  memory->grow(latsites, nlatsites + nlatghosts, 3, "ComputeFrenkel:sites2");
  memory->grow(site_tag, nlatsites + nlatghosts, "ComputeFrenkel:site_tag2");
  memory->grow(noccupants, nlatsites + nlatghosts, "ComputeFrenkel:noccupants2");
  if (clusterID) memory->grow(clusterID, nlatsites + nlatghosts, "ComputeFrenkel:clusterID2");

  // Update latsites, site_tag, clusterID, and nlatghosts
  int kk = nlatsites;
  for (int p = 0; p < nprocs; p++) {
    for (int i = 0; i < n_recv[p]; i++) {
      site_tag[kk] = tag_recv[p][i];
      latsites[kk][0] = x_recv[p][i][0];
      latsites[kk][1] = x_recv[p][i][1];
      latsites[kk][2] = x_recv[p][i][2];
      noccupants[kk] = occup_recv[p][i];
      if (clusterID) clusterID[kk] = clusterID_recv[p][i];
      ghost_index[site_tag[kk]] = kk;
      kk++;
    }
  }

  // Clean up
  memory->destroy(tag_send);
  memory->destroy(x_send);
  memory->destroy(occup_send);
  memory->destroy(clusterID_send);
  memory->destroy(tag_recv);
  memory->destroy(x_recv);
  memory->destroy(occup_recv);
  memory->destroy(clusterID_recv);
}

/****************************************************************************/

// Reverse of the cluster-ID part of exchange_lattice_ghosts(): every ghost site
// sends the (possibly lowered) cluster ID it now holds back to the rank that
// owns it, which keeps the smaller of its stored value and the incoming one.
// This is the ghost -> owner half of the union-find label propagation: without
// it a relabelling discovered on the receiver side of a subdomain boundary
// would be discarded by the next exchange instead of reaching the owner.
// Returns nonzero on this rank if any owned label was lowered.
int ComputeFrenkel::push_ghost_labels_to_owners()
{
  const int me = comm->me;
  const int nprocs = comm->nprocs;
  int changed = 0;

  // Slice the ghost range back into the per-proc groups it arrived in: the same
  // ascending proc order used when the ghosts were appended, ghost_recv_counts[p]
  // entries each (the p == me group carries periodic self-images).
  std::vector<std::vector<tagint>> tag_back(nprocs), cid_back(nprocs);
  int kk = nlatsites;
  for (int p = 0; p < nprocs; p++) {
    tag_back[p].reserve(ghost_recv_counts[p]);
    cid_back[p].reserve(ghost_recv_counts[p]);
    for (int i = 0; i < ghost_recv_counts[p]; i++, kk++) {
      tag_back[p].push_back(site_tag[kk]);
      cid_back[p].push_back(clusterID[kk]);
    }
  }

  // self group: this rank is the owner, so apply the MIN locally
  for (std::size_t i = 0; i < tag_back[me].size(); i++) {
    int idx = site_tag2index(tag_back[me][i]);
    if (idx >= 0 && idx < nlatsites && cid_back[me][i] < clusterID[idx]) {
      clusterID[idx] = cid_back[me][i];
      changed = 1;
    }
  }
  if (nprocs == 1) return changed;

  // Reverse exchange: send each held ghost's label back to its owner.  The send
  // counts are ghost_recv_counts (ghosts we hold from p) and the receive counts
  // are ghost_send_counts (our owned sites mirrored onto p, now coming back).
  std::vector<std::vector<tagint>> rtag(nprocs), rcid(nprocs);
  for (int p = 0; p < nprocs; p++) {
    rtag[p].assign(ghost_send_counts[p], 0);
    rcid[p].assign(ghost_send_counts[p], 0);
  }
  exchange_one(
      [&](int p) {
        return (void *) tag_back[p].data();
      },
      [&](int p) {
        return (void *) rtag[p].data();
      },
      ghost_recv_counts, ghost_send_counts, 1, MPI_LMP_TAGINT, TAG_RTAG);
  exchange_one(
      [&](int p) {
        return (void *) cid_back[p].data();
      },
      [&](int p) {
        return (void *) rcid[p].data();
      },
      ghost_recv_counts, ghost_send_counts, 1, MPI_LMP_TAGINT, TAG_RCID);

  // owner side: keep the smaller of the stored and incoming cluster IDs
  for (int p = 0; p < nprocs; p++) {
    if (p == me) continue;
    for (int i = 0; i < ghost_send_counts[p]; i++) {
      int idx = site_tag2index(rtag[p][i]);
      if (idx >= 0 && idx < nlatsites && rcid[p][i] < clusterID[idx]) {
        clusterID[idx] = rcid[p][i];
        changed = 1;
      }
    }
  }
  return changed;
}

/****************************************************************************/

void ComputeFrenkel::construct_site_nlists()
{
  // Builds neighbor lists for every site in the subdomain

  double dx, dy, dz, drsq;
  double cutsq = cutoff * cutoff;
  int ix, iy, iz, i, j, k;
  nlist.assign(nlatsites + nlatghosts, {});
  for (int n = 0; n < nlatsites + nlatghosts; n++) {
    // Loop over all sites in this site's bin and those adjacent
    find_closest_bin(latsites[n], i, j, k);
    for (int ii = i - 1; ii <= i + 1; ii++) {
      ix = ii;
      for (int jj = j - 1; jj <= j + 1; jj++) {
        iy = jj;
        for (int kk = k - 1; kk <= k + 1; kk++) {
          iz = kk;
          bin_pbc(ix, iy, iz);
          for (int l = 0; l < nlatbins[3]; l++) {
            int m = latbins[ix][iy][iz][l];
            if (m < 0) break;
            if (m == n) continue;
            dx = latsites[n][0] - latsites[m][0];
            dy = latsites[n][1] - latsites[m][1];
            dz = latsites[n][2] - latsites[m][2];
            domain->minimum_image(FLERR, dx, dy, dz);
            drsq = dx * dx + dy * dy + dz * dz;
            if (drsq <= cutsq) {
              // Add site m to neighbor list of site n
              nlist[n].push_back(site_tag[m]);
            }
          }
        }
      }
    }
    // Sort neighbor list and remove duplicates
    std::sort(nlist[n].begin(), nlist[n].end());
    nlist[n].erase(std::unique(nlist[n].begin(), nlist[n].end()), nlist[n].end());
  }
}

/****************************************************************************/

void ComputeFrenkel::put_sites_in_bins()
{
  // Each site is assigned to a unique bin of width cutoff.
  nlatbins[0] = ceil((domain->subhi[0] - domain->sublo[0]) / binwidth);
  nlatbins[1] = ceil((domain->subhi[1] - domain->sublo[1]) / binwidth);
  nlatbins[2] = ceil((domain->subhi[2] - domain->sublo[2]) / binwidth);
  nlatbins[3] = domain->lattice->nbasis * BIN_GROW_SIZE;    // can grow

  // Temporary array to store the index we're currently on
  int ***binindex =
      memory->create(binindex, nlatbins[0], nlatbins[1], nlatbins[2], "ComputeFrenkel:binindex");

  // Create the bins
  memory->destroy(latbins);
  memory->create(latbins, nlatbins[0], nlatbins[1], nlatbins[2], nlatbins[3],
                 "ComputeFrenkel:nlatbins");

  // Initialize the bin contents
  int *iptr = binindex[0][0];
  for (int i = 0; i < nlatbins[0] * nlatbins[1] * nlatbins[2]; i++) iptr[i] = 0;
  iptr = latbins[0][0][0];
  for (int i = 0; i < nlatbins[0] * nlatbins[1] * nlatbins[2] * nlatbins[3]; i++) iptr[i] = -1;

  // Store indices of all lattice points in appropriate bins
  for (int n = 0; n < nlatsites + nlatghosts; n++) {
    int i, j, k, l;
    find_closest_bin(latsites[n], i, j, k);
    l = binindex[i][j][k];
    latbins[i][j][k][l] = n;
    binindex[i][j][k] += 1;
    if (binindex[i][j][k] >= nlatbins[3]) {
      int startval = nlatbins[3];
      nlatbins[3] += BIN_GROW_SIZE;
      memory->grow(latbins, nlatbins[0], nlatbins[1], nlatbins[2], nlatbins[3],
                   "ComputeFrenkel:realloc-latbins");
      for (int ii = 0; ii < nlatbins[0]; ii++)
        for (int jj = 0; jj < nlatbins[1]; jj++)
          for (int kk = 0; kk < nlatbins[2]; kk++)
            for (int ll = startval; ll < nlatbins[3]; ll++) latbins[ii][jj][kk][ll] = -1;
    }
  }

  // Deallocate our counting array
  memory->destroy(binindex);
}

/****************************************************************************/

void ComputeFrenkel::find_closest_bin(double *r, int &i, int &j, int &k)
{
  const int me = comm->me;

  double x = r[0] - domain->sublo[0];
  double y = r[1] - domain->sublo[1];
  double z = r[2] - domain->sublo[2];
  i = std::lround(x / binwidth);
  j = std::lround(y / binwidth);
  k = std::lround(z / binwidth);

  // Apply periodic boundaries, if appropriate
  if (domain->xperiodic && comm->procneigh[0][0] == me)
    i = (nlatbins[0] + i) % nlatbins[0];
  else {
    i = MIN(i, nlatbins[0] - 1);
    i = MAX(i, 0);
  }
  if (domain->yperiodic && comm->procneigh[1][0] == me)
    j = (nlatbins[1] + j) % nlatbins[1];
  else {
    j = MIN(j, nlatbins[1] - 1);
    j = MAX(j, 0);
  }
  if (domain->zperiodic && comm->procneigh[2][0] == me)
    k = (nlatbins[2] + k) % nlatbins[2];
  else {
    k = MIN(k, nlatbins[2] - 1);
    k = MAX(k, 0);
  }
  // (i,j,k) is now the bin on this process whose center is closest to r
}

/****************************************************************************/

void ComputeFrenkel::bin_pbc(int &i, int &j, int &k)
{
  const int me = comm->me;

  if (i < 0) {
    if (domain->xperiodic && comm->procneigh[0][0] == me)
      i = nlatbins[0] - 1;
    else
      i = 0;
  } else if (i >= nlatbins[0]) {
    if (domain->xperiodic && comm->procneigh[0][1] == me)
      i = 0;
    else
      i = nlatbins[0] - 1;
  }
  if (j < 0) {
    if (domain->yperiodic && comm->procneigh[1][0] == me)
      j = nlatbins[1] - 1;
    else
      j = 0;
  } else if (j >= nlatbins[1]) {
    if (domain->yperiodic && comm->procneigh[1][1] == me)
      j = 0;
    else
      j = nlatbins[1] - 1;
  }
  if (k < 0) {
    if (domain->zperiodic && comm->procneigh[2][0] == me)
      k = nlatbins[2] - 1;
    else
      k = 0;
  } else if (k >= nlatbins[2]) {
    if (domain->zperiodic && comm->procneigh[2][1] == me)
      k = 0;
    else
      k = nlatbins[2] - 1;
  }
}

/****************************************************************************/

void ComputeFrenkel::construct_WS_cell()
{
  invoked_construct_WS_cell = update->ntimestep;

  // Constructs the Wigner-Seitz cell (in lattice units) from the basis
  // vectors and lattice directions.

  nnormal = 27 * domain->lattice->nbasis - 1;
  memory->destroy(normal);
  memory->create(normal, nnormal, 3, "ComputeFrenkel:normal");
  double **basis = domain->lattice->basis;
  double *a1 = domain->lattice->a1;
  double *a2 = domain->lattice->a2;
  double *a3 = domain->lattice->a3;

  for (int i = 0; i < nnormal * 3; i++) normal[0][i] = 0.0;

  int n = 0;
  for (int j = 0; j < domain->lattice->nbasis; j++)
    for (int n1 = -1; n1 <= 1; n1++)
      for (int n2 = -1; n2 <= 1; n2++)
        for (int n3 = -1; n3 <= 1; n3++) {
          // Skip basis vector zero; that's the origin!
          if (j == 0 && n1 == 0 && n2 == 0 && n3 == 0) continue;
          // The normal vector is a vector passing through lattice point
          // zero (arbitrary) and another lattice point either in this unit
          // cell or in an adjacent one.  The factor of 0.5 means it is also
          // the midpoint of that line segment.
          normal[n][0] = 0.5 * (basis[j][0] + n1 * a1[0] + n2 * a2[0] + n3 * a3[0] - basis[0][0]);
          normal[n][1] = 0.5 * (basis[j][1] + n1 * a1[1] + n2 * a2[1] + n3 * a3[1] - basis[0][1]);
          normal[n][2] = 0.5 * (basis[j][2] + n1 * a1[2] + n2 * a2[2] + n3 * a3[2] - basis[0][2]);
          n = n + 1;
        }
}

/****************************************************************************/

bool ComputeFrenkel::inside_WS_cell(int n, int k)
{    // Returns true if atom n is inside Wigner-Seitz cell k, false if not
  // This check is not necessary for atoms not near the edge of a process
  // domain; we COULD check that condition and speed this up.  However, we
  // would also have to check against (arbitrary) user-defined regions, so
  // that isn't as trivial as it sounds.

  // First, adjust the coordinates for periodic boundary conditions
  double r[3];
  domain->closest_image(latsites[k], atom->x[n], r);

  // Easy check:  if it's more than half a lattice unit from the center of the
  // Wigner-Seitz cell, it's outside the cell!
  if (fabs(r[0] - latsites[k][0]) > 0.5 * domain->lattice->xlattice) return false;
  if (fabs(r[1] - latsites[k][1]) > 0.5 * domain->lattice->ylattice) return false;
  if (fabs(r[2] - latsites[k][2]) > 0.5 * domain->lattice->zlattice) return false;

  // If we haven't constructed the W-S cell, do it now
  if (invoked_construct_WS_cell != update->ntimestep) construct_WS_cell();

  // It's inside the unit cell drawn with this site at the center, so we
  // now determine whether it's inside the Wigner-Seitz cell for this
  // lattice, centered at this particular lattice site.
  double x, y, z, sx, sy, sz;
  x = r[0];
  y = r[1];
  z = r[2];
  domain->lattice->box2lattice(x, y, z);
  sx = latsites[k][0];
  sy = latsites[k][1];
  sz = latsites[k][2];
  domain->lattice->box2lattice(sx, sy, sz);
  x = x - sx;
  y = y - sy;
  z = z - sz;
  // The coordinates (x, y, z) are now in lattice units with an origin at the
  // center of this Wigner-Seitz cell.

  // The atom is inside the Wigner-Seitz cell if
  // normal[j] . (x,y,z) / ||normal[j]||**2 <= 1 for all j.

  // That is, find whether the projection of (x,y,z) onto each normal vector
  // extends beyond the normal vector itself.
  // I use the equivalent relation here,
  // normal[j] . (x,y,z) - ||normal[j]||**2 <= 0 for all j.
  double norm2sq, diff;
  for (int i = 0; i < nnormal; i++) {
    // Skip zero-length normal vectors (this should never happen)
    if (fabs(normal[i][0]) < SMALL && fabs(normal[i][1]) < SMALL && fabs(normal[i][2]) < SMALL)
      continue;
    norm2sq =
        normal[i][0] * normal[i][0] + normal[i][1] * normal[i][1] + normal[i][2] * normal[i][2];
    diff = x * normal[i][0] + y * normal[i][1] + z * normal[i][2] - norm2sq;
    if (diff > 0.0) return false;
  }

  return true;
}

/****************************************************************************/

bool ComputeFrenkel::tag_is_already_in_occupancy_list(tagint tag, int site)
{

  for (int i = 0; (i < MAX_OCCUPANTS) && (occupant_tag[site][i] > 0); i++)
    if (occupant_tag[site][i] == tag) return true;
  return false;
}

/****************************************************************************/

int ComputeFrenkel::next_free_occupant_tag_index(int s, int linenum)
{
  // Returns the index of the next available occupant tag for site s
  int i;
  for (i = 0; i < MAX_OCCUPANTS && occupant_tag[s][i] > 0;
       i++);    // that is, just increment i until we reach the next available index
  if (i >= MAX_OCCUPANTS)
    error->one(__FILE__, linenum, Error::NOLASTLINE,
               "Greater than {} atoms near a site; incorrect lattice, perhaps?", MAX_OCCUPANTS);
  return i;
}

/****************************************************************************/

void ComputeFrenkel::rescale_lattice_sites()
{
  if (!rescale) return;

  double xprd0, yprd0, zprd0;
  xprd0 = (old_boxhi[0] - old_boxlo[0]);
  yprd0 = (old_boxhi[1] - old_boxlo[1]);
  zprd0 = (old_boxhi[2] - old_boxlo[2]);

  double xcontract, ycontract, zcontract;
  xcontract = domain->xprd / xprd0;
  ycontract = domain->yprd / yprd0;
  zcontract = domain->zprd / zprd0;

  for (int k = 0; k < nlatsites; k++) {
    latsites[k][0] = (latsites0[k][0] - old_boxlo[0]) * xcontract + domain->boxlo[0];
    latsites[k][1] = (latsites0[k][1] - old_boxlo[1]) * ycontract + domain->boxlo[1];
    latsites[k][2] = (latsites0[k][2] - old_boxlo[2]) * zcontract + domain->boxlo[2];
  }
}

/****************************************************************************/

template <typename TYPE>
void ComputeFrenkel::reallocate_array(TYPE **&array, int x1, int y1, int x2, int y2)
{
  // Reallocates array[x1][y1] to array[x2][y2].  Both x2 and y2 can differ
  // from x1 and y1.  Note that if y1 == y2, you should use memory->grow.
  TYPE **arr2 = array;    // Now we just need a dee2
  memory->create(array, x2, y2, "ComputeFrenkel:reallocate2");
  for (int i = 0; i < x2; i++)
    for (int j = 0; j < y2; j++)
      if (i < x1 && j < y1)
        array[i][j] = arr2[i][j];
      else
        array[i][j] = TYPE();
  memory->destroy(arr2);
}

/****************************************************************************/

template <typename TYPE>
void ComputeFrenkel::reallocate_array(TYPE ***&array, int x1, int y1, int z1, int x2, int y2,
                                      int z2)
{

  // Reallocates array[x1][y1][z1] to array[x2][y2][z2].  x2, y2, and z2 can
  // differ from x1, y1, and z1.  Note that if y1 == y2 and z1 == z2, you
  // should use memory->grow.
  TYPE ***arr2 = array;    // Now we just need a dee2
  memory->create(array, x2, y2, z2, "ComputeFrenkel:reallocate3");
  for (int i = 0; i < x2; i++)
    for (int j = 0; j < y2; j++)
      for (int k = 0; k < z2; k++)
        if (i < x1 && j < y1 && k < z1)
          array[i][j][k] = arr2[i][j][k];
        else
          array[i][j][k] = TYPE();
  memory->destroy(arr2);
}

/****************************************************************************/

int ComputeFrenkel::process_neighbor(int x, int y, int z)
{
  // Returns the rank of the process with relative coordinates (x,y,z).
  // This process is (0,0,0).
  int rank = MPI_PROC_NULL;
  if (x == 0 && y == 0 && z == 0) return comm->me;

  int a = comm->myloc[0], b = comm->myloc[1], c = comm->myloc[2];
  if (domain->xperiodic) a = (a + comm->procgrid[0] + x) % comm->procgrid[0];
  if (domain->yperiodic) b = (b + comm->procgrid[1] + y) % comm->procgrid[1];
  if (domain->zperiodic) c = (c + comm->procgrid[2] + z) % comm->procgrid[2];
  if (a < 0 || b < 0 || c < 0) {
    rank = MPI_PROC_NULL;
    error->warning(FLERR, "Domain is inconsistent (got MPI_PROC_NULL next door)");
  } else
    rank = comm->grid2proc[a][b][c];
  return rank;
}

/****************************************************************************/

double ComputeFrenkel::memory_usage()
{
  double nbytes = 0.0;
  if (nlatsites == 0) return nbytes;
  nbytes += (atom->nmax) * sizeof(mindist[0]);
  nbytes += nlatsites * sizeof(site_mindist[0]);
  nbytes += (nlatsites + nlatghosts) * sizeof(noccupants[0]);
  nbytes += nlatsites * MAX_OCCUPANTS * sizeof(occupant_tag[0][0]);
  nbytes += nnormal * 3 * sizeof(normal[0][0]);
  nbytes += (nlatsites + nlatghosts) * 3 * sizeof(latsites[0][0]);
  nbytes += (nlatsites + nlatghosts) * sizeof(site_tag[0]);
  nbytes += nlatbins[0] * nlatbins[1] * nlatbins[2] * nlatbins[3] * sizeof(latbins[0][0][0][0]);
  nbytes += (nlatsites + nlatghosts) * (nlist.empty() ? 0 : nlist[0].size()) * sizeof(tagint);
  nbytes += nlatsites * sizeof(tagint);        // clusterID
  nbytes += 2 * noccupied * sizeof(int);       // cluster_size, cluster_nsites
  nbytes += 3 * noccupied * sizeof(double);    // cluster_center
  nbytes += noccupied * sizeof(tagint);        // occupied_cluster_ID;
  return nbytes;
}
