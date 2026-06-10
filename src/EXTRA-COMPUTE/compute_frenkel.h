/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(frenkel, ComputeFrenkel)

#else

#ifndef LMP_COMPUTE_FRENKEL_H
#define LMP_COMPUTE_FRENKEL_H

#include "compute.h"
#include "region.h"

#include <list>

namespace LAMMPS_NS {

class ComputeFrenkel : public Compute {
 public:
  ComputeFrenkel(class LAMMPS *, int, char **);
  ~ComputeFrenkel() override;
  int modify_param(int, char **) override;

  void init() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  void compute_vector() override;
  void compute_array() override;
  void compute_peratom() override;
  void compute_local() override;
  int compute_image(int *&, double **&) override;
  double memory_usage() override;

 private:
  // graphics-object output for "dump image ... compute <ID> ...":
  // one sphere per defect cluster (vacancy or interstitial)
  int *image_objvec;
  double **image_objarray;
  int image_nmax;

  Region *region;
  std::string sitefile;
  int ifgroup, fgroupbit;
  bool rescale;

  double *mindist;
  double *site_mindist;
  int *noccupants;
  tagint **occupant_tag;
  int nnormal;
  double **normal;
  double cut_vac, cut_int, cutoff, binwidth;
  int nlatsites, nlatghosts;
  double **latsites;
  double **latsites0;
  tagint *site_tag;
  tagint first_local_tag;
  int nlatbins[4];
  int ****latbins;
  std::vector<std::list<tagint>> nlist;    // per-site neighbor lists
  tagint *clusterID;                       // Per-site vector
  std::vector<int> cluster_size;           // Negative => vacancy; positive => interstitial
  std::vector<int> cluster_nsites;         // Number of sites involved in cluster
  double **cluster_center;                 // Geometric center of cluster in x,y,z
  int noccupied;
  tagint *occupied_cluster_ID;    // Per-cluster vector, length noccupied
  double old_boxlo[3], old_boxhi[3];

  bigint invoked_find_defects;
  bigint invoked_find_clusters;
  bigint invoked_construct_WS_cell;

  void create_lattice_sites();
  void put_sites_in_bins();
  int site_tag2index(tagint);
  void exchange_lattice_ghosts();
  void construct_WS_cell();
  bool inside_WS_cell(int, int);
  void find_defects();
  void find_clusters();
  void find_occupied_clusters();
  int clusterID2occupied_index(int);
  void find_closest_bin(double *, int &, int &, int &);
  void bin_pbc(int &, int &, int &);
  bool tag_is_already_in_occupancy_list(tagint, int);
  int next_free_occupant_tag_index(int, int);
  template <typename TYPE> void reallocate_array(TYPE **&, int, int, int, int);
  template <typename TYPE> void reallocate_array(TYPE ***&, int, int, int, int, int, int);
  void construct_site_nlists();
  void rescale_lattice_sites();

  int process_neighbor(int, int, int);

  static int compareIDs(const void *, const void *);
};    // end class ComputeFrenkel
}    // namespace LAMMPS_NS

#endif
#endif
