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

#ifndef LMP_DOMAIN_H
#define LMP_DOMAIN_H

#include "pointers.h"

#include <cmath>
#include <map>
#include <unordered_set>

namespace LAMMPS_NS {
class Region;

class Domain : protected Pointers {
 public:
  int box_exist;                          // 0 = not yet created, 1 = exists
  int dimension;                          // 2 = 2d, 3 = 3d
  int nonperiodic;                        // 0 = periodic in all 3 dims
                                          // 1 = periodic or fixed in all 6
                                          // 2 = shrink-wrap in any of 6
  int xperiodic, yperiodic, zperiodic;    // 0 = non-periodic, 1 = periodic
  int periodicity[3];                     // xyz periodicity as array

  int boundary[3][2];    // settings for 6 boundaries
                         // 0 = periodic
                         // 1 = fixed non-periodic
                         // 2 = shrink-wrap non-periodic
                         // 3 = shrink-wrap non-per w/ min

  int triclinic;          // 0 = orthog box, 1 = triclinic (restricted or general)
  int triclinic_general;  // 1 if general <-> restricted tri mapping is stored, 0 if not

  // orthogonal box

  double xprd, yprd, zprd;                   // global box dimensions
  double xprd_half, yprd_half, zprd_half;    // half dimensions
  double prd[3];                             // array form of dimensions
  double prd_half[3];                        // array form of half dimensions

  // restricted triclinic box
  // xyz prd,xyz prd_half and prd,prd_half = same as if not tilted

  double prd_lamda[3];         // lamda box = (1,1,1)
  double prd_half_lamda[3];    // lamda half box = (0.5,0.5,0.5)

  // orthogonal box global bounds

  double boxlo[3], boxhi[3];

  // restricted triclinic box
  // boxlo/hi = same as if not tilted

  double boxlo_lamda[3], boxhi_lamda[3];    // lamda box = (0,1)
  double boxlo_bound[3], boxhi_bound[3];    // bounding box of tilted domain
  double corners[8][3];                     // 8 corner points

  // orthogonal box & restricted triclinic box

  double minxlo, minxhi;    // minimum size of global box
  double minylo, minyhi;    // when shrink-wrapping
  double minzlo, minzhi;    // tri only possible for non-skew dims

  // orthogonal box

  double sublo[3], subhi[3];    // sub-box bounds on this proc

  // restricted triclinic box
  // sublo/hi = undefined

  double sublo_lamda[3], subhi_lamda[3];    // bounds of subbox in lamda

  // restricted triclinic box

  double xy, xz, yz;                // 3 tilt factors
  double h[6], h_inv[6];            // shape matrix in Voigt ordering
                                    // Voigt = xx,yy,zz,yz,xz,xy
  double h_rate[6], h_ratelo[3];    // rate of box size/shape change

  // general triclinic box
  // boxlo = lower left corner

  double avec[3], bvec[3], cvec[3];  // ABC edge vectors of general triclinic box
  double rotate_g2r[3][3];           // rotation matrix from general --> restricted tri
  double rotate_r2g[3][3];           // rotation matrix from restricted --> general tri

  // box flags

  int box_change;           // 1 if any of next 3 flags are set, else 0
  int box_change_size;      // 1 if box size changes, 0 if not
  int box_change_shape;     // 1 if box shape changes, 0 if not
  int box_change_domain;    // 1 if proc sub-domains change, 0 if not

  int deform_flag;        // 1 if fix deform exist, else 0
  int deform_vremap;      // 1 if fix deform remaps v, else 0
  int deform_groupbit;    // atom group to perform v remap for

  class Lattice *lattice;    // user-defined lattice

  int copymode;
  enum { NO_REMAP, X_REMAP, V_REMAP };

  typedef Region *(*RegionCreator)(LAMMPS *, int, char **);
  typedef std::map<std::string, RegionCreator> RegionCreatorMap;
  RegionCreatorMap *region_map;

  Domain(class LAMMPS *);
  ~Domain() override;
  virtual void init();
  void set_initial_box(int expandflag = 1);
  virtual void set_global_box();
  virtual void set_lamda_box();
  virtual void set_local_box();
  virtual void reset_box();
  virtual void pbc();
  void image_check();
  void box_too_small_check();
  void subbox_too_small_check(double);
  void minimum_image(double &, double &, double &) const;
  void minimum_image(double *delta) const { minimum_image(delta[0], delta[1], delta[2]); }
  void minimum_image_big(double &, double &, double &) const;
  void minimum_image_big(double *delta) const { minimum_image_big(delta[0], delta[1], delta[2]); }
  int closest_image(int, int);
  int closest_image(const double *const, int);
  void closest_image(const double *const, const double *const, double *const);
  void remap(double *, imageint &);
  void remap(double *);
  void remap_near(double *, double *);
  void unmap_inv(double *x, imageint);
  void unmap(double *, imageint);
  void unmap(const double *, imageint, double *);
  void image_flip(int, int, int);
  int ownatom(int, double *, imageint *, int);

  void define_general_triclinic(double *, double *, double *, double *);
  void general_to_restricted_rotation(double *, double *, double *,
                                      double [3][3],
                                      double *, double *, double *);
  void general_to_restricted_coords(double *);
  void restricted_to_general_coords(double *);
  void restricted_to_general_coords(double *, double *);
  void general_to_restricted_vector(double *);
  void restricted_to_general_vector(double *);
  void restricted_to_general_vector(double *, double *x);

  void set_lattice(int, char **);
  void add_region(int, char **);
  void delete_region(Region *);
  void delete_region(const std::string &);
  Region *get_region_by_id(const std::string &) const;
  const std::vector<Region *> get_region_by_style(const std::string &) const;
  const std::vector<Region *> get_region_list();
  void set_boundary(int, char **, int);
  void print_box(const std::string &);
  void boundary_string(char *);

  virtual void lamda2x(int);
  virtual void x2lamda(int);
  virtual void lamda2x(double *, double *);
  virtual void x2lamda(double *, double *);
  int inside(double *);
  int inside_nonperiodic(double *);
  void x2lamda(double *, double *, double *, double *);
  void bbox(double *, double *, double *, double *);
  void box_corners();
  void subbox_corners();
  void lamda_box_corners(double *, double *);

  // minimum image convention check
  // return 1 if any distance > 1/2 of box size
  // indicates a special neighbor is actually not in a bond,
  //   but is a far-away image that should be treated as an unbonded neighbor
  // inline since called from neighbor build inner loop

  inline int minimum_image_check(double dx, double dy, double dz)
  {
    if (xperiodic && fabs(dx) > xprd_half) return 1;
    if (yperiodic && fabs(dy) > yprd_half) return 1;
    if (zperiodic && fabs(dz) > zprd_half) return 1;
    return 0;
  }

 protected:
  double small[3];    // fractions of box lengths
  std::unordered_set<Region *> regions;
};

}    // namespace LAMMPS_NS

#endif
