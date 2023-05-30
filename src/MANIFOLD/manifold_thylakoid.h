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

#ifndef LMP_MANIFOLD_THYLAKOID_H
#define LMP_MANIFOLD_THYLAKOID_H

#include "manifold.h"
#include <vector>

namespace LAMMPS_NS {

namespace user_manifold {
  struct thyla_part;

  class manifold_thylakoid : public manifold {
   public:
    enum { NPARAMS = 3 };
    manifold_thylakoid(LAMMPS *lmp, int, char **);
    ~manifold_thylakoid() override;

    double g(const double *x) override;
    void n(const double *x, double *n) override;

    static const char *type() { return "thylakoid"; }
    const char *id() override { return type(); }
    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }

    void post_param_init() override;

   private:
    void init_domains();

    thyla_part *get_thyla_part(const double *x, int *err_flag, std::size_t *idx = nullptr);
    int is_in_domain(thyla_part *p, const double *x);
    void check_overlap();
    std::vector<thyla_part *> parts;

    thyla_part *make_plane_part(double a, double b, double c, const std::vector<double> &pt);
    thyla_part *make_cyl_part(double a, double b, double c, const std::vector<double> &pt,
                              double R);
    thyla_part *make_sphere_part(const std::vector<double> &pt, double R);
    thyla_part *make_cyl_to_plane_part(double X0, double R0, double R, double s,
                                       const std::vector<double> &pt);

    void set_domain(thyla_part *p, const std::vector<double> &lo, const std::vector<double> &hi);

    // Coefficients for the thylakoid model. At the moment it is just
    // a cylinder, we slowly expand it.
    double pad;    // Padding to make sure periodic images are mapped back properly.
    double LB, lT, lB, wB, LT;

    // Domain size:
    double x0, y0, z0;
    double x1, y1, z1;
    double Lx, Ly, Lz;
  };

}    // namespace user_manifold

}    // namespace LAMMPS_NS

#endif    // LMP_MANIFOLD_THYLAKOID_H
