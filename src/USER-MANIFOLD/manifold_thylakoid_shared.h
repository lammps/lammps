#ifndef MANIFOLD_THYLAKOID_SHARED_H
#define MANIFOLD_THYLAKOID_SHARED_H


#include "lmptype.h"
#include <vector>

namespace LAMMPS_NS {

namespace user_manifold {


  // The thylakoid is composed of many parts
  struct thyla_part {
    enum thyla_types {
      THYLA_TYPE_PLANE,
      THYLA_TYPE_SPHERE,
      THYLA_TYPE_CYL,
      THYLA_TYPE_CYL_TO_PLANE
    };

    thyla_part( int type, double *args, double xlo, double ylo, double zlo,
                double xhi, double yhi, double zhi );
    thyla_part() : type(-1), x0(-1337), y0(-1337), z0(-1337){}
    ~thyla_part();

    double g( const double *x );
    void   n( const double *x, double *n );

    int type;
    double params[7];
    double tol;
    int maxit;

    int err_flag;
    tagint a_id;

    double xlo, xhi, ylo, yhi, zlo, zhi;
    double x0, y0, z0;

  }; // struct thyla_part



  struct thyla_part_geom {
    thyla_part_geom() : pt(3), lo(3), hi(3){}
    std::vector<double> pt, lo, hi;

    // Function for mirroring thyla_geoms:
    enum DIRS { DIR_X, DIR_Y, DIR_Z };
    static void mirror( unsigned int axis, thyla_part_geom *m,
                        const thyla_part_geom *o );

  }; // struct thyla_part_geom


} // namespace user_manifold


} // namespace LAMMPS_NS


#endif // MANIFOLD_THYLAKOID_SHARED_H
