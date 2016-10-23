#ifdef FIX_CLASS

FixStyle(wall/asphere,FixWallAsphere)

#else

#ifndef LMP_FIX_WALL_ASPHERE_H
#define LMP_FIX_WALL_ASPHERE_H

#include "fix_wall.h"

namespace LAMMPS_NS {

class FixWallAsphere : public FixWall {
 public:
  FixWallAsphere(class LAMMPS *, int, char **);
  void precompute(int);
  void wall_particle(int, int, double);

 private:
  double coeff1[6],coeff2[6],coeff3[6],coeff4[6],coeff5[6],coeff6[6];
  class AtomVecEllipsoid *avec;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Particle on or inside fix wall surface

Particles must be "exterior" to the wall in order for energy/force to
be calculated.

*/
