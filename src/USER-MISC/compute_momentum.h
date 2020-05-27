#ifdef COMPUTE_CLASS

ComputeStyle(momentum,ComputeMomentum)

#else

#ifndef LMP_COMPUTE_MOMENTUM_H
#define LMP_COMPUTE_MOMENTUM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeMomentum : public Compute {
 public:
  ComputeMomentum(class LAMMPS *, int, char **);
  virtual ~ComputeMomentum();

  virtual void init();
  virtual void compute_vector();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
