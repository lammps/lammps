#ifdef COMPUTE_CLASS

ComputeStyle(deeptensor/atom,ComputeDeeptensorAtom)

#else

#ifndef LMP_COMPUTE_DEEPTENSOR_ATOM_H
#define LMP_COMPUTE_DEEPTENSOR_ATOM_H

#include "compute.h"
#include "pair_deepmd.h"
#ifdef LMPPLUGIN
#include "DeepTensor.h"
#else
#include "deepmd/DeepTensor.h"
#endif

namespace LAMMPS_NS {

class ComputeDeeptensorAtom : public Compute {
 public:
  ComputeDeeptensorAtom(class LAMMPS *, int, char **);
  ~ComputeDeeptensorAtom();
  void init();
  void compute_peratom();
  double memory_usage();
  void init_list(int, class NeighList *);

 private:
  int nmax;
  double **tensor;
  PairDeepMD dp;
  class NeighList *list;
  deepmd::DeepTensor dt;
  std::vector<int > sel_types;
};

}

#endif
#endif

