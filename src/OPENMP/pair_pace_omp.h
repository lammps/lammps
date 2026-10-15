#ifdef PAIR_CLASS
// clang-format off
PairStyle(pace/omp, PairPACEOMP);
// clang-format on
#else

#ifndef LMP_PAIR_PACE_OMP_H
#define LMP_PAIR_PACE_OMP_H

#include "pair_pace.h"
#include "ace-evaluator/ace_recursive.h"
#include "thr_omp.h"

#include <vector>

namespace LAMMPS_NS {

// inherit from PairPACE and ThrOMP
class PairPACEOMP : public PairPACE, public ThrOMP {
public:
PairPACEOMP(class LAMMPS *);
~PairPACEOMP() override;

void compute(int, int) override;
void coeff(int, char **) override;
double memory_usage() override;

private:
// ACEImpl contains two separate objects that are not thread-safe by default, allocate copies on each thread
// One ACE evaluator per thread
std::vector<ACERecursiveEvaluator *> ace_thr;
// One basis per thread
std::vector<ACECTildeBasisSet *> basis_thr;
// Helper to create one ACERecursiveEvaluator per thread
void build_thread_evaluators();
};

} // namespace LAMMPS_NS

#endif
#endif
