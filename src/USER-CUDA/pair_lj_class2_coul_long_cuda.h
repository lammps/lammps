/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator 

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov 

   See the README file in the top-level LAMMPS directory. 

   ----------------------------------------------------------------------- 

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/ 

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany 

   See the README file in the USER-CUDA directory. 

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/class2/coul/long/cuda,PairLJClass2CoulLongCuda)

#else

#ifndef LMP_PAIR_LJ_CLASS2_COUL_LONG_CUDA_H
#define LMP_PAIR_LJ_CLASS2_COUL_LONG_CUDA_H

#include "pair_lj_class2_coul_long.h"

namespace LAMMPS_NS {

class PairLJClass2CoulLongCuda : public PairLJClass2CoulLong
{
	public:
		PairLJClass2CoulLongCuda(class LAMMPS *);
		void compute(int, int);
		void settings(int, char **);
		void coeff(int, char **);
		void init_list(int, class NeighList *);
		void init_style();
		void ev_setup(int eflag, int vflag);
	protected:
		class Cuda *cuda;
		void allocate();
		bool allocated2;
		class CudaNeighList* cuda_neigh_list;
};

}

#endif
#endif
