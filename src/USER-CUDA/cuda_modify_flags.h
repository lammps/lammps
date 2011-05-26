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

#ifndef CUDA_MODIFY_FLAGS_H
#define CUDA_MODIFY_FLAGS_H

#define INITIAL_INTEGRATE_CUDA  (1 << 16)
#define POST_INTEGRATE_CUDA     (1 << 17)
#define PRE_EXCHANGE_CUDA       (1 << 18)
#define PRE_NEIGHBOR_CUDA       (1 << 19)
#define PRE_FORCE_CUDA          (1 << 20)
#define POST_FORCE_CUDA         (1 << 21)
#define FINAL_INTEGRATE_CUDA    (1 << 22)
#define END_OF_STEP_CUDA        (1 << 23)
#define THERMO_ENERGY_CUDA      (1 << 24)
#define MIN_POST_FORCE_CUDA      (1 << 25)
// remember not to shift over 31 bits

#endif // CUDA_MODIFY_FLAGS_H
