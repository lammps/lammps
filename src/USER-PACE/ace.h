//
// Created by lysogy36 on 31.01.20.
//

#ifndef ACE_H
#define ACE_H

#include "ace_types.h"
#include "atomic_environment.h"
#include "basis.h"
#include "complex.h"
#include "timing.h"

#include "USER-PACE/include/multiarray/multiarray.h"

class ACE {
    // f_ji vector
    DOUBLE_TYPE f_ji[3];

    //cache for xi:fixed, A(xj,n)=A2DLM(l,m)
    Array2D<DOUBLE_TYPE> A_rank1 = Array2D<DOUBLE_TYPE>("A_rank1");
    Array4DLM <Complex> A_dense = Array4DLM<Complex>("A_dense");

    Array2D<DOUBLE_TYPE> weights_rank1 = Array2D<DOUBLE_TYPE>("weights_rank1");
    Array4DLM <Complex> weights_dense = Array4DLM<Complex>("weights_dense");

    //cache for grads: grad_phi(jj,n)=A2DLM(l,m)
    //(neigh_jj)(n=0..nr-1)
    Array2D<DOUBLE_TYPE> DG_dense_cache = Array2D<DOUBLE_TYPE>("DG_dense_cache");
    // (neigh_jj)(n=0..nr-1,l)
    Array3D<DOUBLE_TYPE> R_dense_cache = Array3D<DOUBLE_TYPE>("R_dense_cache");
    Array3D<DOUBLE_TYPE> DR_dense_cache = Array3D<DOUBLE_TYPE>("DR_dense_cache");
    // (neigh_jj)(l,m)
    Array3DLM <Complex> Y_dense_cache = Array3DLM<Complex>("Y_dense_cache");
    Array3DLM <Dycomponent> DY_dense_cache = Array3DLM<Dycomponent>("dY_dense_cache");

    //densities rho(ndensity), index  = 0 .. ndensity-1
    Array1D<DOUBLE_TYPE> rhos = Array1D<DOUBLE_TYPE>("rhos");
    //derivatives wrt. densities  index = 0 .. ndensity-1
    Array1D<DOUBLE_TYPE> dF_dRho = Array1D<DOUBLE_TYPE>("dF_dRho");

    Array1D<Complex> dB_array = Array1D<Complex>("dB_array");

    //pointer to the ACEBasisSet object
    ACEBasisSet *basis_set;

    void init();

public:
    // set of timers for code profiling
    DEFINE_TIMER(loop_over_neighbour);
    DEFINE_TIMER(forces_calc_loop);
    DEFINE_TIMER(forces_calc_neighbour);
    DEFINE_TIMER(phi_calc);
    DEFINE_TIMER(phi_recalc);
    DEFINE_TIMER(energy_calc);
    DEFINE_TIMER(bond_calc);
    DEFINE_TIMER(A_calc);
    DEFINE_TIMER(per_atom_calc);
    DEFINE_TIMER(basis_func_calc);
    DEFINE_TIMER(total_time_calc);
    //TODO: integrate with lammps atoms mapping
    int map_lammps_at_type_to_element[6] = {0, 1, 2, 3, 4, 5}; // mapping from atom types to elements

    //total energy of AtomicEnvironment
    DOUBLE_TYPE total_energy = 0;

    //temporary array for the pair forces between current atom_i and its neighbours atom_k
    //pair_forces_ji(k,3)
    //k = 0..num_of_neighbours(atom_i)-1
    Array2D<DOUBLE_TYPE> pair_forces_ji = Array2D<DOUBLE_TYPE>("pair_forces_ji");

    //total forces array
    //forces(i,3), i = 0..num_of_atoms-1
    Array2D<DOUBLE_TYPE> forces = Array2D<DOUBLE_TYPE>("forces");

    //Per-atom energies
    //energies(i), i = 0..num_of_atoms-1
    Array1D<DOUBLE_TYPE> energies = Array1D<DOUBLE_TYPE>("energies");


    ACE() = default;


    //set the basis function to the ACE evaluator
    void set_basis(ACEBasisSet &bas);

    //compute the energies and forces for each atoms in atomic_environment
    //results are stored in forces and energies arrays
    void compute(AtomicEnvironment &atomic_environment, bool verbose = false);

    //compute the energy and forces for atom_i
    //x - atomic positions [atom_ind][3]
    //type - atomic types [atom_ind]
    //jnum - number of neighbours of atom_i
    //jlist - list of neighbour indices. Indices are for arrays a and type
    void compute_atom(int i, DOUBLE_TYPE **x, const SPECIES_TYPE *type, int jnum, const int *jlist);

};

#endif //ACE_ACE_H

