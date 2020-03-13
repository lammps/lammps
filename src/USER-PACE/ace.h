//
// Created by Yury Lysogorskiy on 31.01.20.
//

#ifndef ACE_H
#define ACE_H


#include "ace_arraynd.h"
#include "ace_array2dlm.h"
#include "ace_c_basis.h"
#include "ace_complex.h"
#include "ace_timing.h"
#include "ace_types.h"

class ACE {

    // f_ji vector
    DOUBLE_TYPE f_ji[3];

    //cache for xi:fixed, A(xj,n)=A2DLM(l,m)
    Array2D<DOUBLE_TYPE> A_rank1 = Array2D<DOUBLE_TYPE>("A_rank1");
    Array4DLM<ACEComplex> A = Array4DLM<ACEComplex>("A");

    Array2D<DOUBLE_TYPE> weights_rank1 = Array2D<DOUBLE_TYPE>("weights_rank1");
    Array4DLM<ACEComplex> weights = Array4DLM<ACEComplex>("weights");

    //cache for grads: grad_phi(jj,n)=A2DLM(l,m)
    //(neigh_jj)(n=0..nr-1)
    Array2D<DOUBLE_TYPE> DG_cache = Array2D<DOUBLE_TYPE>("DG_cache");
    // (neigh_jj)(n=0..nr-1,l)
    Array3D<DOUBLE_TYPE> R_cache = Array3D<DOUBLE_TYPE>("R_cache");
    Array3D<DOUBLE_TYPE> DR_cache = Array3D<DOUBLE_TYPE>("DR_cache");
    // (neigh_jj)(l,m)
    Array3DLM<ACEComplex> Y_cache = Array3DLM<ACEComplex>("Y_cache");
    Array3DLM<Dycomponent> DY_cache = Array3DLM<Dycomponent>("dY_dense_cache");

    //hard-core repulsion
    //(neigh_jj)
    Array1D<DOUBLE_TYPE> DCR_cache = Array1D<DOUBLE_TYPE>("DCR_cache");

    //densities rho(ndensity), index  = 0 .. ndensity-1
    Array1D<DOUBLE_TYPE> rhos = Array1D<DOUBLE_TYPE>("rhos");
    //derivatives wrt. densities  index = 0 .. ndensity-1
    Array1D<DOUBLE_TYPE> dF_drho = Array1D<DOUBLE_TYPE>("dF_drho");

    Array1D<ACEComplex> dB_flatten = Array1D<ACEComplex>("dB_flatten");

    //pointer to the ACEBasisSet object
    ACEBasisSet *basis_set;

    void init();

public:
    // set of timers for code profiling
    ACE_DEFINE_TIMER(loop_over_neighbour);
    ACE_DEFINE_TIMER(forces_calc_loop);
    ACE_DEFINE_TIMER(forces_calc_neighbour);
    ACE_DEFINE_TIMER(phi_calc);
    ACE_DEFINE_TIMER(phi_recalc);
    ACE_DEFINE_TIMER(energy_calc);
    ACE_DEFINE_TIMER(bond_calc);
    ACE_DEFINE_TIMER(A_calc);
    ACE_DEFINE_TIMER(per_atom_calc);
    ACE_DEFINE_TIMER(basis_func_calc);
    ACE_DEFINE_TIMER(total_time_calc);

    //TODO: integrate with lammps atoms mapping
    int map_lammps_at_type_to_element[6] = {0, 0, 1, 2, 3, 4}; // mapping from atom types to elements

    //total energy of ACEAtomicEnvironment
    DOUBLE_TYPE energy = 0;

    //temporary array for the pair forces between current atom_i and its neighbours atom_k
    //neighbours_forces(k,3)
    //k = 0..num_of_neighbours(atom_i)-1
    Array2D<DOUBLE_TYPE> neighbours_forces = Array2D<DOUBLE_TYPE>("neighbours_forces");

    //total forces array
    //forces(i,3), i = 0..num_of_atoms-1
    Array2D<DOUBLE_TYPE> forces = Array2D<DOUBLE_TYPE>("forces");

    //Per-atom energies
    //energies(i), i = 0..num_of_atoms-1
    Array1D<DOUBLE_TYPE> energies = Array1D<DOUBLE_TYPE>("energies");

    ACE() = default;

    //set the basis function to the ACE evaluator
    void set_basis(ACEBasisSet &bas);

    //compute the energy and forces for atom_i
    //x - atomic positions [atom_ind][3]
    //type - atomic types [atom_ind]
    //jnum - number of neighbours of atom_i
    //jlist - list of neighbour indices. Indices are for arrays a and type
    void compute_atom(int i, DOUBLE_TYPE **x, const SPECIES_TYPE *type, int jnum, const int *jlist);

    void resize_neighbours_cache(int max_jnum);
};


#endif //ACE_ACE_H