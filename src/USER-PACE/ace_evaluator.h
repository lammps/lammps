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

class ACEEvaluator {
protected:
    //cache for xi:fixed, A(xj,n)=A2DLM(l,m)
    Array2D<DOUBLE_TYPE> A_rank1 = Array2D<DOUBLE_TYPE>("A_rank1");
    Array4DLM<ACEComplex> A = Array4DLM<ACEComplex>("A");

    //densities rho(ndensity), index  = 0 .. ndensity-1
    Array1D<DOUBLE_TYPE> rhos = Array1D<DOUBLE_TYPE>("rhos");
    //derivatives wrt. densities  index = 0 .. ndensity-1
    Array1D<DOUBLE_TYPE> dF_drho = Array1D<DOUBLE_TYPE>("dF_drho");

    void init(ACEAbstractBasisSet *basis_set);

public:
    // set of timers for code profiling
    ACETimer loop_over_neighbour_timer;
    ACETimer per_atom_calc_timer;


    ACETimer forces_calc_loop_timer;
    ACETimer forces_calc_neighbour_timer;

    ACETimer phi_calc_timer;
    ACETimer phi_recalc_timer;
    ACETimer energy_calc_timer;
    ACETimer bond_calc_timer;
    ACETimer A_calc_timer;
    ACETimer basis_func_calc_timer;
    ACETimer total_time_calc_timer;

    void init_timers();

    //TODO: integrate with lammps atoms mapping
    int map_lammps_at_type_to_element[6] = {0, 0, 1, 2, 3, 4}; // mapping from atom types to elements


    DOUBLE_TYPE e_atom = 0;
    //temporary array for the pair forces between current atom_i and its neighbours atom_k
    //neighbours_forces(k,3)
    //k = 0..num_of_neighbours(atom_i)-1
    Array2D<DOUBLE_TYPE> neighbours_forces = Array2D<DOUBLE_TYPE>("neighbours_forces");

    ACEEvaluator() = default;

    virtual ~ACEEvaluator() = default;

    //compute the energy and forces for atom_i
    //x - atomic positions [atom_ind][3]
    //type - atomic types [atom_ind]
    //jnum - number of neighbours of atom_i
    //jlist - list of neighbour indices. Indices are for arrays a and type
    //this will also update the energies(i) and neighbours_forces(jj, alpha) arrays
    virtual void compute_atom(int i, DOUBLE_TYPE **x, const SPECIES_TYPE *type, int jnum, const int *jlist) = 0;

    virtual void resize_neighbours_cache(int max_jnum) = 0;
};

class ACECTildeEvaluator : public ACEEvaluator {

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


    Array1D<ACEComplex> dB_flatten = Array1D<ACEComplex>("dB_flatten");

    //pointer to the ACEBasisSet object
    ACECTildeBasisSet *basis_set;

    void init(ACECTildeBasisSet *basis_set);

public:

    ACECTildeEvaluator() = default;

    //set the basis function to the ACE evaluator
    void set_basis(ACECTildeBasisSet &bas);

    //compute the energy and forces for atom_i
    //x - atomic positions [atom_ind][3]
    //type - atomic types [atom_ind]
    //jnum - number of neighbours of atom_i
    //jlist - list of neighbour indices. Indices are for arrays a and type
    //this will also update the energies(i) and neighbours_forces(jj, alpha) arrays
    void compute_atom(int i, DOUBLE_TYPE **x, const SPECIES_TYPE *type, int jnum, const int *jlist) override;

    void resize_neighbours_cache(int max_jnum) override;
};


#endif //ACE_ACE_H