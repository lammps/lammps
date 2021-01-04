//
// Created by Yury Lysogorskiy on 31.01.20.
//

#ifndef ACE_RECURSIVE_H
#define ACE_RECURSIVE_H

#include "ace_abstract_basis.h"
#include "ace_arraynd.h"
#include "ace_array2dlm.h"
#include "ace_c_basis.h"
#include "ace_complex.h"
#include "ace_timing.h"
#include "ace_types.h"
#include "ace_evaluator.h"

#include <list>
#include <utility>
#include <algorithm>
#include <map>
#include <vector>

using namespace std;

typedef list <pair<vector<int>, vector<int> >> TPARTITIONS;
typedef pair<vector<int>, vector<int> > TPARTITION;
typedef map<vector<int>, int> TDAGMAP;

class ACEDAG {


    TPARTITIONS find_2partitions(vector<int> v);

    void insert_node(TDAGMAP &dagmap,
                     vector<int> node,
                     vector<ACEComplex> c);

    int dag_idx;    // current index for dag construction


public:

    ACEDAG() = default;

    void init(Array2D<int> Aspec, Array2D<int> AAspec,
              Array1D<int> orders, Array2D<ACEComplex> coeffs);

    Array1D<ACEComplex> AAbuf;
    Array1D<ACEComplex> w;

    Array2D<int> nodes;    // TODO: split into dependent and independent nodes
    Array2D<int> Aspec;
    Array2D<ACEComplex> coeffs;


    int get_num1() { return Aspec.get_dim(0); };

    int get_num2() { return num_nodes - get_num1(); };


    int num_nodes;  // store number of nodes in dag 

    void print();

};


/**
 * Recursive Variant of the ACETildeEvaluator; should be 100% compatible
 */
class ACERecursiveEvaluator : public ACEEvaluator {

    /**
     * Weights \f$ \omega_{i \mu n 0 0} \f$ for rank = 1, see Eq.(10) from implementation notes,
     * 'i' is fixed for the current atom, shape: [nelements][nradbase]
     */
    Array2D<DOUBLE_TYPE> weights_rank1 = Array2D<DOUBLE_TYPE>("weights_rank1");

    /**
     * Weights \f$ \omega_{i \mu n l m} \f$ for rank > 1, see Eq.(10) from implementation notes,
     * 'i' is fixed for the current atom, shape: [nelements][nradbase][l=0..lmax, m]
     */
    Array4DLM<ACEComplex> weights = Array4DLM<ACEComplex>("weights");

    /**
     * cache for gradients of \f$ g(r)\f$: grad_phi(jj,n)=A2DLM(l,m)
     * shape:[max_jnum][nradbase]
     */
    Array2D<DOUBLE_TYPE> DG_cache = Array2D<DOUBLE_TYPE>("DG_cache");


    /**
     * cache for \f$ R_{nl}(r)\f$
     * shape:[max_jnum][nradbase][0..lmax]
     */
    Array3D<DOUBLE_TYPE> R_cache = Array3D<DOUBLE_TYPE>("R_cache");
    /**
     * cache for derivatives of \f$ R_{nl}(r)\f$
     * shape:[max_jnum][nradbase][0..lmax]
     */
    Array3D<DOUBLE_TYPE> DR_cache = Array3D<DOUBLE_TYPE>("DR_cache");
    /**
     * cache for \f$ Y_{lm}(\hat{r})\f$
     * shape:[max_jnum][0..lmax][m]
     */
    Array3DLM<ACEComplex> Y_cache = Array3DLM<ACEComplex>("Y_cache");
    /**
     * cache for \f$ \nabla Y_{lm}(\hat{r})\f$
     * shape:[max_jnum][0..lmax][m]
     */
    Array3DLM<ACEDYcomponent> DY_cache = Array3DLM<ACEDYcomponent>("dY_dense_cache");

    /**
     * cache for derivatives of hard-core repulsion
     * shape:[max_jnum]
     */
    Array1D<DOUBLE_TYPE> DCR_cache = Array1D<DOUBLE_TYPE>("DCR_cache");

    /**
    * Partial derivatives \f$ dB_{i \mu n l m t}^{(r)} \f$  with sequential numbering over [func_ind][ms_ind][r],
    * shape:[func_ms_r_ind]
    */
    Array1D<ACEComplex> dB_flatten = Array1D<ACEComplex>("dB_flatten");

    /**
     * pointer to the ACEBasisSet object
     */
    ACECTildeBasisSet *basis_set = nullptr;

    /**
     * Initialize internal arrays according to basis set sizes
     * @param basis_set
     */
    void init(ACECTildeBasisSet *basis_set);

    /* convert the PACE to the ACE.jl format to prepare for DAG construction*/
    Array2D<int> jl_Aspec;
    Array2D<int> jl_AAspec;
    Array1D<int> jl_orders;
    Array2D<ACEComplex> jl_coeffs;

    void acejlformat();

    /* the main event : the computational graph */
    ACEDAG dag;

public:

    /******* debugging codes...********/
    void print_dag() { dag.print(); }

    void test_acejlformat();

    /********************************/

    ACERecursiveEvaluator() = default;

    explicit ACERecursiveEvaluator(ACECTildeBasisSet &bas) {
        set_basis(bas);
    }

    /**
     * set the basis function to the ACE evaluator
     * @param bas
     */
    void set_basis(ACECTildeBasisSet &bas);

    /**
     * The key method to compute energy and forces for atom 'i'.
     * Method will update the  "e_atom" variable and "neighbours_forces(jj, alpha)" array
     *
     * @param i atom index
     * @param x atomic positions array of the real and ghost atoms, shape: [atom_ind][3]
     * @param type  atomic types array of the real and ghost atoms, shape: [atom_ind]
     * @param jnum  number of neighbours of atom_i
     * @param jlist array of neighbour indices, shape: [jnum]
     */
    void compute_atom(int i, DOUBLE_TYPE **x, const SPECIES_TYPE *type, const int jnum, const int *jlist) override;

    /**
     * Resize all caches over neighbours atoms
     * @param max_jnum  maximum number of neighbours
     */
    void resize_neighbours_cache(int max_jnum) override;
};


#endif //ACE_RECURSIVE_H