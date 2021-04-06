/*
 * Performant implementation of atomic cluster expansion and interface to LAMMPS
 *
 * Copyright 2021  (c) Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 * Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 * Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1
 *
 * ^1: Ruhr-University Bochum, Bochum, Germany
 * ^2: University of Cambridge, Cambridge, United Kingdom
 * ^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
 * ^4: University of British Columbia, Vancouver, BC, Canada
 *
 *
 * See the LICENSE file.
 * This FILENAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


// Created by Christoph Ortner on 20.12.2020

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


typedef pair<vector<int>, vector<int> > TPARTITION;
typedef list<TPARTITION> TPARTITIONS;

typedef map<vector<int>, int> TDAGMAP;

class ACEDAG {

    TPARTITIONS find_2partitions(vector<int> v);

    void insert_node(TDAGMAP &dagmap,
                     vector<int> node,
                     vector<DOUBLE_TYPE> c);

    // the following fields are used only for *construction*, not evaluation
    int dag_idx;     // current index of dag node 
    Array2D<int> nodes_pre; //TODO: YL: better to use vector<>
    Array2D<DOUBLE_TYPE> coeffs_pre; //TODO: YL: better to use vector<>
    Array1D<bool> haschild; //TODO: YL: better to use vector<>

    /* which heuristic to choose for DAG construction? 
     *   0 : the simple original heuristic
     *   1 : prioritize 2-correlation nodes and build the rest from those
     */
    int heuristic = 0;

public:

    ACEDAG() = default;

    void init(Array2D<int> Aspec,  Array2D<int> AAspec, 
              Array1D<int> orders, Array2D<DOUBLE_TYPE> coeffs, 
              int heuristic );

    Array1D<ACEComplex> AAbuf; 
    Array1D<ACEComplex> w; 
    
    Array2D<int> Aspec; 

    // nodes in the graph 
    Array2D<int> nodes;    
    Array2D<DOUBLE_TYPE> coeffs; 

    // total number of nodes in the dag
    int num_nodes;
    // number of interior nodes (with children)
    int num2_int;
    // number of leaf nodes (nc = no child)
    int num2_leaf;


    // number of 1-particle basis functions 
    // (these will be stored in the first num1 entries of AAbuf)
    int get_num1() { return Aspec.get_dim(0); };
    // total number of n-correlation basis functions n > 1.
    int get_num2() { return num_nodes - get_num1(); }; 
    int get_num2_int() { return num2_int; };   // with children
    int get_num2_leaf() { return num2_leaf; };     // without children

    // debugging tool
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
    void init(ACECTildeBasisSet *basis_set, int heuristic);

    /* convert the PACE to the ACE.jl format to prepare for DAG construction*/
    Array2D<int> jl_Aspec;
    Array2D<int> jl_AAspec;
    Array1D<int> jl_AAspec_flat;
    Array1D<int> jl_orders;
    Array2D<DOUBLE_TYPE> jl_coeffs; 
    void acejlformat();

    /* the main event : the computational graph */
    ACEDAG dag; 

    bool recursive = true;

public:


    ACERecursiveEvaluator() = default;

    explicit ACERecursiveEvaluator(ACECTildeBasisSet &bas, 
                                   bool recursive = true) {
        set_recursive(recursive);
        set_basis(bas);
    }

    /**
     * set the basis function to the ACE evaluator
     * @param bas
     */
    void set_basis(ACECTildeBasisSet &bas, int heuristic = 0);

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

    /******* public functions related to recursive evaluator ********/
    
    // print out the DAG for visual inspection
    void print_dag() {dag.print();}
    
    // print out the jl format for visual inspection
    // should be converted into a proper test 
    void test_acejlformat(); 

    void set_recursive(bool tf) { recursive = tf; }

    /********************************/

};


#endif //ACE_RECURSIVE_H